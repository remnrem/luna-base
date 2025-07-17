
//    --------------------------------------------------------------------
//
//    This file is part of Luna.
//
//    LUNA is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    Luna is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with Luna. If not, see <http://www.gnu.org/licenses/>.
//
//    Please see LICENSE.txt for more details.
//
//    --------------------------------------------------------------------

#include "edf/sedf.h"
#include "edf/edf.h"
#include "edf/slice.h"

extern logger_t logger;

sedf_t::sedf_t( edf_t & edf , param_t & param )
{
  
  // which signals to extract? (only data channels)
  const bool no_annotations = true;
  signal_list_t signals = edf.header.signal_list( param.requires( "sig" ) , no_annotations );
  
  
  // get file name for new SEDF
  std::string filename = edf.filename;

  // strip existing extension
  if ( Helper::file_extension( filename, "edf" ) )
    filename = filename.substr(0 , filename.size() - 4 );
  else if ( Helper::file_extension( filename, "edfz" ) )       
    filename = filename.substr(0 , filename.size() - 5 );

  // for now, use .sedf extension for these summary EDFs
  filename += ".sedf";

  // optionally, allow directory change
  if ( param.has( "sedf-dir" ) )
    {
      const std::string outdir = Helper::expand( param.value("sedf-dir") );
      
      if ( outdir[ outdir.size() - 1 ] != globals::folder_delimiter ) 
	Helper::halt("sedf-dir value must end in '" 
		     + std::string(1,globals::folder_delimiter) 
		     + "' to specify a folder" );
      int p=filename.size()-1;
      int v = 0;
      for (int j=p;j>=0;j--)
	{
	  if ( filename[j] == globals::folder_delimiter ) { v=j+1; break; }
	}      
      filename = outdir + filename.substr( v );                  
      // create folder if it does not exist       
      std::string syscmd = globals::mkdir_command + " " + param.value( "sedf-dir" );      
      int retval = system( syscmd.c_str() );      
    }





  //
  // Get statistics: Hjorth or means
  //

  std::map<std::string,std::vector<double> > stats;

  // set epochs 
  const int ne = edf.timeline.first_epoch();
  const int ns = signals.size();
  
  for (int s=0; s<ns; s++)
    {

      edf.timeline.first_epoch();

      channel_type_t ctype = globals::map_channel( signals.label(s) );
      
      if ( ctype == IGNORE_SIGNAL ) continue;
      
      // only for these channel types should we write Hjorth parameters, otherwise just means
      bool hjorth = ctype == EEG || ctype == REF || ctype == EMG || ctype == EOG || ctype == ECG;
      
      const std::string ch  = signals.label( s );

      const std::string ch1 = ( hjorth ? "H1_" : "M_" ) + ch ; // M mean
      const std::string ch2 = ( hjorth ? "H2_" : "L_" ) + ch ; // L min
      const std::string ch3 = ( hjorth ? "H3_" : "U_" ) + ch ; // U max

      logger << "  extracting summary statistics for " << ch << "\n";

      // iterate over epochs
      
      while ( 1 ) 
	{
	  
	  // get epoch
	  int epoch = edf.timeline.next_epoch();      	  
	  if ( epoch == -1 ) break;

	  // get data
          interval_t interval = edf.timeline.epoch( epoch );
          slice_t slice( edf , signals(s) , interval );
	  std::vector<double> * d = slice.nonconst_pdata();

	  // get Hjorth parameters
	  if ( hjorth ) 
	    {
	      double activity = 0 , mobility = 0 , complexity = 0;
	      MiscMath::hjorth( d , &activity , &mobility , &complexity , ! globals::legacy_hjorth );
	      stats[ ch1 ].push_back( activity );
	      stats[ ch2 ].push_back( mobility );
	      stats[ ch3 ].push_back( complexity );
	    }
	  else // otherwise mean/min/max
	    {
	      double min, max;
	      MiscMath::minmax( *d , &min, &max );
	      stats[ ch1 ].push_back( MiscMath::mean( *d ) );
	      stats[ ch2 ].push_back( min );
	      stats[ ch3 ].push_back( max );
	    }
	  
	  // next epoch
	}
  
      // next signal
    }
  
  //
  // Create new SEDF object 
  //

  logger << "  writing SEDF to " << filename << "\n";

  const int ns_summ = stats.size();

  // for a 30 second epoch, set SEDF record duration to 30
  // set numbers of samples per record to 1
  
  const int recdur_summ = edf.timeline.epoch_length();
  const int nr_summ = ne;

  //
  // Create the new SEDF object 
  //

  annotation_set_t annotations;

  edf_t sedf( &annotations );
  
  //
  // Set header
  //

  sedf.header.version = edf.header.version;
  sedf.header.patient_id = edf.header.patient_id;
  sedf.header.recording_info = edf.header.recording_info;
  sedf.header.startdate = edf.header.startdate;
  sedf.header.starttime = edf.header.starttime;
  sedf.header.nbytes_header = 256 + ns_summ * 256;
  sedf.header.ns = 0; // these will be added by add_signal()
  sedf.header.ns_all = ns; // check this... should only matter for EDF access, so okay... 
  sedf.header.nr = edf.header.nr_all = nr_summ;  // likewise, value of nr_all should not matter, but set anyway
  sedf.header.record_duration = recdur_summ; // i.e. epoch length
  sedf.header.record_duration_tp = edf.header.record_duration * globals::tp_1sec;

  //
  // create a (continuous) timeline  [ TODO: see about SEDF+ for masked EDF ] 
  //

  logger << " adding timeline\n";

  sedf.set_edf();
  sedf.set_continuous();
  sedf.timeline.init_timeline();

  //
  // resize data[][], by adding empty records (one per SEDF record == EDF epoch )
  //

  logger << " adding records\n";

  for (int r=0;r<ne;r++)
    {
      edf_record_t record( &sedf ); 
      sedf.records.insert( std::map<int,edf_record_t>::value_type( r , record ) );
    }

  logger << " adding signals\n";

  //
  // add signals (this populates channel-specific 
  //

  std::map<std::string,std::vector<double> >::const_iterator ss = stats.begin();
  while ( ss != stats.end() )
    {

      // -1 implies 1 sample per record, i.e. if positive would be the SR
      // for that signal, but this allows slower channels, by directly specifying
      // the samples per record, (rather than samples per second)
      
      sedf.add_signal( ss->first , -1  , ss->second );

      ++ss;
    }

  // // arbitrary, but need epochs if not set it seems
  // if ( !edf.timeline.epoched() ) 
  //   edf.timeline.set_epoch( 30 , 30 );

  // // if a mask has been set, this will restructure the mask
  // edf.restructure(); 

  //
  // Save this SEDF
  //

  bool saved = sedf.write( filename );

  if ( ! saved ) Helper::halt( "problem trying to write " + filename );

  //
  // Sample list?
  //

  if ( param.has("sample-list") )
    {	  
      std::string file = param.value("sample-list");
      logger << " appending " << filename << " to sample-list " << file << "\n";      
      std::ofstream FL( file.c_str() , std::ios_base::app );
      FL << edf.id << "\t" << filename << "\n";
      FL.close();
    }
  

}

