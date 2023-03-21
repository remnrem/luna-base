
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


#include "helper/helper.h"
#include "db/db.h"

#include "defs/defs.h"
#include "edf/edf.h"
#include "edf/slice.h"

#include <iostream>

extern writer_t writer;

bool identical_headers( const edf_header_t & h1 , const edf_header_t & h2 ) ;

void Helper::merge_EDFs( const std::vector<std::string> & tok )
{

  // by default, get start times from each EDF:
  //  a) check they line up,
  //  b) if there are gaps, add as an EDF+D (i.e. add explicit time track)
  //  c) if there are overlapps, give an error
  
  // otherwise -- ignore times except for the earliest/first file, which is taken as the start
  //  all other files are appended in that order
  // will give a warning if the first file is not the actual first

  
  std::vector<edf_t*> edfs;
  
  std::string id = "merged1";
  std::string filename = "merged.edf";
  std::string slist = "";
  bool use_fixed_order = false; // "fixed"  
    

    
  // expecting a list of file names 
  // but can also include key=value pairs, with keys:
  //   id
  //   edf
  //   sample-list
  //   fixed=T
    
  for (int i=0; i<tok.size(); i++)
    {

      logger << "------------------------------------------------------------\n"
	     << "processing [" << tok[i] << "]\n";
      
      std::vector<std::string> tok2 = Helper::quoted_parse( tok[i] , "=" );
      if ( tok2.size() == 2 )
	{
	  if ( tok2[0] == "id" ) id = tok2[1];
	  else if ( tok2[0] == "edf" ) filename = tok2[1];
	  else if ( tok2[0] == "sample-list" ) slist = tok2[1];
	  else if ( tok2[0] == "fixed" ) use_fixed_order = Helper::yesno( tok2[1] );
	  continue;
	}
      
      const std::string fname = Helper::expand( tok[i] );
      
      if ( ! Helper::fileExists( fname ) )
	{
	  logger << "  ** warning: could not attach " << fname << "\n";
	  continue;
	}
      
      edf_t * edf = new edf_t;
      
      const std::string id = "id" + Helper::int2str( (int)( edfs.size() + 1 ) ) ;
      
      bool okay = edf->attach( fname , id );
      
      if ( ! okay )
	{
	  logger << " ** could not attach " << filename	<< "\n";
          continue;	  
	}

      // only allow starndard EDFs to be merged right now

      if ( edf->header.edfplus ) 
	Helper::halt( "cannot merged EDF+ files : " + fname
		      + "\n (this constraint can be relaxed in future)");
      
      logger << "  attached component EDF " << fname << "\n";

      edfs.push_back( edf ) ;

      logger << "\n";
    }

  const int nf = edfs.size();

   
  //
  // Some outputs
  //

  logger << "------------------------------------------------------------\n"
	 << "\n  attached " << nf << " EDFs\n";  

  // for (int i=0; i<nf; i++)
  //   logger << "  EDF " << i+1 << "\t"
  // 	   << edfs[i]->id << "\t"
  // 	   << edfs[i]->header.startdate << "\t"
  // 	   << edfs[i]->header.starttime << "\n";

  logger  << "  writing merged data:\n"
	  << "     ID           : " << id << "\n"
	  << "     EDF filename : " << filename << "\n";

  if ( slist != "" )
   logger << "     sample-list  : " << slist << "\n"; 
  

  //
  // Get the ordering; use midnight (start) of first day as the primary
  // anchor for times (to the second)
  //

  // encode time as seconds past 1/1/95 00:00
  
  std::map<uint64_t,int> time2edf;

  bool gapped = false;
  
  if ( ! use_fixed_order )
    {
      logger << " ------------------------------------------------------------\n"
	     << " extracting start times from EDFs --> seconds since 1/1/85 00:00\n";
      
      for (int i=0; i<nf; i++)
	{
	  date_t date( edfs[i]->header.startdate );
	  clocktime_t clock( edfs[i]->header.starttime );
	  
	  // seconds past 1/1/85 00:00
	  uint64_t days = date_t::count( date );
	  uint64_t secs = days * 24 * 60 * 60 + clock.rounded_seconds();
	  
	  logger << "  EDF " << i 
		 << "  date: " << edfs[i]->header.startdate
		 << " time: " << edfs[i]->header.starttime
		 << " days: " << days
		 << " --> secs:" << secs << "\n";
	  
	  // check this time does not already exist
	  if ( time2edf.find( secs ) != time2edf.end() )
	    Helper::halt( "EDFs with identical start times specified: " + edfs[i]->header.startdate + " " + edfs[i]->header.starttime );
	  
	  // track
	  time2edf[ secs ] = i;
	  
	}

      //
      // go through order and figure out whether we have a) overlap and/or b) gaps
      //

      bool overlap = false;
      
      // currently, only allow merging of standard EDFs (i.e. simple integer start times)
      std::map<uint64_t,int>::const_iterator ss = time2edf.begin();

      std::map<uint64_t,int>::const_iterator pp = ss;
      ++ss;

      
      logger << "------------------------------------------------------------\n"
	     << "  ordered EDFs (seconds past 1/1/1985 00:00:00)\n";
	
	while ( ss != time2edf.end() )
	{
	  // current : ss
	  // previous : pp
	  
	  // diff between starts
	  double diff = ss->first - pp->first;
	  
	  // implied duration of previous
	  double dur = edfs[ pp->second ]->header.nr * edfs[ pp->second ]->header.record_duration;
	  
	  logger << "  ordered EDFs : prev start = " << pp->first << " ; this start = " << ss->first << "\n"; 
	  logger << "    implied duration of previous based on # records         = " << diff << "\n";
	  logger << "    implied duration of previous based on current EDF start = " << dur << "\n"; 
	  
	  if ( dur - diff > 0.5 )
	    {
	      overlap = true;
	      logger << "  *** warning -- overlapping EDFs implied\n";
	    }
	  else if ( diff - dur > 0.5 )
	    {
	      gapped = true;
	      logger << "  implies gap between previous and current - will output an EDF+D\n";
	    }
	  else
	    logger << "   implies exactly contiguous EDFs\n";
	  
	  logger << "\n";
	  
	  ++ss;
	  ++pp;
	}
	
	if ( overlap )
	  Helper::halt( "found overlapping EDFs -- bailing, cannot merge" );
	
	if ( gapped )
	  logger << " found gaps between EDFs - will generate an EDF+D\n";
	else
	  logger << " no gaps found between EDFs - will generate an EDF (or EDF+C)\n";
	
    }
  
  
  //
  // Check that all headers are compatible
  //   -- initially, headers must be identical
  //   -- this can be relaxed in the future.
  //

  for (int i=1; i<nf; i++)
    if ( ! identical_headers( edfs[0]->header , edfs[i]->header ) )
      Helper::halt( "headers incompatible:" + edfs[0]->filename + " " + edfs[i]->filename );
  
  logger  << "  good, all EDFs have merge-compatible headers\n";
  

  //
  // Get total implied NR for new EDF
  //
  
  int nr = 0;
  
  for (int i=0; i<nf; i++)
    nr += edfs[i]->header.nr;
  
  logger  << "  expecting " << nr
	  << " records (each of "
	  << edfs[0]->header.record_duration << " sec) in the new EDF\n";
  
    
  //
  // Create the new merged EDF
  //
  
  edf_t medf;
  
  int first_edf = use_fixed_order ? 0 : time2edf.begin()->second;
  
  //
  // Set header
  //

  medf.header.version = edfs[ first_edf ]->header.version;
  medf.header.patient_id = id;
  medf.header.recording_info = edfs[ first_edf ]->header.recording_info;
  medf.header.startdate = edfs[ first_edf ]->header.startdate;
  medf.header.starttime = edfs[ first_edf ]->header.starttime;
  medf.header.nbytes_header = edfs[ first_edf ]->header.nbytes_header;
  medf.header.ns = 0; // this is popuilated by edf_t::add_signal()
  medf.header.ns_all = 0; // this is popuilated by edf_t::add_signal()
  medf.header.nr = medf.header.nr_all = nr; // this is sum across all EDFs computed above
  medf.header.record_duration = edfs[ first_edf ]->header.record_duration;
  medf.header.record_duration_tp = edfs[ first_edf ]->header.record_duration_tp;

  
  //
  // create a (continuous/discontinuous) timeline  
  //

  logger << "  adding timeline\n";

  // initially, create the new dataset as a standard EDF
  // i.e. which assumes that all component files are also all EDF
  
  medf.set_edf();
  medf.set_continuous();
  medf.timeline.init_timeline();

 
  //
  // resize data[][], by adding empty records 
  //

  logger << "  adding " << nr << " empty records...\n";

  for (int r=0;r<nr;r++)
    {
      edf_record_t record( &medf ); 
      medf.records.insert( std::map<int,edf_record_t>::value_type( r , record ) );
    }
  
  
  //
  // add signals 
  //

  const int ns = edfs[0]->header.ns ;
  
  for (int s=0; s<ns; s++)
    {

      // skip annotations
      if ( edfs[0]->header.is_annotation_channel( s ) ) continue;
      
      logger << "  compiling channel " << edfs[0]->header.label[s] << "\n";

      std::vector<double> dt;

      for (int j=0; j<nf; j++)
	{

	  edf_t & edf = *(edfs[j]);
	  
	  // get whole signal
	  slice_t slice( edf , s , edf.timeline.wholetrace() );

	  std::vector<double> * d = slice.nonconst_pdata();

	  for (int k=0; k<d->size(); k++) dt.push_back( (*d)[k] ) ;
	  
	} // next EDF

      const int np_obs = dt.size();
      const int np_exp = nr * edfs[0]->header.sampling_freq(s);
      
      if ( np_obs != np_exp )
	Helper::halt( "expected and observed number of sample points did not align" );
      
      // add the signal to the merged EDF
      medf.add_signal( edfs[0]->header.label[ s ] , edfs[0]->header.sampling_freq( s ) , dt );
            
    } // next channel
  

  

  //
  // Add EDF Annotations channel w/ EDF+D times, if needed
  //

  if ( gapped )
    {

      // make a new implied timeline
      // all in tp units past the start of the first EDF
      // [ only populated if gapped == true ]
      std::vector<uint64_t> tps;

      // absolute seconds value (past 1/1/85)
      // i.e. to be subtracted from all subsequent, so
      // first rec starts at 0 seconds;
      uint64_t first_edf_secs = time2edf.begin()->first;
      
      std::map<uint64_t,int>::const_iterator ss = time2edf.begin();
      while ( ss != time2edf.end() )
	{
	  int e = ss->second;
	  uint64_t tp_past_edf_start = globals::tp_1sec * ( ss->first - first_edf_secs );
	  int nr1 = edfs[e]->header.nr;
	  uint64_t tp = edfs[e]->header.record_duration_tp;
	  
	  for (int r=0; r<nr1; r++)
	    {
	      tps.push_back( tp_past_edf_start );
	      tp_past_edf_start += tp;
	    }
	  // next EDF
	  ++ss;
	}

      // update in-memory EDF w/ this information
      
      medf.set_edfplus(); // make EDF+ (adds a time-track and at)
      medf.set_discontinuous(); // set as EDF+D explicitly

      // now update in-memory time-track
      medf.timeline.create_discontinuous_timeline( tps );
      
      // now add EDF annotations w/ explcitly calculated tps
      medf.add_time_track( &tps );


      //
      // Output tidy table w/ information for annotation offsets      
      //

      ss = time2edf.begin();
      while ( ss != time2edf.end() )
        {
          int e = ss->second;
          double sec_past_edf_start = ss->first - first_edf_secs;

	  std::cout << "ANNOT-OFFSET" << "\t"
		    << edfs[e]->filename << "\t"  
		    << sec_past_edf_start << "\n";
          ++ss;
        }
      
    }
  
  //
  // Save this merged EDF
  //

  logger << "  writing merged EDF as " << filename << "\n";
  
  bool saved = medf.write( filename );
  
  if ( ! saved ) Helper::halt( "problem trying to write " + filename );
  
  //
  // Sample list?
  //

  if ( slist != "" )
    {	  
      logger << " appending " << filename << " to sample-list " << slist << "\n";
      std::ofstream FL( slist.c_str() , std::ios_base::app );
      FL << medf.id << "\t" << filename << "\n";
      FL.close();
    }
  
  //
  // check record table
  //

  if ( 0 )
    {
      param_t param1;
      medf.record_table( param1 );
    }



  //
  // Clean up
  //
  
  for (int i=0; i<nf; i++)
    delete edfs[i];      
  
}


bool identical_headers( const edf_header_t & h1 , const edf_header_t & h2 )
{

  // compatible for merging?  currently needs to have *identical*
  // headers (for signal #, SR and EDF record size)

  if ( h1.version != h2.version ) return false;
  
  if ( h1.ns != h2.ns ) return false;
  
  if ( h1.record_duration_tp != h2.record_duration_tp ) return false;

  for ( int s = 0 ; s < h1.ns; s++ )
    {
      if ( h1.label[s] != h2.label[s] ) return false;
      if ( h1.n_samples[s] != h2.n_samples[s] ) return false;
    }
  
  return true;
  
}
