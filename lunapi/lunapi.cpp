
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

#include "lunapi/lunapi.h"

extern logger_t logger;

extern globals global;

void lunapi_bail_function( const std::string & msg )
{
  logger << "*** [lunapi] error: " << msg << "\n";
  return;
}
  

void lunapi_t::init()
{

  global.init_defs();

  globals::bail_function = &lunapi_bail_function;

  //globals::logger_function = &lunapi_message_function;

  // turn off external output
  global.api();
  
  logger << "** lunapi " 
	 << globals::version 
	 << " " << globals::date
	 << "\n";
  
}


void lunapi_t::reset()
{
  globals::problem = false;
  globals::empty = false;
}



bool lunapi_t::attach_edf( const std::string & filename )
{
  if ( ! Helper::fileExists( filename ) ) 
    Helper::halt( "cannot find " + filename );

  bool okay = edf.attach( filename , id , NULL );

  if ( ! okay )
    {
      state = -1;
      return false;
    }

  // EDF+ annotations?
  if ( edf.header.edfplus )
    {
      // must read if EDF+D (but only the time-track will be taken in)                                          
      // if EDF+C, then look at 'skip-edf-annots' flag                                                          
      if ( edf.header.continuous && ! globals::skip_edf_annots )
        edf.timeline.annotations.from_EDF( edf , edf.edfz_ptr() );
      else if ( ! edf.header.continuous )
        edf.timeline.annotations.from_EDF( edf , edf.edfz_ptr() );
    }

  cmd_t::define_channel_type_variables( edf );

  state = 1 ; 
  
  return true;
  
}


bool lunapi_t::attach_annot( const std::string & annotfile )
{

  if ( annotfile.size() == 0 ) return false;

  // is 'annotfile' in fact a folder (i.e. ending in '/') ?
  
  if ( annotfile[ annotfile.size() - 1 ] == globals::folder_delimiter )
    {
      
      // this means we are specifying a folder, in which case search for all files that
      // start id_<ID>_* and attach those

      DIR * dir;
      struct dirent *ent;
      if ( (dir = opendir ( annotfile.c_str() ) ) != NULL )
        {
          /* print all the files and directories within directory */
          while ((ent = readdir (dir)) != NULL)
            {
              std::string fname = ent->d_name;
	      
              if ( Helper::file_extension( fname , "ftr" ) ||
                   Helper::file_extension( fname , "xml" ) ||
                   Helper::file_extension( fname , "eannot" ) ||
                   Helper::file_extension( fname , "annot" ) )
                {
                  edf.load_annotations( annotfile + fname );
                }
            }
          closedir (dir);
        }
      else
        {
          Helper::halt( "could not open folder " + annotfile );
        }
    }
  
  //
  // else a single file, load it
  //
  else
    {
      edf.load_annotations( annotfile );
    }
  
  return true;

}

rtables_t lunapi_t::eval( const std::string & cmdstr )
{

  //
  // set up retval_t mechanism to catch outputs
  //

  retval_t ret;
  
  writer.clear();
  writer.set_types(); // not sure this is needed now...
  writer.use_retval( &ret );

  //
  // set ID 
  //

  writer.id( id , filename );

  logger << "evaluating...\n";

  // track errors? (for moonlight) 
  // R_last_eval_failed = false;
  // R_last_eval_errmsg = "";

  //
  // set command string
  //
  
  cmd_t cmd( cmdstr );
  
  //
  // replace any variables (or @includes, conditionals,etc) into command
  //
  
  cmd.replace_wildcards( id );

  //
  // eval on the current EDF
  //
  
  cmd.eval( edf );
  
  
  //
  // switch off the retval stream (which is local to this function and
  // so will be deleted when leaving this scope) and clear the writer
  // (ensures prior strata not applied to next run)
  //

  writer.use_retval( NULL );


  writer.clear();
  writer.set_types();

  //
  // was a problem flag set?
  //

  if ( globals::problem ) 
    Helper::halt( "problem flag set: likely no unmasked records left?" );

  
  //
  // all done
  //


  rtables_t rtables( ret );

  return rtables;
  
}




//
// fetch signal data : given either list of epochs or intervals
//

Eigen::MatrixXd lunapi_t::slice_epochs( const std::vector<int> & epochs , 					
					const std::string & chstr ,
					const std::string & anstr ,
					std::vector<std::string> * columns )
{
  
  if ( state != 1 ) 
    return Eigen::MatrixXd::Zero(0,0);

  columns->resize(1,"T");


  // Annotations: 0 not found, 1 interval, 2 epoch
  std::vector<std::string> ans = Helper::parse( anstr , "," );
  std::map<std::string,int> atype;   
  for (int i=0;i<ans.size();i++)
    {
      if ( edf.timeline.annotations( ans[i] ) != NULL ) // is this an interval annotation? 
	atype[ ans[i] ] = 1;
      else if ( edf.timeline.epoch_annotation( ans[i] ) ) // or an epoch-annotation?
	atype[ ans[i] ] = 2;
      else
	atype[ ans[i] ] = 0;	  
    }
  
  // alphabetical order of annots:
  std::map<std::string,int>::const_iterator aa = atype.begin();
  while ( aa != atype.end() ) { columns->push_back( aa->first ) ; ++aa; } 

  // get signals  
  signal_list_t signals = edf.header.signal_list( chstr );

  // check similar SRs  
  int fs = -1; 
  for (int s=0; s< signals.size(); s++) 
    {      
      if ( edf.header.is_data_channel( signals(s) ) )
	{
	  columns->push_back( signals.label(s) );
	  if ( fs < 0 ) fs = edf.header.sampling_freq( signals(s) );
	  else if ( edf.header.sampling_freq( signals(s) ) != fs ) 
	    Helper::halt( "requires uniform sampling rate across signals" );	
	}
    }

  // Epochs: TODO, add bound check
  edf.timeline.ensure_epoched();  
  int total_epochs = edf.timeline.num_total_epochs();
  std::vector<interval_t> epoch_intervals;
  for (int epoch=0;epoch<epochs.size();epoch++)
    {
      int epoch0 = epochs[epoch] - 1;
      interval_t interval = edf.timeline.epoch( epoch0 );       
      epoch_intervals.push_back( interval );
    }
    
  // Generate and return data.frame 
  return matrix_internal( epoch_intervals, &epochs , signals, atype );
  
}



Eigen::MatrixXd lunapi_t::slice_intervals( const std::vector<double> & secints ,
					   const std::string & chstr ,
					   const std::string & anstr ,
					   std::vector<std::string> * columns )
{
    
  if ( state != 1 ) 
    return Eigen::MatrixXd::Zero(0,0);

  columns->resize(1,"T");

  // Annotations: 0 not found, 1 interval, 2 epoch
  std::vector<std::string> ans = Helper::parse( anstr , "," );
  std::map<std::string,int> atype;   
  for (int i=0;i<ans.size();i++)
    {
      if ( edf.timeline.annotations( ans[i] ) != NULL ) // is this an interval annotation? 
	atype[ ans[i] ] = 1;
      else if ( edf.timeline.epoch_annotation( ans[i] ) ) // or an epoch-annotation?
	atype[ ans[i] ] = 2;
      else
	atype[ ans[i] ] = 0;	  
    }

  // alphabetical order of annots:
  std::map<std::string,int>::const_iterator aa = atype.begin();
  while ( aa != atype.end() ) { columns->push_back( aa->first ) ; ++aa; } 
  
  // get signals  
  signal_list_t signals = edf.header.signal_list( chstr );

  // check similar SRs  
  int fs = -1; 
  for (int s=0; s< signals.size(); s++) 
    {      
      if ( edf.header.is_data_channel( signals(s) ) )
	{
	  columns->push_back( signals.label(s) );
	  if ( fs < 0 ) fs = edf.header.sampling_freq( signals(s) );
	  else if ( edf.header.sampling_freq( signals(s) ) != fs ) 
	    Helper::halt( "requires uniform sampling rate across signals" );	
	}
    }

  // Intervals: expects a simple dbl vector where START1 STOP1 START2 STOP2 etc
  //   z <- as.numeric(rbind( i$START , i$STOP )
  // in SECONDS

  std::vector<interval_t> intervals;
  
  if ( secints.size() % 2 ) 
    Helper::halt( "internal error, expecting an even sized list");
  
  for (int e=0;e<secints.size();e+=2)
    {
      if ( secints[e+1] < secints[e] ) 
	Helper::halt( "internal error, expecting an even sized list"); 

      // as here the intervals are coming *from* R / the user, we can 
      // assume they will be in one-past-the-end format already
      interval_t interval( secints[e] * globals::tp_1sec , 
			   secints[e+1] * globals::tp_1sec );
      
      intervals.push_back( interval );
    }

  // Generate matrix
  return matrix_internal( intervals, NULL , signals, atype );
  
}


Eigen::MatrixXd lunapi_t::matrix_internal( const std::vector<interval_t> & intervals , 
					   const std::vector<int> * epoch_numbers , 
					   const signal_list_t & signals , 
					   const std::map<std::string,int> & atype )
  

{

  // this supports *either* epochs or generic intervals

  // If epoch_numbers is defined, this it must be the same length
  // as intervals
  
  const bool emode = epoch_numbers != NULL;
  
  if ( emode && intervals.size() != epoch_numbers->size() )
    Helper::halt( "internal error in matrix_internal" );
  
  const int ni = intervals.size();
  const int na = atype.size();
  
  int ns = 0;
  for (int s = 0 ; s < signals.size() ; s++ )
    if ( edf.header.is_data_channel( signals(s) ) ) ++ns;

  // no signals:
  if ( ns == 0 ) 
    Helper::halt( "requires at least one channel/data signal" );
  
  // Point to first epoch
  if ( emode && ! edf.timeline.epoched() ) 
    {
      int n = edf.timeline.set_epoch( globals::default_epoch_len , globals::default_epoch_len );
      logger << " set epochs to default " << globals::default_epoch_len << " seconds, " << n << " epochs\n";
      edf.timeline.first_epoch();
    }
  

  // # of columns: T + NS + NA   
  const int ncols = 1 + ns + na;
  
  // # of rows... pull records 
  int nrows = 0;    
  for (int i=0;i<ni;i++)
    {      
      // Interval (convert to 0-base)
      const interval_t & interval = intervals[i];
      
      // arbitary: first signal
      slice_t slice( edf , signals(0) , interval );      
      const std::vector<uint64_t> * tp = slice.ptimepoints();
      nrows += tp->size();

    }


  // allocate matrix
  Eigen::MatrixXd X = Eigen::MatrixXd::Zero( nrows , ncols );
  

  // first signal starts after T and annotations
  int s_col = 1 + na;
  bool first = true;
  
  //
  // Iterate over signals
  //

  for (int s=0; s<ns; s++) 
    {
      
      uint64_t row = 0;
      
      // Consider each interval
      
      for (int i=0;i<ni;i++)
	{
	  
	  int epoch0 = emode ? (*epoch_numbers)[i] - 1 : 0 ; // used???
	  
	  const interval_t & interval = intervals[i];
	  
	  // Get data
	  
	  slice_t slice( edf , signals(s) , interval );
	  
	  const std::vector<double> * data = slice.pdata();
	  
	  const std::vector<uint64_t> * tp = slice.ptimepoints();

	  int nrows_per_interval = tp->size();
 
	  // Populate signals
	  
	  for (int r=0;r<nrows_per_interval;r++)
	    {
	      
	      // Only add T and annotations once
	      if ( first )
		{		    
		  
		  // elapsed time in seconds
		  X(row,0) = (*tp)[r] * globals::tp_duration ;
		  
		  // Annotations (0/1) E,S 		  
		  int a_col = 1;		  
		  std::map<std::string,int>::const_iterator aa = atype.begin();
		  while ( aa != atype.end() )
		    {		      
		      if ( aa->second == 0 )
			X(row,a_col) = std::numeric_limits<double>::quiet_NaN();
		      else if ( aa->second == 1 )
			{
			  // get exact point      
			  interval_t interval2 = interval_t( (*tp)[r] , (*tp)[r] + 1LLU );
			  annot_t * annot = edf.timeline.annotations( aa->first );
			  annot_map_t events = annot->extract( interval2 );
			  bool has_annot = events.size() ;
			  X(row,a_col) = (int)(has_annot ? 1 : 0 );
			}
		      else if ( aa->second == 2 )
			X(row,a_col) = (int)( edf.timeline.epoch_annotation( aa->first , epoch0 ) ? 1 : 0 ) ;
		      
		      // next annotation
		      ++a_col;
		      ++aa;
		    }
		
		} // end of special case (T/ANNOTS)
	      
	      
	      // Signal data	      
	      X(row,s_col) = (*data)[r];
	      
	      // next row 	      
	      ++row;
	    }

	  // Next interval
	}

      // no need to add T/ANNOTS again if looping back
      first = false;

      // advance to next signal
      ++s_col;
      
      // use Helper::sanitize() to make col names
    }
  
  return X;
  
}




std::vector<std::string> lunapi_t::channels()
{
  std::vector<std::string> chs;
  if ( state != 1 ) return chs;
  signal_list_t signals = edf.header.signal_list( "*" );
  const int ns = signals.size();
  for (int s=0;s<ns;s++) 
    if ( edf.header.is_data_channel( signals(s) ) ) 
      chs.push_back( signals.label(s) );
  return chs; 
}  

std::vector<std::string> lunapi_t::annots() const
{
  if ( state != 1 ) return std::vector<std::string>(0);
  return edf.timeline.annotations.names();
}

