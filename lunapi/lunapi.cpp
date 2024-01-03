
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


// 1) upload signal
// 2) db -> rrtables 
// 3) error/handling/checkpoints

#include "lunapi/lunapi.h"

extern logger_t logger;

extern globals global;

void lunapi_bail_function( const std::string & msg )
{
  logger << "*** [lunapi] error: " << msg << "\n";
  return;
}
  
void lunapi_msg_function( const std::string & msg )
{
  std::cerr << " [lunapi] :: " << msg << "\n";
  return;
}

void lunapi_t::silence( const bool b)
{
  global.cache_log = b;
}

void lunapi_t::init()
{

  global.init_defs();

  globals::bail_function = &lunapi_bail_function;
  
  globals::bail_on_fail = false;
  
  global.R( 1 ); // 1 means to cache the log
  
  writer.nodb();
  
  logger << "** luna " 
   	 << globals::version 
   	 << " " << globals::date
   	 << "\n";

  logger.print_buffer();
  
}

std::string lunapi_t::get_id() const
{
  return id;
}

int lunapi_t::get_state() const
{
  // 0 empty; +1 attached okay, -1 problem
  return state;
}

std::string lunapi_t::get_edf_file() const
{
  return edf_filename;
}

std::string lunapi_t::get_annot_files() const
{
  return Helper::stringize( annot_filenames );
}

void lunapi_t::var( const std::string & key , const std::string & value )
{
  
  // special treatment for `sig`.   This will just append
  // to an existing signallist;  as we don't always want to have
  // to lreset(), make this one case so that sig clears signlist prior 
  // to setting if the signal list is "." 
  
  if ( key == "sig" && value == "." )
    cmd_t::signallist.clear();
  else
    cmd_t::parse_special( key , value );
  
}


void lunapi_t::dropvars( const std::vector<std::string> & keys )
{
  for (int i=0; i<keys.size(); i++) dropvar( keys[i] );
}

void lunapi_t::dropvar( const std::string & key )
{
  std::map<std::string,std::string>::iterator ii = cmd_t::vars.find( key );
  if ( ii != cmd_t::vars.end() )
    cmd_t::vars.erase( ii );
  return;  
}

std::map<std::string,std::variant<std::monostate,std::string> > lunapi_t::vars( const std::vector<std::string> & keys )
{
  std::map<std::string,std::variant<std::monostate,std::string> > r;
  for (int i=0; i<keys.size(); i++) 
    r[ keys[i] ] = var( keys[i] );
  return r;
}

std::variant<std::monostate,std::string> lunapi_t::var( const std::string & key )
{  
  if ( cmd_t::vars.find( key ) == cmd_t::vars.end() ) return std::monostate{};
  return cmd_t::vars[ key ];
}


void lunapi_t::reset()
{
  globals::problem = false;  
  // globals::empty = false;
}


void lunapi_t::refresh()
{
  
  if ( state != -1 ) 
    {
      Helper::halt( "lunapi_t::refresh(): no attached EDF" );
      return;      
    }
    
  // drop edf_t  
  edf.init();
  
  // clear problem flags
  reset();

  // reattach EDF
  attach_edf( edf_filename );
  
  if ( state != 1 ) 
    {
      Helper::halt( "lunapi_t::refresh(): problem reattaching EDF" );
      return;
    }
  
  // reload annotations
  std::set<std::string>::const_iterator aa = annot_filenames.begin();
  while ( aa != annot_filenames.end() )
    {
      edf.load_annotations( *aa );
      ++aa;
    }
  
}

void lunapi_t::drop()
{
  edf.init();
  state = 0;
  edf_filename = "";
  annot_filenames.clear();
}

void lunapi_t::clear()
{
  // clear all variables (both user-defined and special)
  // but do not alter the EDF attachment

  // clear all user-defined variables, signal lists and aliases  
  cmd_t::clear_static_members();
  
  // also reset global variables that may have been changed since
  global.init_defs();
  
  // but need to re-indicate that we are running inside R with no log as default
  global.R( 0 ); // 0 means no log mirroring                                                                                                                        

}

std::map<std::string,datum_t> lunapi_t::status() const 
{
  // datum_t: std::variant<std::monostate{} , double , int , std::string, std::vector<double> , std::vector<int> , std::vector<std::string> >

  std::map<std::string,datum_t> r;

  r[ "state" ]  = state;
  
  if ( state != 1 ) return r;
  
  r[ "edf_file" ] = edf_filename;
  
  r[ "annotation_files" ] = Helper::stringize( annot_filenames );

  int n_data_channels = 0;
  int n_annot_channels = 0;
  for (int i=0;i<edf.header.ns;i++)
    {
      if ( edf.header.is_data_channel( i ) ) ++n_data_channels;
      else ++n_annot_channels;
    }

  r[ "id" ] = edf.id;
  r[ "ns" ] = n_data_channels;
  r[ "nt" ] = edf.header.ns_all;
  r[ "na" ] = (int)edf.timeline.annotations.names().size();
  
  // Record duration, as hh:mm:ss string                                                                                                                            
  uint64_t duration_tp = globals::tp_1sec
    * (uint64_t)edf.header.nr
    * edf.header.record_duration;  
  std::string total_duration_hms = Helper::timestring( duration_tp );
  
  r [ "duration" ] = total_duration_hms;

  // epoch/mask info
  if ( edf.timeline.epoched() )
    {
      r[ "ne" ] = edf.timeline.num_epochs();
      r[ "elen" ] = edf.timeline.epoch_length();
      r[ "nem" ] = edf.timeline.num_total_epochs() - edf.timeline.num_epochs();
    }

  return r;
  
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

  edf_filename = filename;

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
		  annot_filenames.insert( annotfile + fname );
                }
            }
          closedir (dir);
        }
      else
        {
          Helper::halt( "could not open folder " + annotfile );
	  return false;
        }
    }
  
  //
  // else a single file, load it
  //
  else
    {
      edf.load_annotations( annotfile );
      annot_filenames.insert( annotfile );        
    }
  
  return true;

}


// #1 eval returning all output to caller 
std::tuple<std::string,rtables_return_t> lunapi_t::eval_return_data( const std::string & cmdstr )
{  
  const std::string s = eval( cmdstr );
  return std::make_tuple( s , rtables.data() );
}



// #2: eval, but not returning outputs to caller (stores in lunapi_t::rtables)
std::string lunapi_t::eval( const std::string & cmdstr )
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

  writer.id( id , edf_filename );

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
  // get any results
  //

  rtables = rtables_t( ret );

  //
  // was a problem flag set?
  //

  if ( globals::problem ) 
    Helper::halt( "problem flag set: likely no unmasked records left?" );
  
  //
  // all done
  //
  
  return logger.print_buffer();;
  
}


rtable_t lunapi_t::table( const std::string & cmd , const std::string & faclvl ) const
{
  return rtables.table( cmd , faclvl );
}

std::vector<std::string> lunapi_t::variables( const std::string & cmd , const std::string & faclvl ) const
{
  return rtables.table( cmd , faclvl ).cols;
}

rtable_return_t lunapi_t::results( const std::string & cmd , const std::string & faclvl ) const
{
  return rtables.data( cmd , faclvl );
}

rtables_return_t lunapi_t::results() const
{
  return rtables.data() ;
}

//
// fetch signal data :
//   - given either list of epochs or intervals
//   - either combining all data into a single frame, or keeping separate


lint_t lunapi_t::epochs2intervals( const std::vector<int> & epochs )
{

  lint_t r;

  if ( state != 1 )
    return r;
  
  edf.timeline.ensure_epoched();  
  
  int total_epochs = edf.timeline.num_total_epochs();
  
  for (int epoch=0;epoch<epochs.size();epoch++)
    {
      // passed 1-based e-counts
      if ( epochs[epoch] < 1 || epochs[epoch] > total_epochs ) continue;

      // internally, 0-based
      int epoch0 = epochs[epoch] - 1;
      
      interval_t interval = edf.timeline.epoch( epoch0 );

      r.push_back( std::make_tuple( interval.start , interval.stop ) );
    }
  
  return r;
}
  
lint_t lunapi_t::seconds2intervals( const std::vector<std::tuple<double,double> > & s )
{
  lint_t r;
  
  for (int i=0;i<s.size();i++)
    r.push_back( std::make_tuple( std::get<0>(s[i]) * globals::tp_1sec ,
				  std::get<1>(s[i]) * globals::tp_1sec ) );
  
  return r;
}
  


bool lunapi_t::proc_channots( const std::string & chstr ,
			      const std::string & anstr ,
			      std::vector<std::string> * columns,
			      signal_list_t * signals , 
			      std::map<std::string,int> * atype )
{
  
  if ( state != 1 )
    return false;

  // Annotations: 0 not found, 1 interval, 2 epoch
  //  do not support epoch-annots right now
  std::vector<std::string> ans = Helper::parse( anstr , "," );
     
  for (int i=0;i<ans.size();i++)
    {
      if ( edf.timeline.annotations( ans[i] ) != NULL ) // is this an interval annotation? 
	(*atype)[ ans[i] ] = 1;
      else
	(*atype)[ ans[i] ] = 0;	  
    }
  
  // alphabetical order of annots:
  std::map<std::string,int>::const_iterator aa = atype->begin();
  while ( aa != atype->end() ) { columns->push_back( aa->first ) ; ++aa; } 
  
  // get signals  
  *signals = edf.header.signal_list( chstr );
  
  // check similar SRs  
  int fs = -1; 
  for (int s=0; s< signals->size(); s++) 
    {      
      if ( edf.header.is_data_channel( (*signals)(s) ) )
	{
	  columns->push_back( signals->label(s) );
	  if ( fs < 0 ) fs = edf.header.sampling_freq( (*signals)(s) );
	  else if ( edf.header.sampling_freq( (*signals)(s) ) != fs ) 
	    Helper::halt( "requires uniform sampling rate across signals" );	
	}
    }
  return true;
}


ldat_t lunapi_t::data( const std::vector<std::string> & chs ,			
		       const std::vector<std::string> & anns ,
		       const bool time_track )
{
  const interval_t whole = edf.timeline.wholetrace();
  lint_t w;
  w.push_back( std::make_tuple( whole.start , whole.stop ) );  
  return slice( w , chs , anns , time_track ); 
}


ldat_t lunapi_t::slice( const lint_t & intervals , 
			const std::vector<std::string> & chs ,			
			const std::vector<std::string> & anns ,
			const bool time_track )
{
  
  if ( state != 1 ) 
    return std::make_tuple( std::vector<std::string>(0), Eigen::MatrixXd::Zero(0,0) );

  const std::string chstr = Helper::stringize( chs );
  const std::string anstr = Helper::stringize( anns );
  
  // labels
  std::vector<std::string> columns;
  if ( time_track ) columns.push_back( "T" );  
  std::map<std::string,int> atype;     
  signal_list_t signals;
  
  // proc channels/annots
  if ( ! proc_channots( chstr , anstr , &columns , &signals, &atype ) )
    return std::make_tuple( std::vector<std::string>(0), Eigen::MatrixXd::Zero(0,0) );
  
  // pull data
  return std::make_tuple( columns, matrix_internal( intervals, signals, atype , time_track ) );
  
}


ldats_t lunapi_t::slices( const lint_t & intervals , 
			  const std::vector<std::string> & chs ,
			  const std::vector<std::string> & anns ,
			  const bool time_track )
{
  
  
  if ( state != 1 ) 
    return std::make_tuple( std::vector<std::string>(0), std::vector<Eigen::MatrixXd>(0) );

  const std::string chstr = Helper::stringize( chs );
  const std::string anstr = Helper::stringize( anns );

  std::vector<std::string> columns;
  if ( time_track ) columns.push_back( "T" );
  std::map<std::string,int> atype;
  signal_list_t signals;
  
  // get/check channel labels etc
  if ( ! proc_channots( chstr , anstr , &columns , &signals, &atype ) )
    return std::make_tuple( std::vector<std::string>(0), std::vector<Eigen::MatrixXd>(0) );
  
  // iterate over each interval
  std::vector<Eigen::MatrixXd> data;
  for ( int i=0; i<intervals.size(); i++)
    {
      lint_t i1( 1, intervals[i] );
      data.push_back( matrix_internal( i1 , signals, atype , time_track ) );
    }

  // return all 
  return std::make_tuple( columns , data );
  
}



Eigen::MatrixXd lunapi_t::matrix_internal( const lint_t & intervals , 
					   const signal_list_t & signals , 
					   const std::map<std::string,int> & atype ,
					   const bool time_track )
  

{
  
  const int ni = intervals.size();

  const int na = atype.size();

  // count signals
  int ns = 0;
  for (int s = 0 ; s < signals.size() ; s++ )
    if ( edf.header.is_data_channel( signals(s) ) ) ++ns;  
  if ( ns == 0 ) 
    Helper::halt( "requires at least one channel/data signal" );
    
  // # of columns: (T) + NS + NA   
  const int ncols = (int)time_track + ns + na;
  
  // # of rows... pull records 
  int nrows = 0;    
  for (int i=0;i<ni;i++)
    {      
      // Interval 
      const interval_t interval( std::get<0>(intervals[i]) , std::get<1>(intervals[i]) );
      
      // arbitary: first signal
      slice_t slice( edf , signals(0) , interval );      
      const std::vector<uint64_t> * tp = slice.ptimepoints();
      nrows += tp->size();      
    }

  // allocate matrix
  Eigen::MatrixXd X = Eigen::MatrixXd::Zero( nrows , ncols );  

  // first signal starts (after T and) annotations
  int s_col = (int)time_track + na;
  bool first = true;
  
  // Iterate over signals
  for (int s=0; s<ns; s++) 
    {
      
      uint64_t row = 0;
      
      // Consider each interval
      
      for (int i=0;i<ni;i++)
	{
	  
	  const interval_t interval( std::get<0>(intervals[i]) , std::get<1>(intervals[i]) );
	  
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
		  if (time_track) 
		    X(row,0) = (*tp)[r] * globals::tp_duration ;
		  
		  // Annotations (0/1) E,S 		  
		  int a_col = (int)time_track;	// start at 1 or 0
		  
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
		      // else if ( aa->second == 2 )
		      // 	X(row,a_col) = (int)( edf.timeline.epoch_annotation( aa->first , epoch0 ) ? 1 : 0 ) ;
		      
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



std::vector<std::string> lunapi_t::import_db( const std::string & dbfile )
{
  std::set<std::string> ids0;
  return import_db( dbfile , ids0 );
}


std::vector<std::string> lunapi_t::import_db( const std::string & dbfile , const std::set<std::string> & ids )
{

  // this gets populated by the IDs actually read
  std::vector<std::string> obs_ids;

  if ( ! Helper::fileExists( dbfile ) ) return obs_ids;
    
  retval_t ret = writer_t::dump_to_retval( dbfile , &ids , &obs_ids );
  
  logger << "  read data on " << obs_ids.size() << " individuals from " << dbfile << "\n" ;

  // store in the internal rtables cache
  rtables = rtables_t( ret );
  
  return obs_ids;

}


//
// Inserters 
//

bool lunapi_t::update_signal( const std::string & label , const std::vector<double> & x )
{
  if ( state != 1 ) return false;
  if ( ! edf.header.has_signal( label ) ) return false;
  const int slot = edf.header.signal( label );
  
  // void update_signal_retain_range( int s , const std::vector<double> * );

  // void update_signal( int s , const std::vector<double> * , int16_t * dmin = NULL , int16_t * dmax = NULL , 
  // 		      const double * pmin = NULL , const double * pmax = NULL );
  
  edf.update_signal( slot , &x ); 
  
  return true;
}


bool lunapi_t::insert_signal( const std::string & label , const std::vector<double> & x , const int sr )
{
  if ( state != 1 ) return false;

  // edf.add_signal( const std::string & label , const int n_samples , const std::vector<double> & data ,
  //                  double pmin = 0 , double pmax = 0 ,
  //                  int16_t dmin = 0 , int16_t dmax = 0 );
  
  edf.add_signal( label , sr , x );
    
  return true;
}


bool lunapi_t::insert_annotation( const std::string & class_label , const std::vector<std::tuple<double,double > > & x , const bool durcol2 )
{

  if ( state != 1 ) return false;
  if ( x.size() == 0 ) return false;
  if ( class_label == "" ) return false;
  
  const int n = x.size();

  // okay if class_label already exists, this will append new intervals
  annot_t * annot = edf.timeline.annotations.add( class_label );
  
  for (int i=0; i<n; i++ )
    {
      // skip bad elements
      if ( std::get<0>( x[i] ) < 0 || std::get<1>(x[i]) < 0 ) continue;
      
      const uint64_t start = std::get<0>(x[i]) * globals::tp_1sec ;
      const uint64_t stop  = std::get<1>(x[i]) * globals::tp_1sec + ( durcol2 ? start : 0LLU ) ; 
      
      annot->add( "." , // dummy instance ID
                  interval_t( start, stop ) ,
		  "." );  // channel ID dummy
    }

  return true;
}
