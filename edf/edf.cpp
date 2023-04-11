
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

#include "edf.h"
#include "defs/defs.h"
#include "helper/helper.h"
#include "helper/logger.h"
#include "miscmath/miscmath.h"
#include "db/db.h"
#include "edfz/edfz.h"
#include "dsp/resample.h"

#include "slice.h"
#include "tal.h"
#include "eval.h"
#include "clocs/clocs.h"
#include "timeline/timeline.h"

//#include <ftw.h>

#include <iostream>
#include <fstream>

extern writer_t writer;
extern logger_t logger;

void writestring( const std::string & s , int n , FILE * file )
{
  std::string c = s;
  c.resize(n,' ');
  fwrite( c.data() , 1 , n , file );
}

void writestring( const int & s , int n , FILE * file )
{
  std::string c = Helper::int2str(s);
  c.resize(n,' ');
  fwrite( c.data() , 1 , n , file );
}

void writestring( const double & s , int n , FILE * file )
{
  std::string c = Helper::dbl2str_fixed(s,n);
  c.resize(n,' ');
  fwrite( c.data() , 1 , n , file );
}

edf_t::endian_t edf_t::endian = edf_t::MACHINE_LITTLE_ENDIAN;

uint64_t edf_t::get_filesize(FILE *file)
{  
  uint64_t lCurPos, lEndPos;
  lCurPos = ftell(file);
  fseek(file, 0, 2);
  lEndPos = ftell(file);
  fseek(file, lCurPos, 0);
  return lEndPos;
}

int edf_t::get_int( byte_t ** p , int sz )
{
  std::string s = edf_t::get_string( p , sz );
  int t = 0;
  if ( ! Helper::str2int( s , &t ) ) 
    Helper::halt( "problem converting to an integer value: [" + s + "]"  );
  return t;
}

double edf_t::get_double( byte_t ** p , int sz )
{
  std::string s = edf_t::get_string( p , sz );

  double t = 0;
  if ( s == "" ) return -1;
  
  if ( ! Helper::from_string<double>( t , s , std::dec ) ) 
    {     
      logger << "returning -1: [" << s << "] is not a valid real number\n";
      return -1;
    }
  return t;
}

std::string edf_t::get_string( byte_t ** p , int sz )
{
  // only US-ASCII printable characters allowed: 32 .. 126 
  // other characters mapped to '?'
  std::vector<char> buf(sz+1);
  for (int i=0;i<sz;i++)
    {
      buf[i] = **p;
      if ( buf[i] < 32 || buf[i] > 126 ) buf[i] = 63; // '?'
      ++(*p);      
    }
  buf[sz] = '\0';
  std::string str = &buf[0];
  // trim trailing whitespace 
  // (when writing header back out, we expand whitespace to desired length)
  Helper::rtrim(str);
  return str;
}

void edf_t::skip( byte_t ** p , int sz ) 
{
  (*p) += sz;
}

std::vector<char> edf_t::get_bytes( byte_t ** p , int sz )
{
  std::vector<char> buf(sz);
  for (int i=0;i<sz;i++)
    {
      buf[i] = **p;
      ++(*p);
    }
  return buf;
}


inline double edf_record_t::dig2phys( int16_t d , double bv , double offset )
{  
  return bv * ( offset + d ) ; 
}

inline int16_t edf_record_t::phys2dig( double d , double bv , double offset )
{
  return d / bv - offset ; 
}


inline int16_t edf_record_t::tc2dec( char a , char b )
{        
  union 
  {
    int16_t       one[1];    
    unsigned char two[2];      
  } buffer;    
  
  if ( edf_t::endian == edf_t::MACHINE_LITTLE_ENDIAN )
    {
      buffer.two[0] = a;
      buffer.two[1] = b;
    }
  else
    {
      buffer.two[0] = b;
      buffer.two[1] = a;
    }
  return buffer.one[0];   
}


inline void edf_record_t::dec2tc( int16_t x , char * a , char * b )
{
  union 
  {
    int16_t       one[1];    
    unsigned char two[2];      
  } buffer;    

  buffer.one[0] = x;
  
  if ( edf_t::endian == edf_t::MACHINE_LITTLE_ENDIAN )
    {
      *a = buffer.two[0];
      *b = buffer.two[1];
    }
  else
    {
      *b = buffer.two[0];
      *a = buffer.two[1];
    }
  
}



std::string edf_header_t::summary() const
{

  std::stringstream ss;

  ss << "Patient ID     : " << patient_id << "\n"
     << "Recording info : " << recording_info << "\n"
     << "Start date     : " << startdate << "\n"
     << "Start time     : " << starttime << "\n"
     << "\n"
     << "# signals      : " << ns << "\n"
     << "# records      : " << nr << "\n"
     << "Rec. dur. (s)  : " << record_duration << "\n\n";
  
  for (int s=0;s<ns;s++)
    {
      
      ss << "Signal " << (s+1) << " : [" << label[s] << "]\n";
      
      std::string primary = label[s];
      
      // is alias? ( will have been mapped already )
      if ( cmd_t::primary_alias.find( primary ) != cmd_t::primary_alias.end() )
	{
	  std::string aliases = Helper::stringize( cmd_t::primary_alias[ primary ] , " | " );
	  ss << "\taliased from         : " << aliases << "\n"; 
	}
      
      if ( is_annotation_channel( s ) ) 
	ss << "\tannotation channel\n";
      else
	ss << "\tsampling rate        : " << n_samples[s] / (double)record_duration << " Hz\n"
	   << "\t# samples per record : " << n_samples[s] << "\n"	
	   << "\ttransducer type      : " << transducer_type[s] << "\n"
	   << "\tphysical dimension   : " << phys_dimension[s] << "\n"
	   << "\tmin/max (phys)       : " << physical_min[s] << "/" << physical_max[s] << "\n"
	   << "\tEDF min/max (phys)   : " << orig_physical_min[s] << "/" << orig_physical_max[s] << "\n"
	   << "\tmin/max (digital)    : " << digital_min[s] << "/" << digital_max[s] << "\n"
	   << "\tEDF min/max (digital): " << orig_digital_min[s] << "/" << orig_digital_max[s] << "\n"
	   << "\tpre-filtering        : " << prefiltering[s] << "\n\n";
      
    }

  return ss.str();
}

void edf_t::description( const param_t & param ) 
{

  signal_list_t signals = header.signal_list( param.requires( "sig" ) );
  
  bool channel_list = param.has( "channels" );
  
  if ( channel_list )
    {
      for (int s=0;s<signals.size();s++) 
	{
	  if ( header.is_data_channel( signals(s) ) )
	    std::cout << signals.label(s) << "\n";
	}
      return;
    }
    
  uint64_t duration_tp = globals::tp_1sec * (uint64_t)header.nr * header.record_duration ;

  int n_data_channels = 0 , n_annot_channels = 0;
  int n_data_channels_sel = 0 , n_annot_channels_sel = 0;
  
  for (int s=0;s<header.ns;s++) 
    {
      if ( header.is_data_channel(s) )
	++n_data_channels;
      else 
	++n_annot_channels;
    }

  for (int s=0;s<signals.size(); s++) 
    {
      if ( header.is_data_channel( signals(s) ) )
	++n_data_channels_sel;
      else 
	++n_annot_channels_sel;
    }

  clocktime_t et( header.startdate, header.starttime );
  if ( et.valid )
    {
      // go to next time point /after/ end
      double time_sec = ( timeline.last_time_point_tp+1LLU ) * globals::tp_duration ; 
      et.advance_seconds( time_sec );
    }
  
  std::cout << "EDF filename      : " << filename << "\n"
	    << "ID                : " << id << "\n";

  if ( header.edfplus ) 
    std::cout << "Header start time : " << header.starttime << "\n"
	      << "Last observed time: " << et.as_string() << "\n";
  else 
    std::cout << "Clock time        : " << header.starttime << " - " << et.as_string() << "\n";

  std::cout << "Duration          : "
	    << Helper::timestring( duration_tp , ':' , false )
	    << "  " 
	    << header.nr * header.record_duration << " sec" 
	    << "\n"; // not fractional 

  if ( header.edfplus && ! header.continuous )
    {
      clocktime_t st( header.startdate , header.starttime ); // include dates
      double diff_secs = clocktime_t::ordered_difference_seconds( st , et );
      clocktime_t ot( "00.00.00" );
      ot.advance_seconds( diff_secs );
      std::cout << "Duration (w/ gaps): "
		<< ot.as_string() << "  " << diff_secs << " sec\n";
    }
  
  if ( n_data_channels_sel < n_data_channels )
    std::cout << "# signals         : " << n_data_channels_sel << " selected (of " << n_data_channels << ")\n";
  else
    std::cout << "# signals         : " << n_data_channels << "\n";
    
  if ( n_annot_channels > 0 )
    {
      if ( n_annot_channels_sel < n_annot_channels )
	std::cout << "# EDF annotations : " << n_annot_channels_sel << " selected (of " << n_annot_channels << ")\n";
      else
	std::cout << "# EDF annotations : " << n_annot_channels << "\n";
    }
  
  std::cout << "Signals           :";

  int cnt=0;
  for (int s=0;s<signals.size();s++) 
    {
      if ( header.is_data_channel( signals(s) ) )
	std::cout << " " 
		  << signals.label(s) 
		  << "[" << header.sampling_freq( signals(s) ) << "]";
      if ( ++cnt >= 6 ) { cnt=0; std::cout << "\n                   "; } 
    }
  std::cout << "\n\n";
  
  
}


void edf_t::report_aliases() const
{
  // annotations
  std::map<std::string,std::string>::const_iterator aa = timeline.annotations.aliasing.begin();
  while ( aa != timeline.annotations.aliasing.end() )
    {
      writer.level( aa->first , globals::annot_strat );
      writer.value( "ORIG" , aa->second );
      ++aa;
    }
  writer.unlevel( globals::annot_strat );
  
  // channels
  std::map<std::string,std::string>::const_iterator cc = header.aliasing.begin();
  while ( cc != header.aliasing.end() )
    {
      writer.level( cc->first , globals::signal_strat );
      writer.value( "ORIG" , cc->second );
      ++cc;
    }
  writer.unlevel( globals::signal_strat );

}

void edf_t::terse_summary( param_t & param )
{

  // only non-annot signals here
  const bool NO_ANNOTS = true; 

  signal_list_t signals = header.signal_list( param.value( "sig" ) , NO_ANNOTS );
  
  const int ns1 = signals.size();

  const bool write_signals = param.has( "signals" );

  // write EDF type
  std::string edf_type = "EDF";
  if ( header.edfplus ) edf_type = header.continuous ? "EDF+C" : "EDF+D" ; 
  writer.value( "EDF_TYPE" , edf_type ); 

  // write output
  writer.value( "NS_ALL" , header.ns );
  writer.value( "NS" , ns1 );

  writer.value( "NR" , header.nr );
  writer.value( "REC_DUR" , header.record_duration );

  // total duration in TP units
  uint64_t duration_tp = globals::tp_1sec * (uint64_t)header.nr * header.record_duration ;
  std::string total_duration_hms = Helper::timestring( duration_tp , ':' , false );
  writer.value( "TOT_DUR_SEC" , header.nr * header.record_duration );
  writer.value( "TOT_DUR_HMS" , total_duration_hms );

  const std::string pat_id = Helper::trim( header.patient_id ) ;
  writer.value( "EDF_ID" , pat_id == "" ? "." : pat_id );
  writer.value( "START_TIME" , Helper::trim( header.starttime ) );
  writer.value( "START_DATE" , Helper::trim( header.startdate ) );

  // stop time
  clocktime_t et( header.starttime );
  if ( et.valid )
    {
      double time_sec = ( timeline.last_time_point_tp+1LLU) * globals::tp_duration ;
      et.advance_seconds( time_sec );
      writer.value( "STOP_TIME" , et.as_string() );
    }
  
  if ( write_signals ) 
    {
      std::vector<std::string> chs;
      for (int s=0;s<ns1;s++) chs.push_back( signals.label(s) );
      writer.value( "SIGNALS" , Helper::stringize<std::vector<std::string> >( chs ) ); 
      //writer.value( "SIGNALS" , Helper::stringize<std::vector<std::string> >( header.label ) ); 
    }
  

  //for (int s=0;s<header.ns;s++)
  for (int s1=0;s1<ns1;s1++)
    {

      const int s = signals(s1);
     
      // channel name
      writer.level( header.label[s] , globals::signal_strat );

      // channel type
      writer.value( "TYPE" , globals::map_channel_label(  header.label[s]  )  );
      
      // number of samples
      writer.value( "SR" , header.n_samples[s] / (double)header.record_duration );
      
      // physical dimension
      std::string pdim = Helper::trim( header.phys_dimension[s]  );
      writer.value( "PDIM" , pdim != "" ? pdim : "." );

      // transducer type
      std::string transtype = Helper::trim( header.transducer_type[s] );
      writer.value( "TRANS" , transtype != "" ? transtype : "." );

      // physical min/max
      writer.value( "PMIN" , header.physical_min[s] );
      writer.value( "PMAX" , header.physical_max[s] );
      
      // digital min/max
      writer.value( "DMIN" , header.digital_min[s] );
      writer.value( "DMAX" , header.digital_max[s] );

      // sensitivity (unit per bit)
      writer.value( "SENS" , ( header.physical_max[s] - header.physical_min[s] ) / (double)(  header.digital_max[s] -  header.digital_min[s] ) ); 

      // position in (in-memory) EDF
      writer.value( "POS" , s+1 );
            
    }
  
  writer.unlevel( globals::signal_strat );

}




std::set<int> edf_header_t::read( FILE * file , edfz_t * edfz , const std::set<std::string> * inp_signals )
{

  // must be *either* EDF or EDFZ
  
  if ( file != NULL && edfz != NULL ) 
    Helper::halt( "internal error in edf_header_t::read(), unclear whether EDF or EDFZ" );
  
  // Fixed buffer size for header
  // Total header = 256 + ns*256
 
  const int hdrSz = 256; 
  
  // Allocate space in the buffer for the header only
  byte_t * q = new byte_t[ hdrSz ];
  byte_t * q0 = q;

  //
  // Read start of header into the buffer
  //

  size_t rdsz;

  if ( file ) {
    rdsz = fread( q , 1, hdrSz , file);
  } else {
    rdsz = edfz->read( q , hdrSz );
  }
  
  std::set<int> channels;
  
  version        = edf_t::get_string( &q , 8 );
  patient_id     = edf_t::get_string( &q , 80 );
  recording_info = edf_t::get_string( &q , 80 );
  startdate      = edf_t::get_string( &q , 8 );
  starttime      = edf_t::get_string( &q , 8 );
  nbytes_header  = edf_t::get_int( &q , 8 );  
  reserved       = edf_t::get_bytes( &q , 44 );

  // enforce check that reserevd field contains only US-ASCII characters 32-126
  // not clear this is needed, but other software seems to prefer this

  Helper::ascii7( &reserved , ' ' );

  //
  // ensure starttime is in the PM, i.e. 07:00 --> 19:00
  // unless we've otherwise been instructed to respect
  // AM start-times (assume-pm-start=0);  but going to bed at midnight or 
  // 1am should be fine... so 

  //    6am ....  12pm   .... 6pm .... 12am .... 6am 
  //                    |4pm     
  
  //
  // assumes typical sleep onset 
  //
  
  if ( globals::assume_pm_starttime )
    {
      clocktime_t st( starttime );
      if ( st.valid ) 
	{
	  if      ( st.h >= globals::assume_pm_starttime_hour && st.h < 12 ) st.h += 12;
	  else if ( st.h == 12 ) st.h = 0; 
	  starttime = st.as_string();	  
	}
    }

  
  // EDF+C  continuous EDF 
  // EDF+D  discontinuous EDF+
  
  if (    reserved[0] == 'E' 
       && reserved[1] == 'D' 
       && reserved[2] == 'F'
       && reserved[3] == '+' )
    {
      
      if ( reserved[4] == 'C' )
	{
	  edfplus = true;
	  continuous = true;
	}
      else if ( reserved[4] == 'D' )
	{
	  edfplus = true;
	  continuous = false;
	}
    }
  else
    {
      edfplus = false;
      continuous = true;
    }
  
  // check whether we are forcing EDF format
  if ( globals::force_edf )
    {
      logger << "  forcing read as EDF [else remove force-edf=1]\n";

      edfplus = false;
      continuous = true;  
      reserved[0] = ' '; 
      reserved[1] = ' ';
      reserved[2] = ' ';
      reserved[3] = ' ';
      reserved[4] = ' ';

    }

    
  // Number and direction of records/signals
  
  nr                   = edf_t::get_int( &q , 8 );

  // store copy as 'original file' (i.e. if the edf_t is restructured, then
  // nr will be smaller, but we need the same nr_all for remaining file access

  nr_all               = nr; 

  record_duration      = edf_t::get_double( &q , 8 );

  record_duration_tp   = record_duration * globals::tp_1sec;

  ns_all               = edf_t::get_int( &q , 4 );


  // Free buffer
  delete [] q0;


  //
  // Per-signal header information
  //

  
  // read next 256 bytes per signal, i.e. overwriting existing buffer
  byte_t * p = new byte_t[ hdrSz * ns_all ]; 
  byte_t * p0 = p;

  if ( file ) 
    rdsz = fread( p , 1, hdrSz * ns_all , file);      
  else
    rdsz = edfz->read( p , hdrSz * ns_all );

  // for each of 'ns_all' signals
  
  ns = 0; // actual number of important signals

  std::vector<std::string> tlabels;
  std::set<std::string> slabels;
  
  for (int s=0;s<ns_all;s++)
    {
      // signal label, trim leading/trailing spaces
      std::string l = Helper::trim( edf_t::get_string( &p , 16 ) );

      // swap internal spaces? (not for special EDF Annotations channel)
      bool annotation = Helper::imatch( l , "EDF Annotation" , 14 ) ;
      if ( globals::replace_channel_spaces && ! annotation )
	l = Helper::search_replace( l , ' ' , globals::space_replacement );

      // global sanitization of channel labels?
      // but, if allowing spaces, then make these exempt
      // if either 'sanitize_everything' then retrim w/ underscore
      if ( globals::sanitize_everything && ! annotation )
	{
	  if ( globals::replace_channel_spaces )
	    l = Helper::trim( Helper::sanitize( l ) , '_' ) ;
	  else // allow spaces in a sanitized version still
	    l = Helper::trim( Helper::sanitize( l , ' ' ) , '_' ) ;
	}

      // make all data-channels upper case?
      if ( globals::uppercase_channels && ! annotation )
	l = Helper::toupper( l );

      // key on UC version
      std::string uc_l = Helper::toupper( l );
      
      // does this exist already? if so, uniqify 
      if ( slabels.find( uc_l ) != slabels.end() )
	{
	  int inc = 1;
	  while ( 1 ) 
	    {
	      // new unique label?
	      if ( slabels.find( uc_l + "." + Helper::int2str( inc )  ) == slabels.end() )
		{
		  logger << " uniquifying " << l ;
		  l = l + "." + Helper::int2str( inc );
		  uc_l = uc_l + "." + Helper::int2str( inc );
		  logger << " to " << l << "\n";
		  break;
		}
	      else // keep trying
		++inc;
	    }
	}
      
      // store temporary
      tlabels.push_back( l );
      slabels.insert( uc_l );

      // track original LABEL position
      label_all[ uc_l ] = s ;
      
    }

  
  // for each signal, does it match?
  // (and if so, change this to "standard" form)
  
  for (int s=0;s<ns_all;s++)
    {
      
      // retrieve temp label
      std::string l = tlabels[s];
      
      // this match function will change 'l' to match any primary aliase
      // it does a case-insensitive match, but returns the correct (preferred-case) version
      
      bool include = inp_signals == NULL || signal_list_t::match( inp_signals , &l , slabels );
      //std::cout << " l = " << include << " " << l << "\n";
	    
      // imatch allows for case-insensitive match of 'edf annotation*'  (i.e. 14 chars)
      bool annotation = Helper::imatch( l , "EDF Annotation" , 14 ) ;

      // optionally, skip all EDF annotation channels?
      // if this is EDF+C, we can just skip altogether;  otherwise,
      // we need to read the EDF+D time-track (but not other annots)
      
      if ( annotation )
	{	  
	  if ( globals::force_edf ) 
	    include = false;

	  // for EDF+D, will read time-tracks only
	  // for EDF+C, can skip the whole thing
	  if ( globals::skip_edf_annots && continuous )
	    include = false;
	}
    
      //
      // add this channel in 
      //

      
      if ( include ) 
	{

	  channels.insert(s);
	  
	  annotation_channel.push_back( annotation );
	  
	  if ( annotation && ! edfplus ) 
	    {
	      //Helper::halt( "file must be EDF+ to support annotations" );
	      logger << " detected an annotation channel in EDF: will treat as EDF+\n";
	      edfplus = true;
	    }
	  
	  // first annotation channel is time-track
	  if ( annotation && t_track == -1 ) 
	    t_track = label.size();


	  // label mapping only to non-annotation channels
	  if ( ! annotation ) 
	    label2header[ Helper::toupper( l ) ] = label.size(); 
	  
	  label.push_back( l );	  
	  
	  ++ns;

	}
    }

  // transducer type
  for (int s=0;s<ns_all;s++)
    {
      if ( channels.find(s) != channels.end() ) 
	transducer_type.push_back( Helper::trim( edf_t::get_string( &p , 80 ) ) );
      else 
	edf_t::skip( &p , 80 );
    }

  // physical dimension
  for (int s=0;s<ns_all;s++)
    {
      if ( channels.find(s) != channels.end() ) 
	phys_dimension.push_back( Helper::trim( edf_t::get_string( &p , 8 ) ) );
      else
	edf_t::skip( &p , 8 );
    }

  // physical min
  for (int s=0;s<ns_all;s++)
    {
      if ( channels.find(s) != channels.end() ) 
	physical_min.push_back( edf_t::get_double( &p , 8 ) );
      else
	edf_t::skip( &p , 8 );
    }

  // physical max
  for (int s=0;s<ns_all;s++)
    {
      if ( channels.find(s) != channels.end() )       
	physical_max.push_back( edf_t::get_double( &p , 8 ) );
      else
	edf_t::skip( &p , 8 );
    }

  // digital min  
  for (int s=0;s<ns_all;s++)
    {
      if ( channels.find(s) != channels.end() )       
	digital_min.push_back( edf_t::get_int( &p , 8 ) );
      else
	edf_t::skip( &p , 8 );      
    }

  // digital max
  for (int s=0;s<ns_all;s++)
    {
      if ( channels.find(s) != channels.end() )
	digital_max.push_back( edf_t::get_int( &p , 8 ) );
      else
	edf_t::skip( &p , 8 );      
    }


  // prefiltering information
  for (int s=0;s<ns_all;s++)
    {
      if ( channels.find(s) != channels.end() )
	prefiltering.push_back( edf_t::get_string( &p , 80 ) );
      else
	edf_t::skip( &p , 80 );
    }
  
  // number of samples per record
  for (int s=0;s<ns_all;s++)
    {
      int x = edf_t::get_int( &p , 8 );

      // SR == 0 Hz ? 
      if ( x == 0 )
	logger << "  *** warning, " << s << " has SR of 0 and should be dropped\n";
      
      // non-integer SR ? 
      double sr = x / record_duration;
      if ( fabs( trunc(sr) - sr ) > 1e-8 )
	logger << "  *** warning, signal " << s << " has a non-integer SR - advise to RESAMPLE\n";
      
      if ( channels.find(s) != channels.end() )
	n_samples.push_back( x );      
      n_samples_all.push_back( x );
    }
  
  // reserved field
  for (int s=0;s<ns_all;s++)
    {
      if ( channels.find(s) != channels.end() )
	signal_reserved.push_back( edf_t::get_string( &p , 32 ) );
      else
	edf_t::skip( &p , 32 );
    }

  //
  // time-track absolute offset in record (we only care about this 
  //  when reading from disk)
  //

  if ( t_track != -1 )
    {
      t_track_edf_offset = 0;
      for (int ss=0;ss<t_track;ss++)
	t_track_edf_offset += 2 * n_samples_all[ss];	
    }

  //
  // derived values: note, here 'ns' not 'ns_all'
  //

  orig_physical_min = physical_min;
  orig_physical_max = physical_max;

  orig_digital_min = digital_min;
  orig_digital_max = digital_max;
  
  for (int s=0;s<ns;s++)
    {
      double bv = ( physical_max[s] - physical_min[s] ) / (double)( digital_max[s] - digital_min[s] ) ;
      bitvalue.push_back( bv );
      offset.push_back( ( physical_max[s] / bv ) - digital_max[s] ) ;
    }  
  
  
  // clean up buffer
  delete [] p0 ;


  // return mapping of imported channel numbers
  return channels;
  
}




bool edf_record_t::read( int r )
{
  
  // bound checking on 'r' already done, via edf_t::read_record();
  
  // skip if already loaded?
  if ( edf->loaded( r ) ) return false;
  
  // allocate space in the buffer for a single record, and read from file
  
  byte_t * p = new byte_t[ edf->record_size ];
  
  byte_t * p0 = p;

  // EDF?
  if ( edf->file ) 
    {
      
      // determine offset into EDF
      uint64_t offset = edf->header_size + (uint64_t)(edf->record_size) * r;
      
      // find the appropriate record
      fseek( edf->file , offset , SEEK_SET );
      
      // and read it
      size_t rdsz = fread( p , 1, edf->record_size , edf->file );
    }
  else // EDFZ
    {
      
      if ( ! edf->edfz->read_record( r , p , edf->record_size ) ) 
	Helper::halt( "corrupt .edfz or .idx" );      

    }

  // which signals/channels do we actually want to read?
  // header : 0..(ns-1)
  // from record data : 0..(ns_all-1), from which we pick the 'ns' entries is 'channels'
  // data[] is already created for 'ns' signals
  
  // for convenience, use name 'channels' below
  std::set<int> & channels = edf->inp_signals_n;
    
  int s = 0;

  for (int s0=0; s0<edf->header.ns_all; s0++)
    {

      // need to EDF-based header, i.e. if skipped signal still need size to skip
      const int nsamples = edf->header.n_samples_all[s0];
      
      //
      // skip this signal?
      //
      
      if ( channels.find( s0 ) == channels.end() )
	{	  
	  p += 2 * nsamples;
	  continue;
	}
      

      //
      // Data or annotation channel? (note: lookup is based on 's' not
      // 's0', i.w. loaded channels, not all EDF channels
      //
      
      bool annotation = edf->header.is_annotation_channel( s );

      //
      // s0 : actual signal in EDF
      // s  : where this signal will land in edf_t
      //
      
      if ( ! annotation ) 
	{
	  
	  for (int j=0; j < nsamples ; j++)
	    {

	      //int d = tc2dec( **p ,  *((*p)+1)  ); 
	      int16_t d = tc2dec( *p ,  *(p+1)  ); 
	      
	      // advance pointer
	      p += 2;

	      // store digital data-point
	      data[s][j] = d;
	      
	      // physically-scaled data-point	  
	      if ( false )
		if ( d < edf->header.orig_digital_min[s] || d > edf->header.orig_digital_max[s] ) 
		  {	
		    
		  std::cout << "OUT-OF-BOUNDS" << "\t"
			    << edf->id << "\t"
			    << "[" << globals::current_tag << "]\t"
			    << edf->header.label[s] << "\t"
			    << "digt: " << d << "\t"
			    << edf->header.orig_digital_min[s] << " .. " 
			    << edf->header.orig_digital_max[s] << "\t"
			    << "phys: " << edf->header.bitvalue[s] * ( edf->header.offset[s] + d ) << "\t"
			    << edf->header.orig_physical_min[s] << " .. " 
			    << edf->header.orig_physical_max[s] << "\n"; 

		  if (  d < edf->header.orig_digital_min[s] ) 
		    d = edf->header.orig_digital_min[s];
		  else 
		    d = edf->header.orig_digital_max[s];
		  
		  }
	      
	      // concert to physical scale
	      //pdata[s][j] = edf->header.bitvalue[s] * ( edf->header.offset[s] + d );
	      //pdata[s][j] = dig2phys( d , s ) ;
	      
	    }
	}
      else // read as a ANNOTATION
	{
	  
	  // Note, because for a normal signal, each sample takes 2 bytes,
	  // here we read twice the number of datapoints
	  
	  for (int j=0; j < 2 * nsamples; j++)
	    {
	      
	      // store digital data-point
	      data[s][j] = *p;
	      
	      // advance pointer
	      p++;
	      
	    }	  
	  
	}
      

      // next signal

      ++s;

    }


  //
  // Clean up
  //

  delete [] p0;
  
  return true;

}



bool edf_t::read_records( int r1 , int r2 )
{

  // This only tries to load records that are 'retained' and 
  // not already in memory

  if ( r1 < 0 ) r1 = 0;
  if ( r1 > header.nr_all ) r1 = header.nr_all - 1;
  
  if ( r2 < r1 ) r2 = r1;
  if ( r2 > header.nr_all ) r2 = header.nr_all - 1;

  //  std::cerr << "edf_t::read_records :: scanning ... r1, r2 " << r1 << "\t" << r2 << "\n";
  
  for (int r=r1;r<=r2;r++)
    {

      // if ( ! timeline.retained(r) ) 
      // 	std::cerr << "NOT retained " << r << " " << timeline.retained(r) << "\n";
      
      if ( timeline.retained(r) )
	{
	  if ( ! loaded( r ) ) 
	    {
	      edf_record_t record( this ); 
	      record.read( r );
	      records.insert( std::map<int,edf_record_t>::value_type( r , record ) );	      
	    }
	}
    }
  return true;
}


bool edf_t::init_empty( const std::string & i ,
			const int nr ,
			const int rs ,
			const std::string & startdate ,
			const std::string & starttime )
{

  if ( nr == 0 || rs == 0 ) return false;

  id = i;
  
  //
  // Set header
  //

  header.version = "0";
  header.patient_id = id;
  header.recording_info = "";
  header.startdate = startdate;
  header.starttime = starttime;
  header.nbytes_header = 256 + 0 * 256;  // i.e. no signals
  header.ns = 0; // these will be added by add_signal()
  header.ns_all = 0; // check this... should only matter for EDF access, so okay... 
  header.nr = header.nr_all = nr;  // likewise, value of nr_all should not matter, but set anyway
  header.record_duration = rs;
  header.record_duration_tp = header.record_duration * globals::tp_1sec;
  
  //
  // create a timeline
  //

  set_edf();

  set_continuous();

  timeline.init_timeline();


  //
  // resize data[][], by adding empty records
  //

  for (int r=0;r<nr;r++)
    {
      edf_record_t record( this ); 
      records.insert( std::map<int,edf_record_t>::value_type( r , record ) );
    }

  logger << "  created an empty EDF of duration " << rs * nr << " seconds\n";
  
  return true;

}


bool edf_t::read_from_ascii( const std::string & f , // filename
			     const std::string & i , // id
			     const int Fs , // fixed Fs for all signals
			     const std::vector<std::string> & labels0 ,  // if null, look for header
			     const std::string & startdate , 
			     const std::string & starttime 
			     ) 
{
  
  filename = Helper::expand( f );
  
  id = i;
  
  bool has_arg_labels = labels0.size() > 0;

  bool has_header_labels = false;

  std::vector<std::string> labels;

  if ( has_arg_labels ) labels = labels0;

  if ( ! Helper::fileExists( filename ) ) 
    Helper::halt( "could not read " + filename );

  bool compressed = Helper::file_extension( filename , "gz" );
  
  std::ifstream IN1( filename.c_str() , std::ios::in );

  gzifstream ZIN1;
  
  if ( compressed  ) 
    ZIN1.open( filename.c_str() );
  else
    IN1.open( filename.c_str() , std::ios::in );    
  
  std::string line;

  if ( compressed ) 
    {
      Helper::zsafe_getline( ZIN1 , line );
      if ( ZIN1.eof() || line == "" ) Helper::halt( "problem reading from " + filename + ", empty?" );
    }
  else
    {
      Helper::safe_getline( IN1 , line );
      if ( IN1.eof() || line == "" ) Helper::halt( "problem reading from " + filename + ", empty?" );
    }
  
      
  // has a header row (whether we want to use it or not)

  if ( line[0] == '#' ) 
    {
      has_header_labels = true;
      if ( has_arg_labels ) 
	logger << "  ignoring header row in " << filename << " as channel labels specified with --chs\n" ;
      else
	{
	  line = line.substr(1);
	  labels = Helper::parse( line , "\t ," );
	}
    } 

  // if no arg or header labels, we need to make something up 
  else if ( ! has_arg_labels ) 
    {
      std::vector<std::string> tok = Helper::parse( line , "\t ," );
      labels.resize( tok.size() );
      for (int l=0;l<labels.size();l++) labels[l] = "S" + Helper::int2str(l+1) ;
    }

  // and rewind file to start from the beginning, if needed
  
  if ( ! has_header_labels ) 
    {
      if ( compressed )
	{
	  ZIN1.clear();
	  ZIN1.seekg(0, std::ios::beg);
	}
      else
	{
	  IN1.clear();
	  IN1.seekg(0, std::ios::beg);
	}
    }

  
  const int ns = labels.size();

  //
  // Scan file to get number of records
  //
  
  int np = 0;
  while ( !IN1.eof() ) 
    {
      std::string line;

      if ( compressed )
	{	  
	  Helper::zsafe_getline( ZIN1 , line );	  
	  if ( ZIN1.eof() ) break;
	}
      else
	{
	  Helper::safe_getline( IN1 , line );
	  if ( IN1.eof() ) break;
	}
      
	  
      if ( line == "" ) continue;
      ++np;
    }

  // will ignore any partial records at the end of the file
  int nr = np / Fs;
  np = nr * Fs;

  if ( compressed ) 
    IN1.close();
  else
    IN1.close();

  // re-read
  std::ifstream IN2;
  gzifstream ZIN2;
  
  if ( compressed ) 
    ZIN2.open( filename.c_str() );
  else
    IN2.open( filename.c_str() , std::ios::in );

  // skip header?
  if ( has_header_labels ) 
    {
      std::string dummy;
      
      if ( compressed )	
	std::getline( ZIN2 , dummy );
      else
	Helper::safe_getline( IN2 , dummy );	
    }
  

  //
  // Set header
  //

  header.version = "0";
  header.patient_id = id;
  header.recording_info = "";
  header.startdate = startdate;
  header.starttime = starttime;
  header.nbytes_header = 256 + ns * 256;
  header.ns = 0; // these will be added by add_signal()
  header.ns_all = ns; // check this... should only matter for EDF access, so okay... 
  header.nr = header.nr_all = nr;  // likewise, value of nr_all should not matter, but set anyway
  header.record_duration = 1;
  header.record_duration_tp = header.record_duration * globals::tp_1sec;

  

  //
  // create a timeline
  //

  set_edf();

  set_continuous();

  timeline.init_timeline();

  //
  // read data
  //

  logger << "  reading " << ns << " signals, " 
	 << nr << " seconds ("
	 << np << " samples " << Fs << " Hz) from " << filename << "\n";

  Data::Matrix<double> data( np , ns );
  
  for (int p=0;p<np;p++)
    for (int s=0;s<ns;s++)
      {

	if ( compressed )
	  ZIN2 >> data(p,s);
	else
	  IN2 >> data(p,s);
	
	if ( IN2.eof() ) 
	  Helper::halt( filename + " does not contain enough data-points given parameters\n" );
      }
  
  double dd;

  if ( compressed )
    {
      ZIN2 >> dd;
      if ( ! ZIN2.eof() ) 
	logger << " ** warning, truncating potential trailing sample points (<1 second) from end of input\n";
    }
  else
    {
      IN2 >> dd;
      if ( ! IN2.eof() ) 
	logger << " ** warning, truncating potential trailing sample points (<1 second) from end of input\n";
    }
  
  // should now be end of file...

  if ( compressed )
    ZIN2.close();
  else
    IN2.close();


  //
  // resize data[][], by adding empty records
  //

  for (int r=0;r<nr;r++)
    {
      edf_record_t record( this ); 
      records.insert( std::map<int,edf_record_t>::value_type( r , record ) );
    }

  //
  // add signals (this populates channel-specific 
  //
  
  for (int s=0;s<ns;s++)
    add_signal( labels[s] , Fs , *data.col(s).data_pointer() );

  return true;
}



bool edf_t::attach( const std::string & f , 
		    const std::string & i , 
		    const std::set<std::string> * inp_signals ,
		    const bool silent )
{
  
  //
  // Store filename and ID
  //
  
  // expand() expands out any ~/ notation to full path
  filename = Helper::expand( f ) ;

  id = i; 

  //
  // EDF or EDFZ?
  //
  
  file = NULL; 

  edfz = NULL;

  bool edfz_mode = Helper::file_extension( filename , "edfz" )
    || Helper::file_extension( filename , "edf.gz" ); 

  //
  // Attach the file
  //
  
  if ( ! edfz_mode ) 
    {
      if ( ( file = fopen( filename.c_str() , "rb" ) ) == NULL )
	{      
	  file = NULL;
	  logger << " PROBLEM: could not open specified EDF: " << filename << "\n";
	  globals::problem = true;
	  return false;
	}
    }
  else
    {
      edfz = new edfz_t;
      
      // this also looks for the .idx, which sets the record size
      if ( ! edfz->open_for_reading( filename ) ) 
	{
	  delete edfz;
	  edfz = NULL;
	  logger << " PROBLEM: could not open specified .edfz (or .edfz.idx) " << filename << "\n";
	  globals::problem = true;
	  return false;
	}

    }

  
  //
  // Does this look like a valid EDF (i.e. at least contains a header?)
  //

  uint64_t fileSize = 0 ; 

  // for EDF
  if ( file ) 
    {

      fileSize = edf_t::get_filesize( file );
      
      if ( fileSize < 256 ) 
	{
	  logger << " PROBLEM: corrupt EDF, file < header size (256 bytes): " << filename << "\n";
	  globals::problem = true;
	  return false;
	}
    }
  else
    {
      // TODO... need to check EDFZ file. e.g. try reading the last record?
      //
    }

  //
  // Read and parse the EDF header (from either EDF or EDFZ)
  //
  
  // Parse the header and extract signal codes 
  // store so we know how to read records
  
  inp_signals_n = header.read( file , edfz , inp_signals );
  

  //
  // anon header info?
  //

  if ( globals::anon )
    {
      // ID, recording info and startdate --> NULL 
      header.patient_id = header.edfplus ? "X X X X" : ".";
      header.recording_info = header.edfplus ? "Startdate X X X X" : ".";
      header.startdate = "01.01.85";
    }

  //
  // force EDF start times/dates
  //

  if ( globals::force_starttime != "" )
    {
      header.starttime = globals::force_starttime;
      logger << "  forced start-time to " << header.starttime << "\n";
    }
  
  if ( globals::force_startdate != "" )
    {
      header.startdate = globals::force_startdate;
      logger << "  forced start-date to " << header.startdate << "\n";
    }

  //
  // Swap out any signal label aliases at this point
  //
  
  swap_in_aliases();
  
  //
  // EDF+ requires a time-track
  //

  if ( header.edfplus && header.time_track() == -1 ) 
    {
      if ( !header.continuous ) 
	Helper::halt( "EDF+D with no time track" );

      logger << " EDF+ [" << filename << "] did not contain any time-track: adding...\n";

      add_time_track();

    }

  
  //
  // Record details about byte-size of header/records
  //
  
  header_size = 256 + header.ns_all * 256;
  
  record_size = 0;
  
  for (int s=0;s<header.ns_all;s++)
    record_size += 2 * header.n_samples_all[s] ; // 2 bytes each

  
  if ( edfz ) 
    {
      if ( record_size != edfz->record_size )
	{
	  logger << "  EDFZ idx record size = " << edfz->record_size << "\n"
		 << "  EDF record size = " << record_size << "\n";
	  Helper::halt( "internal error, different record size in EDFZ header versus index" );
	}
    }
  


  //
  // Check remaining file size, based on header information
  //    
  
  if ( file ) 
    {
      uint64_t implied = (uint64_t)header_size + (uint64_t)header.nr_all * record_size;
      
      if ( fileSize != implied ) 
	{

	  std::stringstream msg;

	  msg << "details:\n"
	      << "  header size ( = 256 + # signals * 256 ) = " << header_size << "\n"
	      << "  num signals = " << header.ns_all << "\n"	      
	      << "  record size = " << record_size << "\n"
	      << "  number of records = " << header.nr_all << "\n"
	      << "  implied EDF size from header = "
	      << header_size << " + " << record_size << " * " << header.nr_all << " = " << implied << "\n\n"

	      << "  assuming header correct, implies the file has "
	      <<  (double)(fileSize-header_size)/(double)record_size - (double)(implied-header_size)/(double)record_size 
	      << " records too many\n"
	      << "  (where one record is " << header.record_duration << " seconds)\n";
	  
	  if ( ! globals::autofix_edf )
	    {
	      msg << "\nIF you're confident about the remaining data you can add the option:\n\n"
		  << "    luna s.lst fix-edf=T ... \n\n"
		  << "  to attempt to fix this.  This may be appropriate under some circumstances, e.g.\n"
		  << "  if just the last one or two records were clipped.  However, if other EDF header\n"
		  << "  information is incorrect (e.g. number of signals, sample rates), then you'll be\n"
		  << "  dealing with GIGO... so be sure to carefully check all signals for expected properties;\n"
		  << "  really you should try to determine why the EDF was invalid in the first instance, though\n";

	      Helper::halt( "corrupt EDF: expecting " + Helper::int2str(implied) 
			    + " but observed " + Helper::int2str( fileSize) + " bytes" + "\n" + msg.str() );
	    }
	  else
	    {
	      logger << "  warning: EDF has incorrect file size given header information:\n"
		     << msg.str() << "\n";
	      
	      int nr_from_data = floor( (fileSize-header_size)/(double)record_size );
	      
	      logger << "  attempting to fix this, changing the header number of records from " << header.nr_all
		     << " to " << nr_from_data << " ... good luck!\n";
	      
	      // update EDF header internally.
	      header.nr_all = header.nr = nr_from_data;
	      	      
	    }

	  
	}
    }
  
  //
  // Create timeline (relates time-points to records and vice-versa)
  // Here we assume a continuous EDF, but timeline is set up so that 
  // this need not be the case
  //

  timeline.init_timeline();


  //
  // Output some basic information
  //

  if ( ! silent ) 
    {
      
      logger << " duration " << Helper::timestring( timeline.total_duration_tp , '.' , false )  // not fractional
	     << ", " << timeline.total_duration_tp * globals::tp_duration << "s";
      
      clocktime_t et( header.starttime );
      
      if ( et.valid )
	{
	  // nb. going to one past end:
	  double time_sec = ( (timeline.last_time_point_tp+1LLU) * globals::tp_duration ) ;
	  et.advance_seconds( time_sec );
	  logger << " | time " << header.starttime << " - " << et.as_string() ;
	}

      logger << " | date " << header.startdate;
      logger << "\n";
      
      //	 << " hms, last time-point " << Helper::timestring( ++timeline.last_time_point_tp ) << " hms after start\n";
      
      if ( globals::verbose )
	logger << "  " << header.nr_all  << " records, each of " << header.record_duration << " second(s)\n";
      
      logger << "\n signals: " << header.ns << " (of " << header.ns_all << ") selected ";
      
      if ( header.edfplus & header.continuous ) 
	logger << "in an EDF+C file" ;
      else if ( header.edfplus & ! header.continuous ) 
	logger << "in an EDF+D file" ;
      else
	logger << "in a standard EDF file" ;

      for (int s=0;s<header.ns;s++) 
	logger << ( s % 8 == 0 ? "\n  " : " | " ) << header.label[s]; 
      logger << "\n";
    }
  
  return true;
  
}



void edf_t::swap_in_aliases()
{

  // simply get a wildcard-ed signal_list_t
  // as this process of searching for all signals also 
  // swaps in the alias and updates the EDF header
  
  signal_list_t dummy = header.signal_list( "*" );
  
}


std::vector<double> edf_t::fixedrate_signal( uint64_t start , 
					     uint64_t stop , 
					     const int signal , 
					     const int downsample ,
					     std::vector<uint64_t> * tp , 
					     std::vector<int> * rec , 
					     std::vector<int16_t> * ddata ) 
{
  
  std::vector<double> ret;

  if ( tp != NULL ) 
    tp->clear();

  if ( rec != NULL ) 
    rec->clear();

  if ( ddata != NULL )
    ddata->clear();

  //
  // Ensure we are within bounds
  //
  
  if ( stop > timeline.last_time_point_tp + 1 )
    stop = timeline.last_time_point_tp + 1 ;      
  
  //
  // First, determine which records are being requested?
  //
  
  const uint64_t n_samples_per_record = header.n_samples[signal];
  
  // std::cerr << "signal = " << signal << "\t" << header.n_samples.size() << "\t" << header.n_samples_all.size() << "\n";
  // std::cerr << "SR " << n_samples_per_record << "\n";
 
  int start_record, stop_record;
  int start_sample, stop_sample;

  //  std::cerr << "looking for " << start << " to " << stop << "\n";
  
  bool okay = timeline.interval2records( interval_t( start , stop ) , 
					 n_samples_per_record , 
					 &start_record, &start_sample , 
					 &stop_record, &stop_sample );
  
  
  // std::cerr << "records start = " << start_record << " .. " << start_sample << "\n";
  // std::cerr << "records stop  = " << stop_record << " .. " << stop_sample << "\n";
  
  //
  // If the interval is too small (or is applied to a signal with a low sampling rate)
  // we might not find any sample-points in this region.   Not an error per se, but flag
  // (we have to check that all downstream functons will play nicely with an empty set being returned)
  //
  
  if ( ! okay ) 
    {
      logger << " ** warning ... empty intervals returned (check intervals/sampling rates)\n";
      return ret; // i.e. empty
    }

  
  //
  // Ensure that these records are loaded into memory
  // (if they are already, they will not be re-read)
  //

  bool retval = read_records( start_record , stop_record );
  
  //
  // Copy data into a single vector
  //
 
  double bitvalue = header.bitvalue[ signal ];
  double offset   = header.offset[ signal ];

  int r = start_record;

  while ( r <= stop_record )
    {
      //std::cerr << "rec " << r << "\n";

      // std::cout << records.size() << " is REC SIZE\n";
      // std::cout << "foudn " << ( records.find( r ) != records.end() ? " FOUND " : "NOWHERE" ) << "\n";

      const edf_record_t * record = &(records.find( r )->second);

      //std::cerr << " test for NULL " << ( record == NULL ? "NULL" : "OK" ) << "\n";
      
      const int start = r == start_record ? start_sample : 0 ;
      const int stop  = r == stop_record  ? stop_sample  : n_samples_per_record - 1;

      // std::cerr << " start, stop = " << start << "   " << stop << "\n";
      
      // std::cerr << "OUT\t"
      // 		<< record->data.size() << " "
      // 		<< signal << " " 
      // 		<< header.ns << "\n";
      
      for (int s=start;s<=stop;s+=downsample)
	{
	  
	  if ( tp != NULL ) 
	    tp->push_back( timeline.timepoint( r , s , n_samples_per_record ) );
	  if ( rec != NULL ) 
	    rec->push_back( r );

	  // just return digital values...
	  if ( ddata != NULL )
	    ddata->push_back( record->data[ signal ][ s ] );
	  else // ... or convert from digital to physical on-the-fly? (the default)
	    ret.push_back( edf_record_t::dig2phys( record->data[ signal ][ s ] , bitvalue , offset ) );

	}
      
      r = timeline.next_record(r);
      if ( r == -1 ) break;
    }
  
  // will be length==0 if digital == T 
  return ret;  
}



//
// Functions to write an EDF
//


bool edf_header_t::write( FILE * file , const std::vector<int> & ch2slot )
{

  // new number of channels (might be less than original)

  const int ns2 = ch2slot.size();
  
  // regarding the nbytes_header variable, although we don't really
  // use it, still ensure that it is properly set (i.e. we may have
  // added/removed signals, so we need to update before making the EDF)
  
  const int nbytes_header2 = 256 + ns2 * 256;
  
  writestring( version , 8 , file );
  writestring( patient_id , 80 , file );
  writestring( recording_info , 80 , file );
  writestring( startdate , 8 , file );
  writestring( starttime , 8 , file );
  writestring( nbytes_header2 , 8 , file );
  fwrite( reserved.data() , 1 , 44 , file );
  writestring( nr , 8 , file );
  writestring( record_duration , 8 , file );
  writestring( ns2 , 4 , file );

  // for each of 'ns2' signals
  
  for (int s=0;s<ns2;s++)
    writestring( label[ ch2slot[s] ], 16, file );
  
  for (int s=0;s<ns2;s++)
    writestring( transducer_type[ch2slot[s]], 80, file );

  for (int s=0;s<ns2;s++)
    writestring( phys_dimension[ch2slot[s]], 8, file );

  for (int s=0;s<ns2;s++)
    writestring( physical_min[ch2slot[s]], 8, file );

  for (int s=0;s<ns2;s++)
    writestring( physical_max[ch2slot[s]], 8, file );

  for (int s=0;s<ns2;s++)
    writestring( digital_min[ch2slot[s]], 8, file );

  for (int s=0;s<ns2;s++)
    writestring( digital_max[ch2slot[s]], 8, file );

  for (int s=0;s<ns2;s++)
    writestring( prefiltering[ch2slot[s]], 80, file );

  for (int s=0;s<ns2;s++)
    writestring( n_samples[ch2slot[s]], 8, file );
  
  for (int s=0;s<ns2;s++)
    writestring( signal_reserved[ch2slot[s]], 32, file );
  
  return true;
}



bool edf_header_t::write( edfz_t * edfz , const std::vector<int> & ch2slot )
{

  // new number of channels (might be less than original)

  const int ns2 = ch2slot.size();
  
  // regarding the nbytes_header variable, although we don't really
  // use it, still ensure that it is properly set (i.e. we may have
  // added/removed signals, so we need to update before making the EDF
  const int nbytes_header2 = 256 + ns2 * 256;
  
  edfz->writestring( version , 8 );
  edfz->writestring( patient_id , 80 );
  edfz->writestring( recording_info , 80 );
  edfz->writestring( startdate , 8 );
  edfz->writestring( starttime , 8 );
  edfz->writestring( nbytes_header2 , 8 );
  edfz->write( (byte_t*)reserved.data() , 44 );
  edfz->writestring( nr , 8 );
  edfz->writestring( record_duration , 8 );
  edfz->writestring( ns2 , 4 ); // ns2

  // for each of 'ns2' signals
  
  for (int s=0;s<ns2;s++)
    edfz->writestring( label[ch2slot[s]], 16 );
  
  for (int s=0;s<ns2;s++)
    edfz->writestring( transducer_type[ch2slot[s]], 80 );

  for (int s=0;s<ns2;s++)
    edfz->writestring( phys_dimension[ch2slot[s]], 8 );

  for (int s=0;s<ns2;s++)
    edfz->writestring( physical_min[ch2slot[s]], 8 );

  for (int s=0;s<ns2;s++)
    edfz->writestring( physical_max[ch2slot[s]], 8 );

  for (int s=0;s<ns2;s++)
    edfz->writestring( digital_min[ch2slot[s]], 8  );

  for (int s=0;s<ns2;s++)
    edfz->writestring( digital_max[ch2slot[s]], 8 );

  for (int s=0;s<ns2;s++)
    edfz->writestring( prefiltering[ch2slot[s]], 80 );

  for (int s=0;s<ns2;s++)
    edfz->writestring( n_samples[ch2slot[s]], 8 );
  
  for (int s=0;s<ns2;s++)
    edfz->writestring( signal_reserved[ch2slot[s]], 32 );
  
  return true;
}





bool edf_record_t::write( FILE * file , const std::vector<int> & ch2slot )
{

  const int ns2 = ch2slot.size();

  for (int s2=0; s2 < ns2; s2++)
    {
      // get actual from-slot 
      const int s = ch2slot[s2];

      const int nsamples = edf->header.n_samples[s];

      //
      // Normal data channel
      //

      if ( edf->header.is_data_channel(s) )
	{
	  for (int j=0;j<nsamples;j++)
	    {	  
	      char a , b;
	      dec2tc( data[s][j] , &a, &b );	  
	      fputc( a , file );
	      fputc( b , file );
	    }
	}
      
      //
      // EDF Annotations channel
      //
      
      if ( edf->header.is_annotation_channel(s) )
	{
	  //std::cout << " ANNOT WRT ";
	  for (int j=0;j< 2*nsamples;j++)
	    {	  	      
	      char a = j >= data[s].size() ? '\x00' : data[s][j];	      
	      //std::cout << a ;
	      fputc( a , file );
	    }
	  //std::cout << "\n";
	}
    
    }

  return true;
}


bool edf_record_t::write( edfz_t * edfz , const std::vector<int> & ch2slot )
{
  
  const int ns2 = ch2slot.size();
  
  for (int s2=0; s2 < ns2; s2++)
    {
      
      const int s = ch2slot[s2];

      const int nsamples = edf->header.n_samples[s];

      //
      // Normal data channel
      //

      if ( edf->header.is_data_channel(s) )
	{  
	  std::vector<char> d( 2 * nsamples );
	  
	  for (int j=0;j<nsamples;j++)
	    dec2tc( data[s][j] , &(d)[2*j], &(d)[2*j+1] );	  

	  edfz->write( (byte_t*)&(d)[0] , 2 * nsamples );
	  
	}
      
      //
      // EDF Annotations channel
      //
      
      if ( edf->header.is_annotation_channel(s) )
	{      	  	  

	  std::vector<char> d( 2 * nsamples );

	  for (int j=0;j< 2*nsamples;j++)
	    {	  	      
	      char a = j >= data[s].size() ? '\x00' : data[s][j];	      
	      d[j] = a;
	    }
	  
	  edfz->write( (byte_t*)&(d)[0] , 2 * nsamples ); 
	  
	}
    
    }

  return true;
}


bool  edf_t::is_actually_standard_edf()
{
  if ( ! header.edfplus ) return true;
  
  // EDF Annotations (other than time track)?
  if ( has_edf_annots ) return false;

  // for (int s=0;s<header.ns;s++)
  //   {
  //     if ( ! header.is_data_channel(s) )
  // 	if ( s != header.t_track ) return false;
  //   }

  // discontinuous?
  if ( is_actually_discontinuous() ) return false;
  
  return true;
}


bool edf_t::is_actually_discontinuous() 
{

  // definitely continuous
  if ( header.continuous ) return false;
  
  // otherwise, check whether any gaps actually present 
  // (i.e. versus start/end missing, which Luna will still treat
  // as 'discontinuous' for internal reasons
  
  int num_segments = 0;

  int r = timeline.first_record();
  
  uint64_t tp0 = timeline.rec2tp[r];

  uint64_t tp_start = tp0;  

  while ( r != -1 )
    {
      
      // next record
      r = timeline.next_record( r );
      
      // start of this next record
      uint64_t tp;

      bool segend = false;
      
      // end?
      if ( r == -1 )
	{
	  // make this the 'previous'
	  tp0 = tp;
	  segend = true;
	}
      else
	{
	  tp = timeline.rec2tp[r] ;
	  // discontinuity / end of segment?
	  segend = tp - tp0 != header.record_duration_tp ;
	}
      
      // record this segment 
	   
      if ( segend )
	{
	  ++num_segments ;	  
	  // current point becomes start of the next segment
	  tp_start = tp;	  
	}      
      // current point becomes the last one, for next lookup
      tp0 = tp;
    }
  
  // is this discontinuous?
  return num_segments > 1;
      
}


bool edf_t::write( const std::string & f , bool as_edfz , int write_as_edf , bool always_edfd , const std::vector<int> * p_ch2slot )
{

  // write_as_edf 0   -- no, do not force as EDF
  //              1   -- yes, force as EDF but reset start time
  //              2   -- yes, force as EDF and set starttime to NULL (00.00.00) w/ message
  
  //
  // Is this EDF+ truly discontinuous?  i.e. a discontinuous flag is set after any RESTRUCTURE
  // We'll keep this as is here, but for the purpose of WRITE-ing an EDF+ (only), we'll first check 
  // whether it is truly discontinuous (i.e. versus only the start/end was removed).
  //  If always_edfd == T, do not make this change from EDF+D --> EDF+C when writing however.
  //  i.e. this may be desirable if we want offset times always from the EDF+D start
  //

  bool actually_EDFD = is_actually_discontinuous();
  
  bool make_EDFC = (!always_edfd) && (!header.continuous) && (! actually_EDFD ); 

  bool actually_EDF = is_actually_standard_edf();

  if ( actually_EDF && actually_EDFD )
    Helper::halt( "internal error in write() when determining EDF type" );

  if ( actually_EDFD )
    logger << "  data are truly discontinuous\n";
  else
    logger << "  data are not truly discontinuous\n";
  
  //
  // Reset start-time to NULL (i.e. to writing as standard EDF but is actually discontinuous, then)
  // clocktimes will not make sense
  //
  
  bool null_starttime = write_as_edf == 2 && actually_EDFD; 
  
  //
  // Force as standard EDF? 
  //

  if ( ( write_as_edf || actually_EDF ) && ! always_edfd ) 
    {
      logger << "  writing as a standard EDF\n";
      set_edf();
    }
  
  //
  // Deal with start time?  If writing as a truly discontinuous EDF+D, then
  // keep start-time as is (i.e. first epoch might not be 0 seconds).  But if
  // writing as a EDF+C, then make start time == first record (or even if EDF
  // unless we are told otherwise) 
  //
    
  if ( null_starttime ) 
    {
      logger << "  setting EDF starttime to null (00.00.00)\n";
      header.starttime = "00.00.00";
    }
  else if ( write_as_edf == 1 || make_EDFC )  // no changes for EDF+D
    {      
      reset_start_time();
    }
  
  
  //
  // By default, ch2slot will be 0,1,2,...,ns-1   i.e. a straight mapping/ordering of all channels
  //   when writing header and records, we iterate over ch2slot rather than 0..ns-1 however, 
  //   i.e. to allow scenario where we want a reduced set, or a different ordering
  //   this is passed in via the `channels` option of WRITE
  //    where proc_write() will first have checked that all channels actually existed
  //    and will have mapped to the slot numbers.  Therefore, at this point we can
  //    always assume that this will be valid
  //
  
  std::vector<int> ch2slot;

  // a pre-specified channel list?
  if ( p_ch2slot != NULL ) 
    ch2slot = *p_ch2slot;
  else
    {
      // if channels not explicitly specified
      // just take all, in the order they are already in
      for (int s=0; s<header.ns; s++)
	ch2slot.push_back(s);

    }
  
  const int ns2 = ch2slot.size();

  if ( ns2 == 0 )
    {
      logger << "  *** no channels to write to a new EDF... bailing\n";
      return false;
    }
  else
    logger << "  writing " << ns2 << " channels\n";
  
  
  //
  // Write to file
  //

  if ( f == filename )
    Helper::halt( "cannot overwrite an existing file: " + filename );

  filename = f;

  if ( ! as_edfz ) 
    {

      FILE * outfile = NULL;
      
      if ( ( outfile = fopen( filename.c_str() , "wb" ) ) == NULL )      
	{
	  logger << " ** could not open " << filename << " for writing **\n";
	  return false;
	}

      // temporarily change, just for benefit of written header 
      if ( make_EDFC ) set_continuous();

      // write header
      header.write( outfile , ch2slot );

      // change back if needed, as subsequent commands after will be happier
      if ( make_EDFC ) set_discontinuous();

      int r = timeline.first_record();
      while ( r != -1 ) 
	{
	  // we may need to load this record, before we can write it
	  if ( ! loaded( r ) )
	    {
	      edf_record_t record( this ); 
	      record.read( r );
	      records.insert( std::map<int,edf_record_t>::value_type( r , record ) );	      
	    }
	  
	  records.find(r)->second.write( outfile , ch2slot );
	  
	  r = timeline.next_record(r);
	}
      
      fclose(outfile);
    }

  //
  // .edfz and .edfz.idx
  //

  else 
    {

      edfz_t edfz;

      if ( ! edfz.open_for_writing( filename ) )
	{
	  logger << " ** could not open " << filename << " for writing **\n";
	  return false;
	}

      
      if ( make_EDFC ) set_continuous();

      // write header (as EDFZ)
      header.write( &edfz , ch2slot );

      if ( make_EDFC ) set_discontinuous();

      
      int r = timeline.first_record();
      while ( r != -1 ) 
	{
	  
	  // we may need to load this record, before we can write it
	  if ( ! loaded( r ) )
	    {
	      edf_record_t record( this ); 
	      record.read( r );
	      records.insert( std::map<int,edf_record_t>::value_type( r , record ) );	      
	    }
	  
	  
	  // set index :
	  // record -> offset into EDFZ and time-point	  
	  //        -> string representation of EDF Annots

	  // offset into file
	  int64_t offset = edfz.tell();	  

	  // time-point offset
	  uint64_t tp = timeline.timepoint( r );
	
	  // any annots
	  const std::string edf_annot_str = edf_annots[ r ] == "" ? "." : edf_annots[ r ] ;
	  
	  // write to the index
	  edfz.add_index( r , offset , timeline.timepoint( r ) , edf_annot_str  );
	  
	  // now write to the .edfz
	  records.find(r)->second.write( &edfz , ch2slot );
	  
	  // next record
	  r = timeline.next_record(r);
	}
      

      //
      // Write .idx
      //
      
      logger << "  writing EDFZ index to " << filename << ".idx\n";
      
      // update record_size (e.g. if channels dropped)
      // at this point, we will have read all information in from
      // the existing 

      int new_record_size = 0;

      // old
      // for (int s=0;s<header.ns;s++)
      // 	new_record_size += 2 * header.n_samples[s] ; // 2 bytes each                                                       
      
      // now allowing for dropped channels
      for (int s2=0; s2<ns2; s2++)
	{
	  const int s = ch2slot[s2];
	  new_record_size += 2 * header.n_samples[s] ; // 2 bytes each                                                       
	}
      
      edfz.write_index( new_record_size );

      
      //
      // All done
      //

      edfz.close();


    }
  
  logger << "  saved new EDF" 
	 << ( header.edfplus ? ( make_EDFC ? "+C" : "+D" ) : "" )        
	 << ", " << filename << "\n";
    
  return true;
}




void edf_t::drop_signal( const int s )
{

  if ( s < 0 || s >= header.ns ) return;  
  --header.ns;

  // need to track whether this signal was in the list of signals to be read from the original file
  //  -- it needn't be, i.e. if a new channel has been created
  //  -- but if it is, we need to know this when we subsequently read in
  //     new data from disk

  // i.e. it is not whether it was on disk per se, it is whether it would have been included
  //      via an initial sig= specification, i.e. inp_signals for edf_t::attach()
  //      so then we can remove it from edf.inp_signals_n[]
  

  // get original signal slot number (-1 if not present)
  int os = header.original_signal( header.label[ s ] ) ;

  // alter header
  header.label.erase( header.label.begin() + s );
  header.annotation_channel.erase( header.annotation_channel.begin() + s );
  header.transducer_type.erase( header.transducer_type.begin() + s );
  header.phys_dimension.erase( header.phys_dimension.begin() + s );
  header.physical_min.erase( header.physical_min.begin() + s );
  header.physical_max.erase( header.physical_max.begin() + s );
  header.digital_min.erase( header.digital_min.begin() + s );
  header.digital_max.erase( header.digital_max.begin() + s );
  header.orig_physical_min.erase( header.orig_physical_min.begin() + s );
  header.orig_physical_max.erase( header.orig_physical_max.begin() + s );
  header.orig_digital_min.erase( header.orig_digital_min.begin() + s );
  header.orig_digital_max.erase( header.orig_digital_max.begin() + s );
  header.prefiltering.erase( header.prefiltering.begin() + s );
  header.n_samples.erase( header.n_samples.begin() + s );
  header.signal_reserved.erase( header.signal_reserved.begin() + s );
  header.bitvalue.erase( header.bitvalue.begin() + s );
  header.offset.erase( header.offset.begin() + s );
  
  // remove from 'primary input' list (i.e. which is used
  // when reading a new record;  these signal numbers
  // are in the original (EDF-based) counting scheme
  
  if ( os != -1 ) // i.e. present in original signal list
    {
      inp_signals_n.erase( inp_signals_n.find(os) );
    }
  
  // need to remake label2header
  header.label2header.clear();
  for (int l=0;l<header.label.size();l++)     
    if ( header.is_data_channel(l) ) 
      header.label2header[ Helper::toupper( header.label[l] ) ] = l;      
  
  // records
  int r = timeline.first_record();
  while ( r != -1 )
    {
      if ( records.find(r) != records.end() ) 
	records.find(r)->second.drop(s);
      r = timeline.next_record(r);
    }
  
  
}

void edf_record_t::drop( const int s )
{
  data[ s ].clear();
  data.erase( data.begin() + s );
//   pdata[ s ].clear();
//   pdata.erase( pdata.begin() + s );
}

void edf_t::add_signal( const std::string & label ,
			const int Fs ,
			const std::vector<double> & data ,
			double pmin , double pmax ,
			int16_t dmin , int16_t dmax )
{

  const int ndata = data.size();

  // normally, n_samples is Fs * record length.

  // *however*, as we are currently otherwise enforcing that sample rate must be an integer 
  //   we've also added a backdoor for sedf_t creation here, to allow for
  //   has sample rate < 1 Hz and very long records (e.g. 30 seconds): namely,
  //   if Fs is negative, assume this directly encode the n_samples (negative of)
  //   rather than the sample rate per say
  
  const int n_samples = Fs < 0 ? -Fs : Fs * header.record_duration ;

  if ( ndata == 0 ) 
    {
      logger << " **empty EDF, not going to add channel " << label << " **\n";
      return;
    }

  //  std::cout << "nd = " << ndata << " " << header.nr << " " << n_samples << "\n";

  // sanity check -- ie. require that the data is an appropriate length
  if ( ndata != header.nr * n_samples ) 
    {
      logger << " observed n = " << ndata << " but expected = " << header.nr << " * " <<  n_samples << " = " << header.nr * n_samples << "\n";
      Helper::halt( "internal error: problem with length of input data" );  
    }

  //
  // if not othewise specified, get physical signal min/max to determine scaling
  //
  
  if ( pmin == pmax ) 
    {      
      pmin = data[0];
      pmax = data[0];
  
      for (int i=1;i<ndata;i++) 
	{
	  if      ( data[i] < pmin ) pmin = data[i];
	  else if ( data[i] > pmax ) pmax = data[i];
	}
    }

  //
  // if no variation set arbitrary pmin / pmax 
  //
  
  if ( fabs( pmin - pmax ) <= 1e-6 )
    {
      pmin -= 1.0;
      pmax += 1.0;
    }
  
  //
  // determine bitvalue and offset
  //
  
  // if not otherwise specified, set dmin/dmax 
  if ( dmax == dmin ) // i.e. 0 == 0 if not set 
    {
      dmax = 32767;
      dmin = -32768;
    }
  
  double bv = ( pmax - pmin ) / (double)( dmax - dmin );
  double os = ( pmax / bv ) - dmax;

  // store (after converting to digital form)
  
  int c = 0;
  int r = timeline.first_record();
  
  while ( r != -1 ) 
    {
      
      ensure_loaded( r );
      
      std::vector<int16_t> t(n_samples);
      
      for (int i=0;i<n_samples;i++) 
	t[i] = edf_record_t::phys2dig( data[c++] , bv , os );

      records.find(r)->second.add_data(t);

      r = timeline.next_record(r);

    }

  // add to header
  ++header.ns;
    
  header.bitvalue.push_back( bv );
  header.offset.push_back( os );
  
  header.label.push_back( label );
  
  if ( ! Helper::imatch( label , "EDF Annotation" , 14 ) )
    header.label2header[ Helper::toupper( label ) ] = header.label.size()-1;     
  
  header.annotation_channel.push_back( ( header.edfplus ? 
					 Helper::imatch( label , "EDF Annotation" , 14 ) :
					 false ) ) ;

  header.transducer_type.push_back( "n/a" );
  header.phys_dimension.push_back( "n/a" );
  header.physical_min.push_back( pmin );
  header.physical_max.push_back( pmax );
  header.digital_min.push_back( dmin );
  header.digital_max.push_back( dmax );
  header.orig_physical_min.push_back( pmin );
  header.orig_physical_max.push_back( pmax );
  header.orig_digital_min.push_back( dmin );
  header.orig_digital_max.push_back( dmax );
  header.prefiltering.push_back( "n/a" );
  header.n_samples.push_back( n_samples );  
  header.signal_reserved.push_back( "" );  

  // add to TYPES, by recallig this
  cmd_t::define_channel_type_variables( *this );
  
}


std::vector<double> edf_record_t::get_pdata( const int s )
{
  const double & bv     = edf->header.bitvalue[s];
  const double & offset = edf->header.offset[s];
  const int n = data[s].size();
  std::vector<double> r( n );
  for ( int i = 0 ; i < n ; i++ ) r[i] = dig2phys( data[s][i] , bv , offset );
  return r;
}

void edf_record_t::add_data( const std::vector<int16_t> & d )
{
  // store
  data.push_back( d );
}

void edf_record_t::add_annot( const std::string & str )
{
  // create a new data slot
  std::vector<int16_t> dummy; data.push_back(dummy);
  //std::vector<double> pdummy;   pdata.push_back(pdummy);
  // add this to the end
  add_annot( str , data.size()-1 );
}

void edf_record_t::add_annot( const std::string & str , const int signal )
{
  
  if ( signal < 0 || signal >= data.size() ) 
    Helper::halt( "internal error in add_annot()" );

  // convert text to int16_t encoding
  data[ signal ].resize( str.size() );
  for (int s=0;s<str.size();s++) 
    data[signal][s] = (char)str[s];

}

// now redundant
// void edf_record_t::calc_data( double bitvalue , double offset  )
// {
  
//   // convert to physical scale
//   // pdata[s][j] = header->bitvalue[s] * ( header->offset[s] + d );
  
//   const std::vector<double> & pd = pdata[ pdata.size() - 1 ];
//   const int n = pd.size();
//   std::vector<int> d( n );
//   for (int i=0;i<n;i++) d[i] = pd[i]/bitvalue - offset;      
  
//   // create data given min/max etc
//   data.push_back( d );
// }


void edf_t::reset_record_size( const double new_record_duration )
{

  if ( ! header.continuous )
    Helper::halt( "can only change record size for EDF, not EDF+, currently" );

  // this changes the in-memory representation;
  // naturally, new data cannot easily be loaded from disk, so 
  // this command should always write a new EDF and then quit
  
  // original record size (seconds) , and derived value in tp
  // double record_duration;
  // uint64_t record_duration_tp;

  // required : new_record_duration

  // nothing to do?
  if ( header.record_duration == new_record_duration ) return;
  
  std::vector<int> new_nsamples;

  int new_record_size = 0;

  // check that all signals can fit evenly into the new record size
  for (int s=0;s<header.ns;s++)
    {

      if ( header.is_annotation_channel(s) )
	Helper::halt( "cannot change record size for EDF annotations: drop this signal first" );

      int    nsamples = header.n_samples[s];
      double fs = (double)nsamples / header.record_duration;
      
      // does new record size contain an integer number of sample points?

      double implied = new_record_duration * fs;
      
      int   new_nsamples1 = implied;

      if ( fabs( (double)new_nsamples1 - implied ) > 0 ) 
	Helper::halt( "signal " + header.label[s] + " has sample rate " + Helper::int2str( nsamples ) + " per record, "
		      + "\n which cannot be represented in a record of " + Helper::dbl2str( new_record_duration ) );
      
      new_nsamples.push_back( new_nsamples1 );

      // track for record size of the new EDF 
      new_record_size += 2 * new_nsamples1 ; 
      
    }
  
  // buffer for new records
  edf_record_t new_record( this );

  std::map<int,edf_record_t> new_records;
  
  // manually change size of the new buffer record
  for (int s = 0 ; s < header.ns ; s++)
    new_record.data[s].resize( new_nsamples[s] , 0 );
    

  // get implied number of new records (truncate if this goes over)
  int new_nr = floor( header.nr * header.record_duration ) / (double) new_record_duration ;

  for (int r=0;r<new_nr;r++) 
    new_records.insert( std::map<int,edf_record_t>::value_type( r , new_record ) );
  
  // process one signal at a time
  std::vector<int> new_rec_cnt( header.ns , 0 );
  std::vector<int> new_smp_cnt( header.ns , 0 );
  
  int r = timeline.first_record();
  while ( r != -1 ) 
    {
  
      ensure_loaded( r );

      edf_record_t & record = records.find(r)->second;

      for (int s = 0 ; s < header.ns ; s++ )
	{

	  const int n = header.n_samples[s];

	  for (int i = 0 ; i < n ; i++ )
	    {
	      
	      if ( new_smp_cnt[s] == new_nsamples[s] )
		{
		  ++new_rec_cnt[s];
		  new_smp_cnt[s] = 0;
		}
	      
	      if ( new_rec_cnt[s] < new_nr )
		{
		  std::map<int,edf_record_t>::iterator rr = new_records.find( new_rec_cnt[s] );
		  if ( rr == new_records.end() ) Helper::halt( "internal error" );
		  edf_record_t & new_record = rr->second;

//  		  std::cout << "setting " << new_rec_cnt[s] << "\t" << new_smp_cnt[s] << " = " << r << " " << i << "\n";
//  		  std::cout << " sz = " << new_record.data[ s ].size() << " " << record.data[ s ].size() << "\n";
		  new_record.data[ s ][ new_smp_cnt[ s ] ] = record.data[ s ][ i ];
		  
		  ++new_smp_cnt[ s ];
		}

	    } // next sample point
	  
	} // next signal

      r = timeline.next_record(r);

    } // next record
  

  //
  // copy over
  //
  
  records = new_records;
  new_records.clear();

  //
  // and update EDF header
  //

  header.nr = new_nr;
  header.n_samples = new_nsamples;
  header.record_duration = new_record_duration;
  header.record_duration_tp = header.record_duration * globals::tp_1sec;

  // also, update edf.record_size : we won't be reading anything else from the original 
  // EDF, but if we are writing an EDFZ, then edf_t needs the /new/ record size for the
  // index

  record_size = new_record_size ; 

  // make a new timeline 
  timeline.re_init_timeline();

  // all done

}


void edf_t::reference_and_scale( const int s , const int r , const double rescale )
{
  
  //
  // reference and/or rescale 
  //
  
  if ( s < 0 || s >= header.ns ) Helper::halt( "incorrectly specified signal" );
  
  bool hasref = r != -1;
  if ( r < -1 || r >= header.ns || r == s ) Helper::halt( "incorrectly specified reference" );
  
  //
  // check comparable sampling rate
  //
  
  if ( hasref && header.n_samples[ s ] != header.n_samples[ r ] ) 
    Helper::halt( "reference must have similar sampling rate" );
  
  const int ns = header.n_samples[ s ];
  
  
  //
  // for every record (masked or otherwise), 
  // subtract out the reference and rescale (e.g. mV -> uV)
  // 
  
  std::vector<double> d; 
  
  int rec = timeline.first_record();
  while ( rec != -1 )
    {

      ensure_loaded( rec );
      
      edf_record_t & record = records.find(rec)->second;
      
      if ( hasref )
	{
	  std::vector<double> pdata_sig = record.get_pdata(s);
	  std::vector<double> pdata_ref = record.get_pdata(r);
	  
	  for (int i=0;i<ns;i++) d.push_back( ( pdata_sig[i] - pdata_ref[i] ) * rescale );

	}
      else	
	{
	  std::vector<double> pdata_sig = record.get_pdata(s);
	  for (int i=0;i<ns;i++) 
	    {
	      d.push_back( pdata_sig[i] * rescale );
	      //std::cout << "rescale " << pdata_sig[i] << " " << pdata_sig[i] * rescale << "\n";
	    }
	}
      
      rec = timeline.next_record(rec);
      
    }
  
  // update signal
  update_signal( s , &d );
  
}


void edf_t::pairwise_reference( const signal_list_t & signals ,
				const signal_list_t & refs ,
				const bool make_new ,
				const std::vector<std::string> & new_channels ,
				const int new_sr ,
				bool dereference ,
				const bool verbose )
{
  const int ns = signals.size();
  const int nr = refs.size();
  const int nw = new_channels.size();
  
  if ( ns != nr ) Helper::halt( "sig and ref must be same size with 'pairwise' " );
  if ( make_new && nw != ns ) Helper::halt( "sig and new must be same size with 'pairwise' " );
  for (int s=0; s<ns; s++)
    {
      signal_list_t s1 = header.signal_list( signals.label( s ) ) ;
      signal_list_t s2 = header.signal_list( refs.label( s ) ) ;
      
      reference( s1, s2, make_new, new_channels[s], new_sr, dereference, verbose );      
    }
}


void edf_t::reference( const signal_list_t & signals0 ,
		       const signal_list_t & refs ,
		       const bool make_new ,
		       const std::string & new_channel , 
		       const int new_sr ,
		       bool dereference ,
		       const bool verbose )
{

  // copy as we may modify this
  signal_list_t signals = signals0;
  
  const int ns = signals.size();
  const int nr = refs.size();

  // need at least one channel specified

  if ( ns == 0 ) 
    Helper::halt( "must specify sig={ ... }" );
    
  //
  // Create a new channel?
  //

  if ( make_new && ns > 1 )
    Helper::halt( "can only re-reference a single channel if 'new' is specified" );

  std::string ch_label = "" ;
    
  if ( make_new )
    {
      // retain original label for output below
      ch_label = signals.label(0);
      
      // make copy
      copy_signal( header.label[ signals(0) ] , new_channel );
      
      // switch to re-reference this copy now
      signals = header.signal_list( new_channel );

      // do we need to resample? 
    
      const int sig_sr = (int)header.sampling_freq( signals(0) );
     
      // resample sig, if needed (only one)
      // this slot is the 'new' one, so original signal untouched
      if ( new_sr != 0 && sig_sr != new_sr )
	dsptools::resample_channel( *this , signals(0) , new_sr );

      
      // if the reference needs resampling, we need to copy a new
      // channel and do the re-sampling (i.e. to leave the original
      // untouched).   Do this downstream on the 'final' reference
      // (as this may involve multiple channels that have different 
      // SRs.
      
    }

  
  //
  // if nr size is 0, means leave as is
  // if we've requested a new channel, we need to make this still
  //
  
  if ( nr == 0 )
    {
      if ( ! make_new ) Helper::halt( "no valid ref channels specified" );

      // else...
      return;
    }

  //
  // Console logging 
  //

  if ( verbose && nr > 0 )
    {
      logger << ( dereference ? "  dereferencing" : "  referencing" );
      for (int s=0;s<ns;s++) logger << " " << ( make_new ? ch_label : header.label[ signals(s) ] ) ;
      logger << " with respect to";
      if ( nr > 1 ) logger << " the average of";      
      for (int r=0;r<nr;r++) logger << " " << header.label[ refs(r) ];      
      if ( make_new ) logger << " --> " << header.label[ signals(0) ];
      logger << "\n";
    }



  //
  // check SR for all channels  
  //
  
  int np_sig = header.n_samples[ signals(0) ];

  if ( (!make_new ) || ( make_new && new_sr == 0 ) ) 
    {      
      
      for (int s=0;s<ns;s++) 
	if ( header.n_samples[ signals(s) ] != np_sig ) 
	  Helper::halt( "all signals/references must have similar sampling rates" );
      
      for (int r=0;r<nr;r++) 		
	if ( header.n_samples[ refs(r) ] != np_sig ) 
	  Helper::halt( "all signals/references must have similar sampling rates" );
      
    }
  else
    {
      // here we are fixing SR, and we've already done this for the 
      // signal;  we'll do it later for REF, but need to check they
      // all (if >1 ref, i.e. average ref) match

      int np_ref = header.n_samples[ refs(0) ];
      
      for (int r=0;r<nr;r++) 		
	if ( header.n_samples[ refs(r) ] != np_ref ) 
	  Helper::halt( "all references must have similar sampling rates" );

    }


  //
  // Build reference once
  //
  
  std::vector<double> reference;

  // number of samples points per record for reference
  const int np_ref = header.n_samples[ refs(0) ];

  int rec = timeline.first_record();
  while ( rec != -1 )
    {
      ensure_loaded( rec );

      edf_record_t & record = records.find(rec)->second;
      
      std::vector<std::vector<double> > refdata;

      // get data
      for (int r=0;r<nr;r++) 		
	refdata.push_back( record.get_pdata( refs(r) ) );
      
      // average
      for (int i=0;i<np_ref;i++) 
	{
	  double avg = 0;
	  for (int r=0;r<nr;r++) avg += refdata[r][i];
	  if ( nr != 1 ) avg /= (double)nr;
	  reference.push_back( avg );
	}      

      // next record
      rec = timeline.next_record(rec);       
    }


  //
  // Need to resample reference?
  //

  if ( make_new && new_sr != 0 )
    {
      const int ref_sr = (int)header.sampling_freq( refs(0) );
      if ( ref_sr != new_sr )
	{
	  const int refsize = reference.size();
	  
	  reference = dsptools::resample( &reference , ref_sr , new_sr );
	  
	  // ensure exact length... pad if needed
	  if ( reference.size() != refsize )
	    reference.resize( refsize );
	}
    }
  
  //
  // transform signals one at a time, now we have reference in 'reference'
  //
  
  for (int s=0;s<signals.size();s++) 
    {
      
      // do not reference to self
      if ( nr == 1 && signals(s) == refs(0) ) 
	{
	  if ( verbose )
	    logger << " skipping " << refs.label(0) << " to not re-reference to self\n"; 
	  continue;
	}
      
      // transformed signal      
      std::vector<double> d;
      int cc = 0;

      //
      // iterate over records
      // 
      
      int rec = timeline.first_record();
      while ( rec != -1 )
	{
	  
	  ensure_loaded( rec );
	  
	  // now we can access
	  edf_record_t & record = records.find(rec)->second;
	  
	  std::vector<double> d0 = record.get_pdata( signals(s) );
	  
	  if ( dereference ) 
	    for (int i=0;i<np_sig;i++) d.push_back( d0[i] + reference[cc++] );
	  else	    
	    for (int i=0;i<np_sig;i++) d.push_back( d0[i] - reference[cc++] );
	  
	  // next record
	  rec = timeline.next_record(rec); 
	  
	}
      
      // update signal
      update_signal( signals(s) , &d );
      
      // next signal to re-reference
    }

}


bool edf_t::load_annotations( const std::string & f0 )
{

  //
  // parse annotation filename
  //

  const std::string f = Helper::expand( f0 );


  // allow wildcards
    
  if ( ! Helper::fileExists( f ) ) 
    Helper::halt( "annotation file " + f + " does not exist for EDF " + filename );

  //
  // store filename (if needed to be output in a WRITE to the sample-list)
  //
  
  annot_files.push_back( f );

  //
  // Type of input?
  //

  bool xml_mode = Helper::file_extension( f , "xml" );
  
  bool feature_list_mode = Helper::file_extension( f , "ftr" );
  
  //
  // XML files (NSRR, Profusion or Luna formats)
  //
  
  if ( xml_mode ) 
    {
      annot_t::loadxml( f , this );
      return true;
    }
  

  //
  // Feature lists
  //
  
  if ( feature_list_mode && globals::read_ftr )
    {
      
      std::vector<std::string> tok = Helper::parse( f , "/" );

      std::string file_name = tok[ tok.size()-1];	
    
      // filename should be id_<ID>_feature_<FEATURE>.ftr
      int pos = file_name.find( "_feature_" );
      
      if ( pos == std::string::npos || file_name.substr(0,3) != "id_" )  
	Helper::halt( "bad format for feature list file name: id_<ID>_feature_<FEATURE>.ftr" );
      
      std::string id_name = file_name.substr(3,pos-3);

      if ( id_name != id )
	{
	  Helper::warn( ".ftr file id_{ID} does not match EDF ID : [" + id_name + "] vs [" + id + "]" );
	  return false;
	}
      
      std::string feature_name = file_name.substr( pos+9 , file_name.size() - 4 - pos - 9 );
  
      // are we checking whether to add this file or no? 
      
      if ( globals::specified_annots.size() > 0 && 
	   globals::specified_annots.find( feature_name ) == globals::specified_annots.end() ) return false;
      
      // create and load annotation
      
      annot_t * a = timeline.annotations.add( feature_name );
      
      a->name = feature_name;
      a->description = "feature-list";
      a->file = file_name;

      // load features, and track how many
      aoccur[ feature_name ] = a->load_features( f  );
      
      return true;
    }


  //
  // Otherwise, process as an .annot or .eannot file
  //
  
  return annot_t::load( f , *this );
  
}


int  edf_header_t::signal( const std::string & s , bool silent )
{  
  signal_list_t slist = signal_list(s);
  if ( slist.size() != 1 ) 
    {
      if ( ! silent ) 
	logger << " ** could not find signal [" << s << "] of " << label2header.size() << " signals **\n";
      return -1;
    }  
  return slist(0);
}


bool  edf_header_t::has_signal( const std::string & s )
{
  std::vector<std::string> tok = Helper::parse( s , "|" );    
  for (int t=0;t<tok.size();t++)
    {
      // primary name (that might be an alias)?
      if ( label2header.find( Helper::toupper( tok[t] ) ) != label2header.end() )
	return true;
      
      // using aliased (i.e. original) name?
      if ( cmd_t::label_aliases.find( Helper::toupper( tok[t] ) ) != cmd_t::label_aliases.end() )
	return true;
    }
  return false;
}


int  edf_header_t::original_signal_no_aliasing( const std::string & s  )
{  
  std::map<std::string,int>::const_iterator ff = label_all.find( Helper::toupper( s ) );
  if ( ff != label_all.end() ) return ff->second;
  return -1;
}


int  edf_header_t::original_signal( const std::string & s  )
{  

  // look up, with aliases, in original
  // label_all[ ]

  const std::string uc_s = Helper::toupper( s );
  
  std::map<std::string,int>::const_iterator ff = label_all.find( uc_s );
  
  if ( ff != label_all.end() ) return ff->second;
  
  // otherwise, consider if we have aliases
  if ( cmd_t::label_aliases.find( uc_s ) != cmd_t::label_aliases.end() )
    {
      const std::string & s2 = cmd_t::label_aliases[ uc_s ];

      ff = label_all.find( Helper::toupper( s2 ) );

      if ( ff != label_all.end() ) return ff->second;
    }

  // otherwise, look to a primary term
  
  if ( cmd_t::primary_upper2orig.find( uc_s ) != cmd_t::primary_upper2orig.end() )
    {
      // swap PRIMARY -> Primary, and then pull all aliases
      // this returns ALIASES , so we can use w/ label_all[] directly
      std::vector<std::string> & a = cmd_t::primary_alias.find( cmd_t::primary_upper2orig[ uc_s ] )->second;
      for (int i=0;i<a.size();i++)
	{
	  ff = label_all.find( a[i] );
	  if ( ff != label_all.end() ) return ff->second;
	}
    }
  
  return -1;
  
}


signal_list_t edf_header_t::signal_list( const std::string & s , bool no_annotation_channels , bool show_warnings )
{

  signal_list_t r;
  
  // wildcard means all signals '*'
  
  if ( s == "*" )
    {
      for (int s=0;s<label.size();s++)
	{
	  
	  // ? only consider data tracks
	  
	  if ( no_annotation_channels  && 
	       is_annotation_channel( s ) ) continue;
	  
	  std::string lb = label[s];
	  
	  std::string uppercase_lb = Helper::toupper( lb );
	  
	  // swap in alias? [ aliases are always stored as UPPERCASE ]
	  if ( cmd_t::label_aliases.find( uppercase_lb ) != cmd_t::label_aliases.end() ) 
	    {
	      // track
	      aliasing[ cmd_t::label_aliases[ uppercase_lb ] ] = lb;
	      
	      // swap in the primary
	      lb = cmd_t::label_aliases[ uppercase_lb ];
	      label2header[ Helper::toupper( lb ) ] = s;
	      label[s] = lb;
	      
	    }
	  
	  r.add( s, lb );
	}
    }

  //
  // comma-delimited; but within a single signal specification,
  // we allow a pipe-delimited list, where we pick the first that matches
  //  
  //  A,B|C,D|E|F,G  - mean A and ( B or C ) and ( D or E or F ) and G 
  //  

  std::vector<std::string> tok = Helper::quoted_parse( s , "," );    
  
  for (int t=0;t<tok.size();t++)    
    {

      std::vector<std::string> tok2_ = Helper::quoted_parse( tok[t] , "|" );    
      
      // first swap in any aliases, and place those at the front of the list
      // then continue as before

      // swap in alias first? -- this may double alias, but fine.
      
      // e.g. 
      // alias   sigX|sigY|sigZ
      // signal  sigY|sigZ|sig0
      // will make all --> sigX (which is correct)

      std::string alias = "";

      for (int t2=0;t2<tok2_.size();t2++)    
	{
	  const std::string uc_lb = Helper::toupper( tok2_[t2] );
	  
	  if ( cmd_t::primary_upper2orig.find( uc_lb ) != cmd_t::primary_upper2orig.end() )
	    {
	      if ( cmd_t::primary_alias.find( cmd_t::primary_upper2orig[ uc_lb ] ) != cmd_t::primary_alias.end() )
		{
		  //std::cout << "tok2_ " << tok2_[t2] << "and alias ["<<alias<< "]\n";
		  if ( alias == "" ) alias = cmd_t::primary_upper2orig[ uc_lb ] ;  //  OLD = tok2_[t2];
		  else if ( ! Helper::iequals( alias , uc_lb ) )
		    Helper::halt( "more than one alias implied" );
		  //std::cout << "tok2_ " << tok2_[t2] << "and alias ["<<alias<< "]\n";
		}
	    }
	  else if ( cmd_t::label_aliases.find( uc_lb ) != cmd_t::label_aliases.end() ) 
	    {
	      if ( alias == "" ) 
		alias = cmd_t::label_aliases[ uc_lb ];	  
	      else if ( ! Helper::iequals( alias , cmd_t::label_aliases[ uc_lb ] ) )
		Helper::halt( "more than one alias implied" );
	    }
	}
      
      //
      // update list if needed
      //
      // std::cout << " alias = [" << alias << "]\n";
      // std::cout << " primary alias size = " <<  cmd_t::primary_alias.size() << "\n";
      
      std::vector<std::string> tok2;
      if ( alias != "" ) 
	{
	  tok2.push_back( alias );
	  const std::vector<std::string> & avec = cmd_t::primary_alias.find( alias )->second;
	  std::vector<std::string>::const_iterator aa = avec.begin();
	  while ( aa != avec.end() )
	    {	      
	      tok2.push_back( *aa );
	      ++aa;
	    }	  
	  for (int t2=0;t2<tok2_.size();t2++) 
	    {
	      if ( tok2_[t2] != alias ) 
		tok2.push_back( tok2_[t2] );   
	    }
	}
      else
	tok2 = tok2_;
      
      std::set<int> added;

      
      //
      // proceed as before
      //

      //      std::cout << " tok2 = " << tok2.size() << "\n";

      for (int t2=0;t2<tok2.size();t2++)    
	{
	  
	  // std::cout << "t2 = " << t2 << "\t" << tok2[t2] << "\n";
	  // std::cout << "label2header.find size " << label2header.size() << "\n";
	  
	  // add first match found 
	  if ( label2header.find( Helper::toupper( tok2[t2] ) )  != label2header.end() ) 
	    {
	      
	      const int l = label2header[ Helper::toupper( tok2[t2] ) ];
	      
	      //	      std::cout << "found match " << l << "\n";
	      
	      if ( t2 > 0 ) // relabel if wasn't first choice?
		{
		  label2header[ Helper::toupper( tok2[0] ) ] = l;
		  //label[l] = tok2[0];
		}	  
	      
	      //std::cout << "adding N " << label2header[ Helper::toupper( tok2[0] ) ]  << "\n";
	      
	      const int l0 = label2header[ Helper::toupper( tok2[0] ) ] ;
	      
	      if ( added.find( l0 ) == added.end() )
		{
		  r.add( l0 , label[l] );
		  added.insert( l0 );
		}
	      
	      break;
	    }
	}
     
    }

  return r;
}

void edf_header_t::rename_channel( const std::string & old_label , const std::string & new_label )
{
  // expects exact match (i.e. this only called from XML <Signals> / <CanonicalLabel> information  
  // also by SIGNALS pick/rename
  for (int s=0;s<label.size();s++) if ( label[s] == old_label ) label[s] = new_label;
  label_all[ Helper::toupper( new_label ) ] = label_all[ Helper::toupper( old_label ) ];
  label2header[ Helper::toupper( new_label ) ] = label2header[ Helper::toupper( old_label ) ];
  
}


double edf_header_t::sampling_freq( const int s ) const
{  
  if ( s < 0 || s >= n_samples.size() ) return -1;
  return n_samples[ s ] / record_duration;
}

std::vector<double> edf_header_t::sampling_freq( const signal_list_t & signals ) const
{
  const int n = signals.size();
  std::vector<double> fs( n );
  for (int s=0;s<n;s++)
    fs[s] = n_samples[ signals.signals[s] ] / record_duration;
  
  return fs;
}


void edf_header_t::check_channels()
{
  // when loading EDF, we would have made unique  (a, a.1, a.2, etc) any identical channel
  // names;  here we also need to check that aliases aren't making non-unique labels
  // e.g. 
  //   A|B,C
  // but EDF has both "B" and "C" 

  // therefore, for each cmd_t::primary_alias[term].vector[] we need to make sure that 
  // we do not see more than one instance

  bool okay = true;

  std::map<std::string,std::vector<std::string> >::const_iterator ii = cmd_t::primary_alias.begin();
  while ( ii != cmd_t::primary_alias.end() )
    {
      std::set<std::string> obs;
      std::vector<std::string>::const_iterator jj = ii->second.begin();
      while ( jj != ii->second.end() )
	{
	  if ( original_signal_no_aliasing( *jj ) != -1 ) obs.insert( *jj );
	  ++jj;
	}
      if ( obs.size() > 1 ) 
	{
	  okay = false;
	  logger << " different channels map to the same alias term: "
		 << ii->first << " <- " << Helper::stringize( obs , " | " ) << "\n";
	}
      ++ii;
    }
  
  if ( ! okay ) 
    Helper::halt( "problem: different channels present in the EDF are mapped to the same alias" );

  
}


bool edf_t::restructure()
{
  
  //
  // Map back onto original epochs
  //
  
  timeline.set_epoch_mapping();

  // Check that we have anything to do
  
  if ( ! timeline.is_epoch_mask_set() ) 
    {
      logger << "  no epoch mask set, no restructuring needed\n";
      
      writer.value( "NR1" , header.nr );
      writer.value( "NR2" , header.nr );
      
      writer.value( "DUR1" , header.nr * header.record_duration );
      writer.value( "DUR2" , header.nr * header.record_duration );

      return false;
    }
  

  bool any_records_dropped = false;
  int cnt = 0;
  int r = timeline.first_record();
  while ( r != -1 ) 
    {
      if ( timeline.masked_record( r ) )
	{
	  any_records_dropped = true;
	  break;
	}
      ++cnt;
      r = timeline.next_record(r);
    }

  // nothing to do...
  if ( ! any_records_dropped ) 
    {

      writer.value( "NR1" , cnt );
      writer.value( "NR2" , cnt );

      writer.value( "DUR1" , cnt * header.record_duration );
      writer.value( "DUR2" , cnt * header.record_duration );

      return false;

    }


  //
  // We now will have a discontinuous EDF+
  //

  if ( ! header.edfplus ) 
    {
      logger << "  restructuring as an EDF+:";
      set_edfplus();
    }
      
  set_discontinuous();


  //
  // First, check that all necessary records have actually been loaded
  // This will not reload records, and will only load 'retained' records
  //
  
  // Alternative??  if not found below, just read
  //  read_records( 0 , header.nr_all - 1 );

  //
  // Ensure is loaded if we need it
  //

  std::set<int> include;

  for (int r = 0 ; r < header.nr_all; r++)
    {
      
      bool found     = records.find(r) != records.end();
      bool retained  = timeline.retained(r);
      bool unmasked  = !timeline.masked_record(r);
      
      if ( retained )
	if ( unmasked ) 
	  {
	    if ( ! found ) read_records( r, r );	      
	    include.insert( r );
	  }    
    }

  
  //
  // Remove records based on epoch-mask
  //

  std::map<int,edf_record_t> copy = records;
  records.clear();
  
  //
  // Copy back, but now use iterator instead
  //

  std::set<int>::const_iterator ii = include.begin();
  while ( ii != include.end() )
    {
      records.insert( std::map<int,edf_record_t>::value_type( *ii , copy.find( *ii )->second ) );
      ++ii;
    }


  if ( 0 ) 
    {
      for (int r = 0 ; r < header.nr_all; r++)
	{
	  
	  bool found     = copy.find(r) != copy.end();
	  bool retained  = timeline.retained(r);
	  bool unmasked  = !timeline.masked_record(r);
	  
	  if ( retained )
	    {
	      if ( unmasked ) 
		{	      
		  if ( ! found ) Helper::halt( "internal error in restructure()");
		  records.insert( std::map<int,edf_record_t>::value_type( r , copy.find(r)->second ) );
		  include.insert( r );
		}
	    }
	}
    }

  
      
  // set warning flags, if not enough data left
  
  if ( records.size() == 0 ) globals::empty = true;
    

  logger << " keeping " 
	 << records.size() << " records of " 
	 << copy.size() << ", resetting mask\n";
  
  writer.value( "NR1" , (int)copy.size() );
  writer.value( "NR2" , (int)records.size() );
  
  writer.value( "DUR1" , copy.size() * header.record_duration );
  writer.value( "DUR2" , records.size() * header.record_duration );

  int n_data_channels = 0 , n_annot_channels = 0;

  for (int s=0;s<header.ns;s++)
    {
      if ( header.is_data_channel(s) )
        ++n_data_channels;
      else
        ++n_annot_channels;
    }

  // signal info -- total number of channels (data / annot )
  writer.value( "NS" , n_data_channels );
  writer.value( "NA" , n_annot_channels );
  

  // update EDF header
  // nb. header.nr_all stays the same, reflecting the 
  // original file which has not changed

  header.nr = records.size();

  // adjust timeline (now will be a discontinuous track)
  
  timeline.restructure( include );

  return true;

}



void edf_t::update_physical_minmax( const int s )
{
  
  interval_t interval = timeline.wholetrace();  
  slice_t slice( *this , s , interval );
  const std::vector<double> * d = slice.pdata();
  const int n = d->size();

  double pmin = (*d)[0];
  double pmax = (*d)[0];
  
  for (int i=1;i<n;i++)
    {
      if      ( (*d)[i] < pmin ) pmin = (*d)[i];
      else if ( (*d)[i] > pmax ) pmax = (*d)[i];
    }
  
  header.physical_min[s] = pmin;
  header.physical_max[s] = pmax;  
  
  // update bitvalue/offset also

  header.bitvalue[s] = ( pmax - pmin ) / (double)( header.digital_max[s] - header.digital_min[s] );
  header.offset[s] = ( pmax / header.bitvalue[s] ) - header.digital_max[s] ;

}


void edf_t::shift( int s , int shift_sp , bool wrap )
{

  if ( shift_sp == 0 ) return;

  // i.e. parameterize as +ve means to push the series forward
  shift_sp = -shift_sp;
  
  // get data : note, this ignores EDF discontinuities
  
  slice_t slice( *this , s , timeline.wholetrace() );

  const std::vector<double> * d = slice.pdata();
  
  const int np = d->size();

  if ( np <= shift_sp ) return;
  
  std::vector<double> d2( np , 0 );
  
  for (int i=0;i<np;i++)
    {

      int j = i - shift_sp;
      
      if ( j < 0 )
	{
	  if ( wrap ) 
	    {
	      j = np - shift_sp + i;
	      d2[j] = (*d)[i];
	    }
	}
      else if ( j >= np ) 
	{ 
	  if ( wrap ) 
	    {
	      j = j - np;
	      d2[j] = (*d)[i];
	    }
	}
      else
	{
	  d2[j] = (*d)[i];
	}
    }
  
  update_signal( s , &d2 );

}


void edf_t::set_order( param_t & param )
{
  
  

}


void edf_t::copy_signal( const std::string & from_label , const std::string & to_label )
{
  
  const int s1 = header.signal( from_label );
  
  if ( s1 == -1 ) 
    Helper::halt( "could not find signal " + from_label );
  
  if ( header.has_signal( to_label ) ) 
    Helper::halt( to_label + " already exists in the EDF" );
  

  //
  // get data
  //
  
  interval_t interval = timeline.wholetrace();  
  slice_t slice( *this , s1 , interval );
  const std::vector<double> * d = slice.pdata();
  
  
  //
  // add signal (w/ same pmin/pmax and dmin/dmax)
  //
  
  add_signal( to_label , header.sampling_freq(s1) , *d ,
	      header.physical_min[s1] , header.physical_max[s1] ,
	      header.digital_min[s1] , header.digital_max[s1] 
	      );

  
  //
  // and copy the header values that would not have been properly set by add_signal()
  //

  const int s2 = header.signal( to_label );
  
  if ( s2 == -1 ) 
    Helper::halt( "problem with COPY: could not find new signal " + to_label );
  
  header.transducer_type[s2] = header.transducer_type[s1];
  header.phys_dimension[s2] = header.phys_dimension[s1];
  header.prefiltering[s2] = header.prefiltering[s1];
  
}


void edf_t::update_records( int a , int b , int s , const std::vector<double> * d )
{

  if ( header.is_annotation_channel(s) ) 
    Helper::halt( "edf_t:: internal error, cannot update an annotation channel" );

  // keep digital min/max scale as is.

  // for signal s, place back data in 'd' into EDF record structure
  // and update the physical min/max

  const int points_per_record = header.n_samples[s];
  const int n_records = b - a + 1 ;

  //std::cerr << "a,b = " << a << " " << b << " " << header.nr << "\n";

  if ( a < 0 || b < 0 || n_records <= 0 || a >= header.nr_all || b >= header.nr_all )    
    Helper::halt( "bad record specification in edf_t::update_records()" );
  const int n = d->size();
  
  if ( n != n_records * points_per_record )
    Helper::halt( "internal error in update_records()" );
  
  // use existing digital/physical min/max encoding
  // but will need to make sure we stay within digital min/max
  
  const int16_t dmin = header.digital_min[s];
  const int16_t dmax = header.digital_max[s];
  
  double pmin = header.physical_min[s];
  double pmax = header.physical_max[s];
  
  double bv = header.bitvalue[s];
  double os = header.offset[s];
  
  int cnt = 0;
  
  // assume records have already been read in... if they have, this function
  // automatically returns so okay to call just in case

  read_records( a , b );

  for ( int r = a ; r <= b ; r++ ) 
    {
      
      // find records      
      std::vector<int16_t>    & data  = records.find(r)->second.data[ s ];
      
      // check that we did not change sample rate      
      if ( data.size() != points_per_record ) 
	Helper::halt( "changed sample rate, cannot update record" );
      
      for (int p=0;p<points_per_record;p++)
	{
	  double x = (*d)[cnt];
	  if ( x < pmin ) x = pmin;
	  else if ( x > pmax ) x = pmax;

// 	  std::cout << "edit\t" << edf_record_t::dig2phys( data[p] , bv , os ) 
// 		    << "\t"
// 		    << (*d)[cnt] << "\n";
	  data[p] = edf_record_t::phys2dig( (*d)[cnt] , bv , os );
	  ++cnt;	  
	}
    }
}

  
void edf_t::update_signal_retain_range( int s , const std::vector<double> * d )
{
  if ( s < 0 || s >= header.ns ) Helper::halt( "bad 's' value in update_signal_retain_range()" );

  int16_t dmin = header.digital_min[s];
  int16_t dmax = header.digital_max[s];
  double pmin = header.physical_min[s];
  double pmax = header.physical_max[s];

  update_signal( s , d , &dmin, &dmax, &pmin, &pmax );
}

void edf_t::update_signal( int s , const std::vector<double> * d , int16_t * dmin_ , int16_t * dmax_ , double * pmin_ , double * pmax_ )
{

  const bool debug = false;
  
  // if non-null, use these dmax/dmin pmax/pmin values to update signal 
  
  if ( header.is_annotation_channel(s) ) 
    Helper::halt( "edf_t:: internal error, cannot update an annotation channel" );
  
  // for signal s, place back data in 'd' into EDF record structure
  // and update the physical min/max

  const int points_per_record = header.n_samples[s];
  const int n = d->size();

  if ( n != header.nr * points_per_record )
    Helper::halt( "internal error in update_signal()" );

  if ( debug ) std::cout << " n = " << n << "\n";
  
  bool set_minmax = dmin_ != NULL ; 
  
  // use full digital min/max scale if not otherwise specified
  int16_t dmin = -32768;
  int16_t dmax = 32767;  

  double pmin = (*d)[0];
  double pmax = (*d)[0];

  if ( set_minmax )
    {
      pmin = *pmin_;
      pmax = *pmax_;
      dmin = *dmin_;
      dmax = *dmax_;

      if ( dmin == dmax )
	{
	  dmin = -32768;
	  dmax = 32767;
	}
      else if ( dmin > dmax )
	{
	  dmin = *dmax_;
	  dmax = *dmin_;
	}

      if ( pmin == pmax )
        {
          pmin--;
          pmax++;
        }
      else if ( pmin > pmax )
        {
          pmin = *pmax_;
          pmax = *pmin_;
        }

    }
  else
    {
      // empirically find physical min/max for this signal
      for (int i=0;i<n;i++)
	{
	  if      ( (*d)[i] < pmin ) pmin = (*d)[i];
	  else if ( (*d)[i] > pmax ) pmax = (*d)[i];
	}

      // exapand range as needed
      if ( fabs( pmin - pmax ) < 1e-6 )
	{
	  pmin--;
	  pmax++;
	}
      
    }

  if ( debug )
    {
      std::cout << " pmin, pmax = " << pmin << " " << pmax << "\n";
      std::cout << " dmin, dmax = " << dmin << " " << dmax << "\n";
    }

  
  //
  // update header min/max (but leave orig_physical_min/max unchanged)
  //
  
  header.digital_min[s] = dmin;
  header.digital_max[s] = dmax;

  header.physical_min[s] = pmin;
  header.physical_max[s] = pmax;

  
  double bv = ( pmax - pmin ) / (double)( dmax - dmin );
  double os = ( pmax / bv ) - dmax;

  // TODO: we should flag in the header that this signal has been transformed
  // (?always) we expect that all data will have been read, so no more will be read
  // from disk;  if it were, these could be the incorrect BV/OS values, of course. 
  // We should set a flag that header has changed, and bail out if subsequent reads are made

  header.bitvalue[s] = bv;
  header.offset[s] = os;
  
  int cnt = 0;

  if ( debug ) std::cout << " records[] size = " << records.size() << "\n";
  
  int r = timeline.first_record();
  while ( r != -1 ) 
    {

      if ( debug )
	{
	  std::cout << " r = " << r << "\n";
	  if ( records.find(r) == records.end() ) std::cout << " could not find record\n";
	  std::cout << records.find(r)->second.data.size() << " is data[] size\n";
	  std::cout << " s = " << s << "\n";
	}
		  
      // find records
      
      //      std::vector<double> & pdata = records.find(r)->second.pdata[ s ];
      std::vector<int16_t>    & data  = records.find(r)->second.data[ s ];
      
      // check that we did not change sample rate
      
      if ( data.size() != points_per_record ) 
	{
	  //pdata.resize( points_per_record , 0 );
	  data.resize( points_per_record , 0 );
	}

      
      for (int p=0;p<points_per_record;p++)
	{
	  
	  // check that does not exceed range
	  double x = (*d)[cnt];
	  if ( x < pmin ) x = pmin;
	  if ( x > pmax ) x = pmax;
	  
	  // also need to convert to bit-value
	  // pdata[s][j] = header->bitvalue[s] * ( header->offset[s] + d );	  
	  // reverse digital --> physical scaling
	  
	  //data[p] = (*d)[cnt] / bv - os;
	  data[p] = edf_record_t::phys2dig( x , bv , os );
	
	  ++cnt;	  
	}

      r = timeline.next_record(r);
    }
}




edf_record_t::edf_record_t( edf_t * e ) 
{    

  edf = e;

  // only store digital value, convert on-the-fly
  data.resize( edf->header.ns );
  //pdata.resize( edf->header.ns );
  
  for (int s = 0 ; s < edf->header.ns ; s++)
    {
      if ( edf->header.is_annotation_channel(s) )
	{
	  data[s].resize( 2 * edf->header.n_samples[s] , 0 );
	  //pdata[s].resize( edf->header.n_samples[s] * fac , 0 );
	}
      else
	{
	  data[s].resize( edf->header.n_samples[s] , 0 );
	  //pdata[s].resize( edf->header.n_samples[s] , 0 );
	}
    }
  
}


void edf_t::reset_start_time()
{

  // get time of first record
  int r = timeline.first_record();
  if ( r == -1 ) return;

  interval_t interval = timeline.record2interval(r);

  // keep as is?
  if ( interval.start == 0 ) return;
  
  // interval for this record
  logger << "  setting EDF start time from " << header.starttime ;

  clocktime_t et( header.starttime );
  if ( et.valid )
    {
      double time_sec = ( interval.start * globals::tp_duration ) ; 
      et.advance_seconds( time_sec );
    }

  header.starttime = et.as_string();
  logger << " to " << header.starttime  << "\n"; 
}


void edf_t::set_continuous()
{
  if ( ! header.edfplus ) return;
  header.continuous = true;
  header.reserved[0] = 'E';
  header.reserved[1] = 'D';
  header.reserved[2] = 'F';
  header.reserved[3] = '+';
  header.reserved[4] = 'C';

}

void edf_t::set_discontinuous()
{
  if ( ! header.edfplus ) return;
  header.continuous = false;
  header.reserved[0] = 'E';
  header.reserved[1] = 'D';
  header.reserved[2] = 'F';
  header.reserved[3] = '+';
  header.reserved[4] = 'D';
}


void edf_t::set_edfplus()
{
  if ( header.edfplus ) return;
  header.edfplus = true;
  header.continuous = true;    
  set_continuous(); // this sets reversed field EDF+C
  add_time_track();
}

void edf_t::set_edf()
{
  if ( ! header.edfplus ) return;
  header.edfplus = false;
  header.continuous = true;  
  header.reserved[0] = ' '; 
  header.reserved[1] = ' ';
  header.reserved[2] = ' ';
  header.reserved[3] = ' ';
  header.reserved[4] = ' ';

  set_continuous(); 
  drop_time_track();
  drop_annots();
}

void edf_t::drop_annots()
{
  // drop all 'EDF Annot' signals from in-memory EDF
  // (i.e. part of making an EDF from an EDF+
  
  for (int s=0;s<header.ns;s++)
    if ( header.is_annotation_channel(s) )
      drop_signal( s );

  has_edf_annots = false;
  
}

void edf_t::drop_time_track()
{
  // means that the EDF will become 'continuous'
  set_continuous();

  // no TT in any case?
  if ( header.time_track() == -1 ) return;
  drop_signal( header.time_track() );
  
}


int edf_t::add_time_track( const std::vector<uint64_t> * tps )
{
  
  // if tps == NULL, this implies a continuous record
  //   - this will be the typical case -- i.e. if it is
  //     an EDF+D/discontinuous, then (by definition) we will
  //     have read in a time-track

  // however, one exception to this is when merging standard EDFs
  // to make a new EDF (--merge) and if there are gaps between files.
  // here we need to set the EDF+D time-track explicitly- which is
  // done by calling this function by having tps != NULL but a vector
  // of time-points for each record

  const bool contin = tps == NULL;
  
  if ( contin && ! header.continuous ) 
    return header.time_track();
  
  if ( ! header.edfplus )
    set_edfplus();

  
  // time-track already set?
  if ( contin && header.time_track() != -1 )
    return header.time_track();  
  
  // check EDF+D time-track size, if specified
  if ( ! contin )
    if ( tps->size() != header.nr )
      Helper::halt( "internal error: expecting " + Helper::int2str( header.nr )
		    + " records but given time-track for " + Helper::int2str( (int)tps->size() ) );


  // add a new time-track?
  
  if ( header.time_track() == -1 ) 
    {
      
      // update header
      ++header.ns;
      
      // set t_track channel
      header.t_track  = header.ns - 1;
      header.t_track_edf_offset = record_size; // i.e. at end of record
      
      const int16_t dmax = 32767;
      const int16_t dmin = -32768;
      
      // need to set a record size -- this should be enough?
      // default (defs/defs.cpp) is currently 15 (i.e. 30 chars)
      const int n_samples = globals::edf_timetrack_size;
      
      // how many existing 'EDF Annotations' tracks?
      int annot_tracks = 0 ; 
      
      std::map<std::string,int>::const_iterator jj = header.label_all.begin();
      while ( jj != header.label_all.end() )
	{
	  if ( Helper::imatch( jj->first  , "EDF Annotation" , 14 ) ) 
	    annot_tracks++;                     
	  ++jj;
	}
      
      header.label.push_back( "EDF Annotations" + ( annot_tracks > 0 ? Helper::int2str( annot_tracks ) : "" ) );
      header.annotation_channel.push_back( true );

      // note: annot, so not added to header/record signal map label2header

      header.transducer_type.push_back( "" );
      header.phys_dimension.push_back( "" );
      
      header.physical_min.push_back( 0 ); // ignored
      header.physical_max.push_back( 1 ); // ignored
      header.digital_min.push_back( dmin );
      header.digital_max.push_back( dmax );
      
      header.orig_physical_min.push_back( 0 ); // ignored
      header.orig_physical_max.push_back( 1 ); // ignored
      header.orig_digital_min.push_back( dmin );
      header.orig_digital_max.push_back( dmax );
      
      header.prefiltering.push_back( "" );
      header.n_samples.push_back( n_samples );
      header.signal_reserved.push_back( "" );  
      header.bitvalue.push_back( 1 ); // ignored
      header.offset.push_back( 0 );   // ignored
    }
  
  // create each 'TAL' timestamp, and add to record
  double dur_sec = header.record_duration;
  double onset = 0; // start at T=0 [ for EDF+C ], else uses tps[] below

  // uint64_t onset_tp = 0;
  // uint64_t dur_tp = header.record_duration_tp;

  // for each record
  int r = timeline.first_record();

  // counter (EDF+D only, to index tps[]) 
  int rc = 0;

  while ( r != -1 ) 
    {

      // either EDF+C or EDF+D times
      double tsec = contin ? onset : ( (*tps)[rc] / (double)globals::tp_1sec ) ;
      
      std::string ts = "+" + Helper::dbl2str( tsec ) + "\x14\x14\x00";

      // need to make sure that the record (i.e. other signals) 
      // are first loaded into memory...
      
      bool record_in_memory = loaded(r);
      
      if ( ! record_in_memory )
	{
	  
	  // this will be created with ns+1 slots (i.e. 
	  // already with space for the new timetrack, 
	  // so we can add directly)
	  
	  edf_record_t record( this ); 
	  
	  record.read( r );

	  records.insert( std::map<int,edf_record_t>::value_type ( r , record ) );

	}

	  
      //
      // Add the time-stamp as the new track (i.e. if we write as EDF+)
      //

      if ( contin )
	{
      
	  if ( ! record_in_memory ) // record structure already 'updated' from above
	    records.find(r)->second.add_annot( ts , header.t_track );
	  else // push_back on end of record
	    records.find(r)->second.add_annot( ts );
	}
      else
	{
	  // different logic for the EDF+D / --merge case

	  // here, we are adding a EDF+D time-track (from --merge)
	  // there will already be a time-track
	  records.find(r)->second.add_annot( ts , header.t_track );
	  
	}
      
      //
      // And mark the actual record directy (i.e. if this is used in memory)
      // for EDF+C  (if EDF+D, this does not matter)
      
      onset += dur_sec;
      //onset_tp += dur_tp;
      
      // next record [ used for EDF+D ]
      ++rc;
      
      r = timeline.next_record(r);
    }

//   std::cout << "DET1 " << records.begin()->second.pdata.size() << "\t"
// 	    << records.begin()->second.data.size() << "\n";


  return header.time_track();

}



uint64_t edf_t::timepoint_from_EDF( int r )
{

  //
  // for EDFZ, this will be stored in the .idx
  //

  if ( file == NULL )
    return edfz->get_tindex( r );

  //
  // Read this is called when constructing a time-series for 
  // an existing EDF+D, only
  //
  
  if ( ! header.edfplus ) Helper::halt( "should not call timepoint_from_EDF for basic EDF");
  if (   header.continuous ) Helper::halt( "should not call timepoint_from_EDF for EDF+C");
  if (   header.time_track() == -1 ) Helper::halt( "internal error: no EDF+D time-track" );
  
  // allocate buffer space
  int ttsize = 2 * globals::edf_timetrack_size;  
  byte_t * p = new byte_t[ ttsize ];
  byte_t * p0 = p;
  
  // determine offset into EDF
  uint64_t offset = header_size + (uint64_t)(record_size) * r;      
  offset += header.time_track_offset(); 
  
  // time-track is record : edf->header.time_track 
  // find the appropriate record
  fseek( file , offset , SEEK_SET );
      
  // and read only time-track (all of it)
  size_t rdsz = fread( p , 1, ttsize , file );
  
  std::string tt( ttsize , '\x00' );

  int e = 0;
  for (int j=0; j < ttsize; j++)
    {      
      tt[j] = *p;
      if ( tt[j] == '\x14' || tt[j] == '\x15' ) break;
      ++p;
      ++e;
    }
  
  double tt_sec = 0;

  if ( ! Helper::str2dbl( tt.substr(0,e) , &tt_sec ) ) 
    Helper::halt( "problem converting time-track in EDF+" );
  
  delete [] p0;
  
  uint64_t tp = globals::tp_1sec * tt_sec;

  return tp; 
  
}
  
void edf_t::flip( const int s )
{
  if ( header.is_annotation_channel(s) ) return;
  logger << "  flipping polarity of " << header.label[s] << "\n";

  // get all data
  interval_t interval = timeline.wholetrace();
  slice_t slice( *this , s , interval );
  const std::vector<double> * d = slice.pdata();
  std::vector<double> rescaled( d->size() );
  
  for (int i=0;i<d->size();i++)  rescaled[i] = - (*d)[i];

  // update signal (and min/max in header)
  update_signal( s , &rescaled );
  
}

void edf_t::reverse( const int s )
{
  if ( s < 0 || s >= header.ns ) return;
  
  if ( header.is_annotation_channel(s) ) return;
  logger << "  reversing  " << header.label[s] << "\n";

  // get all data
  interval_t interval = timeline.wholetrace();
  slice_t slice( *this , s , interval );
  const std::vector<double> * d = slice.pdata();
  const int np = d->size();
  std::vector<double> reversed( np );
  for (int i=0;i<np;i++)
    {
      reversed[i] = (*d)[np-i-1];
      //      std::cout << " reversed[i]  = " << i << "\t" << reversed[i]  << "\n";
    } 
  update_signal_retain_range( s , &reversed );  
}



void edf_t::rescale( const int s , const std::string & sc , const bool quietly )
{
  
  if ( header.is_annotation_channel(s) ) return;

  bool is_mV = header.phys_dimension[s] == "mV";
  bool is_uV = header.phys_dimension[s] == "uV";
  bool is_V  = header.phys_dimension[s] == "V";

  bool rescale_from_mV_to_uV = is_mV && sc == "uV"; // *1000
  bool rescale_from_uV_to_mV = is_uV && sc == "mV"; // /1000
  
  bool rescale_from_V_to_uV = is_V && sc == "uV"; // * 1e6
  bool rescale_from_V_to_mV = is_V && sc == "mV"; // * 1e3

  if ( ! ( rescale_from_mV_to_uV || rescale_from_uV_to_mV 
	   || rescale_from_V_to_uV || rescale_from_V_to_mV ) ) 
    {
      //logger << " no rescaling needed\n";
      return;
    }

  // get all data
  interval_t interval = timeline.wholetrace();
  slice_t slice( *this , s , interval );
  const std::vector<double> * d = slice.pdata();
  std::vector<double> rescaled( d->size() );

  // get rescaling factor
  double fac = 1;
  if      ( rescale_from_uV_to_mV ) fac = 1.0/1000.0;
  else if ( rescale_from_mV_to_uV ) fac = 1000;
  else if ( rescale_from_V_to_mV )  fac = 1000;
  else if ( rescale_from_V_to_uV )  fac = 1000000;
  
  // rescale
  for (int i=0;i<d->size();i++) 
    {
      rescaled[i] = (*d)[i] * fac;
    }

  // update signal (and min/max in header)
  update_signal( s , &rescaled );

  // update headers
  if ( rescale_from_mV_to_uV || rescale_from_V_to_uV ) 
    {
      if ( ! quietly )
	logger << "  rescaled " << header.label[s] << " to uV\n";
      header.phys_dimension[s] = "uV";     
    }
  
  if ( rescale_from_uV_to_mV || rescale_from_V_to_mV ) 
    {
      if ( ! quietly )
	logger << "  rescaled " << header.label[s] << " to mV\n";
      header.phys_dimension[s] = "mV";
    }
}



void edf_t::minmax( signal_list_t & signals )
{

  int16_t dmax = 0;
  int16_t dmin = 0;
  double pmin = 0 , pmax = 0;
  
  bool any_set = false;

  const int ns = signals.size();

  for (int s=0; s < ns; s++)
    {

      if ( ! header.is_data_channel( signals(s) ) ) continue;

      if ( ! any_set )
	{
	  pmin = header.physical_min[ signals(s) ] ;
	  pmax = header.physical_max[ signals(s) ] ;	  
	  dmin = header.digital_min[ signals(s) ] ;
	  dmax = header.digital_max[ signals(s) ] ;
	  any_set = true;
	}
      else
	{
	  if ( header.physical_min[ signals(s) ] < pmin ) pmin = header.physical_min[ signals(s) ];
	  if ( header.physical_max[ signals(s) ] > pmax ) pmax = header.physical_max[ signals(s) ];
	  if ( header.digital_min[ signals(s) ] < dmin ) dmin = header.digital_min[ signals(s) ];
	  if ( header.digital_max[ signals(s) ] > dmax ) dmax = header.digital_max[ signals(s) ];
	}
    }

  //
  // now rescale each channel to these identical EDF scales
  //

  interval_t interval = timeline.wholetrace();
 
  for (int s=0; s < ns; s++)
    {
      if ( ! header.is_data_channel( signals(s) ) ) continue;
      
      slice_t slice( *this , signals(s) , interval );

      const std::vector<double> * d = slice.pdata();
      
      update_signal( signals(s) , d , &dmin, &dmax , &pmin , &pmax );
    }
  
  
}


bool edf_t::basic_stats( param_t & param )
{
  
  // Run through each record
  // Get min/max
  // Calculate RMS for each signal
  // Get mean/median/SD and skewness/kurtosis
  // optinoally, display a histogram of observed values (and figure out range)
  
  std::string signal_label = param.requires( "sig" );  

  signal_list_t signals = header.signal_list( signal_label );

  std::vector<double> Fs = header.sampling_freq( signals );
  
  bool by_epoch = param.has( "epoch" );

  bool hist = param.has( "encoding" );
  
  const int ns = signals.size();
  
  bool calc_median = true;
  
  int required_sr = param.has( "sr-under" ) ? param.requires_int( "sr-under" ) : 0 ; 

  const bool minimal = param.has( "min" ) || param.has( "minimal" );

  const bool run_pcts = param.has( "pct" ) ? param.yesno("pct") : true ;
  
  for (int s=0; s<ns; s++)
    {
      
      //
      // skip annotation channels
      //
      
      if ( header.is_annotation_channel( signals(s) ) ) continue;

      //
      // SR requirements?
      //

      if ( required_sr != 0 && header.sampling_freq( signals(s) ) > required_sr ) continue;

      if ( header.sampling_freq( signals(s) ) == 0 ) continue;
	   
      //
      // Output signal
      //
      
      writer.level( header.label[ signals(s) ] , globals::signal_strat );

      //
      // Mean, variance, skewmess/kurtosis, RMS, min, max based on per-epoch stats
      //
            
      std::vector<double> e_mean, e_median , e_sd, e_rms, e_skew, e_kurt;
      
      double t_min = 0 , t_max = 0;
      
      logger << " processing " << header.label[ signals(s) ] << " ...\n";
      
      
      //
      // EPOCH-level statistics first
      //
      
      if ( by_epoch ) 
	{

	  timeline.first_epoch();  
	  	  
	  //
	  // Iterate over epcohs
	  //
	  
	  while ( 1 ) 
	    {
	      
	      int epoch = timeline.next_epoch();      
	      
	      if ( epoch == -1 ) break;
	      
	      interval_t interval = timeline.epoch( epoch );
	      
	      //
	      // Get data 
	      //
	      
	      slice_t slice( *this , signals(s) , interval );
	      
	      const std::vector<double> * d = slice.pdata();
	      
	      const int n = d->size();
	      
	      if ( n == 0 ) { continue; } 
	      
	      
	      //
	      // Filter data
	      //
	      
	      double mean   = MiscMath::mean( *d );
	      double median = calc_median ? MiscMath::median( *d ) : 0;
	      
	      double sd     = minimal ? 0 : MiscMath::sdev( *d , mean );
	      double rms    = minimal ? 0 : MiscMath::rms( *d );
	      double skew   = minimal ? 0 : MiscMath::skewness( *d , mean , sd );
	      double kurt   = minimal ? 0 : MiscMath::kurtosis( *d , mean );
	      
	      double min = (*d)[0];
	      double max = (*d)[0];
	      if ( ! minimal ) 
		for (int i = 0 ; i < n ; i++ )
		  {
		    if ( (*d)[i] < min ) min = (*d)[i];
		    if ( (*d)[i] > max ) max = (*d)[i];
		  }
	      
	      std::map<int,double> pct;
	      if ( run_pcts )
		{
		  //pct[ -1 ]  = MiscMath::percentile( *d , 0.001 );
		  pct[ 1 ]  = MiscMath::percentile( *d , 0.01 );
		  pct[ 2 ]  = MiscMath::percentile( *d , 0.02 );
		  pct[ 5 ]  = MiscMath::percentile( *d , 0.05 );
		  pct[ 95 ] = MiscMath::percentile( *d , 0.95 );
		  pct[ 98 ] = MiscMath::percentile( *d , 0.98 );
		  pct[ 99 ] = MiscMath::percentile( *d , 0.99 );
		  //pct[ 999 ]  = MiscMath::percentile( *d , 0.999 );
		  for (int pp=0;pp<9;pp++)
		    pct[ 10 + pp * 10 ] = MiscMath::percentile( *d , 0.1 + pp*0.1 );
		}
	      
	      //
	      // Output
	      //
	      
	      writer.epoch( timeline.display_epoch( epoch ) );
	      
	      writer.value( "MEAN" , mean );

	      if ( calc_median ) 
		writer.value( "MEDIAN" , median );	      

	      if ( ! minimal )
		{
		  writer.value( "MAX"  , max  );
		  writer.value( "MIN"  , min  );	      
		  
		  if ( Helper::realnum( skew ) )
		    writer.value( "SKEW" , skew );
		  
		  if ( Helper::realnum( kurt ) )
		    writer.value( "KURT" , kurt );
		  
		  writer.value( "RMS"  , rms  );

		  if ( run_pcts )
		    {
		      std::map<int,double>::const_iterator pp = pct.begin() ;
		      while ( pp != pct.end() )
			{
			  if ( pp->first == -1 )
			    writer.value( ( "P001" ) + Helper::int2str( pp->first ) , pp->second ) ;
			  else		    
			    writer.value( ( pp->first < 10 ? "P0" : "P" ) + Helper::int2str( pp->first ) , pp->second ) ;
			  ++pp;
			}
		    }
		}
	      
	      
	      //
	      // Record
	      //
	      
	      e_mean.push_back( mean );

	      if ( calc_median ) 
		e_median.push_back( median );

	      if ( ! minimal )
		{
		  if ( t_min == 0 && t_max == 0 ) 
		    { 
		      t_min = min; 
		      t_max = max; 
		    } 
		  
		  if ( min < t_min ) t_min = min;
		  if ( max > t_max ) t_max = max;
		  
		  e_sd.push_back( sd );
		  e_rms.push_back( rms );	  
		  e_skew.push_back( skew );
		  e_kurt.push_back( kurt );
		}
	    }
	  
	  writer.unepoch();
	  
	  
	}
      
      
      //
      // Whole-signal level output
      //
      
      interval_t interval = timeline.wholetrace();
      
      slice_t slice( *this , signals(s) , interval );
	  
      const std::vector<double> * d = slice.pdata();
      
      const int n = d->size();

      if ( n == 0 ) { continue; } 
 
      double mean = MiscMath::mean( *d );
      //double median = calc_median ? MiscMath::median( *d ) : 0 ;

      writer.value( "MEAN" , mean );
      
      if ( ! minimal )
	{
	  double rms  = MiscMath::rms( *d );
	  double sd = MiscMath::sdev( *d );
	  double skew = MiscMath::skewness( *d , mean , sd );
	  double kurt = MiscMath::kurtosis( *d , mean );
	  double min = (*d)[0];
	  double max = (*d)[0];

      
	  for (int i = 0 ; i < n ; i++ )
	    {
	      if ( (*d)[i] < min ) min = (*d)[i];
	      if ( (*d)[i] > max ) max = (*d)[i];
	    }


	  std::map<int,double> pct;
	  if ( run_pcts )
	    {	      
	      pct[ 1 ]  = MiscMath::percentile( *d , 0.01 );
	      pct[ 2 ]  = MiscMath::percentile( *d , 0.02 );
	      pct[ 5 ]  = MiscMath::percentile( *d , 0.05 );
	      pct[ 95 ] = MiscMath::percentile( *d , 0.95 );
	      pct[ 98 ] = MiscMath::percentile( *d , 0.98 );
	      pct[ 99 ] = MiscMath::percentile( *d , 0.99 );
	      for (int pp=0;pp<9;pp++)
		pct[ 10 + pp * 10 ] = MiscMath::percentile( *d , 0.1 + pp*0.1 );
	    }
	  
	  //
	  // Output
	  //
	  	  
	  writer.value( "MAX"  , max  );
	  writer.value( "MIN"  , min  );      	  
	  writer.value( "SKEW" , skew );
	  writer.value( "KURT" , kurt );
	  
	  //if ( calc_median ) writer.value( "MEDIAN" , median );
	  
	  writer.value( "RMS"  , rms  );
	  writer.value( "SD"  , sd  );

	  if ( run_pcts )
	    {
	      std::map<int,double>::const_iterator pp = pct.begin() ;
	      while ( pp != pct.end() ) 
		{
		  writer.value( ( pp->first < 10 ? "P0" : "P" ) + Helper::int2str( pp->first ) , pp->second ) ;
		  ++pp;
		}
	    }
	}
      
      //
      // Also, same strata:  summaries of epoch-level statistics
      //
      
      if ( by_epoch && e_mean.size() > 0 )
	{
	  const int ne = e_mean.size(); 
	  double med_mean  = median_destroy( &e_mean[0] , ne );
	  double med_median  = calc_median ? median_destroy( &e_median[0] , ne ) : 0 ;  
	  writer.value( "MEDIAN.MEAN" , med_mean );
	  if ( calc_median )
	    writer.value( "MEDIAN.MEDIAN" , med_median );

	  writer.value( "NE" , timeline.num_total_epochs() );	  
	  writer.value( "NE1" , ne );
	  
	  if ( ! minimal )
	    {
	  
	      double med_rms  = median_destroy( &e_rms[0] , ne );
	      double med_skew = median_destroy( &e_skew[0] , ne );
	      double med_kurt = median_destroy( &e_kurt[0] , ne );
	      
	      writer.value( "MEDIAN.RMS"  , med_rms );
	      writer.value( "MEDIAN.SKEW" , med_skew );
	      writer.value( "MEDIAN.KURT" , med_kurt );
	    }
	}

      
      //
      // Optional, encoding 
      //
      
      // verbose output: every unique value / count 
      if ( hist )
	{
	  
	  std::map<double,int> counts;
	  for (int i = 0 ; i < n ; i++ )
	    counts[ (*d)[i] ]++;
	  
	  writer.value( "OBS_ENCODING" , (int)counts.size() );
	  
	  // largest possible EDF digital span
	  
	  int span_obs = header.digital_max[ signals(s) ] - header.digital_min[ signals(s) ] + 1;
	  
	  int zero_cells = span_obs - counts.size();
	  
	  writer.value( "MAX_ENCODING" , span_obs );
	  writer.value( "PCT_ENCODING" , counts.size() / (double)span_obs );
	  
	  std::map<double,int>::const_iterator ii = counts.begin();
	  while ( ii != counts.end() )
	    {
	      writer.level( ii->first , globals::value_strat );
	      writer.value( "CNT" , ii->second );
	      ++ii;
	    }
	  writer.unlevel( "VAL" );

	}
      
      //
      // Next channel
      //
            
    }

  //
  // All done
  //

  writer.unlevel( globals::signal_strat );
  
  return true;
  
}




bool signal_list_t::match( const std::set<std::string> * inp_signals ,
			   std::string * l ,
			   const std::set<std::string> & slabels )
{

  //  std::cout << " IN signal_list_t::match() \n";
  
  // inp_signals : list of input signals (EDF or subset of) / any-CASE
  // l           : label to match : any-case, and we want to preserve this, but matching is done in case-insensitive manner
  // slabels     : 
  
  // exact match? (i.e. no "|" alternatives specified)
  // old:   if ( inp_signals->find(*l) != inp_signals->end() ) return true; 
  // new: for case-insensitive match, we need to explicitly consider each
  std::set<std::string>::const_iterator cc = inp_signals->begin();
  while ( cc != inp_signals->end() )
    {
      if ( Helper::iequals( *l , *cc ) ) return true;
      ++cc;
    }
  
  // alternatively, as an alias?  
  if ( cmd_t::label_aliases.find( Helper::toupper( *l ) ) != cmd_t::label_aliases.end() )
    {
      *l = cmd_t::label_aliases[ Helper::toupper( *l ) ];
      // now, does this matc
      // old: return inp_signals->find(*l) != inp_signals->end() ;
      // new: loop over each (as above)
      std::set<std::string>::const_iterator cc = inp_signals->begin();
      while ( cc != inp_signals->end() )
	{
	  if ( Helper::iequals( *l , *cc ) ) return true;
	  ++cc;
	}
      return false;
    }
  
  // subset match (i.e. one of x|y|z)
  // if both 'x' and 'y' exist, always pick 'x' first
  
  std::set<std::string>::const_iterator ii = inp_signals->begin();
  while ( ii != inp_signals->end() )
    {
      std::vector<std::string> tok = Helper::parse( *ii , "|" );
      for (int i=0;i<tok.size();i++) 
	{
	  
	  // if gone preferred value exists in some other slot, then this is not a match
	  // i.e. only include one selection, the preferred one
	  if ( i>0 && slabels.find( tok[0] ) != slabels.end() )  break; 
	    
	  if ( *l == tok[i] ) 
	    {
	      // swap in 'preferred' name
	      if ( i>0 ) *l = tok[0];
	      return true;
	    }
	}
      ++ii;
    }    
  return false;
}








bool edf_t::append( const std::string & filename , 
		    const std::vector<std::string> & channels , 
		    const std::vector<std::vector<std::vector<double> > > & data ) 
{

  // data[rec][channels][samples]
  
  if ( channels.size() == 0 ) return false;
  if ( data.size() == 0 ) return false;

  //
  // Read header of the original (base) EDF
  //

  edf_t base;

  base.attach( filename , "." , NULL , true ); // true implies silent mode (no console logs)
  
  //
  // Check this is not EDFZ
  //
  
  if ( base.edfz )
    Helper::halt( "cannot append to EDFZ" );

  
  //
  // key measures
  //
  //  base.header.ns      # of signals
  //  base.header.nr      # of records (will be increased)
  //  base.header.recdur  # irrelevant

  //  base.header.label[i]       channel label, should match EXACTLY data[r].first
  //  base.header.n_samples[i]   # samples per reccord for signal i, should match data[r][s].size() 
  
  const int n_new_records = data.size();

  //
  // Check all original channels exist, and have correct n_samples[]
  //  - allow for different order
  //  - ignore channels in new data but not present in the original
  
  std::map<std::string,int> ch2slot;
  for (int s=0; s<channels.size(); s++)
    ch2slot[ channels[s] ] = s;

  // consider each original channel - all must be present in the new 
  // data -- but now the order and # does not otherwise have to align
  // we will use ch2ch[] below to select the correct channel from the new 
  // data 
  
  std::vector<int> ch2ch( base.header.ns );
  for (int s=0; s<base.header.ns; s++)
    {
      if ( ch2slot.find( base.header.label[s] ) == ch2slot.end() )
	Helper::halt( "could not find " + base.header.label[s] + " in the to-be-appended data" );
      ch2ch[s] = ch2slot[ base.header.label[s] ];
    }
  
  // if ( base.header.ns != channels.size() )
  //   Helper::halt( "must have exactly same set of channels to append to an EDF" );
  
  // for (int s=0; s<base.header.ns; s++)
  //   if ( channels[s] != base.header.label[s] )
  //     Helper::halt( "must have exactly the same order of channels to append to an EDF" );
  


  //
  // Just check first record, but **assume** all records have the same length (will be checked when writing)
  //
  
  // not required any more: allowing more channels in new data (which will be ignored)
  // if ( data[0].size() != base.header.ns )
  //   Helper::halt( "data[] must exactly match EDF # of channels" );
  
  // nb ch2ch[] mapping from original to new slots
  for (int s=0; s < base.header.ns; s++)
    if ( data[0][ ch2ch[s] ].size() != base.header.n_samples[s] )
      Helper::halt( "data[] must have exactly the same # of samples per record to append" );
  
  //
  // Store key values : these necessarily match the original/base
  //

  const int orig_nr = base.header.nr;
  const int orig_ns = base.header.ns;
  const std::vector<std::string> orig_label = base.header.label;
  const std::vector<double> orig_physical_min = base.header.physical_min;
  const std::vector<double> orig_physical_max = base.header.physical_max;
  const std::vector<int> orig_digital_min = base.header.digital_min;
  const std::vector<int> orig_digital_max = base.header.digital_max;
  const std::vector<int> orig_nsamples = base.header.n_samples;

  //
  // All looks okay, so close original
  //

  base.init();

  //
  // Re-open for reading and writing 
  //

  FILE * mergefile = NULL;
    
  if ( ( mergefile = fopen( filename.c_str() , "rb+" ) ) == NULL )
    Helper::halt( "problem opening " + filename + " to edit header" );

  //
  // Update NR in headar
  //

  std::string c = Helper::int2str( n_new_records + orig_nr );
  c.resize(8,' ');
  fseek( mergefile , 236 , SEEK_SET );  
  fwrite( c.data() , 1 , 8 , mergefile );
  
  //
  // Go to end of file 
  //

  fseek( mergefile , 0, SEEK_END); 

  const int nr = data.size();

  //
  // Precompute EDF offset/bv
  //

  std::vector<double> bv( orig_ns );
  std::vector<double> os( orig_ns );

  for (int s=0; s<orig_ns; s++)
    {      
      bv[s] = ( orig_physical_max[s] - orig_physical_min[s] ) / (double)( orig_digital_max[s] - orig_digital_min[s] );
      os[s] = ( orig_physical_max[s] / bv[s] ) - orig_digital_max[s];
    }
  
  //
  // Iterate over records
  //
  std::map<std::string,int> clippings;
  
  for (int r = 0 ; r < nr ; r++ )
    {
      
      for (int s=0;s<orig_ns;s++)
	{

	  const int nsamples = orig_nsamples[s];

	  // nb. select the correct slot from the new data using ch2ch[]
	  // otherwise, below use s rather than ch2ch[s] as we are 
	  // placing into slot s
	  const std::vector<double> & d = data[r][ ch2ch[s] ];
	  
	  if ( d.size() != nsamples ) 
	    Helper::halt( "hmm... internal error in append()" );
	  
	  const double pmin = orig_physical_min[s] ;
	  const double pmax = orig_physical_max[s] ;
	  
	  for (int j=0;j<nsamples;j++)
	    {

	      double pvalue = d[j];

	      //
	      // range checking
	      //
	      
	      if ( pvalue < pmin )
		{
		  ++clippings[ orig_label[s] ];
		  pvalue = pmin;
		}
	      else if ( pvalue > pmax )
		{
		  ++clippings[ orig_label[s] ];
		  pvalue = pmax;
		}

	      //
	      // physical --> digital scaling [ double --> int16_t ] 
	      //
	      
	      int16_t dvalue = edf_record_t::phys2dig( pvalue , bv[s] , os[s] );

	      //
	      // write in little-endian
	      //
	      
	      char a , b;
              edf_record_t::dec2tc( dvalue , &a, &b );	      
              fputc( a , mergefile );
              fputc( b , mergefile );


	    } // next sample
	} // next channel
    } // next record


  //
  // warnings?
  //

  if ( clippings.size() > 0 )
    {
      logger << "  *** warning: physical values outside of EDF-header specified physical min/max ranges:\n";
      std::map<std::string,int>::const_iterator cc = clippings.begin();
      while ( cc != clippings.end() )
	{
	  logger << "   " << cc->second << " samples for " << cc->first << "\n";
	  ++cc;
	}
    }
  
  //
  // all done, close up the file 
  //

  fclose( mergefile );

  logger << "  appended " << n_new_records << " to " << filename << ", and updated header\n";

  return true;
}


void edf_t::set_headers( param_t & param )
{

  
  if ( param.has( "id" ) )
    {
      header.patient_id = param.value( "id" ) ;
      // also update edf.id
      id = header.patient_id;
      // this will not change any DB output though
      logger << "  set 'id' to " << header.patient_id << "\n";
      if ( header.patient_id.size() > 80 )
	logger << "  *** warning - 'id' will be truncated to 80 characters if saved as EDF\n";
    }
  
  if ( param.has( "recording-info" ) )
    {
      header.recording_info = param.value( "recording-info" ) ;
      logger << "  set 'recording-info' to " << header.recording_info << "\n";
      if ( header.recording_info.size() > 80 )
        logger << "  *** warning - 'recording-info' will be truncated to 80 characters if saved as EDF\n";
    }
  
  if ( param.has( "start-date" ) )
    {
      header.startdate = param.value( "start-date" );
      logger << "  set 'start-date' to " << header.startdate << "\n";
      if ( header.startdate.size() > 8 )
        logger << "  *** warning - 'start-date' will be truncated to 8 characters if saved as EDF\n";
    }

  if ( param.has( "start-time" ) )
    {
      header.starttime = param.value( "start-time" );
      logger << "  set 'start-time' to " << header.starttime << "\n";
      if ( header.starttime.size() > 8 )
        logger << "  *** warning - 'start-time' will be truncated to 8 characters if saved as EDF\n";
    }

  const bool no_annot_channels = true;
  
  signal_list_t signals = header.signal_list( param.value( "sig" ) , no_annot_channels );  
  
  const int ns = signals.size();

  for (int s=0;s<ns;s++)
    {

      const int slot = signals(s);
      
      // transducer_type
      // phys_dimension
      // prefiltering

      if ( param.has( "transducer" ) )
	{
	  header.transducer_type[ slot ] = param.value( "transducer" );
	  logger << "  set " << signals.label(s) << " 'transducer' to " << header.transducer_type[ slot ] << "\n";
	  if ( s == 0 && header.transducer_type[ slot ].size() > 80 )
	    logger << "  *** warning - 'transducer' will be truncated to 80 characters if saved as EDF\n";	  
	}

      if ( param.has( "physical-dimension" ) )
        {
          header.phys_dimension[ slot ] = param.value( "physical-dimension" );
	  logger << "  set " << signals.label(s) << " 'physical-dimension' to " << header.phys_dimension[ slot ] << "\n";
          if ( s == 0 && header.phys_dimension[ slot ].size() > 8 )
            logger << "  *** warning - 'physical-dimension' will be truncated to 8 characters if saved as EDF\n";
        }
      else if ( param.has( "unit" ) )
	{
          header.phys_dimension[ slot ] = param.value( "unit" );
	  logger << "  set " << signals.label(s) << " 'unit' to " << header.phys_dimension[ slot ] << "\n";
	  if ( s == 0 && header.phys_dimension[ slot ].size() > 8 )
            logger << "  *** warning - 'unit' will be truncated to 8 characters if saved as EDF\n";
        }

      if ( param.has( "prefiltering" ) )
        {
          header.prefiltering[ slot ] = param.value( "prefiltering" );
	  logger << "  set " << signals.label(s) << " 'prefiltering' to " << header.prefiltering[ slot ] << "\n";
          if ( s == 0 && header.prefiltering[ slot ].size() > 80 )
            logger << "  *** warning - 'prefiltering' will be truncated to 80 characters if saved as EDF\n";
        }
      
    }
  
}

// int edf_t::read_all()
// {
  
//   for (int r = 0 ; r < header.nr_all; r++)
//     {	
//       bool found     = records.find(r) != records.end();
//       bool retained  = timeline.retained(r);
//       bool unmasked  = !timeline.masked_record(r);
      
//       if ( retained )
// 	if ( unmasked ) 
// 	  if ( ! found )
// 	    {
// 	      read_records( r, r );
// 		++rex;
// 	    }
//     }
//   return rex;
// }


void edf_t::update_edf_pointers( edf_t * p )
{
  for (int r = 0 ; r < header.nr_all; r++)
    {
      bool found = records.find(r) != records.end();
      records.find(r)->second.edf = p; 
    }
}
