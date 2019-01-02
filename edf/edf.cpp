
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
#include "miscmath/miscmath.h"
#include "db/db.h"

#include <iostream>
#include <fstream>

extern writer_t writer;


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

long long edf_t::get_filesize(FILE *file)
{  
  long lCurPos, lEndPos;
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
      if ( ! globals::silent ) 
	std::cerr << "returning -1: [" << s << "] is not a valid real number\n";
      return -1;
    }
  return t;
}

std::string edf_t::get_string( byte_t ** p , int sz )
{
  std::vector<char> buf(sz+1);
  for (int i=0;i<sz;i++)
    {
      buf[i] = **p;
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
     << "Duration       : " << record_duration << "\n\n";
  
  for (int s=0;s<ns;s++)
    {
      ss << "Signal " << (s+1) << " : [" << label[s] << "]\n"
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

void edf_t::description() const
{
  
  uint64_t duration_tp = globals::tp_1sec * (uint64_t)header.nr * header.record_duration ;
  
  std::cout << "EDF file          : " << filename << "\n"
	    << "ID                : " << id << "\n"
	    << "Duration          : " << Helper::timestring( duration_tp ) << "\n"
	    << "Number of signals : " << header.ns << "\n"
	    << "Signals           :";
  
  for (int s=0;s<header.ns;s++) std::cout << " " << header.label[s];
  std::cout << "\n\n";
  
  
}


void edf_t::terse_summary() const
{
  
  writer.unlevel();

  // variable definitions
  writer.var( "NS" , "Number of signals" );
  writer.var( "NR" , "Number of records" ); 
  writer.var( "REC.DUR" , "Record duration (sec)" );
  writer.var( "TOT.DUR.SEC" , "Total recording duration (sec)" );
  writer.var( "TOT.DUR.HMS" , "Total recording duration (hh:mm:ss)" );
  writer.var( "TOT.DUR.TP"  , "Total recording duration (time-points, 10^(-12) secs)" );

  writer.var( "SR" , "Sampling race (points per second)" );
  writer.var( "PDIM" , "Physical dimension/units" );
  writer.var( "PMIN" , "Physical minimum" );
  writer.var( "PMAX" , "Physical maximum" );

  writer.var( "DMIN" , "Digital minimum" );
  writer.var( "DMAX" , "Digital maximum" );

  // write output
  writer.value( "NS" , header.ns );
  writer.value( "NR" , header.nr );
  writer.value( "REC.DUR" , header.record_duration );

  // total duration in TP units
  uint64_t duration_tp = globals::tp_1sec * (uint64_t)header.nr * header.record_duration ;
  std::string total_duration_hms = Helper::timestring( duration_tp );
  writer.value( "TOT.DUR.SEC" , header.nr * header.record_duration );
  writer.value( "TOT.DUR.HMS" , total_duration_hms );
  writer.value( "TOT.DUR.TP" , Helper::int2str( duration_tp ) );

  for (int s=0;s<header.ns;s++)
    {
      // channel name
      writer.level( header.label[s] , globals::signal_strat );
      
      // number of samples
      writer.value( "SR" , header.n_samples[s] / (double)header.record_duration );
      
      // physical dimension
      writer.value( "PDIM" , header.phys_dimension[s] );

      // physical min/max
      writer.value( "PMIN" , header.physical_min[s] );
      writer.value( "PMAX" , header.physical_max[s] );
      
      // digital min/max
      writer.value( "DMIN" , header.digital_min[s] );
      writer.value( "DMAX" , header.digital_max[s] );

    }
  
  writer.unlevel( globals::signal_strat );

}



std::set<int> edf_header_t::read( FILE * file , const std::set<std::string> * inp_signals )
{

  // Fixed buffer size for header
  // Total header = 256 + ns*256
 
  const int hdrSz = 256; 
  
  // Allocate space in the buffer for the header only
  byte_t * q = new byte_t[ hdrSz ];
  byte_t * q0 = q;

  // Read start of header into the buffer
  size_t rdsz = fread( q , 1, hdrSz , file);

  std::set<int> channels;
  
  version        = edf_t::get_string( &q , 8 );
  patient_id     = edf_t::get_string( &q , 80 );
  recording_info = edf_t::get_string( &q , 80 );
  startdate      = edf_t::get_string( &q , 8 );
  starttime      = edf_t::get_string( &q , 8 );
  nbytes_header  = edf_t::get_int( &q , 8 );  
  reserved       = edf_t::get_bytes( &q , 44 );
    
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
  if ( globals::skip_edf_annots )
    {
      if ( ! globals::silent ) 
	std::cerr << " forcing read as EDF\n";
      edfplus = false;
      continuous = true;
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

  rdsz = fread( p , 1, hdrSz * ns_all , file);      

  // for each of 'ns_all' signals
  
  ns = 0; // actual number of important signals

  std::vector<std::string> tlabels;
  std::set<std::string> slabels;
  
  for (int s=0;s<ns_all;s++)
    {
      
      // signal label
      std::string l = edf_t::get_string( &p , 16 );
      
      // does this exist already? if so, uniqify 
      if ( slabels.find( l ) != slabels.end() )
	{
	  int inc = 1;
	  while ( 1 ) 
	    {
	      // new unique label?
	      if ( slabels.find( l + "." + Helper::int2str( inc )  ) == slabels.end() )
		{
		  l = l + "." + Helper::int2str( inc );
		  break;
		}
	      else // keep trying
		++inc;
	    }
	}
      
      // store temporary
      tlabels.push_back( l );
      slabels.insert( l );

      // track original label position
      label_all[ l ] = s ;

    }
  
  // for each signal, does it match?
  // (and if so, change this to "standard" form)
  
  for (int s=0;s<ns_all;s++)
    {
      
      // retrieve temp label
      std::string l = tlabels[s];
      
      bool include = inp_signals == NULL || signal_list_t::match( inp_signals , &l , slabels );

      // imatch allows for case-insensitive match of 'edf annotation*'
      bool annotation = Helper::imatch( l , "EDF Annotation" ) ;

      // optionally skip all EDF annotation channels?
      if ( annotation && globals::skip_edf_annots ) 
	{	  
	  include = false;
	}

      if ( include ) 
	{
	  
	  channels.insert(s);
	  
	  annotation_channel.push_back( annotation );
	  
	  if ( annotation && ! edfplus ) 
	    {
	      //Helper::halt( "file must be EDF+ to support annotations" );
	      if ( ! globals::silent ) 
		std::cerr << " detected an annotation channel in EDF: will treat as EDF+\n";
	      edfplus = true;
	    }
	  
	  // first annotation channel is time-track
	  if ( annotation && t_track == -1 ) 
	    {
	      t_track = label.size(); 
	    }

	  // label mapping only to non-annotation channels
	  if ( ! annotation ) 
	    label2header[ l ] = label.size(); 
	  
	  label.push_back( l );	  
	  
	  ++ns;

	}
    }

  // transducer type
  for (int s=0;s<ns_all;s++)
    {
      if ( channels.find(s) != channels.end() ) 
	transducer_type.push_back( edf_t::get_string( &p , 80 ) );
      else 
	edf_t::skip( &p , 80 );
    }

  // physical dimension
  for (int s=0;s<ns_all;s++)
    {
      if ( channels.find(s) != channels.end() ) 
	phys_dimension.push_back( edf_t::get_string( &p , 8 ) );
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
  // time-track absolute offset in record
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




bool edf_record_t::read( FILE * file , int r )
{
  
  // bound checking on 'r' already done, via edf_t::read_record();
  
  // skip if already loaded?
  if ( edf->loaded( r ) ) return false;
  
  // determine offset into EDF
  long int offset = edf->header_size + (long int)(edf->record_size) * r;
  
  // find the appropriate record
  fseek( file , offset , SEEK_SET );
  
  // allocate space in the buffer for a single record, and read from file
  byte_t * p = new byte_t[ edf->record_size ];
  byte_t * p0 = p;

  // and read it
  size_t rdsz = fread( p , 1, edf->record_size , edf->file );


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

  //std::cerr << "edf_t::read_records :: scanning ... r1, r2 " << r1 << "\t" << r2 << "\n";
  
  for (int r=r1;r<=r2;r++)
    {

//       if ( ! timeline.retained(r) ) 
// 	std::cerr << "NOT retained " << r << " " << timeline.retained(r) << "\n";
      
      if ( timeline.retained(r) )
	{
	  if ( ! loaded( r ) ) 
	    {
	      edf_record_t record( this ); 
	      record.read( file , r );
	      records.insert( std::map<int,edf_record_t>::value_type( r , record ) );	      
	    }
	}
    }
  return true;
}



bool edf_t::attach( const std::string & f , 
		    const std::string & i , 
		    const std::set<std::string> * inp_signals )
{
  
  //
  // Store filename and ID
  //

  filename = f;

  id = i; 


  //
  // Attach the file
  //
  
  if ( ( file = fopen( filename.c_str() , "rb" ) ) == NULL )
    {      
      file = NULL;
      if ( ! globals::silent ) 
	std::cerr << " PROBLEM: could not open specified EDF: " << filename << "\n";
      globals::problem = true;
      return false;
    }
  
  
  //
  // Does this look like a valid EDF (i.e. at least contains a header?)
  //

  long fileSize = edf_t::get_filesize( file );
  
  if ( fileSize < 256 ) 
    {
      if ( ! globals::silent ) 
	std::cerr << " PROBLEM: corrupt EDF, file < header size (256 bytes): " << filename << "\n";
      globals::problem = true;
      return false;
    }


  //
  // Read and parse the EDF header
  //
  
  // Parse the header and extract signal codes 
  // store so we know how to read records

  inp_signals_n = header.read( file , inp_signals );


  //
  // EDF+ requires a time-track
  //
  
  if ( header.edfplus && header.time_track() == -1 ) 
    {
      if ( !header.continuous ) 
	Helper::halt( "EDF+D with no time track" );

      if ( ! globals::silent ) 
	std::cerr << " EDF+ [" << filename << "] did not contain any time-track: adding...\n";

      add_continuous_time_track();

    }

  
  //
  // Record details about byte-size of header/records
  //
  
  header_size = 256 + header.ns_all * 256;
  
  record_size = 0;
  
  for (int s=0;s<header.ns_all;s++)
    record_size += 2 * header.n_samples_all[s] ; // 2 bytes each
  

  //
  // Create timeline (relates time-points to records and vice-versa)
  // Here we assume a continuous EDF, but timeline is set up so that 
  // this need not be the case
  //

  timeline.init_timeline();


  //
  // Check remaining file size, based on header information
  //    
  

  long implied = header_size + header.nr_all * record_size;
  
  if ( fileSize != implied ) 
    {
      Helper::halt( "corrupt EDF: expecting " + Helper::int2str(implied) 
		    + " but observed " + Helper::int2str( fileSize) + " bytes" );
    }

  //
  // Output some basic information
  //

  if ( ! globals::silent ) 
    {
      std::cerr << " total duration " << Helper::timestring( timeline.total_duration_tp ) 
		<< ", with last time-point at " << Helper::timestring( ++timeline.last_time_point_tp ) << "\n";
      std::cerr << " " << header.nr_all  << " records, each of " << header.record_duration << " second(s)\n";
      std::cerr << " " << header.ns << " (of " << header.ns_all << ") signals selected ";
      std::cerr << "in " << ( header.edfplus ? "an EDF+" : "a standard EDF" ) << " file:" ;
      for (int s=0;s<header.ns;s++) 
	std::cerr << ( s % 8 == 0 ? "\n  " : " | " ) << header.label[s]; 
      std::cerr << "\n";
    }

  return true;

}



std::vector<double> edf_t::fixedrate_signal( uint64_t start , 
					     uint64_t stop , 
					     const int signal , 
					     const int downsample ,
					     std::vector<uint64_t> * tp , 
					     std::vector<int> * rec ) 
{


  std::vector<double> ret;
  
  tp->clear();

  rec->clear();

  //
  // ensure we are within bounds
  //
  
  // std::cerr << "start, stop = " << start << "\t" << stop  << "\n";
  // std::cerr << timeline.total_duration_tp << "\n";
  // std::cerr << timeline.last_time_point_tp << "\n";
  
  // Final timepoint; for EDF+ w/ gaps, then                                                                                                                
  // total_duration_tp <= last_time_point_tp + 1                                                                                                             

  uint64_t  last_time_point_tp;
  
  if ( stop >= timeline.last_time_point_tp )
    stop = timeline.last_time_point_tp;      
  
  //
  // Figure out which records are being requested 
  //
  
  //
  // we know total tp-duration of file
  // we here are given start and stop both in tp
  // 

  const uint64_t n_samples_per_record = header.n_samples[signal];
  
  // std::cerr << "signal = " << signal << "\t" << header.n_samples.size() << "\t" << header.n_samples_all.size() << "\n";
  // std::cerr << "FR " << n_samples_per_record << "\n";

  
  int start_record, stop_record;
  int start_sample, stop_sample;
  
  bool okay = timeline.interval2records( interval_t( start , stop ) , 
					 n_samples_per_record , 
					 &start_record, &start_sample , 
					 &stop_record, &stop_sample );

  
  if ( ! okay ) 
    {
      if ( ! globals::silent ) 
	std::cerr << " ** warning ... null record set found... ** \n";
      return ret; // ie empty
    }


  // std::cerr << "records start = " << start_record << " .. " << start_sample << "\n";
  // std::cerr << "records stop  = " << stop_record << " .. " << stop_sample << "\n";
  
  
  // Ensure that these records are loaded into memory
  // (if they are already, they will not be re-read)
  
  bool retval = read_records( start_record , stop_record );
  

  //
  // Copy data into a single vector
  //
 
  double bitvalue = header.bitvalue[ signal ];
  double offset   = header.offset[ signal ];

  int r = start_record;

  while ( r <= stop_record )
    {

      const edf_record_t * record = &(records.find( r )->second);

      const int start = r == start_record ? start_sample : 0 ;
      const int stop  = r == stop_record  ? stop_sample  : n_samples_per_record - 1;
      
       // std::cerr << "OUT\t"
       // 		<< record->pdata.size() << " " 
       // 		<< record->data.size() << " "
       // 		<< signal << " " 
       // 		<< header.ns << "\n";
	      
      for (int s=start;s<=stop;s+=downsample)
	{
	  // convert from digital to physical on-the-fly
	  ret.push_back( edf_record_t::dig2phys( record->data[ signal ][ s ] , bitvalue , offset ) );
	  tp->push_back( timeline.timepoint( r , s , n_samples_per_record ) );
	  rec->push_back( r );
	}
      
      r = timeline.next_record(r);
      if ( r == -1 ) break;
    }

  return ret;  
}



//
// Functions to write an EDF
//


bool edf_header_t::write( FILE * file )
{
  
  // For now, we can only write a EDF+C (i.e. standard EDF, not EDF plus)

  if ( edfplus && ! globals::silent ) 
    std::cerr << " ** warning... can only write as standard EDF (not EDF+) currently... **\n";

  // regarding the nbytes_header variable, although we don't really
  // use it, still ensure that it is properly set (i.e. we may have
  // added/removed signals, so we need to update before making the EDF
  nbytes_header = 256 + ns * 256;
  
  writestring( version , 8 , file );
  writestring( patient_id , 80 , file );
  writestring( recording_info , 80 , file );
  writestring( startdate , 8 , file );
  writestring( starttime , 8 , file );
  writestring( nbytes_header , 8 , file );
  fwrite( reserved.data() , 1 , 44 , file );
  writestring( nr , 8 , file );
  writestring( record_duration , 8 , file );
  writestring( ns , 4 , file );

  // for each of 'ns' signals
  
  for (int s=0;s<ns;s++)
    writestring( label[s], 16, file );
  
  for (int s=0;s<ns;s++)
    writestring( transducer_type[s], 80, file );

  for (int s=0;s<ns;s++)
    writestring( phys_dimension[s], 8, file );

  for (int s=0;s<ns;s++)
    writestring( physical_min[s], 8, file );

  for (int s=0;s<ns;s++)
    writestring( physical_max[s], 8, file );

  for (int s=0;s<ns;s++)
    writestring( digital_min[s], 8, file );

  for (int s=0;s<ns;s++)
    writestring( digital_max[s], 8, file );

  for (int s=0;s<ns;s++)
    writestring( prefiltering[s], 80, file );

  for (int s=0;s<ns;s++)
    writestring( n_samples[s], 8, file );
  
  for (int s=0;s<ns;s++)
    writestring( signal_reserved[s], 32, file );
  
  return true;
}


bool edf_record_t::write( FILE * file )
{

  // check if this has been read?
  
  for (int s=0;s<edf->header.ns;s++)
    {
      
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
	  for (int j=0;j< 2*nsamples;j++)
	    {	  	      
	      char a = j >= data[s].size() ? '\x00' : data[s][j];	      
	      fputc( a , file );	      
	    }
	}
    
    }
  return true;
}


bool edf_t::write( const std::string & f )
{
  
  // Note: currently, we can only write a EDF, not EDF+
  // will be easy to change, but right now, if EDF+, the time-track is
  // not properly written to disk for some reason.

  // ensure is a standard EDF

  set_edf();

  filename = f;

  FILE * outfile = NULL;

  if ( ( outfile = fopen( filename.c_str() , "wb" ) ) == NULL )      
    {
      if ( ! globals::silent ) 
	std::cerr << " ** could not open " << filename << " for writing **\n";
      return false;
    }

  header.write( outfile );

  int r = timeline.first_record();
  while ( r != -1 ) 
    {

      // we may need to load this record, before we can write it
      if ( ! loaded( r ) )
 	{
 	  edf_record_t record( this ); 
 	  record.read( file , r );
 	  records.insert( std::map<int,edf_record_t>::value_type( r , record ) );	      
 	}
      
      records.find(r)->second.write( outfile );
      r = timeline.next_record(r);
    }

  fclose(outfile);

  return true;
}


void edf_t::drop_signal( const int s )
{

  if ( s < 0 || s >= header.ns ) return;  
  --header.ns;

  // was this signal in the original EDF file? 
  bool present_in_EDF_file = header.label_all.find( header.label[s] ) != header.label_all.end() ;
  int os = present_in_EDF_file ? header.label_all[ header.label[ s ] ] : -1 ;

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

  if ( present_in_EDF_file )
    inp_signals_n.erase( inp_signals_n.find(os) );

  // need to remake label2header
  header.label2header.clear();
  for (int l=0;l<header.label.size();l++)     
    if ( header.is_data_channel(l) ) 
      header.label2header[ header.label[l] ] = l;      
  
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

void edf_t::add_signal( const std::string & label , const int Fs , const std::vector<double> & data )
{
  const int ndata = data.size();

  const int n_samples = Fs * header.record_duration ;
  
  if ( ndata == 0 ) 
    {
      if ( ! globals::silent ) 
	std::cerr << " **empty EDF, not going to add channel " << label << " **\n";
      return;
    }

  // sanity check -- ie. require that the data is an appropriate length
  if ( ndata != header.nr * n_samples ) 
    Helper::halt( "internal error: problem with length of input data" );  

  // get physical signal min/max to determine scaling
  
  double pmin = data[0];
  double pmax = data[0];
  
  for (int i=1;i<ndata;i++) 
    {
      if      ( data[i] < pmin ) pmin = data[i];
      else if ( data[i] > pmax ) pmax = data[i];
    }

  // determine bitvalue and offset

  //header
  const int16_t dmax = 32767;
  const int16_t dmin = -32768;

  double bv = ( pmax - pmin ) / (double)( dmax - dmin );
  double os = ( pmax / bv ) - dmax;

  // store (after converting to digital form)
  
  int c = 0;
  int r = timeline.first_record();

  while ( r != -1 ) 
    {

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
  
  if ( ! Helper::imatch( label , "EDF Annotation" ) )
    header.label2header[label] = header.label.size()-1;     

  header.annotation_channel.push_back( ( header.edfplus ? 
					 Helper::imatch( label , "EDF Annotation" ) :
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


  // Redundant now

  
  // create digital version of data in each record given overall scaling

//   r = timeline.first_record();
//   while ( r != -1 )
//     {
//       records.find(r)->second.calc_data( bv , os );
//       r = timeline.next_record(r);
//     }
  
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
    {
      data[signal][s] = (char)str[s];
    } 
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
	Helper::halt( "bad value of ns" );
      
      new_nsamples.push_back( new_nsamples1 );

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
	  //	  std::cerr << "i to n " << n << " becomes " << new_nsamples[s] << "\n";
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


void edf_t::reference( const signal_list_t & signals , const signal_list_t & refs )
{

  if ( signals.size() == 0 || refs.size() == 0 ) Helper::halt( "must specify signal={ ... } and ref={ ... }" );
  const int ns = signals.size();
  const int nr = refs.size();

  if ( ! globals::silent ) 
    {
      std::cerr << " referencing";
      for (int s=0;s<ns;s++) std::cerr << " " << header.label[ signals(s) ];
      std::cerr << " with respect to";
      if ( nr > 1 ) std::cerr << " the average of";
      for (int r=0;r<nr;r++) std::cerr << " " << header.label[ refs(r) ];
      std::cerr << "\n";
    }

  // check SR for all channels  
  int np = header.n_samples[ signals(0) ];

  //
  // Build reference once
  //
  
  std::vector<double> reference;

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
      for (int i=0;i<np;i++) 
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
  // transform signals one at a time, now we have reference in 'reference'
  //
  
  for (int s=0;s<signals.size();s++) 
    {
      
      if ( header.n_samples[ signals(s) ] != np ) 
	Helper::halt( "all signals/references must have similar sampling rates" );
      
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
	  
	  for (int i=0;i<np;i++) d.push_back( d0[i] - reference[cc++] );
	  
	  // next record
	  rec = timeline.next_record(rec); 
	  
	}
      
      // update signal
      update_signal( signals(s) , &d );
      
      // next signal to re-reference
    }

}


bool edf_t::populate_alist( const std::string & f )
{
  
  //
  // peek into each annotation file just to get a list of the
  // available annoations, do not load at this point
  //
  
  // Allow wildcards

  
  if ( ! Helper::fileExists( f ) ) 
    Helper::halt( "annotation file " + f + " does not exist for EDF " + filename );
  
  bool xml_mode = Helper::file_extension( f , "xml" );
  
  bool feature_list_mode = Helper::file_extension( f , "ftr" );
  
  //std::cout << "pop alist [" << f << "]\n";

  //
  // For XML files, we have to parse everything, so do this just once and store
  //
  
  if ( xml_mode ) 
    {
      annot_t::loadxml( f , this );
      return true;
    }
  
  //
  // For feature lists, load now
  //
  
  if ( feature_list_mode )
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
      
      logger  << " extracting [" << feature_name << "] for [" << id_name << "] from " << f << "\n";
      
      // create and load annotation
      
      annot_t * a = timeline.annotations.add( feature_name );
      
      // always qualifies as single TEXTUAL event
      a->name = feature_name;
      a->description = "feature-list";
      a->type.resize(1);
      a->type[0] = ATYPE_TEXTUAL;
      a->cols.resize(1);
      a->cols[0] = ".";

      // map annotations to this file
      alist[ feature_name ] = f;
      flist[f] = feature_name;

      // load features, and track how many
      aoccur[ feature_name ] = a->load_features( f );
      
      return true;
    }


  //
  // For basic text annotation files, only peek at first few lines
  //
  
  std::ifstream FIN( f.c_str(), std::ios::in );
  
  if ( FIN.bad() ) Helper::halt( "annotation file " + f + " cannot be read for EDF " + filename );
  
  while ( ! FIN.eof()  ) 
    {
      std::string line;      
      std::getline( FIN , line );      
      if ( FIN.eof() ) continue;
      if ( line == "" ) continue;
      std::vector<std::string> tok = Helper::parse( line , "\t" );
      const int n = tok.size();
      if ( n == 0 ) continue;
      if ( Helper::iequals( tok[0] , "NAME" ) )
	{	  
	  if ( tok.size() != 2 ) Helper::halt( "problem with NAME format in " + f );
	  alist[ tok[1] ] = f;
	  flist[ f ] = tok[1];
	  //std::cout << "adding a/flist <" << tok[1] << ">  <" << f << ">\n";	
	}
      // if we hit actual data, drop-out
      if ( Helper::iequals( tok[0] , "E" ) || Helper::iequals( tok[0] , "I" ) ) break;
    }
  FIN.close();
  return true;
}


int  edf_header_t::signal( const std::string & s )
{  
  signal_list_t slist = signal_list(s);
  if ( slist.size() != 1 ) 
    {
      if ( ! globals::silent ) 
	std::cerr << " ** could not find signal [" << s << "] of " << label2header.size() << " signals **\n";
      return -1;
    }  
  return slist(0);
}

bool  edf_header_t::has_signal( const std::string & s )
{
  std::vector<std::string> tok = Helper::parse( s , "|" );    
  for (int t=0;t<tok.size();t++)
    if ( label2header.find(tok[t]) != label2header.end() ) return true;
  return false;
}



signal_list_t edf_header_t::signal_list( const std::string & s )
{
  
  signal_list_t r;
  
  // wildcard means all signals '*'

  if ( s == "*" )
    {
      for (int s=0;s<label.size();s++)
	{
	  std::string lb = label[s];

	  // swap in alias?
	  if ( cmd_t::label_aliases.find( lb ) != cmd_t::label_aliases.end() ) 
	    {
	      lb = cmd_t::label_aliases[ lb ];
	      label2header[ lb ] = s;
	      label[s] = lb;
	    }
	  
	  r.add( s, lb );
	}
    }
  
  // comma-delimited; but within a signal, can have options
  // that are pipe-delimited  
  
  std::vector<std::string> tok = Helper::quoted_parse( s , "," );    
  for (int t=0;t<tok.size();t++)    
    {

      std::vector<std::string> tok2_ = Helper::quoted_parse( tok[t] , "|" );    

      // first swap in any aliases, and place those at the front of the list
      // then continue as before

      // swap in alias first? -- this may double alias, but fine.
      
      // eig. 
      // alias   sigX|sigY|sigZ
      // signal  sigY|sigZ|sig0
      // will make all --> sigX (which is correct)

      std::string alias = "";
      for (int t2=0;t2<tok2_.size();t2++)    
	{
	  if ( cmd_t::primary_alias.find( tok2_[t2] ) != cmd_t::primary_alias.end() )
	    {
	      //	      std::cout << "tok2_ " << tok2_[t2] << "and alias ["<<alias<< "]\n";
	      if ( alias == "" ) alias = tok2_[t2];
	      else if ( alias != tok2_[t2] )
		Helper::halt( "more than one alias implied" );
	      //std::cout << "tok2_ " << tok2_[t2] << "and alias ["<<alias<< "]\n";
	    }
	  else if ( cmd_t::label_aliases.find( tok2_[t2] ) != cmd_t::label_aliases.end() ) 
	    {
	      if ( alias == "" ) 
		alias = cmd_t::label_aliases[ tok2_[t2] ];	  
	      else if ( alias != cmd_t::label_aliases[ tok2_[t2] ] )
		Helper::halt( "more than one alias implied" );
	    }
	}
      
      // update list if needed
      std::vector<std::string> tok2;
      if ( alias != "" ) 
	{
	  tok2.push_back( alias );
	  const std::vector<std::string> & avec = cmd_t::primary_alias.find( alias )->second;
	  std::vector<std::string>::const_iterator aa = avec.begin();
	  while ( aa != avec.end() )
	    {
	      //std::cout << "adding " << *aa << "\n";
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
      
      // proceed as before
      for (int t2=0;t2<tok2.size();t2++)    
	{
// 	  std::cout << "t2 = " << t2 << "\t" << tok2[t2] << "\n";
// 	  std::cout << "label2header.find size " << label2header.size() << "\n";

	  // add first match found 
	  if ( label2header.find(tok2[t2]) != label2header.end() ) 
	    {

	      const int l = label2header[tok2[t2]];
	      
	      //std::cout << "found match " << l << "\n";

	      if ( t2 > 0 ) // relabel if wasn't first choice?
		{
		  label2header[ tok2[0] ] = l;
		  label[l] = tok2[0];
		}	  
	      
	      //std::cout << "adding N " << label2header[tok2[0]]  << "\n";
	      if ( added.find( label2header[tok2[0]] ) == added.end() )
		{
		  r.add( label2header[tok2[0]] , tok2[0] ); 
		  added.insert( label2header[tok2[0]] ) ;
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
  
  for (int s=0;s<label.size();s++) if ( label[s] == old_label ) label[s] = new_label;
  label_all[ new_label ] = label_all[ old_label ];
  label2header[ new_label ] = label2header[ old_label ];
  
}


double edf_header_t::sampling_freq( const int s ) const
{

  //   for (int ss=0;ss<n_samples.size();ss++)
  //     std::cout << "GET " << ss << "\t" << n_samples[s] << "\n";
  //   for (int ss=0;ss<n_samples_all.size();ss++)
  //     std::cout << "GET_ALL " << ss << "\t" << n_samples_all[s] << "\n";
  
  if ( s < 0 || s >= n_samples.size() ) return -1;
  return n_samples[ s ] / record_duration;
}

std::vector<double> edf_header_t::sampling_freq( const signal_list_t & signals ) const
{
  const int n = signals.size();
  std::vector<double> fs( n );
  for (int s=0;s<n;s++)
    {
//       std::cout << "SF " << s << " " 
// 		<< signals.signals[s] << " "
// 		<< signals(s) << " "
// 		<< n_samples[ signals.signals[s] ] << " "
// 		<< n_samples[ signals(s) ] << "\n";

      fs[s] = n_samples[ signals.signals[s] ] / record_duration;
    }
  
//   for (int jj=0;jj<n_samples.size();jj++)
//     std::cout << "SF jj " << n_samples[ jj ] << "\n";

  return fs;
}


bool edf_t::restructure()
{
  
  //
  // Map back onto original epochs
  //
  
  timeline.set_epoch_mapping();

  // output headers
  writer.var( "NR1" , "Number of records prior to restructuring" );
  writer.var( "NR2" , "Number of records after restructuring" );
  writer.var( "DUR1" , "Duration (sec) prior to restructuring" );
  writer.var( "DUR2" , "Duration (sec) after restructuring" );

  // Check that we have anything to do
  
  if ( ! timeline.is_epoch_mask_set() ) 
    {
      
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

  if ( ! globals::silent )   
    std::cerr << " restructuring as an EDF+ : ";
  
  set_edfplus();

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
  
  if ( records.size() == 0 ) globals::problem = true;
    

  if ( ! globals::silent ) 
    std::cerr << "keeping " 
	      << records.size() << " records of " 
	      << copy.size() << ", resetting mask\n";


  writer.value( "NR1" , copy.size() );
  writer.value( "NR2" , records.size() );
  
  writer.value( "DUR1" , copy.size() * header.record_duration );
  writer.value( "DUR2" , records.size() * header.record_duration );

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


void edf_t::copy_signal( const std::string & from_label , const std::string & to_label )
{
  
  const int s1 = header.signal( from_label );
  
  if ( s1 == -1 ) 
    {
      if ( ! globals::silent ) 
	std::cerr << "**error** could not find signal " << from_label << "\n";
      return;
    }

  if ( header.has_signal( to_label ) ) 
    {
      if ( ! globals::silent ) 
	std::cerr << "**error** to-signal already exists in copy_signal()\n";
      return;
    }

  //
  // get data
  //

  interval_t interval = timeline.wholetrace();  
  slice_t slice( *this , s1 , interval );
  const std::vector<double> * d = slice.pdata();
  
  //
  // add signal
  //

  add_signal( to_label , header.sampling_freq(s1) , *d );
  
}


void edf_t::update_signal( const int s , std::vector<double> * d )
{
  
  if ( header.is_annotation_channel(s) ) 
    Helper::halt( "edf_t:: internal error, cannot update an annotation channel" );
  
  // for signal s, place back data in 'd' into EDF record structure
  // and update the physical min/max

  const int points_per_record = header.n_samples[s];
  const int n = d->size();
  
  if ( n != header.nr * points_per_record )
    Helper::halt( "internal error in update_signal()" );

  // use full digital min/max scale
  const int16_t dmin = -32768;
  const int16_t dmax = 32767;  
  header.digital_min[s] = dmin;
  header.digital_max[s] = dmax;

  double pmin = (*d)[0];
  double pmax = (*d)[0];
  
  for (int i=0;i<n;i++)
    {
      if      ( (*d)[i] < pmin ) pmin = (*d)[i];
      else if ( (*d)[i] > pmax ) pmax = (*d)[i];
    }
  
  // update physical min/max (but leave orig_physical_min/max unchanged)
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
  
  int r = timeline.first_record();
  while ( r != -1 ) 
    {
      
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
	  
	  //pdata[p] = (*d)[cnt]; 
	  
	  // also need to convert to bit-value
	  // pdata[s][j] = header->bitvalue[s] * ( header->offset[s] + d );	  
	  // reverse digital --> physical scaling
	  
	  //data[p] = (*d)[cnt] / bv - os;
	  data[p] = edf_record_t::phys2dig( (*d)[cnt] , bv , os );
	
	  ++cnt;	  
	}

      r = timeline.next_record(r);
    }
}




edf_record_t::edf_record_t( edf_t * e ) 
{    

  edf = e;

  // just for now, store both digital and physical values
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


void edf_t::set_continuous()
{
  if ( ! header.edfplus ) return;
  header.continuous = true;
  header.reserved[4] = 'C';    

}

void edf_t::set_discontinuous()
{
  if ( ! header.edfplus ) return;
  header.continuous = false;
  header.reserved[4] = 'D';    
}


void edf_t::set_edfplus()
{
  if ( header.edfplus ) return;
  header.edfplus = true;
  header.continuous = true;  
  header.reserved[0] = 'E'; 
  header.reserved[1] = 'D';
  header.reserved[2] = 'F';
  header.reserved[3] = '+';
  set_continuous();
  add_continuous_time_track();
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

  set_continuous(); 
  drop_time_track();
}

void edf_t::drop_time_track()
{
  // means that the EDF will become 'continuous'
  set_continuous();

  // no TT in any case?
  if ( header.time_track() == -1 ) return;
  drop_signal( header.time_track() );
  
}


int edf_t::add_continuous_time_track()
{

  
  // this can only add a time-track to a continuous record
  // i.e. if discontinuous, it must already (by definition) 
  // have a time track

  if ( ! header.continuous ) 
    return header.time_track();

  if ( ! header.edfplus ) set_edfplus();

  // time-track already set?
  if ( header.time_track() != -1 ) return header.time_track();

  // update header
  ++header.ns;

  // set t_track channel
  header.t_track  = header.ns - 1;
  header.t_track_edf_offset = record_size; // i.e. at end of record

  const int16_t dmax = 32767;
  const int16_t dmin = -32768;
  
  // need to set a record size -- this should be enough?
  const int n_samples = globals::edf_timetrack_size;

  // how many existing 'EDF Annotations' tracks?
  int annot_tracks = 0 ; 

  std::map<std::string,int>::const_iterator jj = header.label_all.begin();
  while ( jj != header.label_all.end() )
    {
      if ( Helper::imatch( jj->first  , "EDF Annotation" ) ) 
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
  
  // create each 'TAL' timestamp, and add to record
  double dur_sec = header.record_duration;
  double onset = 0; // start at T=0

  uint64_t onset_tp = 0;
  uint64_t dur_tp = header.record_duration_tp;

  // for each record
  int r = timeline.first_record();
  
  while ( r != -1 ) 
    {

      std::string ts = "+" + Helper::dbl2str( onset ) + "\x14\x14\x00";
      
      // need to make sure that the record (i.e. other signals) 
      // are first loaded into memory...
      
      bool record_in_memory = loaded(r);
      
      if ( ! record_in_memory )
	{

	  // this will be created with ns+1 slots (i.e. 
	  // already with space for the new timetrack, 
	  // so we can add directly)
	  
	  edf_record_t record( this ); 
	  
	  record.read( file , r );

	  records.insert( std::map<int,edf_record_t>::value_type ( r , record ) );

	}

	  
      //
      // Add the time-stamp as the new track (i.e. if we write as EDF+)
      //
      
      if ( ! record_in_memory ) // record already 'updated'
	records.find(r)->second.add_annot( ts , header.t_track );
      else // push_back on end of record
	records.find(r)->second.add_annot( ts );
      
      //
      // And mark the actual record directy (i.e. if this is used in memory)
      //
      
      onset += dur_sec;
      onset_tp += dur_tp;
     
      r = timeline.next_record(r);
    }

//   std::cout << "DET1 " << records.begin()->second.pdata.size() << "\t"
// 	    << records.begin()->second.data.size() << "\n";


  return header.time_track();

}




uint64_t edf_t::timepoint_from_EDF( int r )
{

  //
  // Read this is called when constructing a time-series for 
  // an existing EDF+D, only
  //

  if ( ! header.edfplus ) Helper::halt( "should not call timepoint_from_EDF for basic EDF");
  if (   header.continuous ) Helper::halt( "should not call timepoint_from_EDF for EDF+C");
  if (   header.time_track() == -1 ) Helper::halt( "internal error: no EDF+D time-track" );
  
  // determine offset into EDF
  long int offset = header_size + (long int)(record_size) * r;

  offset += header.time_track_offset(); 

  // time-track is record : edf->header.time_track 
  // find the appropriate record
  fseek( file , offset , SEEK_SET );
  
  int ttsize = 2 * globals::edf_timetrack_size;
  
  // allocate space in the buffer for a single record, and read from file
  byte_t * p = new byte_t[ ttsize ];
  byte_t * p0 = p;

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

  // get all data
  interval_t interval = timeline.wholetrace();
  slice_t slice( *this , s , interval );
  const std::vector<double> * d = slice.pdata();
  std::vector<double> rescaled( d->size() );
  
  for (int i=0;i<d->size();i++)  rescaled[i] = - (*d)[i];

  // update signal (and min/max in header)
  update_signal( s , &rescaled );
  
}


void edf_t::rescale( const int s , const std::string & sc )
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
      if ( ! globals::silent ) 
	std::cerr << " no rescaling needed\n";
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
      if ( ! globals::silent ) 
	std::cerr << " rescaled " << header.label[s] << " to uV\n";
      header.phys_dimension[s] = "uV";     
    }
  
  if ( rescale_from_uV_to_mV || rescale_from_V_to_mV ) 
    {
      if ( ! globals::silent ) 
	std::cerr << " rescaled " << header.label[s] << " to mV\n";
      header.phys_dimension[s] = "mV";
    }
}



bool edf_t::validate( param_t & param )
{
  
  // Run through each record
  // Get min/max
  // Calculate RMS for each signal
  
  std::string signal_label = param.requires( "signal" );  
  signal_list_t signals = header.signal_list( signal_label );
  std::vector<double> Fs = header.sampling_freq( signals );
  
  bool by_epoch = timeline.epoched();

  if ( ! param.has( "epoch" ) ) by_epoch = false;
  
  const int ns = signals.size();

  bool calc_median = param.has( "median" );

  for (int s=0; s<ns; s++)
    {
               
      //
      // skip annotation channels
      //

      if ( header.is_annotation_channel( signals(s) ) ) continue;

      //
      // Whole-signal level output
      //

      interval_t interval = timeline.wholetrace();
      
      slice_t slice( *this , signals(s) , interval );
	  
      const std::vector<double> * d = slice.pdata();
      
      const int n = d->size();

      if ( n == 0 ) { continue; } 
      

      //
      // Output signal
      //
      
      std::string unit = header.phys_dimension[ signals(s) ];

      writer.level( header.label[ signals(s) ] , globals::signal_strat );
      
      //
      // Calculate whole-signal stats
      //

      double mean = MiscMath::mean( *d );
      double median = calc_median ? MiscMath::median( *d ) : 0 ;
      
      double sd   = MiscMath::sdev( *d , mean );
      double rms  = MiscMath::rms( *d );
	  
      double min = (*d)[0];
      double max = (*d)[0];
      
      for (int i = 0 ; i < n ; i++ )
	{
	  if ( (*d)[i] < min ) min = (*d)[i];
	  if ( (*d)[i] > max ) max = (*d)[i];
	}
      
      //
      // Output
      //
	  
      writer.var( "MIN"  , "Signal minimum" );
      writer.var( "MAX"  , "Signal maximum" );
      writer.var( "MEAN" , "Signal mean" );
      if ( calc_median ) 
      writer.var( "MEDIAN" , "Signal median" );
      writer.var( "SD"   , "Signal SD" );
      writer.var( "RMS"  , "Signal RMS" );

      writer.value( "MAX"  , max  );
      writer.value( "MIN"  , min  );
      
      writer.value( "MEAN" , mean );
      if ( calc_median ) 
	writer.value( "MEDIAN" , median );

      writer.value( "SD"   , sd   );
      writer.value( "RMS"  , rms  );
            

      //
      // Mean, variance, RMS, min, max
      //
            
      std::vector<double> e_mean, e_median , e_sd, e_rms;
      
      double t_min = 0 , t_max = 0;
      
      if ( ! globals::silent ) 
	std::cerr << " processing " << header.label[ signals(s) ] << " ...\n";

      
      //
      // Perform EPOCH-level analyses?
      //
      
      if ( ! by_epoch ) 
	{
	  writer.unlevel( globals::signal_strat );
	  continue;
	}

    
      //
      // Each EPOCH
      //

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
	  	 
	  double mean = MiscMath::mean( *d );
	  double median = calc_median ? MiscMath::median( *d ) : 0;
	  double sd   = MiscMath::sdev( *d , mean );
	  double rms  = MiscMath::rms( *d );
	  
	  double min = (*d)[0];
	  double max = (*d)[0];
	  
	  for (int i = 0 ; i < n ; i++ )
	    {
	      if ( (*d)[i] < min ) min = (*d)[i];
	      if ( (*d)[i] > max ) max = (*d)[i];
	    }
	  

	  //
	  // Output
	  //
	  
	  writer.epoch( timeline.display_epoch( epoch ) );

	  writer.value( "MAX"  , max  );
	  writer.value( "MIN"  , min  );

	  writer.value( "MEAN" , mean );
	  if ( calc_median ) 
	    writer.value( "MEDIAN" , median );
	  writer.value( "SD"   , sd   );
	  writer.value( "RMS"  , rms  );

	  //
	  // Record
	  //
	  
	  if ( t_min == 0 && t_max == 0 ) 
	    { 
	      t_min = min; 
	      t_max = max; 
	    } 

	  if ( min < t_min ) t_min = min;
	  if ( max > t_max ) t_max = max;
	  
	  e_mean.push_back( mean );
	  if ( calc_median ) 
	    e_median.push_back( median );
	  e_sd.push_back( sd );
	  e_rms.push_back( rms );	  
	  
	}

      // 
      // Any data?
      //

      if ( e_mean.size() == 0 ) continue;

      //
      // Get median mean/max, etc
      //
      
      const int ne = e_mean.size(); 
      double med_mean  = median_destroy( &e_mean[0] , ne );
      double med_median  = calc_median ? median_destroy( &e_median[0] , ne ) : 0 ;  
      double med_sd   = median_destroy( &e_sd[0] , ne );
      double med_rms  = median_destroy( &e_rms[0] , ne );
      
      //
      // Whole channel
      //

      writer.unepoch();
      
      writer.var( "DMIN"  , "Digital minimum from EDF header" );
      writer.var( "DMAX"  , "Digital maximum from EDF header" );
      writer.var( "UNIT"  , "Physical unit from EDF header" );
      writer.var( "PMIN"  , "Physical minimum from EDF header" );
      writer.var( "PMAX"  , "Physical maximum from EDF header" );

      writer.var( "OMIN" , "Observed physical minimum" );
      writer.var( "OMAX" , "Observed physical maximum" );

      writer.var( "MEDIAN.MEAN" , "Median of per-epoch means" );
      if ( calc_median ) 
	writer.var( "MEDIAN.MEDIAN" , "Median of per-epoch median" );

      writer.var( "MEDIAN.SD"   , "Median of per-epoch SDs" );
      writer.var( "MEDIAN.RMS"  , "Median of per-epoch RMSs" );
    	


      writer.value( "DMIN" , header.digital_min[ signals(s) ] );
      writer.value( "DMAX" , header.digital_max[ signals(s) ] );
      writer.value( "UNIT" , unit );
      writer.value( "PMIN" , header.physical_min[ signals(s) ] );
      writer.value( "PMAX" , header.physical_max[ signals(s) ] );

      writer.value( "OMIN" , t_min );
      writer.value( "OMAX" , t_max );

      writer.value( "MEDIAN.MEAN" , med_mean );
      if ( calc_median )
	writer.value( "MEDIAN.MEDIAN" , med_median );
      writer.value( "MEDIAN.SD"   , med_sd  );
      writer.value( "MEDIAN.RMS"  , med_rms );

    }

  return true;
  
}



