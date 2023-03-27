
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

#ifndef __PSGLIB_EDF_H__
#define __PSGLIB_EDF_H__

#include "timeline/timeline.h"
#include "clocs/clocs.h"
#include "tal.h"
#include "edfz/edfz.h"
#include "edf/signal-list.h"

#include <iostream>
#include <vector>
#include <fstream>
#include <map>
#include <set>
#include <stdint.h>

typedef unsigned char byte_t;

struct param_t;

struct edf_t;

struct cansigs_t;

struct edf_header_t
{
  
  // ASCII values
  
  // char[8]
  std::string    version;
  
  // char[80]
  std::string    patient_id;
  
  // char[80]
  std::string    recording_info;
  
  // dd mm yy    char[8]
  std::string    startdate;
  
  // hh mm ss    char[8]
  std::string    starttime;
  
  // char[8]
  int            nbytes_header;
  
  // Reserved field (used to indicate EDF+ status)
  std::vector<char>    reserved;  
  
  // Number of records (-1 if unknown) | char[8]
  int            nr;     // max in memory, possibly after restructuring (i.e. nr <= nr_all) 
  int            nr_all; // total in file

  // Duration in seconds in EDF / in tp (time-point) intervals internally| char[8]
  double          record_duration;
  uint64_t        record_duration_tp;
  
  // Number of signals | char[8]
  int            ns;       // number finally selected and input
  int            ns_all;   // total number in file
  
  
  //
  // For each of 'ns' signals (i.e. selected ones only)
  //

  std::vector<std::string> label;           // e.g. actual (case-sensitive) label EEG or BodyTemp
  std::map<std::string,int> label_all;      // still need to track orig. for read() records ; use UPPER

  std::vector<std::string> transducer_type; // e.g. AgAgCl electrode
  
  std::vector<std::string> phys_dimension;  // e.g. uV or degreeC
  
  std::vector<double>      physical_min;    // can be changed (i.e. after update_signal() )
  std::vector<double>      physical_max; 
  
  std::vector<double>      orig_physical_min;  // always refer to EDF
  std::vector<double>      orig_physical_max; 

  std::vector<int>         digital_min;  // can be changed after update_signal
  std::vector<int>         digital_max; 

  std::vector<int>         orig_digital_min;  // always refer to EDF
  std::vector<int>         orig_digital_max; 

  std::vector<std::string> prefiltering;
  
  // for each record, for each signal
  std::vector<int>         n_samples;  
  std::vector<int>         n_samples_all; // if skipping signals, still need to track size  
  
  std::vector<std::string> signal_reserved;

  // derived values
  std::vector<double>      bitvalue;

  std::vector<double>      offset;

  // signal label --> edf_header_t slot 
  // i.e. may include annotations

  std::map<std::string,int> label2header;
  
  // EDF+ -- annotation versus data channel
  std::vector<bool>        annotation_channel;
  
  // EDF+ -- time-track annotation channel (or -1)
  int t_track;  // track in loaded (selected) internal space
  int t_track_edf_offset; // track in original EDF space 

  // track alias swaps
  std::map<std::string,std::string> aliasing;
  
  edf_header_t() 
  {
    init();
  }
  
  bool continuous;
  
  bool edfplus;

  void init()
  {
    edfplus    = false;
    continuous = true;
    // original: reserved area filled with NULL
    // reserved.resize( 44 , '\0' );
    reserved.resize( 44 , ' ' );
    ns = ns_all = 0;
    t_track = -1;
    
    label.clear();
    label_all.clear();
    label2header.clear();    
    transducer_type.clear(); 
    phys_dimension.clear(); 
    physical_min.clear();
    physical_max.clear();
    digital_min.clear();
    digital_max.clear();
    orig_physical_min.clear();
    orig_physical_max.clear();
    orig_digital_min.clear();
    orig_digital_max.clear();
    prefiltering.clear();
    n_samples.clear();
    n_samples_all.clear();
    signal_reserved.clear();
    bitvalue.clear();
    offset.clear();
    annotation_channel.clear();
  }

  std::string summary() const;

  std::set<int> read( FILE * file, edfz_t * edfz , const std::set<std::string> * inp_signals );

  bool write( FILE * file , const std::vector<int> & ch2slot );
  
  bool write( edfz_t * edfz , const std::vector<int> & ch2slot );

  int  signal( const std::string & s , bool silent = false );

  bool  has_signal( const std::string & s );

  int  original_signal( const std::string & s );

  int  original_signal_no_aliasing( const std::string & s );

  signal_list_t signal_list( const std::string & s , bool no_annotation_channels = false , bool show_warnings = true );
  
  void signal_alias( const std::string & s );

  double sampling_freq( const int s ) const;
  
  bool is_data_channel( const int s ) const
  {
    if ( s < 0 || s > ns ) return false;
    return !annotation_channel[s];
  }

  bool is_annotation_channel( const int s ) const 
  {
    if ( s < 0 || s > ns ) return false;
    return annotation_channel[s];
  }
  
  int time_track() const { return t_track; }
  int time_track_offset() const { return t_track_edf_offset; }

  std::vector<double> sampling_freq( const signal_list_t & ) const;
  
  void rename_channel( const std::string & o , const std::string & n );

  void check_channels();
  
  void drop_annots_from_signal_list( signal_list_t * slist )
  {
    
    std::vector<int> signals;  
    std::vector<std::string> signal_labels;  
    
    for (int s=0;s<slist->size();s++)
      {
	if ( is_annotation_channel( (*slist)(s) ) ) continue;
	signals.push_back( (*slist)(s) );
	signal_labels.push_back( slist->label(s) );
      }
    
    slist->signals = signals;
    slist->signal_labels = signal_labels;    
  }

 
};



struct edf_record_t
{

  friend struct edf_header_t;
  friend struct edf_t;
  
  //
  // contains all samples for all signals for a single
  // record/time-interval
  //

 public:
  
  edf_record_t( edf_t * e );

  // 
  // Main I/O functions
  //

  // from either EDF or EDFZ, it will determine given edf_t * parent
  bool read( int r ); 

  // for writing, split out into two separate functions (no particular reason for the differences...)
  bool write( FILE * file , const std::vector<int> & ch2slot );

  bool write( edfz_t * , const std::vector<int> & ch2slot );

  
  void add_data( const std::vector<int16_t> & );
  
  std::vector<double> get_pdata( const int signal );
  
  // here we know which slot to add to
  void add_annot( const std::string & , int signal );

  // push back on end
  void add_annot( const std::string & );

  // now redundant
  //  void calc_data( double bitvalue , double offset  );

  void drop( const int s );
  
 private:


  //
  // Primary data store ( signal x samples per record )
  //
  
  edf_t * edf;
  
  std::vector<std::vector<int16_t> >    data;
  
  // compute this on-the-fly
  //  std::vector<std::vector<double> > pdata;  // physically-scaled
  
 public:

  //
  // Helper functions
  //
    
  // 2-byte two's complement integers (little-endian format)
  
  inline static int16_t tc2dec( char a , char b );
  inline static void dec2tc( int16_t x , char * a , char * b );

  // convert between physical (double) and digitial (int16_t) scales
  // directly give bit-value and offset
  inline static double dig2phys( int16_t , double , double );
  inline static int16_t phys2dig( double , double , double );
  
};





struct edf_t
{
  
  friend struct edf_header_t;
  friend struct edf_record_t;

public:

  //
  // Core data
  //

  std::string                filename;

  std::string                id;

  std::vector<std::string>   annot_files; // for WRITE --> sample-list

  edf_header_t               header;
  
  std::map<int,edf_record_t> records;

  std::set<int>              inp_signals_n; // read these signals
  
  int                        record_size;   // bytes per record (for ns_all signals)
  int                        header_size;   // bytes for entire header

  //
  // EDF+ Annotations store (populate on reading, but only used
  // for WRITE edfz , i.e. to store these in the index)
  //

  std::map<int,std::string> edf_annots;
  bool has_edf_annots;
  
  //
  // Primary timeline for masking, annotations, etc
  //

  timeline_t timeline;

  uint64_t timepoint_from_EDF( int r );


  //
  // Channel location data 
  //

  clocs_t clocs;


  //
  // EDFZ
  //
  
  edfz_t * edfz_ptr() const
  {
    return edfz;
  }

  //
  // Annotations
  //

/*   std::map<std::string,std::string> alist; // map annoations -> annotation files (but do not load) */
/*   std::map<std::string,std::string> flist; // annotation files -> annoation name  */

  std::map<std::string,int> aoccur;        // map annoations -> # of occurences

  //
  // Data access
  //

  int records_loaded() const { return records.size(); } 

  
  // has this record already been loaded?

  bool loaded( const int r ) const { return records.find(r) != records.end(); } 

  // load if not loaded
  void ensure_loaded( const int rec )
  {
    // we may need to load this record first, before we can edit it
    if ( ! loaded( rec ) )
      {
	edf_record_t record( this ); 
	record.read( rec );
	records.insert( std::map<int,edf_record_t>::value_type( rec , record ) );	      
      }
  }

  
  std::vector<double> fixedrate_signal( uint64_t start , 
					uint64_t stop , 
					const int signal , 
					const int downsample , 
					std::vector<uint64_t> * tp , 
					std::vector<int> * rec , 
					std::vector<int16_t> * ddata = NULL );
  
  
  tal_t tal( const int signal , const int rec );

  //
  // Manipulate signal data
  //

  void swap_in_aliases();
  
  void drop_signal( const int s );

  void add_signal( const std::string & label , const int n_samples , const std::vector<double> & data ,
		   double pmin = 0 , double pmax = 0 ,
		   int16_t dmin = 0 , int16_t dmax = 0 );
  
  void copy_signal( const std::string & from_label , const std::string & to_label );

  void set_order( param_t & param );

  void set_headers( param_t & param );
  
  void update_signal_retain_range( int s , const std::vector<double> * );

  void update_signal( int s , const std::vector<double> * , int16_t * dmin = NULL , int16_t * dmax = NULL , double * pmin = NULL , double * pmax = NULL );

  void shift( int s , int tp , bool wrap = true ); 

  void update_records( int a , int b , int s , const std::vector<double> * );

  void data_dumper( const std::string & , const param_t & );
  
  void tabulate( param_t & param );

  void seg_dumper( param_t & param );
  
  void record_dumper( param_t & param );
  
  void record_table( param_t & param );

  void data_epoch_dumper( param_t & param , std::set<std::string> * = NULL );

  void epoch_matrix_dumper( param_t & param );

  void head_matrix_dumper( param_t & param );

  // redundant
  void reference_and_scale( const int s , const int r , double rescale = 1 ); // perform single channel referencing

  void pairwise_reference( const signal_list_t & signals ,
			   const signal_list_t & refs ,
			   const bool make_new ,
			   const std::vector<std::string> & new_channels ,
			   const int new_channel_sr , // ignored if not make_new
			   const bool dereference = false ,
			   const bool verbose = true );

  void reference( const signal_list_t & signals ,
		  const signal_list_t & refs ,
		  const bool make_new ,
		  const std::string & new_channel ,
		  const int new_channel_sr , // ignored if not make_new
		  const bool dereference = false ,
		  const bool verbose = true );

  
  void guess_canonicals( param_t & , bool make_signals );

  cansigs_t make_canonicals( const std::vector<std::string> & files,
			     const std::string &  group ,
			     const bool make_signals ,
			     const bool drop_originals , 
			     const std::string & prefix ,
			     const std::set<std::string> * cs = NULL ,
			     const bool only_check_labels = false );
    
  void flip( const int s );

  void reverse( const int s ); 

  void minmax( signal_list_t & );

  void reset_record_size( const double );

  void set_id( const std::string & s ) { id=s; } 



  //
  // Basic stats
  //

  void covar( const std::string & s1, const std::string & s2 );
	      
  //
  // Header reporting 
  //

  void terse_summary( param_t & param );

  void report_aliases() const;
  
  void description( const param_t & );
  

  //
  // Units
  //

  void rescale( const int s , const std::string & sc , const bool quietly = false );

private:
  
  //
  // File buffer for standard EDF
  //
  
  FILE * file;


  //
  // Alternate buffer for EDFZ
  //

  edfz_t * edfz;
  
  
  //
  // Endianness
  // 

  enum endian_t {     
    MACHINE_LITTLE_ENDIAN = 0 ,
    MACHINE_BIG_ENDIAN = 1 };
  
  static endian_t endian;

  endian_t determine_endian() 
  {
    int i = 1;
    char *p = (char *)&i;
    return p[0] == 1 ? MACHINE_LITTLE_ENDIAN : MACHINE_BIG_ENDIAN ;
  }

  

public:
  

  edf_t() : timeline( this )
  {
    endian = determine_endian();    
    file = NULL;
    edfz = NULL;
    init();
  } 


  ~edf_t() 
  {
    init();
  }

  void init()
  {
    if ( file != NULL ) 
      fclose(file);
    file = NULL;

    if ( edfz != NULL ) 
      {
	edfz->close();
	delete edfz;
      }
    edfz = NULL;
    
    header.init();
    records.clear();    
    inp_signals_n.clear();
    has_edf_annots = false;
  }
  
  
  //
  // Primary read modes
  //

  bool attach( const std::string & f , const std::string & id , const std::set<std::string> * inp_signals = NULL , const bool silent = false );
  
  bool read_records( int r , int r2 );

  bool basic_stats( param_t & );

  //
  // Create EDF struct from ASCII
  //
  
  bool read_from_ascii( const std::string & f , // filename
			const std::string & i , // id
			const int Fs , // fixed Fs for all signals
			const std::vector<std::string> & labels0 , // if null, look for header
			const std::string & startdate = "01.01.00", 
			const std::string & starttime = "00.00.00" 
			) ;

  //
  // Initiate empty (but fixed duration) EDF
  //

  bool init_empty( const std::string & id , const int nr , const int rs ,
		   const std::string & startdate ,
		   const std::string & starttime );
  
  //
  // Utility to append to an EDF (static function)
  //

  
  static bool append( const std::string & filename ,
		      const std::vector<std::string> & channels ,
		      const std::vector<std::vector<std::vector<double> > > & data );


  //
  // EDF vs EDF+
  //

  void set_edfplus();

  void set_edf();
  
  void reset_start_time();

  void set_continuous();
  
  void set_discontinuous();

  bool is_actually_discontinuous(); // asks 'but is this /really/ discontinuous?'...

  bool is_actually_standard_edf(); // i.e. effectively continuous AND no EDF Annots (other than time-track)
  
  int add_time_track( const std::vector<uint64_t> * tps = NULL ); // add EDF+C (or EDF+D) track to the EDF (as Annot channel)
  
  void drop_time_track();

  void drop_annots();

  //
  // Annotations
  //

  bool load_annotations( const std::string & f );
  

  //
  // Write EDF(Z) back to file
  //
  
  bool write( const std::string & f , 
	      bool edfz = false , 
	      int null_starttime = 0 , 
	      bool always_edfd = false , 
	      const std::vector<int> * p_ch2slot = NULL );
  
    
  // given a mask, change the representation in memory, 
  // [downstream, we can add new things, like changing the record structure, etc]
  // needs to update header as well (and timeline)
  
  bool restructure();


  //
  // for COPY of edf_ts: read up all records, i.e. make entirely in-memory
  //  
  //  int edf_t::read_all();

  //
  // after copying shallow copy of an edf_t, need to fix things:
  //

  void update_edf_pointers( edf_t * );


  
  //
  // extract and realign samples based on a set of annotations
  //
  
  bool align( const std::vector<std::string> & annots ); 


  //
  // edf-minus & set-timestamps
  //

  bool edf_minus();

  void set_timestamps( param_t & );
  
  // empirical recalculate the physical min/max from data (i.e. 
  // instead of relying on the EDF header;  replace the header values
  
  void update_physical_minmax(int);

  // given a list of (non-epoch) intervals, slice the EDF
  
  void slicer( const std::set<interval_t> & , param_t & , int );
  

  //
  // Helper functions
  //

  static uint64_t    get_filesize(FILE *file);
  static int         get_int( byte_t ** p , int sz );
  static double      get_double( byte_t ** p , int sz );
  static void        skip( byte_t ** p , int sz );
  static std::string get_string( byte_t ** p , int sz );
  static std::vector<char> get_bytes( byte_t ** p , int sz );

  
};


// misc.

void dump_intervals( const std::string & edfs , 
		     const std::string & ints );



#endif
