
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

#ifndef __DEFS_H__
#define __DEFS_H__

#include <map>
#include <set>
#include <string>
#include <complex>
#include <stdint.h>
#include <vector>

#include "cmddefs.h"

struct param_t; 

struct sample_list_t { 
  std::string id;
  std::string edf;
  std::set<std::string> annots;
};

typedef std::complex<double> dcomp;

enum window_function_t
  { 
    WINDOW_NONE = 0 , 
    WINDOW_HAMMING,    
    WINDOW_TUKEY50,
    WINDOW_HANN  
  };

enum frequency_band_t 
  { 
    SLOW, 
    DELTA, 
    THETA, 
    ALPHA, 
    SIGMA, 
    LOW_SIGMA,
    HIGH_SIGMA,
    BETA, 
    GAMMA,
    TOTAL 
  };

enum fft_t 
{
  FFT_FORWARD,
  FFT_INVERSE
};

enum sleep_stage_t
  {
    WAKE     , 
    NREM1    , 
    NREM2    , 
    NREM3    , 
    NREM4    ,    
    REM      ,
    UNSCORED ,
    MOVEMENT ,
    ARTIFACT ,
    LIGHTS_ON ,
    UNKNOWN //  i.e. null marker / not a sleep stage
  };

typedef std::map<sleep_stage_t,std::string> sleep_stage_label_t;
typedef std::map<std::string,sleep_stage_t> sleep_stage_label_lookup_t;

typedef std::pair<double,double> freq_range_t;

struct globals
{
  
  static std::string version;
  static std::string date;

  // command definitions

  static cmddefs_t cmddefs;
  
  // global variables
  static std::map<frequency_band_t,freq_range_t> freq_band;
  static sleep_stage_label_t sleep_stage;
  
  static sleep_stage_label_lookup_t sleep_stage_labels;
  
  static char folder_delimiter;

  static std::string project_path; 

  // in -t output mode:   folder/indiv-id/{value}COMMAND-F{value}.txt{.gz}
  static std::string txt_table_prepend;
  static std::string txt_table_append;

  static bool assume_pm_starttime;
  static int assume_pm_starttime_hour;
  
  static bool force_edf;
  static bool skip_edf_annots;
  static bool skip_nonedf_annots;

  static bool set_annot_inst2hms;
  static bool set_annot_inst2hms_force;

  //
  // Annotation types stored here statically;  these can be properties of both 
  // annot_t, in which case they provide a guide for all instance_t.data fields
  // alternatively, they are specified for every a_var_t, 
  //
  
  enum atype_t { A_NULL_T ,  // i.e. not found, or for annot_t means that no one type is specified
		 A_FLAG_T , 
		 A_MASK_T , 
		 A_BOOL_T , 
		 A_INT_T , 
		 A_DBL_T , 
		 A_TXT_T , 
		 A_BOOLVEC_T , 
		 A_INTVEC_T , 
		 A_DBLVEC_T , 
		 A_TXTVEC_T };
  
  static std::map<atype_t,std::string> type_name;

  static std::map<std::string,atype_t> name_type;
  
  static bool read_ftr;

  static std::string indiv_wildcard;

  static std::string current_tag;

  // drop these EDFs
  static std::set<std::string> id_excludes;

  // only include these EDFs
  static std::set<std::string> id_includes;

  // used by the --build command only
  static std::set<std::string> sl_annot_extensions;
  static std::map<std::string,sample_list_t> sl_data;
  static bool sl_visit_edf;
  static bool sl_link_across_folders;
  
  static int sample_list_min;
  static int sample_list_max;
  static std::string sample_list_id;
  
  static bool write_naughty_list;
  static std::string naughty_list;
  
  // enforce or not the 30-second epoch check
  static bool enforce_epoch_check;

  static int default_epoch_len;

  // output common stratifier labels
  static std::string epoch_strat;
  static std::string time_strat;
  static std::string freq_strat;
  static std::string signal_strat;
  static std::string stage_strat;
  static std::string cycle_strat;
  static std::string band_strat;
  static std::string annot_strat;          // annot class
  static std::string annot_instance_strat; // annot instance ID
  static std::string annot_meta_strat;     // annot instance meta variable label
  static std::string count_strat;
  static std::string sample_strat;
  static std::string cluster_strat;

  // database variables
  static std::string & SQLITE_SCRATCH_FOLDER();  
  
  static std::string print( const freq_range_t & );
    
  // function to bail to if needed
  static void (*bail_function) ( const std::string & msg );
  
  // in CGI mode, set this to T
  static bool silent;

  // is in R mode
  static bool Rmode;
  static bool Rdisp;

  // generic global parameters
  static param_t param;

  static bool problem;

  static bool bail_on_fail;

  // global functions: primary initiation of all globals
  void init_defs();
  
  // modes
  void api();

  void R( bool );


  // default annotation folder (i.e. added to each record in sample-list implicitly)
  static std::string annot_folder;

  // additional annot files 
  static std::vector<std::string> annot_files;
  
  // specified annots
  static std::set<std::string> specified_annots;

  // helper functions to pull out global values
  static std::string band( frequency_band_t b );

  static std::string stage( sleep_stage_t );

  static std::string stage( int );

  static sleep_stage_t stage( const std::string &  );


  static double band_width( frequency_band_t b );
  
  // time track for EDF
  static std::string edf_timetrack_label;
  static int edf_timetrack_size;
  static std::string edf_annot_label;

  // time-units
  static uint64_t tp_1sec;
  static uint64_t tp_1000thsec;
  static double tp_duration;   
};

#endif
