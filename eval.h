
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

#ifndef __LUNA_EVAL_H__
#define __LUNA_EVAL_H__

#include <string>
#include <map>
#include <set>
#include <vector>
#include <sstream>

struct edf_t;


//
// Helper to parse command syntax
//

struct param_t
{

 public:
  
  void add( const std::string & option , const std::string & value = "" ); 

  void add_hidden( const std::string & option , const std::string & value = "" );

  int size() const;
  
  void parse( const std::string & s );
  
  void update( const std::string & id , const std::string & wc );

  void clear();
  
  bool has(const std::string & s ) const;

  std::string value( const std::string & s ) const;
 
  bool single() const;  

  std::string single_value() const ;
  
  std::string requires( const std::string & s ) const;
  
  int requires_int( const std::string & s ) const;
  
  double requires_dbl( const std::string & s ) const;

  std::string dump( const std::string & indent = "  ", const std::string & delim = "\n" ) const;

  std::set<std::string> strset( const std::string & k , const std::string delim = "," ) const;
  
  std::vector<std::string> strvector( const std::string & k , const std::string delim = "," ) const;

  std::vector<double> dblvector( const std::string & k , const std::string delim = "," ) const;

  std::vector<int> intvector( const std::string & k , const std::string delim = "," ) const;

  std::set<std::string> keys() const;

private:

  std::map<std::string,std::string> opt;
  std::set<std::string> hidden;


};


class cmd_t 
{
  
 public:
  
  cmd_t();

  cmd_t( const std::string & str );

  static void add_cmdline_cmd( const std::string & c ) ;
  
  static void parse_special( const std::string & s , const std::string & delim );
  
  // by default (str == NULL) read from command line arguments, 
  // then from STDIN;  otherwise, parse the commands in the 
  // text string 

  bool read( const std::string * str = NULL , 
	     bool silent = false ); 
  
  void reset() ;

  void replace_wildcards( const std::string & id );

  bool eval( edf_t & ) ;
  
  bool empty() const ;

  bool valid() const ;

  bool badline() const ;

  std::string offending() const ;
  
  int num_cmds() const ;
  
  std::string cmd(const int i) ;
    
  param_t & param(const int i) ;
  
  bool process_edfs() const ;
  
  bool is( const int n , const std::string & s ) const;
  
  std::string data() const ;

  bool quit() const ;

  void quit(bool b) ;

  //
  // Static members (i.e. from command line)
  //

  static std::string                        input;  
  static std::map<std::string,std::string>  vars;

  static std::string                        cmdline_cmds;
  static std::string                        stout_file;
  static bool                               append_stout_file;

  // command-specific parameters (i.e. from command-file)
  static std::set<std::string>    signallist;
  
  static std::map<std::string,std::string> label_aliases;
  static std::map<std::string,std::vector<std::string> > primary_alias;

  static void signal_alias( const std::string & s );

  static const std::set<std::string> & signals() ;

  static void clear_signals() ;

  static std::string signal_string() ;

  static void clear_static_members() ;

 private:
  
  // a command is a single input (0+ EDFs) linked to 1 or more
  // commands, that are performed sequentially on each EDF

  std::string line; // raw input from STDIN, stored to report in case of error

  bool error;
  bool will_quit;
      
  std::vector<std::string> cmds;
  std::vector<param_t>     params;
  std::vector<param_t>     original_params;
  
  static std::set<std::string> commands;
  
  static void populate_commands() ;

};



//
// Primary commands
//

void proc_summaries( edf_t & , param_t & );
void proc_headers( edf_t & , param_t & );
void proc_validate( edf_t & , param_t & );
void proc_desc( edf_t & , param_t & );

void proc_stats( edf_t & , param_t & );
void proc_list_annots( edf_t & , param_t & );
void proc_list_all_annots( edf_t & , param_t & );

void proc_tag( param_t & );
void set_tag( const std::string & t = "." );
void proc_anon( edf_t & , param_t & );
void proc_write( edf_t & , param_t & );
void proc_restructure( edf_t & , param_t & );
void proc_drop_signals( edf_t & , param_t & );
void proc_copy_signal( edf_t & , param_t & );
void proc_scale( edf_t & , param_t & , const std::string & s );
void proc_flip( edf_t & , param_t & );
void proc_reference( edf_t & , param_t & );
void proc_rerecord( edf_t & edf , param_t & param );

void proc_dump( edf_t & , param_t & );
void proc_dump_mask( edf_t & , param_t & );
void proc_file_mask( edf_t & , param_t & );
void proc_file_annot( edf_t & , param_t & );
void proc_sleep_stage( edf_t & , param_t & , bool verbose = false );

void proc_record_dump( edf_t & , param_t & );
void proc_intervals( param_t & , const std::string & );
void proc_epoch_dump( edf_t & , param_t & );
void proc_epoch_matrix( edf_t & , param_t & );
void proc_epoch_mask( edf_t & , param_t & );
void proc_mask( edf_t & , param_t & );
void proc_eval( edf_t & , param_t & );
void proc_epoch( edf_t & , param_t & );
void proc_slice( edf_t & , param_t & , int );

void proc_timetrack( edf_t & , param_t & );
void proc_continuous( edf_t & , param_t & );

void proc_covar( edf_t & , param_t & );

void proc_artifacts( edf_t & , param_t & );
void proc_rms( edf_t & , param_t & );
void proc_mse( edf_t & , param_t & );
void proc_lzw( edf_t & , param_t & );
void proc_zratio( edf_t & , param_t & );

void proc_resample( edf_t & , param_t & );
void proc_filter( edf_t & , param_t & );
void proc_filter_legacy( edf_t & , param_t & );
void proc_filter_design( edf_t & , param_t & );
void proc_filter_design_cmdline();
void proc_cwt_design( edf_t & , param_t & );
void proc_cwt_design_cmdline();
void proc_psd( edf_t & , param_t & );
void proc_mtm( edf_t & , param_t & );
void proc_1overf_norm( edf_t & , param_t & );
void proc_tv_denoise( edf_t & , param_t & );
void proc_cwt( edf_t & , param_t & );
void proc_hilbert( edf_t & , param_t & );


void proc_fiplot( edf_t & , param_t & );
void proc_mi( edf_t & , param_t & );
void proc_ica( edf_t & , param_t & );
void proc_emd( edf_t & , param_t & );
void proc_leave_one_out( edf_t & , param_t & );
void proc_correl( edf_t & , param_t & );
void proc_coh( edf_t & , param_t & );
void proc_coh_legacy( edf_t & , param_t & );
void proc_elec_distance( edf_t & , param_t & );

void proc_pac( edf_t & , param_t & );
void proc_cfc( edf_t & , param_t & );

void proc_ecgsuppression( edf_t & , param_t & );
void proc_bpm( edf_t & , param_t & );

void proc_spindles( edf_t & , param_t & );
void proc_polarity( edf_t & , param_t & );
void proc_slowwaves( edf_t & , param_t & );
void proc_cwt( edf_t & , param_t & );
void proc_spike( edf_t & , param_t & );



#endif
