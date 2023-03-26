
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
  
  bool empty(const std::string & s ) const;
  
  // if ! has(X) return F
  // else return yesno(value(X))
  bool yesno(const std::string & s ) const;

  std::string value( const std::string & s , const bool uppercase = false ) const;
 
  bool single() const;  

  std::string single_value() const ;

  std::string single_pair(std::string * ) const ;
  
  std::string requires( const std::string & s , const bool uppercase = false ) const;
  
  int requires_int( const std::string & s ) const;
  
  double requires_dbl( const std::string & s ) const;

  std::string dump( const std::string & indent = "  ", const std::string & delim = "\n" ) const;

  std::set<std::string> strset( const std::string & k , const std::string delim = "," , const bool uppercase = false ) const;
  
  std::vector<std::string> strvector( const std::string & k , const std::string delim = "," , const bool uppercase = false ) const;

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

  static void attach_ivars( const std::string & file );

  static void attach_idmapper( const std::string & file );

  static std::string remap_id( const std::string & id )
  {
    if ( idmapper.find( id ) == idmapper.end() ) return id;
    return idmapper.find( id )->second; 
  }
  
  // these cannot be used as variable names in scripts
  static void register_specials();

  // by default (str == NULL) read from command line arguments, 
  // then from STDIN;  otherwise, parse the commands in the 
  // text string 

  bool read( const std::string * str = NULL , 
	     bool silent = false ); 
  
  void reset() ;

  static void define_channel_type_variables( edf_t & );
  
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

  // project-level variables
  static std::map<std::string,std::string>  vars;

  // individual-specific vars
  static std::map<std::string,std::map<std::string,std::string> >  ivars;

  static std::map<std::string,int> pull_ivar( const std::vector<std::string> & ids , 
					      const std::string & phe );

  static bool pull_ivar( const std::string & id , const std::string & phe , double * x );

  static bool pull_ivar_bool( const std::string & id , const std::string & phe );
  
  // id-mapper
  static std::map<std::string,std::string>  idmapper;
  
  static std::string                        cmdline_cmds;
  static std::string                        stout_file;
  static std::string                        stout_template;
  static bool                               append_stout_file;
  static bool                               plaintext_mode;
  static std::string                        plaintext_root;
  static bool                               has_indiv_wildcard;
  static std::string resolved_outdb( const std::string & id , const std::string & str );
  
  // command-specific parameters (i.e. from command-file)
  static std::set<std::string>    signallist;
  
  static std::map<std::string,std::string> label_aliases;
  static std::map<std::string,std::vector<std::string> > primary_alias;
  static std::map<std::string,std::string> primary_upper2orig;

  static void signal_alias( const std::string & s );

  static const std::set<std::string> & signals() ;

  static void clear_signals() ;

  static std::string signal_string() ;

  static void clear_static_members() ;

  static bool is_special( const std::string & v ) 
  { return specials.find( v ) != specials.end() ; } 

 private:
  
  // a command is a single input (0+ EDFs) linked to 1 or more
  // commands, that are performed sequentially on each EDF

  std::string line; // raw input from STDIN, stored to report in case of error

  bool error;
  bool will_quit;
      
  std::vector<std::string> cmds;
  std::vector<param_t>     params;
  
  static std::set<std::string> commands;
  
  static void populate_commands() ;

  static std::set<std::string> specials;


};



//
// Primary commands
//

void proc_summaries( edf_t & , param_t & );
void proc_aliases( edf_t & , param_t & );
void proc_headers( edf_t & , param_t & );
void proc_set_headers( edf_t & , param_t & );
void proc_set_ivar( edf_t & , param_t & );
void proc_validate( edf_t & , param_t & );
void proc_desc( edf_t & , param_t & );
void proc_dump_vars( edf_t & , param_t & );
void proc_show_channel_map();
void proc_has_signals( edf_t & , param_t & );
void proc_stats( edf_t & , param_t & );
void proc_dupes( edf_t & , param_t & );

void proc_list_annots( edf_t & , param_t & );
void proc_write_annots( edf_t & , param_t & );
void proc_extend_annots( edf_t & , param_t & );
void proc_list_all_annots( edf_t & , param_t & );
void proc_list_spanning_annots( edf_t & , param_t & );
void proc_force_edf( edf_t & , param_t & );
void proc_edf_minus( edf_t & , param_t & );
void proc_set_timestamps( edf_t & , param_t & );

void proc_tag( param_t & );
void set_tag( const std::string & t = "." );
void proc_anon( edf_t & , param_t & );
void proc_write( edf_t & , param_t & );
void proc_restructure( edf_t & , param_t & );
void proc_drop_signals( edf_t & , param_t & );
void proc_rename( edf_t & , param_t & );
void proc_enforce_signals( edf_t & , param_t & );
void proc_copy_signal( edf_t & , param_t & );
void proc_order_signals( edf_t & , param_t & );
void proc_scale( edf_t & , param_t & , const std::string & s );
void proc_minmax( edf_t & , param_t & );
void proc_standardize( edf_t & , param_t & );
void proc_flip( edf_t & , param_t & );
void proc_rectify( edf_t & , param_t & );
void proc_reverse( edf_t & , param_t & );
void proc_reference( edf_t & , param_t & );
void proc_dereference( edf_t & , param_t & );
void proc_adjust( edf_t & , param_t & );
void proc_rerecord( edf_t & edf , param_t & param );
void proc_canonical( edf_t & edf , param_t & param );
void proc_remap_annots( edf_t & edf , param_t & param );

void proc_dump( edf_t & , param_t & );
void proc_dump_mask( edf_t & , param_t & );
void proc_annot_mask( edf_t & , param_t & );
//void proc_chep( edf_t & , param_t & ); // now in timeline_t
void proc_chep_mask( edf_t & , param_t & );
void proc_file_mask( edf_t & , param_t & );
void proc_file_annot( edf_t & , param_t & );
void proc_sleep_stage( edf_t & , param_t & , bool verbose = false );
void proc_suds( edf_t & , param_t & );
void proc_make_suds( edf_t & , param_t & );
void proc_self_suds( edf_t & , param_t & );
void proc_resoap( edf_t & , param_t & );
void proc_rebase_soap( edf_t & , param_t & );
void proc_place_soap( edf_t & , param_t & );

void proc_pops( edf_t & , param_t & );
void proc_eval_stages( edf_t & , param_t & );

void proc_copy_suds_cmdline();
void proc_combine_suds_cmdline();

void proc_hypoxic_burden( edf_t & , param_t & );

void proc_annotate( edf_t & , param_t & );
void proc_annot2signal( edf_t & , param_t & );
void proc_signal2annot( edf_t & , param_t & );
void proc_sig_annot_mean( edf_t & , param_t & );
void proc_annot2cache( edf_t & , param_t & );
void proc_sig_tabulate( edf_t & edf , param_t & param );

void proc_record_dump( edf_t & , param_t & );
void proc_record_table( edf_t & , param_t & );
void proc_dump_segs( edf_t & , param_t & );

void proc_intervals( param_t & , const std::string & );
void proc_align( edf_t & edf , param_t & param );
void proc_epoch_dump( edf_t & , param_t & );
void proc_epoch_matrix( edf_t & , param_t & );
void proc_head_matrix( edf_t & , param_t & );
void proc_epoch_mask( edf_t & , param_t & );
void proc_mask( edf_t & , param_t & );
void proc_eval( edf_t & , param_t & );
void proc_epoch( edf_t & , param_t & );
void proc_slice( edf_t & , param_t & , int );

void proc_freeze( edf_t & , param_t & );
void proc_thaw( edf_t & , param_t & );

void proc_trans( edf_t & , param_t & );

void proc_timetrack( edf_t & , param_t & );
void proc_continuous( edf_t & , param_t & );

void proc_covar( edf_t & , param_t & );


void proc_artifacts( edf_t & , param_t & );
void proc_rms( edf_t & , param_t & );
void proc_mse( edf_t & , param_t & );
void proc_lzw( edf_t & , param_t & );
void proc_zratio( edf_t & , param_t & );
void proc_correct( edf_t & , param_t & );

void proc_resample( edf_t & , param_t & );
void proc_zoh( edf_t & , param_t & );
void proc_moving_average( edf_t & , param_t & );
void proc_filter( edf_t & , param_t & );
void proc_filter_legacy( edf_t & , param_t & );
void proc_filter_design( edf_t & , param_t & );
void proc_filter_design_cmdline();
void proc_cwt_design( edf_t & , param_t & );
void proc_cwt_design_cmdline();
void proc_psd( edf_t & , param_t & );
void proc_fft( edf_t & , param_t & );
void proc_mtm( edf_t & , param_t & );
void proc_irasa( edf_t & , param_t & );
void proc_1overf_norm( edf_t & , param_t & );
void proc_tv_denoise( edf_t & , param_t & );
void proc_otsu( edf_t & , param_t & );
void proc_cwt( edf_t & , param_t & );
void proc_hilbert( edf_t & , param_t & );
void proc_sync(edf_t & , param_t & );
void proc_tsync(edf_t & , param_t & );
void proc_psc( edf_t & , param_t & );
void proc_microstates( edf_t & , param_t & );

void proc_asymm( edf_t & , param_t & );
void proc_tlock( edf_t & , param_t & );
void proc_peaks( edf_t & , param_t & );
void proc_zpeaks( edf_t & , param_t & );
void proc_sedf( edf_t & , param_t & );
void proc_tclst( edf_t & , param_t & );

void proc_fiplot( edf_t & , param_t & );
void proc_mi( edf_t & , param_t & );
void proc_ica( edf_t & , param_t & );
void proc_emd( edf_t & , param_t & );
void proc_dfa( edf_t & , param_t & );

void proc_attach_clocs( edf_t & , param_t & );
void proc_surface_laplacian( edf_t & , param_t & );
void proc_leave_one_out( edf_t & , param_t & );
void proc_chep_based_interpolation( edf_t & , param_t & );
void proc_correl( edf_t & , param_t & );
void proc_coh( edf_t & , param_t & );
void proc_elec_distance( edf_t & , param_t & );
void proc_acf( edf_t & , param_t & );
void proc_psi( edf_t & , param_t & );

void proc_conncoupl( edf_t & , param_t & );
void proc_pac( edf_t & , param_t & );
void proc_cfc( edf_t & , param_t & );

void proc_ecgsuppression( edf_t & , param_t & );
void proc_bpm( edf_t & , param_t & );

void proc_spindles( edf_t & , param_t & );
void proc_coupling( edf_t & , param_t & );
void proc_ripples( edf_t & , param_t & );

void proc_polarity( edf_t & , param_t & );
void proc_slowwaves( edf_t & , param_t & );
void proc_cwt( edf_t & , param_t & );
void proc_rems( edf_t & , param_t & );

void proc_dump_cache( edf_t & , param_t & );

void proc_siggen( edf_t & , param_t & );
void proc_simul( edf_t & , param_t & );
void proc_spike( edf_t & , param_t & );
void proc_shift( edf_t & , param_t & );



#endif
