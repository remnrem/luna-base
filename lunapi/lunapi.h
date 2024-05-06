
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

#ifndef __LUNAPI_H__
#define __LUNAPI_H__

#include "luna.h"
#include "lunapi/rtables.h"
#include "stats/Eigen/Dense"
#include <variant>
#include <optional>
#include <memory>

struct segsrv_t;

typedef std::tuple<std::vector<std::string>,Eigen::MatrixXd> ldat_t;
typedef std::tuple<std::vector<std::string>,std::vector<Eigen::MatrixXd> > ldats_t;

typedef std::vector<std::tuple<uint64_t,uint64_t> > lint_t;
typedef std::vector<std::tuple<std::string,double,double> > lannot_t;
typedef std::vector<std::tuple<std::string,std::string,std::string,std::string,double,double> > lannot_full_t;

typedef std::variant<std::monostate, double , int , std::string, std::vector<double>, std::vector<int>, std::vector<std::string> > datum_t;
typedef std::vector<std::tuple<std::string,std::string,std::set<std::string> > > slist_t;

struct lunapi_inst_t;
typedef std::shared_ptr<lunapi_inst_t> lunapi_inst_ptr;

struct lunapi_t {

public:
  
  lunapi_t( const lunapi_t & obs ) = delete ;
  
  static lunapi_t * inaugurate()
  {
    if ( p_instance == NULL )
      {
	p_instance = new lunapi_t();
	p_instance->init();
	//std::cout << "initiated a new lunapi_t engine (" << globals::version << " " << globals::date << ")\n";	
      }	
    return p_instance;
  }

  static void retire()
  {
    if ( p_instance != NULL )
      {
	//std::cout << "retiring lunapi_t\n";
	delete p_instance;
	p_instance = NULL;
      }
  }

  static std::string version()
  {
    return globals::version + " (" + globals::date + ")";
  }
  
  
private:

  static lunapi_t * p_instance;
  
  lunapi_t() { } 

  // set up Luna library
  void init();
  
  // after importing a db
  rtables_t rtables;

public:

  //
  // Access store (e.g. after eval() or import_db()
  //

  std::vector<std::string> commands() const
  { return rtables.commands(); } 
    
  std::vector<std::pair<std::string,std::string> > strata() const
  { return rtables.list(); } 

  rtable_t table( const std::string & cmd , const std::string & faclvl ) const;

  std::vector<std::string> variables( const std::string & cmd , const std::string & faclvl ) const;
  
  rtable_return_t results( const std::string & cmd , const std::string & faclvl ) const; 

  rtables_return_t results() const;


  //
  // Import/read helper functions
  //

  static std::string cmdfile( const std::string & );

  int includefile( const std::string & );
  
  std::vector<std::string> import_db( const std::string & filename );

  std::vector<std::string> import_db( const std::string & filename , const std::set<std::string> & ids );

  std::vector<std::vector<std::string> > aliases() const;
  
  //
  // Sample list functions
  //
  
  int read_sample_list( const std::string & file );
  
  int build_sample_list( const std::vector<std::string> & toks );
  
  slist_t sample_list() const;
  
  void insert_inst( const std::string & id ,
		    const std::string & edf ,
		    const std::set<std::string> & annots );
  
  int nobs() const;
  
  void clear();

  // generators for new instance types: will always return something (even if no attached EDF/annots)
  lunapi_inst_ptr inst( const std::string & id = "id1" ) const;
  lunapi_inst_ptr inst( const std::string & id , const std::string & edf ) const;
  lunapi_inst_ptr inst( const std::string & id , const std::string & edf , const std::string & annot ) const;
  lunapi_inst_ptr inst( const std::string & id , const std::string & edf , const std::set<std::string> & annots ) const;
  
  // from sample list - may return an empty type 
  std::optional<lunapi_inst_ptr> inst( const int i ) const;

  // working w/ the sample list
  std::optional<int> get_n( const std::string & ) const;
  std::optional<std::string> get_id( const int i ) const;
  std::string get_edf( const int i ) const;
  std::set<std::string> get_annot( const int i ) const;
  
  //
  // Environment variables
  //

  void silence( const bool b);

  bool is_silenced() const;
  
  // set 
  void var( const std::string & key , const std::string & value );
  
  // get
  std::variant<std::monostate,std::string> var( const std::string & key ) const;
  std::map<std::string,std::variant<std::monostate,std::string> > vars( const std::vector<std::string> & keys ) const;
  std::map<std::string,std::string> vars( ) const;
  
  // drop
  void dropvar( const std::string & key );
  void dropvars( const std::vector<std::string> & keys );
  void dropallvars();
  
  // includes resetting all prior attributes
  void re_init();
  
  void clear_ivars(); // clear all ivars for all indivs

  // reset (global) problem/empty flags
  void reset() const;

  // flush any log bufer
  void flush();

  // desc table
  std::vector<std::vector<std::string> > desc();


  
  //
  // Command evaluation
  //

  rtables_return_t eval( const std::string & );
  
private:
  
  std::map<std::string,std::string> edfs;
  std::map<std::string,std::set<std::string> > annots;
  std::map<int,std::string> n2id;
  std::map<std::string,int> id2n;

};


struct lunapi_inst_t {

  //  ~lunapi_inst_t() { std::cout << "destructor lunapi_inst_t " << id << "\n"; }
  ~lunapi_inst_t() { }

  lunapi_inst_t(const lunapi_inst_t&) = delete;

  //
  // new instance (i.e. only from lunapi_t) w/ private ctor
  //

private:

  lunapi_inst_t( const std::string & id ) : id(id) 
  {
    state = 0;
  }

public:

  friend lunapi_t;
  friend segsrv_t;

  //
  // attach data 
  //

  bool attach_edf( const std::string & filename );

  bool attach_annot( const std::string & filename );

  //
  // drop/reset 
  //
  
  // reload an EDF
  void refresh();

  // drop an EDF, all annotations, etc
  void drop();
  
  void ivar( const std::string & key , const std::string & value ); 
  
  std::variant<std::monostate,std::string> ivar( const std::string & key ) const;

  std::map<std::string,std::variant<std::monostate,std::string> > ivars() const;

  void clear_ivar(); // clear for this indiv

  void clear_selected_ivar( const std::set<std::string> & keys ); // selected vars for this indiv
  
    
  //
  // insert signals, annotations
  //

  bool insert_signal( const std::string & label , const std::vector<double> & x , const int sr );

  bool update_signal( const std::string & label , const std::vector<double> & x );
  
  bool insert_annotation( const std::string & label , const std::vector<std::tuple<double,double> > & x , const bool durcol2 = false );
    

  //
  // basic reports
  //

  std::map<std::string,datum_t> status() const; 

  std::vector<std::string> channels();

  std::vector<bool> has_channels( const std::vector<std::string> & );
  
  std::vector<std::string> annots() const;

  std::vector<bool> has_annots( const std::vector<std::string> & );

  bool has_staging();
  
  std::vector<std::string> desc();
  
  //
  // data slices
  //
  
  lint_t epochs2intervals( const std::vector<int> & e );
  
  lint_t seconds2intervals( const std::vector<std::tuple<double,double> > & s );

  ldat_t  data( const std::vector<std::string> & chs , 
		const std::vector<std::string> & anns ,
		const bool time_track );
  
  ldat_t  slice( const lint_t & intervals , 
		 const std::vector<std::string> & chs ,
		 const std::vector<std::string> & anns ,
		 const bool time_track );
  
  ldats_t  slices( const lint_t & intervals , 
		   const std::vector<std::string> & chs ,
		   const std::vector<std::string> & anns ,
		   const bool time_track );
  
  //
  // pull annotations
  //

  lannot_t fetch_annots( const std::vector<std::string> & anns ) const;

  lannot_full_t fetch_full_annots( const std::vector<std::string> & anns ) const;

  
  //
  // Luna commands
  //

  std::string eval( const std::string & );
  std::string eval_project( const std::string & , retval_t * accumulator );
  std::string eval1( const std::string & , retval_t * accumulator );
  
  std::tuple<std::string,rtables_return_t > eval_return_data( const std::string & );
  

  //
  // last output
  //

  rtables_t rtables;

  std::vector<std::string> commands() const
  { return rtables.commands(); } 
    

  std::vector<std::pair<std::string,std::string> > strata() const
  { return rtables.list(); } 

  rtable_t table( const std::string & cmd , const std::string & faclvl ) const;

  std::vector<std::string> variables( const std::string & cmd , const std::string & faclvl ) const;
  
  rtable_return_t results( const std::string & cmd , const std::string & faclvl ) const; 

  rtables_return_t results() const;
  
  //
  // helpers
  //
    
  std::string get_id() const;

  std::string get_edf_file() const;

  std::string get_annot_files() const;

  int get_state() const;

  double last_sec() const;
  
private:

  // 0 empty; +1 attached okay, -1 problem 
  int state;
  
  // store so we can perform a refresh() if requested
  std::string id;
  
  std::string edf_filename;

  std::set<std::string> annot_filenames;
  
  // pointer to the actual data store
  edf_t edf;
  
  // helper functions
  Eigen::MatrixXd matrix_internal( const lint_t & intervals ,				   
				   const signal_list_t & signals ,
				   const std::map<std::string,int> & atype ,
				   const bool time_track );

  bool proc_channots( const std::string & chstr ,
		      const std::string & anstr ,
		      std::vector<std::string> * channels ,
		      signal_list_t * signals , 
		      std::map<std::string,int> * atype );

  
  
};

#endif
