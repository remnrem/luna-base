
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


typedef std::tuple<Eigen::MatrixXd,std::vector<std::string> > ldat_t;
typedef std::tuple<std::vector<Eigen::MatrixXd>,std::vector<std::string> > ldats_t;
typedef std::vector<std::tuple<uint64_t,uint64_t> > lint_t;
typedef std::variant<std::monostate, double , int , std::string, std::vector<double>, std::vector<int>, std::vector<std::string> > datum_t;

struct lunapi_t {

  //
  // mandatory initialization
  //
  
  static void init(); 

  //
  // global variables 
  //

  // set 
  static void var( const std::string & key , const std::string & value );
  
  // get
  static std::variant<std::monostate,std::string> var( const std::string & key );
  static std::map<std::string,std::variant<std::monostate,std::string> > vars( const std::vector<std::string> & keys );

  // drop
  static void dropvar( const std::string & key );
  static void dropvars( const std::vector<std::string> & keys );


  //
  // new instance
  //

  lunapi_t( const std::string & id ) : id(id) 
  {
    reset();  
  }

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
  
  // clear all variables
  void clear();


  //
  // individual variables
  //
  
  /* void ivar( const std::string & key , const std::string & value ); */
  
  /* std::string ivar( const std::string & key ) const; */


  //
  // basic reports
  //

  std::map<std::string,datum_t> status() const; 

  std::vector<std::string> channels();

  std::vector<std::string> annots() const;

  //
  // data slices
  //
  
  lint_t epochs2intervals( const std::vector<int> & e );
  
  lint_t seconds2intervals( const std::vector<std::tuple<double,double> > & s );
  
  ldat_t  slice( const lint_t & intervals , 
		 const std::string & chstr ,
		 const std::string & anstr );
  
  ldats_t  slices( const lint_t & intervals , 
		   const std::string & chstr ,
		   const std::string & anstr );
    
  
  //
  // Luna commands
  //
  
  bool eval( const std::string & );


  //
  // last output
  //

  rtables_t rtables;

  std::vector<std::string> commands() const
  { return rtables.commands(); } 
    

  std::vector<std::pair<std::string,std::string> > strata() const
  { return rtables.list(); } 

  rtable_t table( const std::string & cmd , const std::string & faclvl ) const;

  rtable_data_t data( const std::string & cmd , const std::string & faclvl ) const; 
  

  //
  // helpers
  //
  
  // reset problem/empty flags
  void reset();
  
  bool read_db( const std::string & filename );
  
  
private:

  // 0 empty; +1 attached okay, -1 problem 
  int state;
  
  // store so we can perform a refresh() if requested
  std::string id;
  
  std::string edf_filename;

  std::set<std::string> annot_filenames;

  // the actual data store
  edf_t edf;

  
  // helper functions
  Eigen::MatrixXd matrix_internal( const lint_t & intervals ,				   
				   const signal_list_t & signals ,
				   const std::map<std::string,int> & atype );

  bool proc_channots( const std::string & chstr ,
		      const std::string & anstr ,
		      std::vector<std::string> * channels ,
		      signal_list_t * signals , 
		      std::map<std::string,int> * atype );

  
  
};

#endif
