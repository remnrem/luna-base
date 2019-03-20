
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


#ifndef __ANNOT_H__
#define __ANNOT_H__

#include <vector>
#include <map>
#include <string>
#include <map>
#include <set>
#include <iostream>

// a single 'annotation' (that has to be attached to a 'timeline' and
// therefore a single EDF)

// has one or more events, where an event is localized in time

// events can be loaded as 'epochs' but are internally converted to the 
// interval scale; i.e. we'll only ever have a single 'epoch' definition at
// any one time, and that is in timeline_t

// events can have artibitrary real and textual key-value pairs attached

#include "intervals/intervals.h"
#include "miscmath/miscmath.h"
#include "helper/helper.h"


// annotation_set_t  is the top-level set of annotations

// annot_t    a single coherent annot class (e.g. one FTR)
//    has interval_t --> set of avar_t map 

// a_instance_t is an instance of an annotation 
//   has an id (may or may not be unique) 
//   has map of var -> a_data_t 

//   annot_t : interval_t->ainst_t : set of avar_t : 

//
// for one instance of an annotation, a particular instance is uniquely defined by an 
// interval and an ID
//

struct instance_idx_t;
struct instance_t;
struct avar_t;
struct edf_t;

typedef std::map<instance_idx_t,instance_t*> annot_map_t;
typedef std::map<std::string,avar_t*> instance_table_t;

struct annot_t
{
  
  friend struct timeline_t;

  friend struct annotation_set_t;
  
  

  //
  // Main data members
  //

  std::string name;
  
  // here an annot_t can have a 'type';  in practice, the instances can 
  // have artbitrary meta-data of any type;  however, based on the definition
  // of the annot class, this type will be useful to help clients know what to 
  // expect for this annotation. 

  globals::atype_t type;

  std::string file;

  std::string description;

  // info on types must be specified, so we know how to encode an
  // event (instance of the annotation with fields)
  
  std::map<std::string, globals::atype_t> types;
  
  // main data store
  
  annot_map_t interval_events;
  
  // for clean-up

  std::set<instance_t*> all_instances;

  

  //
  // Constructor/destructor
  //

  annot_t( const std::string & n )  : name(n) 
  { 
    file = description = "";
    type = globals::A_NULL_T;
    types.clear();
  }
  

 ~annot_t() 
   {
     wipe();
   }


  //
  // Primary loaders/savers
  //

  static bool load( const std::string & , edf_t & edf );  

  static bool map_epoch_annotations( edf_t & parent_edf , 
				     const std::vector<std::string> & e , 
				     const std::string & filename , 
				     uint64_t elen , 
				     uint64_t einc );
  
  int  load_features( const std::string & );    
  
  static bool loadxml( const std::string & , edf_t * );

  bool save( const std::string & );  
  
  bool savexml( const std::string & );  
  
  static void dumpxml( const std::string & , bool );

  int num_interval_events() const { return interval_events.size(); } 

  
  //
  // key function to add a new instance to this annotation
  //
  
  instance_t * add( const std::string & id , const interval_t & interval );

  void remove( const std::string & id , const interval_t & interval );
   
  uint64_t minimum_tp() const;  
  
  uint64_t maximum_tp() const;


  //
  // main output function: return pointers to all annotations in a given window
  //
  
  annot_map_t extract( const interval_t & window );
  
  
  std::set<std::string> instance_ids() const;

 private:

  void wipe();
  
  void reset()
  {        
    name = "";
    file = "";
    description = "";
    types.clear();
    interval_events.clear();
    wipe();
  }


  // helper functions 
  
public:

  static std::vector<bool> as_bool_vec( const std::vector<int> & x ) ; 
  static std::vector<bool> as_bool_vec( const std::vector<double> & x ) ;
  static std::vector<bool> as_bool_vec( const std::vector<std::string> & x ) ; 

  static std::vector<int> as_int_vec( const std::vector<bool> & x ) ;
  static std::vector<int> as_int_vec( const std::vector<double> & x ) ;
  static std::vector<int> as_int_vec( const std::vector<std::string> & x ) ;

  static std::vector<double> as_dbl_vec( const std::vector<bool> & x ) ;
  static std::vector<double> as_dbl_vec( const std::vector<int> & x ) ;
  static std::vector<double> as_dbl_vec( const std::vector<std::string> & x ) ;

  static std::vector<std::string> as_txt_vec( const std::vector<bool> & x ) ;
  static std::vector<std::string> as_txt_vec( const std::vector<int> & x ) ;
  static std::vector<std::string> as_txt_vec( const std::vector<double> & x ) ;


};





struct instance_idx_t { 
  
  instance_idx_t( const annot_t * p , 
		  const interval_t & i , 
		  const std::string & s ) 
  : parent(p) , interval(i) , id(s) { } 
  
  const annot_t *   parent;

  interval_t  interval;
  
  std::string id;

  bool operator< ( const instance_idx_t & rhs ) const 
  {
    if ( interval < rhs.interval ) return true;
    if ( interval > rhs.interval ) return false;
    if ( parent->name < rhs.parent->name ) return true;
    if ( parent->name > rhs.parent->name ) return false;
    return id < rhs.id;
  }

};



struct instance_t {   
  
  // an instance then has 0 or more variable/value pairs
  
  std::map<std::string,avar_t*> data;
  
  // for clean-up
  
  std::set<avar_t*> tracker;


  //
  // In/out functions
  //

  //
  // return the type of the stored variable
  //

  globals::atype_t type( const std::string & s ) const;
  
  avar_t * find( const std::string & name ) 
  { 
    std::map<std::string,avar_t*>::iterator aa = data.find( name );
    if ( aa == data.end() ) return NULL;
    return aa->second;
  } 
  
  bool empty() const { return data.size() == 0; } 

  bool single( const std::string * n , const avar_t * d ) const
  {
    n = NULL;
    d = NULL;
    if ( data.size() != 1 ) return false;
    std::map<std::string,avar_t*>::const_iterator aa = data.begin();
    n = &(aa->first);
    d = aa->second;
    return true;
  }

  void check( const std::string & name );
  
  // add flag 
  void set( const std::string & name );

  // add mask
  void set_mask( const std::string & name , const bool b );

  // add integer 
  void set( const std::string & name , const int i );

  // add string
  void set( const std::string & name , const std::string & s );

  // add bool
  void set( const std::string & name , const bool b );
  
  // add double
  void set( const std::string & name , const double d );

  // add integer vec
  void set( const std::string & name , const std::vector<int> &  i );

  // add string vec
  void set( const std::string & name , const std::vector<std::string> & s );

  // add bool vec
  void set( const std::string & name , const std::vector<bool> & b );
  
  // add double vec
  void set( const std::string & name , const std::vector<double> & d );

  // convenience function to add FTR metadate (i.e. str->str key/value pairs)
  void add( const std::map<std::string,std::string> & d )
  {
    std::map<std::string,std::string>::const_iterator ii = d.begin();
    while ( ii != d.end() )
      {
	set( ii->first , ii->second );
	++ii;
      }
  }

  std::string print( const std::string & delim = ";" , const std::string & prelim = "" ) const;
  
  //
  // Misc helper functions
  //
  
  // the instance controls adding data-points, so it is also responsible for clean-up

  ~instance_t();



};




//
// Generic 'value' class
//

struct avar_t
{
  
  // abstract base class for variable -> value mapping
  // avar_t belong to ainst_t, instances of a particular annotation
  
  virtual ~avar_t() { } ;  // virtual destructor for ABC
  
  virtual bool        bool_value() const = 0;
  virtual int         int_value() const = 0;
  virtual double      double_value() const = 0;
  virtual std::string text_value() const = 0;

  virtual std::vector<bool>        bool_vector() const = 0;
  virtual std::vector<int>         int_vector() const = 0;
  virtual std::vector<double>      double_vector() const = 0;
  virtual std::vector<std::string> text_vector() const = 0;

  virtual globals::atype_t     atype() const = 0;  
  
  virtual bool is_vector() const = 0 ; 

  virtual int  size() const = 0 ; 
  
  virtual avar_t *    clone() const = 0;
  
  bool has_value;
  
  bool is_missing() const { return ! has_value; } 
  
  friend std::ostream & operator<<( std::ostream & out , const avar_t & e );
  
};


struct flag_avar_t : public avar_t 
{
 public:
  flag_avar_t() { }  // empty
  ~flag_avar_t() { } 
  flag_avar_t *  clone() const { return new flag_avar_t(*this); }
  void set() { has_value=false; }

  bool bool_value() const { return true; } 
  int int_value() const { return 1; }
  double double_value() const { return 1; } 
  std::string text_value() const { return "."; } 

  std::vector<bool>        bool_vector()   const { return std::vector<bool>(0); } 
  std::vector<int>         int_vector()    const { return std::vector<int>(0); } 
  std::vector<double>      double_vector() const { return std::vector<double>(0); } 
  std::vector<std::string> text_vector()   const { return std::vector<std::string>(0); } 
  bool is_vector() const { return false; } 
  int size() const { return 0; } 

  globals::atype_t atype() const { return globals::A_FLAG_T; }  
};


//
// MASK -- just a special BOOL
//

struct mask_avar_t : public avar_t 
{
 public:
  mask_avar_t( bool b ) { set(b); }
  ~mask_avar_t() { } 
  mask_avar_t *  clone() const { return new mask_avar_t(*this); }
  bool value;
  void set( bool b ) { has_value=true; value=b; } 
  bool bool_value() const { return value; }
  int int_value() const { return value; }
  double double_value() const { return value; }
  std::string text_value() const { return has_value ? ( value ? "true" : "false" ) : "." ; } 
  globals::atype_t atype() const { return globals::A_MASK_T; }

  std::vector<bool>        bool_vector()   const { return std::vector<bool>(1,bool_value()); } 
  std::vector<int>         int_vector()    const { return std::vector<int>(1,int_value()); } 
  std::vector<double>      double_vector() const { return std::vector<double>(1,double_value()); } 
  std::vector<std::string> text_vector()   const { return std::vector<std::string>(1,text_value()); } 
  bool is_vector() const { return false; } 
  int size() const { return 1; } 
  
};



//
// Scalars
//

struct bool_avar_t : public avar_t 
{
 public:
  bool_avar_t( bool b ) { set(b); }
  ~bool_avar_t() { } 
  bool_avar_t *  clone() const { return new bool_avar_t(*this); }
  bool value;
  void set( bool b ) { has_value=true; value=b; } 
  bool bool_value() const { return value; }
  int int_value() const { return value; }
  double double_value() const { return value; }
  std::string text_value() const { return has_value ? ( value ? "true" : "false" ) : "." ; } 
  globals::atype_t atype() const { return globals::A_BOOL_T; }

  std::vector<bool>        bool_vector()   const { return std::vector<bool>(1,bool_value()); } 
  std::vector<int>         int_vector()    const { return std::vector<int>(1,int_value()); } 
  std::vector<double>      double_vector() const { return std::vector<double>(1,double_value()); } 
  std::vector<std::string> text_vector()   const { return std::vector<std::string>(1,text_value()); } 
  bool is_vector() const { return false; } 
  int size() const { return 1; } 
  
};


struct int_avar_t : public avar_t 
{
 public:
  int_avar_t( int i ) { set(i); }
  ~int_avar_t() { } 
  int_avar_t *  clone() const { return new int_avar_t(*this); }
  int value;
  void set( int i ) { has_value=true; value=i; } 
  bool bool_value() const { return value!=0; }
  int int_value() const { return value; }
  double double_value() const { return value; }
  std::string text_value() const { return has_value ? Helper::int2str( value ) : "."; } 
  globals::atype_t atype() const { return globals::A_INT_T; }

  std::vector<bool>        bool_vector()   const { return std::vector<bool>(1,bool_value()); } 
  std::vector<int>         int_vector()    const { return std::vector<int>(1,int_value()); } 
  std::vector<double>      double_vector() const { return std::vector<double>(1,double_value()); } 
  std::vector<std::string> text_vector()   const { return std::vector<std::string>(1,text_value()); } 
  bool is_vector() const { return false; } 
  int size() const { return 1; } 
  
};

struct double_avar_t : public avar_t
{
 public:
  double_avar_t( double d ) { set(d); }
  ~double_avar_t() { } 
  double_avar_t *  clone() const { return new double_avar_t(*this); }
  double value;
  void set( double d ) { has_value=true; value=d; } 
  bool bool_value() const { return value!=0; } // ignores floating-point issues...
  int int_value() const { return value; }
  double double_value() const { return value; }
  std::string text_value() const { return has_value ? Helper::dbl2str( value ) : "."; } 
  globals::atype_t atype() const { return globals::A_DBL_T; }

  std::vector<bool>        bool_vector()   const { return std::vector<bool>(1,bool_value()); } 
  std::vector<int>         int_vector()    const { return std::vector<int>(1,int_value()); } 
  std::vector<double>      double_vector() const { return std::vector<double>(1,double_value()); } 
  std::vector<std::string> text_vector()   const { return std::vector<std::string>(1,text_value()); } 
  bool is_vector() const { return false; } 
  int size() const { return 1; } 
  
};

struct text_avar_t : public avar_t 
{
 public:
  text_avar_t( const std::string & t ) { set(t); }
  ~text_avar_t() { } 
  text_avar_t *  clone() const { return new text_avar_t(*this); }
  std::string value;
  void set( const std::string & t ) { has_value=true; value=t; } 
  bool bool_value() const { return value!="0" && value != "false" ; }
  int int_value() const 
  { 
    if ( ! has_value ) return 0;
    int i=0;
    if ( ! Helper::str2int( value , &i ) ) return 0;
    return i;
  }  
  double double_value() const 
  { 
    if ( ! has_value ) return 0;
    double d = 0;
    if ( ! Helper::str2dbl( value , &d ) ) return 0;
    return d;
  }

  std::string text_value() const { return has_value ? value : "."; } 
  globals::atype_t atype() const { return globals::A_TXT_T; }  

  std::vector<bool>        bool_vector()   const { return std::vector<bool>(1,bool_value()); } 
  std::vector<int>         int_vector()    const { return std::vector<int>(1,int_value()); } 
  std::vector<double>      double_vector() const { return std::vector<double>(1,double_value()); } 
  std::vector<std::string> text_vector()   const { return std::vector<std::string>(1,text_value()); } 
  bool is_vector() const { return false; } 
  int size() const { return 1; } 

};


//
// Vectors
//

struct boolvec_avar_t : public avar_t 
{
 public:
  boolvec_avar_t( const std::vector<bool> & b ) { set(b); }
  ~boolvec_avar_t() { } 
  boolvec_avar_t *  clone() const { return new boolvec_avar_t(*this); }
  std::vector<bool> value;
  void set( const std::vector<bool> & b ) { has_value=true; value=b; } 

  bool bool_value() const { return value.size() > 0 ; }
  int int_value() const { return value.size(); }
  double double_value() const { return value.size(); }
  std::string text_value() const { return Helper::int2str( (int)value.size() ); }
  
  globals::atype_t atype() const { return globals::A_BOOLVEC_T; }
  
  std::vector<bool>        bool_vector()   const { return value; } 
  std::vector<int>         int_vector()    const { return annot_t::as_int_vec( value ); } 
  std::vector<double>      double_vector() const { return annot_t::as_dbl_vec( value ); }
  std::vector<std::string> text_vector()   const { return annot_t::as_txt_vec( value ); }
  bool is_vector() const { return true; } 
  int size() const { return value.size(); } 
  
};


struct intvec_avar_t : public avar_t 
{
 public:
  intvec_avar_t( const std::vector<int> & i ) { set(i); }
  ~intvec_avar_t() { } 
  intvec_avar_t *  clone() const { return new intvec_avar_t(*this); }
  std::vector<int> value;
  void set( const std::vector<int> & i ) { has_value=true; value=i; } 

  bool bool_value() const { return value.size() > 0 ; }
  int int_value() const { return value.size(); }
  double double_value() const { return value.size(); }
  std::string text_value() const { return Helper::int2str( (int)value.size() ); }
  
  globals::atype_t atype() const { return globals::A_INTVEC_T; }
  
  std::vector<bool>        bool_vector()   const { return annot_t::as_bool_vec( value ); } 
  std::vector<int>         int_vector()    const { return value ; } 
  std::vector<double>      double_vector() const { return annot_t::as_dbl_vec( value ); }
  std::vector<std::string> text_vector()   const { return annot_t::as_txt_vec( value ); }
  bool is_vector() const { return true; } 
  int size() const { return value.size(); } 

};

struct doublevec_avar_t : public avar_t
{
 public:
  doublevec_avar_t( const std::vector<double> & d ) { set(d); }
  ~doublevec_avar_t() { } 
  doublevec_avar_t *  clone() const { return new doublevec_avar_t(*this); }
  std::vector<double> value;
  void set( const std::vector<double> & d ) { has_value=true; value=d; } 

  bool bool_value() const { return value.size() > 0 ; }
  int int_value() const { return value.size(); }
  double double_value() const { return value.size(); }
  std::string text_value() const { return Helper::int2str( (int)value.size() ); }
  
  globals::atype_t atype() const { return globals::A_DBLVEC_T; }
  
  std::vector<bool>        bool_vector()   const { return annot_t::as_bool_vec( value ); } 
  std::vector<int>         int_vector()    const { return annot_t::as_int_vec( value ); } 
  std::vector<double>      double_vector() const { return value ; }
  std::vector<std::string> text_vector()   const { return annot_t::as_txt_vec( value ); }
  bool is_vector() const { return true; } 
  int size() const { return value.size(); } 
  
};

struct textvec_avar_t : public avar_t 
{
 public:
  textvec_avar_t( const std::vector<std::string> & t ) { set(t); }
  ~textvec_avar_t() { } 
  textvec_avar_t *  clone() const { return new textvec_avar_t(*this); }
  std::vector<std::string> value;
  void set( const std::vector<std::string> & t ) { has_value=true; value=t; } 

  bool bool_value() const { return value.size() > 0 ; }
  int int_value() const { return value.size(); }
  double double_value() const { return value.size(); }
  std::string text_value() const { return Helper::int2str( (int)value.size() ); }
  
  globals::atype_t atype() const { return globals::A_TXTVEC_T; }
  
  std::vector<bool>        bool_vector()   const { return annot_t::as_bool_vec( value ); } 
  std::vector<int>         int_vector()    const { return annot_t::as_int_vec( value ); } 
  std::vector<double>      double_vector() const { return annot_t::as_dbl_vec( value ); }
  std::vector<std::string> text_vector()   const { return value; }

  bool is_vector() const { return true; } 
  int size() const { return value.size(); } 

};








struct annotation_set_t;
struct edf_t;
struct param_t;

void summarize_annotations( edf_t & edf , param_t & param );





struct annotation_set_t
{
  
  ~annotation_set_t()
  {
    clear();
  }
  
  std::map<std::string,annot_t*> annots;
  
  annot_t * add( const std::string & name ) 
  { 
    
    if ( annots.find( name ) != annots.end() )
      return annots[name];
    
    annot_t * a = new annot_t( name ); 
    annots[ name ] = a;
    return a;
  }
  
  annot_t * operator()( const std::string & name ) 
  {
    return find( name );
  }

  annot_t * find( const std::string & name ) 
  {     
    if ( annots.find( name ) != annots.end() )
      return annots[name];
    else 
      return NULL;
  }

  void clear() 
  { 
    std::map<std::string,annot_t*>::iterator ii = annots.begin();
    while ( ii != annots.end() ) 
      {
	delete ii->second;
	++ii;
      }
    annots.clear(); 
  }
  
  void clear( const std::string & name )
  {
    std::map<std::string,annot_t*>::iterator ii = annots.find( name );
    if ( ii != annots.end() )
      {
	delete ii->second;
	annots.erase( ii ); 
      }
  }
  
  std::vector<std::string> names() const 
  {
    std::vector<std::string> n;
    std::map<std::string,annot_t*>::const_iterator ii = annots.begin();
    while ( ii != annots.end() )
      {
	n.push_back( ii->first );
	++ii;
      }
    return n;
  }
  

  // Attempt to create a single SLEEP STAGE annotation from multiple
  // other annotations
  
  bool make_sleep_stage( const std::string & a_wake = "" , 
			 const std::string & a_n1 = "", 
			 const std::string & a_n2 = "", 
			 const std::string & a_n3 = "", 
			 const std::string & a_n4 = "", 
			 const std::string & a_rem = "",
			 const std::string & a_other = "");

};


#endif
