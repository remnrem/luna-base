
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

struct edf_t;

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

enum etype_t { EVT_BOOL , EVT_INT , EVT_DBL , EVT_TXT };

class event_t
{
 public:
  event_t() { }
  event_t( const std::string & label ) : label(label) { }  

  virtual ~event_t() { } ;  // virtual destructor for ABC
  
  std::string label;
  bool        has_value;
  std::map<std::string,std::string> metadata;
  
  virtual bool       bool_value() const = 0;
  virtual int        int_value() const = 0;
  virtual double     double_value() const = 0;
  virtual std::string text_value() const = 0;
  virtual etype_t    event_type() const = 0;
  
  virtual event_t *  clone() const = 0;
  
  bool operator<( const event_t & rhs ) const
  {
    // simple sorting only on label here
    return label < rhs.label ;
  }    

  friend std::ostream & operator<<( std::ostream & out , const event_t & e );

};

class bool_event_t : public event_t 
{
 public:
  bool_event_t( const std::string & l , bool b ) { set(l,b); }
  ~bool_event_t() { } 
  bool_event_t *  clone() const { return new bool_event_t(*this); }
  bool value;
  void set( const std::string & l , bool b ) { has_value=true; label=l; value=b; } 
  bool bool_value() const { return value; }
  int int_value() const { return value; }
  double double_value() const { return value; }
  std::string text_value() const { return has_value ? ( value ? "1" : "0" ) : "." ; } 
  etype_t event_type() const { return EVT_BOOL; }
  
};

class int_event_t : public event_t 
{
 public:
  int_event_t( const std::string & l , int i ) { set(l,i); }
  ~int_event_t() { } 
  int_event_t *  clone() const { return new int_event_t(*this); }
  int value;
  void set( const std::string & l , int i ) { has_value=true; label=l; value=i; } 
  bool bool_value() const { return value!=0; }
  int int_value() const { return value; }
  double double_value() const { return value; }
  std::string text_value() const { return has_value ? Helper::int2str( value ) : "."; } 
  etype_t event_type() const { return EVT_INT; }
  
};

class double_event_t : public event_t
{
 public:
  double_event_t( const std::string & l , double d ) { set(l,d); }
  ~double_event_t() { } 
  double_event_t *  clone() const { return new double_event_t(*this); }
  double value;
  void set( const std::string & l , double d ) { has_value=true; label=l; value=d; } 
  bool bool_value() const { return value!=0; } // ignore floating-point issue
  int int_value() const { return value; }
  double double_value() const { return value; }
  std::string text_value() const { return has_value ? Helper::dbl2str( value ) : "."; } 
  etype_t event_type() const { return EVT_DBL; }
  
};

class text_event_t : public event_t 
{
 public:
  text_event_t( const std::string & l , const std::string & t ) { set(l,t); }
  ~text_event_t() { } 
  text_event_t *  clone() const { return new text_event_t(*this); }
  std::string value;
  void set( const std::string & l , const std::string & t ) 
  { has_value=true; label=l; value=t; } 
  bool bool_value() const { return value!="0" && value != "F" ; }

  int int_value() const 
  { 
    int i=0;
    if ( Helper::str2int( value , &i ) ) return 0;
    return i;
  }
  
  double double_value() const 
  { 
    double d = 0;
    if ( Helper::str2dbl( value , &d ) ) return 0;
    return d;
  }

  std::string text_value() const { return has_value ? value : "."; } 
  etype_t event_type() const { return EVT_TXT; }
  
};


struct annotation_set_t;

enum ANNOTATION
  {
    ATYPE_UNKNOWN ,     // not allowed
    ATYPE_BINARY ,
    ATYPE_INTEGER , 
    ATYPE_QUANTITATIVE ,
    ATYPE_TEXTUAL ,     
    ATYPE_COLCAT ,       // colour categories
    ATYPE_SLEEP_STAGE ,  // special case of INTEGER
    ATYPE_MASK
  };


struct edf_t;
struct param_t;

void summarize_annotations( edf_t & edf , param_t & param );


typedef std::vector<const event_t*> evt_table_t;
typedef std::map<interval_t,std::vector<const event_t*> > interval_evt_map_t;


class annot_t
{

  friend struct timeline_t;
  friend struct annotation_set_t;

 public:
  
  std::map<ANNOTATION,std::string> type_name;
 
 private: 

 annot_t() { set_type_names(); } 

 annot_t( const std::string & n ) 
   : name(n) 
  { 
    file="";
    cols.clear();
    type.clear();
    min = 0;
    max = 0;
    has_range = false;
    set_type_names(); 
 }

 public:

 ~annot_t() 
   {
     std::set<event_t *>::iterator ii = all_events.begin();
     while ( ii != all_events.end() )
       {	 
	 delete *ii;
	 ++ii;
       }
   }
 
  bool load( const std::string & );  
  int  load_features( const std::string & );    
  bool save( const std::string & );  
  
  static bool loadxml( const std::string & , edf_t * );
  bool savexml( const std::string & );  
  static void dumpxml( const std::string & , bool );
  
  void set_description( const std::string & d ) { description = d; } 
  void add_col( const std::string & c ) { cols.push_back(c); }
  void add_col( const std::string & c , ANNOTATION t ) 
  { cols.push_back(c); type.push_back(t);}

  void add( interval_t interval , const event_t * e ) 
  {
    event_t * ne = e->clone();
    interval_events[ interval ].push_back( ne );
    all_events.insert( (event_t*)ne );
  }

  int num_interval_events() const { return interval_events.size(); } 
  
  void set_type_names()
  { 
    type_name[ ATYPE_UNKNOWN ]      = ".";
    type_name[ ATYPE_BINARY ]       = "BINARY";
    type_name[ ATYPE_INTEGER ]      = "INTEGER";
    type_name[ ATYPE_SLEEP_STAGE ]  = "SLEEP_STAGE";
    type_name[ ATYPE_TEXTUAL ]      = "TEXTUAL";
    type_name[ ATYPE_MASK ]         = "MASK";
    type_name[ ATYPE_QUANTITATIVE ] = "QUANTITATIVE";
  }

  std::string type_string(const int f ) const { return type_name.find( type[f] )->second; } 

  uint64_t minimum_tp() const;  
  uint64_t maximum_tp() const;

  interval_evt_map_t extract( const interval_t & window );

  std::string                              name;
  std::string                              file;
  std::string                              description;

  std::vector<std::string>                 cols;
  std::vector<ANNOTATION>                  type;

  //
  // optionally, a RANGE specified for this annotation
  // (can be used in plotting, etc)
  //

  bool has_range;
  double min;
  double max;
  
  bool univariate() const { return cols.size() == 0; }
  bool multivariate() const { return cols.size() > 1; }

  
  

 private:

  interval_evt_map_t  interval_events;

  // just used to clean up 
  std::set<event_t*> all_events;

  void reset()
  {        
    name = "";
    file = "";
    description = "";
    cols.clear();
    has_range = false;
    min = max = 0;
    interval_events.clear();
  }

};


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
  
  std::vector<std::string> display_names() const 
  {
    // return actual names in the file, i.e. not the handle (map key) 
    // which is often the file name, if the .annot was loaded in...
    std::vector<std::string> n;
    std::map<std::string,annot_t*>::const_iterator ii = annots.begin();
    while ( ii != annots.end() )
      {
	n.push_back( ii->second->name );
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
