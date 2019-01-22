
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


#ifndef __PDC_H__
#define __PDC_H__

#include <vector>
#include <map>
#include <set>

#include "helper/helper.h"
#include "stats/matrix.h"

struct edf_t;

struct param_t;

struct pdc_obs_t { 
	       
  pdc_obs_t( const int q ) { init(q); }

  // observation ID
  std::string id;
  
  // has channel?
  std::vector<bool> ch;

  // (optional) time-series data (per channel)
  std::vector<std::vector<double> > ts;
  
  // permutation-distribution data (per channel) 
  std::vector<std::vector<double> > pd;
  
  // aux. info: primary label
  std::string label;
  
  // aux. info: all other info, as key=value pairs
  std::map<std::string,std::string> aux;

  //
  // Functions
  //

  // encode TS(s) as PD(s)
  void encode( int , int );
  
  // add (combine PD counts; *assumes* channels, m/t etc are the same...
  void add( const pdc_obs_t & rhs )
  {
    if ( pd.size() != rhs.pd.size() ) Helper::halt( "cannot add pdc_obs_t" );
    for (int i=0;i<pd.size();i++)
      {
	if ( pd[i].size() == 0 ) pd[i] = rhs.pd[i];
	else
	  {
	    //std::cout << "Sz = " << pd[i].size() << " " << rhs.pd[i].size() << "\n";
	    if ( pd[i].size() != rhs.pd[i].size() ) Helper::halt( "internal pdc_obs_t prob" );
	    for (int j=0;j<pd[i].size();j++) pd[i][j] += rhs.pd[i][j];
	  }
      }
  }
  
  void norm( double s )
  {
    for (int i=0;i<pd.size();i++)      
      for (int j=0;j<pd[i].size();j++) pd[i][j] /= s;
    
  }


  // get entropy of current PD(s)
  std::vector<double> entropy() const;
  
  void init(int q)
  {
    if ( q == 0 ) Helper::halt( "must set channel space before adding observations");
    id = "";
    label = "";
    aux.clear();
    
    ch.resize( q , false );
    ts.clear(); pd.clear();
    ts.resize( q ); pd.resize( q );

  }
  
  
};


//
// Helper struct for ordering
//

struct pd_dist_t { 
  pd_dist_t( double d , int ix ) : d(d) , ix(ix) { } 
  double d;
  int ix;
  bool operator<(const pd_dist_t & rhs ) const {
    if ( d < rhs.d ) return true;
    if ( d > rhs.d ) return false;
    return ix < rhs.ix;
  }
};



struct pdc_t { 

  friend struct pdc_obs_t;
  
  pdc_t( const bool b = true ) 
  { 
    clear();
    store_timeseries(b);
  }

  static void add( const pdc_obs_t & ob );
  
  void store_timeseries( const bool b ) { store_ts = b; } 
  
  bool store_timeseries() const { return store_ts; } 

  //
  // add an observation
  //

  

  //
  // summarize all channels availability for all observations
  //

  static void channel_check();
  
  
  //
  // set primary parameters
  //

  static void set_param( int _m , int _t )
  {
    m = _m;
    t = _t;
  }


  //
  // encode all observations given current m and t
  //
  
  static void encode_ts();
  
  
  //
  // determine optimal m and t, given entropy heuristic
  // potentially stratified by label
  //

  static void entropy_heuristic_wrapper( param_t & );

  static void entropy_heuristic( int min_m = 2 , int max_m = 7 , int min_t = 1 , int max_t = 5 , bool by_cat = false );

  
  //
  // tidy up
  //

  static void clear()
  {
    obs.clear();
    labels.clear();
    label_count.clear();
    q=0;
    channels.clear();
  }

  
  //
  // Helper test function
  //

  void test();
  
  
  //
  // Read/write current observations 
  //
  
  void write_pdlib( const std::string & );
  
  static void read_pdlib( const std::string & , const std::set<std::string> * incl_chs = NULL );
  
  static void read_tslib( const std::string & );

  //
  // Wrapper to get TS-LIB from EDFs
  //

  static void construct_tslib( edf_t & edf , param_t & param );

  
  //
  // Wrapper to generate an epoch-by-epoch ("ExE") matrix
  //

  static void similarity_matrix( edf_t & edf , param_t & param );
  

  //
  // Helper function to make a PD-LIB from a TS-LIB
  //

  static void construct_pdlib( param_t & param );

  //
  // Add channels
  //

  static int add_channel( const std::string & c )
  {
    std::map<std::string,int>::const_iterator cc = channels.find(c);
    if ( cc != channels.end() ) return cc->second;
    q = channels.size() + 1;
    channels[ c ] = q-1;
    return q-1;
  }
  
  static bool has_channel( const std::string & c ) 
  { 
    return channels.find(c) != channels.end(); 
  } 
  
  static int channel( const std::string & c ) 
  {
    std::map<std::string,int>::const_iterator cc = channels.find( c );
    if ( cc == channels.end() ) return -1;
    return cc->second;
  }

  //
  // Primary routine to autoscore, assuming a library has been attached
  //

  static void simple_sleep_scorer( edf_t & edf , param_t & param );

  //
  // Write results of the above to XML (NSRR format)
  //

  static void write_xml( const std::string & filename , const std::vector<std::string> & stges );
  
  //
  // Get distance between two observations, from all channels 
  //
  
  static double distance( const pdc_obs_t & a, const pdc_obs_t & b );  
  
  //
  // Get distance between two observations, from a subset of channels
  //

  static double distance( const pdc_obs_t & a, const pdc_obs_t & b , const std::vector<int> & chs );  
  
  
  //
  // Generate an all-by-all distance matrix
  //

  static Data::Matrix<double> all_by_all();


  //
  // For a single observation, find the best nmatches in terms of 'label'
  //
  
  static std::set<pd_dist_t> match( const pdc_obs_t & target , const int nbest = 10 );

  static std::map<std::string,double> summarize( const std::set<pd_dist_t> & matches , std::string * cat , double * conf );
    
  
private: 


  
  //
  // Members
  //

  
  static int m;  // embedding dimension (default = 5 range 2 .. 7 )
  
  static int t;  // time delay (default = 1, range 1 .. 10 )

  static bool store_ts;


  // for each observation, stored PD, and optionally, the original TS 
  // the main label, and any aux data

  static std::vector<pdc_obs_t> obs;

  //
  // track channels (these will be similar for all pdc_obs_t)
  //

  // total number of channels
  static int q; 

  // all channels in PDLIB
  static std::map<std::string,int> channels;


  //
  // track labels
  //

  static std::set<std::string> labels;
  static std::map<std::string,int> label_count;
  
  //
  // Helper functions
  //
  

  // primary function to get the PD from a particular time-series
  // if sum set to 0, store integers (unnormalized) 
  // if sum !=0 , then automatically normalize (i.e. return vector with sum of 1.0)
  static std::vector<double> calc_pd( const std::vector<double> & x , int m , int t , int * sum );

  static int num_pd(int m);
  
  static double entropy( const std::vector<double> & );
  
  static double hellinger( const std::vector<double> & , const std::vector<double> & );

  static double squared_hellinger( const std::vector<double> & , const std::vector<double> & );
  
  static double symmetricAlphaDivergence( const std::vector<double> & , const std::vector<double> & );
  
  static int codebook( const std::vector<double> & , int, int, int );
  

};


#endif
