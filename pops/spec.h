
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


#ifdef HAS_LGBM

#ifndef __LUNA_POPS_SPEC_H__
#define __LUNA_POPS_SPEC_H__

#include <map>
#include <set>
#include <cstring>
#include <vector>

enum pops_feature_t
  {
    POPS_LOGPSD ,   // LWR , UPR (fixed 0.25 Hz intervals) 
    POPS_RELPSD ,   // LWR , UPR , NORM-LWR , NORM-UPR ( scale PSD by sum of NORM )  
    POPS_CVPSD ,    // LWR , UPR 
    POPS_BANDS ,    // fixed power bands
    POPS_RBANDS ,   // relative bands
    POPS_VBANDS ,   // variance 
    POPS_COH ,      // two-channel coherence
    POPS_SLOPE ,    // Fixed 30-45 Hz ; fixed other param
    POPS_SKEW , 
    POPS_KURTOSIS ,
    POPS_HJORTH ,
    POPS_FD ,
    POPS_PE ,
    POPS_MEAN ,
    POPS_COVAR ,     // from ivars
    
    // level 2 features
    POPS_TIME ,          // level 2 features
    POPS_SMOOTH , 
    POPS_DENOISE ,  
    POPS_SVD ,     
    POPS_NORM ,    
    POPS_RESCALE, // -1, +1 sigmoid encoding
    POPS_CUMUL, 
    POPS_DERIV, 
    
    // not a feature, but a rule to remove epochs
    POPS_EPOCH_OUTLIER 
    
  };



struct pops_spec_t {  
  
  std::string block;
  
  pops_feature_t ftr;
  
  std::string ch; 

  int size;
  
  std::map<std::string,std::string> arg;

  bool has( const std::string & key ) const
  {
    return arg.find( key ) != arg.end() ;
  }

  double narg( const std::string & key ) const
  {
    if ( ! has( key ) ) return 0;
    double d;
    if ( ! Helper::str2dbl( arg.find( key )->second , &d ) )
      Helper::halt( "problem converting string -> numeric: " + key );
    return d;
  }
  
  int cols(int *) ;
  
};


struct pops_channel_t {

  pops_channel_t( const std::string & ch , 
		  const std::set<std::string> & aliases , 
		  const int sr , 
		  const std::string unit )
  : ch(ch), aliases(aliases), sr(sr) , unit(unit)
  { } 
  
  pops_channel_t() {
    ch = "";
    sr = 0;
    aliases.clear();
    unit = "uV";
  } 

  // main label
  std::string ch;
  // other aliases
  std::set<std::string> aliases;
  int sr;
  std::string unit;

  bool match( const std::string & s , std::string * label ) const
  {
    if ( s == ch ) 
      {
	*label = ch;
	return true;
      }
    
    std::set<std::string>::const_iterator ii = aliases.find( s );
    if ( ii == aliases.end() )
      return false;
    *label = *ii;
    return true;
    
  }


  
};


struct pops_specs_t {
  
  // read a specification file
  void read( const std::string & f );
  
  // get mappings
  void init();
  static std::map<std::string,pops_feature_t> lab2ftr;
  static std::map<pops_feature_t,std::string> ftr2lab;
  static std::set<std::string> lvl2;

  // defaults
  static std::vector<std::string> defaults;
  void init_default();
  
  // specs attached?
  bool loaded() const { return specs.size() != 0 ; } 

  // give columns for a spec/channel level1 feature
  bool has( pops_feature_t , const std::string & ch );
  std::vector<int> cols( pops_feature_t , const std::string & ch );
  std::vector<int> block_cols( const std::string & , int );  
  std::map<pops_feature_t,std::map<std::string,std::vector<int> > > ftr2ch2col; // track which l1-features/cols selected
  void build_colmap(); // function to construct ftr2ch2col, called after read()
  void check_args();
  
  // track level 1 feature/channel pairings
  std::map<pops_feature_t,std::map<std::string,pops_spec_t> > fcmap;
  
  // channel info
  std::map<std::string,pops_channel_t> chs;
  
  // all features (in order)
  std::vector<pops_spec_t> specs;
  
  // blocks --> size
  static std::map<std::string,int> blocksize;
  
  // SELECT <blocks>
  std::set<std::string> selected;

  // DROP <vars>
  std::set<std::string> dropped;
  
  // expanded feature lists
  std::vector<std::string> col_block;
  std::vector<std::string> col_label;
  std::vector<std::string> col_original_label;
  std::vector<std::string> col_root;
  std::vector<bool>        col_select;
  std::vector<int>         col_level;  
  std::map<int,int>        orig2final;
  std::map<int,int>        final2orig;
  
  
  // number of features
  int n1; // level 1 only (all)
  int na; // all features (l1+l2)
  int nf; // final (selected) number of features
  int ns; // number of signals in total (whether used or not) 
  
  int total_cols() const;
  int select_cols() const;

  // get feature labels
  //  std::vector<std::string> total_labels();
  std::vector<std::string> select_labels();

  // prior to any replace=X,Y -- so that any .ranges files
  // will still map to this channel
  std::vector<std::string> select_original_labels();

  // var roots only (for inc-vars and exc-vars) [ OLD ] 
  std::vector<std::string> select_roots();

  // blocks only (for inc-vars and exc-vars)
  std::vector<std::string> select_blocks();
  
  
};

#endif
#endif
