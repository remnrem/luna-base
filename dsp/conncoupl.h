
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

#ifndef __CONNCOUPL_H__
#define __CONNCOUPL_H__

#include <vector>
#include "defs/defs.h"
#include "stats/matrix.h"

struct edf_t;

struct param_t;



struct conncoupl_res_t {

  conncoupl_res_t() { } 

  conncoupl_res_t( const int ne , const int nt )
  {
    stats.resize(ne,nt);
    emp_z.resize(ne,nt);
  }
  
  Data::Matrix<double> stats;
  Data::Matrix<double> emp_z;
  
};


struct conncoupl_t
{

  //
  // filter-Hilbert 
  //
  
  conncoupl_t( edf_t & edf ,
	       signal_list_t & signals ,
	       const int sr , 		  
	       const std::vector<freq_range_t> & fint1 , 
	       const std::vector<freq_range_t> & fint2 ,    		  
	       const double ripple = 0.01 ,  
	       const double tw = 0.5 ,
	       const int nreps = 1000 ,  
	       const int es = 0 ,
	       const bool xpac = false )
    : edf(edf) , signals(signals) , sr(sr) , fint1(fint1) , fint2(fint2) , nreps(nreps), tw(tw) , ripple(ripple) , es(es) , xpac(xpac)  
  {
    
    use_hilbert = true;
    
    setup();
    
    pre_calc();

    calc();

  }

  
  //
  // wavelets
  //

  conncoupl_t( edf_t & edf ,
	       signal_list_t & signals ,
	       const int sr , 		  

	       const std::vector<double> & fc1 ,   // center frequency (min/max)
	       const std::vector<double> & fwhm1 , // 
	       const int num1 ,  // number of inc

	       const std::vector<double> & fc2 ,   // center frequency (min/max)
	       const std::vector<double> & fwhm2 , // 
	       const int num2 ,  // number of inc

	       const double tlen , // length of wavelets (in seconds)
	       
	       const int nreps = 1000 ,
	       const int es = 0 ,
	       const bool xpac = false ,
	       const bool dump_wavelets = false ) :
  edf(edf) , signals(signals) , sr(sr) ,
    fc1(fc1) , fwhm1( fwhm1 ) , num1( num1 ),
    fc2(fc2) , fwhm2( fwhm2 ) , num2( num2 ),
    tlen( tlen ), 
    nreps(nreps) , es(es) , xpac( xpac ) ,
    dump_wavelets( dump_wavelets )
  {

    use_hilbert = false;
    
    setup();
    
    pre_calc();
    
    calc();

  }
  
	

  
  //
  // core functions
  //

  void setup();

  void pre_calc();

  void calc();


  
private:
  
  //  
  // members
  //
  
  edf_t & edf;

  signal_list_t & signals;

  int sr;


  //
  // Aggregate results
  //
  
  std::map<std::string,conncoupl_res_t> results;

  //
  // mode (filter-Hilbert vs wavelets) 
  //
  
  bool use_hilbert;

  

  // wavelets
  
  std::vector<double> fc1, fwhm1;

  std::vector<double> fc2, fwhm2; 

  int num1 , num2 ;
  
  double tlen; 

  bool dump_wavelets;
  
  // hilbert

  std::vector<freq_range_t> fint1; 
  
  std::vector<freq_range_t> fint2; 
  
  double tw;

  double ripple;
    
  // permutations
  
  int nreps;

  int es;  

  int es_pts;
  
  std::vector<int> offset;

  //
  // Contrasts
  //

  // within-channel, different freq. (PAC)    [ Y if second set of freqs given ] 
  // across-channek, within freq. (wPLI, coherence, PLV)  [ Y if multiple channels given ] 
  // cross-channel PAC (cross-frequency) [ Y if xpac option is set ] 
  
  bool xpac;
  
  
  //
  // Core transformed complex signal(s) 
  // populated by pre_calc()
  //  
  // epoch x channel x frequency x sample-points (whole trace)
  //
  
  std::vector<std::vector<std::map<std::string, std::vector<dcomp> > > > a;
  std::vector<std::vector<std::map<std::string, std::vector<dcomp> > > > a_conj;

  // contrasts

  std::map<std::string,freq_range_t> fmap;  
  std::vector<int> s1, s2;
  std::vector<std::string> f1, f2;
  std::vector<double> disp_f1, disp_f2;
  std::vector<bool> cfc, xch;
  

  //
  // helpers
  //

  std::string str( const freq_range_t & f ) const {
    return Helper::dbl2str( f.first ) + ".." + Helper::dbl2str( f.second ) ;
  }

  // misc/old
  
  /* bool bin( double d , int * b , const std::vector<double> & th , const int nbins ); */
  /* double test_uniform( const std::vector<std::vector<double> > & m ); */

  
};



namespace dsptools 
{  
  // wrapper
  void connectivity_coupling( edf_t & , param_t & );
}

#endif
