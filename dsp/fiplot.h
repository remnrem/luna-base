#ifndef __FIPLOT_H__
#define __FIPLOT_H__

#include <vector>
#include <map>

#include "../edf/edf.h"
#include <string>
#include "../db/db.h"


void fiplot_wrapper( edf_t & edf , const param_t & param );


struct fipair_t 
{ 
  fipair_t(double w, double n) : w(w), n(n) { } 
  fipair_t()  { w=0; n=0; } 
  double w; double n; 
};

struct fibin_t { std::map<double,fipair_t> r; };

struct fipoint_t 
{ 
  fipoint_t( int _i, int _j , double _h )
  {
    i = _i; j = _j; h = _h; t = j - i + 1 ;
  }

  int i; 
  int j;
  int t;
  double h; 

  bool operator< ( const fipoint_t & rhs ) const 
  {
    // sort by duration, longest first
    if ( t > rhs.t ) return true;
    if ( t < rhs.t ) return false;
    return i < rhs.i;
  }
}; 

struct fiplot_t
{
  fiplot_t( const std::vector<double> & x , const std::vector<uint64_t> * tp , const int _fs  ,
	    double t_lwr , double t_upr , double t_inc , bool cycles , 
	    double f_lwr , double f_upr , double f_inc , bool logspace = true )
  {
    fs = _fs;
    set_t( t_lwr, t_upr, t_inc , cycles );
    set_f( f_lwr, f_upr, f_inc , logspace );
    proc(x,tp,fs);
  }
  
  void set_t( double lwr , double upr , double inc , bool cyc ) 
  {
    t_lwr = lwr;
    t_upr = upr;
    t_inc = inc;
    cycles = cyc;
  }

  void set_f( double lwr , double upr , double inc , bool logspace ) 
  {
    
    frqs.clear();

    f_lwr = lwr;
    f_upr = upr;
    f_inc = inc; // # of inc if logspace == T
    
    if ( ! logspace )
      {
	for (double f = f_lwr ; f <= f_upr ; f += f_inc ) frqs.push_back( f );
      }
    else
      {
	frqs = MiscMath::logspace( f_lwr , f_upr , f_inc );
      }
  }

  void proc( const std::vector<double> & x , 
	     const std::vector<uint64_t> * tp , 
	     const int fs ); 
    
  int nt; // number of time points
  int nf; // number of frequencies
  
  int fs; // sample rate

  std::vector<double> frqs;
  double f_lwr, f_upr, f_inc;
  
  double t_lwr, t_upr, t_inc;
  bool cycles;

  std::vector<double> cwt( const std::vector<double> & x , const int fs, const double fc , const int num_cycles );
  
  fibin_t intervalize( const std::vector<double> & x , 
		       const std::vector<uint64_t> * tp , 
		       const int fs , 
		       const double t_lwr , 
		       const double t_upr , 
		       const double t_inc , 
		       const bool cycles , 
		       const double fc );
  
  
  
};


#endif
