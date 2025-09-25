
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

#ifndef __GPA_H__
#define __GPA_H__

struct param_t;

#include "stats/Eigen/Dense"
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include "helper/helper.h"


struct bfile_t {

  bfile_t( const std::string & name ) 
    : name( name ) { } 
  
  bool write( const std::vector<std::string> & ids ,
	      const std::vector<std::string> & vars ,
	      const std::map<std::string,std::string> & vargroup,
	      const std::map<std::string,std::string> & basevar,	      
	      const std::map<std::string,std::map<std::string,std::string> > & faclvl , 
	      const Eigen::MatrixXd & X ) ; 
  
  bool read( const std::set<std::string> & incvars ,
	     const std::set<std::string> & excvars ,
	     const std::set<std::string> & inclvars,
	     const std::set<std::string> & exclvars,
	     const std::vector<std::pair<int,int> > & incnums,
	     const std::vector<std::pair<int,int> > & excnums, 
	     const std::set<std::string> & incfacs,
	     const std::set<std::string> & excfacs,
	     const std::map<std::string,std::set<std::string> > & incfaclvls,
	     const std::map<std::string,std::set<std::string> > & excfaclvls,
	     const std::set<std::string> & incgrps,
             const std::set<std::string> & excgrps,	     
	     const std::set<std::string> & xz_incvars,
             const std::set<std::string> & xz_incgrps,	     
	     const std::set<std::string> & y_incvars,
             const std::set<std::string> & y_incgrps,	     
	     std::vector<std::string> * ids ,
	     std::vector<std::string> * vars ,
	     std::map<std::string,std::string> * grps ,
	     std::map<std::string,std::string> * basevar,	     
	     std::map<std::string,std::map<std::string,std::string> > * faclvl ,
	     Eigen::MatrixXd * X );

  int rows() const { return ni; }

  int cols() const { return nv; }
  
  
private:

  std::string name;
  int ni, nv;

  inline static void bwrite( std::ofstream & O , const std::string & s ) 
  {
    uint8_t l = s.size();    
    O.write( (char*)( &l ), sizeof(uint8_t) );
    if ( l ) O.write( s.c_str(), l );      
  }

  inline static void bwrite( std::ofstream & O , int i ) 
  {
    O.write( (char*)( &i ), sizeof(int) );
  }
  
  inline static void bwrite( std::ofstream & O , double d ) 
  {
    O.write( (char*)( &d ), sizeof(double) );
  }
  
  inline static std::string bread_str( std::ifstream & I )
  {
    uint8_t len;
    I.read( (char*)( &len ), sizeof(uint8_t) );
    if ( len == 0 ) return "";
    std::vector<char> b( len );
    I.read( &b[0] , len );
    std::string s( b.begin() , b.end() );
    return s;
  }
  
  inline static int bread_int( std::ifstream & I )
  {
    int i;
    I.read( (char*)( &i ), sizeof(int) );
    return i;
  }
  
  inline static double bread_dbl( std::ifstream & I )
  {
    double d;
    I.read( (char*)( &d ), sizeof(double) );
    return d;
  }
  
  inline static void bskip_dbl( std::ifstream & I , const int n )
  {
    std::vector<double> dummy( n ) ;
    I.read( (char*)( &dummy[0] ), n * sizeof(double) );
  }
  
  inline static void bskip_int( std::ifstream & I , const int n )
  {
    std::vector<double> dummy( n ) ;
    I.read( (char*)( &dummy[0] ), n * sizeof(int) );
  }
  
  
};



struct gpa_t { 

  gpa_t( param_t & param , const bool prep_mode );

  // prep (write bfile)
  void prep();

  // read inputs spec from file
  void parse( const std::string & pfile );
	     
  // read bfile
  void read();

  // subset rows
  void subset( const std::set<int> & , const std::map<int,bool> & );

  // drop NA cols
  void drop_null_columns();
  
  // QC matrix
  void qc( const double winsor , const bool stats_mode );

  // kNN imputation of missing points
  void knn_imputation( param_t & );
  
  // dump binary as text (to stdout)
  void dump();

  // info on faclvl
  void manifest();  

  // summary of manifest
  void summarize();

  // dump means/SDs/Ns, optionally by strata
  void stats( const std::vector<std::string> * s = NULL );

  // make strata (for stats())

  std::vector<std::map<std::string,std::string> > make_strata( const std::vector<std::string> & ,
							       std::vector<std::vector<int>> * x );

  bool next_combo(std::vector<int>& a, const std::vector<int>& L);
  
    
  // run
  void run();   // correct for all X considered
  void run1X(); // correction w/in X 

    
private:

  // ids
  std::vector<std::string> ids;

  // var labels
  std::vector<std::string> vars;
  
  // track factors: var -> factor -> level[str] 
  std::map<std::string,std::map<std::string,std::string> > faclvl;
  
  // track base var
  std::map<std::string,std::string> basevar;

  // track group
  std::map<std::string,std::string> var2group;
  std::map<std::string,std::string> file2group;
  
  // data
  Eigen::MatrixXd X;
  
  // include d-vars, i-vars and co-vars
  std::vector<int> dvs; // Y
  std::vector<int> ivs; // X
  std::vector<int> cvs; // Z

      
  //
  // options  
  //

  // [prep+run] binary matrix file 
  std::string bfile;  
  
  // [prep] input files (and stratifying factors)
  std::map<std::string,std::set<std::string> > infiles;

  // file vars (if defined, restrict to only these basevars from a file) 
  std::map<std::string,std::string> file2basevars;

  // aliases
  std::map<std::string,std::map<std::string,std::string> > file2var2alias;

  // re-mappings ( file -> var -> level -> num-value )
  std::map<std::string,std::map<std::string,std::map<std::string,double> > > file2var2mapping;

  
  // [prep/read] variables in included/exclude [on base vars] 
  std::set<std::string> incvars, excvars;
  std::set<std::string> incfacs, excfacs;
  std::set<std::string> incgrps, excgrps;

  // XYZ specific includes
  std::set<std::string> xz_incvars, xz_incgrps;
  std::set<std::string> y_incvars, y_incgrps;

  // file-specific incvars
  std::map<std::string,std::set<std::string> > file2incvars;
  std::map<std::string,std::set<std::string> > file2excvars;
  
  // file-specific fixed-factors ("fix" in JSON)
  std::map<std::string,std::map<std::string,std::string> > file2fixed;

  // inc/exc bsed on FAC/LVL pairs
  std::map<std::string,std::set<std::string> > incfaclvls;
  std::map<std::string,std::set<std::string> > excfaclvls;
  
  // [read/run] nums (on full manifest) 
  std::vector<std::pair<int,int> > incnums, excnums;

  // [read/run] long-vars (on full manifest) e.g. for dump
  std::set<std::string> inclvars, exclvars;  

  // opts
  bool verbose;
  bool show_xfacs;
  
  // [run] number of permutations;
  //       can be 0, means asymptotic only
  int nreps; 
  bool show_progress;
  
  // adjustments
  bool adj_any;
  bool adj_bonf;
  bool adj_holm;
  bool adj_fdr_bh;
  bool adj_fdr_by;
  
  bool dump_file; // in prep-mode, instead of manifest
  
  // require at least N non-missing obs
  int n_req;
  double n_prop;

  // retain invariants (for dump/manifest mode only)
  bool retain_rows;
  bool retain_cols;
  
  // do perm-corrections independently for each X
  bool correct_all_X;

  
  
};


struct linmod_results_t {  
  // point-wise (empirical) p-values
  std::map<std::string,std::map<std::string,double> > beta;
  std::map<std::string,std::map<std::string,double> > t;
  std::map<std::string,std::map<std::string,double> > p;  
  std::map<std::string,std::map<std::string,double> > emp;
  std::map<std::string,std::map<std::string,double> > emp_corrected;    

  // corrected results
  Eigen::MatrixXd corr;

  // index
  std::map<std::string,std::map<std::string,int> > index;

  // generate corrected results
  void make_corrected( const std::vector<std::string> & xvars ,
		       const std::vector<std::string> & yvars );
		       

  // return corrected results
  double fdr_bh(const std::string & , const std::string & );
  double fdr_by(const std::string & , const std::string & );
  double holm(const std::string & , const std::string & );
  double bonf(const std::string & , const std::string & );

};



struct linmod_t {

  linmod_t()
  {
    ni = ny = nx = nz = 0;
  }
  
  linmod_t( const Eigen::MatrixXd & Y_ ,
	    const std::vector<std::string> & yvars,
	    const Eigen::MatrixXd & X_ ,
	    const std::vector<std::string> & xvars,
	    const Eigen::MatrixXd & Z_ )
  {
    ni = ny = nz = 0;
    set_DV( Y_ );
    set_IV( X_ );
    vname = yvars;
    xname = xvars;
    set_Z( Z_ );
  }
  
  void set_DV( const Eigen::MatrixXd & Y_ );  
  void set_IV( const Eigen::MatrixXd & X_ );
  void set_IV( const Eigen::VectorXd & X_ , const std::string & );
  void set_Z( const Eigen::MatrixXd & Z_ );
  
  //
  // Do the work 
  //

  linmod_results_t run( const int , const bool show_progress = true ) ;

  Eigen::VectorXd get_tstats( const Eigen::VectorXd & B ,
			      const Eigen::MatrixXd & Yres ,
			      const double vx ,
			      const int denom ,
			      Eigen::VectorXd * pvalues = NULL );
  			      
  //
  // Main members
  //

  int ni;
  int ny; // number of DVs (sleep measures, with adjacencies) 
  int nx; // number of X vars
  int nz; 

  std::vector<std::string> vname; // Y vars only
  std::vector<std::string> xname; // X vars

  //
  // Data 
  //  

  // i.e. pairwise single X and Y, but covarying for all Z
  // model:  Y(i) ~ X(j) + Z + e 
  
  // DVs: sleep metrics
  Eigen::MatrixXd Y;
  
  // IVs: (e.g. disease state, 0/1, but tested one at a time) 
  Eigen::MatrixXd X;

  // Covariates
  Eigen::MatrixXd Z;
  
  // permutation matrix
  Eigen::MatrixXi P;

  //
  // helper
  //

  static Eigen::MatrixXd correct( const Eigen::VectorXd & );

  
};



#endif 
