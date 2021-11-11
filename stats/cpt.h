
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

#ifndef __CPT_H__
#define __CPT_H__

struct param_t;
struct clocs_t;

#include "stats/Eigen/Dense"
#include <vector>
#include <map>
#include <set>

void cpt_wrapper( param_t & ) ;


struct cpt_results_t {

  // point-wise empirical p-values

  
  std::map<std::string,double> beta;
  std::map<std::string,double> t;
  std::map<std::string,double> emp;
  std::map<std::string,double> emp_corrected;  

  // point is in a cluster?
  std::map<std::string,int> inclst;
  
  std::map<std::string,double> cluster_emp; // seed -> p-value
  std::map<std::string,std::set<std::string> > cluster_members; // seed-> members
  
  
};


struct cpt_cluster_t {
  
  double stat;
  int seed; 
  std::set<int> members;

  int emp; // p-value 
  
  bool operator< (const cpt_cluster_t & rhs ) const
  {
    if ( stat < rhs.stat ) return true;
    if ( stat > rhs.stat ) return false;
    return seed < rhs.seed;
  }
  
};

struct cpt_clusters_t {

  cpt_clusters_t( const Eigen::VectorXd & T ,
		  double threshold ,
		  const std::map<int,std::set<int> > & adj , 
		  bool two_sided = true , 
		  bool verbose = false , 
		  const std::vector<std::string> * labels = NULL 
		  ); 

  
  void update( double pt )
  {
    // always iterate over same observed set, so can use 'perm'
    // as below and this will track with the same cluster...
   
    int cnt = 0;
    std::set<cpt_cluster_t>::iterator cc = clusters.begin();
    while ( cc != clusters.end() )
      {
	if ( pt >= cc->stat ) ++perm[cnt];
	++cc;
	++cnt;
      }    
  }
  
  // populate
  double max_stat;
  std::set<cpt_cluster_t> clusters;
  std::vector<double> perm;
  
};

struct cpt_t {   

  cpt_t()
  {
    ni = ny = nz = 0;
  }

  cpt_t( const Eigen::MatrixXd & Y_ ,
	 const Eigen::VectorXd & X_ ,
	 const Eigen::MatrixXd & Z_ )
  {
    ni = ny = nz = 0;
    set_DV( Y_ );
    set_IV( X_ );
    set_Z( Z_ );
  }
    
  void set_DV( const Eigen::MatrixXd & Y_ );
  
  void set_IV( const Eigen::VectorXd & X_ );

  void set_Z( const Eigen::MatrixXd & Z_ );


  //
  // Adjacencies
  //


  void calc_adjacencies( const std::vector<std::string> & vname ,
			 const std::map<std::string,std::string> & col2var,
			 const std::map<std::string,double> & col2f,
			 const std::map<std::string,double> & col2t,
			 const std::map<std::string,std::string> & col2ch1,
			 const std::map<std::string,std::string> & col2ch2,
			 double fth ,
			 double tth , 
			 clocs_t * clocs ,
			 double sth , 
			 bool dump_adj );
  
  

  
  //
  // Do the work 
  //

  cpt_results_t run( int , double , bool , bool ) ;

  Eigen::VectorXd get_tstats( const Eigen::VectorXd & B ,
			      const Eigen::MatrixXd & Yres ,
			      const double vx ,
			      const int denom );
  			      
  
  //
  // Main members
  //

  int ni;
  int ny; // number of DVs (sleep measures, with adjacencies) 
  int nz; // 

  //
  // Data 
  //

  // model:  Y ~ X + Z
    
  // DVs: sleep metrics
  Eigen::MatrixXd Y;
  
  // Single IV (e.g. disease state, 0/1) 
  Eigen::VectorXd X;

  // Covariates
  Eigen::MatrixXd Z;
  
  // permutation matrix
  Eigen::MatrixXi P;

  

  //
  // Helpers
  //

  //  void set_P();
  

  //
  // Cluster based 
  //

  std::vector<std::string> vname;
  
  // adjacenies [ columns of X ]
  std::map<int,std::set<int> > adjacencies;

  

};

#endif 
