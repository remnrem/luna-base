
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

#ifndef __AROUSALS_H__
#define __AROUSALS_H__

#include "stats/Eigen/Dense"
#include "intervals/intervals.h"

#include <map>
#include <set>

struct edf_t;
struct param_t;
struct signal_list_t;
struct gaussian_hmm_t;

struct arousals_t {
  
  arousals_t( edf_t & edf , param_t & param );
  
  bool hjorth( const Eigen::VectorXd & x ,
	       double * , double * , double * ,
	       const bool mean_center = false ) const;
  
private:

  void build_ftr_matrix( edf_t & edf ,
			 const signal_list_t & eeg_signals ,
			 const signal_list_t & emg_signals ,
			 Eigen::MatrixXd & Xeeg ,
			 Eigen::MatrixXd & Xemg , 
			 std::vector<int> * state ,
			 std::vector<int> * sequence ,
			 std::vector<double> * sec
			 );

  Eigen::VectorXd calc_eeg_ftrs( const Eigen::MatrixXd & X );
  
  Eigen::VectorXd calc_emg_ftrs( const Eigen::MatrixXd & X , const std::vector<double> &  );

  Eigen::MatrixXd process_ftr_matrix( Eigen::MatrixXd * Xeeg , Eigen::MatrixXd * Xemg , const std::vector<int> & st );
  
  void dump( const std::vector<std::vector<std::vector<Eigen::VectorXd> > > & X ,
	     const std::vector<std::vector<std::vector<double> > > & tt ) const;
  
  std::vector<std::vector<std::vector<Eigen::VectorXd> > > assemble( const Eigen::MatrixXd & Xftr , 
								     const std::vector<int> & state ,
								     const std::vector<int> & seq ,
								     const std::vector<double> & sec ,
								     std::vector<std::vector<std::vector<double> > > * tt );

  std::vector<std::vector<std::vector<Eigen::VectorXd> > >
  extract( const std::vector<std::vector<std::vector<Eigen::VectorXd> > > & X , const std::vector<int> & ex );

  // helper
  Eigen::VectorXd robust_mad_norm( const Eigen::VectorXd & x , const std::vector<int> & st );

  std::vector<std::pair<int,int>>
  merge_events_with_gap_sorted(const std::vector<std::pair<int,int>> &events,
			       int max_gap);

  std::vector<std::pair<int,int>> mask_to_intervals(const std::vector<bool> &mask);
  
  
  std::map<std::string,std::set<interval_t> > event_heuristic( const std::vector<std::vector<std::vector<Eigen::VectorXd> > > & X ,
							       const std::vector<std::vector<std::vector<double> > > & tt );
  


  
  // map states
  void map_states( const Eigen::MatrixXd & mu , int*,int*,int*, bool nrem) const;
  
  // given HMM results, make annotations / channels
  int annotate( const int state ,
		std::vector<std::vector<double> > & tt ,
		std::vector<std::vector<Eigen::VectorXd> > & X , 
		std::vector<std::vector<int> > & paths ,
		std::vector<Eigen::MatrixXd> & posteriors ,
		int,int,int,
		const std::string & aname );

  void add_annot( int idx ,   			  
		  std::vector<std::vector<int> > & paths ,
		  std::vector<std::vector<double> > & tt ,			    
		  const std::string & class_label ,
		  const std::string & inst_label );

  void add_channels( const std::vector<std::vector<Eigen::VectorXd> > & X ,
		     const std::vector<std::vector<double> > & tt ,
		     const std::string & ch_prefix );
    
    
  // k-means based seed for HMM
  void init_kmeans_hmm( gaussian_hmm_t & hmm ,
			const std::vector<std::vector<Eigen::VectorXd>> & sequences );
  
  
  // parent EDF
  edf_t * parent;
  
  // assume single, fixed SR for EEG
  int sr;
  
  // for each contiguous segment
  std::vector<Eigen::VectorXd> nrem_ftr, rem_ftr;
  std::vector<uint64_t> nrem_times, rem_times;


  inline Eigen::MatrixXd stack_sequences( const std::vector<std::vector<Eigen::VectorXd>> &sequences,
					  int dim)
  {
    // Count total samples
    std::size_t total = 0;
    for (const auto &seq : sequences)
      total += seq.size();
    
    Eigen::MatrixXd X(static_cast<int>(total), dim);
    int row = 0;
    for (const auto &seq : sequences) {
      for (const auto &x : seq) {
	if (x.size() != dim)
	  throw std::runtime_error("stack_sequences: dim mismatch");
	X.row(row++) = x.transpose();
      }
    }
    return X;
  }
  
  
};


#endif
