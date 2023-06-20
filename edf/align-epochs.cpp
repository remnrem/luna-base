
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

#include "edf/align-epochs.h"

#include "edf/edf.h"
#include "edf/slice.h"
#include "db/db.h"
#include "helper/logger.h"
#include "helper/helper.h"
#include "stats/eigen_ops.h"

extern writer_t writer;
extern logger_t logger;

align_epochs_t::align_epochs_t( edf_t & edf , param_t & param )
{
  
  // reduced EDF to align
  const std::string edffile2 = param.requires( "edf" );
  
  edf_t edf2;
  
  bool okay = edf2.attach( edffile2 , "." );

  if ( ! okay ) Helper::halt( "could not attach " + edffile2 ) ;

  //
  // signal to use (pull from attached, initial EDF)
  //
  
  const std::string signal_label = param.value( "sig" );

  const bool no_annotations = true;
  
  signal_list_t signals = edf.header.signal_list( signal_label , no_annotations);
  
  ns = signals.size();
  
  if ( ns == 0 ) return;
  
  std::vector<double> Fs = edf.header.sampling_freq( signals );

  //
  // Check all exist, with similar SR, in the second EDF
  //

  for (int s=0; s<ns; s++)
    {
      const std::string slab = signals.label(s);
      if ( ! edf2.header.has_signal( slab ) )
	Helper::halt( "could not find " + slab + " in " + edffile2 );

      const double sr2 = edf2.header.sampling_freq( signals(s) );

      if ( fabs( sr2 - Fs[s] ) > 0.01 )
	Helper::halt( "different sample rate for " + slab + " between EDFs" );
      
    }
  
  //
  // threshold for similarity
  //

  th = param.has( "th" ) ? param.requires_dbl( "th" ) : 2.0 ; 
  if ( th <= 0 ) Helper::halt( "expecting positive value for th" );

  th2 = param.has( "th2" ) ? param.requires_dbl( "th2" ) : 0.2 ; 
  if ( th2 <= 0 ) Helper::halt( "expecting positive value for th2" );

  logger << "  matching epochs when best match is " << th << " SD units below the mean, and " << th2 << " SD units better than next best\n";
    
  //
  // Assume same order of epochs 
  //

  assume_order = param.has( "ordered" ) ? param.yesno( "ordered" ) : true;

  //
  // Attempt to fix, i.e. if abother epoch matches better, take second best
  //
  
  resolve_order = param.has( "resolve" ) ? param.yesno( "resolve" ) : true;
  resolve_order = false;
  
  //
  // verbose output for a single E2 epoch
  //

  verbose = param.has( "verbose" ) ? param.requires_int( "verbose" ) : -1 ; 

  //
  // Looks all good to align
  //
  
  edf.timeline.ensure_epoched();
  edf2.timeline.ensure_epoched();

  ne = edf.timeline.first_epoch();
  ne2 = edf2.timeline.first_epoch();

  logger << "  aligning " << ne2 << " epochs from " << edffile2 << " to the in-memory EDF containing " << ne << " epochs\n";

  if ( ne2 >= ne ) Helper::halt( "expecting " + edffile2 + " to have fewer epochs than the primary attached EDF" );

  // to track discontinuities
  std::vector<interval_t> intervals;

  E1.clear();
  E2.clear();
  
  // get data
  while (1)
    {
      int epoch = edf.timeline.next_epoch();
      
      if ( epoch == -1 ) break;

      E1.push_back( epoch );
      
      interval_t interval = edf.timeline.epoch( epoch );

      intervals.push_back( interval );
      
      eigen_matslice_t mslice( edf , signals , interval );
      
      Eigen::MatrixXd X = mslice.data_ref();
      
      // T = center , F = winsor , 0 = W(winsor) , F = second_rescale , T = ignore invariants (set to 0/VAR1)
      eigen_ops::robust_scale( X , true , false , 0 , false , true );

      X1[ epoch ] = X;
      
    }

  logger << "  extracted " << ns << " signals for " << ne << " epochs from main EDF\n";
  
  // get data from second EDF

  while (1)
    {
      int epoch = edf2.timeline.next_epoch();
      
      if ( epoch == -1 ) break;
      
      E2.push_back( epoch );

      interval_t interval = edf2.timeline.epoch( epoch );
      
      eigen_matslice_t mslice( edf2 , signals , interval );

      Eigen::MatrixXd X = mslice.data_ref();

      // T = center , F = winsor , 0 = W(winsor) , F = second_rescale , T = ignore invariants (set to 0/VAR1)
      eigen_ops::robust_scale( X , true , false , 0 , false , true );
      
      X2[ epoch ] = X;

    }

  logger << "  extracted " << ns << " signals for " << ne2 << " epochs from " << edffile2 << "\n";
    
  mapping.clear();
  mapping2.clear();
  
  if ( E1.size() != ne  ) Helper::halt( "internal error constructing epoch list" );
  if ( E2.size() != ne2 ) Helper::halt( "internal error constructing epoch list for " + edffile2 );

  //
  // Do mapping
  //
  
  std::map<int,double> scs, nxts;
  std::map<int,std::set<int> > rmapping;
  
  // consider each E2 epoch and find best match
  // best_match returns -1 if there is not a single good match

  for (int e2idx = 0; e2idx < ne2 ; e2idx++ )
    {
      int e2 = E2[e2idx];
      double sc, nxt;
      int nxte;

      int e1 = best_match( e2 , &sc, &nxt , &nxte );
      mapping[ e2 ] = e1;
      mapping2[ e2 ] = nxte;
      scs[ e2 ] = sc;
      nxts[ e2 ] = nxt;
      
      // reverse mapping for check below (resolve)
      if ( e1 != -1 )
	{
	  rmapping[ e1 ].insert( e2 );
	  if ( rmapping[ e1 ].size() > 1 )
	    logger << "  *** warning: epoch " << edf.timeline.display_epoch( e1 ) << " preliminarily mapped to multiple epochs in " << edffile2 << "\n";
	}
    }
  
  //
  // Optionally, have a quick go at resolving any clear ambiguities
  //   - if same E1 epoch is assigned to >1, take second best for lower score

  if ( resolve_order )
    {
      int resolved = 0;

      for (int e2idx=0; e2idx < ne2; e2idx++)
	{
	  // E2
	  const int e2 = E2[ e2idx ];

	  // proposed mapping
	  const int e1 = mapping[ e2 ] ;
	  
	  // was this E1 also mapped to other E2 epochs?
	  if ( rmapping[ e1 ].size() > 1 )
	    {
	      // which mapped better? just get absolute distances between E1 and all E2 candidates

	      logger << "  attempting to resolve multiple mappings of " << edf.timeline.display_epoch( e1 ) << "\n";
	      
	      const std::set<int> e2s = rmapping[ e1 ];
	      
	      std::vector<double> dists;
	      double dmin = 9999;
	      int imin = -1;
	      
	      std::set<int>::const_iterator ee = e2s.begin();
	      while ( ee != e2s.end() )
		{
		  const double d = dist( e1 , *ee ) ;
		  if ( imin == -1 || d < dmin )
		    {
		      imin = *ee;
		      dmin = d;
		    }

		  dists.push_back( d );

		  // second best for this e2
		  const int e22 = mapping2[ *ee ];
		  const double d2 = dist( e22 , *ee );

		  int de1 = edf.timeline.display_epoch( e1 );
		  int de2 = edf2.timeline.display_epoch( *ee );
		  int de22 = edf2.timeline.display_epoch( e22 );
		  
		  logger << "   d( " << de1 << ", " << de2 << " ) = " << d << "; "
			 << "  second-best " << de2 << " d( " << de22 << ", " << de2 << " ) = " << d2 << "\n";
		    
		  
		  
		  ++ee;
		}
	      
	    }
	}
    }
  
  
  //
  // check that we have a 1:1 ordering
  //
  
  bool order_okay = true;
  int last_epoch = -2;
  int last_e2 = -1;
  
  std::map<int,int>::const_iterator ii = mapping.begin();
  while ( ii != mapping.end() )
    {
      // only check if mapping
      if ( ii->second != -1 )
	{
	  if ( ii->second <= last_epoch )
	    {
	      order_okay = false;
	      int ee2 = edf2.timeline.display_epoch( ii->first ) ;
	      int ee1 = edf.timeline.display_epoch( ii->second ) ;

	      int eprior2 = edf2.timeline.display_epoch( last_e2 );
	      int eprior = edf.timeline.display_epoch( last_epoch );

	      double sc = scs[ ii->first ] ;
	      double nxt = nxts[ ii->first ] ;
	      int nxte = mapping2[ ii->first ] != -1 ? edf.timeline.display_epoch( mapping2[ ii->first ] ) : -1;
		
	      double sc2 = scs[ last_e2 ] ;
	      double nxt2 = nxts[ last_e2 ] ;
	      int nxte2 = mapping2[ last_e2 ] != -1 ? edf.timeline.display_epoch( mapping2[ last_e2 ] ) : -1;
	      
	      logger << "  epochs aligned out-of-order:\n"
		     << "   " << eprior2 << " in " << edffile2 << " mapped to " << eprior 
		     << " [ score = " << sc2 << " and next-best = " << nxt2 << ", epoch " << nxte << "]\n"
		     << "   " << ee2 << " in " << edffile2 << " mapped to " << ee1 
		     << " [ score = " << sc << " and next-best = " << nxt << ", epoch " << nxte2 << "]\n";
	      
	      if ( assume_order ) 
		Helper::halt( "alignment violated ordering assumption (use ordered=F for ignore)");
	    }
	  last_epoch = ii->second;
	  last_e2 = ii->first ; 
	}
      ++ii;
    }
  
  if ( order_okay ) 
    logger << "  all epochs aligned in correct, increasing order\n";
     
  //
  // output
  //
  
  int count = 0 ;
  int count1 = 0 , count2 = 0;
  
  std::map<int,int>::const_iterator ee = mapping.begin();
  while ( ee != mapping.end() )
    {
      writer.level( edf2.timeline.display_epoch( ee->first ) , "E2" );
      if ( ee->second != -1 ) 
	{
	  writer.value( "E1" , edf.timeline.display_epoch( ee->second )  );
	  ++count;
	}
      else
	{
	  if ( scs[ ee->first ] > -th ) ++count1;
	  if ( nxts[ ee->first ] < th2 ) ++count2;
	}

      if ( mapping2[ ee->first ] != -1 )
	writer.value( "NEXT_E1" , edf.timeline.display_epoch( mapping2[ ee->first ] )  );
      
      writer.value( "D" , scs[ ee->first ] );
      writer.value( "NEXT" , nxts[ ee->first ] ) ;
      ++ee;
    }
  writer.unlevel( "E2" );

  logger << "  aligned " << count << " of " << ne2 << " epochs confidently\n";
  logger << "    " << count1 << " epochs (" << ceil((count1/(double)ne2)*100) << "%) failed based on th = " << th << "\n"
	 << "    " << count2 << " epochs (" << ceil((count2/(double)ne2)*100) << "%) failed based on th2 = " << th2 << "\n";
    
  //
  // copy annotations from original into the new space
  //
  
  
  
}


int align_epochs_t::best_match( const int e2 , double * sc , double * nxt , int * nxte ) const
{

  Eigen::VectorXd D = Eigen::VectorXd::Zero( ne );

  for (int e1idx=0; e1idx<ne; e1idx++)
    {
      const int e1 = E1[e1idx];
      D[ e1idx ] = dist( e1 , e2 );      
    }
  
  int idx = 0;
  const double dmin = D.minCoeff(&idx);

  const double dmean = D.mean();
  
  const double sd = sqrt( ( D.array() - dmean ).square().sum() / (D.size() - 1 ) );
  
  const double dt = dmean - th * sd ;
  
  // get second lowest
  int prior_idx = -1;
  const double prior_min = idx != 0 ? D.head( idx ).minCoeff( &prior_idx ) : 999999;
  // get ne - idx --> take 1 from end, ne = 10 , idx = 8, 
  int past_idx;
  const double past_min = idx != ne - 1 ? D.tail( ne - idx - 1 ).minCoeff( &past_idx ) : 999999;
  
  const double second_best = prior_min < past_min ? prior_min : past_min ;
  const int second_best_epoch = prior_min < past_min ? prior_idx : past_idx;  
   

  // save (next best is in SD units)
  *sc = ( dmin - dmean ) / sd ;
  *nxt = ( second_best - dmean ) / sd ;  
  // express as how much worse next-best fit was (in SD units)
  *nxt = *nxt - *sc;

  // matched?
  const bool matched = dmin <= dt && *nxt >= th2;  
  
  *nxte = matched ? second_best_epoch : -1;
  return matched ? idx : -1;

}

double align_epochs_t::dist( const int e1, const int e2 ) const
{
  double d = 0;

  if ( X1.find(e1) == X1.end() ) Helper::halt( "could not find " + Helper::int2str( e1 ) + " in dataset #1" );
  if ( X2.find(e2) == X2.end() ) Helper::halt( "could not find " + Helper::int2str( e2 ) + " in dataset #2" );
    
  const Eigen::MatrixXd & X = X1.find(e1)->second;
  const Eigen::MatrixXd & Y = X2.find(e2)->second;

  for (int s=0; s<ns; s++)
    d += (X.col(s) - Y.col(s)).squaredNorm();

  if ( d == 0 ) d = DEPS;
  
  return log( d ) ;
}

