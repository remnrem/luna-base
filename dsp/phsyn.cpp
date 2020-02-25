
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

#include "phsyn.h"
#include "edf/edf.h"
#include "eval.h"

#include "db/db.h"
#include "helper/logger.h"
#include "dsp/hilbert.h"
#include "edf/slice.h"
#include "miscmath/crandom.h"

extern writer_t writer;

extern logger_t logger;

void dsptools::phsyn( edf_t & edf , param_t & param )
{
  
  //
  // Signals
  //

  std::string signal_label = param.requires( "sig" );
  signal_list_t signals = edf.header.signal_list( signal_label );  
  const int ns = signals.size();

  
  //
  // Frequencies
  //
  
  std::vector<double> lf0 = param.dblvector( "lwr" );
  std::vector<double> uf0 = param.dblvector( "upr" );
  
  std::vector<freq_range_t> lf;
  std::vector<freq_range_t> uf;

  double w = param.requires_dbl( "w" );

  if ( lf0.size() < 1 || lf0.size() > 3 ) Helper::halt( "expecting lwr to have 1,2 or 3 values: start,end,step" );
  if ( uf0.size() < 1 || uf0.size() > 3 ) Helper::halt( "expecting upr to have 1,2 or 3 values: start,end,step" );
  
  if ( lf0.size() == 1 ) lf0.push_back( lf0[0] );
  if ( lf0.size() == 2 ) lf0.push_back( 1 );
  if ( uf0.size() == 1 ) uf0.push_back( uf0[0] );
  if ( uf0.size() == 2 ) uf0.push_back( 1 );
  if ( lf0[0] > lf0[1] ) Helper::halt( "bad format for lwr" );
  if ( uf0[0] > uf0[1] ) Helper::halt( "bad format for upr" );
  
  if ( lf0[0] - w < 0 ) Helper::halt( "bad format for lwr, lower value too low given w" );
  if ( uf0[0] - w < 0 ) Helper::halt( "bad format for upr, lower value too low given w" );

  for (double v = lf0[0] ; v <= lf0[1] ; v += lf0[2] )
    lf.push_back( freq_range_t( v-w , v+w ) );

  for (double v = uf0[0] ; v <= uf0[1] ; v += uf0[2] )
    uf.push_back( freq_range_t( v-w , v+w ) );

  //
  // Sampling rates
  //

  std::vector<double> Fs = edf.header.sampling_freq( signals );

  //
  // Filter
  // 

  double ripple = param.has( "ripple" ) ? param.requires_dbl( "ripple" ) : 0.05;

  double tw = param.has( "tw" ) ? param.requires_dbl( "tw" ) : 2 ;

  
  //
  // Number of bins, rplicates
  //

  const int nbins = param.has( "nbin" ) ? param.requires_int( "nbin" ) : 20;
  const int nreps = param.has( "nrep" ) ? param.requires_int( "nrep" ) : 1000;


  //
  // requires epochs 
  //

  if ( ! edf.timeline.epoched() ) 
    Helper::halt( "requires EPOCH'ed data" );
  
  const double epoch_sec = edf.timeline.epoch_length();

  bool no_epoch_perm = param.has( "no-epoch-perm" );
  

  //
  // function only applies to entire current trace
  //
  
  interval_t interval = edf.timeline.wholetrace();

  
  //
  // Interate over each signal
  //

  for (int s=0;s<ns;s++)
    {
      
      //
      // Only data-channels
      //

      if ( edf.header.is_annotation_channel( signals(s) ) ) continue;	  
      

      // epoch size in sample points (if using epoch-perms)
      
      const int es = no_epoch_perm ? 0 : Fs[s] * epoch_sec;


      //
      // Get signal
      //
      
      slice_t slice( edf , s , interval );
      
      const std::vector<double> * d = slice.pdata();


      //
      // Calculate
      //

      phsyn_t phsyn( *d , Fs[s] , lf , uf , nbins , nreps , ripple , tw , es );
      
      phsyn.calc();


      //
      // Report
      //

      //writer.level( signals.label(s) , globals::signal_strat );

      //writer.value( "PHSYN" , phsyn.
      
    }
  
  //writer.unlevel( globals::signal_strat );


}


void phsyn_t::calc() 
{
  

  //
  // bin boundaries (0..360)
  //
  
  double bs = 360 / (double)nbins;  
  
  std::vector<double> bb( nbins );
  
  bb[0] = 0;
  
  for (int i=1;i<nbins;i++) bb[i] = bb[i-1] + bs;
    
  //
  // size output grids
  //
  
  obs.resize( nbins );
  perm.resize( nbins );
  pv.resize( nbins );
  z.resize( nbins );
  z2.resize( nbins );
  
  for (int b=0;b<nbins;b++) 
    {      
      obs[b].resize( nbins , 0 );
      perm[b].resize( nbins , 0 );
      pv[b].resize( nbins , 0 );
      z[b].resize( nbins , 0 );
      z2[b].resize( nbins , 0 );
    }
  

  //
  // get phase for each f1 and f2 frequency band  
  //
  
  std::map<freq_range_t,std::vector<double> > ph;  
  
  std::set<freq_range_t> frqs;
  for (int f=0;f<f1.size();f++) frqs.insert( f1[f] );
  for (int f=0;f<f2.size();f++) frqs.insert( f2[f] );
  
  std::set<freq_range_t>::const_iterator ff = frqs.begin();
  
  logger << "  considering " << frqs.size() << " frequency bands\n";

  //
  // get phase angles at each time point
  //  

  while ( ff != frqs.end() )
    {
      
      logger << "  hilbert " << ff->first << "-" << ff->second << "\n";

      // filter-Hilbert
      hilbert_t hilbert( x , sr , ff->first , ff->second , ripple , tw );
      
      // convert to degrees with 0 as pos-to-neg crossing
      std::vector<double> angle = *hilbert.phase();
      for (int i=0;i<angle.size();i++) angle[i] = MiscMath::as_angle_0_pos2neg( angle[i] );
      // store
      ph[ *ff ] = angle;            
      ++ff;
    }
  

  const int npoints = x.size();
  
  //
  // create observed frequency bin
  //
        
  for (int i1=0;i1<f1.size();i1++) 
    for (int i2=0;i2<f2.size();i2++) 
      {

	logger << f1[i1].first << " and " << f2[i2].first << "\n";

	//
	// Clear working count matrix
	//
	
	for (int b1=0;b1<nbins;b1++)
	  for (int b2=0;b2<nbins;b2++)
	    obs[b1][b2] = perm[b1][b2] = z[b1][b2] = z2[b1][b2] = pv[b1][b2] = 0 ;
	
	const std::vector<double> * ph1 = &ph[ f1[i1] ];
	const std::vector<double> * ph2 = &ph[ f2[i2] ];
	
	int b1 = 0, b2 = 0; 
	for (int i=0;i<npoints;i++)
	  {
	    bin( (*ph1)[i] , &b1 , bb , nbins );
	    bin( (*ph2)[i] , &b2 , bb , nbins );

	    // increment by one unit (default, unweighted)
	    obs[ b1 ][ b2 ]++;
	    
	  }
		

	//
	// Observed summary statistic
	//

	double obs_stat = test_uniform( obs );

	
	//
	// Replicates
	//
	
	int emp_stat = 0;
	std::vector<double> perm_stats;

	for (int r = 0; r < nreps; r++) 
	  {

// 	    if ( r % 10 == 0 ) logger << ".";
// 	    if ( r % 100 == 0 ) logger << "\n";
	    
	    // clear working perm matrix
	    
	    for (int b1=0;b1<nbins;b1++)
	      for (int b2=0;b2<nbins;b2++)
		perm[b1][b2] = 0;

	    // initially, whole trace perm
	    
	    const std::vector<double> * ph1 = &ph[ f1[i1] ];
	    const std::vector<double> * ph2 = &ph[ f2[i2] ];
	
	    int b1 = 0, b2 = 0; 

	    // within-epoch perm
	    
	    std::vector<int> permpos( npoints );
	    int j = 0;

	    if ( es ) 
	      {
		int num_epoch = npoints / es ; 
		std::vector<int> shuffle;
		for (int e=0;e<num_epoch;e++) shuffle.push_back( CRandom::rand( es ) );
		//std::cout << "fig " << num_epoch << " epochs\n";
		for (int i=0;i<npoints;i++)
		  {		    
		    int cur_epoch = i / es;
		    if ( cur_epoch >= num_epoch ) continue; // do not permute last epoch, if partial
		    int off = i - cur_epoch * es; 
		    if ( off + shuffle[ cur_epoch ] >= es ) permpos[i] = i - es + shuffle[cur_epoch ];
		    else permpos[i] = i + shuffle[cur_epoch ];
		    //std::cout << "pp\t" << cur_epoch << "\t" << i << "\t" << permpos[i] << "\t" << off << " " << shuffle[cur_epoch ] << "\n";
		  }
	      }
	    else  // otherwise whole trace perm
	      j = CRandom::rand( npoints );
	    

	    for (int i=0;i<npoints;i++)
	      {
		
		bin( (*ph1)[i] , &b1 , bb , nbins );
		
		if ( es ) // within-epoch shuffle
		  j = permpos[i];
		else // whole trace shuffle
		  {
		    ++j;
		    if ( j == npoints ) j = 0;
		  }

		bin( (*ph2)[j] , &b2 , bb , nbins );

		// increment by one unit (default, unweighted)
		perm[ b1 ][ b2 ]++;
	      }

	    

	    // 
	    // summary stat for surrogate 
	    //

	    double perm_stat = test_uniform( perm );
	    
	    perm_stats.push_back( perm_stat );

	    if ( perm_stat >= obs_stat ) ++emp_stat;

	    //
	    // accumulate point level pv/z file
	    //

	    for (int b1=0;b1<nbins;b1++)
	      for (int b2=0;b2<nbins;b2++)
		{
		  if ( perm[b1][b2] >= obs[b1][b2] ) pv[b1][b2]++;
		  z[b1][b2] += perm[b1][b2];
		  z2[b1][b2] += perm[b1][b2] * perm[b1][b2];
		}
	    
	    // next replicate
	  
	  }
	

	double z_stat_mean = MiscMath::mean( perm_stats );
	double z_stat_sd   = MiscMath::sdev( perm_stats );
	double z_stat = ( obs_stat - z_stat_mean ) / z_stat_sd ; 
	
	std::cout << f1[i1].first << "-" << f1[i1].second << "\t"
		  << f2[i2].first << "-" << f2[i2].second << "\t"
	  //<< obs_stat << "\t"
		  << z_stat << "\t"
		  << (emp_stat+1)/double(nreps+1) << "\n";



	if ( true )
	  {
	    
	    for (int b1=0;b1<nbins;b1++)
	      for (int b2=0;b2<nbins;b2++)
		{
		  //std::cout << "ha " << z[b1][b2] << "  " << (double)nreps << "\n";
		  double zmean = z[b1][b2]/(double)nreps;
		  double zsd   = sqrt(   z2[b1][b2]/(double)nreps  - zmean * zmean );
		  
		  std::cout << "res " << b1 << " " << b2 << " " 
			    << obs[b1][b2] << " " 		
			    << (pv[b1][b2] + 1 )/double(nreps+1) << " " 
			    << zmean << " " 
			    << zsd << " " 
			    << (obs[b1][b2] - zmean ) / zsd << "\n";
		}
	  }



	//
	// Next pair of frequencies
	//

      }
  
}

bool phsyn_t::bin( double d , int * b , const std::vector<double> & th , const int nbins )
{
  // expect 0 to 360 
  if ( d < 0 || d > 360 ) return false;
  if ( *b < 0 || *b >= nbins ) return false;
  
  // bin 0 is 0 degree
  // last bin is 360 
  
  while ( 1 ) 
    {      
      //      std::cerr << "b = " << *b << " " << th[*b] << " " << th[ (*b ) + 1]  << "\t" << d << "\n";

      if ( *b == nbins - 1 )
	{
	  if ( d >= th[ *b ] ) return true;
	  *b = 0;
	}
      
      if ( d >= th[*b] && d < th[ (*b ) + 1] ) return true;      

      ++(*b);

      if ( *b == nbins ) *b = 0;
    }
  
  return false;
}


double phsyn_t::test_uniform( const std::vector<std::vector<double> > & m )
{
  // stat = ( O - E )^2 

  const int bs = m.size();
  
  std::vector<double> rows( bs , 0 );
  std::vector<double> cols( bs , 0 );
  double tot = 0;
  for (int b1=0;b1<bs;b1++)
    for (int b2=0;b2<bs;b2++)
      {
	rows[b1] += m[b1][b2];
	cols[b2] += m[b1][b2];
	tot += m[b1][b2];
      }

  double stat = 0; 
    
  for (int b1=0;b1<bs;b1++)
    for (int b2=0;b2<bs;b2++)
      {
	double exp = ( rows[b1] * cols[b2] ) / tot;
	stat += ( m[b1][b2] - exp ) * ( m[b1][b2] - exp );
      }
  return stat;
    
}