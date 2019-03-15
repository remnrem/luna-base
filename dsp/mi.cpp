
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

#include "mi.h"

#include "helper/helper.h"
#include "eval.h"
#include "miscmath/miscmath.h"
#include "miscmath/crandom.h"
#include "edf/edf.h"
#include "edf/slice.h"
#include "db/db.h"

#include <cmath>



extern writer_t writer;


mi_t::mi_t( const std::vector<double> & a , const std::vector<double> & b )
{    
  eps = 1e-60;
  if ( b.size() != a.size() ) Helper::halt("unequal sequence length in MI");
  n = a.size();
  da = a; db = b;    
}

// use previous calculated bins (i.e. done on whole datasets
// so that same grid can be used per epoch

void mi_t::force_thresholds( const std::vector<double> & _tha, 
			     const std::vector<double> & _thb )     
{
  if ( _tha.size() != _thb.size() ) Helper::halt("problem in mi_t::force_thresholds()");
  tha = _tha;
  thb = _thb;
  nbins = tha.size();     
  bin_data();
}


void dsptools::compute_mi( edf_t & edf , param_t & param )
{
  
  
  //
  // Signals
  //
  
  std::string signal_label = param.requires( "sig" );
  
  signal_list_t signals = edf.header.signal_list( signal_label );  
  
  const int ns = signals.size();
  
  
  
  
  //
  // Epochs or whole signal?
  //

  bool epoched = edf.timeline.epoched() && param.has("epoch") ;

  //
  // # of bin rules
  //

  int bin_rule = 1; // FD
  if ( param.has( "scott" ) ) bin_rule = 2;
  else if ( param.has("sturges" ) ) bin_rule = 3;

  int fixed_nbins = 0;
  if ( param.has( "nbins" ) ) fixed_nbins = param.requires_int( "nbins" );


  //
  // Permuation?
  //
  
  int nperms = param.has("permute") ? param.requires_int( "permute" ) : 0 ;


  //
  // Key output variables
  //

  writer.var( "MI" , "Mutual information" );
  
  if ( nperms )
    writer.var( "PMI" , "Empirical p-value, MI" );


  //
  // Ensure similar sampling rates
  //
  
  for (int i=0;i<ns-1;i++)
    {
      
      if ( edf.header.is_annotation_channel( signals(i) ) ) continue;
      
      for (int j=i+1;j<ns;j++)
	{
	  
	  if ( edf.header.is_annotation_channel( signals(j) ) ) continue;
	  
	  //
	  // MI 
	  //
	  
	  const int sr1 = edf.header.sampling_freq( signals(i) );
	  const int sr2 = edf.header.sampling_freq( signals(j) );	      
	  if ( sr1 != sr2 ) Helper::halt( "MI requires similiar sampling rates" );
	  

	  //
	  // stratify output by SIGNALS
	  //
	  
	  writer.level( signals.label(i) + "x" + signals.label(j) , "CHS" );
	  
	  //
	  // First, MI for entire duration
	  //
	  	  
	  interval_t interval = edf.timeline.wholetrace();
	  
	  slice_t slice1( edf , signals(i) , interval );
	  slice_t slice2( edf , signals(j) , interval );
	  
	  const std::vector<double> * d1 = slice1.pdata();
	  const std::vector<double> * d2 = slice2.pdata();

	  mi_t mi( *d1 , *d2 );

	  if ( fixed_nbins )
	    mi.set_nbins( fixed_nbins );
	  else
	    {
	      if      ( bin_rule == 1 ) mi.set_nbins_fd();
	      else if ( bin_rule == 2 ) mi.set_nbins_scott();
	      else if ( bin_rule == 3 ) mi.set_nbins_sturges();
	    }

	  // given # of bins, set thresholds
	  mi.set_thresholds( bin_rule );

	  // calculate MI
	  mi.calc();
	  
	  // empirical p-values
	  if ( nperms ) 
	    {
	      double pemp, pz;
	      mi.permute( nperms , &pemp , &pz );
	      writer.value( "EMP" , pemp );
	      writer.value( "Z" , pz );
	    }

	  
	  //
	  // Output
	  //
	  
	  writer.value( "MI" , mi.mutinf );
	  writer.value( "JINF" , mi.jointinf );
	  writer.value( "TOTCORR" , mi.total_corr );
	  writer.value( "DTOTCORR" , mi.dual_total_corr );
	  writer.value( "INFA" , mi.infa );
	  writer.value( "INFB" , mi.infb );
	  writer.value( "NBINS" , mi.nbins );

	  
	  //
	  // Epochs ? 
	  //
	  
	  if ( epoched ) 
	    {
	      
	      // get thresholds from total analysis
	      const std::vector<double> & tha = mi.tha;
	      const std::vector<double> & thb = mi.tha;
	      
	      int ne = edf.timeline.first_epoch();      
	      
	      while ( 1 ) 
		{

		  //
		  // Stratify by epoch 
		  //

		  int epoch = edf.timeline.next_epoch();      
		  if ( epoch == -1 ) break;
		  
		  
		  interval_t interval = edf.timeline.epoch( epoch );
		  
		  writer.epoch( edf.timeline.display_epoch( epoch ) );		  				
		  
		  //
		  // Get data
		  //

 		  slice_t slice1( edf , signals(i) , interval );
 		  slice_t slice2( edf , signals(j) , interval );
		  
 		  const std::vector<double> * d1 = slice1.pdata();
		  const std::vector<double> * d2 = slice2.pdata();

		  //
		  // Calculate MI for this epoch
		  // but force to whole-trace thresholds
		  //
		  
		  mi_t emi( *d1 , *d2 );		  

		  emi.force_thresholds( tha , thb );

		  emi.calc();

		  // empirical p-values
		  if ( nperms ) 
		    {
		      double pemp, pz;
		      emi.permute( nperms , &pemp , &pz );
		      writer.value( "EMP" , pemp );
		      writer.value( "Z" , pz );
		    }
		  
		  
		  //
		  // Output
		  //
		  
		  writer.value( "MI" , emi.mutinf );
		  writer.value( "JINF" , emi.jointinf );
		  writer.value( "TOTCORR" , emi.total_corr );
		  writer.value( "DTOTCORR" , emi.dual_total_corr );
		  writer.value( "INFA" , emi.infa );
		  writer.value( "INFB" , emi.infb );		  
		  
		  // next epoch 
		}
	      
	      writer.unepoch();
	      
	    }  // end of 'if epoched?'
	
	} // next pair .. 
    } // .. of channels
      

  // 
  // All done 
  //
  
  writer.unlevel( "CHS" );
  
}




void mi_t::calc()
{ 

  // assumes bina and binb are set, and nbins
  
  // univar
  std::vector<double> pa( nbins , 0 );
  std::vector<double> pb( nbins , 0 );

  // joint
  std::vector<std::vector<double> > pab( nbins );
  for (int j=0;j<nbins;j++) pab[j].resize( nbins , 0 );

  // count all
  for (int i=0;i<n;i++) 
    {
      ++pa[ bina[i] ];
      ++pb[ binb[i] ];
      ++pab[ bina[i] ][ binb[i] ];
    }
  
  // init
  infa = 0; infb = 0; jointinf = 0; mutinf = 0;
  //  double mutinf2 = 0;

  // as calcs
  for (int j=0;j<nbins;j++)
    {
      pa[j] /= (double)n;
      pb[j] /= (double)n;
      for (int k=0;k<nbins;k++) pab[j][k] /= (double)n;
    }

  // shannonize, calc
  for (int j=0;j<nbins;j++)
    {
      infa -= pa[j] * log2( pa[j] + eps );
      infb -= pb[j] * log2( pb[j] + eps );
            
      for (int k=0;k<nbins;k++)
	{
	  jointinf -= pab[j][k] * log2( pab[j][k] + eps );	  
	  //mutinf2 += pab[j][k] * log2( pab[j][k] / ( pa[j] * pb[j] + eps ) + eps );
	}
      
    }

  mutinf = infa + infb - jointinf;

  // normalized variants

  // https://en.wikipedia.org/wiki/Mutual_information#Normalized_variants  

  total_corr = mutinf / ( infa < infb ? infa : infb );
  dual_total_corr = mutinf / jointinf;



  
  
  
}


//fd_bins      = ceil(maxmin_range/(2.0*iqr(signal1)*n^(-1/3))); % Freedman-Diaconis 
//scott_bins   = ceil(maxmin_range/(3.5*std(signal1)*n^(-1/3))); % Scott
//sturges_bins = ceil(1+log2(n)); % Sturges


int mi_t::set_nbins_fd()
{
 
  // Freedman-Diaconis rule: 
    
  // get min/max
  double maxa, maxb, mina, minb;
  MiscMath::minmax( da , &mina, &maxa );
  MiscMath::minmax( db , &minb, &maxb );
  double rnga = maxa - mina;
  double rngb = maxb - minb;
  
  // get IQR
  double Qa = MiscMath::iqr( da );
  double Qb = MiscMath::iqr( db );
  
  int nbinsa = ceil( rnga / ( 2 * Qa * pow( n , -1/3.0) ) );
  int nbinsb = ceil( rngb / ( 2 * Qb * pow( n , -1/3.0) ) );
  
  nbins = ceil( nbinsa +  nbinsb / 2.0 );
    
  return nbins;
}


int mi_t::set_nbins_scott()
{
  // Scott's rule
    
  // get min/max
  double maxa, maxb, mina, minb;
  MiscMath::minmax( da , &mina, &maxa );
  MiscMath::minmax( db , &minb, &maxb );
  double rnga = maxa - mina;
  double rngb = maxb - minb;
  
  // get SD
  double Sa = MiscMath::sdev( da );
  double Sb = MiscMath::sdev( db );
  
  int nbinsa = ceil( rnga / ( 3.5 * Sa * pow( n , -1/3.0) ) );
  int nbinsb = ceil( rngb / ( 3.5 * Sb * pow( n , -1/3.0) ) );
  
  nbins = ceil( nbinsa +  nbinsb / 2.0 );

  return nbins;
}

int mi_t::set_nbins_sturges()
{
  nbins = ceil( 1 + log2( n ) );
  return nbins;
}
    

void mi_t::set_nbins(const int b )
{
  nbins = b;
}


int mi_t::set_thresholds( const int bin_rule )
{
  // nbins should be have set

  double maxa, maxb, mina, minb;
  MiscMath::minmax( da , &mina, &maxa );
  MiscMath::minmax( db , &minb, &maxb );
  double rnga = maxa - mina;
  double rngb = maxb - minb;
  double inca = rnga / (double)nbins;
  double incb = rngb / (double)nbins;
  tha.resize( nbins );
  thb.resize( nbins );
  
  for (int i = 0 ; i < nbins  ; i++ ) 
    {
      tha.push_back( mina + i * inca );
      thb.push_back( minb + i * incb );
    }
  
  bin_data();

  return nbins;

}

void mi_t::bin_data()
{

  // and bin data (default is largest bin, i.e. not explicitly tested here, as only t, not t+1 breaks
  bina.resize( n , nbins-1 );
  binb.resize( n , nbins-1 );
  
  for (int i=0;i<n;i++)
    {
      for (int j=1;j<nbins;j++)
	if ( da[i] < tha[j] ) { bina[i] = j-1; break; }       
      for (int j=1;j<nbins;j++)
	if ( db[i] < thb[j] ) { binb[i] = j-1; break; }       
    }  

}


void mi_t::permute( const int nrep , double * pemp , double * pz )
{
  double r = 0;
  std::vector<double> stats;

  for (int p=0;p<nrep;p++)
    {
      //      std::cerr << ".";

      // determine random shift
      int shift = CRandom::rand( n );      
      
      // joint distrib.
      std::vector<std::vector<double> > pab( nbins );
      for (int j=0;j<nbins;j++) pab[j].resize( nbins , 0 );
      
      // count all
      for (int i=0;i<n;i++) 
	{
	  int permi = i + shift;
	  if (  permi >= n ) permi -= n;
	  
	  //	  std::cout << "permi " << i <<" " << permi << " " << n << "\n";
	  ++pab[ bina[i] ][ binb[permi] ];
	}
      
      // init
      double pjointinf = 0;
      
      // shannonize, calc
      for (int j=0;j<nbins;j++)
	for (int k=0;k<nbins;k++)
	  {
	    pab[j][k] /= (double)n;
	    pjointinf -= pab[j][k] * log2( pab[j][k] + eps );
	  }
      
      
      double stat = infa + infb - pjointinf;
      
      if ( stat >= mutinf ) ++r;
      
      stats.push_back( stat );
      
    }
  
  *pemp = ( r+1.0 ) / ( nrep+1.0 ); 

  double null_mean = MiscMath::mean( stats );
  double null_sd   = MiscMath::sdev( stats );
  *pz = ( mutinf - null_mean ) / null_sd ;
  
  return ;  
}

