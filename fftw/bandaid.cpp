
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

#include "fftw/bandaid.h"


bandaid_t::bandaid_t()
{
  init();
  
}

void bandaid_t::init()
{

  // nb. this function clears the current tracker, but any changed freq band defs remain

  track_band.clear();
  
  slow = delta = theta = alpha = sigma = beta = gamma = 0;
  low_sigma = high_sigma = 0;
  denom = total = 0;
  
  bands.clear();
  bands.push_back( SLOW );
  bands.push_back( DELTA );
  bands.push_back( THETA );
  bands.push_back( ALPHA );
  bands.push_back( SIGMA );
  //     bands.push_back( LOW_SIGMA );
  //     bands.push_back( HIGH_SIGMA );
  bands.push_back( BETA );
  bands.push_back( GAMMA );
  bands.push_back( DENOM );

}

int bandaid_t::size() const
{
  return bands.size();
}

void bandaid_t::define_bands( param_t & param )
{
  
  if ( param.has( "slow" ) )
    {
      double f0, f1;
      freq_band_settings( param.value( "slow" ) , &f0, &f1 );
      globals::freq_band[ SLOW ] =  freq_range_t( f0 , f1 ) ;
      logger << "  defining slow as " << f0 << " to " << f1 << " Hz\n";
    }
  
  if ( param.has( "delta" ) )
    {
      double f0, f1;
      freq_band_settings( param.value( "delta" ) , &f0, &f1 );
      globals::freq_band[ DELTA ] =  freq_range_t( f0 , f1 ) ;
      logger << "  defining delta as " << f0 << " to " << f1 << " Hz\n";
    }

  if ( param.has( "theta" ) )
    {
      double f0, f1;
      freq_band_settings( param.value( "theta" ) , &f0, &f1 );
      globals::freq_band[ THETA ] =  freq_range_t( f0 , f1 ) ;
      logger << "  defining theta as " << f0 << " to " << f1 << " Hz\n";
    }
  
  if ( param.has( "alpha" ) )
    {
      double f0, f1;
      freq_band_settings( param.value( "alpha" ) , &f0, &f1 );
      globals::freq_band[ ALPHA ] =  freq_range_t( f0 , f1 ) ;
      logger << "  defining alpha as " << f0 << " to " << f1 << " Hz\n";
    }

  if ( param.has( "sigma" ) )
    {
      double f0, f1;
      freq_band_settings( param.value( "sigma" ) , &f0, &f1 );
      globals::freq_band[ SIGMA ] =  freq_range_t( f0 , f1 ) ;
      logger << "  defining sigma as " << f0 << " to " << f1 << " Hz\n";
    }

  if ( param.has( "beta" ) )
    {
      double f0, f1;
      freq_band_settings( param.value( "beta" ) , &f0, &f1 );
      globals::freq_band[ BETA ] =  freq_range_t( f0 , f1 ) ;
      logger << "  defining beta as " << f0 << " to " << f1 << " Hz\n";
    }

  if ( param.has( "gamma" ) )
    {
      double f0, f1;
      freq_band_settings( param.value( "gamma" ) , &f0, &f1 );
      globals::freq_band[ GAMMA ] =  freq_range_t( f0 , f1 ) ;
      logger << "  defining gamma as " << f0 << " to " << f1 << " Hz\n"; 
    }
  

  //
  // User defined 'TOTAL' --> DENOM? 
  //
  
  globals::freq_band[ DENOM ] = globals::freq_band[ TOTAL ];

  if ( param.has( "total" ) )
    {
      double f0, f1;
      freq_band_settings( param.value( "total" ) , &f0, &f1 );
      logger << "  setting total power (denominator for RELPSD) to " << f0 << " to " << f1 << " Hz\n"; 
      globals::freq_band[ DENOM ] = freq_range_t( f0 , f1 ) ;
    }  
}




void bandaid_t::freq_band_settings( const std::string & b , double * r0 , double * r1 )
{  
  std::vector<std::string> f = Helper::parse( b , ",-" );
  if ( f.size() != 2 ) Helper::halt( "expecting band=lower-upper" );
  double f0, f1;
  if ( ! Helper::str2dbl( f[0] , &f0 ) ) Helper::halt( "expecting numeric for power range" );
  if ( ! Helper::str2dbl( f[1] , &f1 ) ) Helper::halt( "expecting numeric for power range" );
  if ( f0 >= f1 ) Helper::halt( "expecting band=lower-upper" );
  if ( f0 < 0 || f1 < 0 ) Helper::halt( "negative frequencies specified" );
  *r0 = f0;
  *r1 = f1;
}

void bandaid_t::track()
{
  // used w/ MTM (or others) -- assumes 'current' values have just been added
  // (via calc_bands() from a spectrum.. store these)
  
  track_band[ SLOW  ].push_back( slow );
  track_band[ DELTA ].push_back( delta );
  track_band[ THETA ].push_back( theta );
  track_band[ ALPHA ].push_back( alpha );
  track_band[ SIGMA ].push_back( sigma );
  track_band[ BETA  ].push_back( beta );
  track_band[ GAMMA ].push_back( gamma );
  track_band[ DENOM ].push_back( denom );
  
  // track_band[ LOW_SIGMA ].push_back( low_sigma );
  // track_band[ HIGH_SIGMA ].push_back( high_sigma );
    
}

void bandaid_t::track_bands_per_epoch( double slow_ , double delta_, double theta_, double alpha_,
				       double sigma_, double low_sigma_ , double high_sigma_, double beta_,
				       double gamma_, double denom_ )
{

  slow = slow_;
  delta = delta_;
  theta = theta_;
  alpha = alpha_;
  sigma = sigma_;
  low_sigma = low_sigma_;
  high_sigma = high_sigma_;
  beta = beta_;
  gamma = gamma_;
  denom = denom_;

  // total power (may differ from denom)
  total = slow + delta + theta + alpha + sigma + beta + gamma ; 

  //
  // track epoch-level band-power statistics
  //
  
  track();
  
}


double bandaid_t::psdsum( const std::vector<double> & f , const std::vector<double> & x, const freq_range_t & b )
{

  const int N = f.size();
  
  // add is l <= x < y 
  double r = 0;
  for (int i=0;i<N;i++) 
    {
      if ( f[i] >= b.second ) break;
      if ( f[i] >= b.first ) r += x[i];
    }
  
  // area under the curve, scale by bin width 
  double fbin = N > 1 ? f[1] - f[0] : 1 ;
  
  return r * fbin;
  
}

void bandaid_t::calc_bandpower( const std::vector<double> & f , const std::vector<double> & x )
{
  slow = psdsum( f , x, globals::freq_band[ SLOW ] );
  delta = psdsum( f , x, globals::freq_band[ DELTA ] );
  theta = psdsum( f , x, globals::freq_band[ THETA ] );
  alpha = psdsum( f , x, globals::freq_band[ ALPHA ] );
  sigma = psdsum( f , x, globals::freq_band[ SIGMA ] );
  beta = psdsum( f , x, globals::freq_band[ BETA ] );
  gamma = psdsum( f , x, globals::freq_band[ GAMMA ] );
  denom = psdsum( f , x, globals::freq_band[ DENOM ] );
  total = psdsum( f , x, globals::freq_band[ TOTAL ] );
  
}

double bandaid_t::fetch( frequency_band_t b ) const
{
  switch( b ) {
  case SLOW : return slow;
  case DELTA : return delta;
  case THETA : return theta;
  case ALPHA : return alpha;
  case SIGMA : return sigma;
  case BETA : return beta;
  case GAMMA : return gamma;
  case DENOM : return denom;
  case TOTAL : return total;
  case LOW_SIGMA : return low_sigma ;
  case HIGH_SIGMA : return high_sigma;
  default : return 0;
  }
  
}

