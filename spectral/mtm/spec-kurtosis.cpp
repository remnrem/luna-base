
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

#include "spectral/mtm/mtm.h"

// take for 1+ channels, power values, for a) segments x b) bins within a band
//  take average over channels
//  return kurtosis per band 

spectral_kurtosis_t::spectral_kurtosis_t( bool k3 )
{
  kurt3 = k3;
  bands.clear();
  bands.insert( SLOW );
  bands.insert( DELTA );
  bands.insert( THETA );
  bands.insert( ALPHA );
  bands.insert( SIGMA );
  bands.insert( BETA );
  bands.insert( GAMMA );  
} 


void spectral_kurtosis_t::setf( const std::vector<double> & f_ )
{
  f = f_;
}

void spectral_kurtosis_t::add( int ch , const std::vector<std::vector<double> > & x )
{
  // x rows = seg x bins

  if ( x.size() == 0 ) return; // empty

  if ( x[0].size() != f.size() )
    Helper::halt( "internal error in MTM/speckurt" );

  ch2segxf[ ch ] = x;
    
}

void spectral_kurtosis_t::average_channels()
{
  int ns = -1;
  if ( ch2segxf.size() == 0 ) return;

  std::map<int,std::vector<std::vector<double> > >::const_iterator ss = ch2segxf.begin();
  while ( ss != ch2segxf.end() )
    {
      if ( ns == -1 ) ns = ss->second.size();
      else if ( ss->second.size() != ns )
	Helper::halt( "internal error in speckurt" );
      ++ss;
    }

  int nf = f.size();
  segxf.resize( ns );
  for (int i=0;i<ns;i++) segxf[i].resize( nf , 0 );

  ss = ch2segxf.begin();
  while ( ss != ch2segxf.end() )
    {
      for (int i=0;i<ns;i++)
	for (int j=0;j<nf;j++)
	  segxf[i][j] += ss->second[i][j];
      ++ss;
    }

  int nch = ch2segxf.size();
  for (int i=0;i<ns;i++)
    for (int j=0;j<nf;j++)
      segxf[i][j] /= (double)nch;
  
}


double spectral_kurtosis_t::kurtosis( frequency_band_t b , double * sd , double * skew )
{
  freq_range_t band = globals::freq_band[ b ];
  std::vector<double> xx;
  int nf=f.size();
  int n_segs=segxf.size();
  
  for (int fi=0;fi<nf;fi++)
    {
      if ( f[fi] >= band.first && f[fi] < band.second )
	{
	  for (int i=0; i<n_segs; i++)
	    xx.push_back( segxf[i][fi] );
	}
    }
  
  if ( xx.size() < 2 )
    return -999;

  double k = MiscMath::kurtosis( xx ) + ( kurt3 ? 3 : 0 ) ;

  if ( sd != NULL )
    {
      std::vector<double> lntracker( xx.size() );
      for (int i=0; i<xx.size(); i++)
	lntracker[i] = log( xx[i] );
      
      // sd of natural log scaled                                                          
      const double xsd = MiscMath::sdev( lntracker );
      
      // CV, using formula for log-normal data                                             
      *sd = sqrt( exp( xsd * xsd ) -1 );

    }

  if ( skew != NULL )
    *skew = MiscMath::skewness( xx );

  //  std::cout << " len1 = " << xx.size() << "\n";

  return k;
  
}


double spectral_kurtosis_t::kurtosis2( frequency_band_t b , double * sd, double * skew )
{
  // alternate version where we sum in the freq domain first
  // and then kurtosis just for the e.g. 29 values
  
  freq_range_t band = globals::freq_band[ b ];
  
  int nf=f.size();
  int n_segs=segxf.size();
  
  std::vector<double> xx;
  for (int i=0; i<n_segs; i++)
    {
      double b = 0;
      for (int fi=0;fi<nf;fi++)
	if ( f[fi] >= band.first && f[fi] < band.second )
	  b += segxf[i][fi];
      xx.push_back( b );
    }
  
  if ( xx.size() < 2 )
    return -999;
  
  double k = MiscMath::kurtosis( xx ) + ( kurt3 ? 3 : 0 ) ;
  
  if ( sd != NULL )
    {
      std::vector<double> lntracker( xx.size() );
      for (int i=0; i<xx.size(); i++)
	lntracker[i] = log( xx[i] );
      
      // sd of natural log scaled                                                          
      const double xsd = MiscMath::sdev( lntracker );
      
      // CV, using formula for log-normal data                                             
      *sd = sqrt( exp( xsd * xsd ) -1 );
    }

  //  std::cout << " len2 = " << xx.size() << "\n";

  if ( skew != NULL ) 
    *skew = MiscMath::skewness( xx );

  return k;
  
}
