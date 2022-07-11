
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

#include "dsp/mt_spectrogram.h"

#include "dsp/mtm/mtm.h"

mt_spectrogram_t::mt_spectrogram_t( const Data::Matrix<double> & X ,
				    const int sr ,
				    const double npi , // nw
				    const int nwin ,   // t
				    const double segment_size_sec ,
				    const double segment_step_sec ,
				    const double min_f , 
				    const double max_f , 
				    const bool dB ,
				    const bool mean_center )
{

  // 1) get PSD for each observation & then average
  // 2) get average time-series, & then get PSD

  // nb. X is time x obs 
  const int nobs = X.dim2();
  const int np   = X.dim1();
  
  // Step size in sample-points
  
  const int segment_size = sr * segment_size_sec;
  const int segment_step = sr * segment_step_sec;
      
  // nw / npi
  // t  / nwin
  
  Z.clear(); Z_median.clear(); ZZ.clear();
  int out_f = 0 , out_t = 0;

  std::vector<std::vector<std::vector<double> > > Z_trk; // for medians
  
  for (int i=0; i<nobs; i++)
    {

      // get row
      const std::vector<double> * d = X.col_pointer(i)->data_pointer();
      
      //      std::cout << " row size = " << d->size() << "\n";
      
      // do MTM
      mtm_t mtm( npi , nwin ) ;      
      mtm.dB = dB;
      mtm.opt_remove_mean = mean_center;
      
      //
      // do MTM   (i==0 means only give verbose output on first row)
      //
    
      mtm.apply( d , sr , segment_size , segment_step , i == 0 );
      
      
      //
      // get dimensions, and size Z & ZZ ( first round only )
      //
      
      if ( out_f == 0 )
	{
	  frq.clear();
	  for ( int i = 0 ; i < mtm.f.size() ; i++ )
	    if ( mtm.f[i] >= min_f && mtm.f[i] <= max_f )
	      frq.push_back( mtm.f[i] );	  

	  out_f = frq.size();
	  out_t = mtm.espec.size();
	  
	  Z.resize( out_f , out_t );
	  Z_median.resize( out_f , out_t );
	  ZZ.resize( out_f , out_t );

	  Z_trk.resize( out_f );
	  for (int f=0; f<out_f; f++)
	    Z_trk[f].resize( out_t );
	  
	}
      
      const int nsegs = mtm.espec.size();
      if ( nsegs != out_t )
	Helper::halt( "internal problem in mt_spectrogram_t() " );
      
      //
      // Extract and aggregate PSD
      //
      
      for ( int j = 0 ; j < nsegs ; j++)
	{
	  int fidx = 0;
	  for ( int i = 0 ; i < mtm.f.size() ; i++ )
	    if ( mtm.f[i] >= min_f && mtm.f[i] <= max_f )
	      {
		Z[fidx][j] += mtm.espec[j][fidx] ;
		Z_trk[fidx][j].push_back( mtm.espec[j][fidx] );
		ZZ[fidx][j] += mtm.espec[j][fidx] * mtm.espec[j][fidx] ; 
		++fidx;
	      }
	}
      
      // next observation
    }
  

  //
  // normalize
  //

  for (int i=0; i<out_f; i++)
    for (int j=0; j<out_t; j++)
      {
	double mean = Z[i][j] / (double)nobs;
	double var = ZZ[i][j] / (double)nobs - mean * mean ;
	Z[i][j] = mean; // mean
	Z_median[i][j] = MiscMath::median( Z_trk[i][j] ); // median
	ZZ[i][j] = sqrt( var ); // SD	
      }

  //
  // also set up t[]
  //


  // assume odd number of t-bins -w , MID , +w
  
  double mid_t = (double)(out_t-1)/2.0 * segment_step_sec + segment_size_sec / 2.0 ; 
  
  t.clear();
  for (int i=0; i<out_t; i++)
    t.push_back( (  i * segment_step_sec + segment_size_sec / 2.0 ) - mid_t );
  
  
}
