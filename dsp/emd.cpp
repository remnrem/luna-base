
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


// EMD based on methods implemented in R 'EMD' package:
// Notes: https://journal.r-project.org/archive/2009-1/RJournal_2009-1_Kim+Oh.pdf

// TODO:  make extrema robust to clipped/repeated values: get interval that spans zero-cross / max/min
//        different stopping rules
//        boundary effects for IMF

// currently set to : boundary condition 'wave'
//                    stopping rule =    'type1'

// https://fr.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/55938/versions/6/previews/Clustering_toolbox/utils/emd.m/index.html?access_key=

#include "emd.h"
#include "spline.h"

#include "miscmath/miscmath.h"
#include "helper/helper.h"

#include "edf/edf.h"
#include "edf/slice.h"

#include "helper/logger.h"

#include <iostream>
#include <set>

extern logger_t logger;

void dsptools::emd_wrapper( edf_t & edf , param_t & param ) 
{
  // add IMF by epochs

  std::string signal_label = param.requires( "sig" );

  const bool no_annotations = true;

  signal_list_t signals = edf.header.signal_list( signal_label , no_annotations );

  const int ns = signals.size();

  // if SIG --> SIG_IMF_1 , SIG_IMF_2, ...
  std::string component_tag = param.has( "tag" ) ? param.value( "tag" ) : "_IMF_";


  const int max_sift = param.has( "sift" ) ? param.requires_int( "sift" ) : 20 ;
  const int max_imf = param.has( "imf" ) ? param.requires_int( "imf" ) : 10 ;
  
  
  //
  // iterate over each signal
  //

  for (int s=0; s<ns; s++)
    {

      //
      // iterate over epochs
      //
      
      // int ne = edf.timeline.first_epoch();

      // const bool by_epoch = false;
      
      // while ( 1 )
      // 	{
	  
      // 	  int epoch = edf.timeline.next_epoch();

      // 	  if ( epoch == -1 ) break;

	  //interval_t interval = by_epoch ? edf.timeline.epoch( epoch ) : 

      interval_t interval = edf.timeline.wholetrace();
	  
      slice_t slice( edf , signals(s) , interval );
      
      const std::vector<double> * d = slice.pdata();
      
      emd_t emd;

      emd.max_sift = max_sift;
      emd.max_imf = max_imf;
      
      logger << "  processing " << signals.label(s) << "... ";

      const int nimf = emd.proc( d );	  
      
      logger << "  adding " << nimf << " IMFs\n";

      //
      // add signals
      //

      for (int c=0;c<nimf;c++)
        {	  
	  std::string imflab = signals.label(s) + component_tag + Helper::int2str( c+1 );
          edf.add_signal( imflab , edf.header.sampling_freq( signals(s) ) , emd.imf[c] );
        }

      // residual 'IMF0'
      std::string imflab = signals.label(s) + component_tag + "0";
      edf.add_signal( imflab , edf.header.sampling_freq( signals(s) ) , emd.residual );
      
    } // next signal

    
}


// return unique list of points 
std::vector<int> extrema_t::maxindex()
{
  std::set<int> s;
  for (int i=0;i<maxindex_start.size();i++)
    {
      s.insert( maxindex_start[i] );
      s.insert( maxindex_stop[i] );
    }
  std::vector<int> res;
  std::set<int>::iterator ss = s.begin();
  while ( ss != s.end() ) 
    {
      res.push_back( *ss );
      ++ss;
    }
  return res;
}

std::vector<int> extrema_t::minindex()
{
  std::set<int> s;
  for (int i=0;i<minindex_start.size();i++)
    {
      s.insert( minindex_start[i] );
      s.insert( minindex_stop[i] );
    }
  std::vector<int> res;
  std::set<int>::iterator ss = s.begin();
  while ( ss != s.end() ) 
    {
      res.push_back( *ss );
      ++ss;
    }
  return res;
}


extrema_t::extrema_t( const std::vector<double> & x )
{
  
  ncross = 0;
  nmin = 0;
  nmax = 0;
  nextrema = 0;
  
  // n-1 vector of differences:  -1, 0, +1
  const int n = x.size();
  std::vector<int> z1;
  std::vector<int> index1;

  bool seen_plus = false, seen_neg = false;

  int last_diff = 0;
  for (int i=0;i<n-1;i++)
    {
      int diff = 0;

      if ( x[i] != x[i+1] ) 
	{
	  // sign is w.r.t. i, i.e. if i is minimum, then set to -1
	  diff = x[i+1] - x[i] > 0 ? - 1 : +1 ; 
	  
	  if ( last_diff != 0 && last_diff != diff ) 
	    {
	      int ii = i;
	      while ( 1 ) 
		{
		  --ii;
		  if ( ii < 0 ) break;
		  if ( x[ii] != x[i] ) { ++ii; break; } 
		}
	      index1.push_back(ii);
	      z1.push_back(diff);	  
	      if ( diff > 0 ) seen_plus = true; else seen_neg = true;
	    }
	  
	  last_diff = diff; 
	}
    }
  
  nextrema = z1.size();
  
  // no local extrema, return
  if ( index1.size() == 0 || (!seen_plus) || (!seen_neg) ) 
    return;

  minindex_start.clear(); minindex_stop.clear();
  maxindex_start.clear(); maxindex_stop.clear();
  
  if ( nextrema >= 2 )
    {
      for (int j=0;j<nextrema-1;j++)
	{
	  int i1 = index1[j];
	  int tmpi2 = index1[j+1]-1;
	  int i2 = i1;
	  
	  for (int k=i1;k<=tmpi2;k++)
	    if ( x[ k ] == x[ i1 ] ) i2 = k;

	  if ( z1[j] > 0  ) 
	    {
	      maxindex_start.push_back( i1 );
	      maxindex_stop.push_back( i2 );
	    }
	  else
	    {
	      minindex_start.push_back( i1 );
	      minindex_stop.push_back( i2 );
	    }

//  	  std::cout << j+1 << "\t" << i1+1 << " " << i2+1 << " " 
//  		    << (  i1!=i2 ? "ZZ" : "" ) << " "  
//  		    << ( z1[j] > 0 ? "MAX" : "MIN" ) << "\n";
	}     
    }


  // final one
  int i1 = index1[ nextrema - 1 ];
  int tmpi2 = n - 1 - 1;
  int i2 = i1;

  for (int k=i1;k<=tmpi2;k++)
    if ( x[ k ] == x[ i1 ] ) i2 = k;
  
  if ( z1[nextrema-1] > 0  ) 
    {
      maxindex_start.push_back( i1 );
      maxindex_stop.push_back( i2 );
    }
  else
    {
      minindex_start.push_back( i1 );
      minindex_stop.push_back( i2 );
    }
  
//   std::cout << nextrema   << "\t" << i1+1 << " " << i2+1 << " " 
// 	    << (  i1!=i2 ? "ZZ" : "" ) << " "  
// 	    << ( z1[nextrema-1] > 0 ? "MAX" : "MIN" ) << "\n";
  

  
  //
  // Zero-crossing count
  //
  
  // '-1' counts as the first sample point
  for (int i=-1;i<nextrema-1;i++)
    {

      // i.e. implicitly insert first point in

      int ii = i == -1 ? 0 : index1[i];
      int jj = index1[i+1];

      // point 1 : ii
      // point 2 : jj

      // does this point itself equal 0?
      if ( x[ii] == 0 ) 
	{
	  int i1 = ii;
	  int i2 = i1;
	  for (int k=i1;k<=jj;k++)
	    if ( x[k] == 0 ) i2 = k;
	  cross_start.push_back(i1);
	  cross_stop.push_back(i2);
	}
      else if ( x[ii] * x[jj] < 0 )
	{
	  int tmpmin = -1;
	  for (int k=ii;k<=jj;k++) { if ( x[ii] * x[k] <= 0 ) { tmpmin = k; break; } }
	  if ( tmpmin == -1 ) Helper::halt("internal error in extrema_t()");

	  if ( x[tmpmin] == 0 )
	    {
	      // search from tmpmin to jj
	      int tmin = -1, tmax = -1;
	      for (int k=tmpmin;k<=jj;k++)
		{
		  if ( x[k] == 0 ) 
		    {
		      if ( tmin == -1 ) tmin = k;
		      tmax = k;
		    }
		}
	      cross_start.push_back(tmin);
	      cross_stop.push_back(tmax);
	    }
	  else
	    {
	      cross_start.push_back( tmpmin - 1 );
	      cross_stop.push_back( tmpmin );
	    }
	}
    } // next extrema

  //
  // end case for zero-crossing
  //

  i1 = index1[ nextrema-1 ];
  int tmpmin = -1;

  for (int k=i1; k < n ; k++ )
    if ( x[i1] * x[k] <= 0 ) { tmpmin = k; break; }
	
  if ( tmpmin != -1 )
    {
      if ( x[tmpmin] == 0 ) 
	{
	  int tmin = -1, tmax = -1;
	  for (int k=tmpmin;k<n;k++)
	    {
	      if ( x[k] == 0 )
		{
		  if ( tmin == -1 ) tmin = k;
		  tmax = k;
		}
	    }
	  cross_start.push_back(tmin);
	  cross_stop.push_back(tmax);
	}
      else
	{
	  cross_start.push_back( tmpmin - 1 );
	  cross_stop.push_back( tmpmin );
	}
      
    }
  

  //
  // All done, summary
  //
  
  ncross = cross_start.size();
  nmin = minindex_start.size();
  nmax = maxindex_start.size();
  nextrema = nmin + nmax;

  //  for (int i=0;i<cross_start.size();i++) std::cerr << "cross " << i+1 << " " << cross_start[i]+1 << " " << cross_stop[i]+1 << "\n";
//    for (int i=0;i<minindex.size();i++) std::cerr << "min " << i+1 << " " << minindex[i]+1 << "\n";
//    for (int i=0;i<maxindex.size();i++) std::cerr << "max " << i+1 << " " << maxindex[i]+1 << "\n";
//      std::cerr << "ncross, nextrema = " << ncross << " " << nextrema << "\n";
//      std::cerr << "nmin, nmax = " << nmin << " " << nmax << "\n";
//  std::cerr << "nextrema, zc " << nextrema << " " << ncross << " " << ( nextrema == ncross || nextrema == ncross + 1 ? "Y" : "." ) << "\n";

}

std::vector<double> emd_t::sift( const std::vector<double> & x )
{

  //
  // Extract an IMF from 'x'
  //
  
  const int n = x.size();

  std::vector<double> h = x;

  //
  // Begin sifting
  //

  int j = 1;

  while ( 1 ) 
    {

      if ( verbose ) 
	std::cerr << " sifting " << j << "\n";;
      
      
      //
      // mean of envelope
      //

      std::vector<double> m = envelope_mean( h );
      

      //
      // require at least 'x' extrema: if this wasn't met, we should have an empty matrix here
      //
      
      if ( m.size() == 0 ) 
	break;
    
      stop_mode = 1;

      std::vector<double> h1 = h;

      for (int i=0;i<n;i++)
	h1[i] -= m[i];

      // return as IMF if more than max sifts

      if ( j >= max_sift ) 
	{ 

	  if ( verbose ) 
	    std::cerr << "required " << j << " sifting iterations\n"; 
	  
	  extrema_t fex( h1 );

	  if ( verbose ) 
	    std::cerr << "H1 nextrema, zc " << fex.nextrema << " " << fex.ncross << " "
		      << ( fex.nextrema == fex.ncross || fex.nextrema == fex.ncross + 1 ? "Y" : "." ) << "\n";

	  return h1; 
	} 
      
      
      // otherwise, consider other stopping rules
      if ( stop_mode == 1 ) 
	{

	  int cnt = 0;
	  double mx = 0;	      
	  for (int i=0;i<n;i++)
	    if ( fabs(m[i]) > mx )
	      {
		mx = fabs(m[i]);
		++cnt;
	      }
	  
	  if ( verbose ) 
	    std::cerr << "tol,mx,cnt  << " << tol << " " << mx << " " << cnt << "\n";
	  
	  if ( mx < tol )
	    {

	      if ( verbose ) 
		std::cerr << "required " << j << " sifting iterations\n"; 
	      
	      extrema_t fex( h1 );

	      if ( verbose )
		std::cerr << "H1 nextrema, zc " << fex.nextrema << " " << fex.ncross
			  << " " << ( fex.nextrema == fex.ncross || fex.nextrema == fex.ncross + 1 ? "Y" : "." ) << "\n";

	      return h1;
	    }
	}
      else if ( stop_mode == 2 && j >= 2 )
	{
	  
	  double sd = 0;
	  for (int i=1;i<n-1;i++) 
	    {
	      sd +=  ( ( h[i] - h1[i] ) * ( h[i] - h1[i] )  ) / (h[i]*h[i]);

	      if ( verbose ) 
		std::cerr << "cum sd " << sd << " " << h[i] << " " << h1[i] << "\n";
	    }

	  if ( verbose ) 
	    std::cerr << "SD= " << sd << " j=" << j << "\n";
      
	  // stop?
	  if ( sd < 0.3 ) 
	    {
	      
	      if ( verbose ) 
		std::cerr << "required " << j << " sifting iterations\n"; 
	      
	      //
	      // check extrema and ZC numbers
	      //
	      
	      extrema_t fex( h1 );

	      std::cerr << "H1 nextrema, zc " << fex.nextrema << " " << fex.ncross << " "
			<< ( fex.nextrema == fex.ncross || fex.nextrema == fex.ncross + 1 ? "Y" : "." ) << "\n";
	      
	      
	      return h1;
	    }

	}
	  
      //
      // continue sifting, store previous 
      //

      if ( verbose ) 
	std::cerr << " going to continue sifting... back for next j\n";

      // set to re-sift results of previous sift

      h = h1;

      ++j;

    }
  
  std::vector<double> dummy;
  return dummy;
  
}

emd_t::emd_t( const bool v )
  : verbose( v )
{
  
  // defaults
  max_sift = 20;
  max_imf  = 10;
  
}

int emd_t::proc( const std::vector<double> * d )
{
  
  std::vector<double> working = *d;

  // default for tolerance, and stop mode (==1) (from EMD R package)  
  tol = MiscMath::sdev( working ) * 0.1*0.1;

  stop_mode = 1;

  //  std::cout << " tol = " << tol << "\n";
  
  const  int n = d->size();
  
  imf.clear();

  // keep track of IMF
  int k = 0;
  
  while ( 1 ) 
    {

      // Get each IMF

      std::vector<double> h = sift( working );

      // not enought extrema on signal/residual: done
      if ( h.size() == 0 )
	break;	
      
      // Store
      imf.push_back( h );
      
      // make residual 
      for (int i=0;i<n;i++)
	working[i] -= h[i];
      
      // next IMF
      ++k;
      
      if ( k > max_imf ) break;
    }

  if ( verbose )
    logger << "  extracted " << k << " IMF\n";

  //
  // final residual
  //
  
  residual = *d;

  for (int i=0;i<n;i++)
    {
      for (int j=0;j<k;j++)
	residual[i] -= imf[j][i];
      
      // if ( verbose ) 
      // 	{
      // 	  std::cout << "EMD: " << i << "\t" << d[i] ;
      // 	  for (int j=0;j<k;j++) std::cout << "\t" << imf[j][i];
      // 	  std::cout << "\t" << residual[i] << "\n";
      // 	}

    }
  
  return imf.size();
}

void emd_t::hht( double Fs )
{

  const int k = imf.size();

  if ( k == 0 ) return;
    
  //
  // Hilbert transform 
  //
  
  // now we have IMFs and residual: apply HHT

  for (int j=0;j<k;j++)   
    {

      hilbert_t hilbert( imf[j] );
      
      std::vector<double> f = hilbert.instantaneous_frequency( Fs );

      for (int k=0;k<f.size();k++)
	std::cout << "IMF " << j << " " << k << " " << f[k] << "\n"; 
      
    }
  
}


std::vector<double> emd_t::envelope_mean( const std::vector<double> & x )
{

  //
  // get extrema
  //
  
  extrema_t extrema( x );

  //
  // check # of extrema (requires at least 2)  
  //
  
  if ( extrema.nextrema <= 2 ) { std::vector<double> dummy; return dummy; } 

  std::vector<int> minindex = extrema.minindex();
  std::vector<int> maxindex = extrema.maxindex();

  if ( verbose ) 
    std::cerr << "n min/max = " << extrema.nextrema << " " << minindex.size() << " " << maxindex.size() << "\n";
  
  //
  // handle boundary ('wave' condition)
  //
  
  int first_pt_idx  = 0;
  int first_min_idx = minindex[0] ;
  int first_max_idx = maxindex[0] ;

  double first_pt = x[ first_pt_idx ];
  double first_min = x[ first_min_idx ];
  double first_max = x[ first_max_idx ];
  
  double d1 = first_max_idx < first_min_idx ? first_max_idx : first_min_idx ;  
  double d2 = first_max_idx < first_min_idx ? first_min_idx - first_max_idx : first_max_idx - first_min_idx ;  

  double wavefreq1 = 0;
  bool add_first_min = false , add_first_max = false , add_last_min = false , add_last_max = false;

  if ( verbose ) 
    std::cerr << " first_min, first_max = " << first_min << " " << first_max << "  " << first_pt << "\n";

  if ( first_pt <= first_min && first_pt <= first_max ) 
    {
      add_first_min = true;
      wavefreq1 = 2 * d1;
    }
  else if ( first_pt >= first_min && first_pt >= first_max )
    {
      add_first_max = true;
      wavefreq1= 2 * d1;
    }
  else if ( first_pt >= (double)( first_min + first_max ) / 2.0 )
    {
      if ( d2 > 2 * d1 ) 
	wavefreq1 = 2 * d2 ;
      else
	wavefreq1 = d2 + 2 * d1 ;
    }
  else
    {
      if ( d2 > round( 1.5 * d1 ) )
	wavefreq1 = 2 * d2;
      else
	wavefreq1 = d2 + round(1.5 * d1 );
    }


  int last_pt_idx  = x.size()-1;
  int last_min_idx = minindex[ minindex.size()-1] ;
  int last_max_idx = maxindex[ maxindex.size()-1] ;

  double last_pt = x[ last_pt_idx ];
  double last_min = x[ last_min_idx ];
  double last_max = x[ last_max_idx ];
  
  d1 = last_max_idx < last_min_idx ? last_min_idx - last_max_idx : last_max_idx - last_min_idx ;  
  d2 = last_max_idx < last_min_idx ? last_pt_idx - last_min_idx : last_pt_idx - last_max_idx ;  
  
  double wavefreq2 = 0;

  if ( last_pt <= last_min && last_pt <= last_max ) 
    {
      add_last_min = true;
      wavefreq2 = 2 * d2;
    }
  else if ( last_pt >= last_min && last_pt >= last_max ) 
    {
      add_last_max = true;
      wavefreq2 = 2 * d2;
    }
  else if ( last_pt >= ( last_min + last_max ) / 2.0 ) 
    {

      if ( d1 > 2 * d2 ) wavefreq2 = 2 * d1 ;
      else wavefreq2 = d1 + 2 * d2 ;
    }
  else
    {
      if ( d1 > round( 1.5 * d2 ) ) wavefreq2 = 2 * d1;
      else wavefreq2 = d1 + round(1.5 * d2 );
    }

  if ( verbose ) 
    std::cerr << "wavefreqs " << wavefreq1 <<" " << wavefreq2 << "\n";

 
  //
  // Set extrema and values for cubic spline
  //

  std::vector<double> e_min_idx, e_max_idx;
  std::vector<double> e_min_val, e_max_val;

  
  // do we need to add new extrema?
  // i.e. if first/last point has become a new local min/max extrema

  if ( add_first_min )
    {
      first_min = first_pt;
      first_min_idx = first_pt_idx;
    }

  if ( add_first_max ) 
    {
      first_max = first_pt;
      first_max_idx = first_pt_idx;
    }

  if ( add_last_min )
    {
      last_min = last_pt;
      last_min_idx = last_pt_idx;
    }

  if ( add_last_max ) 
    {
      last_max = last_pt;
      last_max_idx = last_pt_idx;
    }
  

  //
  // left boundary
  //

  for (int i=4;i>=1;i--)
    {      
      e_min_idx.push_back( first_min_idx - i *  wavefreq1 );
      e_min_val.push_back( first_min );
      e_max_idx.push_back( first_max_idx - i *  wavefreq1 );
      e_max_val.push_back( first_max );
    }

  //
  // possible new extrema
  //

  if ( add_first_min ) 
    {
      e_min_idx.push_back( 0 );
      e_min_val.push_back( first_min );
    }
  
  if ( add_first_max ) 
    {
      e_max_idx.push_back( 0 );
      e_max_val.push_back( first_max );
    }

  //
  // inner (main signal)
  //
  
  for (int i=0;i<extrema.nmin;i++ ) 
    {
      e_min_val.push_back( x[ minindex[ i ] ] );
      e_min_idx.push_back( minindex[ i ] );
    }

  for (int i=0;i<extrema.nmax;i++ ) 
    {
      e_max_val.push_back( x[ maxindex[ i ] ] );
      e_max_idx.push_back( maxindex[ i ] );
    }
  

  //
  // possible new extrema for final point
  //
  
  if ( add_last_min ) 
    {
      e_min_idx.push_back( x.size() );
      e_min_val.push_back( last_min );
    }

  if ( add_last_max ) 
    {
      e_max_idx.push_back( x.size() );
      e_max_val.push_back( last_max );
    }
  
  //
  // right boundary
  //

  for (int i=1;i<=4;i++)
    {
      e_min_idx.push_back( last_min_idx + i *  wavefreq2 );
      e_min_val.push_back( last_min );
      e_max_idx.push_back( last_max_idx + i *  wavefreq2 );
      e_max_val.push_back( last_max );
    }
  

  //
  // get upper spline
  // nb. requires the _idx is pre-sorted
  //

  // std::cout << " max = " << e_max_idx.size() <<" " << e_max_val.size() << "\n";
  // for (int ii=0; ii<e_max_idx.size() ; ii++)
  //   std::cout << " ii= " << e_max_idx[ii] << " " << e_max_val[ii] << "\n";
  
  tk::spline sa;
  sa.set_points( e_max_idx , e_max_val ); 
  
  tk::spline sb;
  sb.set_points( e_min_idx , e_min_val ); 

  
  //
  // get mean envelope 
  //

  const int n = x.size();
  
  std::vector<double> env( n );
  for (int i=0; i<n; i++)
    env[i] = ( sa( i ) + sb( i ) ) / 2.0 ; 
  
  return env;

}

