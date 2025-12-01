
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

#include "lunapi/segsrv.h"
#include "lunapi/lunapi.h"

#include <algorithm>
#include <limits>
#include <cassert>
#include <vector>
#include <cmath>
#include <random>

segsrv_t::segsrv_t( lunapi_inst_ptr inst ) : p( inst ) , sigmod( this )
{
  awin = bwin = 0;
  
  max_samples_in = 200; // default input throttling
  max_samples_out = 0;  // no output throttling by default

  epoch_sec = 30;

  xpixels = 800; // default, but should always be changed by caller
  
}

void segsrv_t::init()
{
  
  awin = bwin = 0;
  aidx.clear();
  bidx.clear();
  tmap.clear();
  evts.clear();
  segments = p->edf.timeline.segments();
  gaps = p->edf.timeline.gaps( segments ); 

  sigmod.clear();
  
  annot_format6 = true; // plotly
  clip_xaxes = true;
  
  clocktime_t etime( p->edf.header.starttime );
  if ( ! etime.valid )
    {
      logger << "*** invalid EDF start time - setting to 00.00.00***\n";
      clocktime_t mtime( "00.00.00" );
      edf_start = mtime;
    }
  else
    edf_start = etime;
  
  valid_window = false; // this is set by set_window()
  
  // epoch_sec is 30 by default, but could have been changed (min 4s)
  // prior to init() being called
  
  // define epochs based on simple clock time, inc. by epoch_sec, but
  // then flag the epochs that map (fully) into selected signal space
  //  -->  std::map<int,int> clock2signal_emap;
    
  epoch_num = 0;
  clock_epoch_num = 0;
  epoch_sec_starts.clear();
  clk2sig_emap.clear();


  //  double cumul_sec = 0;
  double last_sec = 0;
  std::set<interval_t>::const_iterator ss = segments.begin();
  while ( ss != segments.end() )
    {
      //      cumul_sec += ss->duration_sec();
      last_sec = ss->stop_sec();
      ++ss;
    }
  
  epoch_sec = 30;
  
  uint64_t etp = 0LLU;
  uint64_t epoch_tp = epoch_sec * globals::tp_1sec;    
  uint64_t max_tp = get_total_sec_original() * globals::tp_1sec;
  double esec = 0; // clock time
  double cumul_esec = 0; // actual signal time (used to make band epoch summs)
  uint64_t last_epoch_tp = 0LLU;
  
  int clk_idx = 0;
  int sig_idx = 0;
  
  while ( 1 )
    {
      // allow last partial epoch (for clock time)
      uint64_t e2tp = etp + epoch_tp;
      
      // past end?
      if ( e2tp > max_tp ) break;
            
      // is this a valid & full epoch? (--> sets valid_window F if all missing)
      set_window(  esec , esec + epoch_sec );
      
      uint64_t good_tps = 0LLU;
      const bool gapped = has_gaps( etp , e2tp , &good_tps );
            
      const bool okay = valid_window &&  ( ! gapped ) && ( esec < last_sec );; 
      
      // will be calculating epoch-level stats for this epoch
      if ( okay )
	{
	  clk2sig_emap[ clk_idx ] = sig_idx;
	  epoch_sec_starts.push_back( cumul_esec );
	  
	  // counts as a signal epoch
	  ++sig_idx;
	}
      
      // advance 
      ++clk_idx;
      ++clock_epoch_num;
      esec += epoch_sec;
      etp += epoch_tp;      
      cumul_esec += good_tps * globals::tp_duration;

    }
  
  // define epochs based on simple clock time, inc. by epoch_sec, but
  // then flag the epochs that map (fully) into selected signal space
  
  epoch_num = epoch_sec_starts.size();
  
  // reset
  awin = bwin = 0;

}


void segsrv_t::calc_bands( const std::vector<std::string> & chs )
{
  bands.clear();
  for (int s=0; s<chs.size(); s++)
    bands[ chs[s] ] = Eigen::MatrixXf::Zero(0,0);
}

void segsrv_t::calc_hjorths( const std::vector<std::string> & chs )
{
  hjorth.clear();
  for (int s=0; s<chs.size(); s++)
    hjorth[ chs[s] ] = Eigen::MatrixXf::Zero(0,0);
}

void segsrv_t::do_summaries( const std::string & ch , const int sr , const std::vector<double> * data ,
			     const bool do_band, const bool do_hjorth )
{

  if ( ! ( do_band || do_hjorth ) ) return ;
  
  // this is called per-epoch from populate() only 
  double fft_segment_size = 4;
  double fft_segment_overlap = 2;
  window_function_t window_function = WINDOW_TUKEY50;	   
  
  const double overlap_sec = fft_segment_overlap;
  const double segment_sec  = fft_segment_size;
  const int total_points = epoch_sec * sr;
  const int segment_points = segment_sec * sr;
  const int noverlap_points  = overlap_sec * sr;
  
  // implied number of segments
  int noverlap_segments = floor( ( total_points - noverlap_points) 
				 / (double)( segment_points - noverlap_points ) );
  
  Eigen::MatrixXf X = Eigen::MatrixXf::Zero( epoch_num , 6 );
  Eigen::MatrixXf H = Eigen::MatrixXf::Zero( epoch_num , 3 );
  Eigen::MatrixXf HH = Eigen::MatrixXf::Constant( epoch_num , 101 , std::numeric_limits<float>::quiet_NaN() );

  // for each clock-epoch
  for (int e=0; e<clock_epoch_num; e++)
    {

      // does this map to a full, observed epoch?
      if ( clk2sig_emap.find( e ) != clk2sig_emap.end() )
	{
	  
	  const int e2 = clk2sig_emap[ e ];
	  
	  // start in seconds [ based on extracted signal ] -> samples
	  //   (can do this, as this called prior to decimation
	  //    and so 'sr' is always valid)
	  const int idx = epoch_sec_starts[ e2 ] * sr ;
	  
	  // epoch is 'total_points' from index idx
	  
	  // copy part (can check later, but this should always fit in data[])
	  if ( idx + total_points > data->size() )
	    {
	      std::cerr << " sr = " << sr << " for " << ch << "\n";
	      std::cerr << " data->size()  = " << data->size()  << "\n"
			<< " idx " << idx << " " << total_points << " " << idx + total_points << "\n";
	      Helper::halt( "internal error in do_summaries()");
	    }
	  
	  std::vector<double> edata( total_points , 0 );
	  for (int i=0; i<total_points; i++)
	    edata[i] = (*data)[ idx + i ];
	  	  
	  if ( do_band )
	    {
	      
	      PWELCH pwelch( edata , 
			     sr , 
			     segment_sec , 
			     noverlap_segments , 
			     window_function ,
			     true , // use median
			     false ); // no SD calculated
	      
	      X(e2,0) = log10( pwelch.psdsum( SLOW ) + 0.00001 );
	      X(e2,1) = log10( pwelch.psdsum( DELTA ) + 0.00001);
	      X(e2,2) = log10( pwelch.psdsum( THETA ) + 0.00001);
	      X(e2,3) = log10( pwelch.psdsum( ALPHA ) + 0.00001);
	      X(e2,4) = log10( pwelch.psdsum( SIGMA ) + 0.00001);
	      X(e2,5) = log10( pwelch.psdsum( BETA ) + 0.00001);
	      
	    }
	  
	  if ( do_hjorth )
	    {
	      
	      double activity = 0 , mobility = 0 , complexity = 0;
	      MiscMath::hjorth( &edata , &activity , &mobility , &complexity , ! globals::legacy_hjorth );
	      H(e2,0) = log1p( activity );
	      H(e2,1) = mobility;
	      H(e2,2) = complexity;
	    }
	}
      
      // next epoch 
    }

  //
  // Normalize over epochs (Z-score, winsorized)
  //
  
  const float max_z = 2;
  
  if ( do_band )
    {
      
      // by band
      
      for (int j=0;j<6;j++)
       	{
       	  Eigen::VectorXf T = X.col(j);
       	  const int n = T.size();
	  
       	  double mean = T.mean();
       	  double std_dev = sqrt((T.array() - mean).square().sum() / (n - 1));	  
       	  T = ( T.array() - mean ) / std_dev ;	  	  
       	  X.col(j) = T; 
       	}
      
      // winsorize
      X = X.unaryExpr([&max_z](float x) { return x < -max_z ? -max_z : ( x > max_z ? max_z : x ); } );    

      // second-round norm (whole X)      
      double meanx = X.mean();
      double std_devx = sqrt((X.array() - meanx).square().sum() / (X.size() - 1));
      int nx = X.size() ;      
      X = ( X.array() - meanx ) / std_devx ; 
      X = X.unaryExpr([&max_z](float x) { return x < -max_z ? -max_z : ( x > max_z ? max_z : x ); } );
      
    }
  

  if ( do_hjorth )
    {
      // --> populate HH
      
      // normalize H2 and H3 only
      for (int j=0;j<3;j++)
	{
	  Eigen::VectorXf T = H.col(j);
	  const int n = T.size();
	  double mean = T.mean();
	  double std_dev = sqrt((T.array() - mean).square().sum() / (n - 1));
	  T = ( T.array() - mean ) / std_dev ;
	  
	  for (int i=0; i<n; i++)
	    if ( T[i] < -max_z ) T[i] = -max_z;
	    else if ( T[i] > max_z ) T[i] = max_z;
	  
	  H.col(j) = T; 
	}
      
      float amin = H.col(0).minCoeff();
      float amax = H.col(0).maxCoeff();
      Eigen::VectorXf A = H.col(0);
      float arng = amax - amin;

      if ( arng <= 0 ) arng = 1;
      A = ( A.array() - amin ) / ( amax - amin );
      A = A.array() * 50;
      
      for (int e=0; e<A.size(); e++)
	{
	  int h = A[e] ;
	  for (int y=1;y<h;y++)
	    {
	      HH(e, 50 + y ) = H(e,1);
	      HH(e, 50 - y ) = H(e,2);
	    }
	}
            
    }
  

  
  //
  // Splice into full epoch set 
  //

  Eigen::MatrixXf X2 = Eigen::MatrixXf::Constant( clock_epoch_num , 6   , std::numeric_limits<float>::quiet_NaN() );
  Eigen::MatrixXf H2 = Eigen::MatrixXf::Constant( clock_epoch_num , 101 , std::numeric_limits<float>::quiet_NaN() );
  
  if ( clock_epoch_num == epoch_num )
    {
      X2 = X;
      H2 = HH;
    }
  else
    {

      for (int e=0; e<clock_epoch_num; e++)
	if ( clk2sig_emap.find( e ) != clk2sig_emap.end() )
	  {	    
	    const int e2 = clk2sig_emap[ e ];
	    X2.row(e) = X.row(e2);
	    H2.row(e) = HH.row(e2);
	  }

    }
   
  //
  // Store
  //

  if ( do_band ) bands[ ch ] = X2;
  if ( do_hjorth ) hjorth[ ch ] = H2;

  // all done

}


Eigen::MatrixXf segsrv_t::get_bands( const std::string & ch )
{
  // slow, delta, theta, alpha, sigma, beta
  std::map<std::string,Eigen::MatrixXf>::const_iterator bb = bands.find( ch );
  if ( bb != bands.end() ) return bb->second;
  // else, just return an empty matrix
  return Eigen::MatrixXf::Zero( clock_epoch_num , 6 );  
}

Eigen::MatrixXf segsrv_t::get_hjorths( const std::string & ch )
{
  // h1,h2,h3
  std::map<std::string,Eigen::MatrixXf>::const_iterator hh = hjorth.find( ch );
  if ( hh != hjorth.end() ) return hh->second;
  // else, just return an empty matrix
  return Eigen::MatrixXf::Zero( clock_epoch_num , 3 );  
}


//
// Pull in signals
//

int segsrv_t::populate( const std::vector<std::string> & chs , const std::vector<std::string> & anns )
{
  
  // for now, can only init this once
  //  (i.e. do a single data pull - this avoids if the internal object is changed
  //   after calling, but the view is still open)
  //  that is, we only touch the EDF/instance on creation, or after update w/ 
  
  init();
  
  int count = 0;

  // add channels

  for (int i=0; i<chs.size(); i++)
    if ( add_channel( chs[i] ) ) ++count;
  
  // add annotations
  
  for (int i=0; i<anns.size(); i++)
    if ( add_annot( anns[i] ) ) ++count;

  // build interval tree
  etree.build( evts.begin() , evts.end() );
  
  // make sure we have some sensible default scaling
  
  set_scaling( count ,
	       anns.size() ,
	       1 , // yscale
	       1 , // ygroup
	       0.05, // yheader
	       0.05, // yfooter
	       count ? 0.10 : 1 ,
	       false ); // scaling_fixed_annot
  
  return count;
}



//
// get overall scale
//

std::vector<std::pair<double,double> > segsrv_t::get_time_scale() const
{

  std::vector<std::pair<double,double> > r;

  // 0..1 : viz dist (i.e. plots)    max = segsrv_t::get_ungapped_total_sec()
  // 0..1 : clock dist (i.e. time)   double segsrv_t::get_total_sec_original()

  
  uint64_t viz_min = 0LLU, viz_max = 0LLU;
  uint64_t clk_min = 0LLU, clk_max = 0LLU;

  std::set<interval_t>::const_iterator ss = segments.begin();

  viz_min = 0LLU;
  clk_min = ss->start;
  
  while ( ss != segments.end() )
    {
      viz_max += ss->duration();
      ++ss;
    }

  clk_max = get_total_sec_original() * globals::tp_1sec;

  
  uint64_t curr = 0LLU;
  ss = segments.begin();
  while	( ss != segments.end() )
    {
      // first point
      double pviz = curr / (double)viz_max;
      double pclk = ss->start / (double)clk_max;
      r.push_back( std::pair<double,double>( pviz , pclk ) );

      // end point
      curr += ss->duration();
      pviz = curr / (double)viz_max;
      pclk = ss->stop / (double)clk_max;
      r.push_back( std::pair<double,double>( pviz , pclk ) );
      
      ++ss;
    }
    
  return r;
}



std::set<std::pair<double,double> > segsrv_t::get_gaps() const
{

  // current window = awin, bwin
  // segments included
  
  uint64_t atp = awin * globals::tp_1sec ;
  uint64_t btp = bwin * globals::tp_1sec ;
  
  std::set<std::pair<double,double> > g;      
  std::set<interval_t>::const_iterator gg = gaps.begin();
  while ( gg != gaps.end() )
    {
      // std::cout << "\n considering gap: " << gg->as_string() << "\t"
      // 	    << gg->start << " " << gg->stop << "\n";
      // std::cout << " window = " << atp << " " << btp << "\n";
      
      // does any of the window fall in this gap?
      if ( btp > gg->start && atp < gg->stop )
	{
	  //std::cout << " this gap in window\n";
	  uint64_t start1 = atp > gg->start ? atp : gg->start;
	  uint64_t stop1 = btp < gg->stop ? btp : gg->stop;
	  g.insert( std::pair<double,double>( start1 * globals::tp_duration , stop1 * globals::tp_duration ) );
	  
	  // std::cout << "  --> add(tp) " << start1 << " " << stop1 << "\n";
	  // std::cout << "  --> add(s) " << start1 * globals::tp_duration << " " << stop1 * globals::tp_duration <
	  // < "\n";
	}
      ++gg;
    }
  
  return g;
}


bool segsrv_t::has_gaps( uint64_t atp, uint64_t btp , uint64_t * ungapped ) const
{
  bool has_gaps = false;
  uint64_t t = 0LLU;
  std::set<std::pair<double,double> > g;      
  std::set<interval_t>::const_iterator gg = gaps.begin();
  while ( gg != gaps.end() )
    {
      if ( btp > gg->start && atp < gg->stop )
	{
	  has_gaps = true;

	  uint64_t m1 = atp > gg->start ? atp : gg->start ;
	  uint64_t m2 = btp < gg->stop ? btp : gg->stop ;
	  t += m2 - m1;	  
	}
      ++gg;
    }

  if ( ungapped != NULL )
    *ungapped = ( btp - atp ) - t;
  return has_gaps;
}


double segsrv_t::get_total_sec() const
{
  return p->last_sec();
}

double segsrv_t::get_total_sec_original() const
{
  return p->last_sec_original();
}

Eigen::VectorXf segsrv_t::decimate( const Eigen::VectorXf & x0 , const int sr, const int q )
{
  if ( x0.size() == 0 ) return x0;
  
  // new sample rate (min 1 Hz) 
  const double sr2 = sr / (double)q ;
  const double fc = sr2 * 0.5 ;

  // forward-backward low-pass filtering
  iir_t iir1;
  iir1.init( BUTTERWORTH_LOWPASS , 2 , sr, fc ); 
  Eigen::VectorXf x1 = iir1.apply_bwlp_f( x0.reverse() ).reverse();
  
  iir_t iir2;
  iir2.init( BUTTERWORTH_LOWPASS , 2 , sr, fc );
  Eigen::VectorXf x2 = iir2.apply_bwlp_f( x1 );
  
  // decimate
  return x2(Eigen::seq(0,Eigen::last,q));
  
}
  

bool segsrv_t::add_channel( const std::string & ch )
{

  const int slot = p->edf.header.signal( ch );
  if ( slot == -1 ) return false;
    
  // (original) sample rate
  int sr = p->edf.header.sampling_freq( slot );
  
  // decimate?
  int decimation_fac = 1;
  if ( sr > max_samples_in )
    decimation_fac = sr / max_samples_in ;

  // get all data
  slice_t slice( p->edf , slot , p->edf.timeline.wholetrace() );
  const std::vector<double> * data = slice.pdata();
  const int n = data->size();

  // get some sensible ranges
  const double scaling_plwr = 0.05;
  const double scaling_pupr = 0.95;
  bool is_discrete;
  set_empirical_phys_ranges( ch, data->data() , n, scaling_plwr , scaling_pupr , &is_discrete );

  // show whether channel is discrete or continuous (for downsampling treatment)
  discrete[ ch ] = is_discrete;

  //  std::cout << " is dis = " << is_discrete << " " << ch << "\n";
  
  // do means, min/max & SD, as well as spectral/hjorth summaries? (on original data)
  do_summaries( ch, sr, data , bands.find( ch ) != bands.end() , hjorth.find( ch ) != hjorth.end() );  

  // get signal (copying nfull samples, so may contain zero-padded partial last records)
  // bur only copy n smaples over into the nfull space
  Eigen::VectorXf d = Eigen::VectorXf::Zero( n );
  for (int i=0; i<n; i++) d[i] = (*data)[i];

  // decimate?
  if ( decimation_fac > 1 )
    d = decimate( d , sr , decimation_fac );  

  // store
  sigmap[ ch ] = d;

  // store SR - note, just use original integer SR as a label
  //  for lookup - so it doesn't matter if we've decimated
  srmap[ ch ] = sr;

  // store new SR post any decimation
  decimated_srmap[ sr ] = sr / (double)decimation_fac;

  //  std::cout << "adding::: ch = " << ch << " " << sr << " " << sr / (double)decimation_fac << "\n";
  
  // do we already have a time-track?
  if ( tmap.find( sr ) == tmap.end() )
    {
      // get time-stamps
      const std::vector<uint64_t> * tp = slice.ptimepoints();
      
      // clock
      Eigen::VectorXf ts = Eigen::VectorXf::Zero( n );
      
      // for quick lookup: map of time-in-sec --> idx-int
      // scan original non-decimated tp-vector (length n)
      for (int i=0; i<n; i++)
	ts[ i ] = (*tp)[i] * globals::tp_duration;

      // decimate? (nb. need eval() to avoid aliasing issues)
      if ( decimation_fac > 1 )
      	ts = ts( Eigen::seq(0,Eigen::last,decimation_fac) ).eval(); 
            
      // and the time-series (sample # -> time)
      tmap[ sr ] = ts;
      
    }

  return true;
}


// set window (and track win valid_window whether this is non-null)
//  windows that fall entirely in a gap are invalid

bool segsrv_t::set_window( double a , double b )
{

  // max time (seconds, 1-tp-unit past end) 
  const double tmax = p->last_sec();
  
  // store seconds 
  awin = a < 0 ? 0 : ( a > tmax ? tmax : a ) ;
  bwin = b < 0 ? 0 : ( b > tmax ? tmax : b ) ;

  if ( awin > bwin )
    {
      double tmp = bwin;
      bwin = awin;
      awin = tmp;
    }

  
  std::set<int> srs;
  std::map<std::string,int>::const_iterator cc = srmap.begin();
  while ( cc != srmap.end() )
    {
      srs.insert( cc->second );
      ++cc;
    }

  bool all_okay = true;
  
  std::set<int>::const_iterator ss = srs.begin();
  while ( ss != srs.end() )
    {
      int aa = 0, bb = 0;
      const bool okay = get_tidx( awin, bwin , *ss , &aa , &bb );

      if ( okay )
	{
	  aidx[ *ss ] = aa;
	  bidx[ *ss ] = bb;
	}
      else
	{
	  all_okay = false;
	}
      ++ss;
    }
  
  
  // track
  valid_window = all_okay;
  
  return all_okay;
}


std::string segsrv_t::get_window_left_hms() const
{
  return get_hms( awin );  
}


std::string segsrv_t::get_window_right_hms() const
{
  return get_hms( bwin );
}


std::string segsrv_t::get_hms( const double s ) const
{
  clocktime_t t1 = edf_start;
  t1.advance_seconds( s );
  return t1.as_string( ':' );
}

std::map<double,std::string> segsrv_t::get_hour_ticks() const
{
  std::map<double,std::string> t;
  // take the entire record, and get the positions of the on-the-hour marks
  // we have edf_start

  clocktime_t t1 = edf_start;
  
  // whole record spans s0 to s1 seconds (units = sec past epoch)
  const double s0 = t1.seconds();
  const double s1 = s0 + get_total_sec_original();
  
  // is EDF on the hour? if so include?
  if ( t1.m == 0 && t1.s == 0 )
    t[ 0 ] = "| " + t1.as_string( ':' );
  
  // else, advance through the night
  while ( 1 ) 
    {
      // advance to top of the hour
      t1.advance_next_hr();

      // past end of recording?
      const double s = t1.seconds();      
      if ( s >= s1 ) break;

      // add, as fraction
      t[ s - s0 ] = "| " + t1.as_string( ':' );
      
    }
  
  return t;
}


std::map<double,std::string> segsrv_t::get_clock_ticks(const int n) const
{
  std::map<double,std::string> t;
  //  if ( ! valid_window ) return t;
  if ( n < 1 || n > 100 ) return t;
  
  // for the current window  
  // try to get up to than 'n' ticks at sensible values,
  // key = seconds (within curr window) 

  // determine size
  double sz = bwin - awin;

  // use this per-second count, starting at 0
  int per = sz / n;
  int aint = awin;
  if ( aint < awin ) ++aint;
  
  // start at (awin)  
  clocktime_t t1 = edf_start;
  t1.advance_seconds( aint );  
  t[ aint ] = "| " + t1.as_string( ':' );

  for (int p=1;p<n;p++)
    {
      aint += per;      
      t1.advance_seconds( per );
      t[ aint ] = "| " + t1.as_string( ':' );
    }

  return t;
  
}


Eigen::VectorXf segsrv_t::get_timetrack( const std::string & ch ) const
{
  std::map<std::string,int>::const_iterator ss = srmap.find( ch );
  if ( ss == srmap.end() || ! valid_window )
    {
      Eigen::VectorXf empty = Eigen::VectorXf::Zero(0);
      return empty;
    }
  int sr = ss->second;

  // return SR-specific time-track
  // (which may contain gaps etc so need whole thing)
  
  const Eigen::VectorXf & tt = tmap.find( sr )->second;
  const int aa = aidx.find(sr)->second;
  const int bb = bidx.find(sr)->second;  

  // throttle?
  if ( max_samples_out && ( bb - aa ) > max_samples_out )
    {

      // for discrete signals, decimate
      if ( is_discrete( ch ) )
	{
	  const int original_length = bb - aa;
	  const int reduction_factor = original_length / max_samples_out;      
	  Eigen::VectorXf f1 = tt.segment( aa , bb - aa );
	  return f1( Eigen::seq( 0 , Eigen::last , reduction_factor ) );
	}
      else // for contiuous signals, get envelope (tt-version)
	{
	  return envelope_timetrack( tt.segment( aa , bb - aa ) , xpixels );	  
	}
    }

  return tt.segment( aa , bb - aa );

}


// set scaling params                                                                                                                   
void segsrv_t::set_scaling( const int nchs , const int nanns ,
			    const double yscale , const double ygroup ,
			    const double yheader , const double yfooter ,
			    const double fixed_annot ,
			    const bool clip )
{

  // these parameters adjust the results of get_scaled_signal()
  //  i.e. we will know the signal slot for the requested signal from get_scaled_signal() n1

  // here,
  //   nchs   = total number of channels
  //   nanns  = total number of annotations (probably just keep as a single 'channel')
  
  //   always assume annots have a fixed height if >1 annot (i.e. always use the same scales)
  
  //   if n channels, by default we cut the y-axes into n bands,
  //      have the plot normalized within each

  //   ygroup = 0 all channels are at single midpoint
  //            1 default, as above

  //   yscale = 1 default, scaling within each 'bin'
  //            >1 or <1 (>0) implies multiplying signal by this factor
  //    (if >1, means it can overlap w/ other bins)

  //  yheader, yfooter - allow blank space at top, bottom
   
  scaling_nchs = nchs;
  scaling_nanns = nanns;
  scaling_yscale = yscale;
  scaling_ygroup = ygroup;
  scaling_yheader = yheader;
  scaling_yfooter = yfooter;
  scaling_fixed_annot = fixed_annot;
  scaling_clip = clip;
  
  if ( scaling_yheader < 0 ) scaling_yheader = 0;
  if ( scaling_yheader > 1 ) scaling_yheader = 1;
  if ( scaling_yfooter < 0 ) scaling_yfooter = 0;
  if ( scaling_yfooter > 1 ) scaling_yfooter = 1;
  if ( scaling_yheader + scaling_yfooter > 0.5 ) {
    scaling_yheader = 0;
    scaling_yfooter = 0;
  }
  
  if ( scaling_nchs < 0 ) scaling_nchs = 0;
  if ( scaling_nanns < 0 ) scaling_nanns = 0;
  if ( scaling_yscale < 0 ) scaling_yscale = 1;
  if ( scaling_ygroup < 0 ) scaling_ygroup = 0;
  if ( scaling_ygroup > 1 ) scaling_ygroup = 1;

  if ( scaling_fixed_annot < 0 ) scaling_fixed_annot = 0;
  if ( scaling_fixed_annot > 1 ) scaling_fixed_annot = 1;

  //
  // derive and store channel locs
  //

  // useable bin space
  double y1 = 1 - scaling_yheader - scaling_yfooter;
  
  // adjust for annotations (Y/N = 1+ vs 0)
  if ( scaling_nanns )
    y1 -= scaling_fixed_annot;
  
  // get spacing of signals within 'y1'
  //  ygroup = 1  then space is y1 / nsig
  //         = 0  then space is 0 
  // with all centered around y1/0.5

  if ( scaling_nchs )
    {
      std::vector<double> ymids( scaling_nchs );

      // initial gaps between y-midpoints of each signal
      double ygaps = scaling_ygroup / (double)scaling_nchs ;

      double ymean = 0;
      for (int i=0; i<scaling_nchs; i++)
	{
	  
	  //ymids[i] = i * ygaps ; 
	  // for reverse order, go (scaling_nchs-i+1)
	  ymids[i] = ( scaling_nchs - i + 1 ) * ygaps ; 
	  ymean += ymids[i];
	}
      ymean /= (double)scaling_nchs;
      
      // diff between midpoint of y1
      double delta = 0.5 - ymean;
      
      // adjust, and map to y1
      for (int i=0; i<scaling_nchs; i++)
	ymids[i] = ( ymids[i] + delta ) * y1;
      
      // prior to yscale (magnification) between whole y1 (if ygroup = 0 )
      // and y1/n (if ygroup = 1 )
      // 
      double yheight = (1-scaling_ygroup) * y1 + scaling_ygroup * ( y1 / (double)scaling_nchs ); 

      // apply magnification
      yheight *= scaling_yscale;

      scaling_lwr.resize( scaling_nchs );
      scaling_upr.resize( scaling_nchs );

      // get bands, adding in footer also
      for (int i=0; i<scaling_nchs; i++)
        {
	  scaling_lwr[i] = ymids[i] - 0.5 * yheight + yfooter ;
	  scaling_upr[i] = ymids[i] + 0.5 * yheight + yfooter ;
	}
                 
    }

}


// set physical limits (mins) for a given channel
void segsrv_t::fix_physical_scale( const std::string & ch , const double lwr, const double upr )
{
  std::pair<double,double> lu( lwr < upr ? lwr : upr , lwr < upr ? upr : lwr );
  phys_ranges[ ch ] = lu;    
}

// set empirical precalculated (percentile-based) scale
void segsrv_t::empirical_physical_scale( const std::string & ch )
{
  phys_ranges[ ch ] = empirical_phys_ranges[ ch ] ;
}


void segsrv_t::free_physical_scale( const std::string & ch )
{
  //   std::map<std::string, std::pair<double,double> > phys_ranges;
  std::map<std::string, std::pair<double,double> >::const_iterator pp = phys_ranges.find( ch );
  if ( pp == phys_ranges.end() ) return;
  // erase
  phys_ranges.erase( pp );
}


// calc. robust reasonable ranges
template<typename T>
void segsrv_t::set_empirical_phys_ranges( const std::string & ch , const T * data , const int n , 
					  const double plwr , const double pupr ,
					  bool * is_discrete )
{

  //  std::cout << " n = " << n << "\n";
  // for (int i=0; i<n; i++)
  //   {
  //     double d = static_cast<double>(data[i]);  
  //     std::cout << "  " << d << "\n";
  //   }
      
  // fixed at 9/95 for now 9i.e. this ignores plwr/pupr
  axis_stats_t stats = compute_axis_stats( data , n );

  //  std::cout << " is_discrete = " << stats.is_discrete << "\n";
  
  // track back?
  if ( is_discrete != NULL )
    *is_discrete = stats.is_discrete;

  //  std::cout << "stats.min_val, stats.max_val  = " << stats.min_val << " " <<  stats.max_val  << "\n";

  // for discrete signals, we want to keep the original values
  // treat as 'categorical' (e.g. 0/1 flag) if 10 or fewer discrete values
  if ( stats.is_discrete )
    {
      empirical_phys_ranges[ ch ] = std::pair<double,double>( stats.min_val, stats.max_val );
      return;
    }
  
  // else percentiles
  empirical_phys_ranges[ ch ] = std::pair<double,double>( stats.p5 , stats.p95 );

  //  std::cout << " setting " <<  stats.p5 << " " << stats.p95 << "\n";

}

bool segsrv_t::get_yscale_signal( const int n1 , double * lwr, double * upr ) const
{
  // for a signal, given it's place in the slots, return the scaled (within 0.1)
  // values for this plot  
  if ( n1 < 0 || n1 >= scaling_nchs ) return false;
  *lwr = scaling_lwr[n1];
  *upr = scaling_upr[n1];  
  return true;
}



// get scaled signals (0-1 give other annots)                                
Eigen::VectorXf segsrv_t::get_scaled_signal( const std::string & ch , const int n1 )
{

  if ( ! valid_window ) return Eigen::VectorXf::Zero(0);
  
  // allow to specify n1 separate from actual channel list (i.e. may be
  // showing a subset of all channes, and so we would have recalled 
  // set_scaling() multiple times after populate() [ which is only called once ]
  
  Eigen::VectorXf s = get_signal( ch );
  
  // fixed physical scaling, or auto-scaling?
  // if fixed scaling, optionally, clip if above/below
  
  float smin, smax; 

  std::map<std::string, std::pair<double,double> >::const_iterator pp =	phys_ranges.find( ch );
  if ( pp == phys_ranges.end() )
    {
      // auto-scaling  
      smin = min_skip_nan( s );
      smax = max_skip_nan( s ); // s.maxCoeff();

      // store
      window_phys_range[ ch ] = std::pair<double,double>( smin, smax );
    }
  else
    {
      // use pre-specified values
      smin = pp->second.first;
      smax = pp->second.second;

      // optionally, clip?
      if ( scaling_clip )
	{
	  
	  // compare to observed
	  float omin = min_skip_nan(s); // s.minCoeff();
	  float omax = max_skip_nan(s); // s.maxCoeff();
	  
	  const int n = s.size();      
	  const bool fix_lwr = omin < smin ;
	  const bool fix_upr = omax > smax ;
	  
	  if ( fix_lwr || fix_upr )
	    {
	      for (int i=0;i<n;i++)
		{
		  if ( fix_lwr && s[i] < smin ) s[i] = smin;
		  if ( fix_upr && s[i] > smax ) s[i] = smax;
		}
	    }
	  
	  // store (observed values)
	  window_phys_range[ ch ] = std::pair<double,double>( fix_lwr ? smin : omin , fix_upr ? smax : omax );
	}
      else // store just fixed values
	window_phys_range[ ch ] = std::pair<double,double>( smin , smax );

    }

  const double srange = smax - smin;
  // special case of smin == smax --> set whole signal to 0.5 
  if ( srange < 1e-6f )
    s = Eigen::VectorXf::Zero( s.size() ).array() + 0.5 ;
  else // normalize to 0..1  X = ( X - min ) / ( max - min )    
    s = ( s.array() - smin ) / (float)srange;
  
  // rescale to display units
  double lwr = 0,  upr = 1;
  const bool okay = get_yscale_signal( n1 , &lwr, &upr );

  // to clip between 0 and 1:
  //     vec = vec.cwiseMax(0.0).cwiseMin(1.0);
    
  if ( okay )
    {
      // make signal
      s = s.array() * ( upr - lwr ) + lwr;

      // also track to reconstruct particular y-values from get_scaled_y()
      // to be called **only after this whole-window function**
      track_ylwr[ ch ] = lwr;
      track_yupr[ ch ] = upr;
      track_smin[ ch ] = smin;
      track_smax[ ch ] = smax;
    }
  
  // return
  return s;
      
}

float segsrv_t::get_scaled_y( const std::string & ch , const float y ) const
{
  std::map<std::string,float>::const_iterator ylwr = track_ylwr.find( ch );
  if ( ylwr == track_ylwr.end() ) return -9;
  std::map<std::string,float>::const_iterator yupr = track_yupr.find( ch );
  std::map<std::string,float>::const_iterator smin = track_smin.find( ch );
  std::map<std::string,float>::const_iterator smax = track_smax.find( ch );
    
  double t = ( y - smin->second ) / (float)( smax->second - smin->second );
  return t * ( yupr->second - ylwr->second ) + ylwr->second;

}


Eigen::VectorXf segsrv_t::get_signal( const std::string & ch ) const 
{
  
  std::map<std::string,int>::const_iterator ss = srmap.find( ch );

  if ( ss == srmap.end() || ! valid_window )
    {
      Eigen::VectorXf empty = Eigen::VectorXf::Zero(0);
      return empty;
    }
  
  int sr = ss->second;

  // filtered or original signal?
  bool is_filtered = filtered.find( ch ) != filtered.end();
      
  // return data (signal) - optionally filtered
  const Eigen::VectorXf & data = is_filtered ? sigmap_f.find( ch )->second : sigmap.find( ch )->second;
  const int aa = aidx.find(sr)->second;
  const int bb = bidx.find(sr)->second;  


  // throttle?
  if ( max_samples_out && ( bb - aa ) > max_samples_out )
    {
      const int original_length = bb - aa;
      const int reduction_factor = original_length / max_samples_out;

      // std::cout << " reduction_factor = " << reduction_factor << "\n";
      // int spp = std::max(1, original_length / xpixels );  // samples per pixel
      // std::cout << " spp = " << spp << "\n";


      // for discrete signals, decimate                                                                                                       
      if ( is_discrete( ch ) )
	{
	      
	  // decimation here on second level - but note that the SR may
	  // have been adjusted from initial decimation; note also that SR
	  // is cast to an int here, but for this single purpose (of
	  // designing a low-pass anti-aliasing filter) this should not
	  // matter.
	  
	  return decimate( data.segment( aa , bb - aa ) , decimated_srmap.find( sr )->second , reduction_factor );
	}
      else // for continuous signals
	{
	  return envelope_signal( data.segment( aa , bb - aa ) , xpixels );
	}
    }

  // else return full
  return data.segment( aa , bb - aa );
    
}


// given two times and a sample rate, get indices
bool segsrv_t::get_tidx(double a, double b, int sr, int *aa, int *bb) const
{
  auto it = tmap.find(sr);
  if (it == tmap.end()) return false;
  
  const Eigen::VectorXf &ts = it->second;
  const int n = ts.size();
  if (n == 0) return false;
  
  const float af = static_cast<float>(a);
  const float bf = static_cast<float>(b);
  
  const float *beg = ts.data();
  const float *end = beg + n;
  
  // first sample with time >= a
  const float *pa = std::lower_bound(beg, end, af);
  if (pa == end) return false;
  
  // first sample with time >= b
  const float *pb = std::lower_bound(beg, end, bf);
  bool at_end = false;
  if (pb == end) {
    // move to last valid sample, then we'll add +1 below
    --pb;
    at_end = true;
  }
  
  // gap case: window [a,b) lies entirely between two samples
  if (pa == pb) return false;
  
  const int ia = int(pa - beg);
  const int ib = int(pb - beg);
  
  *aa = ia;
  *bb = at_end ? ib + 1 : ib;
  return true;
}


Eigen::MatrixXf segsrv_t::get_summary_stats( const std::string & ch )
{
  Eigen::MatrixXf r = Eigen::MatrixXf::Zero(0,0);
  return r;
}

Eigen::VectorXf segsrv_t::get_summary_timetrack( const std::string & ch ) const
{
  Eigen::VectorXf t = Eigen::MatrixXf::Zero(0,0);
  return t;
}



//
// Annotations
//

bool segsrv_t::add_annot( const std::string & ch )
{
  
  annot_t * annot = p->edf.annotations->find( ch );
  if ( annot == NULL ) return false;
  if ( annot->interval_events.size() == 0 ) return false;
  
  annot_map_t::const_iterator ii = annot->interval_events.begin();   
  while ( ii != annot->interval_events.end() )
    {
      const instance_idx_t & instance_idx = ii->first;

      bool has_id = ii->first.id != "" && ii->first.id != ".";
      bool has_ch = ii->first.ch_str != "" && ii->first.ch_str != ".";

      std::string meta; 

      if ( has_id )
	{
	  meta = ii->first.id;
	  if ( has_ch )
	    meta += "; " + ii->first.ch_str;
	}
      else
	{
	  if ( has_ch )
	    meta = ii->first.ch_str;
	  else
	    meta = ".";
	}
      
      evt_t evt( ii->first.interval , ii->first.parent->name , meta );

      evts.insert( evt );

      ++ii;
    }
  
  return true;
}





// get events from the current window
std::map<std::string,std::vector<std::pair<double,double> > > segsrv_t::fetch_evts() const
{

  std::map<std::string,std::vector<std::pair<double,double> > > r;
  
  // current window
  uint64_t atp = awin * globals::tp_1sec;
  uint64_t btp = bwin * globals::tp_1sec;


  // legacy search
  if ( 0 )
    {
      uint64_t endtp = btp == 0LLU ? 0LLU : btp - 1LLU;
      interval_t win( atp , btp );  
  
      // set 0-dur marker at end of window (at last included point, end-1)
      evt_t win_end( interval_t( endtp, endtp ) , "__luna#window#marker__" ) ;
      
      // get upper bound of this (i.e. next event that is /not/ in the window) 
      std::set<evt_t>::const_iterator ee = evts.upper_bound( win_end );
      
      while ( 1 )
	{
	  // check
	  if ( ee == evts.begin() ) break;
	  
	  // step back
	  --ee;
	  
	  // done?
	  if ( ee->interval.is_before( win ) ) break;
	  
	  // overlap?
	  if ( ee->interval.overlaps( win ) )
	    r[ ee->name ].push_back( std::pair<double,double>( ee->interval.start_sec() , ee->interval.stop_sec() ) ) ; 
	  
	  // next event      
	}
    }


  //
  // interval tree implementation
  //
  
  auto hits = etree.query_ptrs( atp , btp );
  
  for ( const auto & p : hits)
      r[ p->name ].push_back( std::pair<double,double>( p->interval.start_sec() , p->interval.stop_sec() ) );
  
  return r;
}


// for lunascope in particular (include instance IDs)
std::vector<std::vector<std::string> >
segsrv_t::fetch_all_evts_with_inst_ids( const std::vector<std::string> & avec ,
					const bool hms ) const
{

  // !hms : return annot | inst | startsec | stopsec
  //  hms : return annot | inst | hh:mm:ss | startsec | +duration
  //               0       1      2          3          4

  const int cols = hms ? 5 : 4;
    
  // hms uses clocktime_t edf_start, if it is valid
    
  std::vector<std::vector<std::string> > r;
  
  std::set<std::string> aset = Helper::vec2set( avec );
  
  std::set<evt_t>::const_iterator ee = evts.begin();
  
  while ( ee != evts.end() )
    {
      if ( aset.find( ee->name ) != aset.end() )
	{
	  
	  std::vector<std::string> row( cols );

	  row[0] = ee->name;
	  row[1] = ee->meta;
	  
	  std::string str = ee->name + " | " + ee->interval.as_string(3,"-") ; // 3dp
	  
	  if ( hms )
	    {

	      // clock time (start)
	      if  (edf_start.valid )
		{
		  clocktime_t t = edf_start;
		  t.advance_tp( ee->interval.start );
		  row[2] = t.as_string( ':' , false ); // F = do not include any fractional seconds
		}
	      else
		row[2] = "?";

	      // start in seconds
	      row[3] = Helper::dbl2str( ee->interval.start_sec() , 3 );
	      
	      // duration in seconds	      	      
	      row[4] = Helper::dbl2str( ee->interval.duration_sec() , 3 );
	      
	    }
	  else
	    {
	      // !hms : return annot | inst | startsec | stopsec
	      //               0       1      2          3
	      
	      // start in seconds
	      row[2] = Helper::dbl2str( ee->interval.start_sec() , 3 );
	      
	      // start in seconds
	      row[3] = Helper::dbl2str( ee->interval.stop_sec() , 3 );
	      
	    }

	  r.push_back( row );
	  
	}
      ++ee;
    }
  
  return r;

}


// for selection window
std::vector<std::string> segsrv_t::fetch_all_evts( const std::vector<std::string> & avec , const bool hms ) const
{

  // !hms : return annot | startsec-stopsec
  //  hms : return annot | startsec-stopsec | hh:mm:ss | +duration

  // hms uses clocktime_t edf_start, if it is valid
    
  std::vector<std::string> r;
  std::set<std::string> aset = Helper::vec2set( avec );
  
  std::set<evt_t>::const_iterator ee = evts.begin();

  
  while ( ee != evts.end() )
    {
      if ( aset.find( ee->name ) != aset.end() )
	{

	  std::string str = ee->name + " | " + ee->interval.as_string(3,"-") ; // 3dp
	  
	  if ( hms )
	    {
	      // duration in seconds
	      const std::string dur = Helper::dbl2str( ee->interval.duration_sec() , 3 ); // 3 dp

	      if  (edf_start.valid )
		{
		  
		  // start hh:mm:ss
		  clocktime_t t = edf_start;
		  t.advance_tp( ee->interval.start );
		  std::string clock = t.as_string( ':' , false ); // F = do not include any fractional seconds
		  
		  str += " | " + clock + " | " + dur ;
		}
	      else
		str += " | ? | " + dur ;
	    }

	  r.push_back( str );
	  
	}
      ++ee;
    }
  
  return r;
}

// set format for annotation plots
void segsrv_t::set_annot_format6( const bool b )
{
  annot_format6 = b;
}

// compile a set of selected events for the current window
void segsrv_t::compile_evts( const std::vector<std::string> & anns )
{
  // clear window
  
  // ultimate (in 6-val format or standard) for plotting
  compiled_annots_times.clear();
  compiled_annots_end_times.clear();
  compiled_annots_stacks.clear();
  compiled_annots_end_stacks.clear();

  // for fixed position ordering: total number of annots (whether shown or not in this window)
  const int na = anns.size();
  
  // working intermediate
  std::map<std::string,std::vector<std::pair<double,double> > > annots_times;
  std::map<std::string,std::vector<std::pair<double,double> > > annots_stacks;
  
  // get all events that overlap this window
  std::map<std::string,std::vector<std::pair<double,double> > > wevts = fetch_evts();
  
  // make uniform time-line: extracted from in-window events
  std::set<fevt_t> xevts;
  
  // count # of unique annot types in this window
  std::map<std::string,int> amap;
  
  for (int a=0;a<anns.size();a++)
    {
      // not present?
      if ( wevts.find( anns[a] ) == wevts.end() ) continue;

      const std::string & aname = anns[a];
      
      // track position for y
//      amap[ aname ] = amap.size() ; 
      // alternate: simple scaling (i.e. fixed w.r.t to /all/ annots, not jsut those in the window)
      amap[ aname ] = a;
      
      // events
      const std::vector<std::pair<double,double> > & e = wevts.find( aname )->second;

      for (int i=0;i<e.size();i++)
	xevts.insert( fevt_t( e[i].first , e[i].second , aname ) );
      
    }

  // empty? all done
  if ( xevts.size() == 0 ) return;


  //
  // Simple y-axis positons
  // 

  int maxd = amap.size();
  
  std::set<fevt_t>::const_iterator xx = xevts.begin();
  while ( xx != xevts.end() )
    {
      	  
      // add times (optionally, clipping at window- boundaries) 
      const std::pair<double,double> p2(
					clip_xaxes ? std::max(xx->start, awin) : xx->start,
					clip_xaxes ? std::min(xx->stop,  bwin) : xx->stop
					);
      
      annots_times[ xx->name ].push_back( p2 );
      	  
      // determine depth (N-1) [ scale after when we know max depth ]  
      annots_stacks[ xx->name ].push_back( std::pair<double,double>( amap[ xx->name ] , 0 ) );
      
      // next annot
      ++xx;
    }
  
  // scale depth into plotting values:
  //  annots go from ( 1 - scaling_yheader ) to ( 1 - scaling_fixed_annot - scaling_yheader )
  //   (top to bottom, going down)
  const double abase = 1.0 - scaling_yheader;
  //const double abase = 1.0 - scaling_fixed_annot - scaling_yheader;
  
  // maxd +1 i.e to allow for height of top annot
// OLD:   const double a1height = scaling_fixed_annot / (double)(maxd);
  // alternate: simpler y-axis scaling 
  const double a1height = scaling_fixed_annot / (double)(na);
  
  std::map<std::string,std::vector<std::pair<double,double> > >::iterator cc = annots_stacks.begin();
  while ( cc != annots_stacks.end() )
    {
      
      std::vector<std::pair<double,double> > & ds = cc->second;
      
      for (int i=0; i<ds.size(); i++)
	{
	  double y = ds[i].first; // scale from 1... maxd
	  y = abase - y * a1height; // lower y;
	  ds[i].first = y;
	  ds[i].second = y - a1height;
	}      
      ++cc;
    }

  //
  // for annot 'i' we now have (x1,x2) and (y1,y2) 


  //
  //  if annot_format6 is True, then for easier plotting, we want to
  //  have six entires:
  //  
  //     (x1,x2,x2,x1,x1,NA)
  //     (y1,y1,y2,y2,y1,NA)
  //  i.e. for plotly to make the traces and fill in a rectange
  //

  //
  // otherwise, we want to pass two separate n-sized vectors: starts, stops, and y-offsets
  //
  
  
  cc = annots_times.begin();
  while ( cc != annots_times.end() )
    {
      
      std::vector<std::pair<double,double> > & xs = cc->second;
      std::vector<std::pair<double,double> > & ys = annots_stacks[ cc->first ];

      const int n = xs.size();

      //
      // format6: 
      //

      if ( annot_format6 )
	{

	  std::vector<float> xx( 6 * n );
	  std::vector<float> yy( 6 * n );

	  int xidx = 0;
	  int yidx = 0;
	  for (int i=0; i<n; i++)
	    {
	      xx[xidx++] = xs[i].first ;
	      xx[xidx++] = xs[i].second ;
	      xx[xidx++] = xs[i].second ;
	      xx[xidx++] = xs[i].first ;
	      xx[xidx++] = xs[i].first ;
	      xx[xidx++] = std::numeric_limits<float>::quiet_NaN() ;
	      
	      yy[yidx++] = ys[i].first ;
	      yy[yidx++] = ys[i].first ;
	      yy[yidx++] = ys[i].second ;
	      yy[yidx++] = ys[i].second ;
	      yy[yidx++] = ys[i].first ;
	      yy[yidx++] = std::numeric_limits<float>::quiet_NaN() ;
	      
	    }

	  // store
	  compiled_annots_times[ cc->first ] = xx;
	  compiled_annots_stacks[ cc->first ] = yy;
	}
      else
	{

	  //
	  // or not format6
	  //
	  
	  std::vector<float> xx1( n );
	  std::vector<float> xx2( n );
	  std::vector<float> yy1( n );
	  std::vector<float> yy2( n );

	  for (int i=0; i<n; i++)
	    {
	      xx1[i] = xs[i].first;
	      xx2[i] = xs[i].second;

	      yy1[i] = ys[i].first;
	      yy2[i] = ys[i].second;
	    }

	  // store
	  compiled_annots_times[ cc->first ] = xx1;
	  compiled_annots_end_times[ cc->first ] = xx2;

	  compiled_annots_stacks[ cc->first ] = yy1;
	  compiled_annots_end_stacks[ cc->first ] = yy2;

	}
      ++cc;
    }
  
}

// given a compilation (subset of all evts), get evts for a particular class
std::vector<float> segsrv_t::get_evnts_xaxes( const std::string & ann ) const
{
  std::map<std::string,std::vector<float> >::const_iterator aa = compiled_annots_times.find( ann );
  if ( aa != compiled_annots_times.end() ) return aa->second;
  std::vector<float> empty;
  return empty;  
}

// given a compilation (subset of all evts), get evts for a particular class
std::vector<float> segsrv_t::get_evnts_xaxes_ends( const std::string & ann ) const
{
  std::map<std::string,std::vector<float> >::const_iterator aa = compiled_annots_end_times.find( ann );
  if ( aa != compiled_annots_end_times.end() ) return aa->second;
  std::vector<float> empty;
  return empty;  
}

// given a compilation (subset of all evts), get y-axis stacking for a particular class
std::vector<float> segsrv_t::get_evnts_yaxes( const std::string & ann ) const
{
  std::map<std::string,std::vector<float> >::const_iterator aa = compiled_annots_stacks.find( ann );
  if ( aa != compiled_annots_stacks.end() ) return aa->second;
  std::vector<float> empty;
  return empty;	
}

// given a compilation (subset of all evts), get y-axis stacking for a particular class
std::vector<float> segsrv_t::get_evnts_yaxes_ends( const std::string & ann ) const
{
  std::map<std::string,std::vector<float> >::const_iterator aa = compiled_annots_end_stacks.find( ann );
  if ( aa != compiled_annots_end_stacks.end() ) return aa->second;
  std::vector<float> empty;
  return empty;	
}



//
// sigmod_t
//

void sigmod_t::clear_mod( const std::string & mod_label )
{
  if ( mod_bins.find( mod_label ) == mod_bins.end() )
    return;

  // clear
  mod_bins.erase( mod_label );
  mod_tt.erase( mod_label );

}

void sigmod_t::make_mod( const std::string & mod_label ,
			 const std::string & mod_ch ,			 
			 const std::string & type ,
			 const std::vector<double> & sos , // filter if populated
			 const bool ylim , const double ylwr , const double yupr )
{

  // has this already been created?
  if ( mod_bins.find( mod_label ) != mod_bins.end() )
    return;

  // given a filter to apply (presumably for filter-Hilbert but allow for 'raw' amplitude too)
  const bool has_filter = sos.size() != 0 ;
  
  // type of modulation
  const bool mod_raw = type == "raw";  // take signal as is
  const bool mod_amp = type == "amp";  // filter-Hilbert --> amplitude
  const bool mod_pha = type == "phase"; // filter-Hilbert --> magnitude
  if ( ! ( mod_raw || mod_amp || mod_pha ) )
    return;
  
  // get original mod-sig
  std::map<std::string,Eigen::VectorXf>::const_iterator ss = parent->sigmap.find( mod_ch );
  if ( ss == parent->sigmap.end() ) return;

  // get entire signal 
  Eigen::VectorXf s = parent->sigmap[ mod_ch ];
  
  // filter?
  if ( has_filter ) 
    {            
      sos_filter_t f( sos );
      std::size_t M = sos.size() / 6;
      int pad = std::min(std::max(32, (int)(16*M) ), (int)(s.size() /8 ) );
      sos_filter_prime_with_reflection( f , s , pad );
      f.process( s );
    }

  // currently, yscale not used
  
  // upper limit of each bin
  std::vector<double> cuts( nbins );
  
  // Hilbert transform?
  if ( ! mod_raw )
    {
      std::vector<double> d(s.data(), s.data() + s.size());      
      hilbert_t h( d );
      const std::vector<double> * out;
      if ( mod_amp ) out = h.magnitude();
      else if ( mod_pha ) out = h.phase();
      
      if ( out->size() != s.size() )
	Helper::halt( "internal error: hilbert output size" );
      for (int i = 0; i < s.size(); ++i)
	s[i] = static_cast<float>((*out)[i]);
      
    }
  
  // get range
  double smin = s.minCoeff();
  double smax = s.maxCoeff();

  // find segments
  const float span = smax - smin;
  if (span <= 0) return;  
  const float inv_span = nbins / span;
  
  // compute bin indices for all s
  Eigen::ArrayXf t = (s.array() - smin) * inv_span;
  Eigen::ArrayXi b = t.floor().cast<int>();
  b = b.min(nbins - 1); // clamp right edge
  
  // save modulator: a) need bins, b) time-track, both for entire recording
  mod_bins[ mod_label ] = b;
  
  //  pull entire signal
  mod_tt[ mod_label ] = parent->tmap.find( parent->srmap[ mod_ch ] )->second; 

}


void sigmod_t::apply_mod( const std::string & mod_label ,
			  const std::string & ch ,
			  const int slot )
  
{

  // did we have this modulator made correctly?
  status = mod_tt.find( mod_label ) != mod_tt.end() ;

  if ( ! status ) return;
  
  Eigen::VectorXf tX = parent->get_timetrack( ch );

  Eigen::VectorXf X = parent->get_scaled_signal( ch , slot ); // will set y0 also
  
  // // Inputs: tS (N_S), bS (N_S), tX (N_X), X (N_X)
  segments = bin_X_by_Sbins( mod_tt[ mod_label ] , // tS
			     mod_bins[ mod_label ] , // bS
			     parent->get_timetrack( ch ) , // tX
			     parent->get_scaled_signal( ch , slot ) ); // X  
  
}

 
Eigen::VectorXf sigmod_t::get_timetrack( const int bin ) const
{
  if ( ! status ) return Eigen::VectorXf::Zero(0);
  return segments.t[ bin ];
}

Eigen::VectorXf sigmod_t::get_scaled_signal( const int bin ) const
{
  if ( ! status ) return Eigen::VectorXf::Zero(0);
  return segments.x[ bin ];
}





// Left-constant over [tS[i], tS[i+1)).
// Out-of-range policy: returns -1 for tx < tS[0] or tx >= tS.back().
Eigen::ArrayXi sigmod_t::bins_from_Sbins_at_X(const Eigen::VectorXf& tS,
                                              const Eigen::ArrayXi&  bS,
                                              const Eigen::VectorXf& tX)
{
    assert(tS.size() == bS.size());
    const int nS = static_cast<int>(tS.size());
    const int nX = static_cast<int>(tX.size());

    Eigen::ArrayXi bX = Eigen::ArrayXi::Constant(nX, -1);
    if (nS == 0 || nX == 0) return bX;

    // Start i at the last index with tS[i] <= tX[0] using binary search
    int i = 0;
    {
        const float* s0 = tS.data();
        const float* sN = s0 + nS;
        const float* it = std::upper_bound(s0, sN, tX[0]); // first > tX[0]
        i = static_cast<int>(std::max<std::ptrdiff_t>(0, (it - s0) - 1));
    }

    for (int j = 0; j < nX; ++j) {
        const float tx = tX[j];

        // Out-of-range to -1; change to (tx < tS[0]) only if you want last bin to hold to +inf
        if (tx < tS[0] || tx >= tS[nS - 1]) { bX[j] = -1; continue; }

        while (i + 1 < nS && tS[i + 1] <= tx) ++i;
        bX[j] = bS(i);
    }
    return bX;
}

// Build per-bin (tX, X) arrays with NaN separators so pyqtgraph connect='finite' draws segments only.
sigmod_segment_t sigmod_t::segments_from_bins_dual(const Eigen::VectorXf& tX,
                                                   const Eigen::VectorXf& X,
                                                   const Eigen::ArrayXi&  bX)
{
    assert(tX.size() == X.size());
    assert(X.size() == bX.size());
    assert(nbins > 0);

    const int n = static_cast<int>(X.size());
    const float NaN = std::numeric_limits<float>::quiet_NaN();

    // First pass: count points per bin and number of segment separators
    std::vector<int> counts(nbins, 0), seps(nbins, 0);
    int prev = -1;
    for (int i = 0; i < n; ++i) {
        const int bi = bX(i);
        if (bi != prev && prev >= 0)
	  {
	    counts[prev]++;  
	    seps[prev]++;    // leaving a run -> separator in prev bin
	  }

	if (bi >= 0 && bi < nbins)
	  counts[bi]++;      // accept only valid bins
        prev = (bi >= 0 && bi < nbins) ? bi : -1;
    }
    if (prev >= 0) seps[prev]++;                      // trailing separator for last run
    
    // Allocate output with exact sizes
    sigmod_segment_t out;
    out.t.resize(nbins);
    out.x.resize(nbins);
    for (int k = 0; k < nbins; ++k) {
        const int m = counts[k] + seps[k];                    // one NaN per run plus all points
        out.t[k].resize(m);
        out.x[k].resize(m);
    }

    // Second pass: fill data with NaN separators
    std::vector<int> pos(nbins, 0);
    prev = -1;
    for (int i = 0; i < n; ++i) {
        const int bi = bX(i);

        if (bi != prev && prev >= 0) {
	  // write the extra point + NaN in previous bin
	  int p = pos[prev]++;
	  out.t[prev][p] = tX[i];
	  out.x[prev][p] = X[i];
	  
	  p = pos[prev]++;     // advance pos[prev] again
	  out.t[prev][p] = NaN;
	  out.x[prev][p] = NaN;
	  
        }

        if (bi >= 0 && bi < nbins) {
            int p = pos[bi]++;
            out.t[bi][p] = tX[i];
            out.x[bi][p] = X[i];
        }

        prev = (bi >= 0 && bi < nbins) ? bi : -1;
    }

    if (prev >= 0) {
        int p = pos[prev]++;
        out.t[prev][p] = NaN;
        out.x[prev][p] = NaN;
    }

    return out;
}

// End-to-end: precomputed whole-night S bins -> per-bin segments on X window
sigmod_segment_t sigmod_t::bin_X_by_Sbins(const Eigen::VectorXf& tS,
                                          const Eigen::ArrayXi&  bS,
                                          const Eigen::VectorXf& tX,
                                          const Eigen::VectorXf& X)
{
    assert(tS.size() == bS.size());
    assert(tX.size() == X.size());

    const Eigen::ArrayXi bX = bins_from_Sbins_at_X(tS, bS, tX);
    return segments_from_bins_dual(tX, X, bX);
}


// filtering
void segsrv_t::apply_filter( const std::string & ch , const std::vector<double> & sos )
{    
  // indicate this channel is filtered
  filtered.insert( ch );

  // copy whole signal
  sigmap_f[ ch ] = sigmap[ ch ];
  Eigen::VectorXf & data = sigmap_f[ ch ];
  
  // filter
  sos_filter_t f( sos );
  std::size_t M = sos.size() / 6;
  int pad = std::min(std::max(32, (int)(16*M) ), (int)(data.size() /8 ) );
  sos_filter_prime_with_reflection( f , data, pad );
  f.process( data );
  
  // reset 5/95 percentiles
  empirical_phys_ranges_orig[ ch ] = empirical_phys_ranges[ ch ];
  set_empirical_phys_ranges( ch , data.data(), data.size(), 0.05, 0.95 );
}
  
void segsrv_t::clear_filter( const std::string & ch )
{

  // indicate not filtered
  filtered.erase( ch );
  
  // restore scaling
  empirical_phys_ranges[ ch ] = empirical_phys_ranges_orig[ ch ];
  
  // clear up storage
  sigmap_f[ ch ].resize(0);
  
} 

void segsrv_t::clear_filters()
{
  std::set<std::string> c = filtered;
  std::set<std::string>::const_iterator cc = c.begin();
  while ( cc != c.end() )
    {
      clear_filter( *cc );
      ++cc;
    }    
} 




template<typename T>
axis_stats_t<T> segsrv_t::compute_axis_stats(const T* x,
                                             std::size_t n,
                                             std::size_t max_unique,
                                             std::size_t max_sample)
{
    axis_stats_t<T> out{};
    out.is_discrete = true;
    out.p5  = 0.0;
    out.p95 = 0.0;

    if (n == 0) {
        out.min_val = T{};
        out.max_val = T{};
        return out;
    }

    // Track min/max over entire data
    T min_v = x[0];
    T max_v = x[0];

    // Track unique values (up to max_unique+1)
    std::vector<T> uniques;
    uniques.reserve(max_unique + 1);

    // Reservoir sample (in double) for approximate quantiles
    const std::size_t sample_cap = std::min(max_sample, n);
    std::vector<double> sample;
    sample.reserve(sample_cap);

    // Simple RNG for reservoir sampling
    // (You can move this out / make it deterministic if needed)
    static thread_local std::mt19937 rng{std::random_device{}()};

    std::size_t seen = 0;

    for (std::size_t i = 0; i < n; ++i) {
        T v = x[i];

        // min/max
        if (v < min_v) min_v = v;
        if (v > max_v) max_v = v;

        // Unique tracking, only while we still might be discrete
        if (out.is_discrete) {
            bool found = false;
            for (T u : uniques) {
                if (v == u) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                uniques.push_back(v);
                if (uniques.size() > max_unique) {
                    out.is_discrete = false;
                    uniques.clear();  // no longer needed
                }
            }
        }

        // Reservoir sampling for percentiles (if requested)
        if (sample_cap > 0) {
            const double dv = static_cast<double>(v);

            if (seen < sample_cap) {
                // Fill the reservoir
                sample.push_back(dv);
            } else {
                // Replace elements with decreasing probability
                std::uniform_int_distribution<std::size_t> dist(0, seen);
                std::size_t j = dist(rng);
                if (j < sample_cap) {
                    sample[j] = dv;
                }
            }
            ++seen;
        }
    }

    out.min_val = min_v;
    out.max_val = max_v;

    // Discrete signal: only care about uniques/min/max; p5/p95 unused
    if (out.is_discrete) {
        out.uniques = std::move(uniques);
        return out;
    }

    // Non-discrete: approximate 5th/95th percentiles from the reservoir sample
    if (sample.empty()) {
        out.p5 = out.p95 = 0.0;
        return out;
    }

    const std::size_t m = sample.size();

    auto q_index = [m](double q) -> std::size_t {
        if (m == 1) return 0;
        double idx = q * double(m - 1);
        std::size_t k = static_cast<std::size_t>(idx);
        if (k >= m) k = m - 1;
        return k;
    };

    // 5th percentile
    std::size_t k5 = q_index(0.05);
    std::nth_element(sample.begin(), sample.begin() + k5, sample.end());
    out.p5 = sample[k5];

    // 95th percentile
    std::size_t k95 = q_index(0.95);
    std::nth_element(sample.begin(), sample.begin() + k95, sample.end());
    out.p95 = sample[k95];

    // Debug if needed:
    // std::cout << "k5, k95 = " << k5 << " " << k95 << "\n";

    return out;
}


// template<typename T>
// axis_stats_t<T> segsrv_t::compute_axis_stats(const T* x,
// 					     std::size_t n,
// 					     std::size_t max_unique ,
// 					     std::size_t max_sample )
// {
//   axis_stats_t<T> out{};
//   out.is_discrete = true;
//   out.p5 = 0.0;
//   out.p95 = 0.0;

//   if (n == 0) {
//     // Nothing to do; min/max are left default-initialized
//     out.min_val = T{};
//     out.max_val = T{};
//     return out;
//   }
  
//   // Track min/max over entire data
//   T min_v = x[0];
//   T max_v = x[0];
  
//   // Track unique values (up to max_unique+1)
//   std::vector<T> uniques;
//   uniques.reserve(max_unique + 1);
  
//   // Build a sample (in double) for approximate quantiles
//   const std::size_t sample_size = std::min(max_sample, n);
//   std::vector<double> sample;
//   sample.reserve(sample_size);
  
//   // If n <= sample_size, just take all; else stride
//   std::size_t stride = (sample_size > 0 && sample_size < n)
//     ? n / sample_size
//     : 1;
//   if (stride == 0) stride = 1;
  
//   std::size_t next_sample_idx = 0;
  
//   for (std::size_t i = 0; i < n; ++i) {
//     T v = x[i];
    
//     // min/max
//     if (v < min_v) min_v = v;
//     if (v > max_v) max_v = v;
    
//     // Unique tracking, only while we still might be discrete
//     if (out.is_discrete) {
//       bool found = false;
//       for (T u : uniques) {
// 	if (v == u) {  // equality is fine for "discrete" signals
// 	  found = true;
// 	  break;
// 	}
//       }
//       if (!found) {
// 	uniques.push_back(v);
// 	if (uniques.size() > max_unique) {
// 	  out.is_discrete = false;
// 	  uniques.clear();  // no longer needed
// 	}
//       }
//     }
    
//     // Sampling for percentiles
//     if (stride == 1) {
//       // small vector: just take all
//       sample.push_back(static_cast<double>(v));
//     } else {
//       if (i == next_sample_idx && sample.size() < sample_size) {
// 	sample.push_back(static_cast<double>(v));
// 	next_sample_idx += stride;
//       }
//     }
//   }
  
//   out.min_val = min_v;
//   out.max_val = max_v;
  
//   if (out.is_discrete) {
//     // Discrete signal: we only care about uniques + min/max
//     out.uniques = std::move(uniques);
//     // p5/p95 left as 0.0 (unused in this case)
//     return out;
//   }
  
//   // Non-discrete: approximate 5th/95th percentiles from the sample
//   if (sample.empty()) {
//     out.p5 = out.p95 = 0.0;
//     return out;
//   }
  
//   std::size_t m = sample.size();
//   auto q_index = [m](double q) -> std::size_t {
//     if (m == 1) return 0;
//     double idx = q * double(m - 1);
//     std::size_t k = static_cast<std::size_t>(idx);
//     if (k >= m) k = m - 1;
//     return k;
//   };
  
//   // 5th percentile
//   std::size_t k5 = q_index(0.05);
//   std::nth_element(sample.begin(), sample.begin() + k5, sample.end());
//   out.p5 = sample[k5];
  
//   // 95th percentile
//   std::size_t k95 = q_index(0.95);
//   std::nth_element(sample.begin(), sample.begin() + k95, sample.end());
//   out.p95 = sample[k95];

//   std::cout << " k5, k95  = " << k5 << " " << k95 << "\n";
//   return out;
// }



Eigen::VectorXf segsrv_t::envelope_timetrack( const Eigen::VectorXf & x , const int nx ) const
{

  // t1, t2, t3, ...                                (per original) 
  // --> t1, t1, nan, t2, t2, nan, t3, t3, nan, ... (per n bin)

  const int n0 = x.size();

  Eigen::VectorXf x_center_out = Eigen::VectorXf::Constant(nx,
							   std::numeric_limits<float>::quiet_NaN());
  
  if (n0  == 0 || nx <= 0)
    return x_center_out;
  
  Eigen::VectorXf x_min = Eigen::VectorXf::Constant(nx,
						    std::numeric_limits<float>::infinity());
  Eigen::VectorXf x_max = Eigen::VectorXf::Constant(nx,
						    -std::numeric_limits<float>::infinity());

  // fill x_min/x_max
  for (int i = 0; i < n0; ++i) {
    int px = (int)((int64_t)i * nx / n0);
    if (px < 0) px = 0;
    if (px >= nx) px = nx - 1;
    
    double xi = x[i];
    if (xi < x_min[px]) x_min[px] = xi;
    if (xi > x_max[px]) x_max[px] = xi;
  }
  
  // Convert to bin centers
  for (int px = 0; px < nx; ++px) {
    if (std::isfinite(x_min[px]) && std::isfinite(x_max[px])) {
      x_center_out[px] = 0.5 * (x_min[px] + x_max[px]);
    }
  }

  Eigen::VectorXf tt( 3 * nx );
  int c = 0;
  for (int i=0; i<nx; i++)
    {
      tt[c++] = x_center_out[i];
      tt[c++] = x_center_out[i];
      tt[c++] = std::numeric_limits<float>::quiet_NaN();
    }
  return tt;
}




Eigen::VectorXf segsrv_t::envelope_signal_iqr( const Eigen::VectorXf & y , const int nx ) const
{

  // nb. not actually IQR, but based on robust SD

  const int n0 = y.size();
  if (n0 == 0 || nx <= 0) {
    return Eigen::VectorXf();  // empty
  }
  
  // Per-bin accumulators
  Eigen::VectorXf sum   = Eigen::VectorXf::Zero(nx);
  Eigen::VectorXf sumsq = Eigen::VectorXf::Zero(nx);
  Eigen::VectorXi count = Eigen::VectorXi::Zero(nx);
  
  // Bin samples
  for (int i = 0; i < n0; ++i) {
    int px = static_cast<int>((static_cast<long long>(i) * nx) / n0);
    if (px < 0)      px = 0;
    else if (px >= nx) px = nx - 1;
    
    const float v = y[i];
    sum[px]   += v;
    sumsq[px] += v * v;
    count[px] += 1;
  }
  
  Eigen::VectorXf out(3 * nx);
  const float NaN = std::numeric_limits<float>::quiet_NaN();
  // For Normal: Q1/Q3 = μ ± 0.67448975 σ
  const float k = 0.67448975f;

  float prior_q25, prior_q75;
  
  int c = 0;
  for (int px = 0; px < nx; ++px) {

    if (count[px] == 0) {
      out[c++] = NaN;
      out[c++] = NaN;
      out[c++] = NaN;
      continue;
    }
    
    const float n    = static_cast<float>(count[px]);
    const float mu   = sum[px] / n;
    const float ex2  = sumsq[px] / n;
    float var        = ex2 - mu * mu;
    if (var < 0.0f) var = 0.0f;  // numeric guard
    const float sd   = std::sqrt(var);
    
    float q25 = mu - k * sd;
    float q75 = mu + k * sd;

    // nudge so adjacent bars adjoin
    if ( px != 0 )
      {
	if ( q25 > prior_q75 )
	  q25 = prior_q75;
	if ( q75 < prior_q25 )
	  q75 = prior_q25;
      }

    prior_q25 = q25;
    prior_q75 = q75;
    
    out[c++] = q25;
    out[c++] = q75;
    out[c++] = NaN;   // separator for pyqtgraph spikes
  }
  
  return out;
}
 

Eigen::VectorXf segsrv_t::envelope_signal( const Eigen::VectorXf & y , const int nx ) const
{
  const int n0 = y.size();

  const int spp = n0 / nx;
  if ( spp > 100 )
    return envelope_signal_iqr( y , nx );
  
  Eigen::VectorXf y_min_out = Eigen::VectorXf::Constant(nx,
					std::numeric_limits<float>::infinity());
  Eigen::VectorXf y_max_out = Eigen::VectorXf::Constant(nx,
					-std::numeric_limits<float>::infinity());
  
  if (n0 == 0 || nx <= 0)
    return Eigen::VectorXf::Constant(nx*3,std::numeric_limits<float>::quiet_NaN()) ;
  
  for (int i = 0; i < n0; ++i) {
    int px = (int)((int64_t)i * nx / n0);
    if (px < 0) px = 0;
    if (px >= nx) px = nx - 1;

    float yi = y[i];
    if (yi < y_min_out[px]) y_min_out[px] = yi;
    if (yi > y_max_out[px]) y_max_out[px] = yi;
  }
  
  Eigen::VectorXf yy( 3 * nx );
  int c = 0;
  for (int i=0; i<nx; i++)
    {
      yy[c++] = y_min_out[i];
      yy[c++] = y_max_out[i];
      yy[c++] = std::numeric_limits<float>::quiet_NaN();
    }
  return yy;

}

float segsrv_t::min_skip_nan(const Eigen::VectorXf &v) {
  float m = std::numeric_limits<float>::infinity();
  for (int i = 0; i < v.size(); ++i) {
    float x = v[i];
    if (std::isfinite(x) && x < m)
      m = x;
  }
  return m;
}

float segsrv_t::max_skip_nan(const Eigen::VectorXf &v) {
  float m = -std::numeric_limits<float>::infinity();
  for (int i = 0; i < v.size(); ++i) {
    float x = v[i];
    if (std::isfinite(x) && x > m)
      m = x;
  }
  return m;
}
