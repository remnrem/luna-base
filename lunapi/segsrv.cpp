
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

segsrv_t::segsrv_t( lunapi_inst_ptr inst ) : p( inst ) 
{
  awin = bwin = 0;
  
  max_samples_in = 200; // default input throttling
  max_samples_out = 0;  // no output throttling by default

  epoch_sec = 30;
  
}

void segsrv_t::init()
{

  awin = bwin = 0;
  aidx.clear();
  bidx.clear();
  evts.clear();
  segments = p->edf.timeline.segments();
  gaps = p->edf.timeline.gaps( segments ); 

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
	      std::cout << " data->size()  = " << data->size()  << "\n"
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
	      MiscMath::hjorth( &edata , &activity , &mobility , &complexity );
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

  // make sure we have some sensible default scaling
  
  set_scaling( count ,
	       anns.size() ,
	       1 , // yscale
	       1 , // ygroup
	       0.05, // yheader
	       0.05, // yfooter
	       0.10 ); // scaling_fixed_annot
  
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
  
  // new sample rate
  const int sr2 = sr / q ;
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

  // do we already have a time-track?
  if ( tidx.find( sr ) == tidx.end() )
    {      
      // get time-stampls
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
      
      // lookup index - note, here ts.size() as we may have decimtaed
      std::map<double,int> tt;
      for (int i=0; i<ts.size(); i++)
	tt[ ts[i] ] = i;

      
      // store lookup index (time->sample #)
      tidx[ sr ] = tt;
      
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
      const int original_length = bb - aa;
      const int reduction_factor = original_length / max_samples_out;      
      Eigen::VectorXf f1 = tt.segment( aa , bb - aa );
      return f1( Eigen::seq( 0 , Eigen::last , reduction_factor ) );
    }

  return tt.segment( aa , bb - aa );

}


// set scaling params                                                                                                                   
void segsrv_t::set_scaling( const int nchs , const int nanns ,
			    const double yscale , const double ygroup ,
			    const double yheader , const double yfooter ,
			    const double fixed_annot )
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
  if ( scaling_fixed_annot > 0.5 ) scaling_fixed_annot = 0.5;

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

void segsrv_t::free_physical_scale( const std::string & ch )
{
  //   std::map<std::string, std::pair<double,double> > phys_ranges;
  std::map<std::string, std::pair<double,double> >::const_iterator pp = phys_ranges.find( ch );
  if ( pp == phys_ranges.end() ) return;
  // erase
  phys_ranges.erase( pp );
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
  // if fixed scaling, then clip if above/below
  
  float smin, smax; 

  std::map<std::string, std::pair<double,double> >::const_iterator pp =	phys_ranges.find( ch );
  if ( pp == phys_ranges.end() )
    {
      // auto-scaling  
      smin = s.minCoeff();
      smax = s.maxCoeff();

      // store
      window_phys_range[ ch ] = std::pair<double,double>( smin, smax );
    }
  else
    {
      // use pre-specified values
      smin = pp->second.first;
      smax = pp->second.second;
      
      // but compare to observed
      float omin = s.minCoeff();
      float omax = s.maxCoeff();
      
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
  
  const double srange = smax - smin;
  // special case of smin == smax --> set whole signal to 0.5 
  if ( srange < 1e-3 )
    s = Eigen::VectorXf::Zero( s.size() ).array() + 0.5 ;
  else // normalize to 0..1  X = ( X - min ) / ( max - min )    
    s = ( s.array() - smin ) / (float)srange;

  // rescale to display units
  double lwr = 0,  upr = 1;
  const bool okay = get_yscale_signal( n1 , &lwr, &upr );

  if ( okay )
    s = s.array() * ( upr - lwr ) + lwr;
  
  // return
  return s;
      
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

  // return data (signal)
  const Eigen::VectorXf & data = sigmap.find( ch )->second;
  const int aa = aidx.find(sr)->second;
  const int bb = bidx.find(sr)->second;  

  // throttle?
  if ( max_samples_out && ( bb - aa ) > max_samples_out )
    {
      const int original_length = bb - aa;
      const int reduction_factor = original_length / max_samples_out;

      // decimation here on second level - but note that the SR may
      // have been adjusted from initial decimation; note also that SR
      // is cast to an int here, but for this single purpose (of
      // designing a low-pass anti-aliasing filter) this should not
      // matter.
      
      return decimate( data.segment( aa , bb - aa ) , decimated_srmap.find( sr )->second , reduction_factor );
      
    }
  
  // else return full
  return data.segment( aa , bb - aa );
    
}


// given two times and a sample rate, get indices
bool segsrv_t::get_tidx( double a, double b , int sr , int * aa, int *bb ) const
{

  if ( tidx.find( sr ) == tidx.end() ) return false;

  const std::map<double,int> & ts = tidx.find( sr )->second ;

  // iterator equal/greater than start
  std::map<double,int>::const_iterator abound = ts.lower_bound( a );
  if ( abound == ts.end() ) return false;
 
  // one-past the end
  std::map<double,int>::const_iterator bbound = ts.lower_bound( b );
  if ( bbound == ts.end() ) return false;
  
  // if we are in a gap, then both abound and bbound will point to the
  // same element (i.e. if not end(), then both the same next segment
  // index;  this means the window is not valid

  if ( abound == bbound ) return false;
  
  *aa = abound->second;
  *bb = bbound->second;

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

  annot_t * annot = p->edf.timeline.annotations.find( ch );
  if ( annot == NULL ) return false;
  if ( annot->interval_events.size() == 0 ) return false;
  
  annot_map_t::const_iterator ii = annot->interval_events.begin();   
  while ( ii != annot->interval_events.end() )
    {
      const instance_idx_t & instance_idx = ii->first;
      
      evt_t evt( ii->first.interval , ii->first.parent->name );

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
      
  return r;
}


// for selection window
std::vector<std::string> segsrv_t::fetch_all_evts( const std::vector<std::string> & avec ) const
{
  std::vector<std::string> r;
  std::set<std::string> aset = Helper::vec2set( avec );
  
  std::set<evt_t>::const_iterator ee = evts.begin();
  while ( ee != evts.end() )
    {
      if ( aset.find( ee->name ) != aset.end() )
	r.push_back( ee->name + " | " + ee->interval.as_string(1,"-") );
      ++ee;
    }
  return r;
}


// compile a set of selected events for the current window
void segsrv_t::compile_evts( const std::vector<std::string> & anns )
{
  // clear window

  // ultimate (in 6-val format for plottin)
  compiled_annots_times.clear();
  compiled_annots_stacks.clear();

  // working intermediate
  std::map<std::string,std::vector<std::pair<double,double> > > annots_times;
  std::map<std::string,std::vector<std::pair<double,double> > > annots_stacks;

  // get all events that overlap this window
  std::map<std::string,std::vector<std::pair<double,double> > > wevts = fetch_evts();

  // make uniform time-line: extracted from in-window events
  std::set<fevt_t> xevts;
  
  for (int a=0;a<anns.size();a++)
    {
      // not present?
      if ( wevts.find( anns[a] ) == wevts.end() ) continue;

      const std::string & aname = anns[a];
      
      // events
      const std::vector<std::pair<double,double> > & e = wevts.find( aname )->second;

      for (int i=0;i<e.size();i++)
	xevts.insert( fevt_t( e[i].first , e[i].second , aname ) );
      
    }

  // empty? all done
  if ( xevts.size() == 0 ) return;

  // determine stacking
  std::set<fevt_t> pool;
  int maxd = 0;
  
  std::set<fevt_t>::const_iterator xx = xevts.begin();
  while ( xx != xevts.end() )
    {
      
      // remove any events (that end before start of this one)
      const double c = xx->start;

      std::set<fevt_t> pool1;
      std::set<fevt_t>::const_iterator pp = pool.begin();
      while ( pp != pool.end() )
	{
	  if ( pp->stop > c ) pool1.insert( *pp );
	  ++pp;
	}

      // copy over
      pool = pool1;
      
      // add next event to pool
      pool.insert( *xx );
      
      // add times (w/ clipping at window- boundaries) 
      const std::pair<double,double> p2( xx->start < awin ? awin : xx->start ,
					 xx->stop  > bwin ? bwin : xx->stop );
      
      annots_times[ xx->name ].push_back( p2 );
						    
      
      // determine depth (N-1) [ scale after when we know max depth ]  
      const int d = pool.size();
      annots_stacks[ xx->name ].push_back( std::pair<double,double>( d , 0 ) );

      // track max depth
      if ( d > maxd ) maxd = d;

      // next annot
      ++xx;
    }
  
  // scale depth into plotting values:
  //  annots go from ( 1 - scaling_fixed_annot - scaling_yheader ) to ( 1 - scaling_yheader )
  const double abase = 1.0 - scaling_fixed_annot - scaling_yheader;

  // maxd +1 i.e to allow for height of top annot
  const double a1height = scaling_fixed_annot / (double)(maxd);
  
  std::map<std::string,std::vector<std::pair<double,double> > >::iterator cc = annots_stacks.begin();
  while ( cc != annots_stacks.end() )
    {

      std::vector<std::pair<double,double> > & ds = cc->second;

      for (int i=0; i<ds.size(); i++)
	{
	  double y = ds[i].first; // scale from 1... maxd
	  y = abase + (y-1) * a1height; // lower y;
	  ds[i].first = y;
	  ds[i].second = y + a1height;
	}

      
      ++cc;
    }


  //
  // for annot 'i' we now have (x1,x2) and (y1,y2) 
  //  for easier plotting, we want to have six entires:
  //     (x1,x2,x2,x1,x1,NA)
  //     (y1,y1,y2,y2,y1,NA)
  //  i.e. for plotly to make the traces and fill in a rectange
  //
  
  cc = annots_times.begin();
  while ( cc != annots_times.end() )
    {
      
      std::vector<std::pair<double,double> > & xs = cc->second;
      std::vector<std::pair<double,double> > & ys = annots_stacks[ cc->first ];

      const int n = xs.size();

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

// given a compilation (subset of all evts), get y-axis stacking for a particular class
std::vector<float> segsrv_t::get_evnts_yaxes( const std::string & ann ) const
{
  std::map<std::string,std::vector<float> >::const_iterator aa = compiled_annots_stacks.find( ann );
  if ( aa != compiled_annots_stacks.end() ) return aa->second;
  std::vector<float> empty;
  return empty;	
}
