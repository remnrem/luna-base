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

#include "resp/resp_desat_link.h"

#include "param.h"
#include "helper/helper.h"
#include "helper/logger.h"
#include "db/db.h"
#include "edf/edf.h"
#include "edf/slice.h"
#include "annot/annot.h"
#include "defs/defs.h"
#include "dsp/iir.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <vector>

extern writer_t writer;
extern logger_t logger;

namespace {

struct sleep_iv_t {
  uint64_t start;
  uint64_t stop;
};

struct resp_evt_t {
  interval_t iv;
  std::string type;
  bool in_sleep;
  instance_t * inst;
};

struct arousal_evt_t {
  interval_t iv;
  std::string type;
};

struct template_t {
  bool learned;
  int n_used;
  double mean_resp_dur;
  double baseline_t;
  double nadir_t;
  double recovery_t;
  double drop;
  std::vector<double> waveform;
  std::vector<double> waveform_smooth;
};

struct segment_mapper_t {
  struct seg_t {
    interval_t iv;
    int start_idx;
    int n_samples;
  };

  std::vector<seg_t> segs;
  int total_n;

  segment_mapper_t() : total_n(0) {}

  bool locate( const uint64_t tp , const int fs , int * idx , int * seg_start , int * seg_stop ) const
  {
    for ( int i = 0 ; i < segs.size() ; ++i )
      {
        if ( tp < segs[i].iv.start || tp >= segs[i].iv.stop ) continue;
        const double rel_sec = (double)( tp - segs[i].iv.start ) / (double)globals::tp_1sec;
        int rel_idx = (int)std::llround( rel_sec * fs );
        if ( rel_idx < 0 ) rel_idx = 0;
        if ( rel_idx >= segs[i].n_samples ) rel_idx = segs[i].n_samples - 1;
        *idx = segs[i].start_idx + rel_idx;
        *seg_start = segs[i].start_idx;
        *seg_stop = segs[i].start_idx + segs[i].n_samples - 1;
        return true;
      }
    return false;
  }
};

static std::string lower( const std::string & s )
{
  std::string r = s;
  std::transform( r.begin() , r.end() , r.begin() ,
                  []( unsigned char c ) { return std::tolower( c ); } );
  return r;
}

static std::vector<std::string> resolve_resp_labels( annotation_set_t * annots ,
                                                     const std::vector<std::string> & tokens )
{
  std::vector<std::string> names = annots->names();
  std::set<std::string> resolved;

  for ( int t = 0 ; t < tokens.size() ; ++t )
    {
      const std::string tok = lower( tokens[t] );
      for ( int i = 0 ; i < names.size() ; ++i )
        if ( lower( names[i] ) == tok )
          resolved.insert( names[i] );
    }

  return std::vector<std::string>( resolved.begin() , resolved.end() );
}

static std::vector<sleep_iv_t> sleep_intervals( edf_t & edf )
{
  std::vector<sleep_iv_t> ivs;
  const std::vector<std::string> labels = { "N1" , "N2" , "N3" , "N4" , "R" , "REM" };

  for ( int i = 0 ; i < labels.size() ; ++i )
    {
      annot_t * a = edf.annotations->find( labels[i] );
      if ( a == NULL ) continue;
      annot_map_t evts = a->extract( edf.timeline.wholetrace() );
      for ( annot_map_t::const_iterator aa = evts.begin() ; aa != evts.end() ; ++aa )
        ivs.push_back( { aa->first.interval.start , aa->first.interval.stop } );
    }

  std::sort( ivs.begin() , ivs.end() ,
             []( const sleep_iv_t & a , const sleep_iv_t & b )
             {
               if ( a.start != b.start ) return a.start < b.start;
               return a.stop < b.stop;
             } );

  return ivs;
}

static bool in_sleep( const std::vector<sleep_iv_t> & ivs , const uint64_t tp )
{
  for ( int i = 0 ; i < ivs.size() ; ++i )
    if ( tp >= ivs[i].start && tp < ivs[i].stop )
      return true;
  return false;
}

static double median_copy( std::vector<double> v )
{
  if ( v.empty() ) return std::numeric_limits<double>::quiet_NaN();
  std::sort( v.begin() , v.end() );
  const int n = v.size();
  if ( n % 2 ) return v[ n / 2 ];
  return 0.5 * ( v[ n / 2 - 1 ] + v[ n / 2 ] );
}

static double vector_mean( const std::vector<double> & v )
{
  if ( v.empty() ) return std::numeric_limits<double>::quiet_NaN();
  double s = 0.0;
  for ( int i = 0 ; i < v.size() ; ++i ) s += v[i];
  return s / v.size();
}

static std::vector<double> smooth3_nan( const std::vector<double> & x , const int passes )
{
  std::vector<double> y = x;
  for ( int p = 0 ; p < passes ; ++p )
    {
      std::vector<double> z = y;
      for ( int i = 0 ; i < y.size() ; ++i )
        {
          double s = 0.0;
          int n = 0;
          for ( int j = i - 1 ; j <= i + 1 ; ++j )
            {
              if ( j < 0 || j >= y.size() ) continue;
              if ( std::isnan( y[j] ) ) continue;
              s += y[j];
              ++n;
            }
          z[i] = n ? s / n : std::numeric_limits<double>::quiet_NaN();
        }
      y.swap( z );
    }
  return y;
}

static int clampi( const int x , const int lo , const int hi )
{
  if ( x < lo ) return lo;
  if ( x > hi ) return hi;
  return x;
}

static int argmin_range( const std::vector<double> & x , int a , int b )
{
  int best = -1;
  double bestv = std::numeric_limits<double>::infinity();
  a = std::max( a , 0 );
  b = std::min( b , (int)x.size() - 1 );
  for ( int i = a ; i <= b ; ++i )
    if ( ! std::isnan( x[i] ) && ( best == -1 || x[i] < bestv ) )
      {
        best = i;
        bestv = x[i];
      }
  return best;
}

static int argmax_range( const std::vector<double> & x , int a , int b )
{
  int best = -1;
  double bestv = -std::numeric_limits<double>::infinity();
  a = std::max( a , 0 );
  b = std::min( b , (int)x.size() - 1 );
  for ( int i = a ; i <= b ; ++i )
    if ( ! std::isnan( x[i] ) && ( best == -1 || x[i] > bestv ) )
      {
        best = i;
        bestv = x[i];
      }
  return best;
}

static int rightmost_near_max( const std::vector<double> & x , int a , int b , const double tol )
{
  a = std::max( a , 0 );
  b = std::min( b , (int)x.size() - 1 );
  if ( a > b ) return -1;
  int imax = argmax_range( x , a , b );
  if ( imax == -1 ) return -1;
  const double vmax = x[imax];
  int best = imax;
  for ( int i = a ; i <= b ; ++i )
    if ( ! std::isnan( x[i] ) && x[i] >= vmax - tol ) best = i;
  return best;
}

static int center_of_nadir_plateau( const std::vector<double> & x , int a , int b , const double tol )
{
  a = std::max( a , 0 );
  b = std::min( b , (int)x.size() - 1 );
  if ( a > b ) return -1;
  int imin = argmin_range( x , a , b );
  if ( imin == -1 ) return -1;
  const double vmin = x[imin];
  int l = imin;
  while ( l - 1 >= a && ! std::isnan( x[l-1] ) && x[l-1] <= vmin + tol ) --l;
  int r = imin;
  while ( r + 1 <= b && ! std::isnan( x[r+1] ) && x[r+1] <= vmin + tol ) ++r;
  return ( l + r ) / 2;
}

static double trimmed_mean_copy( std::vector<double> v , const double trim_frac )
{
  if ( v.empty() ) return std::numeric_limits<double>::quiet_NaN();
  std::sort( v.begin() , v.end() );
  const int n = v.size();
  int k = (int)std::floor( trim_frac * n );
  if ( 2 * k >= n ) k = std::max( 0 , ( n - 1 ) / 2 );
  double s = 0.0;
  int m = 0;
  for ( int i = k ; i < n - k ; ++i )
    {
      s += v[i];
      ++m;
    }
  return m ? s / m : std::numeric_limits<double>::quiet_NaN();
}

static std::vector<double> robust_mean_waveform( const std::vector< std::vector<double> > & windows ,
                                                 const int win_n ,
                                                 const double trim_frac )
{
  std::vector<double> out( win_n , std::numeric_limits<double>::quiet_NaN() );
  for ( int i = 0 ; i < win_n ; ++i )
    {
      std::vector<double> vals;
      for ( int r = 0 ; r < windows.size() ; ++r )
        if ( i < windows[r].size() && ! std::isnan( windows[r][i] ) )
          vals.push_back( windows[r][i] );
      if ( ! vals.empty() ) out[i] = trimmed_mean_copy( vals , trim_frac );
    }
  return out;
}

static std::vector<bool> make_artifact_mask( const std::vector<double> & x ,
                                             const double low_art ,
                                             const double spike_th ,
                                             const int expand_n )
{
  const int n = x.size();
  std::vector<bool> bad( n , false );

  for ( int i = 0 ; i < n ; ++i )
    if ( std::isnan( x[i] ) || x[i] < low_art || x[i] > 100.0 )
      bad[i] = true;

  int prev = -1;
  for ( int i = 0 ; i < n ; ++i )
    {
      if ( bad[i] ) continue;
      if ( prev != -1 && std::fabs( x[i] - x[prev] ) > spike_th )
        {
          bad[i] = true;
          bad[prev] = true;
        }
      prev = i;
    }

  if ( expand_n > 0 )
    {
      std::vector<bool> out = bad;
      for ( int i = 0 ; i < n ; ++i )
        if ( bad[i] )
          {
            const int a = std::max( 0 , i - expand_n );
            const int b = std::min( n - 1 , i + expand_n );
            for ( int j = a ; j <= b ; ++j ) out[j] = true;
          }
      bad.swap( out );
    }

  return bad;
}

static std::vector<double> filtered_signal( edf_t & edf ,
                                            const int sig_n ,
                                            const int fs ,
                                            const double lp_cutoff ,
                                            const int lp_order ,
                                            const bool do_round )
{
  slice_t slice( edf , sig_n , edf.timeline.wholetrace() );
  const std::vector<double> * pdata = slice.pdata();
  if ( pdata == NULL ) return std::vector<double>();

  std::vector<double> x = *pdata;
  if ( x.empty() ) return x;

  if ( lp_cutoff > 0 )
    {
      iir_t filt;
      filt.init( BUTTERWORTH_LOWPASS , lp_order , (double)fs , lp_cutoff , 0.0 );
      x = filt.apply( x );
    }

  if ( do_round )
    for ( int i = 0 ; i < x.size() ; ++i )
      if ( ! std::isnan( x[i] ) ) x[i] = std::round( x[i] );

  return x;
}

static segment_mapper_t build_segment_mapper( edf_t & edf ,
                                              const int sig_n )
{
  segment_mapper_t mapper;
  std::set<interval_t> segs = edf.timeline.segments();
  int offset = 0;

  for ( std::set<interval_t>::const_iterator ss = segs.begin() ; ss != segs.end() ; ++ss )
    {
      slice_t slice( edf , sig_n , *ss );
      const std::vector<double> * pdata = slice.pdata();
      const int n = pdata == NULL ? 0 : (int)pdata->size();
      if ( n <= 0 ) continue;
      segment_mapper_t::seg_t s;
      s.iv = *ss;
      s.start_idx = offset;
      s.n_samples = n;
      mapper.segs.push_back( s );
      offset += n;
    }

  mapper.total_n = offset;
  return mapper;
}

static bool interval_valid_fraction( const std::vector<bool> & bad , int a , int b , double * frac )
{
  a = std::max( a , 0 );
  b = std::min( b , (int)bad.size() - 1 );
  if ( a > b ) { *frac = 0.0; return false; }
  int good = 0;
  for ( int i = a ; i <= b ; ++i )
    if ( ! bad[i] ) ++good;
  *frac = (double)good / (double)( b - a + 1 );
  return good > 0;
}

static template_t build_template( const std::vector<resp_evt_t> & events ,
                                  const std::vector<double> & x ,
                                  const std::vector<bool> & bad ,
                                  const segment_mapper_t & mapper ,
                                  const int fs ,
                                  const int win_pre_s ,
                                  const int win_post_s ,
                                  const int min_template_events ,
                                  const int hard_baseline_pre_s ,
                                  const int hard_nadir_min_s ,
                                  const int hard_nadir_max_s ,
                                  const int hard_recovery_max_s ,
                                  const double min_post_valid_frac )
{
  template_t tpl;
  tpl.learned = false;
  tpl.n_used = 0;
  tpl.mean_resp_dur = std::numeric_limits<double>::quiet_NaN();
  tpl.baseline_t = -7.0; // Raichel default HBBaselineTime
  tpl.nadir_t = 22.0;    // Raichel default delaySpO2
  tpl.recovery_t = 34.0; // Raichel default HBRecoveryTime
  tpl.drop = std::numeric_limits<double>::quiet_NaN();

  const int win_pre_n = win_pre_s * fs;
  const int win_post_n = win_post_s * fs;
  const int win_n = win_pre_n + win_post_n + 1;

  std::vector< std::vector<double> > windows;
  std::vector<double> durs;
  int reject_no_segment = 0;
  int reject_window_edge = 0;
  int reject_post_valid = 0;
  int reject_baseline = 0;

  for ( int e = 0 ; e < events.size() ; ++e )
    {
      int end_idx = -1, seg_start = -1, seg_stop = -1;
      if ( ! mapper.locate( events[e].iv.stop , fs , &end_idx , &seg_start , &seg_stop ) )
        {
          ++reject_no_segment;
          continue;
        }

      const int a = end_idx - win_pre_n;
      const int b = end_idx + win_post_n;
      if ( a < seg_start || b > seg_stop )
        {
          ++reject_window_edge;
          continue;
        }

      double frac_post = 0.0;
      if ( ! interval_valid_fraction( bad ,
                                      end_idx + hard_nadir_min_s * fs ,
                                      end_idx + hard_nadir_max_s * fs ,
                                      &frac_post ) )
        {
          ++reject_post_valid;
          continue;
        }
      if ( frac_post < min_post_valid_frac )
        {
          ++reject_post_valid;
          continue;
        }

      std::vector<double> baseline_vals;
      const int bs_a = std::max( a , end_idx - hard_baseline_pre_s * fs );
      const int bs_b = std::min( b , end_idx );
      for ( int i = bs_a ; i <= bs_b ; ++i )
        if ( ! bad[i] ) baseline_vals.push_back( x[i] );
      if ( baseline_vals.size() < fs * 5 )
        {
          ++reject_baseline;
          continue;
        }

      const double baseline0 = median_copy( baseline_vals );
      if ( std::isnan( baseline0 ) )
        {
          ++reject_baseline;
          continue;
        }

      std::vector<double> w( win_n , std::numeric_limits<double>::quiet_NaN() );
      for ( int i = 0 ; i < win_n ; ++i )
        {
          const int idx = a + i;
          if ( idx < 0 || idx >= x.size() || bad[idx] ) continue;
          w[i] = x[idx] - baseline0;
        }

      windows.push_back( w );
      durs.push_back( events[e].iv.duration_sec() );
    }

  tpl.n_used = windows.size();
  logger << "    template rejects   : no-seg=" << reject_no_segment
         << " edge=" << reject_window_edge
         << " post-valid=" << reject_post_valid
         << " baseline=" << reject_baseline << "\n";
  if ( tpl.n_used < min_template_events ) return tpl;

  tpl.learned = true;
  tpl.mean_resp_dur = vector_mean( durs );
  const double trim_frac = windows.size() >= 10 ? 0.10 : 0.0;
  tpl.waveform = robust_mean_waveform( windows , win_n , trim_frac );
  tpl.waveform_smooth = smooth3_nan( tpl.waveform , fs <= 2 ? 2 : 1 );

  const int zero_idx = win_pre_n;
  const int nadir_a = zero_idx + hard_nadir_min_s * fs;
  const int nadir_b = zero_idx + hard_nadir_max_s * fs;
  const int nadir_i0 = argmin_range( tpl.waveform_smooth , nadir_a , nadir_b );
  if ( nadir_i0 == -1 )
    {
      tpl.learned = false;
      return tpl;
    }

  const int base_peak_i0 = argmax_range( tpl.waveform_smooth ,
                                         zero_idx - hard_baseline_pre_s * fs ,
                                         zero_idx );
  if ( base_peak_i0 == -1 )
    {
      tpl.learned = false;
      return tpl;
    }

  const double amp0 = tpl.waveform_smooth[ base_peak_i0 ] - tpl.waveform_smooth[ nadir_i0 ];
  const double plat_tol = std::max( 0.20 , 0.20 * std::fabs( amp0 ) );
  const int nadir_i = center_of_nadir_plateau( tpl.waveform_smooth , nadir_a , nadir_b , plat_tol );
  if ( nadir_i == -1 )
    {
      tpl.learned = false;
      return tpl;
    }

  int base_a = zero_idx - hard_baseline_pre_s * fs;
  int base_b = zero_idx;
  if ( ! std::isnan( tpl.mean_resp_dur ) )
    base_a = std::max( base_a , nadir_i - (int)std::llround( 2.0 * tpl.mean_resp_dur * fs ) );
  const int base_i = rightmost_near_max( tpl.waveform_smooth , base_a , base_b , plat_tol );
  if ( base_i == -1 )
    {
      tpl.learned = false;
      return tpl;
    }

  const int rec_a = nadir_i;
  const int rec_b = std::min( win_n - 1 , zero_idx + hard_recovery_max_s * fs );
  int rec_i = -1;
  const double rec_thresh = tpl.waveform_smooth[base_i] - std::max( 0.20 , 0.25 * std::fabs( amp0 ) );
  for ( int i = rec_a ; i <= rec_b ; ++i )
    if ( ! std::isnan( tpl.waveform_smooth[i] ) && tpl.waveform_smooth[i] >= rec_thresh )
      {
        rec_i = i;
        break;
      }
  if ( rec_i == -1 )
    rec_i = argmax_range( tpl.waveform_smooth , rec_a , rec_b );
  if ( rec_i == -1 )
    {
      tpl.learned = false;
      return tpl;
    }

  tpl.baseline_t = (double)( base_i - zero_idx ) / fs;
  tpl.nadir_t = (double)( nadir_i - zero_idx ) / fs;
  tpl.recovery_t = (double)( rec_i - zero_idx ) / fs;
  tpl.drop = tpl.waveform_smooth[ base_i ] - tpl.waveform_smooth[ nadir_i ];

  // Raichel suspicious-window spirit: if the learned template is too small, fall back.
  if ( std::isnan( tpl.drop ) || tpl.drop < 0.75 )
    tpl.learned = false;

  return tpl;
}

static int first_recovery_idx( const std::vector<double> & x ,
                               int nadir_i ,
                               int stop_i ,
                               const double baseline ,
                               const double nadir ,
                               const double recovery_frac )
{
  if ( stop_i < nadir_i ) return -1;
  const double target = nadir + ( baseline - nadir ) * recovery_frac;
  for ( int i = nadir_i ; i <= stop_i ; ++i )
    if ( ! std::isnan( x[i] ) && x[i] >= target )
      return i;
  return -1;
}

} // namespace

resp_desat_link_t::resp_desat_link_t( edf_t & edf , param_t & param )
{
  const std::string sig_label = param.has( "sig" ) ? param.value( "sig" ) : "SpO2";
  const std::vector<std::string> resp_tokens = param.has( "resp" )
    ? param.strvector( "resp" )
    : std::vector<std::string>{ "apnea" , "hypopnea" };
  const std::vector<std::string> arousal_tokens = param.has( "arousal" )
    ? param.strvector( "arousal" )
    : std::vector<std::string>{ "arousal" };
  const double arousal_post_sec = 5.0;

  const bool sleep_only = param.yesno( "sleep" , true , true );
  const double low_art = param.has( "low" ) ? param.requires_dbl( "low" ) : 40.0;        // Raichel
  const double spike_th = param.has( "spike" ) ? param.requires_dbl( "spike" ) : 10.0;
  int neighbor_s = param.has( "neighbor" ) ? param.requires_int( "neighbor" ) : 2; // Raichel
  double lp_cutoff = param.has( "lp" ) ? param.requires_dbl( "lp" ) : 1.0;         // Raichel
  const int lp_order = param.has( "lp-order" ) ? param.requires_int( "lp-order" ) : 2;   // Raichel
  const bool do_round = ! param.has( "no-round" );

  const int template_pre_s = param.has( "template-pre" ) ? param.requires_int( "template-pre" ) : 120; // Raichel
  const int template_post_s = param.has( "template-post" ) ? param.requires_int( "template-post" ) : 120; // Raichel
  const int hard_baseline_pre_s = param.has( "baseline-pre" ) ? param.requires_int( "baseline-pre" ) : 30;
  const int hard_nadir_min_s = param.has( "nadir-min" ) ? param.requires_int( "nadir-min" ) : 5;
  const int hard_nadir_max_s = param.has( "nadir-max" ) ? param.requires_int( "nadir-max" ) : 45; // Raichel max
  const int hard_recovery_max_s = param.has( "recovery-max" ) ? param.requires_int( "recovery-max" ) : 60;
  const int min_template_events = param.has( "min-template" ) ? param.requires_int( "min-template" ) : 8; // Raichel
  const double drop_th = param.has( "drop" ) ? param.requires_dbl( "drop" ) : 2.0;
  const double recovery_frac = param.has( "recovery" ) ? param.requires_dbl( "recovery" ) : 0.75;
  const int baseline_hw_s = param.has( "baseline-hw" ) ? param.requires_int( "baseline-hw" ) : 8;
  const int nadir_hw_s = param.has( "nadir-hw" ) ? param.requires_int( "nadir-hw" ) : 12;
  const int recovery_hw_s = param.has( "recovery-hw" ) ? param.requires_int( "recovery-hw" ) : 15;
  const int gap_tol_s = param.has( "gap-tol" ) ? param.requires_int( "gap-tol" ) : 2;    // Raichel spirit
  double min_post_valid_frac = param.has( "min-post-valid" ) ? param.requires_dbl( "min-post-valid" ) : 0.70;

  logger << "\n  RESP-LINK : respiratory-event response linker\n";
  logger << "    signal             : " << sig_label << "\n";
  logger << "    respiratory tokens : ";
  for ( int i = 0 ; i < resp_tokens.size() ; ++i )
    logger << ( i ? "," : "" ) << resp_tokens[i];
  logger << "\n";
  logger << "    arousal tokens     : ";
  for ( int i = 0 ; i < arousal_tokens.size() ; ++i )
    logger << ( i ? "," : "" ) << arousal_tokens[i];
  logger << "\n";
  logger << "    arousal rule       : onset >= resp_start and onset < resp_stop + "
         << arousal_post_sec << " sec\n";
  logger << "    sleep-only         : " << ( sleep_only ? "yes" : "no" ) << "\n";
  logger << "    low artifact       : " << low_art << " %\n";
  logger << "    spike threshold    : " << spike_th << " SpO2 units\n";
  logger << "    artifact expand    : +/-" << neighbor_s << " sec\n";
  logger << "    low-pass           : " << lp_cutoff << " Hz (order " << lp_order << ")\n";
  logger << "    rounded            : " << ( do_round ? "yes" : "no" ) << "\n";
  logger << "    template window    : [" << -template_pre_s << "," << template_post_s << "] sec around event end\n";
  logger << "    hard nadir window  : [" << hard_nadir_min_s << "," << hard_nadir_max_s << "] sec after event end\n";
  logger << "    hard recovery stop : +" << hard_recovery_max_s << " sec\n";
  logger << "    min template N     : " << min_template_events << "\n";
  int sig_n = edf.header.signal( sig_label );
  if ( sig_n == -1 )
    {
      logger << "  ** could not find signal '" << sig_label << "', skipping RESP-LINK\n";
      return;
    }

  const int fs = edf.header.sampling_freq( sig_n );
  if ( fs <= 0 )
    {
      logger << "  ** signal has zero sample rate, skipping RESP-LINK\n";
      return;
    }

  const bool auto_low_sr = fs <= 2;
  if ( auto_low_sr )
    {
      if ( ! param.has( "lp" ) ) lp_cutoff = 0.0;
      if ( ! param.has( "neighbor" ) ) neighbor_s = 0;
      if ( ! param.has( "min-post-valid" ) ) min_post_valid_frac = 0.25;
    }

  if ( auto_low_sr )
    {
      logger << "    effective expand   : +/-" << neighbor_s << " sec\n";
      logger << "    effective low-pass : " << lp_cutoff << " Hz\n";
    }

  logger << "    linked drop >=     : " << drop_th << " %\n";
  logger << "    min post-valid     : " << min_post_valid_frac << "\n";
  logger << "    link annotation    : existing respiratory events\n";

  std::vector<double> x = filtered_signal( edf , sig_n , fs , lp_cutoff , lp_order , do_round );
  if ( x.empty() )
    {
      logger << "  ** could not read signal data, skipping RESP-LINK\n";
      return;
    }

  const segment_mapper_t mapper = build_segment_mapper( edf , sig_n );
  if ( mapper.total_n != (int)x.size() )
    logger << "  ** warning: segment/sample map mismatch " << mapper.total_n
           << " vs " << x.size() << "\n";

  const std::vector<bool> bad = make_artifact_mask( x , low_art , spike_th , neighbor_s * fs );
  std::vector<double> x_search = x;
  for ( int i = 0 ; i < x_search.size() ; ++i )
    if ( bad[i] ) x_search[i] = std::numeric_limits<double>::quiet_NaN();
  x_search = smooth3_nan( x_search , auto_low_sr ? 2 : 1 );
  const std::vector<sleep_iv_t> slp = sleep_intervals( edf );
  const bool has_staging = ! slp.empty();

  std::vector<std::string> event_labels = resolve_resp_labels( edf.annotations , resp_tokens );
  if ( event_labels.empty() )
    {
      logger << "  ** no annotation labels matched resp= list, skipping RESP-LINK\n";
      return;
    }

  std::vector<std::string> arousal_labels = resolve_resp_labels( edf.annotations , arousal_tokens );

  logger << "    matched labels     : ";
  for ( int i = 0 ; i < event_labels.size() ; ++i )
    logger << ( i ? "," : "" ) << event_labels[i];
  logger << "\n";
  logger << "    matched arousals   : ";
  if ( arousal_labels.empty() )
    logger << "(none)";
  else
    for ( int i = 0 ; i < arousal_labels.size() ; ++i )
      logger << ( i ? "," : "" ) << arousal_labels[i];
  logger << "\n";
  logger << "    staging            : " << ( has_staging ? "found" : "not found" ) << "\n";
  if ( auto_low_sr )
    logger << "    low-SR mode        : adaptive defaults enabled\n";

  std::map<instance_idx_t,instance_t*> events_map;
  for ( int e = 0 ; e < event_labels.size() ; ++e )
    {
      annot_t * annot = edf.annotations->find( event_labels[e] );
      if ( annot == NULL ) continue;
      annot_map_t evts = annot->extract( edf.timeline.wholetrace() );
      for ( annot_map_t::const_iterator aa = evts.begin() ; aa != evts.end() ; ++aa )
        events_map.insert( *aa );
    }

  if ( events_map.empty() )
    {
      logger << "  ** matched labels had no events, skipping RESP-LINK\n";
      return;
    }

  std::map<instance_idx_t,instance_t*> arousal_map;
  for ( int a = 0 ; a < arousal_labels.size() ; ++a )
    {
      annot_t * annot = edf.annotations->find( arousal_labels[a] );
      if ( annot == NULL ) continue;
      annot_map_t evts = annot->extract( edf.timeline.wholetrace() );
      for ( annot_map_t::const_iterator aa = evts.begin() ; aa != evts.end() ; ++aa )
        arousal_map.insert( *aa );
    }

  std::vector<arousal_evt_t> arousals;
  for ( std::map<instance_idx_t,instance_t*>::const_iterator ii = arousal_map.begin() ; ii != arousal_map.end() ; ++ii )
    {
      arousal_evt_t ar;
      ar.iv = ii->first.interval;
      ar.type = ii->first.parent->name;
      arousals.push_back( ar );
    }

  std::vector<resp_evt_t> events;
  for ( std::map<instance_idx_t,instance_t*>::const_iterator ii = events_map.begin() ; ii != events_map.end() ; ++ii )
    {
      const interval_t & iv = ii->first.interval;
      const uint64_t mid = iv.start + ( iv.stop - iv.start ) / 2;
      const bool sl = has_staging ? in_sleep( slp , mid ) : true;
      if ( sleep_only && has_staging && ! sl ) continue;
      resp_evt_t ev;
      ev.iv = iv;
      ev.type = ii->first.parent->name;
      ev.in_sleep = sl;
      ev.inst = ii->second;
      events.push_back( ev );
    }

  if ( events.empty() )
    {
      logger << "  ** no events remained after sleep filtering, skipping RESP-LINK\n";
      return;
    }

  logger << "    selected events    : " << events.size() << "\n";
  logger << "    arousal events     : " << arousals.size() << "\n";

  const template_t tpl = build_template( events , x , bad , mapper , fs ,
                                         template_pre_s , template_post_s ,
                                         min_template_events ,
                                         hard_baseline_pre_s ,
                                         hard_nadir_min_s ,
                                         hard_nadir_max_s ,
                                         hard_recovery_max_s ,
                                         min_post_valid_frac );

  const bool use_learned = tpl.learned;
  const double mean_resp_dur = use_learned && ! std::isnan( tpl.mean_resp_dur )
    ? tpl.mean_resp_dur
    : 0.0;

  logger << "    template mode      : " << ( use_learned ? "learned" : "fallback" ) << "\n";
  logger << "    template N         : " << tpl.n_used << "\n";
  logger << "    template baseline  : " << tpl.baseline_t << " sec\n";
  logger << "    template nadir     : " << tpl.nadir_t << " sec\n";
  logger << "    template recovery  : " << tpl.recovery_t << " sec\n";
  if ( ! std::isnan( tpl.drop ) )
    logger << "    template drop      : " << tpl.drop << " %\n";

  int n_linked = 0;
  double sum_drop = 0.0;
  double sum_delay = 0.0;
  std::map<std::string,int> reject_counts;
  std::vector<resp_desat_link_event_t> out;

  for ( int e = 0 ; e < events.size() ; ++e )
    {
      const double resp_start_sec = events[e].iv.start / (double)globals::tp_1sec;
      const double resp_stop_sec = events[e].iv.stop / (double)globals::tp_1sec;
      const double resp_dur_sec = events[e].iv.duration_sec();

      resp_desat_link_event_t row;
      row.resp_start_sec = resp_start_sec;
      row.resp_stop_sec = resp_stop_sec;
      row.resp_dur_sec = resp_dur_sec;
      row.resp_type = events[e].type;
      row.in_sleep = events[e].in_sleep;
      row.linked = false;
      row.reject_reason = "UNSET";
      row.arousal_linked = false;
      row.arousal_reject_reason = arousals.empty() ? "NO_AROUSAL_ANNOT" : "NO_AROUSAL";
      row.arousal_type = ".";
      row.arousal_start_sec = std::numeric_limits<double>::quiet_NaN();
      row.arousal_stop_sec = std::numeric_limits<double>::quiet_NaN();
      row.arousal_dur_sec = std::numeric_limits<double>::quiet_NaN();
      row.arousal_onset_rel_start_sec = std::numeric_limits<double>::quiet_NaN();
      row.arousal_onset_rel_stop_sec = std::numeric_limits<double>::quiet_NaN();
      row.baseline_center_t = std::numeric_limits<double>::quiet_NaN();
      row.nadir_center_t = std::numeric_limits<double>::quiet_NaN();
      row.recovery_center_t = std::numeric_limits<double>::quiet_NaN();
      row.prev_rel_stop_sec = std::numeric_limits<double>::quiet_NaN();
      row.next_rel_start_sec = std::numeric_limits<double>::quiet_NaN();
      row.base_a_sec = row.base_b_sec = std::numeric_limits<double>::quiet_NaN();
      row.nad_a_sec = row.nad_b_sec = std::numeric_limits<double>::quiet_NaN();
      row.baseline_sec = row.baseline_val = std::numeric_limits<double>::quiet_NaN();
      row.nadir_sec = row.nadir_val = std::numeric_limits<double>::quiet_NaN();
      row.recovery_sec = row.recovery_val = std::numeric_limits<double>::quiet_NaN();
      row.drop = std::numeric_limits<double>::quiet_NaN();

      const double arousal_latest_sec = resp_stop_sec + arousal_post_sec;
      for ( int a = 0 ; a < arousals.size() ; ++a )
        {
          const double arousal_start_sec = arousals[a].iv.start / (double)globals::tp_1sec;
          if ( arousal_start_sec < resp_start_sec ) continue;
          if ( arousal_start_sec >= arousal_latest_sec ) break;
          row.arousal_linked = true;
          row.arousal_reject_reason = ".";
          row.arousal_type = arousals[a].type;
          row.arousal_start_sec = arousal_start_sec;
          row.arousal_stop_sec = arousals[a].iv.stop / (double)globals::tp_1sec;
          row.arousal_dur_sec = arousals[a].iv.duration_sec();
          row.arousal_onset_rel_start_sec = row.arousal_start_sec - resp_start_sec;
          row.arousal_onset_rel_stop_sec = row.arousal_start_sec - resp_stop_sec;
          break;
        }

      if ( events[e].inst != NULL )
        {
          events[e].inst->set( "arousal_linked" , row.arousal_linked ? 1 : 0 );
          events[e].inst->set( "arousal_reject" , row.arousal_reject_reason );
          if ( row.arousal_linked )
            {
              events[e].inst->set( "arousal_type" , row.arousal_type );
              events[e].inst->set( "arousal_start_t" , row.arousal_start_sec - row.resp_stop_sec );
              events[e].inst->set( "arousal_start" , row.arousal_start_sec );
              events[e].inst->set( "arousal_stop" , row.arousal_stop_sec );
              events[e].inst->set( "arousal_dur" , row.arousal_dur_sec );
              events[e].inst->set( "arousal_onset_from_start" , row.arousal_onset_rel_start_sec );
            }
        }

      int end_idx = -1, seg_start = -1, seg_stop = -1;
      if ( ! mapper.locate( events[e].iv.stop , fs , &end_idx , &seg_start , &seg_stop ) )
        {
          row.reject_reason = "NO_SEG";
          ++reject_counts[ row.reject_reason ];
          if ( events[e].inst != NULL )
            {
              events[e].inst->set( "desat_linked" , 0 );
              events[e].inst->set( "desat_reject" , row.reject_reason );
            }
          out.push_back( row );
          continue;
        }

      const double d_dur = resp_dur_sec - mean_resp_dur;
      const double baseline_center_t = tpl.baseline_t - d_dur;
      const double nadir_center_t = tpl.nadir_t;
      const double recovery_center_t = tpl.recovery_t;
      row.baseline_center_t = baseline_center_t;
      row.nadir_center_t = nadir_center_t;
      row.recovery_center_t = recovery_center_t;

      int prev_rel_stop = -hard_baseline_pre_s * fs;
      if ( e > 0 )
        {
          const int64_t dt_tp = (int64_t)events[e-1].iv.stop - (int64_t)events[e].iv.stop;
          prev_rel_stop = (int)std::llround( (double)dt_tp / (double)globals::tp_1sec * fs );
        }
      row.prev_rel_stop_sec = (double)prev_rel_stop / fs;

      int next_rel_start = hard_recovery_max_s * fs;
      if ( e + 1 < events.size() )
        {
          const int64_t dt_tp = (int64_t)events[e+1].iv.start - (int64_t)events[e].iv.stop;
          next_rel_start = (int)std::llround( (double)dt_tp / (double)globals::tp_1sec * fs );
        }
      row.next_rel_start_sec = (double)next_rel_start / fs;

      int base_a = end_idx + clampi( (int)std::llround( ( baseline_center_t - baseline_hw_s ) * fs ) ,
                                     -hard_baseline_pre_s * fs ,
                                     5 * fs );
      int base_b = end_idx + clampi( (int)std::llround( ( baseline_center_t + baseline_hw_s ) * fs ) ,
                                     -hard_baseline_pre_s * fs ,
                                     5 * fs );
      base_a = std::max( base_a , end_idx + prev_rel_stop + gap_tol_s * fs );
      base_b = std::min( base_b , end_idx + hard_nadir_min_s * fs );

      int nad_a = end_idx + clampi( (int)std::llround( ( nadir_center_t - nadir_hw_s ) * fs ) ,
                                    hard_nadir_min_s * fs ,
                                    hard_nadir_max_s * fs );
      int nad_b = end_idx + clampi( (int)std::llround( ( nadir_center_t + nadir_hw_s ) * fs ) ,
                                    hard_nadir_min_s * fs ,
                                    hard_nadir_max_s * fs );
      nad_b = std::min( nad_b , end_idx + next_rel_start + gap_tol_s * fs );

      if ( base_a < seg_start ) base_a = seg_start;
      if ( base_b > seg_stop ) base_b = seg_stop;
      if ( nad_a < seg_start ) nad_a = seg_start;
      if ( nad_b > seg_stop ) nad_b = seg_stop;
      row.base_a_sec = (double)base_a / fs;
      row.base_b_sec = (double)base_b / fs;
      row.nad_a_sec = (double)nad_a / fs;
      row.nad_b_sec = (double)nad_b / fs;

      if ( base_a > base_b )
        {
          row.reject_reason = "BAD_BASE_WINDOW";
          ++reject_counts[ row.reject_reason ];
          if ( events[e].inst != NULL )
            {
              events[e].inst->set( "desat_linked" , 0 );
              events[e].inst->set( "desat_reject" , row.reject_reason );
            }
          out.push_back( row );
          continue;
        }

      if ( nad_a > nad_b )
        {
          row.reject_reason = "BAD_NADIR_WINDOW";
          ++reject_counts[ row.reject_reason ];
          if ( events[e].inst != NULL )
            {
              events[e].inst->set( "desat_linked" , 0 );
              events[e].inst->set( "desat_reject" , row.reject_reason );
            }
          out.push_back( row );
          continue;
        }

      const int base_i = argmax_range( x_search , base_a , base_b );
      if ( base_i == -1 )
        {
          row.reject_reason = "NO_BASELINE";
          ++reject_counts[ row.reject_reason ];
          if ( events[e].inst != NULL )
            {
              events[e].inst->set( "desat_linked" , 0 );
              events[e].inst->set( "desat_reject" , row.reject_reason );
            }
          out.push_back( row );
          continue;
        }

      row.baseline_sec = (double)base_i / fs;
      row.baseline_val = x[base_i];

      const int nad_i = argmin_range( x_search , std::max( nad_a , base_i ) , nad_b );

      if ( nad_i == -1 )
        {
          row.reject_reason = "NO_NADIR";
          ++reject_counts[ row.reject_reason ];
          if ( events[e].inst != NULL )
            {
              events[e].inst->set( "desat_baseline_t" , row.baseline_sec - row.resp_stop_sec );
              events[e].inst->set( "desat_baseline" , row.baseline_val );
              events[e].inst->set( "desat_linked" , 0 );
              events[e].inst->set( "desat_reject" , row.reject_reason );
            }
          out.push_back( row );
          continue;
        }

      row.nadir_sec = (double)nad_i / fs;
      row.nadir_val = x[nad_i];
      row.drop = row.baseline_val - row.nadir_val;

      int rec_a = nad_i;
      int rec_b = end_idx + clampi( (int)std::llround( ( recovery_center_t + recovery_hw_s ) * fs ) ,
                                    hard_nadir_min_s * fs ,
                                    hard_recovery_max_s * fs );
      rec_b = std::min( rec_b , end_idx + next_rel_start + gap_tol_s * fs );
      rec_b = std::min( rec_b , seg_stop );

      int rec_i = first_recovery_idx( x_search , rec_a , rec_b , x_search[base_i] , x_search[nad_i] , recovery_frac );
      if ( rec_i == -1 )
        rec_i = argmax_range( x_search , rec_a , rec_b );

      if ( rec_i == -1 )
        {
          row.reject_reason = "NO_RECOVERY";
          ++reject_counts[ row.reject_reason ];
          if ( events[e].inst != NULL )
            {
              events[e].inst->set( "desat_baseline_t" , row.baseline_sec - row.resp_stop_sec );
              events[e].inst->set( "desat_baseline" , row.baseline_val );
              events[e].inst->set( "desat_nadir_t" , row.nadir_sec - row.resp_stop_sec );
              events[e].inst->set( "desat_nadir" , row.nadir_val );
              events[e].inst->set( "desat_drop" , row.drop );
              events[e].inst->set( "desat_delay" , row.nadir_sec - row.resp_stop_sec );
              events[e].inst->set( "desat_linked" , 0 );
              events[e].inst->set( "desat_reject" , row.reject_reason );
            }
          out.push_back( row );
          continue;
        }

      bool bad_between = false;
      for ( int i = base_i ; i <= rec_i ; ++i )
        if ( bad[i] ) { bad_between = true; break; }

      row.recovery_sec = (double)rec_i / fs;
      row.recovery_val = x[rec_i];

      if ( events[e].inst != NULL )
        {
          events[e].inst->set( "desat_baseline_t" , row.baseline_sec - row.resp_stop_sec );
          events[e].inst->set( "desat_baseline" , row.baseline_val );
          events[e].inst->set( "desat_nadir_t" , row.nadir_sec - row.resp_stop_sec );
          events[e].inst->set( "desat_nadir" , row.nadir_val );
          events[e].inst->set( "desat_drop" , row.drop );
          events[e].inst->set( "desat_delay" , row.nadir_sec - row.resp_stop_sec );
          events[e].inst->set( "desat_recovery_t" , row.recovery_sec - row.resp_stop_sec );
          events[e].inst->set( "desat_recovery" , row.recovery_val );
        }

      if ( ! bad_between && row.drop >= drop_th && row.nadir_sec > resp_stop_sec )
        {
          row.linked = true;
          row.reject_reason = ".";
          ++n_linked;
          sum_drop += row.drop;
          sum_delay += row.nadir_sec - resp_stop_sec;

          if ( events[e].inst != NULL )
            {
              events[e].inst->set( "desat_linked" , 1 );
              events[e].inst->set( "desat_reject" , row.reject_reason );
            }
        }
      else if ( events[e].inst != NULL )
        {
          if ( bad_between )
            row.reject_reason = "ARTIFACT_BETWEEN";
          else if ( row.drop < drop_th )
            row.reject_reason = "DROP_LT_TH";
          else if ( row.nadir_sec <= resp_stop_sec )
            row.reject_reason = "NADIR_PRE_END";
          else
            row.reject_reason = "UNLINKED";
          ++reject_counts[ row.reject_reason ];
          events[e].inst->set( "desat_linked" , 0 );
          events[e].inst->set( "desat_reject" , row.reject_reason );
        }

      out.push_back( row );
    }

  writer.value( "N_RESP" , (int)out.size() );
  writer.value( "N_TEMPLATE" , tpl.n_used );
  writer.value( "N_LINKED" , n_linked );
  writer.value( "PCT_LINKED" , out.empty() ? 0.0 : 100.0 * n_linked / (double)out.size() );
  writer.value( "TEMPLATE_LEARNED" , (int)use_learned );
  writer.value( "TEMPLATE_BASE_T" , tpl.baseline_t );
  writer.value( "TEMPLATE_NADIR_T" , tpl.nadir_t );
  writer.value( "TEMPLATE_REC_T" , tpl.recovery_t );
  if ( ! std::isnan( tpl.drop ) ) writer.value( "TEMPLATE_DROP" , tpl.drop );
  writer.value( "DROP_MEAN" , n_linked ? sum_drop / n_linked : 0.0 );
  writer.value( "DELAY_MEAN" , n_linked ? sum_delay / n_linked : 0.0 );
  for ( const auto & kv : reject_counts )
    writer.value( "REJECT_" + kv.first , kv.second );

  if ( ! tpl.waveform.empty() )
    {
      const int zero_idx = template_pre_s * fs;
      for ( int i = 0 ; i < tpl.waveform.size() ; ++i )
        {
          writer.level( i + 1 , "TPL" );
          writer.value( "T" , (double)( i - zero_idx ) / fs );
          if ( ! std::isnan( tpl.waveform[i] ) ) writer.value( "D_SPO2" , tpl.waveform[i] );
          if ( i < tpl.waveform_smooth.size() && ! std::isnan( tpl.waveform_smooth[i] ) ) writer.value( "D_SPO2_S" , tpl.waveform_smooth[i] );
          writer.unlevel( "TPL" );
        }
    }

  for ( int i = 0 ; i < out.size() ; ++i )
    {
      writer.level( i + 1 , "RDL" );
      writer.value( "RESP_START" , out[i].resp_start_sec );
      writer.value( "RESP_STOP" , out[i].resp_stop_sec );
      writer.value( "RESP_DUR" , out[i].resp_dur_sec );
      writer.value( "TYPE" , out[i].resp_type );
      writer.value( "SLEEP" , (int)out[i].in_sleep );
      writer.value( "LINKED" , (int)out[i].linked );
      writer.value( "REJECT" , out[i].reject_reason );
      writer.value( "AROUSAL_LINKED" , (int)out[i].arousal_linked );
      writer.value( "AROUSAL_REJECT" , out[i].arousal_reject_reason );
      writer.value( "AROUSAL_TYPE" , out[i].arousal_type );
      if ( ! std::isnan( out[i].arousal_start_sec ) ) writer.value( "AROUSAL_START" , out[i].arousal_start_sec );
      if ( ! std::isnan( out[i].arousal_stop_sec ) ) writer.value( "AROUSAL_STOP" , out[i].arousal_stop_sec );
      if ( ! std::isnan( out[i].arousal_dur_sec ) ) writer.value( "AROUSAL_DUR" , out[i].arousal_dur_sec );
      if ( ! std::isnan( out[i].arousal_onset_rel_start_sec ) ) writer.value( "AROUSAL_ONSET_FROM_START" , out[i].arousal_onset_rel_start_sec );
      if ( ! std::isnan( out[i].arousal_onset_rel_stop_sec ) ) writer.value( "AROUSAL_ONSET_FROM_STOP" , out[i].arousal_onset_rel_stop_sec );
      if ( ! std::isnan( out[i].baseline_center_t ) ) writer.value( "BASE_CENTER_T" , out[i].baseline_center_t );
      if ( ! std::isnan( out[i].nadir_center_t ) ) writer.value( "NADIR_CENTER_T" , out[i].nadir_center_t );
      if ( ! std::isnan( out[i].recovery_center_t ) ) writer.value( "REC_CENTER_T" , out[i].recovery_center_t );
      if ( ! std::isnan( out[i].prev_rel_stop_sec ) ) writer.value( "PREV_REL_STOP" , out[i].prev_rel_stop_sec );
      if ( ! std::isnan( out[i].next_rel_start_sec ) ) writer.value( "NEXT_REL_START" , out[i].next_rel_start_sec );
      if ( ! std::isnan( out[i].base_a_sec ) ) writer.value( "BASE_A" , out[i].base_a_sec );
      if ( ! std::isnan( out[i].base_b_sec ) ) writer.value( "BASE_B" , out[i].base_b_sec );
      if ( ! std::isnan( out[i].nad_a_sec ) ) writer.value( "NAD_A" , out[i].nad_a_sec );
      if ( ! std::isnan( out[i].nad_b_sec ) ) writer.value( "NAD_B" , out[i].nad_b_sec );
      if ( ! std::isnan( out[i].baseline_sec ) ) writer.value( "BASE_T" , out[i].baseline_sec );
      if ( ! std::isnan( out[i].baseline_val ) ) writer.value( "BASELINE" , out[i].baseline_val );
      if ( ! std::isnan( out[i].nadir_sec ) ) writer.value( "NADIR_T" , out[i].nadir_sec );
      if ( ! std::isnan( out[i].nadir_val ) ) writer.value( "NADIR" , out[i].nadir_val );
      if ( ! std::isnan( out[i].recovery_sec ) ) writer.value( "REC_T" , out[i].recovery_sec );
      if ( ! std::isnan( out[i].recovery_val ) ) writer.value( "RECOVERY" , out[i].recovery_val );
      if ( ! std::isnan( out[i].drop ) ) writer.value( "DROP" , out[i].drop );
      writer.unlevel( "RDL" );
    }
}
