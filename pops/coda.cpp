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

#ifdef HAS_LGBM

#include "pops/coda.h"
#include "pops/pops.h"
#include "pops/posteriors.h"
#include "pops/options.h"
#include "lgbm/lgbm.h"
#include "miscmath/crandom.h"
#include "miscmath/miscmath.h"

#include "helper/helper.h"
#include "helper/logger.h"
#include "defs/defs.h"
#include "db/db.h"
#include "edf/edf.h"

#include <cmath>
#include <algorithm>
#include <array>
#include <map>
#include <numeric>
#include <fstream>
#include <limits>

extern logger_t logger;
extern writer_t writer;

namespace {

const char * pops_coda_soft_fail_msg = "POPS-CODA: no usable posterior rows remain after screening";

std::string coda_feature_name_file( const std::string & model_file )
{
  if ( model_file.size() >= 4 &&
       model_file.substr( model_file.size() - 4 ) == ".mod" )
    return model_file.substr( 0 , model_file.size() - 4 ) + ".fnames";
  return model_file + ".fnames";
}

bool coda_subject_in_holdouts( const std::set<std::string> & holdouts ,
                               const std::string & subject_id )
{
  return subject_id != "" && holdouts.find( subject_id ) != holdouts.end();
}

void coda_print_confusion_matrix_3class( const std::vector<int> & obs ,
                                         const std::vector<int> & pred )
{
  const int n = obs.size();
  if ( n != pred.size() )
    Helper::halt( "internal error: unequal vectors in coda_print_confusion_matrix_3class()" );

  std::map<int,std::map<int,int> > res;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      res[i][j] = 0;

  for (int i = 0; i < n; i++)
    res[ obs[i] ][ pred[i] ]++;

  std::map<int,double> rows, cols;
  double tot = 0;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      {
        rows[i] += res[i][j];
        cols[j] += res[i][j];
        tot += res[i][j];
      }

  const std::vector<std::string> labs = { "W" , "R" , "NR" };

  logger << "\t Pred:";
  for (int j = 0; j < 3; j++)
    logger << "\t" << labs[j];
  logger << "\tTot\n";

  logger << "  Obs:";
  for (int i = 0; i < 3; i++)
    {
      logger << "\t" << labs[i];
      for (int j = 0; j < 3; j++)
        logger << "\t" << res[i][j];
      logger << "\t" << Helper::pp( tot > 0 ? rows[i] / tot : 0.0 );
      logger << "\n";
    }

  logger << "\tTot:";
  for (int j = 0; j < 3; j++)
    logger << "\t" << Helper::pp( tot > 0 ? cols[j] / tot : 0.0 );
  logger << "\t1.00\n\n";
}

int coda_effective_label( int label , const bool three_state )
{
  if ( three_state && ( label == POPS_N2 || label == POPS_N3 ) )
    return POPS_N1;
  return label;
}

int count_flagged_rows( const std::vector<bool> & flagged )
{
  int n = 0;
  for (int e = 0; e < (int)flagged.size(); e++)
    if ( flagged[e] ) ++n;
  return n;
}

void emit_coda_soft_fail( const std::string & msg ,
                          const int ne_flagged ,
                          const bool stratify_coda )
{
  logger << "  " << msg << "\n";
  if ( stratify_coda )
    writer.level( "1" , "CODA" );
  writer.value( "NE_OKAY" , 0 );
  writer.value( "NE_FLAGGED" , ne_flagged );
  writer.value( "OK" , 0 );
  writer.value( "MSG" , msg );
  if ( stratify_coda )
    writer.unlevel( "CODA" );
}

std::string coda_bool( const bool b )
{
  return b ? "true" : "false";
}

std::string coda_dbl_or_off( const double x )
{
  return x >= 0 ? Helper::dbl2str( x ) : "off";
}

bool coda_parse_bool_header( const std::string & s )
{
  std::string x = s;
  std::transform( x.begin() , x.end() , x.begin() ,
                  []( unsigned char c ) { return std::tolower( c ); } );
  if ( x == "true" || x == "t" || x == "1" ) return true;
  if ( x == "false" || x == "f" || x == "0" ) return false;
  Helper::halt( "POPS-CODA: could not parse boolean header value: " + s );
  return false;
}

struct coda_file_subject_t {
  std::string id;
  Eigen::MatrixXd P;
  std::vector<int> S;
  std::vector<int> E;
  std::vector<std::string> start;
  std::vector<std::string> stop;
  bool has_prior = false;
  bool has_start = false;
  bool has_stop = false;
};

void load_coda_posteriors_file( const std::string & filename ,
                                bool require_prior ,
                                pops_coda_t::options_t * opt ,
                                std::vector<coda_file_subject_t> * out )
{
  out->clear();

  const std::string expanded_filename = Helper::expand( filename );
  if ( !Helper::fileExists( expanded_filename ) )
    Helper::halt( "POPS-CODA: posteriors file not found: " + filename );

  std::ifstream IN = LunaIO::open_ifstream( expanded_filename );
  if ( !IN.good() )
    Helper::halt( "POPS-CODA: could not open: " + filename );

  std::string hdr;
  Helper::safe_getline( IN , hdr );
  if ( hdr.empty() )
    Helper::halt( "POPS-CODA: empty header in " + filename );

  std::vector<std::string> cols = Helper::parse( hdr , "\t " );

  int col_id = -1 , col_e = -1 , col_prior = -1;
  int col_ppW = -1 , col_ppR = -1;
  int col_ppN1 = -1 , col_ppN2 = -1 , col_ppN3 = -1 , col_ppNR = -1;
  int col_start = -1 , col_stop = -1;

  for (int c = 0; c < (int)cols.size(); c++)
    {
      const std::string & s = cols[c];
      if      ( s == "ID"    ) col_id    = c;
      else if ( s == "E"     ) col_e     = c;
      else if ( s == "PRIOR" || s == "PRIOR30" ) col_prior = c;
      else if ( s == "PP_W"  ) col_ppW   = c;
      else if ( s == "PP_R"  ) col_ppR   = c;
      else if ( s == "PP_N1" ) col_ppN1  = c;
      else if ( s == "PP_N2" ) col_ppN2  = c;
      else if ( s == "PP_N3" ) col_ppN3  = c;
      else if ( s == "PP_NR" ) col_ppNR  = c;
      else if ( s == "START" ) col_start = c;
      else if ( s == "STOP"  ) col_stop  = c;
    }

  if ( col_id  < 0 ) Helper::halt( "POPS-CODA: missing ID column in " + filename );
  if ( col_e   < 0 ) Helper::halt( "POPS-CODA: missing E column in " + filename );
  if ( require_prior && col_prior < 0 )
    Helper::halt( "POPS-CODA: missing PRIOR/PRIOR30 column in " + filename );
  if ( col_ppW < 0 ) Helper::halt( "POPS-CODA: missing PP_W column in " + filename );
  if ( col_ppR < 0 ) Helper::halt( "POPS-CODA: missing PP_R column in " + filename );

  const bool file_5state = ( col_ppN1 >= 0 && col_ppN2 >= 0 && col_ppN3 >= 0 );
  const bool file_3state = ( col_ppNR >= 0 && !file_5state );

  if ( !file_5state && !file_3state )
    Helper::halt( "POPS-CODA: file must have PP_N1+PP_N2+PP_N3 (5-state) "
                  "or PP_NR (3-state) columns; found neither" );

  if ( file_3state && opt != NULL )
    opt->three_state = true;

  const int ns_file = file_5state ? 5 : 3;

  struct Row {
    int e;
    std::array<double,5> pp;
    int prior;
    std::string start;
    std::string stop;
  };

  std::map< std::string , std::vector<Row> > subjects;

  while ( true )
    {
      std::string line;
      Helper::safe_getline( IN , line );
      if ( IN.eof() || IN.bad() ) break;
      if ( line.empty() ) continue;

      std::vector<std::string> tok = Helper::parse( line , "\t " );
      const int nt = (int)tok.size();

      auto get_dbl = [&]( int c ) -> double {
        if ( c < 0 || c >= nt ) return std::numeric_limits<double>::quiet_NaN();
        double d;
        return Helper::str2dbl( tok[c] , &d ) ? d : std::numeric_limits<double>::quiet_NaN();
      };

      auto get_str = [&]( int c ) -> std::string {
        return c >= 0 && c < nt ? tok[c] : "";
      };

      const std::string id = get_str( col_id );
      if ( id.empty() ) continue;

      int ei = -1;
      if ( col_e < nt )
        {
          Helper::str2int( tok[col_e] , &ei );
          ei -= 1;
        }
      if ( ei < 0 ) continue;

      Row r;
      r.e      = ei;
      r.pp[0]  = get_dbl( col_ppW );
      r.pp[1]  = get_dbl( col_ppR );
      r.pp[2]  = file_5state ? get_dbl( col_ppN1 ) : get_dbl( col_ppNR );
      r.pp[3]  = file_5state ? get_dbl( col_ppN2 ) : 0.0;
      r.pp[4]  = file_5state ? get_dbl( col_ppN3 ) : 0.0;
      r.start  = get_str( col_start );
      r.stop   = get_str( col_stop );

      const std::string stg = get_str( col_prior );
      if      ( stg == "W"                  ) r.prior = POPS_WAKE;
      else if ( stg == "R"                  ) r.prior = POPS_REM;
      else if ( stg == "N1" || stg == "NR" ) r.prior = POPS_N1;
      else if ( stg == "N2"                 ) r.prior = POPS_N2;
      else if ( stg == "N3"                 ) r.prior = POPS_N3;
      else                                    r.prior = POPS_UNKNOWN;

      subjects[id].push_back( r );
    }
  IN.close();

  if ( subjects.empty() )
    Helper::halt( "POPS-CODA: no data rows parsed from " + filename );

  for ( auto & kv : subjects )
    {
      std::vector<Row> & rows = kv.second;
      if ( rows.empty() ) continue;

      std::sort( rows.begin() , rows.end() ,
                 []( const Row & a , const Row & b ){ return a.e < b.e; } );

      coda_file_subject_t subj;
      subj.id = kv.first;
      subj.P.resize( rows.size() , ns_file );
      subj.S.resize( rows.size() );
      subj.E.resize( rows.size() );
      subj.start.resize( rows.size() );
      subj.stop.resize( rows.size() );
      subj.has_prior = col_prior >= 0;
      subj.has_start = col_start >= 0;
      subj.has_stop  = col_stop >= 0;

      for (int i = 0; i < (int)rows.size(); i++)
        {
          subj.E[i] = rows[i].e;
          subj.S[i] = rows[i].prior;
          subj.start[i] = rows[i].start;
          subj.stop[i] = rows[i].stop;
          for (int s = 0; s < ns_file; s++)
            subj.P(i,s) = rows[i].pp[s];
        }

      out->push_back( subj );
    }

  logger << "  POPS-CODA: read posteriors for " << out->size()
         << " subjects from " << filename << "\n";
}

}


// ============================================================
//  Helpers
// ============================================================

std::vector<std::string> pops_coda_t::stage_labels() const
{
  return opt.three_state ? pops_t::labels3 : pops_t::labels5;
}

std::string pops_coda_t::coda_default_params( int num_class , int n_iterations ) const
{
  std::string p;
  p += "boosting_type=gbdt ";
  p += "objective=multiclass ";
  p += "metric=multi_logloss ";
  p += "num_class="           + Helper::int2str( num_class ) + " ";
  p += "metric_freq=1 ";
  p += "is_training_metric=true ";
  p += "max_bin=255 ";
  p += "num_trees="           + Helper::int2str( n_iterations ) + " ";
  p += "learning_rate=0.05 ";
  p += "num_leaves=16 ";
  p += "max_depth=4 ";
  p += "feature_fraction=0.8 ";
  p += "bagging_fraction=0.8 ";
  p += "bagging_freq=1 ";
  p += "lambda_l2=5.0 ";
  p += "min_data_in_leaf=50";
  return p;
}

int pops_coda_t::scale_context_epochs( int context_30s_epochs , double row_duration_sec )
{
  if ( context_30s_epochs <= 0 ) return 0;
  if ( row_duration_sec <= 0 ) row_duration_sec = 30.0;
  return std::max( 1 , (int)std::lround( context_30s_epochs * 30.0 / row_duration_sec ) );
}

double pops_coda_t::coda_quantile( std::vector<double> v , double q )
{
  // remove NaN
  v.erase( std::remove_if( v.begin() , v.end() ,
                            [](double x){ return std::isnan(x); } ) , v.end() );
  if ( v.empty() ) return std::numeric_limits<double>::quiet_NaN();
  std::sort( v.begin() , v.end() );
  double idx = q * ( (double)v.size() - 1.0 );
  int lo = (int)idx;
  int hi = lo + 1;
  double frac = idx - lo;
  if ( hi >= (int)v.size() ) return v.back();
  return v[lo] * (1.0 - frac) + v[hi] * frac;
}


// ============================================================
//  Feature-name list  (deterministic, options-driven)
// ============================================================

std::vector<std::string> pops_coda_t::feature_names() const
{
  const std::vector<std::string> slabs = stage_labels();
  const int ns = (int)slabs.size();

  std::vector<std::string> fn;

  // A. current posteriors
  for (int s = 0; s < ns; s++) fn.push_back( "p_" + slabs[s] );

  // B. elapsed time
  fn.push_back( "elapsed_sec" );
  fn.push_back( "elapsed_frac" );

  // C. previous and next epoch posteriors
  if ( opt.include_prev_next )
    {
      for (int s = 0; s < ns; s++) fn.push_back( "prev_p_" + slabs[s] );
      for (int s = 0; s < ns; s++) fn.push_back( "next_p_" + slabs[s] );
    }

  // D. posterior quantiles in past context window
  if ( opt.include_quantiles )
    {
      for (int s = 0; s < ns; s++)
        {
          fn.push_back( "past_median_p_" + slabs[s] );
          fn.push_back( "past_q25_p_"    + slabs[s] );
          fn.push_back( "past_q75_p_"    + slabs[s] );
        }
    }

  // E. posterior quantiles in future context window
  if ( opt.include_quantiles && opt.include_future )
    {
      for (int s = 0; s < ns; s++)
        {
          fn.push_back( "future_median_p_" + slabs[s] );
          fn.push_back( "future_q25_p_"    + slabs[s] );
          fn.push_back( "future_q75_p_"    + slabs[s] );
        }
    }

  // F. hard-stage proportions in past context window
  if ( opt.include_hard_props )
    {
      for (int s = 0; s < ns; s++) fn.push_back( "past_prop_" + slabs[s] );
    }

  // G. hard-stage proportions in future context window
  if ( opt.include_hard_props && opt.include_future )
    {
      for (int s = 0; s < ns; s++) fn.push_back( "future_prop_" + slabs[s] );
    }

  // H. entropy / confidence
  if ( opt.include_entropy_margin )
    {
      fn.push_back( "top_prob" );
      fn.push_back( "second_prob" );
      fn.push_back( "margin" );
      fn.push_back( "entropy" );
    }

  // I. global posterior totals / fractions
  if ( opt.include_global )
    {
      for (int s = 0; s < ns; s++) fn.push_back( "total_p_"    + slabs[s] );
      for (int s = 0; s < ns; s++) fn.push_back( "fraction_p_" + slabs[s] );
    }

  // J. relative elapsed posterior progression (always included)
  for (int s = 0; s < ns; s++) fn.push_back( "rel_elapsed_p_" + slabs[s] );

  // J2. explicit physiology-inspired cumulative summaries
  fn.push_back( ns == 5 ? "sleep_pressure_W_minus_4xN3" : "sleep_pressure_W_minus_4xNR" );
  fn.push_back( ns == 5 ? "cum_frac_N3_in_NREM" : "cum_frac_NR_in_NREM" );

  // K. short / medium soft run-lengths (past-only persistence)
  for (int s = 0; s < ns; s++) fn.push_back( "run_short_p_" + slabs[s] );
  for (int s = 0; s < ns; s++) fn.push_back( "run_med_p_"   + slabs[s] );

  // K2. soft NREM-bout phase proxies (posterior-only, no formal cycle rules)
  fn.push_back( "soft_time_since_rem" );
  fn.push_back( "soft_nrem_run" );
  fn.push_back( "soft_nrem_phase" );

  // L. long decayed history (past-only state pressure)
  for (int s = 0; s < ns; s++) fn.push_back( "hist_long_p_" + slabs[s] );

  // M. long-history contrasts: adjacent-stage pairs in sleep architecture
  fn.push_back( "hist_long_R_minus_W" );
  fn.push_back( "hist_long_R_minus_"  + slabs[2] );   // R vs N1/NR
  fn.push_back( "hist_long_W_minus_"  + slabs[2] );   // W vs N1/NR
  if ( ns == 5 )
    {
      fn.push_back( "hist_long_N1_minus_N2" );
      fn.push_back( "hist_long_N2_minus_N3" );
      fn.push_back( "hist_long_R_minus_N2"  );
      // short/medium run-length contrasts for the N2/N3 boundary specifically
      fn.push_back( "run_short_N2_minus_N3" );
      fn.push_back( "run_med_N2_minus_N3"   );
    }

  return fn;
}


// ============================================================
//  Posterior adaptation and validation
// ============================================================

Eigen::MatrixXd pops_coda_t::adapt_posteriors( const Eigen::MatrixXd & P ) const
{
  // If 3-state mode and P has 5 columns, collapse N1+N2+N3 → NR
  if ( opt.three_state && P.cols() == 5 )
    {
      const int ne = P.rows();
      Eigen::MatrixXd P3( ne , 3 );
      for (int e = 0; e < ne; e++)
        {
          P3(e,0) = P(e,0);                       // W
          P3(e,1) = P(e,1);                       // R
          P3(e,2) = P(e,2) + P(e,3) + P(e,4);    // NR
          // renormalise row
          double s = P3(e,0) + P3(e,1) + P3(e,2);
          if ( s > 0 )
            { P3(e,0) /= s; P3(e,1) /= s; P3(e,2) /= s; }
        }
      return P3;
    }
  return P;
}

Eigen::MatrixXd pops_coda_t::validate_and_normalise( const Eigen::MatrixXd & P ) const
{
  const int expected_cols = n_classes();
  if ( P.cols() != expected_cols && !( opt.three_state && P.cols() == 5 ) )
    Helper::halt( "POPS-CODA: expected " + Helper::int2str(expected_cols) +
                  " posterior columns, got " + Helper::int2str(static_cast<long>(P.cols())) );

  Eigen::MatrixXd Pn = adapt_posteriors( P );

  // Check approximate unit-sum for non-NaN rows; normalise with warning if needed
  bool warned = false;
  const double tol = 1e-3;
  const int ne = Pn.rows();
  const int ns = Pn.cols();
  for (int e = 0; e < ne; e++)
    {
      // skip rows with any NaN
      bool has_nan = false;
      for (int s = 0; s < ns; s++)
        if ( std::isnan( Pn(e,s) ) ) { has_nan = true; break; }
      if ( has_nan ) continue;

      double rowsum = Pn.row(e).sum();
      if ( std::fabs( rowsum - 1.0 ) > tol )
        {
          if ( !warned )
            {
              logger << "  ** POPS-CODA: posteriors do not sum to 1 (max dev "
                     << std::fabs( rowsum - 1.0 )
                     << "); normalising per row\n";
              warned = true;
            }
          if ( rowsum > 0 )
            Pn.row(e) /= rowsum;
        }
    }
  return Pn;
}


// ============================================================
//  Core feature computation
// ============================================================

void pops_coda_t::make_features( const Eigen::MatrixXd & P ,
                                  const std::vector<int> & E ,
                                  const std::vector<bool> & flagged ,
                                  Eigen::MatrixXd & X ,
                                  std::vector<std::string> & fnames ) const
{
  const int ne  = P.rows();
  const int ns  = P.cols();   // already adapted (3 or 5)
  const double NaN = std::numeric_limits<double>::quiet_NaN();
  const double eps = 1e-12;

  fnames = feature_names();  // deterministic names
  const int nf = (int)fnames.size();

  // Initialise all features to NaN; flagged epochs remain NaN throughout
  X = Eigen::MatrixXd::Constant( ne , nf , NaN );

  if ( ne == 0 ) return;

  // ---- valid[e] = !flagged[e] ----------------------------------------
  // Build singly-linked lists: prev_valid[e] / next_valid[e]
  // = index of nearest valid epoch strictly before/after e, or -1
  std::vector<int> prev_valid( ne , -1 );
  std::vector<int> next_valid( ne , -1 );
  {
    int last = -1;
    for (int e = 0; e < ne; e++)
      { prev_valid[e] = last; if ( !flagged[e] ) last = e; }
    last = -1;
    for (int e = ne-1; e >= 0; e--)
      { next_valid[e] = last; if ( !flagged[e] ) last = e; }
  }

  // ---- Count valid epochs and span of valid epoch indices ------------
  int n_valid = 0;
  int first_E = (E.size() > 0) ? E[0] : 0;
  int last_E = first_E;
  bool found_first = false;
  for (int e = 0; e < ne; e++)
    {
      if ( !flagged[e] )
        {
          n_valid++;
          if ( !found_first ) { first_E = E[e]; found_first = true; }
          last_E = E[e];
        }
    }

  // ---- Global totals (I) and cumulative sums (J) ---------------------
  std::vector<double> total_p( ns , 0.0 );
  for (int e = 0; e < ne; e++)
    {
      if ( flagged[e] ) continue;
      for (int s = 0; s < ns; s++) total_p[s] += P(e,s);
    }
  std::vector<double> fraction_p( ns , 0.0 );
  if ( n_valid > 0 )
    for (int s = 0; s < ns; s++) fraction_p[s] = total_p[s] / n_valid;

  // Running cumulative sum (over valid epochs only, in E order)
  std::vector<std::vector<double>> cum_p( ne , std::vector<double>( ns , 0.0 ) );
  for (int e = 0; e < ne; e++)
    {
      if ( e > 0 )
        for (int s = 0; s < ns; s++) cum_p[e][s] = cum_p[e-1][s];
      if ( !flagged[e] )
        for (int s = 0; s < ns; s++) cum_p[e][s] += P(e,s);
    }

  // ---- Derived cumulative physiology-inspired summaries --------------
  const int idx_w = POPS_WAKE;
  const int idx_r = POPS_REM;
  const int idx_n1 = POPS_N1;
  const int idx_n2 = ns == 5 ? POPS_N2 : POPS_N1;
  const int idx_n3 = ns == 5 ? POPS_N3 : POPS_N1;

  std::vector<double> cum_nrem( ne , 0.0 );
  std::vector<double> cum_deep( ne , 0.0 );
  std::vector<double> cum_frac_deep_in_nrem( ne , 0.0 );
  std::vector<double> sleep_pressure( ne , 0.0 );

  for (int e = 0; e < ne; e++)
    {
      double cn = cum_p[e][idx_n1];
      if ( ns == 5 )
        cn += cum_p[e][idx_n2] + cum_p[e][idx_n3];
      cum_nrem[e] = cn;
      cum_deep[e] = cum_p[e][idx_n3];
      cum_frac_deep_in_nrem[e] = cn > eps ? cum_deep[e] / cn : 0.0;
      sleep_pressure[e] = cum_p[e][idx_w] - 4.0 * cum_deep[e];
    }

  // ---- Past-only temporal memory features ----------------------------
  // run_* encode fractional persistence of the current state signal:
  // high current p_s grows the run, weak current p_s collapses it.
  // hist_long encodes longer-range recent pressure for each state.
  const double run_short_alpha = 0.70;
  const double run_med_alpha   = 0.85;
  const double hist_long_alpha = 0.97;

  std::vector<std::vector<double>> run_short( ne , std::vector<double>( ns , 0.0 ) );
  std::vector<std::vector<double>> run_med  ( ne , std::vector<double>( ns , 0.0 ) );
  std::vector<std::vector<double>> hist_long( ne , std::vector<double>( ns , 0.0 ) );
  std::vector<double> soft_time_since_rem( ne , 0.0 );
  std::vector<double> soft_nrem_run( ne , 0.0 );
  std::vector<double> soft_nrem_phase( ne , 0.0 );

  const double soft_tsr_alpha = 0.95;
  const double soft_tsr_rem_reset = 0.90;
  const double soft_nrem_alpha = 0.90;
  const double soft_phase_c = 5.0;

  for (int e = 0; e < ne; e++)
    {
      if ( flagged[e] ) continue;

      const int pv = prev_valid[e];
      for (int s = 0; s < ns; s++)
        {
          const double prev_run_short = pv >= 0 ? run_short[pv][s] : 0.0;
          const double prev_run_med   = pv >= 0 ? run_med[pv][s]   : 0.0;
          const double prev_hist_long = pv >= 0 ? hist_long[pv][s] : 0.0;

          run_short[e][s] = P(e,s) * ( 1.0 + run_short_alpha * prev_run_short );
          run_med[e][s]   = P(e,s) * ( 1.0 + run_med_alpha   * prev_run_med   );
          hist_long[e][s] = P(e,s) + hist_long_alpha * prev_hist_long;
        }

      const double p_r = P(e,idx_r);
      double p_nrem = P(e,idx_n1);
      if ( ns == 5 )
        p_nrem += P(e,idx_n2) + P(e,idx_n3);

      const double prev_tsr = pv >= 0 ? soft_time_since_rem[pv] : 0.0;
      const double prev_nrem_run = pv >= 0 ? soft_nrem_run[pv] : 0.0;

      soft_time_since_rem[e] =
        p_nrem * ( 1.0 + soft_tsr_alpha * prev_tsr ) * ( 1.0 - soft_tsr_rem_reset * p_r );

      soft_nrem_run[e] =
        p_nrem + soft_nrem_alpha * prev_nrem_run * std::max( 0.0 , 1.0 - p_r );

      soft_nrem_phase[e] = soft_time_since_rem[e] / ( soft_time_since_rem[e] + soft_phase_c );
    }

  // ---- Fill each non-flagged epoch -----------------------------------
  for (int e = 0; e < ne; e++)
    {
      if ( flagged[e] ) continue;

      int col = 0;

      // A. Current posteriors
      for (int s = 0; s < ns; s++) X(e, col++) = P(e,s);

      // B. Elapsed time
      double esec = ( E[e] - first_E ) * opt.row_duration_sec;
      X(e, col++) = esec;
      X(e, col++) = ( last_E > first_E ) ? (double)( E[e] - first_E ) / (double)( last_E - first_E ) : 0.0;

      // C. Previous / next epoch posteriors
      if ( opt.include_prev_next )
        {
          int pv = prev_valid[e];
          int nv = next_valid[e];
          for (int s = 0; s < ns; s++) X(e, col++) = (pv >= 0) ? P(pv,s) : NaN;
          for (int s = 0; s < ns; s++) X(e, col++) = (nv >= 0) ? P(nv,s) : NaN;
        }

      // Helper: walk the prev/next chain collecting up to context_epochs values
      auto collect_past = [&]( int start , int ctx ) -> std::vector<std::vector<double>>
        {
          std::vector<std::vector<double>> vals( ns );
          int u = start; int cnt = 0;
          while ( u >= 0 && cnt < ctx )
            { for (int s=0;s<ns;s++) vals[s].push_back(P(u,s)); cnt++; u = prev_valid[u]; }
          return vals;
        };
      auto collect_future = [&]( int start , int ctx ) -> std::vector<std::vector<double>>
        {
          std::vector<std::vector<double>> vals( ns );
          int u = start; int cnt = 0;
          while ( u >= 0 && cnt < ctx )
            { for (int s=0;s<ns;s++) vals[s].push_back(P(u,s)); cnt++; u = next_valid[u]; }
          return vals;
        };

      // D. Past posterior quantiles
      if ( opt.include_quantiles )
        {
          auto pvals = collect_past( prev_valid[e] , opt.context_epochs );
          for (int s = 0; s < ns; s++)
            {
              if ( pvals[s].empty() )
                { X(e,col++) = NaN; X(e,col++) = NaN; X(e,col++) = NaN; }
              else
                {
                  X(e, col++) = coda_quantile( pvals[s] , 0.50 );
                  X(e, col++) = coda_quantile( pvals[s] , 0.25 );
                  X(e, col++) = coda_quantile( pvals[s] , 0.75 );
                }
            }
        }

      // E. Future posterior quantiles
      if ( opt.include_quantiles && opt.include_future )
        {
          auto fvals = collect_future( next_valid[e] , opt.context_epochs );
          for (int s = 0; s < ns; s++)
            {
              if ( fvals[s].empty() )
                { X(e,col++) = NaN; X(e,col++) = NaN; X(e,col++) = NaN; }
              else
                {
                  X(e, col++) = coda_quantile( fvals[s] , 0.50 );
                  X(e, col++) = coda_quantile( fvals[s] , 0.25 );
                  X(e, col++) = coda_quantile( fvals[s] , 0.75 );
                }
            }
        }

      // F. Hard-stage proportions in past context
      if ( opt.include_hard_props )
        {
          std::vector<int> phist( ns , 0 );
          int pcnt = 0;
          int u = prev_valid[e]; int cnt2 = 0;
          while ( u >= 0 && cnt2 < opt.context_epochs )
            {
              int hmax = 0;
              for (int s=1;s<ns;s++) if (P(u,s) > P(u,hmax)) hmax = s;
              phist[hmax]++; pcnt++; cnt2++; u = prev_valid[u];
            }
          if ( pcnt == 0 )
            for (int s=0;s<ns;s++) X(e,col++) = NaN;
          else
            for (int s=0;s<ns;s++) X(e,col++) = (double)phist[s] / pcnt;
        }

      // G. Hard-stage proportions in future context
      if ( opt.include_hard_props && opt.include_future )
        {
          std::vector<int> fhist( ns , 0 );
          int fcnt = 0;
          int u = next_valid[e]; int cnt2 = 0;
          while ( u >= 0 && cnt2 < opt.context_epochs )
            {
              int hmax = 0;
              for (int s=1;s<ns;s++) if (P(u,s) > P(u,hmax)) hmax = s;
              fhist[hmax]++; fcnt++; cnt2++; u = next_valid[u];
            }
          if ( fcnt == 0 )
            for (int s=0;s<ns;s++) X(e,col++) = NaN;
          else
            for (int s=0;s<ns;s++) X(e,col++) = (double)fhist[s] / fcnt;
        }

      // H. Entropy / confidence
      if ( opt.include_entropy_margin )
        {
          double p1 = -1.0 , p2 = -1.0;
          for (int s = 0; s < ns; s++)
            {
              double pv = P(e,s);
              if ( pv > p1 ) { p2 = p1; p1 = pv; }
              else if ( pv > p2 ) p2 = pv;
            }
          if ( p2 < 0.0 ) p2 = 0.0;
          double ent = 0.0;
          for (int s = 0; s < ns; s++)
            {
              double pv = P(e,s);
              ent -= pv * std::log( pv + eps );
            }
          X(e, col++) = p1;
          X(e, col++) = p2;
          X(e, col++) = p1 - p2;
          X(e, col++) = ent;
        }

      // I. Global totals / fractions
      if ( opt.include_global )
        {
          for (int s = 0; s < ns; s++) X(e, col++) = total_p[s];
          for (int s = 0; s < ns; s++) X(e, col++) = fraction_p[s];
        }

      // J. Relative elapsed posterior progression
      for (int s = 0; s < ns; s++)
        X(e, col++) = cum_p[e][s] / ( total_p[s] + opt.lambda );

      // J2. Explicit physiology-inspired cumulative context summaries
      X(e, col++) = sleep_pressure[e];
      X(e, col++) = cum_frac_deep_in_nrem[e];

      // K. Short / medium soft run-lengths
      for (int s = 0; s < ns; s++) X(e, col++) = run_short[e][s];
      for (int s = 0; s < ns; s++) X(e, col++) = run_med[e][s];

      // K2. Soft NREM-cycle proxy features
      X(e, col++) = soft_time_since_rem[e];
      X(e, col++) = soft_nrem_run[e];
      X(e, col++) = soft_nrem_phase[e];

      // L. Long decayed history
      for (int s = 0; s < ns; s++) X(e, col++) = hist_long[e][s];

      // M. Long-history contrasts: adjacent-stage pairs in sleep architecture
      X(e, col++) = hist_long[e][POPS_REM]  - hist_long[e][POPS_WAKE];
      X(e, col++) = hist_long[e][POPS_REM]  - hist_long[e][POPS_N1];
      X(e, col++) = hist_long[e][POPS_WAKE] - hist_long[e][POPS_N1];
      if ( ns == 5 )
        {
          X(e, col++) = hist_long[e][POPS_N1] - hist_long[e][POPS_N2];
          X(e, col++) = hist_long[e][POPS_N2] - hist_long[e][POPS_N3];
          X(e, col++) = hist_long[e][POPS_REM] - hist_long[e][POPS_N2];
          X(e, col++) = run_short[e][POPS_N2]  - run_short[e][POPS_N3];
          X(e, col++) = run_med[e][POPS_N2]    - run_med[e][POPS_N3];
        }

      // Sanity: should have filled exactly nf columns
      if ( col != nf )
        Helper::halt( "POPS-CODA internal error: col=" + Helper::int2str(col) +
                      " nf=" + Helper::int2str(nf) );
    }
}


// ============================================================
//  Training
// ============================================================

int pops_coda_t::accumulate( const Eigen::MatrixXd & P_in ,
                              const std::vector<int> & S ,
                              const std::vector<int> & E ,
                              const std::string & subject_id )
{
  const int ne = P_in.rows();
  if ( ne == 0 ) return 0;
  if ( (int)S.size() != ne || (int)E.size() != ne )
    Helper::halt( "POPS-CODA: P/S/E size mismatch in accumulate()" );

  // Adapt and validate posteriors
  Eigen::MatrixXd P = validate_and_normalise( P_in );
  const int ns = P.cols();   // effective classes after adaptation

  // Build flagged[] from missing labels or NaN posteriors
  std::vector<bool> flagged( ne , false );
  for (int e = 0; e < ne; e++)
    {
      if ( S[e] == POPS_UNKNOWN )
        { flagged[e] = true; continue; }
      for (int s = 0; s < ns; s++)
        if ( std::isnan( P(e,s) ) ) { flagged[e] = true; break; }
    }

  int n_valid = 0;
  int longest_run = 0;
  int current_run = 0;
  int prev_valid_epoch = -1;

  for (int e = 0; e < ne; e++)
    {
      if ( flagged[e] )
        {
          current_run = 0;
          prev_valid_epoch = -1;
          continue;
        }

      ++n_valid;

      if ( prev_valid_epoch >= 0 && E[e] == prev_valid_epoch + 1 )
        ++current_run;
      else
        current_run = 1;

      if ( current_run > longest_run )
        longest_run = current_run;

      prev_valid_epoch = E[e];
    }

  if ( opt.min_subject_contiguous_epochs > 0 &&
       longest_run < opt.min_subject_contiguous_epochs )
    {
      logger << "  POPS-CODA: skipping subject"
             << ( subject_id == "" ? "" : " " + subject_id )
             << " (" << n_valid << " valid epochs; longest contiguous run = "
             << longest_run << " < required "
             << opt.min_subject_contiguous_epochs << ")\n";
      return 0;
    }

  // Build features
  Eigen::MatrixXd X;
  std::vector<std::string> fn;
  make_features( P , E , flagged , X , fn );

  std::vector<Eigen::VectorXd> * X_rows = &X_train_rows;
  std::vector<int> * S_rows = &S_train;
  std::vector<float> * W_rows = &W_train;
  if ( coda_subject_in_holdouts( holdouts , subject_id ) )
    {
      X_rows = &X_valid_rows;
      S_rows = &S_valid;
      W_rows = &W_valid;
    }

  const double row_duration_sec = opt.row_duration_sec <= 0 ? 30.0 : opt.row_duration_sec;
  const bool can_mix =
    row_duration_sec < 30.0 &&
    std::fabs( std::round( 30.0 / row_duration_sec ) - ( 30.0 / row_duration_sec ) ) < 1e-6;
  const int rows_per_epoch = can_mix ? (int)std::lround( 30.0 / row_duration_sec ) : 1;
  const int mid_offset = can_mix ? ( rows_per_epoch - 1 ) / 2 : 0;

  std::map<int,int> epoch_labels;
  if ( can_mix )
    {
      for (int e = 0; e < ne; e++)
        {
          if ( flagged[e] ) continue;
          const int obs_epoch = E[e] / rows_per_epoch;
          if ( epoch_labels.find( obs_epoch ) == epoch_labels.end() )
            epoch_labels[ obs_epoch ] = S[e];
        }
    }

  int n_pure = 0;
  int n_mixed_split = 0;
  int n_mixed_same = 0;
  int n_mixed_partial = 0;

  auto add_row = [&]( int e , int label , float weight )
    {
      if ( weight <= 0 ) return;
      X_rows->push_back( X.row(e) );
      S_rows->push_back( label );
      W_rows->push_back( weight );
    };

  // Collect rows: non-flagged only. In sub-30s mode, expand rows that span
  // two 30s labels into up to two weighted training rows.
  int n_added = 0;
  for (int e = 0; e < ne; e++)
    {
      if ( flagged[e] ) continue;

      // Convert label to CODA class index
      const int label = coda_effective_label( S[e] , opt.three_state );
      if ( label < 0 || label >= ns ) continue;

      if ( ! can_mix )
        {
          add_row( e , label , 1.0f );
          ++n_added;
          ++n_pure;
          continue;
        }

      const int window_start_row = E[e] - mid_offset;
      int offset = window_start_row % rows_per_epoch;
      if ( offset < 0 ) offset += rows_per_epoch;

      if ( offset == 0 )
        {
          add_row( e , label , 1.0f );
          ++n_added;
          ++n_pure;
          continue;
        }

      const int obs_epoch = E[e] / rows_per_epoch;
      const std::map<int,int>::const_iterator nn = epoch_labels.find( obs_epoch + 1 );
      const int next_label_raw = nn == epoch_labels.end() ? POPS_UNKNOWN : nn->second;
      const int next_label = coda_effective_label( next_label_raw , opt.three_state );

      const float w_cur = (float)( rows_per_epoch - offset ) / (float)rows_per_epoch;
      const float w_nxt = (float)offset / (float)rows_per_epoch;

      if ( next_label_raw == POPS_UNKNOWN || next_label < 0 || next_label >= ns )
        {
          add_row( e , label , w_cur );
          ++n_added;
          ++n_mixed_partial;
          continue;
        }

      if ( next_label == label )
        {
          add_row( e , label , 1.0f );
          ++n_added;
          ++n_mixed_same;
          continue;
        }

      add_row( e , label , w_cur );
      add_row( e , next_label , w_nxt );
      n_added += 2;
      ++n_mixed_split;
    }

  if ( can_mix )
    logger << "  POPS-CODA: "
           << ( subject_id == "" ? "training rows" : subject_id )
           << " pure=" << n_pure
           << ", mixed-split=" << n_mixed_split
           << ", mixed-same=" << n_mixed_same
           << ", mixed-partial=" << n_mixed_partial << "\n";

  return n_added;
}

void pops_coda_t::load_validation_ids( const std::string & filename )
{
  holdouts.clear();

  const std::string expanded = Helper::expand( filename );
  if ( ! Helper::fileExists( expanded ) )
    Helper::halt( "POPS-CODA: could not open " + filename );

  std::ifstream in = LunaIO::open_ifstream( expanded , std::ios::in );
  while ( ! in.eof() )
    {
      std::string id;
      in >> id;
      if ( id == "" || in.eof() ) break;
      holdouts.insert( id );
    }
  in.close();

  logger << "  POPS-CODA: read " << holdouts.size()
         << " validation dataset individuals from "
         << filename << "\n";
}

void pops_coda_t::load_exclude_ids( const std::string & filename )
{
  exclude_ids.clear();

  const std::string expanded = Helper::expand( filename );
  if ( ! Helper::fileExists( expanded ) )
    Helper::halt( "POPS-CODA: could not open " + filename );

  std::ifstream in = LunaIO::open_ifstream( expanded , std::ios::in );
  while ( ! in.eof() )
    {
      std::string id;
      in >> id;
      if ( id == "" || in.eof() ) break;
      exclude_ids.insert( id );
    }
  in.close();

  logger << "  POPS-CODA: read " << exclude_ids.size()
         << " exclusion individual"
         << ( exclude_ids.size() == 1 ? "" : "s" )
         << " from " << filename << "\n";
}


void pops_coda_t::train_model( const std::string & config_file ,
                                int n_iterations )
{
  const int nrows = (int)X_train_rows.size();
  if ( nrows == 0 )
    Helper::halt( "POPS-CODA: no training data accumulated" );

  const int nf = (int)X_train_rows[0].size();

  // Assemble matrix from buffered rows
  Eigen::MatrixXd X_train( nrows , nf );
  for (int r = 0; r < nrows; r++)
    X_train.row(r) = X_train_rows[r];

  const int nrows_valid = (int)X_valid_rows.size();
  Eigen::MatrixXd X_valid;
  if ( nrows_valid )
    {
      X_valid.resize( nrows_valid , nf );
      for (int r = 0; r < nrows_valid; r++)
        X_valid.row(r) = X_valid_rows[r];
    }

  if ( (int)W_train.size() != nrows )
    Helper::halt( "POPS-CODA: training row weights do not match training rows" );
  if ( nrows_valid && (int)W_valid.size() != nrows_valid )
    Helper::halt( "POPS-CODA: validation row weights do not match validation rows" );

  const int nc = n_classes();

  logger << "  POPS-CODA: training on " << nrows << " epochs, "
         << nf << " features, " << nc << " classes, "
         << n_iterations << " iterations";
  if ( nrows_valid )
    logger << " with " << nrows_valid << " validation epochs";
  logger << "\n";

  // Set up LightGBM config
  if ( config_file == "." || config_file == "" )
    lgbm.params = coda_default_params( nc , n_iterations );
  else
    lgbm.load_config( Helper::expand( config_file ) );

  lgbm.n_iterations = n_iterations;
  lgbm.early_stopping_rounds = 10;

  // Attach training data
  lgbm.attach_training_matrix( X_train );
  lgbm.attach_training_labels( S_train );
  if ( nrows_valid )
    {
      lgbm.attach_validation_matrix( X_valid );
      lgbm.attach_validation_labels( S_valid );
    }

  lgbm.training_weights = W_train;
  if ( nrows_valid )
    lgbm.validation_weights = W_valid;

  // Apply uniform per-class weights (can be extended later)
  std::vector<std::string> clabs = stage_labels();
  std::vector<double>      cwgts;

  if ( opt.class_weights.empty() )
    {
      // Default to modestly upweighting harder / less prevalent sleep stages.
      if ( opt.three_state )
        cwgts = { 1.0 , 1.5 , 1.25 };          // W, R, NR
      else
        cwgts = { 1.0 , 1.5 , 1.5 , 1.0 , 1.0 }; // W, R, N1, N2, N3
    }
  else
    {
      cwgts = opt.class_weights;
    }

  if ( (int)cwgts.size() != nc )
    Helper::halt( "POPS-CODA: expected " + Helper::int2str(nc) +
                  " class weights but got " + Helper::int2str( (int)cwgts.size() ) );

  logger << "  POPS-CODA: class weights";
  for (int i = 0; i < nc; i++)
    logger << " " << clabs[i] << "=" << cwgts[i];
  logger << "\n";

  lgbm_label_t lw( clabs , cwgts );
  lgbm.apply_weights( lgbm.training , &lgbm.training_weights );
  lgbm.add_label_weights( lgbm.training , &lgbm.training_weights , lw );
  lgbm.apply_weights( lgbm.training , &lgbm.training_weights );
  if ( nrows_valid )
    {
      lgbm.apply_weights( lgbm.validation , &lgbm.validation_weights );
      lgbm.add_label_weights( lgbm.validation , &lgbm.validation_weights , lw );
      lgbm.apply_weights( lgbm.validation , &lgbm.validation_weights );
    }

  lgbm.create_booster();

  const int fitted_nc = lgbm_t::classes( lgbm.booster );
  if ( fitted_nc != nc )
    Helper::halt( "POPS-CODA: trained booster reports " + Helper::int2str( fitted_nc ) +
                  " class" + ( fitted_nc == 1 ? "" : "es" ) +
                  " but expected " + Helper::int2str( nc ) +
                  " -- check CODA LightGBM params/objective" );
}


// ============================================================
//  Save / load
// ============================================================

void pops_coda_t::save( const std::string & f )
{
  if ( !lgbm.has_booster )
    Helper::halt( "POPS-CODA: no trained model to save" );

  const std::string model_file = Helper::expand( f );
  lgbm.save_model( model_file );

  // Write feature names alongside model
  std::string fn_file = coda_feature_name_file( model_file );
  std::ofstream FN = LunaIO::open_ofstream( Helper::expand( fn_file ) , std::ios::out );
  if ( !FN.good() )
    Helper::halt( "POPS-CODA: could not write feature names to " + fn_file );

  std::vector<std::string> fn = feature_names();
  const int nc = n_classes();

  FN << "# POPS-CODA feature names\n";
  FN << "# model_file=" << model_file << "\n";
  FN << "# n_classes=" << nc << "\n";
  FN << "# n_features=" << fn.size() << "\n";
  FN << "# three_state=" << coda_bool( opt.three_state ) << "\n";
  FN << "# row_duration_sec=" << Helper::dbl2str( opt.row_duration_sec ) << "\n";
  FN << "# context_epochs=" << opt.context_epochs << "\n";
  FN << "# lambda=" << Helper::dbl2str( opt.lambda ) << "\n";
  FN << "# include_future=" << coda_bool( opt.include_future ) << "\n";
  FN << "# include_global=" << coda_bool( opt.include_global ) << "\n";
  FN << "# include_prev_next=" << coda_bool( opt.include_prev_next ) << "\n";
  FN << "# include_quantiles=" << coda_bool( opt.include_quantiles ) << "\n";
  FN << "# include_hard_props=" << coda_bool( opt.include_hard_props ) << "\n";
  FN << "# include_entropy_margin=" << coda_bool( opt.include_entropy_margin ) << "\n";
  FN << "# require_all_stages=" << coda_bool( opt.require_all_stages ) << "\n";
  FN << "# broad_stage_qc=" << coda_bool( opt.broad_stage_qc ) << "\n";
  FN << "# min_subject_contiguous_epochs=" << opt.min_subject_contiguous_epochs << "\n";
  FN << "# min_stage_minutes=" << Helper::stringize( opt.min_stage_minutes ) << "\n";
  FN << "# min_stage1_kappa=" << coda_dbl_or_off( opt.min_stage1_kappa ) << "\n";
  FN << "# min_stage1_kappa3=" << coda_dbl_or_off( opt.min_stage1_kappa3 ) << "\n";
  FN << "# class_weights=" << ( opt.class_weights.empty() ? std::string("default") : Helper::stringize( opt.class_weights ) ) << "\n";
  FN << "# feature_names_begin\n";
  for (const auto & s : fn) FN << s << "\n";
  FN.close();

  logger << "  POPS-CODA: wrote " << fn.size() << " feature names to " << fn_file << "\n";
}


void pops_coda_t::load( const std::string & f )
{
  const std::string model_file = Helper::expand( f );
  lgbm.load_model( model_file );

  // Load feature names
  std::string fn_file = coda_feature_name_file( model_file );
  const std::string legacy_fn_file = model_file + ".fnames";
  if ( !Helper::fileExists( Helper::expand( fn_file ) ) )
    {
      if ( Helper::fileExists( Helper::expand( legacy_fn_file ) ) )
        {
          fn_file = legacy_fn_file;
        }
      else
        {
          logger << "  ** POPS-CODA: feature-name file not found (" << fn_file
                 << "); skipping validation\n";
          return;
        }
    }

  loaded_fnames.clear();
  bool have_row_duration_sec = false;
  double loaded_row_duration_sec = 0.0;
  bool have_context_epochs = false;
  int loaded_context_epochs = 0;
  bool have_three_state = false;
  bool loaded_three_state = false;
  bool have_lambda = false;
  double loaded_lambda = 0.0;
  bool have_include_future = false;
  bool loaded_include_future = false;
  bool have_include_global = false;
  bool loaded_include_global = false;
  bool have_include_prev_next = false;
  bool loaded_include_prev_next = false;
  bool have_include_quantiles = false;
  bool loaded_include_quantiles = false;
  bool have_include_hard_props = false;
  bool loaded_include_hard_props = false;
  bool have_include_entropy_margin = false;
  bool loaded_include_entropy_margin = false;
  std::ifstream FN = LunaIO::open_ifstream( Helper::expand( fn_file ) , std::ios::in );
  while ( true )
    {
      std::string line;
      Helper::safe_getline( FN , line );
      if ( FN.eof() || FN.bad() ) break;
      if ( line == "" ) continue;
      if ( line[0] == '#' )
        {
          if ( line.find( "# row_duration_sec=" ) == 0 )
            {
              have_row_duration_sec = true;
              loaded_row_duration_sec = std::stod( line.substr( 19 ) );
            }
          else if ( line.find( "# context_epochs=" ) == 0 )
            {
              have_context_epochs = true;
              loaded_context_epochs = std::stoi( line.substr( 17 ) );
            }
          else if ( line.find( "# three_state=" ) == 0 )
            {
              have_three_state = true;
              loaded_three_state = coda_parse_bool_header( line.substr( 14 ) );
            }
          else if ( line.find( "# lambda=" ) == 0 )
            {
              have_lambda = true;
              loaded_lambda = std::stod( line.substr( 9 ) );
            }
          else if ( line.find( "# include_future=" ) == 0 )
            {
              have_include_future = true;
              loaded_include_future = coda_parse_bool_header( line.substr( 17 ) );
            }
          else if ( line.find( "# include_global=" ) == 0 )
            {
              have_include_global = true;
              loaded_include_global = coda_parse_bool_header( line.substr( 17 ) );
            }
          else if ( line.find( "# include_prev_next=" ) == 0 )
            {
              have_include_prev_next = true;
              loaded_include_prev_next = coda_parse_bool_header( line.substr( 20 ) );
            }
          else if ( line.find( "# include_quantiles=" ) == 0 )
            {
              have_include_quantiles = true;
              loaded_include_quantiles = coda_parse_bool_header( line.substr( 20 ) );
            }
          else if ( line.find( "# include_hard_props=" ) == 0 )
            {
              have_include_hard_props = true;
              loaded_include_hard_props = coda_parse_bool_header( line.substr( 21 ) );
            }
          else if ( line.find( "# include_entropy_margin=" ) == 0 )
            {
              have_include_entropy_margin = true;
              loaded_include_entropy_margin = coda_parse_bool_header( line.substr( 25 ) );
            }
          continue;
        }
      loaded_fnames.push_back( line );
    }
  FN.close();

  // Prediction-time feature semantics come from the model metadata, not CLI flags.
  if ( have_row_duration_sec ) opt.row_duration_sec = loaded_row_duration_sec;
  if ( have_context_epochs ) opt.context_epochs = loaded_context_epochs;
  if ( have_three_state ) opt.three_state = loaded_three_state;
  if ( have_lambda ) opt.lambda = loaded_lambda;
  if ( have_include_future ) opt.include_future = loaded_include_future;
  if ( have_include_global ) opt.include_global = loaded_include_global;
  if ( have_include_prev_next ) opt.include_prev_next = loaded_include_prev_next;
  if ( have_include_quantiles ) opt.include_quantiles = loaded_include_quantiles;
  if ( have_include_hard_props ) opt.include_hard_props = loaded_include_hard_props;
  if ( have_include_entropy_margin ) opt.include_entropy_margin = loaded_include_entropy_margin;

  // Validate against current options
  std::vector<std::string> expected = feature_names();
  if ( loaded_fnames.size() != expected.size() )
    Helper::halt( "POPS-CODA: feature count mismatch: loaded=" +
                  Helper::int2str( (int)loaded_fnames.size() ) +
                  " expected=" + Helper::int2str( (int)expected.size() ) +
                  " for " + model_file );
  for (int i = 0; i < (int)expected.size(); i++)
    if ( loaded_fnames[i] != expected[i] )
      Helper::halt( "POPS-CODA: feature-name mismatch at column " +
                    Helper::int2str( i + 1 ) +
                    " for " + model_file +
                    " (loaded=" + loaded_fnames[i] +
                    ", expected=" + expected[i] + ")" );

  logger << "  POPS-CODA: loaded " << loaded_fnames.size()
         << " feature names from " << fn_file << "\n";
}


// ============================================================
//  Prediction
// ============================================================

void pops_coda_t::predict( const Eigen::MatrixXd & P_in ,
                            const std::vector<int> & E ,
                            const std::vector<bool> & flagged ,
                            int ne_total ,
                            bool has_staging ,
                            const std::vector<int> & S ,
                            edf_t * pedf ,
                            const std::vector<std::string> * start ,
                            const std::vector<std::string> * stop ,
                            bool emit_stage1 )
{
  const bool stratify_coda = opt.output_both_stages;
  if ( !lgbm.has_booster )
    Helper::halt( "POPS-CODA: no model loaded" );

  const int ne = P_in.rows();
  if ( ne == 0 )
    {
      emit_coda_soft_fail( pops_coda_soft_fail_msg , 0 , stratify_coda );
      return;
    }

  const int n_flagged = count_flagged_rows( flagged );
  const int n_usable = ne - n_flagged;
  if ( n_usable == 0 )
    {
      emit_coda_soft_fail( pops_coda_soft_fail_msg , ne , stratify_coda );
      return;
    }

  Eigen::MatrixXd P_sanitised = P_in;
  const double uniform = 1.0 / (double)P_in.cols();
  for (int e = 0; e < ne; e++)
    if ( flagged[e] )
      for (int j = 0; j < P_sanitised.cols(); j++)
        P_sanitised(e,j) = uniform;

  // Adapt and validate posteriors
  Eigen::MatrixXd P = validate_and_normalise( P_sanitised );
  const int ns = P.cols();

  const int nc = n_classes();
  if ( ns != nc )
    Helper::halt( "POPS-CODA: posterior columns (" + Helper::int2str(ns) +
                  ") do not match model classes (" + Helper::int2str(nc) + ")" );

  // Count flagged epochs for logging
  logger << "  POPS-CODA: " << nc << "-class prediction, "
         << ne << " valid epochs (" << n_flagged << " flagged), "
         << n_features() << " features\n";

  // Build features
  Eigen::MatrixXd X;
  std::vector<std::string> fn;
  make_features( P , E , flagged , X , fn );

  // Apply CODA model: get posteriors
  Eigen::MatrixXd P_coda = lgbm.predict( X );   // ne × nc

  if ( opt.do_SHAP )
    SHAP( X , fn , E , flagged );

  if ( (int)P_coda.rows() != ne || (int)P_coda.cols() != nc )
    Helper::halt( "POPS-CODA: predict returned " + Helper::int2str(static_cast<long>(P_coda.rows())) +
                  " x " + Helper::int2str(static_cast<long>(P_coda.cols())) +
                  " but expected " + Helper::int2str(ne) + " x " + Helper::int2str(nc) );

  // Build epoch-to-index map (shared by both output blocks)
  std::map<int,int> e2e;
  for (int e = 0; e < ne; e++) e2e[ E[e] ] = e;

  const std::vector<std::string> & slabs = stage_labels();
  const double epoch_to_min = opt.row_duration_sec / 60.0;
  std::vector<int> pred_stage( ne , POPS_UNKNOWN );
  std::vector<double> pred_conf( ne , std::numeric_limits<double>::quiet_NaN() );
  std::vector<double> dur_pred1( nc , 0.0 );
  std::vector<double> dur_predf( nc , 0.0 );
  std::vector<double> dur_obs( nc , 0.0 );
  std::vector<double> dur_obs_orig( nc , 0.0 );

  int slp_lat_obs = -1 , slp_lat_prd = -1;
  int rem_lat_obs = -1 , rem_lat_prd = -1;
  double mean_conf = 0.0;
  int n_pred_valid = 0;

  std::vector<int> obs, pred;
  obs.reserve( ne );
  pred.reserve( ne );

  for (int e = 0; e < ne; e++)
    {
      if ( flagged[e] ) continue;

      int predx = 0;
      double pmax = P_coda.row(e).maxCoeff( &predx );
      pred_stage[e] = predx;
      pred_conf[e] = pmax;
      mean_conf += pmax;
      ++n_pred_valid;

      if ( predx >= 0 && predx < nc ) dur_pred1[predx]++;
      for (int ss = 0; ss < nc; ss++) dur_predf[ss] += P_coda(e,ss);

      if ( slp_lat_prd == -1 && predx != POPS_WAKE )
        slp_lat_prd = E[e];
      if ( rem_lat_prd == -1 && predx == POPS_REM )
        rem_lat_prd = E[e] - slp_lat_prd;

      if ( ! has_staging || S[e] == POPS_UNKNOWN ) continue;

      const int obsx = coda_effective_label( S[e] , opt.three_state );
      if ( obsx < 0 || obsx >= nc ) continue;

      obs.push_back( obsx );
      pred.push_back( predx );
      dur_obs[obsx]++;
      dur_obs_orig[obsx]++;

      if ( slp_lat_obs == -1 && S[e] != POPS_WAKE )
        slp_lat_obs = E[e];
      if ( rem_lat_obs == -1 && S[e] == POPS_REM )
        rem_lat_obs = E[e] - slp_lat_obs;
    }

  if ( n_pred_valid > 0 ) mean_conf /= (double)n_pred_valid;

  clocktime_t starttime( pedf ? pedf->header.starttime : "." );
  const bool file_has_times = start != NULL || stop != NULL;
  bool hms = pedf != NULL && ! file_has_times;
  if ( pedf && hms && ! starttime.valid )
    {
      logger << " ** could not find valid start-time in EDF header **\n";
      hms = false;
    }

  auto emit_times = [&]( int epoch , int e ) {
    if ( start != NULL && e >= 0 && e < (int)start->size() && (*start)[e] != "" )
      writer.value( "START" , (*start)[e] );
    if ( stop != NULL && e >= 0 && e < (int)stop->size() && (*stop)[e] != "" )
      writer.value( "STOP" , (*stop)[e] );
    if ( ! hms ) return;

    interval_t interval = pedf->timeline.epoch( epoch );

    double tp1_sec = interval.start / (double)globals::tp_1sec;
    clocktime_t present1 = starttime;
    present1.advance_seconds( tp1_sec );
    double tp1_extra = tp1_sec - (long)tp1_sec;

    double tp2_sec = interval.stop / (double)globals::tp_1sec;
    clocktime_t present2 = starttime;
    present2.advance_seconds( tp2_sec );
    double tp2_extra = tp2_sec - (long)tp2_sec;

    writer.value( "START" , present1.as_string(':') +
                  Helper::dbl2str_fixed( tp1_extra , globals::time_format_dp ).substr(1) );
    writer.value( "STOP" , present2.as_string(':') +
                  Helper::dbl2str_fixed( tp2_extra , globals::time_format_dp ).substr(1) );
  };

  auto write_epoch_block = [&]( const std::string & coda_flag ,
                                const bool level_coda ,
                                const Eigen::MatrixXd & post ,
                                const std::vector<int> & pred_vec ,
                                const std::vector<double> & conf_vec ) {
    if ( level_coda )
      writer.level( coda_flag , "CODA" );
    for (int epoch = 0; epoch < ne_total; epoch++)
      {
        writer.epoch( epoch + 1 );

        const std::map<int,int>::const_iterator ii = e2e.find( epoch );
        if ( ii == e2e.end() )
          {
            writer.value( "FLAG" , -1 );
            continue;
          }

        const int e = ii->second;
        if ( e < 0 || e >= ne )
          Helper::halt( "POPS-CODA: e2e index out of range: e=" + Helper::int2str(e) +
                        " ne=" + Helper::int2str(ne) + " epoch=" + Helper::int2str(epoch) );

        emit_times( epoch , e );

        if ( flagged[e] )
          {
            writer.value( "FLAG" , -1 );
            continue;
          }

        if ( nc == 3 )
          {
            writer.value( "PP_W"  , post(e,0) );
            writer.value( "PP_R"  , post(e,1) );
            writer.value( "PP_NR" , post(e,2) );
          }
        else
          {
            writer.value( "PP_W"  , post(e,0) );
            writer.value( "PP_R"  , post(e,1) );
            writer.value( "PP_N1" , post(e,2) );
            writer.value( "PP_N2" , post(e,3) );
            writer.value( "PP_N3" , post(e,4) );
          }

        writer.value( "CONF" , conf_vec[e] );
        writer.value( "PRED" , slabs[ pred_vec[e] ] );

        if ( has_staging )
          writer.value( "PRIOR" , pops_t::label( (pops_stage_t)S[e] ) );

        int flag = 0;
        if ( has_staging && S[e] != POPS_UNKNOWN )
          {
            const bool obs_w  = S[e] == POPS_WAKE;
            const bool obs_r  = S[e] == POPS_REM;
            const bool obs_nr = S[e] == POPS_N1 || S[e] == POPS_N2 || S[e] == POPS_N3;
            const bool prd_w  = pred_vec[e] == POPS_WAKE;
            const bool prd_r  = pred_vec[e] == POPS_REM;
            const bool prd_nr = pred_vec[e] == POPS_N1 || pred_vec[e] == POPS_N2 || pred_vec[e] == POPS_N3;
            if ( coda_effective_label( S[e] , opt.three_state ) != pred_vec[e] )
              {
                flag = 1;
                if ( (obs_w && (prd_r || prd_nr)) ||
                     (obs_r && (prd_w || prd_nr)) ||
                     (obs_nr && (prd_w || prd_r)) ) flag = 2;
              }
          }
        writer.value( "FLAG" , flag );
      }
    writer.unepoch();
    if ( level_coda )
      writer.unlevel( "CODA" );
  };

  if ( emit_stage1 )
    {
      std::vector<int> pred_stage1( ne , POPS_UNKNOWN );
      std::vector<double> conf_stage1( ne , std::numeric_limits<double>::quiet_NaN() );
      for (int e = 0; e < ne; e++)
        {
          if ( flagged[e] ) continue;
          int predx = 0;
          conf_stage1[e] = P.row(e).maxCoeff( &predx );
          pred_stage1[e] = predx;
        }
      write_epoch_block( "0" , stratify_coda , P , pred_stage1 , conf_stage1 );

      if ( stratify_coda )
        writer.level( "0" , "CODA" );

      std::vector<double> dur_pred1_stage1( nc , 0.0 );
      std::vector<double> dur_predf_stage1( nc , 0.0 );
      std::map<int,double> dur_obs_stage1, dur_obs_orig_stage1;
      std::vector<int> obs_stage1, pred_stage1_valid;
      double mean_conf_stage1 = 0.0;
      int n_pred_valid_stage1 = 0;
      int slp_lat_obs_stage1 = -1, slp_lat_prd_stage1 = -1;
      int rem_lat_obs_stage1 = -1, rem_lat_prd_stage1 = -1;

      for (int e = 0; e < ne; e++)
        {
          if ( flagged[e] ) continue;

          const int predx = pred_stage1[e];
          if ( predx == POPS_UNKNOWN ) continue;

          const double pmax = conf_stage1[e];
          mean_conf_stage1 += pmax;
          ++n_pred_valid_stage1;

          if ( predx >= 0 && predx < nc ) dur_pred1_stage1[predx]++;
          for (int ss = 0; ss < nc; ss++) dur_predf_stage1[ss] += P(e,ss);

          if ( slp_lat_prd_stage1 == -1 && predx != POPS_WAKE )
            slp_lat_prd_stage1 = E[e];
          if ( rem_lat_prd_stage1 == -1 && predx == POPS_REM )
            rem_lat_prd_stage1 = E[e] - slp_lat_prd_stage1;

          if ( ! has_staging || S[e] == POPS_UNKNOWN ) continue;

          const int obsx = coda_effective_label( S[e] , opt.three_state );
          if ( obsx < 0 || obsx >= nc ) continue;

          obs_stage1.push_back( obsx );
          pred_stage1_valid.push_back( predx );
          dur_obs_stage1[obsx]++;
          dur_obs_orig_stage1[obsx]++;

          if ( slp_lat_obs_stage1 == -1 && S[e] != POPS_WAKE )
            slp_lat_obs_stage1 = E[e];
          if ( rem_lat_obs_stage1 == -1 && S[e] == POPS_REM )
            rem_lat_obs_stage1 = E[e] - slp_lat_obs_stage1;
        }

      if ( n_pred_valid_stage1 > 0 ) mean_conf_stage1 /= (double)n_pred_valid_stage1;

      writer.value( "CONF" , mean_conf_stage1 );
      if ( slp_lat_prd_stage1 >= 0 ) writer.value( "SLP_LAT_PRD" , slp_lat_prd_stage1 * epoch_to_min );
      if ( rem_lat_prd_stage1 >= 0 ) writer.value( "REM_LAT_PRD" , rem_lat_prd_stage1 * epoch_to_min );

      if ( ! has_staging || obs_stage1.empty() )
        {
          for (int ss = 0; ss < nc; ss++)
            {
              writer.level( nc == 5 ? pops_t::label5( (pops_stage_t)ss )
                                    : pops_t::label3( (pops_stage_t)ss ) ,
                            globals::stage_strat );
              writer.value( "PRF" , dur_predf_stage1[ss] * epoch_to_min );
              writer.value( "PR1" , dur_pred1_stage1[ss] * epoch_to_min );
            }
          writer.level( pops_t::label( POPS_UNKNOWN ) , globals::stage_strat );
          writer.value( "PRF" , ( ne_total - n_pred_valid_stage1 ) * epoch_to_min );
          writer.value( "PR1" , ( ne_total - n_pred_valid_stage1 ) * epoch_to_min );
          writer.unlevel( globals::stage_strat );
          if ( stratify_coda )
            writer.unlevel( "CODA" );
        }
      else
        {
          if ( slp_lat_obs_stage1 >= 0 ) writer.value( "SLP_LAT_OBS" , slp_lat_obs_stage1 * epoch_to_min );
          if ( rem_lat_obs_stage1 >= 0 ) writer.value( "REM_LAT_OBS" , rem_lat_obs_stage1 * epoch_to_min );

          if ( nc == 5 )
            {
              pops_stats_t stats5_stage1( obs_stage1 , pred_stage1_valid , 5 );
              pops_stats_t stats3_stage1( pops_t::NRW( obs_stage1 ) , pops_t::NRW( pred_stage1_valid ) , 3 );
              writer.value( "K" , stats5_stage1.kappa );
              writer.value( "K3" , stats3_stage1.kappa );
              writer.value( "ACC" , stats5_stage1.acc );
              writer.value( "ACC3" , stats3_stage1.acc );
              writer.value( "MCC" , stats5_stage1.mcc );
              writer.value( "MCC3" , stats3_stage1.mcc );
              writer.value( "F1" , stats5_stage1.macro_f1 );
              writer.value( "PREC" , stats5_stage1.macro_precision );
              writer.value( "RECALL" , stats5_stage1.macro_recall );
              writer.value( "F1_WGT" , stats5_stage1.avg_weighted_f1 );
              writer.value( "PREC_WGT" , stats5_stage1.avg_weighted_precision );
              writer.value( "RECALL_WGT" , stats5_stage1.avg_weighted_recall );
              writer.value( "F13" , stats3_stage1.macro_f1 );
              writer.value( "PREC3" , stats3_stage1.macro_precision );
              writer.value( "RECALL3" , stats3_stage1.macro_recall );

              for (int ss = 0; ss < nc; ss++)
                {
                  writer.level( pops_t::label5( (pops_stage_t)ss ) , globals::stage_strat );
                  writer.value( "F1" , stats5_stage1.f1[ss] );
                  writer.value( "PREC" , stats5_stage1.precision[ss] );
                  writer.value( "RECALL" , stats5_stage1.recall[ss] );
                  writer.value( "OBS" , dur_obs_stage1[ss] * epoch_to_min );
                  writer.value( "ORIG" , dur_obs_orig_stage1[ss] * epoch_to_min );
                  writer.value( "PRF" , dur_predf_stage1[ss] * epoch_to_min );
                  writer.value( "PR1" , dur_pred1_stage1[ss] * epoch_to_min );
                }
              writer.level( pops_t::label( POPS_UNKNOWN ) , globals::stage_strat );
              writer.value( "OBS" , ( ne_total - (int)obs_stage1.size() ) * epoch_to_min );
              writer.value( "ORIG" , ( ne_total - (int)obs_stage1.size() ) * epoch_to_min );
              writer.value( "PRF" , ( ne_total - n_pred_valid_stage1 ) * epoch_to_min );
              writer.value( "PR1" , ( ne_total - n_pred_valid_stage1 ) * epoch_to_min );
              writer.unlevel( globals::stage_strat );
            }
          else
            {
              pops_stats_t stats3_stage1( obs_stage1 , pred_stage1_valid , 3 );
              writer.value( "K" , stats3_stage1.kappa );
              writer.value( "K3" , stats3_stage1.kappa );
              writer.value( "ACC" , stats3_stage1.acc );
              writer.value( "ACC3" , stats3_stage1.acc );
              writer.value( "MCC" , stats3_stage1.mcc );
              writer.value( "MCC3" , stats3_stage1.mcc );
              writer.value( "F1" , stats3_stage1.macro_f1 );
              writer.value( "PREC" , stats3_stage1.macro_precision );
              writer.value( "RECALL" , stats3_stage1.macro_recall );
              writer.value( "F1_WGT" , stats3_stage1.avg_weighted_f1 );
              writer.value( "PREC_WGT" , stats3_stage1.avg_weighted_precision );
              writer.value( "RECALL_WGT" , stats3_stage1.avg_weighted_recall );
              writer.value( "F13" , stats3_stage1.macro_f1 );
              writer.value( "PREC3" , stats3_stage1.macro_precision );
              writer.value( "RECALL3" , stats3_stage1.macro_recall );

              for (int ss = 0; ss < nc; ss++)
                {
                  writer.level( pops_t::label3( (pops_stage_t)ss ) , globals::stage_strat );
                  writer.value( "F1" , stats3_stage1.f1[ss] );
                  writer.value( "PREC" , stats3_stage1.precision[ss] );
                  writer.value( "RECALL" , stats3_stage1.recall[ss] );
                  writer.value( "OBS" , dur_obs_stage1[ss] * epoch_to_min );
                  writer.value( "ORIG" , dur_obs_orig_stage1[ss] * epoch_to_min );
                  writer.value( "PRF" , dur_predf_stage1[ss] * epoch_to_min );
                  writer.value( "PR1" , dur_pred1_stage1[ss] * epoch_to_min );
                }
              writer.level( pops_t::label( POPS_UNKNOWN ) , globals::stage_strat );
              writer.value( "OBS" , ( ne_total - (int)obs_stage1.size() ) * epoch_to_min );
              writer.value( "ORIG" , ( ne_total - (int)obs_stage1.size() ) * epoch_to_min );
              writer.value( "PRF" , ( ne_total - n_pred_valid_stage1 ) * epoch_to_min );
              writer.value( "PR1" , ( ne_total - n_pred_valid_stage1 ) * epoch_to_min );
              writer.unlevel( globals::stage_strat );
            }

          if ( stratify_coda )
            writer.unlevel( "CODA" );
        }
    }

  write_epoch_block( "1" , stratify_coda , P_coda , pred_stage , pred_conf );

  if ( stratify_coda )
    writer.level( "1" , "CODA" );

  if ( pops_opt_t::write_coda_posteriors_to_edf() )
    pops_posteriors::add_edf_channels( pedf ,
                                       P_coda ,
                                       E ,
                                       &flagged ,
                                       ne_total ,
                                       nc == 3 ,
                                       pops_opt_t::posterior_prefix_coda ,
                                       "POPS-CODA" );

  writer.value( "CONF" , mean_conf );

  if ( slp_lat_prd >= 0 ) writer.value( "SLP_LAT_PRD" , slp_lat_prd * epoch_to_min );
  if ( rem_lat_prd >= 0 ) writer.value( "REM_LAT_PRD" , rem_lat_prd * epoch_to_min );

  if ( ! has_staging || obs.empty() )
    {
      for (int ss = 0; ss < nc; ss++)
        {
          writer.level( nc == 5 ? pops_t::label5( (pops_stage_t)ss )
                                : pops_t::label3( (pops_stage_t)ss ) ,
                        globals::stage_strat );
          writer.value( "PRF" , dur_predf[ss] * epoch_to_min );
          writer.value( "PR1" , dur_pred1[ss] * epoch_to_min );
        }
      writer.level( pops_t::label( POPS_UNKNOWN ) , globals::stage_strat );
      writer.value( "PRF" , ( ne_total - n_pred_valid ) * epoch_to_min );
      writer.value( "PR1" , ( ne_total - n_pred_valid ) * epoch_to_min );
      writer.unlevel( globals::stage_strat );
      if ( stratify_coda )
        writer.unlevel( "CODA" );
      return;
    }

  if ( slp_lat_obs >= 0 ) writer.value( "SLP_LAT_OBS" , slp_lat_obs * epoch_to_min );
  if ( rem_lat_obs >= 0 ) writer.value( "REM_LAT_OBS" , rem_lat_obs * epoch_to_min );

  if ( nc == 5 )
    {
      pops_stats_t stats5( obs , pred , 5 );
      pops_stats_t stats3( pops_t::NRW( obs ) , pops_t::NRW( pred ) , 3 );

      writer.value( "K" , stats5.kappa );
      writer.value( "K3" , stats3.kappa );
      writer.value( "ACC" , stats5.acc );
      writer.value( "ACC3" , stats3.acc );
      writer.value( "MCC" , stats5.mcc );
      writer.value( "MCC3" , stats3.mcc );
      writer.value( "F1" , stats5.macro_f1 );
      writer.value( "PREC" , stats5.macro_precision );
      writer.value( "RECALL" , stats5.macro_recall );
      writer.value( "F1_WGT" , stats5.avg_weighted_f1 );
      writer.value( "PREC_WGT" , stats5.avg_weighted_precision );
      writer.value( "RECALL_WGT" , stats5.avg_weighted_recall );
      writer.value( "F13" , stats3.macro_f1 );
      writer.value( "PREC3" , stats3.macro_precision );
      writer.value( "RECALL3" , stats3.macro_recall );

      pops_stats_t stats_AAA( obs , pred , 5 , 1 );
      pops_stats_t stats_AAX( obs , pred , 5 , 2 );
      pops_stats_t stats_XAA( obs , pred , 5 , 3 );
      pops_stats_t stats_XAX( obs , pred , 5 , 4 );
      pops_stats_t stats_TRN( obs , pred , 5 , 5 );

      bool set_etype = false;
      if ( stats_AAA.nobs > 10 )
        {
          set_etype = true;
          writer.level( "AAA" , "ETYPE" );
          writer.value( "N" , stats_AAA.nobs );
          writer.value( "ACC" , stats_AAA.acc );
        }
      if ( stats_AAX.nobs > 10 )
        {
          set_etype = true;
          writer.level( "AAX" , "ETYPE" );
          writer.value( "N" , stats_AAX.nobs );
          writer.value( "ACC" , stats_AAX.acc );
        }
      if ( stats_XAA.nobs > 10 )
        {
          set_etype = true;
          writer.level( "XAA" , "ETYPE" );
          writer.value( "N" , stats_XAA.nobs );
          writer.value( "ACC" , stats_XAA.acc );
        }
      if ( stats_XAX.nobs > 10 )
        {
          set_etype = true;
          writer.level( "XAX" , "ETYPE" );
          writer.value( "N" , stats_XAX.nobs );
          writer.value( "ACC" , stats_XAX.acc );
        }
      if ( stats_TRN.nobs > 10 )
        {
          set_etype = true;
          writer.level( "TRN" , "ETYPE" );
          writer.value( "N" , stats_TRN.nobs );
          writer.value( "ACC" , stats_TRN.acc );
        }
      if ( set_etype )
        writer.unlevel( "ETYPE" );

      for (int ss = 0; ss < nc; ss++)
        {
          writer.level( pops_t::label5( (pops_stage_t)ss ) , globals::stage_strat );

          pops_stats_t stats_OAO( obs , pred , 5 , 0 , ss );
          pops_stats_t stats_ss_AAA( obs , pred , 5 , 1 , ss );
          pops_stats_t stats_ss_AAX( obs , pred , 5 , 2 , ss );
          pops_stats_t stats_ss_XAA( obs , pred , 5 , 3 , ss );
          pops_stats_t stats_ss_XAX( obs , pred , 5 , 4 , ss );
          pops_stats_t stats_ss_TRN( obs , pred , 5 , 5 , ss );

          bool set_ss_etype = false;
          if ( stats_OAO.nobs > 10 )
            {
              set_ss_etype = true;
              writer.level( "OAO" , "ETYPE" );
              writer.value( "N" , stats_OAO.nobs );
              writer.value( "ACC" , stats_OAO.acc );
            }
          if ( stats_ss_AAA.nobs > 10 )
            {
              set_ss_etype = true;
              writer.level( "AAA" , "ETYPE" );
              writer.value( "N" , stats_ss_AAA.nobs );
              writer.value( "ACC" , stats_ss_AAA.acc );
            }
          if ( stats_ss_AAX.nobs > 10 )
            {
              set_ss_etype = true;
              writer.level( "AAX" , "ETYPE" );
              writer.value( "N" , stats_ss_AAX.nobs );
              writer.value( "ACC" , stats_ss_AAX.acc );
            }
          if ( stats_ss_XAA.nobs > 10 )
            {
              set_ss_etype = true;
              writer.level( "XAA" , "ETYPE" );
              writer.value( "N" , stats_ss_XAA.nobs );
              writer.value( "ACC" , stats_ss_XAA.acc );
            }
          if ( stats_ss_XAX.nobs > 10 )
            {
              set_ss_etype = true;
              writer.level( "XAX" , "ETYPE" );
              writer.value( "N" , stats_ss_XAX.nobs );
              writer.value( "ACC" , stats_ss_XAX.acc );
            }
          if ( stats_ss_TRN.nobs > 10 )
            {
              set_ss_etype = true;
              writer.level( "TRN" , "ETYPE" );
              writer.value( "N" , stats_ss_TRN.nobs );
              writer.value( "ACC" , stats_ss_TRN.acc );
            }
          if ( set_ss_etype )
            writer.unlevel( "ETYPE" );

          writer.unlevel( globals::stage_strat );
        }

      for (int ss = 0; ss < nc; ss++)
        {
          writer.level( pops_t::label5( (pops_stage_t)ss ) , globals::stage_strat );
          writer.value( "F1" , stats5.f1[ss] );
          writer.value( "PREC" , stats5.precision[ss] );
          writer.value( "RECALL" , stats5.recall[ss] );
          writer.value( "OBS" , dur_obs[ss] * epoch_to_min );
          writer.value( "ORIG" , dur_obs_orig[ss] * epoch_to_min );
          writer.value( "PRF" , dur_predf[ss] * epoch_to_min );
          writer.value( "PR1" , dur_pred1[ss] * epoch_to_min );
        }
      writer.level( pops_t::label( POPS_UNKNOWN ) , globals::stage_strat );
      writer.value( "OBS" , ( ne_total - (int)obs.size() ) * epoch_to_min );
      writer.value( "ORIG" , ( ne_total - (int)obs.size() ) * epoch_to_min );
      writer.value( "PRF" , ( ne_total - n_pred_valid ) * epoch_to_min );
      writer.value( "PR1" , ( ne_total - n_pred_valid ) * epoch_to_min );
      writer.unlevel( globals::stage_strat );

      logger << "  kappa = " << stats5.kappa
             << "; 3-class kappa = " << stats3.kappa
             << " (n = " << obs.size() << " epochs)\n";
      logger << "  Confusion matrix:\n";
      pops_t::tabulate( obs , pred , true );
      logger << "\n";
      logger << "  3-class confusion matrix:\n";
      coda_print_confusion_matrix_3class( pops_t::NRW( obs ) , pops_t::NRW( pred ) );
      logger << "\n";
    }
  else
    {
      pops_stats_t stats3( obs , pred , 3 );

      writer.value( "K" , stats3.kappa );
      writer.value( "K3" , stats3.kappa );
      writer.value( "ACC" , stats3.acc );
      writer.value( "ACC3" , stats3.acc );
      writer.value( "MCC" , stats3.mcc );
      writer.value( "MCC3" , stats3.mcc );
      writer.value( "F1" , stats3.macro_f1 );
      writer.value( "PREC" , stats3.macro_precision );
      writer.value( "RECALL" , stats3.macro_recall );
      writer.value( "F1_WGT" , stats3.avg_weighted_f1 );
      writer.value( "PREC_WGT" , stats3.avg_weighted_precision );
      writer.value( "RECALL_WGT" , stats3.avg_weighted_recall );
      writer.value( "F13" , stats3.macro_f1 );
      writer.value( "PREC3" , stats3.macro_precision );
      writer.value( "RECALL3" , stats3.macro_recall );

      pops_stats_t stats_AAA( obs , pred , 3 , 1 );
      pops_stats_t stats_AAX( obs , pred , 3 , 2 );
      pops_stats_t stats_XAA( obs , pred , 3 , 3 );
      pops_stats_t stats_XAX( obs , pred , 3 , 4 );
      pops_stats_t stats_TRN( obs , pred , 3 , 5 );

      bool set_etype = false;
      if ( stats_AAA.nobs > 10 )
        {
          set_etype = true;
          writer.level( "AAA" , "ETYPE" );
          writer.value( "N" , stats_AAA.nobs );
          writer.value( "ACC" , stats_AAA.acc );
        }
      if ( stats_AAX.nobs > 10 )
        {
          set_etype = true;
          writer.level( "AAX" , "ETYPE" );
          writer.value( "N" , stats_AAX.nobs );
          writer.value( "ACC" , stats_AAX.acc );
        }
      if ( stats_XAA.nobs > 10 )
        {
          set_etype = true;
          writer.level( "XAA" , "ETYPE" );
          writer.value( "N" , stats_XAA.nobs );
          writer.value( "ACC" , stats_XAA.acc );
        }
      if ( stats_XAX.nobs > 10 )
        {
          set_etype = true;
          writer.level( "XAX" , "ETYPE" );
          writer.value( "N" , stats_XAX.nobs );
          writer.value( "ACC" , stats_XAX.acc );
        }
      if ( stats_TRN.nobs > 10 )
        {
          set_etype = true;
          writer.level( "TRN" , "ETYPE" );
          writer.value( "N" , stats_TRN.nobs );
          writer.value( "ACC" , stats_TRN.acc );
        }
      if ( set_etype )
        writer.unlevel( "ETYPE" );

      for (int ss = 0; ss < nc; ss++)
        {
          writer.level( pops_t::label3( (pops_stage_t)ss ) , globals::stage_strat );

          pops_stats_t stats_OAO( obs , pred , 3 , 0 , ss );
          pops_stats_t stats_ss_AAA( obs , pred , 3 , 1 , ss );
          pops_stats_t stats_ss_AAX( obs , pred , 3 , 2 , ss );
          pops_stats_t stats_ss_XAA( obs , pred , 3 , 3 , ss );
          pops_stats_t stats_ss_XAX( obs , pred , 3 , 4 , ss );
          pops_stats_t stats_ss_TRN( obs , pred , 3 , 5 , ss );

          bool set_ss_etype = false;
          if ( stats_OAO.nobs > 10 )
            {
              set_ss_etype = true;
              writer.level( "OAO" , "ETYPE" );
              writer.value( "N" , stats_OAO.nobs );
              writer.value( "ACC" , stats_OAO.acc );
            }
          if ( stats_ss_AAA.nobs > 10 )
            {
              set_ss_etype = true;
              writer.level( "AAA" , "ETYPE" );
              writer.value( "N" , stats_ss_AAA.nobs );
              writer.value( "ACC" , stats_ss_AAA.acc );
            }
          if ( stats_ss_AAX.nobs > 10 )
            {
              set_ss_etype = true;
              writer.level( "AAX" , "ETYPE" );
              writer.value( "N" , stats_ss_AAX.nobs );
              writer.value( "ACC" , stats_ss_AAX.acc );
            }
          if ( stats_ss_XAA.nobs > 10 )
            {
              set_ss_etype = true;
              writer.level( "XAA" , "ETYPE" );
              writer.value( "N" , stats_ss_XAA.nobs );
              writer.value( "ACC" , stats_ss_XAA.acc );
            }
          if ( stats_ss_XAX.nobs > 10 )
            {
              set_ss_etype = true;
              writer.level( "XAX" , "ETYPE" );
              writer.value( "N" , stats_ss_XAX.nobs );
              writer.value( "ACC" , stats_ss_XAX.acc );
            }
          if ( stats_ss_TRN.nobs > 10 )
            {
              set_ss_etype = true;
              writer.level( "TRN" , "ETYPE" );
              writer.value( "N" , stats_ss_TRN.nobs );
              writer.value( "ACC" , stats_ss_TRN.acc );
            }
          if ( set_ss_etype )
            writer.unlevel( "ETYPE" );

          writer.unlevel( globals::stage_strat );
        }

      for (int ss = 0; ss < nc; ss++)
        {
          writer.level( pops_t::label3( (pops_stage_t)ss ) , globals::stage_strat );
          writer.value( "F1" , stats3.f1[ss] );
          writer.value( "PREC" , stats3.precision[ss] );
          writer.value( "RECALL" , stats3.recall[ss] );
          writer.value( "OBS" , dur_obs[ss] * epoch_to_min );
          writer.value( "ORIG" , dur_obs_orig[ss] * epoch_to_min );
          writer.value( "PRF" , dur_predf[ss] * epoch_to_min );
          writer.value( "PR1" , dur_pred1[ss] * epoch_to_min );
        }
      writer.level( pops_t::label( POPS_UNKNOWN ) , globals::stage_strat );
      writer.value( "OBS" , ( ne_total - (int)obs.size() ) * epoch_to_min );
      writer.value( "ORIG" , ( ne_total - (int)obs.size() ) * epoch_to_min );
      writer.value( "PRF" , ( ne_total - n_pred_valid ) * epoch_to_min );
      writer.value( "PR1" , ( ne_total - n_pred_valid ) * epoch_to_min );
      writer.unlevel( globals::stage_strat );

      logger << "  kappa = " << stats3.kappa
             << "; 3-class kappa = " << stats3.kappa
             << " (n = " << obs.size() << " epochs)\n";
      logger << "  Confusion matrix:\n";
      pops_t::tabulate( obs , pred , true );
      logger << "\n";
    }

  if ( stratify_coda )
    writer.unlevel( "CODA" );
}

Eigen::MatrixXd pops_coda_t::rescore( const Eigen::MatrixXd & P_in ,
                                      const std::vector<int> & E ,
                                      const std::vector<bool> & flagged )
{
  if ( !lgbm.has_booster )
    Helper::halt( "POPS-CODA: no model loaded" );

  const int ne = P_in.rows();
  if ( ne == 0 ) return P_in;

  const int n_flagged = count_flagged_rows( flagged );
  if ( ne == n_flagged ) return P_in;

  Eigen::MatrixXd P = validate_and_normalise( P_in );
  const int ns = P.cols();
  const int nc = n_classes();

  if ( ns != nc )
    Helper::halt( "POPS-CODA: posterior columns (" + Helper::int2str(ns) +
                  ") do not match model classes (" + Helper::int2str(nc) + ")" );

  logger << "  POPS-CODA: " << nc << "-class prediction, "
         << ne << " valid epochs (" << n_flagged << " flagged), "
         << n_features() << " features\n";

  Eigen::MatrixXd X;
  std::vector<std::string> fn;
  make_features( P , E , flagged , X , fn );

  Eigen::MatrixXd P_coda = lgbm.predict( X );

  if ( opt.do_SHAP )
    SHAP( X , fn , E , flagged );

  if ( (int)P_coda.rows() != ne || (int)P_coda.cols() != nc )
    Helper::halt( "POPS-CODA: predict returned " + Helper::int2str(static_cast<long>(P_coda.rows())) +
                  " x " + Helper::int2str(static_cast<long>(P_coda.cols())) +
                  " but expected " + Helper::int2str(ne) + " x " + Helper::int2str(nc) );

  return P_coda;
}

void pops_coda_t::SHAP( const Eigen::MatrixXd & X ,
                        const std::vector<std::string> & fnames ,
                        const std::vector<int> & E ,
                        const std::vector<bool> & flagged )
{
  logger << "  calculating SHAP values...\n";
  Eigen::MatrixXd shap = lgbm.SHAP_values( X );

  const int nc = n_classes();
  const int n_features = fnames.size();

  int p = 0;
  for (int c = 0 ; c < nc ; c++)
    {
      if ( nc == 5 )
        writer.level( pops_t::labels5[ c ] , globals::stage_strat );
      else
        writer.level( pops_t::labels3[ c ] , globals::stage_strat );

      for (int r = 0; r < n_features ; r++)
        {
          writer.level( fnames[r] , globals::feature_strat );
          writer.value( "SHAP" , shap.col(p++).cwiseAbs().mean() );
        }
      writer.unlevel( globals::feature_strat );
    }
  writer.unlevel( globals::stage_strat );

  if ( pops_opt_t::epoch_level_SHAP )
    {
      logger << "  reporting epoch-level SHAP values...\n";

      int p = 0;
      for (int c = 0 ; c < nc ; c++)
        {
          if ( nc == 5 )
            writer.level( pops_t::labels5[ c ] , globals::stage_strat );
          else
            writer.level( pops_t::labels3[ c ] , globals::stage_strat );

          for (int r = 0; r < n_features ; r++)
            {
              writer.level( fnames[r] , globals::feature_strat );

              for (int e = 0; e < X.rows() ; e++)
                {
                  if ( flagged[e] ) continue;
                  writer.epoch( E[e] + 1 );
                  writer.value( "SHAP" , shap(e,p) );
                }
              writer.unepoch();

              p++;
            }
          writer.unlevel( globals::feature_strat );
        }

      writer.unlevel( globals::stage_strat );
    }
}


// ============================================================
//  Standalone training from a posteriors file
// ============================================================

void pops_coda_t::train_from_posteriors_file( const std::string & filename ,
                                              const std::string & config_file ,
                                              int n_iterations )
{
  X_train_rows.clear();
  S_train.clear();
  W_train.clear();
  X_valid_rows.clear();
  S_valid.clear();
  W_valid.clear();
  std::vector<coda_file_subject_t> subjects;
  load_coda_posteriors_file( filename , true , &opt , &subjects );

  if ( ! exclude_ids.empty() )
    {
      std::vector<coda_file_subject_t> filtered;
      filtered.reserve( subjects.size() );
      int n_excluded = 0;
      for (int i = 0; i < (int)subjects.size(); i++)
        {
          if ( exclude_ids.find( subjects[i].id ) != exclude_ids.end() )
            {
              ++n_excluded;
              continue;
            }
          filtered.push_back( subjects[i] );
        }
      subjects.swap( filtered );

      logger << "  POPS-CODA: excluded " << n_excluded
             << " individual" << ( n_excluded == 1 ? "" : "s" )
             << " from posteriors file " << filename << "\n";
    }

  // Build per-stage epoch thresholds (W=0,R=1,N1=2,N2=3,N3=4); size-1 applies to all
  const double row_duration_sec = opt.row_duration_sec <= 0 ? 30.0 : opt.row_duration_sec;
  auto mins_to_epochs = [row_duration_sec]( int m ) -> int {
    if ( m <= 0 ) return 0;
    const double rows = ( m * 60.0 ) / row_duration_sec;
    return (int)std::ceil( rows - 1e-9 );
  };
  std::array<int,5> min_stage_epochs;
  if ( opt.min_stage_minutes.size() == 5 )
    for (int s = 0; s < 5; s++) min_stage_epochs[s] = mins_to_epochs( opt.min_stage_minutes[s] );
  else
    min_stage_epochs.fill( mins_to_epochs( opt.min_stage_minutes[0] ) );
  const bool any_min_stage = std::any_of( min_stage_epochs.begin(), min_stage_epochs.end(),
                                          []( int v ){ return v > 0; } );

  struct ScreenCounts {
    int insufficient_contiguous = 0;
    int missing_stages = 0;
    int insufficient_stage_minutes = 0;
    int low_kappa = 0;
    int low_kappa5 = 0;
    int low_kappa3 = 0;
    int dropped_any = 0;
  } screen_counts;

  struct SubjectData {
    std::string id;
    Eigen::MatrixXd P;
    std::vector<int> S;
    std::vector<int> E;
  };

  std::vector<SubjectData> good_subjects;
  good_subjects.reserve( subjects.size() );

  for ( int subj_idx = 0 ; subj_idx < (int)subjects.size() ; subj_idx++ )
    {
      const coda_file_subject_t & src = subjects[subj_idx];
      const int ne = src.P.rows();
      SubjectData subj;
      subj.id = src.id;
      subj.P = src.P;
      subj.S = src.S;
      subj.E = src.E;

      Eigen::MatrixXd P = validate_and_normalise( subj.P );
      const int ns = P.cols();

      std::vector<bool> flagged( ne , false );
      for (int e = 0; e < ne; e++)
        {
          if ( subj.S[e] == POPS_UNKNOWN )
            { flagged[e] = true; continue; }
          for (int s = 0; s < ns; s++)
            if ( std::isnan( P(e,s) ) ) { flagged[e] = true; break; }
        }

      int longest_run = 0;
      int current_run = 0;
      int prev_valid_epoch = -1;
      std::vector<int> stage_counts( ns , 0 );
      std::vector<int> prior;
      std::vector<int> pred;
      std::vector<int> prior5;
      std::vector<int> pred5;

      for (int e = 0; e < ne; e++)
        {
          if ( flagged[e] )
            {
              current_run = 0;
              prev_valid_epoch = -1;
              continue;
            }

          if ( prev_valid_epoch >= 0 && subj.E[e] == prev_valid_epoch + 1 )
            ++current_run;
          else
            current_run = 1;

          if ( current_run > longest_run )
            longest_run = current_run;

          prev_valid_epoch = subj.E[e];

          const int label = coda_effective_label( subj.S[e] , opt.three_state );
          if ( label >= 0 && label < ns )
            stage_counts[ label ]++;

          prior.push_back( label );
          Eigen::Index predx = 0;
          P.row(e).maxCoeff( &predx );
          pred.push_back( (int)predx );

          if ( ! opt.three_state )
            {
              prior5.push_back( subj.S[e] );
              pred5.push_back( (int)predx );
            }
        }

      bool fail_contiguous = false;
      if ( opt.min_subject_contiguous_epochs > 0 &&
           longest_run < opt.min_subject_contiguous_epochs )
        {
          fail_contiguous = true;
          ++screen_counts.insufficient_contiguous;
        }

      // In broad_stage_qc mode (5-state only): collapse N1+N2+N3 → NREM for QC checks.
      // Per-stage min_stage_minutes (size==5) always overrides broad mode for duration checks.
      const bool do_broad = opt.broad_stage_qc && !opt.three_state && ns == 5;
      const bool do_broad_minutes = do_broad && ( opt.min_stage_minutes.size() != 5 );
      const int broad_nrem = do_broad
        ? stage_counts[POPS_N1] + stage_counts[POPS_N2] + stage_counts[POPS_N3]
        : 0;

      // Presence check: broad collapses N1+N2+N3 → requires any NREM > 0
      bool fail_missing_stages = false;
      if ( opt.require_all_stages )
        {
          if ( do_broad )
            {
              if ( stage_counts[POPS_WAKE] == 0 || stage_counts[POPS_REM] == 0 || broad_nrem == 0 )
                { fail_missing_stages = true; ++screen_counts.missing_stages; }
            }
          else
            {
              for (int s = 0; s < ns; s++)
                if ( stage_counts[s] == 0 )
                  { fail_missing_stages = true; ++screen_counts.missing_stages; break; }
            }
        }

      // Duration check: broad collapses NREM only when no per-stage thresholds given
      bool fail_stage_minutes = false;
      if ( any_min_stage )
        {
          if ( do_broad_minutes )
            {
              const int nrem_thresh = std::max( { min_stage_epochs[2], min_stage_epochs[3], min_stage_epochs[4] } );
              if ( stage_counts[POPS_WAKE] < min_stage_epochs[POPS_WAKE] ||
                   stage_counts[POPS_REM]  < min_stage_epochs[POPS_REM]  ||
                   broad_nrem               < nrem_thresh )
                { fail_stage_minutes = true; ++screen_counts.insufficient_stage_minutes; }
            }
          else
            {
              for (int s = 0; s < ns; s++)
                if ( min_stage_epochs[s] > 0 && stage_counts[s] < min_stage_epochs[s] )
                  { fail_stage_minutes = true; ++screen_counts.insufficient_stage_minutes; break; }
            }
        }

      const double kappa =
        opt.three_state
        ? -1.0
        : ( prior5.empty() ? -1.0 : MiscMath::kappa( prior5 , pred5 , POPS_UNKNOWN ) );
      const double kappa3 =
        prior.empty() ? -1.0 : MiscMath::kappa( pops_t::NRW( prior ) , pops_t::NRW( pred ) , POPS_UNKNOWN );

      const bool check_kappa5 = ! opt.three_state && opt.min_stage1_kappa >= 0;
      const bool check_kappa3 = opt.min_stage1_kappa3 >= 0;
      bool fail_kappa = false;
      if ( check_kappa5 && kappa < opt.min_stage1_kappa )
        {
          fail_kappa = true;
          ++screen_counts.low_kappa5;
        }
      if ( check_kappa3 && kappa3 < opt.min_stage1_kappa3 )
        {
          fail_kappa = true;
          ++screen_counts.low_kappa3;
        }
      if ( fail_kappa )
        {
          ++screen_counts.low_kappa;
        }

      if ( fail_contiguous || fail_missing_stages || fail_stage_minutes || fail_kappa )
        {
          ++screen_counts.dropped_any;
          continue;
        }

      good_subjects.push_back( subj );
    }

  logger << "  POPS-CODA: screening summary: "
         << screen_counts.dropped_any << " dropped, "
         << good_subjects.size() << " retained\n";
  logger << "    stage-1 kappa thresholds:";
  if ( opt.three_state )
    logger << " K=off";
  else
    logger << " K=" << ( opt.min_stage1_kappa >= 0 ? Helper::dbl2str( opt.min_stage1_kappa ) : "off" );
  logger << ", K3=" << ( opt.min_stage1_kappa3 >= 0 ? Helper::dbl2str( opt.min_stage1_kappa3 ) : "off" ) << "\n";
  logger << "    low stage-1 kappa: " << screen_counts.low_kappa;
  if ( ! opt.three_state )
    logger << " (K: " << screen_counts.low_kappa5;
  else
    logger << " (K: off";
  logger << ", K3: " << screen_counts.low_kappa3 << ")\n";
  logger << "    missing required stages: " << screen_counts.missing_stages << "\n";
  logger << "    below min-stage-minutes threshold: "
         << screen_counts.insufficient_stage_minutes << "\n";
  logger << "    short contiguous run: " << screen_counts.insufficient_contiguous << "\n";

  if ( good_subjects.empty() )
    Helper::halt( "POPS-CODA: no subjects remain after screening" );

  if ( opt.random_validation_subjects > 0 )
    {
      std::vector<std::string> candidates;
      for (const auto & subj : good_subjects)
        if ( ! coda_subject_in_holdouts( holdouts , subj.id ) )
          candidates.push_back( subj.id );

      const int n_available = (int)candidates.size();
      const int n_select = opt.random_validation_subjects < n_available
        ? opt.random_validation_subjects
        : n_available;

      std::vector<int> order( n_available );
      CRandom::random_draw( order );
      for (int i = 0; i < n_select; i++)
        holdouts.insert( candidates[ order[i] ] );

      logger << "  POPS-CODA: randomly selected " << n_select
             << " additional validation subject"
             << ( n_select == 1 ? "" : "s" );
      if ( n_select < opt.random_validation_subjects )
        logger << " (requested " << opt.random_validation_subjects
               << ", available " << n_available << ")";
      logger << "\n";
    }

  // Accumulate screened-good subjects into training / validation sets
  int n_subj_training = 0;
  int n_subj_validation = 0;
  for (const auto & subj : good_subjects)
    {
      const bool is_validation = coda_subject_in_holdouts( holdouts , subj.id );
      if ( accumulate( subj.P , subj.S , subj.E , subj.id ) > 0 )
        {
          if ( is_validation )
            ++n_subj_validation;
          else
            ++n_subj_training;
        }
    }

  logger << "  POPS-CODA: accumulated training data from " << n_subj_training << " subjects ("
         << (int)X_train_rows.size() << " valid epochs)";
  if ( ! holdouts.empty() )
    logger << "; validation data from " << n_subj_validation << " subjects ("
           << (int)X_valid_rows.size() << " valid epochs)";
  logger << "\n";

  if ( n_subj_training == 0 )
    Helper::halt( "POPS-CODA: no training subjects remain after screening and validation selection" );

  train_model( config_file , n_iterations );
}

void pops_coda_t::predict_from_posteriors_file( const std::string & filename )
{
  std::vector<coda_file_subject_t> subjects;
  load_coda_posteriors_file( filename , false , &opt , &subjects );

  for ( int subj_idx = 0 ; subj_idx < (int)subjects.size() ; subj_idx++ )
    {
      const coda_file_subject_t & subj = subjects[subj_idx];
      const int ne = subj.P.rows();
      if ( ne == 0 ) continue;

      logger << "  --------------------------------------------------------------------------\n";
      logger << "  POPS-CODA: subject " << subj.id
             << " (" << ( subj_idx + 1 ) << "/" << subjects.size() << ")\n";

      std::vector<bool> flagged( ne , false );
      for (int e = 0; e < ne; e++)
        for (int j = 0; j < subj.P.cols(); j++)
          if ( std::isnan( subj.P(e,j) ) ) { flagged[e] = true; break; }

      int ne_total = 0;
      for (int e = 0; e < ne; e++)
        if ( subj.E[e] + 1 > ne_total ) ne_total = subj.E[e] + 1;

      writer.id( subj.id , filename );
      predict( subj.P ,
               subj.E ,
               flagged ,
               ne_total ,
               subj.has_prior ,
               subj.S ,
               NULL ,
               subj.has_start ? &subj.start : NULL ,
               subj.has_stop ? &subj.stop : NULL ,
               opt.output_both_stages );
    }
}

#endif
