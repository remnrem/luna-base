
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

int coda_effective_label( int label , const bool three_state )
{
  if ( three_state && ( label == POPS_N2 || label == POPS_N3 ) )
    return POPS_N1;
  return label;
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

  std::ifstream IN( expanded_filename.c_str() );
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
      else if ( s == "PRIOR" ) col_prior = c;
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
    Helper::halt( "POPS-CODA: missing PRIOR column in " + filename );
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

  // K. short / medium soft run-lengths (past-only persistence)
  for (int s = 0; s < ns; s++) fn.push_back( "run_short_p_" + slabs[s] );
  for (int s = 0; s < ns; s++) fn.push_back( "run_med_p_"   + slabs[s] );

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
    }

  // ---- Fill each non-flagged epoch -----------------------------------
  for (int e = 0; e < ne; e++)
    {
      if ( flagged[e] ) continue;

      int col = 0;

      // A. Current posteriors
      for (int s = 0; s < ns; s++) X(e, col++) = P(e,s);

      // B. Elapsed time
      double esec = ( E[e] - first_E ) * pops_opt_t::epoch_len;
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

      // K. Short / medium soft run-lengths
      for (int s = 0; s < ns; s++) X(e, col++) = run_short[e][s];
      for (int s = 0; s < ns; s++) X(e, col++) = run_med[e][s];

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
  if ( coda_subject_in_holdouts( holdouts , subject_id ) )
    {
      X_rows = &X_valid_rows;
      S_rows = &S_valid;
    }

  // Collect rows: non-flagged only
  int n_added = 0;
  for (int e = 0; e < ne; e++)
    {
      if ( flagged[e] ) continue;

      // Convert label to CODA class index
      const int label = coda_effective_label( S[e] , opt.three_state );

      X_rows->push_back( X.row(e) );
      S_rows->push_back( label );
      ++n_added;
    }

  return n_added;
}

void pops_coda_t::load_validation_ids( const std::string & filename )
{
  holdouts.clear();

  const std::string expanded = Helper::expand( filename );
  if ( ! Helper::fileExists( expanded ) )
    Helper::halt( "POPS-CODA: could not open " + filename );

  std::ifstream in( expanded.c_str() , std::ios::in );
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
  lgbm.add_label_weights( lgbm.training , &lgbm.training_weights , lw );
  lgbm.apply_weights( lgbm.training , &lgbm.training_weights );
  if ( nrows_valid )
    {
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
  std::ofstream FN( Helper::expand( fn_file ).c_str() , std::ios::out );
  if ( !FN.good() )
    Helper::halt( "POPS-CODA: could not write feature names to " + fn_file );

  std::vector<std::string> fn = feature_names();
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
  std::ifstream FN( Helper::expand( fn_file ).c_str() , std::ios::in );
  while ( true )
    {
      std::string line;
      Helper::safe_getline( FN , line );
      if ( FN.eof() || FN.bad() ) break;
      if ( line == "" ) continue;
      loaded_fnames.push_back( line );
    }
  FN.close();

  // Validate against current options
  std::vector<std::string> expected = feature_names();
  if ( loaded_fnames.size() != expected.size() )
    logger << "  ** POPS-CODA: feature count mismatch: loaded="
           << loaded_fnames.size() << " expected=" << expected.size()
           << " — check that CODA options (3-state, context, etc.) match the model\n";

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
  if ( !lgbm.has_booster )
    Helper::halt( "POPS-CODA: no model loaded" );

  const int ne = P_in.rows();
  if ( ne == 0 ) return;

  // Adapt and validate posteriors
  Eigen::MatrixXd P = validate_and_normalise( P_in );
  const int ns = P.cols();

  const int nc = n_classes();
  if ( ns != nc )
    Helper::halt( "POPS-CODA: posterior columns (" + Helper::int2str(ns) +
                  ") do not match model classes (" + Helper::int2str(nc) + ")" );

  // Count flagged epochs for logging
  int n_flagged = 0;
  for (int e = 0; e < ne; e++) if ( flagged[e] ) n_flagged++;

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
  const double epoch_to_min = pops_opt_t::epoch_len / 60.0;
  const bool use_mdl_strata = emit_stage1;

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

  auto write_epoch_block = [&]( const std::string & mdl ,
                                const Eigen::MatrixXd & post ,
                                const std::vector<int> & pred_vec ,
                                const std::vector<double> & conf_vec ) {
    if ( use_mdl_strata )
      writer.level( mdl , "MDL" );
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
      write_epoch_block( "E" , P , pred_stage1 , conf_stage1 );
      if ( use_mdl_strata )
        writer.unlevel( "MDL" );
    }

  write_epoch_block( "CODA" , P_coda , pred_stage , pred_conf );

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
      if ( use_mdl_strata )
        writer.unlevel( "MDL" );
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
      pops_t::tabulate( pops_t::NRW( obs ) , pops_t::NRW( pred ) , true );
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

  if ( use_mdl_strata )
    writer.unlevel( "MDL" );
}

Eigen::MatrixXd pops_coda_t::rescore( const Eigen::MatrixXd & P_in ,
                                      const std::vector<int> & E ,
                                      const std::vector<bool> & flagged )
{
  if ( !lgbm.has_booster )
    Helper::halt( "POPS-CODA: no model loaded" );

  const int ne = P_in.rows();
  if ( ne == 0 ) return P_in;

  Eigen::MatrixXd P = validate_and_normalise( P_in );
  const int ns = P.cols();
  const int nc = n_classes();

  if ( ns != nc )
    Helper::halt( "POPS-CODA: posterior columns (" + Helper::int2str(ns) +
                  ") do not match model classes (" + Helper::int2str(nc) + ")" );

  int n_flagged = 0;
  for (int e = 0; e < ne; e++) if ( flagged[e] ) n_flagged++;

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
  X_valid_rows.clear();
  S_valid.clear();
  std::vector<coda_file_subject_t> subjects;
  load_coda_posteriors_file( filename , true , &opt , &subjects );

  // Build per-stage epoch thresholds (W=0,R=1,N1=2,N2=3,N3=4); size-1 applies to all
  auto mins_to_epochs = []( int m ) -> int {
    return m <= 0 ? 0 : ( m * 60 + 29 ) / 30;
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

      const double kappa = prior.empty() ? -1.0 : MiscMath::kappa( prior , pred , POPS_UNKNOWN );
      bool fail_kappa = false;
      if ( kappa < opt.min_stage1_kappa )
        {
          fail_kappa = true;
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
  logger << "    low stage-1 kappa: " << screen_counts.low_kappa << "\n";
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
               false );
    }
}

#endif
