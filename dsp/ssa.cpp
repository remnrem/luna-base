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


#include "dsp/ssa.h"

#include "param.h"
#include "edf/edf.h"
#include "edf/slice.h"
#include "helper/helper.h"
#include "helper/logger.h"
#include "db/db.h"
#include "stats/eigen_ops.h"
#include "miscmath/miscmath.h"

#include <algorithm>
#include <cmath>
#include <set>

extern logger_t logger;
extern writer_t writer;

namespace {

bool finite_vector( const Eigen::VectorXd & x )
{
  for ( int i = 0 ; i < x.size() ; i++ )
    if ( ! std::isfinite( x[i] ) ) return false;
  return true;
}

std::vector<int> parse_component_list( const std::string & s , const int d )
{
  std::vector<int> idx;
  if ( s == "" ) return idx;

  std::vector<std::string> tok = Helper::quoted_parse( s , "," );
  std::set<int> seen;

  for ( int i = 0 ; i < tok.size() ; i++ )
    {
      int x = 0;
      if ( ! Helper::str2int( Helper::unquote( tok[i] ) , &x ) )
        Helper::halt( "SSA component lists require integer values" );
      if ( x < 1 || x > d )
        Helper::halt( "SSA component index out of range" );
      if ( seen.insert( x - 1 ).second )
        idx.push_back( x - 1 );
    }

  return idx;
}

std::vector< std::vector<int> > parse_component_groups( param_t & param , const int d )
{
  std::vector< std::vector<int> > groups;
  if ( ! param.has( "grp" ) ) return groups;

  const std::vector<std::string> raw = param.strvector( "grp" , ";" );
  for ( int i = 0 ; i < raw.size() ; i++ )
    {
      const std::vector<int> idx = parse_component_list( raw[i] , d );
      if ( idx.size() == 0 )
        Helper::halt( "SSA grp requires one or more component indices in each group" );
      groups.push_back( idx );
    }

  return groups;
}

void maybe_add_signal( edf_t & edf ,
                       const std::string & label ,
                       const int sr ,
                       const Eigen::VectorXd & x )
{
  std::vector<double> copy = eigen_ops::copy_vector( x );
  edf.add_signal( label , sr , copy );
}

}


ssa_t::ssa_t()
{
  n = 0;
  l = 0;
  k = 0;
  d = 0;
}


ssa_t::ssa_t( const std::vector<double> * x , const int l1 )
{
  const int n1 = x->size();
  Eigen::VectorXd t = Eigen::VectorXd::Zero( n1 );
  for ( int i = 0 ; i < n1 ; i++ ) t[i] = (*x)[i];
  fit( t , l1 );
}


ssa_t::ssa_t( const Eigen::VectorXd & t , const int l1 )
{
  fit( t , l1 );
}


Eigen::VectorXd ssa_t::diagonal_average( const Eigen::MatrixXd & Y )
{
  const int L = Y.rows();
  const int K = Y.cols();
  const int N = L + K - 1;

  Eigen::VectorXd out = Eigen::VectorXd::Zero( N );
  Eigen::VectorXd counts = Eigen::VectorXd::Zero( N );

  for ( int i = 0 ; i < L ; i++ )
    for ( int j = 0 ; j < K ; j++ )
      {
        const int m = i + j;
        out[m] += Y(i,j);
        counts[m] += 1.0;
      }

  for ( int i = 0 ; i < N ; i++ )
    if ( counts[i] > 0 ) out[i] /= counts[i];

  return out;
}


void ssa_t::fit( const Eigen::VectorXd & t , const int l1 )
{
  original = t;
  n = t.size();

  if ( n < 4 ) Helper::halt( "SSA requires at least 4 samples" );
  if ( l1 < 2 || l1 > n / 2 )
    Helper::halt( "SSA window length L must be between 2 and n/2" );

  if ( ! finite_vector( t ) )
    Helper::halt( "SSA input contains non-finite values" );

  l = l1;
  k = n - l + 1;

  // Rough memory guard for the main working arrays used below.
  const double approx_elems = (double)l * k + (double)n * l + 2.0 * l * l + (double)k * l;
  const double approx_gb = approx_elems * sizeof(double) / 1e9;
  if ( approx_gb > 2.0 )
    Helper::halt( "SSA working set is too large; reduce L, shorten the signal, or resample first" );

  X = Eigen::MatrixXd::Zero( l , k );
  for ( int i = 0 ; i < k ; i++ )
    X.col(i) = t.segment( i , l );

  Eigen::BDCSVD<Eigen::MatrixXd> svd( X , Eigen::ComputeThinU | Eigen::ComputeThinV );
  U = svd.matrixU();
  V = svd.matrixV();
  sigma = svd.singularValues();
  d = sigma.size();
  lambda = sigma.array().square().matrix();

  TS_comps = Eigen::MatrixXd::Zero( n , d );
  for ( int i = 0 ; i < d ; i++ )
    {
      const Eigen::MatrixXd elem = sigma[i] * U.col(i) * V.col(i).transpose();
      TS_comps.col(i) = diagonal_average( elem );
    }
}


Eigen::VectorXd ssa_t::reconstruct( const std::vector<int> & indices ) const
{
  if ( indices.size() == 0 ) return Eigen::VectorXd::Zero( n );

  Eigen::VectorXd out = Eigen::VectorXd::Zero( n );
  for ( int i = 0 ; i < indices.size() ; i++ )
    {
      const int idx = indices[i];
      if ( idx < 0 || idx >= d )
        Helper::halt( "SSA reconstruction index out of range" );
      out.array() += TS_comps.col(idx).array();
    }
  return out;
}


Eigen::VectorXd ssa_t::reconstruct( const int index ) const
{
  std::vector<int> idx( 1 , index );
  return reconstruct( idx );
}


Eigen::MatrixXd ssa_t::calc_wcorr() const
{
  Eigen::VectorXd w = Eigen::VectorXd::Zero( n );

  for ( int i = 0 ; i < n ; i++ )
    {
      if ( i < l ) w[i] = i + 1;
      else if ( i < k ) w[i] = l;
      else w[i] = n - i;
    }

  Eigen::VectorXd inv_norm = Eigen::VectorXd::Zero( d );
  for ( int i = 0 ; i < d ; i++ )
    {
      const double wi = ( w.array() * TS_comps.col(i).array().square() ).sum();
      inv_norm[i] = wi > 0 ? 1.0 / std::sqrt( wi ) : 0.0;
    }

  Eigen::MatrixXd wc = Eigen::MatrixXd::Identity( d , d );
  for ( int i = 0 ; i < d ; i++ )
    for ( int j = i + 1 ; j < d ; j++ )
      {
        const double wij = ( w.array() * TS_comps.col(i).array() * TS_comps.col(j).array() ).sum();
        const double v = std::fabs( wij * inv_norm[i] * inv_norm[j] );
        wc(i,j) = v;
        wc(j,i) = v;
      }

  return wc;
}


void dsptools::ssa_wrapper( edf_t & edf , param_t & param )
{
  const std::string signal_label = param.requires( "sig" );
  const bool no_annotations = true;
  signal_list_t signals = edf.header.signal_list( signal_label , no_annotations );
  const int ns = signals.size();

  if ( ns == 0 ) Helper::halt( "no signals found matching " + signal_label );

  const bool center = param.has( "center" ) ? param.yesno( "center" ) : true;
  const bool normalize = param.has( "norm" ) ? param.yesno( "norm" ) : false;
  const bool detrend = param.yesno( "detrend" );
  const bool do_winsor = param.has( "winsor" );
  const double winsor_q = do_winsor ? param.requires_dbl( "winsor" ) : 0.0;
  const bool no_new_channels = param.has( "no-new-channels" );
  const bool want_wcorr = param.yesno( "wcorr" );
  const int wcorr_n = param.has( "wcorr-n" ) ? param.requires_int( "wcorr-n" ) : 10;
  const std::string tag = param.has( "tag" ) ? param.value( "tag" ) : "SSA_";

  if ( do_winsor && ( winsor_q < 0 || winsor_q > 0.5 ) )
    Helper::halt( "SSA winsor must be between 0 and 0.5" );
  if ( wcorr_n < 1 )
    Helper::halt( "SSA wcorr-n must be at least 1" );

  logger << "  SSA";
  if ( center ) logger << " center";
  if ( normalize ) logger << " norm";
  if ( detrend ) logger << " detrend";
  if ( do_winsor ) logger << " winsor=" << winsor_q;
  logger << "\n";

  for ( int s = 0 ; s < ns ; s++ )
    {
      const int sr = edf.header.sampling_freq( signals(s) );

      writer.level( signals.label(s) , globals::signal_strat );

      slice_t slice( edf , signals(s) , edf.timeline.wholetrace() );
      const std::vector<double> * d = slice.pdata();
      const int L = param.has( "L" ) ? param.requires_int( "L" ) :
                    param.has( "sec" ) ? (int)std::round( param.requires_dbl( "sec" ) * sr ) :
                    std::max( 2 , std::min( sr , (int)d->size() / 10 ) );
      Eigen::MatrixXd X = Eigen::MatrixXd::Zero( d->size() , 1 );
      for ( int i = 0 ; i < d->size() ; i++ ) X(i,0) = (*d)[i];

      if ( ! finite_vector( X.col(0) ) )
        Helper::halt( "SSA cannot process non-finite values in signal " + signals.label(s) );

      if ( detrend ) eigen_ops::detrend( X );
      if ( center || normalize || do_winsor )
        eigen_ops::robust_scale( X , center , normalize , winsor_q );

      ssa_t ssa( X.col(0) , L );

      writer.value( "N" , ssa.n );
      writer.value( "L" , ssa.l );
      writer.value( "K" , ssa.k );
      writer.value( "D" , ssa.d );

      const int nc = param.has( "nc" ) ? param.requires_int( "nc" ) : std::min( 10 , ssa.d );
      if ( nc < 1 ) Helper::halt( "SSA nc must be at least 1" );
      const int keep = std::min( nc , ssa.d );

      const double lambda_sum = ssa.lambda.sum();
      double cum = 0.0;

      for ( int i = 0 ; i < ssa.d ; i++ )
        {
          writer.level( i + 1 , "C" );
          writer.value( "SIGMA" , ssa.sigma[i] );
          writer.value( "LAMBDA" , ssa.lambda[i] );
          writer.value( "RC_SD" , eigen_ops::sdev( ssa.TS_comps.col(i) ) );
          writer.value( "INC" , (int)( i < keep ) );
          const double ve = lambda_sum > 0 ? ssa.lambda[i] / lambda_sum : 0.0;
          cum += ve;
          writer.value( "VE" , ve );
          writer.value( "CVE" , cum );
        }
      writer.unlevel( "C" );

      if ( want_wcorr )
        {
          const Eigen::MatrixXd wc = ssa.calc_wcorr();
          const int m = std::min( keep , std::min( wcorr_n , ssa.d ) );
          for ( int i = 0 ; i < m ; i++ )
            {
              writer.level( i + 1 , "C1" );
              for ( int j = 0 ; j < m ; j++ )
                {
                  writer.level( j + 1 , "C2" );
                  writer.value( "WCORR" , wc(i,j) );
                }
              writer.unlevel( "C2" );
            }
          writer.unlevel( "C1" );
        }

      const std::vector< std::vector<int> > groups = parse_component_groups( param , ssa.d );

      if ( ! no_new_channels )
        {
          if ( groups.size() == 0 )
            {
              logger << "  adding SSA component channels for " << signals.label(s) << ":";
              for ( int i = 0 ; i < keep ; i++ )
                {
                  const std::string ch = tag + signals.label(s) + "_" + Helper::int2str( i + 1 );
                  maybe_add_signal( edf , ch , sr , ssa.TS_comps.col(i) );
                  logger << " " << ch;
                }
              logger << "\n";
            }
          else
            {
              logger << "  adding SSA grouped reconstructions for " << signals.label(s) << ":";
              for ( int i = 0 ; i < groups.size() ; i++ )
                {
                  const std::string ch = tag + signals.label(s) + "_G" + Helper::int2str( i + 1 );
                  maybe_add_signal( edf , ch , sr , ssa.reconstruct( groups[i] ) );
                  writer.level( i + 1 , "G" );
                  writer.value( "NCOMP" , (int)groups[i].size() );
                  writer.value( "LABEL" , ch );
                  writer.unlevel( "G" );
                  logger << " " << ch;
                }
              logger << "\n";
            }
        }
      else if ( groups.size() > 0 )
        {
          for ( int i = 0 ; i < groups.size() ; i++ )
            {
              writer.level( i + 1 , "G" );
              writer.value( "NCOMP" , (int)groups[i].size() );
              writer.value( "SD" , eigen_ops::sdev( ssa.reconstruct( groups[i] ) ) );
              writer.unlevel( "G" );
            }
        }

      writer.unlevel( globals::signal_strat );
    }
}
