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

#include "stats/dpp-fit.h"

#include "stats/dpp-spec.h"
#include "edf/edf.h"
#include "param.h"
#include "eval.h"
#include "helper/helper.h"
#include "helper/logger.h"
#include "helper/luna_io.h"

#include <fstream>
#include <limits>
#include <cmath>

extern logger_t logger;

//
// dpp_fit_t : cohort-level (--dpp-fit) trainer for a single pooled
// (non-stage-conditioned) LightGBM regression booster. See the
// implementation plan's "Stage 3" section for the full design and
// reused-primitive citations (lgbm_t, dpp_io, cmd_t::pull_ivar,
// pops_t::attach_indiv_weights's add_block_weights idiom, CODA's .fnames
// sidecar pattern).
//

std::vector<std::string> dpp_fit::feature_labels( const dpp_specs_t & specs )
{
  std::vector<std::string> labels;
  for (int i=0; i<(int)specs.specs.size(); i++)
    {
      const dpp_spec_t & s = specs.specs[i];
      const std::string root = s.label_root() + ".w" + Helper::dbl2str( s.window_sec );
      const int ncol = s.cols();
      for (int k=0; k<ncol; k++)
	labels.push_back( root + ".V" + ( ncol == 1 ? "" : Helper::int2str( k + 1 ) ) );
    }
  return labels;
}

dpp_fit_t::dpp_fit_t( param_t & param )
  : param(param) , n_features(0)
{
  out_root = param.requires( "out" );
  phe_label = param.requires( "phe" );
}

void dpp_fit_t::fit()
{
  load_corpus();
  build_feature_labels();
  attach_phenotypes();
  flatten_and_split();
  train_booster();
  save_bundle();
}

void dpp_fit_t::load_corpus()
{
  std::vector<std::string> files;
  if ( param.has( "files" ) ) files = param.strvector( "files" );
  else files.push_back( param.requires( "data" ) );

  individuals = dpp_io::load_files( files , -1 );
  if ( individuals.empty() ) Helper::halt( "no individuals loaded from data=/files=" );

  n_features = -1;
  for (int i=0; i<(int)individuals.size(); i++)
    for (int r=0; r<(int)individuals[i].X.size(); r++)
      {
	const int nf = (int)individuals[i].X[r].size();
	if ( n_features == -1 ) n_features = nf;
	else if ( nf != n_features )
	  Helper::halt( "inconsistent feature count within loaded corpus" );
      }
  if ( n_features <= 0 ) Helper::halt( "no valid feature rows found in the loaded corpus" );

  logger << "  loaded " << individuals.size() << " individuals, " << n_features << " features\n";
}

void dpp_fit_t::build_feature_labels()
{
  dpp_specs_t specs;

  if ( param.has( "spec" ) )
    specs.read( param.value( "spec" ) );
  else
    {
      // cohort-level: no attached EDF to resolve sig= against, so take the
      // channel names directly (dpp_specs_t::init_default() is EDF-independent)
      std::vector<std::string> channels = param.strvector( "sig" );
      specs.init_default( channels );
      specs.apply_inline_overrides( param.has( "windows" )   ? param.value( "windows" )   : "" ,
				    param.has( "filters" )   ? param.value( "filters" )   : "" ,
				    param.has( "features" )  ? param.value( "features" )  : "" ,
				    param.has( "prefilter" ) ? param.value( "prefilter" ) : "" );
    }

  full_labels = dpp_fit::feature_labels( specs );

  if ( (int)full_labels.size() != n_features )
    Helper::halt( "spec= feature count (" + Helper::int2str( (int)full_labels.size() ) +
		 ") does not match the loaded training corpus (" + Helper::int2str( n_features ) + ")" );
}

void dpp_fit_t::attach_phenotypes()
{
  for (int i=0; i<(int)individuals.size(); i++)
    {
      double y;
      if ( ! cmd_t::pull_ivar( individuals[i].id , phe_label , &y ) )
	Helper::halt( "no phenotype " + phe_label + " found for individual " + individuals[i].id );
      phenotype[ individuals[i].id ] = y;
    }
}

namespace {
  bool row_is_usable( const std::vector<double> & row )
  {
    for (int i=0; i<(int)row.size(); i++) if ( Helper::realnum( row[i] ) ) return true;
    return false;
  }
}

void dpp_fit_t::flatten_and_split()
{
  if ( param.has( "validation" ) )
    {
      std::ifstream IN1( Helper::expand( param.value( "validation" ) ).c_str() );
      if ( ! IN1.good() ) Helper::halt( "could not open " + param.value( "validation" ) );
      while ( ! IN1.eof() )
	{
	  std::string id;
	  IN1 >> id;
	  if ( IN1.eof() ) break;
	  if ( id != "" ) validation_ids.insert( id );
	}
      IN1.close();
    }

  // pass 1: count usable (not-all-NaN) rows per individual/split, so the
  // Eigen matrices can be pre-sized rather than dynamically resized
  std::vector<int> kept_rows( individuals.size() , 0 );
  int n_train_rows = 0, n_valid_rows = 0;

  for (int i=0; i<(int)individuals.size(); i++)
    {
      const bool is_valid = validation_ids.find( individuals[i].id ) != validation_ids.end();
      int cnt = 0;
      for (int r=0; r<(int)individuals[i].X.size(); r++)
	if ( row_is_usable( individuals[i].X[r] ) ) ++cnt;
      kept_rows[i] = cnt;
      if ( is_valid ) n_valid_rows += cnt; else n_train_rows += cnt;
    }

  if ( n_train_rows == 0 ) Helper::halt( "no usable training rows in the loaded corpus" );

  const double NaN_value = std::numeric_limits<double>::quiet_NaN();
  Xtrain = Eigen::MatrixXd::Constant( n_train_rows , n_features , NaN_value );
  ytrain.assign( n_train_rows , 0.0 );
  if ( n_valid_rows > 0 )
    {
      Xvalid = Eigen::MatrixXd::Constant( n_valid_rows , n_features , NaN_value );
      yvalid.assign( n_valid_rows , 0.0 );
    }

  // pass 2: fill, tracking per-individual block starts for row-weighting
  int train_row = 0, valid_row = 0;
  for (int i=0; i<(int)individuals.size(); i++)
    {
      if ( kept_rows[i] == 0 ) continue;

      const bool is_valid = validation_ids.find( individuals[i].id ) != validation_ids.end();
      const double y = phenotype[ individuals[i].id ];
      const uint64_t start = is_valid ? (uint64_t)valid_row : (uint64_t)train_row;

      for (int r=0; r<(int)individuals[i].X.size(); r++)
	{
	  if ( ! row_is_usable( individuals[i].X[r] ) ) continue;

	  if ( is_valid )
	    {
	      for (int c=0; c<n_features; c++) Xvalid( valid_row , c ) = individuals[i].X[r][c];
	      yvalid[ valid_row ] = y;
	      ++valid_row;
	    }
	  else
	    {
	      for (int c=0; c<n_features; c++) Xtrain( train_row , c ) = individuals[i].X[r][c];
	      ytrain[ train_row ] = y;
	      ++train_row;
	    }
	}

      const float w = 1.0f / (float)kept_rows[i];
      if ( is_valid ) { istart_valid.push_back( start ); wtable_valid[ start ] = w; }
      else            { istart_train.push_back( start ); wtable_train[ start ] = w; }
    }

  logger << "  " << n_train_rows << " training rows"
	 << ( n_valid_rows > 0 ? ( ", " + Helper::int2str( n_valid_rows ) + " validation rows" ) : std::string("") )
	 << " from " << individuals.size() << " individuals ("
	 << ( validation_ids.size() ) << " held out)\n";
}

void dpp_fit_t::apply_row_weights()
{
  // lgbm.attach_training_matrix()/attach_validation_matrix() already reset
  // the corresponding weight vector to 1.0 internally -- add_block_weights
  // multiplies in the per-individual 1/n_i factor, apply_weights flushes
  // the result to the dataset (lgbm/lgbm.cpp:386-406,1006-1043,970-985)
  lgbm.add_block_weights( lgbm.training , &lgbm.training_weights , istart_train , wtable_train );
  lgbm.apply_weights( lgbm.training , &lgbm.training_weights );

  if ( yvalid.size() > 0 )
    {
      lgbm.add_block_weights( lgbm.validation , &lgbm.validation_weights , istart_valid , wtable_valid );
      lgbm.apply_weights( lgbm.validation , &lgbm.validation_weights );
    }
}

void dpp_fit_t::train_booster()
{
  if ( param.has( "config" ) ) lgbm.load_config( param.value( "config" ) );

  // regression only, this stage -- must be set before create_booster() and
  // stay set through any later predict() (a mismatch found, but not
  // copied, in assoc/massoc.cpp during this stage's research)
  lgbm.qt_mode = true;

  lgbm.attach_training_matrix( Xtrain );
  lgbm.attach_training_qts( ytrain );

  if ( yvalid.size() > 0 )
    {
      lgbm.attach_validation_matrix( Xvalid );
      lgbm.attach_validation_qts( yvalid );
    }

  apply_row_weights();

  if ( param.has( "iterations" ) ) lgbm.n_iterations = param.requires_int( "iterations" );

  const bool verbose = param.has( "verbose" ) ? param.yesno( "verbose" ) : false;
  lgbm.create_booster( verbose );
}

void dpp_fit_t::save_bundle()
{
  const std::string model_file = Helper::expand( out_root + ".mod" );
  lgbm.save_model( model_file );

  const std::string manifest_file = Helper::expand( out_root + ".dpp" );
  std::ofstream O1( manifest_file.c_str() );
  if ( ! O1.good() ) Helper::halt( "could not write " + manifest_file );

  O1 << "# DPP model manifest\n";
  O1 << "# model_file=" << model_file << "\n";
  O1 << "# phe=" << phe_label << "\n";
  O1 << "# mode=regression\n";
  O1 << "# n_features=" << full_labels.size() << "\n";
  O1 << "# feature_names_begin\n";
  for (int i=0; i<(int)full_labels.size(); i++) O1 << full_labels[i] << "\n";
  O1.close();

  logger << "  wrote model to " << model_file << " and manifest to " << manifest_file << "\n";
}

//
// dpp_fit::apply() : per-EDF apply path, called from dsptools::dpp()
// (stats/dpp.cpp) when model= is set
//

void dpp_fit::apply( edf_t & edf , param_t & param , const dpp_specs_t & specs , const dpp_matrix_t & mat )
{
  const std::string root = param.requires( "model" );
  const std::string model_file = Helper::expand( root + ".mod" );
  const std::string manifest_file = Helper::expand( root + ".dpp" );

  if ( ! Helper::fileExists( model_file ) ) Helper::halt( "could not find " + model_file );
  if ( ! Helper::fileExists( manifest_file ) ) Helper::halt( "could not find " + manifest_file );

  // parse the manifest's feature-name list ('#'-prefixed lines are metadata,
  // matching CODA's .fnames convention, pops/coda.cpp)
  std::vector<std::string> loaded_labels;
  std::ifstream IN1( manifest_file.c_str() );
  while ( 1 )
    {
      std::string line;
      Helper::safe_getline( IN1 , line );
      if ( IN1.eof() || IN1.bad() ) break;
      if ( line == "" ) continue;
      if ( line[0] == '#' ) continue;
      loaded_labels.push_back( line );
    }
  IN1.close();

  // validate against the caller's own (re-supplied) spec -- both to
  // actually recompute features on the new recording, and as the
  // cross-check that it matches what the model was trained on
  std::vector<std::string> expected = dpp_fit::feature_labels( specs );
  if ( loaded_labels.size() != expected.size() )
    Helper::halt( "DPP model=: feature count mismatch: loaded=" + Helper::int2str( (int)loaded_labels.size() ) +
		 " expected=" + Helper::int2str( (int)expected.size() ) );
  for (int i=0; i<(int)expected.size(); i++)
    if ( loaded_labels[i] != expected[i] )
      Helper::halt( "DPP model=: feature-name mismatch at column " + Helper::int2str( i + 1 ) +
		   " (loaded=" + loaded_labels[i] + ", expected=" + expected[i] + ")" );

  const int n_features = (int)expected.size();

  const int nr_out = (int)mat.time_sec.size();
  if ( nr_out == 0 ) { logger << "  *** no DPP output rows to apply model to\n"; return; }

  lgbm_t lgbm;
  lgbm.qt_mode = true;
  lgbm.load_model( model_file );

  Eigen::MatrixXd X( nr_out , n_features );
  for (int r=0; r<nr_out; r++)
    for (int c=0; c<n_features; c++)
      X(r,c) = c < (int)mat.X[r].size() ? mat.X[r][c] : std::numeric_limits<double>::quiet_NaN();

  Eigen::MatrixXd Y = lgbm.predict( X );   // nr_out x 1, since qt_mode == true

  if ( edf.is_actually_discontinuous() )
    Helper::halt( "DPP model=: cannot attach a signal to a discontinuous EDF" );

  const double step_sec = param.has( "step" ) ? param.requires_dbl( "step" ) : 30;
  const double n_spr_d = edf.header.record_duration / step_sec;
  const int n_spr = (int) std::lround( n_spr_d );
  if ( std::fabs( n_spr_d - n_spr ) > 1e-6 || n_spr < 1 )
    // add_signal() needs an integer number of samples per record, i.e.
    // record_duration must itself be an integer multiple of step_sec --
    // for the common case of a sparse step= (larger than the recording's
    // native record duration), that means growing the record size to
    // match, not shrinking it
    Helper::halt( "DPP model=: record duration (" + Helper::dbl2str( edf.header.record_duration ) +
		 "s) is not an integer multiple of step=" + Helper::dbl2str( step_sec ) +
		 " -- try RECORD-SIZE dur=" + Helper::dbl2str( step_sec ) + " first" );

  const int ne_total = edf.header.nr * n_spr;

  // pmin/pmax computed from real (non-sentinel) predictions only, then a
  // detectable sentinel placed just outside that range -- passing explicit
  // bounds to add_signal (rather than letting it auto-scan, which would
  // otherwise span down to the sentinel and crush the 16-bit resolution
  // available to genuine values). Mirrors pops_indiv_t::add_edf_signals's
  // approach (pops/posteriors.cpp), generalized from POPS's fixed [0,1]
  // posterior range to an arbitrary regression output range.
  double pmin = std::numeric_limits<double>::infinity();
  double pmax = -std::numeric_limits<double>::infinity();
  for (int r=0; r<nr_out; r++)
    {
      const double v = Y(r,0);
      if ( v < pmin ) pmin = v;
      if ( v > pmax ) pmax = v;
    }
  if ( pmin > pmax ) { pmin = 0; pmax = 1; }
  const double span = pmax - pmin;
  const double SENTINEL = pmin - ( span > 0 ? span : 1.0 ) - 1.0;

  std::vector<double> data( ne_total , SENTINEL );
  for (int r=0; r<nr_out; r++)
    {
      const int slot = (int) std::lround( mat.time_sec[r] / step_sec ) - 1;
      if ( slot < 0 || slot >= ne_total ) continue;
      data[ slot ] = Y(r,0);
    }

  const std::string label = param.has( "label" ) ? param.value( "label" ) : "DPP_Z";
  if ( edf.header.has_signal( label ) )
    Helper::halt( "DPP model=: signal " + label + " already exists" );

  const double Fs_code = -(double)n_spr;
  edf.add_signal( label , Fs_code , data , SENTINEL , pmax );

  logger << "  DPP model=: attached signal " << label << " (" << ne_total << " samples, "
	 << nr_out << " real values, sentinel=" << SENTINEL << ")\n";
}

#endif
