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
#include "stats/dpp-hypno.h"
#include "stats/dpp-vector.h"
#include "stats/dpp-twolevel.h"
#include "stats/glm.h"
#include "stats/matrix.h"
#include "edf/edf.h"
#include "param.h"
#include "eval.h"
#include "defs/defs.h"
#include "helper/helper.h"
#include "helper/logger.h"
#include "helper/luna_io.h"

#include "stats/statistics.h"
#include "miscmath/miscmath.h"

#include <fstream>
#include <limits>
#include <cmath>
#include <algorithm>
#include <cctype>
#include <iostream>
#include <sstream>

extern logger_t logger;

//
// dpp_fit_t : cohort-level (--dpp-fit) trainer for a single pooled
// (non-stage-conditioned) LightGBM regression booster. See the
// implementation plan's "Stage 3" section for the full design and
// reused-primitive citations (lgbm_t, dpp_io, cmd_t::pull_ivar,
// pops_t::attach_indiv_weights's add_block_weights idiom, CODA's .fnames
// sidecar pattern).
//

namespace {
  // create out_root's parent directory if it doesn't already exist, so a
  // long real training run never fails only at the final save step --
  // same globals::mkdir_command/folder_delimiter mechanism WRITE/RECORD-
  // SIZE's edf-dir= and SEDF's sedf-dir= already use (eval.cpp)
  void ensure_parent_dir_exists( const std::string & root )
  {
    const std::string expanded = Helper::expand( root );
    const size_t p = expanded.find_last_of( globals::folder_delimiter );
    if ( p == std::string::npos ) return;
    const std::string dir = expanded.substr( 0 , p );
    if ( dir == "" ) return;
    const std::string syscmd = globals::mkdir_command + " " + dir;
    system( syscmd.c_str() );
  }
}

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
  : param(param) , n_features(0) , model_set_mode(false)
{
  out_root = param.requires( "out" );
  ensure_parent_dir_exists( out_root );
  phe_label = param.requires( "phe" );
  if ( param.has( "covar" ) )
    {
      covariate_labels = Helper::parse( param.value( "covar" ) , "," );
      for (int i=0; i<(int)covariate_labels.size(); i++)
	{
	  while ( ! covariate_labels[i].empty() &&
		  std::isspace( (unsigned char)covariate_labels[i].front() ) )
	    covariate_labels[i].erase( covariate_labels[i].begin() );
	  while ( ! covariate_labels[i].empty() &&
		  std::isspace( (unsigned char)covariate_labels[i].back() ) )
	    covariate_labels[i].pop_back();
	  if ( covariate_labels[i].empty() )
	    Helper::halt( "covar= contains an empty covariate label" );
	  if ( std::find( covariate_labels.begin() , covariate_labels.begin()+i ,
			  covariate_labels[i] ) != covariate_labels.begin()+i )
	    Helper::halt( "covar= contains duplicate covariate label " + covariate_labels[i] );
	}
    }
  if ( param.has( "fit-spec" ) )
    {
      if ( param.has( "covar" ) )
	Helper::halt( "covar= and fit-spec= are mutually exclusive" );
      model_set_mode = true;
      parse_fit_spec();
    }
  else if ( param.has( "covar" ) )
    {
      model_set_mode = true;
      model_spec_t m;
      m.id = "adjusted";
      m.covariates = covariate_labels;
      m.sleep = true;
      model_specs.push_back( m );
    }
  vector_mode = dpp_vector::enabled( param );
  two_level = dpp_twolevel::enabled( param );
  if ( two_level && ! vector_mode ) Helper::halt( "DPP two-level requires vector=T" );
  stage_conditioned = param.has( "hypno" ) || param.has( "hypno-files" );
  if ( (vector_mode || two_level) && stage_conditioned )
    Helper::halt( "DPP vector mode currently uses context features in the vector row and cannot use stage-conditioned hypno corpus fitting" );
  hypno_three_state = param.has( "hypno-three-state" ) ? param.yesno( "hypno-three-state" ) : false;

  n_folds = param.has( "folds" ) ? param.requires_int( "folds" ) : ( model_set_mode ? 5 : 0 );
  save_fold_models = param.has( "folds-save" ) ? param.yesno( "folds-save" ) : false;

  if ( n_folds > 0 && param.has( "validation" ) )
    Helper::halt( "folds= and validation= are mutually exclusive -- CV already "
			 "produces a held-out prediction for every individual" );
  if ( model_set_mode && param.has( "validation" ) )
    Helper::halt( "fit-spec=/covar= model sets require grouped folds and do not use validation=" );
  if ( n_folds == 1 ) Helper::halt( "folds= must be at least 2" );

  // only meaningful during CV folds (each has its own validation set
  // attached); a harmless no-op for the final, no-validation fit
  early_stopping_rounds = param.has( "early-stopping-rounds" ) ? param.requires_int( "early-stopping-rounds" ) :
    ( n_folds > 0 ? 20 : 0 );

  forced_n_iterations = 0;   // set later by cross_validate(), if folds> 0
}

// Explicit, out-of-line destructor (member cleanup is otherwise identical
// to the implicit one -- every member destructs itself). This class is
// defined entirely in the header with no other out-of-line special member
// function, and holds an lgbm_t (RAII-owning raw LightGBM handles) both
// directly and via std::vector<lgbm_t> (stage_lgbm) -- an implicitly
// generated (and therefore implicitly *inline*, separately emitted per
// translation unit) destructor for a class this shape was observed to
// corrupt the heap intermittently (reproduced via dpp-fit's own test
// group: a later, unrelated dpp_fit_t's destructor would abort with
// "pointer being freed was not allocated"). Forcing a single, non-inline
// definition here resolved it. See also lgbm_t's own deleted copy/defined
// move constructors (lgbm/lgbm.h) -- a related, independently confirmed
// hazard (an implicit shallow copy of two lgbm_t's sharing one handle
// would double-free it) fixed alongside this.
dpp_fit_t::~dpp_fit_t() { }

void dpp_fit_t::parse_fit_spec()
{
  const std::string file = Helper::expand( param.requires( "fit-spec" ) );
  std::ifstream IN( file.c_str() );
  if ( ! IN.good() ) Helper::halt( "could not open fit-spec= " + file );

  auto labels = []( const std::string & raw ) {
    std::vector<std::string> x = Helper::parse( raw , "," );
    for (std::string & s : x)
      {
	while ( ! s.empty() && std::isspace((unsigned char)s.front()) ) s.erase(s.begin());
	while ( ! s.empty() && std::isspace((unsigned char)s.back()) ) s.pop_back();
	if ( s.empty() ) Helper::halt( "fit-spec= contains an empty covariate label" );
      }
    return x;
  };

  std::string line;
  int line_no = 0;
  while ( std::getline( IN , line ) )
    {
      ++line_no;
      const size_t comment = line.find('%');
      if ( comment != std::string::npos ) line.erase(comment);
      std::istringstream ss(line);
      std::string directive, id;
      if ( ! ( ss >> directive ) ) continue;
      if ( ! ( ss >> id ) ) Helper::halt( "fit-spec= line " + Helper::int2str(line_no) + ": missing identifier" );

      if ( directive == "MODEL" )
	{
	  model_spec_t m; m.id = id; m.sleep = false;
	  std::string tok;
	  while ( ss >> tok )
	    {
	      const size_t eq = tok.find('=');
	      if ( eq == std::string::npos ) Helper::halt( "fit-spec= line " + Helper::int2str(line_no) + ": expected key=value" );
	      const std::string key = tok.substr(0,eq), val = tok.substr(eq+1);
	      if ( key == "covar" ) m.covariates = labels(val);
	      else if ( key == "sleep" )
		{
		  if ( val == "T" || val == "1" ) m.sleep = true;
		  else if ( val == "F" || val == "0" ) m.sleep = false;
		  else Helper::halt( "fit-spec= line " + Helper::int2str(line_no) + ": sleep= must be T/F or 1/0" );
		}
	      else Helper::halt( "fit-spec= line " + Helper::int2str(line_no) + ": unknown MODEL option " + key );
	    }
	  for (const model_spec_t & z : model_specs)
	    if ( z.id == m.id ) Helper::halt( "fit-spec= duplicate MODEL " + m.id );
	  if ( ! m.sleep && m.covariates.empty() ) Helper::halt( "fit-spec= MODEL " + m.id + " has neither covariates nor sleep features" );
	  model_specs.push_back(m);
	}
      else if ( directive == "COMPARE" )
	{
	  contrast_spec_t c; c.id = id; c.base = ""; c.add = "";
	  std::string tok;
	  while ( ss >> tok )
	    {
	      const size_t eq = tok.find('=');
	      if ( eq == std::string::npos ) Helper::halt( "fit-spec= line " + Helper::int2str(line_no) + ": expected key=value" );
	      const std::string key = tok.substr(0,eq), val = tok.substr(eq+1);
	      if ( key == "base" ) c.base = val;
	      else if ( key == "add" ) c.add = val;
	      else Helper::halt( "fit-spec= line " + Helper::int2str(line_no) + ": unknown COMPARE option " + key );
	    }
	  if ( c.base == "" || c.add == "" || c.base == c.add ) Helper::halt( "fit-spec= line " + Helper::int2str(line_no) + ": invalid COMPARE models" );
	  for (const contrast_spec_t & z : contrast_specs)
	    if ( z.id == c.id ) Helper::halt( "fit-spec= duplicate COMPARE " + c.id );
	  contrast_specs.push_back(c);
	}
      else Helper::halt( "fit-spec= line " + Helper::int2str(line_no) + ": expected MODEL or COMPARE" );
    }
  IN.close();
  if ( model_specs.empty() ) Helper::halt( "fit-spec= contains no MODEL lines" );
  for (const model_spec_t & m : model_specs)
    for (const std::string & c : m.covariates)
      if ( std::find(covariate_labels.begin(),covariate_labels.end(),c) == covariate_labels.end() )
	covariate_labels.push_back(c);
  for (const contrast_spec_t & c : contrast_specs)
    {
      bool base=false, add=false;
      for (const model_spec_t & m : model_specs) { if (m.id==c.base) base=true; if (m.id==c.add) add=true; }
      if ( !base || !add ) Helper::halt( "fit-spec= COMPARE " + c.id + " references an unknown model" );
    }
}

void dpp_fit_t::fit()
{
  if ( two_level )
    {
      if ( model_set_mode )
	Helper::halt( "fit-spec=/covar= model sets are currently supported only for pooled DPP fitting" );
      dpp_twolevel::fit( param );
      return;
    }
  load_corpus();
  build_feature_labels();
  attach_phenotypes();
  attach_covariates();
  if ( stage_conditioned ) load_hypno_corpus();

  if ( model_set_mode )
    {
      fit_model_set();
      return;
    }

  if ( n_folds > 0 )
    {
      assign_folds();
      cross_validate();
    }

  // final, deployed model: always trained on the entire corpus (validation_ids
  // is empty here -- either the user never set it, or cross_validate() just
  // restored it after its own per-fold use)
  flatten_and_split();
  if ( stage_conditioned )
    {
      train_stage_boosters();

      // train_stage_boosters()'s own calibrate() calls just returned
      // {1,0} for every stage (yvalid is empty for this, the final,
      // full-corpus fit -- calibrate()'s own early-exit) -- override with
      // the CV-derived mean calibration, if available, rather than
      // shipping an uncalibrated bundle. Every subject still trained the
      // boosters above; this only changes the affine blending correction.
      if ( n_folds > 0 && ! cv_calib_a.empty() )
	{
	  calib_a = cv_calib_a;
	  calib_b = cv_calib_b;
	  logger << "  applying CV-derived calibration to the final bundle\n";
	}
    }
  else train_booster();
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
  if ( vector_mode )
    {
      full_labels.clear();
      for (int i=0; i<n_features; i++)
        full_labels.push_back( "VEC.F" + Helper::int2str(i+1) );
      return;
    }

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

void dpp_fit_t::attach_covariates()
{
  covariate.clear();
  if ( covariate_labels.empty() ) return;

  int n_numeric = 0;
  for (int i=0; i<(int)individuals.size(); i++)
    for (int c=0; c<(int)covariate_labels.size(); c++)
      {
	double x = std::numeric_limits<double>::quiet_NaN();
	if ( cmd_t::pull_ivar( individuals[i].id , covariate_labels[c] , &x ) )
	  { covariate[ individuals[i].id ][ covariate_labels[c] ] = x; ++n_numeric; }
	else
	  covariate[ individuals[i].id ][ covariate_labels[c] ] =
	    std::numeric_limits<double>::quiet_NaN();
      }

  logger << "  loaded " << covariate_labels.size() << " covariate labels ("
	 << n_numeric << " numeric subject values; missing/non-numeric retained as missing)\n";
}

void dpp_fit_t::fit_model_set()
{
  if ( stage_conditioned ) Helper::halt( "fit-spec= model sets do not support stage-conditioned DPP" );
  if ( n_folds < 2 ) Helper::halt( "fit-spec= requires folds>=2" );

  std::map<std::string,int> model_index;
  for (int i=0; i<(int)model_specs.size(); i++) model_index[ model_specs[i].id ] = i;

  struct result_t { std::map<std::string,double> pred; std::map<std::string,int> n; };
  std::vector<result_t> results( model_specs.size() );

  auto make_matrix = [&]( const model_spec_t & m , const Eigen::MatrixXd & base ,
				 const std::vector<std::pair<std::string,double> > & keys ,
				 std::vector<std::string> * labels )
    {
      const int nc = m.sleep ? base.cols() : 0;
      const int nv = m.covariates.size();
      Eigen::MatrixXd out( base.rows() , nc + nv );
      labels->clear();
      if ( m.sleep ) *labels = full_labels;
      for (const std::string & c : m.covariates)
	labels->push_back("COVAR." + c);
      for (int r=0; r<base.rows(); r++)
	{
	  for (int c=0; c<nc; c++) out(r,c) = base(r,c);
	  const std::string & id = keys[r].first;
	  for (int j=0; j<nv; j++)
	    {
	      double x = std::numeric_limits<double>::quiet_NaN();
	      std::map<std::string,std::map<std::string,double> >::const_iterator ii = covariate.find(id);
	      if ( ii != covariate.end() )
		{
		  std::map<std::string,double>::const_iterator jj = ii->second.find(m.covariates[j]);
		  if ( jj != ii->second.end() ) x = jj->second;
		}
	      // LightGBM handles NaN natively.  Do not impute from the full
	      // cohort (or even from a training fold) and do not add a synthetic
	      // missingness column here; the tree learner can choose its default
	      // direction for each split.
	      out(r,nc+j) = x;
	    }
	}
      return out;
    };

  auto train = [&]( const model_spec_t & m , const Eigen::MatrixXd & tr , const Eigen::MatrixXd & va ,
		    const std::vector<double> & yt , const std::vector<double> & yv ,
		    const std::vector<float> & wt , const std::vector<float> & wv ,
		    const std::vector<std::string> & labels )
    {
      std::unique_ptr<lgbm_t> lg(new lgbm_t());
      if ( param.has("config") ) lg->load_config(param.value("config"));
      lg->qt_mode = true;
      lg->attach_training_matrix(tr);
      if ( labels.size() == (size_t)tr.cols() ) lg->set_feature_names(labels);
      lg->attach_training_qts(yt);
      if ( va.rows() > 0 ) { lg->attach_validation_matrix(va); lg->attach_validation_qts(yv); }
      lg->training_weights = wt; lg->apply_weights(lg->training,&lg->training_weights);
      if ( va.rows() > 0 ) { lg->validation_weights = wv; lg->apply_weights(lg->validation,&lg->validation_weights); }
      if ( param.has("iterations") ) lg->n_iterations = param.requires_int("iterations");
      lg->early_stopping_rounds = early_stopping_rounds;
      lg->create_booster(param.has("verbose") && param.yesno("verbose"));
      return lg;
    };

  assign_folds();
  for (int f=0; f<n_folds; f++)
    {
      validation_ids.clear();
      for (std::map<std::string,int>::const_iterator ii=fold_assignment.begin(); ii!=fold_assignment.end(); ++ii)
	if(ii->second==f) validation_ids.insert(ii->first);
      if ( validation_ids.empty() ) continue;
      flatten_and_split();
      for (int mi=0; mi<(int)model_specs.size(); mi++)
	{
	  std::vector<std::string> labels;
	  Eigen::MatrixXd tr=make_matrix(model_specs[mi],Xtrain,train_row_key,&labels);
	  Eigen::MatrixXd va=make_matrix(model_specs[mi],Xvalid,valid_row_key,&labels);
	  std::unique_ptr<lgbm_t> lg=train(model_specs[mi],tr,va,ytrain,yvalid,indiv_weight_train,indiv_weight_valid,labels);
	  Eigen::MatrixXd p=lg->predict(va);
	  for(int r=0;r<p.rows();r++){const std::string&id=valid_row_key[r].first;results[mi].pred[id]+=p(r,0);results[mi].n[id]++;}
	}
    }

  std::ofstream CO(Helper::expand(out_root+".contrasts").c_str());
  if(!CO.good()) Helper::halt("could not write "+out_root+".contrasts");
  CO << "# DPP model-set contrasts\nNAME\tBASE\tADD\tN\tBASE_RMSE\tADD_RMSE\tDELTA_RMSE\tBASE_R2\tADD_R2\tDELTA_R2\n";
  auto metrics = [&](int mi, const std::set<std::string> & ids) {
    double sy=0,syy=0,sse=0; int n=0;
    for(const std::string&id:ids){ if(!results[mi].n.count(id))continue; double y=phenotype[id],p=results[mi].pred[id]/results[mi].n[id]; sy+=y;syy+=y*y;sse+=(y-p)*(y-p);++n; }
    const double rmse=n?std::sqrt(sse/n):std::numeric_limits<double>::quiet_NaN();
    const double mean=n?sy/n:0, tss=syy-n*mean*mean;
    const double r2=n&&tss>0?1-sse/tss:std::numeric_limits<double>::quiet_NaN();
    return std::make_pair(rmse,r2);
  };
  for(const contrast_spec_t& c:contrast_specs){std::set<std::string>ids;for(auto&x:results[model_index[c.base]].pred)if(results[model_index[c.add]].pred.count(x.first))ids.insert(x.first);auto a=metrics(model_index[c.base],ids),b=metrics(model_index[c.add],ids);CO<<c.id<<"\t"<<c.base<<"\t"<<c.add<<"\t"<<ids.size()<<"\t"<<a.first<<"\t"<<b.first<<"\t"<<a.first-b.first<<"\t"<<a.second<<"\t"<<b.second<<"\t"<<b.second-a.second<<"\n";std::ofstream PO(Helper::expand(out_root+".contrast."+c.id+".subject").c_str());if(!PO.good())Helper::halt("could not write "+out_root+".contrast."+c.id+".subject");PO<<"ID\tY_TRUE\tP_BASE\tP_ADD\tSE_BASE\tSE_ADD\tSE_BASE_MINUS_ADD\n";for(const std::string&id:ids){double y=phenotype[id],pb=results[model_index[c.base]].pred[id]/results[model_index[c.base]].n[id],pa=results[model_index[c.add]].pred[id]/results[model_index[c.add]].n[id],eb=(y-pb)*(y-pb),ea=(y-pa)*(y-pa);PO<<id<<"\t"<<y<<"\t"<<pb<<"\t"<<pa<<"\t"<<eb<<"\t"<<ea<<"\t"<<eb-ea<<"\n";}PO.close();}
  CO.close();

  for (int mi=0; mi<(int)model_specs.size(); mi++)
    {
      std::ofstream OO(Helper::expand(out_root+"."+model_specs[mi].id+".oof.subject").c_str());
      if ( !OO.good() ) Helper::halt("could not write "+out_root+"."+model_specs[mi].id+".oof.subject");
      OO << "# DPP model-set subject-level out-of-fold predictions\n"
         << "# model_id=" << model_specs[mi].id << "\n"
         << "ID\tY_TRUE\tY_PRED\n";
      for (std::map<std::string,double>::const_iterator ii=results[mi].pred.begin(); ii!=results[mi].pred.end(); ++ii)
	OO << ii->first << "\t" << phenotype[ii->first] << "\t"
	   << ii->second / results[mi].n[ii->first] << "\n";
      OO.close();
    }

  validation_ids.clear(); flatten_and_split();
  for(int mi=0;mi<(int)model_specs.size();mi++){
    std::vector<std::string> labels; Eigen::MatrixXd tr=make_matrix(model_specs[mi],Xtrain,train_row_key,&labels); Eigen::MatrixXd empty(0,tr.cols()); std::vector<double> emptyy; std::vector<float> emptyw; std::unique_ptr<lgbm_t>lg=train(model_specs[mi],tr,empty,ytrain,emptyy,indiv_weight_train,emptyw,labels); std::string root=Helper::expand(out_root+"."+model_specs[mi].id);lg->save_model(root+".mod");std::ofstream M(root+".dpp");M<<"# DPP model manifest\n# model_id="<<model_specs[mi].id<<"\n# sleep="<<(model_specs[mi].sleep?"T":"F")<<"\n# covariates=";for(int j=0;j<(int)model_specs[mi].covariates.size();j++){if(j)M<<",";M<<model_specs[mi].covariates[j];}M<<"\n# missing=LIGHTGBM_NATIVE\n# feature_names_begin\n";for(auto&x:labels)M<<x<<"\n";M.close();}
  logger << "  wrote " << model_specs.size() << " DPP model(s) and " << contrast_specs.size() << " contrast(s)\n";
}

void dpp_fit_t::load_hypno_corpus()
{
  stage_labels = dpp_hypno::stage_labels( hypno_three_state );

  std::vector<std::string> files;
  if ( param.has( "hypno-files" ) ) files = param.strvector( "hypno-files" );
  else files.push_back( param.requires( "hypno" ) );

  hypno_individuals = dpp_io::load_files( files , (int)stage_labels.size() );

  // validate alignment against the already-loaded feature corpus: same
  // individuals (by ID), same row count each -- both corpora are written
  // from the identical per-EDF output-time loop (stats/dpp.cpp), so a
  // mismatch here means the wrong hypno=/hypno-files= was supplied, not
  // something to silently paper over
  std::map<std::string,int> hypno_idx;
  for (int i=0; i<(int)hypno_individuals.size(); i++) hypno_idx[ hypno_individuals[i].id ] = i;

  for (int i=0; i<(int)individuals.size(); i++)
    {
      std::map<std::string,int>::const_iterator ii = hypno_idx.find( individuals[i].id );
      if ( ii == hypno_idx.end() )
	Helper::halt( "no hypnodensity corpus entry for individual " + individuals[i].id );
      const dpp_matrix_t & hmat = hypno_individuals[ ii->second ];
      if ( hmat.X.size() != individuals[i].X.size() )
	Helper::halt( "hypnodensity/feature row-count mismatch for individual " + individuals[i].id );

      // row count alone doesn't guarantee row *alignment* -- both corpora
      // are written from the same per-EDF output-time loop under normal
      // use, so their time_sec[] should match exactly row-for-row; if the
      // wrong hypno=/hypno-files= was supplied (e.g. from a different
      // step= run), the row count could still coincidentally match while
      // pairing each row against the wrong stage-probability vector, silently
      // corrupting stage-conditioned training -- catch that explicitly here
      for (int r=0; r<(int)individuals[i].X.size(); r++)
	if ( std::fabs( hmat.time_sec[r] - individuals[i].time_sec[r] ) > 1e-6 )
	  Helper::halt( "hypnodensity/feature time misalignment for individual " + individuals[i].id +
		       " at row " + Helper::int2str( r ) + " (" + Helper::dbl2str( hmat.time_sec[r] ) +
		       " vs " + Helper::dbl2str( individuals[i].time_sec[r] ) + ") -- wrong hypno=/hypno-files=?" );
    }

  logger << "  loaded hypnodensity corpus for " << stage_labels.size() << " stages ("
	 << Helper::stringize( stage_labels ) << ")\n";
}

void dpp_fit_t::assign_folds()
{
  fold_assignment.clear();

  if ( param.has( "fold-file" ) )
    {
      const std::string f = Helper::expand( param.value( "fold-file" ) );
      std::ifstream IN1( f.c_str() );
      if ( ! IN1.good() ) Helper::halt( "could not open fold-file= " + f );
      while ( ! IN1.eof() )
	{
	  std::string id;
	  int fold;
	  IN1 >> id;
	  if ( IN1.eof() ) break;
	  IN1 >> fold;
	  if ( id == "" ) continue;
	  if ( fold < 0 || fold >= n_folds )
	    Helper::halt( "fold-file= " + f + ": fold index " + Helper::int2str( fold ) +
			 " out of range for folds=" + Helper::int2str( n_folds ) + " (individual " + id + ")" );
	  fold_assignment[ id ] = fold;
	}
      IN1.close();

      for (int i=0; i<(int)individuals.size(); i++)
	if ( fold_assignment.find( individuals[i].id ) == fold_assignment.end() )
	  Helper::halt( "fold-file= " + f + ": no fold assignment for individual " + individuals[i].id );

      logger << "  read fold assignment for " << fold_assignment.size() << " individuals from " << f << "\n";
      return;
    }

  // default: deterministic sorted-round-robin, no hashing (portable,
  // trivially reproducible, easy to assert against in a test)
  std::vector<std::string> ids;
  for (int i=0; i<(int)individuals.size(); i++) ids.push_back( individuals[i].id );
  std::sort( ids.begin() , ids.end() );

  for (int i=0; i<(int)ids.size(); i++) fold_assignment[ ids[i] ] = i % n_folds;

  logger << "  assigned " << ids.size() << " individuals to " << n_folds << " folds (sorted round-robin)\n";
}

namespace {
  bool row_is_usable( const std::vector<double> & row )
  {
    for (int i=0; i<(int)row.size(); i++) if ( Helper::realnum( row[i] ) ) return true;
    return false;
  }

  // fractional (average) ranks, standard tie handling -- no generic rank
  // primitive exists elsewhere in this codebase to reuse; Spearman's rho
  // is just Statistics::correlation() on the two rank vectors
  std::vector<double> rank_transform( const std::vector<double> & x )
  {
    const int n = (int)x.size();
    std::vector<int> idx( n );
    for (int i=0; i<n; i++) idx[i] = i;
    std::sort( idx.begin() , idx.end() , [&x]( int a , int b ) { return x[a] < x[b]; } );

    std::vector<double> r( n );
    int i = 0;
    while ( i < n )
      {
	int j = i;
	while ( j + 1 < n && x[ idx[j+1] ] == x[ idx[i] ] ) ++j;
	const double avg_rank = ( (i+1) + (j+1) ) / 2.0;   // 1-based, averaged over the tied block
	for (int k=i; k<=j; k++) r[ idx[k] ] = avg_rank;
	i = j + 1;
      }
    return r;
  }

  double spearman( const std::vector<double> & a , const std::vector<double> & b )
  {
    return Statistics::correlation( rank_transform( a ) , rank_transform( b ) );
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

  // flatten_and_split() is re-invoked once per fold by cross_validate()
  // (stage 5), plus once more for the final full-corpus fit -- these
  // accumulators must not carry rows/boundaries over from a previous call,
  // or add_block_weights() (lgbm.cpp) would pair a stale, out-of-range
  // block boundary against the current (differently-sized) weight vector
  istart_train.clear(); istart_valid.clear();
  wtable_train.clear(); wtable_valid.clear();

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
  indiv_weight_train.assign( n_train_rows , 0.0f );
  train_row_key.resize( n_train_rows );
  // Xvalid/yvalid/indiv_weight_valid/Hvalid must be reset to empty even
  // when n_valid_rows==0 (the common final-fit case, validation_ids
  // empty) -- flatten_and_split() is re-invoked once per cross_validate()
  // fold plus once more for the final full-corpus fit, and without an
  // unconditional reset here, the final fit would silently retain the
  // *previous* (last CV fold's) validation set -- same stale-K-fold-state
  // class of bug already fixed for istart_valid/wtable_valid above, missed
  // here since valid_row_key.clear() was already unconditional but these
  // weren't
  valid_row_key.clear();
  Xvalid = Eigen::MatrixXd( 0 , n_features );
  yvalid.clear();
  indiv_weight_valid.clear();
  if ( n_valid_rows > 0 )
    {
      Xvalid = Eigen::MatrixXd::Constant( n_valid_rows , n_features , NaN_value );
      yvalid.assign( n_valid_rows , 0.0 );
      indiv_weight_valid.assign( n_valid_rows , 0.0f );
      valid_row_key.resize( n_valid_rows );
    }

  const int n_stages = (int)stage_labels.size();
  Hvalid = Eigen::MatrixXd( 0 , n_stages );
  if ( stage_conditioned )
    {
      Htrain = Eigen::MatrixXd::Constant( n_train_rows , n_stages , NaN_value );
      if ( n_valid_rows > 0 ) Hvalid = Eigen::MatrixXd::Constant( n_valid_rows , n_stages , NaN_value );
    }

  // ID -> index into hypno_individuals, for the row-aligned pull below
  // (both corpora are written from the same per-EDF loop, so individual i's
  // row r in 'individuals' corresponds exactly to row r in the matching
  // hypno_individuals entry)
  std::map<std::string,int> hypno_idx;
  if ( stage_conditioned )
    for (int i=0; i<(int)hypno_individuals.size(); i++) hypno_idx[ hypno_individuals[i].id ] = i;

  // pass 2: fill, tracking per-individual block starts for row-weighting
  int train_row = 0, valid_row = 0;
  for (int i=0; i<(int)individuals.size(); i++)
    {
      if ( kept_rows[i] == 0 ) continue;

      const bool is_valid = validation_ids.find( individuals[i].id ) != validation_ids.end();
      const double y = phenotype[ individuals[i].id ];
      const uint64_t start = is_valid ? (uint64_t)valid_row : (uint64_t)train_row;
      const float w = 1.0f / (float)kept_rows[i];

      const dpp_matrix_t * hmat = stage_conditioned ? & hypno_individuals[ hypno_idx[ individuals[i].id ] ] : NULL;

      for (int r=0; r<(int)individuals[i].X.size(); r++)
	{
	  if ( ! row_is_usable( individuals[i].X[r] ) ) continue;

	  if ( is_valid )
	    {
	      for (int c=0; c<n_features; c++) Xvalid( valid_row , c ) = individuals[i].X[r][c];
	      yvalid[ valid_row ] = y;
	      indiv_weight_valid[ valid_row ] = w;
	      valid_row_key[ valid_row ] = std::make_pair( individuals[i].id , individuals[i].time_sec[r] );
	      if ( stage_conditioned )
		for (int s=0; s<n_stages; s++) Hvalid( valid_row , s ) = hmat->X[r][s];
	      ++valid_row;
	    }
	  else
	    {
	      for (int c=0; c<n_features; c++) Xtrain( train_row , c ) = individuals[i].X[r][c];
	      ytrain[ train_row ] = y;
	      indiv_weight_train[ train_row ] = w;
	      train_row_key[ train_row ] = std::make_pair( individuals[i].id , individuals[i].time_sec[r] );
	      if ( stage_conditioned )
		for (int s=0; s<n_stages; s++) Htrain( train_row , s ) = hmat->X[r][s];
	      ++train_row;
	    }
	}

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
  // full_labels may be unpopulated/mismatched here when train_booster() is
  // called directly against a hand-built Xtrain (as some tests do, bypassing
  // load_corpus()/build_feature_labels()) -- skip naming rather than halt,
  // since the real fit() path already enforces full_labels.size()==n_features
  if ( (int)full_labels.size() == Xtrain.cols() ) lgbm.set_feature_names( full_labels );
  lgbm.attach_training_qts( ytrain );

  if ( yvalid.size() > 0 )
    {
      lgbm.attach_validation_matrix( Xvalid );
      lgbm.attach_validation_qts( yvalid );
    }

  apply_row_weights();

  // forced_n_iterations (set by cross_validate(), only for the final,
  // post-CV fit) takes precedence over iterations=/default; during CV
  // folds themselves it's still unset (0), so folds behave exactly as
  // before, governed by iterations=/default + their own early stopping
  if ( forced_n_iterations > 0 ) lgbm.n_iterations = forced_n_iterations;
  else if ( param.has( "iterations" ) ) lgbm.n_iterations = param.requires_int( "iterations" );

  // a no-op unless a validation set is attached above (true for CV folds,
  // false for the final fit) -- see lgbm_t::create_booster()
  lgbm.early_stopping_rounds = early_stopping_rounds;

  const bool verbose = param.has( "verbose" ) ? param.yesno( "verbose" ) : false;
  lgbm.create_booster( verbose );
}

void dpp_fit_t::train_stage_boosters()
{
  const bool verbose = param.has( "verbose" ) ? param.yesno( "verbose" ) : false;
  const int n_stages = (int)stage_labels.size();

  stage_lgbm.clear();
  stage_lgbm.resize( n_stages );
  calib_a.assign( n_stages , 1.0 );
  calib_b.assign( n_stages , 0.0 );

  for (int s=0; s<n_stages; s++)
    {
      // LightGBM binds one dataset to one booster -- five independent
      // stage boosters need five independent lgbm_t instances (some
      // recomputation cost vs. reusing one dataset accepted; simplest
      // correct approach given the wrapper's API shape)
      lgbm_t & lg = stage_lgbm[s];
      if ( param.has( "config" ) ) lg.load_config( param.value( "config" ) );
      lg.qt_mode = true;

      lg.attach_training_matrix( Xtrain );
      // see train_booster()'s matching comment: skip naming rather than
      // halt when full_labels doesn't match Xtrain (direct-unit-test paths)
      if ( (int)full_labels.size() == Xtrain.cols() ) lg.set_feature_names( full_labels );
      lg.attach_training_qts( ytrain );
      if ( yvalid.size() > 0 )
	{
	  lg.attach_validation_matrix( Xvalid );
	  lg.attach_validation_qts( yvalid );
	}

      // per-row weight = (1/n_i) * PP_s(row) -- built directly (bypasses
      // add_block_weights' per-individual-block-constant semantics, since
      // PP_s varies row-to-row within an individual, unlike stage 3's
      // constant-per-individual phenotype weighting). A missing/NaN PP_s
      // (e.g. an INVALID-flagged epoch) gets weight 0 -- excluded from
      // this stage's loss without needing a separate, smaller matrix.
      lg.training_weights.assign( Xtrain.rows() , 0.0f );
      for (int r=0; r<Xtrain.rows(); r++)
	{
	  const double pp = Htrain( r , s );
	  lg.training_weights[r] = Helper::realnum( pp ) ? (float)( indiv_weight_train[r] * pp ) : 0.0f;
	}
      lg.apply_weights( lg.training , &lg.training_weights );

      if ( yvalid.size() > 0 )
	{
	  lg.validation_weights.assign( Xvalid.rows() , 0.0f );
	  for (int r=0; r<Xvalid.rows(); r++)
	    {
	      const double pp = Hvalid( r , s );
	      lg.validation_weights[r] = Helper::realnum( pp ) ? (float)( indiv_weight_valid[r] * pp ) : 0.0f;
	    }
	  lg.apply_weights( lg.validation , &lg.validation_weights );
	}

      // forced_stage_n_iterations (set by cross_validate(), only for the
      // final, post-CV fit) takes precedence -- see train_booster()'s
      // matching comment for the pooled path
      if ( ! forced_stage_n_iterations.empty() ) lg.n_iterations = forced_stage_n_iterations[s];
      else if ( param.has( "iterations" ) ) lg.n_iterations = param.requires_int( "iterations" );

      lg.early_stopping_rounds = early_stopping_rounds;

      lg.create_booster( verbose );

      std::pair<double,double> ab = calibrate( s );
      calib_a[s] = ab.first;
      calib_b[s] = ab.second;

      logger << "  stage " << stage_labels[s] << ": trained, calibration a=" << calib_a[s]
	     << " b=" << calib_b[s] << "\n";
    }
}

std::pair<double,double> dpp_fit_t::calibrate( int stage_idx )
{
  if ( yvalid.size() == 0 ) return std::make_pair( 1.0 , 0.0 );

  Eigen::MatrixXd Zvalid = stage_lgbm[stage_idx].predict( Xvalid );   // nrows x 1

  std::vector<double> yy, zz;
  for (int r=0; r<(int)yvalid.size(); r++)
    {
      const double pp = Hvalid( r , stage_idx );
      if ( ! Helper::realnum( pp ) || pp <= 0.5 ) continue;
      yy.push_back( yvalid[r] );
      zz.push_back( Zvalid(r,0) );
    }

  if ( (int)yy.size() < 5 ) return std::make_pair( 1.0 , 0.0 );

  Data::Vector<double> Y( (int)yy.size() );
  Data::Matrix<double> X( (int)yy.size() , 2 );
  for (int i=0; i<(int)yy.size(); i++) { Y[i] = yy[i]; X(i,0) = 1.0; X(i,1) = zz[i]; }

  GLM g( GLM::LINEAR );
  g.set( Y , X );
  if ( ! g.fit() || ! g.valid() ) return std::make_pair( 1.0 , 0.0 );

  // display()'s per-coefficient "okay" gate (glm.cpp) requires a non-
  // degenerate *standard error*, not just a defined coefficient -- an
  // (almost never seen with real, noisy validation data) exactly-zero-
  // residual fit leaves display() reporting !all_okay, with coef left at
  // its resized-but-unset default of {0,0}. Using {0,0} directly would
  // silently zero out this stage's entire contribution to the blend
  // (a_s=0) rather than falling back to a harmless uncalibrated pass-
  // through -- check the return value and fall back to {1,0} explicitly,
  // matching every other fallback path in this function
  Data::Vector<double> coef;
  if ( ! g.display( &coef ) ) return std::make_pair( 1.0 , 0.0 );
  return std::make_pair( coef[1] , coef[0] );   // {a_s, b_s} = {slope, intercept}
}

//
// stage 5: grouped K-fold CV/OOF -- purely an evaluation layer. Reuses
// flatten_and_split()/train_booster()/train_stage_boosters() unmodified,
// simply re-driving them once per fold with validation_ids set to that
// fold's held-out individuals. See the implementation plan's "Stage 5"
// section for the full design.
//

void dpp_fit_t::cross_validate()
{
  oof_rows.clear();

  const std::set<std::string> saved_validation_ids = validation_ids;

  // per-fold outcomes, aggregated (median iteration count, mean
  // calibration) after the loop and applied to the *final*, full-corpus
  // fit -- see fit()/train_booster()/train_stage_boosters(). No fresh
  // data is held out to produce these; every fold already trains on
  // everyone but itself, so this is free signal, not a new exclusion.
  std::vector<double> fold_best_iters;
  const int n_stages_cv = stage_conditioned ? (int)stage_labels.size() : 0;
  std::vector<std::vector<double> > fold_stage_best_iters( n_stages_cv );
  std::vector<std::vector<double> > fold_calib_a( n_stages_cv ), fold_calib_b( n_stages_cv );

  for (int f=0; f<n_folds; f++)
    {
      validation_ids.clear();
      for (std::map<std::string,int>::const_iterator ii = fold_assignment.begin(); ii != fold_assignment.end(); ++ii)
	if ( ii->second == f ) validation_ids.insert( ii->first );

      if ( validation_ids.empty() )
	{
	  logger << "  fold " << f << ": no individuals assigned, skipping\n";
	  continue;
	}

      flatten_and_split();

      if ( yvalid.size() == 0 )
	{
	  logger << "  fold " << f << ": no usable held-out rows, skipping\n";
	  continue;
	}

      std::vector<double> Z( yvalid.size() , 0.0 );
      std::vector<bool> row_ok( yvalid.size() , true );

      if ( stage_conditioned )
	{
	  train_stage_boosters();   // also populates calib_a/calib_b, calibrated against this fold

	  const int n_stages = (int)stage_labels.size();
	  for (int s=0; s<n_stages; s++)
	    {
	      // best_iteration stays at its 0 default when early stopping
	      // never fired this fold -- fall back to the full configured
	      // count in that case, since 0 itself is not a meaningful
	      // iteration count to record
	      const int bi = stage_lgbm[s].best_iteration > 0 ? stage_lgbm[s].best_iteration : stage_lgbm[s].n_iterations;
	      fold_stage_best_iters[s].push_back( (double) bi );
	      fold_calib_a[s].push_back( calib_a[s] );
	      fold_calib_b[s].push_back( calib_b[s] );
	    }

	  std::vector<Eigen::MatrixXd> Ys( n_stages );
	  for (int s=0; s<n_stages; s++) Ys[s] = stage_lgbm[s].predict( Xvalid );

	  for (int r=0; r<Xvalid.rows(); r++)
	    {
	      double sum_pp = 0 , z = 0;
	      for (int s=0; s<n_stages; s++)
		{
		  const double pp = Hvalid( r , s );
		  if ( ! Helper::realnum( pp ) ) continue;
		  sum_pp += pp;
		  z += pp * ( calib_a[s] * Ys[s](r,0) + calib_b[s] );
		}
	      if ( sum_pp < 0.01 ) { row_ok[r] = false; continue; }
	      Z[r] = z;
	    }
	}
      else
	{
	  train_booster();
	  const int bi = lgbm.best_iteration > 0 ? lgbm.best_iteration : lgbm.n_iterations;
	  fold_best_iters.push_back( (double) bi );
	  Eigen::MatrixXd Y = lgbm.predict( Xvalid );
	  for (int r=0; r<Xvalid.rows(); r++) Z[r] = Y(r,0);
	}

      for (int r=0; r<Xvalid.rows(); r++)
	{
	  if ( ! row_ok[r] ) continue;
	  oof_row_t row;
	  row.id = valid_row_key[r].first;
	  row.time_sec = valid_row_key[r].second;
	  row.y_true = yvalid[r];
	  row.y_pred = Z[r];
	  row.fold = f;
	  oof_rows.push_back( row );
	}

      if ( save_fold_models )
	{
	  if ( stage_conditioned )
	    for (int s=0; s<(int)stage_labels.size(); s++)
	      stage_lgbm[s].save_model( Helper::expand( out_root + ".fold" + Helper::int2str( f ) + "." + stage_labels[s] + ".mod" ) );
	  else
	    lgbm.save_model( Helper::expand( out_root + ".fold" + Helper::int2str( f ) + ".mod" ) );
	}

      logger << "  fold " << f << ": trained on " << Xtrain.rows() << " rows, "
	     << "predicted " << Xvalid.rows() << " held-out rows\n";
    }

  validation_ids = saved_validation_ids;

  // aggregate per-fold outcomes for the final, full-corpus fit (fit(),
  // called after cross_validate() returns) -- median iteration count
  // (robust to one fold's outlier run), mean calibration (a plain
  // ensemble average over independent fold-level estimates)
  if ( ! fold_best_iters.empty() )
    {
      forced_n_iterations = (int) std::lround( MiscMath::median( fold_best_iters ) );
      logger << "  CV-derived n_iterations=" << forced_n_iterations << " (median across " << fold_best_iters.size() << " folds)\n";
    }

  if ( n_stages_cv > 0 && ! fold_stage_best_iters[0].empty() )
    {
      forced_stage_n_iterations.assign( n_stages_cv , 0 );
      cv_calib_a.assign( n_stages_cv , 1.0 );
      cv_calib_b.assign( n_stages_cv , 0.0 );
      for (int s=0; s<n_stages_cv; s++)
	{
	  forced_stage_n_iterations[s] = (int) std::lround( MiscMath::median( fold_stage_best_iters[s] ) );

	  double sum_a = 0 , sum_b = 0;
	  for (int i=0; i<(int)fold_calib_a[s].size(); i++) { sum_a += fold_calib_a[s][i]; sum_b += fold_calib_b[s][i]; }
	  cv_calib_a[s] = sum_a / (double)fold_calib_a[s].size();
	  cv_calib_b[s] = sum_b / (double)fold_calib_b[s].size();

	  logger << "  stage " << stage_labels[s] << ": CV-derived n_iterations=" << forced_stage_n_iterations[s]
		 << ", calibration a=" << cv_calib_a[s] << " b=" << cv_calib_b[s]
		 << " (mean across " << fold_calib_a[s].size() << " folds)\n";
	}
    }

  if ( oof_rows.empty() )
    {
      logger << "  *** no out-of-fold predictions produced\n";
      return;
    }

  std::vector<double> yy( oof_rows.size() ), zz( oof_rows.size() );
  double sse = 0;
  for (int i=0; i<(int)oof_rows.size(); i++)
    {
      yy[i] = oof_rows[i].y_true;
      zz[i] = oof_rows[i].y_pred;
      const double d = yy[i] - zz[i];
      sse += d * d;
    }
  const double r = Statistics::correlation( yy , zz );
  const double rmse = std::sqrt( sse / (double)oof_rows.size() );

  logger << "  cross-validation (row level): " << oof_rows.size() << " out-of-fold predictions, "
	 << "r=" << r << " rmse=" << rmse << "\n";

  write_oof();

  aggregate_oof_by_subject();
  write_oof_subject();
}

void dpp_fit_t::write_oof()
{
  const std::string oof_file = Helper::expand( out_root + ".oof" );
  std::ofstream O1( oof_file.c_str() );
  if ( ! O1.good() ) Helper::halt( "could not write " + oof_file );

  O1 << "# DPP cross-validation out-of-fold predictions\n";
  O1 << "# phe=" << phe_label << "\n";
  O1 << "# folds=" << n_folds << "\n";
  O1 << "ID\tTIME_SEC\tY_TRUE\tY_PRED\tFOLD\n";
  for (int i=0; i<(int)oof_rows.size(); i++)
    O1 << oof_rows[i].id << "\t" << oof_rows[i].time_sec << "\t"
       << oof_rows[i].y_true << "\t" << oof_rows[i].y_pred << "\t" << oof_rows[i].fold << "\n";
  O1.close();

  logger << "  wrote " << oof_rows.size() << " out-of-fold predictions to " << oof_file << "\n";
}

void dpp_fit_t::aggregate_oof_by_subject()
{
  oof_subjects.clear();

  std::map<std::string,double> sum_pred;
  std::map<std::string,int> n_rows;
  std::map<std::string,double> y_true;   // constant per individual; last write wins (all identical)

  for (int i=0; i<(int)oof_rows.size(); i++)
    {
      const std::string & id = oof_rows[i].id;
      sum_pred[ id ] += oof_rows[i].y_pred;
      n_rows[ id ] += 1;
      y_true[ id ] = oof_rows[i].y_true;
    }

  for ( std::map<std::string,double>::const_iterator ii = sum_pred.begin(); ii != sum_pred.end(); ++ii )
    {
      oof_subject_t s;
      s.id = ii->first;
      s.n_rows = n_rows[ ii->first ];
      s.y_pred_mean = ii->second / (double)s.n_rows;
      s.y_true = y_true[ ii->first ];
      oof_subjects.push_back( s );
    }
}

void dpp_fit_t::write_oof_subject()
{
  if ( oof_subjects.empty() ) return;

  // the evaluation that actually matters: the phenotype is per-subject,
  // not per-row, and the row-level r/rmse above mixes in within-subject
  // autocorrelation and unequal per-subject row counts
  std::vector<double> yy( oof_subjects.size() ), zz( oof_subjects.size() );
  double sse = 0;
  for (int i=0; i<(int)oof_subjects.size(); i++)
    {
      yy[i] = oof_subjects[i].y_true;
      zz[i] = oof_subjects[i].y_pred_mean;
      const double d = yy[i] - zz[i];
      sse += d * d;
    }
  const double r_subj   = Statistics::correlation( yy , zz );
  const double rho_subj = spearman( yy , zz );
  const double rmse_subj = std::sqrt( sse / (double)oof_subjects.size() );

  logger << "  cross-validation (subject level, N=" << oof_subjects.size() << "): "
	 << "r=" << r_subj << " rho=" << rho_subj << " rmse=" << rmse_subj << "\n";

  const std::string f = Helper::expand( out_root + ".oof.subject" );
  std::ofstream O1( f.c_str() );
  if ( ! O1.good() ) Helper::halt( "could not write " + f );

  O1 << "# DPP cross-validation out-of-fold predictions, aggregated by subject\n"
     << "# (mean Y_PRED across all of that subject's rows, vs. the true phenotype)\n"
     << "# phe=" << phe_label << "\n"
     << "# n_subjects=" << oof_subjects.size() << "\n"
     << "# pearson_r=" << r_subj << "\n"
     << "# spearman_rho=" << rho_subj << "\n"
     << "# rmse=" << rmse_subj << "\n"
     << "ID\tY_TRUE\tY_PRED_MEAN\tN_ROWS\n";
  for (int i=0; i<(int)oof_subjects.size(); i++)
    O1 << oof_subjects[i].id << "\t" << oof_subjects[i].y_true << "\t"
       << oof_subjects[i].y_pred_mean << "\t" << oof_subjects[i].n_rows << "\n";
  O1.close();

  logger << "  wrote " << oof_subjects.size() << " subject-level out-of-fold summaries to " << f << "\n";
}

namespace {
  // shared by write_importance_table() (single booster) and
  // write_aggregate_importance_table() (summed across stage boosters)
  void write_importance_rows( std::ofstream & O1 ,
			       const std::vector<std::string> & labels ,
			       const std::vector<double> & gain ,
			       const std::vector<double> & split )
  {
    std::vector<int> idx( labels.size() );
    for (int i=0; i<(int)idx.size(); i++) idx[i] = i;

    std::sort( idx.begin() , idx.end() ,
	       [&gain]( int a , int b ) { return gain[a] > gain[b]; } );

    O1 << "RANK\tFEATURE\tGAIN\tSPLIT\n";
    for (int r=0; r<(int)idx.size(); r++)
      O1 << (r+1) << "\t" << labels[ idx[r] ] << "\t" << gain[ idx[r] ] << "\t" << (int)split[ idx[r] ] << "\n";
  }
}

void dpp_fit_t::write_importance_table( lgbm_t & lg , const std::string & file )
{
  const std::vector<double> gain  = lg.feature_importance( 1 ); // C_API_FEATURE_IMPORTANCE_GAIN
  const std::vector<double> split = lg.feature_importance( 0 ); // C_API_FEATURE_IMPORTANCE_SPLIT

  if ( (int)gain.size() != (int)full_labels.size() )
    Helper::halt( "internal error: feature_importance() length does not match full_labels" );

  std::ofstream O1( file.c_str() );
  if ( ! O1.good() ) Helper::halt( "could not write " + file );

  O1 << "# DPP feature importance : " << file << "\n";
  O1 << "# sorted by GAIN (total loss-reduction attributable to splits on that feature), descending\n";
  O1 << "# SPLIT = number of times the feature was chosen as a split variable across all trees\n";
  write_importance_rows( O1 , full_labels , gain , split );
  O1.close();

  logger << "  wrote feature importance table to " << file << "\n";
}

void dpp_fit_t::write_aggregate_importance_table( const std::string & file )
{
  // stage-conditioned bundles have no single shared model -- this combines
  // the five independent stage boosters' own importances (summed gain and
  // summed split, per feature) into one overall ranking, so "which features
  // matter across the whole bundle" has a direct answer alongside the
  // per-stage <root>.<stage>.importance files
  const int n = (int)full_labels.size();
  std::vector<double> agg_gain( n , 0.0 ) , agg_split( n , 0.0 );

  for (int s=0; s<(int)stage_labels.size(); s++)
    {
      const std::vector<double> gain  = stage_lgbm[s].feature_importance( 1 );
      const std::vector<double> split = stage_lgbm[s].feature_importance( 0 );
      if ( (int)gain.size() != n )
	Helper::halt( "internal error: feature_importance() length does not match full_labels" );
      for (int i=0; i<n; i++) { agg_gain[i] += gain[i]; agg_split[i] += split[i]; }
    }

  std::ofstream O1( file.c_str() );
  if ( ! O1.good() ) Helper::halt( "could not write " + file );

  O1 << "# DPP feature importance (aggregate across stage boosters " << Helper::stringize( stage_labels ) << ") : " << file << "\n";
  O1 << "# GAIN/SPLIT summed across all stage boosters -- see <root>.<stage>.importance for per-stage values\n";
  O1 << "# sorted by summed GAIN, descending\n";
  write_importance_rows( O1 , full_labels , agg_gain , agg_split );
  O1.close();

  logger << "  wrote aggregate (cross-stage) feature importance table to " << file << "\n";
}


void dpp_fit_t::save_bundle()
{
  if ( ! stage_conditioned )
    {
      const std::string model_file = Helper::expand( out_root + ".mod" );
      lgbm.save_model( model_file );
      write_importance_table( lgbm , Helper::expand( out_root + ".importance" ) );

      const std::string manifest_file = Helper::expand( out_root + ".dpp" );
      std::ofstream O1( manifest_file.c_str() );
      if ( ! O1.good() ) Helper::halt( "could not write " + manifest_file );

      O1 << "# DPP model manifest\n";
      O1 << "# model_file=" << model_file << "\n";
      O1 << "# phe=" << phe_label << "\n";
      if ( ! covariate_labels.empty() )
	{
	  O1 << "# covariates=";
	  for (int i=0; i<(int)covariate_labels.size(); i++)
	    { if ( i ) O1 << ","; O1 << covariate_labels[i]; }
	  O1 << "\n";
	}
      O1 << "# mode=regression\n";
      if ( vector_mode )
	{
	  std::string vector_time = param.has("vector-time") ? param.value("vector-time") : "RELATIVE";
	  for (char & c : vector_time) c = (char)std::toupper((unsigned char)c);
	  O1 << "# vector=T\n";
	  O1 << "# vector-time=" << vector_time << "\n";
	}
      O1 << "# n_features=" << full_labels.size() << "\n";
      O1 << "# feature_names_begin\n";
      for (int i=0; i<(int)full_labels.size(); i++) O1 << full_labels[i] << "\n";
      O1.close();

      logger << "  wrote model to " << model_file << " and manifest to " << manifest_file << "\n";
      return;
    }

  // stage-conditioned: one .mod per stage, e.g. <root>.W.mod, <root>.N1.mod
  // (POPS's own stage label strings, so "R" not "REM" -- consistency with
  // posterior_labels(), pops/posteriors.cpp)
  for (int s=0; s<(int)stage_labels.size(); s++)
    {
      stage_lgbm[s].save_model( Helper::expand( out_root + "." + stage_labels[s] + ".mod" ) );
      write_importance_table( stage_lgbm[s] , Helper::expand( out_root + "." + stage_labels[s] + ".importance" ) );
    }

  write_aggregate_importance_table( Helper::expand( out_root + ".importance" ) );

  const std::string manifest_file = Helper::expand( out_root + ".dpp" );
  std::ofstream O1( manifest_file.c_str() );
  if ( ! O1.good() ) Helper::halt( "could not write " + manifest_file );

  O1 << "# DPP model manifest\n";
  O1 << "# phe=" << phe_label << "\n";
  if ( ! covariate_labels.empty() )
    {
      O1 << "# covariates=";
      for (int i=0; i<(int)covariate_labels.size(); i++)
	{ if ( i ) O1 << ","; O1 << covariate_labels[i]; }
      O1 << "\n";
    }
  O1 << "# mode=stage-conditioned\n";
  O1 << "# stages=" << Helper::stringize( stage_labels ) << "\n";
  O1 << "# n_features=" << full_labels.size() << "\n";
  for (int s=0; s<(int)stage_labels.size(); s++)
    O1 << "# calib." << stage_labels[s] << ".a=" << calib_a[s] << "\n"
       << "# calib." << stage_labels[s] << ".b=" << calib_b[s] << "\n";
  O1 << "# feature_names_begin\n";
  for (int i=0; i<(int)full_labels.size(); i++) O1 << full_labels[i] << "\n";
  O1.close();

  logger << "  wrote " << stage_labels.size() << " stage models and manifest to " << manifest_file << "\n";
}

//
// dpp_fit::apply() : per-EDF apply path, called from dsptools::dpp()
// (stats/dpp.cpp) when model= is set
//

namespace {

  // manifest fields beyond the feature-name list ('#'-prefixed lines,
  // CODA's .fnames convention). mode defaults to "regression" (stage 3
  // bundles never wrote a mode= line at all in earlier testing, but the
  // current writer always does -- default kept for robustness)
  struct manifest_t
  {
    std::string mode;
    std::string vector_time;
    std::vector<std::string> stages;
    std::map<std::string,double> calib_a, calib_b;
    std::vector<std::string> feature_labels;
  };

  manifest_t read_manifest( const std::string & manifest_file )
  {
    manifest_t m;
    m.mode = "regression";
    m.vector_time = "RELATIVE";

    std::ifstream IN1( manifest_file.c_str() );
    while ( 1 )
      {
	std::string line;
	Helper::safe_getline( IN1 , line );
	if ( IN1.eof() || IN1.bad() ) break;
	if ( line == "" ) continue;

	if ( line[0] == '#' )
	  {
	    if ( line.find( "# mode=" ) == 0 ) m.mode = line.substr( 7 );
	    else if ( line.find( "# vector-time=" ) == 0 ) m.vector_time = line.substr( 14 );
	    else if ( line.find( "# stages=" ) == 0 ) m.stages = Helper::parse( line.substr( 9 ) , "," );
	    else if ( line.find( "# calib." ) == 0 )
	      {
		// "# calib.<stage>.a=<value>" / "# calib.<stage>.b=<value>"
		const std::string rest = line.substr( 8 );
		const size_t eq = rest.find( '=' );
		if ( eq != std::string::npos )
		  {
		    const std::string key = rest.substr( 0 , eq );
		    double val = 0;
		    Helper::str2dbl( rest.substr( eq + 1 ) , &val );
		    if ( key.size() > 2 && key.substr( key.size() - 2 ) == ".a" )
		      m.calib_a[ key.substr( 0 , key.size() - 2 ) ] = val;
		    else if ( key.size() > 2 && key.substr( key.size() - 2 ) == ".b" )
		      m.calib_b[ key.substr( 0 , key.size() - 2 ) ] = val;
		  }
	      }
	    continue;
	  }

	m.feature_labels.push_back( line );
      }
    IN1.close();
    return m;
  }

}

void dpp_fit::apply( edf_t & edf , param_t & param , const dpp_specs_t & specs ,
		     const dpp_matrix_t & mat , const dpp_matrix_t * hmat )
{
  const std::string root = param.requires( "model" );
  const std::string model_file = Helper::expand( root + ".mod" );
  const std::string manifest_file = Helper::expand( root + ".dpp" );

  if ( ! Helper::fileExists( manifest_file ) ) Helper::halt( "could not find " + manifest_file );

  const manifest_t manifest = read_manifest( manifest_file );
  if ( param.has("vector") && param.yesno("vector") )
    {
      std::string requested = param.has("vector-time") ? param.value("vector-time") : "RELATIVE";
      for (char & c : requested) c = (char)std::toupper((unsigned char)c);
      if ( requested != "RELATIVE" && requested != "ONSET" && requested != "EDF" )
	Helper::halt( "vector-time= must be RELATIVE, ONSET, or EDF" );
      if ( requested != manifest.vector_time )
	Helper::halt( "DPP model=: vector-time mismatch: model=" + manifest.vector_time +
		      " requested=" + requested );
    }
  if ( manifest.mode == "two-level" )
    {
      dpp_twolevel::apply( edf , param , mat );
      return;
    }
  const bool stage_conditioned = manifest.mode == "stage-conditioned";

  if ( ( ! stage_conditioned ) && ! Helper::fileExists( model_file ) )
    Helper::halt( "could not find " + model_file );

  // validate feature labels against the caller's own (re-supplied) spec --
  // both to actually recompute features on the new recording, and as the
  // cross-check that it matches what the model was trained on
  std::vector<std::string> expected;
  if ( param.has("vector") && param.yesno("vector") )
    {
      expected.resize( mat.X.empty() ? 0 : mat.X[0].size() );
      for (int i=0; i<(int)expected.size(); i++) expected[i] = "VEC.F" + Helper::int2str(i+1);
    }
  else expected = dpp_fit::feature_labels( specs );
  if ( manifest.feature_labels.size() != expected.size() )
    Helper::halt( "DPP model=: feature count mismatch: loaded=" + Helper::int2str( (int)manifest.feature_labels.size() ) +
		 " expected=" + Helper::int2str( (int)expected.size() ) );
  for (int i=0; i<(int)expected.size(); i++)
    if ( manifest.feature_labels[i] != expected[i] )
      Helper::halt( "DPP model=: feature-name mismatch at column " + Helper::int2str( i + 1 ) +
		   " (loaded=" + manifest.feature_labels[i] + ", expected=" + expected[i] + ")" );

  const int n_features = (int)expected.size();

  const int nr_out = (int)mat.time_sec.size();
  if ( nr_out == 0 ) { logger << "  *** no DPP output rows to apply model to\n"; return; }

  Eigen::MatrixXd X( nr_out , n_features );
  for (int r=0; r<nr_out; r++)
    for (int c=0; c<n_features; c++)
      X(r,c) = c < (int)mat.X[r].size() ? mat.X[r][c] : std::numeric_limits<double>::quiet_NaN();

  // per-row output value and validity (a row can be marked invalid below --
  // stage-conditioned only, when the hypnodensity sum is ~0 -- meaning
  // "no stage information", not "predict zero")
  std::vector<double> Z( nr_out , 0.0 );
  std::vector<bool> row_ok( nr_out , true );

  if ( ! stage_conditioned )
    {
      lgbm_t lgbm;
      lgbm.qt_mode = true;
      lgbm.load_model( model_file );
      Eigen::MatrixXd Y = lgbm.predict( X );   // nr_out x 1, since qt_mode == true
      for (int r=0; r<nr_out; r++) Z[r] = Y(r,0);
    }
  else
    {
      if ( hmat == NULL )
	Helper::halt( "DPP model=: a stage-conditioned model requires hypno-prefix= to be supplied" );
      if ( (int)hmat->time_sec.size() != nr_out )
	Helper::halt( "DPP model=: hypnodensity row count does not match the feature output row count" );

      const int n_stages = (int)manifest.stages.size();
      if ( n_stages == 0 ) Helper::halt( "DPP model=: stage-conditioned manifest has no stages= list" );

      std::vector<lgbm_t> stage_lgbm( n_stages );
      std::vector<double> a( n_stages , 1.0 ) , b( n_stages , 0.0 );
      for (int s=0; s<n_stages; s++)
	{
	  const std::string & stg = manifest.stages[s];
	  const std::string stage_model_file = Helper::expand( root + "." + stg + ".mod" );
	  if ( ! Helper::fileExists( stage_model_file ) ) Helper::halt( "could not find " + stage_model_file );
	  stage_lgbm[s].qt_mode = true;
	  stage_lgbm[s].load_model( stage_model_file );

	  std::map<std::string,double>::const_iterator aa = manifest.calib_a.find( stg );
	  std::map<std::string,double>::const_iterator bb = manifest.calib_b.find( stg );
	  if ( aa != manifest.calib_a.end() ) a[s] = aa->second;
	  if ( bb != manifest.calib_b.end() ) b[s] = bb->second;
	}

      std::vector<Eigen::MatrixXd> Ys( n_stages );
      for (int s=0; s<n_stages; s++) Ys[s] = stage_lgbm[s].predict( X );

      for (int r=0; r<nr_out; r++)
	{
	  double sum_pp = 0 , z = 0;
	  for (int s=0; s<n_stages; s++)
	    {
	      const double pp = s < (int)hmat->X[r].size() ? hmat->X[r][s] : std::numeric_limits<double>::quiet_NaN();
	      if ( ! Helper::realnum( pp ) ) continue;
	      sum_pp += pp;
	      z += pp * ( a[s] * Ys[s](r,0) + b[s] );
	    }
	  if ( sum_pp < 0.01 ) { row_ok[r] = false; continue; }   // "no stage information", not "predict zero"
	  Z[r] = z;
	}
    }

  if ( edf.is_actually_discontinuous() )
    Helper::halt( "DPP model=: cannot attach a signal to a discontinuous EDF" );

  double step_sec = param.has( "step" ) ? param.requires_dbl( "step" ) : 30;
  if ( param.has("vector") && param.yesno("vector") )
    {
      // Normally infer the vector spacing from adjacent observations.  A
      // long-record EDF may legitimately contain only one vector sample per
      // record; in that case the record duration is the only available
      // spacing (and is exactly the intended one-sample-per-record layout).
      if ( ! param.has("step") )
        {
          if ( nr_out >= 2 ) step_sec = mat.time_sec[1] - mat.time_sec[0];
          else step_sec = edf.header.record_duration;
        }
      if ( step_sec <= 0 ) Helper::halt( "DPP vector model application found invalid time points" );
    }
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

  // pmin/pmax computed from real (non-sentinel, row_ok) values only, then a
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
      if ( ! row_ok[r] ) continue;
      if ( Z[r] < pmin ) pmin = Z[r];
      if ( Z[r] > pmax ) pmax = Z[r];
    }
  if ( pmin > pmax ) { pmin = 0; pmax = 1; }
  const double span = pmax - pmin;
  const double SENTINEL = pmin - ( span > 0 ? span : 1.0 ) - 1.0;

  std::vector<double> data( ne_total , SENTINEL );
  for (int r=0; r<nr_out; r++)
    {
      if ( ! row_ok[r] ) continue;
      const int slot = (int) std::lround( mat.time_sec[r] / step_sec ) - 1;
      if ( slot < 0 || slot >= ne_total ) continue;
      data[ slot ] = Z[r];
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
