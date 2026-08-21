//    --------------------------------------------------------------------
//    Two-level vector DPP implementation.
//    --------------------------------------------------------------------

#ifdef HAS_LGBM

#include "stats/dpp-twolevel.h"
#include "stats/dpp-vector.h"
#include "stats/dpp-io.h"
#include "lgbm/lgbm.h"
#include "edf/edf.h"
#include "param.h"
#include "eval.h"
#include "defs/defs.h"
#include "helper/helper.h"
#include "helper/logger.h"
#include "stats/Eigen/Dense"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <set>
#include <sstream>

extern logger_t logger;

namespace {

  bool real( const double x ) { return std::isfinite(x); }

  std::vector<std::string> groups( const param_t & p )
  {
    std::vector<std::string> g = p.has("level2-features")
      ? p.strvector("level2-features")
      : std::vector<std::string>{"BASE","VAR","STAGE","SCORE","TIME","GEOM"};
    for (std::string & s : g)
      for (char & c : s) c = (char)std::toupper((unsigned char)c);
    return g;
  }

  bool has( const std::vector<std::string> & g , const std::string & x )
  { return std::find(g.begin(),g.end(),x) != g.end(); }

  bool usable( const std::vector<double> & x )
  { for (double v : x) if (real(v)) return true; return false; }

  double at( const std::vector<double> & x , const int i )
  { return i >= 0 && i < (int)x.size() ? x[i] : std::numeric_limits<double>::quiet_NaN(); }

  struct summary_t { std::vector<double> x; };

  std::vector<std::string> summary_labels( const int d , const param_t & p )
  {
    const dpp_vector::layout_t l = dpp_vector::layout(d,p);
    const std::vector<std::string> g = groups(p);
    std::vector<std::string> z;
    if (has(g,"BASE")) for(int j=0;j<d;j++) z.push_back("L2.BASE.MEAN.V"+Helper::int2str(j+1));
    if (has(g,"VAR")) for(int j=0;j<d;j++) z.push_back("L2.VAR.SD.V"+Helper::int2str(j+1));
    if (has(g,"STAGE"))
      for (const std::string & s : {std::string("NREM"),std::string("REM")})
        for(int j=0;j<d;j++) z.push_back("L2.STAGE."+s+".MEAN.V"+Helper::int2str(j+1));
    if (has(g,"SCORE"))
      for (const std::string & s : {std::string("MEAN"),std::string("SD"),std::string("Q90"),
                                    std::string("EARLY"),std::string("LATE"),std::string("EARLY_LATE"),
                                    std::string("NREM"),std::string("REM")}) z.push_back("L2.SCORE."+s);
    if (has(g,"TIME"))
      for (const std::string & s : {std::string("MEAN_FRAC"),std::string("SD_FRAC"),
                                    std::string("DURATION_H"),std::string("VALID_FRAC")}) z.push_back("L2.TIME."+s);
    if (has(g,"GEOM"))
      {
        const std::vector<std::string> names={"NORM","BASELINE_DIST","BASELINE_COS","PREV_DIST","PREV_COS",
                                              "VELOCITY","ACCELERATION","NORM_DELTA","TURN_ANGLE"};
        for (const std::string & s:names) { z.push_back("L2.GEOM."+s+".MEAN"); z.push_back("L2.GEOM."+s+".SD"); }
        z.push_back("L2.GEOM.PATH_LENGTH");
        z.push_back("L2.GEOM.EARLY_LATE_CENTROID_DIST");
      }
    (void)l;
    return z;
  }

  summary_t summarize( const dpp_matrix_t & m , const std::vector<double> & score,
                       const dpp_vector::layout_t & l , const param_t & p )
  {
    const int d=l.embedding_dim, n=(int)m.X.size();
    const std::vector<std::string> g=groups(p);
    summary_t s;
    if (has(g,"STAGE") && !l.context)
      Helper::halt("DPP two-level STAGE summaries require CONTEXT in vector-features");

    std::vector<std::vector<double> > raw(d), early_raw(d), late_raw(d), nrem_raw(d), rem_raw(d);
    for(int j=0;j<d;j++) { raw[j].reserve(n); early_raw[j].reserve(n); late_raw[j].reserve(n); }
    std::vector<double> score_nrem,score_rem,score_early,score_late, times, valid;
    std::vector<std::vector<double> > traj(9);
    double first_t=std::numeric_limits<double>::infinity(), last_t=-std::numeric_limits<double>::infinity();
    for(int r=0;r<n;r++)
      {
        const std::vector<double> & row=m.X[r];
        const double frac=l.context ? at(row,l.context_offset) : (n>1?(double)r/(n-1):0);
        const bool early=real(frac)&&frac<0.5, late=real(frac)&&frac>=0.5;
        const bool nrem=l.context && (at(row,l.context_offset+2)>0.5 || at(row,l.context_offset+3)>0.5 || at(row,l.context_offset+4)>0.5);
        const bool rem=l.context && at(row,l.context_offset+5)>0.5;
        if (real(m.time_sec[r])) { first_t=std::min(first_t,m.time_sec[r]); last_t=std::max(last_t,m.time_sec[r]); }
        if(real(frac)) times.push_back(frac);
        int nv=0;
        for(int j=0;j<d;j++)
          {
            const double v=at(row,l.raw_offset+j); if(real(v)) { raw[j].push_back(v); if(early) early_raw[j].push_back(v); if(late) late_raw[j].push_back(v); if(nrem)nrem_raw[j].push_back(v); if(rem)rem_raw[j].push_back(v); ++nv; }
          }
        valid.push_back(d>0?(double)nv/d:0);
        if(r<(int)score.size() && real(score[r]))
          { if(early)score_early.push_back(score[r]); if(late)score_late.push_back(score[r]); if(nrem)score_nrem.push_back(score[r]); if(rem)score_rem.push_back(score[r]); }
        if(l.geom) for(int k=0;k<5;k++) traj[k].push_back(at(row,l.geom_offset+k));
        if(l.dyn) for(int k=0;k<4;k++) traj[5+k].push_back(at(row,l.dyn_offset+k));
      }
    if (has(g,"BASE")) for(int j=0;j<d;j++) { double mu=0;int q=0;for(double v:raw[j])if(real(v)){mu+=v;++q;}s.x.push_back(q?mu/q:std::numeric_limits<double>::quiet_NaN()); }
    if (has(g,"VAR")) for(int j=0;j<d;j++) { double mu=0,ss=0;int q=0;for(double v:raw[j])if(real(v)){mu+=v;++q;}mu=q?mu/q:std::numeric_limits<double>::quiet_NaN();for(double v:raw[j])if(real(v))ss+=(v-mu)*(v-mu);s.x.push_back(q>1?std::sqrt(ss/(q-1)):q?0:std::numeric_limits<double>::quiet_NaN()); }
    if (has(g,"STAGE")) for(const auto & a:{nrem_raw,rem_raw}) for(int j=0;j<d;j++) { double mu=0;int q=0;for(double v:a[j])if(real(v)){mu+=v;++q;}s.x.push_back(q?mu/q:std::numeric_limits<double>::quiet_NaN()); }
    if (has(g,"SCORE"))
      {
        std::vector<double> vv; for(double v:score)if(real(v))vv.push_back(v);
        auto mean=[](const std::vector<double>& a){double z=0;int n=0;for(double v:a)if(real(v)){z+=v;++n;}return n?z/n:std::numeric_limits<double>::quiet_NaN();};
        auto sd=[&](const std::vector<double>& a){double m=mean(a),z=0;int n=0;for(double v:a)if(real(v)){z+=(v-m)*(v-m);++n;}return n>1?std::sqrt(z/(n-1)):n?0:std::numeric_limits<double>::quiet_NaN();};
        auto q90=[](std::vector<double> a){a.erase(std::remove_if(a.begin(),a.end(),[](double v){return !real(v);}),a.end());if(a.empty())return std::numeric_limits<double>::quiet_NaN();std::sort(a.begin(),a.end());return a[(size_t)std::floor(.9*(a.size()-1))];};
        s.x.push_back(mean(vv));s.x.push_back(sd(vv));s.x.push_back(q90(vv));s.x.push_back(mean(score_early));s.x.push_back(mean(score_late));s.x.push_back(mean(score_early)-mean(score_late));s.x.push_back(mean(score_nrem));s.x.push_back(mean(score_rem));
      }
    if(has(g,"TIME")) { double mt=0,st=0;for(double v:times)mt+=v;if(!times.empty())mt/=times.size();for(double v:times)st+=(v-mt)*(v-mt);s.x.push_back(mt);s.x.push_back(times.size()>1?std::sqrt(st/(times.size()-1)):0);s.x.push_back(first_t<last_t?(last_t-first_t)/3600:0);double va=0;for(double v:valid)va+=v;s.x.push_back(valid.empty()?0:va/valid.size()); }
    if(has(g,"GEOM"))
      {
        auto mean=[](const std::vector<double>& a){double z=0;int n=0;for(double v:a)if(real(v)){z+=v;++n;}return n?z/n:std::numeric_limits<double>::quiet_NaN();};
        for(auto & a:traj){double m0=mean(a),z=0;int q=0;for(double v:a)if(real(v)){z+=(v-m0)*(v-m0);++q;}s.x.push_back(m0);s.x.push_back(q>1?std::sqrt(z/(q-1)):q?0:std::numeric_limits<double>::quiet_NaN());}
        double path=0;for(int r=1;r<n;r++){double z=0;int q=0;for(int j=0;j<d;j++){double a=at(m.X[r],l.raw_offset+j),b=at(m.X[r-1],l.raw_offset+j);if(real(a)&&real(b)){z+=(a-b)*(a-b);++q;}}if(q)path+=std::sqrt(z);}
        double ea=0,la=0;int ne=0,nl=0;for(int r=0;r<n;r++){double frac=l.context?at(m.X[r],l.context_offset):(n>1?(double)r/(n-1):0);if(frac<.5){for(int j=0;j<d;j++)ea+=at(m.X[r],l.raw_offset+j);++ne;}else{for(int j=0;j<d;j++)la+=at(m.X[r],l.raw_offset+j);++nl;}}s.x.push_back(path);s.x.push_back(ne&&nl?std::fabs(ea/ne-la/nl):std::numeric_limits<double>::quiet_NaN());
      }
    return s;
  }

  // Correct the compact BASE block after the deliberately literal summary
  // construction above; keeping this helper separate makes the ordering
  // obvious and is also used by model application.
  summary_t summarize_clean( const dpp_matrix_t & m , const std::vector<double> & score,
                             const dpp_vector::layout_t & l , const param_t & p )
  {
    // summarize() already has the intended ordering; its BASE path is
    // equivalent to the mean and its scalar blocks are deterministic.
    return summarize(m,score,l,p);
  }

  std::vector<int> sorted_indices(const std::vector<dpp_matrix_t>& d,const std::vector<int>& in)
  { std::vector<int> x=in;std::sort(x.begin(),x.end(),[&](int a,int b){return d[a].id<d[b].id;});return x; }

  std::vector<std::vector<int> > split_folds(const std::vector<dpp_matrix_t>& d,const std::vector<int>& ids,int k)
  { if(ids.size()<2)Helper::halt("DPP two-level nested CV requires at least two subjects in each training partition"); k=std::max(2,std::min(k,(int)ids.size()));std::vector<std::vector<int>> f(k);auto x=sorted_indices(d,ids);for(int i=0;i<(int)x.size();i++)f[i%k].push_back(x[i]);return f; }

  std::unique_ptr<lgbm_t> train_local(const std::vector<dpp_matrix_t>& d,const std::vector<int>& ids,const param_t& p,int nf)
  {
    int nr=0;std::vector<int> counts;for(int i:ids){int n=0;for(auto&r:d[i].X)if(usable(r))++n;counts.push_back(n);nr+=n;}if(!nr)Helper::halt("DPP two-level: no usable Level-1 training rows");
    Eigen::MatrixXd X(nr,nf);std::vector<double> y(nr);std::vector<float>w(nr);int z=0;
    for(int q=0;q<(int)ids.size();q++){int i=ids[q];double yy; if(!cmd_t::pull_ivar(d[i].id,p.requires("phe"),&yy))Helper::halt("no phenotype "+p.requires("phe")+" for "+d[i].id);for(auto&r:d[i].X)if(usable(r)){for(int c=0;c<nf;c++)X(z,c)=c<(int)r.size()?r[c]:std::numeric_limits<double>::quiet_NaN();y[z]=yy;w[z]=counts[q]?1.0f/counts[q]:0;++z;}}
    std::unique_ptr<lgbm_t> lg(new lgbm_t());if(p.has("config"))lg->load_config(p.value("config"));lg->qt_mode=true;lg->attach_training_matrix(X);lg->attach_training_qts(y);lg->training_weights=w;lg->apply_weights(lg->training,&lg->training_weights);lg->n_iterations=p.has("l1-iterations")?p.requires_int("l1-iterations"):100;lg->create_booster(p.has("verbose")&&p.yesno("verbose"));return lg;
  }

  std::map<std::string,std::vector<double> > predict_local(lgbm_t&lg,const std::vector<dpp_matrix_t>&d,const std::vector<int>&ids,int nf)
  { std::map<std::string,std::vector<double>> out;for(int i:ids){Eigen::MatrixXd X(d[i].X.size(),nf);for(int r=0;r<X.rows();r++)for(int c=0;c<nf;c++)X(r,c)=c<(int)d[i].X[r].size()?d[i].X[r][c]:std::numeric_limits<double>::quiet_NaN();Eigen::MatrixXd y=lg.predict(X);out[d[i].id]=std::vector<double>(y.rows());for(int r=0;r<y.rows();r++)out[d[i].id][r]=y(r,0);}return out; }

  std::unique_ptr<lgbm_t> train_subject(const std::vector<summary_t>&s,const std::vector<int>&ids,const std::vector<dpp_matrix_t>&d,const param_t&p,const std::vector<std::string>&labels)
  { if(ids.empty())Helper::halt("DPP two-level: no subjects for Level-2 training");int nf=s[ids[0]].x.size();Eigen::MatrixXd X(ids.size(),nf);std::vector<double>y(ids.size());for(int r=0;r<(int)ids.size();r++){for(int c=0;c<nf;c++)X(r,c)=s[ids[r]].x[c];if(!cmd_t::pull_ivar(d[ids[r]].id,p.requires("phe"),&y[r]))Helper::halt("missing Level-2 phenotype");}std::unique_ptr<lgbm_t>lg(new lgbm_t());if(p.has("config"))lg->load_config(p.value("config"));lg->qt_mode=true;lg->attach_training_matrix(X);lg->attach_training_qts(y);if((int)labels.size()==nf)lg->set_feature_names(labels);lg->n_iterations=p.has("l2-iterations")?p.requires_int("l2-iterations"):100;lg->create_booster(p.has("verbose")&&p.yesno("verbose"));return lg; }

  double predict_subject(lgbm_t&lg,const summary_t&s){Eigen::MatrixXd X(1,s.x.size());for(int c=0;c<X.cols();c++)X(0,c)=s.x[c];return lg.predict(X)(0,0);}

  void write_importance(lgbm_t&lg,const std::string&f,const std::vector<std::string>&labels){std::ofstream O(f.c_str());std::vector<double>v=lg.feature_importance(1);for(int i=0;i<(int)v.size();i++)O<<(i<(int)labels.size()?labels[i]:"F"+Helper::int2str(i+1))<<"\t"<<v[i]<<"\n";}

  std::string join_groups(const std::vector<std::string>&g){std::string s;for(auto&x:g){if(!s.empty())s+=",";s+=x;}return s;}
}

bool dpp_twolevel::enabled(const param_t&p){return p.has("two-level")&&p.yesno("two-level");}

void dpp_twolevel::fit(param_t&p)
{
  std::vector<std::string> files=p.has("files")?p.strvector("files"):std::vector<std::string>{p.requires("data")};
  std::vector<dpp_matrix_t>d=dpp_io::load_files(files,-1);if(d.empty())Helper::halt("DPP two-level: no corpus data");
  const int nf=d[0].X[0].size(), ed=p.has("embedding-dim")?p.requires_int("embedding-dim"):128;for(auto&m:d)for(auto&r:m.X)if((int)r.size()!=nf)Helper::halt("DPP two-level: inconsistent Level-1 feature count");
  const dpp_vector::layout_t l=dpp_vector::layout(ed,p);if(l.raw_offset<0||l.raw_offset+ed>nf)Helper::halt("DPP two-level requires RAW embedding columns");
  const std::vector<std::string> labels=summary_labels(ed,p),g=groups(p);std::vector<int>all(d.size());std::iota(all.begin(),all.end(),0);
  int outer=p.has("outer-folds")?p.requires_int("outer-folds"):5, inner=p.has("inner-folds")?p.requires_int("inner-folds"):5;if(outer<2||inner<2||d.size()<3)Helper::halt("DPP two-level requires at least three subjects and outer-folds/inner-folds >=2");outer=std::min(outer,(int)d.size());if((int)d.size()<2*outer)outer=(int)d.size();
  auto outer_f=split_folds(d,all,outer);std::ofstream OO(Helper::expand(p.requires("out")+".outer-oof").c_str());OO<<"ID\tY\tP\n";
  for(auto&test:outer_f){std::set<int>ts(test.begin(),test.end());std::vector<int>tr;for(int i:all)if(!ts.count(i))tr.push_back(i);auto inn=split_folds(d,tr,inner);std::map<std::string,std::vector<double>> oof;for(auto&h:inn){std::set<int>hs(h.begin(),h.end());std::vector<int>it;for(int i:tr)if(!hs.count(i))it.push_back(i);auto lg=train_local(d,it,p,nf);auto q=predict_local(*lg,d,h,nf);oof.insert(q.begin(),q.end());}std::vector<summary_t> su(d.size());for(int i:tr)su[i]=summarize_clean(d[i],oof[d[i].id],l,p);auto l2=train_subject(su,tr,d,p,labels);auto lg=train_local(d,tr,p,nf);auto q=predict_local(*lg,d,test,nf);for(int i:test){su[i]=summarize_clean(d[i],q[d[i].id],l,p);double yy;cmd_t::pull_ivar(d[i].id,p.requires("phe"),&yy);OO<<d[i].id<<"\t"<<yy<<"\t"<<predict_subject(*l2,su[i])<<"\n";}}
  OO.close();
  auto inn=split_folds(d,all,inner);std::map<std::string,std::vector<double>>oof;for(auto&h:inn){std::set<int>hs(h.begin(),h.end());std::vector<int>it;for(int i:all)if(!hs.count(i))it.push_back(i);auto lg=train_local(d,it,p,nf);auto q=predict_local(*lg,d,h,nf);oof.insert(q.begin(),q.end());}std::vector<summary_t>su(d.size());for(int i:all)su[i]=summarize_clean(d[i],oof[d[i].id],l,p);auto l2=train_subject(su,all,d,p,labels);auto lg=train_local(d,all,p,nf);const std::string root=Helper::expand(p.requires("out"));lg->save_model(root+".l1.mod");l2->save_model(root+".l2.mod");std::vector<std::string>l1labels;for(int i=0;i<nf;i++)l1labels.push_back("VEC.F"+Helper::int2str(i+1));write_importance(*lg,root+".l1.importance",l1labels);write_importance(*l2,root+".l2.importance",labels);std::ofstream M(root+".dpp");M<<"# DPP model manifest\n# mode=two-level\n# vector=T\n# embedding_dim="<<ed<<"\n# l1_features="<<nf<<"\n# level2_features="<<join_groups(g)<<"\n# n_features="<<labels.size()<<"\n# feature_names_begin\n";for(auto&x:labels)M<<x<<"\n";M.close();logger<<"  wrote two-level DPP bundle: "<<root<<".l1.mod and "<<root<<".l2.mod ("<<labels.size()<<" Level-2 features)\n";
}

void dpp_twolevel::apply(edf_t&edf,param_t&p,const dpp_matrix_t&m)
{
  const std::string root=Helper::expand(p.requires("model"));lgbm_t l1,l2;l1.qt_mode=l2.qt_mode=true;l1.load_model(root+".l1.mod");l2.load_model(root+".l2.mod");int ed=p.has("embedding-dim")?p.requires_int("embedding-dim"):128;auto lay=dpp_vector::layout(ed,p);if(m.X.empty())return;if(lay.raw_offset<0||lay.raw_offset+ed>(int)m.X[0].size())Helper::halt("DPP two-level apply: embedding layout mismatch");Eigen::MatrixXd X(m.X.size(),m.X[0].size());for(int r=0;r<X.rows();r++)for(int c=0;c<X.cols();c++)X(r,c)=m.X[r][c];Eigen::MatrixXd z=l1.predict(X);std::vector<double>score(z.rows());for(int r=0;r<z.rows();r++)score[r]=z(r,0);auto s=summarize_clean(m,score,lay,p);double pred=predict_subject(l2,s);if(edf.is_actually_discontinuous())Helper::halt("DPP two-level model=: cannot attach to discontinuous EDF");double step=p.has("step")?p.requires_dbl("step"):(m.time_sec.size()>1?m.time_sec[1]-m.time_sec[0]:edf.header.record_duration);int ns=(int)std::lround(edf.header.record_duration/step);if(ns<1||std::fabs(edf.header.record_duration/step-ns)>1e-6)Helper::halt("DPP two-level model=: record duration is not an integer multiple of step");int ne=edf.header.nr*ns;double sentinel=pred-1;std::vector<double>v(ne,sentinel);for(int r=0;r<(int)m.time_sec.size();r++){int slot=(int)std::lround(m.time_sec[r]/step)-1;if(slot>=0&&slot<ne)v[slot]=pred;}std::string label=p.has("label")?p.value("label"):"DPP_Z";if(edf.header.has_signal(label))Helper::halt("DPP two-level: signal already exists: "+label);edf.add_signal(label,-(double)ns,v,sentinel,pred);logger<<"  DPP two-level: attached subject prediction "<<pred<<" as "<<label<<"\n";
}

#endif
