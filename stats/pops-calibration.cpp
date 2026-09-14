// POPS posterior calibration and optional global weight/temperature scaling.
// This module deliberately has no LightGBM dependency.

#include "param.h"
#include "db/db.h"
#include "helper/helper.h"
#include "helper/logger.h"
#include "helper/zfile.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <fstream>
#include <limits>
#include <set>
#include <vector>

extern writer_t writer;
extern logger_t logger;

namespace {

const int NS = 5;
const double EPS = 1e-12;
const double RAW_LOGLOSS_EPS = 1e-15; // preserves pre-temperature behavior
const double POSTERIOR_SUM_WARN_TOL = 0.01;
const std::array<std::string,NS> labels = {{ "W", "R", "N1", "N2", "N3" }};

struct row_t {
  std::string id, epoch, observed;
  int y;
  std::array<double,NS> raw;
};
struct qc_t {
  int read = 0, used = 0, bad_stage = 0, bad_posterior = 0, bad_sum = 0;
  int sum_rows = 0;
  int subjects = 0;
};
struct bin_t { int index; double lo, hi; int n; double pred, obs; };
struct stage_t {
  int n = 0;
  double prevalence = 0, mean_p = 0, brier = 0, ece = 0;
  std::vector<bin_t> bins;
};
struct metrics_t {
  std::array<stage_t,NS> stages;
  double logloss = 0, multiclass_brier = 0, ece = 0;
};

int stage_index( const std::string & s ) {
  for (int i = 0; i < NS; ++i) if ( s == labels[i] ) return i;
  return -1;
}
// destrat's tabular writer uses tabs.  Preserve those field positions exactly:
// splitting tabular rows again on spaces can shift columns when an identifier
// or another non-posterior field contains a space.
std::vector<std::string> split_destrat_row( const std::string & line ) {
  return line.find('\t') == std::string::npos
    ? Helper::parse(line," ")
    : Helper::parse(line,'\t');
}
bool valid_positive( double x ) { return std::isfinite(x) && x > 0; }
int argmax( const std::array<double,NS> & p ) {
  int best = 0;
  for (int k = 1; k < NS; ++k) if ( p[k] > p[best] ) best = k;
  return best;
}
double logsumexp( const std::array<double,NS> & z ) {
  double m = z[0];
  for (int k = 1; k < NS; ++k) m = std::max(m, z[k]);
  double total = 0;
  for (int k = 0; k < NS; ++k) total += std::exp(z[k] - m);
  return m + std::log(total);
}

// Required order: remove the known loss weight in log space, then apply T.
std::array<double,NS> transform( const std::array<double,NS> & raw,
                                  const std::array<double,NS> * weights,
                                  double temperature ) {
  std::array<double,NS> z, out;
  for (int k = 0; k < NS; ++k)
    z[k] = ( std::log(std::max(raw[k], EPS))
             - ( weights ? std::log((*weights)[k]) : 0.0 ) ) / temperature;
  const double lse = logsumexp(z);
  for (int k = 0; k < NS; ++k) out[k] = std::exp(z[k] - lse);
  return out;
}
double nll( const std::vector<row_t> & rows,
            const std::array<double,NS> * weights, double temperature ) {
  double total = 0;
  for (const row_t & row : rows) {
    std::array<double,NS> z;
    for (int k = 0; k < NS; ++k)
      z[k] = (std::log(std::max(row.raw[k], EPS))
              - (weights ? std::log((*weights)[k]) : 0.0)) / temperature;
    total += logsumexp(z) - z[row.y];
  }
  return total / rows.size();
}
double fit_temperature( const std::vector<row_t> & rows,
                        const std::array<double,NS> * weights, bool * boundary ) {
  const double lo0 = .05, hi0 = 20.0, ratio = .6180339887498948482;
  double lo = lo0, hi = hi0;
  double x1 = hi - ratio * (hi - lo), x2 = lo + ratio * (hi - lo);
  double f1 = nll(rows,weights,x1), f2 = nll(rows,weights,x2);
  for (int it = 0; it < 180 && hi-lo > 1e-6 * std::max(1.0,(hi+lo)*.5); ++it)
    if (f1 > f2) {
      lo=x1; x1=x2; f1=f2; x2=lo+ratio*(hi-lo); f2=nll(rows,weights,x2);
    } else {
      hi=x2; x2=x1; f2=f1; x1=hi-ratio*(hi-lo); f1=nll(rows,weights,x1);
    }
  const double out = (lo+hi)*.5;
  if(boundary) *boundary = out-lo0 < .001 || hi0-out < .001;
  return out;
}

std::vector<double> edges_for( const std::vector<std::array<double,NS> > & p,
                               int stage, int bins, bool quantile ) {
  std::vector<double> edge(bins+1);
  if (!quantile)
    for (int b=0;b<=bins;++b) edge[b]=b/(double)bins;
  else {
    std::vector<double> x;
    for(const auto & row:p) x.push_back(row[stage]);
    std::sort(x.begin(),x.end());
    edge[0]=0; edge[bins]=1;
    for(int b=1;b<bins;++b) edge[b]=x[(int)std::ceil(b*x.size()/(double)bins)-1];
  }
  return edge;
}
metrics_t metrics( const std::vector<row_t> & rows,
                   const std::vector<std::array<double,NS> > & p,
                   int bins, bool quantile ) {
  metrics_t m;
  for(int s=0;s<NS;++s) {
    stage_t & st=m.stages[s]; st.n=rows.size();
    const std::vector<double> edge=edges_for(p,s,bins,quantile);
    struct a_t { int n=0; double p=0,y=0; }; std::vector<a_t>a(bins);
    for(int i=0;i<(int)rows.size();++i) {
      const double ps=p[i][s], y=rows[i].y==s?1.0:0.0;
      const int b=quantile ? std::min(bins-1,(int)(std::upper_bound(edge.begin(),edge.end(),ps)-edge.begin())-1)
                           : std::min(bins-1,(int)std::floor(ps*bins));
      ++a[b].n; a[b].p+=ps; a[b].y+=y;
      st.mean_p+=ps; st.prevalence+=y; st.brier+=(ps-y)*(ps-y);
    }
    st.mean_p/=st.n; st.prevalence/=st.n; st.brier/=st.n;
    for(int b=0;b<bins;++b) if(a[b].n) {
      bin_t out={b+1,edge[b],edge[b+1],a[b].n,a[b].p/a[b].n,a[b].y/a[b].n};
      st.ece+=out.n/(double)st.n*std::fabs(out.obs-out.pred); st.bins.push_back(out);
    }
    m.ece+=st.ece; m.multiclass_brier+=st.brier;
  }
  m.ece/=NS;
  for(int i=0;i<(int)rows.size();++i) m.logloss-=std::log(std::max(p[i][rows[i].y],RAW_LOGLOSS_EPS));
  m.logloss/=rows.size();
  return m;
}

std::vector<row_t> read_rows( const std::string & filename, qc_t * qc ) {
  const std::string expanded=Helper::expand(filename);
  if(!Helper::fileExists(expanded)) Helper::halt("Luna calibration: posteriors file not found: "+filename);
  std::ifstream in=LunaIO::open_ifstream(expanded);
  if(!in.good()) Helper::halt("Luna calibration: could not open: "+filename);
  std::string header; Helper::safe_getline(in,header);
  if(header.empty()) Helper::halt("Luna calibration: empty header in "+filename);
  const std::vector<std::string> col=split_destrat_row(header);
  int id=-1,e=-1,prior=-1; std::array<int,NS> pp; pp.fill(-1);
  for(int c=0;c<(int)col.size();++c) {
    if(col[c]=="ID") id=c; else if(col[c]=="E") e=c;
    else if(col[c]=="PRIOR"||col[c]=="PRIOR30") prior=c;
    else if(col[c]=="PP_W")pp[0]=c; else if(col[c]=="PP_R")pp[1]=c;
    else if(col[c]=="PP_N1")pp[2]=c; else if(col[c]=="PP_N2")pp[3]=c;
    else if(col[c]=="PP_N3")pp[4]=c;
  }
  if(id<0||e<0||prior<0) Helper::halt("Luna calibration: expected ID, E, and PRIOR/PRIOR30 columns");
  for(int k=0;k<NS;++k) if(pp[k]<0) Helper::halt("Luna calibration: missing PP_"+labels[k]+" column");
  std::vector<row_t> rows; std::set<std::string> subjects; std::string line;
  while(Helper::safe_getline(in,line)) {
    if(line.empty()) continue; ++qc->read;
    const std::vector<std::string> tok=split_destrat_row(line); const int n=tok.size();
    row_t row; row.id=id<n?tok[id]:""; row.epoch=e<n?tok[e]:""; row.observed=prior<n?tok[prior]:""; row.y=stage_index(row.observed);
    if(!row.id.empty()) subjects.insert(row.id);
    bool okay=true;
    for(int k=0;k<NS;++k) {
      double x=std::numeric_limits<double>::quiet_NaN();
      // str2dbl() may leave a finite value on a failed conversion (e.g. NA),
      // so its return value is required for PP-column validity.
      const bool parsed = pp[k]<n && Helper::str2dbl(tok[pp[k]],&x);
      row.raw[k]=x;
      if(!parsed || !std::isfinite(x)||x<0||x>1) okay=false;
    }
    if(!okay) { ++qc->bad_posterior; continue; }
    // Sum the five raw, directly parsed PP columns before any correction,
    // temperature scaling, binning, or other derived calculation.
    const double posterior_sum =
      row.raw[2] + row.raw[3] + row.raw[4] + row.raw[1] + row.raw[0];
    const double sum_deviation = std::fabs(posterior_sum - 1.0);
    ++qc->sum_rows;
    if(sum_deviation > POSTERIOR_SUM_WARN_TOL) ++qc->bad_sum;
    if(row.y<0) { ++qc->bad_stage; continue; }
    rows.push_back(row);
  }
  qc->used=rows.size(); qc->subjects=subjects.size();
  if(rows.empty()) Helper::halt("Luna calibration: no usable rows remain after QC");
  return rows;
}

bool parse_weights( const std::string & text, std::array<double,NS> * weights ) {
  const std::vector<std::string> tokens=Helper::parse(text,",");
  if(tokens.size()!=NS) return false;
  for(int k=0;k<NS;++k)
    if(!Helper::str2dbl(tokens[k],&(*weights)[k]) || !valid_positive((*weights)[k])) return false;
  return true;
}
void save_temperature( const std::string & filename, double temperature,
                       const std::array<double,NS> * weights ) {
  std::ofstream out(Helper::expand(filename).c_str());
  if(!out.good()) Helper::halt("Luna calibration: could not write temperature file: "+filename);
  out.precision(17); out<<"temperature "<<temperature<<"\nweights ";
  if(!weights) out<<"none\n";
  else { for(int k=0;k<NS;++k) { if(k) out<<" "; out<<(*weights)[k]; } out<<"\n"; }
}
double load_temperature( const std::string & filename, const std::array<double,NS> * supplied ) {
  std::ifstream in=LunaIO::open_ifstream(Helper::expand(filename));
  if(!in.good()) Helper::halt("Luna calibration: could not read temperature file: "+filename);
  bool got=false, metadata=false, none=false; double temperature=0; std::array<double,NS> saved;
  std::string line;
  while(Helper::safe_getline(in,line)) {
    const std::vector<std::string>x=Helper::parse(line," \t");
    if(x.empty()||x[0][0]=='#') continue;
    if(x[0]=="temperature") { if(x.size()!=2||!Helper::str2dbl(x[1],&temperature)) Helper::halt("Luna calibration: invalid temperature file"); got=true; }
    else if(x[0]=="weights") {
      metadata=true;
      if(x.size()==2&&x[1]=="none") none=true;
      else { if(x.size()!=NS+1) Helper::halt("Luna calibration: invalid weight metadata");
        for(int k=0;k<NS;++k) if(!Helper::str2dbl(x[k+1],&saved[k])||!valid_positive(saved[k])) Helper::halt("Luna calibration: invalid weight metadata"); }
    } else if(!got&&x.size()==1&&Helper::str2dbl(x[0],&temperature)) got=true;
    else Helper::halt("Luna calibration: unrecognized temperature file content");
  }
  if(!got||!valid_positive(temperature)) Helper::halt("Luna calibration: invalid saved temperature");
  if(metadata) {
    if(none!=(supplied==NULL)) Helper::halt("Luna calibration: saved temperature weight convention does not match training-weights=");
    if(!none) for(int k=0;k<NS;++k)
      if(std::fabs(saved[k]-(*supplied)[k])>1e-12*std::max(1.0,std::fabs(saved[k])))
        Helper::halt("Luna calibration: saved temperature weights do not match training-weights=");
  }
  return temperature;
}

void emit_metrics( const std::string & scale, const metrics_t & m, bool use_scale ) {
  if(use_scale) writer.level(scale,"SCALE");
  for(int s=0;s<NS;++s) for(const bin_t &b:m.stages[s].bins) {
    writer.level(labels[s],"STAGE"); writer.level(b.index,"BIN");
    writer.value("P_LO",b.lo); writer.value("P_HI",b.hi); writer.value("N",b.n);
    writer.value("PRED",b.pred); writer.value("OBS",b.obs);
  }
  // Keep the private command stratum set by main.cpp.  A bare unlevel()
  // clears it, leaving subsequent summary/QC rows with command = NA.
  writer.unlevel("BIN");
  writer.unlevel("STAGE");
  if(use_scale) writer.level(scale,"SCALE");
  for(int s=0;s<NS;++s) {
    const stage_t&st=m.stages[s]; writer.level(labels[s],"STAGE");
    writer.value("PREV",st.prevalence); writer.value("MEAN_P",st.mean_p);
    writer.value("BRIER",st.brier); writer.value("ECE",st.ece);
  }
  writer.unlevel("STAGE"); writer.value("LOGLOSS",m.logloss);
  if(use_scale) writer.unlevel("SCALE");
}
void write_posteriors( const std::string & filename, const std::vector<row_t> & rows,
                       const std::vector<std::pair<std::string,std::vector<std::array<double,NS> > > > & scales ) {
  std::ofstream out(Helper::expand(filename).c_str());
  if(!out.good()) Helper::halt("Luna calibration: could not write posterior-out file: "+filename);
  out<<"ID\tE\tPRIOR\tPP_W\tPP_R\tPP_N1\tPP_N2\tPP_N3";
  for(const auto&scale:scales) for(int k=0;k<NS;++k) out<<"\tPP_"<<labels[k]<<"_"<<scale.first;
  out<<"\n"; out.precision(17);
  for(int i=0;i<(int)rows.size();++i) {
    const row_t&r=rows[i]; out<<r.id<<"\t"<<r.epoch<<"\t"<<r.observed;
    for(int k=0;k<NS;++k) out<<"\t"<<r.raw[k];
    for(const auto&scale:scales) for(int k=0;k<NS;++k) out<<"\t"<<scale.second[i][k];
    out<<"\n";
  }
}

} // namespace

namespace pops_calibration {

void run( param_t & param ) {
  const bool fit=param.has("fit-temperature"), has_value=param.has("temperature"), has_file=param.has("temperature-file");
  if((fit?1:0)+(has_value?1:0)+(has_file?1:0)>1)
    Helper::halt("Luna calibration: fit-temperature, temperature=, and temperature-file= are mutually exclusive");
  const bool has_weights=param.has("training-weights"); std::array<double,NS> weights;
  if(has_weights&&!parse_weights(param.value("training-weights"),&weights))
    Helper::halt("Luna calibration: training-weights= requires five finite positive W,R,N1,N2,N3 values");
  bool use_temperature=fit||has_value||has_file; double temperature=1;
  if(has_value) { temperature=param.requires_dbl("temperature"); if(!valid_positive(temperature)) Helper::halt("Luna calibration: temperature= must be finite and > 0"); }
  if(has_file) temperature=load_temperature(param.requires("temperature-file"),has_weights?&weights:NULL);
  if(param.has("temperature-out")&&!use_temperature) Helper::halt("Luna calibration: temperature-out= requires a fitted or supplied temperature");
  const int bins=param.has("bins")?param.requires_int("bins"):10;
  if(bins<1) Helper::halt("Luna calibration: bins must be at least 1");
  std::string binning=param.has("binning")?param.value("binning"):"fixed";
  std::transform(binning.begin(),binning.end(),binning.begin(),[](unsigned char c){return std::tolower(c);});
  if(binning!="fixed"&&binning!="quantile") Helper::halt("Luna calibration: binning must be fixed or quantile");

  logger << "  CALIBRATION options:\n"
         << "    file=<destrat posterior table>  (required)\n"
         << "    bins=N                          (default 10)\n"
         << "    binning=fixed|quantile          (default fixed)\n"
         << "    fit-temperature\n"
         << "    temperature=T\n"
         << "    temperature-file=FILE\n"
         << "    training-weights=W,R,N1,N2,N3\n"
         << "    temperature-out=FILE\n"
         << "    posterior-out=FILE\n"
         << "  CALIBRATION settings:\n"
         << "    file=" << param.requires("file") << "\n"
         << "    bins=" << bins << "\n"
         << "    binning=" << binning << "\n"
         << "    temperature=" << ( fit ? "fit" : has_value ? param.value("temperature")
                               : has_file ? "file:" + param.value("temperature-file") : "none" ) << "\n"
         << "    training-weights=" << ( has_weights ? param.value("training-weights") : "none" ) << "\n"
         << "    temperature-out=" << ( param.has("temperature-out") ? param.value("temperature-out") : "none" ) << "\n"
         << "    posterior-out=" << ( param.has("posterior-out") ? param.value("posterior-out") : "none" )
         << "\n";

  qc_t qc; const std::vector<row_t> rows=read_rows(param.requires("file"),&qc);
  if(fit) { bool boundary=false; temperature=fit_temperature(rows,has_weights?&weights:NULL,&boundary);
    if(boundary) logger<<"  Luna calibration: fitted temperature reached search boundary; consider wider bounds\n"; }
  std::vector<std::array<double,NS> > raw; for(const row_t&r:rows) raw.push_back(r.raw);
  const metrics_t raw_metrics=metrics(rows,raw,bins,binning=="quantile");
  const bool transformed=has_weights||use_temperature;
  std::vector<std::pair<std::string,std::vector<std::array<double,NS> > > > scales;
  if(transformed) scales.push_back({"RAW",raw});
  if(has_weights) { std::vector<std::array<double,NS> > p; for(const row_t&r:rows)p.push_back(transform(r.raw,&weights,1)); scales.push_back({"WEIGHT",p}); }
  if(use_temperature) { std::vector<std::array<double,NS> > p; for(const row_t&r:rows)p.push_back(transform(r.raw,has_weights?&weights:NULL,temperature)); scales.push_back({has_weights?"WEIGHT_TEMP":"TEMP",p}); }
  int changed_weight=0, changed_temp=0, nontrivial_temp=0;
  if(has_weights) for(int i=0;i<(int)rows.size();++i) if(argmax(raw[i])!=argmax(scales[1].second[i])) ++changed_weight;
  if(use_temperature) {
    const std::vector<std::array<double,NS> >&base=has_weights?scales[1].second:raw;
    const std::vector<std::array<double,NS> >&final=scales.back().second;
    for(int i=0;i<(int)rows.size();++i) if(argmax(base[i])!=argmax(final[i])) {
      ++changed_temp; std::array<double,NS>x=base[i]; std::sort(x.begin(),x.end(),std::greater<double>());
      if(x[0]-x[1]>1e-10) ++nontrivial_temp;
    }
  }
  logger<<"  CALIBRATION QC:\n"
        <<"    usable epochs="<<qc.used<<"\n"
        <<"    subjects="<<qc.subjects<<"\n"
        <<"    excluded overall="<<qc.read-qc.used<<"\n"
        <<"      invalid observed stage="<<qc.bad_stage<<"\n"
        <<"      missing/non-numeric posterior="<<qc.bad_posterior<<"\n";
  if(qc.bad_sum)
    logger<<"  Luna calibration warning: "<<qc.bad_sum
          <<" raw posterior rows have abs(sum(PP_*) - 1) > "
          <<POSTERIOR_SUM_WARN_TOL<<"; raw values were not renormalized\n";
  if(nontrivial_temp) logger<<"  Luna calibration warning: temperature scaling changed "<<nontrivial_temp<<" non-tied argmax calls\n";
  if(!transformed) emit_metrics("",raw_metrics,false);
  else {
    for(const auto&scale:scales) emit_metrics(scale.first,metrics(rows,scale.second,bins,binning=="quantile"),true);
    writer.value("ARGMAX_CHANGED_WEIGHT",changed_weight);
    writer.value("ARGMAX_CHANGED_TEMP",changed_temp);
    if(use_temperature) {
      const metrics_t final=metrics(rows,scales.back().second,bins,binning=="quantile");
      writer.value("TEMPERATURE",temperature);
      writer.value("LOGLOSS_RAW",raw_metrics.logloss); writer.value("LOGLOSS_TEMP",final.logloss);
      writer.value("BRIER_RAW",raw_metrics.multiclass_brier); writer.value("BRIER_TEMP",final.multiclass_brier);
      writer.value("ECE_RAW",raw_metrics.ece); writer.value("ECE_TEMP",final.ece);
    }
  }
  // Aggregate QC shares the unstratified row with LOGLOSS.
  writer.value("N",qc.used);
  writer.value("N_READ",qc.read); writer.value("N_USED",qc.used);
  writer.value("N_BAD_STAGE",qc.bad_stage); writer.value("N_BAD_POSTERIOR",qc.bad_posterior);
  writer.value("N_SUM_ROWS",qc.sum_rows);
  writer.value("N_SUM_NOT_1",qc.bad_sum); // compatibility alias
  writer.value("N_SUM_GT_0_01",qc.bad_sum);
  writer.value("N_SUBJECTS",qc.subjects);
  writer.value("N_EXCLUDED",qc.read-qc.used);
  if(param.has("temperature-out")) save_temperature(param.value("temperature-out"),temperature,has_weights?&weights:NULL);
  if(param.has("posterior-out")) {
    if(!transformed) Helper::halt("Luna calibration: posterior-out= requires training-weights= or a temperature option");
    write_posteriors(param.value("posterior-out"),rows,scales);
  }
}

bool self_test( std::string * message ) {
  const std::array<double,NS> raw={{.8,.1,.05,.03,.02}}, one={{1,1,1,1,1}}, weights={{1,1.5,1.5,1,1.5}};
  const auto t1=transform(raw,NULL,1), soft=transform(raw,NULL,2), sharp=transform(raw,NULL,.5), equal=transform(raw,&one,1);
  double sum=0; for(double x:soft) sum+=x;
  bool okay=std::fabs(sum-1)<1e-12&&argmax(raw)==argmax(t1)&&argmax(raw)==argmax(soft)&&argmax(raw)==argmax(sharp)&&soft[0]<raw[0]&&sharp[0]>raw[0];
  for(int k=0;k<NS;++k) okay=okay&&std::fabs(t1[k]-raw[k])<1e-12&&std::fabs(equal[k]-raw[k])<1e-12;
  std::array<double,NS> weighted; double d=0; for(int k=0;k<NS;++k)d+=raw[k]*weights[k];
  for(int k=0;k<NS;++k) weighted[k]=raw[k]*weights[k]/d;
  const auto recovered=transform(weighted,&weights,1);
  for(int k=0;k<NS;++k) okay=okay&&std::fabs(recovered[k]-raw[k])<1e-12;
  // Deliberately mismatched confidence: the one-dimensional fit must soften
  // overconfident rows and sharpen underconfident rows, never worsening NLL.
  std::vector<row_t> over, under;
  for(int i=0;i<100;++i) {
    row_t a; a.raw={{.9,.1,0,0,0}}; a.y=i<60?0:1; over.push_back(a);
    row_t b; b.raw={{.6,.4,0,0,0}}; b.y=i<90?0:1; under.push_back(b);
  }
  bool at_boundary=false;
  const double over_t=fit_temperature(over,NULL,&at_boundary);
  const double under_t=fit_temperature(under,NULL,&at_boundary);
  okay=okay&&over_t>1&&under_t<1
       &&nll(over,NULL,over_t)<=nll(over,NULL,1)+1e-10
       &&nll(under,NULL,under_t)<=nll(under,NULL,1)+1e-10;
  const double nan=std::numeric_limits<double>::quiet_NaN();
  okay=okay&&!valid_positive(0)&&!valid_positive(-1)&&!valid_positive(nan)
       &&!valid_positive(std::numeric_limits<double>::infinity());
  if(message) *message=okay?"T=1, normalization, argmax, sharpness, and weight-correction checks passed":"temperature/weight transform check failed";
  return okay;
}

} // namespace pops_calibration
