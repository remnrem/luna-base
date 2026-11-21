// Part of Luna; GPL-3.0-or-later (see LICENSE.txt).
#include "dsp/project.h"
#include "helper/helper.h"
#include "helper/logger.h"
#include "edf/edf.h"
#include "param.h"
#include "db/db.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
extern logger_t logger;
extern writer_t writer;
namespace project {
namespace {
void require(bool ok, const std::string & message) {
  if (!ok) throw std::runtime_error(message);
}
std::string read(const std::string & file) {
  std::ifstream in(file, std::ios::binary);
  require(bool(in), "cannot read " + file);
  std::ostringstream s; s << in.rdbuf();
  require(!in.bad(), "read failed: " + file);
  return s.str();
}
std::vector<std::string> split(const std::string & line) {
  std::vector<std::string> v;
  size_t a=0, b;
  do {
    b=line.find('\t',a);
    v.push_back(Helper::trim(line.substr(a,b==std::string::npos ? b : b-a)));
    a=b+1;
  } while (b!=std::string::npos);
  return v;
}
int index(const std::string & s) {
  require(!s.empty() && s.find_first_not_of("0123456789")==std::string::npos, "invalid zero-based index: " + s);
  size_t used=0; long long n=std::stoll(s,&used);
  require(n<=std::numeric_limits<int>::max(), "index too large");
  return int(n);
}
void label(const std::string & s) {
  require(!s.empty() && s==Helper::trim(s) && s.find(',')==std::string::npos &&
          s.find_first_of("\r\n\t")==std::string::npos &&
          !Helper::imatch(s,"EDF Annotation",14), "illegal signal label: " + s);
}
void labels(const std::string & text, int n, bool parcels,
            std::vector<std::string> & names, std::vector<std::string> & ids) {
  std::istringstream in(text); std::string line;
  require(bool(std::getline(in,line)), "empty label TSV");
  auto h=split(line); std::map<std::string,int> columns;
  for (int i=0;i<h.size();++i) require(columns.emplace(h[i],i).second,"duplicate TSV header");
  const std::string key=parcels ? "row" : "column", name=parcels ? "luna_signal" : "channel";
  require(columns.count(key) && columns.count(name) && (!parcels || columns.count("parcel_id")), "missing TSV columns");
  names.resize(n); ids.resize(n); std::vector<bool> seen(n,false); std::set<std::string> unique, unique_ids;
  while (std::getline(in,line)) {
    auto v=split(line);
    require(v.size()==h.size(),"incorrect TSV field count");
    int i=index(v[columns[key]]);
    require(i<n && !seen[i],"duplicate or out-of-range TSV index"); seen[i]=true;
    names[i]=v[columns[name]]; label(names[i]);
    require(unique.insert(normalize(names[i])).second,"duplicate signal name: "+names[i]);
    if (parcels) {
      ids[i]=v[columns["parcel_id"]]; index(ids[i]);
      require(unique_ids.insert(std::to_string(index(ids[i]))).second,"duplicate parcel ID");
    }
  }
  require(std::find(seen.begin(),seen.end(),false)==seen.end(),"gapped TSV indices");
}
std::string relative(const std::string & base, const std::string & file) {
  namespace fs=std::filesystem;
  fs::path p(file);
  require(!file.empty() && !p.is_absolute(),"bundle path must be relative: "+file);
  auto root=fs::canonical(base), resolved=fs::weakly_canonical(root/p);
  auto r=resolved.lexically_relative(root);
  require(!r.empty() && *r.begin()!="..", "bundle path escapes directory: "+file);
  return resolved.string();
}
}
std::string normalize(const std::string & s) { return Helper::toupper(Helper::trim(s)); }
double voltage_scale(const std::string & unit) {
  const std::string s=normalize(unit);
  if (s=="V") return 1;
  if (s=="MV") return 1e-3;
  if (s=="UV" || s=="µV" || s=="μV") return 1e-6;
  throw std::runtime_error("unsupported voltage unit: "+unit+" (use input-unit if necessary)");
}
// Original portable SHA-256 implementation, GPL-3.0-or-later, like this file.
// The algorithm's constants and operations follow FIPS 180-4; no library dependency.
std::string sha256(const std::string & bytes) {
  static const uint32_t k[64]={
    0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
    0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
    0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
    0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
    0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
    0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
    0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
    0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2};
  auto ro=[](uint32_t x, int n){ return (x>>n)|(x<<(32-n)); };
  std::vector<uint8_t> data(bytes.begin(),bytes.end());
  uint64_t bits=uint64_t(data.size())*8;
  data.push_back(0x80); while(data.size()%64!=56) data.push_back(0);
  for(int i=7;i>=0;--i) data.push_back(uint8_t(bits>>(8*i)));
  std::array<uint32_t,8> h={0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19};
  for(size_t offset=0;offset<data.size();offset+=64) {
    uint32_t w[64];
    for(int i=0;i<16;++i) w[i]=(uint32_t(data[offset+4*i])<<24)|(uint32_t(data[offset+4*i+1])<<16)|(uint32_t(data[offset+4*i+2])<<8)|data[offset+4*i+3];
    for(int i=16;i<64;++i) w[i]=w[i-16]+(ro(w[i-15],7)^ro(w[i-15],18)^(w[i-15]>>3))+w[i-7]+(ro(w[i-2],17)^ro(w[i-2],19)^(w[i-2]>>10));
    auto a=h;
    for(int i=0;i<64;++i) {
      uint32_t t=a[7]+(ro(a[4],6)^ro(a[4],11)^ro(a[4],25))+((a[4]&a[5])^(~a[4]&a[6]))+k[i]+w[i];
      uint32_t u=(ro(a[0],2)^ro(a[0],13)^ro(a[0],22))+((a[0]&a[1])^(a[0]&a[2])^(a[1]&a[2]));
      for(int j=7;j>0;--j) a[j]=a[j-1];
      a[4]+=t; a[0]=t+u;
    }
    for(int i=0;i<8;++i) h[i]+=a[i];
  }
  std::ostringstream out; out<<std::hex<<std::setfill('0');
  for(auto v:h) out<<std::setw(8)<<v;
  return out.str();
}
operator_t load(paths_t p) {
  operator_t op;
  if (!p.bundle.empty()) {
    require(p.matrix.empty() && p.channels.empty() && p.parcels.empty() && p.metadata.empty(),"cannot mix bundle and explicit files");
    p.metadata=relative(p.bundle,"metadata.json");
    p.channels=relative(p.bundle,"channels.tsv"); p.parcels=relative(p.bundle,"parcels.tsv");
    p.manifest=relative(p.bundle,"manifest.json");
  }
  const auto metadata_text=read(p.metadata);
  op.metadata=nlohmann::json::parse(metadata_text);
  const auto & m=op.metadata;
  require(m.at("schema_version").is_number_integer() && m.at("schema_version")==1,"unsupported metadata schema_version");
  op.name=m.at("bundle_name").get<std::string>(); require(!op.name.empty(),"empty bundle_name");
  const auto & mx=m.at("matrix"); const auto & shape=mx.at("shape");
  require(shape.is_array() && shape.size()==2 && shape[0].is_number_integer() && shape[1].is_number_integer(),"matrix shape must contain two integers");
  require(shape[0]>0 && shape[1]>0 && shape[0]<=std::numeric_limits<int>::max() && shape[1]<=std::numeric_limits<int>::max(),"invalid matrix dimensions");
  int nr=shape[0].get<int>(), nc=shape[1].get<int>();
  if (!p.bundle.empty()) p.matrix=relative(p.bundle,mx.at("file").get<std::string>());
  else require(!p.matrix.empty() && !p.channels.empty() && !p.parcels.empty(),"matrix, channels, parcels and metadata are required");
  require(mx.at("format")=="headerless tab-separated float64, rows=parcels, columns=channels","unsupported matrix format");
  op.input_unit=mx.at("input_unit").get<std::string>(); voltage_scale(op.input_unit);
  op.output_unit=mx.at("output_unit").get<std::string>();
  require(!op.output_unit.empty() && op.output_unit.size()<=8 && op.output_unit.find_first_of("\t\r\n")==std::string::npos,"invalid output unit");
  const auto & ref=m.at("reference");
  op.reference=ref.at("required").get<std::string>();
  op.embedded_projector=ref.at("embedded_projector").get<bool>();
  op.average_reference=op.reference=="average" || op.reference=="average across all restored "+std::to_string(nc)+" scalp channels";
  require(op.average_reference || op.reference=="none","unsupported reference: "+op.reference);
  require(op.average_reference || !op.embedded_projector,"no-reference metadata cannot embed a reference projector");
  std::map<std::string,std::string> texts={{"matrix",read(p.matrix)},{"channels",read(p.channels)},{"parcels",read(p.parcels)},{"metadata",metadata_text}};
  std::map<std::string,std::string> files={{"matrix",p.matrix},{"channels",p.channels},{"parcels",p.parcels},{"metadata",p.metadata}};
  op.files=files;
  for (auto & v:texts) op.hashes[v.first]=sha256(v.second);
  // Both offline and compact manifests expose the same sha256 map. Only consumed files matter.
  if (p.manifest.empty()) p.manifest=(std::filesystem::path(p.metadata).parent_path()/"manifest.json").string();
  if (std::filesystem::exists(p.manifest)) {
    auto manifest=nlohmann::json::parse(read(p.manifest));
    const auto & hashes=manifest.at("sha256"); require(hashes.is_object(),"invalid manifest sha256 map");
    auto base=std::filesystem::path(p.manifest).parent_path();
    for (auto & f:files) {
      auto key=std::filesystem::absolute(f.second).lexically_normal().lexically_relative(std::filesystem::absolute(base).lexically_normal()).generic_string();
      require(hashes.contains(key),"manifest lacks consumed file: "+key);
      const std::string expected=hashes.at(key).get<std::string>();
      require(expected.size()==64 && normalize(expected)==normalize(op.hashes.at(f.first)),"checksum mismatch: "+f.second);
    }
  }
  std::vector<std::string> unused;
  labels(texts["channels"],nc,false,op.inputs,unused);
  labels(texts["parcels"],nr,true,op.outputs,op.parcel_ids);
  op.weights.resize(nr,nc); std::istringstream in(texts["matrix"]); std::string line; int row=0;
  while (std::getline(in,line)) {
    require(row<nr,"too many matrix rows"); auto v=split(line); require(v.size()==nc,"ragged matrix");
    for(int c=0;c<nc;++c) {
      size_t used=0; double x=std::stod(v[c],&used);
      require(used==v[c].size() && std::isfinite(x),"invalid/nonfinite matrix cell"); op.weights(row,c)=x;
    }
    ++row;
  }
  require(row==nr,"incorrect matrix row count");
  return op;
}
Eigen::MatrixXd multiply(const Eigen::MatrixXd & w, const Eigen::MatrixXd & x) {
  require(w.cols()==x.rows() && w.rows()>0 && w.cols()>0 && x.cols()>0,"kernel dimension mismatch");
  require(w.allFinite() && x.allFinite(),"nonfinite kernel input");
  Eigen::MatrixXd y=w*x; require(y.allFinite(),"nonfinite projection result"); return y;
}
namespace {
// Explicit range setup avoids add_signal's small-span expansion. Keep this
// local: projection needs fixed ranges before streaming records into signals.
struct range_t { double lo, hi, step, offset; };
range_t storage_range(double lo, double hi) {
  const double span=hi-lo, peak=std::max(std::fabs(lo),std::fabs(hi));
  require(std::isfinite(span),"output range overflow");
  // A constant zero needs no amplitude resolution; this range also round-trips
  // exactly through EDF's limited-width physical range fields.
  if (lo==0 && hi==0) return {-32768,32767,1,0};
  double pad=span>0 ? span/65530 : peak*1e-4;
  pad=std::max(pad,peak*16*std::numeric_limits<double>::epsilon());
  double a=lo-pad, b=hi+pad, step=(b-a)/65535;
  require(std::isfinite(a) && std::isfinite(b) && step>0 && std::isfinite(step),"unrepresentable output storage range");
  require(span==0 || step<=span/1000,"output variation cannot be resolved in 16-bit storage");
  require(step<=peak/1000,"output amplitude cannot be resolved in 16-bit storage");
  double offset=b/step-32767;
  require(std::isfinite(offset),"unrepresentable output storage offset");
  return {a,b,step,offset};
}
void run_impl(edf_t & edf, param_t & param) {
  const bool verbose = param.has("verbose");
  paths_t paths;
  for(auto pair: {std::make_pair("matrix",&paths.matrix),{"channels",&paths.channels},{"parcels",&paths.parcels},{"metadata",&paths.metadata},{"bundle",&paths.bundle}})
    if(param.has(pair.first)) *pair.second=param.requires(pair.first);
  operator_t op=load(paths);
  auto & h=edf.header;
  require(h.record_duration>0 && std::isfinite(h.record_duration),"unsupported record duration");
  std::map<std::string,std::vector<int>> available;
  for(int s=0;s<h.ns;++s) available[normalize(h.label[s])].push_back(s);
  std::vector<int> inputs;
  std::vector<double> scales;
  const int ni=op.inputs.size(), no=op.outputs.size();
  const double native_input_scale=voltage_scale(op.input_unit);
  for(auto & name:op.inputs) {
    auto match=available.find(normalize(name));
    require(match!=available.end(),"missing required channel: "+name);
    require(match->second.size()==1,"duplicate required channel: "+name);
    int s=match->second[0]; require(!h.is_annotation_channel(s),"annotation input: "+name);
    inputs.push_back(s);
    scales.push_back(voltage_scale(param.has("input-unit") ? param.requires("input-unit") : h.phys_dimension[s]));
  }
  const int samples=h.n_samples[inputs[0]];
  const double rate=h.sampling_freq(inputs[0]);
  require(samples>0 && rate>0 && std::isfinite(rate),"invalid sample count/rate");
  for(auto s:inputs) require(h.n_samples[s]==samples && h.sampling_freq(s)==rate,"input sample counts/rates do not match; resample first");
  // CHEP describes known bad channel/epoch combinations, not interpolation history.
  for(const auto & epoch:edf.timeline.make_chep_copy())
    for(const auto & bad:epoch.second)
      for(const auto & name:op.inputs) require(normalize(bad)!=normalize(name),"required channel has a CHEP bad-channel mask: "+name);
  std::string prefix=param.has("prefix") ? param.value("prefix") : "";
  std::set<std::string> outputs;
  for(auto & name:op.outputs) {
    name=prefix+name; label(name);
    require(outputs.insert(normalize(name)).second && !available.count(normalize(name)),"output label collision: "+name);
  }
  const std::string stored=param.has("output-unit") ? param.requires("output-unit") : op.output_unit=="A m" ? "nAm" : op.output_unit;
  double output_scale=1;
  if (stored!=op.output_unit) {
    if(op.output_unit=="A m" && stored=="nAm") output_scale=1e9;
    else if(op.output_unit=="nAm" && stored=="A m") output_scale=1e-9;
    else output_scale=voltage_scale(op.output_unit)/voltage_scale(stored);
  }
  require(!stored.empty() && stored.size()<=8,"invalid stored unit");
  std::vector<int> records;
  uint64_t previous_stop=0;
  for(int r=edf.timeline.first_record();r!=-1;r=edf.timeline.next_record(r)) {
    require(r>=0 && r<h.nr_all,"unsupported record index");
    // record2interval exposes inclusive end points in Luna.
    auto interval=edf.timeline.record2interval(r);
    require(interval.stop>=interval.start && interval.duration()==h.record_duration_tp-1 && (records.empty() || interval.start>previous_stop),"unsupported overlapping/inconsistent record timeline");
    previous_stop=interval.stop; records.push_back(r);
  }
  require(!records.empty() && records.size()==h.nr,"empty or inconsistent retained record timeline");
  auto block=[&](int r) {
    edf.ensure_loaded(r);
    Eigen::MatrixXd x(ni,samples);
    for(int c=0;c<ni;++c) {
      auto data=edf.records.at(r).get_pdata(inputs[c]);
      require(data.size()==samples,"inconsistent input record length");
      for(int t=0;t<samples;++t) x(c,t)=data[t]*scales[c];
    }
    require(x.allFinite(),"nonfinite input samples");
    if(op.average_reference) {
      Eigen::RowVectorXd mean=x.colwise().mean();
      x.rowwise()-=mean;
    }
    x/=native_input_scale;
    Eigen::MatrixXd y=multiply(op.weights,x); y*=output_scale;
    require(y.allFinite(),"output unit conversion overflow"); return y;
  };
  std::vector<double> lo(no,std::numeric_limits<double>::infinity()), hi(no,-std::numeric_limits<double>::infinity());
  for(auto r:records) {
    Eigen::MatrixXd y=block(r);
    for(int row=0;row<no;++row) { lo[row]=std::min(lo[row],y.row(row).minCoeff()); hi[row]=std::max(hi[row],y.row(row).maxCoeff()); }
  }
  std::vector<range_t> ranges;
  for(int row=0;row<no;++row) ranges.push_back(storage_range(lo[row],hi[row]));
  const double total=double(records.size())*samples;
  logger << "  PROJECT " << op.name << ": " << ni << " inputs, " << no << " outputs, " << rate << " Hz, " << records.size() << " records, " << total << " samples\n"
         << "  reference: " << op.reference << "; embedded projector: " << (op.embedded_projector ? "yes" : "no") << "\n"
         << "  native output unit " << op.output_unit << "; stored " << stored << "; scale " << output_scale << "\n";
  if(verbose) logger << "  model metadata: " << op.metadata.dump() << "\n"
                     << "  channel quality/interpolation history cannot be inferred from channel availability\n";
  if(param.has("input-unit")) logger << "  overriding all recording input units with " << param.value("input-unit") << "\n";
  writer.value("BUNDLE",op.name); writer.value("MODEL",op.metadata.dump());
  writer.value("N_INPUT",ni); writer.value("N_OUTPUT",no); writer.value("SR",rate);
  writer.value("N_RECORD",double(records.size())); writer.value("N_SAMPLE",total);
  writer.value("REFERENCE",op.reference); writer.value("EMBEDDED",int(op.embedded_projector));
  writer.value("CHECK",int(param.has("check")));
  for(auto & f:op.hashes) {
    if(verbose) logger << "  " << f.first << " " << op.files.at(f.first) << " SHA256 " << f.second << "\n";
    writer.level(f.first,"FILE"); writer.value("SHA256",f.second); writer.value("PATH",op.files.at(f.first)); writer.unlevel("FILE");
  }
  for(int c=0;c<ni;++c) {
    int s=inputs[c];
    if(verbose) logger << "  input " << c << ": " << op.inputs[c] << " -> [" << s << "] " << h.label[s] << "; header unit " << h.phys_dimension[s] << ", to V " << scales[c] << "\n";
    writer.level(c,"INPUT"); writer.value("REQUIRED",op.inputs[c]); writer.value("INDEX",s); writer.value("LABEL",h.label[s]);
    writer.value("HEADER_UNIT",h.phys_dimension[s]); writer.value("TO_V",scales[c]);
    writer.value("UNIT_OVERRIDE",param.has("input-unit") ? param.value("input-unit") : "none"); writer.unlevel("INPUT");
  }
  for(int row=0;row<no;++row) {
    writer.level(row,"OUTPUT"); writer.value("PARCEL_ID",op.parcel_ids[row]); writer.value("LABEL",op.outputs[row]);
    writer.value("NATIVE_UNIT",op.output_unit); writer.value("STORED_UNIT",stored); writer.value("SCALE",output_scale);
    writer.value("PMIN",ranges[row].lo); writer.value("PMAX",ranges[row].hi); writer.value("STEP",ranges[row].step); writer.unlevel("OUTPUT");
    if(verbose) logger << "  output " << row << " parcel " << op.parcel_ids[row] << " -> " << op.outputs[row] << "; " << stored << "; step " << ranges[row].step << "\n";
  }
  if(param.has("check")) { logger << "  PROJECT check complete; EDF unchanged\n"; return; }
  // No expected validation failures beyond this point. init_signal allocates
  // ordinary int16 records; explicit sample count avoids floating point truncation.
  const int first=h.ns;
  for(int row=0;row<no;++row) {
    require(edf.init_signal(op.outputs[row],rate,samples),"failed to initialize output");
    int s=first+row; auto range=ranges[row];
    h.phys_dimension[s]=stored;
    h.physical_min[s]=h.orig_physical_min[s]=range.lo;
    h.physical_max[s]=h.orig_physical_max[s]=range.hi;
    h.bitvalue[s]=range.step; h.offset[s]=range.offset;
  }
  for(auto r:records) {
    Eigen::MatrixXd y=block(r);
    for(int row=0;row<no;++row) {
      std::vector<double> data(samples);
      for(int t=0;t<samples;++t) data[t]=y(row,t);
      edf.update_records(r,r,first+row,&data);
    }
  }
}
}
void run(edf_t & edf, param_t & param) {
  try { run_impl(edf,param); }
  catch(const std::exception & e) { Helper::halt(std::string("PROJECT: ")+e.what()); }
}
}
