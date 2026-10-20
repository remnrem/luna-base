//    --------------------------------------------------------------------
//    Native, text-first reference PCA for HDSTATS feature lines.
//    --------------------------------------------------------------------

#include "hdstats-pc.h"

#include "param.h"
#include "helper/helper.h"
#include "helper/logger.h"
#include "helper/luna_io.h"
#include "db/db.h"
#include "stats/Eigen/Dense"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <mutex>
#include <set>
#include <sstream>

extern logger_t logger;
extern writer_t writer;

namespace {

struct format_t {
  const char * model_magic;
  const char * feature_magic;
  const char * extension;
  const char * label;
  const char * score_prefix;
  const char * command;
};

const format_t hdstats_format = {
  "__HDSTATS_PC_V2", "__HDSTATS_FEATURES_V2", ".hdpc", "HDSTATS-PC", "HDPC", "--hdstats-pc"
};
const format_t hypno_format = {
  "__HYPNO_PC_V1", "__HYPNO_FEATURES_V1", ".hypc", "HYPNO-PC", "HYPC", "--hypno-pc"
};

struct row_t {
  std::string id;
  hdstats_pc::feature_map_t values;
};

struct model_t {
  std::vector<std::string> names;
  Eigen::VectorXd means;
  Eigen::VectorXd sds;
  Eigen::MatrixXd loadings;
  Eigen::VectorXd singular_values;
};

struct panel_def_t {
  std::string name;
  std::string option;
  std::vector<std::string> names;
};

struct fitted_panel_t {
  std::string name;
  std::vector<std::string> declared;
  model_t model;
  int n_input = 0;
  int n_fit = 0;
  double total_ss = 0.0;
};

// This is the executable schema. pops/hdstats.pc is its human-readable guide,
// never an input to fitting or projection.
const std::vector<panel_def_t> & panel_definitions()
{
  static const std::vector<panel_def_t> panels = [] {
    std::vector<panel_def_t> out;
    auto panel = [&]( const std::string & name, const std::string & option ) -> panel_def_t & {
      out.push_back( { name, option, {} } );
      return out.back();
    };
    auto add = []( panel_def_t & p, const std::string & stratum, const std::string & var ) {
      p.names.push_back( hdstats_pc::key( var, stratum ) );
    };
    const std::vector<std::string> nss = { "5", "3" };
    panel_def_t & ambiguity = panel( "AMBIGUITY", "nc-ambiguity" );
    for ( const auto & n : nss )
      for ( const auto & v : { "H", "H_P90", "CONF", "FRAC_LOW_CONF", "MG" } )
        add( ambiguity, "NSS/" + n, v );
    for ( const auto & mix : std::vector<std::string>{ "W_N1", "W_N2", "W_N3", "W_R",
          "N1_N2", "N1_N3", "N1_R", "N2_N3", "N2_R", "N3_R", "N1_N2_N3", "W_N1_R" } )
      for ( const auto & v : { "HD_MINS", "HD_DENS" } ) add( ambiguity, "NSS/5;MIX/" + mix, v );
    for ( const auto & mix : { "W_NR", "W_R", "NR_R", "W_NR_R" } )
      for ( const auto & v : { "HD_MINS", "HD_DENS" } ) add( ambiguity, "NSS/3;MIX/" + std::string(mix), v );

    panel_def_t & instability = panel( "INSTABILITY_STABILITY", "nc-instability-stability" );
    for ( const auto & n : nss )
      {
        const std::string s = "NSS/" + n;
        for ( const auto & v : { "PS", "PS_P90", "PS_LAG", "CORR_H_PS", "PCT_STABLE",
              "H_RATIO_UN_ST", "CONF_DIFF_UN_ST", "PS_RATIO_UN_ST" } ) add( instability, s, v );
        for ( const auto & stable : { "1", "0" } )
          for ( const auto & v : { "H", "CONF", "PS" } ) add( instability, s + ";STABLE/" + stable, v );
      }

    panel_def_t & transitions = panel( "TRANSITIONS", "nc-transitions" );
    const std::vector<std::string> trans5 = { "W", "N1", "N2", "N3", "R" };
    const std::vector<std::string> trans3 = { "W", "NR", "R" };
    for ( const auto & n : nss )
      {
        const std::string s = "NSS/" + n;
        for ( const auto & v : { "TRANS_N", "TRANS_DENS", "BOUT_TRANS_N", "BOUT_TRANS_DENS",
              "TRANS_N_MATCHED", "TRANS_PCT_MATCHED", "TRANS_DUR", "TRANS_PEAK_H",
              "TRANS_MIN_CONF", "TRANS_TOT_PS", "BOUT_TRANS_DUR", "BOUT_TRANS_PEAK_H",
              "BOUT_TRANS_MIN_CONF", "BOUT_TRANS_TOT_PS" } ) add( transitions, s, v );
        const auto & stages = n == "5" ? trans5 : trans3;
        for ( const auto & from : stages ) for ( const auto & to : stages ) if ( from != to )
          {
            add( transitions, s + ";TRANS/" + from + "to" + to, "TRANS_DENS" );
            // N1 rarely qualifies as a settled bout. Its directional bout
            // densities were constant or nearly so across three cohorts.
            if ( n != "5" || ( from != "N1" && to != "N1" ) )
              add( transitions, s + ";BOUT_TRANS/" + from + "to" + to, "BOUT_TRANS_DENS" );
          }
      }

    panel_def_t & architecture = panel( "SLEEP_ARCHITECTURE", "nc-sleep-architecture" );
    for ( const auto & v : { "HD_SPT", "HD_WASO", "HD_SOL", "HD_SME", "HD_REM_LAT", "HD_REM_LAT2" } )
      add( architecture, ".", v );
    for ( const auto & st : trans5 )
      for ( const auto & v : { "HD_MINS", "HD_PCT", "HD_DENS" } ) add( architecture, "SS/" + st, v );
    const std::vector<std::pair<std::string,std::string>> timing = {
      {"3","W"}, {"5","N1"}, {"5","N2"}, {"5","N3"}, {"3","R"},
      {"5","N3_NREM"}, {"5","R_SLEEP"} };
    for ( const auto & t : timing ) for ( const auto & v : { "T50", "COM" } )
      add( architecture, "NSS/" + t.first + ";LRVAR/" + t.second, v );

    panel_def_t & ultradian = panel( "ULTRADIAN", "nc-ultradian" );
    for ( const auto & v : { "U_MAG", "U_COVERAGE", "U_NCYC" } ) add( ultradian, ".", v );
    for ( const auto & metric : { "CONF", "CONF3", "PS", "PS3" } )
      {
        for ( const auto & q : { "D1", "D2", "L1", "L2" } )
          add( ultradian, "UQ/" + std::string(q) + ";U_METRIC/" + metric, "U_MEAN" );
        for ( const auto & v : { "U_SLOPE", "U_MOD", "U_MOD_COS", "U_MOD_SIN", "U_R2", "U_PARTIAL_R2" } )
          add( ultradian, "UVAR/" + std::string(metric), v );
      }

    panel_def_t & overnight = panel( "OVERNIGHT", "nc-overnight" );
    for ( const auto & n : nss )
      {
        const std::string s = "NSS/" + n;
        for ( const auto & lag : { "600", "1800", "3600", "5400" } )
          for ( const auto & v : { "PS", "SIM" } ) add( overnight, s + ";LRLAG/" + lag, v );
        for ( const auto & metric : { "H", "PS" } )
          {
            add( overnight, s + ";LRVAR/" + metric, "SLOPE" );
            for ( const auto & win : { "1800", "3600" } )
              for ( const auto & v : { "WITHIN_SD", "BETWEEN_SD" } )
                add( overnight, s + ";LRVAR/" + metric + ";LRWIN/" + win, v );
          }
      }
    for ( const auto & t : std::vector<std::pair<std::string,std::string>>{
          {"3","W"},{"5","N1"},{"5","N2"},{"5","N3"},{"3","R"} } )
      add( overnight, "NSS/" + t.first + ";LRVAR/" + t.second, "SLOPE" );

    for ( const auto & t : std::vector<std::pair<std::string,std::string>>{
          {"3","W"},{"5","N1"},{"5","N2"},{"5","N3"},{"3","R"} } )
      {
        panel_def_t & stage = panel( "STAGE_" + t.second, "nc-stage-" + std::string(t.second == "W" ? "w" :
          t.second == "R" ? "r" : t.second == "N1" ? "n1" : t.second == "N2" ? "n2" : "n3") );
        for ( const auto & v : { "HD_MINS", "HD_PCT", "HD_DENS" } ) add( stage, "SS/" + t.second, v );
        // N1/N3 argmax-conditioned summaries are unavailable for many
        // otherwise usable recordings; keep them in HDSTATS output only.
        if ( t.second != "N1" && t.second != "N3" )
          for ( const auto & v : { "P_", "PS" } )
            add( stage, "SS/" + t.second + ";NSS/" + t.first,
                 std::string(v) == "P_" ? "P_" + t.second : v );
        for ( const auto & lag : { "600", "1800", "3600", "5400" } )
          add( stage, "SS/" + t.second + ";NSS/" + t.first + ";LRLAG/" + lag, "ACF" );
        for ( const auto & v : { "SLOPE", "T50", "COM" } )
          add( stage, "NSS/" + t.first + ";LRVAR/" + t.second, v );
        for ( const auto & win : { "1800", "3600" } )
          for ( const auto & v : { "WITHIN_SD", "BETWEEN_SD" } )
            add( stage, "NSS/" + t.first + ";LRVAR/" + t.second + ";LRWIN/" + win, v );
      }
    const std::vector<int> expected = {42,28,72,35,43,41,16,14,16,14,16};
    for ( int i = 0; i < (int)out.size(); ++i )
      if ( (int)out[i].names.size() != expected[i] ||
           std::set<std::string>(out[i].names.begin(),out[i].names.end()).size() != out[i].names.size() )
        Helper::halt( "internal HDSTATS panel schema error: " + out[i].name );
    return out;
  }();
  return panels;
}

std::vector<std::string> split_tab( const std::string & s )
{
  std::vector<std::string> out;
  std::string item;
  std::stringstream ss( s );
  while ( std::getline( ss, item, '\t' ) ) out.push_back( item );
  return out;
}

std::string safe_id( const std::string & id )
{
  std::string out;
  for ( const char c : id )
    out += std::isalnum( static_cast<unsigned char>( c ) ) || c == '.' || c == '_' || c == '-'
      ? c : '_';
  return out.empty() ? "_" : out;
}

std::string fmt( const double x )
{
  if ( ! std::isfinite( x ) ) return "NA";
  std::ostringstream ss;
  ss << std::setprecision( 17 ) << x;
  return ss.str();
}

bool parse_value( const std::string & s, double * x )
{
  if ( s == "NA" || s == "." || s.empty() )
    {
      *x = std::numeric_limits<double>::quiet_NaN();
      return true;
    }
  return Helper::str2dbl( s, x );
}

row_t parse_row( const std::string & line, const std::string & source, const int line_no,
                 const format_t & format )
{
  const std::vector<std::string> tok = split_tab( line );
  if ( tok.size() < 4 || tok.size() % 2 )
    Helper::halt( "bad " + std::string( format.label ) + " feature line in " + source + ":" + Helper::int2str( line_no ) );
  row_t row;
  for ( int i = 0; i < (int)tok.size(); i += 2 )
    {
      const std::string & key = tok[i];
      const std::string & value = tok[i+1];
      if ( key == "ID" ) row.id = value;
      else if ( key == "FORMAT" )
        {
          if ( value != format.feature_magic )
            Helper::halt( "unsupported " + std::string( format.label ) + " feature format in " + source + ":" + Helper::int2str( line_no ) );
        }
      else
        {
          double x;
          if ( ! parse_value( value, &x ) )
            Helper::halt( "non-numeric feature '" + key + "' in " + source + ":" + Helper::int2str( line_no ) );
          if ( row.values.find( key ) != row.values.end() )
            Helper::halt( "duplicate feature '" + key + "' in " + source + ":" + Helper::int2str( line_no ) );
          row.values[key] = x;
        }
    }
  if ( row.id.empty() ) Helper::halt( "missing ID in " + source + ":" + Helper::int2str( line_no ) );
  return row;
}

std::vector<row_t> read_rows( const std::vector<std::string> & files, const format_t & format )
{
  std::vector<row_t> rows;
  std::set<std::string> ids;
  for ( const std::string & file : files )
    {
      std::ifstream in = LunaIO::open_ifstream( file, std::ios::in );
      if ( ! in ) Helper::halt( "could not open " + std::string( format.label ) + " feature data " + file );
      std::string line;
      int line_no = 0;
      while ( std::getline( in, line ) )
        {
          ++line_no;
          if ( line.empty() ) continue;
          row_t row = parse_row( line, file, line_no, format );
          if ( ! ids.insert( row.id ).second ) Helper::halt( "duplicate ID in " + std::string( format.label ) + " input: " + row.id );
          rows.push_back( std::move( row ) );
        }
    }
  if ( rows.empty() ) Helper::halt( "no " + std::string( format.label ) + " feature rows read" );
  return rows;
}

void write_model( const std::string & file, const model_t & m, const format_t & format )
{
  std::ofstream out = LunaIO::open_ofstream( file, std::ios::out );
  if ( ! out ) Helper::halt( "could not write " + std::string( format.label ) + " model " + file );
  out << format.model_magic << "\n";
  out << "NF\t" << m.names.size() << "\n";
  for ( int j = 0; j < (int)m.names.size(); ++j )
    out << "F\t" << m.names[j] << "\t" << fmt( m.means(j) ) << "\t" << fmt( m.sds(j) ) << "\n";
  out << "NC\t" << m.loadings.cols() << "\n";
  out << "W";
  for ( int k = 0; k < m.singular_values.size(); ++k ) out << "\t" << fmt( m.singular_values(k) );
  out << "\n";
  for ( int j = 0; j < (int)m.names.size(); ++j )
    {
      out << "V\t" << m.names[j];
      for ( int k = 0; k < m.loadings.cols(); ++k ) out << "\t" << fmt( m.loadings(j,k) );
      out << "\n";
    }
  if ( ! out ) Helper::halt( "error writing " + std::string( format.label ) + " model " + file );
}

void write_correlation_matrix( const std::string & file,
                               const std::vector<std::string> & names,
                               const Eigen::MatrixXd & z,
                               const format_t & format )
{
  std::ofstream out = LunaIO::open_ofstream( file, std::ios::out );
  if ( ! out ) Helper::halt( "could not write " + std::string( format.label ) + " correlation matrix " + file );

  // z is column-standardized using the sample SD, so this is the Pearson
  // correlation matrix for precisely the people and features retained by PCA.
  Eigen::MatrixXd corr = z.transpose() * z / ( z.rows() - 1 );
  corr.diagonal().setOnes();

  out << "VAR";
  for ( const std::string & name : names ) out << "\t" << name;
  out << "\n";
  for ( int i = 0; i < corr.rows(); ++i )
    {
      out << names[i];
      for ( int j = 0; j < corr.cols(); ++j ) out << "\t" << fmt( corr(i,j) );
      out << "\n";
    }
  if ( ! out ) Helper::halt( "error writing " + std::string( format.label ) + " correlation matrix " + file );
}

model_t read_model( const std::string & file, const format_t & format )
{
  std::ifstream in = LunaIO::open_ifstream( file, std::ios::in );
  if ( ! in ) Helper::halt( "could not open " + std::string( format.label ) + " model " + file );
  std::string line;
  if ( ! std::getline( in, line ) || line != format.model_magic ) Helper::halt( "bad " + std::string( format.label ) + " model " + file );
  if ( ! std::getline( in, line ) ) Helper::halt( "truncated " + std::string( format.label ) + " model " + file );
  std::vector<std::string> tok = split_tab( line );
  if ( tok.size() != 2 || tok[0] != "NF" ) Helper::halt( "bad " + std::string( format.label ) + " feature header " + file );
  int nf = 0;
  if ( ! Helper::str2int( tok[1], &nf ) || nf < 1 ) Helper::halt( "bad " + std::string( format.label ) + " feature count " + file );
  model_t m;
  m.names.resize( nf ); m.means.resize( nf ); m.sds.resize( nf );
  for ( int j = 0; j < nf; ++j )
    {
      if ( ! std::getline( in, line ) ) Helper::halt( "truncated " + std::string( format.label ) + " model " + file );
      tok = split_tab( line );
      if ( tok.size() != 4 || tok[0] != "F" || ! parse_value( tok[2], &m.means(j) ) || ! parse_value( tok[3], &m.sds(j) ) || m.sds(j) <= 0 )
        Helper::halt( "bad " + std::string( format.label ) + " feature row " + file );
      m.names[j] = tok[1];
    }
  if ( ! std::getline( in, line ) ) Helper::halt( "truncated " + std::string( format.label ) + " model " + file );
  tok = split_tab( line );
  int nc = 0;
  if ( tok.size() != 2 || tok[0] != "NC" || ! Helper::str2int( tok[1], &nc ) || nc < 1 ) Helper::halt( "bad " + std::string( format.label ) + " component header " + file );
  if ( ! std::getline( in, line ) ) Helper::halt( "truncated " + std::string( format.label ) + " model " + file );
  tok = split_tab( line );
  if ( (int)tok.size() != nc + 1 || tok[0] != "W" ) Helper::halt( "bad " + std::string( format.label ) + " singular values " + file );
  m.singular_values.resize( nc );
  for ( int k = 0; k < nc; ++k ) if ( ! parse_value( tok[k+1], &m.singular_values(k) ) ) Helper::halt( "bad " + std::string( format.label ) + " singular value " + file );
  m.loadings.resize( nf, nc );
  for ( int j = 0; j < nf; ++j )
    {
      if ( ! std::getline( in, line ) ) Helper::halt( "truncated " + std::string( format.label ) + " model " + file );
      tok = split_tab( line );
      if ( (int)tok.size() != nc + 2 || tok[0] != "V" || tok[1] != m.names[j] ) Helper::halt( "bad " + std::string( format.label ) + " loading row " + file );
      for ( int k = 0; k < nc; ++k ) if ( ! parse_value( tok[k+2], &m.loadings(j,k) ) ) Helper::halt( "bad " + std::string( format.label ) + " loading " + file );
    }
  return m;
}

const model_t & cached_model( const std::string & file, const format_t & format )
{
  // HDSTATS is normally invoked once for every recording in a Luna run.  A
  // fixed model must therefore be parsed once per filename, not once per ID.
  static std::map<std::string,model_t> cache;
  static std::mutex cache_mutex;
  const std::string path = std::string( format.model_magic ) + "\t" + Helper::expand( file );
  std::lock_guard<std::mutex> lock( cache_mutex );
  const auto found = cache.find( path );
  if ( found != cache.end() ) return found->second;
  return cache.emplace( path, read_model( Helper::expand( file ), format ) ).first->second;
}

}

std::string hdstats_pc::key( const std::string & variable, const std::string & stratum )
{
  if ( stratum == "." || stratum.empty() ) return variable;
  std::string out = variable + "__";
  for ( const char c : stratum )
    {
      if ( c == ';' ) out += "__";
      else out += c == '/' ? '_' : c;
    }
  return out;
}

void write_features_impl( const std::string & directory,
                          const std::string & id,
                          const hdstats_pc::feature_map_t & features,
                          const format_t & format )
{
  const std::filesystem::path dir( Helper::expand( directory ) );
  std::error_code ec;
  std::filesystem::create_directories( dir, ec );
  if ( ec ) Helper::halt( "could not create " + std::string( format.label ) + " feature directory " + directory );
  const std::filesystem::path file = dir / ( safe_id( id ) + format.extension );
  std::ofstream out = LunaIO::open_ofstream( file.string(), std::ios::out );
  if ( ! out ) Helper::halt( "could not write " + std::string( format.label ) + " feature file " + file.string() );
  out << "ID\t" << id << "\tFORMAT\t" << format.feature_magic;
  for ( const auto & f : features ) out << "\t" << f.first << "\t" << fmt( f.second );
  out << "\n";
  if ( ! out ) Helper::halt( "error writing " + std::string( format.label ) + " feature file " + file.string() );
}

bool project_impl( const std::string & file,
                   const hdstats_pc::feature_map_t & features,
                   std::vector<double> * scores,
                   int * missing,
                   const format_t & format )
{
  const model_t & m = cached_model( file, format );
  Eigen::VectorXd z( m.names.size() );
  int nmissing = 0;
  for ( int j = 0; j < (int)m.names.size(); ++j )
    {
      const auto f = features.find( m.names[j] );
      if ( f == features.end() || ! std::isfinite( f->second ) ) { ++nmissing; z(j) = 0.0; }
      else z(j) = ( f->second - m.means(j) ) / m.sds(j);
    }
  if ( missing != NULL ) *missing = nmissing;
  if ( nmissing ) return false;
  const Eigen::VectorXd s = m.loadings.transpose() * z;
  if ( scores != NULL ) scores->assign( s.data(), s.data() + s.size() );
  return true;
}

void fit_impl( param_t & param, const format_t & format )
{
  const std::vector<std::string> files = param.strvector( "data" );
  if ( files.empty() ) Helper::halt( std::string( format.command ) + " requires data=<concatenated feature file(s)>" );
  const std::string model_file = param.requires( "model" );
  const std::string corr_file = param.has( "corr-matrix" )
    ? param.requires( "corr-matrix" ) : "";
  const int requested_nc = param.has( "nc" ) ? param.requires_int( "nc" ) : 10;
  if ( requested_nc < 1 ) Helper::halt( "nc must be positive" );
  const std::vector<row_t> rows = read_rows( files, format );

  std::set<std::string> name_set;
  for ( const row_t & row : rows ) for ( const auto & f : row.values ) name_set.insert( f.first );
  std::vector<std::string> names( name_set.begin(), name_set.end() );
  Eigen::MatrixXd raw( rows.size(), names.size() );
  std::vector<int> retained;
  for ( int i = 0; i < (int)rows.size(); ++i )
    {
      bool complete = true;
      for ( int j = 0; j < (int)names.size(); ++j )
        {
          const auto f = rows[i].values.find( names[j] );
          raw(i,j) = f == rows[i].values.end() ? std::numeric_limits<double>::quiet_NaN() : f->second;
          if ( ! std::isfinite( raw(i,j) ) ) complete = false;
        }
      if ( complete ) retained.push_back( i );
    }
  if ( retained.size() < 2 ) Helper::halt( std::string( format.label ) + " needs at least two complete feature rows" );
  Eigen::MatrixXd x( retained.size(), names.size() );
  for ( int i = 0; i < (int)retained.size(); ++i ) x.row(i) = raw.row( retained[i] );
  Eigen::VectorXd means = x.colwise().mean();
  Eigen::VectorXd sds( names.size() );
  for ( int j = 0; j < x.cols(); ++j )
    {
      const double ss = ( x.col(j).array() - means(j) ).square().sum();
      sds(j) = std::sqrt( ss / ( x.rows() - 1 ) );
    }
  std::vector<int> keep;
  for ( int j = 0; j < sds.size(); ++j ) if ( std::isfinite( sds(j) ) && sds(j) > 0 ) keep.push_back( j );
  if ( keep.empty() ) Helper::halt( std::string( format.label ) + " has no non-constant features" );
  Eigen::MatrixXd z( x.rows(), keep.size() );
  model_t m;
  m.names.resize( keep.size() ); m.means.resize( keep.size() ); m.sds.resize( keep.size() );
  for ( int j = 0; j < (int)keep.size(); ++j )
    {
      m.names[j] = names[keep[j]]; m.means(j) = means(keep[j]); m.sds(j) = sds(keep[j]);
      z.col(j) = ( x.col(keep[j]).array() - m.means(j) ) / m.sds(j);
    }
  if ( ! corr_file.empty() )
    {
      write_correlation_matrix( corr_file, m.names, z, format );
      logger << "  " << format.label << ": wrote " << m.names.size()
             << " x " << m.names.size() << " feature correlation matrix for "
             << retained.size() << " people to " << corr_file << "\n";
    }
  const int nc = std::min( requested_nc, std::min( (int)z.cols(), (int)z.rows() - 1 ) );
  if ( nc < 1 ) Helper::halt( "not enough " + std::string( format.label ) + " observations/features" );
  Eigen::BDCSVD<Eigen::MatrixXd> svd( z, Eigen::ComputeThinU | Eigen::ComputeThinV );
  m.loadings = svd.matrixV().leftCols( nc );
  m.singular_values = svd.singularValues().head( nc );
  const double total_ss = svd.singularValues().array().square().sum();
  for ( int k = 0; k < nc; ++k )
    {
      Eigen::Index anchor;
      m.loadings.col(k).cwiseAbs().maxCoeff( &anchor );
      if ( m.loadings(anchor,k) < 0 ) m.loadings.col(k) *= -1.0;
    }
  write_model( model_file, m, format );
  logger << "  " << format.label << ": fitted " << nc << " component(s) from " << retained.size()
         << "/" << rows.size() << " complete people and " << m.names.size() << " features\n";
  writer.value( "N", (int)retained.size() );
  writer.value( "N_INPUT", (int)rows.size() );
  writer.value( "N_FEATURES", (int)m.names.size() );
  writer.value( "N_COMPONENTS", nc );

  for ( int j = 0; j < (int)m.names.size(); ++j )
    {
      writer.level( m.names[j], "FTR" );
      writer.value( "MEAN", m.means(j) );
      writer.value( "SD", m.sds(j) );
      writer.unlevel( "FTR" );
    }

  double cumulative_ve = 0.0;
  for ( int k = 0; k < nc; ++k )
    {
      const std::string pc = std::string( format.score_prefix ) + std::string( k < 9 ? "0" : "" ) + Helper::int2str( k + 1 );
      const double ve = total_ss > 0 ? m.singular_values(k) * m.singular_values(k) / total_ss : 0.0;
      cumulative_ve += ve;
      writer.level( pc, "PC" );
      writer.value( "W", m.singular_values(k) );
      writer.value( "VE", ve );
      writer.value( "CVE", cumulative_ve );
      for ( int j = 0; j < (int)m.names.size(); ++j )
        {
          writer.level( m.names[j], "FTR" );
          writer.value( "LOADING", m.loadings(j,k) );
          writer.unlevel( "FTR" );
        }
      writer.unlevel( "PC" );
    }
}

static fitted_panel_t fit_panel( const panel_def_t & def,
                                 const std::vector<row_t> & rows,
                                 const int requested_nc )
{
  fitted_panel_t panel;
  panel.name = def.name;
  panel.declared = def.names;
  panel.n_input = rows.size();
  Eigen::MatrixXd raw( rows.size(), def.names.size() );
  std::vector<int> retained;
  for ( int i = 0; i < (int)rows.size(); ++i )
    {
      bool complete = true;
      for ( int j = 0; j < (int)def.names.size(); ++j )
        {
          const auto f = rows[i].values.find( def.names[j] );
          raw(i,j) = f == rows[i].values.end() ? std::numeric_limits<double>::quiet_NaN() : f->second;
          if ( !std::isfinite(raw(i,j)) ) complete = false;
        }
      if ( complete ) retained.push_back(i);
    }
  panel.n_fit = retained.size();
  if ( panel.n_fit < 2 )
    Helper::halt( "HDSTATS-PC panel " + def.name + " has fewer than two complete people (" +
                  Helper::int2str(panel.n_fit) + "/" + Helper::int2str(panel.n_input) + ")" );
  Eigen::MatrixXd x( panel.n_fit, def.names.size() );
  for ( int i = 0; i < panel.n_fit; ++i ) x.row(i) = raw.row(retained[i]);
  const Eigen::VectorXd means = x.colwise().mean();
  std::vector<int> keep;
  Eigen::VectorXd sds( def.names.size() );
  for ( int j = 0; j < x.cols(); ++j )
    {
      sds(j) = std::sqrt((x.col(j).array() - means(j)).square().sum() / (x.rows() - 1));
      if ( std::isfinite(sds(j)) && sds(j) > 0 ) keep.push_back(j);
    }
  if ( keep.empty() ) Helper::halt( "HDSTATS-PC panel " + def.name + " has no varying features" );
  model_t & m = panel.model;
  Eigen::MatrixXd z( x.rows(), keep.size() );
  m.names.resize(keep.size()); m.means.resize(keep.size()); m.sds.resize(keep.size());
  for ( int j = 0; j < (int)keep.size(); ++j )
    {
      const int k = keep[j];
      m.names[j] = def.names[k]; m.means(j) = means(k); m.sds(j) = sds(k);
      z.col(j) = (x.col(k).array() - means(k)) / sds(k);
    }
  const int max_nc = std::min((int)z.cols(),(int)z.rows()-1);
  if ( requested_nc > max_nc )
    Helper::halt( "HDSTATS-PC panel " + def.name + " requested " + Helper::int2str(requested_nc) +
                  " components but supports at most " + Helper::int2str(max_nc) );
  Eigen::BDCSVD<Eigen::MatrixXd> svd(z,Eigen::ComputeThinU | Eigen::ComputeThinV);
  const Eigen::VectorXd singular = svd.singularValues();
  const double tolerance = singular.size() ? singular(0) * 1e-10 : 0.0;
  if ( singular(requested_nc-1) <= tolerance )
    Helper::halt( "HDSTATS-PC panel " + def.name + " has insufficient numerical rank for nc=" +
                  Helper::int2str(requested_nc) );
  m.loadings = svd.matrixV().leftCols(requested_nc);
  m.singular_values = singular.head(requested_nc);
  panel.total_ss = singular.array().square().sum();
  for ( int k = 0; k < requested_nc; ++k )
    {
      Eigen::Index anchor;
      m.loadings.col(k).cwiseAbs().maxCoeff(&anchor);
      if ( m.loadings(anchor,k) < 0 ) m.loadings.col(k) *= -1.0;
    }
  return panel;
}

static void write_panels( const std::string & file, const std::vector<fitted_panel_t> & panels )
{
  std::ofstream out = LunaIO::open_ofstream(file,std::ios::out);
  if ( !out ) Helper::halt("could not write HDSTATS-PC model " + file);
  out << "__HDSTATS_PC_V3\nNP\t" << panels.size() << "\n";
  for ( const auto & p : panels )
    {
      const model_t & m = p.model;
      out << "PANEL\t" << p.name << "\t" << p.n_input << "\t" << p.n_fit << "\t"
          << p.declared.size() << "\t" << m.names.size() << "\t" << m.loadings.cols()
          << "\t" << fmt(p.total_ss) << "\n";
      for ( const auto & name : p.declared ) out << "I\t" << name << "\n";
      for ( int j = 0; j < (int)m.names.size(); ++j )
        out << "F\t" << m.names[j] << "\t" << fmt(m.means(j)) << "\t" << fmt(m.sds(j)) << "\n";
      out << "W";
      for ( int k = 0; k < m.singular_values.size(); ++k ) out << "\t" << fmt(m.singular_values(k));
      out << "\n";
      for ( int j = 0; j < (int)m.names.size(); ++j )
        {
          out << "V\t" << m.names[j];
          for ( int k = 0; k < m.loadings.cols(); ++k ) out << "\t" << fmt(m.loadings(j,k));
          out << "\n";
        }
      out << "END\n";
    }
  if ( !out ) Helper::halt("error writing HDSTATS-PC model " + file);
}

static std::vector<fitted_panel_t> read_panels( const std::string & file )
{
  std::ifstream in = LunaIO::open_ifstream(file,std::ios::in);
  std::string line;
  if ( !std::getline(in,line) || line != "__HDSTATS_PC_V3" ) Helper::halt("bad HDSTATS-PC V3 model " + file);
  if ( !std::getline(in,line) ) Helper::halt("truncated HDSTATS-PC model " + file);
  auto tok = split_tab(line);
  int np;
  if ( tok.size()!=2 || tok[0]!="NP" || !Helper::str2int(tok[1],&np) || np<1 )
    Helper::halt("bad HDSTATS-PC panel count " + file);
  std::vector<fitted_panel_t> panels;
  std::set<std::string> seen_panels;
  for ( int p = 0; p < np; ++p )
    {
      if ( !std::getline(in,line) ) Helper::halt("truncated HDSTATS-PC panel " + file);
      tok = split_tab(line);
      int n_input,n_fit,nd,nf,nc;
      double total_ss;
      if ( tok.size()!=8 || tok[0]!="PANEL" || tok[1].empty() || !seen_panels.insert(tok[1]).second ||
           !Helper::str2int(tok[2],&n_input) || !Helper::str2int(tok[3],&n_fit) ||
           !Helper::str2int(tok[4],&nd) || !Helper::str2int(tok[5],&nf) ||
           !Helper::str2int(tok[6],&nc) || !parse_value(tok[7],&total_ss) ||
           n_input<2 || n_fit<2 || nd<1 || nf<1 || nc<1 || nc>nf ||
           !std::isfinite(total_ss) || total_ss<=0 ) Helper::halt("bad HDSTATS-PC panel header " + file);
      fitted_panel_t panel;
      panel.name=tok[1]; panel.n_input=n_input; panel.n_fit=n_fit; panel.total_ss=total_ss;
      for ( int j=0; j<nd; ++j )
        {
          if ( !std::getline(in,line) ) Helper::halt("truncated HDSTATS-PC inputs " + file);
          tok=split_tab(line);
          if ( tok.size()!=2 || tok[0]!="I" || tok[1].empty() ||
               std::find(panel.declared.begin(),panel.declared.end(),tok[1])!=panel.declared.end() )
            Helper::halt("bad HDSTATS-PC input schema " + file);
          panel.declared.push_back(tok[1]);
        }
      model_t & m=panel.model;
      m.names.resize(nf);m.means.resize(nf);m.sds.resize(nf);
      for ( int j=0; j<nf; ++j )
        {
          if ( !std::getline(in,line) ) Helper::halt("truncated HDSTATS-PC features " + file);
          tok=split_tab(line);
          if ( tok.size()!=4 || tok[0]!="F" || !parse_value(tok[2],&m.means(j)) ||
               !parse_value(tok[3],&m.sds(j)) || !std::isfinite(m.means(j)) ||
               !std::isfinite(m.sds(j)) || m.sds(j)<=0 ) Helper::halt("bad HDSTATS-PC feature " + file);
          m.names[j]=tok[1];
          if ( std::find(panel.declared.begin(),panel.declared.end(),m.names[j])==panel.declared.end() )
            Helper::halt("HDSTATS-PC undeclared feature " + m.names[j]);
        }
      if ( !std::getline(in,line) ) Helper::halt("truncated HDSTATS-PC singular values " + file);
      tok=split_tab(line);m.singular_values.resize(nc);
      if ( (int)tok.size()!=nc+1 || tok[0]!="W" ) Helper::halt("bad HDSTATS-PC singular values " + file);
      for ( int k=0; k<nc; ++k ) if ( !parse_value(tok[k+1],&m.singular_values(k)) ||
          !std::isfinite(m.singular_values(k)) ) Helper::halt("bad HDSTATS-PC singular value " + file);
      m.loadings.resize(nf,nc);
      for ( int j=0; j<nf; ++j )
        {
          if ( !std::getline(in,line) ) Helper::halt("truncated HDSTATS-PC loadings " + file);
          tok=split_tab(line);
          if ( (int)tok.size()!=nc+2 || tok[0]!="V" || tok[1]!=m.names[j] )
            Helper::halt("bad HDSTATS-PC loading row " + file);
          for ( int k=0; k<nc; ++k ) if ( !parse_value(tok[k+2],&m.loadings(j,k)) ||
              !std::isfinite(m.loadings(j,k)) ) Helper::halt("bad HDSTATS-PC loading " + file);
        }
      if ( !std::getline(in,line) || line!="END" ) Helper::halt("bad HDSTATS-PC panel end " + file);
      panels.push_back(std::move(panel));
    }
  if ( std::getline(in,line) && !line.empty() ) Helper::halt("trailing HDSTATS-PC model data " + file);
  return panels;
}

static const std::vector<fitted_panel_t> & cached_panels( const std::string & file )
{
  static std::map<std::string,std::vector<fitted_panel_t>> cache;
  static std::mutex mutex;
  const std::string path=Helper::expand(file);
  std::lock_guard<std::mutex> lock(mutex);
  auto found=cache.find(path);
  if ( found!=cache.end() ) return found->second;
  return cache.emplace(path,read_panels(path)).first->second;
}

static void fit_hdstats_panels( param_t & param )
{
  const auto files=param.strvector("data");
  if ( files.empty() ) Helper::halt("--hdstats-pc requires data=<feature files>");
  const auto rows=read_rows(files,hdstats_format);
  const int default_nc=param.has("nc") ? param.requires_int("nc") : 1;
  if ( default_nc<1 ) Helper::halt("nc must be positive");
  std::vector<fitted_panel_t> panels;
  for ( const auto & def : panel_definitions() )
    {
      const int nc=param.has(def.option) ? param.requires_int(def.option) : default_nc;
      if ( nc<1 ) Helper::halt(def.option + " must be positive");
      panels.push_back(fit_panel(def,rows,nc));
      logger << "  HDSTATS-PC " << def.name << ": " << panels.back().n_fit << "/" << rows.size()
             << " people, " << panels.back().model.names.size() << "/" << def.names.size()
             << " varying inputs, nc=" << nc << "\n";
    }
  if ( param.has("corr-matrix") )
    {
      const std::string prefix=param.requires("corr-matrix");
      for ( const auto & p : panels )
        {
          const model_t & m=p.model;
          Eigen::MatrixXd z(p.n_fit,m.names.size());
          int i=0;
          for ( const auto & row : rows )
            {
              bool complete=true;
              for ( const auto & name : p.declared )
                {
                  const auto f=row.values.find(name);
                  if ( f==row.values.end() || !std::isfinite(f->second) ) { complete=false; break; }
                }
              if ( !complete ) continue;
              for ( int j=0; j<(int)m.names.size(); ++j )
                z(i,j)=(row.values.at(m.names[j])-m.means(j))/m.sds(j);
              ++i;
            }
          write_correlation_matrix(prefix+"."+p.name+".tsv",m.names,z,hdstats_format);
        }
    }
  const std::string model_file=param.requires("model");
  write_panels(model_file,panels);
  writer.value("N_PANELS",(int)panels.size());
  writer.value("N_INPUT",(int)rows.size());
  for ( const auto & p : panels )
    {
      const model_t & m=p.model;
      writer.level(p.name,"PANEL");
      writer.value("N",p.n_fit);
      writer.value("N_INPUT",p.n_input);
      writer.value("N_DECLARED",(int)p.declared.size());
      writer.value("N_FEATURES",(int)m.names.size());
      writer.value("N_COMPONENTS",(int)m.loadings.cols());
      for ( int j=0; j<(int)m.names.size(); ++j )
        {
          writer.level(m.names[j],"FTR");
          writer.value("MEAN",m.means(j));writer.value("SD",m.sds(j));
          writer.unlevel("FTR");
        }
      double cumulative=0;
      for ( int k=0; k<m.loadings.cols(); ++k )
        {
          const double ve=m.singular_values(k)*m.singular_values(k)/p.total_ss;
          cumulative+=ve;
          writer.level("PC"+Helper::int2str(k+1),"PC");
          writer.value("W",m.singular_values(k));writer.value("VE",ve);writer.value("CVE",cumulative);
          for ( int j=0; j<(int)m.names.size(); ++j )
            {
              writer.level(m.names[j],"FTR");writer.value("LOADING",m.loadings(j,k));writer.unlevel("FTR");
            }
          writer.unlevel("PC");
        }
      writer.unlevel("PANEL");
    }
}

void hdstats_pc::write_features( const std::string & directory,
                                 const std::string & id,
                                 const feature_map_t & features )
{
  write_features_impl( directory, id, features, hdstats_format );
}

bool hdstats_pc::project( const std::string & file,
                           const feature_map_t & features,
                           std::vector<double> * scores,
                           int * missing )
{
  return project_impl( file, features, scores, missing, hdstats_format );
}

std::vector<hdstats_pc::panel_scores_t> hdstats_pc::project_panels(
    const std::string & file, const feature_map_t & features )
{
  std::ifstream in=LunaIO::open_ifstream(file,std::ios::in);
  std::string magic;
  if ( !std::getline(in,magic) ) Helper::halt("could not read HDSTATS-PC model " + file);
  if ( magic==hdstats_format.model_magic )
    {
      panel_scores_t result;
      result.panel="GENERIC";
      project_impl(file,features,&result.scores,&result.missing,hdstats_format);
      return {result};
    }
  if ( magic!="__HDSTATS_PC_V3" ) Helper::halt("unsupported HDSTATS-PC model " + file);
  std::vector<panel_scores_t> results;
  for ( const auto & p : cached_panels(file) )
    {
      panel_scores_t result;
      result.panel=p.name;
      const model_t & m=p.model;
      Eigen::VectorXd z(m.names.size());
      for ( int j=0; j<(int)m.names.size(); ++j )
        {
          const auto f=features.find(m.names[j]);
          if ( f==features.end() || !std::isfinite(f->second) )
            { ++result.missing; z(j)=0.0; }
          else z(j)=(f->second-m.means(j))/m.sds(j);
        }
      if ( result.missing==0 )
        {
          const Eigen::VectorXd s=m.loadings.transpose()*z;
          result.scores.assign(s.data(),s.data()+s.size());
        }
      results.push_back(std::move(result));
    }
  return results;
}

void hdstats_pc::fit( param_t & param )
{
  fit_hdstats_panels( param );
}

void hypno_pc::write_features( const std::string & directory,
                               const std::string & id,
                               const feature_map_t & features )
{
  write_features_impl( directory, id, features, hypno_format );
}

bool hypno_pc::project( const std::string & file,
                        const feature_map_t & features,
                        std::vector<double> * scores,
                        int * missing )
{
  return project_impl( file, features, scores, missing, hypno_format );
}

void hypno_pc::fit( param_t & param )
{
  fit_impl( param, hypno_format );
}
