//    --------------------------------------------------------------------
//    Native, text-first reference PCA for HDSTATS feature lines.
//    --------------------------------------------------------------------

#include "hdstats-pc.h"
#include "hdstats-model.h"

#include "param.h"
#include "helper/helper.h"
#include "helper/logger.h"
#include "helper/luna_io.h"
#include "helper/json.h"
#include "db/db.h"
#include "stats/Eigen/Dense"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <regex>
#include <map>
#include <mutex>
#include <set>
#include <sstream>

// The generated registry is deliberately compiled into luna.  The editable
// manifest lives in docs/hdstats/hdstats.pc; it is not a run-time dependency.
namespace hdstats_pc_generated {
struct generated_feature_t { const char * panel; const char * stratum; const char * variable; };
struct generated_component_t { const char * id; const char * panel; int pc; int orientation; const char * label; };
}
#include "hdstats-pc-schema.h"

extern logger_t logger;
extern writer_t writer;

namespace {

using json = nlohmann::json;

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

struct component_def_t {
  std::string id;
  std::string label;
  int orientation = 1;
};

struct layout_panel_t {
  panel_def_t definition;
  int n_components = 0;
};

struct layout_t {
  json document;
  json feature_descriptions;
  std::vector<layout_panel_t> panels;
};

std::string panel_option( std::string panel )
{
  std::transform( panel.begin(), panel.end(), panel.begin(), []( unsigned char c ) {
    return c == '_' ? '-' : std::tolower( c );
  } );
  return "nc-" + panel;
}

// Generated from docs/hdstats/hdstats.pc.  The manifest controls both fitting
// feature slots and the V3 projection registry.
const std::vector<panel_def_t> & panel_definitions()
{
  static const std::vector<panel_def_t> panels = [] {
    std::vector<panel_def_t> out;
    std::map<std::string,int> position;
    for ( int i = 0; i < hdstats_pc_generated::n_features; ++i )
      {
        const auto & f = hdstats_pc_generated::features[i];
        const std::string panel = f.panel;
        if ( position.find(panel) == position.end() )
          {
            position[panel] = (int)out.size();
            out.push_back( { panel, panel_option(panel), {} } );
          }
        out[position[panel]].names.push_back( hdstats_pc::key( f.variable, f.stratum ) );
      }
    for ( const auto & panel : out )
      if ( panel.names.empty() || std::set<std::string>(panel.names.begin(),panel.names.end()).size() != panel.names.size() )
        Helper::halt( "invalid generated HDSTATS panel schema: " + panel.name );
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

std::string fmt_display( const double x )
{
  if ( ! std::isfinite( x ) ) return "NA";
  std::ostringstream ss;
  ss << std::fixed << std::setprecision( 3 ) << x;
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
  std::string out = variable;
  std::stringstream fields( stratum );
  std::string field;
  while ( std::getline( fields, field, ';' ) )
    {
      const size_t slash = field.find( '/' );
      if ( slash == std::string::npos || slash == 0 || slash + 1 == field.size() )
        Helper::halt( "bad HDSTATS-PC stratum: " + stratum );
      out += "." + field.substr( 0, slash ) + "_" + field.substr( slash + 1 );
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
  hdstats_model::model_t source;
  std::string error;
  if ( !hdstats_model::read(file,&source,&error) ) Helper::halt(error);
  std::vector<fitted_panel_t> panels;
  for ( const auto & source_panel : source.panels )
    {
      fitted_panel_t panel;
      panel.name = source_panel.id;
      panel.n_input = source_panel.n_input;
      panel.n_fit = source_panel.n_fit;
      panel.total_ss = source_panel.total_ss;
      panel.declared = source_panel.declared_features;
      model_t & m = panel.model;
      const int nf = (int)source_panel.features.size();
      const int nc = (int)source_panel.singular_values.size();
      m.names.resize(nf); m.means.resize(nf); m.sds.resize(nf);
      m.singular_values.resize(nc); m.loadings.resize(nf,nc);
      for ( int k = 0; k < nc; ++k ) m.singular_values(k) = source_panel.singular_values[k];
      for ( int j = 0; j < nf; ++j )
        {
          const auto & feature = source_panel.features[j];
          m.names[j] = feature.id;
          m.means(j) = feature.mean;
          m.sds(j) = feature.sd;
          for ( int k = 0; k < nc; ++k ) m.loadings(j,k) = feature.loadings[k];
        }
      panels.push_back(std::move(panel));
    }
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

static std::map<std::string,int> registry_component_counts()
{
  std::map<std::string,int> counts;
  std::map<std::string,std::set<int>> local_pcs;
  for ( int i = 0; i < hdstats_pc_generated::n_components; ++i )
    {
      const auto & c = hdstats_pc_generated::components[i];
      if ( !local_pcs[c.panel].insert(c.pc).second )
        Helper::halt( "duplicate local PC in HDSTATS-PC registry: " + std::string(c.panel) );
      counts[c.panel] = std::max( counts[c.panel], c.pc );
    }
  for ( const auto & item : counts )
    for ( int pc = 1; pc <= item.second; ++pc )
      if ( local_pcs[item.first].find(pc) == local_pcs[item.first].end() )
        Helper::halt( "non-contiguous local PCs in HDSTATS-PC registry: " + item.first );
  return counts;
}

static const std::map<std::string,int> & component_counts()
{
  static const std::map<std::string,int> counts = registry_component_counts();
  return counts;
}

static void validate_registry()
{
  static const bool valid = [] {
    std::map<std::string,std::set<std::string>> declared;
    for ( const auto & panel : panel_definitions() )
      declared[panel.name] = std::set<std::string>(panel.names.begin(),panel.names.end());
    const auto & counts = component_counts();
    for ( const auto & item : counts )
      if ( declared.find(item.first) == declared.end() )
        Helper::halt( "HDSTATS-PC registry has a component for unknown panel " + item.first );
    std::set<std::string> ids;
    for ( int i = 0; i < hdstats_pc_generated::n_components; ++i )
      {
        const auto & c = hdstats_pc_generated::components[i];
        if ( !ids.insert(c.id).second )
          Helper::halt( "duplicate HDSTATS-PC component ID " + std::string(c.id) );
      }
    return true;
  }();
  (void)valid;
}

static void validate_projection_model( const std::vector<fitted_panel_t> & fitted )
{
  validate_registry();
  const auto & definitions = panel_definitions();
  const auto & counts = component_counts();
  if ( fitted.size() != definitions.size() )
    Helper::halt( "HDSTATS-PC model does not match the compiled V3 panel registry" );
  std::map<std::string,const panel_def_t *> definitions_by_name;
  for ( const auto & definition : definitions ) definitions_by_name[definition.name] = &definition;
  std::set<std::string> names;
  for ( const auto & panel : fitted )
    {
      const auto def = definitions_by_name.find(panel.name);
      const auto expected = counts.find(panel.name);
      if ( def == definitions_by_name.end() || expected == counts.end() || !names.insert(panel.name).second ||
           panel.declared != def->second->names || panel.model.loadings.cols() != expected->second )
        Helper::halt( "HDSTATS-PC model panel does not match the current registry: " + panel.name );
    }
}

static const std::set<std::string> & known_feature_universe()
{
  // The compiled registry is the current stable public feature universe.  A
  // JSON layout may redistribute these features freely across panels but may
  // not invent an unimplemented HDSTATS feature key.
  static const std::set<std::string> known = [] {
    std::set<std::string> values;
    for ( int i = 0; i < hdstats_pc_generated::n_features; ++i )
      values.insert( hdstats_pc::key( hdstats_pc_generated::features[i].variable,
                                      hdstats_pc_generated::features[i].stratum ) );
    return values;
  }();
  return known;
}

static layout_t read_layout_spec( const std::string & file )
{
  try
    {
      std::ifstream in = LunaIO::open_ifstream(file,std::ios::in);
      if ( ! in ) Helper::halt("could not open HDSTATS-PC JSON specification " + file);
      layout_t layout;
      layout.document = json::parse(in);
      const json & document = layout.document;
      if ( ! document.is_object() || document.value("format","") != "luna-hdstats-pc-spec" ||
           document.value("version",0) != 1 ||
           document.value("feature_universe","") != "HDSTATS-PC-FEATURES-V1" ||
           ! document.contains("feature_descriptions") || ! document.at("feature_descriptions").is_array() ||
           ! document.contains("panels") || ! document.at("panels").is_array() || document.at("panels").empty() )
        Helper::halt("invalid HDSTATS-PC JSON specification " + file + ": require format/version/feature_universe/feature_descriptions/panels");
      const std::string component_id_style = document.value("component_id_style",std::string("legacy"));
      if ( component_id_style != "legacy" && component_id_style != "r_safe" )
        Helper::halt("invalid HDSTATS-PC component_id_style in " + file + ": expected legacy or r_safe");
      layout.feature_descriptions = document.at("feature_descriptions");

      std::map<std::string,std::pair<std::string,std::string>> feature_metadata;
      for ( const auto & item : layout.feature_descriptions )
        {
          if ( ! item.is_object() || ! item.contains("id") || ! item.at("id").is_string() ||
               ! item.contains("label") || ! item.at("label").is_string() ||
               ! item.contains("description") || ! item.at("description").is_string() )
            Helper::halt("invalid HDSTATS-PC feature description in " + file);
          const std::string id = item.at("id").get<std::string>();
          if ( ! known_feature_universe().count(id) || item.at("label").get<std::string>().empty() ||
               item.at("description").get<std::string>().empty() ||
               ! feature_metadata.emplace(id,std::make_pair(item.at("label").get<std::string>(),item.at("description").get<std::string>())).second )
            Helper::halt("unknown, empty, or duplicate HDSTATS-PC feature description '" + id + "'");
        }

      std::set<std::string> panel_ids, used_features;
      for ( const auto & source_panel : document.at("panels") )
        {
          if ( ! source_panel.is_object() || ! source_panel.contains("id") || ! source_panel.at("id").is_string() ||
               source_panel.at("id").get<std::string>().empty() || ! source_panel.contains("features") ||
               ! source_panel.at("features").is_array() || source_panel.at("features").empty() ||
               ! source_panel.contains("n_components") || ! source_panel.at("n_components").is_number_integer() ||
               source_panel.at("n_components").get<int>() < 1 )
            Helper::halt("invalid HDSTATS-PC JSON specification panel in " + file);
          layout_panel_t panel;
          panel.definition.name = source_panel.at("id").get<std::string>();
          if ( component_id_style == "r_safe" &&
               !std::regex_match(panel.definition.name,std::regex("[A-Za-z][A-Za-z0-9_.]*")) )
            Helper::halt("HDSTATS-PC r_safe panel IDs must be R-safe name prefixes: " + panel.definition.name);
          if ( ! panel_ids.insert(panel.definition.name).second )
            Helper::halt("duplicate HDSTATS-PC JSON panel " + panel.definition.name);
          std::set<std::string> panel_features;
          for ( const auto & source_feature : source_panel.at("features") )
            {
              if ( ! source_feature.is_string() || source_feature.get<std::string>().empty() )
                Helper::halt("invalid HDSTATS-PC JSON feature in panel " + panel.definition.name);
              const std::string feature = source_feature.get<std::string>();
              if ( ! known_feature_universe().count(feature) )
                Helper::halt("unknown HDSTATS-PC JSON feature '" + feature + "' in panel " + panel.definition.name);
              if ( ! panel_features.insert(feature).second )
                Helper::halt("duplicate HDSTATS-PC JSON feature '" + feature + "' in panel " + panel.definition.name);
              if ( ! feature_metadata.count(feature) )
                Helper::halt("missing HDSTATS-PC JSON feature description for '" + feature + "'");
              used_features.insert(feature);
              panel.definition.names.push_back(feature);
            }
          panel.n_components = source_panel.at("n_components").get<int>();
          layout.panels.push_back(std::move(panel));
        }
      for ( const auto & item : feature_metadata )
        if ( ! used_features.count(item.first) )
          Helper::halt("unused HDSTATS-PC JSON feature description '" + item.first + "'");
      return layout;
    }
  catch ( const std::exception & e )
    {
      Helper::halt("invalid HDSTATS-PC JSON specification " + file + ": " + e.what());
    }
  return layout_t(); // unreachable after Helper::halt(), keeps all compilers happy
}

static void write_json_panels( const std::string & file, const layout_t & layout,
                               const std::vector<fitted_panel_t> & panels )
{
  json document;
  document["format"] = "luna-hdstats-pc-model";
  document["version"] = 1;
  const size_t slash = file.find_last_of("/\\");
  const std::string model_name = file.substr(slash == std::string::npos ? 0 : slash + 1);
  document["model_id"] = model_name;
  document["spec"] = layout.document;
  document["panels"] = json::array();
  for ( const auto & panel : panels )
    {
      const model_t & model = panel.model;
      json item;
      item["id"] = panel.name;
      item["n_input"] = panel.n_input;
      item["n_fit"] = panel.n_fit;
      item["total_ss"] = panel.total_ss;
      item["declared_features"] = panel.declared;
      item["singular_values"] = json::array();
      for ( int k = 0; k < model.singular_values.size(); ++k ) item["singular_values"].push_back(model.singular_values(k));
      item["features"] = json::array();
      for ( int j = 0; j < (int)model.names.size(); ++j )
        {
          json feature;
          feature["id"] = model.names[j];
          feature["mean"] = model.means(j);
          feature["sd"] = model.sds(j);
          feature["loadings"] = json::array();
          for ( int k = 0; k < model.loadings.cols(); ++k ) feature["loadings"].push_back(model.loadings(j,k));
          item["features"].push_back(std::move(feature));
        }
      document["panels"].push_back(std::move(item));
    }
  std::ofstream out = LunaIO::open_ofstream(file,std::ios::out);
  if ( ! out ) Helper::halt("could not write HDSTATS-PC JSON model " + file);
  out << document.dump(2) << "\n";
  if ( ! out ) Helper::halt("error writing HDSTATS-PC JSON model " + file);
}

static void fit_hdstats_panels( param_t & param )
{
  validate_registry();
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
      const int expected = component_counts().at(def.name);
      if ( nc != expected )
        Helper::halt( "HDSTATS-PC " + def.name + " requires " + Helper::int2str(expected) +
                      " components for the current registry (got " + Helper::int2str(nc) + ")" );
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

static void fit_json_hdstats_panels( param_t & param )
{
  // A JSON layout owns component counts.  Mixing it with the legacy nc knobs
  // would make a model's declared layout ambiguous.
  if ( param.has("nc") ) Helper::halt("HDSTATS-PC spec= cannot be combined with nc=");
  for ( const auto & definition : panel_definitions() )
    if ( param.has(definition.option) )
      Helper::halt("HDSTATS-PC spec= cannot be combined with " + definition.option + "=");

  const layout_t layout = read_layout_spec( Helper::expand( param.requires("spec") ) );
  const auto files = param.strvector("data");
  if ( files.empty() ) Helper::halt("--hdstats-pc requires data=<feature files>");
  const auto rows = read_rows(files,hdstats_format);
  std::vector<fitted_panel_t> panels;
  for ( const auto & definition : layout.panels )
    {
      panels.push_back( fit_panel( definition.definition, rows, definition.n_components ) );
      logger << "  HDSTATS-PC " << definition.definition.name << ": " << panels.back().n_fit << "/" << rows.size()
             << " people, " << panels.back().model.names.size() << "/" << definition.definition.names.size()
             << " varying inputs, nc=" << definition.n_components << "\n";
    }

  if ( param.has("corr-matrix") )
    {
      const std::string prefix = param.requires("corr-matrix");
      for ( const auto & panel : panels )
        {
          const model_t & model = panel.model;
          Eigen::MatrixXd z(panel.n_fit,model.names.size());
          int i = 0;
          for ( const auto & row : rows )
            {
              bool complete = true;
              for ( const auto & name : panel.declared )
                {
                  const auto feature = row.values.find(name);
                  if ( feature == row.values.end() || ! std::isfinite(feature->second) ) { complete = false; break; }
                }
              if ( ! complete ) continue;
              for ( int j = 0; j < (int)model.names.size(); ++j )
                z(i,j) = ( row.values.at(model.names[j]) - model.means(j) ) / model.sds(j);
              ++i;
            }
          write_correlation_matrix(prefix + "." + panel.name + ".tsv",model.names,z,hdstats_format);
        }
    }

  write_json_panels( param.requires("model"), layout, panels );
  writer.value("N_PANELS",(int)panels.size());
  writer.value("N_INPUT",(int)rows.size());
  for ( int p = 0; p < (int)panels.size(); ++p )
    {
      const fitted_panel_t & panel = panels[p];
      const model_t & model = panel.model;
      writer.level(panel.name,"PANEL");
      writer.value("N",panel.n_fit);
      writer.value("N_INPUT",panel.n_input);
      writer.value("N_DECLARED",(int)panel.declared.size());
      writer.value("N_FEATURES",(int)model.names.size());
      writer.value("N_COMPONENTS",(int)model.loadings.cols());
      for ( int j = 0; j < (int)model.names.size(); ++j )
        {
          writer.level(model.names[j],"FTR");
          writer.value("MEAN",model.means(j)); writer.value("SD",model.sds(j));
          writer.unlevel("FTR");
        }
      double cumulative = 0.0;
      for ( int k = 0; k < model.loadings.cols(); ++k )
        {
          const double ve = model.singular_values(k) * model.singular_values(k) / panel.total_ss;
          cumulative += ve;
          writer.level("PC" + Helper::int2str(k+1),"PC");
          writer.value("W",model.singular_values(k)); writer.value("VE",ve); writer.value("CVE",cumulative);
          for ( int j = 0; j < (int)model.names.size(); ++j )
            {
              writer.level(model.names[j],"FTR"); writer.value("LOADING",model.loadings(j,k)); writer.unlevel("FTR");
            }
          writer.unlevel("PC");
        }
      writer.unlevel("PANEL");
    }
}

static bool is_json_hdstats_model( const std::string & file )
{
  std::ifstream in = LunaIO::open_ifstream(file,std::ios::in);
  if ( ! in ) Helper::halt("could not open HDSTATS-PC model " + file);
  in >> std::ws;
  return in.peek() == '{';
}

static hdstats_model::model_t read_json_hdstats_model( const std::string & file )
{
  hdstats_model::model_t model;
  std::string error;
  if ( ! hdstats_model::read(file,&model,&error) ) Helper::halt(error);
  if ( model.specification_json.empty() ) Helper::halt("expected self-contained HDSTATS-PC JSON model " + file);
  return model;
}

static std::vector<hdstats_pc::component_score_t> project_json_components(
    const std::string & file, const hdstats_pc::feature_map_t & features, const double min_loading )
{
  const hdstats_model::model_t source = read_json_hdstats_model(file);
  std::vector<hdstats_pc::component_score_t> results;
  const double nan = std::numeric_limits<double>::quiet_NaN();
  for ( const auto & definition : source.components )
    {
      const hdstats_model::panel_t * panel = source.find_panel(definition.panel_id);
      if ( ! panel || definition.panel_pc < 1 || definition.panel_pc > (int)panel->singular_values.size() )
        Helper::halt("invalid embedded component layout in HDSTATS-PC JSON model " + file);
      hdstats_pc::component_score_t result;
      result.id = definition.id;
      result.panel = definition.panel_id;
      result.panel_pc = definition.panel_pc;
      result.orientation = definition.orientation;
      result.label = definition.label;
      result.compatible = true;
      result.missing = 0;
      result.score = nan;
      double score = 0.0;
      for ( const auto & feature : panel->features )
        {
          const auto raw = features.find(feature.id);
          const bool present = raw != features.end() && std::isfinite(raw->second);
          if ( ! present ) { ++result.missing; result.compatible = false; }
          const double z = present ? ( raw->second - feature.mean ) / feature.sd : 0.0;
          const double loading = definition.orientation * feature.loadings[definition.panel_pc-1];
          score += loading * z;
          if ( std::fabs(loading) >= min_loading )
            result.contributions.push_back( { feature.id, present ? loading * z : nan,
                                               present ? raw->second : nan } );
        }
      if ( result.compatible ) result.score = score;
      results.push_back(std::move(result));
    }
  return results;
}

static void describe_json_model( const std::string & file, std::ostream & out, const double min_loading,
                                const std::string & descriptor )
{
  hdstats_model::model_t source = read_json_hdstats_model(file);
  if ( ! descriptor.empty() )
    {
      std::string error;
      if ( ! hdstats_model::read_descriptor(descriptor,&source,&error) ) Helper::halt(error);
    }
  out << "PC\tPANEL\tMODEL_PC\tLABEL\tORIENTATION\tDESCRIPTION\tFTR_LOADINGS\n";
  for ( const auto & component : source.components )
    {
      const hdstats_model::panel_t * panel = source.find_panel(component.panel_id);
      if ( ! panel || component.panel_pc < 1 || component.panel_pc > (int)panel->singular_values.size() )
        Helper::halt("invalid embedded component layout in HDSTATS-PC JSON model " + file);
      std::vector<std::string> loadings;
      for ( const auto & feature : panel->features )
        {
          const double loading = component.orientation * feature.loadings[component.panel_pc-1];
          if ( std::fabs(loading) >= min_loading ) loadings.push_back(feature.id + "=" + fmt_display(loading));
        }
      const auto panel_description = source.find_description("PANEL",component.panel_id);
      out << component.id << '\t' << component.panel_id << '\t' << component.panel_pc << '\t'
          << component.label << '\t' << component.orientation << '\t'
          << component.description << (panel_description ? " | " + panel_description->text : "") << '\t';
      for ( int i = 0; i < (int)loadings.size(); ++i )
        { if ( i ) out << "; "; out << loadings[i]; }
      out << '\n';
    }
  if ( ! source.descriptions.empty() )
    {
      out << "\n# DESCRIPTIONS\nENTITY_TYPE\tENTITY_ID\tTEXT\n";
      for ( const auto & item : source.descriptions )
        out << item.entity_type << '\t' << item.entity_id << '\t' << item.text << '\n';
    }
  if ( ! source.external_associations.empty() )
    {
      out << "\n# EXTERNAL ASSOCIATIONS\nENTITY_TYPE\tENTITY_ID\tFACTOR\tSCORE\n";
      for ( const auto & item : source.external_associations )
        out << item.entity_type << '\t' << item.entity_id << '\t' << item.factor << '\t' << fmt_display(item.score) << '\n';
    }
  auto print_correlations = [&]( const char * title, const std::vector<hdstats_model::correlation_t> & items ) {
    if ( items.empty() ) return;
    out << "\n# " << title << "\nFIRST\tSECOND\tCORRELATION\n";
    for ( const auto & item : items ) out << item.first << '\t' << item.second << '\t' << fmt_display(item.value) << '\n';
  };
  print_correlations("PC CORRELATIONS",source.pc_correlations);
  print_correlations("FEATURE CORRELATIONS",source.feature_correlations);
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
  const bool json_model = is_json_hdstats_model(file);
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
  if ( magic!="__HDSTATS_PC_V3" && ! json_model ) Helper::halt("unsupported HDSTATS-PC model " + file);
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

std::vector<hdstats_pc::component_score_t> hdstats_pc::project_components(
    const std::string & file, const feature_map_t & features, const double min_loading )
{
  if ( ! std::isfinite(min_loading) || min_loading < 0.0 )
    Helper::halt("pc-min-loading must be a finite non-negative value");
  if ( is_json_hdstats_model(file) ) return project_json_components(file,features,min_loading);
  std::ifstream in=LunaIO::open_ifstream(file,std::ios::in);
  std::string magic;
  if ( !std::getline(in,magic) ) Helper::halt("could not read HDSTATS-PC model " + file);
  if ( magic == hdstats_format.model_magic )
    Helper::halt("HDSTATS-PC global component projection requires a V3 model: " + file);
  if ( magic != "__HDSTATS_PC_V3" ) Helper::halt("unsupported HDSTATS-PC model " + file);

  const auto & fitted = cached_panels(file);
  validate_projection_model(fitted);
  std::map<std::string,const fitted_panel_t *> panels;
  for ( const auto & panel : fitted ) panels[panel.name] = &panel;

  std::vector<component_score_t> results;
  const double nan = std::numeric_limits<double>::quiet_NaN();
  for ( int i = 0; i < hdstats_pc_generated::n_components; ++i )
    {
      const auto & definition = hdstats_pc_generated::components[i];
      const fitted_panel_t & panel = *panels.at(definition.panel);
      const model_t & model = panel.model;
      component_score_t result;
      result.id=definition.id;
      result.panel=definition.panel;
      result.panel_pc=definition.pc;
      result.orientation=definition.orientation;
      result.label=definition.label;
      result.compatible=true;
      result.missing=0;
      result.score=nan;

      Eigen::VectorXd z(model.names.size());
      std::map<std::string,int> feature_index;
      for ( int j=0; j<(int)model.names.size(); ++j )
        {
          feature_index[model.names[j]]=j;
          const auto feature=features.find(model.names[j]);
          if ( feature==features.end() || !std::isfinite(feature->second) )
            { ++result.missing; result.compatible=false; z(j)=0.0; }
          else z(j)=(feature->second-model.means(j))/model.sds(j);
        }
      if ( result.compatible )
        result.score=definition.orientation * model.loadings.col(definition.pc-1).dot(z);

      for ( int j = 0; j < (int)model.names.size(); ++j )
        {
          const double loading = definition.orientation * model.loadings(j,definition.pc-1);
          if ( std::fabs(loading) < min_loading ) continue;
          const auto raw = features.find(model.names[j]);
          result.contributions.push_back( { model.names[j], result.compatible
            ? loading * z(j) : nan,
            raw == features.end() || !std::isfinite(raw->second) ? nan : raw->second } );
        }
      results.push_back(std::move(result));
    }
  return results;
}

void hdstats_pc::describe( const std::string & file, std::ostream & out, const double min_loading,
                           const std::string & descriptor )
{
  if ( ! std::isfinite(min_loading) || min_loading < 0.0 )
    Helper::halt("min-loading must be a finite non-negative value");
  if ( is_json_hdstats_model(file) )
    { describe_json_model(file,out,min_loading,descriptor); return; }
  if ( ! descriptor.empty() ) Helper::halt("HDSTATS descriptors are currently supported for self-contained JSON V1 models only");
  validate_registry();
  const auto & fitted = cached_panels(file);
  validate_projection_model(fitted);
  std::map<std::string,const fitted_panel_t *> panels;
  for ( const auto & panel : fitted ) panels[panel.name] = &panel;
  out << "PC\tDOMAIN\tMODEL_PC\tLABEL\tORIENTATION\tFTR_LOADINGS\n";
  for ( int i = 0; i < hdstats_pc_generated::n_components; ++i )
    {
      const auto & component = hdstats_pc_generated::components[i];
      const model_t & model = panels.at(component.panel)->model;
      std::vector<std::string> exemplar_features;
      for ( int j = 0; j < (int)model.names.size(); ++j )
        {
          const double loading = component.orientation * model.loadings(j,component.pc-1);
          if ( std::fabs(loading) >= min_loading )
            exemplar_features.push_back( model.names[j] + "=" + fmt_display( loading ) );
        }
      out << component.id << '\t'
          << component.panel << '\t'
          << component.pc << '\t'
          << component.label << '\t'
          << component.orientation << '\t';
      for ( int j = 0; j < (int)exemplar_features.size(); ++j )
        {
          if ( j ) out << "; ";
          out << exemplar_features[j];
        }
      out << '\n';
    }
}

void hdstats_pc::fit( param_t & param )
{
  if ( param.has("spec") ) fit_json_hdstats_panels( param );
  else fit_hdstats_panels( param );
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
