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
// JSON specification lives in docs/hdstats/hdstats-pc.spec.json; it is not a
// run-time dependency.
namespace hdstats_pc_generated {
struct generated_panel_t { const char * id; int n_components; };
struct generated_feature_t { const char * panel; const char * weight_group; const char * id; };
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
  Eigen::VectorXd weights;
  Eigen::MatrixXd loadings;
  Eigen::VectorXd singular_values;
};

struct panel_def_t {
  std::string name;
  std::string option;
  std::vector<std::string> names;
  std::vector<double> weights;
  std::vector<std::string> groups; // parallel to names, from the fitting spec
};

struct fitted_panel_t {
  std::string name;
  std::vector<std::string> declared;
  model_t model;
  // One-based axes selected before rotation. For rotated models these index
  // the ROWS of rotation, not the final output columns.
  std::vector<int> source_components;
  // post-rotation selection only: one-based indices into the rotated basis
  // (rotation's ROWS, i.e. the PANEL,PC numbering of the full nc rotated
  // fit) that were kept as this panel's final output COLUMNS. Empty unless
  // pc-select-stage=post-rotation was used.
  std::vector<int> selected_components;
  int n_input = 0;
  int n_fit = 0;
  double total_ss = 0.0;
  Eigen::MatrixXd rotation; // empty for ordinary PCA; final loadings = V * rotation
  Eigen::VectorXd score_sds;
  Eigen::MatrixXd score_correlations;
  std::vector<std::string> group_ids;
  Eigen::MatrixXd group_mass; // rows = groups, columns = final components
  double rotation_before = 0.0, rotation_after = 0.0;
  int rotation_sweeps = 0;
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

using component_selection_t = std::map<std::string,std::vector<int>>;

// Standardized inputs in a component with p measures receive
// 1/sqrt(G*p), where G is the panel's component count.  Each component then
// contributes 1/G to the trace of the weighted input covariance matrix.
static double equal_component_weight( const int n_groups, const int group_size )
{
  if ( n_groups < 1 || group_size < 1 ) Helper::halt("invalid HDSTATS-PC weight group");
  return 1.0 / std::sqrt( (double)n_groups * group_size );
}

std::string panel_option( std::string panel )
{
  std::transform( panel.begin(), panel.end(), panel.begin(), []( unsigned char c ) {
    return c == '_' ? '-' : std::tolower( c );
  } );
  return "nc-" + panel;
}

// Generated from the v2 JSON specification.  The registry controls the
// feature producer contract, while fitted JSON models retain their own layout.
const std::vector<panel_def_t> & panel_definitions()
{
  static const std::vector<panel_def_t> panels = [] {
    std::vector<panel_def_t> out;
    std::map<std::string,int> position;
    std::map<std::string,int> group_sizes;
    std::map<std::string,std::set<std::string>> panel_groups;
    for ( int i = 0; i < hdstats_pc_generated::n_features; ++i )
      {
        const auto & f = hdstats_pc_generated::features[i];
        ++group_sizes[std::string(f.panel) + "\t" + f.weight_group];
        panel_groups[f.panel].insert(f.weight_group);
      }
    for ( int i = 0; i < hdstats_pc_generated::n_features; ++i )
      {
        const auto & f = hdstats_pc_generated::features[i];
        const std::string panel = f.panel;
        if ( position.find(panel) == position.end() )
          {
            position[panel] = (int)out.size();
            out.push_back( { panel, panel_option(panel), {} } );
          }
        out[position[panel]].names.push_back( f.id );
        out[position[panel]].weights.push_back(
          equal_component_weight( (int)panel_groups[panel].size(),
                                  group_sizes[panel + "\t" + f.weight_group] ) );
      }
    for ( const auto & panel : out )
      if ( panel.names.empty() || panel.names.size() != panel.weights.size() ||
           std::set<std::string>(panel.names.begin(),panel.names.end()).size() != panel.names.size() )
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
  m.names.resize( nf ); m.means.resize( nf ); m.sds.resize( nf ); m.weights = Eigen::VectorXd::Ones( nf );
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
      else z(j) = m.weights(j) * ( f->second - m.means(j) ) / m.sds(j);
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
  m.weights = Eigen::VectorXd::Ones( keep.size() );
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
  m.weights.resize(keep.size());
  for ( int j = 0; j < (int)keep.size(); ++j )
    {
      const int k = keep[j];
      m.names[j] = def.names[k]; m.means(j) = means(k); m.sds(j) = sds(k);
      m.weights(j) = def.weights[k];
      z.col(j) = m.weights(j) * (x.col(k).array() - means(k)) / sds(k);
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
      // The FINAL output width is the loadings length, not singular_values'
      // length: for a JSON V3 post-rotation-selected panel, singular_values
      // is the (wider) rotated basis, while loadings/score_sds are the
      // selected K columns actually used for projection.
      const int nc = nf ? (int)source_panel.features.front().loadings.size() : 0;
      m.names.resize(nf); m.means.resize(nf); m.sds.resize(nf); m.weights = Eigen::VectorXd::Ones(nf);
      m.singular_values.resize(nc); m.loadings.resize(nf,nc);
      for ( int k = 0; k < nc && k < (int)source_panel.singular_values.size(); ++k )
        m.singular_values(k) = source_panel.singular_values[k];
      for ( int j = 0; j < nf; ++j )
        {
          const auto & feature = source_panel.features[j];
          m.names[j] = feature.id;
          m.means(j) = feature.mean;
          m.sds(j) = feature.sd;
          m.weights(j) = feature.weight;
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

static void validate_projection_model( const std::vector<fitted_panel_t> & fitted )
{
  std::set<std::string> names;
  for ( const auto & panel : fitted )
    if ( panel.name.empty() || !names.insert(panel.name).second || panel.model.names.empty() ||
         panel.model.weights.size() != panel.model.names.size() || panel.model.loadings.cols() < 1 )
      Helper::halt( "invalid HDSTATS-PC panel model: " + panel.name );
}

static const std::set<std::string> & known_feature_universe()
{
  // The compiled registry is generated directly from the v2 JSON source.
  static const std::set<std::string> known = [] {
    std::set<std::string> values;
    for ( int i = 0; i < hdstats_pc_generated::n_features; ++i )
      values.insert( hdstats_pc_generated::features[i].id );
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
           document.value("version",0) != 2 ||
           document.value("feature_universe","") != "HDSTATS-PC-FEATURES-V2" ||
           ! document.contains("feature_descriptions") || ! document.at("feature_descriptions").is_array() ||
           ! document.contains("panels") || ! document.at("panels").is_array() || document.at("panels").empty() ||
           ! document.contains("derived_features") || ! document.at("derived_features").is_array() )
        Helper::halt("invalid HDSTATS-PC JSON specification " + file + ": require v2 format, feature descriptions, panels, and derived features");
      if ( ! document.contains("weighting") || ! document.at("weighting").is_object() ||
           document.at("weighting").value("method","") != "equal_component_variance" )
        Helper::halt("HDSTATS-PC JSON specification " + file + ": require equal_component_variance weighting");
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

      std::map<std::string,std::vector<std::pair<std::string,std::string>>> compiled;
      for ( int i = 0; i < hdstats_pc_generated::n_features; ++i )
        compiled[hdstats_pc_generated::features[i].panel].push_back(
          { hdstats_pc_generated::features[i].weight_group, hdstats_pc_generated::features[i].id } );
      std::map<std::string,int> compiled_counts;
      for ( int i = 0; i < hdstats_pc_generated::n_panels; ++i )
        compiled_counts[hdstats_pc_generated::panels[i].id] = hdstats_pc_generated::panels[i].n_components;
      std::set<std::string> panel_ids, used_features;
      for ( const auto & source_panel : document.at("panels") )
        {
          if ( ! source_panel.is_object() || ! source_panel.contains("id") || ! source_panel.at("id").is_string() ||
               source_panel.at("id").get<std::string>().empty() || ! source_panel.contains("weight_groups") ||
               ! source_panel.at("weight_groups").is_array() || source_panel.at("weight_groups").empty() ||
               ! source_panel.contains("n_components") || ! source_panel.at("n_components").is_number_integer() ||
               source_panel.at("n_components").get<int>() < 1 )
            Helper::halt("invalid HDSTATS-PC JSON specification panel in " + file);
          layout_panel_t panel;
          panel.definition.name = source_panel.at("id").get<std::string>();
          if ( ! panel_ids.insert(panel.definition.name).second )
            Helper::halt("duplicate HDSTATS-PC JSON panel " + panel.definition.name);
          std::set<std::string> panel_features, group_ids;
          std::vector<std::pair<std::string,std::string>> observed;
          const int n_groups = (int)source_panel.at("weight_groups").size();
          for ( const auto & group : source_panel.at("weight_groups") )
            {
              if ( ! group.is_object() || ! group.contains("id") || !group.at("id").is_string() ||
                   group.at("id").get<std::string>().empty() || !group.contains("features") ||
                   !group.at("features").is_array() || group.at("features").empty() )
                Helper::halt("invalid HDSTATS-PC weight group in panel " + panel.definition.name);
              const std::string group_id = group.at("id").get<std::string>();
              if ( !group_ids.insert(group_id).second ) Helper::halt("duplicate HDSTATS-PC weight group " + group_id);
              const double weight = equal_component_weight( n_groups, (int)group.at("features").size() );
              for ( const auto & source_feature : group.at("features") )
                {
                  if ( !source_feature.is_string() || source_feature.get<std::string>().empty() )
                    Helper::halt("invalid HDSTATS-PC JSON feature in panel " + panel.definition.name);
                  const std::string feature = source_feature.get<std::string>();
                  if ( !known_feature_universe().count(feature) || !feature_metadata.count(feature) ||
                       !panel_features.insert(feature).second )
                    Helper::halt("unknown, undescribed, or duplicate HDSTATS-PC feature '" + feature + "' in panel " + panel.definition.name);
                  used_features.insert(feature);
                  observed.push_back({group_id,feature});
                  panel.definition.names.push_back(feature);
                  panel.definition.weights.push_back(weight);
                  panel.definition.groups.push_back(group_id);
                }
            }
          if ( compiled.find(panel.definition.name) == compiled.end() || observed != compiled.at(panel.definition.name) )
            Helper::halt("HDSTATS-PC JSON panel does not match the generated v2 registry: " + panel.definition.name);
          panel.n_components = source_panel.at("n_components").get<int>();
          if ( panel.n_components != compiled_counts.at(panel.definition.name) )
            Helper::halt("HDSTATS-PC JSON panel count does not match the generated v2 registry: " + panel.definition.name);
          layout.panels.push_back(std::move(panel));
        }
      if ( panel_ids.size() != compiled.size() || compiled.size() != compiled_counts.size() )
        Helper::halt("HDSTATS-PC JSON panels do not match the generated v2 registry");
      if ( (int)document.at("derived_features").size() != 4 )
        Helper::halt("HDSTATS-PC v2 requires its four declared NREM derived features");
      const std::map<std::string,std::pair<std::vector<std::string>,std::string>> expected_derived = {
        { "NREM_ILR_N1_VS_N2N3", {{"HD_MINS.SS_N1","HD_MINS.SS_N2","HD_MINS.SS_N3"}, "sqrt(2/3) * log(N1 / sqrt(N2 * N3))"} },
        { "NREM_ILR_N2_VS_N3", {{"HD_MINS.SS_N1","HD_MINS.SS_N2","HD_MINS.SS_N3"}, "sqrt(1/2) * log(N2 / N3)"} },
        { "NREM_COM_N1_MINUS_N2", {{"COM.NSS_5.VAR_N1","COM.NSS_5.VAR_N2"}, "COM_N1 - COM_N2"} },
        { "NREM_COM_N2_MINUS_N3", {{"COM.NSS_5.VAR_N2","COM.NSS_5.VAR_N3"}, "COM_N2 - COM_N3"} }
      };
      std::set<std::string> derived_ids;
      for ( const auto & item : document.at("derived_features") )
        {
          if ( !item.is_object() || !item.contains("id") || !item.at("id").is_string() ||
               !item.contains("inputs") || !item.at("inputs").is_array() || !item.contains("formula") || !item.at("formula").is_string() )
            Helper::halt("invalid HDSTATS-PC v2 derived feature");
          const std::string id = item.at("id").get<std::string>();
          const auto expected = expected_derived.find(id);
          if ( expected == expected_derived.end() || !derived_ids.insert(id).second ||
               item.at("formula").get<std::string>() != expected->second.second )
            Helper::halt("unsupported HDSTATS-PC v2 derived feature " + id);
          std::vector<std::string> inputs;
          for ( const auto & input : item.at("inputs") )
            {
              if ( !input.is_string() ) Helper::halt("invalid input for HDSTATS-PC derived feature " + id);
              inputs.push_back(input.get<std::string>());
            }
          if ( inputs != expected->second.first ) Helper::halt("unexpected inputs for HDSTATS-PC derived feature " + id);
        }
      if ( derived_ids.size() != expected_derived.size() ) Helper::halt("missing HDSTATS-PC v2 derived feature");
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

static void derive_nrem_features( hdstats_pc::feature_map_t * values )
{
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const auto value = [&]( const char * id ) {
    const auto found = values->find(id);
    return found == values->end() ? nan : found->second;
  };
  const double n1 = value("HD_MINS.SS_N1");
  const double n2 = value("HD_MINS.SS_N2");
  const double n3 = value("HD_MINS.SS_N3");
  // Posterior-minute masses are non-negative.  This fixed positive floor
  // makes the documented ILR transform finite at an empty stage boundary.
  const double floor = 1e-6;
  if ( std::isfinite(n1) && std::isfinite(n2) && std::isfinite(n3) )
    {
      const double a = std::max(n1,floor), b = std::max(n2,floor), c = std::max(n3,floor);
      (*values)["NREM_ILR_N1_VS_N2N3"] = std::sqrt(2.0/3.0) * std::log(a / std::sqrt(b*c));
      (*values)["NREM_ILR_N2_VS_N3"] = std::sqrt(1.0/2.0) * std::log(b / c);
    }
  else
    (*values)["NREM_ILR_N1_VS_N2N3"] = (*values)["NREM_ILR_N2_VS_N3"] = nan;
  const double com_n1 = value("COM.NSS_5.VAR_N1");
  const double com_n2 = value("COM.NSS_5.VAR_N2");
  const double com_n3 = value("COM.NSS_5.VAR_N3");
  (*values)["NREM_COM_N1_MINUS_N2"] = std::isfinite(com_n1) && std::isfinite(com_n2) ? com_n1 - com_n2 : nan;
  (*values)["NREM_COM_N2_MINUS_N3"] = std::isfinite(com_n2) && std::isfinite(com_n3) ? com_n2 - com_n3 : nan;
}

static std::vector<row_t> with_nrem_features( std::vector<row_t> rows )
{
  for ( auto & row : rows ) derive_nrem_features(&row.values);
  return rows;
}

// pc-select= is deliberately a fit-time input rather than a projection-time
// dependency.  The emitted JSON model contains the selected columns, so a
// deployed pc-proj= model is entirely self-contained.
static component_selection_t read_component_selection( const std::string & file,
                                                        const layout_t & layout )
{
  try
    {
      std::ifstream in = LunaIO::open_ifstream(file,std::ios::in);
      if ( ! in ) Helper::halt("could not open HDSTATS-PC selection " + file);
      const json document = json::parse(in);
      if ( ! document.is_object() || document.value("format","") != "luna-hdstats-pc-selection" ||
           document.value("version",0) != 1 || ! document.contains("panels") ||
           ! document.at("panels").is_object() )
        Helper::halt("invalid HDSTATS-PC selection " + file + ": require V1 format and panels object");

      std::set<std::string> known_panels;
      for ( const auto & panel : layout.panels ) known_panels.insert(panel.definition.name);
      component_selection_t selected;
      for ( auto it = document.at("panels").begin(); it != document.at("panels").end(); ++it )
        {
          const std::string panel = it.key();
          if ( ! known_panels.count(panel) || ! it.value().is_array() || it.value().empty() )
            Helper::halt("invalid HDSTATS-PC selection panel " + panel);
          int previous = 0;
          std::vector<int> axes;
          for ( const auto & value : it.value() )
            {
              if ( ! value.is_number_integer() || value.get<int>() <= previous )
                Helper::halt("HDSTATS-PC selection axes must be positive and ascending in panel " + panel);
              previous = value.get<int>();
              axes.push_back(previous);
            }
          selected.emplace(panel,std::move(axes));
        }
      return selected;
    }
  catch ( const std::exception & e )
    {
      Helper::halt("invalid HDSTATS-PC selection " + file + ": " + e.what());
    }
  return component_selection_t(); // unreachable after Helper::halt()
}

static std::vector<int> identity_components( const int n )
{
  std::vector<int> axes;
  for ( int i = 1; i <= n; ++i ) axes.push_back(i);
  return axes;
}

static void select_panel_components( fitted_panel_t * panel, const std::vector<int> & axes )
{
  const model_t & full = panel->model;
  if ( axes.empty() ) Helper::halt("HDSTATS-PC selection cannot remove every component from " + panel->name);
  for ( const int axis : axes )
    if ( axis < 1 || axis > full.loadings.cols() )
      Helper::halt("HDSTATS-PC selection requests PC" + Helper::int2str(axis) +
                   " but only " + Helper::int2str((int)full.loadings.cols()) + " axes were fitted for " + panel->name);

  model_t retained = full;
  retained.loadings.resize(full.loadings.rows(),axes.size());
  retained.singular_values.resize(axes.size());
  for ( int k = 0; k < (int)axes.size(); ++k )
    {
      retained.loadings.col(k) = full.loadings.col(axes[k]-1);
      retained.singular_values(k) = full.singular_values(axes[k]-1);
    }
  panel->model = std::move(retained);
  panel->source_components = axes;
}

// Group-quartimax: maximize sum_{g,k} (sum_{j in g} B(j,k)^2)^2,
// where B = V*R and R is orthogonal. Each unit-length column then favors a
// small number of the spec's groups, without favoring single-feature sparsity.
// There is no requirement for one axis per group, nor a signed target template.
static Eigen::MatrixXd loading_group_mass( const Eigen::MatrixXd & b,
                                           const std::vector<int> & groups,
                                           const int n_groups )
{
  Eigen::MatrixXd mass = Eigen::MatrixXd::Zero(n_groups,b.cols());
  for ( int j = 0; j < b.rows(); ++j ) mass.row(groups[j]) += b.row(j).array().square().matrix();
  return mass;
}

static void rotate_columns( Eigen::MatrixXd * b, const int i, const int j, const double angle )
{
  const Eigen::VectorXd x = b->col(i), y = b->col(j);
  const double c = std::cos(angle), s = std::sin(angle);
  b->col(i) = c*x + s*y;
  b->col(j) = -s*x + c*y;
}

static void rotate_panel( fitted_panel_t * panel, const panel_def_t & definition )
{
  const Eigen::MatrixXd v = panel->model.loadings;
  const int nc = v.cols();
  std::map<std::string,int> feature_group;
  for ( int j = 0; j < (int)definition.names.size(); ++j )
    {
      const auto & id = definition.groups[j];
      auto pos = std::find(panel->group_ids.begin(),panel->group_ids.end(),id);
      if ( pos == panel->group_ids.end() )
        { panel->group_ids.push_back(id); pos = panel->group_ids.end()-1; }
      feature_group[definition.names[j]] = pos-panel->group_ids.begin();
    }
  std::vector<int> groups;
  for ( const auto & id : panel->model.names ) groups.push_back(feature_group.at(id));
  const int ng = panel->group_ids.size();
  const auto objective = [&]( const Eigen::MatrixXd & b ) {
    return loading_group_mass(b,groups,ng).squaredNorm();
  };
  panel->rotation_before = objective(v);
  double best = -1.0;
  Eigen::MatrixXd best_rotation;
  // Fixed starts, fixed pair order, exact two-column updates: no RNG or new
  // user input. Ties retain the earliest start (ordinary PCA is start zero).
  // The finite search is reproducible, but is not a global-optimum guarantee.
  const int starts = nc > 1 && ng > 1 ? 8 : 1;
  const int max_sweeps = 1000;
  const double tolerance = 1e-10;
  for ( int start = 0; start < starts; ++start )
    {
      Eigen::MatrixXd rotation = Eigen::MatrixXd::Identity(nc,nc);
      if ( start )
        for ( int i = 0; i < nc; ++i )
          for ( int j = i+1; j < nc; ++j )
            rotate_columns(&rotation,i,j,0.7*std::sin(1.0*start*(i+1)*(j+2)));
      Eigen::MatrixXd b = v*rotation;
      double value = objective(b);
      bool converged = false;
      int sweeps = 0;
      for ( ; sweeps < max_sweeps; ++sweeps )
        {
          const double before = value;
          for ( int i = 0; i < nc; ++i )
            for ( int j = i+1; j < nc; ++j )
              {
                Eigen::VectorXd a = Eigen::VectorXd::Zero(ng), d = a;
                for ( int f = 0; f < b.rows(); ++f )
                  {
                    a(groups[f]) += b(f,i)*b(f,i)-b(f,j)*b(f,j);
                    d(groups[f]) += 2*b(f,i)*b(f,j);
                  }
                const double aa = a.squaredNorm()-d.squaredNorm(), bb = 2*a.dot(d);
                const double gain = 0.25*(std::hypot(aa,bb)-aa);
                // Avoid arbitrary rotations in an exactly flat group subspace.
                if ( gain <= 1e-14 ) continue;
                const double angle = 0.25*std::atan2(bb,aa);
                rotate_columns(&b,i,j,angle);
                rotate_columns(&rotation,i,j,angle);
              }
          value = objective(b);
          if ( value-before <= tolerance*std::max(1.0,before) )
            { converged = true; ++sweeps; break; }
        }
      if ( converged && value > best + tolerance )
        { best = value; best_rotation = rotation; panel->rotation_sweeps = sweeps; }
    }
  if ( best_rotation.size() == 0 || best + 1e-8 < panel->rotation_before )
    Helper::halt("HDSTATS-PC group rotation did not converge for " + panel->name);

  // Rotating unit eigenvectors preserves the orthonormal basis and its total
  // variance; scores are generally correlated. Their SDs must be recomputed.
  const Eigen::VectorXd ss = panel->model.singular_values.array().square();
  const Eigen::MatrixXd gram = best_rotation.transpose()*ss.asDiagonal()*best_rotation;
  std::vector<int> order;
  for ( int k = 0; k < nc; ++k ) order.push_back(k);
  std::stable_sort(order.begin(),order.end(),[&](int a,int b) { return gram(a,a)>gram(b,b); });
  panel->rotation.resize(nc,nc);
  for ( int k = 0; k < nc; ++k ) panel->rotation.col(k) = best_rotation.col(order[k]);
  panel->model.loadings = v*panel->rotation;
  for ( int k = 0; k < nc; ++k )
    {
      Eigen::Index anchor;
      panel->model.loadings.col(k).cwiseAbs().maxCoeff(&anchor);
      if ( panel->model.loadings(anchor,k) < 0 )
        { panel->model.loadings.col(k) *= -1; panel->rotation.col(k) *= -1; }
    }
  const Eigen::MatrixXd covariance = panel->rotation.transpose()*ss.asDiagonal()*panel->rotation/(panel->n_fit-1.0);
  panel->score_sds = covariance.diagonal().array().sqrt();
  panel->score_correlations = panel->score_sds.cwiseInverse().asDiagonal()*covariance*panel->score_sds.cwiseInverse().asDiagonal();
  panel->group_mass = loading_group_mass(panel->model.loadings,groups,ng);
  panel->rotation_after = panel->group_mass.squaredNorm();
}

// Post-rotation selection: pick a subset of an ALREADY ROTATED panel's
// columns (axes are one-based indices into the rotated PANEL,PC numbering,
// i.e. rotation's columns / model.loadings' columns). Unlike
// select_panel_components(), this never triggers a further rotation: the
// retained axes are exactly the ones the caller chose to keep. The full
// rotated basis width is preserved in model.singular_values and in
// rotation's row count, so the saved model can still report basis
// provenance separately from the final (reduced) score scaling.
static void select_rotated_components( fitted_panel_t * panel, const std::vector<int> & axes )
{
  const int basis = (int)panel->model.loadings.cols();
  if ( axes.empty() ) Helper::halt("HDSTATS-PC post-rotation selection cannot remove every component from " + panel->name);
  for ( const int axis : axes )
    if ( axis < 1 || axis > basis )
      Helper::halt("HDSTATS-PC post-rotation selection requests PC" + Helper::int2str(axis) +
                   " but only " + Helper::int2str(basis) + " rotated axes were fitted for " + panel->name);

  const Eigen::MatrixXd loadings = panel->model.loadings;
  const Eigen::MatrixXd rotation = panel->rotation;
  const Eigen::VectorXd score_sds = panel->score_sds;
  const Eigen::MatrixXd score_correlations = panel->score_correlations;
  const Eigen::MatrixXd group_mass = panel->group_mass;
  const int k = (int)axes.size();

  panel->model.loadings.resize(loadings.rows(),k);
  panel->rotation.resize(rotation.rows(),k);
  panel->score_sds.resize(k);
  panel->score_correlations.resize(k,k);
  panel->group_mass.resize(group_mass.rows(),k);
  for ( int a = 0; a < k; ++a )
    {
      panel->model.loadings.col(a) = loadings.col(axes[a]-1);
      panel->rotation.col(a) = rotation.col(axes[a]-1);
      panel->score_sds(a) = score_sds(axes[a]-1);
      panel->group_mass.col(a) = group_mass.col(axes[a]-1);
      for ( int b = 0; b < k; ++b )
        panel->score_correlations(a,b) = score_correlations(axes[a]-1,axes[b]-1);
    }
  // model.singular_values, rotation_before/after, rotation_sweeps, group_ids,
  // and source_components are deliberately left untouched: they describe
  // the full rotated basis this selection was taken from.
  panel->selected_components = axes;
}

static json matrix_json( const Eigen::MatrixXd & matrix )
{
  json result = json::array();
  for ( int i = 0; i < matrix.rows(); ++i )
    {
      json row = json::array();
      for ( int j = 0; j < matrix.cols(); ++j ) row.push_back(matrix(i,j));
      result.push_back(std::move(row));
    }
  return result;
}

static void write_json_panels( const std::string & file, const layout_t & layout,
                               const std::vector<fitted_panel_t> & panels,
                               const bool post_rotation_selection = false )
{
  json document;
  document["format"] = "luna-hdstats-pc-model";
  document["version"] = post_rotation_selection ? 3 : ( panels.front().rotation.size() ? 2 : 1 );
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
      item["source_components"] = panel.source_components.empty()
        ? identity_components((int)model.singular_values.size()) : panel.source_components;
      item["singular_values"] = json::array();
      for ( int k = 0; k < model.singular_values.size(); ++k ) item["singular_values"].push_back(model.singular_values(k));
      if ( post_rotation_selection )
        {
          item["n_basis"] = (int)model.singular_values.size();
          item["selected_components"] = panel.selected_components.empty()
            ? identity_components(model.loadings.cols()) : panel.selected_components;
        }
      if ( panel.rotation.size() )
        {
          item["rotation"] = {
            {"method","group_quartimax_v1"}, {"matrix",matrix_json(panel.rotation)},
            {"objective_before",panel.rotation_before}, {"objective_after",panel.rotation_after},
            {"starts",panel.rotation.cols()>1 && panel.group_ids.size()>1 ? 8 : 1},
            {"sweeps",panel.rotation_sweeps}, {"converged",true},
            {"group_ids",panel.group_ids}, {"group_mass",matrix_json(panel.group_mass)}
          };
          item["score_sds"] = json::array();
          for ( int k = 0; k < panel.score_sds.size(); ++k ) item["score_sds"].push_back(panel.score_sds(k));
          item["score_correlations"] = matrix_json(panel.score_correlations);
        }
      item["features"] = json::array();
      for ( int j = 0; j < (int)model.names.size(); ++j )
        {
          json feature;
          feature["id"] = model.names[j];
          feature["mean"] = model.means(j);
          feature["sd"] = model.sds(j);
          feature["weight"] = model.weights(j);
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
                z(i,j)=m.weights(j) * (row.values.at(m.names[j])-m.means(j))/m.sds(j);
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
  // `nc=` is an intentional exploratory override.  Persist it in the embedded
  // specification so projection never needs the source specification again.
  for ( const auto & definition : panel_definitions() )
    if ( param.has(definition.option) )
      Helper::halt("HDSTATS-PC spec= cannot be combined with " + definition.option + "=");

  layout_t layout = read_layout_spec( Helper::expand( param.requires("spec") ) );
  if ( param.has("nc") )
    {
      const int nc = param.requires_int("nc");
      if ( nc < 1 ) Helper::halt("nc must be positive");
      for ( int i = 0; i < (int)layout.panels.size(); ++i )
        {
          layout.panels[i].n_components = nc;
          layout.document["panels"][i]["n_components"] = nc;
        }
    }
  const bool rotate = param.yesno("rotate",true);
  const std::string pc_select_stage = param.has("pc-select-stage") ? param.requires("pc-select-stage") : "pre-rotation";
  if ( pc_select_stage != "pre-rotation" && pc_select_stage != "post-rotation" )
    Helper::halt("pc-select-stage must be pre-rotation or post-rotation");
  const bool post_rotation_selection = pc_select_stage == "post-rotation";
  if ( post_rotation_selection && !rotate )
    Helper::halt("pc-select-stage=post-rotation requires rotate=T");
  const component_selection_t selection = param.has("pc-select")
    ? read_component_selection( Helper::expand(param.requires("pc-select")), layout )
    : component_selection_t();
  const auto files = param.strvector("data");
  if ( files.empty() ) Helper::halt("--hdstats-pc requires data=<feature files>");
  const auto rows = with_nrem_features( read_rows(files,hdstats_format) );
  std::vector<fitted_panel_t> panels;
  for ( int p = 0; p < (int)layout.panels.size(); ++p )
    {
      layout_panel_t & definition = layout.panels[p];
      panels.push_back( fit_panel( definition.definition, rows, definition.n_components ) );
      const auto found = selection.find(definition.definition.name);

      if ( post_rotation_selection )
        {
          // Fit and rotate the FULL requested basis first; selection then
          // only drops already-rotated columns, never triggering a further
          // rotation of the hand-picked subset.
          select_panel_components(&panels.back(),identity_components(definition.n_components));
          rotate_panel(&panels.back(),definition.definition);
          const int basis = definition.n_components;
          const std::vector<int> axes = found == selection.end()
            ? identity_components(basis) : found->second;
          select_rotated_components(&panels.back(),axes);
          definition.n_components = (int)axes.size();
          layout.document["panels"][p]["n_components"] = definition.n_components;
          logger << "  HDSTATS-PC " << definition.definition.name << ": " << panels.back().n_fit << "/" << rows.size()
                 << " people, " << panels.back().model.names.size() << "/" << definition.definition.names.size()
                 << " varying inputs, nc=" << definition.n_components << " of " << basis << " rotated PCs";
          if ( found != selection.end() )
            {
              logger << " (post-rotation selection from";
              for ( const int axis : axes ) logger << " PC" << axis;
              logger << ")";
            }
          logger << ", group rotation " << panels.back().rotation_before << " -> " << panels.back().rotation_after << "\n";
          continue;
        }

      const std::vector<int> axes = found == selection.end()
        ? identity_components(definition.n_components) : found->second;
      select_panel_components(&panels.back(),axes);
      if ( rotate ) rotate_panel(&panels.back(),definition.definition);
      // The embedded spec defines the columns in this deployed model, not the
      // larger exploratory PCA used while choosing them.
      definition.n_components = (int)axes.size();
      layout.document["panels"][p]["n_components"] = definition.n_components;
      logger << "  HDSTATS-PC " << definition.definition.name << ": " << panels.back().n_fit << "/" << rows.size()
             << " people, " << panels.back().model.names.size() << "/" << definition.definition.names.size()
             << " varying inputs, nc=" << definition.n_components;
      if ( found != selection.end() )
        {
          logger << " (selected from";
          for ( const int axis : axes ) logger << " PC" << axis;
          logger << ")";
        }
      if ( rotate )
        logger << ", group rotation " << panels.back().rotation_before << " -> " << panels.back().rotation_after;
      logger << "\n";
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
                z(i,j) = model.weights(j) * ( row.values.at(model.names[j]) - model.means(j) ) / model.sds(j);
              ++i;
            }
          write_correlation_matrix(prefix + "." + panel.name + ".tsv",model.names,z,hdstats_format);
        }
    }

  write_json_panels( param.requires("model"), layout, panels, post_rotation_selection );
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
      writer.value("ROTATED",panel.rotation.size() ? 1 : 0);
      for ( int j = 0; j < (int)model.names.size(); ++j )
        {
          writer.level(model.names[j],"FTR");
          writer.value("MEAN",model.means(j)); writer.value("SD",model.sds(j));
          writer.unlevel("FTR");
        }
      double cumulative = 0.0;
      for ( int k = 0; k < model.loadings.cols(); ++k )
        {
          const double score_sd = panel.rotation.size() ? panel.score_sds(k)
            : model.singular_values(k)/std::sqrt(panel.n_fit-1.0);
          const double w = panel.rotation.size() ? score_sd*std::sqrt(panel.n_fit-1.0) : model.singular_values(k);
          const double ve = w*w / panel.total_ss;
          cumulative += ve;
          writer.level("PC" + Helper::int2str(k+1),"PC");
          if ( ! panel.rotation.size() ) writer.value("SOURCE_PC",panel.source_components.empty() ? k+1 : panel.source_components[k]);
          writer.value("W",w); writer.value("SCORE_SD",score_sd); writer.value("VE",ve); writer.value("CVE",cumulative);
          if ( panel.rotation.size() )
            {
              Eigen::Index dominant;
              const double mass = panel.group_mass.col(k).maxCoeff(&dominant);
              writer.value("GROUP",panel.group_ids[dominant]); writer.value("GROUP_MASS",mass);
              if ( !panel.selected_components.empty() ) writer.value("DISCOVERY_PC",panel.selected_components[k]);
              for ( int g = 0; g < (int)panel.group_ids.size(); ++g )
                { writer.level(panel.group_ids[g],"GROUP"); writer.value("MASS",panel.group_mass(g,k)); writer.unlevel("GROUP"); }
              for ( int b = 0; b < panel.rotation.rows(); ++b )
                { writer.level(panel.source_components[b],"BASIS_PC"); writer.value("ROTATION",panel.rotation(b,k)); writer.unlevel("BASIS_PC"); }
            }
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
      const double score_sd = source.score_unit_variance
        ? panel->score_sds[definition.panel_pc-1] : 1.0;
      if ( ! std::isfinite(score_sd) || score_sd <= 0.0 )
        Helper::halt("invalid HDSTATS-PC score scale in " + file);
      double score = 0.0;
      for ( const auto & feature : panel->features )
        {
          const auto raw = features.find(feature.id);
          const bool present = raw != features.end() && std::isfinite(raw->second);
          if ( ! present ) { ++result.missing; result.compatible = false; }
          const double z = present ? feature.weight * ( raw->second - feature.mean ) / feature.sd : 0.0;
          const double loading = definition.orientation * feature.loadings[definition.panel_pc-1];
          score += loading * z;
          if ( std::fabs(loading) >= min_loading )
            result.contributions.push_back( { feature.id, present ? loading * z / score_sd : nan,
                                               present ? raw->second : nan } );
        }
      if ( result.compatible ) result.score = score / score_sd;
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
  if ( !source.panels.empty() && source.panels.front().rotation_method != "none" )
    out << "# Rotated components (group_quartimax_v1); PC indices refer to final model columns.\n";
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
  feature_map_t effective_features = features;
  derive_nrem_features(&effective_features);
  const bool json_model = is_json_hdstats_model(file);
  const hdstats_model::model_t json_source = json_model ? read_json_hdstats_model(file) : hdstats_model::model_t();
  std::ifstream in=LunaIO::open_ifstream(file,std::ios::in);
  std::string magic;
  if ( !std::getline(in,magic) ) Helper::halt("could not read HDSTATS-PC model " + file);
  if ( magic==hdstats_format.model_magic )
    {
      panel_scores_t result;
      result.panel="GENERIC";
      project_impl(file,effective_features,&result.scores,&result.missing,hdstats_format);
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
          const auto f=effective_features.find(m.names[j]);
          if ( f==effective_features.end() || !std::isfinite(f->second) )
            { ++result.missing; z(j)=0.0; }
          else z(j)=m.weights(j) * (f->second-m.means(j))/m.sds(j);
        }
      if ( result.missing==0 )
        {
          const Eigen::VectorXd s=m.loadings.transpose()*z;
          result.scores.resize(s.size());
          const hdstats_model::panel_t * source_panel = json_model ? json_source.find_panel(p.name) : NULL;
          for ( int k = 0; k < s.size(); ++k )
            {
              const double score_sd = json_model && json_source.score_unit_variance
                ? source_panel->score_sds[k] : 1.0;
              if ( ! std::isfinite(score_sd) || score_sd <= 0.0 )
                Helper::halt("invalid HDSTATS-PC score scale in " + file);
              result.scores[k] = s(k) / score_sd;
            }
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
  feature_map_t effective_features = features;
  derive_nrem_features(&effective_features);
  if ( is_json_hdstats_model(file) ) return project_json_components(file,effective_features,min_loading);
  std::ifstream in=LunaIO::open_ifstream(file,std::ios::in);
  std::string magic;
  if ( !std::getline(in,magic) ) Helper::halt("could not read HDSTATS-PC model " + file);
  if ( magic == hdstats_format.model_magic )
    Helper::halt("HDSTATS-PC global component projection requires a V3 model: " + file);
  if ( magic != "__HDSTATS_PC_V3" ) Helper::halt("unsupported HDSTATS-PC model " + file);

  const auto & fitted = cached_panels(file);
  validate_projection_model(fitted);
  std::vector<component_score_t> results;
  const double nan = std::numeric_limits<double>::quiet_NaN();
  for ( const auto & panel : fitted )
    {
      const model_t & model = panel.model;
      Eigen::VectorXd z(model.names.size());
      int missing = 0;
      for ( int j=0; j<(int)model.names.size(); ++j )
        {
          const auto feature=effective_features.find(model.names[j]);
          if ( feature==effective_features.end() || !std::isfinite(feature->second) )
            { ++missing; z(j)=0.0; }
          else z(j)=model.weights(j) * (feature->second-model.means(j))/model.sds(j);
        }
      for ( int pc = 0; pc < model.loadings.cols(); ++pc )
        {
          component_score_t result;
          result.id=panel.name + "_" + Helper::int2str(pc+1);
          result.panel=panel.name;
          result.panel_pc=pc+1;
          result.orientation=1;
          result.compatible=missing == 0;
          result.missing=missing;
          result.score=result.compatible ? model.loadings.col(pc).dot(z) : nan;
          for ( int j = 0; j < (int)model.names.size(); ++j )
            {
              const double loading = model.loadings(j,pc);
              if ( std::fabs(loading) < min_loading ) continue;
              const auto raw = effective_features.find(model.names[j]);
              result.contributions.push_back( { model.names[j], result.compatible ? loading * z(j) : nan,
                raw == effective_features.end() || !std::isfinite(raw->second) ? nan : raw->second } );
            }
          results.push_back(std::move(result));
        }
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
  if ( ! descriptor.empty() ) Helper::halt("HDSTATS descriptors are supported for self-contained JSON models only");
  const auto & fitted = cached_panels(file);
  validate_projection_model(fitted);
  out << "PC\tDOMAIN\tMODEL_PC\tLABEL\tORIENTATION\tFTR_LOADINGS\n";
  for ( const auto & panel : fitted )
    {
      const model_t & model = panel.model;
      for ( int pc = 0; pc < model.loadings.cols(); ++pc )
        {
          std::vector<std::string> exemplar_features;
          for ( int j = 0; j < (int)model.names.size(); ++j )
            if ( std::fabs(model.loadings(j,pc)) >= min_loading )
              exemplar_features.push_back(model.names[j] + "=" + fmt_display(model.loadings(j,pc)));
          out << panel.name << '_' << pc + 1 << '\t' << panel.name << '\t' << pc + 1 << "\t\t1\t";
          for ( int j = 0; j < (int)exemplar_features.size(); ++j )
            { if ( j ) out << "; "; out << exemplar_features[j]; }
          out << '\n';
        }
    }
}

void hdstats_pc::fit( param_t & param )
{
  if ( !param.has("spec") )
    Helper::halt("--hdstats-pc requires spec=docs/hdstats/hdstats-pc.spec.json");
  fit_json_hdstats_panels( param );
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
