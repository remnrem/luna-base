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

void hdstats_pc::fit( param_t & param )
{
  fit_impl( param, hdstats_format );
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
