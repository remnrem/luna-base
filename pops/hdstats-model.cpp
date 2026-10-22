#include "hdstats-model.h"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <set>
#include <sstream>

namespace hdstats_model {
namespace {

std::vector<std::string> split_tab( const std::string & line )
{
  std::vector<std::string> fields;
  std::string field;
  std::istringstream stream(line);
  while ( std::getline(stream,field,'\t') ) fields.push_back(field);
  if ( !line.empty() && line.back() == '\t' ) fields.push_back("");
  return fields;
}

bool parse_int( const std::string & text, int * value )
{
  if ( text.empty() ) return false;
  char * end = nullptr;
  errno = 0;
  const long parsed = std::strtol(text.c_str(),&end,10);
  if ( errno || end != text.c_str() + text.size() ||
       parsed < std::numeric_limits<int>::min() || parsed > std::numeric_limits<int>::max() )
    return false;
  *value = (int)parsed;
  return true;
}

bool parse_number( const std::string & text, double * value )
{
  if ( text.empty() ) return false;
  char * end = nullptr;
  errno = 0;
  const double parsed = std::strtod(text.c_str(),&end);
  if ( errno || end != text.c_str() + text.size() || !std::isfinite(parsed) ) return false;
  *value = parsed;
  return true;
}

bool entity_type( const std::string & type )
{
  return type == "DOMAIN" || type == "PC" || type == "FEATURE";
}

std::string pair_key( const std::string & first, const std::string & second )
{
  return first < second ? first + "\t" + second : second + "\t" + first;
}

bool ignorable_annotation_line( const std::string & line )
{
  return line.empty() || line[0] == '#';
}

struct reader_t {
  explicit reader_t( const std::string & filename ) : file(filename), in(filename.c_str()) {}

  bool next( std::string * line )
  {
    if ( !std::getline(in,*line) ) return false;
    ++line_number;
    return true;
  }

  bool fail( const std::string & message, std::string * error ) const
  {
    if ( error )
      *error = "HDSTATS model " + file + ", line " + std::to_string(line_number) + ": " + message;
    return false;
  }

  std::string file;
  std::ifstream in;
  int line_number = 0;
};

bool read_annotations( reader_t * reader, model_t * model, std::string * error )
{
  std::string line;
  bool found = false;
  while ( reader->next(&line) )
    {
      if ( ignorable_annotation_line(line) ) continue;
      if ( line != "__HDSTATS_PC_ANNOTATIONS_V1" )
        return reader->fail("unexpected trailing model data",error);
      found = true;
      break;
    }
  if ( !found ) return true;

  std::set<std::string> descriptions;
  std::set<std::string> associations;
  std::set<std::string> pc_pairs;
  std::set<std::string> feature_pairs;
  bool ended = false;
  while ( reader->next(&line) )
    {
      if ( ignorable_annotation_line(line) ) continue;
      if ( line == "END_HDSTATS_PC_ANNOTATIONS" )
        { ended = true; break; }
      const std::vector<std::string> fields = split_tab(line);
      if ( fields.empty() ) return reader->fail("empty annotation record",error);
      if ( fields[0] == "DESC" )
        {
          if ( fields.size() != 4 || !entity_type(fields[1]) || fields[2].empty() || fields[3].empty() )
            return reader->fail("bad DESC annotation",error);
          const std::string key = fields[1] + "\t" + fields[2];
          if ( !descriptions.insert(key).second ) return reader->fail("duplicate DESC annotation",error);
          model->descriptions.push_back( { fields[1], fields[2], fields[3] } );
        }
      else if ( fields[0] == "EXT_ASSOC" )
        {
          double score;
          if ( fields.size() != 5 || !entity_type(fields[1]) || fields[2].empty() || fields[3].empty() ||
               !parse_number(fields[4],&score) ) return reader->fail("bad EXT_ASSOC annotation",error);
          const std::string key = fields[1] + "\t" + fields[2] + "\t" + fields[3];
          if ( !associations.insert(key).second ) return reader->fail("duplicate EXT_ASSOC annotation",error);
          model->external_associations.push_back( { fields[1], fields[2], fields[3], score } );
        }
      else if ( fields[0] == "PC_CORR" || fields[0] == "FTR_CORR" )
        {
          double value;
          if ( fields.size() != 4 || fields[1].empty() || fields[2].empty() || fields[1] == fields[2] ||
               !parse_number(fields[3],&value) ) return reader->fail("bad correlation annotation",error);
          std::set<std::string> & pairs = fields[0] == "PC_CORR" ? pc_pairs : feature_pairs;
          if ( !pairs.insert(pair_key(fields[1],fields[2])).second )
            return reader->fail("duplicate unordered correlation annotation",error);
          correlation_t record = { fields[1], fields[2], value };
          if ( fields[0] == "PC_CORR" ) model->pc_correlations.push_back(record);
          else model->feature_correlations.push_back(record);
        }
      else return reader->fail("unknown HDSTATS model annotation record",error);
    }
  if ( !ended ) return reader->fail("unterminated HDSTATS model annotations",error);
  while ( reader->next(&line) )
    if ( !ignorable_annotation_line(line) ) return reader->fail("trailing data after model annotations",error);
  return true;
}

}

const panel_t * model_t::find_panel( const std::string & id ) const
{
  for ( const auto & panel : panels ) if ( panel.id == id ) return &panel;
  return nullptr;
}

const feature_t * model_t::find_feature( const std::string & panel_id,
                                         const std::string & feature_id ) const
{
  const panel_t * panel = find_panel(panel_id);
  if ( !panel ) return nullptr;
  for ( const auto & feature : panel->features ) if ( feature.id == feature_id ) return &feature;
  return nullptr;
}

const description_t * model_t::find_description( const std::string & entity_type,
                                                 const std::string & entity_id ) const
{
  for ( const auto & description : descriptions )
    if ( description.entity_type == entity_type && description.entity_id == entity_id ) return &description;
  return nullptr;
}

bool read( const std::string & file, model_t * model, std::string * error )
{
  if ( error ) error->clear();
  if ( !model )
    {
      if ( error ) *error = "HDSTATS model reader requires a non-null model pointer";
      return false;
    }
  reader_t reader(file);
  if ( !reader.in )
    {
      if ( error ) *error = "could not open HDSTATS model " + file;
      return false;
    }
  model_t parsed;
  std::string line;
  if ( !reader.next(&line) || line != "__HDSTATS_PC_V3" )
    return reader.fail("expected __HDSTATS_PC_V3",error);
  if ( !reader.next(&line) ) return reader.fail("missing panel count",error);
  std::vector<std::string> fields = split_tab(line);
  int n_panels;
  if ( fields.size() != 2 || fields[0] != "NP" || !parse_int(fields[1],&n_panels) || n_panels < 1 )
    return reader.fail("bad panel count",error);

  std::set<std::string> panel_ids;
  for ( int p = 0; p < n_panels; ++p )
    {
      if ( !reader.next(&line) ) return reader.fail("truncated panel header",error);
      fields = split_tab(line);
      int n_input, n_fit, n_declared, n_features, n_components;
      double total_ss;
      if ( fields.size() != 8 || fields[0] != "PANEL" || fields[1].empty() ||
           !panel_ids.insert(fields[1]).second || !parse_int(fields[2],&n_input) ||
           !parse_int(fields[3],&n_fit) || !parse_int(fields[4],&n_declared) ||
           !parse_int(fields[5],&n_features) || !parse_int(fields[6],&n_components) ||
           !parse_number(fields[7],&total_ss) || n_input < 2 || n_fit < 2 || n_declared < 1 ||
           n_features < 1 || n_components < 1 || n_components > n_features || total_ss <= 0.0 )
        return reader.fail("bad panel header",error);
      panel_t panel;
      panel.id = fields[1]; panel.n_input = n_input; panel.n_fit = n_fit; panel.total_ss = total_ss;
      std::set<std::string> declared;
      for ( int j = 0; j < n_declared; ++j )
        {
          if ( !reader.next(&line) ) return reader.fail("truncated declared features",error);
          fields = split_tab(line);
          if ( fields.size() != 2 || fields[0] != "I" || fields[1].empty() || !declared.insert(fields[1]).second )
            return reader.fail("bad declared feature",error);
          panel.declared_features.push_back(fields[1]);
        }
      std::set<std::string> feature_ids;
      for ( int j = 0; j < n_features; ++j )
        {
          if ( !reader.next(&line) ) return reader.fail("truncated feature statistics",error);
          fields = split_tab(line);
          feature_t feature;
          if ( fields.size() != 4 || fields[0] != "F" || fields[1].empty() || !declared.count(fields[1]) ||
               !feature_ids.insert(fields[1]).second || !parse_number(fields[2],&feature.mean) ||
               !parse_number(fields[3],&feature.sd) || feature.sd <= 0.0 )
            return reader.fail("bad feature statistics",error);
          feature.id = fields[1];
          panel.features.push_back(feature);
        }
      if ( !reader.next(&line) ) return reader.fail("truncated singular values",error);
      fields = split_tab(line);
      if ( (int)fields.size() != n_components + 1 || fields[0] != "W" )
        return reader.fail("bad singular values",error);
      for ( int k = 0; k < n_components; ++k )
        {
          double value;
          if ( !parse_number(fields[k+1],&value) ) return reader.fail("bad singular value",error);
          panel.singular_values.push_back(value);
        }
      for ( int j = 0; j < n_features; ++j )
        {
          if ( !reader.next(&line) ) return reader.fail("truncated feature loadings",error);
          fields = split_tab(line);
          if ( (int)fields.size() != n_components + 2 || fields[0] != "V" || fields[1] != panel.features[j].id )
            return reader.fail("bad feature loading",error);
          for ( int k = 0; k < n_components; ++k )
            {
              double value;
              if ( !parse_number(fields[k+2],&value) ) return reader.fail("bad feature loading value",error);
              panel.features[j].loadings.push_back(value);
            }
        }
      if ( !reader.next(&line) || line != "END" ) return reader.fail("bad panel end",error);
      parsed.panels.push_back(panel);
    }
  if ( !read_annotations(&reader,&parsed,error) ) return false;
  *model = parsed;
  return true;
}

}
