#include "hdstats-model.h"

#include "helper/json.h"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <regex>
#include <set>
#include <sstream>

namespace hdstats_model {
namespace {

using json = nlohmann::json;

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

bool json_string( const json & object, const char * key, std::string * value )
{
  if ( ! object.contains(key) || ! object.at(key).is_string() || object.at(key).get<std::string>().empty() ) return false;
  *value = object.at(key).get<std::string>();
  return true;
}

bool json_number( const json & object, const char * key, double * value )
{
  if ( ! object.contains(key) || ! object.at(key).is_number() ) return false;
  *value = object.at(key).get<double>();
  return std::isfinite(*value);
}

bool json_positive_int( const json & object, const char * key, int * value )
{
  if ( ! object.contains(key) || ! object.at(key).is_number_integer() ) return false;
  *value = object.at(key).get<int>();
  return *value > 0;
}

bool r_safe_name( const std::string & value )
{
  static const std::regex pattern("[A-Za-z][A-Za-z0-9_.]*");
  return std::regex_match(value,pattern);
}

bool read_json_model( const std::string & file, model_t * model, std::string * error )
{
  try
    {
      std::ifstream in(file.c_str());
      if ( ! in )
        { if ( error ) *error = "could not open HDSTATS model " + file; return false; }
      const json document = json::parse(in);
      if ( ! document.is_object() || document.value("format","") != "luna-hdstats-pc-model" ||
           document.value("version",0) != 1 || ! document.contains("spec") ||
           ! document.at("spec").is_object() || ! document.contains("panels") || ! document.at("panels").is_array() )
        { if ( error ) *error = "HDSTATS JSON model " + file + ": expected V1 model object"; return false; }

      const std::string component_id_style = document.at("spec").value("component_id_style",std::string("legacy"));
      if ( component_id_style != "legacy" && component_id_style != "r_safe" )
        { if ( error ) *error = "HDSTATS JSON model " + file + ": invalid component_id_style"; return false; }

      model_t parsed;
      if ( ! json_string(document,"model_id",&parsed.model_id) ||
           document.at("spec").value("format","") != "luna-hdstats-pc-spec" ||
           document.at("spec").value("version",0) != 1 ||
           document.at("spec").value("feature_universe","") != "HDSTATS-PC-FEATURES-V1" ||
           ! document.at("spec").contains("feature_descriptions") || ! document.at("spec").at("feature_descriptions").is_array() ||
           ! document.at("spec").contains("panels") || ! document.at("spec").at("panels").is_array() ||
           document.at("spec").at("panels").empty() )
        { if ( error ) *error = "HDSTATS JSON model " + file + ": invalid embedded specification"; return false; }
      if ( document.contains("fingerprint") )
        {
          if ( ! document.at("fingerprint").is_string() || document.at("fingerprint").get<std::string>().empty() )
            { if ( error ) *error = "HDSTATS JSON model " + file + ": invalid fingerprint"; return false; }
          parsed.fingerprint = document.at("fingerprint").get<std::string>();
        }
      parsed.specification_json = document.at("spec").dump();
      parsed.component_id_style = component_id_style;

      std::map<std::string,std::pair<std::string,std::string>> metadata;
      for ( const auto & item : document.at("spec").at("feature_descriptions") )
        {
          std::string id,label,description;
          if ( ! item.is_object() || ! json_string(item,"id",&id) || ! json_string(item,"label",&label) ||
               ! json_string(item,"description",&description) || label.empty() || description.empty() ||
               ! metadata.emplace(id,std::make_pair(label,description)).second )
            { if ( error ) *error = "HDSTATS JSON model " + file + ": invalid embedded feature metadata"; return false; }
        }

      std::map<std::string,const json *> fit_panels;
      for ( const auto & item : document.at("panels") )
        {
          std::string id;
          if ( ! item.is_object() || ! json_string(item,"id",&id) || ! fit_panels.emplace(id,&item).second )
            { if ( error ) *error = "HDSTATS JSON model " + file + ": invalid or duplicate fitted panel"; return false; }
        }

      std::set<std::string> component_ids, used_features;
      for ( const auto & spec_panel : document.at("spec").at("panels") )
        {
          std::string panel_id;
          if ( ! spec_panel.is_object() || ! json_string(spec_panel,"id",&panel_id) ||
               ! spec_panel.contains("features") || ! spec_panel.at("features").is_array() || spec_panel.at("features").empty() ||
               ! spec_panel.contains("n_components") || ! spec_panel.at("n_components").is_number_integer() ||
               spec_panel.at("n_components").get<int>() < 1 )
            { if ( error ) *error = "HDSTATS JSON model " + file + ": invalid specification panel"; return false; }
          if ( component_id_style == "r_safe" && ! r_safe_name(panel_id) )
            { if ( error ) *error = "HDSTATS JSON model " + file + ": panel IDs must be R-safe name prefixes"; return false; }
          const auto fitted = fit_panels.find(panel_id);
          if ( fitted == fit_panels.end() )
            { if ( error ) *error = "HDSTATS JSON model " + file + ": missing fitted panel " + panel_id; return false; }
          const json & source = *fitted->second;
          panel_t panel;
          panel.id = panel_id;
          if ( ! json_positive_int(source,"n_input",&panel.n_input) || ! json_positive_int(source,"n_fit",&panel.n_fit) ||
               panel.n_fit < 2 || ! json_number(source,"total_ss",&panel.total_ss) || panel.total_ss <= 0 ||
               ! source.contains("declared_features") || ! source.at("declared_features").is_array() ||
               ! source.contains("features") || ! source.at("features").is_array() ||
               ! source.contains("singular_values") || ! source.at("singular_values").is_array() )
            { if ( error ) *error = "HDSTATS JSON model " + file + ": invalid fitted panel " + panel_id; return false; }

          std::set<std::string> declared;
          for ( const auto & feature : source.at("declared_features") )
            if ( ! feature.is_string() || feature.get<std::string>().empty() || ! declared.insert(feature.get<std::string>()).second )
              { if ( error ) *error = "HDSTATS JSON model " + file + ": invalid declared feature"; return false; }
            else panel.declared_features.push_back(feature.get<std::string>());
          if ( panel.declared_features.empty() )
            { if ( error ) *error = "HDSTATS JSON model " + file + ": no declared features"; return false; }

          const int n_components = spec_panel.at("n_components").get<int>();
          if ( (int)source.at("singular_values").size() != n_components || n_components > (int)source.at("features").size() )
            { if ( error ) *error = "HDSTATS JSON model " + file + ": component/loading dimensions disagree"; return false; }
          for ( const auto & value : source.at("singular_values") )
            if ( ! value.is_number() || ! std::isfinite(value.get<double>()) )
              { if ( error ) *error = "HDSTATS JSON model " + file + ": invalid singular value"; return false; }
            else panel.singular_values.push_back(value.get<double>());

          std::set<std::string> feature_ids;
          for ( const auto & source_feature : source.at("features") )
            {
              feature_t feature;
              if ( ! source_feature.is_object() || ! json_string(source_feature,"id",&feature.id) ||
                   ! declared.count(feature.id) || ! feature_ids.insert(feature.id).second ||
                   ! json_number(source_feature,"mean",&feature.mean) || ! json_number(source_feature,"sd",&feature.sd) || feature.sd <= 0 ||
                   ! source_feature.contains("loadings") || ! source_feature.at("loadings").is_array() ||
                   (int)source_feature.at("loadings").size() != n_components )
                { if ( error ) *error = "HDSTATS JSON model " + file + ": invalid fitted feature"; return false; }
              for ( const auto & loading : source_feature.at("loadings") )
                if ( ! loading.is_number() || ! std::isfinite(loading.get<double>()) )
                  { if ( error ) *error = "HDSTATS JSON model " + file + ": invalid loading"; return false; }
                else feature.loadings.push_back(loading.get<double>());
              panel.features.push_back(std::move(feature));
            }
          if ( panel.features.empty() )
            { if ( error ) *error = "HDSTATS JSON model " + file + ": no fitted features"; return false; }
          for ( const auto & declared_feature : panel.declared_features )
            {
              if ( ! metadata.count(declared_feature) )
                { if ( error ) *error = "HDSTATS JSON model " + file + ": missing feature metadata for " + declared_feature; return false; }
              used_features.insert(declared_feature);
            }
          parsed.panels.push_back(std::move(panel));

          for ( int pc = 1; pc <= n_components; ++pc )
            {
              component_t component;
              component.id = component_id_style == "r_safe"
                ? panel_id + "_" + std::to_string(pc)
                : panel_id + "/" + std::to_string(pc);
              component.label = component.id;
              component.orientation = 1;
              if ( ! component_ids.insert(component.id).second )
                { if ( error ) *error = "HDSTATS JSON model " + file + ": invalid or duplicate component ID"; return false; }
              component.panel_id = panel_id;
              component.panel_pc = pc;
              parsed.components.push_back(std::move(component));
            }
        }
      if ( used_features.size() != metadata.size() )
        { if ( error ) *error = "HDSTATS JSON model " + file + ": embedded feature metadata does not match specification panels"; return false; }
      if ( parsed.panels.size() != fit_panels.size() )
        { if ( error ) *error = "HDSTATS JSON model " + file + ": fitted panel is absent from specification"; return false; }
      *model = std::move(parsed);
      return true;
    }
  catch ( const std::exception & e )
    {
      if ( error ) *error = "HDSTATS JSON model " + file + ": " + e.what();
      return false;
    }
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
  std::ifstream probe(file.c_str());
  if ( ! probe )
    {
      if ( error ) *error = "could not open HDSTATS model " + file;
      return false;
    }
  probe >> std::ws;
  if ( probe.peek() == '{' ) return read_json_model(file,model,error);

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

bool read_descriptor( const std::string & file, model_t * model, std::string * error )
{
  if ( error ) error->clear();
  if ( ! model || model->specification_json.empty() )
    { if ( error ) *error = "HDSTATS descriptor requires a self-contained JSON V1 model"; return false; }
  try
    {
      std::ifstream in(file.c_str());
      if ( ! in ) { if ( error ) *error = "could not open HDSTATS descriptor " + file; return false; }
      const json document = json::parse(in);
      if ( ! document.is_object() || document.value("format","") != "luna-hdstats-pc-desc" ||
           document.value("version",0) != 1 || ! document.contains("model") || ! document.at("model").is_object() ||
           ! document.contains("components") || ! document.at("components").is_array() )
        { if ( error ) *error = "HDSTATS descriptor " + file + ": expected V1 descriptor object"; return false; }
      std::string model_id;
      if ( ! json_string(document.at("model"),"model_id",&model_id) || model_id != model->model_id )
        { if ( error ) *error = "HDSTATS descriptor " + file + ": model_id does not match fitted model"; return false; }
      if ( document.at("model").contains("fingerprint") )
        {
          std::string fingerprint;
          if ( ! json_string(document.at("model"),"fingerprint",&fingerprint) || model->fingerprint.empty() || fingerprint != model->fingerprint )
            { if ( error ) *error = "HDSTATS descriptor " + file + ": fingerprint does not match fitted model"; return false; }
        }

      model_t parsed = *model;
      std::set<std::string> component_keys, component_ids;
      for ( const auto & item : document.at("components") )
        {
          std::string key,id,label;
          if ( ! item.is_object() || ! json_string(item,"key",&key) || ! json_string(item,"id",&id) ||
               ! json_string(item,"label",&label) || id.empty() || label.empty() ||
               ! item.contains("orientation") || ! item.at("orientation").is_number_integer() ||
               ! component_keys.insert(key).second || ! component_ids.insert(id).second )
            { if ( error ) *error = "HDSTATS descriptor " + file + ": invalid or duplicate component record"; return false; }
          if ( parsed.component_id_style == "r_safe" && ! r_safe_name(id) )
            { if ( error ) *error = "HDSTATS descriptor " + file + ": component IDs must be R-safe names"; return false; }
          const size_t slash = key.rfind('/');
          int pc = 0;
          if ( slash == std::string::npos || slash == 0 || !parse_int(key.substr(slash+1),&pc) || pc < 1 )
            { if ( error ) *error = "HDSTATS descriptor " + file + ": component key must be PANEL/number"; return false; }
          const std::string panel_id = key.substr(0,slash);
          auto component = std::find_if(parsed.components.begin(),parsed.components.end(),[&](const component_t & c) {
            return c.panel_id == panel_id && c.panel_pc == pc;
          });
          if ( component == parsed.components.end() )
            { if ( error ) *error = "HDSTATS descriptor " + file + ": unknown component " + key; return false; }
          const int orientation = item.at("orientation").get<int>();
          if ( orientation != -1 && orientation != 1 )
            { if ( error ) *error = "HDSTATS descriptor " + file + ": orientation must be -1 or 1"; return false; }
          component->id = id;
          component->label = label;
          component->orientation = orientation;
          component->description = item.value("description",std::string());
        }

      auto read_descriptions = [&]( const char * name ) -> bool {
        if ( ! document.contains(name) ) return true;
        if ( ! document.at(name).is_array() ) return false;
        std::set<std::pair<std::string,std::string>> seen;
        for ( const auto & item : document.at(name) )
          {
            description_t value;
            if ( ! item.is_object() || ! json_string(item,"entity_type",&value.entity_type) ||
                 ! json_string(item,"entity_id",&value.entity_id) || ! json_string(item,"text",&value.text) ||
                 ! seen.insert({value.entity_type,value.entity_id}).second ) return false;
            parsed.descriptions.push_back(std::move(value));
          }
        return true;
      };
      if ( ! read_descriptions("descriptions") )
        { if ( error ) *error = "HDSTATS descriptor " + file + ": invalid descriptions"; return false; }

      if ( document.contains("external_associations") )
        {
          if ( ! document.at("external_associations").is_array() )
            { if ( error ) *error = "HDSTATS descriptor " + file + ": invalid external_associations"; return false; }
          for ( const auto & item : document.at("external_associations") )
            {
              external_association_t value;
              if ( ! item.is_object() || ! json_string(item,"entity_type",&value.entity_type) ||
                   ! json_string(item,"entity_id",&value.entity_id) || ! json_string(item,"factor",&value.factor) ||
                   ! item.contains("score") || ! item.at("score").is_number() || !std::isfinite(item.at("score").get<double>()) )
                { if ( error ) *error = "HDSTATS descriptor " + file + ": invalid external association"; return false; }
              value.score = item.at("score").get<double>();
              parsed.external_associations.push_back(std::move(value));
            }
        }
      auto known_entity = [&]( const std::string & type, const std::string & id ) {
        if ( type == "PANEL" ) return model->find_panel(id) != nullptr;
        if ( type == "PC" )
          return std::any_of(parsed.components.begin(),parsed.components.end(),[&](const component_t & c) { return c.id == id; });
        if ( type == "FEATURE" )
          return std::any_of(model->panels.begin(),model->panels.end(),[&](const panel_t & p) {
            return std::any_of(p.declared_features.begin(),p.declared_features.end(),[&](const std::string & f) { return f == id; });
          });
        return false;
      };
      for ( const auto & value : parsed.descriptions )
        if ( ! known_entity(value.entity_type,value.entity_id) )
          { if ( error ) *error = "HDSTATS descriptor " + file + ": description references unknown entity"; return false; }
      for ( const auto & value : parsed.external_associations )
        if ( ! known_entity(value.entity_type,value.entity_id) )
          { if ( error ) *error = "HDSTATS descriptor " + file + ": association references unknown entity"; return false; }

      auto read_correlations = [&]( const char * name, std::vector<correlation_t> * output ) -> bool {
        if ( ! document.contains(name) ) return true;
        if ( ! document.at(name).is_array() ) return false;
        std::set<std::pair<std::string,std::string>> seen;
        for ( const auto & item : document.at(name) )
          {
            correlation_t value;
            if ( ! item.is_object() || ! json_string(item,"first",&value.first) || ! json_string(item,"second",&value.second) ||
                 value.first == value.second || ! item.contains("value") || ! item.at("value").is_number() ||
                 !std::isfinite(item.at("value").get<double>()) ) return false;
            value.value = item.at("value").get<double>();
            const auto pair = value.first < value.second ? std::make_pair(value.first,value.second) : std::make_pair(value.second,value.first);
            if ( ! seen.insert(pair).second ) return false;
            output->push_back(std::move(value));
          }
        return true;
      };
      if ( ! read_correlations("pc_correlations",&parsed.pc_correlations) ||
           ! read_correlations("feature_correlations",&parsed.feature_correlations) )
        { if ( error ) *error = "HDSTATS descriptor " + file + ": invalid correlation records"; return false; }
      for ( const auto & value : parsed.pc_correlations )
        if ( ! known_entity("PC",value.first) || ! known_entity("PC",value.second) )
          { if ( error ) *error = "HDSTATS descriptor " + file + ": PC correlation references unknown component"; return false; }
      for ( const auto & value : parsed.feature_correlations )
        if ( ! known_entity("FEATURE",value.first) || ! known_entity("FEATURE",value.second) )
          { if ( error ) *error = "HDSTATS descriptor " + file + ": feature correlation references unknown feature"; return false; }
      *model = std::move(parsed);
      return true;
    }
  catch ( const std::exception & e )
    {
      if ( error ) *error = "HDSTATS descriptor " + file + ": " + e.what();
      return false;
    }
}

}
