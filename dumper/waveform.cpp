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

#include "dumper/waveform.h"

#include "annot/annotate.h"
#include "luna.h"
#include "lunapi/waveforms.h"

#include <algorithm>
#include <cerrno>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <limits>
#include <map>
#include <set>
#include <vector>
#ifdef WINDOWS
#include <direct.h>
#include <windows.h>
#else
#include <ftw.h>
#include <sys/stat.h>
#include <sys/types.h>
#endif

namespace {

std::vector<std::string> * g_recursive_lwf_files = NULL;

const char * const LWF_MAGIC = "LWF1";
const int LWF_VERSION = 2;

std::string get_lwf_dir_param( param_t & param )
{
  if ( param.has( "lwf-dir" ) ) return param.value( "lwf-dir" );
  if ( param.has( "dir" ) ) return param.value( "dir" );
  return param.requires( "lwf-dir" );
}

struct channel_info_t {
  int slot;
  std::string label;
  std::string unit;
  int sr;
  int min_samples;
  int max_samples;
};

struct sample_block_info_t {
  int n;
  double data_start_sec;
  double data_stop_sec;
  int feature_qc;
  std::vector<double> feature_values;
};

struct wave_index_t {
  std::string annot;
  std::string instance;
  std::string annot_ch;
  double annot_start_sec;
  double annot_stop_sec;
  double anchor_sec;
  double wave_start_sec;
  double wave_stop_sec;
  uint64_t payload_offset;
  std::vector<sample_block_info_t> blocks;
};

struct file_summary_t {
  std::string file;
  std::string id;
  std::string edf;
  std::string startdate;
  std::string starttime;
  std::string tag;
  std::string align;
  std::vector<std::string> def_annots;
  std::vector<std::string> feature_names;
  std::vector<channel_info_t> channels;
  int n_waves;
  std::map<std::string,int> annot_counts;
};


void write_u32( std::ofstream & O , const uint32_t x )
{
  O.write( (char*)( &x ) , sizeof(uint32_t) );
}


uint32_t read_u32( std::ifstream & I )
{
  uint32_t x = 0;
  I.read( (char*)( &x ) , sizeof(uint32_t) );
  return x;
}


void write_i32( std::ofstream & O , const int x )
{
  O.write( (char*)( &x ) , sizeof(int) );
}


int read_i32( std::ifstream & I )
{
  int x = 0;
  I.read( (char*)( &x ) , sizeof(int) );
  return x;
}


void write_f32( std::ofstream & O , const float x )
{
  O.write( (char*)( &x ) , sizeof(float) );
}


void skip_f32( std::ifstream & I , const int n )
{
  if ( n <= 0 ) return;
  I.seekg( (std::streamoff)n * (std::streamoff)sizeof(float) , std::ios::cur );
}


void write_f64( std::ofstream & O , const double x )
{
  O.write( (char*)( &x ) , sizeof(double) );
}


void write_u64( std::ofstream & O , const uint64_t x )
{
  O.write( (char*)( &x ) , sizeof(uint64_t) );
}


double read_f64( std::ifstream & I )
{
  double x = 0;
  I.read( (char*)( &x ) , sizeof(double) );
  return x;
}


uint64_t read_u64( std::ifstream & I )
{
  uint64_t x = 0LLU;
  I.read( (char*)( &x ) , sizeof(uint64_t) );
  return x;
}


void write_string( std::ofstream & O , const std::string & s )
{
  const uint32_t n = s.size();
  write_u32( O , n );
  if ( n ) O.write( s.data() , n );
}


std::string read_string( std::ifstream & I )
{
  const uint32_t n = read_u32( I );
  if ( n == 0 ) return "";
  std::string s( n , '\0' );
  I.read( &s[0] , n );
  return s;
}


std::string norm_dir( std::string dir )
{
  dir = Helper::expand( dir );
  if ( dir == "" ) return dir;
  if ( dir[ dir.size() - 1 ] != globals::folder_delimiter )
    dir += globals::folder_delimiter;
  return dir;
}


bool path_exists( const std::string & p )
{
  return Helper::fileExists( p );
}


void ensure_dir( const std::string & dir )
{
  if ( dir == "" ) Helper::halt( "empty directory path for WAVEFORM" );
  if ( path_exists( dir ) ) return;
  const std::string syscmd = globals::mkdir_command + " \"" + dir + "\"";
  const int r = system( syscmd.c_str() );
  if ( r != 0 && ! path_exists( dir ) )
    Helper::halt( "could not create directory " + dir );
}


std::string sanitize_filename_part( const std::string & s , const std::string & fallback )
{
  std::string t = Helper::sanitize( Helper::trim( s ) );
  if ( t == "" || t == "." ) return fallback;
  return t;
}


std::string join_strings( const std::vector<std::string> & x , const std::string & delim )
{
  std::stringstream ss;
  for (int i=0; i<x.size(); i++)
    {
      if ( i ) ss << delim;
      ss << x[i];
    }
  return ss.str();
}


std::vector<std::string> list_lwf_files( const std::string & dir )
{
  std::vector<std::string> files;

  DIR * dp = opendir( dir.c_str() );
  if ( dp == NULL ) Helper::halt( "could not open directory " + dir );

  struct dirent * ent;
  while ( ( ent = readdir( dp ) ) != NULL )
    {
      const std::string name = ent->d_name;
      if ( name == "." || name == ".." ) continue;
      if ( ! Helper::file_extension( name , "lwf" ) ) continue;
      files.push_back( dir + name );
    }

  closedir( dp );

  std::sort( files.begin() , files.end() );
  return files;
}


#ifndef WINDOWS

int fn_luna_lwf_ftw(const char * fpath, const struct stat *ptr, int type )
{
  if ( type == FTW_F && g_recursive_lwf_files != NULL )
    {
      const std::string file( fpath );
      if ( Helper::file_extension( file , "lwf" ) )
        g_recursive_lwf_files->push_back( file );
    }
  return 0;
}

#else

void fn_luna_lwf_win( const std::string & fpath )
{
  WIN32_FIND_DATA FindFileData;

  HANDLE hFind = FindFirstFile( (fpath+"\\*").c_str() , &FindFileData );
  if ( hFind == INVALID_HANDLE_VALUE )
    Helper::halt( "search failed on " + fpath );

  bool cont = true;
  while ( cont )
    {
      std::string fname = FindFileData.cFileName;
      if ( fname != "." && fname != ".." )
        {
          if ( FindFileData.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY )
            fn_luna_lwf_win( fpath + "\\" + fname );
          else
            {
              const std::string file = fpath + "\\" + fname;
              if ( g_recursive_lwf_files != NULL && Helper::file_extension( file , "lwf" ) )
                g_recursive_lwf_files->push_back( file );
            }
        }
      cont = FindNextFile( hFind , &FindFileData );
    }

  if ( GetLastError() != ERROR_NO_MORE_FILES )
    Helper::halt( "FindNextFile failed" );

  if ( ! FindClose(hFind) ) Helper::halt( "FindClose failed" );
}

#endif


std::vector<std::string> list_lwf_files_recursive( const std::string & dir )
{
  std::vector<std::string> files;
  g_recursive_lwf_files = &files;

#ifndef WINDOWS
  if ( ftw( Helper::expand( dir ).c_str() , fn_luna_lwf_ftw , 10 ) != 0 )
    {
      g_recursive_lwf_files = NULL;
      Helper::halt( "problem traversing folder " + dir );
    }
#else
  fn_luna_lwf_win( dir );
#endif

  g_recursive_lwf_files = NULL;
  std::sort( files.begin() , files.end() );
  return files;
}


std::vector<std::string> get_summary_dirs( param_t & param )
{
  std::vector<std::string> dirs;
  const std::vector<std::string> raw = Helper::parse( get_lwf_dir_param( param ) , "," );
  for (int i=0; i<raw.size(); i++)
    {
      const std::string dir = norm_dir( Helper::trim( raw[i] ) );
      if ( dir != "" ) dirs.push_back( dir );
    }
  if ( dirs.size() == 0 ) Helper::halt( "no valid lwf-dir= specified for --waveforms" );
  return dirs;
}


std::string next_output_file( const std::string & dir , const std::string & id , const std::string & tag )
{
  const std::string sid = sanitize_filename_part( id , "ID" );
  const std::string stag = sanitize_filename_part( tag , "waveform" );
  const std::string prefix = sid + "__" + stag + "__";

  int seq = 1;
  DIR * dp = opendir( dir.c_str() );
  if ( dp != NULL )
    {
      struct dirent * ent;
      while ( ( ent = readdir( dp ) ) != NULL )
        {
          const std::string name = ent->d_name;
          if ( name.find( prefix ) == 0 && Helper::file_extension( name , "lwf" ) )
            ++seq;
        }
      closedir( dp );
    }

  return dir + prefix + Helper::zero_pad( seq , 4 ) + ".lwf";
}


double get_anchor_sec( const interval_t & interval , const std::string & align )
{
  if ( align == "start" ) return interval.start_sec();
  if ( align == "stop" ) return interval.stop_sec();
  return interval.mid_sec();
}


std::string get_align( param_t & param )
{
  std::string align = param.has( "align" ) ? param.value( "align" ) : "mid";
  if ( align == "midpoint" ) align = "mid";
  if ( align != "start" && align != "mid" && align != "stop" )
    Helper::halt( "align must be start, mid or stop" );
  return align;
}


std::string get_require_mode( param_t & param )
{
  const std::string require = param.has( "require" ) ? param.value( "require" ) : "full";
  if ( require == "full" || require == "any" ) return require;
  Helper::halt( "require must be full or any" );
  return "full";
}


bool get_match_annot_channel( param_t & param )
{
  return param.yesno( "annot-ch-match" , true , true );
}


std::vector<channel_info_t> get_channels( edf_t & edf , param_t & param )
{
  signal_list_t signals = edf.header.signal_list( param.requires( "sig" ) , true );
  if ( signals.size() == 0 ) Helper::halt( "no signals specified for WAVEFORM" );

  std::vector<channel_info_t> channels;
  for (int s=0; s<signals.size(); s++)
    {
      if ( edf.header.is_annotation_channel( signals(s) ) ) continue;
      channel_info_t ch;
      ch.slot = signals(s);
      ch.label = signals.label(s);
      ch.unit = edf.header.phys_dimension[ signals(s) ];
      ch.sr = (int)edf.header.sampling_freq( signals(s) );
      ch.min_samples = std::numeric_limits<int>::max();
      ch.max_samples = 0;
      if ( ch.sr <= 0 )
        Helper::halt( "invalid sample rate for signal " + ch.label + " in WAVEFORM" );
      channels.push_back( ch );
    }

  if ( channels.size() == 0 ) Helper::halt( "no data channels left for WAVEFORM" );
  return channels;
}


std::vector<std::string> get_annots( edf_t & edf , param_t & param )
{
  const std::vector<std::string> requested = param.strvector( "annot" );
  if ( requested.size() == 0 ) Helper::halt( "annot= is required for WAVEFORM" );

  std::set<std::string> expanded =
    annotate_t::root_match( param.strset_xsigs( "annot" ) , edf.annotations->names() );

  std::vector<std::string> annots;
  std::set<std::string> seen;
  for (int a=0; a<requested.size(); a++)
    {
      const std::string pattern = requested[a];

      if ( pattern.size() > 0 && pattern[ pattern.size() - 1 ] == '*' )
        {
          const std::string root = pattern.substr( 0 , pattern.size() - 1 );
          bool matched = false;
          std::set<std::string>::const_iterator ee = expanded.begin();
          while ( ee != expanded.end() )
            {
              if ( ee->size() >= root.size() && ee->substr( 0 , root.size() ) == root )
                {
                  matched = true;
                  if ( seen.insert( *ee ).second ) annots.push_back( *ee );
                }
              ++ee;
            }
          if ( ! matched )
            logger << "  could not find annotation " << pattern << " - skipping\n";
          continue;
        }

      if ( expanded.find( pattern ) == expanded.end() || edf.annotations->find( pattern ) == NULL )
        {
          logger << "  could not find annotation " << pattern << " - skipping\n";
          continue;
        }

      if ( seen.insert( pattern ).second ) annots.push_back( pattern );
    }

  if ( annots.size() == 0 )
    Helper::halt( "could not find any requested annotations for WAVEFORM" );

  return annots;
}


double get_wave_left( param_t & param )
{
  double w_left = 0;
  if ( param.has( "w" ) )
    {
      const double w = param.requires_dbl( "w" );
      if ( w < 0 ) Helper::halt( "w cannot be negative for WAVEFORM" );
      w_left = w;
    }
  if ( param.has( "w-left" ) )
    {
      w_left = param.requires_dbl( "w-left" );
      if ( w_left < 0 ) Helper::halt( "w-left cannot be negative for WAVEFORM" );
    }
  return w_left;
}


double get_wave_right( param_t & param )
{
  double w_right = 0;
  if ( param.has( "w" ) )
    {
      const double w = param.requires_dbl( "w" );
      if ( w < 0 ) Helper::halt( "w cannot be negative for WAVEFORM" );
      w_right = w;
    }
  if ( param.has( "w-right" ) )
    {
      w_right = param.requires_dbl( "w-right" );
      if ( w_right < 0 ) Helper::halt( "w-right cannot be negative for WAVEFORM" );
    }
  return w_right;
}


std::vector<std::string> get_channel_labels( const std::vector<channel_info_t> & channels )
{
  std::vector<std::string> labels;
  labels.reserve( channels.size() );
  for (int c=0; c<channels.size(); c++)
    labels.push_back( channels[c].label );
  return labels;
}


waveform_feature_options_t get_feature_options( param_t & param , bool * enabled )
{
  waveform_feature_options_t options;
  options.catch_enabled = param.has( "catch22" ) || param.has( "catch24" );
  options.catch24 = param.has( "catch24" );
  options.basic_stats = param.has( "basic-stats" ) || options.catch_enabled;
  *enabled = options.basic_stats || options.catch_enabled;
  return options;
}


wave_index_t make_wave_index( const waveform_event_t & event )
{
  wave_index_t idx;
  idx.annot = event.annot;
  idx.instance = event.instance;
  idx.annot_ch = event.annot_ch;
  idx.annot_start_sec = event.annot_start_sec;
  idx.annot_stop_sec = event.annot_stop_sec;
  idx.anchor_sec = event.anchor_sec;
  idx.wave_start_sec = event.wave_start_sec;
  idx.wave_stop_sec = event.wave_stop_sec;
  idx.payload_offset = 0LLU;
  for (int c=0; c<event.blocks.size(); c++)
    {
      sample_block_info_t block;
      block.n = event.blocks[c].values.size();
      block.data_start_sec = event.blocks[c].data_start_sec;
      block.data_stop_sec = event.blocks[c].data_stop_sec;
      block.feature_qc = event.blocks[c].feature_qc;
      idx.blocks.push_back( block );
    }
  return idx;
}


uint64_t string_storage_bytes( const std::string & s )
{
  return (uint64_t)sizeof(uint32_t) + (uint64_t)s.size();
}


uint64_t wave_payload_storage_bytes( const waveform_event_t & event , const std::vector<std::string> & feature_names )
{
  uint64_t n = 0LLU;
  n += string_storage_bytes( event.annot );
  n += string_storage_bytes( event.instance );
  n += string_storage_bytes( event.annot_ch );
  n += string_storage_bytes( event.meta );
  n += 5LLU * (uint64_t)sizeof(double);
  n += (uint64_t)sizeof(int);
  for (int c=0; c<event.blocks.size(); c++)
    {
      n += (uint64_t)sizeof(int);
      n += 2LLU * (uint64_t)sizeof(double);
      if ( feature_names.size() )
        {
          n += (uint64_t)sizeof(int);
          n += (uint64_t)feature_names.size() * (uint64_t)sizeof(double);
        }
      n += (uint64_t)event.blocks[c].values.size() * (uint64_t)sizeof(float);
    }
  return n;
}


uint64_t wave_index_storage_bytes( const wave_index_t & idx )
{
  uint64_t n = 0LLU;
  n += string_storage_bytes( idx.annot );
  n += string_storage_bytes( idx.instance );
  n += string_storage_bytes( idx.annot_ch );
  n += 5LLU * (uint64_t)sizeof(double);
  n += (uint64_t)sizeof(uint64_t);
  n += (uint64_t)sizeof(int);
  n += (uint64_t)idx.blocks.size() * ( (uint64_t)sizeof(int) + 2LLU * (uint64_t)sizeof(double) );
  return n;
}


uint64_t header_storage_bytes( edf_t & edf , const std::string & outfile , const std::string & tag , const std::string & align ,
                               const std::vector<std::string> & annots , const std::vector<channel_info_t> & channels ,
                               const std::vector<std::string> & feature_names )
{
  uint64_t n = 0LLU;
  n += 4LLU;
  n += (uint64_t)sizeof(int);
  n += string_storage_bytes( edf.id );
  n += string_storage_bytes( edf.filename );
  n += string_storage_bytes( outfile );
  n += string_storage_bytes( Helper::trim( edf.header.startdate ) );
  n += string_storage_bytes( Helper::trim( edf.header.starttime ) );
  n += string_storage_bytes( tag );
  n += string_storage_bytes( align );
  n += (uint64_t)sizeof(int);
  for (int i=0; i<annots.size(); i++) n += string_storage_bytes( annots[i] );
  n += (uint64_t)sizeof(int);
  for (int i=0; i<channels.size(); i++)
    {
      n += string_storage_bytes( channels[i].label );
      n += string_storage_bytes( channels[i].unit );
      n += (uint64_t)sizeof(int);
    }
  n += (uint64_t)sizeof(int);
  for (int i=0; i<feature_names.size(); i++) n += string_storage_bytes( feature_names[i] );
  n += (uint64_t)sizeof(int);
  return n;
}


void write_header( std::ofstream & O , edf_t & edf , const std::string & outfile , const std::string & tag , const std::string & align ,
                   const std::vector<std::string> & annots , const std::vector<channel_info_t> & channels ,
                   const std::vector<std::string> & feature_names , const int n_waves )
{
  O.write( LWF_MAGIC , 4 );
  write_i32( O , LWF_VERSION );
  write_string( O , edf.id );
  write_string( O , edf.filename );
  write_string( O , outfile );
  write_string( O , Helper::trim( edf.header.startdate ) );
  write_string( O , Helper::trim( edf.header.starttime ) );
  write_string( O , tag );
  write_string( O , align );

  write_i32( O , annots.size() );
  for (int i=0; i<annots.size(); i++) write_string( O , annots[i] );

  write_i32( O , channels.size() );
  for (int i=0; i<channels.size(); i++)
    {
      write_string( O , channels[i].label );
      write_string( O , channels[i].unit );
      write_i32( O , channels[i].sr );
    }

  write_i32( O , feature_names.size() );
  for (int i=0; i<feature_names.size(); i++) write_string( O , feature_names[i] );

  write_i32( O , n_waves );
}


void write_index_entry( std::ofstream & O , const wave_index_t & idx )
{
  write_string( O , idx.annot );
  write_string( O , idx.instance );
  write_string( O , idx.annot_ch );
  write_f64( O , idx.annot_start_sec );
  write_f64( O , idx.annot_stop_sec );
  write_f64( O , idx.anchor_sec );
  write_f64( O , idx.wave_start_sec );
  write_f64( O , idx.wave_stop_sec );
  write_u64( O , idx.payload_offset );
  write_i32( O , idx.blocks.size() );
  for (int c=0; c<idx.blocks.size(); c++)
    {
      write_i32( O , idx.blocks[c].n );
      write_f64( O , idx.blocks[c].data_start_sec );
      write_f64( O , idx.blocks[c].data_stop_sec );
    }
}


void write_wave( std::ofstream & O , const waveform_event_t & event , const std::vector<std::string> & feature_names )
{
  write_string( O , event.annot );
  write_string( O , event.instance );
  write_string( O , event.annot_ch );
  write_string( O , event.meta );
  write_f64( O , event.annot_start_sec );
  write_f64( O , event.annot_stop_sec );
  write_f64( O , event.anchor_sec );
  write_f64( O , event.wave_start_sec );
  write_f64( O , event.wave_stop_sec );

  write_i32( O , event.blocks.size() );
  for (int c=0; c<event.blocks.size(); c++)
    {
      write_i32( O , event.blocks[c].values.size() );
      write_f64( O , event.blocks[c].data_start_sec );
      write_f64( O , event.blocks[c].data_stop_sec );
      if ( feature_names.size() )
        {
          write_i32( O , event.blocks[c].feature_qc );
          for (int i=0; i<feature_names.size(); i++)
            {
              std::map<std::string,double>::const_iterator ff = event.blocks[c].features.find( feature_names[i] );
              const double value = ff == event.blocks[c].features.end() ? std::numeric_limits<double>::quiet_NaN() : ff->second;
              write_f64( O , value );
            }
        }
      for (int i=0; i<event.blocks[c].values.size(); i++)
        write_f32( O , (float)event.blocks[c].values[i] );
    }
}


file_summary_t read_summary( const std::string & file )
{
  file_summary_t s;
  s.file = file;
  s.n_waves = 0;

  std::ifstream I( file.c_str() , std::ios::binary | std::ios::in );
  if ( ! I ) Helper::halt( "could not open " + file );

  char magic[4];
  I.read( magic , 4 );
  if ( I.gcount() != 4 || std::memcmp( magic , LWF_MAGIC , 4 ) != 0 )
    Helper::halt( "invalid .lwf file: " + file );

  const int version = read_i32( I );
  if ( version != 1 && version != LWF_VERSION )
    Helper::halt( "unsupported .lwf version in " + file );

  s.id = read_string( I );
  s.edf = read_string( I );
  (void)read_string( I ); // stored output filename
  s.startdate = read_string( I );
  s.starttime = read_string( I );
  s.tag = read_string( I );
  s.align = read_string( I );

  const int n_annots = read_i32( I );
  for (int i=0; i<n_annots; i++) s.def_annots.push_back( read_string( I ) );

  const int n_channels = read_i32( I );
  s.channels.resize( n_channels );
  for (int c=0; c<n_channels; c++)
    {
      s.channels[c].slot = -1;
      s.channels[c].label = read_string( I );
      s.channels[c].unit = read_string( I );
      s.channels[c].sr = read_i32( I );
      s.channels[c].min_samples = std::numeric_limits<int>::max();
      s.channels[c].max_samples = 0;
    }

  if ( version >= 2 )
    {
      const int n_features = read_i32( I );
      for (int i=0; i<n_features; i++) s.feature_names.push_back( read_string( I ) );
    }

  s.n_waves = read_i32( I );

  for (int w=0; w<s.n_waves; w++)
    {
      const std::string annot = read_string( I );
      ++s.annot_counts[ annot ];
      (void)read_string( I );
      (void)read_string( I );
      (void)read_f64( I );
      (void)read_f64( I );
      (void)read_f64( I );
      (void)read_f64( I );
      (void)read_f64( I );
      (void)read_u64( I );

      const int nwc = read_i32( I );
      if ( nwc != n_channels )
        Helper::halt( "malformed .lwf file (channel count mismatch): " + file );

      for (int c=0; c<nwc; c++)
        {
          const int n = read_i32( I );
          (void)read_f64( I );
          (void)read_f64( I );
          if ( n < s.channels[c].min_samples ) s.channels[c].min_samples = n;
          if ( n > s.channels[c].max_samples ) s.channels[c].max_samples = n;
        }
    }

  for (int c=0; c<s.channels.size(); c++)
    if ( s.channels[c].min_samples == std::numeric_limits<int>::max() )
      s.channels[c].min_samples = 0;

  return s;
}


void writer_output_for_file( const std::string & cmd , const file_summary_t & s )
{
  writer.level( s.file , "FILE" );
  writer.value( "EDF_ID" , s.id );
  writer.value( "EDF" , s.edf );
  writer.value( "TAG" , s.tag == "" ? "." : s.tag );
  writer.value( "ALIGN" , s.align );
  writer.value( "START_DATE" , s.startdate );
  writer.value( "START_TIME" , s.starttime );
  writer.value( "N_ANNOTS" , (int)s.def_annots.size() );
  writer.value( "ANNOTS" , join_strings( s.def_annots , "," ) );
  writer.value( "N_WAVES" , s.n_waves );
  writer.value( "N_CH" , (int)s.channels.size() );
  writer.value( "N_FEATURES" , (int)s.feature_names.size() );
  writer.value( "FEATURES" , s.feature_names.size() ? join_strings( s.feature_names , "," ) : "." );

  for (int c=0; c<s.channels.size(); c++)
    {
      writer.level( s.channels[c].label , "CH" );
      writer.value( "SR" , s.channels[c].sr );
      writer.value( "UNIT" , s.channels[c].unit );
      writer.value( "MIN_SAMPLES" , s.channels[c].min_samples );
      writer.value( "MAX_SAMPLES" , s.channels[c].max_samples );
      writer.unlevel( "CH" );
    }

  writer.unlevel( "FILE" );
}

} // namespace


void dsptools::waveform( edf_t & edf , param_t & param )
{
  const std::string dir = norm_dir( get_lwf_dir_param( param ) );
  const std::string tag = param.has( "tag" ) ? param.value( "tag" ) : "";
  const std::string align = get_align( param );
  const std::string require = get_require_mode( param );
  const bool match_annot_channel = get_match_annot_channel( param );
  const double w_left = get_wave_left( param );
  const double w_right = get_wave_right( param );

  ensure_dir( dir );

  std::vector<channel_info_t> channels = get_channels( edf , param );
  const std::vector<std::string> annots = get_annots( edf , param );
  const std::vector<std::string> channel_labels = get_channel_labels( channels );
  bool features_enabled = false;
  const waveform_feature_options_t feature_options = get_feature_options( param , &features_enabled );

  waveform_extract_result_t extracted =
    waveform_core::extract_annotation_window_waveforms( edf , annots , channel_labels , w_left , w_right , align , require , match_annot_channel );

  if ( extracted.total_events == 0 )
    Helper::halt( "no annotation events available for WAVEFORM" );

  logger << "  considering " << extracted.total_events << " annotation-defined waveform interval(s)\n";

  if ( features_enabled )
    waveform_core::compute_features( &extracted , feature_options );

  std::map<std::string,int>::const_iterator dd = extracted.dropped.begin();
  while ( dd != extracted.dropped.end() )
    {
      if ( dd->second > 0 )
        logger << "  dropped " << dd->second << " event(s) due to " << dd->first << "\n";
      ++dd;
    }

  for (int e=0; e<extracted.events.size(); e++)
    {
      const waveform_event_t & wave = extracted.events[e];
      for (int c=0; c<wave.blocks.size() && c<channels.size(); c++)
        {
          const int n = wave.blocks[c].values.size();
          if ( n < channels[c].min_samples ) channels[c].min_samples = n;
          if ( n > channels[c].max_samples ) channels[c].max_samples = n;
        }
    }

  if ( extracted.events.size() == 0 )
    Helper::halt( "no valid waveforms could be extracted for WAVEFORM" );

  std::vector<wave_index_t> valid_index;
  valid_index.reserve( extracted.events.size() );
  for (int i=0; i<extracted.events.size(); i++)
    valid_index.push_back( make_wave_index( extracted.events[i] ) );

  const std::string outfile = next_output_file( dir , edf.id , tag );

  logger << "  writing " << extracted.events.size() << " waveform(s) to " << outfile << "\n";

  std::ofstream O( outfile.c_str() , std::ios::binary | std::ios::out );
  if ( ! O ) Helper::halt( "could not open output file " + outfile );

  uint64_t payload_offset = header_storage_bytes( edf , outfile , tag , align , annots , channels , extracted.feature_names );
  for (int i=0; i<valid_index.size(); i++) payload_offset += wave_index_storage_bytes( valid_index[i] );
  for (int i=0; i<valid_index.size(); i++)
    {
      valid_index[i].payload_offset = payload_offset;
      payload_offset += wave_payload_storage_bytes( extracted.events[i] , extracted.feature_names );
    }

  write_header( O , edf , outfile , tag , align , annots , channels , extracted.feature_names , extracted.events.size() );
  for (int i=0; i<valid_index.size(); i++) write_index_entry( O , valid_index[i] );

  for (int i=0; i<extracted.events.size(); i++) write_wave( O , extracted.events[i] , extracted.feature_names );

  O.close();

  file_summary_t s;
  s.file = outfile;
  s.id = edf.id;
  s.edf = edf.filename;
  s.startdate = Helper::trim( edf.header.startdate );
  s.starttime = Helper::trim( edf.header.starttime );
  s.tag = tag;
  s.align = align;
  s.def_annots = annots;
  s.feature_names = extracted.feature_names;
  s.channels = channels;
  s.n_waves = extracted.events.size();
  writer_output_for_file( "WAVEFORMS" , s );
}


void dsptools::waveform_summary( param_t & param )
{
  const std::vector<std::string> dirs = get_summary_dirs( param );
  const bool recur = param.has( "recur" );
  std::vector<std::string> files;
  for (int i=0; i<dirs.size(); i++)
    {
      const std::vector<std::string> dfiles = recur ? list_lwf_files_recursive( dirs[i] ) : list_lwf_files( dirs[i] );
      files.insert( files.end() , dfiles.begin() , dfiles.end() );
    }
  if ( files.size() == 0 )
    Helper::halt( "no .lwf files found in requested folder(s)" );

  std::vector<file_summary_t> summaries;
  summaries.reserve( files.size() );

  int total_waves = 0;
  std::set<std::string> ids;
  std::set<std::string> tags;
  std::set<std::string> annots;
  std::map<std::string,std::set<int> > ch2sr;
  std::map<std::string,int> event_counts;

  for (int i=0; i<files.size(); i++)
    {
      file_summary_t s = read_summary( files[i] );
      summaries.push_back( s );
      total_waves += s.n_waves;
      ids.insert( s.id );
      if ( s.tag != "" ) tags.insert( s.tag );
      for (int a=0; a<s.def_annots.size(); a++) annots.insert( s.def_annots[a] );
      for (int c=0; c<s.channels.size(); c++) ch2sr[ s.channels[c].label ].insert( s.channels[c].sr );
      std::map<std::string,int>::const_iterator aa = s.annot_counts.begin();
      while ( aa != s.annot_counts.end() )
        {
          for (int c=0; c<s.channels.size(); c++)
            {
              const std::string key = s.id + "\t" + s.channels[c].label + "\t" + ( s.tag == "" ? "." : s.tag ) + "\t" + aa->first;
              event_counts[ key ] += aa->second;
            }
          ++aa;
        }
      writer_output_for_file( "WAVEFORMS" , s );
    }

  logger << "\n"
         << "  Summary for waveform folder(s) " << join_strings( dirs , ", " ) << "\n\n"
         << "  # files          = " << summaries.size() << "\n"
         << "  # individuals    = " << ids.size() << "\n"
         << "  # waveforms      = " << total_waves << "\n"
         << "  # unique annots  = " << annots.size() << "\n"
         << "  # unique tags    = " << tags.size() << "\n";

  if ( annots.size() )
    logger << "  annotations      = " << join_strings( std::vector<std::string>( annots.begin() , annots.end() ) , ", " ) << "\n";

  if ( tags.size() )
    logger << "  tags             = " << join_strings( std::vector<std::string>( tags.begin() , tags.end() ) , ", " ) << "\n";

  std::map<std::string,std::set<int> >::const_iterator cc = ch2sr.begin();
  while ( cc != ch2sr.end() )
    {
      if ( cc->second.size() > 1 )
        logger << "  warning: channel " << cc->first << " appears with multiple sample rates\n";
      ++cc;
    }

  std::map<std::string,int>::const_iterator ee = event_counts.begin();
  while ( ee != event_counts.end() )
    {
      std::vector<std::string> tok = Helper::parse( ee->first , "\t" );
      if ( tok.size() == 4 )
        {
          writer.level( tok[0] , "LWF_ID" );
          writer.level( tok[1] , "CH" );
          writer.level( tok[2] , "TAG" );
          writer.level( tok[3] , "ANNOT" );
          writer.value( "N" , ee->second );
          writer.unlevel( "ANNOT" );
          writer.unlevel( "TAG" );
          writer.unlevel( "CH" );
          writer.unlevel( "LWF_ID" );
        }
      ++ee;
    }

}
