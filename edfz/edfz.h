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

#ifndef __LUNA_EDFZ_H__
#define __LUNA_EDFZ_H__

#include <cstdint>
#include <cstdio>
#include <map>
#include <string>
#include <vector>

#include "zlib-1.3.1/zlib.h"

typedef unsigned char byte_t;

struct edfz_t
{
  edfz_t();
  ~edfz_t();

  bool open_for_reading( const std::string& fn );
  bool open_for_writing( const std::string& fn );
  void close();

  size_t read( byte_t* p, const int n );
  int64_t write( byte_t* p, const int n );
  void writestring( const std::string& s, int n );
  void writestring( const int& s, int n );
  void writestring( const double& s, int n );

  bool read_record( int r, byte_t* p, const int n );
  bool read_offset( int64_t offset, byte_t* p, const int n );
  bool is_attached();
  int64_t tell();
  bool seek( int64_t offset );
  bool eof();

  void clear_index();
  void add_index( int r, int64_t offset, uint64_t tp,
                  const std::vector<std::string>& a );
  int64_t get_index( int r );
  uint64_t get_tindex( int r );
  std::vector<std::string> get_annots( int r );

  // Compatibility names; the index is now embedded in the EDFZ file.
  bool read_index();
  bool write_index( const int rs );

  void set_record_size( int rs ) { record_size = rs; }
  void set_chunk_records( int n ) { chunk_records = n > 0 ? n : 64; }

  std::string filename;
  int mode; // 0 closed, -1 read, +1 write
  int record_size;
  int chunk_records;
  std::map<int, int64_t> index;      // record -> chunk file offset
  std::map<int, uint64_t> tindex;    // record -> EDF time point
  std::map<int, std::vector<std::string>> annotation_bytes;

private:
  struct entry_t
  {
    uint64_t chunk_offset;
    uint32_t record_in_chunk;
    uint64_t timepoint;
    uint64_t annot_offset;
    uint32_t annot_length;
  };

  FILE* file;
  uint64_t header_size, nr_records, data_offset;
  uint64_t index_offset, index_size, virtual_pos;
  std::vector<entry_t> entries;
  std::vector<char> annotation_blob;
  std::vector<byte_t> header_buffer, write_buffer, read_buffer;
  uint64_t read_buffer_first_record;
  uint32_t read_buffer_records;
  bool header_complete, preamble_written, index_loaded;

  bool read_preamble();
  bool write_preamble();
  bool patch_preamble();
  bool finish_header();
  bool parse_header_properties();
  bool load_chunk( uint64_t r );
  bool flush_chunk( uint32_t nrec );
  bool read_embedded_index();
  bool ensure_index();
  bool write_embedded_index();
  bool read_exact( void* p, size_t n );
  bool write_exact( const void* p, size_t n );
  uint64_t file_size();
};

#endif
