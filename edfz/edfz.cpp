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

#include "edfz/edfz.h"

#include <algorithm>
#include <cstring>

#include "helper/helper.h"

namespace
{
const unsigned char PREAMBLE_MAGIC[8] = { 'L', 'U', 'N', 'A', 'E', 'D', 'Z', '1' };
const unsigned char CHUNK_MAGIC[4] = { 'C', 'H', 'N', 'K' };
const unsigned char INDEX_MAGIC[4] = { 'I', 'N', 'D', 'X' };
const unsigned char TRAILER_MAGIC[4] = { 'E', 'N', 'D', 'Z' };
const uint32_t FORMAT_VERSION = 3;
const uint32_t PREAMBLE_SIZE = 72;
const uint32_t CHUNK_HEADER_SIZE = 40;
const uint32_t INDEX_ENTRY_SIZE = 40;
const uint32_t TRAILER_SIZE = 32;
const uint64_t TARGET_CHUNKS = 100;
const uint64_t MIN_CHUNK_BYTES = 1ULL << 20;
const uint64_t MAX_CHUNK_BYTES = 16ULL << 20;

void put32( unsigned char* p, uint32_t x )
{
  for ( int i = 0; i < 4; ++i )
    p[i] = (unsigned char)( ( x >> ( 8 * i ) ) & 255 );
}

void put64( unsigned char* p, uint64_t x )
{
  for ( int i = 0; i < 8; ++i )
    p[i] = (unsigned char)( ( x >> ( 8 * i ) ) & 255 );
}

uint32_t get32( const unsigned char* p )
{
  uint32_t x = 0;
  for ( int i = 0; i < 4; ++i )
    x |= ( (uint32_t)p[i] ) << ( 8 * i );
  return x;
}

uint64_t get64( const unsigned char* p )
{
  uint64_t x = 0;
  for ( int i = 0; i < 8; ++i )
    x |= ( (uint64_t)p[i] ) << ( 8 * i );
  return x;
}

bool ascii_int( const byte_t* p, int n, uint64_t* v )
{
  std::string s( (const char*)p, n );
  Helper::rtrim( s );
  int64_t x = 0;
  if ( !Helper::str2signed_int64( s, &x ) || x < 0 )
    return false;
  *v = (uint64_t)x;
  return true;
}
} // namespace

edfz_t::edfz_t()
    : filename(), mode( 0 ), record_size( 0 ), chunk_records( 64 ), file( NULL ), header_size( 0 ),
      nr_records( 0 ), data_offset( 0 ), index_offset( 0 ), index_size( 0 ), virtual_pos( 0 ),
      read_buffer_first_record( 0 ), read_buffer_records( 0 ), header_complete( false ),
      preamble_written( false ), index_loaded( false )
{
}

edfz_t::~edfz_t()
{
  close();
}

bool edfz_t::write_exact( const void* p, size_t n )
{
  return file != NULL && fwrite( p, 1, n, file ) == n;
}

bool edfz_t::read_exact( void* p, size_t n )
{
  return file != NULL && fread( p, 1, n, file ) == n;
}

uint64_t edfz_t::file_size()
{
  if ( file == NULL )
    return 0;
  const long cur = ftell( file );
  fseek( file, 0, SEEK_END );
  const long end = ftell( file );
  fseek( file, cur, SEEK_SET );
  return end < 0 ? 0 : (uint64_t)end;
}

bool edfz_t::write_preamble()
{
  unsigned char p[PREAMBLE_SIZE];
  memset( p, 0, sizeof( p ) );
  memcpy( p, PREAMBLE_MAGIC, 8 );
  put32( p + 8, FORMAT_VERSION );
  put32( p + 12, PREAMBLE_SIZE );
  return write_exact( p, sizeof( p ) );
}

bool edfz_t::read_preamble()
{
  unsigned char p[PREAMBLE_SIZE];
  if ( !read_exact( p, sizeof( p ) ) || memcmp( p, PREAMBLE_MAGIC, 8 ) != 0 ||
       get32( p + 8 ) != FORMAT_VERSION || get32( p + 12 ) != PREAMBLE_SIZE )
    return false;
  header_size = get64( p + 16 );
  record_size = (int)get64( p + 24 );
  nr_records = get64( p + 32 );
  chunk_records = (int)get32( p + 40 );
  data_offset = get64( p + 44 );
  index_offset = get64( p + 52 );
  index_size = get64( p + 60 );
  return header_size >= 256 && record_size > 0 && nr_records >= 0 && chunk_records > 0 &&
         data_offset >= PREAMBLE_SIZE + header_size;
}

bool edfz_t::patch_preamble()
{
  unsigned char p[PREAMBLE_SIZE];
  memset( p, 0, sizeof( p ) );
  memcpy( p, PREAMBLE_MAGIC, 8 );
  put32( p + 8, FORMAT_VERSION );
  put32( p + 12, PREAMBLE_SIZE );
  put64( p + 16, header_size );
  put64( p + 24, (uint64_t)record_size );
  put64( p + 32, nr_records );
  put32( p + 40, (uint32_t)chunk_records );
  put64( p + 44, data_offset );
  put64( p + 52, index_offset );
  put64( p + 60, index_size );
  return fseek( file, 0, SEEK_SET ) == 0 && write_exact( p, sizeof( p ) ) &&
         fseek( file, 0, SEEK_END ) == 0;
}

bool edfz_t::open_for_writing( const std::string& fn )
{
  close();
  filename = fn;
  file = LunaIO::fopen_utf8( fn, "wb+" );
  if ( file == NULL )
    return false;
  mode = 1;
  preamble_written = write_preamble();
  return preamble_written;
}

bool edfz_t::open_for_reading( const std::string& fn )
{
  close();
  filename = fn;
  file = LunaIO::fopen_utf8( fn, "rb" );
  if ( file == NULL )
    return false;
  mode = -1;
  if ( !read_preamble() )
  {
    close();
    return false;
  }
  virtual_pos = 0;
  read_buffer.clear();
  read_buffer_records = 0;
  return fseek( file, PREAMBLE_SIZE, SEEK_SET ) == 0;
}

void edfz_t::close()
{
  if ( file == NULL )
    return;
  if ( mode == 1 )
  {
    if ( !finish_header() )
      Helper::halt( "incomplete EDFZ header" );
    while ( write_buffer.size() >= (size_t)record_size * chunk_records )
      if ( !flush_chunk( chunk_records ) )
        Helper::halt( "problem writing EDFZ chunk" );
    if ( record_size > 0 && !write_buffer.empty() )
    {
      if ( write_buffer.size() % record_size != 0 )
        Helper::halt( "incomplete EDFZ record data" );
      const uint32_t nrec = write_buffer.size() / record_size;
      if ( nrec && !flush_chunk( nrec ) )
        Helper::halt( "problem writing final EDFZ chunk" );
    }
    if ( entries.size() != nr_records )
      Helper::halt( "EDFZ record count does not match header" );
    if ( !write_embedded_index() || !patch_preamble() )
      Helper::halt( "problem finalizing EDFZ index" );
  }
  fclose( file );
  file = NULL;
  mode = 0;
}

bool edfz_t::parse_header_properties()
{
  if ( header_buffer.size() < 256 )
    return false;
  uint64_t x = 0, ns = 0;
  if ( !ascii_int( &header_buffer[184], 8, &x ) || x < 256 ||
       !ascii_int( &header_buffer[236], 8, &nr_records ) ||
       !ascii_int( &header_buffer[252], 4, &ns ) )
    return false;
  header_size = x;
  if ( header_size != 256 + 256 * ns || header_buffer.size() < header_size )
    return false;
  record_size = 0;
  const uint64_t samples_offset = 256 + 216 * ns;
  for ( uint64_t s = 0; s < ns; ++s )
  {
    uint64_t np = 0;
    if ( !ascii_int( &header_buffer[samples_offset + 8 * s], 8, &np ) )
      return false;
    record_size += (int)( 2 * np );
  }
  // Scale the chunk size to the recording while bounding the amount of data
  // inflated for a random read.  Typical EDFs therefore have about 100
  // chunks; very short files use fewer and very large files are capped at
  // roughly 16 MB per chunk.
  const uint64_t total_data_bytes = nr_records * (uint64_t)record_size;
  uint64_t target_chunk_bytes = total_data_bytes / TARGET_CHUNKS;
  if ( target_chunk_bytes < MIN_CHUNK_BYTES )
    target_chunk_bytes = MIN_CHUNK_BYTES;
  if ( target_chunk_bytes > MAX_CHUNK_BYTES )
    target_chunk_bytes = MAX_CHUNK_BYTES;
  chunk_records =
      (int)std::max<uint64_t>( 1, ( target_chunk_bytes + record_size - 1 ) / record_size );
  return record_size > 0;
}

bool edfz_t::finish_header()
{
  if ( header_complete )
    return true;
  if ( !parse_header_properties() )
    return false;
  if ( fseek( file, PREAMBLE_SIZE, SEEK_SET ) != 0 ||
       !write_exact( header_buffer.data(), header_buffer.size() ) )
    return false;
  data_offset = PREAMBLE_SIZE + header_size;
  header_complete = true;
  return fseek( file, 0, SEEK_END ) == 0;
}

int64_t edfz_t::write( byte_t* p, const int n )
{
  if ( mode != 1 || p == NULL || n < 0 )
    return -1;
  int used = 0;
  while ( used < n )
  {
    if ( !header_complete )
    {
      size_t take = (size_t)n - used;
      if ( header_buffer.size() < 256 )
        take = std::min( take, 256 - header_buffer.size() );
      else if ( header_size > header_buffer.size() )
        take = std::min( take, (size_t)( header_size - header_buffer.size() ) );
      header_buffer.insert( header_buffer.end(), p + used, p + used + take );
      used += (int)take;
      if ( header_buffer.size() >= 256 && header_size == 0 )
      {
        uint64_t x = 0;
        if ( ascii_int( &header_buffer[184], 8, &x ) )
          header_size = x;
      }
      if ( header_size > 0 && header_buffer.size() >= header_size && !finish_header() )
        return -1;
      continue;
    }
    write_buffer.insert( write_buffer.end(), p + used, p + n );
    used = n;
    while ( write_buffer.size() >= (size_t)record_size * chunk_records )
      if ( !flush_chunk( chunk_records ) )
        return -1;
  }
  return n;
}

bool edfz_t::flush_chunk( uint32_t nrec )
{
  const uLong source_len = (uLong)nrec * record_size;
  if ( nrec == 0 || record_size <= 0 || write_buffer.size() < (size_t)source_len )
    return false;
  uLongf compressed_len = compressBound( source_len );
  std::vector<byte_t> compressed( compressed_len );
  if ( compress2( compressed.data(), &compressed_len, write_buffer.data(), source_len,
                  Z_DEFAULT_COMPRESSION ) != Z_OK )
    return false;
  const uint64_t offset = (uint64_t)ftell( file ), first = entries.size();
  unsigned char h[CHUNK_HEADER_SIZE];
  memset( h, 0, sizeof( h ) );
  memcpy( h, CHUNK_MAGIC, 4 );
  put64( h + 4, first );
  put32( h + 12, nrec );
  put64( h + 16, source_len );
  put64( h + 24, compressed_len );
  put32( h + 32, crc32( 0L, write_buffer.data(), source_len ) );
  if ( !write_exact( h, sizeof( h ) ) || !write_exact( compressed.data(), compressed_len ) )
    return false;
  for ( uint32_t i = 0; i < nrec; ++i )
  {
    entry_t e;
    e.chunk_offset = offset;
    e.record_in_chunk = i;
    e.timepoint = tindex[(int)( first + i )];
    e.annot_offset = 0;
    e.annot_length = 0;
    entries.push_back( e );
    index[(int)( first + i )] = offset;
  }
  write_buffer.erase( write_buffer.begin(), write_buffer.begin() + source_len );
  return true;
}

bool edfz_t::load_chunk( uint64_t r )
{
  if ( !ensure_index() || r >= entries.size() )
    return false;
  const entry_t& e = entries[r];
  if ( read_buffer_records && r >= read_buffer_first_record &&
       r < read_buffer_first_record + read_buffer_records )
    return true;
  if ( fseek( file, (long)e.chunk_offset, SEEK_SET ) != 0 )
    return false;
  unsigned char h[CHUNK_HEADER_SIZE];
  if ( !read_exact( h, sizeof( h ) ) || memcmp( h, CHUNK_MAGIC, 4 ) != 0 )
    return false;
  const uint64_t first = get64( h + 4 ), unc = get64( h + 16 ), comp = get64( h + 24 );
  const uint32_t nrec = get32( h + 12 );
  if ( first + nrec > entries.size() || unc != (uint64_t)nrec * record_size || comp > file_size() )
    return false;
  std::vector<byte_t> compressed( comp );
  if ( !read_exact( compressed.data(), comp ) )
    return false;
  read_buffer.resize( unc );
  uLongf out = unc;
  if ( uncompress( read_buffer.data(), &out, compressed.data(), comp ) != Z_OK || out != unc ||
       crc32( 0L, read_buffer.data(), unc ) != get32( h + 32 ) )
    return false;
  read_buffer_first_record = first;
  read_buffer_records = nrec;
  return true;
}

size_t edfz_t::read( byte_t* p, const int n )
{
  if ( mode != -1 || p == NULL || n <= 0 )
    return 0;
  size_t done = 0;
  const uint64_t total = header_size + nr_records * (uint64_t)record_size;
  while ( done < (size_t)n && virtual_pos < total )
  {
    if ( virtual_pos < header_size )
    {
      const size_t take = std::min( (uint64_t)n - done, header_size - virtual_pos );
      if ( fseek( file, PREAMBLE_SIZE + virtual_pos, SEEK_SET ) != 0 ||
           !read_exact( p + done, take ) )
        break;
      virtual_pos += take;
      done += take;
      continue;
    }
    const uint64_t pos = virtual_pos - header_size, r = pos / record_size,
                   inside = pos % record_size;
    if ( !load_chunk( r ) )
      break;
    const uint64_t off = ( r - read_buffer_first_record ) * (uint64_t)record_size + inside;
    const size_t take = std::min( (uint64_t)n - done, (uint64_t)read_buffer.size() - off );
    memcpy( p + done, read_buffer.data() + off, take );
    virtual_pos += take;
    done += take;
  }
  return done;
}

bool edfz_t::read_record( int r, byte_t* p, const int n )
{
  if ( mode != -1 || r < 0 || n != record_size || p == NULL || !load_chunk( r ) )
    return false;
  memcpy( p, read_buffer.data() + ( r - read_buffer_first_record ) * (uint64_t)record_size,
          record_size );
  return true;
}

bool edfz_t::read_offset( int64_t offset, byte_t* p, const int n )
{
  return offset >= 0 && seek( offset ) && read( p, n ) == (size_t)n;
}

bool edfz_t::is_attached()
{
  return file != NULL;
}

int64_t edfz_t::tell()
{
  return (int64_t)virtual_pos;
}

bool edfz_t::seek( int64_t offset )
{
  const uint64_t total = header_size + nr_records * (uint64_t)record_size;
  if ( mode != -1 || offset < 0 || (uint64_t)offset > total )
    return false;
  virtual_pos = offset;
  read_buffer_records = 0;
  return true;
}

bool edfz_t::eof()
{
  return virtual_pos >= header_size + nr_records * (uint64_t)record_size;
}

void edfz_t::writestring( const std::string& s, int n )
{
  std::string c = s;
  c.resize( n, ' ' );
  write( (byte_t*)c.data(), n );
}

void edfz_t::writestring( const int& s, int n )
{
  writestring( Helper::int2str( s ), n );
}

void edfz_t::writestring( const double& s, int n )
{
  writestring( Helper::dbl2str_fixed( s, n ), n );
}

void edfz_t::clear_index()
{
  index.clear();
  tindex.clear();
  annotation_bytes.clear();
  entries.clear();
  annotation_blob.clear();
  index_loaded = false;
}

void edfz_t::add_index( int r, int64_t offset, uint64_t tp, const std::vector<std::string>& a )
{
  index[r] = offset;
  tindex[r] = tp;
  annotation_bytes[r] = a;
}

int64_t edfz_t::get_index( int r )
{
  if ( mode == -1 && !ensure_index() )
    return -1;
  return index.find( r ) == index.end() ? -1 : index[r];
}

uint64_t edfz_t::get_tindex( int r )
{
  if ( mode == -1 && !ensure_index() )
    return 0;
  return tindex.find( r ) == tindex.end() ? 0 : tindex[r];
}

std::vector<std::string> edfz_t::get_annots( int r )
{
  if ( mode == -1 && !ensure_index() )
    return std::vector<std::string>();
  const std::map<int, std::vector<std::string>>::const_iterator i = annotation_bytes.find( r );
  return i == annotation_bytes.end() ? std::vector<std::string>() : i->second;
}

bool edfz_t::ensure_index()
{
  return index_loaded || ( mode == -1 && read_embedded_index() );
}

bool edfz_t::read_embedded_index()
{
  if ( index_loaded )
    return true;

  const uint64_t sz = file_size();

  if ( sz < TRAILER_SIZE || fseek( file, sz - TRAILER_SIZE, SEEK_SET ) != 0 )
    return false;

  unsigned char t[TRAILER_SIZE];

  if ( !read_exact( t, sizeof( t ) ) || memcmp( t, TRAILER_MAGIC, 4 ) != 0 )
    return false;

  index_offset = get64( t + 4 );

  index_size = get64( t + 12 );

  if ( index_offset + index_size > sz - TRAILER_SIZE || index_size < 16 ||
       fseek( file, index_offset, SEEK_SET ) != 0 )
    return false;

  unsigned char ih[16];

  if ( !read_exact( ih, sizeof( ih ) ) || memcmp( ih, INDEX_MAGIC, 4 ) != 0 ||
       get32( ih + 4 ) != FORMAT_VERSION )
    return false;

  const uint64_t n = get64( ih + 8 );

  if ( n != nr_records || index_size < 16 + n * INDEX_ENTRY_SIZE )
    return false;

  entries.resize( n );

  annotation_blob.resize( index_size - 16 - n * INDEX_ENTRY_SIZE );

  for ( uint64_t i = 0; i < n; ++i )
  {
    unsigned char e[INDEX_ENTRY_SIZE];

    if ( !read_exact( e, sizeof( e ) ) )
      return false;

    entries[i].chunk_offset = get64( e );
    entries[i].record_in_chunk = get32( e + 8 );
    entries[i].timepoint = get64( e + 12 );
    entries[i].annot_offset = get64( e + 20 );
    entries[i].annot_length = get32( e + 28 );
    index[(int)i] = entries[i].chunk_offset;
    tindex[(int)i] = entries[i].timepoint;
  }

  if ( !annotation_blob.empty() && !read_exact( annotation_blob.data(), annotation_blob.size() ) )
    return false;

  for ( uint64_t i = 0; i < n; ++i )
    if ( entries[i].annot_length )
    {
      if ( entries[i].annot_offset + entries[i].annot_length > annotation_blob.size() )
        return false;

      const unsigned char* p =
          (const unsigned char*)annotation_blob.data() + entries[i].annot_offset;
      const unsigned char* end = p + entries[i].annot_length;

      if ( end - p < 4 )
        return false;

      const uint32_t na = get32( p );
      p += 4;

      std::vector<std::string> values;
      values.reserve( na );

      for ( uint32_t j = 0; j < na; ++j )
      {
        if ( end - p < 4 )
          return false;

        const uint32_t nbytes = get32( p );
        p += 4;

        if ( (uint64_t)( end - p ) < nbytes )
          return false;

        values.push_back( std::string( (const char*)p, nbytes ) );
        p += nbytes;
      }

      if ( p != end )
        return false;

      annotation_bytes[(int)i] = values;
    }
    else
      annotation_bytes[(int)i] = std::vector<std::string>();
  index_loaded = true;
  return true;
}

bool edfz_t::write_embedded_index()
{
  index_offset = (uint64_t)ftell( file );

  const uint64_t n = entries.size();

  annotation_blob.clear();

  for ( uint64_t i = 0; i < n; ++i )
  {
    entries[i].timepoint = tindex[(int)i];

    const std::vector<std::string> a = get_annots( (int)i );

    if ( !a.empty() )
    {
      entries[i].annot_offset = annotation_blob.size();

      const size_t start = annotation_blob.size();
      unsigned char count[4];

      put32( count, a.size() );
      annotation_blob.insert( annotation_blob.end(), count, count + 4 );

      for ( std::vector<std::string>::const_iterator j = a.begin(); j != a.end(); ++j )
      {
        unsigned char length[4];

        put32( length, j->size() );
        annotation_blob.insert( annotation_blob.end(), length, length + 4 );
        annotation_blob.insert( annotation_blob.end(), j->begin(), j->end() );
      }

      entries[i].annot_length = annotation_blob.size() - start;
    }
    else
    {
      entries[i].annot_offset = 0;
      entries[i].annot_length = 0;
    }
  }

  unsigned char ih[16];

  memset( ih, 0, sizeof( ih ) );
  memcpy( ih, INDEX_MAGIC, 4 );
  put32( ih + 4, FORMAT_VERSION );
  put64( ih + 8, n );

  if ( !write_exact( ih, sizeof( ih ) ) )
    return false;

  for ( uint64_t i = 0; i < n; ++i )
  {
    unsigned char e[INDEX_ENTRY_SIZE];

    memset( e, 0, sizeof( e ) );
    put64( e, entries[i].chunk_offset );
    put32( e + 8, entries[i].record_in_chunk );
    put64( e + 12, entries[i].timepoint );
    put64( e + 20, entries[i].annot_offset );
    put32( e + 28, entries[i].annot_length );

    if ( !write_exact( e, sizeof( e ) ) )
      return false;
  }

  if ( !annotation_blob.empty() && !write_exact( annotation_blob.data(), annotation_blob.size() ) )
    return false;

  index_size = (uint64_t)ftell( file ) - index_offset;

  unsigned char t[TRAILER_SIZE];

  memset( t, 0, sizeof( t ) );
  memcpy( t, TRAILER_MAGIC, 4 );
  put64( t + 4, index_offset );
  put64( t + 12, index_size );
  put32( t + 20, FORMAT_VERSION );
  return write_exact( t, sizeof( t ) );
}

bool edfz_t::read_index()
{
  return ensure_index();
}

bool edfz_t::write_index( const int rs )
{
  record_size = rs;
  return true;
}
