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

#include "sstore.h"
#include "db/sqlwrap.h"
#include "helper/helper.h"

sstore_t::sstore_t( const std::string & f1  ) 
{
  
  std::string f = Helper::expand( f1 );

  if ( attached() ) dettach();
  
  if ( f == "-" || f == "." ) { dettach(); } 
  
  sql.open(f); 
  
  sql.synchronous(false);
  
  filename = f;

  sql.query(" CREATE TABLE IF NOT EXISTS base ("
            "   ch   VARCHAR(2) , "
            "   id   VARCHAR(8) NOT NULL , "
	    "   lvl  VARCHAR(8) , "
            "   n    INTEGER , "
	    "   val  VARCHAR(20) );" );

  sql.query(" CREATE TABLE IF NOT EXISTS epochs ("
	    "   epoch INTEGER NOT NULL , "
            "   ch   VARCHAR(2) , "
            "   id   VARCHAR(8) NOT NULL , "
	    "   lvl  VARCHAR(8) , "
            "   n    INTEGER , "
	    "   val  VARCHAR(20) ); " );

  sql.query(" CREATE TABLE IF NOT EXISTS intervals ("
	    "   start UNSIGNED BIG INT NOT NULL , "
	    "   stop  UNSIGNED BIG INT NOT NULL , "	    
            "   ch   VARCHAR(2) , "
            "   id   VARCHAR(8) NOT NULL , "
	    "   lvl  VARCHAR(8) , "
            "   n    INTEGER , "
	    "   val  VARCHAR(20) );" );

  init();
  
}

bool sstore_t::init()
{
  
  // in all cases, ch can be null
  
  // sets
  stmt_insert_base      = sql.prepare(" INSERT OR REPLACE INTO base ( ch , id , lvl , n , val ) values( :ch, :id, :lvl , :n , :val ); " );
  stmt_insert_epoch     = sql.prepare(" INSERT OR REPLACE INTO epochs ( epoch , ch , id , lvl , n , val ) values( :epoch, :ch, :id, :lvl , :n , :val ); " );
  stmt_insert_interval  = sql.prepare(" INSERT OR REPLACE INTO intervals ( start , stop , ch , id , lvl , n , val ) values( :start , :stop, :ch, :id, :lvl , :n , :val ); " );

  // gets
  stmt_fetch_base = sql.prepare( "SELECT * FROM base;" );
  
  stmt_fetch_epoch = sql.prepare( "SELECT * FROM epochs WHERE epoch == :epoch ;" );
  stmt_fetch_all_epochs = sql.prepare( "SELECT * FROM epochs ;" );

  stmt_fetch_interval = sql.prepare( "SELECT * FROM intervals WHERE start BETWEEN :a AND :b " );
  stmt_fetch_all_intervals = sql.prepare( "SELECT * FROM intervals; " );
   
  stmt_fetch_keys = sql.prepare( "SELECT id, ch, lvl , COUNT(1) FROM base GROUP BY id, ch, lvl ;" );
  stmt_fetch_keys_epochs = sql.prepare( "SELECT id, ch, lvl , COUNT(1) FROM epochs GROUP BY id, ch, lvl ;" );
  stmt_fetch_keys_intervals = sql.prepare( "SELECT id, ch, lvl , COUNT(1) FROM intervals GROUP BY id, ch, lvl ;" );

  return true;
}


bool sstore_t::release()
{
  
  sql.finalise( stmt_insert_base );
  sql.finalise( stmt_insert_epoch );
  sql.finalise( stmt_insert_interval );
  
  sql.finalise( stmt_fetch_base );
  sql.finalise( stmt_fetch_epoch );
  sql.finalise( stmt_fetch_all_epochs );
  sql.finalise( stmt_fetch_interval );
  sql.finalise( stmt_fetch_all_intervals );

  sql.finalise( stmt_fetch_keys );
  sql.finalise( stmt_fetch_keys_epochs );
  sql.finalise( stmt_fetch_keys_intervals );

  return true;
}


bool sstore_t::index()
{
  if ( ! attached() ) return false;
  
  sql.query( "CREATE INDEX IF NOT EXISTS e_idx ON epochs( epoch ); " );
  sql.query( "CREATE INDEX IF NOT EXISTS i_idx ON intervals( start , stop ); " );
  
  // schema changed, so update prepared queries
  release();
  init();  
  return true;
}


bool sstore_t::drop_index()
{
  if ( ! attached() ) return false;
  sql.query( "DROP INDEX IF EXISTS e_idx;" );
  sql.query( "DROP INDEX IF EXISTS i_idx;" );
  // schema changed, so update prepared queries
  release();
  init(); 
  return true;
}


bool sstore_t::dettach()
{
  release();
  sql.close();
  return true;
}



void sstore_t::insert_base( const std::string & id , const std::string & value , const std::string * ch , const std::string * lvl )
{

  // 0 = string
  // 1 = double
  // 1+ vector of length n
  
  sql.bind_text( stmt_insert_base , ":id" , id );

  if ( lvl == NULL )
    sql.bind_null( stmt_insert_base , ":lvl" );
  else
    sql.bind_text( stmt_insert_base , ":lvl" , *lvl );

  sql.bind_int( stmt_insert_base , ":n" ,  0 );

  sql.bind_text( stmt_insert_base , ":val" ,  value );

  if ( ch == NULL ) 
    sql.bind_null( stmt_insert_base , ":ch" );
  else
    sql.bind_text( stmt_insert_base , ":ch" , *ch );

  sql.step( stmt_insert_base );
  sql.reset( stmt_insert_base );
  
  
}

void sstore_t::insert_base( const std::string & id , const double      & value , const std::string * ch , const std::string * lvl )
{

  sql.bind_text( stmt_insert_base , ":id" , id );
  sql.bind_int( stmt_insert_base , ":n" ,  1 );  // single double
  sql.bind_double( stmt_insert_base , ":val" ,  value );

  if ( lvl == NULL ) 
    sql.bind_null( stmt_insert_base , ":lvl" );
  else
    sql.bind_text( stmt_insert_base , ":lvl" , *lvl );

  if ( ch == NULL ) 
    sql.bind_null( stmt_insert_base , ":ch" );
  else
    sql.bind_text( stmt_insert_base , ":ch" , *ch );

  sql.step( stmt_insert_base );
  sql.reset( stmt_insert_base );
  

}

void sstore_t::insert_base( const std::string & id , const std::vector<double> & value , const std::string * ch , const std::string * lvl )
{

  const int n = value.size();
  if ( n == 1 ) insert_base( id , value[0] , ch );

  sql.bind_text( stmt_insert_base , ":id" , id );
  sql.bind_int( stmt_insert_base , ":n" ,  n );  // vector of doubles

  if ( ch == NULL ) 
    sql.bind_null( stmt_insert_base , ":ch" );
  else
    sql.bind_text( stmt_insert_base , ":ch" , *ch );
  
  if ( lvl == NULL ) 
    sql.bind_null( stmt_insert_base , ":lvl" );
  else
    sql.bind_text( stmt_insert_base , ":lvl" , *lvl );

  // not portable... whole array as a blob...
  
  sqlite3_bind_blob( stmt_insert_base , 
		     sqlite3_bind_parameter_index( stmt_insert_base , ":val" ) , 
		     &(value[0]) , 
		     n * sizeof(double) , 
		     0 );
 
  sql.step( stmt_insert_base );
  sql.reset( stmt_insert_base );
  

}




sstore_data_t sstore_t::fetch_base()
{
  sstore_data_t data;
  
  while ( sql.step( stmt_fetch_base ) )
    {
      // 0  ch
      // 1  id
      // 2  lvl
      // 3  n
      // 4  value
      
      sstore_key_t key;
      sstore_value_t val;
      
      bool has_channel = ! sql.is_null( stmt_fetch_base , 0 );
      key.ch = has_channel ? sql.get_text( stmt_fetch_base , 0 ) : "" ;
      key.id = sql.get_text( stmt_fetch_base , 1 );

      bool has_level = ! sql.is_null( stmt_fetch_base , 2 );
      key.lvl = has_level ? sql.get_text( stmt_fetch_base , 2 ) : "" ; 
      
      
      int n = sql.get_int( stmt_fetch_base , 3 ) ;
      
      if ( n == 0 ) // text
	{
	  val.is_text = true;
	  val.str_value = sql.get_text( stmt_fetch_base , 4 );
	}
      else if ( n == 1 ) // double 
	{
	  val.is_double = true;
	  val.dbl_value = sql.get_double( stmt_fetch_base , 4 );	  
	}
      else // vector
	{
	  val.is_vector = true;
	  val.vec_value.resize( n );
	  
	  const double *pBuffer = reinterpret_cast<const double*>( sqlite3_column_blob( stmt_fetch_base, 4) );
	  
	  std::copy( pBuffer, pBuffer + val.vec_value.size(), &val.vec_value[0] );
	  //std::cout << " read " << val.vec_value.size() << " vec items\n";
	}
      
      // store;
      data.data[ key ] = val;
     
    }
  sql.reset( stmt_fetch_base );
  
  return data;
}



//
// Epoch and interval level inserts
//


void sstore_t::insert_epoch( const int e , const std::string & id , const std::string & value , const std::string * ch , const std::string * lvl )
{
  sql.bind_int( stmt_insert_epoch , ":epoch" ,  e );
  sql.bind_text( stmt_insert_epoch , ":id" , id );

  if ( lvl == NULL ) 
    sql.bind_null( stmt_insert_epoch , ":lvl" );
  else
    sql.bind_text( stmt_insert_epoch , ":lvl" , *lvl );

  sql.bind_int( stmt_insert_epoch , ":n" ,  0 );
  sql.bind_text( stmt_insert_epoch , ":val" ,  value );
  if ( ch == NULL ) 
    sql.bind_null( stmt_insert_epoch , ":ch" );
  else
    sql.bind_text( stmt_insert_epoch , ":ch" , *ch );

  sql.step( stmt_insert_epoch );
  sql.reset( stmt_insert_epoch );
}


void sstore_t::insert_epoch( const int e , const std::string & id , const double      & value , const std::string * ch , const std::string * lvl  )
{
  sql.bind_int( stmt_insert_epoch , ":epoch" ,  e );
  sql.bind_text( stmt_insert_epoch , ":id" , id );

  if ( lvl == NULL ) 
    sql.bind_null( stmt_insert_epoch , ":lvl" );
  else
    sql.bind_text( stmt_insert_epoch , ":lvl" , *lvl );

  sql.bind_int( stmt_insert_epoch , ":n" ,  1 );
  sql.bind_double( stmt_insert_epoch , ":val" ,  value );
  if ( ch == NULL ) 
    sql.bind_null( stmt_insert_epoch , ":ch" );
  else
    sql.bind_text( stmt_insert_epoch , ":ch" , *ch );
  sql.step( stmt_insert_epoch );
  sql.reset( stmt_insert_epoch );
}

void sstore_t::insert_epoch( const int e , const std::string & id , const std::vector<double> & value , const std::string * ch , const std::string * lvl )
{

  const int n = value.size();
  if ( n == 1 ) insert_epoch( e , id , value[0] , ch );
  
  sql.bind_int( stmt_insert_epoch , ":epoch" ,  e );  
  sql.bind_text( stmt_insert_epoch , ":id" , id );
  sql.bind_int( stmt_insert_epoch , ":n" ,  n );  // vector of doubles

  if ( lvl == NULL ) 
    sql.bind_null( stmt_insert_epoch , ":lvl" );
  else
    sql.bind_text( stmt_insert_epoch , ":lvl" , *lvl );

  if ( ch == NULL ) 
    sql.bind_null( stmt_insert_epoch , ":ch" );
  else
    sql.bind_text( stmt_insert_epoch , ":ch" , *ch );
  
  // not portable... whole array as a blob...
  
  sqlite3_bind_blob( stmt_insert_epoch , 
		     sqlite3_bind_parameter_index( stmt_insert_epoch , ":val" ) , 
		     &(value[0]) , 
		     n * sizeof(double) , 
		     0 );
 
  sql.step( stmt_insert_epoch );
  sql.reset( stmt_insert_epoch );

}

void sstore_t::insert_interval( const uint64_t a , const uint64_t b , const std::string & id , const std::string & value , const std::string * ch , const std::string * lvl )
{
  sql.bind_uint64( stmt_insert_interval , ":start" ,  a );
  sql.bind_uint64( stmt_insert_interval , ":stop" ,  b );
  sql.bind_text( stmt_insert_interval , ":id" , id );

  if ( lvl == NULL ) 
    sql.bind_null( stmt_insert_interval , ":lvl" );
  else
    sql.bind_text( stmt_insert_interval , ":lvl" , *lvl );

  sql.bind_int( stmt_insert_interval , ":n" ,  0 );
  sql.bind_text( stmt_insert_interval , ":val" ,  value );
  if ( ch == NULL ) 
    sql.bind_null( stmt_insert_interval , ":ch" );
  else
    sql.bind_text( stmt_insert_interval , ":ch" , *ch );

  sql.step( stmt_insert_interval );
  sql.reset( stmt_insert_interval );

} 

void sstore_t::insert_interval( const uint64_t a , const uint64_t b , const std::string & id , const double      & value , const std::string * ch , const std::string * lvl )
{
  sql.bind_uint64( stmt_insert_interval , ":start" ,  a );
  sql.bind_uint64( stmt_insert_interval , ":stop" ,  b );
  sql.bind_text( stmt_insert_interval , ":id" , id );

  if ( lvl == NULL ) 
    sql.bind_null( stmt_insert_interval , ":lvl" );
  else
    sql.bind_text( stmt_insert_interval , ":lvl" , *lvl );

  sql.bind_int( stmt_insert_interval , ":n" ,  1 );
  sql.bind_double( stmt_insert_interval , ":val" ,  value );
  if ( ch == NULL ) 
    sql.bind_null( stmt_insert_interval , ":ch" );
  else
    sql.bind_text( stmt_insert_interval , ":ch" , *ch );

  sql.step( stmt_insert_interval );
  sql.reset( stmt_insert_interval );

}  

void sstore_t::insert_interval( const uint64_t a , const uint64_t b , const std::string & id , const std::vector<double> & value , const std::string * ch , const std::string * lvl )
{
  const int n = value.size();
  if ( n == 1 ) insert_interval( a , b , id , value[0] , ch );
  
  sql.bind_uint64( stmt_insert_interval , ":start" ,  a );  
  sql.bind_uint64( stmt_insert_interval , ":stop" ,  b );  
  sql.bind_text( stmt_insert_interval , ":id" , id );

  if ( lvl == NULL ) 
    sql.bind_null( stmt_insert_interval , ":lvl" );
  else
    sql.bind_text( stmt_insert_interval , ":lvl" , *lvl );

  sql.bind_int( stmt_insert_interval , ":n" ,  n );  // vector of doubles

  if ( ch == NULL ) 
    sql.bind_null( stmt_insert_interval , ":ch" );
  else
    sql.bind_text( stmt_insert_interval , ":ch" , *ch );
  
  // not portable... whole array as a blob...
  
  sqlite3_bind_blob( stmt_insert_interval , 
		     sqlite3_bind_parameter_index( stmt_insert_interval , ":val" ) , 
		     &(value[0]) , 
		     n * sizeof(double) , 
		     0 );
 
  sql.step( stmt_insert_interval );
  sql.reset( stmt_insert_interval );

}



//
// Epoch and interval level fetches
//
  
sstore_data_t sstore_t::fetch_epoch( const int epoch )
{
  sstore_data_t data;
  
  // select particular epoch
  sql.bind_int( stmt_fetch_epoch , ":epoch" , epoch );
  
  while ( sql.step( stmt_fetch_epoch ) )
    {
      
      // 0  epoch
      // 1  ch
      // 2  id
      // 3  lvl 
      // 4  n
      // 5  value
      
      sstore_key_t key;
      sstore_value_t val;
      
      bool has_channel = ! sql.is_null( stmt_fetch_epoch , 1 );
      key.ch = has_channel ? sql.get_text( stmt_fetch_epoch , 1 ) : "" ;
      key.id = sql.get_text( stmt_fetch_epoch , 2 );

      bool has_level = ! sql.is_null( stmt_fetch_epoch , 3 );
      key.lvl = has_level ? sql.get_text( stmt_fetch_epoch , 3 ) : "" ;


      int n = sql.get_int( stmt_fetch_epoch , 4 ) ;
      
      if ( n == 0 ) // text
	{
	  val.is_text = true;
	  val.str_value = sql.get_text( stmt_fetch_epoch , 5 );
	}
      else if ( n == 1 ) // double 
	{
	  val.is_double = true;
	  val.dbl_value = sql.get_double( stmt_fetch_epoch , 5 );	  
	}
      else // vector
	{
	  val.is_vector = true;
	  val.vec_value.resize( n );
	  
	  const double *pBuffer = reinterpret_cast<const double*>( sqlite3_column_blob( stmt_fetch_epoch, 5) );
	  
	  std::copy( pBuffer, pBuffer + val.vec_value.size(), &val.vec_value[0] );

	}
      
      // store;
      data.data[ key ] = val;
     
    }
  sql.reset( stmt_fetch_epoch );
  
  return data;
  
}

std::map<int,sstore_data_t> sstore_t::fetch_epochs()
{
  
  std::map<int,sstore_data_t> data;

  while ( sql.step( stmt_fetch_all_epochs ) )
    {
      
      // 0  epoch
      // 1  ch
      // 2  id
      // 3  lvl
      // 4  n
      // 5  value
      
      sstore_key_t key;
      sstore_value_t val;
      
      int epoch = sql.get_int( stmt_fetch_all_epochs , 0 );

      bool has_channel = ! sql.is_null( stmt_fetch_all_epochs , 1 );
      key.ch = has_channel ? sql.get_text( stmt_fetch_all_epochs , 1 ) : "" ;
      key.id = sql.get_text( stmt_fetch_all_epochs , 2 );

      bool has_level = ! sql.is_null( stmt_fetch_all_epochs , 3 );
      key.lvl = has_level ? sql.get_text( stmt_fetch_all_epochs , 3 ) : "" ;

      int n = sql.get_int( stmt_fetch_all_epochs , 4 ) ;
      
      if ( n == 0 ) // text
	{
	  val.is_text = true;
	  val.str_value = sql.get_text( stmt_fetch_all_epochs , 5 );
	}
      else if ( n == 1 ) // double 
	{
	  val.is_double = true;
	  val.dbl_value = sql.get_double( stmt_fetch_all_epochs , 5 );	  
	}
      else // vector
	{
	  val.is_vector = true;
	  val.vec_value.resize( n );
	  
	  const double *pBuffer = reinterpret_cast<const double*>( sqlite3_column_blob( stmt_fetch_all_epochs, 5) );
	  
	  std::copy( pBuffer, pBuffer + val.vec_value.size(), &val.vec_value[0] );

	}
      
      // store;
      data[ epoch ].data[ key ] = val;
     
    }
  sql.reset( stmt_fetch_all_epochs );
  
  return data;

}

sstore_data_t sstore_t::fetch_interval( const interval_t & interval )
{
  sstore_data_t data;

  // select particular interval
  sql.bind_uint64( stmt_fetch_interval , ":start" , interval.start );
  sql.bind_uint64( stmt_fetch_interval , ":stop" , interval.stop );
  
  while ( sql.step( stmt_fetch_interval ) )
    {
      
      // 0  start
      // 1  stop
      // 2  ch
      // 3  id
      // 4  lvl
      // 5  n
      // 6  value
      
      sstore_key_t key;
      sstore_value_t val;
      
      bool has_channel = ! sql.is_null( stmt_fetch_interval , 2 );
      key.ch = has_channel ? sql.get_text( stmt_fetch_interval , 2 ) : "" ;
      key.id = sql.get_text( stmt_fetch_interval , 3 );

      bool has_level = ! sql.is_null( stmt_fetch_interval , 4 );
      key.lvl = has_level ? sql.get_text( stmt_fetch_interval , 4 ) : "" ;

      int n = sql.get_int( stmt_fetch_interval , 5 ) ;
      
      if ( n == 0 ) // text
	{
	  val.is_text = true;
	  val.str_value = sql.get_text( stmt_fetch_interval , 6 );
	}
      else if ( n == 1 ) // double 
	{
	  val.is_double = true;
	  val.dbl_value = sql.get_double( stmt_fetch_interval , 6 );	  
	}
      else // vector
	{
	  val.is_vector = true;
	  val.vec_value.resize( n );
	  
	  const double *pBuffer = reinterpret_cast<const double*>( sqlite3_column_blob( stmt_fetch_interval, 6) );
	  
	  std::copy( pBuffer, pBuffer + val.vec_value.size(), &val.vec_value[0] );

	}
      
      // store;
      data.data[ key ] = val;
     
    }
  sql.reset( stmt_fetch_interval );


  return data;

}


std::map<interval_t,sstore_data_t> sstore_t::fetch_intervals()
{

  std::map<interval_t,sstore_data_t> data;
  
  while ( sql.step( stmt_fetch_all_intervals ) )
    {
      
      // 0  start
      // 1  stop
      // 2  ch
      // 3  id
      // 4  lvl
      // 5  n
      // 6  value
      
      sstore_key_t key;
      sstore_value_t val;

      interval_t interval( sql.get_uint64( stmt_fetch_all_intervals , 0 ) , 
			   sql.get_uint64( stmt_fetch_all_intervals , 1 ) );
			         
      bool has_channel = ! sql.is_null( stmt_fetch_all_intervals , 2 );
      key.ch = has_channel ? sql.get_text( stmt_fetch_all_intervals , 2 ) : "" ;
      key.id = sql.get_text( stmt_fetch_all_intervals , 3 );
      
      bool has_level = ! sql.is_null( stmt_fetch_all_intervals , 4 );
      key.lvl = has_level ? sql.get_text( stmt_fetch_all_intervals , 4 ) : "" ;

      int n = sql.get_int( stmt_fetch_all_intervals , 5 ) ;
      
      if ( n == 0 ) // text
	{
	  val.is_text = true;
	  val.str_value = sql.get_text( stmt_fetch_all_intervals , 6 );
	}
      else if ( n == 1 ) // double 
	{
	  val.is_double = true;
	  val.dbl_value = sql.get_double( stmt_fetch_all_intervals , 6 );	  
	}
      else // vector
	{
	  val.is_vector = true;
	  val.vec_value.resize( n );
	  
	  const double *pBuffer = reinterpret_cast<const double*>( sqlite3_column_blob( stmt_fetch_all_intervals, 6) );
	  
	  std::copy( pBuffer, pBuffer + val.vec_value.size(), &val.vec_value[0] );

	}
      
      // store;
      data[ interval ].data[ key ] = val;
     
    }
  sql.reset( stmt_fetch_all_intervals );

  
  return data;
  
}


std::map<sstore_key_t,int> sstore_t::keys()
{
  std::map<sstore_key_t,int> keys;
  while ( sql.step( stmt_fetch_keys ) )
    {
      sstore_key_t key;
      key.id  = sql.get_text( stmt_fetch_keys , 0 );
      key.ch  = sql.get_text( stmt_fetch_keys , 1 );
      key.lvl = sql.get_text( stmt_fetch_keys , 2 );      
      if ( key.ch == "" ) key.ch = ".";
      if ( key.lvl == "" ) key.lvl = ".";
      keys[ key ] += sql.get_int( stmt_fetch_keys , 3 );
    }
  sql.reset( stmt_fetch_keys );
  return keys;
}

std::map<sstore_key_t,int> sstore_t::keys_epoch()
{
  std::map<sstore_key_t,int> keys;
  while ( sql.step( stmt_fetch_keys_epochs ) )
    {
      sstore_key_t key;
      key.id  = sql.get_text( stmt_fetch_keys_epochs , 0 );
      key.ch  = sql.get_text( stmt_fetch_keys_epochs , 1 );
      key.lvl = sql.get_text( stmt_fetch_keys_epochs , 2 );      
      if ( key.ch == "" ) key.ch = ".";
      if ( key.lvl == "" ) key.lvl = ".";
      keys[ key ] += sql.get_int( stmt_fetch_keys_epochs , 3 );
    }
  sql.reset( stmt_fetch_keys_epochs );
  return keys;
}

std::map<sstore_key_t,int> sstore_t::keys_interval()
{
  std::map<sstore_key_t,int> keys;
  while ( sql.step( stmt_fetch_keys_intervals ) )
    {
      sstore_key_t key;
      key.id  = sql.get_text( stmt_fetch_keys_intervals , 0 );
      key.ch  = sql.get_text( stmt_fetch_keys_intervals , 1 );
      key.lvl = sql.get_text( stmt_fetch_keys_intervals , 2 );            
      if ( key.ch == "" ) key.ch = ".";
      if ( key.lvl == "" ) key.lvl = ".";
      keys[ key ] += sql.get_int( stmt_fetch_keys_intervals , 3 );

    }
  sql.reset( stmt_fetch_keys_intervals );
  return keys;

}

