
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


#include "timeline/timeline.h"

#include "edf/edf.h"
#include "edf/slice.h"
#include "db/db.h"
#include "helper/logger.h"

extern writer_t writer;
extern logger_t logger;


  //
  // Channel-specific epoch masks (ch/ep mask)
  //
  
static void proc_chep( edf_t & edf , param_t & param );

bool timeline_t::is_chep_mask_set() const
{
  return chep.size() != 0;
} 

void timeline_t::clear_chep_mask()
{
  chep.clear();
} 

std::map<int,std::set<std::string> > timeline_t::make_chep_copy() const
{
  return chep;
}

void timeline_t::set_chep_mask( const int e , const std::string & s )
{
  chep[ display_epoch( e ) ].insert(s);
} 

void timeline_t::merge_chep_mask( const std::map<int,std::set<std::string> > & m ) 
{
  if ( chep.size() == 0 ) { chep = m ; return; } 
  std::map<int,std::set<std::string> >::const_iterator ii = m.begin();
  while ( ii != m.end() )
    {
      std::set<std::string>::const_iterator jj = ii->second.begin();
      while ( jj != ii->second.end() )
	{
	  chep[ ii->first ].insert( *jj );
	  ++jj;
	}
      ++ii;
    }
}


bool timeline_t::unset_chep_mask( const int e , const std::string & s ) 
{ 
  // return T if anything removed
  int e1 = display_epoch( e );    
  std::map<int,std::set<std::string> >::iterator ii = chep.find( e1 ) ;
  if ( ii == chep.end() ) return false;
  std::set<std::string>::iterator jj = ii->second.find( s );
  if ( jj == ii->second.end() ) return false;
  ii->second.erase( jj );
  return true; 
} 

bool timeline_t::masked( const int e , const std::string & s ) const 
{
  std::map<int,std::set<std::string> >::const_iterator ee = chep.find( display_epoch( e ) );
  if ( ee == chep.end() ) return false;
  return ee->second.find( s ) != ee->second.end() ;
}

// save/load cheps

signal_list_t timeline_t::masked_channels_sl( const int e0 , const signal_list_t & signals ) const
{
  const bool silent_mode = true;
  int e = display_epoch( e0 );
  //  std::cout << "e in TLM = " << e0 << " " << e << "\n";
  signal_list_t msigs;
  std::vector<std::string> m = masked_channels( e0 , signals );
  for (int i=0;i<m.size();i++) 
    {
      int chn = edf->header.signal( m[i] , silent_mode );
      if ( chn != -1 ) 
	msigs.add( chn , m[i] );
    }
  return msigs;
}

signal_list_t timeline_t::unmasked_channels_sl( const int e0 , const signal_list_t & signals ) const
{
  const bool silent_mode = true;
  signal_list_t usigs;
  int e = display_epoch( e0 );
  if ( e == -1 ) return usigs;
  std::vector<std::string> u = unmasked_channels( e0 , signals );
  for (int i=0;i<u.size();i++) 
    {
      int chn = edf->header.signal( u[i] , silent_mode );
      if ( chn != -1 ) 
	usigs.add( chn , u[i] );
    }
  return usigs;
}


std::vector<std::string> timeline_t::masked_channels( const int e0 , const signal_list_t & signals ) const
{
  int e = display_epoch( e0 );
  //  std::cerr << " e , e0 = " << e << " " << e0 << "\n";
  std::vector<std::string> m;
  const int ns = signals.size();
  bool any_masked = chep.find( e ) != chep.end() ;
  if ( ! any_masked ) return m; // all good

  const std::set<std::string> & masked = chep.find(e)->second ;
  for (int s=0; s<ns; s++) 
    {
      if ( masked.find( signals.label(s) ) != masked.end() )
	m.push_back( signals.label(s) );
    }
  return m;
}


std::vector<std::string> timeline_t::unmasked_channels( const int e0 , const signal_list_t & signals ) const
{

  int e = display_epoch( e0 );

  std::vector<std::string> u;
  const int ns = signals.size();
  bool any_masked = chep.find( e ) != chep.end() ;
  if ( ! any_masked ) 
    {
      // all good
      for (int s=0; s<ns; s++) u.push_back( signals.label(s) );
      return u; 
    }

  const std::set<std::string> & masked = chep.find(e)->second ;
  for (int s=0; s<ns; s++) 
    {
      if ( masked.find( signals.label(s) ) == masked.end() )
	u.push_back( signals.label(s) );
    }
  return u;
}

void timeline_t::collapse_chep2epoch( signal_list_t signals , const double pct , const int k )
{
  
  // drop any non-data channels (modifies 'signals')
  edf->header.drop_annots_from_signal_list( &signals );

  logger << "  masking epochs";
  if ( k ) logger << " with " << k << " or more masked channels";
  if ( pct < 1 ) logger << ( k ? ", or " : " with >" ) << pct * 100 << "% masked channels: ";

  // if k or more channels are masked --> set epoch mask  [ default 1 ] 
  // if more than pct channels are masked --> set epoch mask [ default 0 ]
  // automatically adjust main 'mask' (using set_mask(), i.e. respecting mask_mode etc)
  
  int epoch = 0; 
  int masked = 0;
  int masks_set = 0;

  std::map<int,std::set<std::string> >::iterator ee = chep.begin();
  while ( ee != chep.end() )
    {
      // **assume** same signals overlap

      int sz = ee->second.size();

      int epoch = ee->first;

      if ( ( k != 0 && sz >= k ) || 
	   ( sz / (double)signals.size() > pct ) ) 
	{
	  // change main epoch mask
	  int epoch0 = display2curr_epoch( epoch );

	  // if this epoch is still present in current file, set mask
	  if ( epoch0 != -1 ) 
	    if ( set_epoch_mask( epoch0 ) ) ++masked;
	  
	  // and also set all CHEP masks (to signals) for this epoch
	  for (int s=0;s<signals.size();s++) ee->second.insert( signals.label(s) );
	  
	}
      
      if ( mask[epoch] ) ++masks_set; 
    
      ++ee;
    }
  
  logger << masked << " epochs\n";
  
}


signal_list_t timeline_t::collapse_chep2ch( signal_list_t signals , 
					    const double pct , const int k  , 
					    bool bad_set_all_bad ,
					    bool good_set_all_good )
{

  // identify which channels (from the set 'signals') have more than 'k' (or more than pct %) of epochs masked
  // return this as a set of 'bad channels'

  // if k or more channels are masked --> set epoch mask  [ default 1 ] 
  // if more than pct channels are masked --> set epoch mask [ default 0 ]
  
  // pct is defined with the denominator as the # of data channels in 'signals'

  // drop any non-data channels (modifies 'signals')
  edf->header.drop_annots_from_signal_list( &signals );

  logger << "  masking channels";
  if ( k ) logger << " with " << k << " or more masked epochs";
  if ( pct < 1 ) logger << (k?", or " : " with > " ) << pct *100 << "% masked epochs:";

  // count of bad epochs per channel
  std::map<std::string,int> c;
  int ns = signals.size();
  int ne = num_epochs();
  for (int i=0; i<ns; i++) c[ signals.label(i) ] = 0;
  
  // get channel slots lookup-table 
  std::map<std::string,int> l2s;
  for (int i=0; i<ns; i++) 
    l2s[ signals.label(i) ] = signals(i);


  std::map<int,std::set<std::string> >::const_iterator ee = chep.begin();
  while ( ee != chep.end() )
    {
      std::set<std::string>::const_iterator ss = ee->second.begin();
      while ( ss != ee->second.end() )
	{
	  if ( c.find( *ss ) != c.end() ) c[ *ss ]++;
	  ++ss;
	}
      ++ee;
    }
  
  signal_list_t good_signals;
  signal_list_t bad_signals;

  const bool silent_mode = true;

  std::map<std::string,int>::const_iterator cc = c.begin(); 
  while ( cc != c.end() ) 
    {
      if ( l2s.find( cc->first ) != l2s.end() ) 
	{
	  if ( ! ( ( k != 0 && cc->second >= k ) || 
		   ( cc->second/(double)ne > pct ) ) )
	    good_signals.add( l2s[ cc->first ] , cc->first );
	  else
	    bad_signals.add( l2s[ cc->first ] , cc->first );
	}
      ++cc;
    }
  
 
  std::set<std::string> good_sigs;
  for (int i=0;i<good_signals.size();i++) 
    good_sigs.insert( good_signals.label(i) );
 
  // set all epochs as masked for a 'bad channel'?
  if ( bad_set_all_bad ) 
    {
      for (int i=0; i<ns; i++) 
	{
	  const std::string label = signals.label(i);
	  if ( good_sigs.find( label ) == good_sigs.end() ) 
	    {
	      logger << " " << label;
	      for (int e=0;e<ne;e++) chep[ display_epoch( e ) ].insert( label );
	    }
	}
    }
      
  // set all eoochs as unmasked for a 'good channel'?
  if ( good_set_all_good )
    {
      for (int i=0; i<ns; i++) 
	{
	  const std::string label = signals.label(i);
	  if ( good_sigs.find( label ) != good_sigs.end() ) 
	    for (int e=0;e<ne;e++) 
	      {
		int e1 = display_epoch( e );
		std::set<std::string>::iterator ii = chep[ e1 ].find( label );
		if ( ii != chep[e1].end() ) chep[ e1 ].erase( ii );
	      }
	}
    }

  logger << "\n";
  
  return bad_signals;
}


void timeline_t::dump_chep_mask( signal_list_t signals , bool write_out )
{

  const int ne = first_epoch();

  // for report to console
  int total_masked = 0;
  int total_total = 0;
  std::map<int,int> track_epochs;
  std::map<std::string,int> track_channels;
  
  // use signals list to restrict summary to channels of interest only
  edf->header.drop_annots_from_signal_list( &signals );
  const int ns = signals.size();
  
  std::map<std::string,int> chtots;
  
  while ( 1 ) 
    {
      
      int e = next_epoch_ignoring_mask();      
      
      if ( e == -1 ) break;
      
      int eptot = 0;

      interval_t interval = epoch( e );
      
      int depoch = display_epoch( e );
      
      if ( write_out )
	writer.epoch( depoch );

      if ( chep.find( depoch ) == chep.end() )
	{
	  for (int s=0;s<ns;s++)     
	    {

	      ++total_total;

	      if ( write_out )
		{
		  writer.level( signals.label(s) , globals::signal_strat );
		  writer.value( "CHEP" , false );
		}
	    }
	  
	  if ( write_out )
	    writer.unlevel( globals::signal_strat );
	}
      else
	{

	  // this epoch has 1+ channel masked 
	  
	  track_epochs[ depoch ]++;

	  const std::set<std::string> & ss = chep.find( depoch )->second;

	  for (int s=0;s<ns;s++)     
	    {
	      
	      const std::string label = signals.label(s);

	      // track total
	      ++total_total;
	      
	      bool masked = ss.find( label ) != ss.end() ;
		  
	      if ( write_out )
		{
		  writer.level( label , globals::signal_strat );
		  writer.value( "CHEP" , masked );
		}
	      
	      if ( masked ) 
		{
		  track_channels[ label ]++;
		  ++total_masked;
		  chtots[ label ]++;
		  eptot++;
		}
	      
	    }
	  
	  if ( write_out )
	    writer.unlevel( globals::signal_strat );	  
	}
      
      if ( write_out )
	writer.value( "CHEP" , eptot );
      
    } // next epoch
  
  if ( write_out )
    writer.unepoch();
  
  // ch totals
  if ( write_out )
    {
      for (int s=0;s<ns;s++)     
	{
	  writer.level( signals.label(s) , globals::signal_strat );
	  writer.value( "CHEP" , chtots[ signals.label(s) ] );
	}    
      writer.unlevel( globals::signal_strat );
    }

  //
  // report to console
  //

  int partial_masks = 0; // ignores if epoch or channels complelety masked
  int partial_masks2 = 0; // should match the above..

  int epochs_totally_masked = 0 , channels_totally_masked = 0;
  std::map<int,int>::const_iterator ii = track_epochs.begin();
  while ( ii != track_epochs.end() )
    {
      if ( ii->second == ns ) ++epochs_totally_masked;
      else partial_masks += ii->second;
      ++ii;
    }

  std::map<std::string,int>::const_iterator jj = track_channels.begin();
  while ( jj != track_channels.end() )
    {
      if ( jj->second == ne ) ++channels_totally_masked;
      else partial_masks2 += jj->second ;
      ++jj;
    }

  // ne is the current number of unmasked epochs

  // total epoch count
  const int ne_all = num_total_epochs();

  logger << "  CHEP summary:\n"
	 << "   " << total_masked << " of " << total_total << " channel/epoch pairs masked (" 
	 << round( 100 * ( total_masked / double(total_total ) ) ) << "%)\n"
	 
	 << "   " << track_epochs.size() << " of " << ne_all << " epochs with 1+ masked channel, " 
	 << epochs_totally_masked << " with all channels masked\n"

	 << "   " << track_channels.size() << " of " << ns << " channels with 1+ masked epoch, " 
	 << channels_totally_masked << " with all epochs masked\n";

  // Hmmm... need to revisit this, not clear   
  // ignore totally masked channels/epochs.. how much is left?
  // logger << "   " << partial_masks  << " (" << round( 100 * ( partial_masks / double(total_total ) ) )  << "%) channels/epochs partially masked\n";
  // logger << "   " << partial_masks2 << " (" << round( 100 * ( partial_masks2 / double(total_total ) ) )  << "%) channels/epochs partially masked\n";

}


void timeline_t::read_chep_file( const std::string & f , bool reset )
{

  if ( reset ) clear_chep_mask();
  
  if ( ! Helper::fileExists( f ) ) Helper::halt( f + " does not exist" );

  std::ifstream FIN( f.c_str() , std::ios::in );
  
  // **assumes** same epoch count, channel names
  // user's responsibility to keep that as it should be

  bool silent_mode = true;

  while ( 1 ) 
    {
      std::string ch;
      int e;      
      FIN >> e >> ch ; 
      if ( FIN.eof() ) break;
      if ( ch == "" ) break;      
      int chn = edf->header.signal( ch , silent_mode );      
      if ( chn != -1 ) chep[ e ].insert( ch );  // i.e. expecting display epoch encoding (1-based)
    }
  
  FIN.close();
}

void timeline_t::write_chep_file( const std::string & f ) const
{
  std::ofstream FOUT( f.c_str() , std::ios::out );
  if ( FOUT.bad() ) Helper::halt( "could not open " + f );
  std::map<int,std::set<std::string> >::const_iterator ee = chep.begin();
  while ( ee != chep.end() )
    {
      const std::set<std::string> & chs = ee->second;
      std::set<std::string>::const_iterator cc = chs.begin();
      while ( cc != chs.end() )
	{
	  FOUT << ee->first << "\t" 
	       << *cc << "\n";
	  //	       << edf->header.label[ *cc ] << "\n";
	  ++cc;
	}
      ++ee;
    }
  FOUT.close();
}

