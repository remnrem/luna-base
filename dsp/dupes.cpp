
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

#include "dsp/dupes.h"

#include "param.h"
#include "stats/matrix.h"
#include "edf/edf.h"
#include "edf/slice.h"

#include "helper/helper.h"
#include "helper/logger.h"
#include "db/db.h"

extern writer_t writer;
extern logger_t logger;
  
void dsptools::dupes( edf_t & edf , param_t & param )
{

  //
  // find duplicates across all signals; also, flat signals
  //
  
  signal_list_t signals = edf.header.signal_list( param.requires( "sig" ) );

  const int ns = signals.size();
  
  const bool physical_check = param.yesno( "physical" );

  // for physical checks
  const double peps = param.has( "eps" ) ? param.requires_dbl( "eps" ) : 0.01 ;  

  // to be different, at least 10% of epoch must be discordant at 'eps'
  const double pdur = param.has( "prop" ) ? param.requires_dbl( "prop" ) : 0.1 ; 

  const bool pdur0 = pdur == 0 ;   
  
  if ( peps < 0 || pdur < 0 )
    Helper::halt( "eps and prop arguments require positive values\n" );
  
  logger << "  checking signal duplicates/flat signals based on "
	 << ( physical_check ? "physical" : "digital" ) << " values\n";
  
  if ( physical_check )
    {
      logger << "  using epsilon ('eps') = " << peps << "\n";
      if ( ! pdur0 ) 
	logger << "  flagging if at least " << pdur
	       << " proportion of an epoch is discordant based on eps\n";
      else
	logger << "  flagging if any sample-point is discordant based on eps\n";
    }
  
  
  
  //
  // does this channel show variability?
  //

  std::set<int> variable;

  //
  // does this channel have illegal phys min/max ranges?
  //

  std::set<int> rangeless;  
  
  for (int s=0; s<ns; s++)
    {
      const int slot = signals(s);
      
      if ( edf.header.digital_min[ slot ] == edf.header.digital_max[ slot ] )
	rangeless.insert( s );
    }
  
  if ( physical_check )
    {
      for (int s=0; s<ns; s++)
	{
	  const int slot = signals(s);
	  
	  if ( fabs( edf.header.physical_min[ slot ] - edf.header.physical_max[ slot ] ) < peps )
	    rangeless.insert( s );
	  
	}
    }
  

  //
  // expect number of (valid) variable and divergent signals/pairs
  //
  
  const int exp_var = ns - rangeless.size();
  
  const int exp_div = exp_var > 1 ? ( exp_var * (exp_var-1) ) / 2 : 0 ; 
  
  
  //
  // does this pair (i,j) of channels show divergence (where i < j)
  //
  
  std::map<int,std::set<int> > divergent;
  
  //
  // different Fs implies divergent
  //

  for (int i=0; i<ns-1; i++)
    for (int j=i+1; j<ns; j++)
      if ( edf.header.n_samples[i] != edf.header.n_samples[j] )
	divergent[ i ].insert( j );	  
  

  //
  // Do epoch-wise (but giving up on a channel pair as needed)
  //
  
  int ne = edf.timeline.first_epoch();

  bool redundancy = true;

  int checked = 0 ; 
  
  while ( 1 )
    {
      
      int epoch = edf.timeline.next_epoch();
      
      if ( epoch == -1 ) break;
      
      interval_t interval = edf.timeline.epoch( epoch );

      ++ checked;
      
      for ( int s=0; s<ns; s++ )
	{

	  //
	  // is this a valid signal? if not, skip
	  //
	  
	  if ( rangeless.find( s ) != rangeless.end() )
	    continue;
	  
	  //
	  // get digital signal
	  //
	  
	  // 1 = no downsampling; ! physical_check --> get digital int16_t values
	  slice_t slice( edf , signals(s) , interval , 1 , ! physical_check );
	  
	  const std::vector<int16_t> * d = physical_check ? NULL : slice.ddata();

	  const std::vector<double> * p = physical_check ? slice.pdata() : NULL ;

	  const int np = physical_check ? p->size() : d->size();
	  
	  //
	  // variability?
	  //
	  
	  const bool test_flat = variable.find( s ) == variable.end();	  

	  int cnt = 0;
	  
	  if ( test_flat )
	    {
	      if ( physical_check )
		{
		  for (int i=1; i<np; i++)
		    if ( fabs( (*p)[i] - (*p)[i-1] ) > peps )
                      {

			++cnt;
			
			if ( pdur0 || ( ( cnt / (double)np ) >= pdur ) ) 
			  {
			    variable.insert( s );
			    break;
			  }
		      }
		}
	      else
		{
		  for (int i=1; i<np; i++)
		    if ( (*d)[i] != (*d)[i-1] )
		      {
			variable.insert( s );
			break;
		      }
		}
	      
	    }

	  
	  //
	  // Pairs
	  //
	  
	  for (int s2=s+1; s2<ns; s2++)
	    {
	      
	      //std::cout << "  CONSIDER " << signals.label(s) << " x " << signals.label(s2) << "\n";

	      // skip if an invalid signal
	      
	      if ( rangeless.find( s2 ) != rangeless.end() )
		continue;
	      
	      const bool test_divergent = ! ( divergent.find( s ) != divergent.end() &&
					      divergent[ s ].find( s2 ) != divergent[ s ].end() ) ;  
	      
	      if ( test_divergent )
		{
		  //std::cout << " testing " << signals.label(s) << " x " << signals.label(s2) << "\n";
		  
		  slice_t slice2( edf , signals(s2) , interval , 1 , ! physical_check );
		  
		  const std::vector<int16_t> * d2 = physical_check ? NULL : slice2.ddata();
		  
		  const std::vector<double> * p2 = physical_check ? slice2.pdata() : NULL ;
		  
		  const int np2 = physical_check ? p2->size() : d2->size();
		  
		  if ( np2 != np )
		    Helper::halt( "internal error in dupes() " );

		  int cnt = 0;
		  
		  if ( physical_check )
		    {
		      for (int i=0; i<np; i++)
			{
			  // std::cout << " i/np " << i << "/" << np << "  "
			  //  	    << " (*p)[i] " << (*p)[i] << "\t" << (*p2)[i] << "\t"
			  // 	    << peps << "\t"
			  // 	    << cnt / (double)np			    
			  // 	    << "\n" ;
			  
			  if ( fabs( (*p)[i] - (*p2)[i] ) > peps )
			    {
			      
			      ++cnt;
			      //std::cout <<"  -- diverging...\n";
			      
			      if ( pdur0 || ( cnt / (double)np ) >= pdur )
				{				  
				  //std::cout <<" DIVERGED\n";
				  divergent[ s ].insert( s2 );
				  break;
				}
			    }
			}
		    }
		  else
		    {		      
		      for (int i=0; i<np; i++)
			{
			  if ( (*d)[i] != (*d2)[i] )
			    {
			      divergent[ s ].insert( s2 );
			      break;
			    }
			}
		    }
		}
	      
	    } // next pair
	  
	} // next signal
      
      //
      // early stopping?
      //
      
      // if   a) all signals are variable
      // and  b) all pairs are divergent
      
      const int n_var = variable.size();
      int n_div = 0;
      std::map<int,std::set<int> >::const_iterator dd = divergent.begin();
      while ( dd != divergent.end() )
	{
	  n_div += dd->second.size();
	  ++dd;
	}

      //std::cout << " ep " << epoch << " " << n_var << " " << ns << " " << n_div <<" " << exp_div << "\n";
      
      // expecting (n(n-1))/2      
      
      if ( exp_div == n_div && n_var == ns ) 
	{
	  redundancy = false;
	  logger << "  no duplicates or flat signals,"
		 << " early stopping after " << checked << " epochs\n";
	  break;
	}
      
    } // next epoch
  
  
  //
  // perhaps finished early, as found all signals/pairs varied
  //
  
  if ( ! redundancy )
    {
      writer.value( "FLAT" , 0 );
      writer.value( "DUPES" , 0 );      
      return;
    }  
    
  
  int n_div = 0;
  std::map<int,std::set<int> >::const_iterator dd = divergent.begin();
  while ( dd != divergent.end() )
    {
      const std::set<int> & set2 = dd->second;      
      n_div += set2.size();      
      ++dd;
    }
  
  const bool has_flat = exp_var - variable.size() ;
  const bool has_dupes = exp_div - n_div ;

  writer.value( "INVALID" , (int)(rangeless.size() ) );
  writer.value( "FLAT" , (int)(exp_var - variable.size() ));
  writer.value( "DUPES" , exp_div - n_div );

  
  
  // which signals?

  if ( rangeless.size() )
    {
      std::set<int>::const_iterator rr = rangeless.begin();
      while ( rr != rangeless.end() )
	{
	  writer.level( signals.label(*rr) , globals::signal_strat );
	  writer.value( "INVALID" , 1 );	  
	  ++rr;
	}
      writer.unlevel( globals::signal_strat );
    }
  


  std::set<int> indupe;
  
  if ( has_dupes )
    {
      for (int s=0; s<ns-1; s++)
	{

	  if ( rangeless.find( s ) != rangeless.end() )
	    continue;
	  
	  for (int s2=s+1; s2<ns; s2++)
	    {

	      if ( rangeless.find( s2 ) != rangeless.end() )
		continue;
	      
	      if ( ! ( divergent.find( s ) != divergent.end()
		       && divergent[ s ].find( s2 ) != divergent[ s ].end() ) )
		{
		  indupe.insert( s );
		  indupe.insert( s2 );
		  writer.level( signals.label(s) + "," + signals.label(s2) , "CHS" );
		  writer.value( "DUPE" , 1 );
		}
	    }
	}
      writer.unlevel( "CHS" );
    }


  if ( has_flat || has_dupes )
    {
      for (int s=0; s<ns; s++)
	{
	  const bool is1 = variable.find( s ) == variable.end() ;
	  const bool is2 = indupe.find( s ) != indupe.end();
	  
	  if ( is1 || is2 )
	    {
	      writer.level( signals.label(s) , globals::signal_strat );	      
	      writer.value( "FLAT" , (int)is1 );
	      writer.value( "DUPE" , (int)is2 );	      
	    }
	}
      writer.unlevel( globals::signal_strat );
    }

  
  logger << "  found " << rangeless.size() << " signals with invalid (empty) ranges\n";
  
  logger << "  found " << exp_var - variable.size()  << " flat signals\n";

  if ( has_dupes )
    logger << "  found " << exp_div - n_div << " duplicated pairs, involving "
	   << indupe.size() << " unique channels\n";
  else
    logger << "  found " << exp_div - n_div << " duplicated pairs\n";
  
}


