
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

#include "interpolate.h"

#include "edf/edf.h"
#include "edf/slice.h"

#include "clocs/clocs.h"
#include "clocs/topo.h"
#include "stats/statistics.h"
#include "pwl_interp_2d_scattered.h"

#include "eval.h"
#include "db/db.h"
#include "helper/helper.h"
#include "helper/logger.h"


extern logger_t logger;

extern writer_t writer;

void dsptools::chep_based_interpolation( edf_t & edf , param_t & param )
{

  // check that data are epoched 
  if ( ! edf.timeline.epoched() ) 
    Helper::halt( "requires epoch'ed data" ) ;
  
  // we require that a chep-mask has been set (or else nothing to do)
  if ( ! edf.timeline.is_chep_mask_set() ) 
    {
      logger << "  leaving interpolate... either CHEP not set, or no bad channel/epoch pairs\n";
      return;
    }

  // load channel locations
  std::string cloc_file = param.requires( "clocs" );  
  clocs_t clocs;
  clocs.load_cart( cloc_file );
  
  // get signals, dropping any non-data channels
  std::string signal_label = param.requires( "sig" );
  signal_list_t signals = edf.header.signal_list( signal_label );  
  edf.header.drop_annots_from_signal_list( &signals );
  const int ns = signals.size();

  // no valid signals? quit
  if ( ns == 0 ) 
    { 
      logger << "  no signals to interpolate, leaving\n"; 
      return; 
    } 

   
  // check that all signals have the same sample rate
  int sr = 0;
  for (int s=0;s<signals.size();s++)
    {
      if ( sr == 0 ) sr = edf.header.sampling_freq( signals(s) ) ;
      if ( edf.header.sampling_freq( signals(s) ) != sr ) 
	Helper::halt( "requires all signals to have similar sampling rate, see RESAMPLE" );      
    }



  //
  // Step through each epoch/channel
  //

  
  int ne = edf.timeline.first_epoch();

  logger << " now interpolating " << ne << " epochs\n";

  int cnt = 0;

  // if no good signals, set epoch-level mask
  int cnt_masked = 0 , cnt_noaction = 0;
  int cnt_interpolated_epochs = 0 , cnt_interpolated_cheps = 0 ; 

  while ( 1 ) 
    {
      
      int epoch = edf.timeline.next_epoch();      

      std::cout << "eee = " << epoch << "\n";

      if ( epoch == -1 ) break;
            
      logger << ".";

      if ( ++cnt % 50 == 0 ) logger << " " << cnt << " epochs\n";

      writer.epoch( edf.timeline.display_epoch( epoch ) );

      interval_t interval = edf.timeline.epoch( epoch );
      
      //
      // get signal data 
      //

      mslice_t mslice( edf , signals , interval );
      
      Data::Matrix<double> D = mslice.extract();



      //
      // Get good/bad channel lists from chep
      //
      
      signal_list_t good_signals = edf.timeline.unmasked_channels_sl( epoch , signals );

      signal_list_t bad_signals = edf.timeline.masked_channels_sl( epoch , signals );

      // index good signals, within context of 'signals' i.e. of 0..ns-1  where ns = signals.size()
      // (i.e. which might be smaller than total number of EDF signals)
      
      std::vector<int> good_signals_idx;
      for (int s=0;s<signals.size();s++) 
	if ( ! edf.timeline.masked( epoch , signals(s) ) ) 
	  good_signals_idx.push_back( s ); // i.e. different encoding, relative to signals()

      std::cerr << "epoch e " << epoch << " " << good_signals.size() << " " << bad_signals.size() << "\n";
      

      // nothing to do
      if ( bad_signals.size() == 0 ) 
	{
	  ++cnt_noaction;
	  continue;
	}

      // hopeless case: set epoch mask
      if ( good_signals.size() == 0 ) 
	{
	  int mc = edf.timeline.set_epoch_mask( epoch );
	  if ( mc == 1 ) ++cnt_masked;
	  continue;
	}

      
      //
      // Construct interpolation matrices
      //

      Data::Matrix<double> invG;

      Data::Matrix<double> Gi;
      
      clocs.make_interpolation_matrices( good_signals , bad_signals , &invG , &Gi );

      
      //
      // interpolate
      //
      
      Data::Matrix<double> I = clocs.interpolate( D , good_signals_idx , invG, Gi );


      //
      // Update original signal
      //
      
      int a , b ; 
      
      if ( ! edf.timeline.epoch_records( epoch , &a , &b ) )
	Helper::halt( "internal error in interpolate()... are non-overlappnig epochs correctly set?" );
	   
      std::cerr << " epoch UPDTE = " << epoch << " --> " << a << " " << b << "\n";

      for (int s=0;s<bad_signals.size();s++)	
	edf.update_records( a , b , bad_signals(s) , I.col(s).data_pointer() );
      
      
      cnt_interpolated_epochs++;
      cnt_interpolated_cheps += bad_signals.size();

      //
      // Next epoch...
      //
    }
  
  logger << " all done\n";

  logger << " set mask for " << cnt_masked << " epochs without any good channels\n";
  logger << " skipped " << cnt_noaction << " epochs without any bad channels\n";
  logger << " interpolated " << cnt_interpolated_epochs << " epochs, for " << cnt_interpolated_cheps << " ch/epoch pairs\n";

}



void dsptools::leave_one_out( edf_t & edf , param_t & param )
{
  
  // expect following input / args
  //   clocs file 
  //   

  // clocs file
  std::string cloc_file = param.requires( "clocs" );
  
  clocs_t clocs;

  clocs.load_cart( cloc_file );
  
  // signal list
  std::string signal_label = param.requires( "sig" );
  signal_list_t signals = edf.header.signal_list( signal_label );  

  // drop any non-data channels (modifies 'signals')
  edf.header.drop_annots_from_signal_list( &signals );

  const int ns = signals.size();
  if ( ns == 0 ) return;

  // that data are epoched 
  bool epoch = edf.timeline.epoched();
  if ( ! epoch ) Helper::halt( "requires epoch'ed data" ) ;
  
  // check all same SR
  int sr = 0;
  for (int s=0;s<signals.size();s++)
    {
      if ( sr == 0 ) sr = edf.header.sampling_freq( signals(s) ) ;
      if ( edf.header.sampling_freq( signals(s) ) != sr ) 
	Helper::halt( "requires all signals to have similar sampling rate" );      
    }


  //
  // for each channel, assume it is bad, and calculate G and Gi based on all other channels
  //
  
  std::vector<Data::Matrix<double> > invG;
  std::vector<Data::Matrix<double> > Gi;
  std::vector<std::vector<int> > good_channels;
  
  logger << " generating leave-one-out G matrices for " << signals.size() << " signals\n";
  
  for (int s=0;s<ns;s++)
    {

      signal_list_t good_signals;
      signal_list_t bad_signals;
      
      // generate 'good channel' list (i.e. all other than 's'
      std::vector<int> gc;
      for (int s2=0;s2<ns;s2++) 
	{
	  if ( s != s2 ) 
	    {
	      gc.push_back( s2 );
	      good_signals.add( signals(s2) , signals.label(s2) );
	    }
	  else
	    {
	      bad_signals.add( signals(s2) , signals.label(s2) );
	    }
	}
      
      //
      // make matrices
      //

      Data::Matrix<double> _invG;
      Data::Matrix<double> _Gi;
      
      clocs.make_interpolation_matrices( good_signals , bad_signals , &_invG , &_Gi );
      
      //
      // save all
      //

      invG.push_back( _invG );
      Gi.push_back( _Gi );
      good_channels.push_back( gc );      
      
    }
  

  //
  // Step through each epoch/channel
  //

  
  int ne = edf.timeline.first_epoch();

  logger << " now iterating through " << ne << " epochs\n";
	
  while ( 1 ) 
    {
      
      int epoch = edf.timeline.next_epoch();      
       
      if ( epoch == -1 ) break;
      
      writer.epoch( edf.timeline.display_epoch( epoch ) );

      interval_t interval = edf.timeline.epoch( epoch );
      
      mslice_t mslice( edf , signals , interval );
      
      Data::Matrix<double> D = mslice.extract();
 
      //
      // interpole for each channel
      //
      
      for (int s=0;s<ns;s++)
	{

	  Data::Matrix<double> & _invG = invG[s];
	  Data::Matrix<double> & _Gi = Gi[s];
	  std::vector<int> & _good_channels = good_channels[s];

	  // interpolate
	  Data::Matrix<double> I = clocs.interpolate( D , _good_channels , _invG, _Gi );
	  
	  // calculate error 
// 	  logger << "X " << s << "\t" 
// 		    << signals.label(s) << "\t" ;
	  
// 	  double error = 0;
	  
// 	  const int nr = I.dim1();

// 	  for (int i=0;i<nr;i++)
// 	    {
// // 	      std::cout << s << "\t"
// // 			<< epoch << "\t"
// // 			<< I[i][0] << "\t"
// // 			<< D[i][s] << "\n";

// 	      error += ( I[i][0] - D[i][s] ) * ( I[i][0] - D[i][s] ) ;
// 	    }
	  // normalize
	  //	  error /= (double)nr;


	  
	  // correlation
	  double r = Statistics::correlation( *I.col_pointer(0)->data_pointer() , 
					      *D.col_pointer(s)->data_pointer() );

	  
	  // write
	  writer.level( signals.label(s) , globals::signal_strat );
	  writer.value( "R" , r ); 
	  
	}
      writer.unlevel( globals::signal_strat );
            
    }
    writer.unepoch();

}









Data::Matrix<double> dsptools::interpolate2D( const std::vector<double> & x , 
					      const std::vector<double> & y , 
					      const std::vector<double> & z , // values 
					      const double xmin , 
					      const double xmax ,
					      const int    nx , 
					      const double ymin , 
					      const double ymax ,
					      const int    ny ) 
{

  // 2D interpolation of scattered points to a uniform 2d grid, i.e. for topoplot()
  // see: pwl_interp_2d_scattered.cpp
  // which uses code from:
  // https://people.sc.fsu.edu/~jburkardt/cpp_src/pwl_interp_2d_scattered/pwl_interp_2d_scattered.html

  // number of nodes
  const int n = x.size();

  // input (x,y) coordinates, interleaved in a 2n vector
  std::vector<double> node_xy( 2 * n );
  int j=0;
  for (int i=0;i<2*n;i+=2) { 
    node_xy[i] = x[j]; 
    node_xy[i+1] = y[j]; 
    ++j; 
  } 

  //
  //  Set up the Delaunay triangulation
  //

  int num_triangles;
  std::vector<int> triangle_node( 2 * 3 * n );  // 3n is max size
  std::vector<int> neighbor( 2 * 3 * n );
    
  r8tris2 ( n , &(node_xy[0]), num_triangles, &(triangle_node[0]), &(neighbor[0]) );

  // r8tris2 ( int          number of nodes
  //           double *     node co-ordinates (2n) in (x,y) pairs
  // output
  //           int & t      number of triangles
  //           int * [3*t]  triangile nodes
  //           int * [3*t]  triangle neighbors
  
  for ( int j = 0; j < num_triangles; j++ )
    for ( int i = 0; i < 3; i++ )
      if ( 0 < neighbor[i+j*3] )
	neighbor[i+j*3] = neighbor[i+j*3] - 1;


  // print some output (SKIP)
   triangulation_order3_print ( n , num_triangles , &(node_xy[0]), &(triangle_node[0]), &(neighbor[0]) );

  //
  // Evaluate interpolants
  //
  
  const int nxy = nx * ny;
  std::vector<double> xyi( 2 * nxy ); 
 
  double xstep = ( xmax - xmin ) / (double)nx;
  double ystep = ( ymax - ymin ) / (double)ny;
  j = 0;
  for (int xi=0;xi<nx;xi++)
    {
      double xp = xmin + xi * xstep;	
      for (int yi=0;yi<ny;yi++)
 	{
 	  xyi[ j ] = xp;
 	  xyi[ j+1 ] = ymin + yi * ystep;	  
 	  j += 2;
 	}
    }
   
 
  // copy input to honor const status  
  std::vector<double> inp = z;
 
  // perform interpolation
  
  double * zi = pwl_interp_2d_scattered_value ( n,                    // number of points
						&(node_xy[0]),        // x,y co-oords
						&(inp[0]),            // &(z[0]) 
						num_triangles,       
						&(triangle_node[0]), 
						&(neighbor[0]), 
						nxy,                  // number of points to interpolate
						&(xyi[0]) );          // co-ords for interpolation points (2*ni)

  Data::Matrix<double> Z( nx, ny );
  int k = 0;
  for ( int i = 0; i < nx; i++ )      
    for ( int j = 0; j < ny; j++ )
      Z[i][j] = zi[k++];

  delete [] zi ;

  return Z;
  
}




void dsptools::interpolate2D( topo_t * topo , const std::vector<double> & z )
{

  // alternate interface for above
  
  // input (x,y) coordinates, interleaved in a 2n vector, contain ed in topo->xy
  
  //
  //  Set up the Delaunay triangulation
  //

  int num_triangles;
  std::vector<int> triangle_node( 2 * 3 * topo->inp_n );  // 3n is max size
  std::vector<int> neighbor( 2 * 3 * topo->inp_n );

//   std::cout << "topo->inp_n = " << topo->inp_n << "\n";
//   std::cout << "topo->inp_xy " << topo->inp_xy.size() << "\n";
//   std::cout << "z = " << z.size() << "\n";

  r8tris2 ( topo->inp_n , 
 	    &(topo->inp_xy[0]), 
 	    num_triangles, 
 	    &(triangle_node[0]), 
 	    &(neighbor[0]) );
   
  for ( int j = 0; j < num_triangles; j++ )
    for ( int i = 0; i < 3; i++ )
      if ( 0 < neighbor[i+j*3] )
	neighbor[i+j*3] = neighbor[i+j*3] - 1;

//   triangulation_order3_print ( topo->inp_n , num_triangles , 
// 			       &(topo->inp_xy[0]), &(triangle_node[0]), &(neighbor[0]) );



  //
  // Evaluate interpolants
  //
  
  double * zi = pwl_interp_2d_scattered_value ( topo->inp_n, 
						&(topo->inp_xy[0]),  
						&(z[0]),          
						num_triangles,       
						&(triangle_node[0]), 
						&(neighbor[0]), 
						topo->out_n , 
						&(topo->out_xy[0]) ); 
  
   topo->out_z.resize( topo->out_n );
   for ( int k = 0; k < topo->out_n; k++ ) 
     {
       topo->out_z[k] = zi[k];
     }

  delete [] zi ;

  return;
}

