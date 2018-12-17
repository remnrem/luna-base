
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


#include "graphics/graphics.h"
#include "edf/edf.h"
#include "cwt/cwt.h"
#include "spindles.h"
#include "mspindles.h"

void draw_spindles( edf_t & edf , 
		    param_t & param , 
		    const std::string & filename ,
		    int s , 
		    const std::vector<spindle_t> & spindles , 		    
		    const std::map<uint64_t, double> * avgmap )
{
  
  bool spectrogram = param.has( "heatmap" );

  spectrogram = true;

  pdf_t pdf;
  
  // hard code # of boxes per page
  const int nx = 2;
  const int ny = 4;
  
  pdf.add_page( nx , ny );
  
  pdf.set_line_width( 0.1 );
  pdf.set_grayscale_fill( 0.95 );
  pdf.set_grayscale( 0.2 );  

  const int nspindles = spindles.size();
  
  // for scale, get total min/max for CWT signal
  double avg_min , avg_max ; 
  std::vector<double> averaged;
  std::map<uint64_t,double>::const_iterator ii = avgmap->begin();
  while ( ii != avgmap->end() ) { averaged.push_back( ii->second ); ++ii; }
  MiscMath::minmax( averaged , &avg_min , &avg_max );
  
  int curr = 0;
  for (int i = 0 ; i < nspindles ; i++ )
    {

      // next box/or page
      if ( ! pdf.set_grid( curr++ ) )
       	{
      	  pdf.add_page( nx , ny );

	  pdf.set_line_width( 0.1 );
	  pdf.set_grayscale_fill( 0.05 );
	  pdf.set_grayscale( 0.2 );  

	  curr = 0;
	  pdf.set_grid( curr++ );
	}

      // bounding rectangle
      pdf.set_line_width( 0.1 );
      pdf.set_grayscale_fill( 0.95 );
      pdf.set_grayscale( 0.2 );  
      pdf.rectangle( 0, 0, 1, 1 );
      pdf.stroke_fill();

      // ID and spindle count
      pdf.set_font_color( "black" );
      pdf.set_fontsize( 6 );
      pdf.text( 0.05,0.1, "ID: " + edf.id + " | " 
		+ globals::current_tag 
		+ " | Spindle: " + Helper::int2str( i+1 ) 
		+ " | Q = " + Helper::dbl2str( spindles[i].qual ).substr(0,4) );
    
      
      //
      // Overall time-line for this spindle?
      //
      
      interval_t interval = spindles[i].tp;
      uint64_t   mx = edf.timeline.last_time_point_tp;

      double pos = interval.mid();      
      double position = pos / (double)mx;

      double st1 = interval.start;
      double st2 = interval.stop;
      
      pdf.set_color( "black" );
      pdf.set_line_width( 0.2 );      
      pdf.move( 0.05 , 0.15 );
      pdf.line( 0.95 , 0.15 );
      pdf.stroke();
      pdf.move( 0.05 , 0.13 );
      pdf.line( 0.05 , 0.17 );
      pdf.stroke();
      pdf.move( 0.95 , 0.13 );
      pdf.line( 0.95 , 0.17 );
      pdf.stroke();

      pdf.set_color( "silver" );
      pdf.set_line_width( 2 );      
      pdf.move( 0.05 + 0.9 * position , 0.13 );
      pdf.line( 0.05 + 0.9 * position , 0.17 );
      pdf.stroke();

           
      //
      // Extract raw data
      //
      
      const double flanking = 5; // 5seconds either side
      const uint64_t flanking_msec = flanking * globals::tp_1sec;
      
      uint64_t spindle_start = interval.start;
      uint64_t spindle_stop  = interval.stop;

      // get signal up to 5sec around the spindle
      interval.expand( flanking_msec );
      
      //      std::cout << "looking for " << interval.start << " .. " << interval.stop << "\n";
      
      slice_t slice( edf , s , interval );
      
      std::vector<double> d = *slice.pdata();
      std::vector<uint64_t> tp = *slice.ptimepoints();
      
      // what we actually get may be less than 
      // asked for (i.e. discontinuous EDF, edges), so track
      
      interval_t actual_slice = slice.duration();
      double actual_sec = actual_slice.duration() * globals::tp_duration;

      pdf.set_font_color( "black" );
      pdf.set_fontsize( 6 );
      pdf.text( 0.75,0.1, "window(sec): " + Helper::dbl2str( actual_sec ) );
      
      // location in time
      
      std::string timestr = Helper::timestring( edf.header.starttime , actual_slice );

//       std::stringstream ss;
//       ss << " ( " << spindle_start << " - " << spindle_stop << " )";
//       timestr.append( ss.str() );
      
      pdf.set_font_color( "black" );
      pdf.set_fontsize( 8 );
      pdf.text( 0.45,0.1, timestr );
      
      const int npoints = d.size();

      // signal min/max 
      double min, max;
      MiscMath::minmax( d , &min , &max );
      
      double plot_left = 0.1;
      double plot_right = 0.9;
      double plot_xinc = ( plot_right - plot_left ) / (double)npoints ;

      double plot_top = 0.3;
      double plot_bottom = 0.5;
      double plot_yrange = plot_bottom - plot_top;

      // scale
      pdf.set_color( "black" );
      pdf.set_line_width( 0.2 );      
      pdf.move( plot_left - 0.005 , plot_top );
      pdf.line( plot_left - 0.005 , plot_bottom );
      pdf.stroke();

      pdf.set_color( "black" );
      pdf.set_line_width( 0.05 );      

      double xpoint = plot_left;
      double ypoint = plot_bottom - plot_yrange * ( ( d[0] - min ) / ( max - min )  ) ;
      pdf.move( xpoint , ypoint );
      
      for (int j=1;j<npoints;j++)
	{
	  double xpoint = plot_left + j * plot_xinc;
	  double ypoint = plot_bottom - plot_yrange * ( ( d[j] - min ) / ( max - min )  ) ;
	  pdf.line( xpoint , ypoint );	  
	}
      pdf.stroke();
      

      //
      // Underlne spindle(s) [ at top and also bottom ]
      //
      
      pdf.set_color( "blue" );
      pdf.set_line_width( 2 );      
      
      uint64_t first_point = actual_slice.start;
      uint64_t last_point  = actual_slice.stop;
      
      int spindle_start_point = npoints * ( spindle_start - first_point ) / double( last_point - first_point );  
      int spindle_end_point = npoints * ( spindle_stop - first_point ) / double( last_point - first_point );  
      
      pdf.move( plot_left + spindle_start_point * plot_xinc , 0.23 );
      pdf.line( plot_left + spindle_end_point * plot_xinc , 0.23 );
      pdf.stroke();
      pdf.text( plot_left + spindle_start_point * plot_xinc , 0.22 , Helper::int2str(i+1) );

      if ( spectrogram )
	{
	  pdf.move( plot_left + spindle_start_point * plot_xinc , 0.93 );
	  pdf.line( plot_left + spindle_end_point * plot_xinc , 0.93 );
	  pdf.stroke();
	}
      
      // prior / after ? 
      int k = i;
      while ( 1 ) 
	{
	  --k;
	  if ( k < 0 ) break;
	  if ( spindles[k].tp.stop < first_point ) break;
	  // here, we have an inrange spindle
	  uint64_t this_start = spindles[k].tp.start;
	  uint64_t this_stop  = spindles[k].tp.stop;
	  if ( this_start < first_point ) this_start = first_point;
	  
	  int spindle_start_point = npoints * ( this_start - first_point ) / double( last_point - first_point );  
	  int spindle_end_point = npoints * ( this_stop - first_point ) / double( last_point - first_point );  	  
	  pdf.set_color( "olive" );
	  pdf.set_line_width( 2 );      
	  pdf.move( plot_left + spindle_start_point * plot_xinc , 0.25 );
	  pdf.line( plot_left + spindle_end_point * plot_xinc , 0.25 );
	  pdf.stroke();	  
	  pdf.text( plot_left + spindle_start_point * plot_xinc , 0.22 , Helper::int2str(k+1) );
	  if ( spectrogram )
	    {
	      pdf.move( plot_left + spindle_start_point * plot_xinc , 0.95 );
	      pdf.line( plot_left + spindle_end_point * plot_xinc , 0.95 );
	      pdf.stroke();	  
	    }
	}

      k = i;
      while ( 1 ) 
	{
	  ++k;
	  if ( k == spindles.size() ) break;
	  if ( spindles[k].tp.start > last_point ) break;
	  // here, we have an inrange spindle
	  uint64_t this_start = spindles[k].tp.start;
	  uint64_t this_stop  = spindles[k].tp.stop;
	  if ( this_start < first_point ) this_start = first_point;
	  
	  int spindle_start_point = npoints * ( this_start - first_point ) / double( last_point - first_point );  
	  int spindle_end_point = npoints * ( this_stop - first_point ) / double( last_point - first_point );  	  
	  pdf.set_color( "olive" );
	  pdf.set_line_width( 2 );      
	  pdf.move( plot_left + spindle_start_point * plot_xinc , 0.25 );
	  pdf.line( plot_left + spindle_end_point * plot_xinc , 0.25 );
	  pdf.stroke();	  
	  pdf.text( plot_left + spindle_start_point * plot_xinc , 0.22 , Helper::int2str(k+1) );
	  if ( spectrogram )
	    {
	      pdf.move( plot_left + spindle_start_point * plot_xinc , 0.95 );
	      pdf.line( plot_left + spindle_end_point * plot_xinc , 0.95 );
	      pdf.stroke();	  
	    }	  
	}
      
	
      //
      // Extract CWT signal.
      //

      // shift plot down
      plot_top = 0.6;
      plot_bottom = 0.7;
      plot_yrange = plot_bottom - plot_top;

      // scale
      pdf.set_color( "black" );
      pdf.set_line_width( 0.2 );      
      pdf.move( plot_left - 0.005 , plot_top );
      pdf.line( plot_left - 0.005 , plot_bottom );
      pdf.stroke();


      const int Fs = edf.header.sampling_freq( s );      
      const uint64_t period = 1.0/(double)Fs * globals::tp_1sec;      
      
      pdf.set_color( "green" );
      pdf.set_line_width( 0.1 );      
      
      // 5second window in sample-points
      
      std::vector<double> extract(npoints,0);      
      for (int j=0; j<npoints; j++)
	{
	  std::map<uint64_t,double>::const_iterator jj = avgmap->find( tp[j] );
	  if ( jj != avgmap->end() )
	    extract[j] = jj->second ;	  
	}

      // get range again: use full average min/max range

      // MiscMath::minmax( extract , &min , &max );      
      // std::stringstream dsp_min;
      // dsp_min.precision(2);
      // dsp_min << max / avg_max;
      // pdf.textbox( plot_left - 0.05 , plot_top , plot_left, plot_top+0.1, Helper::dbl2str( max ) , HPDF_TALIGN_RIGHT );
      // pdf.textbox( plot_left - 0.05 , plot_bottom , plot_left , plot_bottom+0.1 , Helper::dbl2str( min ) , HPDF_TALIGN_RIGHT );
      // pdf.text( plot_left - 0.05 , plot_bottom , dsp_min.str() );
      
      pdf.move( plot_left , plot_bottom - plot_yrange * ( ( extract[0] - avg_min ) / ( avg_max - avg_min )  ) );    
		
      for (int j = 1 ; j < npoints ; j++ )
	{
	  double value = extract[j];
	  pdf.line( plot_left + j * plot_xinc , plot_bottom - plot_yrange * ( ( value - avg_min ) / ( avg_max - avg_min )  ) );    		    
	}
      pdf.stroke();
      

      //
      // Finally, create a heat-map based on the full 5-20 Hz range (0.5Hz intervals)
      //
      
      plot_top = 0.75;
      plot_bottom = 0.9;
      plot_yrange = plot_bottom - plot_top;

      CWT cwt;
      
      cwt.set_sampling_rate( Fs );
      const double lwr = 8;
      const double upr = 18;
      const double inter = 0.5;
      std::vector<double> fx;
      const int num_cycles = 12; // i.e. higher freq. resolution
      for (double f=lwr; f<=upr; f += inter )
	{
	  cwt.add_wavelet( f , num_cycles );  // f( Fc , number of cycles ) 
	  fx.push_back(f);
	}
      cwt.load( &d );
      cwt.run();
      
      // either use result() or raw_result()
      
      double cwt_min = cwt.raw_result(0,0);
      double cwt_max = cwt.raw_result(0,0);
      int np = cwt.points();
      int nf = cwt.freqs();
      std::vector<std::vector<double> > hm( nf );
      for (int fi=0;fi<nf;fi++) hm[fi].resize( np );
      
      //
      // normalize heatmap
      //
      
      for (int fi=0; fi<nf; fi++)
	for (int ti=0; ti<np; ti++)
	  {
	    double c = cwt.raw_result(fi,ti);
	    if      ( c < cwt_min ) cwt_min = c;
	    else if ( c > cwt_max ) cwt_max = c;
	  }
      
      for (int fi=0; fi<nf; fi++)
	for (int ti=0; ti<np; ti++)
	  {
	    const double c = cwt.raw_result(fi,ti) ;
	    hm[fi][ti] = ( c - cwt_min ) / ( cwt_max - cwt_min );	    
	  }
      
      
      pdf.heatmap( plot_left , plot_bottom , plot_right , plot_top , hm , fx );
      
      
      // next spindle
    }
  
  
  //
  // write PDF
  //
  
  pdf.write( filename );
  
}



void draw_mspindles( edf_t & edf , 
		     param_t & param , 
		     const std::string & filename ,
		     std::vector<int> s , 
		     const std::vector<mspindle_t> & spindles )
{  
  // do something  
}

