
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


#include "cfc.h"
#include <cmath>

#include "edf/edf.h"
#include "hilbert.h"
#include "stats/matrix.h"
#include "stats/glm.h"

#include "helper/helper.h"
#include "helper/logger.h"
#include "db/db.h"
#include "eval.h"

#include "edf/slice.h"

extern writer_t writer;

extern logger_t logger;


cfc_t::cfc_t( const std::vector<double> & d , 
	      const double a1,
	      const double a2,
	      const double b1,
	      const double b2 , 
	      const double sr , 
	      const double tw ,
	      const double ripple )
  : d(d), a1(a1), a2(a2), b1(b1), b2(b2), sr(sr), tw(tw), ripple(ripple)
{ 
  if ( a2 <= a1 ) Helper::halt("cfc: invalid lower frequency band");
  if ( b2 <= b1 ) Helper::halt("cfc: invalid upper frequency band");
  if ( a2 >= b1 ) Helper::halt("cfc: invalid lower/upper frequency band combination");

}


void dsptools::cfc( edf_t & edf , param_t & param )
{
  
  //
  // Extract lower and upper frequency bands
  //

  if ( ! param.has( "a" ) ) Helper::halt( "CFC requires a=lwr,upr b=lwr,upr" );
  if ( ! param.has( "b" ) ) Helper::halt( "CFC requires a=lwr,upr b=lwr,upr" );
  std::vector<double> fa = param.dblvector("a");
  std::vector<double> fb = param.dblvector("b");
  if ( fa.size() != 2 || fb.size() !=2 ) Helper::halt( "CFC requires a=lwr,upr b=lwr,upr" );


  //
  // Get signals
  //
  
  std::string signal_label = param.requires( "sig" );
  signal_list_t signals = edf.header.signal_list( signal_label );    
  const int ns = signals.size();

  if ( ns == 0 )
    logger << "  no valid signals specified for CFC\n";

  
  //
  // Output
  //
  
  writer.var( "FA" , "Lower CFC frequency range" ); 
  writer.var( "FB" , "Upper CFC frequency range" ); 
  writer.var( "OKAY" , "Valid CFC results returned" );
  writer.var( "R2_PAC" , "Phase-amplitude coupling [0,1]" );
  writer.var( "C_AMP" , "Amplitude-amplitude coupling (AAC) correlation [-1,+1]" );
  writer.var( "Z_AMP" , "Standardized AAC correlation" );
  writer.var( "R2_TOT" , "Total CFC R-squared (phase and amplitude)" );
  

  
  //
  // using epochs or entire time-line?
  //
  
  bool epoched = param.has( "epoch" ) && edf.timeline.epoched();


  //
  // Only one set of frequencies per run of this function, but still
  // useful to track as a level
  //  
  

  const std::string level = 
    Helper::dbl2str( fa[0] ) + "-" + Helper::dbl2str( fa[1] )  
    + "x" 
    + Helper::dbl2str( fb[0] ) + "-" + Helper::dbl2str( fb[1] ) ;

  writer.level( level , "FRQS" );

     

  
  //
  // for each signal
  //

  for (int s=0;s<ns;s++)
    {
      
      logger << " glm method CFC (" << level << ") for " << signals.label(s) << "\n";
   
      //
      // Output stratifier
      //

      writer.level( signals.label(s) , globals::signal_strat );

      const int srate = edf.header.sampling_freq( signals(s) ) ; 
      
      if ( epoched ) edf.timeline.first_epoch();
      
      // either for each epoch, or for entire trace
      
      while ( 1 ) 
	{
	  
	  //
	  // fetch either an EPOCH or the entire timeline
	  //
	  
	  interval_t interval;

	  int epoch = -1;
	  
	  if ( epoched )
	    {
	      epoch = edf.timeline.next_epoch();
	      if ( epoch == -1 ) break;
	      interval = edf.timeline.epoch( epoch );
	    }
	  else
	    {
	      interval = edf.timeline.wholetrace();
	    }
	  
	  //
	  // Fetch data slice
	  //
	  
	  slice_t slice( edf , signals(s) , interval );
	  
	  const std::vector<double> * signal = slice.pdata();

	
	  //
	  // Calculate PAC
	  //
	  
	  cfc_t cfc( *signal , fa[0] , fa[1] , fb[0] , fb[1] , srate );
	  
	  bool okay = cfc.glm();
	  
	  if ( ! okay ) Helper::halt( "problem in CFC calculation" );

	  //
	  // Output
	  //

	  if ( epoched ) 
	    writer.epoch( edf.timeline.display_epoch( epoch ) );
	  
	  
	  writer.value( "OKAY" , okay );
	  writer.value( "R2_PAC" , cfc.r_PAC );
	  writer.value( "C_AMP" , cfc.c_AMP );
	  writer.value( "Z_AMP" , cfc.z_AMP );
	  writer.value( "R2_TOT" , cfc.r2_TOT );


	  //
	  // Done
	  //

	  if ( ! epoched ) break;
	  
	} // next epoch
      
      if ( epoched )
	writer.unepoch();
      
    } // next signal

  writer.unlevel( globals::signal_strat );
  
  writer.unlevel( "FRQS" );

}

bool cfc_t::glm()
{
  
  // Q: potentially should trim start and stop of windows,
  // althogh when working w/ large epochs, not really necessary

  // Step 1) filter-Hilbert signal at both bands
  
  hilbert_t ha( d , sr , a1, a2 , ripple , tw );
  hilbert_t hb( d , sr , b1, b2 , ripple , tw );
  
  // Step 2) Obtain amp(a), mod( phase(a), 2PI)  and amp(b)
  
  std::vector<double> ampa = * ha.magnitude();
  std::vector<double> pha  = * ha.phase();
  for (int i=0; i<pha.size(); i++) pha[i] = fmod( pha[i] , 2 * M_PI );
  std::vector<double> ampb = * hb.magnitude();
  
  // Step 3) Normalize 
  
  const int nrow = ampa.size();

  // DV
  ampb  = MiscMath::Z( ampb );

  // Predictors
  std::vector<double> pha_sin( nrow );
  std::vector<double> pha_cos( nrow );
  for (int i=0;i<nrow;i++) 
    {
      pha_sin[i] = sin( pha[i] );
      pha_cos[i] = cos( pha[i] );
    }
  
  pha_sin   = MiscMath::Z( pha_sin );
  pha_cos   = MiscMath::Z( pha_cos );  
  ampa      = MiscMath::Z( ampa );
  
  // Step 4) Generate predictors

  // nb. no intercept, as we performed above normalization
  Data::Matrix<double> x( nrow , 3 , 1 );
  
  for (int i = 0 ; i < nrow ; i++)
    {
      x[i][0] = pha_sin[i];
      x[i][1] = pha_cos[i];
      x[i][2] = ampa[i];
    }

  Data::Vector<double> y(ampb);
  
  //
  // Fit GLM
  //

  GLM glm( GLM::LINEAR );
  
  glm.set( y , x ); 
  
  glm.fit();
      
  bool valid = glm.valid();
 
  std::vector<bool> mask;
  Data::Vector<double> beta;
  Data::Vector<double> se;
  Data::Vector<double> lowci;
  Data::Vector<double> uprci;
  Data::Vector<double> statistic;
  Data::Vector<double> pvalue;
  
  glm.display( &beta, &se, &pvalue , &mask, &lowci, &uprci, &statistic );
  
  //    const int nterms = beta.size();
  //    for (int b = 0 ; b < nterms ; b++)
  //      {
  //        std::cout << "b" << b << "\t" 
  //                  << beta[b] << "\t" 
  //                  << se[b] << "\t" 
  //                  << pvalue[b] << "\t" 
  //                  << lowci[b] << "\t" 
  //                  << uprci[b] << "\n";
  //      }
  
  //
  // Calculate measures
  //

  // phase-ampitude coupling ( 0..1 )  r_PAC = sqrt( b1^2 + b2^2 ) ; report here R^2
  r_PAC = beta[0] * beta[0] + beta[1] * beta[1] ;
  
  // amplitude-amplitude coupling (correl -1 .. +1)
  c_AMP = beta[2];
  z_AMP = 0.5 * log( ( 1 + c_AMP ) / ( 1 - c_AMP ) );

  // total r^2
  r2_TOT = glm.calc_rsqr();
  
  
  return valid;
}
