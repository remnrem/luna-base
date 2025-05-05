
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
#include "annot/annotate.h"  // for root_match()
#include "dsp/spline.h"

extern writer_t writer;
extern logger_t logger;

// implementations of various timeline-based functions - or those
// that work with combinations of annots/signals/epochs/records
// and so do not obviously fit somewhere else
//

//  S2A (two versions; orig, + waveform version)
//  A2S
//  SPANNING
//  ANNOTS
//  MEANS
//  META
//  AXA

// also:
//  internal annot2sp() function (used by spindles analysis)
//  helper function for S2A(waveform)  s2a_helper() 

uint64_t s2a_interp( const std::vector<double> * d ,
                     const std::vector<uint64_t> * tp ,
                     const int i ,
                     const double t );

void timeline_t::annot2signal( const param_t & param )
{

  // create a new signal based on one or more annotations
  if ( ! param.has( "annot" ) ) Helper::halt( "no annotations specified: e.g. annot=A1,A2" );
  std::vector<std::string> anames = param.strvector_xsigs( "annot" );

  // SR of new signals
  const int sr = param.requires_int( "sr" );

  // use instance ID as a numeric value e.g. for NREMC 1 , 2 , 3 
  const bool numeric_instance = param.has( "numeric-inst" );

  // if not otherwise specified, use annot names as new channel labels
  std::vector<std::string> labels = param.has( "label" ) ? param.strvector( "label" ) : anames;

  if ( anames.size() != labels.size() )
    Helper::halt( "label size does not match annot size" );

  // get whole signal size for this SR

  const int np = sr * edf->header.record_duration * edf->header.nr;  

  const uint64_t srtp = (1.0/sr) * globals::tp_1sec;

  // create synthetic signal, 0/1 for presence/absence of the annotation
  for (int a=0; a<anames.size(); a++)
    {
      // does annot exist?
      annot_t * annot = edf->annotations->find( anames[a] );
      if ( annot == NULL ) continue;

      // get all events
      const annot_map_t & events = annot->interval_events;

      // new channel to be populated/added to EDF
      std::vector<double> adat( np , 0 );

      annot_map_t::const_iterator aa = events.begin();
      while ( aa != events.end() )
	{

	  const interval_t & interval = aa->first.interval;
	  
	  // convert from time-points to (nearest) sample-points
	  // (after removing the N+1 end point in annotations)
	  
	  int start = interval.start / srtp;
	  int stop  = (interval.stop-1LLU) / srtp;
	  
	  if ( start < 0 || stop >= np )
	    Helper::halt( "internal error in timeline_t::annot2signal()" );

	  double value = 1;
	  if ( numeric_instance )
	    {
	      if ( aa->first.id == "" || aa->first.id == "." )
		value = 0;
	      else
		{
		  // set value to numeric
		  if ( ! Helper::str2dbl( aa->first.id , &value ) )
		    Helper::halt( "requires numeric instance IDs" ); 
		}
	    }
	  
	  // populate (up to and including the start/stop, as we removed the final +1 TP above)
	  for (int p=start; p<=stop; p++)
	    adat[p] = value;

	  // next annotation
	  ++aa;
	}

      //
      // track total time implicated
      //
   
      int points = 0;
      for (int i=0;i<adat.size();i++) if ( adat[i] > 0 ) ++points;
      double seconds = points / sr;
      int minutes = seconds / 60.0;
      if ( minutes > 0 )
	seconds -= minutes * 60.0;

   
      //
      // write as a new signal
      //
   
      logger << "  adding " << events.size() << " "
	     << anames[a] << " annotations (spanning ";
      if ( minutes > 0 ) logger << minutes << " min " << seconds << " sec)";
      else logger << seconds << " sec)";

      if ( numeric_instance )
	logger << " as numeric instance-ID signal " << labels[a] << "\n";
      else
	logger << " as 0/1 signal " << labels[a] << "\n";
   
      edf->add_signal( labels[a]  , sr , adat );

    }
   
}



void timeline_t::signal2annot( const param_t & param )
{

  //
  // handled separately
  //

  if ( param.has( "waves" ) )
    return signal2annot_cuts( param );
  
  //
  // signal to use
  //

  std::string signal_label = param.requires( "sig" );

  signal_list_t signals = edf->header.signal_list( signal_label );

  if ( signals.size() == 0 ) Helper::halt( "could not find any signals: " + signal_label );

  const int ns = signals.size();


  //
  // S2A encoding
  //

  // encoding=LABEL,lwr,upr
  // encoding=label,val,+win
  // bins=min,max,n
  // q=N
  // pos/neg         --> only make annots for above/below X (abs)
  // pos-pct/neg-pct --> only make annots for above/below X (percentile)
  // encoding=0,360  --> but need to make separate values, i.e. otherwise will be
  //                     one big annot; so cut at zero-crossnnigs (phase 0)
  
  //  VALUE :  X    // --> X+EPS
  //           X-Y
  //           X+Y  // eps

  if ( ! ( param.has( "encoding" )
	   || param.has( "encoding2" )
	   || param.has( "bins" )
	   || param.has( "q")
	   || param.has( "pos" )
	   || param.has( "neg" )
	   || param.has( "pos-pct" )
	   || param.has( "neg-pct" ) )
       ) 
    Helper::halt( "no valid encoding\n"
		  "    encoding=label,value,...\n"
		  " or encoding2=label,value1,value2,...\n"
		  " or bins=min,max,n\n"
		  " or waves/half-waves\n"
		  " or q=n\n"
		  " or pos/neg=value\n"
		  " or pos-pct/neg-pct=pct" );

  bool e2 = param.has( "encoding" );
  bool e3 = param.has( "encoding2" );
  bool eb = param.has( "bins" );
  bool eq = param.has( "q" );
  bool etop = param.has( "pos" );
  bool ebot = param.has( "neg" );
  bool etopp = param.has( "pos-pct" );  
  bool ebotp = param.has( "neg-pct" );

  if ( e2 + e3 + eb + eq + etop + ebot + etopp + ebotp > 1 )
    Helper::halt( "can only specify one of encoding|encoding2|bins|q|waves|half-waves|pos|neg|pos-pct|neg-pct" );

  std::string bin_label = "";
  if ( param.has( "bin-label" ) ) bin_label = param.value( "bin-label" );
  else if ( param.has( "no-bin-label" ) ) bin_label = "";
  else if ( eb ) bin_label = "B";
  else if ( eq ) bin_label = "Q";
  else if ( etop || etopp ) bin_label = "POS";
  else if ( ebot || ebotp ) bin_label = "NEG";

  std::vector<std::string> enc; 
  int nxy = -1;

  const int nq = eq ? param.requires_int( "q" ) : 0 ; 

  if ( eq && ( nq < 1 || nq > 200 ) )
    Helper::halt( "q value must be between 2 and 200" );

  double th = 0;
  if ( etop ) th = param.requires_dbl( "pos" );
  else if ( etopp ) th = param.requires_dbl( "pos-pct" );
  else if ( ebot ) th = param.requires_dbl( "neg" );
  else if ( ebotp ) th = param.requires_dbl( "neg-pct" );

  if ( etopp || ebotp )
    if ( th <= 0 || th >= 1 )
      Helper::halt( "percentile thresholds must be between 0 and 1" );
   
  //
  // get encodings (although Q-encodings and topp/botp are signal specific, so do below)
  //

  if ( e2 ) 
    {
      enc = param.strvector( "encoding" );
      nxy = 2;
    }
  else if ( e3 ) 
    {
      enc = param.strvector( "encoding2" );
      nxy = 3;
    }
  else if ( ! ( eq || etop || ebot || etopp || ebotp ) ) // bins
    {
      // make 'encoding2' style string
      std::vector<double> b = param.dblvector("bins");
      if ( b.size() != 3 ) Helper::halt( "expecting bins=min,max,n" );
      const int n = b[2];
      const double bmin = b[0];
      const double bmax = b[1];
      if ( bmin >= bmax || n == 0 ) Helper::halt( "expecting bins=min,max,m" );
      const double binc = (b[1] - b[0] ) / (double)n;
      logger << "  setting " << n << " bins of interval size " << binc << "\n";

      nxy = 3;
      enc.clear();
      int ndigs = MiscMath::num_digits( n );
      for (int i=0; i<n; i++)
	{
	  enc.push_back( bin_label + Helper::zero_pad( i+1 , ndigs ) ); // annotate label
	  enc.push_back( Helper::dbl2str( bmin + i * binc ) );
	  enc.push_back( Helper::dbl2str( bmin + (i+1) * binc ) );
	}
    }

  if ( enc.size() % nxy != 0 )
    Helper::halt( "requires " + Helper::int2str( nxy) + " args per encoding value" );


  //
  // Either make one annot class (and labels are instances)
  //  or each label --> a distinct class
  //

  const bool use_class = param.has( "class" );

  std::string class_name = use_class ? param.value( "class" ) : "" ; 

  //
  // Append channel name to label as instance ID 
  //

  const bool add_ch_label = param.has( "add-channel-label" ) ? param.yesno( "add-channel-label" ) : false ; 

  //
  // Span EDF discontinuities or no?
  //

  bool span_disc = param.has( "span-gaps" );


  //
  // Parse encodings
  //

  // encode first (label) with XX,1  XX,2  etc to allow
  // duplicate labels; but remove those when printing to annots
  std::map<std::string,std::pair<double,double> > e;

  int num_digs = MiscMath::num_digits( enc.size() );

  for (int i=0; i<enc.size(); i += nxy )
    {

      // make each label unique , i.e. to have a one-to-many mapping
      // of labels to ranges
      std::string label = enc[i] + "," + Helper::zero_pad(i , num_digs) ;
       
      double ex = 0 ;
      if ( ! Helper::str2dbl( enc[i+1] , &ex ) )
	Helper::halt( "bad numeric value for encoding" + enc[i+1] );

      // default window
      bool window = true;
      double ey = 0.05;
   
      // 'encoding2' or 'bins'
      if ( e3 || eb )
	{
	  if ( ! Helper::str2dbl( enc[i+2] , &ey ) )
	    Helper::halt( "bad numeric value for encoding" + enc[i+2] );
	  
	  window = enc[i+2].substr(0,1) == "+" ;
   
	  if ( ! window ) 
	    {
	      if ( ey < ex )
		{
		  double t = ex;
		  ex = ey;
		  ey = t;
		}
	    }
	}

      if ( window )
	{
	  double w = ey;
	  ey = ex + w;
	  ex -= w;	  
	}           
   
      // record
      e[ label ] = std::make_pair( ex , ey );

    }

  // handle abs-threshold (pos/neg) cases (-pct variants below)
  // but add e[] here just so output below is correct (i.e. 1 annot)
  if ( etop || ebot || etopp || ebotp )
    e[ bin_label ] = std::make_pair( th , th );

  logger << "  encoding " << ( eq ? nq : e.size() ) << " annotation instances\n";


  //
  // For each signal
  //

  for (int s=0; s<ns; s++)
    {
   
      if ( edf->header.is_annotation_channel( signals(s) ) )
	continue;
   
      //
      // get signal data
      //

      const std::string ch_label = signals.label(s);

      slice_t slice( *edf , signals(s) , wholetrace() );  

      std::vector<double> * d = slice.nonconst_pdata();  
      
      const std::vector<uint64_t> * tp = slice.ptimepoints();

      //
      // Get quantiles?
      //

      if ( eq )
	{
	  // wipe any current encoding
	  e.clear();

	  int num_digs = MiscMath::num_digits(nq);
	  
	  const double pi = 1 / (double)nq;
	  double p = 0;
	  for (int i=0; i<nq; i++)
	    {
	      const double lwr = MiscMath::percentile( *d , p );
	      const double upr = MiscMath::percentile( *d , i == nq-1 ? 1.0 : p + pi );
	      e[ bin_label + Helper::zero_pad(i+1,num_digs) ] = std::make_pair( lwr , upr );
	      p += pi;	      
	    }
	}

      if ( etop || etopp )
	{
	  // wipe any current encoding
	  e.clear();	  
	  const double upr = etopp ? MiscMath::percentile( *d , 1-th ) : th;
	  e[ bin_label ] = std::make_pair( upr , upr ); // only use 1 val
	}

      if ( ebot || ebotp )
	{
          // wipe any current encoding
	  e.clear();
          const double lwr = ebotp ? MiscMath::percentile( *d , th ) : th;
          e[ bin_label ] = std::make_pair( lwr , lwr ); // only use 1 val	  
	}


      //
      // Naming
      //
      //
      //  e[][int] : key is 'bin' label (e.g. B1, or POS; can be set by bin-label) 
      //    --> by default, this is assigned as the 'class' name
      //    --> if class=XXX then label class name is XXX, and label --> annotation instance ID
      //    --> if add-channel-label=T, then label appended to class label
      //
      //   class  add-ch      class    inst
      //   F      F           label    .
      //   =XX    F           XX       label
      //   F      T           label_CH .
      //   =XX    T           XX_CH    label
      
      //
      // Add annot class? (but not adding a channel label)
      //
      
      int sr = edf->header.sampling_freq( signals(s) );
      uint64_t dt = globals::tp_1sec / sr ;
	
      std::map<std::string,std::pair<double,double> >::const_iterator ee = e.begin();
      
      while ( ee != e.end() )
	{
	  
	  const std::string & label = ee->first;	  

	  // remove xxx,N uniquificaiton
	  const std::vector<std::string> ll = Helper::parse( label , "," );
	  const std::string display_label = ll.size() == 0 ? "." : ll[0] ;
      	    
	  double ex = ee->second.first;
	  double ey = ee->second.second;
	  
	  // get annot_t to add to/
	  //   (note: if exists, then add() returns existing set, so earier
	  //          to use add() rather than find() ) 
	  
	  std::string class_label = use_class ? class_name : display_label ;
	  std::string inst_label  = use_class ? display_label      : "."   ;
	  if ( add_ch_label ) class_label += "_" + ch_label;
	  
	  annot_t * a = annotations->add( class_label );
	  
	  if ( a == NULL ) Helper::halt( "internal error in signal2annot()" );
	  
	  // iterate over signal points, find in-range intervals	  
	  const int n = d->size();
	  if ( n == 0 ) {++ee; continue; }
	  
	  bool in;
	  if      ( etop || etopp ) in = (*d)[0] >= ex;
	  else if ( ebot || ebotp ) in = (*d)[0] <= ex;
	  else                      in = (*d)[0] >= ex && (*d)[0] <= ey ;
	  
	  uint64_t start = (*tp)[0];
	  
	  int cnt = 0;
	  
	  for (int i=0; i<n; i++)
	    {
	      // did we just cross a gap, or is this the last data-point?
	      bool gap = span_disc ? false : ( i != 0 ? discontinuity( *tp , sr , i-1 , i ) : false ) ; 

	      // if ( gap )
	      // 	std::cout << " GAP " << i << " of " << n << "\n";
	      
	      // last observed sample?
	      bool end = i == n - 1;
	      
	      // still in region?
	      bool in1;
	      if      ( etop || etopp ) in1 = (*d)[i] >= ex;
	      else if ( ebot || ebotp ) in1 = (*d)[i] <= ex;
	      else                      in1 = (*d)[i] >= ex && (*d)[i] <= ey;
	      
	      // end of an interval? 
	      if ( in && ( gap || end || ! in1 ) ) 
		{	      

		  // 1-past-end encoding
		  uint64_t stop = end ? last_time_point_tp + 1LLU : (*tp)[i] ;

		  // but adjust for gap (i.e. one sample point from prior point) 
		  if ( gap )
		    {
		      stop = (*tp)[i-1] + dt;
		      //std::cout << "adjusting " << start << " --> " << (*tp)[i] << " becomes " << stop << "\n";
		    }

		  a->add( inst_label , interval_t( start , stop ) , ch_label );
		  
		  // update status (i.e. may still be a new interval after a gap)
		  in = in1;
		  
		  if ( gap && in1 ) 
		    {
		      start = (*tp)[i];
		      // unlikely, but could be gap and then last single sample
		      if ( end )
			{
			  a->add( inst_label , 
				  interval_t( start , last_time_point_tp + 1LLU ) ,
				  ch_label ); 			  
			}
		    }
	      ++cnt;
	    }
	  else if ( in1 && ! in ) // ... or start a new interval
	    {
	      start = (*tp)[i];
	      in = true;
	      if ( i == n - 1 ) // single point interval?
		{
		  a->add( inst_label ,
			  interval_t( start , last_time_point_tp + 1LLU ) ,
			  ch_label );
		}
	    }	      
       }
	  
	  
      logger << "  added " << cnt
	     << " intervals for " << class_label << "/" << inst_label ;

      if ( etop || etopp )
	logger << " based on " 
	       << ch_label << " >= " << ex << "\n";
      else if ( ebot || ebotp )
	logger << " based on " 
	       << ch_label << " <= " << ex << "\n";
      else
	logger << " based on " << ex
	       << " <= " << ch_label << " <= " << ey << "\n";
      
      // next label
      ++ee;
	}
      
      // next signal
    }
  
}


void timeline_t::signal2annot_cuts( const param_t & param )
{
  
  // as signal2annot_cuts() is getting messy, split off cuts
  // functionality here; i.e. expecting angular values, e.g. 0 .. 360
  // , and so a 'cut' of 0 means whenever we 'cross' 0 to start/stop
  // an annotation
  
  
  //
  // labels (always adds signal label X)
  //
  
  if ( param.empty( "waves" ) )
    Helper::halt( "no label specified waves={label}" );
  
  const std::string wave_label = param.value( "waves" ) ;
  
  const bool add_ch_inst_label = param.has( "add-channel-inst-label" ) ? param.yesno( "add-channel-inst-label" ) : false ;

  const bool add_ch_class_label = param.has( "add-channel-class-label" ) ? param.yesno( "add-channel-class-label" ) : false ;

  // pos2neg?
  const bool pos2neg = param.has( "pos2neg" ) ? param.yesno( "pos2neg" ) : false; 


    
  //
  // Selection criteria
  //

  const bool sel_tmin   = param.has( "t-min" );
  const bool sel_tmax   = param.has( "t-max" );

  const double th_tmin = sel_tmin ? param.requires_dbl( "t-min" ) : 0 ;
  const double th_tmax = sel_tmax ? param.requires_dbl( "t-max" ) : 0 ; 

  // also, allow for per phase-bin min/max outlier-detection 
  const bool sel_phbin_tmin   = param.has( "t-min-phbin" );
  const bool sel_phbin_tmax   = param.has( "t-max-phbin" );

  const double th_phbin_tmin = sel_phbin_tmin ? param.requires_dbl( "t-min-phbin" ) : 0 ;
  const double th_phbin_tmax = sel_phbin_tmax ? param.requires_dbl( "t-max-phbin" ) : 0 ;

  
  // percentile-based amg , e.g. mag=20 means in top 20% of average MAG value
  // compared to all other waves
  const bool sel_mag = param.has( "mag-percentile" );
  const double th_mag = sel_mag ? param.requires_dbl( "mag-percentile" ) : 0 ;
  if ( sel_mag && ( th_mag <= 0 || th_mag > 1 ) ) Helper::halt( "mag-percentile must be between 0 and 1" );

  // norm values and only take is above th_magz
  const bool sel_magz = param.has( "mag-z" );
  const double th_magz = sel_magz ? param.requires_dbl( "mag-z" ) : 0 ;

  const bool use_mag = sel_mag || sel_magz ;

  const bool use_mono = param.has( "monotonic" );

  const bool add_slope = param.has( "slope" ); // SLOPE
  const bool add_state = param.has( "state" ); // STATE
  const bool add_bins  = param.has( "bins" );  // BIN (x12 fixed)

  if ( ( add_slope || add_state || add_bins ) && ! use_mono )
    Helper::halt( "requires 'monotonic' flag if using slope, state or bins" );

  if ( ( th_phbin_tmin || th_phbin_tmax ) && !  add_bins )
    Helper::halt( "cannot use t-min-phbin or t-max-phbin without bins" );

  
  //
  // signal(s) to use: assume phase-angles in main signal
  //
  
  //   X --> X_ht_ang   (for defining waves)
  //     --> X_ht_mag   (for any amplitude stuff)
  
  std::string signal_label = param.requires( "sig" );
    
  signal_list_t signals = edf->header.signal_list( signal_label );
  
  if ( signals.size() == 0 )
    Helper::halt( "could not find any signals: " + signal_label );

  const int ns = signals.size();

  
  // allow alternative signal names (if not using HILBERT) 
  
  const std::string ph_ext = param.has( "phase-ext" ) ? param.value( "phase-ext" ) : "_ht_ang";

  const std::string mag_ext = param.has( "mag-ext" ) ? param.value( "mag-ext" ) : "_ht_mag";
  
  //
  // For each signal
  //

  for (int s=0; s<ns; s++)
    {

      int okay_cnt =0;
	  

      const std::string sig_label = signals.label(s);

      const std::string phase_label = sig_label + ph_ext; 
      
      const std::string mag_label   = sig_label + mag_ext; 
      
      const int phase_slot = edf->header.signal( phase_label );

      const int mag_slot = edf->header.signal( mag_label );

      if ( phase_slot == -1 )
	{
	  logger << "  ** could not find " << phase_label << ", skipping...\n";
	  continue;
	}
      
      if ( edf->header.is_annotation_channel( phase_slot ) )
	continue;
      
      if ( use_mag )
	{
	  if ( mag_slot == -1 ) continue;

	  if ( edf->header.is_annotation_channel( mag_slot ) )
	    continue;
	}

      //
      // Annotations
      //
      
      annot_t * a_full = add_ch_class_label ?
	edf->annotations->add( wave_label + "_FULL_" + sig_label ) :
	edf->annotations->add( wave_label + "_FULL" ) ;

      //
      // candidate waves (neg hw, then pos hw)
      //

      std::set<interval_t> fwaves;

      //
      // track QC
      //

      int mono_cnt = 0;
      int mag_cnt = 0;
      int dur_cnt = 0;
      int dur_phbin_cnt = 0;
      int all_cnt = 0;
      
      writer.level( sig_label , globals::signal_strat );

      
      //
      // get signal data: phase
      //

      slice_t slice( *edf , phase_slot , wholetrace() );  
      
      std::vector<double> ph ;
      
      if ( pos2neg )
	{
	  logger << "  shifting phase by 180-degrees to detect positive-to-negative waves\n";
	  // shift phase angles by 180 for pos2neg 
	  std::vector<double> xx = *slice.nonconst_pdata();	  
	  ph.resize( xx.size() );  
	  for (int i=0; i<xx.size(); i++)
	    {
	      if ( xx[i] >= 180 ) ph[i] = xx[i] - 180.0 ;
	      else ph[i] = xx[i] + 180.0 ;
	    }	    
	}
      else
	ph = *slice.nonconst_pdata();  

      const std::vector<uint64_t> * tp = slice.ptimepoints();
      
      const int sr = edf->header.sampling_freq( phase_slot );

      const uint64_t dt = globals::tp_1sec / sr ;
      
      const int n = ph.size();

      if ( n == 0 ) continue;

      // always start not in wave (i.e. have to cross the 0/180/360 points
      // first

      int start = -1 , stop = -1;
      
      //
      // iterate over signal
      //

      // nb. starting at point i=1
      for (int i=1; i<n; i++)
	{
	  
	  // did we just cross a gap? if so, cancel any putative wave
	  bool gap = i != 0 ? discontinuity( *tp , sr , i-1 , i ) : false ; 
	  
	  if ( gap )
	    {	      
	      start = -1;
	      stop = -1;
	      continue;
	    }

	  
	  // did we just cross the key point: because of wrapping
	  //   cross '0-deg'      -->   goes from > 180 to < 180 

	  const bool crosses = ph[i-1] >  180 && ph[i] <  180 ;
	    	  
	  if ( crosses ) 
	    {
	      
	      // close an existing wave:

	      if ( start != -1 )
		{
		  
		  stop = i; // +1 ending
		  
		  //
		  // duration criterion?
		  //
		  
		  bool okay = true;		  
		  
		  if ( sel_tmin || sel_tmax )
		    {
		      const uint64_t start_tp = (*tp)[start];
		      const uint64_t stop_tp  = (*tp)[stop];		      
		      const double dur = ( stop_tp - start_tp ) * globals::tp_duration;
		      if ( sel_tmin && dur < th_tmin ) okay = false;
		      if ( sel_tmax && dur > th_tmax ) okay = false;
		    }

		  if ( okay ) 
		    fwaves.insert( interval_t( (uint64_t)start , (uint64_t)stop ) );
		  else
		    ++dur_cnt; // QC checks

		  
		  ++all_cnt; // all putative waves

		}

	      // and start a new one at same place
	      start = i;
	      stop = -1;
	    }
	  
	} // go to next sample

      
      //
      // Require monotonic? Also that starts/stops in <30 and >330
      //

      if ( use_mono )
	{
	  
	  // select waves
	  std::set<interval_t> fwaves2 = fwaves;
	  fwaves.clear();
	  
	  std::set<interval_t>::const_iterator aa = fwaves2.begin();
          while ( aa != fwaves2.end() )
	    {
	      int start = aa->start;
	      int stop  = aa->stop;

	      bool okay = true;
	      
	      for (int p=start+1;p<stop; p++)
		if ( ph[p] <= ph[p-1] )
		  okay = false;
	      
	      // also flag if not fully spanning the wave-form
	      // i.e. so we can be sure to find all 12 phase bins 
	      if ( ph[start] >= 30 || ph[stop-1] <= 330 )
		okay = false;
	      
	      if ( okay )
		fwaves.insert( *aa );
	      else
		++mono_cnt;
	      
	      ++aa;
	    }

	}

      
      //
      // Require magnitude?
      //
      
      if ( use_mag )
	{
	  slice_t mag_slice( *edf , mag_slot , wholetrace() );
	  
	  std::vector<double> * dm = slice.nonconst_pdata();

	  if ( dm->size() != ph.size() )
	    Helper::halt( "phase and magnitude signals must have the same sample rates" );

	  std::vector<double> v;

	  // get average
	  std::set<interval_t>::const_iterator aa = fwaves.begin();
	  while ( aa != fwaves.end() )
	    {
	      int start = aa->start;
	      int stop  = aa->stop;

	      if ( stop - start == 0 )
		v.push_back( 0 );
	      else
		{
		  double x = 0;
		  for (int p=start; p<stop; p++)
		    x += (*dm)[p] ;
		  v.push_back( x / double( stop - start ) );		  
		}	     
	      ++aa;
	    }

	  // normalize, if needed
	  if ( sel_magz )
	    v = MiscMath::Z( v );
	  
	  // get percentile threshold, if needed
	  double percentile = sel_mag && v.size() > 1 ? MiscMath::percentile( v , th_mag ) : 0 ;
	  
	  // select waves
	  std::set<interval_t> fwaves2 = fwaves;
	  fwaves.clear();
	  
	  int i = 0;
	  aa = fwaves2.begin();
          while ( aa != fwaves2.end() )
	    {
	      bool okay = true;
	      
	      if ( sel_mag && v[i] < percentile ) okay = false;
	      if ( sel_magz && v[i] < th_magz ) okay = false;
	      
	      if ( okay )
		fwaves.insert( *aa );
	      
	      ++i;
	      ++aa;
	    }

	  mag_cnt = fwaves2.size() - fwaves.size();
	  
	}
      
            
      //
      // break down completed cycles into BIN, STATE and SLOPE sub-variables?
      //
      
      
      // between start and stop(-1) we should have a smoothish increasing line of
      // 0 to 360 (assuming reasonable # of samples etc given SR)
      
      // STATE : start .. first past 180
      //         first past 180 .. stop
      
      // SLOPE : 0 to 90 (FALL)
      //         90 to 270 (RISE)
      //         270 to 360 (FALL)

      // track potential exclusions of full waveas based on phbin dur criteria
      std::set<interval_t> fwaves_exclusions;

      if ( add_state || add_slope || add_bins )
	{
	  
	  annot_t * a_pos = NULL;
	  annot_t * a_neg = NULL;
	  annot_t * a_rise = NULL;
	  annot_t * a_fall = NULL;
	  annot_t * a_bins = NULL;
	  
	  if ( add_state )
	    {	      
	      a_pos = add_ch_class_label ?
		edf->annotations->add( wave_label + ( pos2neg ? "_NEG_" : "_POS_" ) + sig_label ) :
		edf->annotations->add( wave_label + ( pos2neg ? "_NEG" : "_POS" ) ) ;
	      
	      a_neg = add_ch_class_label ?
		edf->annotations->add( wave_label + ( pos2neg ? "_POS_" : "_NEG_" ) + sig_label ) :
		edf->annotations->add( wave_label + ( pos2neg ? "_POS" : "_NEG" ) ) ;
	    }

	  if ( add_slope )
	    {	      
	      a_rise = add_ch_class_label ?
		edf->annotations->add( wave_label + ( pos2neg ? "_FALL_" : "_RISE_" ) + sig_label ) :
		edf->annotations->add( wave_label + ( pos2neg ? "_FALL" : "_RISE" ) ) ;
	      
	      a_fall = add_ch_class_label ?
		edf->annotations->add( wave_label + ( pos2neg ? "_RISE_" : "_FALL_" ) + sig_label ) :
		edf->annotations->add( wave_label + ( pos2neg ? "_RISE" : "_FALL" ) ) ;
	    }
	  
	  if ( add_bins )
	    {	      
	      a_bins = add_ch_class_label ?
		edf->annotations->add( wave_label + "_BIN_" + sig_label ) :
		edf->annotations->add( wave_label + "_BIN" ) ;
	    }
	  
	  
	  //
	  // iterate over events, adding if meets criteria
	  //
	  
	  std::set<interval_t>::const_iterator aa = fwaves.begin();
	  
	  while ( aa != fwaves.end() )
	    {
	      
	      int start = aa->start;
	      int stop  = aa->stop;
	      
	      // get adjusted points (interpolate between samples for more
	      // accurate annotation durations)
	      const uint64_t start_tp = s2a_interp( &ph, tp, start, 0 ); 
	      const uint64_t stop_tp  = s2a_interp( &ph, tp, stop, 0 ); 
	      
	      interval_t tinterval = interval_t( start_tp , stop_tp );
	      
	      // find all phase landmarks (30-deg bins)
	      
	      std::vector<uint64_t> bs( 13 , 0LLU );
	      bs[0] = tinterval.start;
	      bs[12] = tinterval.stop;
	      
	      // *because we enfore monotonic* we know that the
	      // wave at least starts/end below 30-deg / after 330-deg,
	      // so we are guaranteed to see all 12 phase bins; but check anyway
	      
	      int bcnt = 0; // should be populating 
	      
	      for (int b=1; b<12; b++)
		{
		  
		  const double phase_angle = b * 30.0;		  
		  
		  for (int p=start;p<stop; p++)
		    {
		      
		      if ( ph[p] >= phase_angle )
			{
			  
			  // get closest to 'phage_angle' (between p and p-1)

			  bs[ b ] = s2a_interp( &ph, tp, p, phase_angle );   
			  
			  // update start to avoid minor retracing of steps
			  // but allow for case where we have a zero-gap (i.e.
			  // make sure each element of bs is filled
			  
			  start = p == start ? start : p - 1;

			  // track that we found this point

			  ++bcnt;

			  // and skip to the next 

			  break;
			}
		    }
		  
		  // next bin
		}


	      //
	      // is this waveform okay?
	      //
	      
	      bool okay = bcnt == 11 ;
	      
	      // additional phase-bin duration criteria?
	      if ( okay && ( sel_phbin_tmin || sel_phbin_tmax ) )
		{
		  
		  for (int b=1; b<13; b++)
		    {
                      const uint64_t tp1 = bs[ b-1 ] ;
                      const uint64_t tp2 = bs[ b ] ;
		      
		      const uint64_t tp = tp2 > tp1 ? tp2 - tp1 : 0LLU ; 
		      const double t = tp * globals::tp_duration ;

		      if ( sel_phbin_tmin && t < th_phbin_tmin ) okay = false;
		      if ( sel_phbin_tmax && t > th_phbin_tmax ) okay = false;
		    }
		}

	      
	      //
	      // skip this interval due to phbin criteria
	      //
	      
	      if ( ! okay )
		{
		  ++dur_phbin_cnt;
		  ++aa;
		  continue;
		}
	      
	      
	      //
	      // add annotations 
	      //
	      
	      a_full->add( add_ch_inst_label ? sig_label : "FULL" , 
			   tinterval , 
			   sig_label );
	      
	      
	      //
	      // add in intervals
	      //

	      if ( add_bins )
		{
		  
		  for (int b=1; b<13; b++)
		    {
		      uint64_t tp1 = bs[ b-1 ] ;
		      uint64_t tp2 = bs[ b ] ;
		      
		      if ( tp2 > tp1 )
			{
			  const std::string ph_label = ( b>9 ? "B" : "B0" ) + Helper::int2str( b );
			  a_bins->add( ph_label , interval_t( tp1 , tp2 ) , sig_label );
			}
		    }
		}
	      
	      
	      if ( add_state )
		{
		  const uint64_t tp180 = bs[ 6 ];

		  // add
		  a_neg->add( add_ch_inst_label ? sig_label : ( pos2neg ? "POS" : "NEG" ) ,
			      interval_t( tinterval.start , tp180 ) , 
			      sig_label );
		  
		  a_pos->add( add_ch_inst_label ? sig_label : ( pos2neg ? "NEG" : "POS" ) ,
			      interval_t( tp180 , tinterval.stop ) , 
			      sig_label );
		  
		}


	      if ( add_slope )
		{
		  const uint64_t tp90  = bs[ 3 ];
		  const uint64_t tp270 = bs[ 9 ];

		  // add
		  a_rise->add( add_ch_inst_label ? sig_label : ( pos2neg ? "FALL" : "RISE" ),
			       interval_t( tp90  , tp270 ) , 
			       sig_label );
		  
		  a_fall->add( add_ch_inst_label ? sig_label : ( pos2neg ? "RISE" : "FALL" ),
			       interval_t( tinterval.start , tp90 ) , 
			       sig_label );
		  
		  a_fall->add( add_ch_inst_label ? sig_label : ( pos2neg ? "RISE" : "FALL" ),
			       interval_t( tp270  , tinterval.stop ) , 
			       sig_label );
		  
		}
	      
	      	     

		  //
		  // next wave
		  //
		  
		  ++aa;
		}
	      
	      

	      logger << "  added " << fwaves.size()
		     << " waves for " << sig_label << "\n";
	      
	      //
	      // QC outputs
	      //
	      
	      writer.value( "EXC1_DUR"  , dur_cnt );
	      writer.value( "EXC2_MONO" , mono_cnt );
	      writer.value( "EXC3_MAG"  , mag_cnt );
	      writer.value( "EXC4_PDUR"  , dur_phbin_cnt );
	      writer.value( "N" , (int)fwaves.size() );
	      writer.value( "N0" , all_cnt );
	      
	      
	}      
      
      // next signal
    }
  
  writer.unlevel( globals::signal_strat );
  
}



void timeline_t::list_spanning_annotations( const param_t & param )
{
    
  if ( mask_set ) 
    Helper::halt( "cannot run SPANNING with a MASK set... use RE" );

  
  // currently, SPANNING only for continuous EDFs
  // if ( ! edf->header.continuous )
  //   Helper::halt( "currently, can only run SPANNING on continuous EDF" );
  
  // given a /set/ of annotations, determine 
  //   - seconds outside of EDF
  //   - total duration of signal covered by these (seconds)
  //   - coverage as a proportion of EDF file
  //   - coverage as a proportion of in-memory representation
  //   - number of contiguous blocks of the requested annotations
  // etc
  

  //
  // which signals: either look at all, or the requested set
  //
  
  std::vector<std::string> requested = param.has( "annot" ) && param.value( "annot" ) != "." 
    ? param.strvector_xsigs( "annot" ) 
    : annotations->names() ;
  

  //
  // Get all annotations (i.e. not stratified by epoch), sort by time and collapse
  //
  

  std::set<instance_idx_t> events;


  //
  // iterate over each annotation
  //

  for (int a = 0 ; a < requested.size() ; a++ ) 
    {
      
      annot_t * annot = annotations->find( requested[a] );
      
      if ( annot == NULL ) continue;
      
      const int num_events = annot->num_interval_events();
      
      //
      // iterator over interval/event map
      //
      
      annot_map_t::const_iterator ii = annot->interval_events.begin();
      while ( ii != annot->interval_events.end() )
	{	  
	  const instance_idx_t & instance_idx = ii->first;	  
	  events.insert( instance_idx );	  
	  ++ii;
	}
      
      
    }


  //
  // track total coverage, etc
  //

  uint64_t total = 0;
  
  uint64_t total_all = 0;

  uint64_t total_collapsed = 0;
  
  uint64_t invalid_tps = 0;

  int over_extended = 0;

  int annot_blocks = 0;

  // keep track of where longest spanning annot reaches to 
  // or 0 if past the previous spanning annot
  
  uint64_t earliest = 0;

  uint64_t furthest = 0;
  
  // TODO: track parts of recording that are not spanned
  //std::set<interval_t> unspanned;

  std::set<instance_idx_t>::const_iterator aa = events.begin();
  while ( aa != events.end() )
    {
      
      const interval_t & interval = aa->interval;
      
      //
      // track total (uncollapsed) duration across all ANNOTs
      // i.e. whether valid or not
      //

      total_all += interval.duration();
      
      // 
      // what overlap, if any?
      //
      
      uint64_t vtp = valid_tps( interval );
      
      bool is_valid = interval.duration() == vtp;

      if ( ! is_valid ) 
	{

	  // duration of annots that do not map to a EDF region
	  invalid_tps +=  interval.duration() - vtp; 
	  
	  // count of intervals that do not perfectly match valid regions
	  ++over_extended;
	  
	  // report
	  writer.level( over_extended , globals::count_strat );
	  
	  writer.value( "ANNOT" , aa->parent->name );
	  writer.value( "INST" , aa->id );
	  writer.value( "START" , interval.start_sec() );
	  writer.value( "STOP" , interval.stop_sec() );
	  writer.value( "DUR" , interval.stop_sec() - interval.start_sec() );
	  writer.unlevel( globals::count_strat );
	}
      
          
      //
      // track collapsed duration, but here only consider completely 'valid' intervals
      //

      if ( is_valid ) 
	{
	  
	  // only count whole annotations for this total 
	  // (i.e. entire annot must be contained in a contiguous segment
	  //   of the record )

	  total += interval.duration();

	  // start of a 'new' region?
	  
	  if ( furthest == 0 ) 
	    {
	      earliest = interval.start;
	      furthest = interval.stop;
	      ++annot_blocks;
	    }
	  else // we already have at least one region counted
	    {
	      
	      // is the old region finished?  if so, add
	      if ( interval.start > furthest ) 
		{
		  total_collapsed += furthest - earliest ;
		  earliest = interval.start;
		  furthest = interval.stop;
		  ++annot_blocks; // track that this starts a new block
		}
	      else // add to current region 
		{
		  if ( interval.stop > furthest )
		    furthest = interval.stop;
		}
	      
	    }
	 	  
	}
	  
      // next segment
      ++aa;
      
    }
  
  // add final interval(s)
  
  total_collapsed +=  furthest - earliest ;
  

  //
  // Report
  //
  

  writer.value( "REC_SEC" , Helper::tp2sec( total_duration_tp )  );
  writer.value( "REC_HMS" , Helper::timestring( total_duration_tp , ':' )  );

  writer.value( "ANNOT_N" , (int)events.size() );
  writer.value( "ANNOT_SEC" , Helper::tp2sec( total )  );
  writer.value( "ANNOT_HMS" , Helper::timestring( total , ':' )  );

  
  // do any (valid) annots overlap each other?
  writer.value( "ANNOT_OVERLAP" , ( total_collapsed < total ? "YES" : "NO" )  );

  // how many annots over-extended beyond range of EDF?
  writer.value( "INVALID_N" , over_extended );
  writer.value( "VALID_N" , (int)events.size() - over_extended );
  
  // number of annotation segments, i.e. annotation-based analog of 
  // the SEGMENTS command  
  writer.value( "NSEGS" , annot_blocks );

  // // total annotation duration, whether overlaping or not
  // writer.value( "ALL_ANNOT_SEC" , Helper::tp2sec( total_all )  );
  
  // extent of this over-extension
  writer.value( "INVALID_SEC" , Helper::tp2sec( invalid_tps ) );
  
  writer.value( "SPANNED_PCT" , 100 * ( Helper::tp2sec( total_collapsed ) / Helper::tp2sec( total_duration_tp ) )  );
  writer.value( "SPANNED_SEC" , Helper::tp2sec( total_collapsed )  );
  writer.value( "SPANNED_HMS" , Helper::timestring( total_collapsed , ':' ) );

  writer.value( "UNSPANNED_SEC" , Helper::tp2sec( total_duration_tp - total_collapsed ) );
  writer.value( "UNSPANNED_PCT" , 100 * ( 1 - Helper::tp2sec( total_collapsed ) / Helper::tp2sec( total_duration_tp ) )  );
  writer.value( "UNSPANNED_HMS" , Helper::timestring( total_duration_tp - total_collapsed , ':' ) ); 
    
}


void timeline_t::list_all_annotations( const param_t & param )
{

  //
  // Options
  //

  // count annotations per epoch
  bool per_epoch = param.has( "epoch" );

  // do this either way, as EDF+ mode requires epochs to locate annots
  if ( ! epoched() ) 
    {
      int ne = set_epoch( globals::default_epoch_len , globals::default_epoch_len );
      logger << "  set epochs to default " << globals::default_epoch_len << " seconds, " << ne << " epochs\n";
    }

  // how to decide whether an interval overlaps a mask or not?
  //  start  -- keep annotations that start in an unmasked region
  //  any    -- keep annotations that have any overlap in an unmasked region
  //  all    -- only keep annotations that are completely within unmasked regions
  
  int keep_mode = 0; 
  if ( param.has( "any" ) ) keep_mode = 0;
  if ( param.has( "all" ) ) keep_mode = 1;
  if ( param.has( "start" ) ) keep_mode = 2;  
  
  logger << "  keeping annotations based on ";
  if ( keep_mode == 0 )      logger << "any overlap with";
  else if ( keep_mode == 1 ) logger << "complete (all) overlap with";
  else if ( keep_mode == 2 ) logger << "starting in";
  logger << " an unmasked region\n";
  
  bool show_masked = param.has("show-masked");
  if ( show_masked )
    logger << "  and also showing masked annotations\n";
  
  // annotation names

  std::vector<std::string> names = annotations->names();

  // restrict to a subset? (allow wildcards here as well as xsigs )
  std::set<std::string> req_annots;

  if ( param.has( "annot" ) )
    req_annots = annotate_t::root_match( param.strset_xsigs( "annot" ) , names );

  const bool restricted = req_annots.size();
  
  //
  // Per epoch summary of all annotations
  //

  if ( per_epoch ) 
    {
      
      first_epoch();
      
      while ( 1 ) 
	{
	  
	  int e = show_masked ? next_epoch_ignoring_mask() : next_epoch();
	  
	  if ( e == -1 ) break;
	  
	  writer.epoch( display_epoch( e ) );

	  interval_t interval = epoch( e );

	  // get each annotations
	  for (int a=0;a<names.size();a++)
	    {
	      
	      // ignore this annot?
	      if ( restricted && req_annots.find( names[a] ) == req_annots.end() )
		continue;
	      
	      annot_t * annot = annotations->find( names[a] );

	      // get overlapping annotations for this epoch
	      annot_map_t events = annot->extract( interval );

	      // list
	      annot_map_t::const_iterator ii = events.begin();
	      while ( ii != events.end() )
		{	  
		  const instance_idx_t & instance_idx = ii->first;

		  const instance_t * instance = ii->second;

		  const interval_t & interval = instance_idx.interval;
		  
		  bool is_masked = false;
		  
		  // keep if any part of A overlaps any unmasked region
		  if      ( keep_mode == 0 ) is_masked = ! interval_overlaps_unmasked_region( interval );
		  
		  // ...or, only if entire A is in unmasked region
		  else if ( keep_mode == 1 ) is_masked = ! interval_is_completely_unmasked( interval );
		  
		  // ...or, if start of A is in an unmasked region
		  else if ( keep_mode == 2 ) is_masked = interval_start_is_masked( interval ) ;

		  // skip?
		  if ( is_masked && ! show_masked ) { ++ii; continue; } 
		  
		  // else display
		  writer.level( instance_idx.id , "INST" );
		  writer.level( interval.as_string() , "INTERVAL" );
		  writer.level( instance_idx.ch_str , globals::signal_strat );

		  writer.value( "EMASK" , masked( e ) );
		  writer.value( "AMASK" , is_masked );
		  
		  
		  ++ii;
		}      

	      writer.unlevel( "INTERVAL" );
	      writer.unlevel( "INST" );
	      writer.unlevel( globals::signal_strat );

	    }
	}
      
      writer.unepoch();

      // all done now for epoch-stratified listing
      return;
    }
  


  //
  // Get all annotations (i.e. not stratified by epoch)
  //

  
  // sort by time, collapse across events
  std::map<instance_idx_t,const instance_t*> events;
  
  // class
  std::map<std::string,int> counts;
  std::map<std::string,double> dur;
  
  // class x inst
  std::map<std::string,std::map<std::string,int> > counts2;
  std::map<std::string,std::map<std::string,double> > dur2;
  
  
  // iterate over each annotation
  for (int a = 0 ; a < names.size() ; a++ ) 
    {
      
      // ignore this annot?                                                                                                                           
      if ( restricted && req_annots.find( names[a] ) == req_annots.end() )
	continue;

      annot_t * annot = annotations->find( names[a] );
      
      if ( annot == NULL ) Helper::halt( "internal problem in list_all_annotations()" );

      const int num_events = annot->num_interval_events();
      
      if ( 0 ) 
	{
	  std::cout << names[a] << "\n";
	  std::cout << " ne = " << num_events << "\n";
	  std::cout << " file = " << annot->file << "\n";
	  
	  const int nf = annot->types.size();
	  std::cout << " fields = " << nf << "\n";

	  std::map<std::string, globals::atype_t>::const_iterator tt = annot->types.begin();
	  while ( tt != annot->types.end() )
	    {
	      std::cout << "  " << tt->first << ", is " << globals::type_name[ tt->second ] << "\n";
	      ++tt;
	    }
	  std::cout << "\n";	  

	}

      
      //
      // iterator over interval/event map
      //

      annot_map_t::const_iterator ii = annot->interval_events.begin();
      while ( ii != annot->interval_events.end() )
	{	  

	  const instance_idx_t & instance_idx = ii->first;
	  const instance_t * instance = ii->second;

	  bool keep_this = false;

	  // allow for 0-duration annots: in EDF+D mode, these
	  // functions (interval_overlaps_unmasked_region() etc) will
	  // not return anything, as they end up calling a function to determine
	  // directly the record count spanned.   As we don't want to mess w/ those
	  // deep functions for now, just make a fix here
	  // all other uses of masked_interval() etc are based on epochs, which
	  // will never have zero-duration
	  
	  // 0-duration time-stamps changed to have an arbitrary 1LLU duration:
	  interval_t search = instance_idx.interval;
	  if ( search.duration() == 0LLU ) search.stop += 1LLU;

	  if ( show_masked )
	    keep_this = true;
	  else
	    {
	      if      ( keep_mode == 0 ) keep_this = interval_overlaps_unmasked_region( search );
	      else if ( keep_mode == 1 ) keep_this = interval_is_completely_unmasked( search );
	      else if ( keep_mode == 2 ) keep_this = ! interval_start_is_masked( search );
	    }

	  if ( keep_this )
	    {      	      
	      events[ instance_idx ] = instance ; 
	      
	      counts[ annot->name ]++;
	      counts2[ annot->name ][ instance_idx.id ] ++;
	      
	      dur[ annot->name ] += instance_idx.interval.duration_sec(); 	      
	      dur2[ annot->name ][ instance_idx.id ] += instance_idx.interval.duration_sec(); 	      
	    }
	      
	  
	  ++ii;
	}
    }

  //
  // report HMS?
  //
  
  clocktime_t starttime( edf->header.starttime );
  bool hms = true;  
  if ( ! starttime.valid )
    {
      logger << "  *** could not find valid start-time in EDF header ***\n";
      hms = false;
    }


  // now print all by time point
  std::map<instance_idx_t,const instance_t*>::const_iterator aa = events.begin();
  while ( aa != events.end() )
    {
      
      const instance_idx_t & instance_idx = aa->first;
      
      const interval_t & interval = instance_idx.interval;
      
      const instance_t * instance = aa->second;

      // stratify output by interval
    
      writer.interval( interval );
      
      writer.level( instance_idx.parent->name , globals::annot_strat );
      
      writer.level( instance_idx.id , globals::annot_instance_strat );	  

      writer.value( "START" , interval.start_sec() );

      // NOTE: not sure why we previously did this... for output only, keep consistent form
      // do not +1 time-unit
      //writer.value( "STOP" , interval.stop_sec_exact() );

      writer.value( "STOP" , interval.stop_sec() );
      
      writer.value( "DUR" , interval.stop_sec() - interval.start_sec() );

      // channel label
      writer.value( "CH" , instance_idx.ch_str );

      // HMS : elapsed
      // HMS : clock

      if ( hms )
	{
	  
	  double tp1_sec = interval.start_sec();
	  clocktime_t present1 = starttime;
	  //present1.advance( tp1_sec / 3600.0 );
	  present1.advance_seconds( tp1_sec );
	  
	  // add down to 1/100th of a second
	  double tp1_extra = tp1_sec - (long)tp1_sec;

	  // Not sure why we used this form previously... to be consistent, stick with STOP being +1 end
	  
	  // stop_sec_exact() return last time point (rather than usual 1-past-the-end)
	  //double tp2_sec = interval.stop_sec_exact();

	  double tp2_sec = interval.stop_sec();

	  clocktime_t present2 = starttime;
	  //present2.advance( tp2_sec / 3600.0 );
	  present2.advance_seconds( tp2_sec );
	  
	  double tp2_extra = tp2_sec - (long)tp2_sec;
	   
	  writer.value( "START_HMS"  , present1.as_string(':') +  Helper::dbl2str_fixed( tp1_extra , globals::time_format_dp  ).substr(1) );
	  writer.value( "STOP_HMS"   , present2.as_string(':') +  Helper::dbl2str_fixed( tp2_extra , globals::time_format_dp  ).substr(1) );

	  // elapsed time (00:00:00 is start of EDF)
	  clocktime_t present3;
	  present3.advance_seconds( tp1_sec );	  	  
	  // add down to 1/100th of a second
	  tp1_extra = tp1_sec - (long)tp1_sec;

	  clocktime_t present4;	  
	  present4.advance_seconds( tp2_sec );	  
	  tp2_extra = tp2_sec - (long)tp2_sec;
	   	  
	  // std::cout << "xx\t" << present3.as_string(':') 
	  // 	    << "\t" 
	  // 	    << tp1_sec << "\t" 
	  // 	    << (long)tp1_sec << "\t"
	  // 	    << tp1_extra << "\t" 
	  // 	    << Helper::dbl2str_fixed( tp1_extra , globals::time_format_dp ) << "\n";

	  writer.value( "START_ELAPSED_HMS"  , present3.as_string(':') +  Helper::dbl2str_fixed( tp1_extra , globals::time_format_dp ).substr(1) );
	  writer.value( "STOP_ELAPSED_HMS"   , present4.as_string(':') +  Helper::dbl2str_fixed( tp2_extra , globals::time_format_dp ).substr(1) );


	}
	

      if ( ! instance->empty() ) 
	{
	  writer.value(  "VAL" , instance->print() );
	}

      if ( show_masked ) 
	{
	  
	  bool start_masked = interval_start_is_masked( interval ) ; 
	  bool some_masked = interval_overlaps_masked_region( interval );
	  bool all_masked = interval_is_completely_masked( interval );
	  bool some_unmasked = interval_overlaps_unmasked_region( interval );
	  bool all_unmasked = interval_is_completely_unmasked( interval );
	  
	  writer.value( "START_MASKED"  , start_masked );
	  writer.value( "SOME_MASKED"   , some_masked );
	  writer.value( "ALL_MASKED"    , all_masked );
	  writer.value( "SOME_UNMASKED" , some_unmasked );
	  writer.value( "ALL_UNMASKED"  , all_unmasked );
	}
      
      writer.unlevel( globals::annot_instance_strat );

      writer.unlevel( globals::annot_strat );
      
      ++aa;
      
    }
  writer.uninterval();


  //
  // final counts, durations by class
  //
  
  std::map<std::string,int>::const_iterator cc = counts.begin();
  while ( cc != counts.end() ) 
    {

      writer.level( cc->first , globals::annot_strat );
      writer.value( "COUNT" , cc->second );      
      writer.value( "DUR" , dur[ cc->first ] );

      if ( counts2[ cc->first ].size() > 0 )
	{
	  std::map<std::string,int>::const_iterator dd = counts2[ cc->first ].begin();
	  while ( dd != counts2[ cc->first ].end() )
	    {	      
	      writer.level( dd->first , globals::annot_instance_strat );	  
	      writer.value( "COUNT" , dd->second );      
	      writer.value( "DUR" , dur2[ cc->first ][ dd->first ] );
	      ++dd;
	    }
	  writer.unlevel( globals::annot_instance_strat );
	}

      ++cc;
    }
  writer.unlevel( globals::annot_strat );
}


int timeline_t::annot2sp( edf_t & edf , const std::string & astr ,
			  bool only_this_channel ,
			  std::vector<interval_t> * sample_points , 
			  std::vector<interval_t> * time_points , 
			  int * orig_n , 
			  std::string ch , int sr )
{
  
  sample_points->clear();
  time_points->clear();
  
  // use "CH" to get SR (unless it is otherwise specified)
  //  but read all annots (irrespective of channel) unless only_this_channel == 1
  //  (in which case, we require a specified CH ratehr than a SR (where we simply
  //   find the first match) 

  if ( only_this_channel && ( ch == "" || ch == "." ) )
    Helper::halt( "require a specified channel for annot2sp() " );
  
  // either, find the SR of the given channel:
  if ( sr == 0 )
    {
      signal_list_t signals = edf.header.signal_list( ch );
      if ( signals.size() == 0 ) return 0;
      if ( signals.size() != 1 ) Helper::halt( "problem matching a single channel" );
      std::vector<double> Fs = edf.header.sampling_freq( signals );      
      sr = Fs[0];
    }
  else
    {      
      signal_list_t signals = edf.header.signal_list( "*" );
      std::vector<double> Fs = edf.header.sampling_freq( signals );
      for (int s=0; s<Fs.size(); s++)
	{
	  if ( (int)(Fs[s]) == sr )
	    {
	      ch = signals.label( s );
	      break;
	    }
	}
    }

  
  if ( sr == 0 || ch == "" || ch == "." )
    Helper::halt( "problem finding a channel w/ SR matching" );


  signal_list_t signals = edf.header.signal_list( ch );
  if ( signals.size() != 1 ) Helper::halt( "problem matching a single channel" );
    
  logger << "  using " << ch << " (SR = " << sr << ") to align annotations to sample-points\n";

  
  // must map to within 1 sample (n.b. if at edge, ignored)
  const double max_diff = 1/(double)sr;
  logger << "  mapping to closest sample-point within " << max_diff << " seconds\n";
  
  
  //
  // get the annotation
  //
  
  annot_t * annot = annotations->find( astr );
  
  if ( annot == NULL )
    Helper::halt( "could not find annotation class " + astr );

  
  
  //
  // get time-points for this SR (via pull of a dummy channel) 
  //

  slice_t slice( edf , signals(0) , edf.timeline.wholetrace() );
  
  const std::vector<uint64_t> * tp = slice.ptimepoints();

  const int np = tp->size();
  
  //
  // Iterate over elements, and build a single ordered table of all times
  //
  
  std::map<uint64_t,int> times;

  *orig_n = 0;
  
  const annot_map_t & events = annot->interval_events;  
  annot_map_t::const_iterator aa = events.begin();
  while ( aa != events.end() )
    {
      // get this annot interval
      const instance_idx_t instance = aa->first;
      
      bool add = ( ! only_this_channel ) || instance.ch_str == ch ; 
      if ( add )
	{
	  // track count
	  (*orig_n)++;

	  times[ instance.interval.start ] = -1;
	  times[ instance.interval.stop ] = -1;
	}
      ++aa;      
    }
  
  //
  // now map all *unique & sorted* times 
  //
  
  // index of tp-map (starts at 1)
  int idx = 1;

  std::map<uint64_t,int>::iterator tt = times.begin();
  while ( tt != times.end() )
    {
      
      const uint64_t & curr = tt->first;
      const uint64_t & prior = (*tp)[idx-1];
      const uint64_t & next = (*tp)[idx];
      
      // shift sample-point window up
      if ( next < curr )
	{	  
	  ++idx;
	  if ( idx == np ) break;
	  continue; // i.e. bounce back but do not update ++tt
	}

      //      std::cout << "\n considering = " << curr << "\n";      
      //std::cout << " idx now = " << idx << " ; curr , prior next = " << curr << " " << prior << " " << next << "\n";
      
      // is in-between these two points?
      if ( curr >= prior && curr <= next )
	{
	  const uint64_t d1 = curr - prior;
	  const uint64_t d2 = next - curr ;
	  
	  const bool first = d1 < d2 ; 
	  
	  //	  std::cout << " closest first = " << first  << "\n";
	  
	  double df = ( first ? d1 : d2 ) * globals::tp_duration ; 
	  
	  //std::cout << "  in nbetween = " << df << "\n";
	  
	  // close enough?
	  if ( df <= max_diff )
	    {
	      int sp = first ? idx-1 : idx ;

	      //std::cout << " storing " << sp << "\n";

	      // store
	      tt->second = sp ;	      
	    }

	}
	  

      // advance to next point 
      
      ++tt;
    }
  
  
  //
  // map back to starts and stops
  //
      
  aa = events.begin();
  while ( aa != events.end() )
    {
      
      int start = -1 , stop = -1;

      const instance_idx_t & instance = aa->first;
      
      bool add = ( ! only_this_channel ) || instance.ch_str == ch ; 
      
      if ( add )
	{
	  if ( times.find( instance.interval.start ) != times.end() )
	    start = times[ instance.interval.start ];
	  
	  if ( times.find( instance.interval.stop ) != times.end() )
	    stop = times[ instance.interval.stop ];
	  
	  // std::cout << " start , stop = " << start << "\t" << stop << " <---- "
	  // 	    << instance.interval.start << "  " << instance.interval.stop << "\n";
	  
	  // add this event (both SP and TP to ensure these are aligned
	  // in the returned value (i.e. to go to makeing a spindle_t ) 
	  if ( start != -1 && stop != -1 )
	    {
	      sample_points->push_back( interval_t( (uint64_t)start , (uint64_t)stop ) );
	      time_points->push_back( instance.interval );
	    }
	}
      
      ++aa;

    }

  
  return sample_points->size();
  
}


void timeline_t::signal_means_by_annot( const param_t & param )
{

  
  //
  // annots
  //

  // for root-match
  std::vector<std::string> names = annotations->names();
  
  if ( ! param.has( "annot" ) ) Helper::halt( "no annotations specified: e.g. annot=A1,A2" );

  std::vector<std::string> anames = Helper::set2vec( annotate_t::root_match( param.strset_xsigs( "annot" ) , names ) );

  
  //
  // ignore annotation instannce IDs?
  //

  const bool ignore_instance_ids = ! param.has( "by-instance" );


  //
  // min-max normalized means 
  //

  const bool norms = param.has( "norm" );

  //
  // flanking windows
  //

  const bool flanking = param.has( "w" );

  const double flanking_tp = flanking ? param.requires_dbl( "w" ) * globals::tp_1sec : 0 ;
  
  //
  // signals
  //
  
  std::string signal_label = param.requires( "sig" );
  const bool no_annotations = true;
  signal_list_t signals = edf->header.signal_list( signal_label , no_annotations );  
  const int ns = signals.size();

  if ( ns == 0 ) return;
  
  const int Fs = edf->header.sampling_freq( signals(0) );
  for (int s=1; s<ns; s++)
    if ( edf->header.sampling_freq( signals(s) ) != Fs )
      Helper::halt( "signals must have similar sampling rates" );
  
    
  //
  // stores
  //

  // class -> [instance] -> N   [ assumes same SR across channel ]
  // class -> [instance] -> channel -> sum 
  
  std::map<std::string,std::map<std::string,int> > an;
  std::map<std::string,std::map<std::string,std::map<int,double> > > ax;

  // flanking : prior
  std::map<std::string,std::map<std::string,int> > left_an;
  std::map<std::string,std::map<std::string,int> > right_an;

  std::map<std::string,std::map<std::string,std::map<int,double> > > left_ax;
  std::map<std::string,std::map<std::string,std::map<int,double> > > right_ax;

  //
  // iterate over annots
  //


  for (int a=0; a<anames.size(); a++)
    {
      
      // does annot exist?
      annot_t * annot = edf->annotations->find( anames[a] );
      if ( annot == NULL ) continue;
      const std::string & class_name = anames[a];
      
      // get all events
      const annot_map_t & events = annot->interval_events;

      //std::cout << " considering " << events.size() << "\n";
      //int idx = 0;
      
      annot_map_t::const_iterator aa = events.begin();
      while ( aa != events.end() )
	{
	  // instance ID (or not)
	  const std::string inst_id = ignore_instance_ids ? "." : aa->first.id ;
	  
	  // ++idx;
	  // std::cout << " idx " << idx << " class_name = " << class_name <<  " " << inst_id << "\n";

	  
	  // get main interval
	  const interval_t & interval = aa->first.interval;
      
	  eigen_matslice_t mslice( *edf , signals , interval );	  
	  
	  const Eigen::MatrixXd & X = mslice.data_ref();
          const int rows = X.rows();
          const int cols = X.cols();

	  // add to count, accumulate mean
	  an[ class_name ][ inst_id ] += rows;	  
	  Eigen::ArrayXd sum = X.array().colwise().sum();	  
	  for (int s=0; s<ns; s++)
	    ax[ class_name ][ inst_id ][ s ] += sum[ s ];

	  // repeat for flanking regions?
	  if ( flanking )
	    {
	      // left
	      interval_t left = interval;
	      left.shift_left( flanking_tp );	      
	      eigen_matslice_t left_mslice( *edf , signals , left );
	      const Eigen::MatrixXd & left_X = left_mslice.data_ref();
	      left_an[ class_name ][ inst_id ] += left_X.rows();
	      Eigen::ArrayXd left_sum = left_X.array().colwise().sum();
	      for (int s=0; s<ns; s++)
		left_ax[ class_name ][ inst_id ][ s ] += left_sum[ s ];

	      // right
	      interval_t right = interval;
	      right.shift_right( flanking_tp );	      
	      eigen_matslice_t right_mslice( *edf , signals , right );
	      const Eigen::MatrixXd & right_X = right_mslice.data_ref();
	      right_an[ class_name ][ inst_id ] += right_X.rows();
	      Eigen::ArrayXd right_sum = right_X.array().colwise().sum();
	      for (int s=0; s<ns; s++)
		right_ax[ class_name ][ inst_id ][ s ] += right_sum[ s ];	      	      
	    }
	  
	  // next annotation
	  ++aa;
	}

    } // next annotation
  
  
  //
  // Report means, along w/ normalized values (optionally) 
  //

 
  // first, by by channel

  for (int s=0; s<ns; s++)
    {

      writer.level( signals.label( s ) , globals::signal_strat );
      
      // class -> [instance] -> N   [ assumes same SR across channel ]
      // class -> [instance] -> channel -> sum 
      // std::map<std::string,std::map<std::string,int> > an;
      // std::map<std::string,std::map<std::string,std::map<int,double> > > ax;  
      //                                                     ^
      //                                                     |
      //                                                 channel slot

      //
      // Build norm tables (min/max ranges, within annot class only)
      //
      
      std::map<std::string,double> ann2inst_min;
      std::map<std::string,double> ann2inst_max;
      
      std::map<std::string,std::map<std::string,int> >::const_iterator ii = an.begin();
      while ( ii != an.end() )
	{
 
	  double minval = 0 , maxval = 0;
	  
	  // by instance ID
	  std::map<std::string,int>::const_iterator jj = ii->second.begin();
          while ( jj != ii->second.end() )
            {
	      
	      const double x = ax[ ii->first ][ jj->first ][ s ] / (double)jj->second;

	      // min/max
	      if ( jj == ii->second.begin() )
		minval = maxval = x;
	      else if ( x < minval )
		minval = x;
	      else if ( x > maxval )
		maxval = x;
	      
	      ++jj;
	    }

	  // get mean for this class
	  ann2inst_min[ ii->first ] = minval;
	  ann2inst_max[ ii->first ] = maxval;
	  
	  ++ii;
	}

      
      //
      // Report outputs
      //
      
      // by annotation class
      ii = an.begin();
      while ( ii != an.end() )
	{
	  
	  writer.level( ii->first , globals::annot_strat );
	  
	  const bool do_inst_norms = ii->second.size() > 2 && ! ignore_instance_ids ;
	  
	  // by instance ID 	  
	  std::map<std::string,int>::const_iterator jj = ii->second.begin();
	  while ( jj != ii->second.end() )
	    {
	      
	      // if ignoring instance IDs, then only a single '.' here, so skip
	      // adding as a factor
	      if ( ! ignore_instance_ids ) 
		writer.level( jj->first , globals::annot_instance_strat );
	      
	      // by channel
	      // ax = ax[ ii->first ][ jj->first ][ s ] 
	      const double x = ax[ ii->first ][ jj->first ][ s ] / (double)jj->second ;
	      
	      // main mean
	      writer.value( "M" , x );
	      writer.value( "S" , jj->second / (double)Fs ); // span in seconds
	      
	      // normed? 
	      if ( do_inst_norms )
		{
		  writer.value( "M1" , ( x - ann2inst_min[ ii->first ] ) / ( ann2inst_max[ ii->first ] - ann2inst_min[ ii->first ] ) );
		}
	      
		
	      // flanking regions?
	      if ( flanking )
		{
		  writer.value( "L" , left_ax[ ii->first ][ jj->first ][ s ] / (double)left_an[ ii->first ][ jj->first ] ) ;
		  writer.value( "R" , right_ax[ ii->first ][ jj->first ][ s ] / (double)right_an[ ii->first ][ jj->first ] ) ;
		}

	      // next instance
	      ++jj;
	    }	  

	  if ( ! ignore_instance_ids )
	    writer.unlevel( globals::annot_instance_strat );
	  
	  // next class
	  ++ii;	  
	}

      writer.unlevel( globals::annot_strat );
      
      // next channel
    }
  
  writer.unlevel( globals::signal_strat );
  
  // all done
}


std::set<interval_t> timeline_t::gaps( const std::set<interval_t> & segs )
{
  std::set<interval_t> g;

  // no segs - implies one big gap
  if ( segs.size() == 0 )
    {
      g.insert( interval_t( 0LLU , last_time_point_tp + 1LLU  ) );
      return g;
    }

  // start/end?
  std::set<interval_t>::const_iterator ss = segs.begin();
  if ( ss->start != 0LLU ) g.insert( interval_t( 0LLU , ss->start ) );
  
  std::set<interval_t>::const_iterator ss1 = segs.begin();
  ++ss1;
  
  while ( ss1 != segs.end() )
    {
      g.insert( interval_t( ss->stop , ss1->start ) );
      ++ss1;
      ++ss;
    }

  // end? back up one
  --ss1;
  if ( ss1->stop != last_time_point_tp + 1LLU )
    g.insert( interval_t( ss1->stop , last_time_point_tp + 1LLU ) );
  return g;
}

std::set<interval_t> timeline_t::segments() 
{
  // return a list of current segments, i.e. mirroring the internal EDF+
  // same logic/code as for SEGMENTS command (edf/dumper.cpp)

  std::set<interval_t> segs;
  
  // we only need to consider this for discontinuous EDF+ 
  
  if ( edf->header.continuous || ! edf->header.edfplus )
    {
      segs.insert( interval_t( 0 , total_duration_tp ) );
      return segs;
    }

  
  // need to query the time-tracks
    
  int r = first_record();
  
  uint64_t tp0 = rec2tp[r];
  
  uint64_t tp_start = tp0;

  while ( r != -1 )
    {
      
      // next record
      r = next_record( r );
      
      // start of this next record
      uint64_t tp;
      
      bool segend = false;
      
      // end?
      if ( r == -1 )
	{
	  // make this the 'previous'
	  tp0 = tp;
	  segend = true;
	}
      else
	{
	  tp = rec2tp[r] ;
	  
	  // discontinuity / end of segment?
	  // allow for minor round, must be within
	  // 1/10,000th of a second
	  // 0.0001 * 1e9 =  1e+05 tps
	  
          uint64_t len = tp - tp0;
	  uint64_t dif = len > edf->header.record_duration_tp ? len - edf->header.record_duration_tp : edf->header.record_duration_tp - len ;
          segend = dif > 10000LLU;
	  
	}
      
      // record this segment 
      
      if ( segend )
	{
	  
	  double secs1 = tp_start * globals::tp_duration ; 	  
	  double secs2 = tp0 * globals::tp_duration + edf->header.record_duration; 

	  // start = tp_start
	  // end   = tp0 + header.record_duration_tp  (i.e. up to and +1 past end of record)
	  segs.insert( interval_t( tp_start , tp0 + edf->header.record_duration_tp ) );
	  
	  // current point becomes start of the next segment
	  tp_start = tp;
	  
	}
      
      // current point becomes the last one, for next lookup
      tp0 = tp;
      
    }

  return segs;
  
}


void timeline_t::set_annot_metadata( const param_t & param )
{

  //  annot  : annotations to add MD to
  //  md     : key name of meta-data
  //  w     : specify a window +/- x seconds around each
  //  w-left / w-right : window only to left or right
  
  //  sig   : name of signal(s)
  //          implies signal mode;
  //          should be a single signal 
  //                   
  
  //  other : other annots
  //          cannot include self
  //          can be multiple (in which case, all annots are effectively pooled)
  
  //
  // functions: 
  //   signal-mode functions: mean, min, max, range, sd
  //     
  //   annot-mode functions:
  //        overlap (0/1) : is this key annot spanned by 1+ other annot?
  //        complete-overlap (0/1)  : is key annot completely spanned by (1+) of other annots?
  //        count (N) : count number of instances of all others
  //        nearest : distance to nearestN (in seconds) [ 0 ]
  //        nearest-midpoint : similar, but based on annot midpoints
  
  
  //
  // annot(s) to add MD to 
  //

  if ( ! param.has( "annot" ) ) Helper::halt( "no annotations specified: e.g. annot=A1,A2" );
  std::vector<std::string> anames = param.strvector_xsigs( "annot" );

	  
  //
  // flanking windows
  //

  const bool flanking = param.has( "w" );
  const double flanking_tp = flanking ? param.requires_dbl( "w" ) * globals::tp_1sec : 0 ;

  const bool left_flanking = param.has( "w-left" );
  const double left_flanking_tp = left_flanking ? param.requires_dbl( "w-left" ) * globals::tp_1sec : 0 ;

  const bool right_flanking = param.has( "w-right" );
  const double right_flanking_tp = right_flanking ? param.requires_dbl( "w-right" ) * globals::tp_1sec : 0 ;

  if ( flanking && ( left_flanking || right_flanking ) )
    Helper::halt( "cannot specify both 'w' and 'left-w' or 'right-w'" );

  //
  // flatten other annots (i.e. to define complete-overlap, etc,  we ignore that we have N2-N2 epoch boundary
  //   but for other contexts, e.g. count, we need to have option to not flatten
  //

  const bool flatten = param.has( "flatten" );
  
  //
  // signals
  //
  
  const bool no_annotations = true;
  signal_list_t signals = edf->header.signal_list( param.value( "sig" ) , no_annotations );  
  const int ns = signals.size() ;
  // if ns > 1 then the channel label gets added to mdtag;

  const bool signal_mode = param.value("sig") != "*" && ns > 0 ; 
  
  if ( signal_mode )
    {
      if ( ns == 0 ) Helper::halt( "no valid signals found" );
      const int Fs = edf->header.sampling_freq( signals(0) );      
      for (int s=1; s<ns; s++)
	if ( edf->header.sampling_freq( signals(s) ) != Fs )
	  Helper::halt( "signals must have similar sampling rates" );
    }

  //
  // special case: add metadata for duration of event
  //

  const bool dur_mode = param.has( "dur" ) ;
  
  //
  // other annots
  //

  const bool other_mode = param.has( "other" ) ;
  
  std::vector<std::string> oanames;
  if ( other_mode ) 
    oanames = param.strvector_xsigs( "other" );      

  
  //
  // can only be in a single mode
  //

  if ( dur_mode + other_mode + signal_mode != 1 )
    Helper::halt( "exactly one of 'other', 'sig' or 'dur' must be specified" );

  //
  // functions
  //

  std::string fn = "";
  
  if ( signal_mode )
    {
      int fc = 0;
      if ( param.has( "mean" ) ) { fn = "mean"; ++fc; } 
      if ( param.has( "min" ) ) { fn = "min"; ++fc; } 
      if ( param.has( "max" ) ) { fn = "max"; ++fc; }
      if ( param.has( "range" ) ) { fn = "range"; ++fc; }
      if ( fc != 1 ) Helper::halt( "must specify exactly one of: mean, min, max, range" );
    }
  else if ( other_mode )
    {
      int fc = 0; 
      if ( param.has( "overlap" ) ) {fn = "overlap"; ++fc; }
      if ( param.has( "complete-overlap" ) ) {fn = "complete-overlap"; ++fc; } // A is completely spanned by O
      if ( param.has( "whole-other" ) ) {fn = "whole-other" ; ++fc; } // O is completely spanned by A 
      if ( param.has( "count" ) ) {fn = "count"; ++fc; }
      if ( param.has( "nearest" )) {fn = "nearest"; ++fc; }
      if ( param.has( "nearest-start" )) {fn = "nearest-start"; ++fc; }
      if ( param.has( "nearest-stop" )) {fn = "nearest-stop"; ++fc; }
      if ( param.has( "nearest-midpoint" )) {fn = "nearest-midpoint"; ++fc; }

      if ( fc != 1 ) Helper::halt( "must specify exactly one of: count, overlap, complete-overlap, whole-other, nearest, nearest-midpoint, nearest-start, nearest-stop" ); 
    }


  //
  // search window for nearest comparisons (1 minute by default) 
  //

  double nearest_search_sec = 60;

  if ( fn == "nearest" && ! param.empty( "nearest" ) )
    nearest_search_sec = param.requires_dbl( "nearest" );

  if ( fn == "nearest-start" && ! param.empty( "nearest-start" ) )
    nearest_search_sec = param.requires_dbl( "nearest-start" );
  
  if ( fn == "nearest-stop" && ! param.empty( "nearest-stop" ) )
    nearest_search_sec = param.requires_dbl( "nearest-stop" );

  if ( fn == "nearest-midpoint" && ! param.empty( "nearest-midpoint" ) )
    nearest_search_sec = param.requires_dbl( "nearest-midpoint" );
  
  if ( nearest_search_sec < 0 ) nearest_search_sec = fabs( nearest_search_sec );
  
  uint64_t nearest_search_tp = nearest_search_sec * globals::tp_1sec;

  const bool nearest_mode = fn.substr(0,7) == "nearest";

  if ( nearest_mode )
    logger << "  using " << fn << " search window of " << nearest_search_sec << " seconds\n";

  //
  // MD tag 
  //

  std::string mdtag = param.requires( "md" );
  
  //
  // general other table (for nearest functions)
  //
  
  std::map<interval_t,std::string> allevs;  

  if ( nearest_mode )
    {
      for (int a=0; a<oanames.size(); a++)
	{      
	  annot_t * a1 = annotations->find( oanames[a] );
	  if ( a1 == NULL ) continue;
	  
	  // get /all/ annotations
	  const annot_map_t & events = a1->interval_events;
	  
	  // list events in this epoch (any span)
	  annot_map_t::const_iterator ii = events.begin();
	  while ( ii != events.end() )
	    {	  	      
	      const instance_idx_t & instance_idx = ii->first;	      
	      const interval_t & an_interval = instance_idx.interval;	  
	      if ( fn == "nearest-midpoint" ) 
		allevs[ interval_t( an_interval.mid() , an_interval.mid() ) ] = oanames[a];
	      else if ( fn == "nearest-start" )
		allevs[ interval_t( an_interval.start , an_interval.start ) ] = oanames[a];
	      else if ( fn == "nearest-stop" )
		allevs[ interval_t( an_interval.stop , an_interval.stop ) ] = oanames[a];
	      else if ( fn == "nearest" )
		allevs[ an_interval ] = oanames[a];
	      ++ii;
	    }
	}
      
      logger << "  built a table of " << allevs.size() << " other events for nearest lookups\n";
      
    }

  
  //
  // iterate over annots
  //

  for (int a=0; a<anames.size(); a++)
    {
      
      // does annot exist?
      annot_t * annot = edf->annotations->find( anames[a] );
      if ( annot == NULL ) continue;
      const std::string & class_name = anames[a];
      
      // get all events
      const annot_map_t & events = annot->interval_events;

      annot_map_t::const_iterator aa = events.begin();
      while ( aa != events.end() )
	{

	  const instance_idx_t & idx = aa->first;
	  instance_t * instance = aa->second;

	  // copy main interval
	  interval_t interval = idx.interval;


	  //
	  // special case: dur 
	  //

	  if ( dur_mode )
	    {
	      // use original interval, i.e. this is done pre-expansions
	      instance->set( mdtag , interval.duration_sec() );

	      // all done now --> to next event
	      ++aa;
	      continue;
	    }
	

	  //
	  // expand?
	  //
	  
	  if ( flanking ) interval.expand( flanking_tp );
	  else if ( left_flanking ) interval.expand_left( left_flanking_tp );
	  else if ( right_flanking ) interval.expand_right( right_flanking_tp );

	  //
	  // signal mode
	  //
	  
	  if ( signal_mode )
	    {	      	      
	      eigen_matslice_t mslice( *edf , signals , interval );	  
	      
	      const Eigen::MatrixXd & X = mslice.data_ref();
	      
	      Eigen::ArrayXd stat;

	      if ( fn == "mean" ) stat = X.array().colwise().mean();
	      else if ( fn == "min" ) stat = X.array().colwise().minCoeff();
	      else if ( fn == "max" ) stat = X.array().colwise().maxCoeff();	  
	      else if ( fn == "range" ) stat = X.array().colwise().maxCoeff() - X.array().colwise().minCoeff();
	      
	      // add as meta-data
	      if ( ns == 1 ) instance->set( mdtag , stat[0] );
	      else for (int s=0; s<ns; s++) instance->set( mdtag + "_" + signals.label(s) , stat[s] );
	    }

	  
	  //
	  // annot-mode
	  //

	  if ( other_mode )
	    {
	      	      
	      
	      //
	      // count/overlap/etc
	      //
	      
	      if ( ! nearest_mode )
		{
	      
		  // build up new, epoch-based annotation map
		  std::set<interval_t> nevs;
		  
		  for (int a=0; a<oanames.size(); a++)
		    {
		      annot_t * a1 = annotations->find( oanames[a] );
		      
		      if ( a1 == NULL ) continue;
	      
		      // get overlapping annotations (spanning this this window) 
		      annot_map_t events = a1->extract( interval );
		      
		      // list events in this epoch (any span)
		      annot_map_t::const_iterator ii = events.begin();
		      while ( ii != events.end() )
			{	  	      
			  const instance_idx_t & instance_idx = ii->first;	      
			  const instance_t * instance = ii->second;
			  const interval_t & an_interval = instance_idx.interval;
			  
			  //std::cout << " an_interval = " << instance_idx.id << "\t" << an_interval.as_string() << "\n";
			  
			  // keep intersection w/ this spanning epoch	      
			  nevs.insert( an_interval );
			  ++ii;
			}
		      
		    }

		  //
		  // flatten other events? (i.e. for better def. of complete-overlap and whole-other ) 
		  //
				
		  if ( flatten )
		    nevs = annotate_t::flatten( nevs );
		  
		  
		  //
		  // we now have a list of all potential events in nevs;  do processing
		  //
		  
		  if ( fn == "count" )
		    {
		      // number of other annots overlapping this 
		      instance->set( mdtag , (int)nevs.size() );
		    }
		  else if ( fn == "overlap" )
		    {
		      // similar to count, but reduced to 0/1
		      instance->set( mdtag , (int)( nevs.size() != 0 ) );
		    }
		  else if ( fn == "complete-overlap" )
		    {
		      // 0 vs 1 : is A completely spanned by at least one O?
		      int x = 0;
		      std::set<interval_t>::const_iterator ii = nevs.begin();
		      while ( ii != nevs.end() )
			{
			  if ( interval.is_completely_spanned_by( *ii ) )
			    {
			      x = 1; 
			      break;
			    }
			  ++ii;
			}
		      instance->set( mdtag ,x );
		    }
		  else if ( fn == "whole-other" )
		    {
		      // 0 vs 1 : is at least one O completely spanned by
		      int x = 0;
		      std::set<interval_t>::const_iterator ii = nevs.begin();
		      while ( ii != nevs.end() )
			{
			  if ( ii->is_completely_spanned_by( interval ) )
			    {
			      x = 1; 
			      break;
			    }
			  ++ii;
			}
		      instance->set( mdtag ,x );
		    }
		}

	      
	      //
	      // nearest search (uses allevs, not nevs)
	      //
	      
	      if ( nearest_mode )
		{

		  // default nearest:  0 if overlaps, else  STOP --> START (-ve)   or STOP-->START (+ve)
		  //   others: point-based , as is
		  
		  // do search: set target 
		  interval_t nidx = interval;
		  if ( fn == "nearest-midpoint" )
		    nidx.start = nidx.stop = interval.mid();
		  else if ( fn == "nearest-start" )
		    nidx.start = nidx.stop = interval.start;
		  else if ( fn == "nearest-stop" )
		    nidx.start = nidx.stop = interval.stop;

		  //		  std::cout << " nidx = " << nidx.as_string() << "\n";
		  // do general lookup:: first event /after/ the query (upper-bound)
		  //  (and above, the allevs set will have been constructed w/ the a
		  //   appropriate type (whole, midpoint, start or stop)
		  std::map<interval_t,std::string>::const_iterator bb = allevs.upper_bound( nidx );
		  
		  // we might not find /any/ other events
		  bool any = false;
		  bool first_comp = true;
		  uint64_t dist = 0;
		  int dsign = 0;
		  std::string matched = ".";
		  const bool full_mode = fn == "nearest" ;
		
		  //		  std::cout << " going into loop " << allevs.size() << "\n";
		  
		  if ( allevs.size() )
		    {
		      
		      while ( 1 )
			{
			  
			  // past end?
			  if ( bb == allevs.end() )
			    {
			      //std::cout << " hit the end...\n";
			      --bb;
			      continue;
			    }
			  
			  // std::cout << " considering " << bb->first.as_string() << " " << bb->second << "\n";
			  
			  // found at least one contender?
			  any = true ;
			  
			  if ( full_mode )
			    {
			      
			      // don't think we'll have any actual overlap here (would have been
			      // caught above... but doen't hurt to add here in any case)
			      
			      const bool before = bb->first.stop  <= nidx.start ;
			      const bool after  = bb->first.start >= nidx.stop ;
			      
			      uint64_t d1 = before ? nidx.start - bb->first.stop :
				( after ? bb->first.start - nidx.stop : 0 ) ;
			      
			      if ( d1 < dist || first_comp )
				{
				  dist = d1;
				  dsign = before ? -1 : ( after ? 1 : 0 ) ;
				  matched = bb->second;
				  first_comp = false;
				}
			      
			    }
			  else
			    {
			      // can use start in all cases, all 0-tp points
			      const bool after = bb->first.start > nidx.start ;
			      uint64_t d1 = after ? bb->first.start - nidx.start : nidx.start - bb->first.start ;
			      if ( d1 < dist || first_comp )
				{				  
				  dist = d1;
				  dsign = d1 == 0 ? 0 : ( after ? 1 : -1 ) ;
				  matched = bb->second;
				  first_comp = false;
				}
			    }

			  // best possible match?
			  if ( dist == 0 ) break;
			      
			  // all done? (nb. the window means the /whole/ event must be in the window
			  //  i.e. as we may have the case of something starting a long, long time before
			  //       but persisting right up to the actual event...

			  if ( nidx.start > bb->first.start && nidx.start - bb->first.start > nearest_search_tp ) break;
			  if ( nidx.start < bb->first.start &&  bb->first.start - nidx.start > nearest_search_tp ) break;
			  
			  // all done?
			  if ( bb == allevs.begin() ) break;
			  
			  // step back in time
			  --bb;
			}
		      
		      // check that the nearest event matches criteria (i.e. if was upper bound
		      // but that came a long time after , and we had no other matches)
		      
		      if ( dist > nearest_search_tp ) any = false;
			  
		      // report
		      if ( any )
			{
			  double sec = dist / (double)globals::tp_1sec;
			  instance->set( mdtag , dsign * sec );
			  instance->set( mdtag + "_id"  , matched );
			}
		      else 
			{
			  //instance->set( mdtag , "NA" ); // missing
			  std::string lab1 = "."; // use so correct type w/ set()
			  instance->set( mdtag + "_id"  , lab1 );
			}
		    }
		} // end of nearest mode
	    }  // end of other mode
	  
	  // next event
	  ++aa;
	}
      
    } // next annotation class
    
  // all done
}




//
// Implements AXA 
//

void timeline_t::annot_crosstabs( const param_t & param )
{

  // for root-match
  std::vector<std::string> names = annotations->names();
  
  // get list of annotations
  std::vector<std::string> requested = param.has( "annot" ) && param.value( "annot" ) != "."
    ? param.strvector_xsigs( "annot" )
    : names ; 
  
  // allow root-matching
  //  requested = Helper::set2vec( annotate_t::root_match( Helper::vec2set( requested ) , names ) );
  
  // group annots by class only, or also by instance IDs?
  const bool by_instance = param.has( "by-instance" );
  
  // only consider cls-cls comparisons where the instance ID matches
  //   e.g. cls = ISO(band)  inst = ISO(phasebin)
  const bool match_instance = param.has( "match-instance" );
  if ( match_instance && ! by_instance ) Helper::halt( "match-instance requires by-instance is set" );

  
  // flatten events?
  const bool flatten = param.has( "flatten" ) ? param.yesno( "flatten" ) : false;

  // event-level output?
  const bool verbose = param.has( "verbose" ) ? param.yesno( "verbose" ) : false;

  // within channel only? 
  //  note: some annotations might not have a channel specified,
  //        in which case the catch all "." counts as a 'channel'
  const bool within_channel = param.has( "within-channel" );

  // anchor (-1,0,+1) for start, mid, stop
  const int anchor = param.has( "start" ) ? -1 : ( param.has( "stop" ) ? +1 : 0 ) ; 
  
  // time limit for match  ( neg means no window) 
  const double window = param.has( "w" ) ? param.requires_dbl( "w" ) : -1 ; 
  
  // count implied annotations 
  if ( requested.size() == 0 )
    Helper::halt( "no annotations" );
  
  //
  // build up table of events
  //

  // ch -> annot -> interval lists
  std::map<std::string,std::map<std::string,std::set<interval_t> > > events;
  
  // track annot instance ID (for match-instance) 
  std::map<std::string,std::string> label2instance;


  //
  // iterate over each annotation
  //
  
  for (int a = 0 ; a < requested.size() ; a++ ) 
    {
      
      annot_t * annot = annotations->find( requested[a] );
      
      if ( annot == NULL ) continue;
      
      const int num_events = annot->num_interval_events();
      
      logger << "  found " << num_events << " instances of " << requested[a] << "\n";
      
      //
      // iterator over interval/event map
      //
            
      const std::string label = requested[a];
      
      annot_map_t::const_iterator ii = annot->interval_events.begin();
      while ( ii != annot->interval_events.end() )
	{	  
	  
	  const instance_idx_t & instance_idx = ii->first;
	  
	  // add to the list

	  const std::string ch_str = within_channel ? instance_idx.ch_str : "." ;
	  const std::string label1  = label + ( by_instance ? "_" + instance_idx.id : "" ); 
	  if ( match_instance ) label2instance[ label1 ] = instance_idx.id; 

	  // record
	  events[ ch_str ][ label1 ].insert( instance_idx.interval );
	  
	  ++ii;
	}
      
    }

  //
  // Flatten all events first?
  //  - if within channel, then flattening only happens w/in channels too
  //

  if ( flatten )
    {
      std::map<std::string,std::map<std::string,std::set<interval_t> > >::iterator aa = events.begin();
      while ( aa != events.end() )
	{
	  std::map<std::string,std::set<interval_t> > & events1 = aa->second;
	  std::map<std::string,std::set<interval_t> >::iterator ee = events1.begin();
	  while ( ee != events1.end() )
	    {	      
	      int n1 = ee->second.size();
	      ee->second = annotate_t::flatten( ee->second );
	      int n2 = ee->second.size();
	      if ( n2 < n1 ) logger << "  reduced " << ee->first << " from " << n1 << " to " << n2 << " events\n";
	      ++ee;
	    }
	  ++aa;
	}
    }
  
  
  //
  // Over each channel
  //

  std::map<std::string,std::map<std::string,std::set<interval_t> > >::const_iterator ch = events.begin();
  while ( ch != events.end() )
    {
      
      // track channel?

      if ( within_channel )
	writer.level( ch->first , globals::signal_strat );


      //
      // Pull out annotations (for this channel) 
      //
      
      const std::map<std::string,std::set<interval_t> > & events1 = ch->second;
      
      //
      // Nothing to do?
      //
      
      if ( events1.size() < 2 )
	{
	  logger << "  *** nothing to do, fewer than two annotation classes found";
	  if ( within_channel ) logger << " for channel " << ch->first ;
	  logger << "\n";
	}
            
      //
      // Consider all pairs of events
      //

      std::map<std::string,std::set<interval_t> >::const_iterator bb = events1.begin();
      while ( bb != events1.end() )
	{
	  
	  writer.level( bb->first , globals::annot_strat );

	  // we want to keep 'a' 'as is' (i.e. might be flattened, but might not be
	  //  but for 'b', here we should flatten to make the looks easier., for the
	  //  overlap analyses;  for the nearest analysis, keep the original set,
	  //  and extract the anchor point here:

	  // for 'nearest' distance
	  std::set<double> b1;
	  std::set<interval_t>::const_iterator bb1 = bb->second.begin();
	  while ( bb1 != bb->second.end() )
	    {
	      if      ( anchor == -1 ) b1.insert( bb1->start_sec() );	      
	      else if ( anchor == +1 ) b1.insert( bb1->stop_sec() );
	      else                     b1.insert( bb1->mid_sec() );
	      ++bb1;
	    }

	  // flatten remaining values
	  const std::set<interval_t> b = annotate_t::flatten( bb->second );
	  
	  std::map<std::string,std::set<interval_t> >::const_iterator aa = events1.begin();
	  while ( aa != events1.end() )
	    {
	      
	      // does instance match?
	      if ( match_instance ) 
		{
		  // skip if not
		  if ( label2instance[ aa->first ] != label2instance[ bb->first ] )
		    {
		      ++aa;
		      continue;
		    }
		}
	      
	      writer.level( aa->first , "SEED" );
	      
	      const std::set<interval_t> & a = aa->second;

	      // std::cout << "\n\n------------------------------"
	      //  		<< aa->first << " .... " << bb->first << "\n";
	      
	      // to track outputs
	      std::vector<double> sav_p, sav_t, sav_n, sav_a;
	      
	      // nearest (within w) [ will only contain instances that were within range ] 
	      std::vector<double> sav_d, sav_dabs;
	      
	      // determine for lists a and b :
	      //   total % of overlap per all 'a'
	      //   mean % of overlap per 'a' event
	      //   num_secs of 'b' given 'a'
	      //   time from seed 'anchor' to nearest 'b' 'anchor' (signed + abs) 
	      //   % of 'a' with at least some 'b' overlap
	      
	      int sidx = 0;

	      // for each 'seed' (i.e. conditioning event)
	      std::set<interval_t>::const_iterator seed = a.begin();
	      while ( seed != a.end() )
		{
		  
		  ++sidx;

		  // skip zero-duration seeds here
		  if ( fabs( seed->duration_sec() ) < 1e-8 )
		    {
		      ++seed;
		      continue;
		    }
		  
		  // which b events span this, if any?		  
		  // find the first annot not before (at or after) the seed

		  //
		  // Find nearest distance (using the unflattened b1 list of anchor positions)
		  //

		  const double seed_sec = anchor == -1 ? seed->start_sec()
		    : ( anchor == 1 ? seed->stop_sec() : seed->mid_sec() );
		  
		  std::set<double>::const_iterator cc = b1.upper_bound( seed_sec );
		  
		  std::set<double>::const_iterator closest = cc;
		  
		  std::vector<double> distances; 
		  
		  while ( 1 )
		    {
		      if ( closest == b1.end() ) break;
		      if ( *closest > seed_sec ) break;
		      ++closest;
		    }

		  // one past
		  if ( closest != b1.end() )
		    {
		      //std::cout << " dst: back-track , considering " << *closest << "\n";  
		      distances.push_back( *closest - seed_sec );
		    }
		  
		  // now count back
		  while ( 1 )
		    {		      
		      
		      if ( closest == b1.begin() ) break;		  
		      --closest;
		      
		      //std::cout << " dst: back-track , considering " << *closest << "\n";
		      
		      distances.push_back( *closest - seed_sec );
		      
		      if ( *closest < seed_sec )
			{
			  //std::cout << " dst: gone past , breakingh\n";
			  break;
			}
		      
		    }

		  // get min distance (Could have done above, but whatev)

		  if ( window > 0 && distances.size() > 0 ) 
		    {
		      
		      //std::cout << sidx << "  DIST " << bb->first << " " << aa->first  << " "  << distances.size() << "\n ";
		      
		      double q1 = window + 10000;
		      double q  = 0;
		      int okay = 0;
		      
		      for (int i=0; i<distances.size(); i++)
			{

			  const double d1 = fabs( distances[i] );

			  //std::cout << " candidate = " << d1 << " " << distances[i] << "\n";

			  if ( d1 < window )
			    {
			      ++okay;

			      //std::cout << " con " << distances[i]  << " " << d1 << "\n";
			      
			      if ( d1 < q1 ) 
				{
				  q1 = d1;
				  q = distances[i];
				}
			    }
			}

		      // if ( okay ) 
		      // 	std::cout << " dst: *** saving " << q << "\n";
		      // else
		      // 	std::cout << " dst:  no within-window match found...\n";
		      
		      // save signed closest annot		      
		      if ( okay )
			{
			  sav_d.push_back( q );
			  sav_dabs.push_back( fabs( q ) );
			}
		    }
		  
		  
		  //
		  // Find overlaps (using the flattened b)
		  //
		  
		  std::set<interval_t>::const_iterator oclosest = b.upper_bound( *seed );

		  std::set<interval_t> overlaps;
		  
		  while ( 1 )
		    {
		      if ( oclosest == b.end() ) break;
		      if ( oclosest->start >= seed->stop ) break;
		      ++oclosest;
		    }
		  
		  // now count back
		  while ( 1 )
		    {		      
		      
		      if ( oclosest == b.begin() ) break;		  
		      --oclosest;
		      
		      //std::cout << " olap: back-track , considering " << oclosest->as_string() << "\n";
		      
		      if ( oclosest->stop <= seed->start )
			{
			  //std::cout << " olpa: gone past , breakingh\n";
			  break;
			}
		      
		      interval_t o( oclosest->start > seed->start ? oclosest->start : seed->start ,
				    oclosest->stop < seed->stop  ? oclosest->stop : seed->stop );
		      
		      overlaps.insert( o );		  
		    }
		  
		  const int n_olap = overlaps.size();
		  
		  double t_olap = 0;
		  
		  std::set<interval_t>::const_iterator oo = overlaps.begin();
		  while ( oo != overlaps.end() )
		    {
		      t_olap += oo->duration_sec();
		      ++oo;
		    }
		  
		  const double p_olap = t_olap / seed->duration_sec();
		  
		  sav_n.push_back( n_olap );
		  sav_t.push_back( t_olap );
		  sav_p.push_back( p_olap );
		  sav_a.push_back( n_olap > 0 );
		  
		  
		  //
		  // Verbose output?
		  //
		  
		  if ( verbose )
		    {
		      std::cout << " olap = " << bb->first << " " << aa->first << " = " << overlaps.size() << "\n";
		      
		    }
		  
		  ++seed;
		}
	      
	      //
	      // Now summarize outputs
	      //
	      
	      // mean (or median) per seed interval
	      writer.value( "P" , MiscMath::mean( sav_p ) );
	      //writer.value( "P_MD" , MiscMath::median( sav_p ) );
	      
	      writer.value( "T" , MiscMath::mean( sav_t ) );
	      //writer.value( "T_MD" , MiscMath::median( sav_t ) );
	      
	      writer.value( "N" , MiscMath::mean( sav_n ) );
	      //writer.value( "N_MD" , MiscMath::median( sav_n ) );
	      
	      if ( sav_d.size() > 0 )
		{
		  writer.value( "D" , MiscMath::mean( sav_d ) );
		  writer.value( "DABS" , MiscMath::mean( sav_dabs ) );
		}
	      writer.value( "D_N", (int)sav_d.size() );

	      writer.value( "A" , MiscMath::mean( sav_a ) );
	      
	      // grand totals
	      writer.value( "TOT_N" , MiscMath::sum( sav_n ) );
	      writer.value( "TOT_T" , MiscMath::sum( sav_t ) );
	      
	      ++aa;
	    }
	  writer.unlevel( "SEED" );
	  
	  ++bb;
	}
      writer.unlevel( globals::annot_strat );

      ++ch;
    }

  if ( within_channel )
    writer.unlevel( globals::signal_strat );
}



uint64_t s2a_interp( const std::vector<double> * d ,     // angles (0-360)
		     const std::vector<uint64_t> * tp ,  // time-points
		     const int i ,  // first point *after* target T
		     const double t ) // target
{

  // typically: e.g. if t = 90
  //     i-1   89.7
  //     T     90.0
  //     i     90.2
  // i.e. this spans the range
  
  // however,

  // scenario1: possible at 0 we wrapped around phase: (e.g. T = 0)
  //   i-1    359.8
  //   T        0.0
  //   i        0.2
  
  // scenario2: unlikely, but possible target is high
  //   i-1    349.8
  //   T      350.0
  //   i      2.1

  // so, test whether 
  const bool scenario1 = (*d)[i-1] > t ;
  const bool scenario2 = (*d)[i] < t ; 
  
  // scenario1 : reduce first by 360
  // scenario2 : increase second by 360

  if ( scenario1 && scenario2 ) // should not happen
    Helper::halt( "internal logic error in s2a_interp()" );
  
  const double first  = scenario1 ? (*d)[i-1] - 360.0 : (*d)[i-1] ;
  const double second = scenario2 ? (*d)[i] + 360.0 : (*d)[i] ;

  if ( first > t || first > second || t > second )
    Helper::halt( "internal logic error in s2a_interp()" );

  // get proportional scaling
  const double p = ( t - first ) / ( second - first );

  // return an adjusted tp (i.e. no longer aligned with a sample point per se)
  const uint64_t d1 = (*tp)[i] - (*tp)[i-1];
  const uint64_t retval = (*tp)[i-1] + p * d1;
  return retval;
}
