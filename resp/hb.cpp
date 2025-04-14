
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


// Functions in this file implement Ali Azarbarzin's Hypoxic Burden method
// Reference: 

#include "resp/hb.h"

#include "helper/helper.h"
#include "helper/logger.h"
#include "db/db.h"
#include "edf/edf.h"
#include "edf/slice.h"
#include "dsp/resample.h"
#include "stats/eigen_ops.h"

#include <deque>

extern writer_t writer;
extern logger_t logger;

hb_t::hb_t( edf_t & edf , param_t & param )
{

  hb_results_t res;

  //
  // Parameters
  //
  
  const double sao2_threshold = param.has( "th" ) ? param.requires_dbl( "th" ) : 40;

  const std::string oxy_label = param.requires( "oxygen" );

  const std::string hr_label = param.requires( "hr" );

  const bool remove_wake_events = param.has( "no-wake" );

  
  //
  // Constants
  //
   
  const int dT = 1;

  const double rangeX = 100;

  const double MagAv_thres = 1.5;

  
  //
  // SaO2 
  //  

  int oxy_n = edf.header.signal( oxy_label );  
  if ( oxy_n == -1 )
    {
      logger << "  could not find oxygen desaturation signal " << oxy_label << "\n";
      res.valid = false;
      return;
    }

  
  //
  // HR
  //
  
  int hr_n = edf.header.signal( hr_label );  
  if ( hr_n == -1 )
    {
      logger << "  could not find heart-rate signal " << hr_label << "\n";
      res.valid = false;
      return;
    }

  
  //
  // Fs: use oxygen signal SR unless otherwise specified:
  //
  
  const int fs = param.has( "sr" ) 
    ? param.requires_int( "sr" )
    : edf.header.n_samples[ oxy_n ] / edf.header.record_duration; 

  const uint64_t fs_tp = globals::tp_1sec / fs;
  
    
  // 
  // Resample as needed (e.g. 32 Hz in example data); use ZOH resampling by default
  //

  int converter = param.has( "method" ) ? dsptools::converter( param.value( "method" ) ) : dsptools::converter( "ZOH" );
  
  

  if ( edf.header.sampling_freq( oxy_n ) != fs ) 
    {
      logger << "  resampling oxygen channel using method '" << dsptools::converter( converter ) << "'\n"; 
      dsptools::resample_channel( edf, oxy_n  , fs , converter );
    }
  
  if ( edf.header.sampling_freq( hr_n ) != fs )
    {
      logger << "  resampling HR channel using method '" << dsptools::converter( converter ) << "'\n";
      dsptools::resample_channel( edf, hr_n  , fs , converter );
    }
  
  //
  // Annotations
  //

  const std::string annot_arousal   = param.has( "arousal" )   ? param.value( "arousal" )   : "arousal_standard" ;
  
  std::vector<std::string> event_labels;

  if ( param.has( "events" ) )
    event_labels = param.strvector( "events" );
  else
    {
      const std::string annot_obs_ap    = param.has( "apnea-obs" ) ? param.value( "apnea-ap" )  : "apnea_obstructive" ;
      const std::string annot_obs_cen   = param.has( "apnea-cen" ) ? param.value( "apnea-cen" ) : "apnea_central";
      const std::string annot_obs_mixed = param.has( "apnea-mix" ) ? param.value( "apeaa-mix" ) : "apnea_mixed";
      const std::string annot_hypop_50  = param.has( "hypopnea-50" ) ? param.value( "hypopnea-50" ) : "hypopnea";
      const std::string annot_hypop_30  = param.has( "hypopnea-30" ) ? param.value( "hypopnea-30" ) : "hypopnea";
      
      event_labels.push_back( annot_obs_ap );
      event_labels.push_back( annot_obs_cen );
      event_labels.push_back( annot_obs_mixed );
      event_labels.push_back( annot_hypop_50 );
      event_labels.push_back( annot_hypop_30 );

    }
  
  
  //
  // Extract sleep staging
  //

  edf.annotations->make_sleep_stage( edf.timeline );
      
  if ( ! edf.timeline.hypnogram.construct( &edf.timeline , param , false ) )
    Helper::halt( "problem extracting stage annotations" );

  const std::vector<sleep_stage_t> & stages = edf.timeline.hypnogram.stages;
  
  const int n_epochs = stages.size();

  // 0=R, 1,2,3,4=NR, 5=R

  const int npe = edf.timeline.epoch_length() * fs;
  const int np  = n_epochs * npe;
  
  std::vector<int> ss;
  
  for (int e = 0; e < n_epochs ; e++ )
    {
      // nb. includes unscored, unknown, movement and LightsOn w/ 'wake'
      int s = 0;  
      if      ( stages[e] == NREM1 ) s = 1;
      else if ( stages[e] == NREM2 ) s = 2;
      else if ( stages[e] == NREM3 ) s = 3;
      else if ( stages[e] == NREM4 ) s = 4;
      else if ( stages[e] == REM )   s = 5;
      for (int i=0; i<npe; i++) ss.push_back( s );
    }
  

  //
  // Expand arousal annotation into a 0/1 binary sample-level vector
  //
  
  annot_t * annot = edf.annotations->find( annot_arousal );

  std::vector<bool> arousals( np , false );

  int a_cnt = 0;
  double a_dur = 0;
  int aa_cnt = 0;
  if ( annot != NULL )
    {
      annot_map_t ars = annot->extract( edf.timeline.wholetrace() );
     	
      annot_map_t::const_iterator aa = ars.begin();
      while ( aa != ars.end() )
	{
	  // get nearest sample points
	  int start = aa->first.interval.start / fs_tp ;
	  int stop  = aa->first.interval.stop / fs_tp ;

	  // track length
	  ++a_cnt;
	  a_dur += aa->first.interval.duration_sec();
	  aa_cnt += stop - start + 1 ; 
	  // goes past end?
	  if ( stop >= np ) stop = np - 1;
	  for (int p=start; p<=stop; p++) arousals[p] = true;
	  ++aa;
	}      
      
    }
  else
    logger << "  no arousal annotation track found\n";
  
  logger << "  " << a_cnt << " arousals found, spanning " << a_dur << " secs (" << a_dur/60.0 << " mins)\n";
  

  
  //
  // Annotations
  //

  // typedef std::map<instance_idx_t,instance_t*> annot_map_t;
  annot_map_t events;

  for (int e = 0 ; e < event_labels.size() ; e++ )
    {      
      annot_t * annot = edf.annotations->find( event_labels[e] );

      if ( annot == NULL ) continue;
			     
      // get all events (nb. probably an easier way to do this, can't recall canonical form)
      annot_map_t evts = annot->extract( edf.timeline.wholetrace() );
      
      annot_map_t::const_iterator aa = evts.begin();
      while ( aa != evts.end() )
	{
	  events.insert( *aa );
	  ++aa;
	}      
      
    }

  //
  // Track number of events
  //

  int ne = events.size();

  logger << "  " << ne << " matching event annotations found\n";
  
  if ( 0 )
    {
      annot_map_t::const_iterator aa = events.begin();
      while ( aa != events.end() )
	{
	  const instance_idx_t & idx = aa->first;
	  std::cerr << "evt " << idx.parent->name << "\t"
		    << idx.id << "\t"
		    << idx.interval.start << " - "
		    << idx.interval.stop << "\n";
	  ++aa;
	}
    }
  
  
  // Event types:
  // OAp        : obstructive apneas
  // CAp        : central apneas
  // MAp        : mixed apneas
  // H30 Hyp    : hypopneas
  // H50 Hyp50  : hypopnea >50% or Unsure events


  //
  // Signals: just pull out entire signals
  //
  
  slice_t slice_oxy( edf , oxy_n , edf.timeline.wholetrace() );

  std::vector<double> SaO2 = *slice_oxy.pdata();

  const int n = SaO2.size();
  
  // nearest-neighbour interpolation of low values
  if ( sao2_threshold > 0 )
    {
      
      for (int i=0; i<n; i++)
	{

	  if ( SaO2[i] >= sao2_threshold ) continue;

	  // search for nearest valid value(s) and fill in everything
	  int lwr_idx = i , upr_idx = i;
	  double lwr_val = SaO2[i];
	  double upr_val = SaO2[i];
	  bool lwr = false , upr = false;
	  
	  while ( 1 )
	    {
	      if ( lwr_idx == 0 ) break;
	      --lwr_idx;
	      if ( SaO2[ lwr_idx ] >= sao2_threshold )
		{
		  lwr_val = SaO2[ lwr_idx ];
		  ++lwr_idx;
		  lwr = true;
		  break;
		}
	    }

	  while ( 1 )
	    {
	      ++upr_idx;
	      if ( upr_idx == n ) { --upr_idx; break; }
	      if ( SaO2[ upr_idx ] >= sao2_threshold )
		{
		  upr_val = SaO2[ upr_idx ];
		  --upr_idx;
		  upr = true;
		  break;
		}	      
	    }

	  if ( ! ( lwr || upr ) )
	    Helper::halt( "bad oxygen channel: all sub-threshold" );

	  double imputed = 0;
	  if ( lwr ) imputed += lwr_val;
	  if ( upr ) imputed += upr_val;
	  if ( lwr && upr ) imputed /= 2.0;

	  for (int ii=lwr_idx; ii<=upr_idx; ii++)
	    SaO2[ii] = imputed;
	  
	  // advance to next point, given we've imputed this block
	  if ( upr ) i = upr_idx;
	}
    }

  
  //
  // Time signal
  //

  std::vector<double> seconds( n );

  const std::vector<uint64_t> * tp = slice_oxy.ptimepoints();

  for (int i=0; i<n; i++)
    seconds[i] = (double)(*tp)[i] * globals::tp_duration;

  
   //
   // HR signals
   //

   slice_t slice_hr( edf , hr_n , edf.timeline.wholetrace() );

   std::vector<double> HR = *slice_hr.pdata();

   
   //
   // Sleep time calculation (in minutes)
   //
   
   const double edur = edf.timeline.epoch_length() / 60.0 ;
   
   for (int e = 0; e < n_epochs ; e++ )
     {
       if      ( stages[e] == NREM1 ) { res.TST_N1 += edur; }
       else if ( stages[e] == NREM2 ) { res.TST_N2 += edur; }
       else if ( stages[e] == NREM3 ) { res.TST_N3 += edur; }
       else if ( stages[e] == NREM4 ) { res.TST_N3 += edur; } // collapse to N3
       else if ( stages[e] == REM )   { res.TST_REM += edur; }
       res.TIB += edur;
     }
   
   res.TST_NREM = res.TST_N1 + res.TST_N2 + res.TST_N3 ;
   res.TST = res.TST_NREM + res.TST_REM;
   
   logger << "  TST: " << res.TST << " (NREM: " << res.TST_NREM << ", REM: " << res.TST_REM << ")\n";
   
   
   //
   // Remove events that start in wake?
   //
  
   if ( remove_wake_events )
     {
       // filter event list (evtSt, evtEnd) to only consider events that happen during sleep
       
       // SleepStageEvent=interp1(TimeSig,SleepStageSig,evtSt,'nearest');
       // evtSt=evtSt(SleepStageEvent>0 & SleepStageEvent<6);
       // evtEnd=evtEnd(SleepStageEvent>0 & SleepStageEvent<6);
       
       annot_map_t all_events = events;
       events.clear();
       
       annot_map_t::const_iterator aa = all_events.begin();
       while ( aa != all_events.end() )
	 {
	   const instance_idx_t & idx = aa->first;
	   int start_sp = idx.interval.start / fs_tp ;
	   if ( start_sp >= np ) continue;
	   if ( ss[ start_sp ] > 0 && ss[ start_sp ] <= 5 ) events.insert( *aa );
	   ++aa;	  
	 }
       
       logger << "  subsetting to " << events.size() << " (of " << all_events.size() << ") events that start during sleep\n";
       
       ne = events.size();
     }
   
   //
   // Create window/time track (in seconds )
   //
  
   std::vector<double> t_sp;
   for (int t=-rangeX ; t<= rangeX ; t+= dT )
     t_sp.push_back(t);
   
   const int nt = t_sp.size();
  
   // Time_new=-rangeX:dT:rangeX;
   
   // rows events
   // cols = time in seconds +/- 100 seconds from the END of each event
   
   // Time_OtherSigs=repmat(Time_new,length(evtEnd),1)+repmat(evtEnd,1,length(Time_new));
   
   // for SaO2, SS, arousal & HR, construct [ event x time ] matrices,
   // which are synced to the event stop, +/- rangeX samples
   
   // pull out all the stops... ha ha
   
   std::vector<int> stop;  
   annot_map_t::const_iterator aa = events.begin();
   while ( aa != events.end() )
     {
       stop.push_back( aa->first.interval.stop / fs_tp );
       ++aa;
     }
   


   //
   // Event times (seconds)
   //
   
   std::vector<double> evtSt, evtEnd;  
   aa = events.begin();
   while ( aa != events.end() )
     {
       evtSt.push_back( aa->first.interval.start_sec() );
       evtEnd.push_back( aa->first.interval.stop_sec() );
       ++aa;
     }
   
   //
   // make windows
   //
   
   Eigen::MatrixXd SaO2Ve = Eigen::MatrixXd::Zero( ne , nt );
   Eigen::MatrixXi SleepStageVe = Eigen::MatrixXi::Zero( ne , nt );
   Eigen::MatrixXd HRVe = Eigen::MatrixXd::Zero( ne , nt );
   Eigen::MatrixXi ArousalVe = Eigen::MatrixXi::Zero( ne , nt );
   
   // checks
   // nb. SaO2 and HR may be slightly longer, as they may include a final, incomplete epoch
   if ( SaO2.size() < np ) Helper::halt( "internal discrepancy, SaO2 size" );
   if ( HR.size() < np ) Helper::halt( "internal discrepancy, HR size" );
   if ( ss.size() != np ) Helper::halt( "internal discrepancy, SS size" );
   if ( arousals.size() != np ) Helper::halt( "internal discrepancy, Arousal size" );
   
   for (int e=0; e<ne; e++)
     {
       int a = stop[e] - rangeX*fs;
       int b = stop[e] + rangeX*fs;
       int t0 = 0;
       for (int t=a; t<=b; t+=fs )
	 {
	   // ensure is within range (i.e. last obs carried forward/backward)
	   int t1 = t < 0 ? 0 : ( t >= np ? np-1 : t );
	   SaO2Ve(e,t0) = SaO2[ t1 ];
	   ArousalVe(e,t0) = arousals[ t1 ] ;
	   SleepStageVe(e,t0) = ss[ t1 ];
	   HRVe(e,t0) = HR[t1];
	   ++t0;
	 }
     }
   
   // ArousalVe=interp1(TimeSig,ArousalSig,Time_OtherSigs,'nearest');
   // SleepStageVe=interp1(TimeSig,SleepStageSig,Time_OtherSigs,'nearest');
   // SaO2Ve=interp1(TimeSig,SaO2Sig,Time_OtherSigs,'nearest');
   // HRVe=interp1(TimeSig,HRSig,Time_OtherSigs,'nearest');
   
   //
   // Metrics: AHI
   //
   
   res.TotalAHI= 60.0 * ne / res.TST ;
  
   
   //
   // No events? ... bail
   //
  
   if ( ne < 1 )
     {
       res.valid = false;
       logger << "  ** no events for hypoxic burden analysis, leaving\n";
       return;
     }
   
   
   //
   // Hypoxic burden
   //
   
   const bool use_modified_search_window = true;
   
   // 80 default;  alternative == 15
   const int MaxWin = param.has( "max-win" ) ? param.requires_int( "max-win" ) : 80;
   
   // for each time-point, get mean SpO2 (i.e. averaged over events/rows of SaO2Ve )
   Eigen::ArrayXd SpO2Mean = SaO2Ve.colwise().mean();
   
   // [HypoxicBurden,BaselineSat,HypoxicBurdenAll,BaselineSatAll,searchWin] =
   // 	FindBurden(SaO2Ve,SpO2Mean,Time_new,TST,MaxWin);
   
   logger << "  estimating hypoxic burden for all events";
   
   hb_find_burden_t burden = find_burden( SaO2Ve , SpO2Mean , t_sp , res.TST , MaxWin );

   //
   // Overall output
   //

   writer.value( "TST" , res.TST );
   writer.value( "AHI" , res.TotalAHI );
   
   if ( burden.valid )
     {
       writer.value( "HB" , burden.HB );
       writer.value( "BLSAT" , burden.BaselineSat );
     }
   
   
   //
   // Most frequency sleep stage for each event, i.e. is the event basically 'REM', 'NREM2', etc etc?
   //
   
   // MostFreqStage=mode(SleepStageVe,2);
   
   std::vector<sleep_stage_t> ss_mode( ne );
   for (int e=0; e<ne; e++)
     ss_mode[e] = modal_stage( SleepStageVe.row(e) );
   
   logger << ", also stratifying by sleep stage\n";

   
   //
   // Get NREM/REM specific means
   //
   
   // =nanmean(SaO2Ve(MostFreqStage>0 & MostFreqStage<5,:),1);
   // SpO2Mean_REM=nanmean(SaO2Ve(MostFreqStage==5,:),1);

   const int nwin = SpO2Mean.size();
   
   Eigen::ArrayXd SpO2Mean_NREM = Eigen::ArrayXd::Zero( nwin );
   Eigen::ArrayXd SpO2Mean_REM  = Eigen::ArrayXd::Zero( nwin );
   
   int denom_rem = 0 , denom_nrem = 0;
   
   for (int e=0; e<ne; e++)
     {
       if ( ss_mode[e] == REM )
	 {
	   ++denom_rem;
	   for (int j=0; j<nwin; j++) SpO2Mean_REM[j] += SaO2Ve(e,j);
	 }
       else
	 {
	   ++denom_nrem;
	   for (int j=0; j<nwin; j++) SpO2Mean_NREM[j] += SaO2Ve(e,j);
	 }
     }
   
   if ( denom_rem > 0 ) 
     SpO2Mean_REM /= denom_rem;
   
   if ( denom_nrem > 0 )
     SpO2Mean_NREM /= denom_nrem;
   

   //
   // Repeat burden analysis conditional on each stage
   //
   
   // [HypoxicBurden_N1,BaselineSat_N1]=FindBurden(SaO2Ve(MostFreqStage==1,:),SpO2Mean_NREM,Time_new,TST_N1,MaxWin);
   // [HypoxicBurden_N2,BaselineSat_N2]=FindBurden(SaO2Ve(MostFreqStage==2,:),SpO2Mean_NREM,Time_new,TST_N2,MaxWin);
   // [HypoxicBurden_N3,BaselineSat_N3]=FindBurden(SaO2Ve(MostFreqStage==3 | MostFreqStage==4,:),SpO2Mean_NREM,Time_new,TST_N3,MaxWin);
   // [HypoxicBurden_REM,BaselineSat_REM]=FindBurden(SaO2Ve(MostFreqStage==5,:),SpO2Mean_REM,Time_new,TST_REM,MaxWin);
   // [HypoxicBurden_NREM,BaselineSat_NREM]=FindBurden(SaO2Ve(MostFreqStage>0 & MostFreqStage<5,:),SpO2Mean_NREM,Time_new,TST_NREM,MaxWin);
   

   // N1
   std::vector<bool> incl = which_events( ss_mode , "N1" );
   hb_find_burden_t hb_n1  = find_burden( SaO2Ve , SpO2Mean_NREM , t_sp , res.TST_N1   , MaxWin , &incl );
   if ( hb_n1.valid )
     {
       writer.level( "N1" , globals::stage_strat );
       writer.value( "HB" ,  hb_n1.HB );
       writer.value( "BLSAT" , hb_n1.BaselineSat );
       writer.value( "NE" , hb_n1.ne );
       writer.value( "MINS" , res.TST_N1 );
       writer.value( "AHI" , 60.0 * hb_n1.ne / res.TST_N1 );
       
     }
   
   // N2
   incl = which_events( ss_mode , "N2" );
   hb_find_burden_t hb_n2  = find_burden( SaO2Ve , SpO2Mean_NREM , t_sp , res.TST_N2   , MaxWin , &incl );
   if ( hb_n2.valid )
     {
       writer.level( "N2" , globals::stage_strat );
       writer.value( "HB" ,  hb_n2.HB );
       writer.value( "BLSAT" , hb_n2.BaselineSat );
       writer.value( "NE" , hb_n2.ne );
       writer.value( "MINS" , res.TST_N2 );
       writer.value( "AHI" , 60.0 * hb_n2.ne / res.TST_N2 );
     }
   
   // N3
   incl = which_events( ss_mode , "N3" );
   hb_find_burden_t hb_n3  = find_burden( SaO2Ve , SpO2Mean_NREM , t_sp , res.TST_N3   , MaxWin , &incl );
   if ( hb_n3.valid )
     {
       writer.level( "N3" , globals::stage_strat );
       writer.value( "HB" ,  hb_n3.HB );
       writer.value( "BLSAT" , hb_n3.BaselineSat );
       writer.value( "NE" , hb_n3.ne );
       writer.value( "MINS" , res.TST_N3 );
       writer.value( "AHI" , 60.0 * hb_n3.ne / res.TST_N3 );
     }
   
   // NREM
   incl = which_events( ss_mode , "NREM" );
   hb_find_burden_t hb_nr  = find_burden( SaO2Ve , SpO2Mean_NREM , t_sp , res.TST_NREM , MaxWin , &incl );
   if ( hb_nr.valid )
     {      
       writer.level( "NREM" , globals::stage_strat );
       writer.value( "HB" ,  hb_nr.HB );
       writer.value( "BLSAT" , hb_n3.BaselineSat );
       writer.value( "NE" , hb_nr.ne );
       writer.value( "MINS" , res.TST_NREM );
       writer.value( "AHI" , 60.0 * hb_nr.ne / res.TST_NREM );
     }
   
   // REM
   if ( denom_rem )
     {
       incl = which_events( ss_mode , "REM" );
       hb_find_burden_t hb_rem = find_burden( SaO2Ve , SpO2Mean_REM  , t_sp , res.TST_REM  , MaxWin , &incl );
       if ( hb_rem.valid )
	 {
	   writer.level( "REM" , globals::stage_strat );
	   writer.value( "HB" ,  hb_rem.HB );
	   writer.value( "BLSAT" , hb_rem.BaselineSat );
	   writer.value( "NE" , hb_rem.ne );
	   writer.value( "MINS" , res.TST_REM );
	   writer.value( "AHI" , 60.0 * hb_rem.ne / res.TST_REM );
	 }
     }
   
   writer.unlevel( globals::stage_strat );
   

  
   //
   // Find desats, for hypoxic burden based on 3% or 4% desats
   //
   
   hb_find_desats_t desats = find_desats( eigen_ops::copy_array( SaO2 ) , fs , MagAv_thres );
    
   const int n_desats = desats.dsatStEnd.rows();

   logger << "  identified " << n_desats << " desats: ";
   
   // desats.dsstStEnd( ,0-2) = start, nadir, end
   const Eigen::ArrayXi & DesatSt = desats.dsatStEnd.col(0);
   const Eigen::ArrayXi & DesatNadir = desats.dsatStEnd.col(1);
   const Eigen::ArrayXi & DesatEnd = desats.dsatStEnd.col(2);
   
  // duration of desats, resats
   Eigen::ArrayXd DesatDur = Eigen::ArrayXd::Zero( n_desats );
   Eigen::ArrayXd ResatDur = Eigen::ArrayXd::Zero( n_desats );
   for (int i=0; i<n_desats; i++)
     {
       DesatDur[i] = seconds[ DesatNadir[i] ] - seconds[ DesatSt[i] ] ;
       ResatDur[i] = seconds[ DesatEnd[i] ] - seconds[ DesatNadir[i] ] ;
     }
   
   // magnitude of desats (start, nadir, stop)
   Eigen::ArrayXd DesatNadirMag = Eigen::ArrayXd::Zero( n_desats );
   Eigen::ArrayXd DesatStartMag = Eigen::ArrayXd::Zero( n_desats );
   Eigen::ArrayXd DesatEndMag   = Eigen::ArrayXd::Zero( n_desats );
   
   Eigen::ArrayXd DesatStartTime = Eigen::ArrayXd::Zero( n_desats );
   Eigen::ArrayXd DesatNadir_t   = Eigen::ArrayXd::Zero( n_desats );
   
   for (int i=0; i<n_desats; i++)
     {      
       DesatNadirMag[i]  = SaO2[ DesatNadir[i] ];
       DesatStartMag[i]  = SaO2[ DesatSt[i] ];
       DesatEndMag[i]    = SaO2[ DesatEnd[i] ];
       DesatStartTime[i] = seconds[ DesatSt[i] ];
       DesatNadir_t[i]   = seconds[ DesatNadir[i] ];
     }
   
   
   //
   // Construct lists of 3% and 4% events
   //
   
   // Events: EvIdx(:,1);
   // Desats linked to Respiratory events: EvIdx(:,2)
   
   // EvIdx=[0 0];
   // DesatMag=[0 0];
   
   std::map<int,int> EvIdx; // event --> linked desat
   std::map<int,double> DesatMag;  // event --> linked desat mag
   int last_desat = -1;
   
   // dummy end on evtSt
   // evtSt(length(evtEnd)+1)=TimeSig(end);
   // ne is number of events:
   
   for (int i=0 ; i<ne ; i++ )
     {

       // does this event a) has a desat nadir in the second half of the event
       //             &   b) that desar nadir is prior to the start of the next event
       //             &   c) within 75 seconds of the end of the event 
       //  if so, link that desat to this event
       
       double e_start = evtSt[i];
       double e_end = evtEnd[i];
       double e_mid = ( e_start + e_end ) / 2.0;
       double e_next = i < ne - 1 ? evtSt[i+1] : seconds[ seconds.size() - 1 ];

       //       std::cerr << " dets " << e_start << " " << e_end << " " << e_mid << " " << e_next << "\n";

       // idx = find(  DesatNadir_t >= (evtSt(ii) + evtEnd(ii) ) / 2
       // 		  &
       // 		  ( DesatNadir_t < evtSt(ii+1) & DesatNadir_t - evtEnd(ii) < 75 ) );
       
       int idx = -1;
       int ni = 0;

       while ( ni < n_desats )
	 {
	   
	   //std::cerr << " comp: desat-nadir " << DesatNadir_t[ ni ] << "\n";
	   
	   if ( DesatNadir_t[ ni ] >= e_mid
		&& DesatNadir_t[ ni ] < e_next
		&& DesatNadir_t[ ni ] - e_end < 75 )
	     {
	       idx = ni;
	       break;
	     }
	   ++ni;	  
	 }
       
       // if we find a desat, and this desat has not been linked to the previous 
       // if ( ~isempty(idx) )
       // 	{

       //std::cerr << "event " << i << " " << idx << "\n";
       
       if ( idx != -1 )
	 {
	  
	  // is the linked-desat already linked to the prior event?
	  //  if so, just ignore this event

	  // if ( idx(1)~=EvIdx(end,2) ) 
	  //   {

	  if ( idx != last_desat ) 
	    {
	      // EvIdx=[EvIdx;[ii idx(1)]];
	      // DesatMag=[DesatMag;[ii MagDown(idx(1))]];	      
	      EvIdx[ i ] = idx; 
	      DesatMag[ i ] = desats.MagDown[ idx ];	      
	      //std::cerr << " adding " << DesatMag[i] << "\n";
	      last_desat = idx;
	    }
	}
    }
  
  // not needed
  // EvIdx(1,:)=[];
  //   DesatMag(1,:)=[];
  
  // %%% List of desats with 3% or 4% desats
  // Idx3p=zeros(size(SaO2Ve,1),1);
  // Idx3p(DesatMag(DesatMag(:,2)>=3,1))=1;
  
  // Idx4p=zeros(size(SaO2Ve,1),1);
  // Idx4p(DesatMag(DesatMag(:,2)>=4,1))=1;


  // indicate which EVENTS have linked desats at 3 or 4%

  std::vector<bool> incl_3pct( ne , false );
  std::vector<bool> incl_4pct( ne , false );
  int cnt3 = 0 , cnt4 = 0;
  std::map<int,int>::const_iterator ee = EvIdx.begin();
  while ( ee != EvIdx.end() )
    {
      // Q. note: as EDF can have lots of 2.99999 etc, use round() here in the comparison...
      double desat = round( DesatMag[ ee->first ] );
      if ( desat >= 3 ) { incl_3pct[ ee->first ] = true; ++cnt3; }
      if ( desat >= 4 ) { incl_4pct[ ee->first ] = true; ++cnt4; }
      ++ee;
    }

  logger << cnt3 << " 3% and " << cnt4 << " 4% desats\n";

  //
  // Repeat HB analyses: total, NREM and REM for either 3% or 4%
  //

  writer.level( 3 , "DESAT" );

    
  // HB3=FindBurden(SaO2Ve(Idx3p==1,:),SpO2Mean,Time_new,TST,MaxWin);
  // HB3_REM=FindBurden(SaO2Ve(Idx3p==1 & MostFreqStage==5,:),SpO2Mean_REM,Time_new,TST_REM,MaxWin);
  // HB3_NREM=FindBurden(SaO2Ve(Idx3p==1 & MostFreqStage>0 & MostFreqStage<5,:),SpO2Mean_NREM,Time_new,TST_NREM,MaxWin);
  
  // nb. which_events() can take 3rd arg to specify baseline set (only T from that included, AND logic)


  logger << "  estimating HB only for events linked to 3%/4% desats\n";
    
  //
  // All, 3%
  //

  if ( enough( incl_3pct ) )
    {
      hb_find_burden_t hb = find_burden( SaO2Ve , SpO2Mean  , t_sp , res.TST  , MaxWin , &incl_3pct );
      if ( hb.valid )
	{
	  writer.value( "HB" ,  hb.HB );
	  writer.value( "NE" , hb.ne );	      
	}
    }
  
  //
  // REM 3%
  //
  
  incl = which_events( ss_mode , "REM" , &incl_3pct );
  if ( enough( incl ) )
    {
      hb_find_burden_t hb = find_burden( SaO2Ve , SpO2Mean_REM  , t_sp , res.TST_REM  , MaxWin , &incl );
      if ( hb.valid )
	{
	  writer.level( "REM" , globals::stage_strat );
	  writer.value( "HB" ,  hb.HB );
	  writer.value( "NE" , hb.ne );	      
	}
    }
  
  //
  // NREM 3%
  //
  
  incl = which_events( ss_mode , "NREM" , &incl_3pct );
  if ( enough( incl ) )
    {
      hb_find_burden_t hb = find_burden( SaO2Ve , SpO2Mean_NREM  , t_sp , res.TST_NREM  , MaxWin , &incl );
      if ( hb.valid )
	{
	  writer.level( "NREM" , globals::stage_strat );
	  writer.value( "HB" ,  hb.HB );
	  writer.value( "NE" , hb.ne );	      
	}
    }
  
  writer.unlevel( globals::stage_strat );


  //
  // 4% events
  //

  writer.level( 4 , "DESAT" );
  
  // HB4=FindBurden(SaO2Ve(Idx4p==1,:),SpO2Mean,Time_new,TST,MaxWin);
  // HB4_REM=FindBurden(SaO2Ve(Idx4p==1 & MostFreqStage==5,:),SpO2Mean_REM,Time_new,TST_REM,MaxWin);
  // HB4_NREM=FindBurden(SaO2Ve(Idx4p==1 & MostFreqStage>0 & MostFreqStage<5,:),SpO2Mean_NREM,Time_new,TST_NREM,MaxWin);


  //
  // All, 3%
  //

  if ( enough( incl_4pct ) )
    {
      hb_find_burden_t hb = find_burden( SaO2Ve , SpO2Mean  , t_sp , res.TST  , MaxWin , &incl_4pct );
      if ( hb.valid )
	{
	  writer.value( "HB" ,  hb.HB );
	  writer.value( "NE" , hb.ne );	      
	}
    }
  
  //
  // REM 4%
  //
  
  incl = which_events( ss_mode , "REM" , &incl_4pct );
  if ( enough( incl ) )
    {
      hb_find_burden_t hb = find_burden( SaO2Ve , SpO2Mean_REM  , t_sp , res.TST_REM  , MaxWin , &incl );
      if ( hb.valid )
	{
	  writer.level( "REM" , globals::stage_strat );
	  writer.value( "HB" ,  hb.HB );
	  writer.value( "NE" , hb.ne );	      
	}
    }
  
  //
  // NREM 4%
  //
  
  incl = which_events( ss_mode , "NREM" , &incl_4pct );
  if ( enough( incl ) )
    {
      hb_find_burden_t hb = find_burden( SaO2Ve , SpO2Mean_NREM  , t_sp , res.TST_NREM  , MaxWin , &incl );
      if ( hb.valid )
	{
	  writer.level( "NREM" , globals::stage_strat );
	  writer.value( "HB" ,  hb.HB );
	  writer.value( "NE" , hb.ne );	      
	}
    }
  
  writer.unlevel( globals::stage_strat );
  
  writer.unlevel( "DESAT" );


  //
  // All desats including non-event ones
  //

  
  // Total Area all desats
  double Dsat_area_REM = 0 , Dsat_area_REM_3 = 0;
  double Dsat_area_NREM = 0 , Dsat_area_NREM_3 = 0;

  // requires at least 3 desats 
  if ( n_desats > 2 )
    {
      for ( int ii=0; ii<n_desats; ii++ )
	{
	  // track 3pct desats separately
	  bool is_3pct = round( desats.MagDown[ ii ] ) >= 3.0;
	  
	  // temporary SpO2 signal during this desat
	  // Dsat_temp=SaO2Sig(dsatStEnd(ii,1):dsatStEnd(ii,3));
	  const int len = DesatEnd[ ii ] - DesatSt[ ii ] + 1;
	  Eigen::ArrayXd Dsat_temp = Eigen::ArrayXd::Zero( len );
	  int cnt = 0;
	  for (int t = DesatSt[ ii ]; t <= DesatEnd[ ii ]; t++)
	    Dsat_temp[cnt++] = SaO2[t];
	  
	  // sum, max, 
	  double sum = Dsat_temp.sum();
	  double max = Dsat_temp.maxCoeff();

	  // check # of missing values..
	  // nb. these are currently not tracked in SaO2
	  //	  if ( sum(isnan(Dsat_temp))/length(Dsat_temp)<=0.5 )

	  if ( true )
	    {

	      // SleepStageTemp = mode( SleepStageSig( dsatStEnd(ii,1) : dsatStEnd(ii,3) ) );
	      
	      Eigen::ArrayXi ss_temp = Eigen::ArrayXi::Zero( len );
	      cnt = 0;
	      for (int t = DesatSt[ ii ]; t <= DesatEnd[ ii ]; t++)
		ss_temp[cnt++] = ss[t];
	      sleep_stage_t mss = modal_stage( ss_temp );
	      
	      // area_temp=nansum(max(Dsat_temp)-Dsat_temp)*(1/Fs);
	      double area_temp = (max - Dsat_temp).sum() * (1.0 / fs );
	      
	      if ( mss == NREM1 || mss == NREM2 || mss == NREM3 || mss == NREM4 )
		{		  
		  Dsat_area_NREM += area_temp;
		  if ( is_3pct ) Dsat_area_NREM_3 += area_temp;
		}
	      else if ( mss == REM )
		{
		  Dsat_area_REM += area_temp;
		  if ( is_3pct ) Dsat_area_REM_3 += area_temp;
		}
	    }
	}  
    }
  
  double Dsat_area_Tot   = Dsat_area_NREM   + Dsat_area_REM;
  double Dsat_area_Tot_3 = Dsat_area_NREM_3 + Dsat_area_REM_3;

  // HB based on all desats 
  writer.value( "HB_TOT" , Dsat_area_Tot / res.TST );
  writer.level( "REM" , globals::stage_strat );
  writer.value( "HB_TOT" , Dsat_area_REM / res.TST_REM );  
  writer.level( "NREM" , globals::stage_strat );
  writer.value( "HB_TOT" , Dsat_area_NREM / res.TST_NREM );
  writer.unlevel( globals::stage_strat );

  // HB based only on 3% + desats
  writer.level( 3 , "DESAT" );
  writer.value( "HB_TOT" , Dsat_area_Tot_3 / res.TST );
  writer.level( "REM" , globals::stage_strat );
  writer.value( "HB_TOT" , Dsat_area_REM_3 / res.TST_REM );  
  writer.level( "NREM" , globals::stage_strat );
  writer.value( "HB_TOT" , Dsat_area_NREM_3 / res.TST_NREM );
  writer.unlevel( globals::stage_strat );
  writer.unlevel( "DESAT" );
  
  
  //
  // delta-HR function
  //

  logger << "  estimating event-related delta HR metrics\n";
  
  std::vector<double> dHR_start( ne );
  std::vector<double> dHR_end( ne );

  // Start=round(evtLevelSummary.evtSt-(evtLevelSummary.evtEnd-rangeX));
  // End=evtLevelSummary.evtEnd-(evtLevelSummary.evtEnd-rangeX);

  // Q. End definition?

  for (int e=0; e<ne; e++)
    {
      dHR_start[e] = evtSt[e] - ( evtEnd[e] - rangeX );
      dHR_end[e] = evtEnd[e] - ( evtEnd[e] - rangeX ); // huh, i.e. rangeX...
    }  

  // [dHR_meanbsline,dHR_minbsline,IndHRRep_mean,IndHRRep_min,N,meanbsline,minbsline]=SummarizeHR_AA(HRVe,Start,End,Time_new);
  
  delta_hr_t dHR = SummarizeHR_AA( HRVe , dHR_start, dHR_end , t_sp );
   
  writer.value( "DHR_MEAN_BL" , dHR.dHR_meanbsline );
  writer.value( "DHR_MIN_BL" , dHR.dHR_minbsline );
  writer.value( "DHR_N" , dHR.ne );
  
  // do we want these??
  // evtLevelSummary.dHR_meanBsln_PU=IndHRRep_mean;
  // evtLevelSummary.dHR_minBsln_PU=IndHRRep_min;
  // evtLevelSummary.meanBslnHR_PU=meanbsline;
  // evtLevelSummary.minBslnHR_PU=minbsline;  
  
}



hb_peakdet_t hb_t::peakdet( const Eigen::ArrayXd & v ,
			    double delta ,
			    bool flip = false )
{
  std::vector<double> t( v.size() );
  for (int i=0; i<t.size(); i++) t[i] = i;
  return peakdet( v , delta , t , flip );
}


hb_peakdet_t hb_t::peakdet( const Eigen::ArrayXd & v ,
			    double delta ,
			    const std::vector<double> & x ,
			    const bool flip = false )
{

  hb_peakdet_t r;

  // flip signal first? min->max
  int sgn = flip ? -1 : +1 ; 

  const int n = v.size();
  
  if ( n != x.size() )
    Helper::halt( "internal error in peakdet()" );
  
  double mn = 0, mnpos = 0;
  double mx = 0, mxpos = 0;
  
  bool lookformax = true;
  bool lookformin = false;

  for (int i=0; i<n; i++)
    {

      double th = sgn * v[i];      

      if ( th > mx ) { mx = th; mxpos = x[i]; }
      if ( th < mn ) { mn = th; mnpos = x[i]; }

      if ( lookformax )
	{
	  if ( th < mx - delta )
	    {
	      r.maxV.push_back( mx );
	      r.maxX.push_back( mxpos );
	      mn = th; mnpos = x[i];
	      lookformax = false;
	    }
	}
      else
	{
	  if ( th > mn + delta )
	    {
	      r.minV.push_back( mn );
	      r.minX.push_back( mnpos );
	      mx = th; mxpos = x[i];
	      lookformax = true;
	    }
	}
    }

  return r;
}

hb_find_burden_t hb_t::find_burden( const Eigen::MatrixXd & SpO2Mtx_orig ,
				    const Eigen::ArrayXd & SpO2Mean_orig ,
				    const std::vector<double> Time ,
				    const double TST , 
				    const int MaxWin ,
				    const std::vector<bool> * incl )
{

  hb_find_burden_t r;
  r.valid = false;
  
  // r.HB
  // r.BaselineSat
  // r.SpO2MtxDiff
  // r.BaselineSatAll ,
  // r.searchWin

  // CHECK: HypoxicBurdenAll == SpO2MtxDiff ?
  
  //
  // take desats as +ve
  //
  
  Eigen::ArrayXd SpO2Mean = 100 - SpO2Mean_orig;


  //
  // ... also, transpose to rows=time, cols=events
  //
  
  Eigen::MatrixXd SpO2Mtx = 100 - SpO2Mtx_orig.transpose().array();

  //
  // Keep means as supplied, but if we have an incl[], then splice out these rows
  //

  if ( incl != NULL )
    {
      if ( incl->size() != SpO2Mtx.cols() )
	Helper::halt( "problem in find_burden()" );

      int nn = 0;
      for (int i=0; i<incl->size(); i++)
	if ( (*incl)[i] ) ++nn;

      // no valid events
      if ( nn == 0 )
	{
	  logger << "  no valid events in find_burden()\n";
	  return r;
	}
      
      const int rows = SpO2Mtx.rows();
      
      Eigen::MatrixXd copy = SpO2Mtx;      
      SpO2Mtx.resize( rows , nn );
      
       nn = 0;
      for (int i=0; i<incl->size(); i++)
	{
	  if ( (*incl)[i] )
	    {
	      for (int r=0 ; r<rows; r++)
		SpO2Mtx(r,nn) = copy(r,i);
	      ++nn;
	    }
	}
    }
  

  // std::cerr << "Mean\n"
  // 	    << SpO2Mean << "\n\n";
  
  // Q. how NaNs get here..
  //    currently, no artifact removal on SaO2
  //    and out-of-range values replaced w/ the last/first value
  //    will need to add a boolean mask to track missing values
  //    ?versus, just remove the whole event?

  // NOTE: [ skipping this check, as currently no NaN possible ]
  // cannot have 5+ NaNs
  //  if [ # of NaNs in SpO2Mean ] >= 5 ... do not run
  
  //
  // Get peaks
  //

  //  [maxS,minS] = peakdet(SpO2Mean,0.1,Time);
  hb_peakdet_t peaks = peakdet( SpO2Mean, 0.1, Time );
  
  // Q. flag here?
  if ( peaks.minX.size() == 0 || peaks.maxX.size() == 0 )
    {
      logger << "  problem finding min/max peaks\n";
      return r;
    }
  
  // Search for minimum SpO2 in [eventEnd-10 to EvenEnd+MaxWin]  
  // max_resp=maxS(maxS(:,1)>=-10 & maxS(:,1)<=MaxWin,:);
  // max_resp=max_resp(max_resp(:,2)==max(max_resp(:,2)),:);

  double max_resp = -1;
  int max_resp_idx = -1;
  for (int i = 0 ; i < peaks.maxV.size() ; i++)
    {      
      if ( peaks.maxX[i] >= -10 && peaks.maxX[i] <= MaxWin )
	if ( peaks.maxV[i] >= max_resp )
	{
	  max_resp = peaks.maxV[i];
	  max_resp_idx = peaks.maxX[i];
	}
    }
  
  // if a minimum SpO2 in the ensemble-averaged signal is found
  // not found, bail
  if ( max_resp_idx == -1 )
    {
      logger << "  no minimum found in SpO2 average\n";
      return r; 
    }

  // find the pre-min maximum in the ensemble-averaged SpO2
  // min_pre=minS(minS(:,1)<max_resp(1,1),:); 
  // min_pre=min_pre(end,:); // if multiple, take last

  double min_pre = -1;
  int min_pre_idx = 0;
  for (int i = 0 ; i < peaks.minX.size() ; i++)
    if ( peaks.minX[i] < max_resp_idx )
      {
	// i.e. if multiple, this takes the last
	min_pre = peaks.minV[i]; 
	min_pre_idx = peaks.minX[i];
      }
  // find the post-min maximum in the ensemble-averaged SpO2
  // min_post=minS(minS(:,1)>max_resp(1,1),:);
  // min_post=min_post(1,:);

  double min_post = -1;
  int min_post_idx = 0;
  for (int i = 0 ; i < peaks.minX.size() ; i++)
    if ( peaks.minX[i] > max_resp_idx )
      {
	min_post = peaks.minV[i]; 
	min_post_idx = peaks.minX[i];
	break; // i.e. if multiple, take first
      }

  // Requires that both pre and post minimum are defined

  if ( min_pre < 0 || min_post < 0 )
    {
      logger << "  requires both pre/post minima are defined\n";
      return r;
    }
  
  // Spo2 Matrix during search window : only consider rows within min_pre/post
  //  SpO2MtxSrchWin=SpO2Mtx(Time>=min_pre(1,1) & Time <= min_post(1,1),:); 

  int rows = 0;
  
  for (int t=0; t<Time.size(); t++)
    if ( Time[t] >= min_pre_idx && Time[t] <= min_post_idx ) ++rows;
    
  // extract matrix w/ only the new data
  
  const int cols = SpO2Mtx.cols();

  Eigen::MatrixXd SpO2MtxSrchWin = Eigen::MatrixXd::Zero( rows , cols );
  
  rows = 0;
  for (int t=0; t<Time.size(); t++)
    if ( Time[t] >= min_pre_idx && Time[t] <= min_post_idx )
      {
	for (int j=0;j<cols;j++)
	  SpO2MtxSrchWin( rows , j ) = SpO2Mtx( t , j );
	++rows;
      }

  // SpO2MtxSrchWin [ search-window/samples x events ] 
  
  // Maximum SpO2 during search window is defined as baseline SpO2
  // for each event

  //  SpO2BaslineMtx=repmat(min(SpO2MtxSrchWin),size(SpO2MtxSrchWin,1),1); 

  Eigen::ArrayXd SpO2BaslineMtx = SpO2MtxSrchWin.colwise().minCoeff(); 
  
  // return matrix window x events  (but all 
  
  // Baseline Sat for each event
  //  r.BaselineSatAll=100-mean(SpO2BaslineMtx);
  r.BaselineSatAll = 100 - SpO2BaslineMtx;

  // Mean Baseline Sat all events
  //  r.BaselineSat=mean(100-mean(SpO2BaslineMtx)); 
  r.BaselineSat = ( 100 - SpO2BaslineMtx ).mean();

  // Remove baseline from SpO2 curve during search window
  //  r.SpO2MtxDiff=SpO2MtxSrchWin-SpO2BaslineMtx; 

  r.SpO2MtxDiff = SpO2MtxSrchWin;
  
  for (int i=0; i < r.SpO2MtxDiff.rows(); i++) 
    for (int j=0; j < r.SpO2MtxDiff.cols() ; j++) 
      r.SpO2MtxDiff(i,j) = r.SpO2MtxDiff(i,j) - SpO2BaslineMtx[j];
  
  // Hypoxic burden
  // Q. where do NAN come in?
  //  r.HB=nansum(nansum(SpO2MtxDiff))/TST; 
  r.HB = r.SpO2MtxDiff.array().sum() / TST;

  // track number of events included here
  r.ne = r.BaselineSatAll.size();
  
  // Area under desat for each event
  r.SpO2MtxDiff = r.SpO2MtxDiff.colwise().sum(); 
  
  // return search window
  r.searchWin_lwr = min_pre_idx ;
  r.searchWin_upr = min_post_idx;

  // all good
  r.valid = true;

  return r;

}


sleep_stage_t hb_t::modal_stage( const Eigen::ArrayXi & d )
{
  // 0=W/other; 1,2,3=NR; 5=R
  // if ties, call W > R > NR
  std::map<int,int> counts;
  const int n = d.size();
  for ( int i=0; i<n; i++)
    ++counts[d[i]];
  
  int max = counts[0];
  for (int i=1;i<=5;i++)
    if ( counts[i] > max ) max = counts[i];
  if ( counts[0] == max ) return WAKE;
  if ( counts[5] == max ) return REM;
  if ( counts[1] == max ) return NREM1;
  if ( counts[2] == max ) return NREM2;
  if ( counts[3] == max ) return NREM3;
  return WAKE;
}


std::vector<bool> hb_t::which_events( const std::vector<sleep_stage_t> & ss , const std::string & s ,
				      const std::vector<bool> * orig )
{
  
  const int n = ss.size();
  
  std::vector<bool> incl( n , false );
  
  if ( s == "N1" )
    for (int i=0; i<n; i++) 
      incl[i] = ss[i] == NREM1;
  
  if ( s == "N2" )
    for (int i=0; i<n; i++) 
      incl[i] = ss[i] == NREM2;
  
  if ( s == "N3" )
    for (int i=0; i<n; i++) 
      incl[i] = ss[i] == NREM3;

  if ( s == "REM" )
    for (int i=0; i<n; i++) 
      incl[i] = ss[i] == REM;
  
  if ( s == "NREM" )
    for (int i=0; i<n; i++) 
      incl[i] = ss[i] == NREM1 || ss[i] == NREM2 || ss[i] == NREM3;

  // additional mask? i.e. if not included here, set to F   
  if ( orig != NULL && orig->size() == n )
    for (int i=0; i<n; i++)
      if ( ! (*orig)[i] ) incl[i] = false;
  
  return incl;

}




// Oxygen desaturation finder [ AA & SS ]

hb_find_desats_t hb_t::find_desats( const Eigen::ArrayXd & s , 
				    int fs ,
				    double MagAv_thres )
{

  const double dt = 1.0 / fs;

  // we assume that we've already filled in missing varaibles
  // prior to calling, i.e. done at start of hb_t()

  // SaO2=fillmissing(SaO2_in,'nearest'); 

  // [SaO2Mins,SaO2Maxs]=peakdet(-SaO2,0.5);

  const bool flip_signal = true; 

  hb_peakdet_t p = peakdet( s , 0.5 , flip_signal );

  // constrain indices 
  // SaO2MaxIdx=SaO2Maxs(:,1);
  // SaO2MinIdx=SaO2Mins(:,1);

  // use deque for popping off front and back
  // **** nb. swaping of order minX --> MaxIdx, and vice versa **** 
  // to mirror  [SaO2Mins,SaO2Maxs] = peakdet() above, which usually
  //  is expected to return max / min 

  std::deque<int> SaO2MaxIdx( p.minX.begin() , p.minX.end() ); // nb. swapping min/max
  std::deque<int> SaO2MinIdx( p.maxX.begin() , p.maxX.end() ); // nb. swapping min/max

  //  std::cerr << SaO2MinIdx.size() << " " << SaO2MaxIdx.size() << "\n";

  // while SaO2MaxIdx(1)<0
  //        SaO2MaxIdx(1)=[];
  //  end

  while ( SaO2MaxIdx.front() < 0 )
    SaO2MaxIdx.pop_front();
    
  // while SaO2MinIdx(end)>=SaO2MaxIdx(end)
  //       SaO2MinIdx(end)=[];
  // end
  
  while ( SaO2MinIdx.back() >= SaO2MaxIdx.back() )
    SaO2MinIdx.pop_back();
  
  // while SaO2MinIdx(1)<SaO2MaxIdx(1)
  //    SaO2MinIdx(1)=[];
  // end

  while ( SaO2MinIdx.front() <= SaO2MaxIdx.front() )
    SaO2MinIdx.pop_front();
 
  // the above ensures that starts and ends with MAX
  //   [MAX] MIN [MAX] MIN [MAX]
  // i.e. we should always have 1 less MIN than MAX
  //  MagDown is 
  
  // get diffs
  // MagDown = SaO2(SaO2MaxIdx(1:end-1)) - SaO2(SaO2MinIdx(1:end));
  // MagUp = SaO2(SaO2MaxIdx(2:end))     - SaO2(SaO2MinIdx(1:end));
  // MagAv = (MagDown+MagUp)/2;

  int n_min = SaO2MinIdx.size();
  int n_max = SaO2MaxIdx.size();

  // n_min should be 1 less than n_max

  //  std::cerr << " n_min, n_max = " << n_min << " " << n_max << "\n";
  if ( n_max - n_min != 1 ) logger << "  *** hmm, warning ... need to check find_desats() implementation ***\n";
  
  // these should be the same; guess that must be
  // guaranteed by peakdet()
  
  std::vector<double> MagDown, MagUp, MagAv;  
  for (int i=0; i<n_min; i++ )    
    {
      double md = s[ SaO2MaxIdx[i] ]   - s[ SaO2MinIdx[i] ] ;
      double mu = s[ SaO2MaxIdx[i+1] ] - s[ SaO2MinIdx[i] ] ;
      MagDown.push_back( md );
      MagUp.push_back( mu );
      MagAv.push_back( md+mu / 2.0 );      
    }

  
  //
  // Iterate
  //

  while ( 1 )
    {
      // [minVTinspexp,minVTinspexp_i] = min(MagDown);
      // [minVTinverted,minVTinverted_i] = min(MagUp);
      // [minVT,VTpattern] = min([minVTinspexp,minVTinverted]);
      
      int minVTinspexp_i = -1;
      double minVTinspexp = min( MagDown , &minVTinspexp_i );
      
      int minVTinverted_i = -1;
      double minVTinverted = min( MagUp , &minVTinverted_i );
      
      double minVT = minVTinspexp <= minVTinverted ? minVTinspexp : minVTinverted;
      int VTpattern = minVTinspexp <= minVTinverted ? 1 : 2 ;

      // all done?
      if ( minVT > MagAv_thres ) break;

      if ( VTpattern == 1 ) // down
	{
	  // i=minVTinspexp_i;
	  // SaO2MaxIdx(i)=[];
	  // SaO2MinIdx(i)=[];
	  
	  SaO2MaxIdx.erase( SaO2MaxIdx.begin() + minVTinspexp_i );
	  SaO2MinIdx.erase( SaO2MinIdx.begin() + minVTinspexp_i );
	}
      else if ( VTpattern == 2 ) // %up
	{
 	  // i=minVTinverted_i;
	  // SaO2MaxIdx(i+1)=[];
	  // SaO2MinIdx(i)=[];

	  SaO2MaxIdx.erase( SaO2MaxIdx.begin()+minVTinverted_i+1 );
	  SaO2MinIdx.erase( SaO2MinIdx.begin()+minVTinverted_i  );

	}

      // recalculate
      n_min = SaO2MinIdx.size();

      MagDown.clear();
      MagUp.clear();
      MagAv.clear();
      for (int i=0; i<n_min; i++ )    
	{
	  double md = s[ SaO2MaxIdx[i] ] - s[ SaO2MinIdx[i] ] ;
	  double mu = s[ SaO2MaxIdx[i+1] ] - s[ SaO2MinIdx[i] ] ;
	  MagDown.push_back( md );
	  MagUp.push_back( mu );
	  MagAv.push_back( md+mu / 2.0 );      
	}

      // if isempty(SaO2MaxIdx) || isempty(SaO2MinIdx)
      //      break
      //  end
      
      if ( SaO2MaxIdx.size() == 0 || SaO2MinIdx.size() == 0 ) break;
      
    }

  
  // MagDown = SaO2(SaO2MaxIdx(1:end-1))-SaO2(SaO2MinIdx(1:end));
  // MagUp = SaO2(SaO2MaxIdx(2:end))-SaO2(SaO2MinIdx(1:end));

  n_min = SaO2MinIdx.size();
  MagDown.clear();
  MagUp.clear();
  MagAv.clear();
  
  for (int i=0; i<n_min; i++ )
    {
      double md = s[ SaO2MaxIdx[i] ] - s[ SaO2MinIdx[i] ] ;
      double mu = s[ SaO2MaxIdx[i+1] ] - s[ SaO2MinIdx[i] ] ;
      MagDown.push_back( md );
      MagUp.push_back( mu );
      MagAv.push_back( md+mu / 2.0 );
    }

  // constants
  
  const int searchmaxleftrange = 120;
  const int searchmaxrightrange = 120;
  
  const int searchmaxleftrangei = round( searchmaxleftrange/dt );
  const int searchmaxrightrangei = round( searchmaxrightrange/dt );
  
  // SpO2prei = SaO2MaxIdx(1:end-1);
  // SpO2posti = SaO2MaxIdx(2:end);
  
  std::deque<int> SpO2prei = SaO2MaxIdx;
  SpO2prei.pop_back();
  
  std::deque<int> SpO2posti = SaO2MaxIdx;
  SpO2posti.pop_front();

  const int n = SaO2MinIdx.size();
  
   for ( int i=0; i<n; i++ )
     {
  
       int ileft = SaO2MinIdx[i]-SaO2MaxIdx[i];
       if ( ileft > searchmaxleftrangei )
	 ileft = searchmaxleftrangei;

       int iright = SaO2MaxIdx[i+1]-SaO2MinIdx[i];
       if ( iright > searchmaxrightrangei )
	 iright = searchmaxrightrangei;

       // left: find index of last max (so use reverse array)
       // [valL,~] = max( SaO2( SaO2MinIdx(i) - ileft : SaO2MinIdx(i) ) );
       // Il = find( SaO2(SaO2MinIdx(i) - ileft : SaO2MinIdx(i) ) ==valL , 1 , 'last' );
       // SpO2prei(i) = SaO2MinIdx(i)-ileft+Il-1;
       
       int from_idx = SaO2MinIdx[i] - ileft;
       int len      = ileft + 1;
       int to_idx   = SaO2MinIdx[i];

      // Q. hmm, need to ensure fits in range??
       if ( from_idx < 0 )
	 {
	   len = to_idx + 1;
	   from_idx = 0;
	 }
       
       int Il;
       double valL = s.segment( from_idx , len ).reverse().maxCoeff(&Il);
       Il = ( to_idx - from_idx + 1 ) - Il - 1 ; // reverse index back
       SpO2prei[i] = SaO2MinIdx[i] - ileft+Il; // -1; 0-based arrats, so skip -1
  
       // right: find index of first max
       // [valR,~]=max(SaO2(SaO2MinIdx(i):SaO2MinIdx(i)+iright));
       // Ir = find(SaO2(SaO2MinIdx(i):SaO2MinIdx(i)+iright)==valR,1,'first');
       //  SpO2posti(i) = SaO2MinIdx(i)+Ir-1;

       from_idx = SaO2MinIdx[i];
       len      = iright + 1 ;
       to_idx   = SaO2MinIdx[i] + len - 1 ;

       // Q. hmm, need to ensure fits in range??       
       if ( to_idx >= s.size() ) len -= to_idx - s.size() + 1 ;
 	 
       int Ir;
       double valR = s.segment( from_idx , len ).maxCoeff(&Ir);
       SpO2posti[i] = SaO2MinIdx[i] + Ir; // -1; 0-based arrays so skip
       
     }

   //
   // Return results
   //

   //     MagDown = SaO2(SpO2prei)-SaO2(SaO2MinIdx);
   //     MagUp = SaO2(SpO2posti)-SaO2(SaO2MinIdx);
   //     Is = [SpO2prei SaO2MinIdx SpO2posti];
   
   hb_find_desats_t r;
   
   r.MagDown.resize(n);
   r.MagUp.resize(n);
   r.dsatStEnd.resize( n , 3 );

   
   for (int i=0; i<n; i++)
     {
       r.MagDown[i] = s[SpO2prei[i]] - s[SaO2MinIdx[i]] ;
       r.MagUp[i] = s[SpO2posti[i]] - s[SaO2MinIdx[i]] ;
       r.dsatStEnd(i,0) = SpO2prei[i];
       r.dsatStEnd(i,1) = SaO2MinIdx[i];
       r.dsatStEnd(i,2) = SpO2posti[i];       
     }

     
  return r;
}


bool hb_t::enough( const std::vector<bool> & x , int th )
{
  int c = 0;
  const int n = x.size();
  for (int i=0; i<n; i++)
    if ( x[i] ) ++c;
  return c >= th;
}

double hb_t::min( const std::vector<double> & s , int * idx )
{
  const int n = s.size();
  double m = s[n-1];
  // go backwards, so we get the first instance of the min if ties
  for (int i=n-1; i>=0; i--)
    if ( s[i] <= m ) { m = s[i]; *idx = i; }
  return m;
}




 delta_hr_t hb_t::SummarizeHR_AA( const Eigen::MatrixXd & HRVe ,
 				 const std::vector<double> & Start ,
 				 const std::vector<double> & End ,
 				 const std::vector<double> & Time )
 {
   
   // check Time_new is SP or SEC?
   
   // get 1% and 99.99% percentiles
   //percntiles = prctile(reshape(HR,[],1),[1 99.99]); %1th and 99.99th percentile
   const int ne = HRVe.rows();
   const int nt = HRVe.cols(); 
   
   std::vector<double> xx;
   for (int i=0; i<ne; i++)
     for (int j=0; j<nt; j++)
       xx.push_back( HRVe(i,j) );

   double p01 = MiscMath::percentile( xx , 0.0100 ); 
   double p99 = MiscMath::percentile( xx , 0.9999 ); 
   xx.clear();

   logger << "  HR percentiles (1st, 99.99th) = " << p01 << ", " << p99 << "\n";

   // remove outlier values
   // or FLAG outliers (i.e. no NaN)
   // HR (HR<percntiles(1) | HR > percntiles(2) ) = NaN;
   
   std::vector<std::vector<bool> > nan( ne );
   for (int i=0; i<ne; i++)
     {
       nan[i].resize( nt , false );
       for (int j=0; j<nt; j++)
	 if ( HRVe(i,j) < p01 || HRVe(i,j) > p99 )
	   nan[i][j] = true;
     }

   // baseline1=nanmean(HR,2); % average of all 200s
   // baseline2=nanmin(HR(:,round(nanmean(Start))+1:round(nanmean(End))+1),[],2); % min value during events

   Eigen::ArrayXd baseline1 = Eigen::ArrayXd::Zero( ne ); // average over time (for each event)
   Eigen::ArrayXd baseline2 = Eigen::ArrayXd::Zero( ne ); // min for each event (During avg event time)
   std::vector<bool>   nan_event(ne,false); // track NaNs
   
   // Q/TODO: confirm +1 needed
   int mean_start = round( MiscMath::mean( Start ) ) + 1 ;
   int mean_end   = round( MiscMath::mean( End ) ) + 1 ;
   
   for (int i=0; i<ne; i++)
     {
       double m = 999999;
       double s = 0;
       int cnt = 0;
       for (int j=0; j<nt; j++)
	 {
	   if ( ! nan[i][j] )
	     {
	       if ( HRVe(i,j) < m )
		 m = HRVe(i,j);
	       ++cnt;
	       s += HRVe(i,j);
	     }

	   if ( cnt == 0 )
	     nan_event[i] = true;
	   else
	     {
	       baseline1[i] = s / (double)(cnt);
	       baseline2[i] = m;	   
	     }
	 }
     }
   
   // Default search window is between [event end - eventDur/2 : event end + 50 seconds]

   // defaultSrchWin= round(nanmean(End))-round(nanmean(End)-nanmean(Start))/2+1 : 151;
   
   int def_sw_start = (double)round(MiscMath::mean(End)) - (double)round( MiscMath::mean(End)-MiscMath::mean(Start))/2.0 + 1 ; 
   int def_sw_end   = 151;
   
   //  Modify the end of search window using average HR
   //  Ensembled avg HR
   //   AvgHR=nanmean(HR,1);
   
   Eigen::ArrayXd AvgHR = Eigen::ArrayXd::Zero( nt );
   std::vector<bool>   nan_AvgHR( nt , false );
   
   for (int j=0; j<nt; j++)
     {
       int cnt = 0;
       double s = 0;
       for (int i=0; i<ne; i++)
	 if ( ! nan[i][j] )
	   {
	     ++cnt;
	     s += HRVe(i,j);
	   }

       if ( cnt == 0 )
	 nan_AvgHR[j] = true;
       else
	 AvgHR[j] = s/(double)cnt;
     }

  
   //
   // Requires atleast 5 non-missing events
   //

   int n_nonmissing = 0;
   for (int i=0; i<ne; i++) if ( ! nan_event[i] ) ++n_nonmissing;
   
   //if sum(~isnan(baseline1))>=5


   if ( n_nonmissing >= 5 )
     {
       
       // [maxS,minS] = peakdet(AvgHR,0.25,Time);

       hb_peakdet_t p = peakdet( AvgHR , 0.25 , Time );
       
       // %%%Max response occurs somewhere between event end -10 and event end
       // 	%%%+50 seconds

       const int MaxWin=50;

       //       if ~isempty(maxS) & ~isempty(minS)
       
       if ( p.minX.size() != 0 && p.maxX.size() != 0 ) 
	 {					

	   //  max_resp=maxS(maxS(:,1)>=-10 & maxS(:,1)<=MaxWin,:);
           //  max_resp=max_resp(max_resp(:,2)==max(max_resp(:,2)),:);

	   double max_resp = -999;
	   double max_resp_x = -999;
	   for (int i=0;i<p.maxX.size();i++)
	     {
	       if ( p.maxX[i] >= -10 && p.maxX[i] <= MaxWin )
		 if ( p.maxV[i] > max_resp )
		   {
		     max_resp = p.maxV[i];
		     max_resp_x = p.maxX[i];
		   }
	     }

	   // did we find a valid max?

	   //if ~isempty(max_resp)

	   if ( max_resp > 0 )
	     {

	       // min_pre=minS(minS(:,1)<max_resp(1,1),:);
	       // if ~isempty(min_pre)
               //      min_pre=min_pre(end,:);
               // end

	       // get last minima before max
	       // i.e. overvwrites, will pick last
	       double min_pre = -99999;
	       for (int i=0; i<p.minX.size(); i++)
		 if ( p.minX[i] < max_resp_x )
		   min_pre = p.minX[i];

	       
	       // min_post=minS(minS(:,1)>max_resp(1,1),:);
	       // if ~isempty(min_post)
               //   min_post=min_post(1,:);
               // end

	       // get first min post max; i.e. track backwards

	       double min_post = -99999;
               for (int i=p.minX.size()-1 ; i>=0 ; i--)
                 if ( p.minX[i] > max_resp_x )
                   min_post = p.minX[i];

	       // 	if ~isempty(min_pre) && ~isempty(min_post)
               //      defaultSrchWin=find(Time>=min_pre(1,1) & Time <= min_post(1,1));
	       // 	end

	       if ( min_pre > -88888 && min_post > -88888 )
		 {
		   int s1 = round(min_pre);
		   int s2 = round(min_post);
		   
		   def_sw_start = -1;
		   def_sw_end = nt-1;

		   for (int t=0; t<nt; t++)
		     {
		       if ( def_sw_start == -1 && Time[t] >= s1 ) def_sw_start = t;
		       if ( Time[t] <= s2 ) def_sw_end = t;
		     }		   
		 }	       
	     }
	 }
     }
   

   // defaultSrchWin=round(defaultSrchWin);
   // (done above)
   
   // get max value during search window
   //  maxHR=nanmax(HR(:,defaultSrchWin),[],2); 

   Eigen::ArrayXd maxHR = Eigen::ArrayXd::Zero( ne );

   for (int i=0; i<ne; i++)
     {
       double mx = -9;
       for (int j=def_sw_start; j<=def_sw_end; j++)
	 if ( ! nan[i][j] )
	   if ( HRVe(i,j) > mx ) mx = HRVe(i,j);
       maxHR[i] = mx;       
     }

   delta_hr_t r;

   r.dHR_meanbsline = 0;
   r.dHR_minbsline = 0;
   
   // dHR_meanbsline=nanmean(maxHR-baseline1);
   //  dHR_minbsline=nanmean(maxHR-baseline2);
   //  IndHRRep_mean=maxHR-baseline1;
   //  IndHRRep_min=maxHR-baseline2;
   //  N=sum(~isnan(maxHR));
   //  meanbsline=baseline1;
   //  minbsline=baseline2;

   r.IndHRRep_mean = Eigen::ArrayXd::Zero( ne );
   r.IndHRRep_min  = Eigen::ArrayXd::Zero( ne );
   
   int cnt = 0;
   for (int i=0; i<ne; i++)
     {
       if ( ! nan_event[i] )
	 {
	   ++cnt;
	   r.dHR_meanbsline += maxHR[i] - baseline1[i] ;
	   r.dHR_minbsline  += maxHR[i] - baseline2[i] ;
	 }
       r.IndHRRep_mean[i] = maxHR[i] - baseline1[i];
       r.IndHRRep_min[i]  = maxHR[i] - baseline2[i];
     }

   if ( cnt > 0 )
     {
       r.dHR_meanbsline /= (double)cnt;
       r.dHR_minbsline /= (double)cnt;
     }

   r.ne = cnt;
   r.meanbsline = baseline1;
   r.minbsline  = baseline2;
   r.nan        = nan_event;

  return r;

 }


