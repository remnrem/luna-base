
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

#include "dsp/arousals.h"
#include "defs/defs.h"
#include "edf/edf.h"
#include "edf/slice.h"
#include "stats/eigen_ops.h"
#include "miscmath/miscmath.h"

#include "helper/logger.h"
#include "db/db.h"

extern logger_t logger;
extern writer_t writer;


arousals_t::arousals_t( edf_t & edf , param_t & param )
{

  // use EPOCH mechanism to add new Hjorth signals
  // requires EPOCH inc is 1 Hz

  if ( ! edf.timeline.epoched() )
    Helper::halt( "currently, must have 1-second epoch increment (use e.g. EPOCH inc=1 len=4)" );
  
  if ( edf.timeline.epoch_increment_tp() != globals::tp_1sec )
    Helper::halt( "currently, must have 1-second epoch increment (use e.g. EPOCH inc=1 len=4)" );

  // also for symmetry, require the EPOCH duration is odd # of whole seconds;
  if ( edf.timeline.epoch_len_tp_uint64_t() % globals::tp_1sec != 0LLU )
    Helper::halt( "currently, epoch duration must be an exact/odd number of seconds" );
  
  int esec = edf.timeline.epoch_len_tp_uint64_t() / globals::tp_1sec;

  if ( esec % 2 == 0 ) Helper::halt( "currently, epoch duration must be an exact/odd number of seconds" );
  
  if ( edf.timeline.generic_epochs() )
    Helper::halt( "can not have generic epochs" );
  
  if ( edf.timeline.epoch_any_offset() )
    Helper::halt( "cannot use with an EPOCH offset (e.g. from align)" );
  
  // to save the signal, the EDF record size must be consistent w/ the
  // new SR of the signal i.e. to have same number of samples per EDF
  // record, cannot have SR = 0.5 Hz with 1-second records, etc
  // perhaps for now, enfore that EDF record size is a multiple of 1
  // second, and always output at 1 Hz

  uint64_t d = edf.header.record_duration_tp % globals::tp_1sec;
  if ( d != 0 )
    Helper::halt( "must have EDF records that are a multiple of 1 second (use RECORD-SIZE)" );
 

  //
  // Options
  //
  
  const std::string new_sig_prefix = param.has( "prefix" ) ? param.value( "prefix" ) : "" ; 

  const bool write_combined_metrics = param.has( "combined" ) ? param.yesno( "combined" ) : true ;
  const bool write_channel_specific_metrics = param.has( "per-channel" ) ? param.yesno( "per-channel" ) : false;
    
  const double lag_sec = param.has( "w" ) ? param.requires_dbl( "w" ) : 10.0;  
  const double threshold = param.has( "th" ) ? param.requires_dbl( "th" ) : 3 ;  
  const double influence = param.has( "influence" ) ? param.requires_dbl( "influence" ) : 0 ;  
  const double mindur = param.has( "sec" ) ? param.requires_dbl( "sec" ) : 3.0; 
  const double maxdur = param.has( "max-sec" ) ? param.requires_dbl( "max-sec" ) : 15.0;
  
  double winsor_th = param.has( "no-winsor" ) ? -9 : 0.005;
  if ( param.has( "winsor" ) ) winsor_th = param.requires_dbl( "winsor" ) ;
  if ( winsor_th > 0.2 ) Helper::halt( "winsor should be less than 0.2" );
				        
  if ( lag_sec < 1 ) Helper::halt( "w must be > 1 second" );
  if ( threshold < 0 ) Helper::halt( "th must be positive" );
  if ( influence < 0 || influence > 1 ) Helper::halt( "influence must be 0 .. 1" );
  if ( mindur < 1 || mindur > 15 ) Helper::halt( "sec (min duration) must be between 1 and 15 seconds" );
  
  const std::string aname = param.has( "annot" ) ? param.value( "annot" ) : "arousal_luna" ;    
  
  //
  // Get signals
  //
  
  const bool has_eeg = param.has( "eeg" );
  const bool has_emg = param.has( "emg" );
  
  const std::string eeg_signal_label = has_eeg ? param.value( "eeg" ) : "" ;
  const std::string emg_signal_label = has_emg ? param.value( "emg" ) : "" ;
    
  const bool NO_ANNOTS = true;   
  signal_list_t eeg_signals = edf.header.signal_list( eeg_signal_label , NO_ANNOTS );
  signal_list_t emg_signals = edf.header.signal_list( emg_signal_label , NO_ANNOTS );  
  
  const int ns_eeg = eeg_signals.size();
  const int ns_emg = emg_signals.size();
  const int ns     = ns_eeg + ns_emg;
  
  logger << "  running AROUSALS for " << ns_eeg << " EEG and " << ns_emg << " EMG signals\n";

  if ( ns == 0 ) return;
  
  
  //
  // Annots
  //
  
  annot_t * a = edf.timeline.annotations.add( aname );
  a->description = "Arousals";
  

  //
  // Obtain sampling freqs (Hz) - these must be similar
  //
  
  std::vector<double> Fs_eeg = edf.header.sampling_freq( eeg_signals );
  std::vector<double> Fs_emg = edf.header.sampling_freq( emg_signals );
  
  const int sr = Fs_eeg.size() ? Fs_eeg[0] : Fs_emg[0] ;

  for (int i=0; i<ns_eeg; i++)
    {
      if ( fabs( Fs_eeg[i] - sr ) > 1e-4 )
	Helper::halt( "all signals must have similar sample rates" );      
    }
  
  for (int i=0; i<ns_emg; i++)
    {
      if ( fabs( Fs_emg[i] - sr ) > 1e-4 )
	Helper::halt( "all signals must have similar sample rates" );      
    }
  
  
  //
  // Figure out alignment/sampling of new signals - and padding needs to happen
  // for each segment;   will be a simple N = (D-1)/2 seconds/samples both a start
  // and end, where D is integer, odd number of seconds;  
  //  
  
  const int pad = (esec-1)/2;
  
  std::cout << " padding " << pad << "s at start/end of each contiguous epoch\n"; 
  
  // expected # of samples given 1 Hz rate
  const int expected = edf.header.nr * edf.header.record_duration ;
  
  // but we have this many epochs, because of gaps
  const int obs = edf.timeline.num_epochs();
  
  // this is the total adjustment needed;  we'll splice in padding
  // along the way; but we need to size for expected overall;
  const int diff = expected - obs;
  
  
  //
  // track Hjorths
  //      
  
  const int ne = edf.timeline.first_epoch();
  
  Eigen::MatrixXd H1 = Eigen::MatrixXd::Zero( expected , ns );
  Eigen::MatrixXd H2 = Eigen::MatrixXd::Zero( expected , ns );
  Eigen::MatrixXd H3 = Eigen::MatrixXd::Zero( expected , ns );
  
  // track 'bad' (i.e. flat) epochs channel-wise
  std::map<int,std::set<int> > bad;
  
  
  //
  // iterate over epochs
  //
  
  // to track gaps -- will ultimately span
    
  std::vector<uint64_t> estarts;
  
  // prior sample
  uint64_t prior = 0LLU;
  
  int row = 0 ;
  
  while ( 1 ) 
    {
      
      int epoch = edf.timeline.next_epoch();      

      // final epoch?
      if ( epoch == -1 )
	{
	  // final 0-padding
	  for (int p=0;p<pad;p++)
	    {
	      for (int s=0;s<ns; s++) bad[ row ].insert( s );		
	      ++row;
	      uint64_t e = prior + globals::tp_1sec * ( p+1 ) ;
	      //std::cout << " pad (end) " << p << " " << e << "\n";
	      estarts.push_back( e ) ;
	    }
	  break;
	}
      
      // get current 
      interval_t interval = edf.timeline.epoch( epoch );

      // crossed gap?
      const bool crossed_gap = row != 0 && interval.start != prior + globals::tp_1sec ; 
      
      // first epoch?
      if ( row == 0 )
	{
	  
	  for (int p=0; p<pad; p++)
	    {	      
	      for (int s=0;s<ns; s++) bad[ row ].insert( s );
	      uint64_t e = interval.start - globals::tp_1sec * ( pad - p ) ;
	      //std::cout << " pad (1st) " << p << " " << e << "\n";
	      estarts.push_back( e ) ;
	      ++row;	      	      
	      // pad:::  0  1  2
	      //   --->  3  2  1  --->  3 - i
	    }

	}
      else
	{
	  
	  // gap between epochs? (i.e. cross-segment epoch vs prior ) 	  
	  if ( crossed_gap )
	    {
	      
	      // 1) end of prior segment
	      for (int p=0; p<pad; p++)
		{		  
		  for (int s=0;s<ns; s++) bad[ row ].insert( s );
		  uint64_t e = prior + globals::tp_1sec * ( p + 1 )  ;
		  //std::cout << " @ " << prior << " pad (prior) " << p << " " << e << "\n";
		  estarts.push_back( e ) ;
		  ++row;
		}
	      
	      // 2) start of current segment
	      for (int p=0; p<pad; p++)
                {		  
                  for (int s=0;s<ns; s++) bad[ row ].insert( s );
		  //std::cout << " @ " << interval.start << " pad (post) " << p << " " << interval.start - globals::tp_1sec * ( pad - p ) << "\n";
		  estarts.push_back( interval.start - globals::tp_1sec * ( pad - p ) ) ;
                  ++row;
                }	      
	      
	    }
	}
      
      // track the last
      prior = interval.start;      
      
      // track epoch positions (i.e. should be position by 1 second if no gap)      
      estarts.push_back( interval.start );

      
      //
      // EEG
      //
      
      eigen_matslice_t eeg( edf, eeg_signals, interval );
      
      // signals
      const Eigen::MatrixXd & X = eeg.data_ref();
      for (int s=0; s<ns_eeg; s++)
	{
	  const bool okay = hjorth( X.col(s) , &H1(row,s), &H2(row,s), &H3(row,s) );	  
	  if ( ! okay ) bad[ row ].insert( s );
	}      
      
      
      //
      // EMG
      //

      eigen_matslice_t emg( edf, emg_signals, interval );

      // signals
      const Eigen::MatrixXd & Y = emg.data_ref();
      for (int s=0; s<ns_emg; s++)
	{
	  const bool okay = hjorth( Y.col(s) , &H1(row,ns_eeg+s), &H2(row,ns_eeg+s), &H3(row,ns_eeg+s) );
	  if ( ! okay ) bad[ row ].insert( ns_eeg+s );
	}

      
      //
      // next epoch
      //
      
      ++row;
      
    }

  if ( bad.size() > 0 ) 
    logger << "  " << bad.size() << " epochs with at least one bad (flat) channels\n";


  // check
  if ( row != expected || estarts.size() != expected )
    Helper::halt( "internal error in AROUSALS" );

  
  std::vector<std::vector<int> > good_epochs(ns), bad_epochs(ns);
  
  for (int r=0; r<expected; r++)
    {
      
      std::map<int,std::set<int> >::const_iterator qq = bad.find( r );
      
      if ( qq == bad.end() )
	{
	  for (int s=0; s<ns; s++) good_epochs[s].push_back( r );
	}
      else
	{	  
	  const std::set<int> & t = qq->second;	  
	  for (int s=0; s<ns; s++)
	    {
	      if ( t.find( s ) == t.end() )
		good_epochs[s].push_back( r );
	      else
		bad_epochs[s].push_back( r );
	    }
	}
    }
  

    
  //
  // normalize
  //
   
  for (int s=0; s<ns; s++)
    {
            
      if ( bad_epochs[s].size() == 0 )
	{
	  
	  Eigen::VectorXd v1 = H1.col(s);
	  Eigen::VectorXd v2 = H2.col(s);
	  Eigen::VectorXd v3 = H3.col(s);
	  eigen_ops::robust_scale( v1 , true, true, winsor_th  , true );
	  eigen_ops::robust_scale( v2 , true, true, winsor_th , true );
	  eigen_ops::robust_scale( v3 , true, true, winsor_th , true );      	  
	  H1.col(s) = v1;
	  H2.col(s) = v2;
	  H3.col(s) = v3;
	  
	}
      else if ( good_epochs[s].size() > 2 )
	{
	  
	  Eigen::VectorXd v1 = H1.col(s)( good_epochs[s] );
	  Eigen::VectorXd v2 = H2.col(s)( good_epochs[s] );
	  Eigen::VectorXd v3 = H3.col(s)( good_epochs[s] );

	  eigen_ops::robust_scale( v1 , true, true, winsor_th , true );
	  eigen_ops::robust_scale( v2 , true, true, winsor_th , true );
	  eigen_ops::robust_scale( v3 , true, true, winsor_th , true );      
	  
	  H1.col(s)( good_epochs[s] ).array() = v1;
	  H2.col(s)( good_epochs[s] ).array() = v2;
	  H3.col(s)( good_epochs[s] ).array() = v3;
	  
	  H1.col(s)( bad_epochs[s] ).array() = 0;
	  H2.col(s)( bad_epochs[s] ).array() = 0;
	  H3.col(s)( bad_epochs[s] ).array() = 0;
	  
	}
    }

  
  
  //
  // Get contiguous regions for arousal detection (i.e. should be sleep-only EDF)
  //

  // account for pad
  uint64_t pad_tp = pad * globals::tp_1sec;
  for (int e=0; e<expected; e++)
    estarts[e] += pad_tp;
  
  // i.e. 5 means gap between element 5 and 6 (0-based)
  std::vector<int> egaps;
  
  for (int e=1; e<expected; e++)
    if ( estarts[e] - estarts[e-1] != globals::tp_1sec )
      egaps.push_back(e-1);
  // & add the last one
  egaps.push_back( expected - 1 );
  
  // for (int e=0; e<expected; e++)
  //   {
  //     const bool g = estarts[e] - estarts[e-1] != globals::tp_1sec ;
  //     std::cout << "estart " << e << " " << estarts[e] << " " << g << "\n";
  //   }
  
  // for (int e=0; e<egaps.size(); e++)
  //   std::cout << " egaps " << e << " "<< egaps[e] << "\n";
  

  //
  // Combine statistics (but keep any EEG and EMG separate) 
  //
  
  Eigen::VectorXd H1_eeg, H2_eeg, H3_eeg;
  Eigen::VectorXd H1_emg, H2_emg, H3_emg;
  
  if ( ns_eeg >= 1 )
    {
      H1_eeg = H1.leftCols(ns_eeg).rowwise().mean();
      H2_eeg = H2.leftCols(ns_eeg).rowwise().mean();
      H3_eeg = H3.leftCols(ns_eeg).rowwise().mean();
    }
  
  if ( ns_emg >= 1 )
    {
      H1_emg = H1.rightCols(ns_emg).rowwise().mean();
      H2_emg = H2.rightCols(ns_emg).rowwise().mean();
      H3_emg = H3.rightCols(ns_emg).rowwise().mean();
    }

  
  //
  // Find peaks : channel specific
  //

  const int lag = lag_sec ; // value in samples, but here we've fixed to 1 Hz outputs
  
  for (int s=0; s<ns; s++)
    {
      
      const bool is_eeg = s < ns_eeg; 
      
      for (int h=1; h<=3; h++)
	{
	  
	  // urgh, need to get rid of all these stupid copies...      
	  std::vector<double> vec;
	  
	  if      ( h == 1 ) vec = eigen_ops::copy_vector( H1.col(s) );
	  else if ( h == 2 ) vec = eigen_ops::copy_vector( H2.col(s) );
	  else if ( h == 3 ) vec = eigen_ops::copy_vector( H3.col(s) );
	  
	  // detect peaks / arousals
	  
	  // const std::vector<double> & x , int lag , double threshold , double influence = 0 , 
	  // 	int mindur = 0 , double max = 0 , 
	  // 	double threshold2 = 0 , int mindur2 = 0 , 
	  // 	bool noneg = false , 
	  // 	std::vector<interval_t> * regions = NULL , 
	  // 	std::vector<int> * top_pks = NULL , 
	  // 	bool verbose = false );

	  std::vector<double> locsmooth;
	  
	  std::vector<int> z = MiscMath::batched_smoothedZ( vec , egaps,
							    lag , threshold , influence , mindur ,
							    0, 0, 0, true , &locsmooth ); // +ve peaks only

	  //std::cout << "locsmooth = " << z.size() << " " << vec.size() << " " << locsmooth.size() << " " << H1.rows() << " "<< H1.cols() << "\n";

	  // update signal w/ smoothed variant	  
	  // & normalize / winsorize again
          if      ( h == 1 )
	    {
	      Eigen::VectorXd vec = eigen_ops::copy_array( locsmooth );
	      eigen_ops::robust_scale( vec , true, true, winsor_th , true );
	      H1.col(s) = vec;
	    }
	  else if ( h == 2 ) 
	    {
	      Eigen::VectorXd vec = eigen_ops::copy_array( locsmooth );
	      eigen_ops::robust_scale( vec , true, true, winsor_th , true );
	      H2.col(s) = vec;
	    }
	  else if ( h == 3 ) 
	    {
	      Eigen::VectorXd vec = eigen_ops::copy_array( locsmooth );
	      eigen_ops::robust_scale( vec , true, true, winsor_th , true );
	      H3.col(s) = vec;
	    }
	  
	  // Compile to intervals      
	  std::vector<interval_t> arousals = combine( z, estarts , mindur, maxdur );
	  
	  const std::string l1 = "H" + Helper::int2str( h );
	  const std::string l2 = is_eeg ? eeg_signals.label(s) : emg_signals.label( s - ns_eeg );
	  
	  logger << "  adding " << arousals.size() << " arousals for " << l2 << " based on H" << h << "\n";
	  
	  for (int i=0; i<arousals.size(); i++)
	    a->add( l1, arousals[i], l2 );
	  
	}
      
    }
  
  //
  // Analysis of combined sigs
  //
  
  if ( ns_eeg )
    {
      
    }
  
  
  //
  // Write new signals
  //

  if ( write_combined_metrics )
    {

      const int sr = 1;
      
      if ( ns_eeg )
	{
	  
	  const std::string slab = new_sig_prefix + "_EEG";
	  
	  std::vector<double> h1 = eigen_ops::copy_vector( H1_eeg );
	  std::vector<double> h2 = eigen_ops::copy_vector( H2_eeg );
	  std::vector<double> h3 = eigen_ops::copy_vector( H3_eeg );
	  
	  // zero-pad as needed
	  // h1.insert( h1.begin(), pad1 , 0 );
	  // h2.insert( h2.begin(), pad1 , 0 );
	  // h3.insert( h3.begin(), pad1 , 0 );
	  // for (int j=0;j<pad2;j++)
	  //   {
	  //     h1.push_back( 0 );
	  //     h2.push_back( 0 );
	  //     h3.push_back( 0 );
	  //   }

	  logger << "  writing EEG trace to " << slab << " (H1, H2 & H3)\n";
	  
	  // write
	  edf.add_signal( slab + "_H1" , sr , h1 );
	  edf.add_signal( slab + "_H2" , sr , h2 );
	  edf.add_signal( slab + "_H3" , sr , h3 );
	}

      if ( ns_emg )
	{
	  
	  const std::string slab = new_sig_prefix + "_EMG";
	  
	  std::vector<double> h1 = eigen_ops::copy_vector( H1_emg );
	  std::vector<double> h2 = eigen_ops::copy_vector( H2_emg );
	  std::vector<double> h3 = eigen_ops::copy_vector( H3_emg );
	  
	  // zero-pad as needed
	  // h1.insert( h1.begin(), pad1 , 0 );
	  // h2.insert( h2.begin(), pad1 , 0 );
	  // h3.insert( h3.begin(), pad1 , 0 );
	  // for (int j=0;j<pad2;j++)
	  //   {
	  //     h1.push_back( 0 );
	  //     h2.push_back( 0 );
	  //     h3.push_back( 0 );
	  //   }
	  
	  logger << "  writing EMG trace to " << slab << " (H1, H2 & H3)\n";
	  
	  // write
	  edf.add_signal( slab + "_H1" , sr , h1 );
	  edf.add_signal( slab + "_H2" , sr , h2 );
	  edf.add_signal( slab + "_H3" , sr , h3 );
	}
      
    }
  

  if ( write_channel_specific_metrics )
    {
      for (int s=0; s<ns; s++)
	{
	  
	  // add signal (will always be 1 Hz)
	  const std::string slab = new_sig_prefix + ( s < ns_eeg ? eeg_signals.label(s) : emg_signals.label( s - ns_eeg ) ) ;
	  const int sr = 1;
	  
	  std::vector<double> h1 = eigen_ops::copy_vector( H1.col(s) );
	  std::vector<double> h2 = eigen_ops::copy_vector( H2.col(s) );
	  std::vector<double> h3 = eigen_ops::copy_vector( H3.col(s) );
	  
	  // zero-pad as needed
	  // h1.insert( h1.begin(), pad1 , 0 );
	  // h2.insert( h2.begin(), pad1 , 0 );
	  // h3.insert( h3.begin(), pad1 , 0 );
	  // for (int j=0;j<pad2;j++)
	  //   {
	  //     h1.push_back( 0 );
	  //     h2.push_back( 0 );
	  //     h3.push_back( 0 );
	  //   }
	  
	  // write
	  edf.add_signal( slab + "_H1" , sr , h1 );
	  edf.add_signal( slab + "_H2" , sr , h2 );
	  edf.add_signal( slab + "_H3" , sr , h3 );
	  
	  //
	  // next signal
	  //
	}
    }
  
   
}



bool arousals_t::hjorth( const Eigen::VectorXd & x , double * activity , double * mobility , double * complexity , const bool mean_center ) const
{

  const int n = x.size();  
  if ( n == 0 ) return false;

  const Eigen::VectorXd dxV = x.tail(n-1) - x.head(n-1);
  const Eigen::VectorXd ddxV = dxV.tail(n-2) - dxV.head(n-2);
  
  const double mx2 = ( mean_center ? (x.array()-x.mean()).matrix().squaredNorm() : x.squaredNorm() ) / double(n);
  const double mdx2 = dxV.squaredNorm() / double(n-1);
  const double mddx2 = ddxV.squaredNorm() / double(n-2);
  
  *activity   = mx2;
  *mobility   = mdx2 / mx2;
  *complexity = sqrt( mddx2 / mdx2 - *mobility );
  *mobility   = sqrt( *mobility );

  if ( ! Helper::realnum( *activity ) ) return false;
  if ( ! Helper::realnum( *mobility ) ) return false;
  if ( ! Helper::realnum( *complexity ) ) return false;
  return true;
}


std::vector<interval_t> arousals_t::combine( const std::vector<int> & x ,
					     const std::vector<uint64_t> & tp ,					     
					     const double mindur ,
					     const double maxdur ) const
{
  
  // tp is from estarts above (i.e. the 1 Hz signal starts)
  //  but is the start of the e.g. 3 second epoch; so to account
  //  for padding, we need to shift everything by 'pad' to get the 'right'
  //  value (i.e. aligned w/ the true clock)
  
  // x is the 0/1 from smoothedZ() peak finder
  
  std::vector<interval_t> ints;

  if ( x.size() != tp.size() ) Helper::halt( "internal error in arousals_t::combine()" );
  
  bool px = x[0];
  uint64_t ptp = tp[0];

  const int n = x.size();
  
  for (int i=1; i<n; i++)
    {

      // end of region , or a gap?
      const bool gap = tp[i] - tp[i-1] != globals::tp_1sec;
      const bool end = i == n - 1 ;

      // close out if gap or end?
      if ( px && ( gap || end ) )
	{	  

	  const double dur = ( tp[i] - ptp ) * globals::tp_duration;
	  if ( dur >= mindur && dur <= maxdur )
	    ints.push_back( interval_t( ptp , tp[i] ) );
	  
	  px = false;
	  // allow for (unlikely) immediate restart)
	  // but if at end, we can quit
	  if (end) continue;
	}
      
      // start of a new arousal
      if ( x[i] == 1 && ! px )
	{
	  px = true;
	  ptp = tp[i];
	}
      else if ( x[i] == 0 && px ) // else end?
	{
	  px = false;

	  const double dur = ( tp[i] - ptp ) * globals::tp_duration;
	  if ( dur >= mindur && dur <= maxdur )
	    ints.push_back( interval_t( ptp , tp[i] ) );
	}      
    }
  
  return ints;
}



