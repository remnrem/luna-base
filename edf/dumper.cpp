
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

#include "edf.h"
#include "tal.h"

#include "annot/annot.h"
#include "helper/helper.h"
#include "main.h"
#include "defs/defs.h"
#include "db/db.h"

extern writer_t writer;

extern annotation_set_t annotations;

void edf_t::record_dumper( param_t & param )
{

  bool show_signals = true;
  bool show_annots  = true;
  
  if ( param.has( "no-signals" ) ) show_signals = false; 
  if ( param.has( "no-annots" ) ) show_annots = false;


  //
  // Annotations present? (i.e. already loaded)
  //
  
  std::vector<std::string> annots = timeline.annotations.names();
  
  int r = timeline.first_record();
  
  while ( r != -1 )
    {
      
      std::cout << "Record " << r+1 << " of " << header.nr_all << " total (" << header.nr << " retained)\n";

      //
      // Interval for this record
      //

      interval_t interval = timeline.record2interval(r); 
      
      //
      // Get annotations
      //
      
      if ( show_annots ) 
	{
	  std::cout << "Generic Annotatons-----------------------\n";
	  
	  for (int a=0;a<annots.size();a++)
	    {
	      
	      annot_t * annot = timeline.annotations( annots[a] );
	      
	      annot_map_t events = annot->extract( interval );
	      
	      annot_map_t::const_iterator ii = events.begin();
	      while ( ii != events.end() )
		{	      
		  
		  const instance_idx_t & instance_idx = ii->first;
		  const instance_t * instance = ii->second;
		  
		  std::cout << annot->name << "\t" 
			    << instance_idx.id << "\t"
			    << instance_idx.interval.as_string() ;
		  
		  instance_table_t::const_iterator jj = instance->data.begin();
		  while ( jj != instance->data.end() )
		    {
		      std::cout << "\t" << jj->first 
				<< "[" << globals::type_name[ jj->second->atype() ] << "]="
				<< jj->second->text_value();
		      ++jj;
		    }

		  std::cout << "\n";
		  ++ii;
		}
	    }
	  
	  //
	  // EDF annotations
	  //
	  
	  std::cout << "EDF Annotations--------------------------\n";
	  
	  for ( int s = 0 ; s < header.ns; s ++ )
	    {
	      
	      if ( header.is_annotation_channel( s ) )
		{	      
		  
		  tal_t t = tal( s , r );
		  
		  std::cout << "Signal " << s+1 << " " 
			    << header.label[s] << "\n"
			    << t << "\n\n";
		  
		}
	      
	      
	    } // next signal
	  
	}

      //
      // Get data 
      //

      if ( show_signals ) 
	{

	  std::cout << "Data signals-----------------------------\n";
	  
	  for ( int s = 0 ; s < header.ns; s ++ )
	    {
	      
	      std::cout << "s = " << s << "\n";
	      
	      if ( header.is_data_channel( s ) )
		{
		  
		  std::cout << "interval = " << interval << "\n";
		  
		  slice_t data( *this , s , interval );
		  
		  const std::vector<double> * d = data.pdata();
		  
		  const std::vector<uint64_t> * tp = data.ptimepoints();
		  
		  const int n = d->size(); 
		  
		  std::cout.precision(8);
		  
		  for (int i=0;i<n;i++)
		    {	      
		      
		      uint64_t sec = (*tp)[i] / globals::tp_1sec;
		      uint64_t rem = (*tp)[i] - ( sec * globals::tp_1sec );
		      double   rem2 = (double)rem / (double) globals::tp_1sec ;
		      double   sec2 = (double)sec + rem2;
		      
		      std::cout << "RECORD-DUMP" << "\t" 
				<< header.label[s] << "\t"
				<< "rec=" << r << "\t"
				<< (i+1) << "/" << n << "\t"
				<< "tp=" << (*tp)[i] << "\t"
				<< "sec=" << (*tp)[i] * globals::tp_duration << "\t"
				<< sec << "\t"
				<< rem << "\t"
				<< rem2 << "\t"
				<< sec2 << "\t"
				<< (*d)[i] << "\n";
		    }
		}
	    } // next signal
	  
	}

      r = timeline.next_record( r );
      
    } // next record
  
}

void edf_t::data_dumper( const std::string & signal_labels , const param_t & param )
{
  

  //
  // Attach signals
  //
  
  signal_list_t signals = header.signal_list( signal_labels );

  if ( signals.size() != 1 ) 
    Helper::halt( "DUMP currently only for single channels; see MATRIX" );
  

  //
  // Options
  //

  bool hms = param.has("hms");
  clocktime_t starttime( header.starttime );
  if ( ! starttime.valid ) hms = false;

  bool sec = param.has("sec");

  bool only_signal = param.has("minimal");

  //
  // What annotations are present? (i.e. already loaded)
  //
  
  std::vector<std::string> annots = timeline.annotations.names();
  

  //
  // Point to first epoch
  //
  
  timeline.first_epoch();
  

  //
  // for each each epoch 
  //
  
  while ( 1 ) 
    {
      
      //
      // Get next epoch
      //
      
      int epoch = timeline.next_epoch();      
	  
      if ( epoch == -1 ) break;
      
      interval_t interval = timeline.epoch( epoch );
      
      //
      // Collate 'header'
      //
 
      // ID and time point

      std::stringstream ss;

      if ( ! only_signal ) 
	ss << "DUMP\t" 
	   << id << "\t"
	   << "epoch=" << epoch + 1 ;  // all input/output is 1-based for epochs
      
      //
      // Get annotations
      //

      if ( ! only_signal ) 
	{
	  std::map<std::string,std::set<std::string> > atxt;
	  
	  for (int a=0;a<annots.size();a++)
	    {
	      
	      annot_t * annot = timeline.annotations( annots[a] );
	      
	      annot_map_t events = annot->extract( interval );
	      
	      // collapse
	      
	      annot_map_t::const_iterator ii = events.begin();
	      while ( ii != events.end() )
		{	      

		  const instance_idx_t & instance_idx = ii->first;
		  const instance_t * instance = ii->second;
		  
		  instance_table_t::const_iterator jj = instance->data.begin();
		  while ( jj != instance->data.end() )
		    {
		      if (  jj->second == NULL ) 
			atxt[ jj->first ].insert( "." );
		      else
			atxt[ jj->first ].insert( jj->second->text_value() );
		      ++jj;
		    }
		  ++ii;
		}
	    }
	  
	  // display
	  
	  ss << "\t";
	  
	  std::map<std::string,std::set<std::string> >::iterator ai = atxt.begin();
	  while ( ai != atxt.end() )
	    {
	      if ( ai == atxt.begin() ) ss << "epoch-ann:"; else ss << ";";
	      ss << ai->first << "=";
	      
	      std::set<std::string>::iterator si = ai->second.begin();
	      while ( si != ai->second.end() )
		{
		  if ( si != ai->second.begin() ) ss << ",";
		  ss << *si;
		  ++si;
		}
	      ++ai;
	    }
	}

    
      //
      // Get data 
      //

      slice_t data( *this , signals(0) , interval );


      //
      // assumes only a single signal... okay for now, 
      // as sampling rate may be different in any case
      //

      const std::vector<double> * d = data.pdata();
      const std::vector<uint64_t> * tp = data.ptimepoints();

      //
      // Now display all data points within this EPOCH
      //
    
      const int n = d->size();  // number of input points
      
      for (int i=0;i<n;i++)
	{	      

	  if ( ! only_signal ) 
	    {
	      std::cout << ss.str() << "\t" 
			<< "tp=" << (*tp)[i] ;
	    
	      if ( sec ) std::cout << "\t" << (*tp)[i] / (double)globals::tp_1sec;
	      if ( hms ) 
		{
		  double tp_sec =  (*tp)[i] / (double)globals::tp_1sec;
		  clocktime_t present = starttime;
		  present.advance( tp_sec / 3600.0 );
		  std::cout << "\t" << present.as_string();
		}
	  
	      // signal 	  
	      std::cout << "\t" << (*d)[i] << "\n";
	    }
	  else
	    {
	      std::cout << (*d)[i] << "\n";
	    } 
	}
      
    } // next epoch
}




void edf_t::data_epoch_dumper( param_t & param , std::set<std::string> * selected_annots )
{
  
  std::cerr << " listing " << timeline.num_total_epochs() << " epochs, of which " 
	    << timeline.num_epochs() << " are unmasked\n";

  bool show_times  = param.has("show-times");
  bool hide_masked = param.has("hide-masked");
  bool indiv_only  = param.has("indiv-only");
  bool show_indiv  = param.has("indiv") || indiv_only;
  bool show_epoch  = ! indiv_only;
  bool hide_false  = param.has("show-all") ? false : true; 


  //
  // What annotations are present? (i.e. already loaded)
  //
  
  std::vector<std::string> annots = timeline.annotations.names();
  
  //
  // Point to first epoch
  //
  
  timeline.first_epoch();
  

  //
  // Summary statistics for this individual
  //
  
  int ecnt = 0; // all epochs

  std::map<std::string,int> ecnts;
  
  //
  // Set up epoch annotations
  //
  
  std::set<std::string> epoch_annotations = timeline.epoch_annotations();

  bool has_epoch_annotations = epoch_annotations.size() > 0 ;
  

  //
  // for each each epoch 
  //
  
  while ( 1 ) 
    {
      
      //
      // Get next epoch
      //
      
      int epoch = timeline.next_epoch_ignoring_mask();      

      if ( epoch == -1 ) break;
      
      if ( hide_masked && timeline.masked_epoch( epoch ) )
	continue;


      interval_t interval = timeline.epoch( epoch );

      //
      // Output 
      //
      


      //
      // Collate 'header'
      //
 
      // ID, epoch #, mask setting and time point
      
      if ( show_epoch )
	{
	  
	  writer.epoch( timeline.display_epoch( epoch ) );
	  
	  writer.var( "E1" , "Epoch number ignoring original structure" );
	  writer.var( "MASK" , "Masked epoch (1=Y)" );
	  writer.var( "INTERVAL" , "Interval start-stop (secs)" );
	  
	  writer.value( "E1" , epoch+1 );
	  writer.value( "MASK" , timeline.masked_epoch( epoch ) ? 1 : 0 );
	  writer.value( "INTERVAL" , interval.as_string() );
	  
	}
      

      //
      // Collapsed 'bool' epoch level annotations
      //
      
      if ( has_epoch_annotations )
	{

	  bool any_annot = false;

	  std::set<std::string>::const_iterator aa = epoch_annotations.begin();
	  while ( aa != epoch_annotations.end() )
	    {	      
	      
	      bool has_annot = timeline.epoch_annotation( *aa , epoch );
	      
	      if ( show_epoch && has_annot )
		{
		  writer.level( *aa , globals::annot_strat );
		  writer.var( "PRESENT" , "Epoch has annotation?" );		      		  
		  writer.value( "PRESENT" , has_annot );
		}
	      
	      if ( has_annot ) 
		{
		  any_annot = true;
		  ecnts[ *aa ]++;
		}

	      ++aa;
	    }
	  
	  if ( any_annot ) ++ecnt;
	  
	}
      

      //
      // Display full (values/times) for annotations
      //

      std::map<std::string,std::set<std::string> > atxt;
      std::map<std::string,std::map<std::string,std::string> > atimes; // time-points

      for (int a=0;a<annots.size();a++)
	{
	  
	  if ( selected_annots != NULL && selected_annots->find( annots[a] ) == selected_annots->end() ) continue;

	  annot_t * annot = timeline.annotations( annots[a] );

	  annot_map_t events = annot->extract( interval );
	  
	  // collapse
	  
	  annot_map_t::const_iterator ii = events.begin();
	  while ( ii != events.end() )
	    {	 

	      const instance_idx_t & instance_idx = ii->first;
	      const instance_t * instance = ii->second;

	      instance_table_t::const_iterator jj = instance->data.begin();
	      while ( jj != instance->data.end() )
		{
		  
		  std::string s = ".";
		  
		  if (  jj->second != NULL ) 
		    s = jj->second->text_value() ;

		  atxt[ jj->first ].insert( s );
		  
		  if ( show_times )
		    {		      
		      std::string & t = atimes[ jj->first ][ s ];
		      if ( t == "" ) t = ii->first.interval.as_string() ;
		      else t += "," + ii->first.interval.as_string() ;
		    }

		  ++jj;
		}
	      ++ii;
	    }
	}


      // display
      
      std::map<std::string,std::set<std::string> >::iterator ai = atxt.begin();
      while ( ai != atxt.end() )
	{
	  
	  if ( show_epoch )
	    {
	      
	      writer.level( ai->first , globals::annot_strat );
	      
	      // optionally, times
	      
	      std::map<std::string,std::string>::const_iterator tii;
	      if ( show_times )
		{
		  std::map<std::string,std::map<std::string,std::string> >::const_iterator ti = atimes.find( ai->first );
		  tii = ti->second.begin();
		}

	      int acnt = 0;

	      std::set<std::string>::iterator si = ai->second.begin();
	      while ( si != ai->second.end() )
		{
		  writer.level( ++acnt , globals::count_strat );
		  
		  writer.var( "ANNOT" , "Annotation" );
		  writer.value( "ANNOT" , *si );
		  
		  if ( show_times )
		    {
		      writer.var( "ANNOT_TIME" , "Annotation timestamp" );
		      writer.value( "ANNOT_TIME" , tii->second );
		    }
		  
		  ++si;
		}
	      writer.unlevel( globals::count_strat );
	    }

	  ++ai;
	}
      
      if ( show_epoch ) 
	writer.unlevel( globals::annot_strat );
      
      
    } // next epoch

  
  if ( show_epoch ) 
    writer.unepoch();



  //
  // summary
  //
  
  if ( show_indiv )
    {

      writer.var( "N" , "Total number of epochs" );
      writer.var( "NE_FLAGGED" , "Total number of flagged epochs" );
            
      writer.value( "N" , timeline.num_epochs() );
      writer.value( "NE_FLAGGED" , ecnt );

      std::map<std::string,int>::const_iterator ee = ecnts.begin();
      while ( ee != ecnts.end() ) 
	{
	  writer.level( ee->first , globals::annot_strat );
	  writer.var( "N_ANNOT" , "Number of annotation instances" );	  
	  writer.value( "N_ANNOT" , ee->second );
	  ++ee;
	}
      writer.unlevel( globals::annot_strat );
    }

}






void edf_t::epoch_matrix_dumper( param_t & param , std::set<std::string> * selected_annots )
{
  
  // dump output to STDOUT rather than a DB;  3 rows per subject

  // ID, # epochs, # epoch length (sec) , sample rate, # annotations , 
  // ID, ne, len, fs, na, ns 
  // na * ne : { all annots ( grouped by epoch ) ... }    
  // ns * ne * len*fs : { all signal points (grouped by signal, then epoch, then time-point } 
  
  std::cerr << " dumping  " << timeline.num_epochs() << " unmasked epochs in matrix-format to stdout\n";
  
  // Standard output is header row: ID epoch elapsed-time(sec) annot1 annot2 ... signal1 signal2 signal3 ... 
  // and then tab-delimited entries below, each line is a single data point.

  // An alternate format is used for reading into R for subsequent processing, making this easier by 
  // splitting elements into spearate rows: 

  // I  ID  #epochs  #length-of-epoch-seconds  SR  #annots  #signals
  // SL signal labels
  // EL epoch labels
  // AH annotation header

  //  epoch level annotations:
  // A { epoch 1, bool for each annot } { epoch 2, bool for each annot per datapoint }

  //  raw signals, ordered by signal #
  // S { signal1 alldatapoints } { signal2 all datapoints } ...


  bool alternative_format = param.has( "format2" );
  
 
  //
  // Get signals
  //  
  
  std::string signal_label = param.requires( "signal" );   
  
  signal_list_t signals = header.signal_list( signal_label );  
  
  const int ns = signals.size();

  if ( ns == 0 ) return;


  //
  // Check FS for all signals
  //

  int fs = header.sampling_freq( signals( 0 ) ) ;
  
  for (int s=1; s<ns; s++) 
    {
      if ( header.sampling_freq( signals(s) ) != fs ) 
	Helper::halt( "EPOCH-MATRIX requires uniform sampling rate across signals" ); 
    }
  
      
  //
  // Which annotations are present? (i.e. already loaded)
  //
  
  
  //
  // Point to first epoch
  //
  
  timeline.first_epoch();

  const int ne = timeline.num_epochs();

  //
  // Set up epoch annotations
  //
  
  std::set<std::string> epoch_annotations = timeline.epoch_annotations();
  
  bool has_epoch_annotations = epoch_annotations.size() > 0 ;
   
  const int na = selected_annots == NULL ? 0 : selected_annots->size();


  if ( alternative_format )
    {
  
      //
      // Output
      //
      
      std::cout << "I" << "\t" 
		<< id << "\t" 
		<< ne << "\t"
		<< timeline.epoch_length() << "\t"
		<< fs << "\t"
		<< na << "\t"
		<< ns << "\n";
      
      //
      // Signal labels
      //
      
      std::cout << "SL";
      for (int s = 0 ; s < ns ; s++ )
	std::cout << "\t" << header.label[ signals(s) ] ;
      std::cout << "\n";

      //
      // Epoch labels
      //
      
      std::cout << "EL";
      timeline.first_epoch();       
      while ( 1 ) 
	{
	  int epoch = timeline.next_epoch();	   
	  if ( epoch == -1 ) break;
	  std::cout << "\t" << timeline.display_epoch( epoch ) ;
	}
      std::cout << "\n";
      

      //
      // Output annotations 
      //
      
      if ( selected_annots != NULL )
	{
	  
	  std::vector<std::string> annots = timeline.annotations.names();
	  
	  std::cout << "AH";
	  
	  std::set<std::string>::const_iterator aa = selected_annots->begin();
	  while ( aa != selected_annots->end() )
	    {
	      std::cout << "\t" << *aa;
	      ++aa;
	    }
	  std::cout << "\nA";
	  
	  timeline.first_epoch();       
  
	  while ( 1 ) 
	    {
	      
	      //
	      // Get next epoch
	      //
	      
	      int epoch = timeline.next_epoch();
	   
	      if ( epoch == -1 ) break;
	   
	      interval_t interval = timeline.epoch( epoch );
	   
	      //
	      // Annotations
	      //
	      
	      std::set<std::string>::const_iterator aa = selected_annots->begin();
	      while ( aa != selected_annots->end() )
		{
		  
		  //std::cout << "annots[a] = " << annots[a] << "\n";
		  
		  annot_t * annot = timeline.annotations( *aa );
		  
		  if ( annot == NULL ) 
		    {
		      std::cout << "\t" << 0;
		      ++aa;
		      continue;
		    }

		  annot_map_t events = annot->extract( interval );
		  
		  bool has_annot = events.size() ;
		  
		  std::cout << "\t" << ( has_annot ? 1 : 0 ) ;
		  
		  ++aa;
		}
	      
	      // next epoch
	    }
	  std::cout << "\n";
	}
      

      //
      // Output signals
      //
      
      std::cout << "S";
   
 
      //
      // For each signal
      //
      
      for (int s = 0 ; s < ns ; s++ )
	{	  
	  
	  timeline.first_epoch();
       
	  while ( 1 ) 
	    {
	      
	      //
	      // Get next epoch
	      //
	      
	      int epoch = timeline.next_epoch();
	      
	      if ( epoch == -1 ) break;
	      
	      interval_t interval = timeline.epoch( epoch );
	      
	      slice_t slice( *this , signals(s) , interval );
	      const std::vector<double> * signal = slice.pdata();
	      const int np = signal->size();
	      for (int i=0;i<np;i++) std::cout << "\t" << (*signal)[i];	  
	      
	    } // Next epoch           
	  
	} // Next signal
      
      // all done
      std::cout << "\n";

      return;
      
    }

  
  //
  // Standard matrix format
  //
  

  //
  // Header
  //

  bool include_hms = param.has("hms");

  std::cout << "ID\tE\tT";
  
  if ( include_hms ) std::cout << "\tHMS";
  
  clocktime_t starttime( header.starttime );
  bool invalid_hms = ! starttime.valid;
  
 
  // Annots

  // annots we actually have for this EDF and so should be checked
  // if NULL, means this annot is not present
  std::vector<annot_t*> annotlist;  

  if ( selected_annots != NULL )
    {
      std::set<std::string>::const_iterator aa = selected_annots->begin();
      while ( aa != selected_annots->end() )
	{
	  // header
	  std::cout << "\t" << *aa;
	  
	  // track
	  annot_t * annot = timeline.annotations( *aa );
	  annotlist.push_back( annot ); // NULL if annot not present
	  
	  ++aa;
	}	    
    }
  

  // Signals
  for (int s = 0 ; s < ns ; s++ )
    std::cout << "\t" << header.label[ signals(s) ] ;
  
  std::cout << "\n";

 
    
  //
  // Iterate over epochs, display  
  //

  
  timeline.first_epoch();       
  while ( 1 ) 
    {

      int epoch = timeline.next_epoch();	   
      if ( epoch == -1 ) break;      

      // get all signals for this epoch
      
      interval_t interval = timeline.epoch( epoch );
      std::vector<std::vector<double> > sigdat( ns );
      
      // track time-points, ie. may be a discontinuous file
      std::vector<uint64_t> tp;

      for (int s = 0 ; s < ns ; s++ )
	{	  	  
	  slice_t slice( *this , signals(s) , interval );
	  const std::vector<double> * signal = slice.pdata();
	  sigdat[s] = *signal;
	  if ( s== 0 ) tp = *slice.ptimepoints();
	}
      
      // now iterate over all time-points
      
      const int np = sigdat[0].size();
      
      for (int t=0;t<np;t++)
	{
	  
	  double tp_sec = tp[t] / (double)globals::tp_1sec; 
	  
	  // output rows
	  std::cout << id << "\t"
		    << timeline.display_epoch( epoch ) << "\t"
		    << tp_sec ;

	  if ( include_hms )
	    {
	      clocktime_t present = starttime;
	      present.advance( tp_sec / 3600.0 );
	      std::cout << "\t" << present.as_string();
	    }

	  // annots
	  
	  if ( selected_annots != NULL )
	    {
	      
	      // get exact point
	      interval_t interval2 = interval_t( tp[t] , tp[t] );
	      
	      for (int a=0;a<annotlist.size();a++)
		{
		  if ( annotlist[a] == NULL ) 
		    std::cout << "\t0";
		  else
		    {
		      annot_map_t events = annotlist[a]->extract( interval2 );
		      bool has_annot = events.size() ;		      
		      std::cout << "\t" << ( has_annot ? 1 : 0 ) ;
		    }		      
		} // next annot
	    } 

	  
	  // signals
	  for (int s=0;s<ns;s++) std::cout << "\t" << sigdat[s][t];

	  
	  // done, next row/time-point
	  std::cout << "\n";

	  
	} // next time-point
      
    } // next epoch
  
  
}

