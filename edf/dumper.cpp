
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
#include "slice.h"
#include "tal.h"

#include "annot/annot.h"
#include "helper/helper.h"
#include "helper/logger.h"

#include "eval.h"
#include "defs/defs.h"
#include "db/db.h"

extern writer_t writer;

extern logger_t logger;

extern annotation_set_t annotations;


void edf_t::record_table( param_t & param )
{
  
  // iterate over each record

  int r = timeline.first_record();
  
  while ( r != -1 )
    {
                      
      // Interval for this record
      
      interval_t interval = timeline.record2interval(r); 

      // basic information

      std::cout << "RECS\t"
		<< id << "\t";

      std::cout << r+1 << "\t" 
		<< header.nr << "/" << header.nr_all ;
      
      std::cout << "\t" << interval.as_string() ;

      // epoch information?
      
      if ( timeline.epoched() )
	{
	  std::cout << "\t";
   	  std::map<int,bool>  epochs = timeline.spanning_epoch_masks( r );
	  std::map<int,bool>::const_iterator ii = epochs.begin();
	  while ( ii != epochs.end() ) 
	    {

	      interval_t epoch_interval = timeline.epoch(ii->first );
	      std::cout << " ";
	      if ( ii->second ) std::cout << "[" ;
	      std::cout << timeline.display_epoch( ii->first ) ;
	      std::cout << ";" << epoch_interval.as_string() ;
	      if ( ii->second ) std::cout << "]" ;
	      ++ii;
	    }

	}
      
      // next record

      std::cout << "\n";

      r = timeline.next_record( r );

    }
  
  
}

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
		      
// 		      uint64_t sec = (*tp)[i] / globals::tp_1sec;
// 		      uint64_t rem = (*tp)[i] - ( sec * globals::tp_1sec );
// 		      double   rem2 = (double)rem / (double) globals::tp_1sec ;
// 		      double   sec2 = (double)sec + rem2;
		      
		      std::cout << "RECORD-DUMP" << "\t" 
				<< header.label[s] << "\t"
				<< "rec=" << r << "\t"
				<< (i+1) << "/" << n << "\t"
				<< (*tp)[i] << "\t"
				<< (*tp)[i] * globals::tp_duration << "\t"
// 				<< sec << "\t"
// 				<< rem << "\t"
// 				<< rem2 << "\t"
// 				<< sec2 << "\t"
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
  // output precision
  //

  std::cout.precision(8);

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
	    
	      if ( sec ) 
		std::cout << "\t" 
			  << (*tp)[i] / (double)globals::tp_1sec;

	      if ( hms ) 
		{
// 		  double tp_sec =  (*tp)[i] / (double)globals::tp_1sec;
 		  clocktime_t present = starttime;
// 		  present.advance( tp_sec / 3600.0 );
// 		  std::cout << "\t" << present.as_string();

		  interval_t now( (*tp)[i] , (*tp)[i]+1LLU );
		  std::string t1, t2;
		  if ( Helper::hhmmss( present , now , &t1,&t2 , 5 ) ) 
		    std::cout << "\t" << t1;
		  else
		    std::cout << "\t.";

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
  
  
  bool show_times  = param.has( "show-times" );
  
  bool hide_masked = param.has( "hide-masked" );
  
  
  //
  // What annotations are present? (i.e. already loaded)
  //
  
  std::vector<std::string> annots = timeline.annotations.names();
  
  //
  // Point to first epoch
  //
  
  timeline.first_epoch();
  
  
  std::cerr << " listing " << timeline.num_total_epochs() << " epochs, of which " 
	    << timeline.num_epochs() << " are unmasked\n";
  
  
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
      
      writer.epoch( timeline.display_epoch( epoch ) );
      
      //	  writer.var( "E1" , "Epoch number ignoring original structure" );
      writer.var( "MASK" , "Masked epoch (1=Y)" );
      writer.var( "INTERVAL" , "Interval start-stop (secs)" );
      
      //writer.value( "E1" , epoch+1 );
      writer.value( "MASK" , timeline.masked_epoch( epoch ) ? 1 : 0 );
      writer.value( "INTERVAL" , interval.as_string() );
	  

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
	      
	      if ( has_annot )
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
  
  writer.unlevel( globals::annot_strat );
  
  writer.unepoch();


  //
  // summary
  //
  
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






void edf_t::epoch_matrix_dumper( param_t & param )
{
  
  //
  // requires output to a specific file
  //
  
  std::string filename = param.requires( "file" );
  
  std::ofstream OUT( filename.c_str() , std::ios::out );
  
  //
  // Minimal output?
  //

  bool minimal = param.has( "min" ) || param.has( "minimal" );

  //
  // Set up annotations: both interval and epoch level
  //
  
  bool show_annots = param.has( "annot" );
  
  std::map<std::string,int> atype; // 0 not found, 1 interval, 2 epoch
   
  int na_int = 0 , na_epoch = 0 , na = 0;

  int ns_data = 0;

  if ( minimal ) show_annots = false;
  
  if ( show_annots ) 
    {
      
      std::vector<std::string> a = param.strvector( "annot" ); // expects comma-delimited list
      
      for (int i=0;i<a.size();i++)
	{
	  
	  if ( timeline.annotations( a[i] ) != NULL ) // is this an interval annotation? 
	    { atype[ a[i] ] = 1; ++na_int; } 
	  else if ( timeline.epoch_annotation( a[i] ) ) // or an epoch-annotation?
	    { atype[ a[i] ] = 2; ++na_epoch; } 
	  else
	    atype[ a[i] ] = 0;
	  
	}
      
      na = na_int + na_epoch;
      
    }
  
  
  
  //
  // Standard format (header row, then columns)
  // 
  
  // Standard output is header row: ID epoch elapsed-time(sec) annot1 annot2 ... signal1 signal2 signal3 ... 
  // and then tab-delimited entries below, each line is a single data point.
  
  
  // ID, # epochs, # epoch length (sec) , sample rate, # annotations , 
  // ID, ne, len, fs, na, ns 
  // na * ne : { all annots ( grouped by epoch ) ... }    
  // ns * ne * len*fs : { all signal points (grouped by signal, then epoch, then time-point } 
  

  //
  // Alternate output format
  //
  //
  // For use with reading into R for subsequent processing, making this easier by 
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

  if ( minimal ) alternative_format = false;
  
  
  //
  // Get signals
  //  
  
  std::string signal_label = param.requires( "sig" );   
  
  signal_list_t signals = header.signal_list( signal_label );  
  
  const int ns = signals.size();
  

  //
  // requires at least one signal
  //

  if ( ns == 0 ) Helper::halt( "no signals specified for MATRIX" );
  

  //
  // Check FS for all signals
  //

  int fs = -1;
  
  for (int s=0; s<ns; s++) 
    {
      
      if ( header.is_data_channel( s ) )
	{
	  if ( fs < 0 ) fs = header.sampling_freq( signals(s) );
	  
	  else if ( header.sampling_freq( signals(s) ) != fs ) 
	    Helper::halt( "MATRIX requires uniform sampling rate across signals" ); 
	  
	  ++ns_data;
	}
    }
  
  
  
  //
  // Point to first epoch
  //
  
  if ( ! timeline.epoched() ) 
    {
      int n = timeline.set_epoch( globals::default_epoch_len , globals::default_epoch_len );
      logger << " set epochs to default " << globals::default_epoch_len << " seconds, " << n << " epochs\n";
    }

  timeline.first_epoch();

  const int ne = timeline.num_epochs();

  
  //
  // Output to log
  //
  
  std::cerr << " dumping " << ne << " unmasked epochs in " ;
  
  if ( minimal ) std::cerr << "minimal";
  else
    std::cerr << ( alternative_format ? "alternative" : "standard" ) ;
  
  std::cerr << " matrix-format to stdout\n";  


  

  if ( alternative_format )
    {
  
      //
      // Output
      //
      
      OUT << "I" << "\t" 
	  << id << "\t" 
	  << ne << "\t"
	  << timeline.epoch_length() << "\t"
	  << fs << "\t"
	  << na << "\t"
	  << ns_data << "\n";
      

      //
      // Epoch labels
      //
      
      OUT << "E";
      timeline.first_epoch();       
      while ( 1 ) 
	{
	  int epoch = timeline.next_epoch();	   
	  if ( epoch == -1 ) break;
	  OUT << "\t" << timeline.display_epoch( epoch ) ;
	}
      OUT << "\n";
      


      //
      // Output annotations (at the epoch level)
      //
      
      if ( show_annots )
	{
	  
	  // each row is a single annotation
	  
	  std::map<std::string,int>::const_iterator aa = atype.begin();
	  while ( aa != atype.end() )
	    {
	      
	      OUT << "A" << "\t" 
		  << aa->first ;


	      timeline.first_epoch();       
	      
	      while ( 1 ) 
		{
		  
		  //
		  // Get next epoch
		  //
		  
		  int epoch = timeline.next_epoch();
		  
		  if ( epoch == -1 ) break;

		  // no annotation: all 0

		  if ( aa->second == 0 ) 
		    {
		      OUT << "\t0";
		      continue;
		    }
		  
		  // interval annotation?

		  if ( aa->second == 1 ) 
		    {
		      
		      interval_t interval = timeline.epoch( epoch );
		      
		      annot_t * annot = timeline.annotations( aa->first );
		      
		      annot_map_t events = annot->extract( interval );
		      
		      bool has_annot = events.size() ;
		      
		      OUT << "\t" << ( has_annot ? 1 : 0 ) ;
		    }
		  
		  // epoch annotation
		  
		  if ( aa->second == 2 ) 
		    {
		      OUT << "\t" << ( timeline.epoch_annotation( aa->first , epoch ) ? 1 : 0 ) ;
		    }
		  

		} // next epoch
	      

	      OUT << "\n";

	      ++aa;
	    }
	 	  
	  
	}
      
     

      
      //
      // Output signals
      //
      
      
      
 
      //
      // For each signal
      //
      
      for (int s = 0 ; s < ns ; s++ )
	{	  
	  
	  // skip non-data channels
	  
	  if ( ! header.is_data_channel( s ) ) continue;
	  
	  OUT << "S" << "\t"	  
	      << header.label[ signals(s) ] ;
	  
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

	      for (int i=0;i<np;i++) OUT << "\t" << (*signal)[i];	  
	      
	    } // Next epoch           

	  OUT << "\n";
	  
	} // Next signal
      
      // all done

      
      OUT.close();

      return;
      
    }


  

  //
  // Standard matrix format
  //
  

  //
  // Header
  //

  bool include_hms = param.has( "hms" ) || param.has( "hms2" );
  bool include_hms2 = param.has( "hms2" );

  if ( ! minimal )
    {
      OUT << "ID\tE\tS\tSP\tT";
      
      if ( include_hms ) OUT << "\tHMS";
    }

  clocktime_t starttime( header.starttime );

  bool invalid_hms = ! starttime.valid;
  

  //
  // Annots header
  //
  
  if ( show_annots ) 
    {
      std::map<std::string,int>::const_iterator aa = atype.begin();
      while ( aa != atype.end() )
	{
	  OUT << "\t" << aa->first;
	  ++aa;
	}      
    }
  

  // Signals header
  
  if ( ! minimal )
    {
      for (int s = 0 ; s < ns ; s++ )
	{
	  if ( header.is_data_channel( s ) )
	    OUT << "\t" << header.label[ signals(s) ] ;
	}
      OUT << "\n";
    }


    
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

      std::vector<std::vector<double> > sigdat( ns_data );
      
      // track time-points, ie. may be a discontinuous file
      std::vector<uint64_t> tp;

      int s2 = 0;
      for (int s = 0 ; s < ns ; s++ )
	{	  	  
	  if ( header.is_data_channel( s ) )
	    {
	      slice_t slice( *this , signals(s) , interval );
	      const std::vector<double> * signal = slice.pdata();
	      sigdat[s2] = *signal;
	      if ( s2 == 0 ) tp = *slice.ptimepoints();
	      ++s2; // next data signal
	    }
	}
      
      // now iterate over all time-points (rows)
      
      const int np = sigdat[0].size();
      
      for (int t=0;t<np;t++)
	{
	  
	  double tp_sec = tp[t] / (double)globals::tp_1sec; 

	  double tp_sec_past_estart = ( tp[t] - interval.start) / (double)globals::tp_1sec; 
	  
	  // output rows
	  if ( ! minimal )
	    {
	      OUT << id << "\t"
		  << timeline.display_epoch( epoch ) << "\t"	    
		  << floor( tp_sec ) << "\t"
		  << t - fs * floor( tp_sec_past_estart )  << "\t"	  
		  << tp_sec;
	      
	      if ( include_hms )
		{
		  clocktime_t present = starttime;
		  
		  if ( include_hms2 ) 
		    {
		      interval_t now( tp[t] , tp[t]+1LLU );
		      std::string t1, t2;
		      if ( Helper::hhmmss( present , now , &t1,&t2 , 5 ) ) 
			OUT << "\t" << t1;
		      else
			OUT << "\t.";
		    }
		  else
		    {
		      present.advance( tp_sec / 3600.0 );
		      OUT << "\t" << present.as_string();
		    }
		}
	    }

	  // annots

	  if ( show_annots )
	    {
	      
	      std::map<std::string,int>::const_iterator aa = atype.begin();
	      while ( aa != atype.end() )
		{

		  if ( aa->second == 0 ) 
		    OUT << "\t0";
		  else if ( aa->second == 1 ) 
		    {		  
		      // get exact point
		      interval_t interval2 = interval_t( tp[t] , tp[t] + 1LLU );
		      annot_t * annot = timeline.annotations( aa->first );
		      annot_map_t events = annot->extract( interval2 );
		      bool has_annot = events.size() ;
		      OUT << "\t" << ( has_annot ? 1 : 0 ) ;
		    }
		  else if ( aa->second == 2 ) 
		    OUT << "\t" << ( timeline.epoch_annotation( aa->first , epoch ) ? 1 : 0 ) ;		  
		  
		  ++aa;
		}
	      
	    } 

	  
	  // signals
	  for (int s=0;s<ns_data;s++) OUT << ( ( ! minimal ) || s>0 ? "\t" : "" ) << sigdat[s][t];
	  
	  // done, next row/time-point
	  OUT << "\n";	  
	  
	} // next time-point
      
    } // next epoch
  
  
  OUT.close();
  
  return;

}

