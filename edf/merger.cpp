
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


#include "helper/helper.h"
#include "db/db.h"

#include "defs/defs.h"
#include "edf/edf.h"
#include "edf/slice.h"

#include <iostream>

extern writer_t writer;

bool identical_headers( const edf_header_t & h1 , const edf_header_t & h2 ) ;

void Helper::merge_EDFs( const std::vector<std::string> & tok )
{

  // ignore times:: assume that first file is the start, then just append all subsequent in order

  const bool use_fixed_order = false;
  
  std::vector<edf_t*> edfs;

  std::string id = "id";
  std::string filename = "merged.edf";
  std::string slist = "";

  // expecting a list of file names 
  // but can also include key=value pairs, with keys:
  //   id
  //   edf
  //   sample-list
  
  for (int i=0; i<tok.size(); i++)
    {

      logger << "processing [" << tok[i] << "]\n";
      
      std::vector<std::string> tok2 = Helper::quoted_parse( tok[i] , "=" );
      if ( tok2.size() == 2 )
	{
	  if ( tok2[0] == "id" ) id = tok2[1];
	  else if ( tok2[0] == "edf" ) filename = tok2[1];
	  else if ( tok2[0] == "sample-list" ) slist = tok2[1];
	  continue;
	}
      
      const std::string fname = Helper::expand( tok[i] );
      
      if ( ! Helper::fileExists( fname ) )
	{
	  logger << "  ** warning: could not attach " << fname << "\n";
	  continue;
	}
      
      edf_t * edf = new edf_t;
      
      const std::string id = "id" + Helper::int2str( (int)( edfs.size() + 1 ) ) ;
      
      
      bool okay = edf->attach( fname , id );
      
      if ( ! okay )
	{
	  logger << " ** could not attach " << filename	<< "\n";
          continue;	  
	}

      logger << "  attached component EDF " << fname << "\n";

      edfs.push_back( edf ) ;

    }

  const int nf = edfs.size();

   
  //
  // Some outputs
  //

  logger << "\n  attached " << nf << " EDFs\n";  

  for (int i=0; i<nf; i++)
    logger << "  EDF " << i+1 << "\t"
	   << edfs[i]->id << "\t"
	   << edfs[i]->header.startdate << "\t"
	   << edfs[i]->header.starttime << "\n";

  logger  << "   writing merged data\n"
	  << "     ID           : " << id << "\n"
	  << "     EDF filename : " << filename << "\n";

  if ( slist != "" )
   logger << "     sample-list  : " << slist << "\n"; 
  
  
  //
  // Check that all headers are compatible
  //   -- initially, headers must be identical
  //   -- this can be relaxed in the future.
  //

  for (int i=1; i<nf; i++)
    if ( ! identical_headers( edfs[0]->header , edfs[i]->header ) )
      Helper::halt( "headers incompatible:" + edfs[0]->filename + " " + edfs[i]->filename );
  
  logger  << "  good, all EDFs have merge-compatible headers\n";
  

  //
  // Get total implied NR for new EDF
  //
  
  int nr = 0;

  for (int i=0; i<nf; i++)
    nr += edfs[i]->header.nr;

  logger  << "  expecting " << nr
	  << " records (each of "
	  << edfs[0]->header.record_duration << " sec) in the new EDF\n";

  //
  // check that all segments are contiguous
  //   -- NOT IMPLEMENTED YET
  //   -- initially, we will not allow overlaps or gaps
  //   -- downstream, can allow gaps and overlaps (i.e. if different channels for same time)
  //

  // <---- TODO ----->

  
  //
  // get earliest start time;
  //
  
  std::string first_date = edfs[0]->header.startdate;
  for (int i=1; i<nf; i++)
    if ( edfs[i]->header.startdate < first_date )
      first_date = edfs[i]->header.startdate;
  
  
  clocktime_t t1( edfs[0]->header.starttime );
  if ( !t1.valid )
    Helper::halt( edfs[0]->filename + " does not have a valid start time" );

  
  for (int i=1; i<nf; i++)
    {
      // only check same day 
      if ( first_date == edfs[i]->header.startdate )
	{
	  clocktime_t t2( edfs[i]->header.starttime );
	  if ( !t2.valid )
	    Helper::halt( edfs[i]->filename + " does not have a valid start time" );
	  
	  if ( clocktime_t::earlier( t1 , t2 ) == 2 )
	    t1 = t2;
	}
    }

  logger << "  first record starts at " << first_date << "  " << t1.as_string() << "\n";
  
  //
  // Create the new merged EDF
  //
  
  edf_t medf;
  
  //
  // Set header
  //

  medf.header.version = edfs[0]->header.version;
  medf.header.patient_id = edfs[0]->header.patient_id;
  medf.header.recording_info = edfs[0]->header.recording_info;
  medf.header.startdate = first_date;
  medf.header.starttime = t1.as_string();
  medf.header.nbytes_header = edfs[0]->header.nbytes_header;
  medf.header.ns = 0; // this is popuilated by edf_t::add_signal()
  medf.header.ns_all = 0; // this is popuilated by edf_t::add_signal()
  medf.header.nr = medf.header.nr_all = nr; // this is sum across all EDFs computed above
  medf.header.record_duration = edfs[0]->header.record_duration;
  medf.header.record_duration_tp = edfs[0]->header.record_duration_tp;

  
  //
  // create a (continuous) timeline  
  //

  logger << " adding timeline\n";

  medf.set_edf();
  medf.set_continuous();
  medf.timeline.init_timeline();


  //
  // resize data[][], by adding empty records 
  //

  logger << " adding " << nr << " empty records...\n";

  for (int r=0;r<nr;r++)
    {
      edf_record_t record( &medf ); 
      medf.records.insert( std::map<int,edf_record_t>::value_type( r , record ) );
    }
  

  //
  // add signals 
  //

  const int ns = edfs[0]->header.ns ;
  
  for (int s=0; s<ns; s++)
    {

      // skip annotations
      if ( edfs[0]->header.is_annotation_channel( s ) ) continue;
      
      logger << "  compiling channel " << edfs[0]->header.label[s] << "\n";

      std::vector<double> dt;

      for (int j=0; j<nf; j++)
	{

	  edf_t & edf = *(edfs[j]);
	  
	  // get whole signal
	  slice_t slice( edf , s , edf.timeline.wholetrace() );

	  std::vector<double> * d = slice.nonconst_pdata();

	  for (int k=0; k<d->size(); k++) dt.push_back( (*d)[k] ) ;
	  
	} // next EDF

      const int np_obs = dt.size();
      const int np_exp = nr * edfs[0]->header.sampling_freq(s);
      
      if ( np_obs != np_exp )
	Helper::halt( "expected and observed number of sample points did not align" );
      
      // add the signal to the merged EDF
      medf.add_signal( edfs[0]->header.label[ s ] , edfs[0]->header.sampling_freq( s ) , dt );
            
    } // next channel
  

  //
  // Save this merged EDF
  //

  logger << "  writing merged EDF as " << filename << "\n";
  
  bool saved = medf.write( filename );
  
  if ( ! saved ) Helper::halt( "problem trying to write " + filename );
  
  //
  // Sample list?
  //

  if ( slist != "" )
    {	  
      logger << " appending " << filename << " to sample-list " << slist << "\n";
      std::ofstream FL( slist.c_str() , std::ios_base::app );
      FL << medf.id << "\t" << filename << "\n";
      FL.close();
    }
  
  //
  // Clean up
  //
  
  for (int i=0; i<nf; i++)
    delete edfs[i];      
  
}


bool identical_headers( const edf_header_t & h1 , const edf_header_t & h2 )
{

  // compatible for merging?  currently needs to have *identical*
  // headers (for signal #, SR and EDF record size)

  if ( h1.version != h2.version ) return false;
  
  if ( h1.ns != h2.ns ) return false;
  
  if ( h1.record_duration_tp != h2.record_duration_tp ) return false;

  for ( int s = 0 ; s < h1.ns; s++ )
    {
      if ( h1.label[s] != h2.label[s] ) return false;
      if ( h1.n_samples[s] != h2.n_samples[s] ) return false;
    }
  
  return true;
  
}
