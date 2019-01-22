
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


#include "pdc.h"

#include "helper/logger.h"
#include "db/db.h"
#include "eval.h"
#include "edf/edf.h"
#include "edf/slice.h"
#include "dsp/resample.h"

#include <string>
#include <iostream>
#include <cmath>
#include <set>



extern writer_t writer;

extern logger_t logger;

int pdc_t::m = 5;
int pdc_t::t = 1;
int pdc_t::q = 0;
bool pdc_t::store_ts = true;
std::vector<pdc_obs_t> pdc_t::obs;
std::set<std::string> pdc_t::labels;
std::map<std::string,int> pdc_t::label_count;
std::map<std::string,int> pdc_t::channels;  

///////////////////////////////////////////////////////////////////////////
// File format notes
//
//  one obs/channel per row
//  all channels must be grouped together, as separate lines 
//  but with the same ts-id (i.e. one pdc_obs_t)
//
//
// ts-library     
//  e.g. file1.tslib
//  ts-id indiv-id ch-id cat-id aux fs sp TS 
//
// pd-library (should match ts-library line for line)
//  PD stored as ints, divide by 'sum' to get vector of probabilities
//  e.g. file.pdlib
//  ts-id indiv-id ch-id cat-id aux fs m t sum PD 
//
///////////////////////////////////////////////////////////////////////////



void pdc_t::construct_tslib( edf_t & edf , param_t & param )
{
  
  // This implements luna command TSLIB:

  // TSLIB  ts-lib=lib/n1  sr=100  label=N1
  
  // ts-lib  sets root name of file
  // sr      sets sample rate
  // label   sets the 'label' for the records
  
  // will generate a series of files:  lib/n1-id001.tslib 
  //                                   lib/n1-id002.tslib 
  // etc

  // Note: this assumes that all relevant filtering and signal
  // selection has already been done via standard Luna scripting/masking
  //  i.e. above example is only for N1 sleep
  
  // per-individual TSLIB files can be easily concatenated across individuals
  // and a single Luna command (--pdlib) used to determine m/t and encode these
  // to create a matching PDLIB file

  
  std::string outfile = param.requires( "ts-lib" ) + "-" + edf.id + ".tslib" ;
  
  std::ofstream OUT1( outfile.c_str() , std::ios::out );
  
  //
  // Signals and sample-rate
  //
  
  std::string signal_label = param.requires( "signal" );   
  
  signal_list_t signals = edf.header.signal_list( signal_label );  

  const int ns = signals.size();

  // desired
  int sr = param.requires_int( "sr" );

  // actual
  std::vector<double> Fs = edf.header.sampling_freq( signals );

  // resampling?

  for (int s=0;s<ns;s++)
    {
      if ( edf.header.is_annotation_channel( signals(s) ) ) continue;
      
      if ( edf.header.sampling_freq( signals(s) ) != sr ) 
	{
	  logger << "resampling channel " << signals.label(s) 
		    << " from " << edf.header.sampling_freq( signals(s) )
		    << " to " << sr << "\n";
	  dsptools::resample_channel( edf, signals(s) , sr );
	}
    }  
  
  //
  // Category ID
  //
  
  std::string cat_id = param.requires( "cat" );
  

  //
  // Assumptions: 
  //   30 second epochs
  //   take only 10 seconds from the middle of each epoch when constructing a TS-LIB
  
  // note: for the test subject, each of the 3 10-sec intervals can be independently tested
  

  if ( edf.timeline.epoch_length() != 30 ) 
    Helper::halt( "TSLIB assumes 30-second epochs" );
  
  
  //
  // Iterate over each epoch
  //
  
  int ne = edf.timeline.first_epoch();
  
  int cnt = 0;
  
  while ( 1 ) 
    {
      
      int epoch = edf.timeline.next_epoch();      
      
      if ( epoch == -1 ) break;
      
      interval_t interval = edf.timeline.epoch( epoch );
      
      ++cnt;

      for ( int s=0; s<ns; s++ )
	{
	  
	  // only consider data tracks
	  
	  if ( edf.header.is_annotation_channel( signals(s) ) ) continue;
	  
	  slice_t slice( edf , signals(s) , interval );
	  
	  std::vector<double> * d = slice.nonconst_pdata();
	  
	  const std::vector<uint64_t> * tp = slice.ptimepoints();

	  // check epoch length is exactly 30s, otherwise skip
	  const int na = d->size();
	  if ( na != sr * 30 ) continue;
	  
	  // select middle 10 seconds
	  int start = 10 * sr - 1 ;
	  int end   = start + 10 * sr - 1 ;
	  int np    = end - start + 1 ; 
	  
	  // header
	  OUT1 << "e-" << epoch << "\t" 
	       << edf.id << "\t"
	       << signals.label(s) << "\t" 
	       << cat_id << "\t"
	       << "." << "\t"   // no auxillary information for now; format key=value;key=value
	       << sr << "\t"
	       << np ;
	  
	  // TS
	  for (int i = start ; i <= end ; i++ ) 
	    OUT1 << "\t" << (*d)[i] ;
	  
	  OUT1 << "\n";

	}
      
    }
  
  logger << " output " << cnt << " epochs for " << ns << " signals to TS-lib " << outfile << "\n";

  OUT1.close();

}



void pdc_t::entropy_heuristic_wrapper( param_t & param )
{

  // load TS-LIB

  std::string tslib = param.requires( "ts-lib" );

  read_tslib( tslib );

  // call entropy_heuristic( ) 
  
  bool by_cat = false;

  int m_min = 2; int m_max = 7;
  int t_min = 1; int t_max = 5;

  if ( param.has( "m" ) ) 
    {
      std::vector<int> p = param.intvector( "m" );
      if ( p.size() != 2 ) Helper::halt( "m=lwr,upr" );
      m_min = p[0]; m_max = p[1];
    }

  if ( param.has( "t" ) ) 
    {
      std::vector<int> p = param.intvector( "t" );
      if ( p.size() != 2 ) Helper::halt( "t=lwr,upr" );
      t_min = p[0]; t_max = p[1];
    }
  
  if ( param.has( "stratify" ) ) by_cat = true;
  
  entropy_heuristic( m_min, m_max, t_min, t_max, by_cat );
  
  
}


void pdc_t::construct_pdlib( param_t & param )
{

  //
  // initialize/clear obs/channel/label lists
  //

  clear();

  
  //
  // Always require that signals are specified 
  //

  if ( ! param.has( "ch" ) ) Helper::halt( "no ch={list} specified" );
  
  std::vector<std::string> c = param.strvector( "ch" ); 
  
  for (int j=0;j<c.size();j++) add_channel( c[j] );
  
  
  //
  // Function called directly from the luna command line, i.e. no EDFs to iterate through
  // Takes a TS-LIB file and creates a corresponding PD-LIB file, given set 'm' and 't' values
  //

  std::string infile  = param.requires( "ts-lib" );
  
  if ( param.has( "entropy" ) ) 
    {
      // automatically set m and t
      entropy_heuristic_wrapper( param );
    }
  else 
    {
      m = param.requires_int( "m" );
      t = param.requires_int( "t" );
    }

  //
  // Write a PDLIB?  If not, (i.e. only looking at min-e, then all done
  //

  if ( ! param.has( "pd-lib" ) ) return;

  //
  // Re-read TSLIB, encode using set values of m and t, and write a PDLIB
  //
  
  std::string outfile = param.requires( "pd-lib" );
  
  logger << "building " << outfile << " from " 
	    << infile << ", with m=" <<m << " and t=" << t << "\n";
    

  //
  // Simple read, line-for-line, except, only take channels specified in ch list
  //
  
  Helper::fileExists( infile );
  
  std::ifstream IN( infile.c_str() , std::ios::in );
  std::ofstream OUT( outfile.c_str() , std::ios::out );
  
  int cnt = 0; 

  while ( ! IN.eof() )
    {
      std::string ts_id, indiv_id, ch_id, cat_id, aux;
      int sr, sp;
      std::vector<double> x;
      IN >> ts_id;
      if ( IN.eof() ) break;
      IN >> indiv_id >> ch_id >> cat_id >> aux;
      IN >> sr >> sp;

      for (int i=0;i<sp;i++) 
	{
	  double t;
	  IN >> t;
	  x.push_back(t);	  
	}
     

      //
      // Keep this channel?
      //

      if ( ! has_channel( ch_id ) ) continue;
      
      
      //
      // Encode PD and write channel to PDLIB
      //
    
      int normalize = 0 ; // do not normalize

      std::vector<double> pd = pdc_t::calc_pd( x , m , t , &normalize );      

      // uniquify the ts_id here
      OUT << ts_id << "-" << indiv_id << "\t" 
	  << indiv_id << "\t"
	  << ch_id << "\t"
	  << cat_id << "\t"
	  << aux << "\t"
	  << sr << "\t"
	  << m << "\t"
	  << t << "\t"
	  << normalize ; 
      
      for (int j=0;j<pd.size();j++) OUT << "\t" << (int)pd[j] ;
      OUT << "\n";
      
      ++cnt;

    }
  
  IN.close();

  OUT.close();
  
  logger << " done.\n";
}



//
// Read a TS-LIB (for entropy heuristic)
//

void pdc_t::read_tslib( const std::string & tslib )
{
  
  if ( ! Helper::fileExists( tslib ) ) 
    Helper::halt( "could not find " + tslib );
  
  std::ifstream IN( tslib.c_str() , std::ios::in );
  
  logger << " reading ts-lib " << tslib << "\n";
  
  std::map<std::string,int> label_count;

  std::string previous_ts_id = ""; 
  
  //
  // observation construct with the correct number of slots always
  //
  
  pdc_obs_t ob(q);
  
  int cnt = 0 ;

  while ( ! IN.eof() ) 
    {
      
      // read next row
      std::string ts_id; 
      IN >> ts_id ;
      
      // last row?
      if ( IN.eof() ) 
	{
	  if ( previous_ts_id != "" ) 
	    {
	      label_count[ ob.label ]++;
	      add( ob );
	    }
	  break;
	}
      
      // line counter
      ++cnt;

      //  ts-id indiv-id ch-id cat-id aux fs sp TS 
    
      // read rest of line
      std::string indiv_id , ch_id , cat_id , aux;
      int sr, sp;
      std::vector<double> x;
      
      IN >> indiv_id >> ch_id >> cat_id >> aux;
      IN >> sr >> sp;
      
      for (int i=0;i<sp;i++) 
 	{
 	  double t;
 	  IN >> t;
 	  x.push_back(t);	  
 	}
      
      // a new interval? if so, save the old one
      
      if ( ts_id != previous_ts_id )
	{
	  // already read one observation?
	  if ( previous_ts_id != "" ) 
	    {
	      label_count[ ob.label ]++;
	      add( ob );
	    }      

	  // start this new observation of q potential channels
	  ob.init(q);
	  
	  ob.id = ts_id;
	  ob.label = cat_id;
	  // ignore aux encouding for now

	  // find channel and add to correct slot
	  int c = channel( ch_id );
	  
	  if ( c >= 0 ) 
	    {
	      ob.ch[c] = true;
	      ob.ts[c] = x;
	    }
	  
	  // track who we're recording
	  previous_ts_id = ts_id;
	}
      
      // otherwise, assume we are adding another channel to this observation
      // i.e. we can keep all other information as is
      
      // find channel and add to correct slot
      int c = channel( ch_id );
      
      if ( c >= 0 ) 
	{
	  ob.ch[c] = true;
	  ob.ts[c] = x;
	}
           
      // next row
    }
  
  IN.close();
 
  logger << " scanned " << cnt << " segments and read " << obs.size() << " observations\n";
  std::map<std::string,int>::const_iterator ii = label_count.begin();
  while ( ii != label_count.end() ) 
    {
      logger << "  " << ii->first << "\t" << ii->second << "\n";
      ++ii;
    }

  
  //
  // check/summarize channels
  //
  
  channel_check();

}





//
// Read a PD-LIB
//

void pdc_t::read_pdlib( const std::string & pdlib , const std::set<std::string> * incl_chs )
{
  
  if ( ! Helper::fileExists( pdlib ) ) 
    Helper::halt( "could not find " + pdlib );
  
  std::ifstream IN( pdlib.c_str() , std::ios::in );
  
  logger << " reading pd-lib " << pdlib << "\n";
  
  std::map<std::string,int> label_count;
  
  std::string previous_ts_id = ""; 
  
  pdc_obs_t ob(q);
  
  int cnt = 0;

  while ( ! IN.eof() ) 
    {
      
      // read next row
      std::string ts_id; 
      IN >> ts_id ;
      
      // last row?
      if ( IN.eof() ) 
	{
	  if ( previous_ts_id != "" ) 
	    {
	      label_count[ ob.label ]++;
	      add( ob );
	    }
	  break;
	}

      
      // read rest of line
      std::string indiv_id , ch_id , cat_id , aux;
      int sr, _m, _t; 
      double _sum;

      IN >> indiv_id >> ch_id >> cat_id >> aux;
      IN >> sr >> _m >> _t >> _sum;
      
      int nm = num_pd( _m );
      if ( nm == -1 ) Helper::halt("internal problem in pdc");
      std::vector<double> pd( nm );
      for (int j=0;j<nm;j++) 
	{
	  // DP will be as int, so normalize by sum upon reading
	  double pdval;
	  IN >> pdval;
	  pd[j] = pdval / _sum;
	}




      //
      // are we including this channel?
      //
      
      if ( incl_chs && incl_chs->find( ch_id ) == incl_chs->end() ) 
	continue;


      //
      // Note we're including this line
      //

      ++cnt; // line counter
     
      //
      // a new individual?  save the old one
      //

      if ( ts_id != previous_ts_id )
	{
	  // already read one observation?
	  if ( previous_ts_id != "" ) 
	    {
	      label_count[ ob.label ]++;
	      add( ob );
	    }      

	  // start this new observation
	  ob.init(q);
	  
	  ob.id = ts_id;
	  ob.label = cat_id;
	  // ignore aux encouding for now
	  
	  // find channel and add to correct slot
	  int c = channel( ch_id );
	  
	  if ( c >= 0 ) 
	    {
	      ob.ch[c] = true;
	      ob.pd[c] = pd;
	    }

	  // track who were are recording
	  previous_ts_id = ts_id;
	}
      
      // otherwise, assume we are adding another channel to this observation
      // i.e. we can keep all other information as is
      
      // find channel and add to correct slot
      int c = channel( ch_id );
      
      if ( c >= 0 ) 
	{
	  ob.ch[c] = true;
	  ob.pd[c] = pd;
	}
    
      // next row
    }
  
  IN.close();
  
  logger << " scanned " << cnt << " lines and read " << obs.size() << " observations\n";
  std::map<std::string,int>::const_iterator ii = label_count.begin();
  while ( ii != label_count.end() ) 
    {
      logger << "  " << ii->first << "\t" << ii->second << "\n";
      ++ii;
    }

  
  //
  // check/summarize channels (optional)
  //
  
  channel_check();

}



void pdc_t::test()
{

  clear();

  // read test.dat, 500x10  (but transpose)
  // times two, for multivariate comparison

  std::vector<std::vector<std::vector<double> > > data( 10 );

  for (int i=0;i<10;i++) 
    {
      data[i].resize( 2 );
      for (int j=0;j<2;j++) 
	data[i][j].resize( 500 , 0 );
    }
  
  for (int j=0;j<500;j++) 
    for (int i=0;i<10;i++)   
      {
	double d;
	std::cin >> d;
	data[i][0][j] = d;
      }

  for (int j=0;j<500;j++) 
    for (int i=0;i<10;i++)   
      {
	double d;
	std::cin >> d;
	data[i][1][j] = d;
      }
  
  
  add_channel( "CH1" );
  add_channel( "CH2" );
  
  for (int i=0;i<10;i++)
    {
      pdc_obs_t obs(q);

      obs.id = "obs-" + Helper::int2str(i+1);
      obs.label = "L1";
      
      int c = channel( "CH1" );
      obs.ch[c] = true;
      obs.ts[c] = data[i][0];

      c = channel( "CH2" );
      obs.ch[c] = true;
      obs.ts[c] = data[i][2];
      
      add( obs );
    }


  // init/check  (sets 'q')

  channel_check();
    

  // get entropy heuristics

  entropy_heuristic();
  
  // encode PD
  
  encode_ts();
    
  // distance matrix

  Data::Matrix<double> D = all_by_all();
  
  for (int i=0;i<10;i++)
    {
      for (int j=0;j<10;j++) logger << "\t" << D[i][j];
      logger << "\n";
    }
    
  return;  



}


Data::Matrix<double> pdc_t::all_by_all()
{


  const int N = obs.size();

  logger << " calculating " << N << "-by-" << N << " distance matrix\n";

  if ( N == 0 ) Helper::halt("internal error: PD not encoded in pdc_t");

  Data::Matrix<double> D( N , N );
//   for (int i=0;i<N;i++) 
//     D[i].resize(N,0);
  
  for (int i=0;i<(N-1);i++)
    for (int j=i+1;j<N;j++)
      D[i][j] = D[j][i] = distance( obs[i] , obs[j] );

  return D;
}




void pdc_t::encode_ts()
{  
  logger << " encoding with m="<<m << ", t=" << t << "\n";  
  // encode each observation in PD space
  const int N = obs.size();
  for (int i=0;i<N;i++) obs[i].encode( m , t );
}  

// encode TS(s) as PD(s)
void pdc_obs_t::encode( int m , int t )
{
  const int ns = ts.size();
  pd.resize(ns);
  int normalize = 1; 
  for (int i=0;i<ns;i++)
    pd[i] = pdc_t::calc_pd( ts[i] , m , t , &normalize );
}
  
// get (per-channel) entropy of current PD(s)
std::vector<double> pdc_obs_t::entropy() const
{
  const int ns = pd.size();
  std::vector<double> e( ns );
  for (int i=0;i<ns;i++)
    e[i] = pdc_t::entropy( pd[i] );
  return e;
}

void pdc_t::entropy_heuristic( int m_min , int m_max , int t_min , int t_max , bool by_cat )
{

  if ( m_min < 2 || m_max > 7 ) Helper::halt("invalid m ranges");
  if ( t_min < 1 || t_max > 5 ) Helper::halt("invalid t ranges");
  
  const int n = obs.size();
  if ( n == 0 ) Helper::halt( "no time series loaded" );

  double min_entropy = 1; 
  
  for (int mi = m_min ; mi <= m_max ; mi++ )
    {

      writer.level( mi, "PDC_M" );

      for (int ti = t_min ; ti <= t_max ; ti++ )
	{
	  writer.level( ti, "PDC_T" );
	  
	  std::vector<double> e;
	  
	  for (int i=0;i<n;i++)
	    {
	      // encode TS and PD with mi/ti 
	      obs[i].encode( mi, ti );
	      // and then get resulting entropy values per channel
	      std::vector<double> e1 = obs[i].entropy( );
	      for (int j=0;j<e1.size();j++) e.push_back( e1[j] );
	    }
	  
	  // nb. mean over all observations and channels
	  double mean_entropy = MiscMath::mean( e );
	  
	  writer.value( "E" , mean_entropy );
	  
	  if ( mean_entropy < min_entropy ) 
	    {
	      set_param( mi, ti );
	      min_entropy = mean_entropy;
	    }
	}
      writer.unlevel( "PDC_T" );
    }
  writer.unlevel( "PDC_M" );
  
  writer.value( "PDC_OPT_M" , m );
  writer.value( "PDC_OPT_T" , t );
  
  logger << " based on min entropy, setting m = " << m << ", t = " << t << "\n";
  
  //
  // As above, but stratify by label
  //  
  
  if ( by_cat && labels.size() > 1 ) 
    {
      
      logger << " additionally, stratifying by " << labels.size() << " distinct labels\n";
      
      std::set<std::string>::const_iterator ii = labels.begin();
      while ( ii != labels.end() )
	{
	  
	  writer.level( *ii, "PDC_LABEL" );
	  
	  double min_entropy = 1; 
	  int best_m = m_min; int best_t = t_min;
	  
	  for (int mi = m_min ; mi <= m_max ; mi++ )
	    {
	      
	      writer.level( mi, "PDC_M" );
	      
	      for (int ti = t_min ; ti <= t_max ; ti++ )
		{
		  
		  writer.level( ti, "PDC_T" );
		  
		  std::vector<double> e;
		  
		  for (int i=0;i<n;i++)
		    {
		      // only consider observations with matching labels
		      if ( obs[i].label != *ii ) continue;
		      
		      // encode TS and PD with mi/ti 
		      obs[i].encode( mi, ti );
		      
		      // and then get resulting entropy values per channel
		      std::vector<double> e1 = obs[i].entropy( );
		      for (int j=0;j<e1.size();j++) e.push_back( e1[j] );
		    }
		  
		  // nb. mean over all observations and channels
		  double mean_entropy = MiscMath::mean( e );
		  
		  writer.value( "E" , mean_entropy );
		  
		  if ( mean_entropy < min_entropy ) 
		    {
		      best_m = mi;
		      best_t = ti;
		      min_entropy = mean_entropy;
		    }
		  
		}

	      writer.unlevel( "PDC_T" );
	      
	    }

	  writer.unlevel( "PDC_M" );
	  
	  writer.value( "PDC_OPT_M" , best_m );
	  writer.value( "PDC_OPT_T" , best_t );	  
	  
	  ++ii;
	}
      
      writer.unlevel( "PDC_LABEL" );
      
    }
    
}


void pdc_t::add( const pdc_obs_t & ob ) {
  obs.push_back( ob );
  labels.insert( ob.label );
  label_count[ ob.label ]++;
  // channels are preset so do not record here
}

double pdc_t::distance( const pdc_obs_t & a, const pdc_obs_t & b )
{

  if ( q == 0 ) 
    return 0;
  
  if ( a.pd[0].size() != b.pd[0].size() ) 
    Helper::halt( "incompatible PD -- check similar m used" );

  
  // univariate
  if ( q == 1 ) 
    return symmetricAlphaDivergence( a.pd[0] , b.pd[0] );
  
  // in multichannel case, define total obs-obs distance as sqrt( sum(d^2) ) 
  double d = 0;
  for (int k=0;k<q;k++) 
    d += MiscMath::sqr( symmetricAlphaDivergence( a.pd[k] , b.pd[k] ) );    
  d = sqrt(d);
  return d;

}


double pdc_t::distance( const pdc_obs_t & a, const pdc_obs_t & b , const std::vector<int> & chs )
{
 
  if ( q == 0 || chs.size() == 0 ) 
    return 0;
  
  if ( a.pd[0].size() != b.pd[0].size() ) 
    Helper::halt( "incompatible PD -- check similar m used" );
  
  // univariate (and includes this channel)
  if ( q == 1 && chs[0] == 0 ) 
    return symmetricAlphaDivergence( a.pd[0] , b.pd[0] );
  
  // in multichannel case, define total obs-obs distance as sqrt( sum(d^2) ) 
  double d = 0;
  for (int k=0;k<chs.size();k++)
    {
      if ( chs[k] >= q ) return 0;
      d += MiscMath::sqr( symmetricAlphaDivergence( a.pd[chs[k]] , b.pd[chs[k]] ) );    
    }
  d = sqrt(d);
  return d;
}



void pdc_t::channel_check()
{

  const int N = obs.size();
  if ( N == 0 ) return;
  
  std::map<std::string,int> chs;
  
  for ( int i = 0 ; i < N ; i++ ) 
    {
      
      std::string c;
      std::map<std::string,int>::const_iterator cc = channels.begin();
      while ( cc != channels.end() )
	{
	  if ( obs[i].ch[ cc->second ] ) 
	    {
	      if ( c == "" ) c = cc->first;
	      else c += "; " + cc->first;
	    }
	  ++cc;
	}

      // record
      chs[ c ]++;

    }

  logger << " of " << N << " observations, following breakdown by available channels:\n";
  std::map<std::string,int>::const_iterator cc = chs.begin();
  while ( cc != chs.end() )
    {
      logger << " " << cc->second << "\t" << cc->first << "\n";
      ++cc;
    }
    
}


