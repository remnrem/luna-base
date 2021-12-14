
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

#include "suds.h"

#include <vector>
#include <map>
#include <set>
#include <iomanip>

#include "helper/helper.h"
#include "helper/logger.h"
#include "db/db.h"

#include "dirent.h"

#include "stats/eigen_ops.h"
#include "stats/lda.h"
#include "stats/statistics.h"
#include "miscmath/miscmath.h"
#include "miscmath/crandom.h"

#include "edf/edf.h"
#include "edf/slice.h"

#include "dsp/resample.h"
#include "fftw/fftwrap.h"
#include "dsp/mtm/mtm.h"
#include "dsp/tv.h"

extern logger_t logger;

extern writer_t writer;

//
// Utilities to construct a SUDS trainer library
//  - libraries can be either text or binary; with utility to copy between them
//  - libraries designed to be concatenated together (either binary or text)
//    when reading back, i.e. so just need a single file to be distributed
//

// Format
//   - version number (SUDS3)
//   - trainer ID
//   - nf  = number of features (expected to match the corresponding model file)
//   - nc  = number of components
//   - nve = number of (valid) epochs
//   - observed stages [ nve ]
//   - feature (X) means (over epochs) [ nf ]
//   - feature (X) SDs (over epochs) [ nf ]
//   - U [ nve x nc ]
//   - D [ nc ]
//   - V [ nc x nc ]
//   - X [ nve x nf ]  :: original features, nb. this only needs to be read for weight-trainers

// Information on the signals (SR, lwr/frq, features, etc) is now in the model-file
// Optionally, there can also be a weights-file that gives weights for each feature


void suds_indiv_t::write( edf_t & edf , param_t & param ) const
{

  // write as a single file in the folder specifed by 'db'

  const std::string folder = Helper::expand( param.requires( "db" ) );
  
  const int ns = suds_t::model.chs.size();
  
  // create output folder if it does not exist
  std::string syscmd = globals::mkdir_command + " " + folder ;
  int retval = system( syscmd.c_str() );

  // for saving trainers: use EDF ID, or a fake ID?  (e.g. 'ids=suds')
  std::string suds_id = suds_t::fake_ids ? suds_t::fake_id_root + "_" + Helper::int2str( suds_t::fake_ids++ ) : edf.id;
  
  std::string filename = folder + globals::folder_delimiter + suds_id;

  logger << "  writing trainer data to " << filename << "\n";
  
  std::ofstream OUT1( filename.c_str() , std::ios::out );

  // file version code == 3
  OUT1 << "SUDS3\n";

  // ID
  OUT1 << suds_id << "\n";

  // NVE, NS, NF, NC
  OUT1 << "% number of 1) valid epochs, 2) signals, 3) features, 4) SVD components\n";
  OUT1 << nve << "\n"
       << ns << "\n"
       << nf << "\n"
       << nc << "\n";
    
  // STAGE COUNTS
  OUT1 << "% number of stages, then # of epochs per stage\n";
  OUT1 << counts.size() << "\n";
  std::map<std::string,int>::const_iterator ss = counts.begin();
  while ( ss != counts.end() )
    {
      OUT1 << ss->first << "\n";
      OUT1 << ss->second << "\n";
      ++ss;
    }
  
  // Stages by epoch-by-epoch
  OUT1 << "% Epoch-wise stage assignments ( epoch # --> stage )\n";
  for (int i=0;i<nve;i++)
    OUT1 << epochs[i] << "\n"
	 << y[i] << "\n";

  // Feature summary statistics
  OUT1 << "% Hjorth parameter summary stats (mean, SD) for H1, H2, H3\n";
  for (int s=0; s<ns; s++)
    OUT1 << mean_h1[s] << "\n"
	 << sd_h1[s] << "\n"
	 << mean_h2[s] << "\n"
	 << sd_h2[s] << "\n"
	 << mean_h3[s] << "\n"
	 << sd_h3[s] << "\n";      

  // PSC components
  OUT1 << "% SVD W diagonals\n"; 
  for (int j=0;j<nc;j++)
    OUT1 << W[j] << "\n";
  
  // V
  OUT1 << "% SVD V matrix (" << nf << " by " << nc << ")\n";
  for (int i=0;i<nf;i++)
    for (int j=0;j<nc;j++)
      OUT1 << V(i,j) << "\n";
  
  // U (to re-estimate LDA model upon loading, i.e
  //  to use lda.predict() on target
  //  only needs to be nc rather than nf
  OUT1 << "% SVD U matric (" << nve << " by " << nc << ")\n";
  for (int i=0;i<nve;i++)
    for (int j=0;j<nc;j++)
      OUT1 << U(i,j) << "\n";

  
  // X, raw feature data, as defined by the model-file
  //  --> if this trainer, this is being used as a 'weight trainer' (re-prediction)
  //  --> i.e. will project this individual's raw data into the target space

  OUT1 << "% Feature matrix X (" << nve << " by " << nf << ")\n";
  for (int i=0;i<nve;i++)
    for (int j=0;j<nf;j++)
      OUT1 << X(i,j) << "\n";
  
  //
  // All done
  //
  
  OUT1.close();
  
}

bool next( std::ifstream & IN1 , std::string * line )
{
  while ( 1 ) 
    {
      Helper::safe_getline( IN1 , *line );
      if ( IN1.eof() ) return false; 
      if ( *line == "" ) continue;
      if ( (*line)[0] == '%' ) continue;
      break;
    }
  return true;
}

void suds_t::text2binary( const std::string & texfile ,
			  const std::string & binfile , 
			  const bool with_features )
  
{
  
  // convert format from text to binary 
  // read text in: note, may be concatenated
  
  if ( ! Helper::fileExists( Helper::expand( texfile ) ) )
    Helper::halt( "could not open " + Helper::expand( texfile ) );
  
  // read text from here...
  std::ifstream IN1( Helper::expand( texfile ).c_str() , std::ios::in );
  
  // and write out binary to here.
  std::ofstream OUT1( Helper::expand( binfile ).c_str() , std::ios::binary | std::ios::out );
  
  logger << "  copying from " << texfile << " to " << binfile << " (text2binary conversion)\n";
  if ( with_features ) logger << "  including feature matrices in final output\n";
  else logger << "  not including feature matrices in final output\n";
  
  
  //
  // Write out
  //

  int p = -1;
  int n_indiv = 0;
  uint64_t ecnt = 0;
  
  while ( 1 )
    {

      int i;
      double d;
      std::string line;
      
      // SUDX code (w/ 'f' suffix for features)
      if ( ! next(IN1 , &line ) ) break;
      
      suds_indiv_t::bwrite( OUT1 , with_features ? line + "f" : line );
      
      // ID
      next(IN1,&line);
      suds_indiv_t::bwrite( OUT1 , line );

      // NVE
      int tnve = 0 , tns = 0, tnf = 0 , tnc = 0;
      next(IN1,&line);
      if ( ! Helper::str2int( line , &i ) )
	Helper::halt( "bad numeric" );
      tnve = i;
      suds_indiv_t::bwrite( OUT1 , i );      
      ecnt += tnve;
      
      // NS
      next(IN1,&line);
      if ( ! Helper::str2int( line , &i ) )
	Helper::halt( "bad numeric" );
      tns = i;
      suds_indiv_t::bwrite( OUT1 , i );      

      // NF
      next(IN1,&line);
      if ( ! Helper::str2int( line , &i ) )
	Helper::halt( "bad numeric" );
      tnf = i;
      suds_indiv_t::bwrite( OUT1 , i );      

      // NC
      next(IN1,&line);
      if ( ! Helper::str2int( line , &i ) )
	Helper::halt( "bad numeric" );
      tnc = i;
      suds_indiv_t::bwrite( OUT1 , i );
      
      // Stage counts
      next(IN1,&line);      
      if ( ! Helper::str2int( line , &i ) )
        Helper::halt( "bad numeric" );      
      int tstages = i;
      suds_indiv_t::bwrite( OUT1 , i );
      
      // Each stage count
      for (int j=0;j<tstages; j++)
	{
	  // stage label
	  next(IN1,&line);
	  suds_indiv_t::bwrite( OUT1 , line );	  

	  // stage count
	  next(IN1,&line);
	  if ( ! Helper::str2int( line , &i ) )
	    Helper::halt( "bad numeric(2)" );
	  suds_indiv_t::bwrite( OUT1 , i );
	  
	}
      
      // Stages epoch-by-epoch
      for (int j=0;j<tnve; j++)
	{
	  // epoch number
          next(IN1,&line);
	  if ( ! Helper::str2int( line , &i ) )
            Helper::halt( "bad numeric(3)" );
          suds_indiv_t::bwrite( OUT1 , i );
          
	  // stage label
          next(IN1,&line);
	  suds_indiv_t::bwrite( OUT1 , line );
          
	}
      
      // Hjorth summary statistics
      for (int j=0;j<tns; j++)
	{
	  for (int h=0; h<3; h++ )
	    {
	      // feature mean (over epochs)
	      next(IN1,&line);
	      if ( ! Helper::str2dbl( line , &d ) )
		Helper::halt( "bad numeric(4)" );
	      suds_indiv_t::bwrite( OUT1 , d );          
	      
	      // feature SD (over epochs)
	      next(IN1,&line);
	      if ( ! Helper::str2dbl( line , &d ) )
		Helper::halt( "bad numeric(5)" );
	      suds_indiv_t::bwrite( OUT1 , d );          
	    }
	}
      
      
      // SVD components: W
      for (int j=0;j<tnc; j++)
	{
          next(IN1,&line);
	  if ( ! Helper::str2dbl( line , &d ) )
            Helper::halt( "bad numeric(6)" );
          suds_indiv_t::bwrite( OUT1 , d );          
	}
      
      // SVD components: V 
      for (int j=0;j<tnf; j++)
        for (int k=0;k<tnc; k++)
	  {
	    next(IN1,&line);
	    if ( ! Helper::str2dbl( line , &d ) )
	      Helper::halt( "bad numeric(7)" );
	    suds_indiv_t::bwrite( OUT1 , d );	    
	  }
      
      // SVD components: U
      for (int i=0;i<tnve;i++)
	for (int j=0;j<tnc;j++)
	  {
	    next(IN1,&line);
	    if ( ! Helper::str2dbl( line , &d ) )
              Helper::halt( "bad numeric(8)" );
            suds_indiv_t::bwrite( OUT1 , d );            
	  }
      
      // Original features: X  (optional)
      for (int i=0;i<tnve;i++)
        for (int j=0;j<tnf;j++)
	  {
	    next(IN1,&line);
	    
	    if ( with_features ) // i.e. might skip in output
	      {
		if ( ! Helper::str2dbl( line , &d ) )
		  Helper::halt( "bad numeric(9)" );
		suds_indiv_t::bwrite( OUT1 , d );            
	      }
	  }

      ++n_indiv;

      logger << "  " << n_indiv << " trainers compiled...\n";
      
    }

  //
  // mark EOF explicitly (i.e. as 'SUDSX' version number for next indiv.
  //
  
  suds_indiv_t::bwrite( OUT1 , "_END_" );
  

  //
  // All done
  //

  OUT1.close();
  IN1.close();

  logger << "  in total, converted " << n_indiv << " trainers (" << ecnt << " epochs)\n";
  
}


void suds_indiv_t::bwrite( std::ofstream & O , const std::string & s ) 
{
  uint8_t l = s.size();
  O.write( (char*)( &l ), sizeof(uint8_t) );
  O.write( s.c_str(), l );
}

void suds_indiv_t::bwrite( std::ofstream & O , int i ) 
{
  O.write( (char*)( &i ), sizeof(int) );
}

void suds_indiv_t::bwrite( std::ofstream & O , double d ) 
{
  O.write( (char*)( &d ), sizeof(double) );
}

std::string suds_indiv_t::bread_str( std::ifstream & I )
{
  uint8_t len;
  I.read( (char*)( &len ), sizeof(uint8_t) );
  std::vector<char> b( len );
  I.read( &b[0] , len );
  std::string s( b.begin() , b.end() );
  return s;
}

int suds_indiv_t::bread_int( std::ifstream & I )
{
  int i;
  I.read( (char*)( &i ), sizeof(int) );
  return i;
}

double suds_indiv_t::bread_dbl( std::ifstream & I )
{
  double d;
  I.read( (char*)( &d ), sizeof(double) );
  return d;
}

void suds_indiv_t::bskip_dbl( std::ifstream & I , const int n )
{
  std::vector<double> dummy( n ) ;
  I.read( (char*)( &dummy[0] ), n * sizeof(double) );
}


std::vector<suds_indiv_t*> suds_t::binary_reload( const std::string & filename , bool load_rawx )
{
  
  if ( ! Helper::fileExists( Helper::expand( filename ) ) )
    Helper::halt( "cannot open " + filename );
  
  std::ifstream IN1( Helper::expand( filename ).c_str() , std::ios::binary | std::ios::in );

  std::vector<suds_indiv_t*> bank;

  // assume this might contain multiple individuals

  while ( 1 )
    {

      // SUDSX magic number
      std::string suds = suds_indiv_t::bread_str( IN1 );

      // all done?
      if ( suds == "_END_" ) break;

      // otherwise, check format 
      if ( suds != "SUDS3" && suds != "SUDS3f" )
	Helper::halt( "bad file format for " + filename );
      
      if ( suds == "SUDS3f" && ! load_rawx ) Helper::halt( "library has features, load as 'wdb' " );
      if ( suds == "SUDS3"  &&   load_rawx ) Helper::halt( "library does not have features, load as 'db' " );

      // adding a new individual
      suds_indiv_t * person = new suds_indiv_t;

      // track progress
      if ( bank.size() % 50 == 0 ) logger << "\n ";
      if ( bank.size() % 10 == 0 ) logger << " ";
      logger << ".";

      // ID
      person->id = suds_indiv_t::bread_str( IN1 );

      // NVE
      person->nve = suds_indiv_t::bread_int( IN1 );

      // NS
      int ns0 = suds_indiv_t::bread_int( IN1 );
      if ( ns0 != suds_t::ns )
	Helper::halt( "different specification of 'ns' " );
      
      // NF
      person->nf = suds_indiv_t::bread_int( IN1 );
      
      // NC (which may be lower than higher bound)
      int this_nc = suds_indiv_t::bread_int( IN1 );
  
      if ( this_nc == 0 )
	Helper::halt( "0 PSCs for " + filename );
      
      // set 'nc' for this individual
      person->nc = this_nc;
  
      // stage summaries
      int nstages = suds_indiv_t::bread_int( IN1 );
      
      for (int i=0; i<nstages; i++)
	{
	  const std::string sname = suds_indiv_t::bread_str( IN1 );
	  const int scnt = suds_indiv_t::bread_int( IN1 );
	  person->counts[ sname ] = scnt;
	}

      // stages
      person->y.resize( person->nve );
      person->epochs.resize( person->nve );
      for (int i=0;i<person->nve;i++)
	{
	  person->epochs[i] = suds_indiv_t::bread_int( IN1 );
	  person->y[i] = suds_indiv_t::bread_str( IN1 );
	}
      person->obs_stage = suds_t::type( person->y );

      // Hjorth summary stats (mean/SD)
      person->mean_h1.resize( suds_t::ns );
      person->sd_h1.resize( suds_t::ns ) ;
      person->mean_h2.resize( suds_t::ns );
      person->sd_h2.resize( suds_t::ns ) ;
      person->mean_h3.resize( suds_t::ns );
      person->sd_h3.resize( suds_t::ns ) ;
    
      for (int s=0;s<suds_t::ns;s++)
	{
	  person->mean_h1[s] = suds_indiv_t::bread_dbl( IN1 );
	  person->sd_h1[s] = suds_indiv_t::bread_dbl( IN1 );
	  person->mean_h2[s] = suds_indiv_t::bread_dbl( IN1 );
	  person->sd_h2[s] = suds_indiv_t::bread_dbl( IN1 );
	  person->mean_h3[s] = suds_indiv_t::bread_dbl( IN1 );
	  person->sd_h3[s] = suds_indiv_t::bread_dbl( IN1 );
	}
      
      // SVD: W [ only nc ]
      person->W.resize( person->nc );
      for (int j=0;j<person->nc;j++)
	person->W[j] = suds_indiv_t::bread_dbl( IN1 );
      
      // V [ only nc cols ] 
      person->V.resize( person->nf, person->nc );
      for (int i=0;i<person->nf;i++)
	for (int j=0;j<person->nc;j++)
	  person->V(i,j) = suds_indiv_t::bread_dbl( IN1 );
      
      // U (to reestimate LDA model upon loading, i.e
      //  to use lda.predict() on target 
      person->U.resize( person->nve , person->nc );
      for (int i=0;i<person->nve;i++)
	for (int j=0;j<person->nc;j++)
	  person->U(i,j) = suds_indiv_t::bread_dbl( IN1 );
      
      // X original feature matrix [ nve x nf ]
      // i.e. if this trainer is being used as a 'weight trainer',
      // i.e. will project this individuals raw data into the target space
      
      if ( load_rawx )
	{            
	  person->X.resize( person->nve , person->nf );
	  for (int i=0;i<person->nve;i++)
	    for (int j=0;j<person->nf;j++)
	      person->X(i,j) = suds_indiv_t::bread_dbl( IN1 );
	}

      // else
      // 	{
      // 	  // we need to skip these elements
      // 	  suds_indiv_t::bskip_dbl( IN1 , person->nve * person->nf );
      // 	  std::cout << " nah\n";
      // 	}

      //
      // add this person
      //
      
      bank.push_back( person );

      //
      // get the next person
      //
    }

  IN1.close();
  
  //
  // All done
  //

  return bank;
}


//
// attach a BINARY singlefile library (plus/minus feature data)
//

void suds_t::attach_db( const std::string & file0 , bool read_db , bool read_wdb )
{

  // already populated?
  if ( read_db && bank.size() > 0 ) return;
  if ( read_wdb && wbank.size() > 0 ) return;
  
  if ( ! ( read_db || read_wdb ) )
    Helper::halt( "bad call to suds_t::attach_db()" );
  
  const std::string filename = Helper::expand( file0 );
  
  if ( ! Helper::fileExists( filename ) )
    Helper::halt( "cannot open " + filename );
  
  // are we populating the 'bank', 'wbank' or both?
  std::map<std::string,suds_indiv_t*> * b1 = read_db  ? &bank  : NULL ;
  std::map<std::string,suds_indiv_t*> * b2 = read_wdb ? &wbank : NULL ;
  
  logger << "  attaching training data from " << filename ;
  
  //
  // Read all data from binary file 
  //

  std::vector<suds_indiv_t*> trainers = binary_reload( filename , read_wdb );
  
  const int nt = trainers.size();
  
  //
  // for primary trainers only (! read_psd ) track Hjorth distributions
  //
  
  Eigen::MatrixXd h1_mean( nt , suds_t::ns );
  Eigen::MatrixXd h1_sd( nt , suds_t::ns );

  Eigen::MatrixXd h2_mean( nt , suds_t::ns );
  Eigen::MatrixXd h2_sd( nt , suds_t::ns );

  Eigen::MatrixXd h3_mean( nt , suds_t::ns );
  Eigen::MatrixXd h3_sd( nt , suds_t::ns );

  //
  // process each
  //
  
  for ( int i=0; i<nt ; i++)
    {
      
      suds_indiv_t * trainer = trainers[i];	  

      trainer->fit_lda();
      
      //
      // store in the relevant bank(s):
      //

      if ( b1 ) (*b1)[ trainer->id ] = trainer;

      if ( b2 ) (*b2)[ trainer->id ] = trainer;

      
      //
      // for primary trainers only, copy over feature summary stats
      //
      
      if ( read_db ) 
	{
	  for (int s=0; s<suds_t::ns; s++)
	    {
	      h1_mean(i,s) = trainer->mean_h1[s] ;
	      h1_sd(i,s) = trainer->sd_h1[s] ;
	      
	      h2_mean(i,s) = trainer->mean_h2[s] ;
	      h2_sd(i,s) = trainer->sd_h2[s] ;
	      
	      h3_mean(i,s) = trainer->mean_h3[s] ;
	      h3_sd(i,s) = trainer->sd_h3[s] ;
	    }
	}
      
    } 
  
  logger << "\n"
	 << "  attached " << trainers.size() << " trainers\n";


  
  //
  // from primary trainers only, track feature-wise 95% CI limits
  //
  

  if ( read_db )
    {
      suds_t::hjorth1_lwr95.resize( suds_t::ns );
      suds_t::hjorth1_upr95.resize( suds_t::ns );      

      suds_t::hjorth2_lwr95.resize( suds_t::ns );
      suds_t::hjorth2_upr95.resize( suds_t::ns );      
      
      suds_t::hjorth3_lwr95.resize( suds_t::ns );
      suds_t::hjorth3_upr95.resize( suds_t::ns );      

      Eigen::ArrayXd mean_h1mean = h1_mean.colwise().mean();
      Eigen::ArrayXd mean_h1sd   = h1_sd.colwise().mean();

      Eigen::ArrayXd mean_h2mean = h2_mean.colwise().mean();
      Eigen::ArrayXd mean_h2sd   = h2_sd.colwise().mean();
      
      Eigen::ArrayXd mean_h3mean = h3_mean.colwise().mean();
      Eigen::ArrayXd mean_h3sd   = h3_sd.colwise().mean();
      
      for (int s=0; s<suds_t::ns; s++)
	{
	  suds_t::hjorth1_lwr95[s] = mean_h1mean[s] - suds_t::hjorth_outlier_th * mean_h1sd[s];
	  suds_t::hjorth1_upr95[s] = mean_h1mean[s] + suds_t::hjorth_outlier_th * mean_h1sd[s];
	  
	  suds_t::hjorth2_lwr95[s] = mean_h2mean[s] - suds_t::hjorth_outlier_th * mean_h2sd[s];
	  suds_t::hjorth2_upr95[s] = mean_h2mean[s] + suds_t::hjorth_outlier_th * mean_h2sd[s];
	  
	  suds_t::hjorth3_lwr95[s] = mean_h3mean[s] - suds_t::hjorth_outlier_th * mean_h3sd[s];
	  suds_t::hjorth3_upr95[s] = mean_h3mean[s] + suds_t::hjorth_outlier_th * mean_h3sd[s];
	}
    }

}

