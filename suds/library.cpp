
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

#include "param.h"
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
#include "spectral/mtm/mtm.h"
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
//   - version number (SUDS1)
//   - trainer ID
//   -  hasX, hasLDA, hasQDA?
//   - nf  = number of features (expected to match the corresponding model file)
//   - nc  = number of components
//   - nve = number of (valid) epochs
//   - observed stages [ nve ]
//   - feature (X) means (over epochs) [ nf ]
//   - feature (X) SDs (over epochs) [ nf ]
//   - D [ nc ]
//   - V [ nc x nc ]
//   - optional: LDA model
//   - optional: QDA model
//   - optional: X [ nve x nf ]  :: original features, nb. this only needs to be read for weight-trainers

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

  const bool output_features = param.yesno( "output-X" );
  const bool output_lda = param.yesno( "output-LDA" ) && lda_model.valid ;
  const bool output_qda = param.yesno( "output-QDA" ) && qda_model.valid ;
  
  std::ofstream OUT1( filename.c_str() , std::ios::out );
  
  // file version code
  OUT1 << suds_t::suds_lib_version << "\n";
  
  // ID
  OUT1 << suds_id << "\n";

  // contents
  OUT1 << "X:" << ( output_features ? "Y" : "N" ) << "\n";
  OUT1 << "LDA:" << ( output_lda ? "Y" : "N" ) << "\n";
  OUT1 << "QDA:" << ( output_qda ? "Y" : "N" ) << "\n";
  
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
    OUT1 //<< epochs[i] << "\n"  --> no need for trainer epoch numbers
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
  // OUT1 << "% SVD U matric (" << nve << " by " << nc << ")\n";
  // for (int i=0;i<nve;i++)
  //   for (int j=0;j<nc;j++)
  //     OUT1 << U(i,j) << "\n";

  
  //
  // LDA model
  //

  if ( output_lda )
    {
      OUT1 << "%LDA model\n";
      
      // number of groups
      OUT1 << lda_model.prior.size() << "\n";
      
      // number of predictors/variables
      OUT1 << lda_model.means.cols() << "\n";
      
      // priors
      for (int i=0;i<lda_model.prior.size();i++)
	OUT1 << lda_model.prior[i] << "\n";
            
      // counts
      std::map<std::string,int>::const_iterator ii = lda_model.counts.begin();
      while ( ii != lda_model.counts.end() ) 
	{
	  OUT1 << ii->first << "\n"
	       << ii->second << "\n";
	  ++ii;
	}

      // means      
      for (int i=0; i<lda_model.means.rows(); i++)
	for (int j=0; j<lda_model.means.cols(); j++)
	  OUT1 << lda_model.means(i,j) << "\n";
      
      // scaling
      OUT1 << lda_model.scaling.rows() << "\n"
	   << lda_model.scaling.cols() << "\n";
      for (int j=0; j < lda_model.scaling.rows(); j++)
	for (int k=0; k < lda_model.scaling.cols(); k++)
	  OUT1 << lda_model.scaling(j,k) << "\n";
	      
      // n
      OUT1 << lda_model.n << "\n";

      // labels
      for (int i=0;i<lda_model.labels.size();i++)
	OUT1 << lda_model.labels[i] << "\n";
      
    }
  
  
  //
  // QDA model
  //
  
  if ( output_qda )
    {
      
      OUT1 << "%QDA model\n";

      // number of groups
      OUT1 << qda_model.prior.size() << "\n";
      
      // number of predictors/variables
      OUT1 << qda_model.means.cols() << "\n";
      
      // priors
      for (int i=0;i<qda_model.prior.size();i++)
	OUT1 << qda_model.prior[i] << "\n";
      
      // rows (redundant, but keep)
      for (int i=0;i<qda_model.rows.size();i++)
	OUT1 << qda_model.rows[i] << "\n";
      
      // counts
      std::map<std::string,int>::const_iterator ii = qda_model.counts.begin();
      while ( ii != qda_model.counts.end() ) 
	{
	  OUT1 << ii->first << "\n"
	       << ii->second << "\n";
	  ++ii;
	}

      // means
      for (int i=0; i<qda_model.means.rows(); i++)
	for (int j=0; j<qda_model.means.cols(); j++)
	  OUT1 << qda_model.means(i,j) << "\n";

      // scaling
      for (int i=0; i<qda_model.scaling.size(); i++)
	{
	  const Eigen::MatrixXd & m = qda_model.scaling[i];
	  for (int j=0; j < m.rows(); j++)
	    for (int k=0; k < m.cols(); k++)
	      OUT1 << m(j,k) << "\n";
	}
      
      // ldet
      for (int i=0;i<qda_model.ldet.size();i++)
	OUT1 << qda_model.ldet[i] << "\n";

      // n
      OUT1 << qda_model.n << "\n";

      // labels
      for (int i=0;i<qda_model.labels.size();i++)
	OUT1 << qda_model.labels[i] << "\n";
      
    }
  
  // 
  // X, raw feature data, as defined by the model-file
  //  --> if this trainer, this is being used as a 'weight trainer' (re-prediction)
  //  --> i.e. will project this individual's raw data into the target space

  if ( output_features )
    {
      OUT1 << "% Feature matrix X (" << nve << " by " << nf << ")\n";
      for (int i=0;i<nve;i++)
	for (int j=0;j<nf;j++)
	  OUT1 << X(i,j) << "\n";
    }
  
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
      //std::cout << " line [" << *line << "]\n";
      if ( *line == "" ) continue;
      if ( (*line)[0] == '%' ) continue;
      
      break;
    }
  return true;
}

void suds_t::text2binary( const std::string & texfile ,
			  const std::string & binfile , 
			  const bool with_features )  // ignored for now
  
{
  
  // convert format from text to binary 
  // read text in: note, may be concatenated

  // old version could drop features from the text version (with_features = F) 
  // this is ignored for now ... i..e. create originals w/ or w/out features as desired
  
  if ( ! Helper::fileExists( Helper::expand( texfile ) ) )
    Helper::halt( "could not open " + Helper::expand( texfile ) );
  
  // read text from here...
  std::ifstream IN1( Helper::expand( texfile ).c_str() , std::ios::in );
  
  // and write out binary to here.
  std::ofstream OUT1( Helper::expand( binfile ).c_str() , std::ios::binary | std::ios::out );
  
  logger << "  copying from " << texfile << " to " << binfile << " (text2binary conversion)\n";
    
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
      
      // SUDX code 
      if ( ! next(IN1 , &line ) ) break;      
      suds_indiv_t::bwrite( OUT1 , line );

      // ID
      next(IN1,&line);
      suds_indiv_t::bwrite( OUT1 , line );

      // Contents:
      // X:?
      next(IN1,&line);
      const bool has_features = line == "X:Y";
      suds_indiv_t::bwrite( OUT1 , line );
      
      // LDA:?
      next(IN1,&line);
      const bool has_lda = line == "LDA:Y";
      suds_indiv_t::bwrite( OUT1 , line );

      // QDA:?
      next(IN1,&line);
      const bool has_qda = line == "QDA:Y";
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

	  // epoch number  [ NOT USED ANY MORE ]
          // next(IN1,&line);
	  // if ( ! Helper::str2int( line , &i ) )
          //   Helper::halt( "bad numeric(3)" );
          // suds_indiv_t::bwrite( OUT1 , i );
          
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
      // for (int i=0;i<tnve;i++)
      // 	for (int j=0;j<tnc;j++)
      // 	  {
      // 	    next(IN1,&line);
      // 	    if ( ! Helper::str2dbl( line , &d ) )
      //         Helper::halt( "bad numeric(8)" );
      //       suds_indiv_t::bwrite( OUT1 , d );            
      // 	  }
      

      if ( has_lda )
	{
	  // number of groups
	  next(IN1,&line);
	  if ( ! Helper::str2int( line , &i ) )
	    Helper::halt( "bad numeric" );
	  int ng = i;
	  suds_indiv_t::bwrite( OUT1 , i );
	  
	  // number of variables
	  next(IN1,&line);
	  if ( ! Helper::str2int( line , &i ) )
	    Helper::halt( "bad numeric" );
	  int nv = i;
	  suds_indiv_t::bwrite( OUT1 , i );
	  
	  // priors
	  for (int k=0; k<ng; k++)
	    {
	      next(IN1,&line);
	      if ( ! Helper::str2dbl( line , &d ) )
		Helper::halt( "bad numeric(6)" );
	      suds_indiv_t::bwrite( OUT1 , d );
	    }
	  
	  // counts (str --> int )
	  for (int k=0;k<ng;k++)
	    {
	      // group label (str)
	      next(IN1,&line);
	      suds_indiv_t::bwrite( OUT1 , line );
	      
	      // count (int)
	      next(IN1,&line);
	      if ( ! Helper::str2int( line , &i ) )
		Helper::halt( "bad numeric(3)" );
	      suds_indiv_t::bwrite( OUT1 , i );	      
	    }
	  
	  // means
	  for (int k=0;k<ng;k++)
	    for (int m=0; m<nv; m++)
	      {		
		next(IN1,&line);
		if ( ! Helper::str2dbl( line , &d ) )
		  Helper::halt( "bad numeric(3)" );
		suds_indiv_t::bwrite( OUT1 , d );
	      }
	  
	  // scaling (get row/col size explicitly)
	  next(IN1,&line);
	  if ( ! Helper::str2int( line , &i ) )
	    Helper::halt( "bad numeric(3)" );
	  if ( i != nv ) Helper::halt( "format problem" );
	  suds_indiv_t::bwrite( OUT1 , i );

	  next(IN1,&line);
	  if ( ! Helper::str2int( line , &i ) )
	    Helper::halt( "bad numeric(3)" );
	  const int nr = i;	  
	  suds_indiv_t::bwrite( OUT1 , i );
	  
	  for (int k=0; k<nv; k++)
	    for (int m=0; m<nr; m++)
	      {		
		next(IN1,&line);
		if ( ! Helper::str2dbl( line , &d ) )
		  Helper::halt( "bad numeric(3)" );
		suds_indiv_t::bwrite( OUT1 , d );
	      }
	  
	  // n
	  next(IN1,&line);
	  if ( ! Helper::str2int( line , &i ) )
	    Helper::halt( "bad numeric(3)" );
	  suds_indiv_t::bwrite( OUT1 , i );

	  // labels
	  for (int k=0;k<ng;k++)
	    {
	      next(IN1,&line);
	      suds_indiv_t::bwrite( OUT1 , line );	      
	    }
	  
	}
      
      
      if ( has_qda )
	{
	  
	  // number of groups
	  next(IN1,&line);
	  if ( ! Helper::str2int( line , &i ) )
	    Helper::halt( "bad numeric" );
	  int ng = i;
	  suds_indiv_t::bwrite( OUT1 , i );
	  
	  // number of variables
	  next(IN1,&line);
	  if ( ! Helper::str2int( line , &i ) )
	    Helper::halt( "bad numeric" );
	  int nv = i;
	  suds_indiv_t::bwrite( OUT1 , i );
	  
	  // priors
	  for (int k=0; k<ng; k++)
	    {
	      next(IN1,&line);
	      if ( ! Helper::str2dbl( line , &d ) )
		Helper::halt( "bad numeric(6)" );
	      suds_indiv_t::bwrite( OUT1 , d );
	    }
	  
	  // rows
	  for (int k=0; k<ng; k++)
	    {
	      next(IN1,&line);
	      if ( ! Helper::str2int( line , &i ) )
		Helper::halt( "bad numeric(6)" );
	      suds_indiv_t::bwrite( OUT1 , i );
	    }

	  // counts (str --> int )
	  for (int k=0;k<ng;k++)
	    {
	      // group label (str)
	      next(IN1,&line);
	      suds_indiv_t::bwrite( OUT1 , line );
	      
	      // count (int)
	      next(IN1,&line);
	      if ( ! Helper::str2int( line , &i ) )
		Helper::halt( "bad numeric(3)" );
	      suds_indiv_t::bwrite( OUT1 , i );	      
	    }
	  
	  // means
	  for (int k=0;k<ng;k++)
	    for (int m=0; m<nv; m++)
	      {		
		next(IN1,&line);
		if ( ! Helper::str2dbl( line , &d ) )
		  Helper::halt( "bad numeric(3)" );
		suds_indiv_t::bwrite( OUT1 , d );
	      }
	  
	  // scaling (get row/col size explicitly)
	  for (int g=0; g<ng; g++)
	    for (int k=0; k<nv; k++)
	      for (int m=0; m<nv; m++)
	      {
		next(IN1,&line);
		if ( ! Helper::str2dbl( line , &d ) )
		  Helper::halt( "bad numeric(3)" );
		suds_indiv_t::bwrite( OUT1 , d );
	      }
	  
	  // ldet
	  for (int k=0; k<ng; k++)
	    {
	      next(IN1,&line);
	      if ( ! Helper::str2dbl( line , &d ) )
		Helper::halt( "bad numeric(6)" );
	      suds_indiv_t::bwrite( OUT1 , d );
	    }
	  
	  // n
	  next(IN1,&line);
	  if ( ! Helper::str2int( line , &i ) )
	    Helper::halt( "bad numeric(3)" );
	  suds_indiv_t::bwrite( OUT1 , i );
	  
	  // labels
	  for (int k=0;k<ng;k++)
	    {
	      next(IN1,&line);
	      suds_indiv_t::bwrite( OUT1 , line );	      
	    }

	}

      
      //
      // Original features: X  (optional)
      //

      if ( has_features )
	{
	  for (int j=0;j<tnve;j++)
	    for (int k=0;k<tnf;k++)
	      {
		next(IN1,&line);
		
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
      const std::string suds = suds_indiv_t::bread_str( IN1 );
      
      // all done?
      if ( suds == "_END_" ) break;
      
      // otherwise, check format version
      if ( suds != suds_t::suds_lib_version ) 
	Helper::halt( "bad file format for " + filename
		      + ", expecting " + suds_t::suds_lib_version
		      + " but found " + suds );

      // track progress
      if ( bank.size() % 50 == 0 ) logger << "\n ";
      if ( bank.size() % 10 == 0 ) logger << " ";
      logger << ".";
      
      // adding a new individual
      suds_indiv_t * person = new suds_indiv_t;

      // ID
      person->id = suds_indiv_t::bread_str( IN1 );

      // get contents::
      //    - features (X) included Y/N
      //    - LDA model included Y/N
      //    - QDA model included Y/N

      const bool has_features = suds_indiv_t::bread_str( IN1 ) == "X:Y";
      const bool has_lda = suds_indiv_t::bread_str( IN1 ) == "LDA:Y";
      const bool has_qda = suds_indiv_t::bread_str( IN1 ) == "QDA:Y";
    
      if ( has_features    && ! load_rawx ) Helper::halt( "library has features, load as 'wdb' " );
      if ( (!has_features) &&   load_rawx ) Helper::halt( "library does not have features, load as 'db' " );
      
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
	  
	  // note: do not read epoch numbers any more, we don't need
	  //  these in non-targets (and so they are not stored)
	  //person->epochs[i] = suds_indiv_t::bread_int( IN1 );
	  person->epochs[i] = i+1; 

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

      // not needed now
      // U (to reestimate LDA model upon loading, i.e
      //  to use lda.predict() on target 
      // person->U.resize( person->nve , person->nc );
      // for (int i=0;i<person->nve;i++)
      // 	for (int j=0;j<person->nc;j++)
      // 	  person->U(i,j) = suds_indiv_t::bread_dbl( IN1 );


      //
      // LDA model?
      //

      if ( has_lda )
	{
	  // number of groups
	  const int ng = suds_indiv_t::bread_int( IN1 );

	  // number of variables
	  const int nv = suds_indiv_t::bread_int( IN1 );
	  
	  // priors
	  person->lda_model.prior.resize( ng );
	  for (int i=0;i<ng;i++)
	    person->lda_model.prior[i] = suds_indiv_t::bread_dbl( IN1 );
	  
	  // counts
	  for (int i=0;i<ng;i++)
	    {
	      const std::string s = suds_indiv_t::bread_str( IN1 );
	      person->lda_model.counts[ s ] = suds_indiv_t::bread_int( IN1 );
	    }
	  
	  // means
	  person->lda_model.means.resize( ng , nv );
	  for (int i=0; i<ng; i++)
	    for (int j=0; j<nv; j++)
	      person->lda_model.means(i,j) = suds_indiv_t::bread_dbl( IN1 );


	  // scaling
	  int s1 = suds_indiv_t::bread_int( IN1 );
	  int s2 = suds_indiv_t::bread_int( IN1 );
	  person->lda_model.scaling.resize( s1 , s2 );
	  for (int j=0; j < s1; j++)
	    for (int k=0; k < s2; k++)
	      person->lda_model.scaling(j,k) = suds_indiv_t::bread_dbl( IN1 );
	  	  
	  // n
	  person->lda_model.n = suds_indiv_t::bread_int( IN1 );
	  
	  // labels
	  person->lda_model.labels.resize( ng );
	  for (int i=0;i<person->lda_model.labels.size();i++)
	    person->lda_model.labels[i] =  suds_indiv_t::bread_str( IN1 );
	  
	}
	  
      //
      // QDA model?
      //
      
      if ( has_qda )
	{

	  // number of groups
	  const int ng = suds_indiv_t::bread_int( IN1 );

	  // number of variables
	  const int nv = suds_indiv_t::bread_int( IN1 );
	  
	  // priors
	  person->qda_model.prior.resize( ng );
	  for (int i=0;i<ng;i++)
	    person->qda_model.prior[i] = suds_indiv_t::bread_dbl( IN1 );
	  
	  // rows (redundant, but keep)
	  person->qda_model.rows.resize( ng );
	  for (int i=0;i<ng;i++)
	    person->qda_model.rows[i] = suds_indiv_t::bread_int( IN1 );
	  
	  // counts
	  for (int i=0;i<ng;i++)
	    {
	      const std::string s = suds_indiv_t::bread_str( IN1 );
	      person->qda_model.counts[ s ] = suds_indiv_t::bread_int( IN1 );
	    }
	  
	  // means
	  person->qda_model.means.resize( ng , nv );
	  for (int i=0; i<ng; i++)
	    for (int j=0; j<nv; j++)
	      person->qda_model.means(i,j) = suds_indiv_t::bread_dbl( IN1 );

	  // scaling
	  person->qda_model.scaling.resize( ng );
	  for (int i=0; i<ng; i++)
	    {
	      person->qda_model.scaling[i].resize( nv , nv );
	      Eigen::MatrixXd & m = person->qda_model.scaling[i];
	      for (int j=0; j < m.rows(); j++)
		for (int k=0; k < m.cols(); k++)
		  m(j,k) = suds_indiv_t::bread_dbl( IN1 );
	    }
	  
	  // ldet
	  person->qda_model.ldet.resize( ng );
	  for (int i=0;i<ng;i++)
	    person->qda_model.ldet[i] = suds_indiv_t::bread_dbl( IN1 );
	  
	  // n
	  person->qda_model.n = suds_indiv_t::bread_int( IN1 );

	  // labels
	  person->qda_model.labels.resize( ng );
	  for (int i=0;i<person->qda_model.labels.size();i++)
	    person->qda_model.labels[i] =  suds_indiv_t::bread_str( IN1 );
	      
	  
	}
      
      
      // Features (X) ? 
      //  - feature matrix [ nve x nf ]
      //    only need this is trainer is to be used as a 'weight trainer'
      //     (i.e. will project this individual's raw data into the target SVD space)
      
      if ( has_features )
	{ 
	  person->X.resize( person->nve , person->nf );
	  for (int i=0;i<person->nve;i++)
	    for (int j=0;j<person->nf;j++)
	      person->X(i,j) = suds_indiv_t::bread_dbl( IN1 );
	}
      
      //
      // add this person 
      //
      
      if ( suds_t::single_trainer == "" )
	bank.push_back( person );
      else
	{
	  if ( person->id == suds_t::single_trainer ) 
	    bank.push_back( person );
	  else
	    delete person;
	}

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
      
      // LDA/QDA models now precomputed and stored in the library
      // no U matrix, so no need/ability to recompte, so comment 
      // out this line:

      //trainer->fit_qlda();
      

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



void suds_t::attach_lib( const std::string & infile ) 
{
  // already populated?
  if ( bank.size() != 0 ) return;
  
  // look for infile.fit, infile.svd and infile.hjorth
  // (can extend this to have multiple .fit and .svd pairs too)
 
  logger << "  attaching pre-fit trainer library " << infile << "\n";

  attach_db_prefit( infile );
  
  attach_hjorth_limits( infile + ".hjorth" );
  
  logger << "  bank size = " << bank.size() << "\n";

}


//
// attach a single prefit trainer model (i.e. no indiv. level data)
//

void suds_t::attach_db_prefit( const std::string & infile )
{
  
  // reads infile.fit
  //       infile.svd
  //  also, function below will read static values of infile.hjorth 
  //  

  //
  // LDA/QDA model (i.e. based on 1 or more real trainers)
  //
  
  suds_indiv_t * trainer = new suds_indiv_t;
  
  trainer->qda_model.read( infile + ".fit" ) ; 
  
  bank[ trainer->id ] = trainer;
  
  //
  // V and W matrices from SVD (for target projection)
  //
  std::string svdfile = Helper::expand( infile + ".svd" ) ;
  if ( ! Helper::fileExists( svdfile ) )
    Helper::halt( "could not find " + svdfile );
  
  std::ifstream I1( svdfile.c_str() , std::ios::in );

  int vw;
  I1 >> vw;
  trainer->W.resize( vw );
  for (int i=0; i<vw; i++) 
    I1 >> trainer->W[i];
  
  trainer->nc = trainer->W.size();

  int vr , vc;
  I1 >> vr >> vc;
  trainer->V.resize( vr , vc );
  for (int i=0; i<vr; i++) 
    for (int j=0; j<vc; j++) 
      I1 >> trainer->V(i,j);
  
  I1.close();
  
}


void suds_t::attach_hjorth_limits( const std::string & hjorthfile )
{

  // this file is created only by --combine-suds

  // i.e.  to read back in pre-fit data, we will call attach_db_prefit() above [ with the QDA/LDA model ]
  //       and then this function to set the lower/upper 95% CI limits

  if ( ! Helper::fileExists( Helper::expand( hjorthfile ) ) )
    Helper::halt( "could not open " + hjorthfile );

  std::ifstream I1( Helper::expand( hjorthfile ).c_str() , std::ios::in ) ;
  
  int ns0;
  I1 >> ns0;
  if ( suds_t::ns != ns0 ) 
    {
      logger << "  expecting " << ns << " signals, but " << hjorthfile << " has " << ns0 << "\n";
      Helper::halt( "bad hjorthfile" );
    }
  
  suds_t::hjorth1_lwr95.resize( suds_t::ns );
  suds_t::hjorth1_upr95.resize( suds_t::ns );      
  
  suds_t::hjorth2_lwr95.resize( suds_t::ns );
  suds_t::hjorth2_upr95.resize( suds_t::ns );      
  
  suds_t::hjorth3_lwr95.resize( suds_t::ns );
  suds_t::hjorth3_upr95.resize( suds_t::ns );      
      
  for (int s=0; s<suds_t::ns; s++)
    {
      double h1_m, h2_m, h3_m;
      double h1_s, h2_s, h3_s;
      I1 >> h1_m >> h1_s 
	 >> h2_m >> h2_s 
	 >> h3_m >> h3_s;

      suds_t::hjorth1_lwr95[s] = h1_m - h1_s * suds_t::hjorth_outlier_th;
      suds_t::hjorth1_upr95[s] = h1_m + h1_s * suds_t::hjorth_outlier_th;

      suds_t::hjorth2_lwr95[s] = h2_m - h2_s * suds_t::hjorth_outlier_th;
      suds_t::hjorth2_upr95[s] = h2_m + h2_s * suds_t::hjorth_outlier_th;

      suds_t::hjorth3_lwr95[s] = h3_m - h3_s * suds_t::hjorth_outlier_th;
      suds_t::hjorth3_upr95[s] = h3_m + h3_s * suds_t::hjorth_outlier_th;

    }
   
  I1.close();
  
}



void suds_t::combine_trainers( param_t & param )
{

  logger << "  combining multiple trainer feature sets...\n";

  suds_t::set_options( param );
  
  // we must have NC explicitly set here (i.e. as we are not reading a model specification at this point)
  suds_t::nc = param.requires_int( "nc" );
  if ( nc < 1 || nc > 50 ) Helper::halt( "bad nc value" );

  std::string infile = param.requires( "from" );
  std::string outfile  = param.requires( "to" );
  
  // convert format from text to binary 
  // read text in: note, may be concatenated
  
  if ( ! Helper::fileExists( Helper::expand( infile ) ) )
    Helper::halt( "could not open " + Helper::expand( infile ) );
  
  // read binary file here (w/ multiple indivs)
  std::ifstream IN1( Helper::expand( infile ).c_str() , std::ios::binary | std::ios::in );


  // only take first n individuals

  const int n_start = param.has( "first" ) ? param.requires_int( "first" ) : -1 ; 
  const int n_stop = param.has( "last" ) ? param.requires_int( "last" ) : 0 ; 



  // start
  int p = -1;
  int n_indiv = 0;
  int added = 0; 

  // epoch-level
  std::vector<std::string> stages;
  
  // indiv-summaries
  std::vector<std::vector<double> > h1_means, h2_means, h3_means;
  std::vector<std::vector<double> > h1_vars, h2_vars, h3_vars;
  
  // check that number of features matches 
  int first_nf = 0;
  int first_ns = 0;

  //
  // A new mega-indiv to be created
  //  
  
  suds_indiv_t mega;
  mega.id = param.requires( "id" );
  mega.trainer = true;
  mega.nve = 0;
  mega.nf = 0;
  mega.nc = suds_t::nc;
  mega.X = Eigen::MatrixXd::Zero( 0 , 0 );
  
  // staging
  std::vector<suds_stage_t> obs_stage; 
  
  //
  // Iterate over files
  //
  
  while ( 1 )
    {

      // all done?
      if ( n_stop != 0 && n_indiv == n_stop ) break;
      
      // add this person, or skip? 
      const bool add_person = n_start == -1 || n_start <= n_indiv ; 

      // SUDSX magic number
      const std::string suds = suds_indiv_t::bread_str( IN1 );
      
      // all done?
      if ( suds == "_END_" ) break;
      
      // otherwise, check format version
      if ( suds != suds_t::suds_lib_version ) 
	Helper::halt( "bad file format for " + infile
		      + ", expecting " + suds_t::suds_lib_version
		      + " but found " + suds );
      
      // ID
      std::string id = suds_indiv_t::bread_str( IN1 );      
      std::cout << n_indiv + 1 << "\t" << id << "\t" << ( add_person ? "added" : "skipped" ) << "\n";
      
      // get contents::
      //    - features (X) included Y/N
      //    - LDA model included Y/N
      //    - QDA model included Y/N

      const bool has_features = suds_indiv_t::bread_str( IN1 ) == "X:Y";
      const bool has_lda = suds_indiv_t::bread_str( IN1 ) == "LDA:Y";
      const bool has_qda = suds_indiv_t::bread_str( IN1 ) == "QDA:Y";
      
      if ( ! has_features ) 
	Helper::halt( "file " + id + " does not contain raw features:: cannot compile into a single trainer\nrun MAKE-SUDS with output-X=T" );
            
      int tnve = suds_indiv_t::bread_int( IN1 );
      int tns = suds_indiv_t::bread_int( IN1 );
      int tnf = suds_indiv_t::bread_int( IN1 );
      int tnc = suds_indiv_t::bread_int( IN1 );
      
      // on first individual, set space
      if ( added == 0 ) 
	{
	  first_ns = tns;
	  h1_means.resize( first_ns );
	  h2_means.resize( first_ns );	  
	  h3_means.resize( first_ns );
	  h1_vars.resize( first_ns );
	  h2_vars.resize( first_ns );	  
	  h3_vars.resize( first_ns );
	}
      else if ( first_ns != tns ) 
	Helper::halt( "all inputs must have same # of signals" );
      

      // Stage counts
      int tstages = suds_indiv_t::bread_int( IN1 );
      
      for (int i=0; i<tstages; i++)
	{
	  const std::string sname = suds_indiv_t::bread_str( IN1 );
	  const int scnt = suds_indiv_t::bread_int( IN1 );	  
	}

      // Stages epoch-by-epoch
      for (int j=0;j<tnve; j++)
	{                    
	  // add stages to mega indiv.
	  mega.obs_stage.push_back( suds_t::type( suds_indiv_t::bread_str( IN1 ) ) );
	}
      
      // Hjorth summary statistics
      for (int j=0;j<tns; j++)
	{
	  for (int h=0; h<3; h++ )
	    {
	      double d = suds_indiv_t::bread_dbl( IN1 ); 
	      
	      if ( h == 0 ) h1_means[j].push_back( d );
	      else if ( h == 1 ) h2_means[j].push_back( d );
	      else h3_means[j].push_back( d );

	      // feature SD (over epochs)
	      d = suds_indiv_t::bread_dbl( IN1 ); 	      
	      if ( h == 0 ) h1_vars[j].push_back( d*d );
	      else if ( h == 1 ) h2_vars[j].push_back( d*d );
	      else h3_vars[j].push_back( d*d );
	    }
	}
      
      // skip SVD components (these will be recalculated) 
      // SVD components: W
      for (int j=0;j<tnc; j++)
	{
	  double d =  suds_indiv_t::bread_dbl( IN1 );
	}
      
      // SVD components: V 
      for (int j=0;j<tnf; j++)
        for (int k=0;k<tnc; k++)
	  {
	    double d =  suds_indiv_t::bread_dbl( IN1 );
	  }      
      
      // SVD components: U
      // for (int i=0;i<tnve;i++)
      // 	for (int j=0;j<tnc;j++)
      // 	  next(IN1,&line);
      
      //
      // LDA : ignore 
      //
      
      if ( has_lda ) 
	{

	  // number of groups
	  const int ng = suds_indiv_t::bread_int( IN1 );
	  
	  // number of variables
	  const int nv = suds_indiv_t::bread_int( IN1 );
	  
	  // priors
	  for (int i=0;i<ng;i++)
	    {
	      double d = suds_indiv_t::bread_dbl( IN1 );
	    }	  

	  // counts
	  for (int i=0;i<ng;i++)
	    {
	      const std::string s = suds_indiv_t::bread_str( IN1 );
	      int j = suds_indiv_t::bread_int( IN1 );
	    }
	  
	  // means
	  for (int i=0; i<ng; i++)
	    for (int j=0; j<nv; j++)
	      {
		double d = suds_indiv_t::bread_dbl( IN1 );
	      }

	  // scaling
	  int s1 = suds_indiv_t::bread_int( IN1 );
	  int s2 = suds_indiv_t::bread_int( IN1 );
	  for (int j=0; j < s1; j++)
	    for (int k=0; k < s2; k++)
	      {
		double d = suds_indiv_t::bread_dbl( IN1 );
	      }  

	  // n
	  int j = suds_indiv_t::bread_int( IN1 );
	  
	  // labels
	  for (int i=0;i<ng;i++)
	    {
	      std::string s = suds_indiv_t::bread_str( IN1 );
	    }	  
	}
      
      //
      // QDA model?
      //
      
      if ( has_qda )
	{
	  
	  // number of groups
	  const int ng = suds_indiv_t::bread_int( IN1 );
	  
	  // number of variables
	  const int nv = suds_indiv_t::bread_int( IN1 );
	  
	  // priors	  
	  for (int i=0;i<ng;i++)
	    {
	      double d = suds_indiv_t::bread_dbl( IN1 );
	    }

	  // rows (redundant, but keep)
	  for (int i=0;i<ng;i++)
	    {
	      int j = suds_indiv_t::bread_int( IN1 );
	    }
	  
	  // counts
	  for (int i=0;i<ng;i++)
	    {
	      const std::string s = suds_indiv_t::bread_str( IN1 );
	      int j = suds_indiv_t::bread_int( IN1 );
	    }
	  
	  // means
	  for (int i=0; i<ng; i++)
	    for (int j=0; j<nv; j++)
	      {
		double d = suds_indiv_t::bread_dbl( IN1 );
	      }

	  // scaling
	  for (int i=0; i<ng; i++)
	    {
	      for (int j=0; j < nv; j++)
		for (int k=0; k < nv; k++)
		  {
		    double d = suds_indiv_t::bread_dbl( IN1 );
		  }
	    }
	  
	  // ldet	  
	  for (int i=0;i<ng;i++)
	    {
	      double d = suds_indiv_t::bread_dbl( IN1 );
	    }	  

	  // n
	  int j = suds_indiv_t::bread_int( IN1 );

	  // labels
	  for (int i=0;i<ng;i++)
	    {
	      std::string l1 = suds_indiv_t::bread_str( IN1 );
	    }
	}

      
      //
      // Original features: X 
      //

      // make space for new data 
      int r = mega.X.rows();      
      int r1 = r;

      if ( added == 0 )  // need to set cols too the first time
	mega.X = Eigen::MatrixXd::Zero( tnve , first_nf );
      else	
	mega.X.conservativeResize( r + tnve , Eigen::NoChange );

      for (int i=0;i<tnve;i++)
	{
	  for (int j=0;j<tnf;j++)
	    {
	      mega.X(r,j) = suds_indiv_t::bread_dbl( IN1 );
	    }
	  ++r; // next row
	}


      //
      // next individual
      //
      
      if ( add_person ) ++added;
      
      ++n_indiv;
      
      logger << "  " << added << " trainer compiled (" << n_indiv << " considered)\n";
      
    }

  IN1.close();


  //
  // Done reading and combining features
  //
  
  int ecnt = mega.obs_stage.size();
  logger << "  read " << ecnt << " epochs " << mega.X.rows() << " x " << mega.X.cols() << "\n";
  
  // set labels 
  mega.y = suds_t::str( mega.obs_stage );
  
  mega.counts.clear();
  for (int i=0;i<mega.y.size();i++) mega.counts[mega.y[i]]++;
  std::map<std::string,int>::const_iterator cc = mega.counts.begin();
  logger << "  epoch counts:";
  while ( cc != mega.counts.end() )
    {
      logger << " " << cc->first << ":" << cc->second ;
      ++cc;
    }
  logger << "\n";


  
  //
  // Perform SVD 
  //
  
  int rc = 0;

  annotation_set_t dummy_annotations;
  edf_t dummy( &dummy_annotations ); // this is not used/needed in any of the proc_() calls below
  // i.e. as the feature matrix and stages are already built/compiled
  mega.nve = ecnt;

  suds_helper_t helper( dummy , param );
  helper.ne = helper.nge = ecnt;
  helper.ns = 0; // should not matter
  helper.has_prior_staging = true;
  helper.retained.resize( ecnt , true );
  helper.valid.resize( ecnt , true );

  logger << "  performing primary SVD " << suds_t::nc << " components\n";
  rc = mega.proc_main_svd( &helper );
  if ( rc == 0 ) Helper::halt( "problem in proc_main_svd()" );

  //
  // drop components that do not track well w/ stage
  //

  logger << "  dropping uninformative columns...\n";
  rc = mega.proc_prune_cols( &helper );
  if ( rc == 0 ) Helper::halt( "problem in proc_prune_cols()" );

  //
  // get class label counts
  //
  logger << "  compiling stage labels...\n";
  rc = mega.proc_class_labels( &helper );
  if ( rc == 0 ) Helper::halt( "problem in proc_class_labels()" );

  //
  // some final metrics
  //

  int ne = mega.proc_coda( &helper );
  
  
  //
  // Fit trainer model
  //
  
  logger << "  fitting QDA model...\n";

  qda_t qda( mega.y , mega.U );

  qda_model_t fit = qda.fit( suds_t::flat_priors );
  
  //
  // Save trainer model
  //

  logger << "  writing model fit to " << outfile << ".fit \n";
  
  fit.write( outfile + ".fit" );


  //
  // Write Hjorth 95% limits for this entire set
  //

  logger << "  writing Hjorth outlier values to " << outfile << ".hjorth\n";

  std::ofstream H1( ( outfile + ".hjorth" ).c_str() , std::ios::out );

  H1 << first_ns << "\n"; 
  
  for (int s=0; s<first_ns; s++)
    {
      double h1_mean = MiscMath::mean( h1_means[s] );
      double h2_mean = MiscMath::mean( h2_means[s] );
      double h3_mean = MiscMath::mean( h3_means[s] );
      
      double h1_sd = sqrt( MiscMath::mean( h1_vars[s] ) );
      double h2_sd = sqrt( MiscMath::mean( h2_vars[s] ) );
      double h3_sd = sqrt( MiscMath::mean( h3_vars[s] ) );
      
      H1 << h1_mean << " " << h1_sd << " "
	 << h2_mean << " " << h2_sd << " " 
	 << h3_mean << " " << h3_sd << "\n";
      
    }
  
  
  H1.close();

  //
  // Write V and W matrices for this SVD 
  //

  logger << "  writing SVD V and W matrices to " << outfile << ".svd\n";

  std::ofstream SVD1( ( outfile + ".svd" ).c_str() , std::ios::out );

  SVD1 << mega.W.size() << "\n"
       << mega.W << "\n"
       << mega.V.rows() << " " << mega.V.cols() << "\n"       
       << mega.V << "\n";
  
  SVD1.close();

  //
  // All done
  //


  logger << "  in total, converted " << n_indiv << " trainers (" << ecnt << " epochs)\n";

  
  
}
