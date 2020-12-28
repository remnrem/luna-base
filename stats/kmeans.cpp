
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


#include <iomanip>

#include "stats/kmeans.h"
#include "stats/statistics.h"
#include "miscmath/crandom.h"

#include "miscmath/crandom.h"

double kmeans_t::randf(double m)
{
  //return m * rand() / (RAND_MAX - 1.);
  return CRandom::rand() * m;
}


double kmeans_t::dist2( const point_t & a, const point_t & b )
{
  double d = 0;
  for (int i=0;i<n;i++) d += ( a.x[i] - b.x[i] ) * ( a.x[i] - b.x[i] ) ; 
  return d;
}

int kmeans_t::nearest( const point_t & pt, // point
		       const std::vector<point_t> & cent,  // current means
		       double * d2 , // return distance
		       int * lim ) // optional, only consider clusters < lim 
{

  // for point 'pt', find the nearest cluster center (in cent),
  // considering only up to *lim, if lim != NULL
  
  double min_d = std::numeric_limits<double>::max();
  int min_i;
  
  int i = 0;
  
  std::vector<point_t>::const_iterator cc = cent.begin();
  while ( cc != cent.end() )
    {

      double d = dist2( *cc, pt );

      if ( d < min_d )
	{
	  min_d = d; 
	  min_i = i;	  
	}
      
      ++i;
      ++cc;

      // stop here? (i.e. if durig seeding only considering a subset of the (currently-initiated) clusters
      if ( lim && i == *lim )
	{
	  //std::cout << "stopping at " << i << " of " << cent.size() << "\n";
	  break;
	}
    }
  
  // also return the value?
  if ( d2 != NULL ) *d2 = min_d;
  //  std::cout << "near grp = " << min_i << "\n";
  // return group index
  return min_i;
}


void kmeans_t::kpp( std::vector<point_t> & pts,  
		    std::vector<point_t> & cent )
{
  
  // #define for_len for (j = 0, p = pts; j < len; j++, p++)
  
  // number of data points
  int len = pts.size();
  
  // number of clusters
  int n_cent = cent.size();
    
  int i, j;
  int n_cluster;
  
  double sum;
  
  std::vector<double> d( len );
  
  // pick an initial seed at random from the data 
  //cent[0] = pts[ rand() % len ];


  int first_random_seed = CRandom::rand( len ) ;
  //  std::cout << " INIT == " << first_random_seed + 1 << "\n";
  cent[0] = pts[ first_random_seed ];
  
  // subsequenrly, pick others based on distances to this (kmeans++ algorithm)
  // where prob. of being selected is proportional to the distance to the nearest
  // of the previously selected clusters
  for (int n_cluster = 1 ; n_cluster < n_cent; n_cluster++) 
    {
      sum = 0;
      
      for (int j = 0 ; j < len ; j++ )
	{
	  double pd;
	  // find nearest neighbor, but only looking up to the n_cluster (not including)
	  int nn = nearest(pts[j], cent, &pd , &n_cluster );
	  d[j] = pd;
	  sum += d[j];
	}

      // std::cout << "n_cluster = " << n_cluster << "\n";
      // std::cout << "pre sum = " << sum <<"\n";

      // select a number between 0 and sum
      sum = randf(sum);

      //      std::cout << "post sum = " << sum <<"\n";

      // for (int j = 0 ; j < len ; j++ )
      // 	std::cout << " all - pnt = " << j << " " << d[j] << "\n";

      
      // select point based on weighted distances (squared euclidian)
      for (int j = 0 ; j < len ; j++ )
	{
	  //std::cout << " pnt = " << j << " " << d[j] << "\n";
	  if ( (sum -= d[j] ) > 0 ) continue;
	  cent[ n_cluster ] = pts[ j ];
	  //	  std::cout << " INIT == " << j+1 << "\n";
	  break;
	}	
    }  

  std::map<int,int> cnts;
  // set class for all points 
  for (int j = 0 ; j < len ; j++ )
    {      
      pts[j].group = nearest( pts[j], cent, NULL );
      cnts[ pts[j].group ]++;
    }

  // std::map<int,int>::const_iterator cc = cnts.begin();
  // while ( cc != cnts.end() )
  //   {
  //     std::cout << "init " << cc->first << " = " << cc->second << "\n";
  //     ++cc;
  //   }
  
}


Data::Matrix<double> kmeans_t::lloyd( const Data::Matrix<double> & X , int nk , std::vector<int> * sol )
{
  const int nr = X.dim1();
  const int nc = X.dim2();
  
  // convert to  std::vector<point_t> 
  std::vector<point_t> d( nr );
  
  for (int r=0; r<nr; r++) 
    d[r] = point_t( X.row(r) );
  
  std::vector<point_t> cent = lloyd( d , nk );
  
  // get centroid means
  Data::Matrix<double> ret( nk , nc );
  for (int k=0; k<nk; k++)
    for (int c=0; c<nc; c++)
      ret(k,c) = cent[k].x[c];

  // get solutions for each observation
  if ( sol != NULL )
    {
      sol->resize( nr );
      for (int r=0; r<nr; r++) (*sol)[r] = d[r].group;
    }

  // get variance explained

  variance_explained( d , cent );
  std::cout << "VE = " << between << " " << within << " B = " <<  ( between / ( between + within ) ) << " W = " << ( within / ( within + between ) ) << "\n";
  
  
  // class means (but transposred to channels x classes)
  return Statistics::transpose(  ret );
  
}


std::vector<kmeans_t::point_t> kmeans_t::lloyd( std::vector<kmeans_t::point_t> & pts, int nk )
{

  if ( pts.size() < 2 ) Helper::halt( "passing only 2 points to lloyd()" );

  //
  // track number of variables, for dist() calculations
  //

  n = pts[0].x.size();

  //
  // cluster means
  //
  
  std::vector<point_t> cent( nk );
  for (int k=0; k<nk; k++) cent[k] = point_t( n );


  
  //
  // Use kmeans++ to initialize 
  //

  kpp( pts , cent );


  //
  // begin k-means iterations
  //

 
  int len = pts.size();


  int changed = 0;
  int niter = 0;

  do {

    // track iterations
    ++niter;

    /* group element for centroids are used as counters */
    
    std::vector<point_t>::iterator cc = cent.begin();
    while ( cc != cent.end() )
      {
	cc->clear();
	++cc;
      }
    
    // for_len
    for (int j = 0; j < len; j++)
      {
	// get current class for this observation
	point_t & c = cent[ pts[j].group ];
	c.add( pts[j] );
      }
    
    // scale
    cc = cent.begin();
    while ( cc != cent.end() )
      {
	cc->scale();
	++cc;
      }
    
    changed = 0;
    
    // find closest centroid of each point
    for (int j = 0; j < len; j++)
      {
	point_t & p = pts[j];
	int min_i = nearest( p, cent, NULL );
	
	if (min_i != p.group) 
	  {
	    changed++;
	    p.group = min_i;
	  }
      }

    //    std::cout << "changed = " << changed << " " << len << " " << ( len >> 10 ) << "\n";
    
  } while ( changed > (len >> 10)); /* stop when 99.9% of points are good */
  
  std::cout << "completed in " << niter << " iterations\n";
  
  // populate class assignments    
  int i = 0;
  std::vector<point_t>::iterator cc = cent.begin();
  while ( cc != cent.end() )
    {
      cc->group = i++;
      ++cc;
    }


  // get variance explained
  
  return cent;
}


Data::Matrix<double> kmeans_t::kmeans( const Data::Matrix<double> & X , const int nk , std::vector<int> * sol )
{  
  return lloyd( X , nk , sol );    
}


void kmeans_t::test2()
{
  
  Data::Matrix<double> X( 100 , 10 );
  for (int i=0;i<50;i++)
    for (int j=0;j<5;j++)
      X(i,j) += 2;
  
  for (int i=0;i<50;i++)
    for (int j=0;j<5;j++)
      X(i,j) += CRandom::rand(10);
  
  Data::Matrix<double> km = lloyd( X , 2 );
  
   std::cout << "KM\n" << km.print() << "\n";
  
}


void kmeans_t::variance_explained( const std::vector<point_t> & pts , const std::vector<point_t> & cent )
{
  
  point_t grand_mean( n );
  
  const int nr = pts.size();
  const int nk = cent.size();
  const int nc = n;

  // get grand mean
  for (int r=0; r<nr; r++)
    for (int c=0; c<nc; c++)
      grand_mean.x[c] += pts[r].x[c];
  for (int c=0; c<nc; c++)
    grand_mean.x[c] /= (double)nr;

  // get total SS
  double tot_ss = 0;
  for (int r=0; r<nr; r++)
    tot_ss += dist2( grand_mean , pts[r] );

  // get within SS
  within_ss.resize( nk );
  Data::Vector<double> counts( nk );

  for (int r=0; r<nr; r++)
    {
      int group = pts[r].group ;
      counts[ group ]++;
      within_ss[ group ] += dist2( pts[r] , cent[ group ] );
    }

  within = 0;
  for (int k=0; k<nk; k++)
    {
      within_ss[k] /= counts[k];
      within += within_ss[k] ;
    }
  
  between = tot_ss - within;
  
}




//
// EEG modified K means
//


modkmeans_all_out_t modkmeans_t::fit( const Data::Matrix<double> & data )
{
  // fix:
  nreps = 1 ;


  
  // input is currently N x C, 
  const int N = data.dim1();
  const int C = data.dim2();

  // copy 
  X.resize( N , C );
  for (int i=0;i<N;i++)
    for (int j=0;j<C;j++)
      X(i,j) = data(i,j);
  
  //
  // normalize data by the average STD of channels?
  // (e.g. can be helpful if pooling across individuals)
  //
  
  if ( normalize )
    {
      Data::Vector<double> sdev = Statistics::sdev( X , Statistics::mean( data ) ) ;
      double fac = Statistics::mean( sdev );
      for (int r=0; r<N; r++)
	for (int c=0;c<C;c++)
	  X(r,c) /= fac;	
    }

  
  //
  // We need channels x samples 
  //

  //X.transposeInPlace();
  X = Statistics::transpose( X );

  
  //
  // constant total SS
  //
    
  double const1 = 0;

  for (int i=0; i<C; i++)
    for (int j=0; j<N; j++)
      const1 += X(i,j) * X(i,j);
  
    
  //
  // Use GEV as the goodness of fit metric
  //


  //
  // Global field power GFP
  //

  Data::Vector<double> GFP = Statistics::sdev( X , Statistics::mean( X ) );
  double GFP_const = Statistics::sum_squares( GFP );
  double GEV_opt = 0;

  
  // else
  // case 'CV'
  // sig2_mcv_opt = inf;
  // case 'dispersion'
  // W_opt = inf;

  //
  // Final results
  //

  modkmeans_all_out_t results;


  bool new_best = false;
  
  //
  // Iterate over all K values
  //

  for (int ki=0; ki<ks.size(); ki++)
    {

      const int K = ks[ki];
      
      // K_ind is just index for --> == ki

      // Finding best fit amongst a given number of restarts based on selected
      // measure of fit

      double GEV_best = 0;

      //
      // For each replicate
      //


      for (int r = 0; r < nreps ; r++)
	{

	  std::cout << "rep = " << r << "\n";
	  

	  //
	  // 1) get segmentation
	  //

	  //            % The original Basic N-Microstate Algorithm (Table I in [1])
	  // [A,L,Z,sig2,R2,MSE,ind] = segmentation(X,K,const1,opts);
	  
	  modkmeans_out_t result = segmentation( X , K , const1 );
	 

	  
	  // ----------------------------------------------------------------------------------------------
	  //
	  // 2) segmentation smoothing?
	  //
	  //             [L,sig2,R2,MSE,ind] = smoothing(X,A,K,const1,opts);

	  // TODO...	  

	  // END (SMOOTHING) -----------------------------------------------------------------------------
	  


	  //
	  // Check for better fit
	  //

	  std::cout << "map_corr ...\n";

	  std::cout << " A = " << result.A.dim1() << " " << result.A.dim2() << "\n";
	  std::cout << " L = " << result.L.size() << "\n";

	  // map_corr = columncorr(X,A(:,L));

	  Data::Vector<double> map_corr( N );

	  for (int j=0;j<N;j++)
	    {
	      double r = Statistics::correlation( X.col(j) , result.A.col( result.L[j] ) );
	      if  ( r < -1 ) Helper::halt( "problem with modkmeans" );
	      map_corr[j] = r;
	    }

	  //GEV = sum((GFP.*map_corr).^2) / GFP_const;
	  double GEV = 0;
	  for (int j=0; j<N; j++)
	    GEV += ( GFP[j] * map_corr[j] ) * ( GFP[j] * map_corr[j] ) ;
	  GEV /= GFP_const;
	  	  
	  if ( GEV > GEV_best )
	    {
	      new_best = true;
	      GEV_best = GEV;
	    }

	  std::cout << "GEV = " << GEV << "\n";
	  
	  // Update if new best found

	  if ( new_best )
	    {
	      std::cout << "is a new best\n";
	      results.kres[K] = result;
	      
	      // final.A[K] = resuA;
	      // res.L_all[K] = L;
	      // res.sig2_all[K] = sig2;
	      // res.Z_all[K] = Z;
	      // res.R2_all[K] = R2;
	      // res.MSE_all[K] = MSE;

	      new_best = false;
	    }
	  
	  std::cout << " next K\n";
	  //
	  // Next replicate for this K
	  //

	} 

      //
      // After finishing all replicates for this K 
      //

      //sig2_mcv[k] = sig2_all[k] * ( (C-1)^-1  *   (C-1-K))^-2;
      results.kres[K].sig2_modk_mcv = results.kres[K].sig2 * pow( pow(C-1,-1) * (C-1-K) , -2 );
      
      //
      // Checking for best fit for different values of K: GEV
      //
      
      if ( GEV_best > GEV_opt )
	{
	  new_best = true;
	  GEV_opt = GEV_best;
	}

      std::cout << "GEV_best = " << GEV_best << "\n";
      
      //
      // Update to track optimal across all K considered
      //
      
      if ( new_best )
	{
	  std::cout << "is a new K best\n";
	  results.A = results.kres[K].A;
	  results.L = results.kres[K].L;
	  results.K = K;
	  new_best = false;
	}
      
    } // next K

  
  // save additional results

  // Res.Z_all = Z_all;
  // Res.A_all = A_all;
  // Res.L_all = L_all;
  // Res.R2 = R2_all;
  // Res.MSE = MSE_all;
  // Res.sig2_modk = sig2_all;
  // Res.sig2_modk_mcv = sig2_mcv;
  // Res.K_act = K_act;

  return results;  
  
}


modkmeans_out_t modkmeans_t::segmentation( const Data::Matrix<double> & X , int K , double const1 )
{	  
  
  max_iterations = 5;
  
  const int C = X.dim1();
  const int N = X.dim2();
  
  // Step 1

  double sig2_old = 0;
  double sig2 = std::numeric_limits<double>::max();

  // Step 2a

  // selecting K random timepoints (0 to N-1)  to use as initial microstate maps
  Data::Matrix<double> A( C , K );
  std::vector<int> L( N );

  
  if ( 0 )
    {
      std::set<int> selected;
      while ( 1 ) {
		int pick = CRandom::rand( N );
		if ( selected.find( pick ) != selected.end() ) continue;
		int j = selected.size();
		for (int i=0; i<C; i++)
		  A(i,j) = X(i,pick);
		selected.insert( pick );
		if ( selected.size() == K ) break; 
	      }
	    }
	  else
	    {
	      // fix 
	      // 105907        4775       41538
	      if ( K == 3 )
		{
		  std::vector<int> picks = { 105907     ,   4775   ,    41538  } ;
		  for (int p=0;p<3;p++)
		    for (int i=0; i<C; i++)
		      A(i,p) = X(i,picks[p]-1);	      
		}
	      else if ( K == 4 )
		{
		  std::cout << " K4!\n";
		  std::vector<int> picks = { 22222 , 105907     ,   4775   ,    41538  } ;
		  for (int p=0;p<4;p++)
		    for (int i=0; i<C; i++)
		      A(i,p) = X(i,picks[p]-1);	      
		}
	      
	    }

	  
	  std::cout << "norm each channel\n";
	  
	  // normalize each channel

	  // A = bsxfun(@rdivide,A,sqrt(diag(A*A')));% normalising
	  // CxK . KxC -> CxC, but only need to evaluate the diagonal elements
	  
	  for (int i=0; i<C; i++)
	    {
	      // for the i'th diagonal element, get the normalization factor (for channel/row of A)
	      double norm = 0;
	      for (int k=0; k<K; k++)
		norm += A(i,k) * A(i,k); 
	      norm = sqrt( norm );
	      for (int k=0;k<K;k++)
		A(i,k) /= norm;
	    }

	  std::cout << "init iterations\n";

	  // initialize
	  int ind = 0; // iteration counter
	
	  // Iterations (step 3 to 6)	  
		  
	  while ( fabs(sig2_old - sig2) >= threshold * sig2 && max_iterations > ind )
	    {
	      
	      ++ind;

	      std::cout << "iteration = " << ind << " ( of max " << max_iterations << "); ";
	      std::cout << " sig2 = " << fabs(sig2_old - sig2) << "\t" << threshold * sig2 << "\n";
	      
	      sig2_old = sig2;
    
	      // Step 3

	      // Z = A'*X;
	      
	      // A'  KxC CxN --> Z = KxN
	      // do transpose in place, by hand rather than two function calls
	      
	      Data::Matrix<double> Z( K , N );
	      for (int i=0; i<K; i++)
		for (int j=0; j<N; j++)
		  for (int k=0; k<C; k++)
		    Z(i,j) += A(k,i) * X(k,j); // nb. transpose of A

	      // std::cout << "SX = " << X.dim1() << " " << X.dim2() << "\n";
	      // std::cout << A.print( "A" , 10 , 0 ) << "\n";
	      // std::cout << X.print( "X" , 10 , 10 ) << "\n";
	      // std::cout << Z.print( "Z" , 0 , 10 ) << "\n";

	      // [~,L] = max(Z.^2);
	      // get max and put in L
	      
	      std::map<int,std::vector<int> > K_idx; // track freq of each class, for below
	      for (int j=0; j<N; j++)
		{
		  int idx = 0;
		  double max = Z(0,j) * Z(0,j);
		  for (int i=1; i<K; i++)
		    {
		      double t = Z(i,j) * Z(i,j);
		      if ( t > max )
			{
			  idx = i;
			  max = t;
			}
		    }
		  L[j] = idx;
		  K_idx[idx].push_back( j );		  
		}


	      // TMP check on L
	      // std::map<int,int> qc;
	      // for (int i=0;i<L.size();i++) qc[L[i]]++;
	      // std::map<int,int>::const_iterator qq = qc.begin();
	      // while ( qq != qc.end() )
	      // 	{
	      // 	  std::cout << "L " << qq->first << " = " << qq->second << "\n";
	      // 	  ++qq;
	      // 	}
	      

	      
	      // Step 4
	      
	      for (int k=0; k<K; k++)
		{
		  // std::cout << "Step4 = k = " << k << " " << K << "\n";
		  
		  // std::cout << " K_idx[k].size()  = " << K_idx[k].size()  << "\n";
		  
		  // A is Cx(N==K)
		  if ( K_idx[k].size() == 0 )
		    {
		      // A =  % no members of this microstate
		      // A(:,k) = 0;
		      for (int i=0;i<C;i++) A(i,k) = 0;
		    }
		  else
		    {
		      // extract channel maps for this K into S matrix, which is CxC
		      Data::Matrix<double> S( C , C );
		      //S = X(:,k_idx)*X(:,k_idx)';
		      // CxNk x NkxC
		      
		      // hmm.  element access via operator() is slow in the inner loop
		      // need to investigate that.. but for now, use a raw vec/vec array
		      
		      const std::vector<int> & kidx = K_idx[k];
		      const int nk = kidx.size();
		      
		      // SLOW....
		      // for (int i=0; i<C; i++)
		      // 	for (int j=0; j<C; j++)
		      // 	  for (int l=0; l<nk; l++)
		      // 	    S(i,j) += X(i,kidx[l]) * X(j,kidx[l]);

		      // make XX transpose for potentially faster maxtrix mult
		      std::vector<std::vector<double> > XX( C );

		      for (int i=0; i<C; i++)
			XX[i].resize( nk );
		      for (int i=0; i<C; i++)
			for (int j=0; j<nk; j++)
			  XX[i][j] = X(i,kidx[j]);
		      
		      for (int i=0; i<C; i++)
			for (int j=0; j<C; j++)
			  {
			    double tmp = 0;
			    for (int l=0; l<nk; l++)
			      tmp += XX[i][l] * XX[j][l];
			    S(i,j) = tmp;
			  }
		      
		      // finding eigenvector with largest value and normalising it
		      // [eVecs,eVals] = eig(S,'vector');
		      // [~,idx] = max(abs(eVals));
		      // A(:,k) = eVecs(:,idx);
		      // A(:,k) = A(:,k)./sqrt(sum(A(:,k).^2));
		      
		      bool okay = true;
		      
		      //std::cout << S.print("S",10,10) << "\n";
		      
		      Statistics::Eigen eigen = Statistics::eigenvectors( S , &okay );

		      if ( ! okay ) Helper::halt( "problem in modkmeans()" );

		      // eigen values:
		      //  eigen.d  (eigenvalues)
		      //  eigen.z  (eigenvectors)
		      
		      int mx_idx = 0;
		      double mx = fabs( eigen.d[0] ) ;
		      const int ne = eigen.d.size();
		      if ( eigen.d.size() != C ) Helper::halt( "problem in modkmeans()" );
		      for (int i=1; i<C; i++)
			{
			  if ( fabs( eigen.d[i] ) > mx )
			    {
			      mx = eigen.d[i];
			      mx_idx = i;
			    }
			}

		      // copy into A, and normalize		      
		      double norm = 0;
                      for (int i=0; i<C; i++)
			{
			  A(i,k) = eigen.z(i,mx_idx);
			  norm += A(i,k) * A(i,k);
			}
		      norm = sqrt(norm);
		      for (int i=0; i<C; i++)
			A(i,k) /= norm;

		    }
		
		} // next 'k' of K
	      
	      // Step 5

	      //sig2 = (const1 - sum( sum( A(:,L). *X ).^2) ) / (N*(C-1));
	      // L contains index of X; and element operations: so
	      
	      Data::Vector<double> g(N) ; // sum( sum( A(:,L). *X ).^2) )
	      
	      for (int j=0; j<N; j++)
		{
		  for ( int i=0; i<C; i++)
		    g[j] += A(i,L[j] ) * X(i,j);
		  g[j] *= g[j];
		}

	      double gsum = Statistics::sum( g );

	      sig2 = ( const1 - gsum ) / (double)(N*(C-1));

	      
	    } // end of iterations
	  
	  //%% Saving solution converged on (step 7 and 8)
	      
	  // Step 7
	  // Z = A'*X; % NOTE, not setting non-activated microstates to zero
	  
	  Data::Matrix<double> Z( K , N );
	  for (int i=0; i<K; i++)
	    for (int j=0; j<N; j++)
	      for (int k=0; k<C; k++)
		Z(i,j) += A(k,i) * X(k,j); // nb. transpose of A                                                                                                                    

	  // [~,L] = max(Z.^2);
	  
	  for (int j=0; j<N; j++)
	    {
	      int idx = 0;
	      double max = Z(0,j) * Z(0,j);
	      for (int i=1; i<K; i++)
		{
		  double t = Z(i,j) * Z(i,j);
		  if ( t > max )
		    {
		      idx = i;
		      max = t;
		    }
		}
	      L[j] = idx;	      
	    }

	  
	  //  Step 8
	  // sig2_D = const1 / (N*(C-1));
	  // R2 = 1 - sig2/sig2_D;
	  // activations = zeros(size(Z));
	  // activations = zeros(size(Z));
	  // for n=1:N; activations(L(n),n) = Z(L(n),n); end % setting to zero
	  // MSE = mean(mean((X-A*activations).^2));
	  // end

	  double sig2_D = const1 / (double)(N*(C-1));
	  
	  double R2 = 1.0 - sig2/sig2_D; 


	  // MSE = mean(mean((X-A*activations).^2));                                                                                                                                             //    X - A * act
	  //    CxN -  CxK * KxN  

	  // but only one row of 'act' is non-zero for a given column;
	  //  therefore, can reduce the matrix multipleication for: A * activations
	  //  and directly go from Z & A 

	  Data::Matrix<double> XX = X;
	  for (int i=0;i<C;i++)
	    for (int j=0;j<N;j++)
	      {
		XX(i,j) -= A(i,L[j]) * Z(L[j],j);
		XX(i,j) *= XX(i,j);
	      }

	  double MSE = Statistics::mean( Statistics::mean( XX ) );

	  //
	  // Package up results	  
	  //
	  
	  modkmeans_out_t result;
	  result.A = A;
	  result.L = L;
	  result.Z = Z;
	  result.R2 = R2;
	  result.sig2 = sig2;
	  result.MSE = MSE;
	  result.iter = ind;
	  return result;
 
}	  





