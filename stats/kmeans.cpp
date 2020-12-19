
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

#include "stats/kmeans.h"

#include "miscmath/crandom.h"

double kmeans_t::randf(double m)
{
  return m * rand() / (RAND_MAX - 1.);
}


double kmeans_t::dist2( const point_t & a, const point_t & b )
{
  double d = 0;
  for (int i=0;i<n;i++) d += ( a.x[i] - b.x[i] ) * ( a.x[i] - b.x[i] ) ; 
  return d;
}

int kmeans_t::nearest( const point_t & pt, // point
		       const std::vector<point_t> & cent,  // current means
		       double * d2 ) // distance
{
  
  int min_d = std::numeric_limits<int>::max();
  
  int min_i = pt.group;
  
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
    }
  
  // also return the value?
  if ( d2 != NULL ) *d2 = min_d;

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
  cent[0] = pts[ rand() % len ];
  
  // subsequenrly, pick others based on distances to this (kmeans++ algorithm)    
  for (int n_cluster = 1 ; n_cluster < n_cent; n_cluster++) 
    {
      sum = 0;
      
      for (int j = 0 ; j < len ; j++ )
	{
	  double pd;
	  int nn = nearest(pts[j], cent, &pd );
	  d[j] = pd;
	  sum += d[j];
	}
      
      sum = randf(sum);
      
      for (int j = 0 ; j < len ; j++ )
	{
	  if ( (sum -= d[j] ) > 0 ) continue;
	  cent[ n_cluster ] = pts[ j ];
	  break;
	}	
    }  
  // set class for all points 
  for (int j = 0 ; j < len ; j++ )
    pts[j].group = nearest( pts[j], cent, NULL );
  
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

  // class means
  return ret;
  
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
  

  // return means
  
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
