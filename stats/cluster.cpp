
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


#include "cluster.h"
#include "helper/logger.h"

extern logger_t logger;

#include <map>

cluster_solution_t cluster_t::build( const Data::Matrix<double> & D , const int preK )
{
 
  bool calc_silhouette = preK == 0 ;
  
  const int ni = D.dim1();
  
  const int max_cluster_N = preK > 0 ? preK : ni ;
  
  const int max_cluster_size = 0; 
  
  //
  // seed initial solution
  //
  
  // cluster --> individuals
  std::vector<std::vector<int> > cl;
  
  for (int i=0;i<ni;i++)
    {
      std::vector<int> t(1);
      t[0] = i;
      cl.push_back(t);
    }

  
  cluster_solution_t final_sol;

  int c=1;
  
  bool done=false;
  
  // Matrix of solutions  
  std::vector< std::vector<int> > sol(ni);
  for (int i=0;i<ni;i++) sol[i].resize(ni);
  
  std::vector<double> hist(1);
  
  // Build solution
  for (int i=0; i<cl.size(); i++)
    for (int j=0; j<cl[i].size(); j++)
      sol[cl[i][j]][0] = i;
  
  // silhouette metric
  int best_sil_idx = 0;
  int best_sil_K = 0;
  double best_sil = -1; // lowest possible value      


  while( ! done )
    {
      
      double dmin = 999;
      
      int imin=-1;
      int jmin=-1;
      
	  
      // 1. Find min/max distance between pairable clusters
      
      for (int i=0; i<cl.size()-1; i++)
	for (int j=i+1; j<cl.size(); j++)
	  {	
	    
	    if ( max_cluster_size==0 ||  
		 (( cl[i].size()+cl[j].size()) <= max_cluster_size) )
	      {
	    
		//double d = groupAvgLink(mdist,cl[i],cl[j]) ;
		double d = cldist( D , cl[i], cl[j]);	    
		
		
		// Are these individuals/clusters more similar?
		
		if ( d < dmin )
		  {
		    imin=i;
		    jmin=j;
		    dmin=d;
		  }
	      }
	  }
  
      if (imin==-1) {
	done=true;
	//printLOG("Cannot make clusters that satisfy constraints at step "+int2str(c)+"\n");	
	goto done_making_clusters;
      }
	  
      // Save merge distance 
      hist.push_back(dmin);
      

      // 2. Join these clusters 
      for(int j=0;j<cl[jmin].size();j++)
	cl[imin].push_back(cl[jmin][j]);
      cl.erase(cl.begin()+jmin);
      if (cl.size()==1 || cl.size()==max_cluster_N) done=true;
      

      // List entire sample
      std::cout << "Merge step " << c << "\t" << hist[c] << "\n";
	  
      // Build solution
      for (int i=0; i<cl.size(); i++)
	for (int j=0; j<cl[i].size(); j++)
	  sol[cl[i][j]][c] = i;
            
      // Calculate average within/between cluster distances
//       double between = 0, within = 0; 
//       int withinN = 0, betweenN = 0;
      
//       for (int j1=0; j1<sol.size(); j1++)
// 	for (int j2=0; j2<sol.size(); j2++)
// 	  {
// 	    if (j1 < j2)
// 	      {
// 		if(sol[j1][c] == sol[j2][c])
// 		  {
// 		    within += D[j2][j1];
// 		    withinN++;		  
// 		  }
// 		else
// 		  {
// 		    between += D[j2][j1];
// 		    betweenN++;		  
// 		  }
// 	      }
// 	  }
      
//       std::cout << "\t" << ( betweenN > 0 ? between/(double)betweenN : 0 )
// 		<< "\t" << ( withinN > 0 ? within/(double)withinN  : 0 )
// 		<< "\t" << ( betweenN > 0 && withinN > 0 ? ( between/(double)betweenN )  / ( within/(double)withinN ) : 0 )
// 		<< "\n";
	  
    

      //
      // silhouette method to determine best # of clusters
      //
    
      if ( calc_silhouette )
	{
	 
	  const int K = cl.size();
	  
	  std::vector<double> sil( ni , 0 );
	  
	  //	  if ( K >= 2 && K < ni  )
	  if ( K >= 2 && K <= 20  )
	    {
	      
	      for (int i=0;i<ni;i++)
		{
		
		  // this individual currently assigned to:
		  const int assign_k = sol[i][c];
		  
		  int n = cl[assign_k].size();
		  
		  // is this cluster size 1? 
		  if ( n == 1 ) sil[i] = 0; 
		  else
		    {
		      
		      // a = average distance to all other members of this cluster
		      double a = 0;
		      
		      for (int j=0;j<n;j++) 
			if ( cl[assign_k][j] != i ) // skip this indiv. 
			  a += D( i , cl[assign_k][j] );
		      a /= (double)(n-1);
		      
		      // b = smallest average distance to all members of another cluster
		      double min_b = 99999;
		      
		      for (int k=0;k<K;k++)
			{
			  if ( k == assign_k ) continue; // skip this cluster
			  double b = 0;
			  int n = cl[k].size();
			  for (int j=0;j<n;j++) b += D( i , cl[k][j] );
			  b /= (double)n;
			  if ( b < min_b ) min_b = b;
			} // next cluster
		      
		      // calc. sil
		      double s = ( min_b - a ) / ( a > min_b ? a : min_b ); 
		      
		      sil[i] = s;
		    }
		  
		  //std::cout << "K,I,S = " << K << "\t" << i << "\t" << assign_k << "\t" << sil[K][i] << "\n";
		}
	      
	      // get average silhouette score
	      
	      double avg_sil = 0;
	      for (int i=0;i<ni;i++) avg_sil += sil[i];
	      avg_sil /= (ni);
	      //std::cout << " avg sil = " << K << "\t" << avg_sil << "\n";
	      
	      if ( avg_sil > best_sil )
		{
		  best_sil = avg_sil;
		  best_sil_idx = c;
		  best_sil_K = K;
		  //std::cerr << "updating " << best_sil << " " << c << " K=" << K << "\n";
		}
	    }

	}

      // Next merge
      c++;
    }
  
 done_making_clusters:

   
     
  //////////////////////////////////
  // Best solution is final solution
  
  int best = cl.size();
  int best_c = hist.size() - 1; // last sol

  // use silhouette?
  if ( calc_silhouette )
    {
      best = best_sil_K;
      best_c = best_sil_idx;
    }

  final_sol.k = best;
  final_sol.best.resize( ni );  
  
  logger << " stopped clustering at K=" << best << "\n";

  // copy best solution
  for (int j=0; j<sol.size(); j++)
    final_sol.best[j] = sol[j][best_c] ;

  
//   for (int i=0; i<cl.size(); i++)
//     {
      
// //       std::cout << "SOL-" << i << "\t"
// //        		<< cl[i].size() << "\t";
      
//       for (int j=0; j<cl[i].size(); j++)
// 	{
// 	  final_sol.best[ cl[i][j] ] = i ; 
// 	  //	  std::cout << " " << cl[i][j] ;
// 	}
//       //      std::cout << "\n";
//     }
  
  //  std::cout << "\n";
  
  
  //   std::cout << "\n";

//   // final
  
//   for (int j=0; j<sol.size(); j++)
//     {
//       // Display...
//       std::cout << j << " | " ;
      
//       for (int i=0; i<sol[0].size(); i++)
// 	std::cout << sol[j][i]+1 << " ";
	      
//       std::cout  << "\n";
//     }
//   std::cout << "\n";
  
  return final_sol;
}


double cluster_t::cldist( const Data::Matrix<double> & D ,			  
			  std::vector<int> & a, 
			  std::vector<int> & b)
{
  double l;

  l = a[0]>b[0] ? D[a[0]][b[0]] : D[b[0]][a[0]]; 
  
  for (int i=0; i<a.size(); i++)
    for (int j=0; j<b.size(); j++)
      {
	if ( a[i] > b[j] )
	  {
	    if ( D[a[i]][b[j]] > l ) l = D[a[i]][b[j]];
	  }
	else
	  {
	    if ( D[b[j]][a[i]] > l ) l = D[b[j]][a[i]];	    
	  }

      }
  return l;
}


double cluster_t::groupAvgLink( const Data::Matrix<double> & D , 
				std::vector<int> & a, 
				std::vector<int> & b)
{
  
  double s = 0;
  
  for (int i=0; i<a.size(); i++)
    for (int j=0; j<b.size(); j++)
      {
	
	if ( a[i] > b[j] )
	  {
	    s += D[a[i]][b[j]];
	  }
	else
	  {
	    s += D[b[j]][a[i]];	    
	  }
	
      }
  
  return 1.0 / ( a.size() * b.size() ) *  s ;

}






