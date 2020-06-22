
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


#include "glm.h"
#include "statistics.h"

#include <cmath>
#include <iostream>

#include "helper/logger.h"

extern logger_t logger;

//
// Pre-processing steps (standardization; checking for multi-collinearity)
//

void GLM::ci(double ci)
{
  ci_zt = Statistics::ltqnorm( 1 - (1 - ci ) / 2.0  );
}

void GLM::vif(double v)
{
  vif_threshold = v;
}

bool GLM::valid() const
{
  return all_valid;
}

void GLM::set( Data::Vector<double> & y , 
	       Data::Matrix<double> & x , 
	       std::vector<int> * cl , 
	       std::vector<bool> * mask )
{

  // need to go over all individuals first, as some may be filtered
  // ultimately, put the filter inthe matrix class, so that we can
  // ignore this
  
  int n1 = mask ? mask->size() : y.size();
  for (int i=0; i < n1; i++)
    {
      if ( (!mask) || (*mask)[i] )
	{
	  if ( model == LOGISTIC )
	    Y.push_back( y[i] ? 1 : 0 );
	  else 
	    Y.push_back( y[i] );
	  Data::Vector<double> r = x.row(i);	  
	  X.add_row( x.row(i) );
	  if ( cl ) clst.push_back( (*cl)[i] );
	}
    }


  nind = Y.size();  
  np = x.dim2();
  
  if ( model == LOGISTIC ) 
    {
      pr.resize(nind,0);
      V.resize(nind,0);
    }
}


void GLM::set_variance()
{
  varY=0;
  meanY = 0;
  int actualN=0;

  for (int i=0; i<nind; i++)
    {
      actualN++;
      meanY += Y[i];
    }

  if (actualN==0)
    {
      varY=0;
      return;
    }
  
  meanY /= (double)actualN;
  
  for (int i=0; i<nind; i++)
    varY += (Y[i] - meanY) * (Y[i] - meanY);      
  varY /= (double)(actualN-1);

  if ( actualN != nind ) 
    Helper::halt("internal error in GLM()");
}

void GLM::standardise()
{
  
  // Get mean and variance for all variable
  double sdY = sqrt(varY);
  for (int i=0; i<nind; i++)
    Y[i] = ( Y[i] - meanY ) / sdY;

  Data::Vector<double> mX(np);
  Data::Vector<double> sdX(np);

  // Standardise all predictors, except intercept
  for (int i=0; i<nind; i++)
    for (int j=1; j<np; j++)
      mX[j] += X(i,j);
  for (int j=1; j<np; j++)
    mX[j] /= nind;
  for (int i=0; i<nind; i++)
    for (int j=1; j<np; j++)
      sdX[j] += ( X(i,j) - mX[j] ) * ( X(i,j) - mX[j] );
  for (int j=1; j<np; j++)
    {
      sdX[j] = sqrt( sdX[j] / (nind-1) );
      if ( sdX[j] == 0 ) sdX[j] = 1;
    }
  for (int i=0; i<nind; i++)
    for (int j=1; j<np; j++)
      X(i,j) = ( X(i,j) - mX[j] ) / sdX[j];
}



bool GLM::check_VIF()
{
  
  valid( false );

  // Calculate correlation matrix for X
  // Skip intercept (first term, X[,0]


  int p = nind;
  int q = np - 1 ;

  if ( p < 2 || q < 2 ) 
    {
      valid(true);
      return true;
    }
  
  Data::Vector<double> m(q);
  Data::Matrix<double> c(q,q);

  for (int i=0; i<p; i++)
    for (int j=0; j<q; j++)
      m(j) += X(i,j+1);

  for (int j=0; j<q; j++)
    m(j) /= (double)p;
  
  for (int i=0; i<p; i++)
    for (int j1=0; j1<q; j1++)
      for (int j2=j1; j2<q; j2++)
	c(j1,j2) += ( X(i,j1+1) - m(j1) ) * ( X(i,j2+1) - m(j2) );  
  
  for (int j1=0; j1<q; j1++)
    for (int j2=j1; j2<q; j2++)
      c(j1,j2) /= (double)(p-1);

  for (int j1=0; j1<q; j1++)
    for (int j2=j1+1; j2<q; j2++)
      {
	c(j1,j2) /= sqrt( c(j1,j1) * c(j2,j2) );
	c(j2,j1) = c(j1,j2);	
	if ( c(j2,j1) > 0.999 ) return false;	
      }         

  // Any item with zero variance?  
  for (int j=0; j<q; j++)
    {
      if ( c(j,j) == 0 || ! Helper::realnum( c(j,j) ) ) return false;
      c(j,j) = 1;
    }

  bool flag = true;
  c = Statistics::inverse( c , &flag);  
  if ( ! flag ) all_valid = false;
    
  double maxVIF = 0;
  for (int j=0;j<q;j++)
    {      
      // r^2 = 1 - 1/x where x is diagonal element of inverted
      // correlation matrix
      // As VIF = 1 / ( 1 - r^2 ) , implies VIF = x      
      if ( c(j,j) > vif_threshold ) return false;
    }

  valid(true);
  return true;
}



bool GLM::fit_linear()
{
  
  // if only intercept + term, use a quicker
  // version

  if ( np == 2 && ! cluster ) return fit_univariate_linear();
  
  all_valid = true;  // replace with checkVIF etc
  
  if ( np == 0 || nind == 0 || ! all_valid )
    {
      all_valid = false;
      return false;
    }


  // Allocate space
 
  coef.resize( np );  
  S.resize( np, np );
  
  // DV 

  set_variance();  
  if ( standard_beta ) standardise();

  sig.resize( nind , sqrt( 1.0 / sqrt( (double)nind) ) );
  
  w.resize( np );
  u.resize( nind , np );
  v.resize( np , np );
   
  //  Perform "svdfit(C,Y,sig,b,u,v,w,chisq,function)"
  
  const double TOL=1.0e-13;
  
  double tmp,thresh,sum;

  Data::Vector<double> b(nind);

  for (int i=0; i < nind; i++) 
    {
      Data::Vector<double> afunc = X.row(i);
      tmp=1.0/sig[i];
      for (int j=0; j<np; j++) u(i,j) = afunc[j] * tmp;
      b[i] = Y[i] * tmp;
    }
  
  if ( ! Statistics::svdcmp(u,w,v) ) 
    {
      all_valid = false;
      return false;
    }
  
  double wmax = 0.0;
  for (int j=0;j<np;j++)
    if (w[j] > wmax) wmax = w[j];
  thresh=TOL*wmax;
  for (int j=0;j<np;j++)
    if (w[j] < thresh) w[j]=0.0;

  Statistics::svbksb(u,w,v,b,coef);
  

//   chisq = 0.0;
//   for (i=0;i<nind;i++) 
//     {
//       afunc = X[i];
//       sum=0.0;
//       for (int j=0; j<np; j++) sum += coef[j] * afunc[j];
//       chisq += ( tmp = (Y(i) - sum) / sig(i), tmp*tmp );
//     }


  /////////////////////////////////////////
  // Obtain covariance matrix of estimates

  // Robust cluster variance estimator
  // V_cluster = (X'X)^-1 * \sum_{j=1}^{n_C} u_{j}' * u_j * (X'X)^-1 
  // where u_j = \sum_j cluster e_i * x_i 

  // Above, e_i is the residual for the ith observation and x_i is a
  // row vector of predictors including the constant.

  // For simplicity, I omitted the multipliers (which are close to 1)
  // from the formulas for Vrob and Vclusters.

  // The formula for the clustered estimator is simply that of the
  // robust (unclustered) estimator with the individual ei*s replaced
  // by their sums over each cluster. 

  // http://www.stata.com/support/faqs/stat/cluster.html
  // SEE http://aje.oxfordjournals.org/cgi/content/full/kwm223v1#APP1

  // Williams, R. L. 2000.  A note on robust variance estimation for
  // cluster-correlated data. Biometrics 56: 64

  //  t ( y - yhat X  ) %*%  ( y - yhat)  / nind - np
  // = variance of residuals 
  // j <- ( t( y- m %*% t(b) ) %*% ( y - m %*% t(b) ) ) / ( N - p ) 
  // print( sqrt(kronecker( solve( t(m) %*% m ) , j )  ))
  

  ////////////////////////////////////////////////
  // OLS variance estimator = s^2 * ( X'X )^-1
  // where s^2 = (1/(N-k)) \sum_i=1^N e_i^2
  

  // 1. Calcuate S = (X'X)^-1
  
  bool okay = true;

  Data::Matrix<double> S0 = Statistics::inverse( Statistics::transpose( X ) * X , & okay );
  
  if ( ! okay ) 
    {
      all_valid = false;
      return false;
    }

  
  ////////////////////////
  // Calculate s^2 (sigma)

  if ( ! cluster )
    {

      double sigma= 0.0;

      for (int i=0; i < nind; i++)
        {
          double partial = 0.0;
          for (int j=0; j<np; j++)
            partial += coef(j) * X(i,j);
          partial -= Y[i];
          sigma += partial * partial;
        }
      sigma /= nind-np; 
      
      for (int i=0; i<np; i++)
        for (int j=0; j<np; j++)
          S(i,j) = S0(i,j) * sigma;
    }
  

  ///////////////////////////
  // Robust-cluster variance

  if ( cluster )
    {
    
      Data::Matrix<double> sc( nc , np );
      
      for (int i=0; i<nind; i++)
	{
	  double partial = 0.0;
	  for (int j=0; j<np; j++)
	    partial += coef[j] * X(i,j);
	  partial -= Y[i];
	  
	  for (int j=0; j<np; j++)
	    sc( clst[i] , j ) += partial * X(i,j);
	}
    
      Data::Matrix<double> meat(np,np);

      for (int k=0; k < nc; k++)
	{      
	  for (int i=0; i<np; i++)
	    for (int j=0; j<np; j++)
	      meat[i][j] += sc[k][i] * sc[k][j];       
	}
      
      // make the sandwich
      S = S0 * meat * S0;
      
    }
  return true;
}


bool GLM::fit_univariate_linear() 
{

  // Speed-up version for univariate case has set coef and S

  if ( np != 2 || nind == 0 )
    {
      all_valid = false;
      return false;
    }
  
  coef.resize(2);
  S.resize(2,2);

  double x_mean=0, x_var=0;
  double y_mean=0, y_var=0;
  double y_x_covar=0;
  
  for (int i=0; i<nind; i++)
    {
      y_mean += Y(i);
      x_mean += X(i,1);
    }

  x_mean /= (double)nind;
  y_mean /= (double)nind;      
        
  for (int i=0; i<nind; i++)
    {
      double ty = Y(i) - y_mean;
      double tx = X(i,1) - x_mean;
      y_var += ty*ty;
      x_var += tx*tx;
      y_x_covar += tx * ty;
    }

  y_var /= (double)nind - 1;
  x_var /= (double)nind - 1;
  y_x_covar /= (double)nind - 1;
  
  // b_1 term and SE
  coef(1) = y_x_covar / x_var;
  S(1,1) = (y_var/x_var - (y_x_covar*y_x_covar)/(x_var*x_var) ) / (double)( nind - 2 ) ;  
  
  // also populate for the intercept
  double sse = 0;
  for (int i=0;i<nind;i++)
    {
      double e = Y(i) - coef(1) * X(i,1) ;
      sse += e * e;
    }
  coef(0) = y_mean - coef(1) * x_mean;
  S(0,0) = sqrt(sse/(double)(nind-2.0)) * sqrt( (1/(double)nind ) + ( ( x_mean*x_mean ) / ( S(1,1) ) ) ) ;

  return true;
}


Data::Vector<double> GLM::get_var()
{  
  double multiplier = 1;  
  // Ignore... ??needs revisiting??
  // if (cluster) multiplier = (double)(nc)/((double)(nc-np));
  
  Data::Vector<double> var(np);
  for (int i=0; i<np; i++) 
    var(i) = multiplier * S(i,i);
  return var;
}


Data::Vector<double> GLM::get_SE()
{
  double multiplier = 1;  
  //  if (cluster) multiplier = (double)(nc)/((double)(nc-np));
  Data::Vector<double> var( np );
  for (int i=0; i<np; i++) 
    var[i] = sqrt ( multiplier * S(i,i) ) ;  
  return var;
}


std::string GLM::summary() 
{

  std::vector<bool> mask;
  Data::Vector<double> beta;
  Data::Vector<double> se;
  Data::Vector<double> lowci;
  Data::Vector<double> uprci;
  Data::Vector<double> statistic;
  Data::Vector<double> pvalue;
  display( &beta, &se, &pvalue , &mask, &lowci, &uprci, &statistic );
  std::stringstream ss;
  for (int i=0; i< mask.size(); i++)
    {
      if ( ! mask[i] ) 
	ss << "NA\tNA\tNA\tNA\tNA\tNA\n";
      else
	{
	  ss << beta[i] << "\t"
	     << se[i] << "\t"
	     << lowci[i] << "\t"
	     << uprci[i] << "\t"
	     << statistic[i] << "\t"
	     << pvalue[i] << "\n";      
	}
    }
  return ss.str();
}

bool GLM::test_valid() const
{
  return test_var() < 1e-20 || ! Helper::realnum( test_var() ) ? false : all_valid;
}

double GLM::test_var() const
{
  // Note:: assume 'multiplier' of 1.0 (which should be fine, as not using this 
  // for HW-estimates now, I recall
  return S(t,t);
}

double GLM::test_coef() const
{
  
  return all_valid ? ( model == LINEAR ? coef[t] : exp( coef[t] ) ) : 0 ;
}

double GLM::test_se() const
{
  return sqrt( test_var() );
}

double GLM::test_pval() const
{
  return all_valid ? ( model == LINEAR ? 
		       Statistics::t_prob( test_statistic() , Y.size()-np ) :
		       Statistics::chi2_prob( Statistics::SQR( test_statistic() ) , 1 ) ) : 1.0 ;
}

double GLM::test_statistic() const
{  
  return all_valid ? coef[t] / test_se() : 0 ;
}

double GLM::test_lower_ci() const
{
  return all_valid ? ( model == LINEAR ? 
		       coef[t] - ci_zt * test_se() :
		       exp( coef[t] - ci_zt * test_se() ) ) : 0 ; 
}

double GLM::test_upper_ci() const
{
  return all_valid ? ( model == LINEAR ? 
		       coef[t] + ci_zt * test_se() :
		       exp( coef[t] + ci_zt * test_se() ) ) : 0 ;
}


bool GLM::display( Data::Vector<double> * beta , 
		   Data::Vector<double> * se , 
		   Data::Vector<double> * pvalue ,
		   std::vector<bool> * mask , 
		   Data::Vector<double> * lowci , 
		   Data::Vector<double> * uprci , 
		   Data::Vector<double> * statistic )
{
  
  Data::Vector<double> var = all_valid ? get_var() : Data::Vector<double>(np) ; 
  
  if ( mask ) mask->resize( np );
  if ( beta ) beta->resize( np );
  if ( se ) se->resize( np );
  if ( lowci ) lowci->resize( np );
  if ( uprci ) uprci->resize( np );
  if ( statistic ) statistic->resize( np );
  if ( pvalue ) pvalue->resize( np );
  
  bool all_okay = true;
  
  for (int p = 0; p < np ; p++) 
    {

      bool okay = var[p] < 1e-20 || ! Helper::realnum( var[p] ) ? false : all_valid;
      
      if ( mask ) (*mask)[p] = okay;
      
      if ( okay )
        {	  
	  double sderr = sqrt( var[p] );

	  double Z = coef[p] / sderr;
	 
	  if ( se ) (*se)(p) = sderr;  
	  if ( statistic ) (*statistic)(p) = Z;
	  
	  if ( model == LINEAR ) 
	    {
	      if ( beta ) (*beta)(p) = coef[p];
	      if ( lowci ) (*lowci)(p) = coef[p] - ci_zt * sderr;
	      if ( uprci ) (*uprci)(p) = coef[p] + ci_zt * sderr;
	      if ( pvalue ) (*pvalue)(p) = Statistics::t_prob( Z, Y.size()-np ) ;
	    }
	  else
	    {
	      if ( beta ) (*beta)(p) = exp( coef[p] );
	      if ( lowci ) (*lowci)(p) = exp( coef[p] - ci_zt * sderr );
	      if ( uprci ) (*uprci)(p) = exp( coef[p] + ci_zt * sderr );
	      if ( pvalue ) (*pvalue)(p) = Statistics::chi2_prob( Z*Z , 1 );
	    }
        }
      else
	all_okay = false;
    }

  return all_okay;

}


double GLM::calc_RSS()
{
  if ( model == LOGISTIC ) return 0;
  
  // Calculate residual sum of squares (RSS)  

  // Might already be calculated?  
  if ( RSS >= 0 ) return RSS;

  // std::cerr << nind << " " << np << "\n";
  // std::cerr << "coef " << coef.size() << " " << X.dim1() << " " << X.dim2() << "\n";
    
  RSS = 0;  
  for (int i=0; i<nind; i++)
    {
      double t = Y[i];      
      for ( int p=0; p<np; p++)
        t -= coef[p] * X(i,p);      
      t *= t;
      RSS += t; 
    }  
  return RSS;
}


double GLM::calc_rsqr()
{
  if ( model == LOGISTIC ) return -1;
  
  // Return coefficient of determination. If not already calculated,
  // get residual sum of squares first (set to -1)
  
  if ( RSS < 0 ) RSS = calc_RSS();
  double SSy = varY * (nind-1);
  double r = ( SSy - RSS ) / SSy;
  return r > 0 ? ( r > 1 ? 1 : r ) : 0;

}

double GLM::calc_adj_rsqr()
{
  if ( model == LOGISTIC ) return -1;
  double ra =  1 - ( (double)(nind-1)/(double)(nind-np-1) ) 
    * ( 1 - calc_rsqr() );  
  return ra > 0 ? ( ra > 1 ? 1 : ra ) : 0; 
}

double GLM::calc_MallowC( GLM * submodel )
{  
  if ( model == LOGISTIC ) return -1;
  // Mallow's C = RSSm / S^2 + 2(m+1)-n
  // where S^2 = RSSk / (n-k-1);
  
  double Sk = calc_RSS() / ( nind - np - 1); 
  return ( submodel->calc_RSS() / Sk ) + 2 * ( submodel->np+1 )-nind;
}


double GLM::calc_FTest(GLM * submodel)
{
  double RSSk = calc_RSS();
  double RSSm = submodel->calc_RSS();  
  return ( ( RSSm - RSSk ) / (double)( np - submodel->np ) )
    / ( RSSk / (double)(nind - np - 1 ) );
}


//
// Logistic functions
//

bool GLM::fit_logistic() 
{

  coef.resize(np);
  S.resize(np,np);
  
  if ( np==0 || nind==0 || ! all_valid )
    { return false; }
  
//   for( int i=0;i<nind; i++)
//     {
//       std::cout << "I" << i << "\t" << Y(i);
//       for (int j=0;j<np; j++) std::cout << "\t" <<  X(i,j) ;
//       std::cout << "\n";
//     }
//   std::cout << "\n";

  // Newton-Raphson to fit logistic model
    
  bool converge = false;
  int it = 0;

  while ( ! converge && it < 20 ) 
    {
        
      // Determine p and V
      for (int i=0; i<nind; i++)
	{
	  double t = 0;
	  for (int j=0; j<np; j++)
	    t += coef[j] * X(i,j);
	  pr[i] = 1/(1+exp(-t));
	  V[i] = pr[i] * (1-pr[i]);
	}
        
      // Update coefficients
      // b <- b +  solve( t(X) %*% V %*% X ) %*% t(X) %*% ( y - p ) 
      
      Data::Matrix<double> T(np,np);

      for (int j=0; j<np; j++)
	for (int k=j; k<np; k++) 
	  {
	    double sum = 0;
	    for (int i=0; i<nind; i++)
	      sum += X(i,j) * V(i) * X(i,k);
	    T(j,k) = T(k,j) = sum;        
	  }

      bool flag = true;
      T = Statistics::inverse( T, &flag);
      if ( ! flag ) 
	{
	  all_valid = false;
	  return false;
	}

      Data::Matrix<double> T2(np,nind);
      
      // note implicit transpose of X
      for (int i=0; i<np; i++)
	for (int j=0; j<nind; j++)
	  for (int k=0; k<np; k++)
	    T2(i,j) += T(i,k) * X(j,k);
      
      Data::Vector<double> t3(nind);
      for (int i=0; i < nind; i++) 
	t3(i) = Y(i) - pr(i);
      
      Data::Vector<double> ncoef(np);
      for (int j=0; j<np; j++) 
	for (int i=0; i<nind; i++) 
	  ncoef(j) += T2(j,i) * t3(i);

      // Update coefficients, and check for 
      // convergence
      double delta = 0;
      for (int j=0; j<np; j++)        
	{
	  delta += fabs( ncoef[j] );
	  coef(j) += ncoef(j);
	}

      //      std::cout << "delta = " << delta << "\n";

      if ( delta < 1e-6 )
	converge = true;

      // Next iteration
      it++;
    }


  /////////////////////////////////////////
  // Obtain covariance matrix of estimates
      
  // S <- solve( t(X) %*% V %*% X )    
      
  // Transpose X and multiple by diagonal V
  Data::Matrix<double> Xt(np,nind);  
  for (int i=0; i<nind; i++)
    for (int j=0; j<np; j++) 
      Xt(j,i) = X(i,j) * V[i];
  
  bool flag = true;
  S = Statistics::inverse( Xt * X , &flag );    
  if ( ! flag ) 
    {
      all_valid = false;
      return false;
    }


  if ( cluster ) HuberWhite();
  
  return true;
  
}


double GLM::get_loglik()
{
  if ( model != LOGISTIC ) return 0;
  
  // Return -2 * sample log-likelihood
  // We assume the model is fit, and all Y's are either 0 or 1

  double lnlk = 0;  
  for (int i=0; i<nind; i++)
    {
      double t = 0;
      for (int j=0; j<np; j++)
	t += coef[j] * X(i,j);
      lnlk += Y[i] == 1 ? log( 1/(1+exp(-t))) : log(1 - (1/(1+exp(-t))) );
    }  
  return -2 * lnlk; 
}


void GLM::HuberWhite()
{  
  // Calculate sandwich variance estimators, potentially allowing for
  // clustered data
  
  // Works to update the S matrix, variance/covariance matrix  
  // Originally, S will contain this, uncorrected  
  // Calcuate S = (XtX)^-1  
  
  Data::Matrix<double> S0 = S;   
  Data::Matrix<double> sc(nc,np);
  
  for (int i=0; i<nind; i++)
    {
      double err = Y[i] - pr[i];      
      for (int j=0; j<np; j++)
	sc( clst[i], j) += err * X(i,j);
    }
  
  Data::Matrix<double> meat(np,np);
  for (int k=0; k<nc; k++)
    {      
      for (int i=0; i<np; i++)
	for (int j=0; j<np; j++)
	  meat(i,j) += sc(k,i) * sc(k,j);  
    }
  S = S0 * meat * S0;  
}


double GLM::statistic() const 
{  
  if ( ! valid() ) return 0;
  // just single parameter 't' for now
  double tmp = coef[t] * coef[t];
  double tmp2 = S(t,t);
  return tmp/tmp2;
}


double GLM::linear_hypothesis( Data::Matrix<double> & H, Data::Vector<double> & h)
{
  

  // (H v - h)' (H S H')^-1 (H b - h) ~ X^2 with j df 
  // where H = constraint matrix (j x (p+1)
  //       h = null coefficient values
  //       S = estimated covariance matrix of coefficients

  //   return ( H*b - h ).transpose() 
  //     * ( H*v*H.transpose() ).inverse()
  //     * ( H*b - h ) ;
  
  // Number of constraints
  int nc = h.size(); // # of constraints
  
  // Hb-h
  Data::Vector<double> outer = ( H * coef ) - h ;
  
  // (HVH')^-1  
  bool okay = true;
  Data::Matrix<double> inner = Statistics::inverse( H * S * Statistics::transpose( H ) , &okay );
  
  if ( ! okay ) 
    {
      logger << "** problem inverting in linear_hypothesis()\n";
      valid( false );
      return 0;
    }

  // (Hb-h)'.(HVH')^-1.(Hb-h) --> X^2 nb. do not need first t()
  return Statistics::matrix_inner_product( outer * inner , outer );
    
}

