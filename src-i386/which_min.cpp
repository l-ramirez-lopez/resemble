#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>    // OpenMP
#endif

using namespace Rcpp;
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]

//' @title A function to compute row-wise index of minimum values of a square distance matrix
//' @usage 
//' which_min(X,cores)
//' @param X a square \code{matrix} of distance
//' @param cores number of cores used to run the computation
//' @return a \code{vector} of the indices of the minimum value in each row of the input \code{matrix}
//' @details Used internally to find the nearest neighbours
//' @keywords internal
//' @useDynLib resemble
//' @author Antoine Stevens 
// [[Rcpp::export]]   
NumericVector which_min(NumericMatrix X, int cores){  
   #ifdef _OPENMP
   omp_set_num_threads(cores);
   #endif
   int nX = X.nrow(), kX = X.ncol();
   arma::mat XX(X.begin(), nX, kX, false); 
   arma::uword  index;
   arma::uvec vindex(nX);   
   #pragma omp parallel for schedule(static) 
   for(int i = 0; i < nX; i++){
    arma::rowvec x = XX.row(i);
    x(i) = arma::datum::nan; // remove diag
    double z = x.min(index);
    vindex[i] = index;
   }
   return wrap(vindex +1);   
}

//' @title A function to compute indices of minimum values of a distance vector
//' @usage 
//' which_minV(X,cores)
//' @param X a \code{vector} of distance (as computed in \code{resemble:::fastDistVV} or \code{base::dist})
//' @param cores number of cores used to run the computation
//' @return a \code{vector} of the indices of the nearest neighbours
//' @details 
//' Used internally to find the nearest neighbours. 
//' It searches in lower (or upper?) trianguular matrix. Therefore this must be the format of the 
//' input data. The piece of code int \code{len = (sqrt(X.size()*8+1)+1)/2} generated an error in CRAN
//' since \code{sqrt} cannot be applied to integers.
//' @keywords internal
//' @useDynLib resemble
//' @author Antoine Stevens 
// [[Rcpp::export]]   
NumericVector which_minV(NumericVector X,int cores){  
  #ifdef _OPENMP
  omp_set_num_threads(cores);
  #endif
  arma::uword  index;
  double vct = (sqrt(((double)X.size())*8.0+1.0)+1.0)/2.0;
  int len = (int)vct;
  // int len = (sqrt(X.size()*8+1)+1)/2;
  arma::uvec vindex(len);    
  int i,j;
  #pragma omp parallel for private(i,j) schedule(dynamic)
  for(i = 0; i < len; i++){       
    arma::vec x(len); 
    for(j = 0; j < i; j++){
      // triangular sequence
      int k = j*len - (j*(j+3)/2) + i - 1;
      x[j] = X(k);        
    }        
    for(j = i+1; j < len; j++){
      // triangular sequence
      int k2 = i*len - (i*(i+3)/2) + j - 1;
      x[j] = X(k2);             
    }
    x[i] = arma::datum::nan; // remove diag
    double z = x.min(index);
    vindex[i] = index;
  }
  return wrap(vindex +1);
}