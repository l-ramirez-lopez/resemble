#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>    // OpenMP
#endif

using namespace Rcpp;
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]

//' @title A function to compute indices of minimum values of a distance vector
//' @description For internal use only
//' @usage 
//' which_min_vector(X)
//' @param X a vector of distances 
//' @return a vector of the indices of the nearest neighbors
//' @details 
//' Used internally to find the nearest neighbors. 
//' It searches in lower (or upper) triangular matrix. Therefore this must be the format of the 
//' input data. The piece of code int \code{len = (sqrt(X.size()*8+1)+1)/2} generated an error in CRAN
//' since \code{sqrt} cannot be applied to integers.
//' @keywords internal
//' @useDynLib resemble
//' @author Antoine Stevens 
// [[Rcpp::export]]  
NumericVector which_min_vector(NumericVector X){  
  arma::uword  index;
  double vct = (sqrt(((double)X.size()) * 8.0 + 1.0) + 1.0) / 2.0;
  int len = (int)vct;
  // int len = (sqrt(X.size()*8+1)+1)/2;
  arma::uvec vindex(len);    
  int i,j;
#if defined(_OPENMP) 
  #pragma omp parallel for private(i,j) schedule(dynamic)
#endif
  for(i = 0; i < len; i++){       
    arma::vec x(len); 
    for(j = 0; j < i; j++){
      // triangular sequence
      int k = j * len - (j * (j + 3) / 2) + i - 1;
      x[j] = X(k);        
    }        
    for(j = i+1; j < len; j++){
      // triangular sequence
      int k2 = i * len - (i * (i + 3) / 2) + j - 1;
      x[j] = X(k2);             
    }
    x[i] = arma::datum::nan; // remove diag
    x.min(index); // don't assign result to a value since we are interested only in the index
    vindex[i] = index;    
  }
  return wrap(vindex +1);
}
