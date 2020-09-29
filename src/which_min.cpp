#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>    // OpenMP
#endif

using namespace Rcpp;
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]

//' @title A function to compute row-wise index of minimum values of a square distance matrix
//' @description For internal use only
//' @usage 
//' which_min(X)
//' @param X a square matrix of distances
//' @return a vector of the indices of the minimum value in each row of the input matrix
//' @details Used internally to find the nearest neighbors
//' @keywords internal
//' @useDynLib resemble
//' @author Antoine Stevens 
// [[Rcpp::export]]  
NumericVector which_min(NumericMatrix X){  
   int nX = X.nrow(), kX = X.ncol();
   arma::mat XX(X.begin(), nX, kX, false); 
   arma::uword  index;
   arma::uvec vindex(nX);
#if defined(_OPENMP) 
   #pragma omp parallel for schedule(static) 
#endif
   for(int i = 0; i < nX; i++){
    arma::rowvec x = XX.row(i);
    x(i) = arma::datum::nan; // remove diag
    x.min(index); // don't assign result to a value since we are interested only in the index
    vindex[i] = index;    
   }
   return wrap(vindex +1);   
}
