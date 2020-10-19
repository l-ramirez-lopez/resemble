#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace RcppArmadillo;


//' @title Moving/rolling correlation distance of two matrices
//' @description Computes a moving window correlation distance between two data matrices
//' @usage 
//' moving_cor_diss(X,Y,w)
//' @param X a matrix
//' @param Y a matrix
//' @param w window size (must be odd)
//' @return a matrix of correlation distance
//' @keywords internal
//' @useDynLib resemble
//' @author Leonardo Ramirez-Lopez and Antoine Stevens
// [[Rcpp::export]]
NumericMatrix moving_cor_diss(arma::mat X, arma::mat Y, int w){  
  arma::mat rmwF = arma::zeros(X.n_rows, Y.n_rows);
  if(!w%2) {
    throw std::invalid_argument("w must be odd");
  }
  int gap = (w - 1) / 2;
  if((Y.n_cols - w ) < 1){
    arma::mat rmwF = arma::cor(trans(X), trans(Y));
  } else {
    int ny = Y.n_cols;
    for(int i = gap; i < ny - gap; i++){
      // sum of the correlations
      rmwF += arma::cor(X.cols(i-gap, i + gap).t(), Y.cols(i - gap, i + gap).t()); 
    }
    // get the average
    rmwF = rmwF / (Y.n_cols - (2*gap)); 
  }
  // get the distance
  arma::mat scmw = (1 - rmwF)/2;  
  return wrap(scmw.t());
}
