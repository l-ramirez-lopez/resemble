#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//' @title Moving/rolling correlation distance of two matrices
//' @description Computes a moving window correlation distance between two data matrices
//' @usage 
//' movcorDist(X,Y,w)
//' @param X a \code{matrix}
//' @param Y a \code{matrix}
//' @param w window size (must be odd)
//' @return a \code{matrix} of correlation distance
//' @keywords internal
//' @useDynLib resemble
//' @author Leonardo Ramirez-Lopez and Antoine Stevens
// [[Rcpp::export]]
NumericMatrix movcorDist(NumericMatrix X, NumericMatrix Y,  int w){  
  arma::mat XX(X.begin(), X.nrow(), X.ncol(), false);
  arma::mat YY(Y.begin(), Y.nrow(), Y.ncol(), false);    
  arma::mat rmwF = arma::zeros(X.nrow(),Y.nrow());
  if(!w%2) throw exception("'w' shoudl be odd") ;
  int gap = (w-1)/2;
  if((Y.ncol() - w ) < 1){
    arma::mat rmwF = arma::cor(XX.t(), YY.t());
  } else {
    for(int i = gap; i < Y.ncol() - gap; i++){
      rmwF += arma::cor(XX.cols(i-gap,i+gap).t(), YY.cols(i-gap,i+gap).t()); // sum of the correlations
    }
    rmwF = rmwF/(Y.ncol()- (2*gap)); // get the average
  }
  arma::mat scmw = (1 - rmwF)/2;  // get the distance
  return wrap(scmw.t());
}
