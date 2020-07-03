#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>    // OpenMP
#endif
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

//' @title A fast distance algorithm for two matrices written in C++ 
//' @description Computes distances between two data matrices using 
//' "euclid", "cor", "cosine" 
//' @usage 
//' fast_diss(X, Y, method)
//' @param X a matrix
//' @param Y a matrix
//' @param method a \code{string} with possible values "euclid", "cor", "cosine"
//' @return a distance matrix
//' @keywords internal
//' @useDynLib resemble
//' @author Antoine Stevens and Leonardo Ramirez-Lopez
// [[Rcpp::export]]   
arma::mat fast_diss(NumericMatrix X, NumericMatrix Y, String method){  
  
  double eps = sqrt(DOUBLE_EPS);
  
  int nX = X.nrow(), kX = X.ncol(), nY = Y.nrow(), kY = Y.ncol();
  arma::mat XX(X.begin(), nX, kX, false); // reuses memory and avoids extra copy
  arma::mat YY(Y.begin(), nY, kY, false); // reuses memory and avoids extra copy 
  if(method == "euclid"){
    arma::mat output = arma::ones(nY,1) * arma::sum(arma::square(XX),1).t() + arma::sum(arma::square(YY),1)  * arma::ones(1,nX) - 2 * YY * XX.t();
    return output;
  }   
  if(method=="cor"){
    arma::mat output = (1 - arma::cor(XX.t(), YY.t()))/2;   
    return wrap(output.t());
  }
  else{ // cosine
    arma::mat numerator = XX * YY.t();
    arma::mat dvsr = arma::sqrt(arma::sum(arma::square(XX),1)) * arma::sqrt(arma::sum(arma::square(YY),1)).t();
    arma::mat output = arma::acos(numerator/dvsr);     
    return output.t();
  }   
}



//' @title A fast algorithm of (squared) Euclidean cross-distance for vectors written in C++ 
//' @description A fast (parallel for linux) algorithm of (squared) Euclidean cross-distance for vectors written in C++ 
//' @usage 
//' fast_diss_vector(X)
//' @param X a vector.
//' @return a vector of distance (lower triangle of the distance matrix, stored by column)
//' @details used internally in ortho_projection
//' @author Antoine Stevens
//' @keywords internal 
//' @useDynLib resemble
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]   
NumericVector fast_diss_vector(NumericVector X){  
  int nX = X.size();
  int n = ((nX*nX)-nX)/2;
  NumericVector output(n);  
#if defined(_OPENMP) 
  #pragma omp parallel for schedule(dynamic)
#endif
  for(int i = 0; i < nX-1; i++)
    for(int j = i+1; j < nX; j++){
      double x = X(j)-X(i);
      output(nX*i - (i * (i + 3) / 2)  + j  - 1) =  x * x; 
    }
    return output;         
}

//' @title A fast (serial) algorithm of Euclidean (non-squared) cross-distance for vectors written in C++ 
//' @description A fast algorithm of (squared) Euclidean cross-distance for vectors written in C++ 
//' @usage 
//' fastDistVVL(X)
//' @param X a vector
//' @return a vector of distance (lower triangle of the distance matrix, stored by column)
//' @details used internally in ortho_projection
//' @author Leo Ramirez-Lopez
//' @keywords internal 
//' @useDynLib resemble
// [[Rcpp::export]] 
NumericVector fastDistVVL(NumericVector X){
  int nX = X.size();
  // Compute the length of the output vector
  // that will contain the lower tringle
  // of the distance matrix...
  int n = ((nX * nX) - nX) / 2;
  // ... and create the vector 
  NumericVector output(n);
  
  int starti = 0;
  int endi = 0;  
  for(int i = 0; i < (nX - 1); i++){
    int ii = i + 1;
    int length_vector;
    starti = (i * (nX - 1)) - (((i * i) - i) / 2);
    endi = ((ii * (nX - 1)) - (((ii * ii) - ii) / 2)) - 1;
    //  create a vector or indices 
    Range indices_1(starti, endi);
    Range indices_2(ii, nX - 1);
    length_vector = indices_1.size();
    //  Create a vector of length length_vector with 
    //  repeated X[i] values
    NumericVector ith_point(length_vector, X[i]);
    output[indices_1] = ith_point - X[indices_2];
  }
  
  // compute power of 2
  // convert to absolute differences (equivalent to euclidean)
  output = pow(output, 2);
  return output;
}


//' @title A function to compute indices of minimum values of a distance vector
//' @description For internal use only
//' @usage 
//' minDissV(X,omp_threads)
//' @param X a vector of distance (as computed in \code{resemble:::fastDistVV} or \code{base::dist}).
//' @return a vector of the indices of the nearest neighbors
//' @details 
//' Used internally to find the nearest neighbors. 
//' It searches in lower (or upper?) trianguular matrix. Therefore this must be the format of the 
//' input data. The piece of code int \code{len = (sqrt(X.size()*8+1)+1)/2} generated an error in CRAN
//' since \code{sqrt} cannot be applied to integers.
//' @keywords internal
//' @useDynLib resemble
//' @author Antoine Stevens 
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]  
NumericVector minDissV(NumericVector X){  
  // For the distances
  int nX = X.size();
  int n = ((nX * nX) - nX) / 2;
  NumericVector doutput(n);   
  
  // For the indices
  arma::uword  index;
  double vct = (sqrt(((double)doutput.size()) * 8.0 + 1.0) + 1.0) / 2.0;
  int len = (int)vct;
  // int len = (sqrt(Ds.size()*8+1)+1)/2;
  arma::uvec vindex(len); 
  int i,j;
#if defined(_OPENMP) 
  #pragma omp parallel for private(i,j) schedule(dynamic)
#endif
  for(int i = 0; i < nX - 1; i++)
    for(int j = i + 1; j < nX; j++){
      double x = X(j) - X(i);
      doutput(nX*i - (i * (i + 3) / 2)  + j  - 1) =  x * x;          
    }
    
    for(i = 0; i < len; i++){       
      arma::vec x(len); 
      for(j = 0; j < i; j++){
        // triangular sequence
        int k = j * len - (j * (j + 3) / 2) + i - 1;
        x[j] = doutput(k);        
      }        
      for(j = i + 1; j < len; j++){
        // triangular sequence
        int k2 = i * len - (i * (i + 3) / 2) + j - 1;
        x[j] = doutput(k2);             
      }
      x[i] = arma::datum::nan; // remove diag
      x.min(index); // don't assign result to a value since we are interested only in the index
      vindex[i] = index;    
    }
    return wrap(vindex + 1);
}
