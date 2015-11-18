#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>    // OpenMP
#endif
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

//' @title A fast distance algorithm for two matrices written in C++ 
//' @description Computes distances between two data matrices using "euclid", "cor", "cosine" 
//' @usage 
//' fastDist(X,Y,method)
//' @param X a \code{matrix}
//' @param Y a \code{matrix}
//' @param method a \code{string} with possible values "euclid", "cor", "cosine"
//' @return a distance \code{matrix}
//' @keywords internal
//' @useDynLib resemble
//' @author Antoine Stevens and Leonardo Ramirez-Lopez
// [[Rcpp::export]]   
arma::mat fastDist(NumericMatrix X, NumericMatrix Y, String method){  
   int nX = X.nrow(), kX = X.ncol(), nY = Y.nrow(), kY = Y.ncol();
   arma::mat XX(X.begin(), nX, kX, false); // reuses memory and avoids extra copy
   arma::mat YY(Y.begin(), nY, kY, false); // reuses memory and avoids extra copy 
   if(method=="euclid"){
     arma::mat output = arma::ones(nY,1) * arma::sum(arma::square(XX),1).t() + arma::sum(arma::square(YY),1)  * arma::ones(1,nX) - 2 * YY * XX.t();
     return output;
   }   
   if(method=="cor"){
     arma::mat output = (1 - arma::cor(XX.t(), YY.t()))/2;   
     return output.t();
   }
   else{ // cosine
      arma::mat numerator = XX * YY.t();
      arma::mat dvsr = arma::sqrt(arma::sum(arma::square(XX),1)) * arma::sqrt(arma::sum(arma::square(YY),1)).t();
      arma::mat output = arma::acos(numerator/dvsr);     
      return output.t();
   }   
}

//' @title A fast distance algorithm for a matrix and a vector written in C++
//' @description A fast distance algorithm for a matrix and a vector written in C++  
//' @usage 
//' fastDistV(X,Y,method)
//' @param X a \code{matrix}
//' @param Y a \code{vector}
//' @param method a \code{string} with possible values "euclid", "cor", "cosine"
//' @return a distance \code{vector}
//' @author Antoine Stevens and Leonardo Ramirez-Lopez
//' @keywords internal 
//' @useDynLib resemble
// [[Rcpp::export]]   
NumericVector fastDistV(NumericMatrix X, NumericVector Y, String method){  
   int nX = X.nrow(), kX = X.ncol(),  kY = Y.size();
   arma::mat XX(X.begin(), nX, kX, false); // reuses memory and avoids extra copy
   arma::rowvec YY(Y.begin(), kY, false); // reuses memory and avoids extra copy 
   typedef std::vector<double> stdvec;
   if(method=="euclid"){
     stdvec output = arma::conv_to<stdvec>::from(arma::sum(arma::square(XX),1).t() + arma::sum(arma::square(YY))  * arma::ones(1,nX) - 2 * YY * XX.t());
     return wrap(output);
   }   
   if(method=="cor"){
    stdvec output = arma::conv_to<stdvec>::from((1 - arma::cor(XX.t(), YY.t()))/2);   
     return wrap(output);
   }
   else{ //cosine
      arma::mat numerator = XX * YY.t();
      arma::mat dvsr = arma::sqrt(arma::sum(arma::square(XX),1)) * arma::sqrt(arma::sum(arma::square(YY),1)).t();
      stdvec output = arma::conv_to<stdvec>::from(arma::acos(numerator/dvsr));     
      return wrap(output);
   }   
}

//' @title A fast (parallel) algorithm of (squared) Euclidean cross-distance for vectors written in C++ 
//' @description A fast (parallel) algorithm of (squared) Euclidean cross-distance for vectors written in C++ 
//' @usage 
//' fastDistVV(X, cores)
//' @param X a \code{vector}
//' @return a \code{vector} of distance (lower triangle of the distance matrix, stored by column)
//' @details used internally in orthoProjection
//' @author Antoine Stevens
//' @keywords internal 
//' @useDynLib resemble
// [[Rcpp::export]]   
NumericVector fastDistVV(NumericVector X, int cores){  
   #ifdef _OPENMP
   omp_set_num_threads(cores);
   #endif
   int nX = X.size();
   int n = ((nX*nX)-nX)/2;
   NumericVector output(n);   
   #pragma omp parallel for schedule(dynamic)
   for(int i = 0; i < nX-1; i++)
     for(int j = i+1; j < nX; j++){
     double x = X(j)-X(i);
     output(nX*i - (i*(i+3)/2)  + j  - 1) =  x*x;          
   }
   return output;         
}