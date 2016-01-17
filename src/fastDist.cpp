#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>    // OpenMP
#endif
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;


//' @title A function to compute indices of minimum values of a distance vector
//' @description For internal use only
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
    x.min(index); // don't assign result to a value since we are interested only in the index
    vindex[i] = index;    
  }
  return wrap(vindex +1);
}


//' @title A fast (parallel for linux) algorithm of (squared) Euclidean cross-distance for vectors written in C++ 
//' @description A fast (parallel for linux) algorithm of (squared) Euclidean cross-distance for vectors written in C++ 
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

//' @title A fast (serial) algorithm of (squared) Euclidean cross-distance for vectors written in C++ 
//' @description A fast (parallel for linux) algorithm of (squared) Euclidean cross-distance for vectors written in C++ 
//' @usage 
//' fastDistVVL(X)
//' @param X a \code{vector}
//' @return a \code{vector} of distance (lower triangle of the distance matrix, stored by column)
//' @details used internally in orthoProjection
//' @author Leo Ramirez-Lopez
//' @keywords internal 
//' @useDynLib resemble
//' 
// [[Rcpp::export]] 
NumericVector fastDistVVL(NumericVector X){  
  int nX = X.size();
  int n = ((nX*nX)-nX)/2;
  NumericVector output(n);
  int lst = nX;
  int ed = 0;
  int st;
  for(int i = 0; i < nX - 1; i++){
    int ii = i + 1;
    st = (i * (nX - 1)) - (((i*i)-i)/2); // remove exesive ()
    ed = (ii * (nX - 1)) - (((ii*ii)-ii)/2); // remove exesive ()
    Range idx(st,ed);
    Range idx2(i+1,nX-1);
    lst = idx.size();
    NumericVector xx(lst, X[i]);
    output[idx] = abs(xx - X[idx2]);
    //output[idx] = xx - X[idx2];
  }
  //output = abs(output);
  return output;         
}



//' @title A function to compute indices of minimum values of a distance vector
//' @description For internal use only
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
NumericVector minDissV(NumericVector X,int cores){  
#ifdef _OPENMP
  omp_set_num_threads(cores);
#endif
  // For the distances
  int nX = X.size();
  int n = ((nX*nX)-nX)/2;
  NumericVector doutput(n);   
  
  
  // For the indices
  arma::uword  index;
  double vct = (sqrt(((double)doutput.size())*8.0+1.0)+1.0)/2.0;
  int len = (int)vct;
  // int len = (sqrt(Ds.size()*8+1)+1)/2;
  arma::uvec vindex(len);    
  int i,j;
  
#pragma omp parallel for private(i,j) schedule(dynamic)
  for(int i = 0; i < nX-1; i++)
    for(int j = i+1; j < nX; j++){
      double x = X(j)-X(i);
      doutput(nX*i - (i*(i+3)/2)  + j  - 1) =  x*x;          
    }
    
    for(i = 0; i < len; i++){       
      arma::vec x(len); 
      for(j = 0; j < i; j++){
        // triangular sequence
        int k = j*len - (j*(j+3)/2) + i - 1;
        x[j] = doutput(k);        
      }        
      for(j = i+1; j < len; j++){
        // triangular sequence
        int k2 = i*len - (i*(i+3)/2) + j - 1;
        x[j] = doutput(k2);             
      }
      x[i] = arma::datum::nan; // remove diag
      x.min(index); // don't assign result to a value since we are interested only in the index
      vindex[i] = index;    
    }
    return wrap(vindex +1);
}

