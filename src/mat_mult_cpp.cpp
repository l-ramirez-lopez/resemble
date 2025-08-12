#include <Rcpp.h>
extern "C" {
#include <R_ext/BLAS.h>   // R's portable BLAS shim
}
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix mat_mult_cpp(const NumericMatrix& A, const NumericMatrix& B) {
  // dimensions
  int m  = A.nrow();
  int k  = A.ncol();
  int kB = B.nrow();
  int n  = B.ncol();
  if (k != kB) stop("non-conformable: ncol(A) != nrow(B)");
  
  NumericMatrix C(m, n);
  

  char transA = 'N', transB = 'N';
  double alpha = 1.0, beta = 0.0;
  int lda = m, ldb = k, ldc = m;  
  
  F77_CALL(dgemm)(
      &transA, &transB,
      &m, &n, &k,
      &alpha,
      A.begin(), &lda,   
      B.begin(), &ldb,
      &beta,
      C.begin(), &ldc
  FCONE FCONE        
  );
  
  return C;
}
