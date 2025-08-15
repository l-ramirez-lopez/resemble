#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <limits>
#include <vector>
#include <algorithm>
#ifdef _OPENMP
 #include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;



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
  
  //double eps = sqrt(DOUBLE_EPS);
  //FIXME check numerical precision in Rcpp
  //in some cases it returns 0s as -1e-14 
  //perhaps due to reuse memory?
  // Option: make the inputs arma::mat X arma::mat Y
  
  int nX = X.nrow(), kX = X.ncol(), nY = Y.nrow(), kY = Y.ncol();
  arma::mat XX(X.begin(), nX, kX, false); // reuses memory and avoids extra copy
  arma::mat YY(Y.begin(), nY, kY, false); // reuses memory and avoids extra copy 
  if(method == "euclid"){
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
// [[Rcpp::export]]   
NumericVector fast_diss_vector(NumericVector X) {  
  int nX = X.size();
  int n = nX * (nX - 1) / 2;
  NumericVector output(n);  
  
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int i = 0; i < nX - 1; i++) {
    for (int j = i + 1; j < nX; j++) {
      double x = X[j] - X[i];
      int idx = nX * i - (i * (i + 3) / 2) + j - 1;
      output[idx] = x * x;
    }
  }
  
  return output;         
}

//' @title Symmetric Euclidean Distance Matrix from a Single Matrix
//' @description Computes the symmetric pairwise Euclidean distance matrix
//' from the rows of a single matrix. Only the upper triangle is computed,
//' then mirrored to the lower triangle to reduce redundant operations.
//'
//' Squared distances are calculated first, then square-rooted. Tiny values
//' below 1e-14 are clamped to zero to ensure numerical stability.
//'
//' @param X An `n x p` numeric matrix where each row is a p-dimensional
//' observation.
//'
//' @return An `n x n` symmetric matrix of Euclidean distances.
//'
//' @details
//' This function is designed for efficient in-memory computation of full
//' distance matrices using only the upper triangle to avoid redundant work.
//' It applies the square root only once per distance and avoids unnecessary
//' writes to the diagonal (which is always 0). The implementation leverages
//' memory locality and avoids data duplication.
//'
//' Suitable for moderate to large `n`. If OpenMP is enabled, the upper-triangle
//' loop can optionally be parallelized for additional speedup on multicore systems.
//' Note that the function assumes that the input matrix `X` is not empty and
// [[Rcpp::export]]
Rcpp::NumericMatrix fast_self_euclid(const arma::mat& X) {
  const arma::uword n = X.n_rows;
  const arma::uword p = X.n_cols;
  const double scale = 1.0 / std::sqrt(static_cast<double>(p));
  const double eps = 1e-14;
  
  arma::mat D(n, n, arma::fill::zeros);
  
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (arma::sword i = 0; i < static_cast<arma::sword>(n - 1); ++i) {
    for (arma::uword j = i + 1; j < n; ++j) {
      const double dist2 = arma::accu(arma::square(X.row(i) - X.row(j)));
      double d = std::sqrt(dist2) * scale;
      if (d < eps) d = 0.0;
      D(i, j) = d;
    }
  }
  
  // Symmetrize: copy upper triangle to lower triangle
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (arma::sword i = 0; i < static_cast<arma::sword>(n - 1); ++i) {
    for (arma::uword j = i + 1; j < n; ++j) {
      D(j, i) = D(i, j);
    }
  }
  
  return Rcpp::wrap(D);
}

// //' @title A fast (serial) algorithm of Euclidean (non-squared) cross-distance for vectors written in C++ 
// //' @description A fast algorithm of (squared) Euclidean cross-distance for vectors written in C++ 
// //' @usage 
// //' fastDistVVL(X)
// //' @param X a vector
// //' @return a vector of distance (lower triangle of the distance matrix, stored by column)
// //' @details used internally in ortho_projection
// //' @author Leo Ramirez-Lopez
// //' @keywords internal 
// //' @useDynLib resemble
// // [[Rcpp::export]] 
// NumericVector fastDistVVL(NumericVector X){
//   int nX = X.size();
//   // Compute the length of the output vector
//   // that will contain the lower tringle
//   // of the distance matrix...
//   int n = ((nX * nX) - nX) / 2;
//   // ... and create the vector 
//   NumericVector output(n);
//   
//   int starti = 0;
//   int endi = 0;  
//   for(int i = 0; i < (nX - 1); i++){
//     int ii = i + 1;
//     int length_vector;
//     starti = (i * (nX - 1)) - (((i * i) - i) / 2);
//     endi = ((ii * (nX - 1)) - (((ii * ii) - ii) / 2)) - 1;
//     //  create a vector or indices 
//     Range indices_1(starti, endi);
//     Range indices_2(ii, nX - 1);
//     length_vector = indices_1.size();
//     //  Create a vector of length length_vector with 
//     //  repeated X[i] values
//     NumericVector ith_point(length_vector, X[i]);
//     output[indices_1] = ith_point - X[indices_2];
//   }
//   
//   // compute power of 2
//   // convert to absolute differences (equivalent to euclidean)
//   output = pow(output, 2);
//   return output;
// }


// //' @title A function to compute indices of minimum values of a distance vector
// //' @description For internal use only
// //' @usage 
// //' minDissV(X)
// //' @param X a vector of distance (as computed in \code{resemble:::fastDistVV} or \code{base::dist}).
// //' @return a vector of the indices of the nearest neighbors
// //' @details 
// //' Used internally to find the nearest neighbors. 
// //' It searches in lower (or upper?) trianguular matrix. Therefore this must be the format of the 
// //' input data. The piece of code int \code{len = (sqrt(X.size()*8+1)+1)/2} generated an error in CRAN
// //' since \code{sqrt} cannot be applied to integers.
// //' @keywords internal
// //' @useDynLib resemble
// //' @author Antoine Stevens 
// // [[Rcpp::plugins(openmp)]]
// // [[Rcpp::export]]  
// NumericVector minDissV(NumericVector X){  
//   // For the distances
//   int nX = X.size();
//   int n = ((nX * nX) - nX) / 2;
//   NumericVector doutput(n);   
//   
//   // For the indices
//   arma::uword  index;
//   double vct = (sqrt(((double)doutput.size()) * 8.0 + 1.0) + 1.0) / 2.0;
//   int len = (int)vct;
//   // int len = (sqrt(Ds.size()*8+1)+1)/2;
//   arma::uvec vindex(len); 
//   int i,j;
// #if defined(_OPENMP) 
// #pragma omp parallel for private(i,j) schedule(dynamic)
// #endif
//   for(int i = 0; i < nX - 1; i++)
//     for(int j = i + 1; j < nX; j++){
//       double x = X(j) - X(i);
//       doutput(nX*i - (i * (i + 3) / 2)  + j  - 1) =  x * x;          
//     }
//     
//     for(i = 0; i < len; i++){       
//       arma::vec x(len); 
//       for(j = 0; j < i; j++){
//         // triangular sequence
//         int k = j * len - (j * (j + 3) / 2) + i - 1;
//         x[j] = doutput(k);        
//       }        
//       for(j = i + 1; j < len; j++){
//         // triangular sequence
//         int k2 = i * len - (i * (i + 3) / 2) + j - 1;
//         x[j] = doutput(k2);             
//       }
//       x[i] = arma::datum::nan; // remove diag
//       x.min(index); // don't assign result to a value since we are interested only in the index
//       vindex[i] = index;    
//     }
//     return wrap(vindex + 1);
// }


// //' @title Moving/rolling correlation distance of two matrices
// //' @description Computes a moving window correlation distance between two data matrices
// //' @usage
// //' moving_cor_diss(X,Y,w)
// //' @param X a matrix
// //' @param Y a matrix
// //' @param w window size (must be odd)
// //' @return a matrix of correlation distance
// //' @keywords internal
// //' @useDynLib resemble
// //' @author Leonardo Ramirez-Lopez and Antoine Stevens
// // [[Rcpp::export]]
// NumericMatrix moving_cor_diss(arma::mat X, arma::mat Y, int w){
//   arma::mat rmwF = arma::zeros(X.n_rows, Y.n_rows);
//   if(!w%2) {
//     throw std::invalid_argument("w must be odd");
//   }
//   int gap = (w - 1) / 2;
//   if((Y.n_cols - w ) < 1){
//     arma::mat rmwF = arma::cor(trans(X), trans(Y));
//   } else {
//     int ny = Y.n_cols;
//     for(int i = gap; i < ny - gap; i++){
//       // sum of the correlations
//       rmwF += arma::cor(X.cols(i-gap, i + gap).t(), Y.cols(i - gap, i + gap).t());
//     }
//     // get the average
//     rmwF = rmwF / (Y.n_cols - (2*gap));
//   }
//   // get the distance
//   arma::mat scmw = (1 - rmwF)/2;
//   return wrap(scmw.t());
// }



















// Minimal helper
static inline void validate_window_or_full(int w, arma::uword T) {
  if (w < 1 || static_cast<arma::uword>(w) > T || (((w % 2) == 0) && w != static_cast<int>(T))) {
    throw std::invalid_argument("w must be odd and ≤ number of columns; alternatively set w == number of columns for a full-window.");
  }
}




// --- Full-window correlation via sums/products (templated) -------------------
template<typename eT>
static inline arma::Mat<eT>
cor_dense_from_stats_xy_t(const arma::Mat<eT> &X, const arma::Mat<eT> &Y) {
  if (X.n_cols != Y.n_cols) {
    throw std::invalid_argument("X and Y must have the same number of columns");
  }
  const eT T  = static_cast<eT>(X.n_cols);
  const eT iw = eT(1) / T;
  const eT ic = eT(1) / (T - eT(1));
  
  arma::Col<eT> Sx  = arma::sum(X, 1);
  arma::Col<eT> Sy  = arma::sum(Y, 1);
  arma::Col<eT> Sxx = arma::sum(arma::square(X), 1);
  arma::Col<eT> Syy = arma::sum(arma::square(Y), 1);
  
  arma::Mat<eT> Sxy = X * Y.t();                        // m×n
  arma::Mat<eT> C   = Sxy - (Sx * Sy.t()) * iw;         // cov*(T-1)
  C *= ic;                                              // cov
  
  arma::Col<eT> varx = (Sxx - (Sx % Sx) * iw) * ic;
  arma::Col<eT> vary = (Syy - (Sy % Sy) * iw) * ic;
  
  arma::Col<eT> inv_sdx = eT(1) / arma::sqrt(arma::clamp(varx, eT(0), std::numeric_limits<eT>::infinity()));
  arma::Col<eT> inv_sdy = eT(1) / arma::sqrt(arma::clamp(vary, eT(0), std::numeric_limits<eT>::infinity()));
  
  C.each_row() %= inv_sdy.t();
  C.each_col() %= inv_sdx;
  return C;                                            // m×n
}

// --- Per-tile rolling corr from stats (templated) ----------------------------
template<typename eT>
static inline void
corr_tile_xy_fast_t(
  arma::Mat<eT> &C,                  // mi×nj (output; reused as cov)
  const arma::Mat<eT> &Sxy,          // mi×nj
  const arma::Col<eT> &Sx_i,         // mi
  const arma::Col<eT> &Sy_j,         // nj
  const arma::Col<eT> &Sxx_i,        // mi
  const arma::Col<eT> &Syy_j,        // nj
  const eT iw,                       // 1/w
  const eT icv,                      // 1/(w-1)
  arma::Col<eT> &var_i,              // mi scratch
  arma::Col<eT> &var_j,              // nj scratch
  arma::Col<eT> &inv_i,              // mi scratch
  arma::Col<eT> &inv_j               // nj scratch
) {
  C = Sxy;
  C -= (Sx_i * Sy_j.t()) * iw;
  C *= icv;
  
  var_i = (Sxx_i - (Sx_i % Sx_i) * iw) * icv;
  var_j = (Syy_j - (Sy_j % Sy_j) * iw) * icv;
  
  inv_i = eT(1) / arma::sqrt(arma::clamp(var_i, eT(0), std::numeric_limits<eT>::infinity()));
  inv_j = eT(1) / arma::sqrt(arma::clamp(var_j, eT(0), std::numeric_limits<eT>::infinity()));
  
  C.each_row() %= inv_j.t();
  C.each_col() %= inv_i;
}

// --- XY rolling distance (templated core, OpenMP tiled; start-index windows) -
template<typename eT>
static inline arma::Mat<eT>
moving_cor_diss_xy_impl(const arma::Mat<eT> &X,
                        const arma::Mat<eT> &Y,
                        int w,
                        arma::uword block_x,
                        arma::uword block_y)
{
  if (X.n_cols != Y.n_cols) {
    throw std::invalid_argument("X and Y must have the same number of columns");
  }
  const arma::uword m = X.n_rows, n = Y.n_rows, T = Y.n_cols;
  
  validate_window_or_full(w, T);
  
  // Full-window path
  if (w == static_cast<int>(T)) {
    arma::Mat<eT> C = cor_dense_from_stats_xy_t<eT>(X, Y);  // m×n
    arma::Mat<eT> D = (eT(1) - C) / eT(2);                  // m×n
    // Clamp distances to [0,1] to remove tiny round-off
    D.transform([](eT d){ return (d < eT(0)) ? eT(0) : (d > eT(1) ? eT(1) : d); });
    return D.t();                                           // n×m
  }
  
  const int  n_centers = static_cast<int>(T) - w + 1;
  const eT   iw = eT(1) / static_cast<eT>(w);
  const eT   icv = eT(1) / static_cast<eT>(w - 1);
  
  arma::Mat<eT> out(n, m, arma::fill::zeros); // n×m
  
  arma::Mat<eT> X_sq = arma::square(X);  // m × T
  arma::Mat<eT> Y_sq = arma::square(Y);  // n × T
  
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(dynamic)
#endif
  for (arma::sword ia = 0; ia < static_cast<arma::sword>(m); ia += static_cast<arma::sword>(block_x)) {
    for (arma::sword ja = 0; ja < static_cast<arma::sword>(n); ja += static_cast<arma::sword>(block_y)) {
      
      const arma::uword i0 = static_cast<arma::uword>(ia);
      const arma::uword j0 = static_cast<arma::uword>(ja);
      const arma::uword i1 = std::min(m, i0 + block_x);
      const arma::uword j1 = std::min(n, j0 + block_y);
      const arma::uword mi = i1 - i0;
      const arma::uword nj = j1 - j0;
      
      auto Xi = X.rows(i0, i1 - 1);  
      auto Yj = Y.rows(j0, j1 - 1);  
      auto Xi_sq = X_sq.rows(i0, i1 - 1);
      auto Yj_sq = Y_sq.rows(j0, j1 - 1);
      
      // First window [0, w-1]
      arma::Col<eT> Sx_i  = arma::sum(Xi.cols(0, w - 1), 1);
      arma::Col<eT> Sy_j  = arma::sum(Yj.cols(0, w - 1), 1);
      arma::Col<eT> Sxx_i = arma::sum(Xi_sq.cols(0, w - 1), 1);
      arma::Col<eT> Syy_j = arma::sum(Yj_sq.cols(0, w - 1), 1);
      arma::Mat<eT> Sxy = Xi.cols(0, w - 1) * Yj.cols(0, w - 1).t(); 
      
      arma::Mat<eT> sumCorr(mi, nj, arma::fill::zeros);
      arma::Mat<eT> C(mi, nj, arma::fill::none);
      arma::Col<eT> var_i(mi, arma::fill::none), var_j(nj, arma::fill::none);
      arma::Col<eT> inv_i(mi, arma::fill::none), inv_j(nj, arma::fill::none);
      
      for (int s = 0; s < n_centers; ++s) {
        corr_tile_xy_fast_t<eT>(
          C, Sxy, Sx_i, Sy_j, Sxx_i, Syy_j, iw, icv, var_i, var_j, inv_i, inv_j
        ); 
        sumCorr += C;
        
        if (s + 1 < n_centers) {
          const int out_idx = s, in_idx = s + w;
          Sx_i  += Xi.col(in_idx)    - Xi.col(out_idx);
          Sy_j  += Yj.col(in_idx)    - Yj.col(out_idx);
          Sxx_i += Xi_sq.col(in_idx) - Xi_sq.col(out_idx);
          Syy_j += Yj_sq.col(in_idx) - Yj_sq.col(out_idx);
          Sxy   += Xi.col(in_idx) * Yj.col(in_idx).t()
            - Xi.col(out_idx) * Yj.col(out_idx).t();
        }
      }
      
      // Optional, cheap: clamp the accumulated correlation sums to [-nc, nc]
      const eT nc = static_cast<eT>(n_centers);
      sumCorr = arma::clamp(sumCorr, -nc, nc);
      
      arma::Mat<eT> Dblk = (eT(1) - (sumCorr / nc)) / eT(2);
      
      // Final safety clamp to [0,1]
      Dblk = arma::clamp(Dblk, eT(0), eT(1));
      
      // disjoint write: Y-rows × X-rows block in output
      out.submat(j0, i0, j1 - 1, i1 - 1) = Dblk.t();
    }
  }
  
  return out;
}


//' @title Rolling correlation distance between X and Y
//' @param X Numeric matrix m×T
//' @param Y Numeric matrix n×T
//' @param w Window size (odd in [1..T] or exactly T)
//' @param block_x Tile size for rows of X (default 1024)
//' @param block_y Tile size for rows of Y (default 1024)
//' @param precision "double" (default) or "float32"/"single"
//' @return n×m distance matrix (R double matrix)
// [[Rcpp::export]]
arma::mat moving_cor_diss_xy(
    const arma::mat &X,
    const arma::mat &Y,
    int w,
    arma::uword block_x = 1024,
    arma::uword block_y = 1024,
    std::string precision = "double"
) {
  for (auto &ch : precision) ch = static_cast<char>(std::tolower(ch));
  
  if (precision == "float32" || precision == "float" || precision == "single") {
    arma::fmat Xf = arma::conv_to<arma::fmat>::from(X);
    arma::fmat Yf = arma::conv_to<arma::fmat>::from(Y);
    arma::fmat Df = moving_cor_diss_xy_impl<float>(Xf, Yf, w, block_x, block_y);
    return arma::conv_to<arma::mat>::from(Df);
  } else {
    return moving_cor_diss_xy_impl<double>(X, Y, w, block_x, block_y);
  }
}























// Compute correlation for one tile in-place using a single large buffer.
// Assumes no zero-variance rows/cols (no masking branches).
static inline void corr_tile_from_stats_serial_fast(
    arma::mat &C,                 // mi x mj (output; reused as "cov" scratch)
    const arma::mat &Sxy,         // mi x mj
    const arma::vec &Sx_i,        // mi
    const arma::vec &Sx_j,        // mj
    const arma::vec &Sxx_i,       // mi
    const arma::vec &Sxx_j,       // mj
    const double iw,              // 1.0 / w
    const double icv,             // 1.0 / (w-1)
    arma::vec &var_i,             // mi scratch
    arma::vec &var_j,             // mj scratch
    arma::vec &inv_i,             // mi scratch
    arma::vec &inv_j              // mj scratch
) {
  // cov -> C (avoid extra big copy)
  C = Sxy;
  C -= (Sx_i * Sx_j.t()) * iw;
  C *= icv;
  
  // sd -> reciprocals (avoid divisions later)
  var_i = (Sxx_i - (Sx_i % Sx_i) * iw) * icv;
  var_j = (Sxx_j - (Sx_j % Sx_j) * iw) * icv;
  
  // clamp tiny negatives from roundoff; still vectorized
  inv_i = 1.0 / arma::sqrt(arma::clamp(var_i, 0.0, std::numeric_limits<double>::infinity()));
  inv_j = 1.0 / arma::sqrt(arma::clamp(var_j, 0.0, std::numeric_limits<double>::infinity()));
  
  // r = cov ./ (sd_i sd_j^T)  -> scale by reciprocals
  C.each_row() %= inv_j.t();
  C.each_col() %= inv_i;
}

//' @title Rolling correlation distance within X (upper-triangle, tiled, serial-optimized)
//' @description Mean rolling-window correlation distance among rows of X.
//'              Serial only; computes upper-triangle tiles densely and mirrors.
//' @param X Numeric matrix (m x T)
//' @param w Odd window size
//' @param block_rows Tile size in rows (default 1024)
//' @return m x m symmetric distance matrix
//' @useDynLib resemble, .registration=TRUE
// [[Rcpp::export]]
arma::mat moving_cor_diss_self_f64(const arma::mat &X, int w,
                                   arma::uword block_rows = 1024) {
  const arma::uword T = X.n_cols;
  validate_window_or_full(w, T);
  const arma::uword m_u = X.n_rows;
  
  // Full-window fallback if T < w
  if (T < static_cast<arma::uword>(w)) {
    arma::mat D = (1.0 - arma::cor(X.t())) / 2.0;
    // Clamp to [0,1] to kill tiny round-off excursions
    D.transform([](double d){ return (d < 0.0) ? 0.0 : (d > 1.0 ? 1.0 : d); });
    // Symmetrise & zero diagonal
    D = 0.5 * (D + D.t());
    D.diag().zeros();
    return D;
  }
  
  const int   m  = static_cast<int>(m_u);
  const int   bs = static_cast<int>(block_rows);
  const int   n_centers = static_cast<int>(T) - w + 1;
  const double iw  = 1.0 / static_cast<double>(w);
  const double icv = 1.0 / static_cast<double>(w - 1);
  
  arma::mat out(m_u, m_u, arma::fill::zeros);
  
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for (int i0 = 0; i0 < m; i0 += bs) {
    for (int j0 = i0; j0 < m; j0 += bs) {
      const int i1 = std::min(m, i0 + bs);
      const int j1 = std::min(m, j0 + bs);
      const int mi = i1 - i0;
      const int mj = j1 - j0;
      
      auto Xi = X.rows(static_cast<arma::uword>(i0),
                       static_cast<arma::uword>(i1 - 1));
      auto Xj = X.rows(static_cast<arma::uword>(j0),
                       static_cast<arma::uword>(j1 - 1));
      
      arma::mat Xi_sq = arma::square(Xi);
      arma::mat Xj_sq = arma::square(Xj);
      arma::vec Sx_i  = arma::sum(Xi.cols(0, w - 1), 1);
      arma::vec Sx_j  = arma::sum(Xj.cols(0, w - 1), 1);
      arma::vec Sxx_i = arma::sum(Xi_sq.cols(0, w - 1), 1);
      arma::vec Sxx_j = arma::sum(Xj_sq.cols(0, w - 1), 1);
      arma::mat Sxy   = Xi.cols(0, w - 1) * Xj.cols(0, w - 1).t();
      
      arma::mat sumCorr(mi, mj, arma::fill::zeros);
      arma::mat C(mi, mj, arma::fill::none);
      arma::vec var_i(mi, arma::fill::none), var_j(mj, arma::fill::none);
      arma::vec inv_i(mi, arma::fill::none), inv_j(mj, arma::fill::none);
      
      for (int s = 0; s < n_centers; ++s) {
        corr_tile_from_stats_serial_fast(
          C, Sxy, Sx_i, Sx_j, Sxx_i, Sxx_j, iw, icv, var_i, var_j, inv_i, inv_j
        );
        sumCorr += C;
        
        if (s + 1 < n_centers) {
          const int out_idx = s, in_idx = s + w;
          Sx_i  += Xi.col(in_idx)    - Xi.col(out_idx);
          Sx_j  += Xj.col(in_idx)    - Xj.col(out_idx);
          Sxx_i += Xi_sq.col(in_idx) - Xi_sq.col(out_idx);
          Sxx_j += Xj_sq.col(in_idx) - Xj_sq.col(out_idx);
          Sxy   += Xi.col(in_idx) * Xj.col(in_idx).t()
            - Xi.col(out_idx) * Xj.col(out_idx).t();
        }
      }
      
      // Optional: clamp average correlation range in one cheap pass
      const double nc = static_cast<double>(n_centers);
      sumCorr.transform([nc](double v){
        return (v < -nc) ? -nc : (v > nc ? nc : v);
      });
      
      arma::mat Dblk = (1.0 - (sumCorr / nc)) / 2.0;
      
      // Final safety clamp to [0,1]
      Dblk.transform([](double d){
        return (d < 0.0) ? 0.0 : (d > 1.0 ? 1.0 : d);
      });
      
      out.submat(static_cast<arma::uword>(i0), static_cast<arma::uword>(j0),
                 static_cast<arma::uword>(i1 - 1),
                 static_cast<arma::uword>(j1 - 1)) = Dblk;
    }
  }
  
  // Single mirror to lower triangle & zero diagonal
  for (int j = 0; j < m; ++j) {
    out(j, j) = 0.0;
    for (int i = j + 1; i < m; ++i) out(i, j) = out(j, i);
  }
  return out;
}




// ─────────────────────────────────────────────────────────────────────────────
// Float32 tile
// ─────────────────────────────────────────────────────────────────────────────
static inline void corr_tile_from_stats_serial_fast_f32(
    fmat       &C,        // mi × mj (output; reused as "cov" scratch)
    const fmat &Sxy,      // mi × mj
    const fvec &Sx_i,     // mi
    const fvec &Sx_j,     // mj
    const fvec &Sxx_i,    // mi
    const fvec &Sxx_j,    // mj
    const float iw,       // 1.0f / w
    const float icv,      // 1.0f / (w-1)
    fvec &var_i,          // mi scratch
    fvec &var_j,          // mj scratch
    fvec &inv_i,          // mi scratch
    fvec &inv_j           // mj scratch
) {
  // cov -> C
  C = Sxy;
  C -= (Sx_i * Sx_j.t()) * iw;    // rank-1 correction
  C *= icv;
  
  // sd (variance → inv sd)
  var_i = (Sxx_i - (Sx_i % Sx_i) * iw) * icv;
  var_j = (Sxx_j - (Sx_j % Sx_j) * iw) * icv;
  
  // clamp tiny negatives from roundoff
  inv_i = 1.0f / sqrt(clamp(var_i, 0.0f, std::numeric_limits<float>::infinity()));
  inv_j = 1.0f / sqrt(clamp(var_j, 0.0f, std::numeric_limits<float>::infinity()));
  
  // r = cov ./ (sd_i sd_j^T)
  C.each_row() %= inv_j.t();
  C.each_col() %= inv_i;
}

// ─────────────────────────────────────────────────────────────────────────────
//' @title Rolling correlation distance within X (float32, tiled)
//' @description Mean rolling-window correlation distance among rows of X (single precision).
//' @param X Numeric matrix (m x T)
//' @param w Odd window size
//' @param block_rows Tile size in rows (default 1024)
//' @return m x m symmetric distance matrix (returned as double for R)
//' @useDynLib resemble, .registration=TRUE
// [[Rcpp::export]]
arma::mat moving_cor_diss_self_f32(const arma::mat &X, int w,
                                   arma::uword block_rows = 1024) {
  const uword T   = X.n_cols;
  validate_window_or_full(w, T);
  const uword m_u = X.n_rows;
  
  // Convert once to float32 (R is double; conversion cost is negligible vs compute)
  fmat Xf = conv_to<fmat>::from(X);
  
  // Full-window fallback if T < w  (use float path)
  if (T < static_cast<uword>(w)) {
    // r = cor(t(X)) in float, then D = (1 - r)/2; finally return double
    fmat Xt = Xf.t();
    fmat Rf = cor(Xt);                   // m×m, float
    fmat Df = 0.5f * (1.0f - Rf);
    Df.transform([](float d){
      return (d < 0.0f) ? 0.0f : (d > 1.0f ? 1.0f : d);
    });
    Df.diag().zeros();
    Df = 0.5f * (Df + Df.t());           // symmetry hygiene
    return conv_to<mat>::from(Df);
  }
  
  const int   m  = static_cast<int>(m_u);
  const int   bs = static_cast<int>(block_rows);
  const int   n_centers = static_cast<int>(T) - w + 1;
  const float iw  = 1.0f / static_cast<float>(w);
  const float icv = 1.0f / static_cast<float>(w - 1);
  
  // Float working buffers; final ‘out’ we keep in double (you can also keep float and cast once at end)
  fmat out_f(m_u, m_u, fill::zeros);
  
  // Parallel upper-triangle tiles
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for (int i0 = 0; i0 < m; i0 += bs) {
    for (int j0 = i0; j0 < m; j0 += bs) {
      const int i1 = std::min(m, i0 + bs);
      const int j1 = std::min(m, j0 + bs);
      const int mi = i1 - i0;
      const int mj = j1 - j0;
      
      auto Xi = Xf.rows(static_cast<uword>(i0), static_cast<uword>(i1 - 1));
      auto Xj = Xf.rows(static_cast<uword>(j0), static_cast<uword>(j1 - 1));
      
      fmat Xi_sq = square(Xi);
      fmat Xj_sq = square(Xj);
      
      // First window [0, w-1]
      fvec Sx_i  = sum(Xi.cols(0, w - 1), 1);
      fvec Sx_j  = sum(Xj.cols(0, w - 1), 1);
      fvec Sxx_i = sum(Xi_sq.cols(0, w - 1), 1);
      fvec Sxx_j = sum(Xj_sq.cols(0, w - 1), 1);
      fmat Sxy   = Xi.cols(0, w - 1) * Xj.cols(0, w - 1).t();   // SGEMM
      
      fmat sumCorr(mi, mj, fill::zeros);
      fmat C(mi, mj, fill::none);
      fvec var_i(mi, fill::none), var_j(mj, fill::none);
      fvec inv_i(mi, fill::none), inv_j(mj, fill::none);
      
      for (int s = 0; s < n_centers; ++s) {
        corr_tile_from_stats_serial_fast_f32(
          C, Sxy, Sx_i, Sx_j, Sxx_i, Sxx_j, iw, icv, var_i, var_j, inv_i, inv_j
        );
        sumCorr += C;
        
        if (s + 1 < n_centers) {
          const int out_idx = s, in_idx = s + w;
          Sx_i  += Xi.col(in_idx)    - Xi.col(out_idx);
          Sx_j  += Xj.col(in_idx)    - Xj.col(out_idx);
          Sxx_i += Xi_sq.col(in_idx) - Xi_sq.col(out_idx);
          Sxx_j += Xj_sq.col(in_idx) - Xj_sq.col(out_idx);
          Sxy   += Xi.col(in_idx) * Xj.col(in_idx).t()
            - Xi.col(out_idx) * Xj.col(out_idx).t();
        }
      }
      
      fmat Dblk = 0.5f * (1.0f - (sumCorr / static_cast<float>(n_centers)));
      
      Dblk.transform([](float d){
        return (d < 0.0f) ? 0.0f : (d > 1.0f ? 1.0f : d);
      });
      
      out_f.submat(static_cast<uword>(i0), static_cast<uword>(j0),
                   static_cast<uword>(i1 - 1), static_cast<uword>(j1 - 1)) = Dblk;
    }
  }
  
  // Mirror once & zero diag (still in float)
  for (int j = 0; j < m; ++j) {
    out_f(j, j) = 0.0f;
    for (int i = j + 1; i < m; ++i) out_f(i, j) = out_f(j, i);
  }
  
  // Return as double (R owns doubles); or change signature to return float-packed file
  return conv_to<mat>::from(out_f);
}

// --- R wrapper: precision selectable ("double" | "float32") ------------------

//' @title Rolling correlation distance within X (templated, OpenMP; precision)
//' @param X Numeric matrix (m×T)
//' @param w Window size (odd in [1..T] or exactly T)
//' @param block_rows Tile size (default 1024)
//' @param precision "double" (default) or "single (i.e."float32")
//' @return m×m distance matrix (double for R)
// [[Rcpp::export]]
arma::mat moving_cor_diss_self(const arma::mat &X,
                               int w,
                               arma::uword block_rows = 1024,
                               std::string precision = "double") {
  if (precision == "double") {
    return moving_cor_diss_self_f64(X, w, block_rows);
  } else if (precision == "single") {
    return moving_cor_diss_self_f32(X, w, block_rows);
  } else {
    // Default to double if unspecified
    return moving_cor_diss_self_f64(X, w, block_rows);
  }
}






























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
  // #if defined(_OPENMP) 
  // #pragma omp parallel for schedule(static) 
  // #endif
  for(int i = 0; i < nX; i++){
    arma::rowvec x = XX.row(i);
    x(i) = arma::datum::nan; // remove diag
    x.min(index); // don't assign result to a value since we are interested only in the index
    vindex[i] = index;    
  }
  return wrap(vindex +1);   
}


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
  int i;
  int j;
  // #if defined(_OPENMP) 
  // #pragma omp parallel for private(i,j) schedule(dynamic)
  // #endif

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
  return wrap(vindex + 1);
}


//' @title Top-k Order Indices for Matrix Columns
//' @description 
//' Internal helper function. For each column in a numeric matrix, 
//' returns the indices of the top-\code{k} smallest values (i.e., closest neighbors).
//'
//' @details 
//' This function is implemented in C++ using \code{Rcpp}. It performs a partial sort 
//' on each column to retrieve the \code{k} smallest elements efficiently. This is 
//' particularly useful when working with large dissimilarity matrices, where 
//' performance is critical.
//'
//' @param mat A numeric matrix (e.g., dissimilarity matrix).
//' @param k Integer. Number of top smallest elements to return for each column.
//'
//' @return An integer matrix of dimensions \code{k × ncol(mat)}. Each column contains 
//' the 1-based row indices of the top-\code{k} smallest values in the corresponding column of \code{mat}.
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
IntegerMatrix top_k_order(
    const arma::mat& mat, 
    int k, 
    IntegerVector skip = IntegerVector::create()
) {
  int nrow = mat.n_rows;
  int ncol = mat.n_cols;
  IntegerMatrix result(k, ncol);
  
  std::unordered_set<int> skip_set(skip.begin(), skip.end());
  
  for (int j = 0; j < ncol; ++j) {
    std::vector<std::pair<double, int>> values;
    
    // Collect finite values with their row indices
    for (int i = 0; i < nrow; ++i) {
      double val = mat(i, j);
      if (std::isnan(val)) continue; // mimic na.last = TRUE
      values.push_back({val, i + 1}); // R uses 1-based indexing
    }
    
    // Sort all values ascending
    std::sort(values.begin(), values.end());
    
    std::vector<int> selected;
    for (size_t i = 0; i < values.size(); ++i) {
      if (i == 0) {
        // Always keep the first (closest)
        selected.push_back(values[i].second);
      } else {
        if (skip_set.find(values[i].second) == skip_set.end()) {
          selected.push_back(values[i].second);
        }
      }
      if (selected.size() == (size_t)k) break;
    }
    
    // Fill column
    for (int i = 0; i < k; ++i) {
      result(i, j) = (i < (int)selected.size()) ? selected[i] : NA_INTEGER;
    }
  }
  
  return result;
}


//' @title Extract Column Elements by Row Indices
//' @description
//' Efficiently extracts specific elements from each column of a numeric matrix 
//' based on provided row indices, using Rcpp for performance.
//'
//' @param mat A numeric matrix from which values will be extracted.
//' @param idx An integer matrix of the same number of columns as \code{mat}. 
//' Each column of \code{idx} contains the row indices to extract from the 
//' corresponding column of \code{mat}.
//'
//' @details
//' For each column \code{j}, this function returns a vector of values 
//' \code{mat[idx[, j], j]}. It is implemented in C++ using Rcpp for speed 
//' and is suitable for large-scale operations (e.g., tens of thousands of rows).
//'
//' Indices in \code{idx} are assumed to be 1-based (R-style). Internally, 
//' they are adjusted for 0-based indexing in C++.
//'
//' @return A numeric matrix of the same dimension as \code{idx}, where each 
//' element is extracted from \code{mat} according to the corresponding 
//' row and column in \code{idx}.
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
NumericMatrix extract_by_index(NumericMatrix mat, IntegerMatrix idx) {
  int nrow = idx.nrow();
  int ncol = idx.ncol();
  NumericMatrix out(nrow, ncol);
  
#pragma omp parallel for
  for (int j = 0; j < ncol; j++) {
    for (int i = 0; i < nrow; i++) {
      out(i, j) = mat(idx(i, j) - 1, j);
    }
  }
  
  return out;
}

//' @title Group Disparity Mask for Nearest Neighbors (C++)
//' @description
//' Internal C++ function that creates a logical matrix indicating whether
//' each neighbor in `kidxmat` belongs to a different group than the sample.
//'
//' @param kidxmat An integer matrix of neighbor indices. Each column
//'        corresponds to a sample, and each row to a ranked neighbor.
//'        Indices refer to rows in the `group` vector.
//' @param group An integer vector representing the group assignment
//'        for each sample. Its length must match the number of rows
//'        in the reference data.
//'
//' @return
//' A logical matrix of the same dimensions as `kidxmat`, where each
//' element is `TRUE` if the corresponding neighbor belongs to a
//' different group than the target sample, and `FALSE` otherwise.
//'
//' @details
//' For each column in `kidxmat`, corresponding to a sample,
//' the function:
//' \itemize{
//'   \item Extracts the group assignment of the sample.
//'   \item Checks whether each neighbor (by index) belongs to the same group.
//'   \item Stores `TRUE` if the neighbor is from a different group,
//'         `FALSE` otherwise.
//' }
//' The result can be used as a mask to filter out same-group neighbors.
//'
//' @note
//' Parallelized using OpenMP for performance.
//' Assumes 1-based indexing from R (internally adjusted).
//' @noRd
//' @keywords internal
//' @author Leonardo Ramirez-Lopez

// [[Rcpp::export]]
Rcpp::LogicalMatrix not_in_same_group(
    const Rcpp::IntegerMatrix& kidxmat,
    const Rcpp::IntegerVector& group
) {
  int nrow = kidxmat.nrow();
  int ncol = kidxmat.ncol();
  Rcpp::LogicalMatrix result(nrow, ncol);
  
#pragma omp parallel for
  for (int j = 0; j < ncol; ++j) {
    int grp_j = group[j];
    for (int i = 0; i < nrow; ++i) {
      int idx = kidxmat(i, j) - 1;
      result(i, j) = (group[idx] != grp_j);
    }
  }
  
  return result;
}


//' @title Quantile Statistics from Nearest Neighbors (C++)
//' @description
//' Internal C++ function to compute quantiles of response values
//' among filtered nearest neighbors for each sample.
//'
//' @param kidxmat An integer matrix of nearest neighbor indices.
//'        Each column corresponds to one observation (e.g., in Xu),
//'        and rows represent the ranks of the neighbors.
//' @param kidxgrop A logical matrix of the same shape as `kidxmat`
//'        indicating which neighbors should be retained (`TRUE`)
//'        based on group membership.
//' @param Yr A numeric vector of reference response values. The
//'        values are indexed by `kidxmat` to compute quantiles.
//' @param k An integer vector of the number of neighbors to use
//'        for each observation (length must match `ncol(kidxmat)`).
//' @param probs A numeric vector of probabilities for which to
//'        compute quantiles (e.g., `c(0, 0.05, 0.5, 0.95, 1)`).
//'
//' @return
//' A numeric matrix of size `(length(k) * ncol(kidxmat)) × length(probs)`
//' where each row contains the quantiles for a particular observation
//' and each column corresponds to a probability in `probs`.
//'
//' @details
//' For each observation, the function:
//' \itemize{
//'   \item Selects its top-k nearest neighbors.
//'   \item Filters neighbors based on the logical group mask.
//'   \item Extracts corresponding values from `Yr`.
//'   \item Computes requested quantiles.
//' }
//' Missing values in `Yr` are ignored. If no valid neighbors are
//' found for a given sample, all quantiles will be set to `NA`.
//'
//' @note
//' This is an internal function intended for performance-critical
//' applications. It should not be exported to the user.
//' @noRd
//' @keywords internal
//' @author Leonardo Ramirez-Lopez
// [[Rcpp::export]]
NumericMatrix compute_nn_quantiles(
    const IntegerMatrix& kidxmat,
    const LogicalMatrix& kidxgrop,
    const NumericVector& Yr,
    const IntegerVector& k,
    const NumericVector& probs
) {
  int n_samples = k.size();
  int n_neighbors = kidxmat.ncol();
  int n_quantiles = probs.size();
  
  NumericMatrix result(n_samples * n_neighbors, n_quantiles);
  
#pragma omp parallel for
  for (int i = 0; i < n_samples; i++) {
    int ik = k[i];
    for (int j = 0; j < n_neighbors; j++) {
      std::vector<double> values;
      for (int n = 0; n < ik; n++) {
        if (kidxgrop(n, j)) {
          int idx = kidxmat(n, j) - 1;
          if (!NumericVector::is_na(Yr[idx])) {
            values.push_back(Yr[idx]);
          }
        }
      }
      
      NumericVector q(n_quantiles, NA_REAL);
      if (!values.empty()) {
        std::sort(values.begin(), values.end());
        for (int qidx = 0; qidx < n_quantiles; qidx++) {
          double pos = probs[qidx] * (values.size() - 1);
          int idx_below = std::floor(pos);
          int idx_above = std::ceil(pos);
          double weight = pos - idx_below;
          if (idx_below == idx_above) {
            q[qidx] = values[idx_below];
          } else {
            q[qidx] = (1 - weight) * values[idx_below] + weight * values[idx_above];
          }
        }
      }
      
      for (int qidx = 0; qidx < n_quantiles; qidx++) {
        result(i * n_neighbors + j, qidx) = q[qidx];
      }
    }
  }
  
  return result;
}