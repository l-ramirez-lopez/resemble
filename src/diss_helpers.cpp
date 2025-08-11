#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <limits>
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
NumericVector fast_diss_vector(NumericVector X){  
  int nX = X.size();
  int n = ((nX*nX)-nX)/2;
  NumericVector output(n);  
  // #if defined(_OPENMP) 
  // #pragma omp parallel for schedule(dynamic)
  // #endif
  for(int i = 0; i < nX-1; i++)
    for(int j = i+1; j < nX; j++){
      double x = X(j)-X(i);
      output(nX*i - (i * (i + 3) / 2)  + j  - 1) =  x * x; 
    }
  return output;         
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
                        arma::uword block_y) {
  if (X.n_cols != Y.n_cols) {
    throw std::invalid_argument("X and Y must have the same number of columns");
  }
  const arma::uword m = X.n_rows, n = Y.n_rows, T = Y.n_cols;
  
  validate_window_or_full(w, T);
  
  // full-window path
  if (w == static_cast<int>(T)) {
    arma::Mat<eT> C = cor_dense_from_stats_xy_t<eT>(X, Y); // m×n
    arma::Mat<eT> D = (eT(1) - C) / eT(2);                 // m×n
    return D.t();                                          // n×m
  }
  
  const int  n_centers = static_cast<int>(T) - w + 1;
  const eT   iw = eT(1) / static_cast<eT>(w);
  const eT   icv = eT(1) / static_cast<eT>(w - 1);
  
  arma::Mat<eT> out(n, m, arma::fill::zeros); // n×m
  
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
      arma::Mat<eT> Xi_sq = arma::square(Xi);                 // precompute per tile
      arma::Mat<eT> Yj_sq = arma::square(Yj);
      
      // first window [0, w-1]
      arma::Col<eT> Sx_i  = arma::sum(Xi.cols(0, w - 1), 1);
      arma::Col<eT> Sy_j  = arma::sum(Yj.cols(0, w - 1), 1);
      arma::Col<eT> Sxx_i = arma::sum(Xi_sq.cols(0, w - 1), 1);
      arma::Col<eT> Syy_j = arma::sum(Yj_sq.cols(0, w - 1), 1);
      arma::Mat<eT> Sxy   = Xi.cols(0, w - 1) * Yj.cols(0, w - 1).t(); // mi×nj
      
      arma::Mat<eT> sumCorr(mi, nj, arma::fill::zeros);
      arma::Mat<eT> C(mi, nj, arma::fill::none);
      arma::Col<eT> var_i(mi, arma::fill::none), var_j(nj, arma::fill::none);
      arma::Col<eT> inv_i(mi, arma::fill::none), inv_j(nj, arma::fill::none);
      
      for (int s = 0; s < n_centers; ++s) {
        corr_tile_xy_fast_t<eT>(C, Sxy, Sx_i, Sy_j, Sxx_i, Syy_j, iw, icv,
                                var_i, var_j, inv_i, inv_j);
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
      
      arma::Mat<eT> Dblk = (eT(1) - (sumCorr / static_cast<eT>(n_centers))) / eT(2); // mi×nj
      // disjoint write: Y-rows × X-rows block in output
      out.submat(j0, i0, j1 - 1, i1 - 1) = Dblk.t();
    }
  }
  
  return out;
}

// --- R wrapper: precision = "double" (default) or "float32" ------------------

//' @title Rolling correlation distance between X and Y (templated, OpenMP)
//' @param X Numeric matrix m×T
//' @param Y Numeric matrix n×T
//' @param w Window size (odd in [1..T] or exactly T)
//' @param block_x Tile size for rows of X (default 1024)
//' @param block_y Tile size for rows of Y (default 1024)
//' @param precision "double" (default) or "float32"
//' @return n×m distance matrix (R double matrix)
// [[Rcpp::export]]
arma::mat moving_cor_diss_xy_prec(const arma::mat &X,
                                 const arma::mat &Y,
                                 int w,
                                 arma::uword block_x = 1024,
                                 arma::uword block_y = 1024,
                                 std::string precision = "double") {
 for (auto &ch : precision) ch = static_cast<char>(std::tolower(ch));
 
 if (precision == "float32" || precision == "float" || precision == "single") {
   arma::fmat Xf = arma::conv_to<arma::fmat>::from(X);
   arma::fmat Yf = arma::conv_to<arma::fmat>::from(Y);
   arma::fmat Df = moving_cor_diss_xy_impl<float>(Xf, Yf, w, block_x, block_y);
   return arma::conv_to<arma::mat>::from(Df);  // <-- fix: proper fmat -> mat conversion
 } else {
   arma::mat  Dd = moving_cor_diss_xy_impl<double>(X, Y, w, block_x, block_y);
   return Dd;
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
   
   // Parallel upper-triangle tiles: j0 starts at i0
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
       
       arma::mat Dblk = (1.0 - (sumCorr / static_cast<double>(n_centers))) / 2.0;
       out.submat(static_cast<arma::uword>(i0), static_cast<arma::uword>(j0),
                  static_cast<arma::uword>(i1 - 1),
                  static_cast<arma::uword>(j1 - 1)) = Dblk;
     }
   }
   
   /* mirror once, single-thread */
   for (int j = 0; j < m; ++j) {
     out(j, j) = 0.0;
     for (int i = j + 1; i < m; ++i) out(i, j) = out(j, i);
   }
   
   // Single-thread mirror to lower triangle & zero diag
   for (int j = 0; j < m; ++j) {
     out(j, j) = 0.0;
     for (int i = j + 1; i < m; ++i) {
       out(i, j) = out(j, i);
     }
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
// Float32 main: rolling correlation distance (upper triangle, tiled, OpenMP)
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
   arma::mat Df = moving_cor_diss_self_f64(X, w, block_rows);
   return Df;
 }
 if (precision == "single") {
   arma::mat Df = moving_cor_diss_self_f32(X, w, block_rows);
   return Df;
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

