#include <RcppArmadillo.h> 
#include <math.h>
#include <iostream>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;


//' @title Function for identifiying the column in a matrix with the largest standard deviation
//' @description Identifies the column with the largest standard deviation. For internal use only!
//' @usage get_col_largest_sd(X)
//' @param X a matrix.
//' @return a value indicating the index of the column with the largest standard deviation. 
//' @author Leonardo Ramirez-Lopez
//' @keywords internal 
//' @useDynLib resemble
//' @noRd
// [[Rcpp::export]]
NumericVector get_col_largest_sd(arma::mat X){
  //arma::mat Xz(X.begin(), X.nrow(), X.ncol(), false); 
  int nx = X.n_cols;
  //int nx = X.ncol();
  int maxsd = 0;
  
  double jstd;
  double jstdn;
  
  // Identify the column with the largest standard deviation
  maxsd = 0;
  for (int j = 1; j < nx; j++) {
    jstd = arma::stddev(X.col(maxsd));
    jstdn = arma::stddev(X.col(j));
    if(jstdn > jstd){
      maxsd = j;
    }
  }
  return Rcpp::wrap(maxsd);
}


//' @title Function for computing the standard deviation of each column in a matrix
//' @description Computes the standard deviation of each column in a matrix. For internal use only!
//' @usage get_column_sds(X)
//' @param X a a matrix.
//' @return a vector of standard deviation values. 
//' @author Leonardo Ramirez-Lopez
//' @keywords internal 
//' @noRd
//' @useDynLib resemble
// [[Rcpp::export]]
NumericVector get_column_sds(arma::mat X){
 arma::mat sds = arma::stddev(X, 0, 0);
 return Rcpp::wrap(sds);
}

//' @title Function for computing the overall variance of a matrix
//' @description Computes the variance of a matrix. For internal use only!
//' @usage overall_var(X)
//' @param X a matrix.
//' @return a vector of standard deviation values. 
//' @author Leonardo Ramirez-Lopez
//' @keywords internal 
//' @noRd
//' @useDynLib resemble
// [[Rcpp::export]]
NumericVector overall_var(arma::mat X){
  double n_rows_x = (double)X.n_rows;
  double ovar = sum(var(X)) * (n_rows_x - 1);
  return Rcpp::wrap(ovar);
}

//' @title Function for computing the mean of each column in a matrix
//' @description Computes the mean of each column in a matrix. For internal use only!
//' @usage get_column_means(X)
//' @param X a a matrix.
//' @return a vector of mean values. 
//' @author Leonardo Ramirez-Lopez
//' @keywords internal 
//' @useDynLib resemble
//' @noRd
// [[Rcpp::export]]
NumericVector get_column_means(arma::mat X){
  arma::mat mn = arma::mean(X, 0);
  return Rcpp::wrap(mn);
}

//' @title Function for computing the maximum of each column
//' @description Computes the max of each column in a matrix. For internal use only!
//' @usage get_column_maxs(X)
//' @param X a a matrix.
//' @return a vector of mean values. 
//' @author Leonardo Ramirez-Lopez
//' @keywords internal 
//' @noRd
//' @useDynLib resemble
// [[Rcpp::export]]
NumericVector get_column_maxs(arma::mat X) {
  arma::rowvec mx = arma::max(X, 0);
  return Rcpp::wrap(mx);
}

//' @title Function for computing sum of each column in a matrix
//' @description Computes the sum of each column in a matrix. For internal use only!
//' @usage get_column_sums(X)
//' @param X a matrix.
//' @return a vector of standard deviation values. 
//' @author Leonardo Ramirez-Lopez
//' @keywords internal 
//' @useDynLib resemble
//' @noRd
// [[Rcpp::export]]
NumericVector get_column_sums(arma::mat X){
  arma::mat sm = arma::sum(X);
  return Rcpp::wrap(sm);
}

//' @title Computes the weights for pls regressions
//' @description
//' This is an internal function that computes the wights required for obtaining
//' each vector of pls scores. Implementation is done in C++ for improved performance.
//' @param X a numeric matrix of spectral data.
//' @param Y a matrix of one column with the response variable.
//' @param algorithm a character string indicating what method to use. Options are:
//' \code{'pls'} for pls (using covariance between X and Y), 
//' \code{'mpls'} for modified pls (using correlation between X and Y as in 
//' Shenk and Westerhaus, 1991; Westerhaus 2014) or
//' \code{'xls'} for extended pls (as implemented in BUCHI NIRWise PLUS software).
//' @param xls_min_w an integer indicating the minimum window size for the "xls"
//' method. Only used if \code{algorithm = 'xls'}. Default is 3 (as in BUCHI NIRWise PLUS software).
//' @param xls_max_w an integer indicating the maximum window size for the "xls"
//' method. Only used if \code{algorithm = 'xls'}. Default is 15 (as in BUCHI NIRWise PLUS software).
//' @author Leonardo Ramirez-Lopez and Claudio Orellano
//' @references
//' Shenk, J. S., & Westerhaus, M. O. (1991). Populations structuring of 
//' near infrared spectra and modified partial least squares regression. 
//' Crop Science, 31(6), 1548-1555.
//' 
//' Westerhaus, M. (2014). Eastern Analytical Symposium Award for outstanding 
//' Wachievements in near infrared spectroscopy: my contributions to 
//' Wnear infrared spectroscopy. NIR news, 25(8), 16-20.
//' @return a `matrix` of one column containing the weights.
//' @keywords internal 
//' @noRd
//' @useDynLib resemble
// [[Rcpp::export]]
arma::mat get_weights(const arma::mat& X, 
                      const arma::mat& Y, 
                      const String& algorithm = "pls", 
                      const int xls_min_w = 3, 
                      const int xls_max_w = 15) {
  arma::mat w;
  
  if (algorithm == "pls") {
    w = trans(X) * Y;
  } 
  else if (algorithm == "mpls") {
    w = arma::cor(X, Y);
  }
  else if (algorithm == "xls") {
    int n_cols_x = X.n_cols;
    w = arma::zeros<arma::mat>(n_cols_x, 1);
    for (int i = 0; i < n_cols_x; i++) {
      int j_max = std::min(i + xls_max_w, n_cols_x - 1);
      for (int j = i + xls_min_w; j <= j_max; j++) {
        double corr_val = arma::as_scalar(arma::cor(Y, X.col(i) - X.col(j)));
        w(i, 0) += corr_val;
        w(j, 0) -= corr_val;
      }
    }
  }
  
  // Normalise to unit length
  double norm_w = arma::norm(w);
  w /= norm_w;
  
  return w;
}

//' @title Internal Cpp function for computing the weights of the PLS components 
//' necessary for weighted average PLS
//' @description For internal use only!. 
//' @usage
//' get_local_pls_weights(projection_mat, 
//'           xloadings, 
//'           coefficients, 
//'           new_x, 
//'           min_component, 
//'           max_component, 
//'           scale, 
//'           Xcenter, 
//'           Xscale)
//' @param projection_mat the projection matrix generated either by the \code{opls} function.
//' @param xloadings .
//' @param coefficients the matrix of regression coefficients.
//' @param new_x a matrix of one new spectra to be predicted.
//' @param min_component an integer indicating the minimum number of pls components.
//' @param max_component an integer indicating the maximum number of pls components.
//' @param scale a logical indicating whether the matrix of predictors used to create the regression model was scaled.
//' @param Xcenter a vector with the values that must be used for centering \code{newdata}.
//' @param Xscale if \code{scale = TRUE} a vector with the values that must be used for scaling \code{newdata}.
//' @return a matrix of one row with the weights for each component between the max. and min. specified. 
//' @author Leonardo Ramirez-Lopez
//' @noRd
//' @keywords internal 
//' @useDynLib resemble
// [[Rcpp::export]]
Rcpp::NumericMatrix get_local_pls_weights(
    arma::mat projection_mat, 
    arma::mat xloadings,
    arma::mat coefficients,
    arma::mat new_x,
    int min_component, 
    int max_component, 
    bool scale,
    arma::rowvec Xcenter,
    arma::rowvec Xscale
) {
   arma::mat Xz = new_x;
   arma::mat whgt;
   
   if(scale){
     Xz = Xz / Xscale;
   }
   
   //Necessary to center
   Xz = Xz - Xcenter;
   
   arma::mat xrmsres = arma::zeros(1, max_component);
   
   arma::mat sc = Xz * projection_mat.cols(0, max_component - 1);
   for(int i = (min_component - 1); i < max_component; i++){
     arma::mat xrec = sc.cols(0,i) * xloadings.rows(0, i);
     xrmsres.col(i) = sqrt(arma::mean(arma::mean(pow(Xz - xrec, 2), 0), 1));
   }
   
   arma::mat rmsb = sqrt(get_column_means(pow(coefficients.cols(0, max_component - 1), 2)));
   arma::mat rmsb_x = trans(rmsb.rows(min_component - 1, max_component - 1)) % xrmsres.cols(min_component - 1, max_component - 1);
   arma::mat whgtn = pow(rmsb_x, -1);
   whgt  = whgtn / arma::repmat(sum(whgtn, 1), 1, whgtn.n_cols);
   return Rcpp::wrap(whgt);
 }

//' @title orthogonal scores algorithn of partial leat squares (opls) projection
//' @description Computes orthogonal socres partial least squares (opls) 
//' projection with the NIPALS algorithm. It allows multiple response variables.
//' Although the main use of the function is for projection, it also retrieves 
//' regression coefficients. NOTE: For internal use only!
//' @usage 
//' opls_for_projection(X, Y, ncomp, scale,
//'                     maxiter, tol,
//'                     pcSelmethod = "var",
//'                     pcSelvalue = 0.01, 
//'                     algorithm = "pls", 
//'                     xls_min_w = 3, 
//'                     xls_max_w = 15)
//' @param X a matrix of predictor variables.
//' @param Y a matrix of either a single or multiple response variables.
//' @param ncomp the number of pls components.
//' @param scale logical indicating whether \code{X} must be scaled.
//' @param maxiter maximum number of iterations.
//' @param tol limit for convergence of the algorithm in the nipals algorithm.
//' @param pcSelmethod if \code{regression = TRUE}, the method for selecting the 
//' number of components. 
//' Options are: \code{'manual'}, \code{'cumvar'} (for selecting the number of 
//' principal components based on a given  cumulative amount of explained 
//' variance) and \code{'var'} (for selecting the number of principal components 
//' based on a given amount of explained variance). Default is \code{'cumvar'}.
//' @param pcSelvalue a numerical value that complements the selected method 
//' (\code{pcSelmethod}). 
//' If \code{'cumvar'} is chosen (default), \code{pcSelvalue} must be a value 
//' (larger than 0 and below 1) indicating the maximum amount of cumulative 
//' variance that the retained components should explain. Default is 0.99. 
//' If \code{'var'} is chosen, \code{pcSelvalue} must be a value (larger than 0 
//' and below 1) indicating that components that explain (individually) 
//' a variance lower than this threshold must be excluded. If \code{'manual'} 
//' is chosen, \code{pcSelvalue} has no effect and the number of components 
//' retrieved are the one specified in \code{ncomp}.
//' @param algorithm (for weights computation) a character string indicating 
//' what method to use. Options are:
//' \code{'pls'} for pls (using covariance between X and Y), 
//' \code{'mpls'} for modified pls (using correlation between X and Y) or
//' \code{'xls'} for extended pls (as implemented in BUCHI NIRWise PLUS software).
//' @param xls_min_w (for weights computation) an integer indicating the minimum window size for the "xls"
//' method. Only used if \code{algorithm = 'xls'}. Default is 3 (as in BUCHI NIRWise PLUS software).
//' @param xls_max_w (for weights computation) an integer indicating the maximum window size for the "xls"
//' method. Only used if \code{algorithm = 'xls'}. Default is 15 (as in BUCHI NIRWise PLUS software).
//' @return a list containing the following elements:
//' \itemize{
//' \item{\code{coefficients}: the matrix of regression coefficients.}
//' \item{\code{bo}: a matrix of one row containing the intercepts for 
//' each component.}
//' \item{\code{scores}: the matrix of scores.}
//' \item{\code{X_loadings}: the matrix of X loadings.}
//' \item{\code{Y_loadings}: the matrix of Y loadings.}
//' \item{\code{projection_mat}: the projection matrix.}
//' \item{\code{Y}: the \code{Y} input.}
//' \item{\code{variance}: a \code{list} conating two objects: \code{x_var} 
//' and \code{y_var}. 
//' These objects contain information on the explained variance for the \code{X} 
//' and \code{Y} matrices respectively.}
//' \item{\code{transf}: a \code{list} conating two objects: \code{Xcenter} 
//' and \code{Xscale}}. 
//' \item{\code{weights}: the matrix of wheights.}
//' }
//' @author Leonardo Ramirez-Lopez
//' @noRd
//' @keywords internal 
//' @useDynLib resemble
// [[Rcpp::export]]
List opls_for_projection(arma::mat X, 
                         arma::mat Y, 
                         int ncomp,
                         bool scale,            
                         double maxiter,
                         double tol,
                         String pcSelmethod = "var",
                         double pcSelvalue = 0.01, 
                         String algorithm = "pls", 
                         const int xls_min_w = 3, 
                         const int xls_max_w = 15) {
  
  int ny = Y.n_cols;
  int nynf = ncomp * Y.n_cols;
  
  arma::mat weights = arma::zeros(ncomp, X.n_cols);
  arma::mat scores = arma::zeros(X.n_rows, ncomp);
  arma::mat Xloadings = arma::zeros(ncomp, X.n_cols);
  arma::mat Yloadings = arma::zeros(ncomp, ny);
  arma::mat coefficients = arma::zeros(X.n_cols, nynf);
  arma::mat bo = arma::zeros(ny, ncomp);
  arma::mat explained_var = arma::zeros(3, ncomp);
  arma::mat yex = arma::zeros(ny, ncomp);
  arma::mat Xscale;
  arma::mat x_scale_vec;
  arma::mat x_center_vec;
  arma::mat Xz = X;
  
  if (scale) {
    Xscale = arma::repmat(Rcpp::as<arma::mat>(get_column_sds(Xz)), Xz.n_rows, 1);
    Xz = Xz / Xscale;
    x_scale_vec = Xscale.row(0);
  } else {
    x_scale_vec = arma::ones<arma::rowvec>(X.n_cols);
  }
  x_center_vec = Rcpp::as<arma::mat>(get_column_means(Xz));
  Xz = Xz - arma::repmat(x_center_vec, Xz.n_rows, 1);
  
  arma::mat Xpls = Xz;
  arma::mat Ypls = Y;
  double xvar;
  
  // Variance of Xpls
  xvar = overall_var(Xpls)(0);
  
  // Matrices to declare
  arma::mat iypls;
  arma::mat Yplsb;
  arma::vec largest_sd_col;
  int j;
  bool keepg;
  arma::mat previous_ts = arma::zeros(Xz.n_rows, 1);
  arma::mat lb;
  arma::mat ts;
  arma::mat w;
  arma::mat p;
  arma::mat q;
  double ireconstructed_var;
  arma::mat cx;
  arma::mat cy;
  arma::mat projection_matrix;
  
  int ith_comp = 0;
  
  // ==========================================================================
  // NIPALS algorithms (pls, mpls, xls)
  // ==========================================================================
  if (algorithm == "pls" || algorithm == "mpls" || algorithm == "xls") {
    for (int i = 0; i < ncomp; i++) {
      Yplsb = Ypls;
      // Select the Y variable with the largest standard deviation
      largest_sd_col = get_col_largest_sd(Ypls);
      iypls = Ypls.col(largest_sd_col[0]);
      previous_ts.fill(0);
      
      j = 0;
      keepg = true;
      
      while (keepg) {
        if (j > 0) {
          previous_ts = ts;
        }
        w = get_weights(Xpls, iypls, algorithm, xls_min_w, xls_max_w);
        ts = Xpls * w;
        p = (trans(Xpls) * ts) / repmat((trans(ts) * ts), Xpls.n_cols, 1);
        q = (trans(Yplsb) * ts) / repmat((trans(ts) * ts), Yplsb.n_cols, 1);
        iypls = (Yplsb * q) / repmat((trans(q) * q), Xpls.n_rows, 1);
        lb = abs(sum((ts - previous_ts) / ts));
        keepg = lb[0] > tol;
        j = j + 1;
        if (maxiter <= j) {
          keepg = false;
        }
      }
      
      // Step 4: The residual matrix of X is computed
      cx = ts * trans(p);
      Xpls = Xpls - cx;
      // The vector of residuals of Y is also computed
      cy = ts * trans(q);
      Ypls = Ypls - cy;
      
      // Save the matrices corresponding to the loadings and scores
      weights.row(i) = trans(w);
      scores.col(i) = ts;
      Xloadings.row(i) = trans(p);
      Yloadings.row(i) = trans(q);
      
      ireconstructed_var = overall_var(cx)(0);
      explained_var(0, i) = ireconstructed_var;
      explained_var(1, i) = explained_var(0, i) / xvar;
      explained_var(2, i) = sum(explained_var.row(0)) / xvar;
      
      ith_comp = ith_comp + 1;
      
      if (pcSelmethod != "manual") {
        if (pcSelmethod == "var" || pcSelmethod == "cumvar") {
          bool chk;
          if (pcSelmethod == "cumvar") {
            chk = explained_var(2, i) > pcSelvalue;
          } else {
            chk = explained_var(1, i) < pcSelvalue;
          }
          if (chk) {
            ncomp = ith_comp - 1;
            ith_comp = ith_comp - 2;
            if (i == 0 && pcSelmethod == "var") {
              throw std::invalid_argument("With the current value in the 'pc_selection' argument, no components are selected. Try another value.");
            }
            break;
          }
        }
      }
    }
  }
  
  // ==========================================================================
  // SIMPLS algorithm
  // ==========================================================================
  if (algorithm == "simpls") {
    // Cross-product matrix
    arma::mat S = trans(Xz) * Y;
    
    // Orthonormal basis for deflation
    arma::mat V(Xz.n_cols, ncomp, arma::fill::zeros);
    
    arma::vec r, t, p_load, v;
    arma::mat q_mat;
    double tt;
    double cumulative_var;
    
    for (int i = 0; i < ncomp; i++) {
      // Weight vector: dominant left singular vector of S
      if (ny == 1) {
        r = S.col(0);
      } else {
        arma::mat U;
        arma::vec s;
        arma::mat Vt;
        arma::svd_econ(U, s, Vt, S, "left");
        r = U.col(0);
      }
      r = r / arma::norm(r);
      
      // Scores
      t = Xz * r;
      tt = arma::as_scalar(trans(t) * t);
      
      // X-loadings (for deflation orthogonalisation)
      p_load = (trans(Xz) * t) / tt;
      
      // Y-loadings
      q_mat = (trans(Y) * t) / tt;
      
      // Modified Gram-Schmidt orthogonalisation
      v = p_load;
      for (int k = 0; k < i; k++) {
        v = v - V.col(k) * arma::as_scalar(trans(V.col(k)) * v);
      }
      v = v / arma::norm(v);
      V.col(i) = v;
      
      // Deflate cross-product matrix
      S = S - v * (trans(v) * S);
      
      // Store results
      weights.row(i) = trans(r);
      scores.col(i) = t;
      Yloadings.row(i) = trans(q_mat);
      
      // Compute explained variance for SIMPLS
      arma::mat scores_so_far = scores.cols(0, i);
      arma::mat TtT = trans(scores_so_far) * scores_so_far;
      arma::mat Xloadings_temp = arma::solve(TtT, trans(scores_so_far) * Xz);
      arma::mat Xrec = scores_so_far * Xloadings_temp;
      cumulative_var = overall_var(Xrec)(0);
      
      if (i == 0) {
        explained_var(0, i) = cumulative_var;
      } else {
        explained_var(0, i) = cumulative_var - sum(explained_var.row(0).cols(0, i - 1));
      }
      explained_var(1, i) = explained_var(0, i) / xvar;
      explained_var(2, i) = cumulative_var / xvar;
      
      ith_comp = ith_comp + 1;
      
      // Early stopping - same logic as NIPALS
      if (pcSelmethod != "manual") {
        if (pcSelmethod == "var" || pcSelmethod == "cumvar") {
          bool chk;
          if (pcSelmethod == "cumvar") {
            chk = explained_var(2, i) > pcSelvalue;
          } else {
            chk = explained_var(1, i) < pcSelvalue;
          }
          if (chk) {
            ncomp = ith_comp - 1;
            ith_comp = ith_comp - 2;
            if (i == 0 && pcSelmethod == "var") {
              throw std::invalid_argument("With the current value in the 'pc_selection' argument, no components are selected. Try another value.");
            }
            break;
          }
        }
      }
    }
    
    // Recompute X-loadings post-hoc for all computed components
    int actual_ncomp;
    if (pcSelmethod != "manual" && (pcSelmethod == "var" || pcSelmethod == "cumvar")) {
      actual_ncomp = ncomp + 1;
      if (actual_ncomp > (int)scores.n_cols) {
        actual_ncomp = scores.n_cols;
      }
    } else {
      actual_ncomp = ncomp;
    }
    
    arma::mat TtT = trans(scores.cols(0, actual_ncomp - 1)) * scores.cols(0, actual_ncomp - 1);
    Xloadings.rows(0, actual_ncomp - 1) = arma::solve(TtT, trans(scores.cols(0, actual_ncomp - 1)) * Xz);
  }
  
  // ==========================================================================
  // Common post-processing (same as original)
  // ==========================================================================
  arma::uvec pc_indices;
  if (pcSelmethod != "manual") {
    if (pcSelmethod == "var" || pcSelmethod == "cumvar") {
      if (pcSelmethod == "var") {
        pc_indices = find(explained_var.row(1) >= pcSelvalue);
      } else {
        pc_indices = find(explained_var.row(2) <= pcSelvalue && explained_var.row(2) > 0);
        pc_indices = pc_indices + 1;
        pc_indices.insert_rows(0, 1);
        ncomp = ncomp + 1;
      }
      weights = weights.rows(pc_indices);
      // Keep all the coefficients for all the Ys
      coefficients = coefficients.cols(0, (ncomp * ny) - 1);
      bo = bo.cols(pc_indices);
      scores = scores.cols(pc_indices);
      Xloadings = Xloadings.rows(pc_indices);
      Yloadings = Yloadings.rows(pc_indices);
      explained_var = explained_var.cols(pc_indices);
      yex = yex.cols(pc_indices);
    }
  }
  
  // Compute projection matrix
  if (algorithm == "simpls") {
    projection_matrix = trans(weights);
  } else {
    projection_matrix = trans(weights) * arma::solve(Xloadings * trans(weights), arma::eye(Xloadings.n_rows, Xloadings.n_rows));
  }
  
  // Compute Y explained variance
  arma::mat yexi;
  arma::mat cop;
  for (int i = 0; i < ny; i++) {
    yexi = scores % arma::repmat(trans(Yloadings.col(i)), scores.n_rows, 1);
    cop = pow(arma::cor(Y.col(i), yexi.col(0)), 2);
    yex(i, 0) = cop(0, 0);
    for (int j = 1; j < ncomp; j++) {
      yexi.col(j) = yexi.col(j - 1) + yexi.col(j);
      cop = arma::cor(Y.col(i), yexi.col(j));
      yex(i, j) = pow(cop(0, 0), 2);
    }
  }
  
  // Compute coefficients and intercepts
  arma::mat ymean = arma::mean(Y);
  arma::vec ymean_vec = arma::vectorise(ymean);
  arma::mat y_hat_mean;
  arma::vec y_hat_mean_vec;
  int idx = 0;
  for (int k = 0; k < ny; k++) {
    arma::mat kth_loading = Yloadings.col(k);
    for (int j = 0; j < ncomp; j++) {
      coefficients.col(idx) = projection_matrix.cols(0, j) * kth_loading.rows(0, j);
      y_hat_mean = x_center_vec * coefficients.col(idx);
      y_hat_mean_vec = arma::vectorise(y_hat_mean);
      bo(k, j) = ymean_vec(k) - y_hat_mean_vec(0);
      idx = idx + 1;
    }
  }
  
  return Rcpp::List::create(
    Rcpp::Named("ncomp") = ncomp,
    Rcpp::Named("coefficients") = coefficients,
    Rcpp::Named("bo") = bo,
    Rcpp::Named("scores") = scores,
    Rcpp::Named("X_loadings") = Xloadings,
    Rcpp::Named("Y_loadings") = Yloadings,
    Rcpp::Named("projection_mat") = projection_matrix,
    Rcpp::Named("Y") = Y,
    Rcpp::Named("variance") = Rcpp::List::create(
      Rcpp::Named("original_x_var") = xvar,
      Rcpp::Named("x_var") = explained_var,
      Rcpp::Named("y_var") = yex
    ),
    Rcpp::Named("transf") = Rcpp::List::create(
      Rcpp::Named("Xcenter") = x_center_vec,
      Rcpp::Named("Xscale") = x_scale_vec
    ),
    Rcpp::Named("weights") = weights
  );
}
//' @title orthogonal scores algorithn of partial leat squares (opls_get_all)
//' @description Computes orthogonal socres partial least squares (opls_get_all) 
//' regressions with the NIPALS algorithm. It retrives a comprehensive set of
//' pls outputs (e.g. vip and sensivity radius). It allows multiple response 
//' variables. NOTE: For internal use only!
//' @usage 
//' opls_get_all(X, 
//'              Y, 
//'              ncomp, 
//'              scale, 
//'              maxiter, 
//'              tol, 
//'              algorithm = "pls", 
//'              xls_min_w = 3, 
//'              xls_max_w = 15)
//' @param X a matrix of predictor variables.
//' @param Y a matrix of either a single or multiple response variables.
//' @param ncomp the number of pls components.
//' @param scale logical indicating whether \code{X} must be scaled.
//' @param maxiter maximum number of iterations.
//' @param tol limit for convergence of the algorithm in the nipals algorithm.
//' @param algorithm (for weights computation) a character string indicating 
//' what method to use. Options are:
//' \code{'pls'} for pls (using covariance between X and Y), 
//' \code{'mpls'} for modified pls (using correlation between X and Y) or
//' \code{'xls'} for extended pls (as implemented in BUCHI NIRWise PLUS software).
//' @param xls_min_w (for weights computation) an integer indicating the minimum window size for the "xls"
//' method. Only used if \code{algorithm = 'xls'}. Default is 3 (as in BUCHI NIRWise PLUS software).
//' @param xls_max_w (for weights computation) an integer indicating the maximum window size for the "xls"
//' method. Only used if \code{algorithm = 'xls'}. Default is 15 (as in BUCHI NIRWise PLUS software).
//' @return a list containing the following elements:
//' \itemize{
//' \item{\code{ncomp}: the number of components used.}
//' \item{\code{coefficients}: the matrix of regression coefficients.}
//' \item{\code{bo}: a matrix of one row containing the intercepts for each component.}
//' \item{\code{scores}: the matrix of scores.}
//' \item{\code{X_loadings}: the matrix of X loadings.}
//' \item{\code{Y_loadings}: the matrix of Y loadings.}
//' \item{\code{vip}: the projection matrix.}
//' \item{\code{selectivity_ratio}: the matrix of selectivity ratio (see Rajalahti, Tarja, et al. 2009).}
//' \item{\code{Y}: the \code{Y} input.}
//' \item{\code{variance}: a \code{list} conating two objects: \code{x_var} and \code{y_var}. 
//' These objects contain information on the explained variance for the \code{X} and \code{Y} matrices respectively.}
//' \item{\code{transf}: a \code{list} conating two objects: \code{Xcenter} and \code{Xscale}}. 
//' \item{\code{weights}: the matrix of wheights.}
//' } 
//' @author Leonardo Ramirez-Lopez
//' @keywords internal 
//' @noRd
//' @useDynLib resemble
// [[Rcpp::export]]
List opls_get_all(
    arma::mat X, 
    arma::mat Y, 
    int ncomp,
    bool scale,            
    double maxiter,
    double tol, 
    String algorithm = "pls", 
    const int xls_min_w = 3, 
    const int xls_max_w = 15
) {
  
  int ny = Y.n_cols;
  int nynf = ncomp * Y.n_cols;
  
  arma::mat weights = arma::zeros(ncomp, X.n_cols);
  arma::mat scores = arma::zeros(X.n_rows, ncomp);
  arma::mat Xloadings = arma::zeros(ncomp, X.n_cols);
  arma::mat Yloadings = arma::zeros(ncomp, ny);
  arma::mat coefficients = arma::zeros(X.n_cols, nynf);
  arma::mat bo = arma::zeros(ny, ncomp);
  arma::mat explained_var = arma::zeros(3, ncomp);
  arma::mat yex = arma::zeros(ny, ncomp);
  arma::mat Xscale;
  arma::mat x_scale_vec;
  arma::mat x_center_vec;
  arma::mat Xz = X;
  
  if (scale) {
    Xscale = arma::repmat(Rcpp::as<arma::mat>(get_column_sds(Xz)), Xz.n_rows, 1);
    Xz = Xz / Xscale;
    x_scale_vec =  Xscale.row(0);
  }
  x_center_vec = Rcpp::as<arma::mat>(get_column_means(Xz));
  Xz = Xz - arma::repmat(x_center_vec, Xz.n_rows, 1);
  
  arma::mat Xpls = Xz;
  arma::mat Ypls = Y;
  
  //variance of Xpls
  double xvar = overall_var(Xpls)(0);
  
  // matrices to declare
  arma::mat iypls;
  arma::mat Yplsb;
  arma::vec lagest_sd_col;
  int j;
  bool keepg;  
  arma::mat previous_ts = arma::zeros(Xz.n_rows, 1);
  arma::mat lb;
  arma::mat cr;
  arma::mat ts;
  arma::mat w;
  arma::mat p;
  arma::mat q;
  double ireconstructed_var;
  arma::mat cx;
  arma::mat cy;
  arma::mat projection_matrix;
  arma::mat sratio = arma::zeros(weights.n_rows, weights.n_cols);
  
  
  // ==========================================================================
  // NIPALS algorithms (pls, mpls, xls)
  // ==========================================================================
  if (algorithm == "pls" || algorithm == "mpls" || algorithm == "xls") {
    for (int i = 0; i < ncomp; i++) {
      Yplsb = Ypls;
      // Select the Y variable with the largest standard deviation
      lagest_sd_col = get_col_largest_sd(Ypls);
      iypls = Ypls.col(lagest_sd_col[0]);
      previous_ts.fill(0);
      
      j = 0;
      keepg = true;
      while (keepg) {
        if (j > 0) {
          previous_ts = ts;
        }
        // //Step 1: Compute a vector of loading weights...
        // // 1.1 Compute the 'scaling factor'
        // cr = sqrt(trans(iypls) * Xpls * trans(Xpls) * iypls);
        // // 1.2 The weights are computed as the cross product of
        // // X0 and Y0 divided by the 'scaling factor'...
        // w = (trans(Xpls) * iypls) / repmat(cr, Xpls.n_cols, 1);
        w = get_weights(Xpls, iypls, algorithm, xls_min_w, xls_max_w);
        // Step 2: Compute the scores...
        ts = Xpls * w;
        // Step 3: Compute the X-loadings (p) and the Y-loadings (q)...
        p = (trans(Xpls) * ts) / repmat((trans(ts) * ts), Xpls.n_cols, 1);
        q = (trans(Yplsb) * ts) / repmat((trans(ts) * ts), Yplsb.n_cols, 1);
        iypls = (Yplsb * q) / repmat((trans(q) * q), Xpls.n_rows, 1) ;
        lb = abs(sum((ts - previous_ts) / ts));
        keepg = lb[0] > tol;
        j = j + 1;
        if(maxiter <= j) {
          keepg = false;
        }
      }
      // Step 4: The residual matrix
      // of X is finally computed and...
      cx = ts * trans(p);
      Xpls = Xpls - cx;
      // ... the vector of residuals of Y is also computed
      cy = ts * trans(q);
      Ypls = Ypls - cy;
      // save the matrices corresponding to the loadings
      // and scores..
      weights.row(i) = trans(w);
      scores.col(i) = ts;
      Xloadings.row(i) = trans(p);
      Yloadings.row(i) = trans(q);
      
      ireconstructed_var = overall_var(cx)(0);
      explained_var(0,i) = ireconstructed_var;
      explained_var(1,i) = explained_var(0,i) / xvar;
      explained_var(2,i) = sum(explained_var.row(0)) / xvar;
    }
  }
  
  // ==========================================================================
  // SIMPLS algorithm
  // ==========================================================================
  if (algorithm == "simpls") {
    // Cross-product matrix
    arma::mat S = trans(Xz) * Y;
    
    // Orthonormal basis for deflation
    arma::mat V(Xz.n_cols, ncomp, arma::fill::zeros);
    
    arma::vec r, t, p_load, v;
    arma::mat q_mat;
    double tt;
    
    // For explained variance calculation in SIMPLS
    arma::mat Xpls_simpls = Xz;
    
    for (int i = 0; i < ncomp; i++) {
      // Weight vector: dominant left singular vector of S
      if (ny == 1) {
        r = S.col(0);
      } else {
        arma::mat U;
        arma::vec s;
        arma::mat Vt;
        arma::svd_econ(U, s, Vt, S, "left");
        r = U.col(0);
      }
      r = r / arma::norm(r);
      
      // Scores
      t = Xz * r;
      tt = arma::as_scalar(trans(t) * t);
      
      // X-loadings (for deflation orthogonalisation)
      p_load = (trans(Xz) * t) / tt;
      
      // Y-loadings
      q_mat = (trans(Y) * t) / tt;
      
      // Modified Gram-Schmidt orthogonalisation
      v = p_load;
      for (int k = 0; k < i; k++) {
        v = v - V.col(k) * arma::as_scalar(trans(V.col(k)) * v);
      }
      v = v / arma::norm(v);
      V.col(i) = v;
      
      // Deflate cross-product matrix
      S = S - v * (trans(v) * S);
      
      // Store results
      weights.row(i) = trans(r);
      scores.col(i) = t;
      Yloadings.row(i) = trans(q_mat);
      
      // Explained variance (compute reconstruction from cumulative scores)
      // For SIMPLS, we compute explained variance differently since X is not deflated
      arma::mat scores_so_far = scores.cols(0, i);
      arma::mat TtT_inv_temp = arma::inv(trans(scores_so_far) * scores_so_far);
      arma::mat Xloadings_temp = trans(trans(Xz) * scores_so_far * TtT_inv_temp);
      arma::mat Xrec = scores_so_far * Xloadings_temp;
      double cumulative_var = overall_var(Xrec)(0);
      
      if (i == 0) {
        explained_var(0, i) = cumulative_var;
      } else {
        explained_var(0, i) = cumulative_var - sum(explained_var.row(0).cols(0, i - 1));
      }
      explained_var(1, i) = explained_var(0, i) / xvar;
      explained_var(2, i) = cumulative_var / xvar;
    }
    
    // Recompute X-loadings for proper reconstruction
    arma::mat TtT_inv = arma::inv(trans(scores) * scores);
    Xloadings = trans(trans(Xz) * scores * TtT_inv);
  }
  
  // ==========================================================================
  // Common post-processing (shared by all algorithms)
  // ==========================================================================
  
  projection_matrix = trans(weights) * arma::solve(Xloadings * trans(weights), arma::eye(Xloadings.n_rows, Xloadings.n_rows));
  
  arma::mat yexi;
  arma::mat cop;
  for (int i = 0; i < ny; i++) {
    yexi = scores % arma::repmat(trans(Yloadings.col(i)), scores.n_rows, 1) ;
    cop = pow(arma::cor(Y.col(i), yexi.col(0)), 2);
    //double(*cop2) = reinterpret_cast <double(*)> (cop); //does not work
    yex(i,0) = cop(0,0);
    for(int j = 1; j < ncomp; j++){
      yexi.col(j) = yexi.col(j-1) + yexi.col(j);
      cop = arma::cor(Y.col(i), yexi.col(j));
      yex(i,j) = pow(cop(0,0), 2);
    }
  }
  
  arma::mat ymean = arma::mean(Y);
  arma::vec ymean_vec = arma::vectorise(ymean);
  arma::mat y_hat_mean;
  arma::mat y_hat_mean_vec;
  arma::mat wtp;
  arma::mat ttp;
  arma::mat ptp;
  arma::mat expvar;
  arma::mat resvar;
  arma::mat bd;
  
  int idx = 0;
  for (int k = 0; k < ny; k++) {
    arma::mat jth_loading = Yloadings.col(k);
    for (int j = 0; j < ncomp; j++) {
      // Estimate coefficients
      coefficients.col(idx) = projection_matrix.cols(0,j) * jth_loading.rows(0,j);
      // computation for target projection and its selectivity ratio
      // Rajalahti, Tarja, et al. 
      // Biomarker discovery in mass spectral profiles by means of selectivity ratio plot.
      // Chemometrics and Intelligent Laboratory Systems 95.1 (2009): 35-48.
      bd = sqrt(trans(coefficients.col(idx)) * coefficients.col(idx));
      wtp = coefficients.col(idx)/repmat(bd, coefficients.n_rows, 1);
      ttp = Xz * wtp;
      ptp = trans(ttp) * Xz;
      ptp = ptp/repmat(trans(ttp) * ttp, 1, coefficients.n_rows);
      expvar = ttp * ptp;
      resvar = arma::var(Xz - expvar, 0, 0);
      expvar = arma::var(expvar, 0, 0);
      sratio.row(j) = expvar/resvar ;
      // compute the intercept
      y_hat_mean = x_center_vec * coefficients.col(idx);
      y_hat_mean_vec = arma::vectorise(y_hat_mean);
      bo(k, j) = ymean_vec(k) - y_hat_mean_vec(0);
      idx = idx + 1;
    }
  }
  arma::mat ss = arma::zeros(Yloadings.n_rows, Yloadings.n_cols);
  arma::mat wss = arma::zeros(Yloadings.n_rows, Yloadings.n_cols);
  arma::mat ssw = arma::zeros(scores.n_rows, Yloadings.n_rows);
  arma::mat vip = arma::zeros(weights.n_rows, weights.n_cols);
  arma::mat ss1 = pow(Yloadings, 2);
  arma::mat ss2 = get_column_sums(pow(scores, 2));
  ss = ss1 % ss2;
  
  arma::mat sqweights = trans(pow(weights, 2));
  wss = get_column_sums(sqweights);
  ssw = sqweights % trans(arma::repmat(ss/wss, 1, weights.n_cols));
  arma::mat cssw = arma::zeros(weights.n_cols, weights.n_rows);
  double sclr = ssw.n_rows;
  int it = ssw.n_rows;
  for(int i = 0; i < it; i++) {
    cssw.row(i) = cumsum(ssw.row(i));
  }
  cssw = cssw * sclr;
  vip = pow(cssw / arma::repmat(trans(cumsum(ss)), weights.n_cols, 1), 0.5);
  //FIXME: For every Y store the coefficients independently
  return Rcpp::List::create(
    Rcpp::Named("ncomp") = ncomp,
    Rcpp::Named("coefficients") = coefficients,
    Rcpp::Named("bo") = bo,
    Rcpp::Named("scores") = scores,
    Rcpp::Named("X_loadings") = Xloadings,
    Rcpp::Named("Y_loadings") = Yloadings,
    Rcpp::Named("projection_mat") = projection_matrix,
    Rcpp::Named("vip") = vip,
    Rcpp::Named("selectivity_ratio") = trans(sratio),
    Rcpp::Named("Y") = Y,
    Rcpp::Named("variance") = Rcpp::List::create(
      Rcpp::Named("original_x_var") = xvar,
      Rcpp::Named("x_var") = explained_var,
      Rcpp::Named("y_var") = yex
    ),
    Rcpp::Named("transf") = Rcpp::List::create(
      Rcpp::Named("Xcenter") = x_center_vec,
      Rcpp::Named("Xscale") = x_scale_vec
    ),
    _["weights"] = weights
  );
}

//' @title orthogonal scores algorithn of partial leat squares (opls)
//' @description Computes orthogonal socres partial least squares (opls) 
//' regressions with the NIPALS algorithm. It allows multiple response variables. 
//' It does not return the variance information of the components. NOTE: For 
//' internal use only!
//' @usage 
//' opls(X, 
//'      Y, 
//'      ncomp, 
//'      scale, 
//'      maxiter, 
//'      tol, 
//'      algorithm = "pls", 
//'      xls_min_w = 3, 
//'      xls_max_w = 15)
//'      
//' @param X a matrix of predictor variables.
//' @param Y a matrix of either a single or multiple response variables.
//' @param ncomp the number of pls components.
//' @param scale logical indicating whether \code{X} must be scaled.
//' @param maxiter maximum number of iterations.
//' @param tol limit for convergence of the algorithm in the nipals algorithm.
//' @param algorithm (for weights computation) a character string indicating 
//' what method to use. Options are:
//' \code{'pls'} for pls (using covariance between X and Y), 
//' \code{'mpls'} for modified pls (using correlation between X and Y) or
//' \code{'xls'} for extended pls (as implemented in BUCHI NIRWise PLUS software).
//' @param xls_min_w (for weights computation) an integer indicating the minimum window size for the "xls"
//' method. Only used if \code{algorithm = 'xls'}. Default is 3 (as in BUCHI NIRWise PLUS software).
//' @param xls_max_w (for weights computation) an integer indicating the maximum window size for the "xls"
//' method. Only used if \code{algorithm = 'xls'}. Default is 15 (as in BUCHI NIRWise PLUS software).
//' @return a list containing the following elements:
//' \itemize{
//' \item{\code{coefficients}: the matrix of regression coefficients.}
//' \item{\code{bo}: a matrix of one row containing the intercepts for each component.}
//' \item{\code{scores}: the matrix of scores.}
//' \item{\code{X_loadings}: the matrix of X loadings.}
//' \item{\code{Y_loadings}: the matrix of Y loadings.}
//' \item{\code{projection_mat}: the projection matrix.}
//' \item{\code{Y}: the \code{Y} input.}
//' \item{\code{transf}: a \code{list} conating two objects: \code{Xcenter} and \code{Xscale}}. 
//' \item{\code{weights}: the matrix of wheights.}} 
//' @author Leonardo Ramirez-Lopez
//' @noRd
//' @keywords internal 
//' @useDynLib resemble
// [[Rcpp::export]]
List opls(arma::mat X, 
          arma::mat Y, 
          int ncomp,
          bool scale,            
          double maxiter,
          double tol, 
          String algorithm = "pls", 
          const int xls_min_w = 3, 
          const int xls_max_w = 15) {
  
  int ny = Y.n_cols;
  int nynf = ncomp * Y.n_cols;
  
  arma::mat weights = arma::zeros(ncomp, X.n_cols);
  arma::mat scores = arma::zeros(X.n_rows, ncomp);
  arma::mat Xloadings = arma::zeros(ncomp, X.n_cols);
  arma::mat Yloadings = arma::zeros(ncomp, ny);
  arma::mat coefficients = arma::zeros(X.n_cols, nynf);
  arma::mat bo = arma::zeros(ny, ncomp);
  arma::mat Xscale;
  arma::mat x_scale_vec;
  arma::mat x_center_vec;
  arma::mat Xz = X;
  
  if (scale) {
    Xscale = arma::repmat(Rcpp::as<arma::mat>(get_column_sds(Xz)), Xz.n_rows, 1);
    Xz = Xz / Xscale;
    x_scale_vec =  Xscale.row(0);
  }
  x_center_vec = Rcpp::as<arma::mat>(get_column_means(Xz));
  Xz = Xz - arma::repmat(x_center_vec, Xz.n_rows, 1);
  
  // matrices to declare
  arma::mat Xpls = Xz;
  arma::mat Ypls = Y;
  arma::mat iypls;
  arma::mat Yplsb;
  arma::vec lagest_sd_col;
  int j;
  bool keepg;
  arma::mat previous_ts = arma::zeros(Xz.n_rows, 1);
  arma::mat lb;
  arma::mat cr;
  arma::mat ts;
  arma::mat w;
  arma::mat p;
  arma::mat q;
  arma::mat cx;
  arma::mat cy;
  arma::mat projection_matrix;
  
  // ==========================================================================
  // NIPALS algorithms (pls, mpls, xls)
  // ==========================================================================
  if (algorithm == "pls" || algorithm == "mpls" || algorithm == "xls") {
    for (int i = 0; i < ncomp; i++) {
      Yplsb = Ypls;
      // Select the Y variable with the largest standard deviation
      lagest_sd_col = get_col_largest_sd(Ypls);
      iypls = Ypls.col(lagest_sd_col[0]);
      previous_ts.fill(0);
      
      j = 0;
      keepg = true;
      
      while (keepg) {
        if (j > 0) {
          previous_ts = ts;
        }
        //Step 1: Compute a vector of loading weights...
        // // 1.1 Compute the 'scaling factor'
        // cr = sqrt(trans(iypls) * Xpls * trans(Xpls) * iypls);
        // // 1.2 The weights are computed as the cross product of
        // // X0 and Y0 divided by the 'scaling factor'...
        // w = (trans(Xpls) * iypls) / repmat(cr, Xpls.n_cols, 1);
        w = get_weights(Xpls, iypls, algorithm, xls_min_w, xls_max_w);
        // Step 2: Compute the scores...
        ts = Xpls * w;
        // Step 3: Compute the X-loadings (p) and the Y-loadings (q)...
        p = (trans(Xpls) * ts) / repmat((trans(ts) * ts), Xpls.n_cols, 1);
        q = (trans(Yplsb) * ts) / repmat((trans(ts) * ts), Yplsb.n_cols, 1);
        iypls = (Yplsb * q) / repmat((trans(q) * q), Xpls.n_rows, 1) ;
        lb = abs(sum((ts - previous_ts) / ts));
        keepg = lb[0] > tol;
        j = j + 1;
        if (maxiter <= j) {
          keepg = false;
        }
      }
      // Step 4: The residual matrix
      // of X is finally computed and...
      cx = ts * trans(p) ;
      Xpls = Xpls - cx;
      // ... the vector of residuals of Y is also computed
      cy = ts * trans(q);
      Ypls = Ypls - cy;
      // save the matrices corresponding to the loadings
      // and scores..
      weights.row(i) = trans(w);
      scores.col(i) = ts;
      Xloadings.row(i) = trans(p);
      Yloadings.row(i) = trans(q);
    }
  }
  
  // ==========================================================================
  // SIMPLS algorithm
  // ==========================================================================
  if (algorithm == "simpls") {
    // Cross-product matrix
    arma::mat S = trans(Xz) * Y;
    
    // Orthonormal basis for deflation
    arma::mat V(Xz.n_cols, ncomp, arma::fill::zeros);
    
    arma::vec r, t, p_load, v;
    arma::mat q_mat;
    double tt;
    
    for (int i = 0; i < ncomp; i++) {
      // Weight vector: dominant left singular vector of S
      if (ny == 1) {
        r = S.col(0);
      } else {
        arma::mat U;
        arma::vec s;
        arma::mat Vt;
        arma::svd_econ(U, s, Vt, S, "left");
        r = U.col(0);
      }
      r = r / arma::norm(r);
      
      // Scores
      t = Xz * r;
      tt = arma::as_scalar(trans(t) * t);
      
      // X-loadings (for deflation orthogonalisation)
      p_load = (trans(Xz) * t) / tt;
      
      // Y-loadings
      q_mat = (trans(Y) * t) / tt;
      
      // Modified Gram-Schmidt orthogonalisation
      v = p_load;
      for (int k = 0; k < i; k++) {
        v = v - V.col(k) * arma::as_scalar(trans(V.col(k)) * v);
      }
      v = v / arma::norm(v);
      V.col(i) = v;
      
      // Deflate cross-product matrix
      S = S - v * (trans(v) * S);
      
      // Store results
      weights.row(i) = trans(r);
      scores.col(i) = t;
      Yloadings.row(i) = trans(q_mat);
    }
    
    // Recompute X-loadings for proper reconstruction
    arma::mat TtT_inv = arma::inv(trans(scores) * scores);
    Xloadings = trans(trans(Xz) * scores * TtT_inv);
  }
  
  // ==========================================================================
  // Common post-processing (shared by all algorithms)
  // ==========================================================================
  projection_matrix = trans(weights) * arma::solve(Xloadings * trans(weights), arma::eye(Xloadings.n_rows, Xloadings.n_rows));
  
  arma::mat ymean = arma::mean(Y);
  arma::vec ymean_vec = arma::vectorise(ymean);
  arma::mat y_hat_mean;
  arma::vec y_hat_mean_vec;
  int idx = 0;
  for (int k = 0; k < ny; k++) {
    arma::mat jth_loading = Yloadings.col(k);
    for (int j = 0; j < ncomp; j++) {
      coefficients.col(idx) = projection_matrix.cols(0, j) * jth_loading.rows(0, j);
      y_hat_mean = x_center_vec * coefficients.col(idx);
      y_hat_mean_vec = arma::vectorise(y_hat_mean);
      bo(k, j) = ymean_vec(k) - y_hat_mean_vec(0);
      idx = idx + 1;
    }
  }
  //FIXME: For every Y store the coefficients independently
  return Rcpp::List::create(
    Rcpp::Named("ncomp") = ncomp,
    Rcpp::Named("coefficients") = coefficients,
    Rcpp::Named("bo") = bo,
    Rcpp::Named("scores") = scores,
    Rcpp::Named("X_loadings") = Xloadings,
    Rcpp::Named("Y_loadings") = Yloadings,
    Rcpp::Named("projection_mat") = projection_matrix,
    Rcpp::Named("Y") = Y,
    Rcpp::Named("transf") = Rcpp::List::create(
      Rcpp::Named("Xcenter") = x_center_vec,
      Rcpp::Named("Xscale") = x_scale_vec
    ),
    _["weights"] = weights
  );
}

//' @title Fast orthogonal scores algorithm of partial least squares (PLS)
//' @description Computes orthogonal scores partial least squares (PLS) 
//' regression using either NIPALS or SIMPLS algorithm. Supports multiple 
//' response variables. In contrast to \code{opls}, this function omits 
//' auxiliary outputs (e.g. scores, explained variance) not required for 
//' local regression. For internal use only.
//' @usage 
//' opls_get_basics(X, Y, ncomp, scale, 
//'                 maxiter, tol, 
//'                 algorithm = "pls", 
//'                 xls_min_w = 3, 
//'                 xls_max_w = 15)
//' @param X a matrix of predictor variables.
//' @param Y a matrix of either a single or multiple response variables.
//' @param ncomp the number of PLS components.
//' @param scale logical indicating whether \code{X} must be scaled.
//' @param maxiter maximum number of iterations (only used for NIPALS-based 
//' algorithms: \code{'pls'}, \code{'mpls'}, \code{'xls'}).
//' @param tol convergence tolerance for the NIPALS algorithm (only used for 
//' NIPALS-based algorithms).
//' @param algorithm a character string indicating the PLS algorithm to use:
//' \itemize{
//'   \item{\code{'pls'}: standard PLS using covariance between X and Y for 
//'     weight computation (NIPALS algorithm).}
//'   \item{\code{'mpls'}: modified PLS using correlation between X and Y for 
//'     weight computation (NIPALS algorithm). See Shenk and Westerhaus (1991).}
//'   \item{\code{'xls'}: extended PLS as implemented in BUCHI NIRWise PLUS 
//'     software (NIPALS algorithm).}
//'   \item{\code{'simpls'}: SIMPLS algorithm (de Jong, 1993). Computationally 
//'     faster as it avoids iterative X deflation. Parameters \code{maxiter}, 
//'     \code{tol}, \code{xls_min_w}, and \code{xls_max_w} are ignored.}
//' }
//' @param xls_min_w an integer indicating the minimum window size for the 
//' \code{'xls'} method. Only used if \code{algorithm = 'xls'}. Default is 3.
//' @param xls_max_w an integer indicating the maximum window size for the 
//' \code{'xls'} method. Only used if \code{algorithm = 'xls'}. Default is 15.
//' @return a list containing:
//' \itemize{
//'   \item{\code{ncomp}: the number of PLS components.}
//'   \item{\code{coefficients}: the matrix of regression coefficients.}
//'   \item{\code{bo}: a matrix containing the intercepts for each component.}
//'   \item{\code{X_loadings}: the matrix of X loadings.}
//'   \item{\code{Y_loadings}: the matrix of Y loadings.}
//'   \item{\code{projection_mat}: the projection matrix for computing scores 
//'     from new data.}
//'   \item{\code{transf}: a list containing:
//'     \itemize{
//'       \item{\code{Xcenter}: row vector of column means used for centering.}
//'       \item{\code{Xscale}: row vector of column standard deviations used 
//'         for scaling (ones if \code{scale = FALSE}).}
//'     }
//'   }
//'   \item{\code{weights}: the matrix of PLS weights.}
//' }
//' @references
//' de Jong, S. (1993). SIMPLS: An alternative approach to partial least 
//' squares regression. Chemometrics and Intelligent Laboratory Systems, 
//' 18(3), 251-263.
//' 
//' Shenk, J.S., & Westerhaus, M.O. (1991). Populations structuring of near 
//' infrared spectra and modified partial least squares regression. Crop 
//' Science, 31(6), 1548-1555.
//' @author Leonardo Ramirez-Lopez
//' @keywords internal 
//' @noRd
//' @useDynLib resemble
// [[Rcpp::export]]
List opls_get_basics(
    arma::mat X, 
    arma::mat Y, 
    int ncomp,
    bool scale,            
    double maxiter,
    double tol, 
    String algorithm = "pls", 
    const int xls_min_w = 3, 
    const int xls_max_w = 15
) {
  int ny = Y.n_cols;
  int nynf = ncomp * ny;
  
  arma::mat weights(ncomp, X.n_cols, arma::fill::zeros);
  arma::mat scores(X.n_rows, ncomp, arma::fill::zeros);
  arma::mat Xloadings(ncomp, X.n_cols, arma::fill::zeros);
  arma::mat Yloadings(ncomp, ny, arma::fill::zeros);
  arma::mat coefficients(X.n_cols, nynf, arma::fill::zeros);
  arma::mat bo(ny, ncomp, arma::fill::zeros);
  arma::rowvec x_scale_vec;
  arma::rowvec x_center_vec;
  arma::mat Xz = X;
  
  // Scale-then-center (original order)
  if (scale) {
    x_scale_vec = arma::stddev(Xz, 0, 0);
    Xz.each_row() /= x_scale_vec;
  } else {
    x_scale_vec = arma::ones<arma::rowvec>(X.n_cols);
  }
  x_center_vec = arma::mean(Xz, 0);
  Xz.each_row() -= x_center_vec;
  
  arma::mat Xpls = Xz;
  arma::mat Ypls = Y;
  arma::mat projection_matrix;
  
  // ==========================================================================
  // NIPALS algorithms (pls, mpls, xls)
  // ==========================================================================
  if (algorithm == "pls" || algorithm == "mpls" || algorithm == "xls") {
    arma::mat ts, w, p;
    double tsts;
    
    if (ny == 1) {
      // Fast path for PLS1 (no iteration needed)
      for (int i = 0; i < ncomp; i++) {
        w = get_weights(Xpls, Ypls.col(0), algorithm, xls_min_w, xls_max_w);
        ts = Xpls * w;
        tsts = arma::as_scalar(trans(ts) * ts);
        p = (trans(Xpls) * ts) / tsts;
        double q_val = arma::as_scalar(trans(Ypls) * ts) / tsts;
        
        Xpls -= ts * trans(p);
        Ypls.col(0) -= ts * q_val;
        
        weights.row(i) = trans(w);
        scores.col(i) = ts;
        Xloadings.row(i) = trans(p);
        Yloadings(i, 0) = q_val;
      }
    } else {
      // General path for PLS2 (requires iteration)
      arma::mat iypls, Yplsb, q;
      arma::vec imsd;
      arma::mat previous_ts(Xz.n_rows, 1, arma::fill::zeros);
      double qtq, lb;
      int j;
      bool keepg;
      
      for (int i = 0; i < ncomp; i++) {
        Yplsb = Ypls;
        imsd = get_col_largest_sd(Ypls);
        iypls = Ypls.col(static_cast<arma::uword>(imsd[0]));
        previous_ts.zeros();
        
        j = 0;
        keepg = true;
        
        while (keepg) {
          if (j > 0) {
            previous_ts = ts;
          }
          w = get_weights(Xpls, iypls, algorithm, xls_min_w, xls_max_w);
          ts = Xpls * w;
          tsts = arma::as_scalar(trans(ts) * ts);
          p = (trans(Xpls) * ts) / tsts;
          q = (trans(Yplsb) * ts) / tsts;
          qtq = arma::as_scalar(trans(q) * q);
          iypls = (Yplsb * q) / qtq;
          lb = arma::as_scalar(arma::abs(arma::sum((ts - previous_ts) / ts)));
          keepg = lb > tol;
          j++;
          if (j >= maxiter) {
            keepg = false;
          }
        }
        
        Xpls -= ts * trans(p);
        Ypls -= ts * trans(q);
        
        weights.row(i) = trans(w);
        scores.col(i) = ts;
        Xloadings.row(i) = trans(p);
        Yloadings.row(i) = trans(q);
      }
    }
  }
  
  // ==========================================================================
  // SIMPLS algorithm
  // ==========================================================================
  // ==========================================================================
  // SIMPLS algorithm
  // ==========================================================================
  if (algorithm == "simpls") {
    arma::mat V(Xz.n_cols, ncomp, arma::fill::zeros);
    arma::vec r, t, p_load, v;
    double tt;
    
    if (ny == 1) {
      // Fast path for single Y
      const arma::vec y = Y.col(0);
      arma::vec S = trans(Xz) * y;
      double q_val, vtS;
      
      for (int i = 0; i < ncomp; i++) {
        r = S / arma::norm(S);
        t = Xz * r;
        tt = arma::as_scalar(trans(t) * t);
        p_load = (trans(Xz) * t) / tt;
        q_val = arma::as_scalar(trans(y) * t) / tt;
        
        // Modified Gram-Schmidt
        v = p_load;
        for (int k = 0; k < i; k++) {
          v -= V.col(k) * arma::as_scalar(trans(V.col(k)) * v);
        }
        v /= arma::norm(v);
        V.col(i) = v;
        
        vtS = arma::as_scalar(trans(v) * S);
        S -= v * vtS;
        
        weights.row(i) = trans(r);
        scores.col(i) = t;
        Yloadings(i, 0) = q_val;
      }
    } else {
      // General path for multiple Y
      arma::mat S = trans(Xz) * Y;
      arma::mat q_mat;
      
      for (int i = 0; i < ncomp; i++) {
        arma::mat U;
        arma::vec s;
        arma::mat Vt;
        arma::svd_econ(U, s, Vt, S, "left");
        r = U.col(0);
        r /= arma::norm(r);
        
        t = Xz * r;
        tt = arma::as_scalar(trans(t) * t);
        p_load = (trans(Xz) * t) / tt;
        q_mat = (trans(Y) * t) / tt;
        
        // Modified Gram-Schmidt
        v = p_load;
        for (int k = 0; k < i; k++) {
          v -= V.col(k) * arma::as_scalar(trans(V.col(k)) * v);
        }
        v /= arma::norm(v);
        V.col(i) = v;
        
        S -= v * (trans(v) * S);
        
        weights.row(i) = trans(r);
        scores.col(i) = t;
        Yloadings.row(i) = trans(q_mat);
      }
    }
    
    arma::mat TtT_inv = arma::inv_sympd(trans(scores) * scores);
    Xloadings = trans(trans(Xz) * scores * TtT_inv);
  }
  
  // ==========================================================================
  // Common post-processing
  // ==========================================================================
  if (algorithm == "simpls") {
    projection_matrix = trans(weights);
  } else {
    projection_matrix = trans(weights) * arma::solve(Xloadings * trans(weights), arma::eye(ncomp, ncomp));
  }
  
  arma::rowvec ymean = arma::mean(Y, 0);
  double y_hat_mean;
  int idx = 0;
  
  for (int k = 0; k < ny; k++) {
    arma::vec jth_loading = Yloadings.col(k);
    for (int j = 0; j < ncomp; j++) {
      coefficients.col(idx) = projection_matrix.cols(0, j) * jth_loading.rows(0, j);
      y_hat_mean = arma::as_scalar(x_center_vec * coefficients.col(idx));
      bo(k, j) = ymean(k) - y_hat_mean;
      idx++;
    }
  }
  
  return Rcpp::List::create(
    Rcpp::Named("ncomp") = ncomp,
    Rcpp::Named("coefficients") = coefficients,
    Rcpp::Named("bo") = bo,
    Rcpp::Named("X_loadings") = Xloadings,
    Rcpp::Named("Y_loadings") = Yloadings,
    Rcpp::Named("projection_mat") = projection_matrix,
    Rcpp::Named("transf") = Rcpp::List::create(
      Rcpp::Named("Xcenter") = x_center_vec,
      Rcpp::Named("Xscale") = x_scale_vec
    ),
    Rcpp::Named("weights") = weights
  );
}

// List opls_get_basics(
//     arma::mat X, 
//     arma::mat Y, 
//     int ncomp,
//     bool scale,            
//     double maxiter,
//     double tol, 
//     String algorithm = "pls", 
//     const int xls_min_w = 3, 
//     const int xls_max_w = 15
// ) {
//   int ny = Y.n_cols;
//   int nynf = ncomp * Y.n_cols;
//   
//   arma::mat weights = arma::zeros(ncomp, X.n_cols);
//   arma::mat scores = arma::zeros(X.n_rows, ncomp);
//   arma::mat Xloadings = arma::zeros(ncomp, X.n_cols);
//   arma::mat Yloadings = arma::zeros(ncomp, ny);
//   arma::mat coefficients = arma::zeros(X.n_cols, nynf);
//   arma::mat bo = arma::zeros(ny, ncomp);
//   arma::mat Xscale;
//   arma::mat x_scale_vec;
//   arma::mat x_center_vec;
//   arma::mat Xz = X;
//   
//   if (scale) {
//     Xscale = arma::repmat(Rcpp::as<arma::mat>(get_column_sds(Xz)), Xz.n_rows, 1);
//     Xz = Xz / Xscale;
//     x_scale_vec =  Xscale.row(0);
//   }
//   x_center_vec = Rcpp::as<arma::mat>(get_column_means(Xz));
//   Xz = Xz - arma::repmat(x_center_vec, Xz.n_rows, 1);
//   
//   // matrices to declare
//   arma::mat Xpls = Xz;
//   arma::mat Ypls = Y;
//   arma::mat iypls;
//   arma::mat Yplsb;
//   arma::vec imsd;
//   int j;
//   bool keepg;
//   arma::mat previous_ts = arma::zeros(Xz.n_rows, 1);
//   arma::mat lb;
//   arma::mat cr;
//   arma::mat ts;
//   arma::mat w;
//   arma::mat p;
//   arma::mat q;
//   arma::mat cx;
//   arma::mat cy;
//   arma::mat projection_matrix;
//   
//   // ==========================================================================
//   // NIPALS algorithms (pls, mpls, xls)
//   // ==========================================================================
//   if (algorithm == "pls" || algorithm == "mpls" || algorithm == "xls") {
//     
//     for (int i = 0; i < ncomp; i++) {
//       Yplsb = Ypls;
//       // Select the Y variable with the largest standard deviation
//       imsd = get_col_largest_sd(Ypls);
//       iypls = Ypls.col(imsd[0]);
//       previous_ts.fill(0);
//       
//       j = 0;
//       keepg = true;
//       
//       while (keepg) {
//         if(j > 0) {
//           previous_ts = ts;
//         }
//         //Step 1: Compute a vector of loading weights...
//         // // 1.1 Compute the 'scaling factor'
//         // cr = sqrt(trans(iypls) * Xpls * trans(Xpls) * iypls);
//         // // 1.2 The weights are computed as the cross product of
//         // // X0 and Y0 divided by the 'scaling factor'...
//         // w = (trans(Xpls) * iypls) / repmat(cr, Xpls.n_cols, 1);
//         w = get_weights(Xpls, iypls, algorithm, xls_min_w, xls_max_w);
//         // Step 2: Compute the scores...
//         ts = Xpls * w;
//         // Step 3: Compute the X-loadings (p) and the Y-loadings (q)...
//         p = (trans(Xpls) * ts) / repmat((trans(ts) * ts), Xpls.n_cols, 1);
//         q = (trans(Yplsb) * ts) / repmat((trans(ts) * ts), Yplsb.n_cols, 1);
//         iypls = (Yplsb * q) / repmat((trans(q) * q), Xpls.n_rows, 1) ;
//         lb = abs(sum((ts - previous_ts) / ts));
//         keepg = lb[0] > tol;
//         j = j + 1;
//         if (maxiter <= j) {
//           keepg = false;
//         }
//       }
//       
//       // Step 4: The residual matrix
//       // of X is finally computed and...
//       cx = ts * trans(p) ;
//       Xpls = Xpls - cx;
//       // ... the vector of residuals of Y is also computed
//       cy = ts * trans(q);
//       Ypls = Ypls - cy;
//       // save the matrices corresponding to the loadings
//       // and scores..
//       weights.row(i) = trans(w);
//       scores.col(i) = ts;
//       Xloadings.row(i) = trans(p);
//       Yloadings.row(i) = trans(q);
//     }
//   }
//   // ==========================================================================
//   // SIMPLS algorithm
//   // ==========================================================================
//   if (algorithm == "simpls") {
//     // Cross-product matrix
//     arma::mat S = trans(Xz) * Y;
//     
//     // Orthonormal basis for deflation
//     arma::mat V(Xz.n_cols, ncomp, arma::fill::zeros);
//     
//     arma::vec r, t, p_load, v;
//     arma::mat q_mat;
//     double tt;
//     
//     for (int i = 0; i < ncomp; i++) {
//       // Weight vector: dominant left singular vector of S
//       if (ny == 1) {
//         r = S.col(0);
//       } else {
//         arma::mat U;
//         arma::vec s;
//         arma::mat Vt;
//         arma::svd_econ(U, s, Vt, S, "left");
//         r = U.col(0);
//       }
//       r = r / arma::norm(r);
//       
//       // Scores
//       t = Xz * r;
//       tt = arma::as_scalar(trans(t) * t);
//       
//       // X-loadings (for deflation orthogonalisation)
//       p_load = (trans(Xz) * t) / tt;
//       
//       // Y-loadings
//       q_mat = (trans(Y) * t) / tt;
//       
//       // Modified Gram-Schmidt orthogonalisation
//       v = p_load;
//       for (int k = 0; k < i; k++) {
//         v = v - V.col(k) * arma::as_scalar(trans(V.col(k)) * v);
//       }
//       v = v / arma::norm(v);
//       V.col(i) = v;
//       
//       // Deflate cross-product matrix
//       S = S - v * (trans(v) * S);
//       
//       // Store results
//       weights.row(i) = trans(r);
//       scores.col(i) = t;
//       Yloadings.row(i) = trans(q_mat);
//     }
//     
//     // Recompute X-loadings for proper reconstruction
//     arma::mat TtT_inv = arma::inv(trans(scores) * scores);
//     Xloadings = trans(trans(Xz) * scores * TtT_inv);
//   }
//   
//   // ==========================================================================
//   // Common post-processing (shared by all algorithms)
//   // ==========================================================================
//   projection_matrix = trans(weights) * arma::solve(Xloadings * trans(weights), arma::eye(Xloadings.n_rows, Xloadings.n_rows));
//   
//   arma::mat ymean = arma::mean(Y);
//   arma::vec ymean_vec = arma::vectorise(ymean);
//   arma::mat y_hat_mean;
//   arma::mat y_hat_mean_vec;
//   int idx = 0;
//   for (int k = 0; k < ny; k++) {
//     arma::mat jth_loading = Yloadings.col(k);
//     for (int j = 0; j < ncomp; j++) {
//       //FIXME: For every Y store the coefficients independently
//       coefficients.col(idx) = projection_matrix.cols(0, j) * jth_loading.rows(0, j);
//       y_hat_mean = x_center_vec * coefficients.col(idx);
//       y_hat_mean_vec = arma::vectorise(y_hat_mean);
//       bo(k,j) = ymean_vec(k) - y_hat_mean_vec(0);
//       idx = idx + 1;
//     }
//   }
//   //FIXME: For every Y store the coefficients independently
//   return Rcpp::List::create(
//     Rcpp::Named("ncomp") = ncomp,
//     Rcpp::Named("coefficients") = coefficients,
//     Rcpp::Named("bo") = bo,
//     Rcpp::Named("X_loadings") = Xloadings,
//     Rcpp::Named("Y_loadings") = Yloadings,
//     Rcpp::Named("projection_mat") = projection_matrix,
//     Rcpp::Named("transf") = Rcpp::List::create(
//       Rcpp::Named("Xcenter") = x_center_vec,
//       Rcpp::Named("Xscale") = x_scale_vec
//     ),
//     _["weights"] = weights
//   );
// }

//' @title Prediction function for the \code{opls} and \code{fopls} functions
//' @description Predicts response values based on a model generated by either by \code{opls} or the \code{fopls} functions. 
//' For internal use only!. 
//' @usage predict_opls(bo, b, ncomp, newdata, scale, Xscale)
//' @param bo a numeric value indicating the intercept.
//' @param b the matrix of regression coefficients.
//' @param ncomp an integer value indicating how may components must be used in the prediction.
//' @param newdata a matrix containing the predictor variables.
//' @param scale a logical indicating whether the matrix of predictors used to create the regression model was scaled.
//' @param Xscale if \code{scale = TRUE} a matrix of one row with the values that must be used for scaling \code{newdata}.
//' @return a matrix of predicted values.
//' @author Leonardo Ramirez-Lopez
//' @keywords internal 
//' @noRd
//' @useDynLib resemble
// [[Rcpp::export]]
Rcpp::NumericMatrix predict_opls(
    arma::mat bo, 
    arma::mat b, 
    int ncomp, 
    arma::mat newdata,
    bool scale,
    arma::mat Xscale
) {
  
  if (scale) {
    // newdata = newdata / arma::repmat(Xscale, newdata.n_rows, 1);
    newdata = newdata.each_row() / Xscale;
  } 
  
  // Not Necessary to center since b0 is used
  // Xz = Xz - arma::repmat(Xcenter, newdata.n_rows, 1);
  
  arma::mat predicted = (newdata * b.cols(0, ncomp - 1)) + arma::repmat(bo.cols(0, ncomp - 1), newdata.n_rows, 1);
  return Rcpp::wrap(predicted);
}

//' @title Projection function for the \code{opls} function
//' @description Projects new spectra onto a PLS space based on a model generated by either by \code{opls} or the \code{opls2} functions. 
//' For internal use only!. 
//' @usage project_opls(projection_mat, ncomp, newdata, scale, Xcenter, Xscale)
//' @param projection_mat the projection matrix generated by the \code{opls} function.
//' @param ncomp an integer value indicating how may components must be used in the prediction.
//' @param newdata a matrix containing the predictor variables.
//' @param scale a logical indicating whether the matrix of predictors used to create the regression model was scaled.
//' @param Xscale if \code{scale = TRUE} a matrix of one row with the values that must be used for scaling \code{newdata}.
//' @param Xcenter a matrix of one row with the values that must be used for centering \code{newdata}.
//' @return a matrix corresponding to the new spectra projected onto the PLS space 
//' @author Leonardo Ramirez-Lopez
//' @keywords internal 
//' @noRd
//' @useDynLib resemble
// [[Rcpp::export]]
Rcpp::NumericMatrix project_opls(
    arma::mat projection_mat, 
    int ncomp, 
    arma::mat newdata,
    bool scale,
    arma::mat Xcenter,
    arma::mat Xscale
){
  
  if(scale){
    // newdata = newdata / arma::repmat(Xscale, newdata.n_rows, 1);
    newdata = newdata.each_row() / Xscale;
  }
  
  //Necessary to center
  newdata = newdata - arma::repmat(Xcenter, newdata.n_rows, 1);
  
  arma::mat proj = newdata * projection_mat.cols(0, ncomp - 1);
  
  return Rcpp::wrap(proj);
}

//' @title Projection to pls and then re-construction
//' @description Projects spectra onto a PLS space and then reconstructs it back.
//' @usage reconstruction_error(x, 
//'                             projection_mat, 
//'                             xloadings, 
//'                             scale, 
//'                             Xcenter, 
//'                             Xscale, 
//'                             scale_back = FALSE)
//' @param x a matrix to project.
//' @param projection_mat the projection matrix generated by the \code{opls_get_basics} function.
//' @param xloadings the loadings matrix generated by the \code{opls_get_basics} function.
//' @param scale logical indicating if scaling is required
//' @param Xcenter a matrix of one row with the centering values
//' @param Xscale a matrix of one row with the scaling values
//' @param scale_back compute the reconstruction error after de-centering the 
//' data and de-scaling it.
//' @return a matrix of 1 row and 1 column.
//' @author Leonardo Ramirez-Lopez
//' @keywords internal 
//' @noRd
//' @useDynLib resemble
// [[Rcpp::export]]
Rcpp::NumericMatrix reconstruction_error(
    arma::mat x, 
    arma::mat projection_mat, 
    arma::mat xloadings,
    bool scale,
    arma::mat Xcenter,
    arma::mat Xscale, 
    bool scale_back = false
) {
  
  if (scale){
    // x = x / arma::repmat(Xscale, x.n_rows, 1);
    x = x.each_row() / Xscale;
  }
  
  // Necessary to center
  x = x.each_row() - Xcenter;
  // x = x - arma::repmat(Xcenter, x.n_rows, 1);
  
  arma::mat xrec = x;
  arma::mat xrmse;
  xrec = x * projection_mat * xloadings; 
  
  if (scale_back) {
    x = x.each_row() + Xcenter;
    xrec = xrec.each_row() + Xcenter;
    if (scale) {
      x = x.each_row() % Xscale;
      xrec = xrec.each_row() % Xscale;
    }
  }
  // if(scale){
  //   xrec = xrec % arma::repmat(Xscale, x.n_rows, 1);
  // }
  // 
  // //Necessary to center
  // xrec = xrec + arma::repmat(Xcenter, newdata.n_rows, 1);
  
  xrmse = arma::mean(sqrt(arma::mean(pow(x - xrec, 2), 0)), 1);
  return Rcpp::wrap(xrmse);
}

//' @title Internal Cpp function for performing leave-group-out cross-validations for pls regression 
//' @description For internal use only!. 
//' @usage opls_cv_cpp(X, Y, scale, method, 
//'                   mindices, pindices, 
//'                   min_component, ncomp, 
//'                   new_x, 
//'                   maxiter, tol, 
//'                   wapls_grid, 
//'                   algorithm, 
//'                   statistics = TRUE)
//' @param X a matrix of predictor variables.
//' @param Y a matrix of a single response variable.
//' @param scale a logical indicating whether the matrix of predictors 
//' (\code{X}) must be scaled.
//' @param method the method used for regression. One of the following options: 
//' \code{'pls'} or \code{'wapls'} or \code{'completewapls1p'}.
//' @param mindices a matrix with \code{n} rows and \code{m} columns where 
//' \code{m} is equivalent to the number of resampling iterations. The elements 
//' of each column indicate the indices of the observations to be used for 
//' modeling at each iteration.
//' @param pindices a matrix with \code{k} rows and \code{m} columns where 
//' \code{m} is equivalent to the number of 
//' resampling iterations. The elements of each column indicate the indices of 
//' the observations to be used for predicting at each iteration.
//' @param min_component an integer indicating the number of minimum pls 
//' components (if the \code{method = 'pls'}).
//' @param ncomp an integer indicating the number of pls components.
//' @param new_x a matrix of one row corresponding to the observation to be 
//' predicted (if the \code{method = 'wapls'}).
//' @param maxiter maximum number of iterations.
//' @param tol limit for convergence of the algorithm in the nipals algorithm.
//' @param wapls_grid the grid on which the search for the best combination of 
//' minimum and maximum pls factors of \code{'wapls'} is based on in case 
//' \code{method = 'completewapls1p'}.
//' @param algorithm either pls (\code{'pls'}) or modified pls (\code{'mpls'}). 
//' See \code{get_weigths} function.
//' @param statistics a logical value indicating whether the precision and 
//' accuracy statistics are to be returned, otherwise the predictions for each 
//' validation segment are retrieved.
//' @return 
//' if \code{statistics = true} a list containing the following one-row matrices:
//' \itemize{
//' \item{\code{rmse_seg}: the RMSEs.}
//' \item{\code{st_rmse_seg}: the standardized RMSEs.}
//' \item{\code{rsq_seg}: the coefficients of determination.}
//' } 
//' 
//' if \code{statistics = false} a list containing the following one-row matrices:
//' \itemize{
//' \item{\code{predictions}: the predictions of each of the validation 
//' segments in \code{pindices}. Each column in \code{pindices} contains the 
//' validation indices of a segment.}
//' \item{\code{st_rmse_seg}: the standardized RMSEs.}
//' \item{\code{rsq_seg}: the coefficients of determination.}
//' } 
//' 
//' If \code{method = "wapls"}, data of the pls weights are output in this 
//' list(\code{compweights}).
//'
//' If \code{method = "completewapls1"}, data of all the combination of 
//' components passed in \code{wapls_grid} are 
//' output in this list(\code{complete_compweights}).
//' 
//' @author Leonardo Ramirez-Lopez
//' @keywords internal 
//' @noRd
//' @useDynLib resemble
// [[Rcpp::export]]
List opls_cv_cpp(
    arma::mat X, 
    arma::mat Y, 
    bool scale,
    String method,
    arma::mat mindices,
    arma::mat pindices,
    int min_component,
    int ncomp,
    arma::mat new_x,
    double maxiter, 
    double tol,
    arma::mat wapls_grid, 
    String algorithm,
    bool statistics = true
){
  // Validate method argument
  if (method != "pls" && method != "wapls" && method != "completewapls1") {
    Rcpp::stop("'method' must be 'pls', 'wapls', or 'completewapls1'");
  }
  
  arma::mat rmseseg;
  arma::mat strmseseg;
  arma::mat rsqseg;
  
  arma::mat compweights;
  arma::mat crcompweights;
  
  arma::mat nypred;
  arma::mat ypred;
  arma::mat predictions;
  int pred_rows;
  if (!statistics) {
    // integer multiplication uses "*" while double uses "%"
    pred_rows = mindices.n_cols * pindices.n_rows;
    predictions = arma::zeros(pred_rows, ncomp);
  }
  int preds_counter = 0;
  
  if (method == "pls") {
    rmseseg = arma::zeros(ncomp, mindices.n_cols);
    strmseseg = arma::zeros(ncomp, mindices.n_cols);
    rsqseg = arma::zeros(ncomp, mindices.n_cols);
    
    List transf;
    
    for(int i = 0; (unsigned)i < mindices.n_cols; i++){
      
      // The subset for fitting the model
      arma::vec irows = mindices.col(i);
      arma::mat xmatslice = arma::zeros(mindices.n_rows, X.n_cols);
      arma::mat ymatslice = arma::zeros(mindices.n_rows, Y.n_cols);
      
      
      for (int j = 0; (unsigned)j < irows.size(); j++) {
        xmatslice.row(j) = X.row(irows(j)-1);
        ymatslice.row(j) = Y.row(irows(j)-1);
      }
      
      
      // The subset for predicting with the model
      arma::vec pirows = pindices.col(i);
      arma::mat pxmatslice = arma::zeros(pindices.n_rows, X.n_cols);
      arma::mat pymatslice = arma::zeros(pindices.n_rows, Y.n_cols);
      
      for (int j = 0; (unsigned)j < pirows.size(); j++) {
        pxmatslice.row(j) = X.row(pirows(j)-1);
        pymatslice.row(j) = Y.row(pirows(j)-1);
      }
      
      arma::mat rpymatslice;
      rpymatslice = arma::repmat(pymatslice, 1, ncomp);
      
      List fit = Rcpp::as<Rcpp::List>(opls_get_basics(xmatslice, ymatslice, ncomp, scale, maxiter, tol, algorithm));
      
      transf = fit["transf"];   
      
      ypred = Rcpp::as<arma::mat>(predict_opls(
        fit["bo"], 
           fit["coefficients"],
              ncomp, 
              pxmatslice,
              scale,
              transf["Xscale"]
      ));
      
      if (!statistics){
        // predictions.row(i) = ypred;
        predictions.rows(preds_counter, preds_counter + ypred.n_rows - 1) = ypred;
        preds_counter = preds_counter +  ypred.n_rows;
      } else {
        arma::mat rdl = sqrt(get_column_means(pow(rpymatslice - ypred, 2)));
        rmseseg.col(i) = rdl;
        arma::mat mimav = arma::zeros(1,1);
        mimav.col(0) = max(pymatslice) - min(pymatslice);
        strmseseg.col(i) = rmseseg.col(i) / arma::repmat(mimav, ncomp, 1);
        rsqseg.col(i) = pow(arma::cor(ypred, pymatslice), 2);
      }
    }
  }
  
  if (method == "wapls"){
    
    rmseseg = arma::zeros(1, mindices.n_cols);
    strmseseg = arma::zeros(1, mindices.n_cols);
    rsqseg = arma::zeros(1, mindices.n_cols);
    
    // define the wapls weights directly here
    List cfit = Rcpp::as<Rcpp::List>(opls(X, Y, ncomp, scale, maxiter, tol, algorithm));
    List ctransf = cfit["transf"];
    
    compweights = arma::zeros(1, ncomp);
    compweights.cols(min_component-1, ncomp-1) =  Rcpp::as<arma::mat>(
      get_local_pls_weights(
        cfit["projection_mat"], 
            cfit["X_loadings"],
                cfit["coefficients"],
                    new_x,
                    min_component, 
                    ncomp, 
                    scale,
                    ctransf["Xcenter"],
                           ctransf["Xscale"]
      )
    );
    
    arma::mat rcompweights = arma::repmat(compweights, pindices.n_rows, 1);
    
    List transf;
    
    for(int i = 0; (unsigned)i < mindices.n_cols; i++){
      
      // The subset for fitting the model
      arma::vec irows = mindices.col(i);
      arma::mat xmatslice = arma::zeros(mindices.n_rows, X.n_cols);
      arma::mat ymatslice = arma::zeros(mindices.n_rows, Y.n_cols);
      
      
      for (int j = 0; (unsigned)j < irows.size(); j++) {
        xmatslice.row(j) = X.row(irows(j)-1);
        ymatslice.row(j) = Y.row(irows(j)-1);
      }
      
      
      // The subset for predicting with the model
      arma::vec pirows = pindices.col(i);
      arma::mat pxmatslice = arma::zeros(pindices.n_rows, X.n_cols);
      arma::mat pymatslice = arma::zeros(pindices.n_rows, Y.n_cols);
      
      for (int j = 0; (unsigned)j < pirows.size(); j++) {
        pxmatslice.row(j) = X.row(pirows(j)-1);
        pymatslice.row(j) = Y.row(pirows(j)-1);
      }
      
      
      List fit = Rcpp::as<Rcpp::List>(opls(xmatslice, ymatslice, ncomp, scale, maxiter, tol, algorithm));
      
      transf = fit["transf"];   
      
      nypred = Rcpp::as<arma::mat>(predict_opls(fit["bo"], 
                                                fit["coefficients"], 
                                                   ncomp, 
                                                   pxmatslice,
                                                   scale,
                                                   transf["Xscale"]));
      ypred = arma::sum(rcompweights % nypred, 1);
      
      if (!statistics){
        // predictions.row(i) = ypred;
        predictions.rows(preds_counter, preds_counter + ypred.n_rows - 1) = ypred;
        preds_counter = preds_counter +  ypred.n_rows;
      } else {
        arma::mat rdl = sqrt(get_column_means(pow(pymatslice - ypred, 2)));
        rmseseg.col(i) = rdl;
        arma::mat mimav = arma::zeros(1, 1);
        mimav.col(0) = max(pymatslice) - min(pymatslice);
        strmseseg.col(i) = rmseseg.col(i) / mimav;
        rsqseg.col(i) = pow(arma::cor(ypred, pymatslice), 2);
      }
    }
  }
  
  if (method == "completewapls1") {
    rmseseg = arma::zeros(wapls_grid.n_rows, mindices.n_cols);
    strmseseg = arma::zeros(wapls_grid.n_rows, mindices.n_cols);
    rsqseg = arma::zeros(wapls_grid.n_rows, mindices.n_cols);
    
    // define the wapls weights directly here
    List cfit = Rcpp::as<Rcpp::List>(opls_get_basics(X, Y, ncomp, scale, maxiter, tol, algorithm));
    List ctransf = cfit["transf"];
    
    compweights = arma::zeros(1, ncomp);
    compweights.cols(min_component-1, ncomp-1) =  Rcpp::as<arma::mat>(
      get_local_pls_weights(
        cfit["projection_mat"], 
            cfit["X_loadings"],
                cfit["coefficients"],
                    new_x,
                    min_component, 
                    ncomp, 
                    scale,
                    ctransf["Xcenter"],
                           ctransf["Xscale"]
      )
    );
    
    crcompweights = arma::zeros(wapls_grid.n_rows, ncomp); 
    for(int i = 0; (unsigned)i < crcompweights.n_rows; i++){
      int minpls = wapls_grid(i,0);
      int maxpls = wapls_grid(i,1);
      arma::mat subw = arma::zeros(1, ncomp); 
      subw.cols(minpls - 1, maxpls - 1) = compweights.cols(minpls - 1, maxpls - 1);
      arma::mat sumsubw = arma::repmat(sum(subw, 1), 1, ncomp);
      crcompweights.row(i) = subw / sumsubw;
    }
    
    arma::mat rcompweights = arma::repmat(compweights, pindices.n_rows, 1);
    
    List transf;
    
    for(int i = 0; (unsigned)i < mindices.n_cols; i++){
      // The subset for fitting the model
      arma::vec irows = mindices.col(i);
      arma::mat xmatslice = arma::zeros(mindices.n_rows, X.n_cols);
      arma::mat ymatslice = arma::zeros(mindices.n_rows, Y.n_cols);
      
      for (int j = 0; (unsigned)j < irows.size(); j++) {
        xmatslice.row(j) = X.row(irows(j)-1);
        ymatslice.row(j) = Y.row(irows(j)-1);
      }

      // The subset for predicting with the model
      arma::vec pirows = pindices.col(i);
      arma::mat pxmatslice = arma::zeros(pindices.n_rows, X.n_cols);
      arma::mat pymatslice = arma::zeros(pindices.n_rows, Y.n_cols);
      
      for (int j = 0; (unsigned)j < pirows.size(); j++) {
        pxmatslice.row(j) = X.row(pirows(j)-1);
        pymatslice.row(j) = Y.row(pirows(j)-1);
      }

      List fit = Rcpp::as<Rcpp::List>(opls_get_basics(xmatslice, ymatslice, ncomp, scale, maxiter, tol, algorithm));
      
      transf = fit["transf"];   
      
      arma::mat nypred;
      arma::mat ypred;
      arma::mat rpymatslice;
      rpymatslice = arma::repmat(pymatslice, 1, wapls_grid.n_rows);
      
      nypred = (Rcpp::as<arma::mat>(predict_opls(
        fit["bo"], 
           fit["coefficients"], 
              ncomp, 
              pxmatslice,
              scale,
              transf["Xscale"]
      )));
      
      ypred = nypred * trans(crcompweights);
      
      //ypred = arma::sum(rcompweights % nypred, 1);
      if (!statistics){
        // predictions.row(i) = ypred;
        predictions.rows(preds_counter, preds_counter + ypred.n_rows - 1) = ypred;
        preds_counter = preds_counter +  ypred.n_rows;
      } else {
        arma::mat rdl = sqrt(get_column_means(pow(rpymatslice - ypred, 2)));
        rmseseg.col(i) = rdl;
        arma::mat mimav = arma::zeros(1, 1);
        mimav.col(0) = max(pymatslice) - min(pymatslice);
        strmseseg.col(i) = rmseseg.col(i) / arma::repmat(mimav, wapls_grid.n_rows, 1);
        rsqseg.col(i) = pow(arma::cor(ypred, pymatslice), 2);
      }
    }
  }
  
  if (!statistics) { 
    return Rcpp::List::create(
      Rcpp::Named("predictions") = predictions,
      Rcpp::Named("compweights") = compweights,
      Rcpp::Named("complete_compweights") = crcompweights
    );
  } else {
    // here all the weights are output from 1 to ncomp (if method == wapls)
    // zeroes are assigned to those which are not selected at the begining
    return Rcpp::List::create(
      Rcpp::Named("rmse_seg") = rmseseg,
      Rcpp::Named("st_rmse_seg") = strmseseg,
      Rcpp::Named("rsq_seg") = rsqseg,
      Rcpp::Named("compweights") = compweights,
      Rcpp::Named("complete_compweights") = crcompweights
    );
  }
}

//' @title Orthogonal scores algorithm of partial least squares for gesearch
//' @description Computes orthogonal scores partial least squares (PLS) 
//' regression using either NIPALS or SIMPLS algorithm. This function is 
//' optimised for the \code{gesearch} evolutionary search and computes only 
//' the outputs required for weakness score evaluation: predictions, 
//' reconstruction error, and score-space dissimilarity.
//' 
//' NOTE: This function supports only a single response variable (PLS1). 
//' For internal use only.
//' @usage 
//' opls_gesearch(Xr, Yr, Xu, ncomp, scale,
//'               response = FALSE, reconstruction = TRUE,
//'               similarity = TRUE, fresponse = TRUE,
//'               algorithm = "pls")
//'         
//' @param Xr a matrix of predictor variables for the reference/training set.
//' @param Yr a single-column matrix of the response variable for the 
//' reference/training set. Only single-response (PLS1) is supported.
//' @param Xu a matrix of predictor variables for the target/test set.
//' @param ncomp the number of PLS components.
//' @param scale logical indicating whether \code{Xr} and \code{Xu} must be 
//' scaled. Centering is always applied using parameters derived from 
//' \code{Xr}.
//' @param response logical indicating whether to compute predictions for 
//' \code{Xu}. Used for the response weakness score (\code{w_r}) in 
//' \code{gesearch}. Default is \code{FALSE}.
//' @param reconstruction logical indicating whether to compute the 
//' reconstruction error of \code{Xu}. Used for the reconstruction weakness 
//' score (\code{w_q}) in \code{gesearch}. Default is \code{TRUE}.
//' @param similarity logical indicating whether to compute the distance 
//' between \code{Xr} and \code{Xu} in the PLS score space. Used for the 
//' similarity weakness score (\code{w_d}) in \code{gesearch}. Default is 
//' \code{TRUE}.
//' @param fresponse logical indicating whether to compute the proportion of 
//' response variance not explained by the model. Default is \code{TRUE}.
//' @param algorithm a character string indicating the PLS algorithm to use:
//' \itemize{
//'   \item{\code{'pls'}: standard PLS using covariance between X and Y for 
//'     weight computation (NIPALS algorithm).}
//'   \item{\code{'mpls'}: modified PLS using correlation between X and Y for 
//'     weight computation (NIPALS algorithm). See Shenk and Westerhaus (1991).}
//'   \item{\code{'simpls'}: SIMPLS algorithm (de Jong, 1993). Computationally 
//'     faster as it avoids iterative X deflation.}
//' }
//' @return a list containing:
//' \itemize{
//'   \item{\code{ncomp}: the number of components used.}
//'   \item{\code{pred_response}: predictions for \code{Xu} (only if 
//'     \code{response = TRUE}).}
//'   \item{\code{rmse_reconstruction}: RMSE of the spectral reconstruction 
//'     for \code{Xu} (only if \code{reconstruction = TRUE}).}
//'   \item{\code{score_dissimilarity}: mean Euclidean distance between 
//'     \code{Xr} and \code{Xu} scores in Mahalanobis-scaled PLS space 
//'     (only if \code{similarity = TRUE}).}
//'   \item{\code{residual_variance}: proportion of response variance not 
//'     explained by the model (only if \code{fresponse = TRUE}).}
//' }
//' @details
//' This function is designed for repeated evaluation within the 
//' \code{gesearch} evolutionary search algorithm, where it may be called 
//' ~10^5 times per run. It computes only the outputs necessary for 
//' calculating weakness scores:
//' \itemize{
//'   \item{Response weakness (\code{w_r}): prediction RMSE on the target, 
//'     requires \code{response = TRUE}.}
//'   \item{Reconstruction weakness (\code{w_q}): spectral reconstruction 
//'     error, requires \code{reconstruction = TRUE}.}
//'   \item{Similarity weakness (\code{w_d}): Mahalanobis distance in score 
//'     space, requires \code{similarity = TRUE}.}
//' }
//' 
//' Preprocessing applies scaling (if requested) followed by centering, 
//' using parameters derived from \code{Xr}. The same transformation is 
//' applied to \code{Xu}.
//' 
//' The \code{'simpls'} algorithm is faster than NIPALS-based methods 
//' (\code{'pls'}, \code{'mpls'}) as it avoids iterative X deflation. 
//' However, reconstruction errors and score-space distances differ 
//' numerically between algorithms (rankings are typically similar). 
//' Do not mix algorithms within a single \code{gesearch} run.
//' 
//' @references
//' de Jong, S. (1993). SIMPLS: An alternative approach to partial least 
//' squares regression. Chemometrics and Intelligent Laboratory Systems, 
//' 18(3), 251-263.
//' 
//' Shenk, J.S., & Westerhaus, M.O. (1991). Populations structuring of 
//' near infrared spectra and modified partial least squares regression. 
//' Crop Science, 31(6), 1548-1555.
//' 
//' @author Leonardo Ramirez-Lopez
//' @keywords internal 
//' @noRd
//' @useDynLib resemble
// [[Rcpp::export]]
List opls_gesearch(
    arma::mat Xr, 
    arma::mat Yr,
    arma::mat Xu, 
    int ncomp,
    bool scale,     
    bool response = false, 
    bool reconstruction = true,
    bool similarity = true,
    bool fresponse = true,
    String algorithm = "pls"
) {
  
  int ny = Yr.n_cols;
  int nynf = ncomp * ny;
  
  arma::mat weights(ncomp, Xr.n_cols, arma::fill::zeros);
  arma::mat scores(Xr.n_rows, ncomp, arma::fill::zeros);
  arma::mat Xrloadings(ncomp, Xr.n_cols, arma::fill::zeros);
  arma::vec Yrloadings(ncomp);
  arma::mat coefficients(Xr.n_cols, nynf, arma::fill::zeros);
  arma::mat bo(ny, ncomp, arma::fill::zeros);
  arma::rowvec Xscale;
  arma::rowvec Xcenter;
  arma::mat Xrz = Xr;
  arma::mat Xuz = Xu;
  
  if (scale) {
    Xscale = arma::stddev(Xrz, 0, 0);
    Xrz.each_row() /= Xscale;
    Xuz.each_row() /= Xscale;
  }
  Xcenter = arma::mean(Xrz, 0);
  Xrz.each_row() -= Xcenter;
  Xuz.each_row() -= Xcenter;
  
  arma::mat Xrpls = Xrz;
  arma::mat Yrpls = Yr;
  arma::mat ts, w, p, q;
  arma::mat projection_matrix;
  double tsts;
  
  // ==========================================================================
  // NIPALS (pls, mpls)
  // ==========================================================================
  if (algorithm == "pls" || algorithm == "mpls") {
    for (int i = 0; i < ncomp; i++) {
      w = get_weights(Xrpls, Yrpls.col(0), algorithm, 0, 0);
      ts = Xrpls * w;
      tsts = arma::as_scalar(trans(ts) * ts);
      p = (trans(Xrpls) * ts) / tsts;
      q = (trans(Yrpls) * ts) / tsts;
      
      Xrpls -= ts * trans(p);
      Yrpls -= ts * trans(q);
      
      weights.row(i) = trans(w);
      scores.col(i) = ts;
      Xrloadings.row(i) = trans(p);
      Yrloadings(i) = q(0, 0);
    }
  }
  
  // ==========================================================================
  // SIMPLS
  // ==========================================================================
  if (algorithm == "simpls") {
    arma::vec S = trans(Xrz) * Yr.col(0);
    arma::mat V(Xrz.n_cols, ncomp, arma::fill::zeros);
    arma::vec r, t, p_load, v;
    double tt, q_val, vtS;
    
    for (int i = 0; i < ncomp; i++) {
      r = S / arma::norm(S);
      t = Xrz * r;
      tt = arma::as_scalar(trans(t) * t);
      p_load = (trans(Xrz) * t) / tt;
      q_val = arma::as_scalar(trans(Yr.col(0)) * t) / tt;
      
      v = p_load;
      for (int k = 0; k < i; k++) {
        v -= V.col(k) * arma::as_scalar(trans(V.col(k)) * v);
      }
      v /= arma::norm(v);
      V.col(i) = v;
      
      vtS = arma::as_scalar(trans(v) * S);
      S -= v * vtS;
      
      weights.row(i) = trans(r);
      scores.col(i) = t;
      Yrloadings(i) = q_val;
    }
    
    // Compute X-loadings only if needed (stable version)
    if (reconstruction || similarity) {
      arma::mat TtT_inv = arma::inv_sympd(trans(scores) * scores);
      Xrloadings = trans(trans(Xrz) * scores * TtT_inv);
    }
  }
  
  // ==========================================================================
  // Projection matrix (only if needed)
  // ==========================================================================
  if (response || reconstruction || similarity) {
    if (algorithm == "simpls") {
      projection_matrix = trans(weights);
    } else {
      projection_matrix = trans(weights) * arma::solve(Xrloadings * trans(weights), arma::eye(ncomp, ncomp));
    }
  }
  
  // ==========================================================================
  // Response prediction
  // ==========================================================================
  arma::mat predicted;
  if (response) {
    double ymean = arma::as_scalar(arma::mean(Yr));
    coefficients = cumsum(projection_matrix.each_row() % Yrloadings.t(), 1);
    predicted = Xuz * coefficients.col(ncomp - 1) + ymean;
  }
  
  // ==========================================================================
  // Residual variance
  // ==========================================================================
  arma::mat residual_variance;
  if (fresponse) {
    residual_variance = pow(arma::cor(Yr, scores * Yrloadings), 2);
    residual_variance(0, 0) = 1.0 - residual_variance(0, 0);
  }
  
  // ==========================================================================
  // Reconstruction error (in standardised space)
  // ==========================================================================
  arma::mat xrmse;
  if (reconstruction) {
    // we don't need to add/subtract Xcenter because it 
    // cancels when computing Xu - xrec
    arma::mat err = Xuz - Xuz * projection_matrix * Xrloadings;
    if (scale) {
      err.each_row() %= Xscale;
    }
    xrmse = arma::mean(sqrt(arma::mean(pow(err, 2), 0)), 1);
  }
  
  // ==========================================================================
  // Similarity (score dissimilarity)
  // ==========================================================================
  arma::mat diss_score;
  if (similarity) {
    arma::mat scores_xu = Xuz * projection_matrix;
    
    if (ncomp > 1) {
      arma::mat scores_joined = join_cols(scores, scores_xu);
      arma::mat covsc = arma::cov(scores_joined);
      arma::mat U;
      arma::vec s;
      arma::mat Vsvd;
      svd_econ(U, s, Vsvd, covsc, "left");
      arma::mat sqrt_sm = arma::inv_sympd(U * diagmat(sqrt(s)) * U.t());
      scores_xu *= sqrt_sm;
      scores *= sqrt_sm;
    } else {
      arma::mat scores_joined = join_cols(scores, scores_xu);
      double scmean = arma::as_scalar(arma::mean(scores_joined));
      double scscale = arma::as_scalar(arma::stddev(scores_joined));
      scores_xu = (scores_xu - scmean) / scscale;
      scores = (scores - scmean) / scscale;
    }
    
    // Pairwise squared Euclidean distance
    arma::vec ss_ref = arma::sum(arma::square(scores), 1);
    arma::vec ss_xu = arma::sum(arma::square(scores_xu), 1);
    arma::mat diss = arma::repmat(ss_ref.t(), scores_xu.n_rows, 1) + 
      arma::repmat(ss_xu, 1, scores.n_rows) - 
      2.0 * scores_xu * scores.t();
    diss /= static_cast<double>(ncomp);
    diss_score = arma::mean(arma::mean(sqrt(diss), 0), 1);
  }
  
  return Rcpp::List::create(
    Rcpp::Named("ncomp") = ncomp,
    Rcpp::Named("pred_response") = predicted,
    Rcpp::Named("rmse_reconstruction") = xrmse,
    Rcpp::Named("score_dissimilarity") = diss_score,
    Rcpp::Named("residual_variance") = residual_variance
  );
}
// List opls_gesearch(
//     arma::mat Xr, 
//     arma::mat Yr,
//     arma::mat Xu, 
//     int ncomp,
//     bool scale,     
//     bool response = false, 
//     bool reconstruction = true,
//     bool similarity = true,
//     bool fresponse = true,
//     String algorithm = "pls"
// ) {
//   
//   int ny = Yr.n_cols;
//   int nynf = ncomp * Yr.n_cols;
//   
//   arma::mat weights = arma::zeros(ncomp, Xr.n_cols);
//   arma::mat scores = arma::zeros(Xr.n_rows, ncomp);
//   arma::mat Xrloadings = arma::zeros(ncomp, Xr.n_cols);
//   arma::vec Yrloadings (ncomp);
//   arma::mat coefficients = arma::zeros(Xr.n_cols, nynf);
//   arma::mat bo = arma::zeros(ny, ncomp);
//   arma::mat Xscale;
//   arma::mat Xcenter;
//   // arma::mat Xr_scale_vec;
//   // arma::mat Xr_center_vec;
//   arma::mat Xrz = Xr;
//   arma::mat Xuz = Xu;
//   if (scale) {
//     Xscale = arma::stddev(Xrz, 0, 0);
//     Xrz = Xrz.each_row() / Xscale;
//     Xuz = Xuz.each_row() / Xscale;
//   }
//   Xcenter = arma::mean(Xrz, 0);
//   Xrz = Xrz.each_row() - Xcenter;
//   Xuz = Xuz.each_row() - Xcenter;
//   
//   // matrices to declare
//   arma::mat Xrpls = Xrz;
//   arma::mat Yrpls = Yr;
//   arma::mat iypls;
//   arma::mat Yrplsb;
//   arma::vec lagest_sd_col;
//   arma::mat lb;
//   arma::mat cr;
//   arma::mat ts;
//   arma::mat w;
//   arma::mat p;
//   arma::mat q;
//   arma::mat cx;
//   arma::mat cy;
//   arma::mat projection_matrix;
//   
//   if (algorithm == "pls" || algorithm == "mpls" || algorithm == "xls") {
//     for (int i = 0; i < ncomp; i++) {
//       Yrplsb = Yrpls;
//       iypls = Yrpls.col(0);
//       w = get_weights(Xrpls, iypls, algorithm, 0, 0);
//       ts = Xrpls * w;
//       p = (trans(Xrpls) * ts) / repmat((trans(ts) * ts), Xrpls.n_cols, 1);
//       q = (trans(Yrplsb) * ts) / repmat((trans(ts) * ts), Yrplsb.n_cols, 1);
//       iypls = (Yrplsb * q) / repmat((trans(q) * q), Xrpls.n_rows, 1) ;
//       cx = ts * trans(p) ;
//       Xrpls = Xrpls - cx;
//       cy = ts * trans(q);
//       Yrpls = Yrpls - cy;
//       weights.row(i) = trans(w);
//       scores.col(i) = ts;
//       Xrloadings.row(i) = trans(p);
//       arma::vec qvec = arma::vectorise(q);
//       Yrloadings(i) = qvec(0);
//     }
//   }
//   
//   if (algorithm == "simpls") {
//     // Cross-product vector (single Y)
//     arma::vec S = trans(Xrz) * Yr.col(0);
//     
//     // Orthonormal basis for deflation
//     arma::mat V(Xrz.n_cols, ncomp, arma::fill::zeros);
//     
//     arma::vec r, t, p_load, v;
//     double tt, q_val;
//     
//     for (int i = 0; i < ncomp; i++) {
//       // Weight vector (just normalise S for single Y)
//       r = S / arma::norm(S);
//       
//       // Scores
//       t = Xrz * r;
//       tt = arma::as_scalar(trans(t) * t);
//       
//       // X-loadings (for deflation orthogonalisation)
//       p_load = (trans(Xrz) * t) / tt;
//       
//       // Y-loading (scalar for single Y)
//       q_val = arma::as_scalar(trans(Yr.col(0)) * t) / tt;
//       
//       // Modified Gram-Schmidt orthogonalisation
//       v = p_load;
//       for (int k = 0; k < i; k++) {
//         v = v - V.col(k) * arma::as_scalar(trans(V.col(k)) * v);
//       }
//       v = v / arma::norm(v);
//       V.col(i) = v;
//       
//       // Deflate cross-product vector
//       S = S - v * arma::as_scalar(trans(v) * S);
//       
//       // Store results
//       weights.row(i) = trans(r);
//       scores.col(i) = t;
//       Yrloadings(i) = q_val;
//     }
//     
//     // Recompute X-loadings for proper reconstruction
//     arma::mat TtT_inv = arma::inv(trans(scores) * scores);
//     Xrloadings = trans(trans(Xrz) * scores * TtT_inv);
//   }
//   
//   if (response | reconstruction | similarity) {
//     projection_matrix = trans(weights) * arma::solve(Xrloadings * trans(weights), arma::eye(Xrloadings.n_rows, Xrloadings.n_rows));
//   }
//   
//   // arma::mat yrmse;
//   arma::mat predicted;
//   if (response) {
//     arma::vec ymean;
//     ymean = arma::mean(Yr);
//     coefficients = cumsum(projection_matrix.each_row() % Yrloadings.t(), 1);
//     predicted = Xuz * coefficients.col(ncomp - 1);
//     predicted.for_each( [ymean](arma::mat::elem_type& val) { val += ymean(0); });
//     // yrmse = sqrt(arma::mean(pow(predicted - Yu, 2), 0));
//   }
//   
//   arma::mat residual_variance;
//   if (fresponse) {
//     // residual_variance = arma::mean(pow(Yr - (scores * Yrloadings), 2));
//     residual_variance = pow(arma::cor(Yr, (scores * Yrloadings)), 2);
//     residual_variance.col(0) = 1 - residual_variance.col(0);
//   }
// 
//   arma::mat xrmse;
//   if (reconstruction) {
//     // we don't need to add/subtract Xcenter because it 
//     // cancels when computing Xu - xrec
//     arma::mat err = Xuz - Xuz * projection_matrix * Xrloadings;
//     if (scale) {
//      err.each_row() %= Xscale;
//     }
//     xrmse = arma::mean(sqrt(arma::mean(pow(err, 2), 0)), 1);
//   }
//   
//   arma::mat diss_score;
//   arma::mat scores_xu;  
//   if (similarity) {
//     arma::mat covsc;
//     arma::mat U;
//     arma::vec s;
//     arma::mat V;
//     arma::mat sqrt_sm;
// 
//     scores_xu = Xuz * projection_matrix.cols(0, ncomp - 1);
//     if (ncomp > 1) {
//       // project onto a Mahalanobis space (since by adding the scores of Xuz
//       // the jonined matrix of the scores of Xrz and Xuz is not exactly orthogonal
//       // its the covariance matrix is not diagonal
//       covsc = arma::cov(join_cols( scores, scores_xu ));
//       svd_econ(U, s, V, covsc, "left");
//       
//       sqrt_sm = arma::solve(U * diagmat(sqrt(s.col(0))) * U.t(), arma::eye(ncomp, ncomp));
//       scores_xu = scores_xu * sqrt_sm;
//       scores = scores * sqrt_sm;
//     } else {
//       arma::vec scmean;
//       scmean = mean(join_cols( scores, scores_xu ));
//       scores.for_each( [scmean](arma::mat::elem_type& val) { val += scmean(0); });
//       scores_xu.for_each( [scmean](arma::mat::elem_type& val) { val += scmean(0); });
//       
//       arma::vec scscale;
//       scscale = Rcpp::as<arma::mat>(get_column_sds((join_cols( scores, scores_xu ))));
//       scores.for_each( [scscale](arma::mat::elem_type& val) { val /= scscale(0); });
//       scores_xu.for_each( [scscale](arma::mat::elem_type& val) { val /= scscale(0); });
//     }
//     // compute the distance matrix
//     arma::mat diss = arma::ones(scores_xu.n_rows, 1) * arma::sum(arma::square(scores), 1).t() + arma::sum(arma::square(scores_xu), 1)  * arma::ones(1, scores.n_rows) - 2 * scores_xu * scores.t();
//     diss.for_each( [ncomp](arma::mat::elem_type& val) { val /= ncomp; });
//     diss_score = arma::mean(arma::mean(sqrt(diss), 0), 1);
//   }
//   
//   return Rcpp::List::create(
//     Rcpp::Named("ncomp") = ncomp,
//     Rcpp::Named("pred_response") = predicted,
//     Rcpp::Named("rmse_reconstruction") = xrmse,
//     Rcpp::Named("score_dissimilarity") = diss_score,
//     Rcpp::Named("residual_variance") = residual_variance
//   );
// }



/// Gaussian process regression with linear kernel
//' @title Gaussian process regression with linear kernel (gaussian_process)
//' @description Carries out a gaussian process regression with a linear kernel (dot product). For internal use only!
//' @usage gaussian_process(X, Y, noisev, scale) 
//' @param X a matrix of predictor variables
//' @param Y a matrix with a single response variable
//' @param noisev a value indicating the variance of the noise for Gaussian process regression. Default is 0.001. a matrix with a single response variable
//' @param scale a logical indicating whether both the predictors 
//' and the response variable must be scaled to zero mean and unit variance.
//' @return a list containing the following elements:
//' \itemize{
//' \item{\code{b}: the regression coefficients.}
//' \item{\code{Xz}: the (final transformed) matrix of predictor variables.}
//' \item{\code{alpha}: the alpha matrix.}
//' \item{\code{is.scaled}: logical indicating whether both the predictors and response variable were scaled to zero mean and unit variance.}
//' \item{\code{Xcenter}: if matrix of predictors was scaled, the centering vector used for \code{X}.}
//' \item{\code{Xscale}: if matrix of predictors was scaled, the scaling vector used for \code{X}.}
//' \item{\code{Ycenter}: if matrix of predictors was scaled, the centering vector used for \code{Y}.}
//' \item{\code{Yscale}: if matrix of predictors was scaled, the scaling vector used for \code{Y}.}
//' }
//' @author Leonardo Ramirez-Lopez
//' @keywords internal 
//' @noRd
//' @useDynLib resemble
// [[Rcpp::export]]
List gaussian_process(arma::mat X, 
                      arma::mat Y, 
                      float noisev = 0.001,
                      bool scale = true
){
  
  // matrices to declare
  arma::mat K;
  arma::mat Xz = X;
  arma::mat Yz = Y;
  
  arma::mat Xscale;
  arma::mat Xcent;
  arma::mat xc;
  arma::mat xs;
  arma::mat Yscale;
  arma::mat Ycent;
  arma::mat yc;
  arma::mat ys;
  if(scale){
    xc = Rcpp::as<arma::mat>(get_column_means(Xz));
    Xcent = arma::repmat(xc, Xz.n_rows, 1);
    Xz = Xz - Xcent;
    xs = Rcpp::as<arma::mat>(get_column_sds(Xz));
    Xscale = arma::repmat(xs, Xz.n_rows, 1);
    Xz = Xz / Xscale;
    
    yc = arma::mean(Yz);
    Ycent = arma::repmat(yc, Yz.n_rows, 1);
    Yz = Yz - Ycent;
    ys = arma::stddev(Yz);
    Yscale = arma::repmat(ys, Yz.n_rows, 1);
    Yz = Yz / Yscale;
  }
  
  K = Xz * trans(Xz);
  
  arma::mat vrnc = arma::zeros(X.n_rows, X.n_rows);
  
  for(int i = 0; (unsigned)i < X.n_rows; i++){
    vrnc(i,i) = noisev;
  }
  
  arma::mat alpha = arma::solve(K + vrnc, arma::eye(X.n_rows, X.n_rows)) * Yz;
  
  arma::mat b = trans(Xz) * alpha;
  
  String method = "gaussian_process";
  
  return Rcpp::List::create(
    Rcpp::Named("b") = b,
    Rcpp::Named("Xz") = Xz,
    Rcpp::Named("alpha") = alpha,
    Rcpp::Named("is_scaled") = scale,
    Rcpp::Named("Xcenter") = xc,
    Rcpp::Named("Xscale") = xs,
    Rcpp::Named("Ycenter") = yc,
    Rcpp::Named("Yscale") = ys, 
    Rcpp::Named("regression") = method
  );
}


//' @title Prediction function for the \code{gaussian_process} function (Gaussian process regression with dot product covariance)
//' @description Predicts response values based on a model generated by the \code{gaussian_process} function (Gaussian process regression with dot product covariance). For internal use only!. 
//' @usage predict_gaussian_process(Xz, alpha, newdata, scale, Xcenter, Xscale, Ycenter, Yscale)
//' @param newdata a matrix containing the predictor variables
//' @param scale a logical indicating whether the matrix of predictors used to create the regression model 
//' (in the \code{gaussian_process} function) was scaled
//' @param Xcenter if \code{center = TRUE} a matrix of one row with the values that must be used for centering \code{newdata}.
//' @param Xscale if \code{scale = TRUE} a matrix of one row with the values that must be used for scaling \code{newdata}.
//' @param Ycenter if \code{center = TRUE} a matrix of one row with the values that must be used for accounting for the centering of the response variable.
//' @param Yscale if \code{scale = TRUE} a matrix of one row with the values that must be used  for accounting for the scaling of the response variable.
//' @return a matrix of predicted values
//' @author Leonardo Ramirez-Lopez
//' @keywords internal 
//' @noRd
//' @useDynLib resemble
// [[Rcpp::export]]
NumericVector predict_gaussian_process(arma::mat Xz, 
                                       arma::mat alpha, 
                                       arma::mat newdata,
                                       bool scale,
                                       arma::mat Xcenter,
                                       arma::mat Xscale,
                                       arma::mat Ycenter,
                                       arma::mat Yscale
){
  
  arma::mat newdatatr = newdata;
  
  if(scale){
    newdatatr = newdatatr.each_row() - Xcenter;
    newdatatr = newdatatr.each_row() / Xscale;
  }
  
  arma::mat predicted = newdatatr * trans(Xz) * alpha;
  
  if(scale){
    predicted = predicted % arma::repmat(Yscale, newdata.n_rows, 1) + arma::repmat(Ycenter, newdata.n_rows, 1);
  }
  
  return Rcpp::wrap(predicted);
}

//' @title Internal Cpp function for performing leave-group-out cross 
//' validations for gaussian process
//' @description For internal use only!. 
//' @usage gaussian_process_cv(X, Y, mindices, pindices, noisev = 0.001,  
//' scale = TRUE, statistics = TRUE)
//' @param X a matrix of predictor variables.
//' @param Y a matrix of a single response variable.
//' @param mindices a matrix with \code{n} rows and \code{m} columns where \code{m} is equivalent to the number of 
//' resampling iterations. The elements of each column indicate the indices of the observations to be used for modeling at each 
//' iteration.
//' @param pindices a matrix with \code{k} rows and \code{m} columns where \code{m} is equivalent to the number of 
//' resampling iterations. The elements of each column indicate the indices of the observations to be used for predicting at each 
//' iteration.
//' @param scale a logical indicating whether both the predictors 
//' and the response variable must be scaled to zero mean and unit variance.
//' @param statistics a logical value indicating whether the precision and 
//' accuracy statistics are to be returned, otherwise the predictions for each 
//' validation segment are retrieved.
//' @return a list containing the following one-row matrices:
//' \itemize{
//' \item{\code{rmse.seg}: the RMSEs.}
//' \item{\code{st.rmse.seg}: the standardized RMSEs.}
//' \item{\code{rsq.seg}: the coefficients of determination.}
//' } 
//' @author Leonardo Ramirez-Lopez
//' @keywords internal 
//' @noRd
//' @useDynLib resemble
// [[Rcpp::export]]
List gaussian_process_cv(arma::mat X, 
                         arma::mat Y, 
                         arma::mat mindices,
                         arma::mat pindices,
                         float noisev = 0.001,
                         bool scale = true, 
                         bool statistics = true
){
  
  arma::mat rmseseg = arma::zeros(1, mindices.n_cols);
  arma::mat strmseseg = arma::zeros(1, mindices.n_cols);
  arma::mat rsqseg = arma::zeros(1, mindices.n_cols);
  
  arma::mat predictions;
  int pred_rows;
  if (!statistics) {
    // integer multiplication uses "*" while double uses "%"
    pred_rows = mindices.n_cols * pindices.n_rows;
    predictions = arma::zeros(pred_rows, 1);
  }
  int preds_counter = 0;
  
  List transf;
  
  for(int i = 0; (unsigned)i < mindices.n_cols; i++){
    
    // The subset for fitting the model
    arma::vec irows = mindices.col(i);
    arma::mat xmatslice = arma::zeros(mindices.n_rows, X.n_cols);
    arma::mat ymatslice = arma::zeros(mindices.n_rows, Y.n_cols);
    
    
    for (int j = 0; (unsigned)j < irows.size(); j++) {
      xmatslice.row(j) = X.row(irows(j)-1);
      ymatslice.row(j) = Y.row(irows(j)-1);
    }
    
    
    // The subset for predicting with the model
    arma::vec pirows = pindices.col(i);
    arma::mat pxmatslice = arma::zeros(pindices.n_rows, X.n_cols);
    arma::mat pymatslice = arma::zeros(pindices.n_rows, Y.n_cols);
    
    for (int j = 0; (unsigned)j < pirows.size(); j++) {
      pxmatslice.row(j) = X.row(pirows(j)-1);
      pymatslice.row(j) = Y.row(pirows(j)-1);
    }
    
    List fit = Rcpp::as<Rcpp::List>(gaussian_process(xmatslice, ymatslice, noisev, scale));
    
    arma::mat ypred;
    
    ypred = Rcpp::as<arma::mat>(predict_gaussian_process(fit["Xz"], fit["alpha"], pxmatslice, scale, fit["Xcenter"], fit["Xscale"], fit["Ycenter"], fit["Yscale"]));
    
    if (!statistics){
      // predictions.row(i) = ypred;
      predictions.rows(preds_counter, preds_counter + ypred.n_rows - 1) = ypred;
      preds_counter = preds_counter +  ypred.n_rows;
    } else {
      arma::mat rdl = sqrt(get_column_means(pow(pymatslice - ypred, 2)));
      rmseseg.col(i) = rdl;
      arma::mat mimav = arma::zeros(1,1);
      mimav.col(0) = max(pymatslice) - min(pymatslice);
      strmseseg.col(i) = rmseseg.col(i) / mimav;
      rsqseg.col(i) = pow(arma::cor(ypred, pymatslice), 2);
    }
  }
  if (!statistics) { 
    return Rcpp::List::create(
      Rcpp::Named("predictions") = predictions
    );
  } else {
    return Rcpp::List::create(
      Rcpp::Named("rmse_seg") = rmseseg,
      Rcpp::Named("st_rmse_seg") = strmseseg,
      Rcpp::Named("rsq_seg") = rsqseg
    );
  }
}

//' @title Principal components based on  the non-linear iterative partial least squares (nipals) algorithm
//' @description Computes orthogonal socres partial least squares (opls) regressions with the NIPALS algorithm. It allows multiple response variables. 
//' For internal use only!
//' @usage 
//' pca_nipals(X, ncomp, center, scale,
//'            maxiter, tol,
//'            pcSelmethod = "var",
//'            pcSelvalue = 0.01)
//' @param X a matrix of predictor variables.
//' @param ncomp the number of pls components.
//' @param scale logical indicating whether \code{X} must be scaled.
//' @param maxiter maximum number of iterations.
//' @param tol limit for convergence of the algorithm in the nipals algorithm.
//' @param pcSelmethod the method for selecting the number of components. 
//' Options are: \code{'cumvar'} (for selecting the number of principal components based on a given 
//' cumulative amount of explained variance) and \code{"var"} (for selecting the number of principal 
//' components based on a given amount of explained variance). Default is \code{'var'}
//' @param pcSelvalue a numerical value that complements the selected method (\code{pcSelmethod}). 
//' If \code{"cumvar"} is chosen, it must be a value (larger than 0 and below 1) indicating the maximum 
//' amount of cumulative variance that the retained components should explain. If \code{"var"} is chosen, 
//' it must be a value (larger than 0 and below 1) indicating that components that explain (individually) 
//' a variance lower than this threshold must be excluded. If \code{"manual"} is chosen, it must be a value 
//' specifying the desired number of principal components to retain. Default is 0.01.
//' @return a list containing the following elements:
//' \itemize{
//' \item{\code{pc_scores}: a matrix of principal component scores.}
//' \item{\code{pc_loadings}: a matrix of of principal component loadings.}
//' \item{\code{variance}: a matrix of the variance of the principal components.} 
//' \item{\code{scale}: a \code{list} conating two objects: \code{center} and \code{scale}, which correspond to the vectors used to center and scale the input matrix.} 
//' } 
//' @author Leonardo Ramirez-Lopez
//' @keywords internal 
//' @noRd
//' @useDynLib resemble
// [[Rcpp::export]]
List pca_nipals(arma::mat X, 
                int ncomp,
                bool center,            
                bool scale,            
                double maxiter,
                double tol,
                String pcSelmethod = "var",
                double pcSelvalue = 0.01
){
  
  arma::mat Xscale;
  arma::mat Xcenter;
  arma::mat x_scale_vec;
  arma::mat x_center_vec;
  arma::mat Xz = X;
  
  if(center){
    Xcenter = Rcpp::as<arma::mat>(get_column_means(Xz));
    Xz = Xz - arma::repmat(Xcenter, Xz.n_rows, 1);
    x_center_vec = Xcenter.row(0);
  }
  
  if(scale){
    Xscale = arma::repmat(Rcpp::as<arma::mat>(get_column_sds(Xz)), Xz.n_rows, 1);
    Xz = Xz / Xscale;
    x_scale_vec =  Xscale.row(0);
  }
  
  arma::mat Xpls = Xz;
  //variance of Xpls
  double xvar = overall_var(Xpls)(0);  
  // int max_pcs = arma::min(arma::size(Xpls));
  
  // matrices to declare
  int iter;
  bool keepg;
  arma::mat pp;
  arma::mat cx;
  double ireconstructed_var;
  arma::mat pp_std;
  arma::mat tt_ith;
  arma::mat val_ith;
  arma::mat pc_scores = arma::zeros(Xpls.n_rows, ncomp);
  arma::mat pc_loadings = arma::zeros(Xpls.n_cols, ncomp);
  arma::mat explained_var = arma::zeros(3, ncomp);
  
  int ith_comp = 0;
  for (int i = 0; i < ncomp; i++){
    arma::mat tt = Xpls.col(0);
    iter = 0;
    keepg = true;
    while(keepg){
      pp = (trans(Xpls) * tt) /  repmat((trans(tt) * tt), Xpls.n_cols, 1);
      pp_std = pp / repmat(sqrt(trans(pp) * pp), pp.n_rows, 1);
      tt_ith = trans(trans(pp_std) * trans(Xpls));
      val_ith = sqrt(trans(tt - tt_ith) * (tt - tt_ith));
      keepg = val_ith[0] > tol;
      iter = iter + 1;
      if(maxiter <= iter){
        keepg = false;
      }
      tt.col(0) = tt_ith.col(0);
    }
    cx = (tt * trans(pp_std));
    Xpls = Xpls - cx;
    pc_scores.col(i) = tt.col(0);
    pc_loadings.col(i) = pp_std.col(0);
    
    ireconstructed_var = overall_var(cx)(0);
    explained_var(0,i) = ireconstructed_var;
    explained_var(1,i) = explained_var(0,i) / xvar;
    explained_var(2,i) = sum(explained_var.row(0)) / xvar;   
    
    ith_comp = ith_comp + 1;
    if(pcSelmethod == "var" || pcSelmethod == "cumvar")
    {
      bool chk;
      if(pcSelmethod == "cumvar"){
        chk = explained_var(2,i) > pcSelvalue;
      }
      else{
        chk = explained_var(1,i) < pcSelvalue;
      }
      if(chk)
      {
        ncomp = ith_comp - 1;
        ith_comp = ith_comp - 2;
        if (i == 0) {
          throw std::invalid_argument("With the current value in the 'pc_selection' argument, no components are selected. Try another value.");
        }
        break;
      }
    }
  }
  
  arma::uvec pc_indices;
  if(pcSelmethod == "var") 
  {
    pc_indices = find(explained_var.row(1) >= pcSelvalue); 
    pc_scores = pc_scores.cols(pc_indices);
    pc_loadings = pc_loadings.cols(pc_indices);
    explained_var = explained_var.cols(pc_indices);
  }
  
  if(pcSelmethod == "cumvar") 
  {
    pc_indices = find(explained_var.row(2) <= pcSelvalue && explained_var.row(2) > 0); 
    pc_indices = pc_indices + 1;
    pc_indices.insert_rows(0, 1);
    // pc_indices = find(explained_var.row(2) > 0); 
    pc_scores = pc_scores.cols(pc_indices);
    pc_loadings = pc_loadings.cols(pc_indices);
    explained_var = explained_var.cols(pc_indices);
  }
  
  return Rcpp::List::create(
    Rcpp::Named("pc_indices") = pc_indices,
    Rcpp::Named("pc_scores") = pc_scores,
    Rcpp::Named("pc_loadings") = trans(pc_loadings),
    Rcpp::Named("original_x_variance") = xvar, 
    Rcpp::Named("pc_variance") = explained_var,
    Rcpp::Named("scale") = Rcpp::List::create(
      Rcpp::Named("center") = x_center_vec,
      Rcpp::Named("scale") = x_scale_vec
    )
  );
}


//////////////////////////////////////////////////////
// ------ FUNCTIONS FOR THE LIBRARY OF MODELS ----- //
/////////////////////////////////////////////////////


//' @title Internal: Fit a local weighted PLS model and predict for a query point
//'
//' @description
//' Fits a local Partial Least Squares (PLS) model using a neighborhood subset
//' and computes a weighted prediction for a target sample. The weighting is
//' done over multiple components using a provided evaluation grid.
//'
//' @param X Numeric matrix of predictors from the local neighborhood
//'   (observations in rows, variables in columns).
//' @param Y Numeric matrix (single column) of corresponding response values.
//' @param xval Numeric matrix (single row) representing the query sample to
//'   predict.
//' @param emgrid Numeric matrix used to weight component-wise predictions.
//' @param ncomp_max Integer. Maximum number of PLS components to fit.
//' @param ncomp_min Integer. Minimum number of PLS components to use in
//'   prediction.
//' @param scale Logical. Whether to scale predictors before PLS fitting.
//' @param max_iter Numeric. Maximum number of iterations for the PLS algorithm.
//' @param tol Numeric. Convergence tolerance for the PLS algorithm.
//' @param algorithm Character. PLS algorithm to use: \code{"mpls"} (default),
//'   \code{"pls"} (nipals), or \code{"simpls"}.
//'
//' @return
//' A numeric vector of weighted predictions (length equal to number of rows
//' in \code{emgrid}).
//'
//' @details
//' The function performs the following steps:
//' \enumerate{
//'   \item Fits a PLS model on \code{X} and \code{Y} with up to
//'         \code{ncomp_max} components.
//'   \item Extracts centering and scaling parameters from the fitted model.
//'   \item Computes component weights for \code{xval} using
//'         \code{get_local_pls_weights()}.
//'   \item Predicts component-wise responses using \code{predict_opls()}.
//'   \item Applies \code{emgrid} weights scaled by component weights.
//'   \item Returns row-normalized weighted average of predictions.
//' }
//'
//' @author Leonardo Ramirez-Lopez
//' @keywords internal
//' @noRd
//' @useDynLib resemble
// [[Rcpp::export]]
Rcpp::NumericVector ith_local_fit(
    arma::mat X, arma::mat Y,
    arma::mat xval, arma::mat emgrid,
    int ncomp_max, int ncomp_min,
    bool scale, double max_iter, double tol, 
    String algorithm = "mpls"
) {
  // Step 1: Compute PLS model
  List ipls = opls_get_basics(X, Y, ncomp_max, scale, max_iter, tol, algorithm);
  
  // Step 2: Extract transformation info
  Rcpp::List transf = ipls["transf"];
  arma::rowvec Xcenter = as<arma::rowvec>(transf["Xcenter"]);
  arma::rowvec Xscale  = as<arma::rowvec>(transf["Xscale"]);
  
  // Step 3: Compute local weights
  arma::rowvec pcweights = as<arma::rowvec>(get_local_pls_weights(
    ipls["projection_mat"],
    ipls["X_loadings"],
    ipls["coefficients"],
    xval,
    ncomp_min,
    ncomp_max,
    scale,
    Xcenter,
    Xscale
  ));
  
  // Step 4: Predict using PLS model
  arma::mat preds = as<arma::mat>(predict_opls(
    ipls["bo"],
    ipls["coefficients"],
    ncomp_max,
    xval,
    scale,
    Xscale
  ));
  
  // Step 5: Multiply emgrid with pcweights (column-wise)
  arma::mat wmgrid = emgrid;
  // uivalent to sweep(..., MARGIN = 2)
  wmgrid.each_row() %= pcweights;  
  
  // Step 6: Normalize rows (equivalent to sweep(..., MARGIN = 1))
  arma::vec rowsums = arma::sum(wmgrid, 1);
  wmgrid.each_col() /= rowsums;
  
  // Step 7: Select and transpose relevant predictions
  arma::mat preds_slice = preds.cols(ncomp_min - 1, ncomp_max - 1).t();
  
  // Step 8: Matrix multiply
  arma::mat ipreds = wmgrid * preds_slice;
  Rcpp::NumericVector result = Rcpp::wrap(ipreds.col(0));
  return result;
}



//' @title Compute Final Local PLS Model Outputs
//' @description
//' This function calculates the final local model coefficients, intercept,
//' scaled variable importance (VIP), and selectivity ratio for a single observation
//' using a weighted combination of PLS components.
//'
//' @param X A matrix of predictor variables used for calibration.
//' @param Y A matrix of response variables used for calibration.
//' @param new_x A single observation (1 x p) of predictor variables to compute weights.
//' @param ncomp_min The minimum number of PLS components to include in the final model.
//' @param ncomp_max The maximum number of PLS components to include in the final model.
//' @param scale Logical indicating whether to scale the data.
//' @param maxiter Maximum number of iterations allowed during the NIPALS algorithm.
//' @param tol Tolerance threshold for convergence in the iterative algorithm.
//' @param algorithm \code{'mpls'} (defalt), \code{'pls'} (nipals), \code{'simpls'}.
//'
//' @return A list with the following elements:
//' \itemize{
//'   \item \code{ib0}: Intercept term, computed as a weighted sum of component-specific intercepts.
//'   \item \code{ibs}: Weighted regression coefficients.
//'   \item \code{ivips}: Weighted variable importance (VIP) scores, scaled by component SDs.
//'   \item \code{isratio}: Weighted selectivity ratios, scaled by component SDs.
//'   \item \code{Xscale}: Scaling vector used to scale \code{X} (if \code{scale = TRUE}).
//' }
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Rcpp::List final_fits_cpp(
    const arma::mat& X,
    const arma::mat& Y,
    const arma::mat& new_x,
    int ncomp_min,
    int ncomp_max,
    bool scale,
    double maxiter,
    double tol, 
    String algorithm = "mpls"
) {
  // 1) Fit iPLS (unchanged contract)
  Rcpp::List ipls = opls_get_all(
    X, Y, ncomp_max, scale, maxiter, tol, algorithm
  );
  
  // 2) Extract transform data safely
  Rcpp::List transf = ipls["transf"];
  
  arma::mat Xcenter_m = Rcpp::as<arma::mat>(transf["Xcenter"]);
  arma::rowvec Xcenter = arma::rowvec(Xcenter_m);     // expect 1 x p
  
  arma::mat Xscale_m = Rcpp::as<arma::mat>(transf["Xscale"]);
  arma::rowvec Xscale;
  if (Xscale_m.n_elem == 0) {
    // When scale == false, some builds return empty. Use ones(1 x p).
    Xscale = arma::rowvec(X.n_cols, arma::fill::ones);
  } else {
    Xscale = arma::rowvec(Xscale_m);                  // 1 x p
  }
  
  // 3) Local PLS weights
  //    get_local_pls_weights returns a 1 x k matrix -> make it (k x 1)
  arma::mat pcw_m = Rcpp::as<arma::mat>(get_local_pls_weights(
    ipls["projection_mat"],
        ipls["X_loadings"],
            ipls["coefficients"],
                new_x,
                ncomp_min,
                ncomp_max,
                scale,
                Xcenter,
                Xscale
  ));
  arma::vec pcweights = arma::conv_to<arma::vec>::from(pcw_m.t()); // (k x 1)
  
  // component indices in 0-based C++
  arma::uvec idx_comp = arma::regspace<arma::uvec>(
    static_cast<unsigned>(ncomp_min - 1),
    static_cast<unsigned>(ncomp_max - 1)
  );
  
  // 4) Final coefficients: (p x k) * (k x 1) -> (p x 1)
  arma::mat coefs = Rcpp::as<arma::mat>(ipls["coefficients"]);
  arma::mat ibs   = coefs.cols(idx_comp) * pcweights;
  
  // 5) Final intercept
  arma::mat  bo_m   = Rcpp::as<arma::mat>(ipls["bo"]); // (ny x ncomp)
  arma::rowvec b0_r = bo_m.row(0);                     // first response
  arma::vec intercepts = arma::conv_to<arma::vec>::from(
    b0_r.cols(idx_comp).t()
  );                                                   // (k x 1)
  double ib0 = arma::dot(intercepts, pcweights);
  
  // 6) Scaled VIPs
  arma::mat vip       = Rcpp::as<arma::mat>(ipls["vip"]); // (p x m)
  arma::mat vip_slice = vip.cols(idx_comp);               // (p x k)
  
  arma::vec vip_sds_c = Rcpp::as<arma::vec>(get_column_sds(vip_slice));
  vip_sds_c.transform(
    [](double v){ return (std::abs(v) < 1e-12) ? 1.0 : v; }
  );
  arma::rowvec vip_sds = vip_sds_c.t();                   // (1 x k)
  
  arma::mat svips = vip_slice;                            // (p x k)
  svips.each_row() /= vip_sds;                            // broadcast
  arma::mat ivips = svips * pcweights;                    // (p x 1)
  
  // 7) Scaled selectivity ratios
  arma::mat sr       = Rcpp::as<arma::mat>(ipls["selectivity_ratio"]);
  arma::mat sr_slice = sr.cols(idx_comp);                 // (p x k)
  
  arma::vec sr_sds_c = Rcpp::as<arma::vec>(get_column_sds(sr_slice));
  sr_sds_c.transform(
    [](double v){ return (std::abs(v) < 1e-12) ? 1.0 : v; }
  );
  arma::rowvec sr_sds = sr_sds_c.t();                     // (1 x k)
  
  arma::mat ssratio = sr_slice;                           // (p x k)
  ssratio.each_row() /= sr_sds;                           // broadcast
  arma::mat isratio = ssratio * pcweights;                // (p x 1)
  
  return Rcpp::List::create(
    Rcpp::Named("ib0") = ib0,
    Rcpp::Named("ibs") = ibs,
    Rcpp::Named("ivips") = ivips,
    Rcpp::Named("isratio") = isratio,
    Rcpp::Named("Xcenter") = Xcenter,
    Rcpp::Named("Xscale") = Xscale
  );
}


//' @title Internal: Predict using local PLS coefficients and optional
//' dissimilarities
//'
//' @description
//' Computes predictions for a new observation using local PLS models
//' represented by coefficients (\code{plslib}). The prediction is based on
//' inverse-scaled feature values. If a dissimilarity vector is provided, it is
//' prepended to the input features before inverse scaling.
//'
//' @param plslib A numeric matrix of PLS model coefficients (n_models × p+1).
//'   First column is the intercept; remaining columns are coefficients for 
//'   scaled features.
//' @param xscale A numeric matrix of scaling values (n_models × p), same 
//'   number of columns as \code{plslib[,-1]}.
//' @param Xu A numeric vector of length p representing the query sample to 
//'   be predicted.
//' @param dxrxu Optional numeric vector of dissimilarities between \code{Xu} 
//'   and reference samples. If provided, prepended to \code{Xu} as additional 
//'   predictive features. Default is \code{R_NilValue} (NULL).
//'
//' @return A numeric vector of length \code{nrow(plslib)} containing the 
//'   predicted response for \code{Xu} from each local model.
//'
//' @details
//' For each local model (row of \code{plslib}), the function:
//' \enumerate{
//'   \item Optionally prepends \code{dxrxu} to \code{Xu}
//'   \item Computes inverse-scaled features: \code{Xu / xscale[i,]}
//'   \item Computes prediction: \code{intercept + sum(coefficients * scaled_features)}
//' }
//'
//' @note
//' This is an internal function optimised for performance in the prediction
//' loop of \code{predict.liblex}.
//'
//' @keywords internal
//' @noRd
//' @author Leonardo Ramirez-Lopez
// [[Rcpp::export]]
NumericVector ith_pred_cpp(
   const NumericMatrix& plslib,
   const NumericMatrix& xscale,
   const NumericVector& Xu,
   Rcpp::Nullable<NumericVector> dxrxu = R_NilValue
) {
 int n_models = plslib.nrow();
 int p = xscale.ncol();
 NumericVector ipred(n_models);
 
 // Build input vector: optionally prepend dissimilarities
 NumericVector dxu;
 if (dxrxu.isNotNull()) {
   NumericVector dxrxu_vec(dxrxu);
   int nd = dxrxu_vec.size();
   dxu = NumericVector(nd + Xu.size());
   for (int i = 0; i < nd; i++) {
     dxu[i] = dxrxu_vec[i];
   }
   for (int i = 0; i < Xu.size(); i++) {
     dxu[nd + i] = Xu[i];
   }
 } else {
   dxu = Xu;
 }
 
 // Compute predictions for each model
 for (int i = 0; i < n_models; i++) {
   double pred = plslib(i, 0);  // intercept
   for (int j = 0; j < p; j++) {
     double scaled = dxu[j] / xscale(i, j);
     pred += plslib(i, j + 1) * scaled;
   }
   ipred[i] = pred;
 }
 
 return ipred;
}

