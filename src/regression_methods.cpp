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
//' @useDynLib resemble
// [[Rcpp::export]]
NumericVector get_column_sds(arma::mat X){
  arma::mat sds = arma::stddev(X, 0, 0);
  return Rcpp::wrap(sds);
}

//' @title Function for computing the mean of each column in a matrix
//' @description Computes the mean of each column in a matrix. For internal use only!
//' @usage get_column_means(X)
//' @param X a a matrix.
//' @return a vector of mean values. 
//' @author Leonardo Ramirez-Lopez
//' @keywords internal 
//' @useDynLib resemble
// [[Rcpp::export]]
NumericVector get_column_means(arma::mat X){
  arma::mat mn = arma::mean(X, 0);
  return Rcpp::wrap(mn);
}

//' @title Function for computing sum of each column in a matrix
//' @description Computes the sum of each column in a matrix. For internal use only!
//' @usage get_column_sums(X)
//' @param X a matrix.
//' @return a vector of standard deviation values. 
//' @author Leonardo Ramirez-Lopez
//' @keywords internal 
//' @useDynLib resemble
// [[Rcpp::export]]
NumericVector get_column_sums(arma::mat X){
  arma::mat sm = arma::sum(X);
  return Rcpp::wrap(sm);
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
//'                     pcSelvalue = 0.01)
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
//' @return a list containing the following elements:
//' \itemize{
//' \item{\code{coefficients}}{ the matrix of regression coefficients.}
//' \item{\code{bo}}{ a matrix of one row containing the intercepts for 
//' each component.}
//' \item{\code{scores}}{ the matrix of scores.}
//' \item{\code{X_loadings}}{ the matrix of X loadings.}
//' \item{\code{Y_loadings}}{ the matrix of Y loadings.}
//' \item{\code{projection_mat}}{ the projection matrix.}
//' \item{\code{Y}}{ the \code{Y} input.}
//' \item{\code{variance}}{ a \code{list} conating two objects: \code{x_var} 
//' and \code{y_var}. 
//' These objects contain information on the explained variance for the \code{X} 
//' and \code{Y} matrices respectively.}
//' \item{\code{transf}}{ a \code{list} conating two objects: \code{Xcenter} 
//' and \code{Xscale}}. 
//' \item{\code{weights}}{ the matrix of wheights.}
//' }
//' @author Leonardo Ramirez-Lopez
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
                         double pcSelvalue = 0.01
){
  
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
  
  if(scale){
    Xscale = arma::repmat(Rcpp::as<arma::mat>(get_column_sds(Xz)), Xz.n_rows, 1);
    Xz = Xz / Xscale;
    x_scale_vec =  Xscale.row(0);
  }
  x_center_vec = Rcpp::as<arma::mat>(get_column_means(Xz));
  Xz = Xz - arma::repmat(x_center_vec, Xz.n_rows, 1);
  
  arma::mat Xpls = Xz;
  arma::mat Ypls = Y;
  
  //variance of Xpls
  double xvar = sum(pow(get_column_sds(Xpls), 2));  
  
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
  arma::mat cx;
  arma::mat cy;
  arma::mat projection_matrix;
  
  arma::mat tsrp_lastith;
  
  int ith_comp = 0;
  for (int i = 0; i < ncomp; i++){
    Yplsb = Ypls;
    // Select the Y variable with the largest standard deviation
    lagest_sd_col = get_col_largest_sd(Ypls);
    iypls = Ypls.col(lagest_sd_col[0]);
    previous_ts.fill(0);
    
    j = 0;
    keepg = true;
    
    while (keepg) {
      if(j > 0) {
        previous_ts = ts;
      }
      //Step 1: Compute a vector of loading weights...
      // 1.1 Compute the 'scaling factor'
      cr = sqrt(trans(iypls) * Xpls * trans(Xpls) * iypls);
      // 1.2 The weights are computed as the cross product of
      // X0 and Y0 divided by the 'scaling factor'...
      w = (trans(Xpls) * iypls) / repmat(cr, Xpls.n_cols, 1);
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
    explained_var(0,i) = arma::var(scores.col(i));
    explained_var(1,i) = explained_var(0,i) / xvar;
    explained_var(2,i) = sum(explained_var.row(0)) / xvar;
    
    
    ith_comp = ith_comp + 1;
    
    if (pcSelmethod != "manual") {
      if (pcSelmethod == "var" || pcSelmethod == "cumvar") {
        bool chk;
        if (pcSelmethod == "cumvar") {
          chk = explained_var(2,i) > pcSelvalue;
        } else {
          chk = explained_var(1,i) < pcSelvalue;
        }
        if (chk) {
          ncomp = ith_comp - 1;
          ith_comp = ith_comp - 2;
          if(i == 0) {
            throw std::invalid_argument("With the current value in the 'pc_selection' argument, no components are selected. Try another value.");
          }
          if (pcSelmethod == "cumvar") {
            ncomp = ncomp + 1;
          }
          break;
        }
      } 
    }
  }
  
 
  
    
  arma::uvec pc_indices;
  if (pcSelmethod != "manual") {
    if (pcSelmethod == "var" || pcSelmethod == "cumvar") {
      if (pcSelmethod == "var") {
        pc_indices = find(explained_var.row(1) >= pcSelvalue); 
      } else {
        //pc_indices = find(explained_var.row(2) <= pcSelvalue && explained_var.row(2) > 0);
        pc_indices = find(explained_var.row(2) > 0); 
      }
      weights = weights.rows(pc_indices);
      coefficients = coefficients.cols(pc_indices);
      bo = bo.cols(pc_indices);
      scores = scores.cols(pc_indices);
      Xloadings = Xloadings.rows(pc_indices);
      Yloadings = Yloadings.rows(pc_indices);
      explained_var = explained_var.cols(pc_indices);
      yex = yex.cols(pc_indices);
    }
  }
  // convert this to standard deviation
  explained_var.row(0) = sqrt(explained_var.row(0));
  
  projection_matrix = trans(weights) * arma::solve(Xloadings * trans(weights), arma::eye(Xloadings.n_rows, Xloadings.n_rows));
  
  arma::mat yexi;
  arma::mat cop;
  for (int i = 0; i < ny; i++) {
    yexi = scores % arma::repmat(trans(Yloadings.col(i)), scores.n_rows, 1) ;
    cop = pow(arma::cor(Y.col(i), yexi.col(0)), 2);
    yex(i,0) = cop(0,0);
    for (int j = 1; j < ncomp; j++) {
      yexi.col(j) = yexi.col(j-1) + yexi.col(j);
      cop = arma::cor(Y.col(i), yexi.col(j));
      yex(i,j) = pow(cop(0,0), 2);
    }
  }
  
  
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
//'              tol)
//' @param X a matrix of predictor variables.
//' @param Y a matrix of either a single or multiple response variables.
//' @param ncomp the number of pls components.
//' @param scale logical indicating whether \code{X} must be scaled.
//' @param maxiter maximum number of iterations.
//' @param tol limit for convergence of the algorithm in the nipals algorithm.
//' @return a list containing the following elements:
//' \itemize{
//' \item{\code{ncomp}}{ the number of components used.}
//' \item{\code{coefficients}}{ the matrix of regression coefficients.}
//' \item{\code{bo}}{ a matrix of one row containing the intercepts for each component.}
//' \item{\code{scores}}{ the matrix of scores.}
//' \item{\code{X_loadings}}{ the matrix of X loadings.}
//' \item{\code{Y_loadings}}{ the matrix of Y loadings.}
//' \item{\code{vip}}{ the projection matrix.}
//' \item{\code{selectivity_ratio}}{ the matrix of selectivity ratio (see Rajalahti, Tarja, et al. 2009).}
//' \item{\code{Y}}{ the \code{Y} input.}
//' \item{\code{variance}}{ a \code{list} conating two objects: \code{x_var} and \code{y_var}. 
//' These objects contain information on the explained variance for the \code{X} and \code{Y} matrices respectively.}
//' \item{\code{transf}}{ a \code{list} conating two objects: \code{Xcenter} and \code{Xscale}}. 
//' \item{\code{weights}}{ the matrix of wheights.}} 
//' @author Leonardo Ramirez-Lopez
//' @keywords internal 
//' @useDynLib resemble
// [[Rcpp::export]]
List opls_get_all(arma::mat X, 
                  arma::mat Y, 
                  int ncomp,
                  bool scale,            
                  double maxiter,
                  double tol){
  
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
  double xvar = sum(pow(get_column_sds(Xpls), 2));  
  
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
  arma::mat cx;
  arma::mat cy;
  arma::mat projection_matrix;
  arma::mat sratio = arma::zeros(weights.n_rows, weights.n_cols);
  
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
      // 1.1 Compute the 'scaling factor'
      cr = sqrt(trans(iypls) * Xpls * trans(Xpls) * iypls);
      // 1.2 The weights are computed as the cross product of
      // X0 and Y0 divided by the 'scaling factor'...
      w = (trans(Xpls) * iypls) / repmat(cr, Xpls.n_cols, 1);
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
    explained_var(0,i) = arma::var(scores.col(i));
    explained_var(1,i) = explained_var(0,i)/xvar;
    explained_var(2,i) = sum(explained_var.row(0)) / xvar;
  }
  // convert this to standard deviation
  explained_var.row(0) = sqrt(explained_var.row(0));
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
//'      tol)
//' @param X a matrix of predictor variables.
//' @param Y a matrix of either a single or multiple response variables.
//' @param ncomp the number of pls components.
//' @param scale logical indicating whether \code{X} must be scaled.
//' @param maxiter maximum number of iterations.
//' @param tol limit for convergence of the algorithm in the nipals algorithm.
//' @return a list containing the following elements:
//' \itemize{
//' \item{\code{coefficients}}{ the matrix of regression coefficients.}
//' \item{\code{bo}}{ a matrix of one row containing the intercepts for each component.}
//' \item{\code{scores}}{ the matrix of scores.}
//' \item{\code{X_loadings}}{ the matrix of X loadings.}
//' \item{\code{Y_loadings}}{ the matrix of Y loadings.}
//' \item{\code{projection_mat}}{ the projection matrix.}
//' \item{\code{Y}}{ the \code{Y} input.}
//' \item{\code{transf}}{ a \code{list} conating two objects: \code{Xcenter} and \code{Xscale}}. 
//' \item{\code{weights}}{ the matrix of wheights.}} 
//' @author Leonardo Ramirez-Lopez
//' @keywords internal 
//' @useDynLib resemble
// [[Rcpp::export]]
List opls(arma::mat X, 
          arma::mat Y, 
          int ncomp,
          bool scale,            
          double maxiter,
          double tol){
  
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
      // 1.1 Compute the 'scaling factor'
      cr = sqrt(trans(iypls) * Xpls * trans(Xpls) * iypls);
      // 1.2 The weights are computed as the cross product of
      // X0 and Y0 divided by the 'scaling factor'...
      w = (trans(Xpls) * iypls) / repmat(cr, Xpls.n_cols, 1);
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

//' @title fast orthogonal scores algorithn of partial leat squares (opls)
//' @description Computes orthogonal socres partial least squares (opls) 
//' regressions with the NIPALS algorithm. It allows multiple response variables. 
//' In contrast to code{opls} function, this one does not compute unnecessary 
//' data for (local) regression.
//' For internal use only!
//' @usage 
//' opls_get_basics(X, Y, ncomp, scale, 
//'                 maxiter, tol)
//' @param X a matrix of predictor variables.
//' @param Y a matrix of either a single or multiple response variables.
//' @param ncomp the number of pls components.
//' @param scale logical indicating whether \code{X} must be scaled.
//' @param maxiter maximum number of iterations.
//' @param tol limit for convergence of the algorithm in the nipals algorithm.
//' @return a list containing the following elements:
//' \itemize{
//' \item{\code{coefficients}}{ the matrix of regression coefficients.}
//' \item{\code{bo}}{ a matrix of one row containing the intercepts for each component.}
//' \item{\code{Y_loadings}}{ the matrix of Y loadings.}
//' \item{\code{projection_mat}}{ the projection matrix.}
//' \item{\code{transf}}{ a \code{list} conating two objects: \code{Xcenter} and \code{Xscale}}. 
//' } 
//' @author Leonardo Ramirez-Lopez
//' @keywords internal 
//' @useDynLib resemble
// [[Rcpp::export]]
List opls_get_basics(arma::mat X, 
                     arma::mat Y, 
                     int ncomp,
                     bool scale,            
                     double maxiter,
                     double tol
){
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
  arma::vec imsd;
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
  
  for (int i = 0; i < ncomp; i++) {
    Yplsb = Ypls;
    // Select the Y variable with the largest standard deviation
    imsd = get_col_largest_sd(Ypls);
    iypls = Ypls.col(imsd[0]);
    previous_ts.fill(0);
    
    j = 0;
    keepg = true;
    
    while (keepg) {
      if(j > 0) {
        previous_ts = ts;
      }
      //Step 1: Compute a vector of loading weights...
      // 1.1 Compute the 'scaling factor'
      cr = sqrt(trans(iypls) * Xpls * trans(Xpls) * iypls);
      // 1.2 The weights are computed as the cross product of
      // X0 and Y0 divided by the 'scaling factor'...
      w = (trans(Xpls) * iypls) / repmat(cr, Xpls.n_cols, 1);
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
  
  projection_matrix = trans(weights) * arma::solve(Xloadings * trans(weights), arma::eye(Xloadings.n_rows, Xloadings.n_rows));
  
  arma::mat ymean = arma::mean(Y);
  arma::vec ymean_vec = arma::vectorise(ymean);
  arma::mat y_hat_mean;
  arma::mat y_hat_mean_vec;
  int idx = 0;
  for (int k = 0; k < ny; k++) {
    arma::mat jth_loading = Yloadings.col(k);
    for (int j = 0; j < ncomp; j++) {
      coefficients.col(idx) = projection_matrix.cols(0, j) * jth_loading.rows(0, j);
      y_hat_mean = x_center_vec * coefficients.col(idx);
      y_hat_mean_vec = arma::vectorise(y_hat_mean);
      bo(k,j) = ymean_vec(k) - y_hat_mean_vec(0);
      idx = idx + 1;
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
    _["weights"] = weights
  );
}


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
//' @useDynLib resemble
// [[Rcpp::export]]
Rcpp::NumericMatrix predict_opls(arma::mat bo, 
                                 arma::mat b, 
                                 int ncomp, 
                                 arma::mat newdata,
                                 bool scale,
                                 arma::mat Xscale
){

  if (scale) {
    newdata = newdata / arma::repmat(Xscale, newdata.n_rows, 1);
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
//' @useDynLib resemble
// [[Rcpp::export]]
Rcpp::NumericMatrix project_opls(arma::mat projection_mat, 
                                 int ncomp, 
                                 arma::mat newdata,
                                 bool scale,
                                 arma::mat Xcenter,
                                 arma::mat Xscale
){
  
  if(scale){
    newdata = newdata / arma::repmat(Xscale, newdata.n_rows, 1);
  }
  
  //Necessary to center
  newdata = newdata - arma::repmat(Xcenter, newdata.n_rows, 1);
  
  arma::mat proj = newdata * projection_mat.cols(0, ncomp - 1);
  
  return Rcpp::wrap(proj);
}

//' @title Projection to pls and then re-construction
//' @description Projects spectra onto a PLS space and then reconstructs it back.
//' @usage project_opls(projection_mat, ncomp, newdata, scale, Xcenter, Xscale)
//' @param x a matrix to project.
//' @param projection_mat the projection matrix generated by the \code{opls_get_basics} function.
//' @param xloadings the loadings matrix generated by the \code{opls_get_basics} function.
//' @return a matrix of 1 row and 1 column.
//' @author Leonardo Ramirez-Lopez
//' @keywords internal 
//' @useDynLib resemble
// [[Rcpp::export]]
Rcpp::NumericMatrix reconstruction_error(arma::mat x, 
                                         arma::mat projection_mat, 
                                         arma::mat xloadings){
  arma::mat xrec = x;
  arma::mat xrmse;
  xrec = x * projection_mat * xloadings; 
  xrmse = sqrt(arma::mean(arma::mean(pow(x - xrec, 2), 0), 1));
  return Rcpp::wrap(xrmse);
}


//' @title Internal Cpp function for computing the weights of the PLS components 
//' necessary for weighted average PLS
//' @description For internal use only!. 
//' @usage
//' get_pls_weights(projection_mat, 
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
//' @param Xcenter a matrix of one row with the values that must be used for centering \code{newdata}.
//' @param Xscale if \code{scale = TRUE} a matrix of one row with the values that must be used for scaling \code{newdata}.
//' @return a matrix of one row with the weights for each component between the max. and min. specified. 
//' @author Leonardo Ramirez-Lopez
//' @keywords internal 
//' @useDynLib resemble
// [[Rcpp::export]]
Rcpp::NumericMatrix get_pls_weights(arma::mat projection_mat, 
                                    arma::mat xloadings,
                                    arma::mat coefficients,
                                    arma::mat new_x,
                                    int min_component, 
                                    int max_component, 
                                    bool scale,
                                    arma::mat Xcenter,
                                    arma::mat Xscale
){
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
    arma::mat xrec = sc.cols(0,i) * xloadings.rows(0,i);
    xrmsres.col(i) = sqrt(arma::mean(arma::mean(pow(Xz - xrec, 2), 0), 1));
  }
  
  arma::mat rmsb = sqrt(get_column_means(pow(coefficients.cols(0, max_component - 1), 2)));
  arma::mat rmsb_x = trans(rmsb.rows(min_component - 1, max_component - 1)) % xrmsres.cols(min_component - 1, max_component - 1);
  arma::mat whgtn = pow(rmsb_x, -1);
  whgt  = whgtn / arma::repmat(sum(whgtn, 1), 1, whgtn.n_cols);
  return Rcpp::wrap(whgt);
}

//' @title Internal Cpp function for performing leave-group-out cross-validations for pls regression 
//' @description For internal use only!. 
//' @usage opls_cv_cpp(X, Y, scale, method, 
//'                   mindices, pindices, 
//'                   min_component, ncomp, 
//'                   new_x, 
//'                   maxiter, tol, 
//'                   wapls_grid)
//' @param X a matrix of predictor variables.
//' @param Y a matrix of a single response variable.
//' @param scale a logical indicating whether the matrix of predictors (\code{X}) must be scaled.
//' @param method the method used for regression. One of the following options: \code{'pls'} or \code{'wapls'} or \code{'completewapls1p'}.
//' @param mindices a matrix with \code{n} rows and \code{m} columns where \code{m} is equivalent to the number of 
//' resampling iterations. The elements of each column indicate the indices of the observations to be used for modeling at each 
//' iteration.
//' @param pindices a matrix with \code{k} rows and \code{m} columns where \code{m} is equivalent to the number of 
//' resampling iterations. The elements of each column indicate the indices of the observations to be used for predicting at each 
//' iteration.
//' @param min_component an integer indicating the number of minimum pls components (if the \code{method = 'pls'}).
//' @param ncomp an integer indicating the number of pls components.
//' @param new_x a matrix of one row corresponding to the observation to be predicted (if the \code{method = 'wapls'}).
//' @param maxiter maximum number of iterations.
//' @param tol limit for convergence of the algorithm in the nipals algorithm.
//' @param wapls_grid the grid on which the search for the best combination of minimum and maximum pls factors of \code{'wapls'} is based on in case \code{method = 'completewapls1p'}.
//' @return a list containing the following one-row matrices:
//' \itemize{
//' \item{\code{rmse_seg}}{ the RMSEs.}
//' \item{\code{st_rmse_seg}}{ the standardized RMSEs.}
//' \item{\code{rsq_seg}}{ the coefficients of determination.}
//' } 
//' @author Leonardo Ramirez-Lopez
//' @keywords internal 
//' @useDynLib resemble
// [[Rcpp::export]]
List opls_cv_cpp(arma::mat X, 
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
                 arma::mat wapls_grid
){
  arma::mat rmseseg;
  arma::mat strmseseg;
  arma::mat rsqseg;
  
  arma::mat compweights;
  arma::mat crcompweights;
  
  if(method == "pls"){
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
      
      List fit = Rcpp::as<Rcpp::List>(opls_get_basics(xmatslice, ymatslice, ncomp, scale, maxiter, tol));
      
      transf = fit["transf"];   
      
      arma::mat ypred;
      
      ypred = Rcpp::as<arma::mat>(predict_opls(fit["bo"], 
                                               fit["coefficients"], 
                                                  ncomp, 
                                                  pxmatslice,
                                                  scale,
                                                  transf["Xscale"]));
      
      arma::mat rdl = sqrt(get_column_means(pow(rpymatslice - ypred, 2)));
      rmseseg.col(i) = rdl;
      arma::mat mimav = arma::zeros(1,1);
      mimav.col(0) = max(pymatslice) - min(pymatslice);
      strmseseg.col(i) = rmseseg.col(i) / arma::repmat(mimav, ncomp, 1);
      rsqseg.col(i) = pow(arma::cor(ypred, pymatslice), 2);
    }
  }
  
  if(method == "wapls"){
    
    rmseseg = arma::zeros(1, mindices.n_cols);
    strmseseg = arma::zeros(1, mindices.n_cols);
    rsqseg = arma::zeros(1, mindices.n_cols);
    
    // define the wapls weights directly here
    List cfit = Rcpp::as<Rcpp::List>(opls(X, Y, ncomp, scale, maxiter, tol));
    List ctransf = cfit["transf"];
    
    compweights = arma::zeros(1, ncomp);
    compweights.cols(min_component-1, ncomp-1) =  Rcpp::as<arma::mat>(get_pls_weights(cfit["projection_mat"], 
                                                                      cfit["X_loadings"],
                                                                          cfit["coefficients"],
                                                                              new_x,
                                                                              min_component, 
                                                                              ncomp, 
                                                                              scale,
                                                                              ctransf["Xcenter"],
                                                                                     ctransf["Xscale"]));
    
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
      
      
      List fit = Rcpp::as<Rcpp::List>(opls(xmatslice, ymatslice, ncomp, scale, maxiter, tol));
      
      transf = fit["transf"];   
      
      arma::mat nypred;
      arma::mat ypred;
      
      nypred = Rcpp::as<arma::mat>(predict_opls(fit["bo"], 
                                                fit["coefficients"], 
                                                   ncomp, 
                                                   pxmatslice,
                                                   scale,
                                                   transf["Xscale"]));
      ypred = arma::sum(rcompweights % nypred, 1);
      
      arma::mat rdl = sqrt(get_column_means(pow(pymatslice - ypred, 2)));
      rmseseg.col(i) = rdl;
      arma::mat mimav = arma::zeros(1,1);
      mimav.col(0) = max(pymatslice) - min(pymatslice);
      strmseseg.col(i) = rmseseg.col(i) / mimav;
      rsqseg.col(i) = pow(arma::cor(ypred, pymatslice), 2);
    }
  }
  
  if(method == "completewapls1"){
    rmseseg = arma::zeros(wapls_grid.n_rows, mindices.n_cols);
    strmseseg = arma::zeros(wapls_grid.n_rows, mindices.n_cols);
    rsqseg = arma::zeros(wapls_grid.n_rows, mindices.n_cols);
    
    // define the wapls weights directly here
    List cfit = Rcpp::as<Rcpp::List>(opls_get_basics(X, Y, ncomp, scale, maxiter, tol));
    List ctransf = cfit["transf"];
    
    compweights = arma::zeros(1, ncomp);
    compweights.cols(min_component-1, ncomp-1) =  Rcpp::as<arma::mat>(get_pls_weights(cfit["projection_mat"], 
                                                                      cfit["X_loadings"],
                                                                          cfit["coefficients"],
                                                                              new_x,
                                                                              min_component, 
                                                                              ncomp, 
                                                                              scale,
                                                                              ctransf["Xcenter"],
                                                                                     ctransf["Xscale"]));
    
    crcompweights = arma::zeros(wapls_grid.n_rows, ncomp); 
    for(int i = 0; (unsigned)i < crcompweights.n_rows; i++){
      int minpls = wapls_grid(i,0);
      int maxpls = wapls_grid(i,1);
      arma::mat subw = arma::zeros(1, ncomp); 
      subw.cols(minpls - 1, maxpls - 1) = compweights.cols(minpls - 1, maxpls - 1);
      arma::mat sumsubw = arma::repmat(sum(subw, 1), 1, ncomp);
      crcompweights.row(i) = subw/sumsubw;
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
      
      
      List fit = Rcpp::as<Rcpp::List>(opls_get_basics(xmatslice, ymatslice, ncomp, scale, maxiter, tol));
      
      transf = fit["transf"];   
      
      arma::mat nypred;
      arma::mat ypred;
      arma::mat rpymatslice;
      rpymatslice = arma::repmat(pymatslice, 1, wapls_grid.n_rows);
      
      
      nypred = (Rcpp::as<arma::mat>(predict_opls(fit["bo"], 
                                                 fit["coefficients"], 
                                                    ncomp, 
                                                    pxmatslice,
                                                    scale,
                                                    transf["Xscale"])));
      
      ypred = nypred * trans(crcompweights);
      
      //ypred = arma::sum(rcompweights % nypred, 1);
      
      arma::mat rdl = sqrt(get_column_means(pow(rpymatslice - ypred, 2)));
      rmseseg.col(i) = rdl;
      arma::mat mimav = arma::zeros(1,1);
      mimav.col(0) = max(pymatslice) - min(pymatslice);
      strmseseg.col(i) = rmseseg.col(i) / arma::repmat(mimav, wapls_grid.n_rows, 1);
      rsqseg.col(i) = pow(arma::cor(ypred, pymatslice), 2);
    }
  }
  
  
  // here all the weights are output from 1 to ncomp (if method == wapls)
  // zeroes are assigned to those which are not selected at the begining
  return Rcpp::List::create(
    Rcpp::Named("rmse_seg") = rmseseg,
    Rcpp::Named("st_rmse_seg") = strmseseg,
    Rcpp::Named("rsq_seg") = rsqseg,
    Rcpp::Named("compweights") = compweights,
    Rcpp::Named("crcompweights") = crcompweights
  );
}


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
//' \item{\code{Xz}}{ the (final transformed) matrix of predictor variables.}
//' \item{\code{alpha}}{ the alpha matrix.}
//' \item{\code{is.scaled}}{ logical indicating whether both the predictors and response variable were scaled to zero mean and unit variance.}
//' \item{\code{Xcenter}}{ if matrix of predictors was scaled, the centering vector used for \code{X}.}
//' \item{\code{Xscale}}{ if matrix of predictors was scaled, the scaling vector used for \code{X}.}
//' \item{\code{Ycenter}}{ if matrix of predictors was scaled, the centering vector used for \code{Y}.}
//' \item{\code{Yscale}}{ if matrix of predictors was scaled, the scaling vector used for \code{Y}.}
//' }
//' @author Leonardo Ramirez-Lopez
//' @keywords internal 
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
  
  return Rcpp::List::create(
    Rcpp::Named("Xz") = Xz,
    Rcpp::Named("alpha") = alpha,
    Rcpp::Named("is_scaled") = scale,
    Rcpp::Named("Xcenter") = xc,
    Rcpp::Named("Xscale") = xs,
    Rcpp::Named("Ycenter") = yc,
    Rcpp::Named("Yscale") = ys
  );
}


//' @title Prediction function for the \code{gaussian_process} function (Gaussian process regression with dot product covariance)
//' @description Predicts response values based on a model generated by the \code{gaussian_process} function (Gaussian process regression with dot product covariance). For internal use only!. 
//' @usage predict_gaussian_process(Xz, alpha, newdata, scale, Xcenter, Xscale, Ycenter, Yscale)
//' @param Xz the final (scaled?) matrix of predictors used to create the regression model in the \code{gaussian_process} function
//' @param alpha the alpha matrix corresponding to the regression model in the \code{gaussian_process} function
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
    newdatatr = newdatatr - arma::repmat(Xcenter, newdata.n_rows, 1);
    newdatatr = newdatatr / arma::repmat(Xscale, newdata.n_rows, 1);
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
//' @usage gaussian_process_cv(X, Y, mindices, pindices, noisev = 0.001, scale)
//' @param X a matrix of predictor variables.
//' @param Y a matrix of a single response variable.
//' @param mindices a matrix with \code{n} rows and \code{m} columns where \code{m} is equivalent to the number of 
//' resampling iterations. The elements of each column indicate the indices of the observations to be used for modeling at each 
//' iteration.
//' @param pindices a matrix with \code{k} rows and \code{m} columns where \code{m} is equivalent to the number of 
//' resampling iterations. The elements of each column indicate the indices of the observations to be used for predicting at each 
//' iteration.
//' @param ncomp an integer indicating the number of pls components.
//' @param scale a logical indicating whether both the predictors 
//' and the response variable must be scaled to zero mean and unit variance.
//' @return a list containing the following one-row matrices:
//' \itemize{
//' \item{\code{rmse.seg}}{ the RMSEs.}
//' \item{\code{st.rmse.seg}}{ the standardized RMSEs.}
//' \item{\code{rsq.seg}}{ the coefficients of determination.}
//' } 
//' @author Leonardo Ramirez-Lopez
//' @keywords internal 
//' @useDynLib resemble
// [[Rcpp::export]]
List gaussian_process_cv(arma::mat X, 
                         arma::mat Y, 
                         arma::mat mindices,
                         arma::mat pindices,
                         float noisev = 0.001,
                         bool scale = true
){
  
  arma::mat rmseseg = arma::zeros(1, mindices.n_cols);
  arma::mat strmseseg = arma::zeros(1, mindices.n_cols);
  arma::mat rsqseg = arma::zeros(1, mindices.n_cols);
  
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
    
    
    arma::mat rdl = sqrt(get_column_means(pow(pymatslice - ypred, 2)));
    rmseseg.col(i) = rdl;
    arma::mat mimav = arma::zeros(1,1);
    mimav.col(0) = max(pymatslice) - min(pymatslice);
    strmseseg.col(i) = rmseseg.col(i) / mimav;
    rsqseg.col(i) = pow(arma::cor(ypred, pymatslice), 2);
  }
  return Rcpp::List::create(
    Rcpp::Named("rmse_seg") = rmseseg,
    Rcpp::Named("st_rmse_seg") = strmseseg,
    Rcpp::Named("rsq_seg") = rsqseg
  );
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
//' @param Y a matrix of either a single or multiple response variables.
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
//' \item{\code{pc_scores}}{ a matrix of principal component scores.}
//' \item{\code{pc_loadings}}{ a matrix of of principal component loadings.}
//' \item{\code{variance}}{ a matrix of the variance of the principal components.} 
//' \item{\code{scale}}{ a \code{list} conating two objects: \code{center} and \code{scale}, which correspond to the vectors used to center and scale the input matrix.} 
//' } 
//' @author Leonardo Ramirez-Lopez
//' @keywords internal 
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
  double xvar = sum(pow(get_column_sds(Xpls), 2));  
  // int max_pcs = arma::min(arma::size(Xpls));
  
  // matrices to declare
  int iter;
  bool keepg;
  arma::mat pp;
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
    Xpls = Xpls - (tt * trans(pp_std));
    pc_scores.col(i) = tt.col(0);
    pc_loadings.col(i) = pp_std.col(0);
    explained_var(0,i) = arma::var(pc_scores.col(i));
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
    //pc_indices = find(explained_var.row(2) <= pcSelvalue && explained_var.row(2) > 0); 
    pc_indices = find(explained_var.row(2) > 0); 
    pc_scores = pc_scores.cols(pc_indices);
    pc_loadings = pc_loadings.cols(pc_indices);
    explained_var = explained_var.cols(pc_indices);
  }
  // convert this to standard deviation
  explained_var.row(0) = sqrt(explained_var.row(0));
  
  return Rcpp::List::create(
    Rcpp::Named("pc_indices") = pc_indices,
    Rcpp::Named("pc_scores") = pc_scores,
    Rcpp::Named("pc_loadings") = trans(pc_loadings),
    Rcpp::Named("pc_variance") = explained_var,
    Rcpp::Named("scale") = Rcpp::List::create(
      Rcpp::Named("center") = x_center_vec,
      Rcpp::Named("scale") = x_scale_vec
    )
  );
}
