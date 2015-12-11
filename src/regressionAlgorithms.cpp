#include <RcppArmadillo.h>
#include <math.h>
#include <iostream>
#ifdef _OPENMP
#include <omp.h>    // OpenMP
#endif

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;


//' @title Function for identifiying the column in a \code{matrix} with the largest standard deviation
//' @description Identifies the column with the largest standard deviation. For internal use only!
//' @usage wcolSds(X)
//' @param X a \code{matrix}.
//' @return a value indicating the index of the column with the largest standard deviation. 
//' @useDynLib resemble
//' @author Leonardo Ramirez-Lopez
//' @keywords internal 
//' @useDynLib resemble
// [[Rcpp::export]]
NumericVector wcolSds(arma::mat X){
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


//' @title Function for computing the standard deviation of each column in a \code{matrix}
//' @description Computes the standard deviation of each column in a \code{matrix}. For internal use only!
//' @usage colSds(X)
//' @param X a a \code{matrix}.
//' @return a vector of standard deviation values. 
//' @useDynLib resemble
//' @author Leonardo Ramirez-Lopez
//' @keywords internal 
//' @useDynLib resemble
// [[Rcpp::export]]
NumericVector colSds(arma::mat X){
  //arma::mat Xz(X.begin(), X.nrow(), X.ncol(), false); 
  //   int nx = X.n_cols;
  //   arma::mat sds = arma::zeros(1, X.n_cols);
  //   
  //   // Compute the standard deviations
  //   for (int k = 0; k < nx; k++) {
  //     sds.col(k) = arma::stddev(X.col(k));
  //   }
  arma::mat sds = arma::stddev(X, 0, 0);
  return Rcpp::wrap(sds);
}

//' @title Function for computing the mean of each column in a \code{matrix}
//' @description Computes the mean of each column in a \code{matrix}. For internal use only!
//' @usage colSds(X)
//' @param X a a \code{matrix}.
//' @return a vector of mean values. 
//' @useDynLib resemble
//' @author Leonardo Ramirez-Lopez
//' @keywords internal 
//' @useDynLib resemble
// [[Rcpp::export]]
NumericVector cms(arma::mat X){
  //arma::mat Xz(X.begin(), X.nrow(), X.ncol(), false); 
  //   int nx = X.n_cols;
  //   arma::mat mn = arma::zeros(1, X.n_cols);
  //   
  //   // Compute the standard deviations
  //   for (int j = 0; j < nx; j++) {
  //     mn.col(j) = arma::mean(X.col(j));
  //   }
  arma::mat mn = arma::mean(X, 0);
  return Rcpp::wrap(mn);
}



//' @title orthogonal scores algorithn of partial leat squares (opls)
//' @description Computes orthogonal socres partial least squares (opls) regressions with the NIPALS algorithm. It allows multiple response variables. 
//' For internal use only!
//' @usage 
//' opls(X, Y, ncomp, scale, 
//'      maxiter, tol, 
//'      regression = TRUE, 
//'      pcSelmethod = "cumvar", 
//'      pcSelvalue = 0.99)
//' @param X a \code{matrix} of predictor variables.
//' @param Y a \code{matrix} of either a single or multiple response variables.
//' @param ncomp the number of pls components.
//' @param scale logical indicating whether \code{X} must be scaled.
//' @param maxiter maximum number of iterations.
//' @param tol limit for convergence of the algorithm in the nipals algorithm.
//' @param regression a logical indicating if the function is being used for regression. Otherwise it is used only for projection. Default is \code{TRUE}.
//' @param pcSelmethod if \code{regression = TRUE}, the method for selecting the number of components. 
//' Options are: \code{'cumvar'} (for selecting the number of principal components based on a given 
//' cumulative amount of explained variance) and \code{"var"} (for selecting the number of principal 
//' components based on a given amount of explained variance). Default is \code{'cumvar'}
//' @param pcSelvalue a numerical value that complements the selected method (\code{pcSelmethod}). 
//' If \code{"cumvar"} is chosen, it must be a value (higher than 0 and lower than 1) indicating the maximum 
//' amount of cumulative variance that the retained components should explain. If \code{"var"} is chosen, 
//' it must be a value (higher than 0 and lower than 1) indicating that components that explain (individually) 
//' a variance lower than this threshold must be excluded. If \code{"manual"} is chosen, it must be a value 
//' specifying the desired number of principal components to retain. Default is 0.99.
//' @return a list containing the following elements:
//' \itemize{
//' \item{\code{coefficients}}{ the \code{matrix} of regression coefficients.}
//' \item{\code{bo}}{ a \code{matrix} of one row containing the intercepts for each component.}
//' \item{\code{scores}}{ the \code{matrix} of scores.}
//' \item{\code{X.loadings}}{ the \code{matrix} of X loadings.}
//' \item{\code{Y.loadings}}{ the \code{matrix} of Y loadings.}
//' \item{\code{projectionM}}{ the projection \code{matrix}.}
//' \item{\code{variance}}{ a \code{list} conating two objects: \code{x.var} and \code{y.var}. 
//' These objects contain information on the explained variance for the \code{X} and \code{Y} matrices respectively.}
//' \item{\code{transf}}{ a \code{list} conating two objects: \code{Xcenter} and \code{Xscale}. 
//' These objects contain information on the explained variance for the \code{X} and \code{Y} matrices respectively.}
//' } 
//' @useDynLib resemble
//' @author Leonardo Ramirez-Lopez
//' @keywords internal 
//' @useDynLib resemble
// [[Rcpp::export]]
List opls(arma::mat X, 
          arma::mat Y, 
          int ncomp,
          bool scale,            
          double maxiter,
          double tol,
          bool regression = true,
          String pcSelmethod = "cumvar",
          double pcSelvalue = 0.99
){
  
  int nPf = ncomp;
  int ny = Y.n_cols;
  int nynf = ncomp * Y.n_cols;
  
  arma::mat weights = arma::zeros(nPf, X.n_cols);
  arma::mat scores = arma::zeros(X.n_rows, nPf);
  arma::mat Xloadings = arma::zeros(nPf, X.n_cols);
  arma::mat Yloadings = arma::zeros(nPf, ny);
  arma::mat coefficients = arma::zeros(X.n_cols, nynf);
  arma::mat bo = arma::zeros(ny, nPf);
  arma::mat exv = arma::zeros(3, nPf);
  arma::mat yex = arma::zeros(ny, nPf);
  
  
  // if false, it reuses memory and avoids extra copy
  // the problem is that the matrix cannot be easily overwriten 
  // retrieves random results
  
  //  arma::mat Xz(X.begin(), X.nrow(), X.ncol(), true);     
  //  arma::mat Yz(Y.begin(), Y.nrow(), Y.ncol() , true); 
  
  arma::mat Xscale;
  arma::mat Xs;
  arma::mat xfcntr;
  arma::mat Xz = X;
  
  if(scale){
    Xscale = arma::repmat(Rcpp::as<arma::mat>(colSds(Xz)), Xz.n_rows, 1);
    Xz = Xz / Xscale;
    Xs =  Xscale.row(0);
  }
  
  xfcntr = Rcpp::as<arma::mat>(cms(Xz));
  
  //Xcenter = arma::repmat(Rcpp::as<arma::mat>(cms(Xz)), Xz.n_rows, 1);
  Xz = Xz - arma::repmat(xfcntr, Xz.n_rows, 1);
  
  
  arma::mat Xpls = Xz;
  arma::mat Ypls = Y;
  
  //variance of Xpls
  double xvar = sum(pow(colSds(Xpls), 2));  
  
  // matrices to declare
  arma::mat iypls;
  arma::mat Yplsb;
  arma::vec imsd;
  int j;
  bool keepg;
  arma::mat lb;
  arma::mat cr;
  arma::mat tsr;
  arma::mat ts;
  arma::mat w;
  arma::mat p;
  arma::mat q;
  arma::mat cx;
  arma::mat cy;
  arma::mat tsrp;
  arma::mat ev;
  arma::mat prjM;
  
  int nff = 0;
  
  for (int i = 0; i < ncomp; i++){
    Yplsb = Ypls;
    Xpls = Xpls;      
    // Select the Y variable with the largest standard deviation
    imsd = wcolSds(Ypls);
    iypls = Ypls.col(imsd[0]);
    arma::mat tsr = arma::zeros(Xpls.n_rows, 1);
    tsrp = tsr + tsr;
    
    j = 0;
    keepg = true;
    
    while(keepg){
      if(j > 0)
      {
        tsrp = ts;
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
      lb = abs(sum((ts - tsrp) / ts));
      keepg = lb[0] > tol;
      j = j + 1;
      if(maxiter <= j){
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
    exv(0,i) = arma::var(scores.col(i));
    exv(1,i) = sum(exv.row(0)) / xvar;
    exv(2,i) = exv(0,i)/xvar;
    
    nff = nff + 1;
    
    if(!regression)
    {
      if(pcSelmethod == "var" || pcSelmethod == "cumvar")
      {
        bool chk;
        if(pcSelmethod == "cumvar"){
          chk = exv(1,i) > pcSelvalue;
        }
        else{
          chk = exv(2,i) < pcSelvalue;
        }
        if(chk)
        {
          ncomp = nff - 1;
          nff = nff - 2;
          if(i == 0) throw exception("With the current value in the 'pcSelection' argument, no components are selected. Try another value.");
          break;
        }
      } 
    }
  }
  
  if(!regression){
    weights = weights.rows(0,nff);
    coefficients = coefficients.cols(0, nff);
    bo = bo.cols(0,nff);
    scores = scores.cols(0,nff);
    Xloadings = Xloadings.rows(0,nff);
    Yloadings = Yloadings.rows(0,nff);
    exv = exv.cols(0,nff);
    yex = yex.cols(0,nff);
  }
  
  
  // convert this to standard deviation
  exv.row(0) = sqrt(exv.row(0));
  
  prjM = trans(weights) * arma::solve(Xloadings * trans(weights), arma::eye(Xloadings.n_rows, Xloadings.n_rows));
  
  arma::mat yexi;
  arma::mat cop;
  for(int i = 0; i < ny; i++){
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
  
  arma::vec ymeanv;
  arma::mat ymean = arma::mean(Y);
  ymeanv = arma::vectorise(ymean);
  
  
  int idx = 0;
  for(int k = 0; k < ny; k++){
    arma::mat jyload = Yloadings.col(k);
    for(int j = 0; j < ncomp; j++){
      coefficients.col(idx) = prjM.cols(0,j) * jyload.rows(0,j);
      arma:: vec leov;
      arma::mat leo = xfcntr * coefficients.col(idx);
      leov = arma::vectorise(leo);
      bo(k,j) = ymeanv(k) - leov(0);
      idx = idx + 1;
    }
  }
  
  return Rcpp::List::create(
    Rcpp::Named("ncomp") = ncomp,
    Rcpp::Named("coefficients") = coefficients,
    Rcpp::Named("bo") = bo,
    Rcpp::Named("scores") = scores,
    Rcpp::Named("X.loadings") = Xloadings,
    Rcpp::Named("Y.loadings") = Yloadings,
    Rcpp::Named("projectionM") = prjM,
    Rcpp::Named("Y") = Y,
    Rcpp::Named("variance") = Rcpp::List::create(
      Rcpp::Named("x.var") = exv,
      Rcpp::Named("y.var") = yex
    ),
    Rcpp::Named("transf") = Rcpp::List::create(
      Rcpp::Named("Xcenter") = xfcntr,
      Rcpp::Named("Xscale") = Xs
    ),
    _["weights"] = weights
  );
}



//' @title Prediction function for the \code{opls} and \code{opls2} functions
//' @description Predicts response values based on a model generated by either by \code{opls} or the \code{opls2} functions. 
//' For internal use only!. 
//' @usage predopls(bo, b, ncomp, newdata, scale, Xscale)
//' @param bo a numeric value indicating the intercept.
//' @param b the \code{matrix} of regression coefficients.
//' @param ncomp an integer value indicating how may components must be used in the prediction.
//' @param newdata a \code{matrix} containing the predictor variables.
//' @param scale a logical indicating whether the matrix of predictors used to create the regression model 
//' (either in \code{opls} or \code{opls2}) was scaled.
//' @param Xscale if \code{scale = TRUE} a \code{matrix} of one row with the values that must be used for scaling \code{newdata}.
//' @return a \code{matrix} of predicted values
//' @useDynLib resemble
//' @author Leonardo Ramirez-Lopez
//' @keywords internal 
//' @useDynLib resemble
// [[Rcpp::export]]
Rcpp::NumericMatrix predopls(arma::mat bo, 
                             arma::mat b, 
                             int ncomp, 
                             arma::mat newdata,
                             bool scale,
                             arma::mat Xscale
){
  //arma::mat Xz(newdata.begin(), newdata.nrow(), newdata.ncol(), true);
  //arma::mat predicted = arma::zeros(Xz.n_rows, ncomp);
  
  arma::mat Xz; 
  
  //Not necessary to center
  //   if(center){
  //     Xz = Xz - arma::repmat(Xcenter, Xz.n_rows, 1);
  //   }
  
  if(scale){
    Xz = newdata / arma::repmat(Xscale, newdata.n_rows, 1);
  }
  else{
    Xz = newdata;
  }
    
  arma::mat predicted = (Xz * b.cols(0, ncomp - 1)) + arma::repmat(bo.cols(0, ncomp - 1), Xz.n_rows, 1);
    
  return Rcpp::wrap(predicted);
}



//' @title Projection function for the \code{opls} function
//' @description Projects new spectra onto a PLS space based on a model generated by either by \code{opls} or the \code{opls2} functions. 
//' For internal use only!. 
//' @usage projectpls(projectionm, ncomp, newdata, scale, Xcenter, Xscale)
//' @param projectionm the projection \code{matrix} generated by the \code{opls} function.
//' @param ncomp an integer value indicating how may components must be used in the prediction.
//' @param newdata a \code{matrix} containing the predictor variables.
//' @param scale a logical indicating whether the matrix of predictors used to create the regression model was scaled.
//' @param Xscale if \code{scale = TRUE} a \code{matrix} of one row with the values that must be used for scaling \code{newdata}.
//' @param Xcenter a \code{matrix} of one row with the values that must be used for centering \code{newdata}.
//' @return a \code{matrix} corresponding to the new spectra projected onto the PLS space 
//' @useDynLib resemble
//' @author Leonardo Ramirez-Lopez
//' @keywords internal 
//' @useDynLib resemble
// [[Rcpp::export]]
Rcpp::NumericMatrix projectpls(arma::mat projectionm, 
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
  
  arma::mat proj = newdata * projectionm.cols(0, ncomp - 1);
  
  return Rcpp::wrap(proj);
}


//' @title Internal Cpp function for computing the weights of the PLS components necessary for weighted average PLS
//' @description For internal use only!. 
//' @usage
//' waplswCpp(projectionm, 
//'           xloadings, 
//'           coefficients, 
//'           newX, 
//'           minF, 
//'           maxF, 
//'           scale, 
//'           Xcenter, 
//'           Xscale)
//' @param projectionm the projection \code{matrix} generated either by the \code{opls} function.
//' @param xloadings .
//' @param coefficients the \code{matrix} of regression coefficients.
//' @param newX a \code{matrix} of one new spectra to be predicted.
//' @param minF an integer indicating the minimum number of pls components.
//' @param maxF an integer indicating the maximum number of pls components.
//' @param scale a logical indicating whether the matrix of predictors used to create the regression model 
//' (either in \code{opls} or \code{opls2}) was scaled.
//' @param Xcenter a \code{matrix} of one row with the values that must be used for centering \code{newdata}.
//' @param Xscale if \code{scale = TRUE} a \code{matrix} of one row with the values that must be used for scaling \code{newdata}.
//' @return a \code{matrix} of one row with the weights for each component between the max. and min. specified. 
//' @useDynLib resemble
//' @author Leonardo Ramirez-Lopez
//' @keywords internal 
//' @useDynLib resemble
// [[Rcpp::export]]
Rcpp::NumericMatrix waplswCpp(arma::mat projectionm, 
                              arma::mat xloadings,
                              arma::mat coefficients,
                              arma::mat newX,
                              int minF, 
                              int maxF, 
                              bool scale,
                              arma::mat Xcenter,
                              arma::mat Xscale
){
  arma::mat Xz = newX;
  arma::mat whgt;
  
  if(scale){
    Xz = Xz / Xscale;
  }
  
  //Necessary to center
  Xz = Xz - Xcenter;
  
  arma::mat xrmsres = arma::zeros(1, maxF);
  
  arma::mat sc = Xz * projectionm.cols(0, maxF - 1);
  for(int i = (minF - 1); i < maxF; i++){
    arma::mat xrec = sc.cols(0,i) * xloadings.rows(0,i);
    xrmsres.col(i) = sqrt(arma::mean(arma::mean(pow(Xz - xrec, 2), 0), 1));
  }
  
  arma::mat rmsb = sqrt(cms(pow(coefficients.cols(0, maxF - 1), 2)));
  arma::mat rmsb_x = trans(rmsb.rows(minF - 1, maxF - 1)) % xrmsres.cols(minF - 1, maxF - 1);
  arma::mat whgtn = pow(rmsb_x, -1);
  whgt  = whgtn / arma::repmat(sum(whgtn, 1), 1, whgtn.n_cols);
  // Another way for computing the weights based on a tricubic function
  // whgt = ((rms.b_x)-min(rms.b_x))/diff(range(rms.b_x))
  // whgt = (1 - (whgt^3))^3
  // whgt[which(whgt == 0)] <- 0.00001
  // whgt = whgt/sum(whgt)
  // here only the selected weights are output from minF to max F
  return Rcpp::wrap(whgt);
}

//' @title Internal Cpp function for performing leave-group-out cross validations for pls regression 
//' @description For internal use only!. 
//' @usage pplscv_cpp(X, Y, scale, method, mindices, pindices, minF, ncomp, newX, maxiter, tol)
//' @param X a \code{matrix} of predictor variables.
//' @param Y a \code{matrix} of a single response variable.
//' @param scale a logical indicating whether the matrix of predictors (\code{X}) must be scaled.
//' @param method the method used for regression. One of the following options: \code{'pls'} or \code{'wapls1'}.
//' @param mindices a \code{matrix} with \code{n} rows and \code{m} columns where \code{m} is equivalent to the number of 
//' resampling iterations. The elements of each column indicate the indices of the samples to be used for modeling at each 
//' iteration.
//' @param pindices a \code{matrix} with \code{k} rows and \code{m} columns where \code{m} is equivalent to the number of 
//' resampling iterations. The elements of each column indicate the indices of the samples to be used for predicting at each 
//' iteration.
//' @param minF an integer indicating the number of minimum pls components (if the \code{method = 'pls'}).
//' @param ncomp an integer indicating the number of pls components.
//' @param newX a \code{matrix} of one row corresponding to the sample to be predicted (if the \code{method = 'wapls1'}).
//' @param maxiter maximum number of iterations.
//' @param tol limit for convergence of the algorithm in the nipals algorithm.
//' @return a list containing the following one-row matrices:
//' \itemize{
//' \item{\code{rmse.seg}}{ the RMSEs.}
//' \item{\code{st.rmse.seg}}{ the standardized RMSEs.}
//' \item{\code{rsq.seg}}{ the coefficients of determination.}
//' } 
//' @useDynLib resemble
//' @author Leonardo Ramirez-Lopez
//' @keywords internal 
//' @useDynLib resemble
// [[Rcpp::export]]
List pplscv_cpp(arma::mat X, 
                arma::mat Y, 
                bool scale,
                String method,
                arma::mat mindices,
                arma::mat pindices,
                int minF,
                int ncomp,
                arma::mat newX,
                double maxiter, 
                double tol
){
  arma::mat rmseseg;
  arma::mat strmseseg;
  arma::mat rsqseg;
  
  arma::mat compweights;
  
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
      
      List fit = Rcpp::as<Rcpp::List>(opls(xmatslice, ymatslice, ncomp, scale, maxiter, tol));
      
      transf = fit["transf"];   
      
      arma::mat ypred;
      
      ypred = Rcpp::as<arma::mat>(predopls(fit["bo"], 
                                            fit["coefficients"], 
                                               ncomp, 
                                               pxmatslice,
                                               scale,
                                               transf["Xscale"]));
      
      arma::mat rdl = sqrt(cms(pow(rpymatslice - ypred, 2)));
      rmseseg.col(i) = rdl;
      arma::mat mimav = arma::zeros(1,1);
      mimav.col(0) = max(pymatslice) - min(pymatslice);
      strmseseg.col(i) = rmseseg.col(i) / arma::repmat(mimav, ncomp, 1);
      rsqseg.col(i) = pow(arma::cor(ypred, pymatslice), 2);
    }
  }
  if(method == "wapls1"){
    
    rmseseg = arma::zeros(1, mindices.n_cols);
    strmseseg = arma::zeros(1, mindices.n_cols);
    rsqseg = arma::zeros(1, mindices.n_cols);
    
    // define the wapls1 weights directly here
    List cfit = Rcpp::as<Rcpp::List>(opls(X, Y, ncomp, scale, maxiter, tol));
    List ctransf = cfit["transf"];
    
    compweights = arma::zeros(1, ncomp);
    compweights.cols(minF-1, ncomp-1) =  Rcpp::as<arma::mat>(waplswCpp(cfit["projectionM"], 
                                                             cfit["X.loadings"],
                                                                 cfit["coefficients"],
                                                                     newX,
                                                                     minF, 
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
      
      nypred = Rcpp::as<arma::mat>(predopls(fit["bo"], 
                                             fit["coefficients"], 
                                                ncomp, 
                                                pxmatslice,
                                                scale,
                                                transf["Xscale"]));
      ypred = arma::sum(rcompweights % nypred, 1);
      
      arma::mat rdl = sqrt(cms(pow(pymatslice - ypred, 2)));
      rmseseg.col(i) = rdl;
      arma::mat mimav = arma::zeros(1,1);
      mimav.col(0) = max(pymatslice) - min(pymatslice);
      strmseseg.col(i) = rmseseg.col(i) / mimav;
      rsqseg.col(i) = pow(arma::cor(ypred, pymatslice), 2);
    }
  }
  
  // here all the weights are output from 1 to ncomp (if method == wapls1)
  return Rcpp::List::create(
    Rcpp::Named("rmse.seg") = rmseseg,
    Rcpp::Named("st.rmse.seg") = strmseseg,
    Rcpp::Named("rsq.seg") = rsqseg,
    Rcpp::Named("compweights") = compweights
  );
}

  
/// Gaussian process regression with linear kernel
//' @title Gaussian process regression with linear kernel (gprdp)
//' @description Carries out a gaussian process regression with a linear kernel (dot product). For internal use only!
//' @usage gprdp(X, Y, noisev, scale) 
//' @param X a matrix of predictor variables
//' @param Y a matrix with a single response variable
//' @param noisev a value indicating the variance of the noise for Gaussian process regression. Default is 0.001. a matrix with a single response variable
//' @param scale a logical indicating whether both the predictors 
//' and the response variable must be scaled to zero mean and unit variance.
//' @return a list containing the following elements:
//' \itemize{
//' \item{\code{Xz}}{ the (final transformed) \code{matrix} of predictor variables.}
//' \item{\code{alpha}}{ the alpha \code{matrix}.}
//' \item{\code{is.scaled}}{ logical indicating whether both the predictors and response variable were scaled to zero mean and unit variance.}
//' \item{\code{Xcenter}}{ if matrix of predictors was scaled, the centering vector used for \code{X}.}
//' \item{\code{Xscale}}{ if matrix of predictors was scaled, the scaling vector used for \code{X}.}
//' \item{\code{Ycenter}}{ if matrix of predictors was scaled, the centering vector used for \code{Y}.}
//' \item{\code{Yscale}}{ if matrix of predictors was scaled, the scaling vector used for \code{Y}.}
//' }
//' @useDynLib resemble
//' @author Leonardo Ramirez-Lopez
//' @keywords internal 
//' @useDynLib resemble
// [[Rcpp::export]]
List gprdp(arma::mat X, 
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
    xc = Rcpp::as<arma::mat>(cms(Xz));
    Xcent = arma::repmat(xc, Xz.n_rows, 1);
    Xz = Xz - Xcent;
    xs = Rcpp::as<arma::mat>(colSds(Xz));
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
    Rcpp::Named("is.scaled") = scale,
    Rcpp::Named("Xcenter") = xc,
    Rcpp::Named("Xscale") = xs,
    Rcpp::Named("Ycenter") = yc,
    Rcpp::Named("Yscale") = ys
  );
}


//' @title Prediction function for the \code{gprdp} function (Gaussian process regression with dot product covariance)
//' @description Predicts response values based on a model generated by the \code{gprdp} function (Gaussian process regression with dot product covariance). For internal use only!. 
//' @usage predgprdp(Xz, alpha, newdata, scale, Xcenter, Xscale, Ycenter, Yscale)
//' @param Xz the final (scaled?) matrix of predictors used to create the regression model in the \code{gprdp} function
//' @param alpha the alpha matrix corresponding to the regression model in the \code{gprdp} function
//' @param newdata a \code{matrix} containing the predictor variables
//' @param scale a logical indicating whether the matrix of predictors used to create the regression model 
//' (in the \code{gprdp} function) was scaled
//' @param Xcenter if \code{center = TRUE} a \code{matrix} of one row with the values that must be used for centering \code{newdata}.
//' @param Xscale if \code{scale = TRUE} a \code{matrix} of one row with the values that must be used for scaling \code{newdata}.
//' @param Ycenter if \code{center = TRUE} a \code{matrix} of one row with the values that must be used for accounting for the centering of the response variable.
//' @param Yscale if \code{scale = TRUE} a \code{matrix} of one row with the values that must be used  for accounting for the scaling of the response variable.
//' @return a \code{matrix} of predicted values
//' @useDynLib resemble
//' @author Leonardo Ramirez-Lopez
//' @keywords internal 
//' @useDynLib resemble
// [[Rcpp::export]]
NumericVector predgprdp(arma::mat Xz, 
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

//' @title Internal Cpp function for performing leave-group-out cross validations for pls regression 
//' @description For internal use only!. 
//' @usage pgpcv_cpp(X, Y, mindices, pindices, noisev = 0.001, scale)
//' @param X a \code{matrix} of predictor variables.
//' @param Y a \code{matrix} of a single response variable.
//' @param mindices a \code{matrix} with \code{n} rows and \code{m} columns where \code{m} is equivalent to the number of 
//' resampling iterations. The elements of each column indicate the indices of the samples to be used for modeling at each 
//' iteration.
//' @param pindices a \code{matrix} with \code{k} rows and \code{m} columns where \code{m} is equivalent to the number of 
//' resampling iterations. The elements of each column indicate the indices of the samples to be used for predicting at each 
//' iteration.
//' @param mindices a \code{matrix} with \code{n} rows and \code{m} columns where \code{m} is equivalent to the number of 
//' resampling iterations. The elements of each column indicate the indices of the samples to be used for modeling at each 
//' iteration.
//' @param pindices a \code{matrix} with \code{k} rows and \code{m} columns where \code{m} is equivalent to the number of 
//' resampling iterations. The elements of each column indicate the indices of the samples to be used for predicting at each 
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
//' @useDynLib resemble
//' @author Leonardo Ramirez-Lopez
//' @keywords internal 
//' @useDynLib resemble
// [[Rcpp::export]]
List pgpcv_cpp(arma::mat X, 
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

    List fit = Rcpp::as<Rcpp::List>(gprdp(xmatslice, ymatslice, noisev, scale));

    arma::mat ypred;
    
    ypred = Rcpp::as<arma::mat>(predgprdp(fit["Xz"], fit["alpha"], pxmatslice, scale, fit["Xcenter"], fit["Xscale"], fit["Ycenter"], fit["Yscale"]));


    arma::mat rdl = sqrt(cms(pow(pymatslice - ypred, 2)));
    rmseseg.col(i) = rdl;
    arma::mat mimav = arma::zeros(1,1);
    mimav.col(0) = max(pymatslice) - min(pymatslice);
    strmseseg.col(i) = rmseseg.col(i) / mimav;
    rsqseg.col(i) = pow(arma::cor(ypred, pymatslice), 2);
  }
  return Rcpp::List::create(
    Rcpp::Named("rmse.seg") = rmseseg,
    Rcpp::Named("st.rmse.seg") = strmseseg,
    Rcpp::Named("rsq.seg") = rsqseg
  );
}
