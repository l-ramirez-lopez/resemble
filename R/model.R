#' @title Global spectral calibration model
#' @name model
#' @aliases model
#' @aliases predict.resemble_model
#' @description
#'
#' \loadmathjax
#'
#' Fits a global calibration model for spectral data using partial least
#' squares (PLS) or Gaussian process regression (GPR). Unlike
#' \code{\link{mbl}}, which builds local models for each prediction,
#' \code{model()} fits a single global model to the reference data.
#'
#' @usage
#' model(Xr, Yr, fit_method, control = model_control(), verbose = TRUE)
#'
#' \method{predict}{resemble_model}(object, newdata, ncomp = NULL, ...)
#'
#' @param Xr A numeric matrix of predictor variables, with observations in rows
#'   and variables in columns. Typically spectral data.
#' @param Yr A numeric matrix with one column containing the response variable
#'   values corresponding to the observations in \code{Xr}.
#' @param fit_method An object created by \code{\link{fit_pls}} or
#'   \code{\link{fit_gpr}} specifying the regression method and its parameters,
#'   including centring and scaling options.
#'   \code{fit_wapls()} is not supported for global models because it requires
#'   target observations to compute weights.
#' @param control A list created by \code{model_control()} specifying the
#'   cross-validation settings. The default is \code{model_control()}.
#' @param verbose Logical indicating whether progress information should be
#'   printed. Default is \code{TRUE}.
#'
#' @param object A \code{resemble_model} object returned by \code{\link{model}}.
#' @param newdata A numeric matrix of new observations with the same number of
#'   columns as the training data.
#' @param ncomp For PLS models, the number of components to use for prediction.
#'   The default is the number used during fitting. Ignored for GPR models. For 
#'   prediction, this can be a vector of integers representing the predictions
#'   coming from models with the requested components. 
#' @param ... Additional arguments, currently unused.
#'
#' @details
#'
#' \code{model()} provides a straightforward interface for fitting global
#' calibration models, in contrast to the local modelling approach implemented
#' in \code{\link{mbl}}. This is useful when:
#' \itemize{
#'   \item{the relationship between spectra and the response is approximately
#'   linear across the full dataset}
#'   \item{a single portable model is needed for deployment}
#'   \item{computational efficiency is prioritised over predictive performance
#'   in heterogeneous datasets}
#' }
#'
#' \subsection{Cross-validation}{
#' When \code{validation_type = "lgo"}, stratified random sampling is used
#' to create training-validation splits that preserve the distribution of the
#' response variable. Cross-validation results include RMSE, R², and
#' standardised RMSE for each component in PLS models, or overall in GPR
#' models.
#' }
#'
#' \subsection{PLS models}{
#' When \code{fit_pls()} is used, the function fits a PLS model with the
#' specified number of components. If cross-validation is enabled, the optimal
#' number of components is selected based on the minimum RMSE.
#' }
#'
#' \subsection{GPR models}{
#' When \code{fit_gpr()} is used, the function fits a Gaussian process
#' regression model with a dot-product covariance function. The noise variance
#' parameter controls regularisation.
#' }
#'
#' @return
#' For \code{model()}, an object of class \code{resemble_model} containing:
#' \itemize{
#'   \item \code{fit_method}: Fitting method object used.
#'   \item \code{control}: Model control object used.
#'   \item \code{model}: Fitted model object.
#'   \item \code{cv_results}: Cross-validation results, if
#'     \code{validation_type = "lgo"}.
#'   \item \code{n_obs}: Number of observations used.
#'   \item \code{n_vars}: Number of predictor variables.
#' }
#'
#' For \code{predict.resemble_model()}, a numeric matrix of predictions. For
#' PLS models, columns correspond to different numbers of components.
#'
#' @author
#' \href{https://orcid.org/0000-0002-5369-5120}{Leonardo Ramirez-Lopez}
#'
#' @seealso
#' \code{\link{fit_pls}} and \code{\link{fit_gpr}} for specifying the fitting
#' method;
#' \code{\link{model_control}} for controlling aspects of the modelling
#' process;
#' \code{\link{mbl}} for memory-based (local) learning.
#'
#' @examples
#' \dontrun{
#' library(prospectr)
#' data(NIRsoil)
#'
#' # Preprocess spectra
#' Xr <- savitzkyGolay(NIRsoil$spc, m = 1, p = 2, w = 11)
#' Yr <- NIRsoil$CEC
#'
#' # Remove missing values
#' ok <- !is.na(Yr)
#' Xr <- Xr[ok, ]
#' Yr <- as.matrix(Yr[ok])
#'
#' # Fit a PLS model with 10 components and cross-validation
#' # Scaling is controlled via fit_pls()
#' pls_mod <- model(
#'   Xr = Xr,
#'   Yr = Yr,
#'   fit_method = fit_pls(ncomp = 10, scale = FALSE),
#'   control = model_control(validation_type = "lgo", number = 10)
#' )
#'
#' # View cross-validation results
#' pls_mod$cv_results
#'
#' # Fit a GPR model (centring/scaling controlled via fit_gpr())
#' gpr_mod <- model(
#'   Xr = Xr,
#'   Yr = Yr,
#'   fit_method = fit_gpr(noise_variance = 0.001, scale = TRUE),
#'   control = model_control(validation_type = "lgo")
#' )
#' }
#'
#' @references
#' de Jong, S. (1993). SIMPLS: An alternative approach to partial least
#' squares regression. \emph{Chemometrics and Intelligent Laboratory Systems},
#' 18(3), 251--263.
#'
#' Rasmussen, C. E., & Williams, C. K. I. (2006). \emph{Gaussian Processes for
#' Machine Learning}. MIT Press.
#'
#' Shenk, J. S., & Westerhaus, M. O. (1991). Population structuring of near
#' infrared spectra and modified partial least squares regression.
#' \emph{Crop Science}, 31(6), 1548--1555.
#' 
# =============================================================================
# model
# =============================================================================

#' @export
model <- function(
    Xr, 
    Yr,
    fit_method,
    control = model_control(),
    verbose = TRUE
) {
  
  call_f <- match.call()
  
  # ---------------------------------------------------------------------------
  # Input validation
  # ---------------------------------------------------------------------------
  if (!is.matrix(Xr)) {
    Xr <- as.matrix(Xr)
  }
  
  if (!is.matrix(Yr)) {
    Yr <- as.matrix(Yr)
  }
  
  if (nrow(Xr) != nrow(Yr)) {
    stop("'Xr' and 'Yr' must have the same number of rows.", call. = FALSE)
  }
  
  if (ncol(Yr) != 1L) {
    stop("'Yr' must have exactly one column.", call. = FALSE)
  }
  
  if (anyNA(Xr)) {
    stop("'Xr' contains missing values.", call. = FALSE)
  }
  
  if (anyNA(Yr)) {
    stop("'Yr' contains missing values.", call. = FALSE)
  }
  
  if (!inherits(fit_method, "fit_method")) {
    stop("'fit_method' must be a fit_*() object (e.g., fit_pls(), fit_gpr()).",
         call. = FALSE)
  }
  
  if (inherits(fit_method, "fit_wapls")) {
    stop(
      "'fit_wapls()' is not supported for global models. ",
      "Use 'fit_pls()' or 'fit_gpr()' instead.",
      call. = FALSE
    )
  }
  
  if (!inherits(control, "model_control")) {
    stop("'control' must be created with model_control().", call. = FALSE)
  }
  
  n_obs  <- nrow(Xr)
  n_vars <- ncol(Xr)
  
  # ---------------------------------------------------------------------------
  # Dispatch to appropriate fitting function
  # ---------------------------------------------------------------------------
  if (inherits(fit_method, "fit_pls")) {
    result <- .fit_model_pls(
      Xr         = Xr,
      Yr         = Yr,
      fit_method = fit_method,
      control    = control,
      verbose    = verbose
    )
  } else if (inherits(fit_method, "fit_gpr")) {
    result <- .fit_model_gpr(
      Xr         = Xr,
      Yr         = Yr,
      fit_method = fit_method,
      control    = control,
      verbose    = verbose
    )
  } else {
    stop("Unsupported 'fit_method' type.", call. = FALSE)
  }
  
  # ---------------------------------------------------------------------------
  # Assemble output
  # ---------------------------------------------------------------------------
  out <- list(
    fit_method = fit_method,
    control    = control,
    model      = result$model,
    cv_results = result$cv_results,
    n_obs      = n_obs,
    n_vars     = n_vars
  )
  attr(out, "call") <- call_f
  class(out) <- "resemble_model"
  out
}


# =============================================================================
# Internal: PLS model fitting
# =============================================================================

.fit_model_pls <- function(
    Xr, Yr, fit_method, control, verbose
) {
  ncomp    <- fit_method$ncomp
  method   <- fit_method$method
  scale    <- fit_method$scale
  max_iter <- fit_method$max_iter
  tol      <- fit_method$tol
  do_cv    <- control$validation_type == "lgo"
  
  cv_results <- NULL
  
  # Cross-validation
  if (do_cv) {
    if (verbose) message("Running cross-validation...")
    
    cv_samples <- sample_stratified(
      y           = Yr,
      p           = control$p,
      number      = control$number,
      group       = NULL,
      replacement = FALSE,
      seed        = NULL
    )
    
    cv_raw <- opls_cv_cpp(
      X             = Xr,
      Y             = Yr,
      scale         = scale,
      method        = "pls",
      mindices      = cv_samples$hold_in,
      pindices      = cv_samples$hold_out,
      min_component = 1L,
      ncomp         = ncomp,
      new_x         = matrix(0, 1, 1),
      maxiter       = max_iter,
      tol           = tol,
      wapls_grid    = matrix(0, 0, 0),
      algorithm     = method,
      statistics    = TRUE
    )
    
    cv_results <- data.frame(
      ncomp      = seq_len(ncomp),
      rmse       = rowMeans(cv_raw$rmse_seg),
      rmse_sd    = apply(cv_raw$rmse_seg, 1, sd),
      st_rmse    = rowMeans(cv_raw$st_rmse_seg),
      r2         = rowMeans(cv_raw$rsq_seg),
      row.names  = NULL
    )
    
    best_ncomp <- which.min(cv_results$rmse)
    cv_results$optimal <- seq_len(ncomp) == best_ncomp
  }
  
  # Fit final model
  if (verbose) message("Fitting model...")
  
  model_fit <- opls_get_all(
    X         = Xr,
    Y         = Yr,
    ncomp     = ncomp,
    scale     = scale,
    maxiter   = max_iter,
    tol       = tol,
    algorithm = method
  )
  
  list(
    model      = model_fit,
    cv_results = cv_results,
    transf     = model_fit$transf
  )
}


# =============================================================================
# Internal: GPR model fitting
# =============================================================================

.fit_model_gpr <- function(
    Xr, Yr, fit_method, control, verbose
) {
  
  noise_variance <- fit_method$noise_variance
  gpr_center     <- fit_method$center
  gpr_scale      <- fit_method$scale
  do_cv          <- control$validation_type == "lgo"
  
  cv_results <- NULL
  
  # Cross-validation
  if (do_cv) {
    if (verbose) message("Running cross-validation...")
    
    cv_samples <- sample_stratified(
      y           = Yr,
      p           = control$p,
      number      = control$number,
      group       = NULL,
      replacement = FALSE,
      seed        = NULL
    )
    
    cv_raw <- gaussian_process_cv(
      X        = Xr,
      Y        = Yr,
      mindices = cv_samples$hold_in,
      pindices = cv_samples$hold_out,
      noisev   = noise_variance,
      scale    = gpr_scale
    )
    
    cv_results <- data.frame(
      rmse      = mean(cv_raw$rmse_seg),
      rmse_sd   = sd(cv_raw$rmse_seg),
      st_rmse   = mean(cv_raw$st_rmse_seg),
      r2        = mean(cv_raw$rsq_seg),
      row.names = NULL
    )
  }
  
  # Fit final model
  if (verbose) message("Fitting model...")
  
  model_fit <- gaussian_process(
    X      = Xr,
    Y      = Yr,
    noisev = noise_variance,
    scale  = gpr_scale
  )
  
  transf <- list(
    Xcenter = model_fit$Xcenter,
    Xscale  = model_fit$Xscale,
    Ycenter = model_fit$Ycenter,
    Yscale  = model_fit$Yscale
  )
  
  list(
    model      = model_fit,
    cv_results = cv_results,
    transf     = transf
  )
}



# =============================================================================
# Predict method
# =============================================================================
#' @export
predict.resemble_model <- function(object, newdata, ncomp = NULL, ...) {
  
  if (!is.matrix(newdata)) {
    newdata <- as.matrix(newdata)
  }
  
  if (ncol(newdata) != object$n_vars) {
    stop(
      "'newdata' must have ", object$n_vars, " columns (same as training data).",
      call. = FALSE
    )
  }
  
  if (inherits(object$fit_method, "fit_pls")) {
    
    if (is.null(ncomp)) {
      ncomp <- object$fit_method$ncomp
      sel_ncomp <- 1:ncomp
    } else {
      sel_ncomp <- ncomp
    }
    
    if (ncomp > object$fit_method$ncomp) {
      stop(
        "'ncomp' (", ncomp, ") exceeds the number of components in the model (",
        object$fit_method$ncomp, ").",
        call. = FALSE
      )
    }
    
    yhat <- predict_opls(
      bo      = object$model$bo,
      b       = object$model$coefficients,
      ncomp   = ncomp,
      newdata = newdata,
      scale   = object$fit_method$scale,
      Xscale  = object$model$transf$Xscale
    )
    colnames(yhat) <- paste0("ncomp", seq_len(ncol(yhat)))
    yhat <- yhat[, sel_ncomp, drop = FALSE]
    

    
  } else if (inherits(object$fit_method, "fit_gpr")) {
    
    yhat <- predict_gaussian_process(
      Xz      = object$model$Xz,
      alpha   = object$model$alpha,
      newdata = newdata,
      scale   = object$model$is_scaled,
      Xcenter = object$model$Xcenter,
      Xscale  = object$model$Xscale,
      Ycenter = object$model$Ycenter,
      Yscale  = object$model$Yscale
    )
    
    yhat <- as.matrix(yhat)
    colnames(yhat) <- "predicted"
    
  } else {
    stop("Unknown model type.", call. = FALSE)
  }
  
  rownames(yhat) <- rownames(newdata)
  yhat
}

