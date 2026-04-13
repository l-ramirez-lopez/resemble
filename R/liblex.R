#' @title Build a precomputed library of localised experts using memory-based learning
#' @name liblex
#' @aliases liblex predict.liblex
#' 
#' @description
#' Constructs a library of local predictive models based on memory-based 
#' learning (MBL). For each anchor observation, a local regression model is 
#' fitted using its nearest neighbors from the reference set. This 
#' implementation is based on the methods proposed in 
#' Ramirez-Lopez et al. (2026b).
#'
#' @usage
#' liblex(Xr, Yr, neighbors,
#'        diss_method = diss_pca(ncomp = ncomp_by_opc()),
#'        fit_method = fit_wapls(min_ncomp = 3, max_ncomp = 15),
#'        anchor_indices = NULL, gh = TRUE, group = NULL,
#'        control = liblex_control(), verbose = TRUE, ...)
#'
#' \method{predict}{liblex}(object, newdata, diss_method = NULL,
#'         weighting = c("gaussian", "tricube", "triweight", "triangular",
#'                       "quartic", "parabolic", "cauchy", "none"),
#'         adaptive_bandwidth = TRUE, reliability_weighting = TRUE,
#'         range_prediction_limits = FALSE, residual_cutoff = NULL,
#'         enforce_indices = NULL, probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
#'         verbose = TRUE, allow_parallel = TRUE, blas_threads = 1L, ...)
#'         
#' \method{plot}{liblex}(x, ...)
#'
#' @param Xr A numeric matrix of predictor variables with dimensions `n × p`
#'   (observations in rows, variables in columns). Column names are required.
#' @param Yr A numeric vector or single-column matrix of length `n` containing
#'   the response variable values corresponding to `Xr`. Missing values are
#'   allowed (see Details).
#' @param neighbors A neighbor selection object specifying how to select
#'   neighbors. Use \code{\link{neighbors_k}()} for fixed k-nearest neighbors
#'   or \code{\link{neighbors_diss}()} for dissimilarity threshold-based
#'   selection.
#' @param diss_method For \code{liblex}: either a `diss_*` object specifying 
#'   the dissimilarity method, or a precomputed numeric dissimilarity matrix.
#'   If a `diss_*` object (e.g., [diss_pca()], [diss_pls()], [diss_correlation()],
#'   [diss_euclidean()]), dissimilarities are computed internally. Additional 
#'   parameters (centering, scaling, number of components) are controlled 
#'   within the object.
#'   If a matrix is provided:
#'   \itemize{
#'     \item When `anchor_indices = NULL`: must be square (`n × n`) with zeros 
#'       on the diagonal.
#'     \item When `anchor_indices` is specified (length `m`): must have 
#'       dimensions `n × m`, with `diss_method[anchor_indices, ]` having zeros 
#'       on the diagonal.
#'   }
#'   Default is `diss_pca(ncomp = ncomp_by_opc())`.
#'
#'   For \code{predict.liblex}: a dissimilarity method object created by 
#'   \code{diss_*()} constructors. If not provided, uses the method stored in
#'   \code{object}. Required if \code{object} was built with a precomputed
#'   dissimilarity matrix.
#' @param fit_method A `local_fit` object specifying the local regression 
#'   method. Currently supported: [fit_pls()] and [fit_wapls()]. 
#'   [fit_gpr()] is not yet supported. Default is 
#'   `fit_wapls(min_ncomp = 3, max_ncomp = 15)`.
#' @param anchor_indices An optional integer vector specifying row indices 
#'   of `Xr` around which to build local models. If `NULL` (default), models 
#'   are built for all observations. See Details.
#' @param gh Logical indicating whether to compute the GH distance 
#'   (Mahalanobis distance in PLS score space) for each anchor observation. 
#'   Default is `TRUE`.
#' @param group An optional factor assigning group labels to observations in 
#'   `Xr`. Used for leave-group-out validation to avoid pseudo-replication. 
#'   When an observation is selected for validation, all observations from 
#'   the same group are excluded from model fitting. Default is `NULL` 
#'   (each observation is its own group).
#' @param control A list of control parameters created by [liblex_control()].
#'   Default is `liblex_control()`.
#' @param object A fitted object of class \code{"liblex"}, created by
#'   \code{\link{liblex}} (for \code{predict}).
#' @param newdata A numeric matrix or data frame containing new predictor 
#'   values. Must include all predictors used in \code{object}.
#' @param weighting Character string specifying the kernel weighting function
#'   applied to neighbours when combining predictions. Options are:
#'   \code{"gaussian"} (default), \code{"tricube"}, \code{"triweight"},
#'   \code{"triangular"}, \code{"quartic"}, \code{"parabolic"}, \code{"cauchy"},
#'   or \code{"none"} (equal weights). See Details for kernel definitions.
#' @param adaptive_bandwidth Logical indicating whether to use adaptive 
#'   bandwidth for kernel weighting. When \code{TRUE} (default), the bandwidth 
#'   is set to the maximum dissimilarity within the neighborhood of each observation,
#'   so weights adapt to local density. When \code{FALSE}, a fixed global 
#'   bandwidth is used across all predictions, which may result in uniform 
#'   weights in sparse regions or overly concentrated weights in dense regions.
#' @param reliability_weighting Logical indicating whether to weight expert 
#'   predictions by their estimated reliability. When \code{TRUE} (default), 
#'   the contribution of each experrt is additionally weighted by the inverse of its 
#'   cross-validation residual variance, giving more influence to models that 
#'   performed well during fitting. When \code{FALSE}, only dissimilarity-based 
#'   kernel weights are used.
#' @param probs A numeric vector of probabilities in \eqn{[0, 1]} for computing
#'   weighted quantiles of expert predictions. Default is
#'   \code{c(0.05, 0.25, 0.5, 0.75, 0.95)}.
#' @param range_prediction_limits Logical. If \code{TRUE}, predictions falling
#'   outside the 5th–95th percentile range of neighbour response values are
#'   clipped to those limits. Default is \code{FALSE}.
#' @param residual_cutoff Numeric threshold for excluding models. Models with
#'   absolute residuals exceeding this value are penalized during neighbour
#'   selection. Default is \code{NULL} (no exclusion).
#' @param enforce_indices Optional integer vector specifying model indices that
#'   must always be included in each prediction neighborhood. These models are
#'   assigned the minimum dissimilarity of the neighborhood to ensure selection.
#'   Default is \code{NULL} (no enforced models).
#' @param allow_parallel Logical indicating whether parallel computation is
#'   permitted if a backend is registered. Default is \code{TRUE}.
#' @param verbose Logical indicating whether to display progress messages. 
#'   Default is `TRUE`.
#' @param blas_threads Integer specifying the number of BLAS threads to use.
#'   Default is \code{1L}. Requires the \pkg{RhpcBLASctl} package for thread
#'   control.
#' @param x An object of class `"liblex"` as returned by [liblex()].
#' @param ... Additional arguments (currently unused).
#' 
#' @details
#' By default, local models are constructed for all `n` observations in the
#' reference set. Alternatively, specify a subset of `m` observations 
#' (`m < n`) via `anchor_indices` to reduce computation.
#'
#' Each local model uses neighbors selected from the full reference set, but
#' models are only built for anchor observations. This is useful for large
#' datasets where building models for all observations is computationally
#' prohibitive.
#'
#' When dissimilarity methods depend on `Yr` (e.g., PLS-based distances), the
#' response values of anchor observations are excluded during dissimilarity
#' computation for efficiency. However, anchor response values are always
#' used when fitting local models.
#'
#' The number of anchors must not exceed 90% of `nrow(Xr)`; to build models 
#' for all observations, use `anchor_indices = NULL`.
#'
#' ## Relationship between anchors and neighborhood size
#'
#' The `neighbors` argument controls the neighborhood size (`k`) used both 
#' for fitting local models and for retrieving experts during prediction. 
#' When `anchor_indices` is specified, the number of available experts equals 
#' the number of anchors. If `max(k)` exceeds the number of anchors and 
#' tuning selects a large optimal `k`, prediction will retrieve fewer experts 
#' than specified. For reliable predictions, ensure the number of anchors is 
#' at least as large as the maximum `k` value being evaluated.
#' 
#' ## Missing values in Yr
#' 
#' Missing values in `Yr` are permitted. Observations with missing response
#' values can still serve as neighbors but are excluded from model fitting
#' as target observations.
#'
#' ## GH distance
#' 
#' The GH distance is computed independently from `diss_method` using a PLS
#' projection with optimized component selection. This provides a measure of
#' how far each observation lies from the center of the reference set in the
#' PLS score space.
#'
#' ## Validation and tuning
#' 
#' When `control$mode = "validate"` or `control$tune = TRUE`, nearest-neighbor
#' cross-validation is performed. For each anchor observation, its nearest
#' neighbor is excluded, a model is fitted on remaining neighbors, and the
#' excluded neighbor's response is predicted. This provides validation 
#' statistics for parameter selection.
#'
#' ## Prediction
#' 
#' For each observation in \code{newdata}, the \code{predict} method:
#' \enumerate{
#'   \item Computes dissimilarities to anchor observations (or their 
#'     neighbourhood centres) stored in \code{object}.
#'   \item Selects the \code{k} nearest neighbours based on the optimal 
#'     \code{k} determined during model fitting.
#'   \item Applies kernel weighting based on dissimilarity.
#'   \item Combines expert predictions using weighted averaging.
#' }
#'
#' ## Kernel weighting functions
#' 
#' The weighting functions follow Cleveland and Devlin (1988). Let \eqn{d} be 
#' the normalised dissimilarity (scaled to \eqn{[0, 1]} within the neighbourhood 
#' when \code{adaptive_bandwidth = TRUE}). The available kernels are:
#' \itemize{
#'   \item \code{"gaussian"}: \eqn{w = \exp(-d^2)}
#'   \item \code{"tricube"}: \eqn{w = (1 - d^3)^3}
#'   \item \code{"triweight"}: \eqn{w = (1 - d^2)^3}
#'   \item \code{"triangular"}: \eqn{w = 1 - d}
#'   \item \code{"quartic"}: \eqn{w = (1 - d^2)^2}
#'   \item \code{"parabolic"}: \eqn{w = 1 - d^2}
#'   \item \code{"cauchy"}: \eqn{w = 1 / (1 + d^2)}
#'   \item \code{"none"}: \eqn{w = 1} (equal weights)
#' }
#'
#' @return 
#' \strong{For \code{liblex}:} A list of class `"liblex"` (when 
#' `control$mode = "build"`) or `"liblex_validation"` (when 
#' `control$mode = "validate"`) containing:
#' \itemize{
#'   \item \code{dissimilarity}: List containing the dissimilarity method 
#'     and matrix.
#'   \item \code{fit_method}: Fit constructor from \code{fit_method}.
#'   \item \code{gh}: If `gh = TRUE`, a list with GH distances and the PLS 
#'     projection.
#'   \item \code{results}: Data frame of validation statistics for each 
#'     parameter combination (if validation was performed).
#'   \item \code{best}: The optimal parameter combination based on 
#'     `control$metric`.
#'   \item \code{optimal_params}: List with optimal `k` and `ncomp` values.
#'   \item \code{residuals}: Residuals from predictions using optimal 
#'     parameters.
#'   \item \code{coefficients}: (Build mode only) List of regression 
#'     coefficients: `B0` (intercepts), `B` (slopes).
#'   \item \code{vips}: (Build mode only) Variable importance in projection 
#'     scores.
#'   \item \code{selectivity_ratios}: (Build mode only) Selectivity ratios 
#'     for each predictor.
#'   \item \code{scaling}: (Build mode only) Centering and scaling vectors 
#'     for prediction.
#'   \item \code{neighborhood_stats}: Statistics (response quantiles) for 
#'     each neighborhood size.
#'   \item \code{anchor_indices}: The anchor indices used.
#'   \item \code{neighbors}: The object passed to \code{neighbors}.
#' }
#'
#' \strong{For \code{predict.liblex}:} A list with the following components:
#' \itemize{
#'   \item \code{predictions}: A data frame containing:
#'     \itemize{
#'       \item \code{pred}: Weighted mean predictions.
#'       \item \code{pred_sd}: Weighted standard deviation of expert 
#'         predictions.
#'       \item \code{q*}: Weighted quantiles at probabilities specified by 
#'         \code{probs}.
#'       \item \code{gh}: Global Mahalanobis distance (if computed during 
#'         fitting).
#'       \item \code{min_yr}: Minimum response value (5th percentile) across 
#'         neighbours.
#'       \item \code{max_yr}: Maximum response value (95th percentile) across 
#'         neighbours.
#'       \item \code{below_min}: Logical indicating prediction below 
#'         \code{min_yr}.
#'       \item \code{above_max}: Logical indicating prediction above 
#'         \code{max_yr}.
#'     }
#'   \item \code{neighbors}: A list with:
#'     \itemize{
#'       \item \code{indices}: Matrix of neighbour indices (models) for each 
#'         observation.
#'       \item \code{dissimilarities}: Matrix of corresponding dissimilarity 
#'         scores.
#'     }
#'   \item \code{expert_predictions}: A list with:
#'     \itemize{
#'       \item \code{weights}: Matrix of kernel weights applied to each expert.
#'       \item \code{predictions}: Matrix of raw predictions from each expert.
#'       \item \code{weighted}: Matrix of weighted predictions from each expert.
#'     }
#' }
#'
#' @references 
#' Cleveland, W. S., & Devlin, S. J. (1988). Locally weighted regression:
#' An approach to regression analysis by local fitting.
#' \emph{Journal of the American Statistical Association}, 83(403), 596–610.
#'
#' Naes, T., Isaksson, T., & Kowalski, B. (1990). Locally weighted regression
#' and scatter correction for near-infrared reflectance data.
#' \emph{Analytical Chemistry}, 62(7), 664–673.
#'
#' Ramirez-Lopez, L., Behrens, T., Schmidt, K., Stevens, A., Dematte, J.A.M.,
#' Scholten, T. (2013). The spectrum-based learner: A new local approach for
#' modeling soil vis-NIR spectra of complex datasets. Geoderma, 195-196, 
#' 268-279.
#' 
#' Ramirez-Lopez, L., Metz, M., Lesnoff, M., Orellano, C.,
#' Perez-Fernandez, E., Plans, M., Breure, T., Behrens, T.,
#' Viscarra Rossel, R., & Peng, Y. (2026b). Rethinking local spectral
#' modelling: From per-query refitting to model libraries. 
#' \emph{Analytica Chimica Acta}, under review.
#'
#' Rajalahti, T., Arneberg, R., Berven, F.S., Myhr, K.M., Ulvik, R.J., 
#' Kvalheim, O.M. (2009). Biomarker discovery in mass spectral profiles by 
#' means of selectivity ratio plot. Chemometrics and Intelligent Laboratory 
#' Systems, 95(1), 35-48.
#'
#' @seealso 
#' [liblex_control()] for control parameters, [neighbors_k()] for neighborhood
#' specification, [diss_pca()], [diss_pls()], [diss_correlation()] for 
#' dissimilarity methods, [fit_pls()], [fit_wapls()] for fitting methods.
#'
#' @author \href{https://orcid.org/0000-0002-5369-5120}{Leonardo Ramirez-Lopez}
#' 
#' @examples
#' \dontrun{
#' library(prospectr)
#' data(NIRsoil)
#'
#' # Preprocess spectra
#' NIRsoil$spc_pr <- savitzkyGolay(
#'   detrend(NIRsoil$spc, wav = as.numeric(colnames(NIRsoil$spc))),
#'   m = 1, p = 1, w = 7
#' )
#'
#' # Missing values in the response are allowed
#' train_x <- NIRsoil$spc_pr[NIRsoil$train == 1, ]
#' train_y <- NIRsoil$Ciso[NIRsoil$train == 1]
#' test_x  <- NIRsoil$spc_pr[NIRsoil$train == 0, ]
#' test_y  <- NIRsoil$Ciso[NIRsoil$train == 0]
#'
#' # Build library
#' model_library <- liblex(
#'   Xr = train_x,
#'   Yr = train_y,
#'   neighbors = neighbors_k(c(30, 40)),
#'   diss_method = diss_correlation(ws = 27, scale = TRUE),
#'   fit_method = fit_wapls(
#'     min_ncomp = 4,
#'     max_ncomp = 17,
#'     scale = FALSE,
#'     method = "mpls"
#'   ),
#'   control = liblex_control(tune = TRUE)
#' )
#'
#' # Visualise neighborhood centroids and samples to predict
#' matplot(
#'   as.numeric(colnames(model_library$scaling$local_x_center)),
#'   t(test_x),
#'   col = rgb(1, 0, 0, 0.3),
#'   lty = 1,
#'   type = "l",
#'   xlab = "Wavelength (nm)",
#'   ylab = "First derivative detrended absorbance"
#' )
#' matlines(
#'   as.numeric(colnames(model_library$scaling$local_x_center)),
#'   t(model_library$scaling$local_x_center),
#'   col = rgb(0, 0, 1, 0.3),
#'   lty = 1,
#'   type = "l"
#' )
#' grid(lty = 1)
#' legend(
#'   "topright",
#'   legend = c("Samples to predict", "Neighborhood centroids"),
#'   col = c(rgb(1, 0, 0, 0.8), rgb(0, 0, 1, 0.8)),
#'   lty = 1,
#'   lwd = 2,
#'   bty = "n"
#' )
#'
#' # Predict new observations
#' y_hat_liblex <- predict(model_library, test_x)
#'
#' # Predicted versus observed values
#' lims <- range(y_hat_liblex$predictions$pred, test_y, na.rm = TRUE)
#' plot(
#'   y_hat_liblex$predictions$pred,
#'   test_y,
#'   pch = 16,
#'   col = rgb(0, 0, 0, 0.5),
#'   xlab = "Predicted",
#'   ylab = "Observed",
#'   xlim = lims,
#'   ylim = lims
#' )
#' abline(a = 0, b = 1, col = "red")
#' grid(lty = 1)
#' 
#' ## run liblex in parallel (requires a parallel backend, e.g., doParallel)
#' library(doParallel)
#' n_cores <- min(2, parallel::detectCores())
#' clust <- makeCluster(n_cores)
#' registerDoParallel(clust)
#' 
#' model_library2 <- liblex(
#'   Xr = train_x,
#'   Yr = train_y,
#'   neighbors = neighbors_k(c(30, 40)),
#'   fit_method = fit_wapls(min_ncomp = 4, max_ncomp = 17, method = "simpls")
#' )
#' 
#' y_hat_liblex2 <- predict(model_library2, test_x)
#' registerDoSEQ()
#' try(stopCluster(clust))
#' }
#' @export
######################################################################
# resemble
# Copyright (C) 2026 Leonardo Ramirez-Lopez and Antoine Stevens
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
######################################################################


liblex <- function(
    Xr, Yr, neighbors,
    diss_method = diss_pca(ncomp = ncomp_by_opc()),
    fit_method = fit_wapls(min_ncomp = 3, max_ncomp = 15),
    anchor_indices = NULL,
    gh = TRUE, 
    group = NULL,
    control = liblex_control(),
    verbose = TRUE,
    ...
) {
  f_call <- match.call()
  
  # --- Xr validation ---
  if (!is.matrix(Xr)) {
    Xr <- as.matrix(Xr)
  }
  if (is.null(colnames(Xr))) {
    stop("'Xr' must have column names", call. = FALSE)
  }
  if (!is.numeric(Xr)) {
    stop("'Xr' must be numeric", call. = FALSE)
  }
  
  # --- Yr validation ---
  
  if (!is.null(ncol(Yr))) {
    if (ncol(Yr) != 1L) {
      stop("'Yr' must be a numeric vector or single-column matrix (only one response variable allowed at the moment)", call. = FALSE)
    }
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("'verbose' must be a single TRUE or FALSE value", call. = FALSE)
  }
  if (length(Yr) != nrow(Xr)) {
    stop("'Yr' must have the same number of observations as 'Xr'", call. = FALSE)
  }
  if (!is.matrix(Yr)) {
    Yr <- matrix(Yr, nrow = nrow(Xr))
  }
  if (!is.numeric(Yr) || ncol(Yr) != 1L) {
    stop("'Yr' must be a numeric vector or single-column matrix", call. = FALSE)
  }
  if (nrow(Yr) != nrow(Xr)) {
    stop("'Yr' must have the same number of observations as 'Xr'", call. = FALSE)
  }
  
  # --- neighbors validation ---
  if (missing(neighbors)) {
    stop("'neighbors' is required", call. = FALSE)
  }
  if (!(inherits(neighbors, "neighbors_k") || inherits(neighbors, "neighbors_diss"))) {
    stop("'neighbors' must be a 'neighbors_k' object", call. = FALSE)
  }

  # --- diss_method validation ---
  if (is.matrix(diss_method)) {
    if (!is.numeric(diss_method)) {
      stop("'diss_method' matrix must be numeric", call. = FALSE)
    }
    if (nrow(diss_method) != ncol(diss_method)) {
      stop("'diss_method' matrix must be square", call. = FALSE)
    }
    if (nrow(diss_method) != nrow(Xr)) {
      stop("'diss_method' matrix dimensions must match nrow(Xr)", call. = FALSE)
    }
    precomputed_diss <- TRUE
  } else if (inherits(diss_method, "diss_method")) {
    precomputed_diss <- FALSE
  } else {
    stop("'diss_method' must be a diss_* object or a dissimilarity matrix",
         call. = FALSE)
  }
  
  # --- fit_method validation ---
  if (!inherits(fit_method, "fit_method")) {
    stop("'fit_method' must be a fit_* object (e.g., fit_wapls(), fit_pls())",
         call. = FALSE)
  }
  
  # --- anchor_indices validation ---
  if (is.null(anchor_indices)) {
    if (verbose) {
      message(
        "No 'anchor_indices' provided; building local models for all observations in 'Xr'."
      )
    }
  } else {
    if (!is.numeric(anchor_indices) || any(is.na(anchor_indices))) {
      stop("'anchor_indices' must be a numeric vector without NA", call. = FALSE)
    }
    anchor_indices <- as.integer(anchor_indices)
    if (any(anchor_indices < 1L) || any(anchor_indices > nrow(Xr))) {
      stop("'anchor_indices' must be between 1 and nrow(Xr)", call. = FALSE)
    }
    if (length(anchor_indices) > floor(nrow(Xr) * 0.90)) {
      stop("'anchor_indices' exceeds 90% of nrow(Xr); ",
           "set anchor_indices = NULL to build models for all observations",
           call. = FALSE)
    }
    if (any(duplicated(anchor_indices))) {
      stop("'anchor_indices' contains duplicate values", call. = FALSE)
    }
    maximum_k <- if (inherits(neighbors, "neighbors_k")) {
      max(neighbors$k)
    } else {
      neighbors$k_max
    }
    
    if (maximum_k > length(anchor_indices)) {
      warning(
        "Max. number of neighbors (", maximum_k, ") exceeds number of anchors (", 
        length(anchor_indices), "); prediction may retrieve fewer experts than ",
        "optimal if tuning selects a large k.",
        call. = FALSE
      )
    }
    
  }
  
  # --- gh validation ---
  if (!is.logical(gh) || length(gh) != 1L || is.na(gh)) {
    stop("'gh' must be TRUE or FALSE", call. = FALSE)
  }
  
  # --- control validation ---
  if (!inherits(control, "liblex_control")) {
    stop("'control' must be created by liblex_control()", call. = FALSE)
  }
  
  # --- chunk_size validation (from control) ---
  min_size <- if (!is.null(anchor_indices)) {
    min(nrow(Xr), length(anchor_indices))
  } else {
    nrow(Xr)
  }
  if (control$chunk_size > min_size) {
    stop("'chunk_size' in control cannot exceed ", min_size, call. = FALSE)
  }
  
  # --- verbose validation ---
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("'verbose' must be TRUE or FALSE", call. = FALSE)
  }
  

  if (!precomputed_diss && inherits(diss_method, c("diss_pca", "diss_pls"))) {
    diss_method$ncomp <- cap_ncomp_to_data(
      diss_method$ncomp,
      n_obs = nrow(Xr),
      n_vars = ncol(Xr),
      verbose = verbose
    )
  }
  
  # --- BLAS thread management ---
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    old_blas_threads <- RhpcBLASctl::blas_get_num_procs()
    if (old_blas_threads != control$blas_threads) {
      RhpcBLASctl::blas_set_num_threads(control$blas_threads)
      on.exit(RhpcBLASctl::blas_set_num_threads(old_blas_threads), add = TRUE)
    }
  }
  
  # ---------------------------------------------------------------------------
  # Parallel setup
  # ---------------------------------------------------------------------------
  "%mydo%" <- get("%do%")
  if (control$allow_parallel && getDoParRegistered()) {
    "%mydo%" <- get("%dopar%")
  }
  
  if (is.matrix(diss_method)) {
    diss_method_type <- "Precomputed dissimilarity matrix"
    
    # --- Validate matrix ---
    if (!is.numeric(diss_method)) {
      stop("'diss_method' matrix must be numeric", call. = FALSE)
    }
    if (nrow(diss_method) != nrow(Xr)) {
      stop("'diss_method' must have nrow equal to nrow(Xr)", call. = FALSE)
    }

    # --- Validate dimensions based on anchor_indices ---
    if (!is.null(anchor_indices)) {
      if (ncol(diss_method) != length(anchor_indices)) {
        stop("'diss_method' must have ncol equal to length(anchor_indices)",
             call. = FALSE)
      }
      if (any(abs(diag(diss_method[anchor_indices, ])) > 1e-8)) {
        stop("'diss_method[anchor_indices, ]' must have zeros on the diagonal",
             call. = FALSE)
      }
      Yr_anchor <- Yr[anchor_indices, , drop = FALSE]
    } else {
      if (ncol(diss_method) != nrow(Xr)) {
        stop("'diss_method' must be square with dimensions equal to nrow(Xr)",
             call. = FALSE)
      }
      if (any(abs(diag(diss_method)) > 1e-8)) {
        stop("'diss_method' diagonal must contain only zeros", call. = FALSE)
      }
      Yr_anchor <- Yr
    }

    # --- Store and free original matrix ---
    dsm <- list(
      diss_method = diss_method_type,
      dissimilarity = diss_method
    )
    rm(diss_method)
    gc()
  } else if (inherits(diss_method, "diss_method")) {
    precomputed_diss <- FALSE
    diss_method_type <- class(diss_method)[1L]

    if (verbose) {
      cat("Computing dissimilarities... \n")
    }

    if (!is.null(anchor_indices)) {
      Yr_anchor <- Yr[anchor_indices, , drop = FALSE]
      if (sum(!is.na(Yr_anchor)) < 3L) {
        stop("At least 3 non-missing values in 'Yr' required for anchor indices",
             call. = FALSE)
      }
      
      # Create copy to avoid modifying user input
      diss_method_internal <- diss_method
      
      # Force projection if needed for GH computation
      if (gh || inherits(diss_method_internal, c("diss_pls")) || inherits(diss_method_internal, c("diss_pca"))) {
        diss_method_internal$return_projection <- TRUE
      }
      
      dsm <- dissimilarity(
        Xr = Xr[-anchor_indices, , drop = FALSE],
        Xu = Xr[anchor_indices, , drop = FALSE],
        diss_method = diss_method_internal,
        Yr = Yr[-anchor_indices, , drop = FALSE]
      )
      
      # Expand dissimilarity matrix to full Xr dimensions
      # (rows for anchor_indices were not computed, fill with NA)
      diss_full <- matrix(NA_real_, nrow(Xr), ncol(dsm$dissimilarity))
      diss_full[-anchor_indices, ] <- dsm$dissimilarity
      dsm$dissimilarity <- diss_full
      rm(diss_full)
      gc()

      # Compute dissimilarities among anchor observations
      if (inherits(diss_method, c("diss_pca", "diss_pls"))) {
        # For projection-based methods, use scores from the existing projection
        n_non_anchor <- nrow(Xr) - length(anchor_indices)
        z_anchor <- dsm$projection$scores[-(1:n_non_anchor), , drop = FALSE]
        z_anchor <- scale(
          z_anchor,
          center = FALSE,
          scale = dsm$projection$scores_sd
        )
        
        # Self-dissimilarity among anchor observations
        dsm$dissimilarity[anchor_indices, ] <- fast_self_euclid(z_anchor)
        rm(z_anchor)
        gc()

        # Reorder projection scores to match original Xr order
        real_order <- order(c(seq_len(nrow(Xr))[-anchor_indices], anchor_indices))
        dsm$projection$scores <- dsm$projection$scores[real_order, , drop = FALSE]
        rownames(dsm$projection$scores) <- paste0("Xr_", seq_len(nrow(Xr)))
        
      } else {
        # For non-projection methods: pre-scale using full Xr, then compute
        # dissimilarity with center/scale disabled
        
        X_anchor <- Xr[anchor_indices, , drop = FALSE]
        
        # Pre-scale using full Xr statistics
        do_center <- isTRUE(diss_method$center)
        do_scale <- isTRUE(diss_method$scale)
        
        if (do_center || do_scale) {
          X_anchor <- scale(
            X_anchor,
            center = if (do_center) get_column_means(Xr) else FALSE,
            scale = if (do_scale) get_column_sds(Xr) else FALSE
          )
        }
        
        # Create modified diss_method with center/scale disabled
        diss_method_anchor <- diss_method
        diss_method_anchor$center <- FALSE
        diss_method_anchor$scale <- FALSE
        
        dsm_anchor <- dissimilarity(
          Xr = X_anchor,
          Xu = NULL,
          diss_method = diss_method_anchor,
          Yr = Yr_anchor
        )
        
        dsm$dissimilarity[anchor_indices, ] <- dsm_anchor$dissimilarity
        rm(dsm_anchor, X_anchor, diss_method_anchor)
        gc()
      }
    } else {
      Yr_anchor <- Yr
      # Create copy to avoid modifying user input
      diss_method_internal <- diss_method
      
      # Force projection for projection-based methods
      if (inherits(diss_method_internal, c("diss_pca", "diss_pls"))) {
        diss_method_internal$return_projection <- TRUE
      }

      dsm <- dissimilarity(
        Xr = Xr,
        Xu = NULL,
        diss_method = diss_method_internal,
        Yr = Yr
      )
    }
    
    sml <- list(diss_method = diss_method)
    dsm <- append(sml, dsm)
  } else {
    stop(
      "'diss_method' must be created by a diss_*() function or be a ",
      "numeric dissimilarity matrix", call. = FALSE
    )
  }
  # --- Compute GH distances (independent of diss_method) ---
  if (gh) {
    # Compute GH distances independently via PLS projection
    # FIXME: this is not ideal, as it requires an additional projection step.
    # in addition the optimisation of components is hard coded
    gh_projection <- ortho_projection(
      Xr = Xr,
      Xu = NULL,
      Yr = Yr,
      ncomp = ncomp_by_opc(min(nrow(Xr), 40L)),
      method = "pls",
      scale = fit_method$scale
    )

    gh_center <- colMeans(gh_projection$scores)
    
    # Use all scores when anchor_indices is NULL
    scores_for_gh <- if (is.null(anchor_indices)) {
      gh_projection$scores
    } else {
      gh_projection$scores[anchor_indices, , drop = FALSE]
    }
    
    gh_all <- as.vector(f_diss(
      scores_for_gh,
      Xu = t(gh_center),
      diss_method = "mahalanobis",
      center = FALSE,
      scale = FALSE
    ))
    
    dsm$gh <- list(
      gh_Xr = gh_all,
      projection = gh_projection
    )
  } else {
    dsm$gh <- NULL
  }

  # Determine maximum neighborhood size
  if (inherits(neighbors, "neighbors_k")) {
    max_k <- max(neighbors$k)
  } else {
    # neighbors_diss
    max_k <- neighbors$k_max
  }

  # Validate against non-missing Yr values
  n_valid_yr <- sum(!is.na(Yr))
  if (n_valid_yr < max_k) {
    stop(
      "Maximum neighborhood size (", max_k, ") exceeds the number of ",
      "non-missing observations in 'Yr' (", n_valid_yr, ")",
      call. = FALSE
    )
  }
  
  if (inherits(neighbors, "neighbors_k")) {
    # kidxmat <- diss_to_neighbors(
    #   diss_matrix = dsm$dissimilarity,
    #   k = max_k,
    #   k_diss = NULL,
    #   k_range = NULL,
    #   spike = -which(is.na(Yr)),
    #   return_dissimilarity = FALSE,
    #   skip_first = FALSE,
    #   keep_self = TRUE
    # )
    k <- neighbors$k
    len_neighbors <- length(k)
    names_nnstats <- paste0(
      "neighbor ", 
      neighbors$method, 
      " index ", 
      seq_along(neighbors$k), 
      " (k = ", neighbors$k, ")"
    )
  }
  
  if (inherits(neighbors, "neighbors_diss")) {
    # kidxmat <- diss_to_neighbors(
    #   diss_matrix = dsm$dissimilarity,
    #   k = NULL,
    #   k_diss = neighbors$threshold,
    #   k_range = c(neighbors$k_min, neighbors$k_max),
    #   spike = -which(is.na(Yr)),
    #   return_dissimilarity = FALSE,
    #   skip_first = FALSE,
    #   keep_self = TRUE
    # )
    len_neighbors <- length(neighbors$threshold)
    
    format_threshold <- function(x, digits = 4) {
      natural <- format(x, scientific = FALSE, trim = TRUE, digits = 8)
      rounded <- formatC(x, format = "f", digits = digits)
      
      decimal_part <- sub("^[^.]*\\.?", "", natural)
      decimal_places <- nchar(decimal_part)
      
      if (decimal_places > digits) {
        paste0(rounded, "...")
      } else {
        natural
      }
    }
    
    names_nnstats <- paste0(
      "neighbor ", 
      neighbors$method, 
      " index ", 
      seq_along(neighbors$threshold), 
      " (diss < ", sapply(neighbors$threshold, format_threshold), ")"
    )
  }

  
  # 
  # 
  # kidxmat2 <- diss_to_neighbors(
  #   diss_matrix = dsm$dissimilarity,
  #   k = NULL, 
  #   k_diss = 0.45, 
  #   k_range = c(5, 50),
  #   spike = -which(is.na(Yr)),
  #   return_dissimilarity = FALSE,
  #   skip_first = FALSE, 
  #   keep_self = TRUE
  # )
  # 
  # sum(abs((kidxmat2$neighbors) - (kidxmat)))
  # 

  # Find indices of the nearest neighbors
  kidxmat <- list(
    neighbors = top_k_order(
      dsm$dissimilarity,
      k = max_k,
      skip = which(is.na(Yr))
    )
  )
  
  # Extract dissimilarities for the selected neighbors
  kidxmat$neighbors_diss <- extract_by_index(dsm$dissimilarity, kidxmat$neighbors)
  
  # kidxmat$neighbors_diss
  # kidxmat$neighbors
  
  # --- group validation ---  
  # For each anchor sample, identify neighbors not in the same group
  # (used for leave-group-out validation)
  if (is.null(group)) {
    group <- as.factor(seq_len(nrow(Xr)))
    kidxgrop <- matrix(
      TRUE, 
      nrow = nrow(kidxmat$neighbors), 
      ncol = ncol(kidxmat$neighbors)
    )
    kidxgrop[1L, ] <- FALSE
  } else {
    if (length(group) != nrow(Xr)) {
      stop("'group' must have length equal to nrow(Xr)", call. = FALSE)
    }
    group <- as.factor(group)
    kidxgrop <- not_in_same_group(kidxmat$neighbors, group = group)
  }
  
  # Compute nearest neighbor statistics (quantiles of Yr in neighborhoods)
  probs <- c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1.00)
  dnms <- list(
    seq_len(ncol(kidxgrop)),
    paste0(probs * 100, "%")
  )

  nnstats <- vector("list", len_neighbors)
  #######################################################################
  #######################################################################
  #######################################################################
  #######################################################################
  # FIXME: THE compute_nn_quantiles() FUNCTION NEEDS TO BE REFACTORED
  # IN A WAY THAT IT CAN TAKE THE DISSIMILARITY MATRIX AND TO EACH COLUMN IN 
  # kidxmat$neighbors, IT APPLIES THE THRESHOLD IN neighbors$threshold
  # AND THEN COMPUTES THE QUANTILES 
  # OR EVEN BETTER IT TAKES THE FOLLOWIN ks OBJECT:
  
  if (inherits(neighbors, "neighbors_diss")) {
    # The list of number of neighbors per each threshold dissimilarity
    # nrow(ks) is equivalent to the number of thresholds
    # ncol(ks) is equivalent to the observations
    ks <- apply(
      kidxmat$neighbors_diss, 
      MARGIN = 2, 
      FUN = function(x, thresholds) {
        ith_ks <- sapply(thresholds, FUN = function(d, thresholds) sum(x < thresholds, na.rm = TRUE), d = x)
        ith_ks
      }, 
      thresholds = neighbors$threshold
    )
    ks[ks < neighbors$k_min] <- neighbors$k_min
    ks[ks > neighbors$k_max] <- neighbors$k_max
    thresholds_ks <- neighbors$threshold
  } else {
    ks <- matrix(k, nrow = length(k), ncol = nrow(Xr)) 
    thresholds_ks <- neighbors$k
  }

  for (ii in seq_len(len_neighbors)) {
    nnstats[[ii]] <- compute_nn_quantiles(
      kidxmat = kidxmat$neighbors,
      kidxgrop = kidxgrop,
      Yr = Yr,
      k = ks[ii, ],
      probs = probs
    )
    dimnames(nnstats[[ii]]) <- dnms
  }
  names(nnstats) <- names_nnstats

  #######################################################################
  #######################################################################
  #######################################################################
  
  
  # RMSE/(q"75%" - q"25%")
  # Similar to the ratio of performance to inter-quartile distance (RPIQ)
  # Bellon-Maurel, V., Fernandez-Ahumada, E., Palagos, B., Roger, J.M., 
  # McBratney, A., 2010. Critical review of chemometric indicators commonly 
  # used for assessing the quality of the prediction of soil attributes by 
  # NIR spectroscopy. TrAC Trends in Analytical Chemistry, 29(9), pp.1073-1081.
  itq <- vapply(nnstats, function(x) {
    x[, "75%"] - x[, "25%"]
  }, numeric(nrow(nnstats[[1L]])))
  
  # Dissimilarity as predictors setup
  # (Note: diss_usage not yet implemented in liblex, placeholder for future)
  dssm <- NULL
  npredictors <- ncol(Xr)

  # Progress bar setup
  addit <- if (control$mode == "build") 1L else 0L
  if (verbose) {
    pb <- txtProgressBar(min = 0, max = len_neighbors + addit, char = "-")
  }
  
  # PLS component range from fit_method
  if (inherits(fit_method, "fit_wapls")) {
    min_ncomp <- fit_method$min_ncomp
    max_ncomp <- fit_method$max_ncomp
  } else if (inherits(fit_method, "fit_pls")) {
    min_ncomp <- fit_method$ncomp
    max_ncomp <- fit_method$ncomp
  } else {
    stop(
      "'fit_method' must be fit_pls() or fit_wapls(); fit_gpr() is not yet supported",
      call. = FALSE
    )
  }
  
  if (control$tune) {
    # Grid of possible min/max PLS component combinations for optimization
    sgrid <- expand.grid(minpls = min_ncomp:max_ncomp, maxpls = min_ncomp:max_ncomp)
    sgrid <- sgrid[sgrid$minpls <= sgrid$maxpls, ]
    rownames(sgrid) <- seq_len(nrow(sgrid))
    
    # Binary grid to switch components on/off for each combination
    emgrid <- t(
      sapply(
        seq_len(nrow(sgrid)),
        FUN = function(q, x, wv) {
          wv[x[q, 1L]:x[q, 2L]] <- 1L
          wv
        },
        x = sgrid,
        wv = rep(0L, max_ncomp)
      )
    )
    emgrid <- emgrid[, min_ncomp:max_ncomp, drop = FALSE]
  } else {
    emgrid <- matrix(1L, nrow = 1L, ncol = length(min_ncomp:max_ncomp))
    sgrid <- data.frame(
      minpls = min_ncomp,
      maxpls = max_ncomp
    )
  }
  
  if (verbose) cat("Fitting models...\n")
  ## perform the nearest neighbor predictions  
  if (control$mode == "validate" || control$tune) {
    nnpreds <- sapply(
      seq_len(len_neighbors),
      FUN = .get_all_fits,
      Xr = Xr,
      Yr = Yr,
      k = ks,
      ncomp_min = min_ncomp,
      ncomp_max = max_ncomp,
      emgrid = emgrid,
      scale = fit_method$scale,
      max_iter = fit_method$max_iter,
      tol = fit_method$tol,
      algorithm = fit_method$method,
      kidxmat = kidxmat$neighbors,
      kidxgrop = kidxgrop,
      dissimilarity_mat = dssm,
      pb = if (verbose) pb else NULL,
      chunk_size = control$chunk_size
    )
    # Organize the results (in nnpreds)
    # Compute prediction statistics: R², RMSE, bias, and RPIQ-like metric
    pparam <- matrix(NA_real_, nrow(emgrid), 4L)
    
    sstats <- function(y, yhat, itqk, pparam) {
      # Prediction errors (observed - predicted)
      me <- sweep(-yhat, MARGIN = 1L, STATS = y, FUN = "+", check.margin = FALSE)
      
      # R² (squared correlation)
      pparam[, 1L] <- cor(y, yhat, use = "complete.obs")^2
      
      # RMSE (root mean squared error)
      pparam[, 2L] <- colMeans(me^2, na.rm = TRUE)^0.5
      
      # Bias (mean error)
      pparam[, 3L] <- colMeans(me, na.rm = TRUE)
      
      
      # RMSE / IQR: similar to RPIQ (ratio of performance to inter-quartile distance)
      # See Bellon-Maurel et al. (2010)
      pparam[, 4L] <- colMeans(
        sweep(me,
              MARGIN = 1L,
              STATS = itqk,
              FUN = "/",
              check.margin = FALSE)^2,
        na.rm = TRUE
      )^0.5
      
      pparam
    }
    
    # Iterator function to subset columns from multiple matrices simultaneously
    # Used for parallel processing of prediction statistics
    isubset_cols <- function(...) {
      sargs <- names(match.call())[-1L]
      nextEl <- function(..ii..) {
        sapply(sargs, FUN = function(x) get(x)[, ..ii..], simplify = FALSE)
      }
      list(nextElem = nextEl)
    }
    
    # Create iterator for nnpreds and itq matrices
    # Iterator call
    itr <- isubset_cols(nnpreds = nnpreds, itq = itq)
    
    # Then in kpredstats:
    kpredstats <- function(..k.., itr, pparam, y) {
      ne <- itr$nextElem(..k..)
      
      sstats(
        y = y,
        yhat = t(matrix(ne$nnpreds, nrow(pparam))),
        itqk = ne$itq,
        pparam = pparam
      )
    }
    
    # Compute prediction statistics for each neighborhood size
    predperformance <- lapply(
      seq_len(ncol(nnpreds)),
      FUN = kpredstats,
      itr = itr,
      pparam = pparam,
      y = Yr_anchor
    )

    # Combine prediction statistics into a data.frame
    predperformance <- data.frame(do.call("rbind", predperformance))
    colnames(predperformance) <- c("r2", "rmse", "me", "st_rmse")
    
    # Add PLS component range and neighborhood size columns
    
    if (inherits(neighbors, "neighbors_k")) {
      predperformance <- data.frame(
        min_ncomp = rep(sgrid$minpls, times = length(thresholds_ks)),
        max_ncomp = rep(sgrid$maxpls, times = length(thresholds_ks)),
        k = rep(thresholds_ks, each = nrow(pparam)),
        predperformance
      )
      ks_param_name <- "k"
    } else {
      predperformance <- data.frame(
        min_ncomp = rep(sgrid$minpls, times = length(thresholds_ks)),
        max_ncomp = rep(sgrid$maxpls, times = length(thresholds_ks)),
        diss_threshold = rep(thresholds_ks, each = nrow(pparam)),
        k_min = rep(neighbors$k_min, each = nrow(pparam)),
        k_max = rep(neighbors$k_max, each = nrow(pparam)),
        predperformance
      )
      ks_param_name <- "diss_threshold"
    }
      
    
    # Find optimal parameters based on selected metric
    if (control$metric == "rmse") {
      bestp <- predperformance[which.min(predperformance[, "rmse"]), , drop = FALSE][1L, ]
    } else {
      # metric == "r2"
      bestp <- predperformance[which.max(predperformance[, "r2"]), , drop = FALSE][1L, ]
    }
    # Extract optimal parameters
    optimal_param <- bestp[[ks_param_name]]
    optimal_min_ncomp <- bestp$min_ncomp
    optimal_max_ncomp <- bestp$max_ncomp
    
    
    # Extract the vector of predictions corresponding to the best predictions
    # (optimal k, optimal pls range)
    pls_item_idx <- which(
      sgrid$minpls == optimal_min_ncomp & sgrid$maxpls == optimal_max_ncomp
    )
    pls_item_idx <- seq(pls_item_idx, by = nrow(pparam), length.out = length(Yr_anchor))
    
    best_preds <- itr$nextElem(which(thresholds_ks == optimal_param))$nnpreds[pls_item_idx]
    best_preds_residuals <- Yr_anchor - best_preds
    
  } else {
    # No validation: use maximum k and full PLS range
    optimal_param <- max(thresholds_ks)
    optimal_min_ncomp <- min_ncomp
    optimal_max_ncomp <- max_ncomp
    best_preds_residuals <- NULL
    bestp <- NULL
    predperformance <- NULL
    
    if (verbose) {
      setTxtProgressBar(pb, 1L)
    }
  }
  
  if (control$mode == "build") {
    
    # final vector of optimal ks
    optimal_ks <- ks[which(thresholds_ks == optimal_param), ]
    
    # Calculate number of variables for library storage
    n_var <- 1L + (5L * ncol(Xr))
    
    # Template matrix for storing PLS library coefficients
    plslib_template <- matrix(
      NA_real_,
      nrow = n_var,
      ncol = control$chunk_size
    )
    
    # Number of iterations for parallel processing
    n_iter <- ceiling(ncol(dsm$dissimilarity) / control$chunk_size)

    # Build PLS library in parallel (or sequentially if no backend registered)
    # Declare foreach variables to avoid R CMD check NOTE
    ksubsets <- ithbarrio <- iset <- NULL
    plslib <- foreach(
      i = seq_len(n_iter),
      .export = c(
        "ith_pred_subsets",
        "ith_subsets_list",
        "ith_subsets_by_group",
        ".get_all_fits",
        "ith_local_fit",
        "final_fits_cpp"
      ),
      ithbarrio = ith_subsets_list(
        x = Xr,
        y = Yr,
        kindx = kidxmat$neighbors[seq_len(max(optimal_ks)), , drop = FALSE],
        neighborhood_sizes = optimal_ks,
        D = dssm,
        chunk_size = control$chunk_size
      )
    ) %mydo% {
      iplslib <- plslib_template
      for (j in seq_along(ithbarrio)) {
        ij_pls <- final_fits_cpp(
          X = ithbarrio[[j]]$x,
          Y = ithbarrio[[j]]$y,
          new_x = ithbarrio[[j]]$x[1L, , drop = FALSE],
          ncomp_min = optimal_min_ncomp,
          ncomp_max = optimal_max_ncomp,
          scale = fit_method$scale,
          maxiter = fit_method$max_iter,
          tol = fit_method$tol, 
          algorithm = fit_method$method
        )
        iplslib[, j] <- unlist(ij_pls, recursive = FALSE, use.names = FALSE)
      }
      
      # Trim unused columns if last chunk is smaller
      if (j < ncol(iplslib)) {
        iplslib <- iplslib[, seq_len(j), drop = FALSE]
      }
      iplslib
    }

    # Combine parallel results
    plslib <- do.call("cbind", plslib)
    
    if (verbose) {
      setTxtProgressBar(pb, length(thresholds_ks) + 1L)
    }
    
    plslib <- t(plslib)
    
    # Setup column names based on whether dissimilarity is used as predictors
    namesk <- NULL
    npredictors <- ncol(Xr)

    # Extract scaling vectors from library matrix
    xscale <- plslib[, -seq_len((4L * npredictors) + 1L), drop = FALSE]
    plslib <- plslib[, seq_len((4L * npredictors) + 1L), drop = FALSE]
    
    xcenter <- plslib[, -seq_len((3L * npredictors) + 1L), drop = FALSE]
    plslib <- plslib[, seq_len((3L * npredictors) + 1L), drop = FALSE]
    
    # Extract VIPs and selectivity ratios
    plsvips <- plslib[, -c(seq_len(npredictors + 1L), 
                           (ncol(plslib) - npredictors + 1L):ncol(plslib)), 
                      drop = FALSE]
    plssratios <- plslib[, (ncol(plslib) - npredictors + 1L):ncol(plslib), drop = FALSE]
    plslib <- plslib[, seq_len(npredictors + 1L), drop = FALSE]
    
    # Assign column names
    colnames(plslib) <- c("b0", namesk, colnames(Xr))
    colnames(xcenter) <- colnames(xscale) <- c(namesk, colnames(Xr))
    colnames(plssratios) <- colnames(plsvips) <- c(namesk, colnames(Xr))
    
    # Extract regression coefficients
    bs <- list(
      B0 = plslib[, 1L],
      B = plslib[, colnames(plslib) %in% colnames(Xr), drop = FALSE]
    )
    
    
    # Global centering/scaling vectors for prediction
    # (used when applying the library to new observations)
    if (precomputed_diss) {
      global_center <- rep(0, ncol(Xr))
      global_scale <- rep(1, ncol(Xr))
    } else {
      global_center <- if (isTRUE(diss_method$center)) get_column_means(Xr) else rep(0, ncol(Xr))
      global_scale <- if (isTRUE(diss_method$scale)) get_column_sds(Xr) else rep(1, ncol(Xr))
    }
    
    # Scaling information for prediction
    xcenter <- xcenter[, colnames(xcenter) %in% colnames(Xr), drop = FALSE]
    xscale <- xscale[, colnames(xscale) %in% colnames(Xr), drop = FALSE]
    iscale <- list(
      center = global_center,
      scale = global_scale,
      # Adjust local centers by local scales for prediction
      # our pls first scales and then centers, so we need to adjust the 
      # local centers accordingly
      local_x_center = xcenter * xscale, 
      local_x_scale = xscale
    )
    
    param_name <- if (inherits(neighbors, "neighbors_k")) "k" else "diss_threshold"
    optimal_params <- c(
      setNames(list(optimal_param), param_name),
      list(ncomp = c(min = optimal_min_ncomp, max = optimal_max_ncomp))
    )
    
    fresults <- list(
      dissimilarity = dsm[!names(dsm) %in% "gh"],
      fit_method = fit_method,
      gh = dsm$gh,
      results = predperformance,
      best = bestp,
      optimal_params = optimal_params,
      
      
      residuals = best_preds_residuals,
      coefficients = bs,
      vips = plsvips,
      selectivity_ratios = plssratios,
      scaling = iscale,
      neighborhood_stats = nnstats
    )
    
    # Row names for output
    if (!is.null(rownames(Xr))) {
      row_names <- rownames(Xr)
    } else {
      row_names <- seq_len(nrow(Xr))
    }
    
    if (!is.null(anchor_indices)) {
      row_names <- row_names[anchor_indices]
    }
    
    # Assign row names to nearest neighbor statistics
    fresults$neighborhood_stats <- lapply(
      fresults$neighborhood_stats,
      FUN = function(x, nms) {
        rownames(x) <- nms
        x
      },
      nms = row_names
    )
    # Assign row names to library components
    rownames(fresults$vips) <-
      rownames(fresults$selectivity_ratios) <-
      rownames(fresults$coefficients$B) <-
      names(fresults$coefficients$B0) <- row_names
    
    if (!is.null(fresults$coefficients$Bk)) {
      rownames(fresults$coefficients$Bk) <- row_names
    }
  } else {
    
    # Validation only (no library built)
    fresults <- list(
      dissimmilarity = dsm[!names(dsm) %in% "gh"],
      fit_method = fit_method,
      gh = dsm$gh,
      results = predperformance,
      best = bestp,
      residuals = best_preds_residuals,
      neighborhood_stats = nnstats
    )
    
    # Row names for output
    if (!is.null(rownames(Xr))) {
      row_names <- rownames(Xr)
    } else {
      row_names <- seq_len(nrow(Xr))
    }
    
    if (!is.null(anchor_indices)) {
      row_names <- row_names[anchor_indices]
    }
    
    # Assign row names to nearest neighbor statistics
    fresults$neighborhood_stats <- lapply(
      fresults$neighborhood_stats,
      FUN = function(x, nms) {
        rownames(x) <- nms
        x
      },
      nms = row_names
    )
  }
  
  if (verbose) close(pb)
  
  fresults$anchor_indices <- if (!is.null(anchor_indices)) anchor_indices else seq_len(nrow(Xr))
  
  if (control$mode == "validate" || control$tune) {
    names(fresults$residuals) <- row_names
  }
  
  # add the original object passed to neighbors
  fresults$neighbors <- neighbors
  
  attr(fresults, "call") <- f_call
  class(fresults) <- c("liblex", "list")
  fresults
}



#' @aliases liblex
#' @export
predict.liblex <- function(
    object,
    newdata,
    diss_method = NULL,
    weighting = c(
      "gaussian", "tricube", "triweight", "triangular",
      "quartic", "parabolic", "cauchy", "none"
    ),
    adaptive_bandwidth = TRUE,
    reliability_weighting = TRUE,
    range_prediction_limits = FALSE,
    residual_cutoff = NULL,
    enforce_indices = NULL,
    probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
    verbose = TRUE,
    allow_parallel = TRUE,
    blas_threads = 1L,
    ...
) {
  
  # --- Validate allow_parallel ---
  if (!is.logical(allow_parallel) || length(allow_parallel) != 1L) {
    stop("'allow_parallel' must be TRUE or FALSE", call. = FALSE)
  }
  
  # ---------------------------------------------------------------------------
  # Parallel setup
  # ---------------------------------------------------------------------------
  "%mydo%" <- get("%do%")
  if (allow_parallel && getDoParRegistered()) {
    "%mydo%" <- get("%dopar%")
  }
  
  # --- Validate blas_threads ---
  if (!is.numeric(blas_threads) || length(blas_threads) != 1L || 
      blas_threads < 1L || blas_threads != as.integer(blas_threads)) {
    stop("'blas_threads' must be a positive integer", call. = FALSE)
  }
  blas_threads <- as.integer(blas_threads)

  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    old_blas_threads <- blas_get_num_procs()
    if (old_blas_threads != blas_threads) {
      blas_set_num_threads(blas_threads)
      on.exit(blas_set_num_threads(old_blas_threads), add = TRUE)
    }
  } else if (Sys.info()["sysname"] == "Linux" && blas_threads == 1L) {
    message(
      "Tip: Install 'RhpcBLASctl' for optimal performance on Linux:\n",
      "  install.packages('RhpcBLASctl')"
    )
  }
  
  # Validate weighting argument
  weighting <- match.arg(weighting)
  
  # --- Validate object ---
  if (!inherits(object, "liblex")) {
    stop("'object' must be of class 'liblex'", call. = FALSE)
  }
  
  if (is.null(object$coefficients)) {
    stop("Cannot predict: object was built with mode = 'validate'", call. = FALSE)
  }
  
  # --- Validate adaptive_bandwidth ---
  if (!is.logical(adaptive_bandwidth) || length(adaptive_bandwidth) != 1L || is.na(adaptive_bandwidth)) {
    stop("'adaptive_bandwidth' must be TRUE or FALSE", call. = FALSE)
  }
  
  # --- Validate reliability_weighting ---
  if (!is.logical(reliability_weighting) || length(reliability_weighting) != 1L || is.na(reliability_weighting)) {
    stop("'reliability_weighting' must be TRUE or FALSE", call. = FALSE)
  }
  
  # --- Validate probs ---
  if (!is.numeric(probs) || any(is.na(probs)) || any(probs < 0) || any(probs > 1)) {
    stop("'probs' must be a numeric vector with values in [0, 1]", call. = FALSE)
  }
  
  # --- Validate range_prediction_limits ---
  if (!is.logical(range_prediction_limits) || length(range_prediction_limits) != 1L || is.na(range_prediction_limits)) {
    stop("'range_prediction_limits' must be TRUE or FALSE", call. = FALSE)
  }
  
  # --- Validate residual_cutoff ---
  if (!is.null(residual_cutoff)) {
    if (!is.numeric(residual_cutoff) || length(residual_cutoff) != 1L || residual_cutoff <= 0) {
      stop("'residual_cutoff' must be a single positive number or NULL", call. = FALSE)
    }
  }
  
  # --- Validate enforce_indices ---
  if (!is.null(enforce_indices)) {
    if (!is.numeric(enforce_indices) || any(enforce_indices < 1L)) {
      stop("'enforce_indices' must be positive integers or NULL", call. = FALSE)
    }
    enforce_indices <- as.integer(enforce_indices)
    n_models <- nrow(object$coefficients$B)
    if (any(enforce_indices > n_models)) {
      stop("'enforce_indices' contains values exceeding number of models (", n_models, ")", 
           call. = FALSE)
    }
  }
  
  # --- Validate verbose ---
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("'verbose' must be TRUE or FALSE", call. = FALSE)
  }
 
  # --- Validate diss_method ---
  stored_diss <- object$dissimilarity$diss_method
  
  if (missing(diss_method)) {
    if (is.character(stored_diss) && stored_diss == "Precomputed dissimilarity matrix") {
      stop(
        "Model built with precomputed dissimilarity matrix; ",
        "'diss_method' argument is required for prediction",
        call. = FALSE
      )
    }
    diss_method <- stored_diss
  } else {
    # Validate user-provided diss_method
    if (!inherits(diss_method, "diss_method")) {
      stop("'diss_method' must be a diss_*() object", call. = FALSE)
    }
    # User provided diss_method — warn if different from stored
    if (inherits(stored_diss, "diss_method") && !identical(diss_method, stored_diss)) {
      warning(
        "Overriding stored 'diss_method' with user-supplied value",
        call. = FALSE
      )
    }
  }

  # --- Validate newdata ---
  required_vars <- colnames(object$coefficients$B)
  
  if (!all(required_vars %in% colnames(newdata))) {
    missing_vars <- setdiff(required_vars, colnames(newdata))
    stop(
      "Missing predictor variables in 'newdata': ",
      paste(missing_vars[seq_len(min(5L, length(missing_vars)))], collapse = ", "),
      if (length(missing_vars) > 5L) paste0(" ... and ", length(missing_vars) - 5L, " more"),
      call. = FALSE
    )
  }
  
  newdata <- as.matrix(newdata[, required_vars, drop = FALSE])
  
  # --- Compute dissimilarities between newdata and anchor observations ---
  if (inherits(diss_method, c("diss_pca", "diss_pls"))) {
    # Projection-based methods: compute distances in score space
    
    # Project newdata onto the stored projection
    scores_newdata <- predict(object$dissimilarity$projection, newdata)
    
    # Get scores for anchor observations
    scores_anchors <- object$dissimilarity$projection$scores[object$anchor_indices, , drop = FALSE]
    
    # Standardize using moments from all training scores (not just anchors)
    score_center <- get_column_means(object$dissimilarity$projection$scores)
    score_scale <- get_column_sds(object$dissimilarity$projection$scores)
    
    scores_anchors_scaled <- scale(scores_anchors, center = score_center, scale = score_scale)
    scores_newdata_scaled <- scale(scores_newdata, center = score_center, scale = score_scale)
    
    # Euclidean distance in standardized score space
    dsmxu <- dissimilarity(
      Xr = scores_anchors_scaled,
      Xu = scores_newdata_scaled,
      diss_method = diss_euclidean(center = FALSE, scale = FALSE)
    )
    
  } else {
    # Non-projection methods: compute distances in original variable space
    # Use stored neighborhood centers as reference points
    
    # Apply stored global centering/scaling
    # Scale the local centers and newdata using the global center and 
    # scaling vectors 
    local_centers_scaled <- scale(
      object$scaling$local_x_center,
      center = object$scaling$center,
      scale = object$scaling$scale
    )
    newdata_scaled <- scale(
      newdata,
      center = object$scaling$center,
      scale = object$scaling$scale
    )
    
    # Copy diss_method to avoid modifying the original object
    # Disable center/scale since already applied above
    diss_method_copy <- diss_method
    diss_method_copy$center <- FALSE
    diss_method_copy$scale <- FALSE
   
    dsmxu <- dissimilarity(
      Xr = local_centers_scaled,
      Xu = newdata_scaled,
      diss_method = diss_method_copy
    )
  }
  
  # --- Compute GH distance for newdata (if available) ---
  gh_newdata <- NULL
  if (!is.null(object$gh)) {
    # Project newdata onto the stored GH projection space
    gh_scores_newdata <- predict(object$gh$projection, newdata)
    
    # GH = distance from each observation to the training centroid
    gh_center <- colMeans(object$gh$projection$scores)
    
    gh_newdata <- as.vector(dissimilarity(
      Xr = gh_scores_newdata,
      Xu = t(gh_center),
      diss_method = diss_euclidean(center = TRUE, scale = TRUE)
    )$dissimilarity)
  }
  
  # --- Progress message for model retrieval ---
  if (verbose) {
    diss_type <- class(diss_method)[1L]
    
    if (inherits(diss_method, "diss_correlation")) {
      if (!is.null(diss_method$ws)) {
        cat("Retrieving models using correlation dissimilarity with window size", 
            diss_method$ws, "...\n")
      } else {
        cat("Retrieving models using correlation dissimilarity with full window...\n")
      }
    } else {
      cat("Retrieving models using", sub("^diss_", "", diss_type), "dissimilarity...\n")
    }
  }
  
  # --- Identify nearest neighbors for each new observation ---
  # Neighbors are sorted by dissimilarity (ascending order)
  
  k_min <- diss_threshold <- NULL
  if (inherits(object$neighbors, "neighbors_diss")) {
    k_min <- object$neighbors$k_min
    k_max <- object$neighbors$k_max
    diss_threshold <- object$optimal_params$diss_threshold
  } else {
    k_min <- k_max <- object$optimal_params$k    
  }
  
  if (k_max > nrow(dsmxu$dissimilarity)) {
    k_max <- nrow(dsmxu$dissimilarity)
    message("Number of available models is", nrow(dsmxu$dissimilarity)) 
  }
  
  if (k_min > nrow(dsmxu$dissimilarity)) {
    k_min <- nrow(dsmxu$dissimilarity)
  }
  
  neighbor_indices <- top_k_neighbors(
    D = dsmxu$dissimilarity,
    k_min = k_min,
    k_max = k_max,
    threshold = diss_threshold
  )
  
  # --- Filter models with high residuals (optional) ---
  # Models exceeding residual_cutoff are penalized with max dissimilarity
  # to prevent their selection as neighbors
  if (!is.null(object$residuals) && !is.null(residual_cutoff)) {
    abs_res <- abs(object$residuals)
    
    # Flag models with residuals exceeding cutoff
    high_residual_flag <- abs_res >= residual_cutoff
    high_residual_flag[is.na(high_residual_flag)] <- TRUE
    
    # Penalize high-residual models by assigning max dissimilarity
    neighbor_diss <- sapply(
      seq_along(neighbor_indices),
      FUN = function(col_idx, diss, nn, flag) {
        knns <- nn[[col_idx]]
        d <- diss[knns, col_idx]
        # Assign max dissimilarity to flagged models
        d[flag[knns]] <- max(d)
        d
      },
      diss = dsmxu$dissimilarity,
      nn = neighbor_indices,
      flag = high_residual_flag
    )
    
  } else {
    # No residual filtering: extract sorted dissimilarities for neighbors
    # neighbor_diss <- apply(
    #   dsmxu$dissimilarity,
    #   MARGIN = 2L,
    #   FUN = sort
    # )[seq_len(k_max), , drop = FALSE]
    
    neighbor_diss <- lapply(
      seq_along(neighbor_indices),
      function(i) dsmxu$dissimilarity[neighbor_indices[[i]], i]
    )
    
    high_residual_flag <- rep(FALSE, nrow(object$coefficients$B))
  }
  
  # --- Enforce specific models into all neighborhoods (optional) ---
  # Enforced models are prepended with minimal dissimilarity to ensure selection
  if (!is.null(enforce_indices)) {
    neighbor_indices <- lapply(neighbor_indices, function(x) c(enforce_indices, x))
    
    # Assign minimal dissimilarity to enforced models
    min_diss_per_obs <- lapply(
      neighbor_diss, 
      FUN = function(x, repetitions) rep(min(x), repetitions), 
      repetitions = length(enforce_indices)
    ) 
    # add the minimal dissimilarity for the enforced models to the existing 
    # neighbor dissimilarities
    neighbor_diss <- Map(c, min_diss_per_obs, neighbor_diss)
  }
  
  # --- Compute neighbor weights ---
  if (weighting == "none") {
    # Equal weights for all neighbors
    dweights <- lapply(neighbor_diss, function(x) rep(1 / length(x), length(x)))
  } else {
    # Normalize dissimilarities to [0, 1] range per observation
    dweights <- lapply(neighbor_diss, function(d) d / max(d))
    
    # Apply kernel weighting function
    sigma <- if (adaptive_bandwidth) NULL else 0.5
    dweights <- compute_pred_weights(dweights, weighting = weighting, sigma = sigma)
  }
  
  # --- Apply reliability weighting (based on model residuals) ---
  if (reliability_weighting && !is.null(object$residuals)) {
    # Extract residuals for neighbor models
    neighbor_residuals <- lapply(neighbor_indices, function(idx) object$residuals[idx])
    
    # Reliability: inverse of squared residual magnitude
    # Scale by median absolute residual to make dimensionless
    eps <- median(abs(object$residuals), na.rm = TRUE) * 0.1
    eps <- max(eps, .Machine$double.eps)  # ensure non-zero
    
    reliability <- lapply(neighbor_residuals, function(r) {
      rel <- 1 / (1 + (r / eps)^2)
      rel[is.na(rel)] <- 0.5  # models with missing residuals get neutral weight
      rel
    })
    
    # Combine distance and reliability weights
    dweights <- Map(`*`, dweights, reliability)
    
    # Re-normalize
    dweights <- lapply(dweights, function(w) {
      s <- sum(w, na.rm = TRUE)
      if (s == 0) s <- 1
      w / s
    })
  }
  
  # --- Prepare coefficient library and scaling parameters ---
  # Standard model without dissimilarity predictors
  coef_library <- cbind(
    object$coefficients$B0,
    object$coefficients$B
  )
  local_scale <- object$scaling$local_x_scale
  
  
  # --- Extract local prediction subsets for each new observation ---
  local_subsets <- ith_pred_subsets(
    plslib = coef_library,
    Xu = newdata,
    xunn = neighbor_indices,
    xscale = local_scale,
    dxrxu = NULL
  )
  
  # --- Compute predictions for each new observation ---
  if (verbose) {
    cat("Computing predictions...\n")
  }
  
  predictions_raw <- foreach(
    i = seq_len(nrow(newdata)),
    iset = local_subsets
  ) %mydo% {
    ith_pred_cpp(
      plslib = iset$iplslib,
      xscale = iset$ixscale,
      Xu = iset$ixu,
      dxrxu = iset$idxrxu
    )
  }
  
  # --- Combine predictions into matrix ---
  max_len <- max(lengths(predictions_raw))
  predictions_raw <- t(do.call("cbind", lapply(predictions_raw, function(x) {
    length(x) <- max_len
    x
  })))
  
  colnames(predictions_raw) <- paste0("expert_", seq_len(ncol(predictions_raw)))
  rownames(predictions_raw) <- rownames(newdata) %||% seq_len(nrow(newdata))
  
  
  
  max_len_weights <- max(lengths(dweights))
  dweights <- t(do.call("cbind", lapply(dweights, function(x) {
    length(x) <- max_len_weights
    x
  })))

  # --- Compute weighted predictions and uncertainty ---
  # Transpose weights to match predictions matrix layout
  rownames(dweights) <- rownames(predictions_raw)
  colnames(dweights) <- colnames(predictions_raw)
  
  # Weighted predictions per expert
  weighted_predictions <- predictions_raw * dweights
  
  # Weighted mean prediction per observation
  pred_mean <- rowSums(weighted_predictions, na.rm = TRUE)
  
  # Weighted standard deviation
  centered_dev <- sweep(predictions_raw, 1L, pred_mean, FUN = "-")^2
  pred_var <- rowSums(centered_dev * dweights, na.rm = TRUE)
  pred_sd <- sqrt(pred_var)

  # Weighted quantiles (exclude last column to avoid edge effects)
  pred_quantiles <- weighted_quantiles(
    predictions_raw,
    weights = dweights,
    probs = probs,
    exclude_last = TRUE
  )
  
  # Assemble predictions data frame
  predictions <- data.frame(
    pred = pred_mean,
    pred_sd = pred_sd,
    pred_quantiles
  )
  
  # --- Add GH distance to predictions ---
  predictions$gh <- gh_newdata
  
  # --- Preserve row names ---
  if (!is.null(rownames(newdata))) {
    rownames(predictions) <- rownames(newdata)
  }
  
  # --- Assemble output list ---
  result <- list(
    predictions = predictions,
    neighbors = list(
      indices = neighbor_indices,
      dissimilarities = neighbor_diss
    ),
    expert_predictions = list(
      weights = dweights,
      predictions = predictions_raw,
      weighted = weighted_predictions
    )
  )

  # --- Compute prediction limits from neighbor statistics ---
  
  if ( inherits(object$neighbors, "neighbors_diss") ) {
    nn_stats_indx <- which(object$neighbors$threshold == object$optimal_params[[1]])
  } else {
    nn_stats_indx <- which(object$neighbors$k == object$optimal_params[[1]])
  }
  
  neighborhood_stats <- object$neighborhood_stats[[nn_stats_indx]]

  # Extract min/max response values across neighbors for each observation
  min_yr <- sapply(
    seq_len(nrow(newdata)),
    FUN = function(i, stats, nn) min(stats[nn[[i]], "5%"]),
    stats = neighborhood_stats,
    nn = neighbor_indices
  )
  max_yr <- sapply(
    seq_len(nrow(newdata)),
    FUN = function(i, stats, nn) max(stats[nn[[i]], "95%"]),
    stats = neighborhood_stats,
    nn = neighbor_indices
  )
  
  result$predictions$min_yr <- min_yr
  result$predictions$max_yr <- max_yr
  
  # Flag predictions outside neighborhood range
  result$predictions$below_min <- result$predictions$pred < min_yr
  result$predictions$above_max <- result$predictions$pred > max_yr
  
  # Clip predictions to neighborhood range (optional)
  if (range_prediction_limits) {
    below <- result$predictions$below_min
    above <- result$predictions$above_max
    result$predictions$pred[below] <- min_yr[below]
    result$predictions$pred[above] <- max_yr[above]
  }
  
  result
}


compute_pred_weights <- function(
    xudss, 
    weighting = "gaussian",
    sigma = NULL,
    gamma = 0.3, 
    normalize = TRUE
) {
  valid_methods <- c(
    "triangular", "quartic", "triweight", "tricube",
    "parabolic", "gaussian", "cauchy"
  )
  weighting <- match.arg(weighting, valid_methods)
  
  # Normalize distances per observation to [0, 1]
  dscaled <- lapply(xudss, function(d) {
    d <- d / max(d)
    pmin(pmax(d, 0), 1)
  })
  
  # Adaptive bandwidth: use median neighbor distance per observation
  if (weighting == "gaussian") {
    if (is.null(sigma)) {
      sigma <- vapply(xudss, function(d) median(d) * 0.5, numeric(1))
      sigma[sigma == 0] <- min(sigma[sigma > 0], na.rm = TRUE)
      sigma[!is.finite(sigma)] <- 0.5
    }
    single_sigma <- length(sigma) == 1L
  }
  
  # Compute kernel weights
  dweights <- switch(
    weighting,
    triangular = lapply(dscaled, function(d) 1 - d),
    parabolic  = lapply(dscaled, function(d) 1 - d^2),
    quartic    = lapply(dscaled, function(d) (1 - d^2)^2),
    triweight  = lapply(dscaled, function(d) (1 - d^2)^3),
    tricube    = lapply(dscaled, function(d) (1 - d^3)^3),
    gaussian   = {
      if (single_sigma) {
        s2 <- 2 * sigma^2
        lapply(dscaled, function(d) exp(-(d^2) / s2))
      } else {
        Map(function(d, s) exp(-(d^2) / (2 * s^2)), dscaled, sigma)
      }
    },
    cauchy = lapply(dscaled, function(d) 1 / (1 + (d / gamma)^2))
  )
  
  dweights <- lapply(dweights, function(w) {
    w[!is.finite(w)] <- 0
    pmax(w, 0)
  })
  
  if (normalize) {
    dweights <- lapply(dweights, function(w) {
      s <- sum(w, na.rm = TRUE)
      if (s == 0) s <- 1
      w / s
    })
  }
  
  attr(dweights, "weighting") <- weighting
  attr(dweights, "sigma") <- sigma
  dweights
}


#' @aliases liblex
#' @export
plot.liblex <- function(x, ...) {
  if (!inherits(x, "liblex")) {
    stop("'x' must be an object of class 'liblex'", call. = FALSE)
  }
  
  
  # Set up 2-panel layout
  
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))
  
  # --- Panel 1: Best RMSE per neighborhood parameter ---
  if (!is.null(x$results)) {
    bresult <- NULL
    
    if (inherits(x$neighbors, "neighbors_k")) {
      for (k in unique(x$results$k)) {
        kth_r <- x$results[x$results$k == k, , drop = FALSE]
        kth_r <- kth_r[which.min(kth_r$rmse), , drop = FALSE]
        bresult <- rbind(bresult, kth_r)
      }
      param <- "k"
      param_name <- "Number of neighbors (k)"
    } else {
      for (thr in unique(x$results$diss_threshold)) {
        kth_r <- x$results[x$results$diss_threshold == thr, , drop = FALSE]
        kth_r <- kth_r[which.min(kth_r$rmse), , drop = FALSE]
        bresult <- rbind(bresult, kth_r)
      }
      param <- "diss_threshold"
      param_name <- "Dissimilarity threshold"
    }
    
    plot(
      bresult[[param]],
      bresult$rmse,
      type = "b",
      pch = 16,
      col = "dodgerblue",
      xlab = param_name,
      ylab = "RMSE",
      main = "Best RMSE at each neighbor selection parameter"
    )
    grid(lty = 1)
    
    # Mark optimal point
    opt_idx <- which.min(bresult$rmse)
    points(
      bresult[[param]][opt_idx],
      bresult$rmse[opt_idx],
      pch = 16,
      col = "red",
      cex = 1.5
    )
  } else {
    plot.new()
    text(0.5, 0.5, "No validation results available\n(tune = FALSE)", cex = 1.2)
  }
  
  # --- Panel 2: Neighborhood centroids ---
  if (!is.null(x$scaling$local_x_center)) {
    centroids <- x$scaling$local_x_center
    coords <- as.numeric(colnames(centroids))
    
    if (any(is.na(coords))) {
      coords <- seq_len(ncol(centroids))
      xlab_text <- "Variable index"
    } else {
      xlab_text <- "Spectral coordinate"
    }
    
    matplot(
      coords,
      t(centroids),
      col = rgb(0, 0, 1, 0.3),
      lty = 1,
      type = "l",
      xlab = xlab_text,
      ylab = "Intensity",
      main = paste0("Neighborhood centroids (n = ", nrow(centroids), ")")
    )
    grid(lty = 1)
  } else {
    plot.new()
    text(0.5, 0.5, "No centroids available\n(mode = 'validate')", cex = 1.2)
  }
  
  invisible(NULL)
}

