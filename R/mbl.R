#' @title Memory-based learning (mbl)
#' @description
#' \loadmathjax
#' Memory-based learning (a.k.a. instance-based learning or local regression)
#' is a non-linear lazy learning approach for predicting a response variable
#' from predictor variables. For each observation in a prediction set, a local
#' regression is fitted using a subset of similar observations (nearest
#' neighbors) from a reference set. This function does not produce a global
#' model.
#'
#' @usage
#' mbl(Xr, Yr, Xu, Yu = NULL,
#'     neighbors,
#'     diss_method = diss_pca(ncomp = ncomp_by_opc()),
#'     diss_usage = c("none", "predictors", "weights"),
#'     fit_method = fit_wapls(min_ncomp = 3, max_ncomp = 15),
#'     spike = NULL, group = NULL,
#'     gh = FALSE,
#'     control = mbl_control(),
#'     verbose = TRUE, seed = NULL,
#'     k, k_diss, k_range, method, pc_selection,
#'     center, scale, documentation, ...)
#'
#' @usage \method{plot}{mbl}(x, what = c("validation", "gh"), metric = "rmse", ncomp = c(1, 2), ...)
#'
#' @param Xr A matrix of predictor variables for the reference data
#'   (observations in rows, variables in columns). Column names are required.
#' @param Yr A numeric vector or single-column matrix of response values
#'   corresponding to \code{Xr}. NA values are not permitted.
#' @param Xu A matrix of predictor variables for the data to be predicted
#'   (observations in rows, variables in columns). Must have the same column
#'   names as \code{Xr}.
#' @param Yu An optional numeric vector or single-column matrix of response
#'   values corresponding to \code{Xu}. Used for computing prediction
#'   statistics. Default is \code{NULL}.
#' @param neighbors A neighbor selection object specifying how to select
#'   neighbors. Use \code{\link{neighbors_k}()} for fixed k-nearest neighbors
#'   or \code{\link{neighbors_diss}()} for dissimilarity threshold-based
#'   selection.
#' @param diss_method A dissimilarity method object or a precomputed
#'   dissimilarity matrix. Available constructors:
#'   \itemize{
#'     \item \code{\link{diss_pca}()}: Mahalanobis distance in PCA score space. 
#'     This is the default where the number of components is optimized using 
#'     side information (see \code{ncomp_by_opc}()).
#'     \item \code{\link{diss_pls}()}: Mahalanobis distance in PLS score space
#'     \item \code{\link{diss_euclidean}()}: Euclidean distance
#'     \item \code{\link{diss_mahalanobis}()}: Mahalanobis distance
#'     \item \code{\link{diss_cosine}()}: Cosine dissimilarity
#'     \item \code{\link{diss_correlation}()}: Correlation-based dissimilarity
#'   }
#'   A precomputed matrix can also be passed. When \code{diss_usage =
#'   "predictors"}, it must be square with dimensions
#'   \code{(nrow(Xr) + nrow(Xu))} and zeros on the diagonal. Otherwise, it must
#'   have \code{nrow(Xr)} rows and \code{nrow(Xu)} columns.
#' @param diss_usage How dissimilarity information is used in local models:
#'   \itemize{
#'     \item \code{"none"} (default): dissimilarities used only for neighbor
#'       selection
#'     \item \code{"predictors"}: local dissimilarity matrix columns added as
#'       predictors
#'     \item \code{"weights"}: neighbors weighted by dissimilarity using a
#'       tricubic function
#'   }
#' @param fit_method A local fitting method object. Available constructors:
#'   \itemize{
#'     \item \code{\link{fit_pls}()}: Partial least squares regression
#'     \item \code{\link{fit_wapls}()}: Weighted average PLS (default)
#'     \item \code{\link{fit_gpr}()}: Gaussian process regression
#'   }
#' @param spike An integer vector indicating indices of observations in
#'   \code{Xr} to force into (positive values) or exclude from (negative values)
#'   all neighborhoods. Default is \code{NULL}. Spiking does not change
#'   neighborhood size; forced observations displace the most distant neighbors.
#' @param gh Logical indicating whether to compute global Mahalanobis (GH)
#'   distances. Default is \code{FALSE}. GH distances measure how far each
#'   observation lies from the center of the reference set in PLS score space.
#'   The computation uses a fixed methodology: PLS projection with the number
#'   of components selected via \code{\link{ncomp_by_opc}()} (capped at 40).
#'   This is independent of the \code{diss_method} argument.
#' @param group An optional factor assigning group labels to \code{Xr}
#'   observations (e.g., measurement batches). Used to avoid pseudo-replication
#'   in cross-validation: when one observation is held out, all observations
#'   from its group are also removed.
#' @param control A list from \code{\link{mbl_control}()} specifying validation
#'   type, tuning options, and other settings.
#' @param verbose Logical indicating whether to display a progress bar.
#'   Default is \code{TRUE}. Not shown during parallel execution.
#' @param seed An integer for random number generation, enabling reproducible
#'   cross-validation results. Default is \code{NULL}.
#'   
#' @param x An object of class \code{mbl} (as returned by \code{\link{mbl}}).
#' @param what Character vector specifying what to plot. Options are
#'   \code{"validation"} (validation statistics) and/or \code{"gh"} (PLS scores
#'   used for GH distance computation). Default is both.
#' @param metric Character string specifying which validation statistic to plot.
#'   Options are \code{"rmse"}, \code{"st_rmse"}, or \code{"r2"}. Only used when
#'   \code{"validation"} is in \code{what}.
#' @param ncomp Integer vector of length 1 or 2 specifying which PLS components
#'   to plot. Default is \code{c(1, 2)}. Only used when \code{"gh"} is in
#'   \code{what}.
#' @param ... Additional arguments passed to \code{\link[graphics]{plot}}. 
#' Currently unused for `mbl()`.
#' 
#' @param k Deprecated.
#' @param k_diss Deprecated.
#' @param k_range Deprecated.
#' @param method Deprecated.
#' @param pc_selection Deprecated.
#' @param center Deprecated.
#' @param scale Deprecated.
#' @param documentation Deprecated.
#' @param ... Additional arguments (currently unused).
#' 
#' @details
#' ## Spiking
#' The \code{spike} argument forces specific reference observations into or out
#' of neighborhoods. Positive indices are always included; negative indices are
#' always excluded. When observations are forced in, the most distant neighbors
#' are displaced to maintain neighborhood size. See Guerrero et al. (2010).
#'
#' ## Dissimilarity usage
#' When \code{diss_usage = "predictors"}, the local dissimilarity matrix columns
#' are appended as additional predictor variables, which can improve predictions
#' (Ramirez-Lopez et al., 2013a).
#'
#' When \code{diss_usage = "weights"}, neighbors are weighted using a tricubic
#' function (Cleveland and Devlin, 1988; Naes et al., 1990):
#'
#' \mjdeqn{W_{j} = (1 - v^{3})^{3}}{W_j = (1 - v^3)^3}
#'
#' where \mjeqn{v = d(xr_i, xu_j) / \max(d)}{v = d(xr_i, xu_j) / max(d)}.
#'
#' ## GH distance
#' The global Mahalanobis distance (GH) measures how far each observation lies
#' from the center of the reference set. It is always computed using a PLS
#' projection with the number of components optimized via
#' \code{\link{ncomp_by_opc}()} (maximum 40 components or \code{nrow(Xr)},
#' whichever is smaller). This methodology is fixed and independent of the
#' \code{diss_method} specified for neighbor selection.
#'
#' GH distances are useful for identifying extrapolation: observations with
#' high GH values lie far from the calibration space and may yield unreliable
#' predictions.
#'
#' ## Grouping
#' The \code{group} argument enables leave-group-out cross-validation. When
#' \code{validation_type = "local_cv"} in \code{\link{mbl_control}()}, the
#' \code{p} parameter refers to the proportion of groups (not observations)
#' retained per iteration.
#'
#' @return A list of class \code{mbl} containing:
#' \itemize{
#'   \item \code{control}: control parameters from \code{control}
#'   \item \code{fit_method}: fit constructor from \code{fit_method} 
#'   \item \code{Xu_neighbors}: list with neighbor indices and dissimilarities
#'   \item \code{dissimilarities}: dissimilarity method and matrix (if
#'     \code{return_dissimilarity = TRUE} in \code{control})
#'   \item \code{n_predictions}: number of predictions made
#'   \item \code{gh}: GH distances for \code{Xr} and \code{Xu} (if
#'     \code{gh = TRUE})
#'   \item \code{validation_results}: validation statistics by method
#'   \item \code{results}: list of data.frame objects with predictions, one per
#'     neighborhood size
#'   \item \code{seed}: the seed value used
#' }
#'
#' Each results table contains:
#' \itemize{
#'   \item \code{o_index}: observation index
#'   \item \code{k}: number of neighbors used
#'   \item \code{k_diss}, \code{k_original}: (\code{neighbors_diss} only)
#'     threshold and original count
#'   \item \code{ncomp}: (\code{fit_pls} only) number of PLS components
#'   \item \code{min_ncomp}, \code{max_ncomp}: (\code{fit_wapls} only)
#'     component range
#'   \item \code{yu_obs}, \code{pred}: observed and predicted values
#'   \item \code{yr_min_obs}, \code{yr_max_obs}: response range in neighborhood
#'   \item \code{index_nearest_in_Xr}, \code{index_farthest_in_Xr}: neighbor
#'     indices
#'   \item \code{y_nearest}, \code{y_farthest}: neighbor response values
#'   \item \code{diss_nearest}, \code{diss_farthest}: neighbor dissimilarities
#'   \item \code{y_nearest_pred}: (NNv validation) leave-one-out prediction
#'   \item \code{loc_rmse_cv}, \code{loc_st_rmse_cv}: (local_cv validation) CV
#'     statistics
#'   \item \code{loc_ncomp}: (local dissimilarity only) components used locally
#' }
#'
#' @author
#' \href{https://orcid.org/0000-0002-5369-5120}{Leonardo Ramirez-Lopez} and
#' Antoine Stevens
#'
#' @references
#' Cleveland, W. S., and Devlin, S. J. 1988. Locally weighted regression: an
#' approach to regression analysis by local fitting. Journal of the American
#' Statistical Association 83:596-610.
#'
#' Guerrero, C., Zornoza, R., Gomez, I., Mataix-Beneyto, J. 2010. Spiking of
#' NIR regional models using observations from target sites: Effect of model
#' size on prediction accuracy. Geoderma 158:66-77.
#'
#' Naes, T., Isaksson, T., Kowalski, B. 1990. Locally weighted regression and
#' scatter correction for near-infrared reflectance data. Analytical Chemistry
#' 62:664-673.
#'
#' Ramirez-Lopez, L., Behrens, T., Schmidt, K., Stevens, A., Dematte, J.A.M.,
#' Scholten, T. 2013a. The spectrum-based learner: A new local approach for
#' modeling soil vis-NIR spectra of complex data sets. Geoderma 195-196:268-279.
#'
#' Ramirez-Lopez, L., Behrens, T., Schmidt, K., Viscarra Rossel, R., Dematte,
#' J.A.M., Scholten, T. 2013b. Distance and similarity-search metrics for use
#' with soil vis-NIR spectra. Geoderma 199:43-53.
#'
#' Rasmussen, C.E., Williams, C.K. 2006. Gaussian Processes for Machine
#' Learning. MIT Press.
#'
#' Shenk, J., Westerhaus, M., Berzaghi, P. 1997. Investigation of a LOCAL
#' calibration procedure for near infrared instruments. Journal of Near
#' Infrared Spectroscopy 5:223-232.
#'
#' @seealso
#' \code{\link{mbl_control}}, \code{\link{neighbors_k}},
#' \code{\link{neighbors_diss}}, \code{\link{diss_pca}}, \code{\link{diss_pls}},
#' \code{\link{fit_pls}}, \code{\link{fit_wapls}}, \code{\link{fit_gpr}},
#' \code{\link{search_neighbors}}
#'
#' @examples
#' \dontrun{
#' library(prospectr)
#' data(NIRsoil)
#'
#' # Preprocess: detrend + first derivative with Savitzky-Golay
#' sg_det <- savitzkyGolay(
#'   detrend(NIRsoil$spc, wav = as.numeric(colnames(NIRsoil$spc))),
#'   m = 1, p = 1, w = 7
#' )
#' NIRsoil$spc_pr <- sg_det
#'
#' # Split data
#' test_x <- NIRsoil$spc_pr[NIRsoil$train == 0 & !is.na(NIRsoil$CEC), ]
#' test_y <- NIRsoil$CEC[NIRsoil$train == 0 & !is.na(NIRsoil$CEC)]
#' train_x <- NIRsoil$spc_pr[NIRsoil$train == 1 & !is.na(NIRsoil$CEC), ]
#' train_y <- NIRsoil$CEC[NIRsoil$train == 1 & !is.na(NIRsoil$CEC)]
#'
#' # Example 1: Spectrum-based learner (Ramirez-Lopez et al., 2013)
#' ctrl <- mbl_control(validation_type = "NNv")
#'
#' sbl <- mbl(
#'   Xr = train_x,
#'   Yr = train_y,
#'   Xu = test_x,
#'   neighbors = neighbors_k(seq(40, 140, by = 20)),
#'   diss_method = diss_pca(ncomp = ncomp_by_opc(40)),
#'   fit_method = fit_gpr(),
#'   control = ctrl
#' )
#' sbl
#' plot(sbl)
#' get_predictions(sbl)
#'
#' # Example 2: With known Yu
#' sbl_2 <- mbl(
#'   Xr = train_x,
#'   Yr = train_y,
#'   Xu = test_x,
#'   Yu = test_y,
#'   neighbors = neighbors_k(seq(40, 140, by = 20)),
#'   fit_method = fit_gpr(),
#'   control = ctrl
#' )
#' plot(sbl_2)
#'
#' # Example 3: LOCAL algorithm (Shenk et al., 1997)
#' local_algo <- mbl(
#'   Xr = train_x,
#'   Yr = train_y,
#'   Xu = test_x,
#'   Yu = test_y,
#'   neighbors = neighbors_k(seq(40, 140, by = 20)),
#'   diss_method = diss_correlation(),
#'   diss_usage = "none",
#'   fit_method = fit_wapls(min_ncomp = 3, max_ncomp = 15),
#'   control = ctrl
#' )
#' plot(local_algo)
#'
#' # Example 4: Using dissimilarity as predictors
#' local_algo_2 <- mbl(
#'   Xr = train_x,
#'   Yr = train_y,
#'   Xu = test_x,
#'   Yu = test_y,
#'   neighbors = neighbors_k(seq(40, 140, by = 20)),
#'   diss_method = diss_pca(ncomp = ncomp_by_opc(40)),
#'   diss_usage = "predictors",
#'   fit_method = fit_wapls(min_ncomp = 3, max_ncomp = 15),
#'   control = ctrl
#' )
#' plot(local_algo_2)
#'
#' # Example 5: Parallel execution
#' library(doParallel)
#' n_cores <- min(2, parallel::detectCores())
#' clust <- makeCluster(n_cores)
#' registerDoParallel(clust)
#'
#' local_algo_par <- mbl(
#'   Xr = train_x,
#'   Yr = train_y,
#'   Xu = test_x,
#'   Yu = test_y,
#'   neighbors = neighbors_k(seq(40, 140, by = 20)),
#'   diss_method = diss_correlation(),
#'   fit_method = fit_wapls(min_ncomp = 3, max_ncomp = 15),
#'   control = ctrl
#' )
#'
#' registerDoSEQ()
#' try(stopCluster(clust))
#' }
#'
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

mbl <- function(
  Xr, Yr, Xu, Yu = NULL,
  neighbors,
  diss_method = diss_pca(ncomp = ncomp_by_opc()),
  diss_usage = c("none", "predictors", "weights"),
  fit_method = fit_wapls(min_ncomp = 3, max_ncomp = 15),
  spike = NULL,
  group = NULL,
  gh = FALSE,
  control = mbl_control(),
  verbose = TRUE,
  seed = NULL,
  k, k_diss, k_range, method, pc_selection,
  center, scale, documentation,
    ...
) {
  f_call <- match.call()
  
  # ---------------------------------------------------------------------------
  # Block removed arguments
  # ---------------------------------------------------------------------------
  if (!missing(k) || !missing(k_diss) || !missing(k_range)) {
    stop(
      "Arguments 'k', 'k_diss', 'k_range' have been removed.\n",
      "Use neighbors_k() or neighbors_diss() instead.\n",
      "Example: mbl(..., neighbors = neighbors_k(c(40, 80, 120)))",
      call. = FALSE
    )
  }
  
  if (!missing(method)) {
    stop(
      "Argument 'method' has been renamed to 'fit_method'.\n",
      "Use fit_pls(), fit_wapls(), or fit_gpr() constructors.\n",
      "Example: mbl(..., fit_method = fit_wapls(min_ncomp = 3, max_ncomp = 15))",
      call. = FALSE
    )
  }
  
  if (!missing(pc_selection)) {
    stop(
      "Argument 'pc_selection' has been removed.\n",
      "Component selection is now specified in diss_*() constructors.\n",
      "Example: mbl(..., diss_method = diss_pca(ncomp = ncomp_by_opc(40)))",
      call. = FALSE
    )
  }
  
  if (!missing(center) || !missing(scale)) {
    stop(
      "Arguments 'center' and 'scale' have been removed.\n",
      "These are now set in diss_*() and fit_*() constructors.\n",
      "Example: diss_pca(center = TRUE, scale = FALSE)",
      call. = FALSE
    )
  }
  
  if (!missing(documentation)) {
    stop(
      "Argument 'documentation' has been removed.",
      call. = FALSE
    )
  }
  
  # ---------------------------------------------------------------------------
  # Set BLAS threads to reduce overhead (restored on exit)
  # ---------------------------------------------------------------------------
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    old_blas_threads <- blas_get_num_procs()
    if (old_blas_threads != control$blas_threads) {
      blas_set_num_threads(control$blas_threads)
      on.exit(blas_set_num_threads(old_blas_threads), add = TRUE)
    }
  } else if (Sys.info()["sysname"] == "Linux" && control$blas_threads == 1L) {
    message(
      "Tip: Install 'RhpcBLASctl' for optimal performance on Linux:\n",
      "  install.packages('RhpcBLASctl')"
    )
  }
  # ---------------------------------------------------------------------------
  # Validate constructor arguments
  # ---------------------------------------------------------------------------
  if (missing(neighbors)) {
    stop("'neighbors' is required. Use neighbors_k() or neighbors_diss().", 
         call. = FALSE)
  }
  
  if (!inherits(neighbors, "neighbors")) {
    stop(
      "'neighbors' must be created by neighbors_k() or neighbors_diss().",
      call. = FALSE
    )
  }
  
  # ---------------------------------------------------------------------------
  # Validate inputs
  # ---------------------------------------------------------------------------
  # --- control validation ---
  if (!inherits(control, "mbl_control")) {
    stop("'control' must be created by mbl_control()", call. = FALSE)
  }
  
  if (!is.logical(verbose) || length(verbose) != 1L) {
    stop("'verbose' must be TRUE or FALSE.", call. = FALSE)
  }
  
  if (!is.matrix(Xr)) {
    Xr <- as.matrix(Xr)
  }
  
  if (!is.matrix(Xu)) {
    Xu <- as.matrix(Xu)
  }
  
  if (!is.null(Yu)) {
    if (!is.matrix(Yu)) {
      Yu <- as.matrix(Yu)
    }
  }
  
  if (!is.matrix(Yr)) {
    Yr <- as.matrix(Yr)
  }
  
  
  
  n_xr <- nrow(Xr)
  n_xu <- nrow(Xu)
  n_total <- n_xr + n_xu
  
  if (ncol(Xr) != ncol(Xu)) {
    stop("ncol(Xr) must equal ncol(Xu).", call. = FALSE)
  }
  
  if (ncol(Xr) < 4L) {
    stop("Xr must have at least 4 columns.", call. = FALSE)
  }
  
  if (NROW(Yr) != nrow(Xr)) {
    stop("length(Yr) must equal nrow(Xr).", call. = FALSE)
  }
  
  if (any(is.na(Yr))) {
    stop("NA values in Yr are not supported.", call. = FALSE)
  }
  
  if (!is.null(Yu) && NROW(Yu) != nrow(Xu)) {
    stop("length(Yu) must equal nrow(Xu).", call. = FALSE)
  }
  
  
  # ---------------------------------------------------------------------------
  # Validate diss_method
  # ---------------------------------------------------------------------------
  if (is.matrix(diss_method)) {
    # User passed a precomputed dissimilarity matrix
    diss_matrix <- diss_method
    diss_method <- NULL
    # center <- FALSE
    # scale <- FALSE
  } else if (!inherits(diss_method, "diss_method")) {
    stop(
      "'diss_method' must be a diss_*() constructor or a numeric matrix.\n",
      "Available constructors: diss_pca(), diss_pls(), diss_correlation(), ",
      "diss_euclidean(), diss_cosine(), diss_mahalanobis()",
      call. = FALSE
    )
  } else {
    # ---------------------------------------------------------------------------
    # Validate ncomp against data dimensions
    # ---------------------------------------------------------------------------
    if (!is.null(diss_method$ncomp)) {
      max_ncomp <- .get_max_ncomp(diss_method$ncomp)
      max_allowed <- min(n_total, ncol(Xr))
      
      if (max_ncomp > max_allowed) {
        warning(
          "Requested ncomp (", max_ncomp, ") exceeds max allowed by data dimensions (",
          max_allowed, ").\nIt will be capped to ", max_allowed, ".",
          call. = FALSE
        )
      }
    }
    diss_matrix <- NULL
    # center <- if (!is.null(diss_method$center)) diss_method$center else TRUE
    # scale <- if (!is.null(diss_method$scale)) diss_method$scale else FALSE
  }
  
  if (!inherits(fit_method, "fit_method")) {
    stop(
      "'fit_method' must be created by a fit_*() constructor.\n",
      "Available: fit_pls(), fit_wapls(), fit_gpr()",
      call. = FALSE
    )
  }
  
  diss_usage <- match.arg(diss_usage)
  
  # ---------------------------------------------------------------------------
  # Convert neighbors to internal format
  # ---------------------------------------------------------------------------
  if (inherits(neighbors, "neighbors_k")) {
    k <- neighbors$k
    k_diss <- NULL
    k_range <- NULL
    k_max <- max(k)
    k_diss_max <- NULL
    
    # Validate k against data
    if (k_max > n_xr) {
      stop(
        "k (", k_max, ") cannot exceed nrow(Xr) (", n_xr, ").",
        call. = FALSE
      )
    }
  } else {
    # neighbors_diss
    k <- NULL
    k_diss <- neighbors$threshold
    k_range <- c(neighbors$k_min, neighbors$k_max)
    k_max <- NULL
    k_diss_max <- max(k_diss)
    
    if (is.infinite(k_range[2])) {
      if (verbose) {
        message(
          "setting 'k_max' (", k_range[2], ") to nrow(Xr) (", n_xr, ").",
          call. = FALSE
        )
      }
      k_range[2] <- n_xr
    }
    
    # Validate k_range against data
    if (k_range[2] > n_xr) {
      stop(
        "k_max (", k_range[2], ") cannot exceed nrow(Xr) (", n_xr, ").",
        call. = FALSE
      )
    }
  }
  
  # ---------------------------------------------------------------------------
  # Extract from fit_method constructor
  # ---------------------------------------------------------------------------
  is_fit_pls <- inherits(fit_method, "fit_pls")
  is_fit_wapls <- inherits(fit_method, "fit_wapls")
  is_fit_gpr <- inherits(fit_method, "fit_gpr")
  
  # For fit_and_predict()
  if (is_fit_pls) {
    pred_method <- "pls"
    pls_c <- fit_method$ncomp
  } else if (is_fit_wapls) {
    pred_method <- "wapls"
    pls_c <- c(fit_method$min_ncomp, fit_method$max_ncomp)
  } else {
    pred_method <- "gpr"
    pls_c <- NULL
  }
  
  noise_variance <- fit_method$noise_variance
  fit_scale <- fit_method$scale
  
  # ---------------------------------------------------------------------------
  # Parallel setup
  # ---------------------------------------------------------------------------
  "%mydo%" <- get("%do%")
  if (control$allow_parallel && getDoParRegistered()) {
    "%mydo%" <- get("%dopar%")
  }
  
  
  # ---------------------------------------------------------------------------
  # Validate and set names
  # ---------------------------------------------------------------------------
  if (is.null(colnames(Xr))) {
    stop("Xr must have column names.", call. = FALSE)
  }
  
  if (is.null(colnames(Xu))) {
    stop("Xu must have column names.", call. = FALSE)
  }
  
  if (!identical(colnames(Xr), colnames(Xu))) {
    stop("Column names of Xr and Xu must be identical.", call. = FALSE)
  }
  
  rownames(Xr) <- seq_len(nrow(Xr))
  rownames(Xu) <- seq_len(nrow(Xu))
  
  # ---------------------------------------------------------------------------
  # Validation type
  # ---------------------------------------------------------------------------
  validation_type <- control$validation_type
  is_local_cv <- "local_cv" %in% validation_type
  is_nnv <- "NNv" %in% validation_type
  
  if (is_local_cv && is_nnv) {
    validation_type <- "both"
  }
  
  # ---------------------------------------------------------------------------
  # Validate group
  # ---------------------------------------------------------------------------
  if (!is.null(group)) {
    if (length(group) != nrow(Xr)) {
      stop("length(group) must equal nrow(Xr).", call. = FALSE)
    }
    group <- as.factor(group)
  }
  
  # ---------------------------------------------------------------------------
  # Setup for Xu
  # ---------------------------------------------------------------------------
  if (is_nnv && n_xu < 3L) {
    stop(
      "Nearest neighbor validation requires at least 3 observations in Xu.",
      call. = FALSE
    )
  }
  
  if (!is.null(Yu)) {
    Yu <- as.matrix(Yu)
    if (nrow(Yu) != n_xu) {
      stop("nrow(Yu) must equal nrow(Xu).", call. = FALSE)
    }
  }
  
  pre_nms_ng <- "Xu_"
  ln <- n_xu
  first_nn <- 1L
  y_output_name <- "yu_obs"
  y_hat_output_name <- "pred"
  val_summary_name <- "Yu_prediction_statistics"
  
  # ---------------------------------------------------------------------------
  # Validate fit_method ncomp against data dimensions
  # ---------------------------------------------------------------------------
  if (is_fit_pls) {
    if (fit_method$ncomp > min(n_xr, ncol(Xr))) {
      warning(
        "Requested ncomp (", fit_method$ncomp, ") in fit_pls() exceeds max allowed (",
        min(n_xr, ncol(Xr)), ").\nIt will be capped internally.",
        call. = FALSE
      )
    }
  } else if (is_fit_wapls) {
    max_allowed <- min(n_xr, ncol(Xr))
    if (fit_method$max_ncomp > max_allowed) {
      warning(
        "Requested max_ncomp (", fit_method$max_ncomp, ") in fit_wapls() exceeds max allowed (",
        max_allowed, ").\nIt will be capped internally.",
        call. = FALSE
      )
    }
  }
  
  has_projection <- FALSE
  # Copy diss_method to avoid modifying user's input
  diss_method_internal <- diss_method
  # ---------------------------------------------------------------------------
  # Compute neighborhoods
  # ---------------------------------------------------------------------------
  if (is.null(diss_matrix)) {
    # Force return_projection = TRUE when needed for diss_usage = "predictors"
    if (diss_usage == "predictors") {
      if (inherits(diss_method_internal, "diss_pca") || 
          inherits(diss_method_internal, "diss_pls")) {
        if (!isTRUE(diss_method_internal$return_projection)) {
          diss_method_internal$return_projection <- TRUE
        }
      }
    }
    # Using diss_method constructor
    neighborhoods <- get_neighbor_info(
      Xr = Xr, 
      Xu = Xu,
      diss_method = diss_method_internal, 
      Yr = Yr,
      neighbors = neighbors,
      spike = spike,
      return_dissimilarity = control$return_dissimilarity,
      diss_usage = diss_usage
    )
    
    diss_xr_xu <- neighborhoods$dissimilarity
    
    if (!is.null(neighborhoods$projection)) {
      diss_projection <- neighborhoods$projection
      has_projection <- TRUE
    } else {
      has_projection <- FALSE
    }
    
  } else {
    # Using precomputed diss_matrix
    has_projection <- FALSE
    diss_xr_xr <- NULL
    dim_diss <- dim(diss_matrix)
    
    if (diss_usage == "predictors") {
      if (dim_diss[1] != n_total || dim_diss[2] != n_total || 
          any(diag(diss_matrix) != 0)) {
        stop(
          "When diss_usage = 'predictors', the dissimilarity matrix must be ",
          "square (", n_total, " x ", n_total, ") with zeros on the diagonal.",
          call. = FALSE
        )
      }
  
      diss_xr_xr <- diss_matrix[1:n_xr, 1:n_xr]
      diss_xr_xu <- diss_matrix[1:n_xr, (n_xr + 1):n_total]
      
    } else {
      # diss_usage %in% c("weights", "none")
      if (dim_diss[1] != n_xr || dim_diss[2] != n_xu) {
        stop(
          "When diss_usage = 'weights' or 'none', the dissimilarity matrix ",
          "must have nrow(Xr) (", n_xr, ") rows and nrow(Xu) (", n_xu, ") columns.",
          call. = FALSE
        )
      }
      diss_xr_xu <- diss_matrix
    }
    
    neighborhoods <- diss_to_neighbors(
      diss_xr_xu,
      k = k_max, 
      k_diss = k_diss_max,
      k_range = k_range,
      spike = spike,
      return_dissimilarity = control$return_dissimilarity,
      skip_first = FALSE
    )
    
    neighborhoods$diss_xr_xr <- diss_xr_xr
  }
  
  # ---------------------------------------------------------------------------
  # Compute GH distances (independent of diss_method)
  # ---------------------------------------------------------------------------
  if (gh) {
    # FIXME: this is not ideal, as it requires an additional projection step.
    # in addition the optimisation of components is hard coded
    gh_projection <- ortho_projection(
      Xr = Xr, 
      Xu = Xu,
      Yr = Yr,
      ncomp = ncomp_by_opc(min(n_xr, 40L)),
      method = "pls",
      scale = fit_scale
    )
    
    gh_center <- colMeans(gh_projection$scores[1:n_xr, , drop = FALSE])
    
    neighborhoods$gh <- list(
      gh_Xr = as.vector(f_diss(
        gh_projection$scores[1:n_xr, , drop = FALSE],
        Xu = t(gh_center),
        diss_method = "mahalanobis",
        center = FALSE, 
        scale = FALSE
      )),
      gh_Xu = as.vector(f_diss(
        gh_projection$scores[(n_xr + 1):n_total, , drop = FALSE],
        Xu = t(gh_center),
        diss_method = "mahalanobis",
        center = FALSE, 
        scale = FALSE
      )),
      projection = gh_projection
    )
  } else {
    neighborhoods$gh <- NULL
  }
  
  # ---------------------------------------------------------------------------
  # Determine smallest neighborhood size
  # ---------------------------------------------------------------------------
  if (inherits(neighbors, "neighbors_k")) {
    min_k <- min(neighbors$k)
    smallest_neighborhood <- neighborhoods$neighbors[1:min_k, , drop = FALSE]
    smallest_n_neighbors <- colSums(!is.na(smallest_neighborhood))
  } else {
    # neighbors_diss
    min_threshold <- min(neighbors$threshold)
    min_diss <- neighborhoods$neighbors_diss <= min_threshold
    
    if (!is.null(spike)) {
      spike_in <- sum(spike > 0)
      if (spike_in > 0) {
        min_diss[1:spike_in, ] <- TRUE
      }
    }
    
    smallest_neighborhood <- neighborhoods$neighbors
    smallest_neighborhood[!min_diss] <- NA
    smallest_n_neighbors <- colSums(!is.na(smallest_neighborhood))
    smallest_n_neighbors[smallest_n_neighbors < neighbors$k_min] <- neighbors$k_min
    smallest_n_neighbors[smallest_n_neighbors > neighbors$k_max] <- neighbors$k_max
  }
  
  if (is_local_cv) {
    min_n_samples <- floor(min(smallest_n_neighbors) * control$p) - 1L
    
    if (inherits(neighbors, "neighbors_k")) {
      min_neighborhood <- min(neighbors$k)
    } else {
      min_neighborhood <- neighbors$k_min
    }
    min_cv_samples <- floor(min_neighborhood * (1 - control$p))
    
    if (min_cv_samples < 3L) {
      stop(
        "Local cross-validation requires at least 3 observations in ",
        "the hold-out set. The current cross-validation parameters ",
        "leave less than 3 observations in some neighborhoods.",
        call. = FALSE
      )
    }
  } else {
    min_n_samples <- smallest_n_neighbors - 1L
  }
  
  # ---------------------------------------------------------------------------
  # Validate PLS components vs neighborhood size
  # ---------------------------------------------------------------------------
  if (is_fit_pls || is_fit_wapls) {
    if (is_fit_pls) {
      max_pls <- fit_method$ncomp
    } else {
      max_pls <- fit_method$max_ncomp
    }
    
    if (any(min_n_samples < max_pls)) {
      stop(
        "More PLS components than observations in some neighborhoods.\n",
        "If 'local_cv' is being used, consider that some observations ",
        "in the neighborhoods are held out for local validation.",
        call. = FALSE
      )
    }
  }
  
  # ---------------------------------------------------------------------------
  # Prepare iteration neighborhoods
  # ---------------------------------------------------------------------------
  is_local_diss <- !is.null(diss_method) && 
    (inherits(diss_method, "diss_local_pca") || 
       inherits(diss_method, "diss_local_pls"))
  
  if (is_local_diss) {
    iter_neighborhoods <- ith_mbl_neighbor(
      Xr = Xr, Xu = Xu, Yr = Yr, Yu = Yu,
      diss_usage = "none",
      neighbor_indices = neighborhoods$neighbors,
      neighbor_diss = neighborhoods$neighbors_diss,
      group = group
    )
  } else {
    iter_neighborhoods <- ith_mbl_neighbor(
      Xr = Xr, Xu = Xu, Yr = Yr, Yu = Yu,
      diss_usage = diss_usage,
      neighbor_indices = neighborhoods$neighbors,
      neighbor_diss = neighborhoods$neighbors_diss,
      diss_xr_xr = neighborhoods$diss_xr_xr,
      group = group
    )
  }
  
  r_fields <- c(
    "o_index", "k_diss", "k_original", "k", "ncomp", "min_ncomp", "max_ncomp",
    y_output_name, y_hat_output_name, "yr_min_obs", "yr_max_obs",
    "index_nearest_in_Xr", "index_farthest_in_Xr",
    "y_nearest", "y_nearest_pred",
    "y_farthest", "diss_nearest", "diss_farthest",
    "loc_rmse_cv", "loc_st_rmse_cv", "loc_ncomp", "rep"
  )
  
  n_ith_result <- ifelse(
    inherits(neighbors, "neighbors_k"), 
    length(neighbors$k), 
    length(neighbors$threshold)
  )
  
  template_pred_results <- as.data.frame(
    matrix(
      NA, n_ith_result, length(r_fields),
      dimnames = list(NULL, r_fields)
    )
  )
  
  template_pred_results$rep[1] <- 0
  
  if (inherits(neighbors, "neighbors_k")) {
    template_pred_results$k <- neighbors$k
  } else {
    template_pred_results$k_diss <- neighbors$threshold
  }
  
  pg_bar_width <- 10
  n_characters <- nchar(n_xu)
  n_iter <- n_xu
  to_erase <- pg_bar_width + (2 * n_characters) + 8
  to_erase <- paste(rep(" ", to_erase), collapse = "")
  
  if (verbose) {
    cat("\033[32m\033[3mPredicting...\n\033[23m\033[39m")
  }
  
  pred_obs <- foreach(
    i = seq_len(n_iter),
    ith_observation = iter_neighborhoods,
    .inorder = FALSE,
    .export = c(
      "ortho_diss", "fit_and_predict", "pls_cv",
      "get_col_sds", "get_wapls_weights"
    ),
    .noexport = c("Xr", "Xu")
  ) %mydo% {
  
    ith_pred_results <- template_pred_results
    additional_results <- NULL
    ith_pred_results$o_index[] <- i
    
    if (is_local_diss) {
      ith_observation <- get_ith_local_neighbors(
        ith_xr = ith_observation$ith_xr,
        ith_xu = ith_observation$ith_xu,
        ith_yr = ith_observation$ith_yr,
        ith_yu = ith_observation$ith_yu,
        diss_usage = diss_usage,
        ith_neig_indices = ith_observation$ith_neig_indices,
        neighbors = neighbors,
        spike = spike,
        diss_method = diss_method,
        ith_group = ith_observation$ith_group,
        mbl_is_parallel = control$allow_parallel
      )
      ith_pred_results$loc_ncomp[] <- ith_observation$ith_ncomp
      additional_results$ith_neig_indices <- ith_observation$ith_neig_indices
      additional_results$ith_neigh_diss <- ith_observation$ith_neigh_diss
    }
    
    if (verbose) {
      cat(paste0("\033[34m\033[3m", i, "/", n_iter, "\033[23m\033[39m"))
      pb <- txtProgressBar(width = pg_bar_width, char = "\033[34m_\033[39m")
    }
    
    if (inherits(neighbors, "neighbors_diss")) {
      ith_diss <- ith_observation$ith_neigh_diss
      if (!is.null(spike)) {
        n_spike_hold <- sum(spike > 0)
        if (n_spike_hold > 0) {
          ith_diss[seq_len(n_spike_hold)] <- 0
        }
      }
      k_diss_vec <- neighbors$threshold
      k_min <- neighbors$k_min
      k_max <- neighbors$k_max
      ith_pred_results$k_original <- sapply(k_diss_vec, FUN = function(x, d) sum(d < x), d = ith_diss)
      ith_pred_results$k <- ith_pred_results$k_original
      ith_pred_results$k[ith_pred_results$k_original < k_min] <- k_min
      ith_pred_results$k[ith_pred_results$k_original > k_max] <- k_max
    } else {
      ith_pred_results$k <- neighbors$k
    }
    
    
    for (kk in seq_len(nrow(ith_pred_results))) {
      if (verbose) {
        setTxtProgressBar(pb, kk / nrow(ith_pred_results))
      }
      
      # If the sample has not been predicted before,
      # then create a model and predict it (useful only when k_diss is used)
      current_k <- ith_pred_results$k[kk]
      
      # Skip refitting if neighborhood size unchanged from previous iteration
      # (optimization for neighbors_diss where different thresholds may prorduce 
      # same k)
      if (current_k != ifelse(kk == 1, 0, ith_pred_results$k[kk - 1])) {
        
        # Subset predictors: when diss_usage == "predictors", include local 
        # dissimilarity columns alongside original spectral variables
        if (diss_usage == "predictors") {
          keep_cols <- c(
            1:current_k,
            (1 + ith_observation$n_k):ncol(ith_observation$ith_xr)
          )
          i_k_xr <- ith_observation$ith_xr[1:current_k, keep_cols]
          i_k_xu <- ith_observation$ith_xu[, keep_cols, drop = FALSE]
        } else {
          i_k_xr <- ith_observation$ith_xr[1:current_k, ]
          i_k_xu <- ith_observation$ith_xu
        }
        
        # for extracting some basic stats
        i_k_yr <- ith_observation$ith_yr[first_nn:current_k, , drop = FALSE]
        
        # i_k_yu <- ith_observation$ith_yu
        kth_diss <- ith_observation$ith_neigh_diss[first_nn:current_k]
        i_idx <- ith_observation$ith_neig_indices[first_nn:current_k]
        
        
        ith_pred_results$rep[kk] <- 0
        ith_yr_range <- range(i_k_yr)
        ith_pred_results$yr_min_obs[kk] <- ith_yr_range[1]
        ith_pred_results$yr_max_obs[kk] <- ith_yr_range[2]
        ith_pred_results$diss_farthest[kk] <- max(kth_diss)
        ith_pred_results$diss_nearest[kk] <- min(kth_diss)
        ith_pred_results$y_farthest[kk] <- i_k_yr[which.max(kth_diss)]
        ith_pred_results$y_nearest[kk] <- i_k_yr[which.min(kth_diss)]
        ith_pred_results$index_nearest_in_Xr[kk] <- i_idx[which.min(kth_diss)]
        ith_pred_results$index_farthest_in_Xr[kk] <- i_idx[which.max(kth_diss)]
        
        # Full neighbor set for model fitting (including spiked observations).
        # Stats above excluded spike via first_nn offset. This set is then
        # subsetted at each k iteration.
        i_k_yr <- ith_observation$ith_yr[1:current_k, , drop = FALSE]
        
        
        if (!is.null(group)) {
          i_k_group <- factor(ith_observation$ith_group[1:current_k])
        } else {
          i_k_group <- NULL
        }
        
        if (diss_usage == "weights") {
          # Weights are defined according to a tricubic function
          # as in Cleveland and Devlin (1988) and Naes and Isaksson (1990).
          std_kth_diss <- kth_diss / max(kth_diss)
          kth_weights <- (1 - (std_kth_diss^3))^3
          kth_weights[which(kth_weights == 0)] <- 1e-04
        } else {
          kth_weights <- rep(1, current_k)
        }

        # Local model fit and prediction
        i_k_pred <- fit_and_predict(
          x = i_k_xr,
          y = i_k_yr,
          pred_method = pred_method,
          scale = fit_scale,
          pls_c = pls_c,
          weights = kth_weights,
          newdata = i_k_xu,
          CV = is_local_cv,
          tune = control$tune_locally,
          group = i_k_group,
          p = control$p,
          number = control$number,
          noise_variance = noise_variance,
          range_prediction_limits = control$range_prediction_limits,
          pls_max_iter = fit_method$max_iter,
          pls_tol = fit_method$tol,
          seed = seed,
          algorithm = fit_method$method
        )
        
        ith_pred_results[[y_hat_output_name]][kk] <- i_k_pred$prediction
        
        # Track selected PLS components for NNv validation
        selected_pls <- NULL
        
        if (is_local_cv) {
          # Select best model: tuned or fixed components
          if (control$tune_locally) {
            best_row <- which.min(i_k_pred$validation$cv_results$rmse_cv)
          } else {
            best_row <- ifelse(is_fit_pls, fit_method$ncomp, 1)
          }
          # Store optimized component counts
          if (is_fit_pls) {
            ith_pred_results$ncomp[kk] <- i_k_pred$validation$cv_results$ncomp[best_row]
            selected_pls <- ith_pred_results$ncomp[kk]
          }
          if (is_fit_wapls) {
            ith_pred_results$min_ncomp[kk] <- i_k_pred$validation$cv_results$min_component[best_row]
            ith_pred_results$max_ncomp[kk] <- i_k_pred$validation$cv_results$max_component[best_row]
            selected_pls <- i_k_pred$validation$cv_results[best_row, 1:2]
          }
          
          # Store local CV statistics
          ith_pred_results$loc_rmse_cv[kk] <- i_k_pred$validation$cv_results$rmse_cv[best_row]
          ith_pred_results$loc_st_rmse_cv[kk] <- i_k_pred$validation$cv_results$st_rmse_cv[best_row]
        } else {
          # No local CV: use fixed component counts from fit_method
          if (is_fit_pls) {
            ith_pred_results$ncomp[kk] <- fit_method$ncomp
            selected_pls <- ith_pred_results$ncomp[kk]
          }
          if (is_fit_wapls) {
            ith_pred_results$min_ncomp[kk] <- fit_method$min_ncomp
            ith_pred_results$max_ncomp[kk] <- fit_method$max_ncomp
            selected_pls <- c(fit_method$min_ncomp, fit_method$max_ncomp)
          }
        }
        
        
        if (is_nnv) {
          # Leave-nearest-neighbor-out validation: exclude nearest neighbor (or its group)
          # and predict its value using the remaining neighbors
          if (!is.null(group)) {
            out_group <- which(i_k_group == i_k_group[[ith_observation$local_index_nearest]])
          } else {
            out_group <- ith_observation$local_index_nearest
          }
          
          nearest_pred <- fit_and_predict(
            x = i_k_xr[-out_group, ],
            y = i_k_yr[-out_group, , drop = FALSE],
            pred_method = pred_method,
            scale = fit_scale,
            pls_c = selected_pls,
            noise_variance = noise_variance,
            newdata = i_k_xr[ith_observation$local_index_nearest, , drop = FALSE],
            CV = FALSE,
            tune = FALSE,
            range_prediction_limits = control$range_prediction_limits,
            pls_max_iter = fit_method$max_iter,
            pls_tol = fit_method$tol,
            seed = seed,
            algorithm = fit_method$method
          )$prediction
          
          ith_pred_results$y_nearest_pred[kk] <- nearest_pred / kth_weights[1]
        }
      } else {
        # Neighborhood size unchanged: copy previous results
        ith_k_diss <- ith_pred_results$k_diss[kk]
        ith_pred_results[kk, ] <- ith_pred_results[kk - 1, ]
        ith_pred_results$rep[kk] <- 1
        ith_pred_results$k_diss[kk] <- ith_k_diss
      }
    }
    
    if (verbose) {
      if (kk == nrow(ith_pred_results) & i != n_iter) {
        cat("\r", to_erase, "\r")
      }
      
      if (i == n_iter) {
        cat("\n")
      }
      # do not use close() (it prints a new line)
      ## close(pb)
    }
    list(
      results = ith_pred_results,
      additional_results = additional_results
    )
  }
  
  # Reorder results by observation index (foreach may return out of order)
  iteration_order <- sapply(
    pred_obs,
    FUN = function(x) x$results$o_index[1]
  )
  pred_obs <- pred_obs[order(iteration_order, decreasing = FALSE)]
  
  # Combine all iteration results into single table
  results_table <- do.call(
    "rbind", lapply(pred_obs, FUN = function(x) x$results)
  )
  
  # No k adjustment needed when Xu is provided (always the case now)
  fix_k <- 0
  
  if (is_local_diss) {
    # Reconstruct local dissimilarity matrix from iteration results
    diss_xr_xu <- do.call(
      "cbind",
      lapply(iteration_order,
             FUN = function(ii, x, m) {
               idc <- x[[ii]]$additional_results$ith_neig_indices
               d <- x[[ii]]$additional_results$ith_neigh_diss
               m[idc] <- d
               m
             },
             x = pred_obs,
             m = matrix(NA, nrow(Xr), 1)
      )
    )
    class(diss_xr_xu) <- c("local_ortho_diss", "matrix")
    dimnames(diss_xr_xu) <- list(
      paste0("Xr_", seq_len(nrow(diss_xr_xu))),
      paste0(pre_nms_ng, seq_len(ncol(diss_xr_xu)))
    )
    
    # Reconstruct neighbor indices from iteration results
    neighborhoods$neighbors <- do.call(
      "cbind", lapply(iteration_order,
                      FUN = function(ii, x, m) {
                        idc <- x[[ii]]$additional_results$ith_neig_indices
                        m[seq_len(length(idc))] <- idc
                        m
                      },
                      x = pred_obs,
                      m = matrix(NA, max(results_table$k), 1)
      )
    )
  }
  
  
  out <- c(
    if (is.null(Yu)) "yu_obs",
    if (all(is.na(results_table$k_original))) "k_original",
    if (!is_nnv) "y_nearest_pred",
    if (!is_fit_wapls) c("min_ncomp", "max_ncomp"),
    if (!is_fit_pls) "ncomp",
    if (!is_local_cv) c("loc_rmse_cv", "loc_st_rmse_cv"),
    if (!is_local_diss) "loc_ncomp",
    "rep"
  )
  
  results_table <- results_table[, !(names(results_table) %in% out), drop = FALSE]
  
  
  if (inherits(neighbors, "neighbors_diss")) {
    # Split results by dissimilarity threshold
    k_diss_vec <- neighbors$threshold
    results_table <- lapply(
      k_diss_vec,
      FUN = function(x, tbl, sel) tbl[tbl[[sel]] == x, ],
      tbl = results_table,
      sel = "k_diss"
    )
    names(results_table) <- paste0("k_diss_", k_diss_vec)
    
    # Compute percentage of observations bounded by k_range
    k_min <- neighbors$k_min
    k_max <- neighbors$k_max
    p_bounded <- sapply(
      results_table,
      FUN = function(x, k_min, k_max) {
        sum(x$k_original <= k_min | x$k_original >= k_max)
      },
      k_min = k_min,
      k_max = k_max
    )
    col_ks <- data.frame(
      k_diss = k_diss_vec,
      p_bounded = paste0(round(100 * p_bounded / n_xu, 3), "%")
    )
  } else {
    # Split results by k
    k_vec <- neighbors$k - fix_k
    results_table <- lapply(
      k_vec,
      FUN = function(x, tbl, sel) tbl[tbl[[sel]] == x, ],
      tbl = results_table,
      sel = "k"
    )
    names(results_table) <- paste0("k_", k_vec)
    col_ks <- data.frame(k = k_vec)
  }
  
  
  if (is_nnv) {
    # Compute nearest neighbor validation statistics
    nn_stats <- function(x) {
      nn_rmse <- (mean((x$y_nearest - x$y_nearest_pred)^2))^0.5
      nn_st_rmse <- nn_rmse / diff(range(x$y_nearest))
      nn_rsq <- (cor(x$y_nearest, x$y_nearest_pred))^2
      c(nn_rmse = nn_rmse, nn_st_rmse = nn_st_rmse, nn_rsq = nn_rsq)
    }
    loc_nn_res <- do.call("rbind", lapply(results_table, FUN = nn_stats))
    loc_nn_res <- cbind(
      col_ks,
      rmse = loc_nn_res[, "nn_rmse"],
      st_rmse = loc_nn_res[, "nn_st_rmse"],
      r2 = loc_nn_res[, "nn_rsq"]
    )
  } else {
    loc_nn_res <- NULL
  }
  
  if (is_local_cv) {
    # Compute mean local cross-validation statistics
    mean_loc_res <- function(x) {
      mean_loc_rmse <- mean(x$loc_rmse_cv)
      mean_loc_st_rmse <- mean(x$loc_st_rmse_cv)
      c(loc_rmse = mean_loc_rmse, loc_st_rmse = mean_loc_st_rmse)
    }
    loc_res <- do.call("rbind", lapply(results_table, mean_loc_res))
    loc_res <- cbind(
      col_ks,
      rmse = loc_res[, "loc_rmse"],
      st_rmse = loc_res[, "loc_st_rmse"]
    )
  } else {
    loc_res <- NULL
  }
  
  if (!is.null(Yu)) {
    # Add observed values to results and compute prediction statistics
    for (i in seq_along(results_table)) {
      results_table[[i]][[y_output_name]] <- Yu
    }
    yu_stats <- function(x, y_hat, y) {
      y_rmse <- mean((x[[y_hat]] - x[[y]])^2, na.rm = TRUE)^0.5
      y_st_rmse <- y_rmse / diff(range(x[[y_hat]]), na.rm = TRUE)
      y_rsq <- cor(x[[y_hat]], x[[y]], use = "complete.obs")^2
      c(rmse = y_rmse, st_rmse = y_st_rmse, rsq = y_rsq)
    }
    pred_res <- do.call(
      "rbind",
      lapply(results_table, yu_stats, y_hat = y_hat_output_name, y = y_output_name)
    )
    pred_res <- cbind(
      col_ks,
      rmse = pred_res[, "rmse"],
      st_rmse = pred_res[, "st_rmse"],
      r2 = pred_res[, "rsq"]
    )
  } else {
    pred_res <- NULL
  }
  
  diss_method_name <- if (is_local_diss) {
    paste0(class(diss_method)[1], " (locally computed)")
  } else {
    class(diss_method)[1]
  }
  
  if (control$return_dissimilarity) {
    diss_list <- list(
      diss_method = diss_method,  
      diss_usage = diss_usage,  
      diss_xr_xu = diss_xr_xu
    )
    if (has_projection) {
      diss_list$projection <- diss_projection
    }
  } else {
    diss_method_return <- NULL
    if (inherits(diss_method, "diss_method")) {
      diss_method_return <- diss_method
    } 
    diss_list <- list(
      diss_method = diss_method_return,  
      diss_usage = diss_usage
    )
  }
  
  
  colnames(neighborhoods$neighbors) <- paste0(pre_nms_ng, seq_len(ln))
  rownames(neighborhoods$neighbors) <- paste0("k_", seq_len(nrow(neighborhoods$neighbors)))
  
  val_list <- structure(
    list(
      loc_res,
      loc_nn_res,
      pred_res
    ),
    names = c(
      "local_cross_validation",
      "nearest_neighbor_validation",
      val_summary_name
    )
  )
  
  results_list <- list(
    control = control,
    diss_usage = diss_usage,
    fit_method = fit_method,
    spike = ifelse(is.null(spike), FALSE, TRUE), 
    dissimilarities = diss_list,
    Xu_neighbors = list(
      neighbors = neighborhoods$neighbors,
      neighbors_diss = neighborhoods$neighbors_diss
    ),
    n_predictions = n_xu,
    gh = neighborhoods$gh,
    validation_results = val_list,
    results = results_table,
    seed = seed
  )
  
  attr(results_list, "call") <- f_call
  class(results_list) <- c("mbl", "list")
  
  results_list
}

#' @aliases mbl
#' @export
plot.mbl <- function(x,
                     what = c("validation", "gh"),
                     metric = "rmse",
                     ncomp = c(1, 2),
                     ...) {
  # Validate arguments
  what <- match.arg(what, c("validation", "gh"), several.ok = TRUE)
  metric <- match.arg(metric, c("rmse", "st_rmse", "r2"))
  
  if (!inherits(x, "mbl")) {
    stop("'x' must be an object of class 'mbl'.", call. = FALSE)
  }
  
  # Save and restore graphics parameters
  
  opar <- par("mfrow", "mar")
  on.exit(par(opar))
  
  # Set up multi-panel if needed
  if (length(what) == 2 && !is.null(x$gh)) {
    par(mfrow = c(1, 2))
  }
  
  
  # Process plot arguments
  plot_args <- .process_plot_args(list(...))
  
  # Plot validation results
  if ("validation" %in% what) {
    .plot_validation(x, metric, plot_args)
  }
  
  # Plot GH scores
  if ("gh" %in% what) {
    if (is.null(x$gh)) {
      message("GH distance not available in this object.")
    } else {
      .plot_gh_scores(x, ncomp, plot_args)
    }
  }
  
  mtext(plot_args$main, outer = TRUE, cex = 2, line = -2)
  invisible(x)
}


# =============================================================================
# Internal plotting helpers
# =============================================================================

#' Process plot arguments with defaults
#' @keywords internal
.process_plot_args <- function(dots) {
  args <- dots
  
  # Extract or set main title
  if ("main" %in% names(args)) {
    main <- args$main
    args$main <- NULL
  } else {
    main <- "Memory-based learning results"
  }
  
  # Set defaults
  if (!"col.axis" %in% names(args)) {
    args$col.axis <- grey(0.3)
  }
  if (!"pch" %in% names(args)) {
    args$pch <- 16
  }
  
  # Remove arguments that are set internally
  internal_args <- c("col", "xlab", "ylab", "type", "ylim", "xlim")
  args <- args[!names(args) %in% internal_args]
  
  args$main <- main
  args
}


#' Plot validation results
#' @keywords internal
.plot_validation <- function(object, metric, plot_args) {
  # Collect validation results
  val_data <- list()
  colors <- character()
  
  if (!is.null(object$validation_results$nearest_neighbor_validation)) {
    val_data$NNv <- cbind(
      object$validation_results$nearest_neighbor_validation,
      val = "NNv"
    )
    colors <- c(colors, "dodgerblue")
  }
  
  if (!is.null(object$validation_results$local_cross_validation)) {
    val_data$local_cv <- cbind(
      object$validation_results$local_cross_validation,
      r2 = NA,
      val = "local_cv"
    )
    colors <- c(colors, "green4")
  }
  
  if (!is.null(object$validation_results$Yu_prediction_statistics)) {
    val_data$Yu <- cbind(
      object$validation_results$Yu_prediction_statistics,
      val = "Yu prediction"
    )
    colors <- c(colors, "red")
  }
  
  if (!is.null(object$validation_results$Yr_fitted_statistics)) {
    val_data$Yr <- cbind(
      object$validation_results$Yr_fitted_statistics,
      val = "Yr fitted"
    )
    colors <- c(colors, "orange")
  }
  
  if (length(val_data) == 0) {
    par(mfrow = c(1, 1))
    message("No validation results to plot.")
    return(invisible(NULL))
  }
  
  tpl <- do.call(rbind, val_data)
  col_names <- colnames(tpl)
  
  # Determine ID variable
  id_var <- if ("k" %in% col_names) "k" else "k_diss"
  
  # Select columns for metric
  keep_cols <- c(id_var, metric, "val")
  keep_cols <- keep_cols[keep_cols %in% col_names]
  tpl_subset <- data.frame(tpl)[, keep_cols, drop = FALSE]
  
  # Reshape to wide format
  to_plot <- reshape(
    tpl_subset,
    timevar = "val",
    idvar = id_var,
    direction = "wide"
  )
  colnames(to_plot) <- gsub("[.]", " ", colnames(to_plot))
  
  # Remove r2 for local_cv (not available)
  if (metric == "r2") {
    drop_col <- grepl("local_cv", colnames(to_plot))
    to_plot <- to_plot[, !drop_col, drop = FALSE]
    colors <- colors[colors != "green4"]
  }
  
  # Plot
  y_range <- range(to_plot[, -1], na.rm = TRUE)
  y_range[2] <- y_range[2] * 1.1
  
  do.call("matplot", c(
    list(
      x = to_plot[, 1],
      y = to_plot[, -1, drop = FALSE],
      type = "b",
      xlab = id_var,
      ylab = metric,
      ylim = y_range,
      col = colors
    ),
    plot_args[!names(plot_args) %in% "main"]
  ))
  
  grid(col = rgb(0.3, 0.3, 0.3, 0.1), lty = 1, lwd = 1)
  mtext("Validation results", col = grey(0.3))
  
  legend(
    "topright",
    legend = colnames(to_plot)[-1],
    col = colors,
    pch = plot_args$pch,
    box.lty = 0,
    bg = NA
  )
  
  invisible(NULL)
}


#' Plot GH distance scores
#' @keywords internal
.plot_gh_scores <- function(object, ncomp, plot_args) {
  xr_scores <- object$gh$projection$scores
  n_components <- object$gh$projection$ncomp
  
  # Transform to Mahalanobis space if multicomponent
  if (n_components > 1) {
    xr_scores <- euclid_to_mahal(xr_scores)
    xr_scores <- sweep(xr_scores, 2, colMeans(xr_scores), "-")
  }
  
  # Split Xr and Xu scores
  xu_idx <- grep("Xu_", rownames(xr_scores))
  xr_idx <- grep("Xr_", rownames(xr_scores))
  
  xu_scores <- xr_scores[xu_idx, , drop = FALSE]
  xr_scores <- xr_scores[xr_idx, , drop = FALSE]
  
  # Colors
  xr_col <- rgb(0, 0, 0.4, 0.5)
  xu_col <- rgb(1, 0, 0, 0.5)
  
  if (n_components == 1) {
    .plot_gh_1d(xr_scores, xu_scores, xr_col, xu_col, plot_args)
  } else {
    .plot_gh_2d(xr_scores, xu_scores, ncomp, xr_col, xu_col, plot_args)
  }
  
  invisible(NULL)
}


#' Plot 1D GH scores
#' @keywords internal
.plot_gh_1d <- function(xr_scores, xu_scores, xr_col, xu_col, plot_args) {
  scores_combined <- c(xr_scores[, 1], xu_scores[, 1])
  set_labels <- c(rep("Xr", nrow(xr_scores)), rep("Xu", nrow(xu_scores)))
  
  df <- data.frame(
    index = seq_along(scores_combined),
    score = scores_combined,
    set = set_labels
  )
  df <- df[order(df$score), ]
  df$index <- seq_len(nrow(df))
  
  rng <- range(scores_combined)
  rng <- 1.2 * c(-max(abs(rng)), max(abs(rng)))
  
  do.call("plot", c(
    list(
      x = df$index[df$set == "Xr"],
      y = df$score[df$set == "Xr"],
      ylim = rng,
      col = xr_col,
      xlab = "Ordered PLS values",
      ylab = "PLS 1"
    ),
    plot_args[!names(plot_args) %in% "main"]
  ))
  
  points(
    df$index[df$set == "Xu"],
    df$score[df$set == "Xu"],
    col = xu_col,
    pch = plot_args$pch
  )
  
  mtext("Partial least squares scores", col = grey(0.3))
  grid(col = rgb(0.3, 0.3, 0.3, 0.1), lty = 1, lwd = 1)
  legend("topleft", legend = c("Xr", "Xu"), col = c(xr_col, xu_col),
         pch = plot_args$pch, cex = 0.8, box.lty = 0, bg = NA)
}


#' Plot 2D GH scores
#' @keywords internal
.plot_gh_2d <- function(xr_scores, xu_scores, ncomp, xr_col, xu_col, plot_args) {
  rng <- 1.2 * range(xr_scores[, ncomp], xu_scores[, ncomp])
  
  xlab <- paste0("PLS ", ncomp[1], " (standardized)")
  ylab <- paste0("PLS ", ncomp[2], " (standardized)")
  
  do.call("plot", c(
    list(
      x = xr_scores[, ncomp[1]],
      y = xr_scores[, ncomp[2]],
      xlab = xlab,
      ylab = ylab,
      xlim = rng,
      ylim = rng,
      col = xr_col
    ),
    plot_args[!names(plot_args) %in% "main"]
  ))
  
  points(xu_scores[, ncomp, drop = FALSE], col = xu_col, pch = plot_args$pch)
  
  mtext("Partial least squares (Mahalanobis space)", col = grey(0.3))
  grid(col = rgb(0.3, 0.3, 0.3, 0.1), lty = 1, lwd = 1)
  
  legend("topright", legend = c("Xr", "Xu"), col = c(xr_col, xu_col),
         pch = plot_args$pch, cex = 0.8, box.lty = 0, bg = NA)
  
  # Draw reference circles
  max_score <- ceiling(max(abs(xr_scores[, ncomp])))
  for (r in seq_len(max_score)) {
    theta <- seq(0, 2 * pi, length.out = 101)
    lines(r * cos(theta), r * sin(theta),
          col = rgb(0.3, 0.3, 0.3, 0.3), lty = 1, lwd = 0.5)
  }
}
