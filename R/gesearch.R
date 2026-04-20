#' @title Evolutionary sample search for context-specific calibrations
#' @name gesearch
#' @aliases gesearch gesearch.default gesearch.formula
#' @aliases predict.gesearch plot.gesearch
#'
#' @description
#' Implements an evolutionary search algorithm that selects a subset from large
#' reference datasets (e.g., spectral libraries) to build context-specific
#' calibrations. The algorithm iteratively removes weak or non-informative
#' samples based on prediction error, spectral reconstruction error, or
#' dissimilarity criteria. This implementation is based on the methods proposed 
#' in Ramirez-Lopez et al. (2026a).
#' @usage
#'
#' \method{gesearch}{default}(Xr, Yr, Xu, Yu = NULL, Yu_lims = NULL,
#'          k, b, retain = 0.95, target_size = k,
#'          fit_method = fit_pls(ncomp = 10),
#'          optimization = "reconstruction",
#'          group = NULL, control = gesearch_control(),
#'          intermediate_models = FALSE,
#'          verbose = TRUE, seed = NULL, pchunks = 1L, ...)
#'
#' \method{gesearch}{formula}(formula, train, test, k, b, target_size, fit_method,
#'          ..., na_action = na.pass)
#'
#' \method{predict}{gesearch}(object, newdata, type = "response",
#'          what = c("final", "all_generations"), ...)
#'
#' \method{plot}{gesearch}(x, which = c("weakness", "removed"), ...)
#'
#' @param Xr A numeric matrix of predictor variables for the reference data
#'   (observations in rows, variables in columns).
#' @param Yr A numeric vector or single-column matrix of response values
#'   corresponding to \code{Xr}. Only one response variable is supported.
#' @param Xu A numeric matrix of predictor variables for target observations
#'   (same structure as \code{Xr}).
#' @param Yu An optional numeric vector or single-column matrix of response 
#'   values for \code{Xu}. Required when \code{optimization} includes 
#'   \code{"response"}. Default is \code{NULL}.
#' @param Yu_lims A numeric vector of length 2 specifying expected response
#'   limits for the target population. Used with \code{optimization = "range"}.
#' @param k An integer specifying the number of samples in each resampling
#'   subset (gene size).
#' @param b An integer specifying the target average number of times each
#'   training sample is evaluated per iteration. Higher values (e.g., >40)
#'   produce more stable results but increase computation time.
#' @param retain A numeric value in (0, 1] specifying the proportion of samples
#'   retained per iteration. Default is 0.95. Values >0.9 are recommended for
#'   stability. See \code{\link{gesearch_control}} for retention strategy.
#' @param target_size An integer specifying the target number of selected
#'   samples (gene pool size). Must be >= \code{k}. Default is \code{k}.
#' @param fit_method A fit method object created with \code{\link{fit_pls}}.
#'   Specifies the regression model and scaling used during the search.
#'   Currently only \code{fit_pls()} is supported.
#' @param optimization A character vector specifying optimization criteria:
#'   \itemize{
#'     \item \code{"reconstruction"}: (default) Retains samples based on
#'       spectral reconstruction error of \code{Xu} in PLS space.
#'     \item \code{"response"}: Retains samples based on RMSE of predicting
#'       \code{Yu}. Requires \code{Yu}.
#'     \item \code{"similarity"}: Retains samples based on Mahalanobis distance
#'       between \code{Xu} and training samples in PLS score space.
#'     \item \code{"range"}: Removes samples producing predictions outside
#'       \code{Yu_lims}.
#'   }
#'   Multiple criteria can be combined, e.g.,
#'   \code{c("reconstruction", "similarity")}.
#' @param group An optional factor assigning group labels to training
#'   observations. Used for leave-group-out cross-validation to avoid
#'   pseudo-replication.
#' @param control A list created with \code{\link{gesearch_control}} containing
#'   additional algorithm parameters.
#' @param intermediate_models A logical indicating whether to store models for
#'   each intermediate generation. Default is \code{FALSE}.
#' @param verbose A logical indicating whether to print progress information.
#'   Default is \code{TRUE}.
#' @param seed An integer for random number generation to ensure
#'   reproducibility. Default is \code{NULL}.
#' @param pchunks An integer specifying the chunk size used for memory-efficient
#' parallel processing. Larger values divide the workload into smaller pieces,
#' which can help reduce memory pressure. Default is 1L.
#' @param formula A \code{\link[stats]{formula}} defining the model.
#' @param train A data.frame containing training data with model variables.
#' @param test A data.frame containing test data with model variables.
#' @param na_action A function for handling missing values in training data.
#'   Default is \code{\link[stats]{na.pass}}.
#' @param object A fitted \code{gesearch} object (for \code{predict}).
#' @param newdata A matrix or data.frame of new observations. For formula-fitted
#'   models, a data.frame containing all predictor variables is accepted. For
#'   non-formula models, a matrix is required.
#' @param type A character string specifying the prediction type. Currently only
#'   \code{"response"} is supported.
#' @param what A character string specifying which models to use for prediction:
#'   \code{"final"} (default) for predictions from final models only, or
#'   \code{"all_generations"} for predictions from all intermediate generations
#'   plus the final models.
#' @param x A \code{gesearch} object (for \code{plot}).
#' @param which Character string specifying what to plot: 
#'   \code{"weakness"} (maximum weakness scores per generation) or 
#'   \code{"removed"} (cumulative samples removed).
#' @param ... Additional arguments passed to methods.
#'
#' @details
#' The \code{gesearch} algorithm requires a large reference dataset (\code{Xr})
#' where the sample search is conducted, target observations (\code{Xu}), and
#' three tuning parameters: \code{k}, \code{b}, and \code{retain}.
#'
#' The target observations (\code{Xu}) should represent the population of
#' interest. These may be selected via algorithms like Kennard-Stone when
#' response values are unavailable.
#'
#' The algorithm iteratively removes weak samples from \code{Xr} based on:
#' \itemize{
#'   \item Increased RMSE when predicting \code{Yu}
#'   \item Increased PLS reconstruction error on \code{Xu}
#'   \item Increased dissimilarity to \code{Xu} in PLS space
#' }
#'
#' A resampling scheme identifies samples that consistently appear in
#' high-error subsets. These are labeled weak and removed. The process
#' continues until approximately \code{target_size} samples remain.
#'
#' The `gesearch()` function also returns a final model fitted on the selected
#' samples, which can be used for prediction. This model is internally validated
#' by cross-validation using only the selected samples from the training/reference
#' set. If `Yu` is available, a model fitted only on the selected reference samples
#' is first used to predict the target samples. The final model is then refitted
#' using both the selected reference samples and the target samples used to guide
#' the search, provided that response values are available for those target samples.
#'
#' ## Parameter guidance
#' \itemize{
#'   \item \code{k}: Number of samples per resampling subset. See Lobsey et
#'     al. (2017) for guidance.
#'   \item \code{b}: Resampling intensity. Higher values increase stability
#'     but computational cost.
#'   \item \code{retain}: Proportion retained per iteration. Values >0.9
#'     recommended.
#' }
#' 
#' ## Prediction
#' The \code{predict} method generates predictions from a fitted 
#' \code{gesearch} object. If the model was fitted with a formula, 
#' \code{newdata} is validated and transformed to the appropriate model matrix.
#' 
#' When \code{what = "all_generations"}, the return value is a named list with
#' one element per generation, where each element contains a prediction 
#' matrix. This option requires \code{intermediate_models = TRUE} during 
#' fitting.
#'
#' @return
#' \strong{For \code{gesearch}:} A list of class \code{"gesearch"} containing:
#' \itemize{
#'   \item \code{x_local}: Matrix of predictors for selected samples.
#'   \item \code{y_local}: Vector of responses for selected samples.
#'   \item \code{indices}: Indices of selected samples from original training set.
#'   \item \code{complete_iter}: Number of completed iterations.
#'   \item \code{iter_weakness}: List with iteration-level weakness statistics.
#'   \item \code{samples}: List of sample indices retained at each iteration.
#'   \item \code{n_removed}: data.frame of samples removed per iteration.
#'   \item \code{control}: Copy of control parameters.
#'   \item \code{fit_method}: Fit constructor from \code{fit_method}.
#'   \item \code{validation_results}: Cross-validation in the training only set 
#'   validation on the test set using models built only with the samples found.
#'   \item \code{final_models}: Final PLS model containing coefficients, loadings, 
#'     scores, VIP, and selectivity ratios.
#'   \item \code{intermediate_models}: List of models per generation (if
#'     \code{intermediate_models = TRUE}).
#'   \item \code{seed}: RNG seed used.
#' }
#' 
#' \strong{For \code{predict.gesearch}:}
#' \itemize{
#'   \item If \code{what = "final"}: a prediction matrix with 
#'     \code{nrow(newdata)} rows and one column per PLS component.
#'   \item If \code{what = "all_generations"}: a named list of generations, 
#'     where each generation contains a prediction matrix as above.
#' }
#'
#' @author
#' \href{https://orcid.org/0000-0002-5369-5120}{Leonardo Ramirez-Lopez},
#' Claudio Orellano,
#' \href{https://orcid.org/0000-0001-5416-4520}{Craig Lobsey},
#' \href{https://orcid.org/0000-0003-1540-4748}{Raphael Viscarra Rossel}
#'
#' @references
#' Lobsey, C.R., Viscarra Rossel, R.A., Roudier, P., Hedley, C.B. 2017.
#' rs-local data-mines information from spectral libraries to improve local
#' calibrations. European Journal of Soil Science 68:840-852.
#'
#' Kennard, R.W., Stone, L.A. 1969. Computer aided design of experiments.
#' Technometrics 11:137-148.
#'
#' Rajalahti, T., Arneberg, R., Berven, F.S., Myhr, K.M., Ulvik, R.J.,
#' Kvalheim, O.M. 2009. Biomarker discovery in mass spectral profiles by means
#' of selectivity ratio plot. Chemometrics and Intelligent Laboratory Systems
#' 95:35-48.
#' 
#' Ramirez-Lopez, L., Viscarra Rossel, R., Behrens, T., Orellano, C.,
#' Perez-Fernandez, E., Kooijman, L., Wadoux, A. M. J.-C., Breure, T.,
#' Summerauer, L., Safanelli, J. L., & Plans, M. (2026a). When spectral
#' libraries are too complex to search: Evolutionary subset selection for
#' domain-adaptive calibration. \emph{Analytica Chimica Acta}, under review.
#'
#' @seealso
#' \code{\link{fit_pls}}, \code{\link{gesearch_control}}, \code{\link{mbl}}
#'
#' @examples
#' \dontrun{
#' library(prospectr)
#' data(NIRsoil)
#'
#' # Preprocess
#' sg_det <- savitzkyGolay(
#'   detrend(NIRsoil$spc, wav = as.numeric(colnames(NIRsoil$spc))),
#'   m = 1, p = 1, w = 7
#' )
#' NIRsoil$spc_pr <- sg_det
#'
#' # Split data
#' train_x <- NIRsoil$spc_pr[NIRsoil$train == 1 & !is.na(NIRsoil$Ciso), ]
#' train_y <- NIRsoil$Ciso[NIRsoil$train == 1 & !is.na(NIRsoil$Ciso)]
#' test_x <- NIRsoil$spc_pr[NIRsoil$train == 0 & !is.na(NIRsoil$Ciso), ]
#' test_y <- NIRsoil$Ciso[NIRsoil$train == 0 & !is.na(NIRsoil$Ciso)]
#'
#' # Basic search with reconstruction and similarity optimizations
#' gs <- gesearch(
#'   Xr = train_x, Yr = train_y,
#'   Xu = test_x, Yu = test_y,
#'   k = 50, b = 100, retain = 0.97,
#'   target_size = 200,
#'   fit_method = fit_pls(ncomp = 15, method = "mpls"),
#'   optimization = c("reconstruction", "similarity"),
#'   control = gesearch_control(retain_by = "probability"),
#'   seed = 42
#' )
#'
#' # Predict
#' preds <- predict(gs, test_x)
#'
#' # Plot progress
#' plot(gs)
#' plot(gs, which = "removed")
#'
#' # With reconstruction and response optimization (requires Yu)
#' gs_response <- gesearch(
#'   Xr = train_x, Yr = train_y,
#'   Xu = test_x, Yu = test_y,
#'   k = 50, b = 100, retain = 0.97,
#'   target_size = 200,
#'   fit_method = fit_pls(ncomp = 15),
#'   optimization = c("reconstruction", "response"),
#'   seed = 42
#' )
#'
#' # Parallel processing
#' library(doParallel)
#' n_cores <- min(2, parallel::detectCores() - 1)
#' cl <- makeCluster(n_cores)
#' registerDoParallel(cl)
#'
#' gs_parallel <- gesearch(
#'   Xr = train_x, Yr = train_y,
#'   Xu = test_x,
#'   k = 50, b = 100, retain = 0.97,
#'   target_size = 200,
#'   fit_method = fit_pls(ncomp = 15),
#'   pchunks = 3,
#'   seed = 42
#' )
#'
#' stopCluster(cl)
#' registerDoSEQ()
#' }
#'
#' @importFrom stats as.formula update na.pass model.frame model.matrix
#' @importFrom stats model.extract terms .MFclass delete.response
#' @export gesearch
#' 
NULL
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

#' @aliases gesearch
"gesearch" <-
  function(...) {
    UseMethod("gesearch")
  }


#' @aliases gesearch
#' @export
#'
## add an argument to initualize with a given amount of observations, e.g. if 
## the library is 30.000 you could initialize with 15.000 observations and to 
## compensate you can increase the number of resampling iterations
## perhaps a number of initial resampling iterations can also work say we 
## start with 10.000 for the first loop and then it goes to what b specifies
gesearch.default <- function(
    Xr,
    Yr,
    Xu,
    Yu = NULL,
    Yu_lims = NULL,
    k,
    b,
    retain = 0.95,
    target_size = k,
    fit_method = fit_pls(ncomp = 10),
    optimization = "reconstruction",
    group = NULL,
    control = gesearch_control(),
    intermediate_models = FALSE,
    verbose = TRUE,
    seed = NULL,
    pchunks = 1L,
    ...
) {
  crossover <- TRUE

  dots <- list(...)
  if ("formula" %in% names(dots)) {
    if (!inherits(dots$formula, "formula")) {
      stop("The 'formula' argument must be a formula object.", call. = FALSE)
    }
  }
  
  
  # --- Input validation ---
  if (!inherits(Xr, c("matrix", "data.frame", "formula"))) {
    stop(
      "You are attempting to pass a ", class(Xr)[1L], " as 'Xr' or 'formula'.\n",
      "Expected a numeric matrix for gesearch.default() or a formula for gesearch.formula().",
      call. = FALSE
    )
  }
  
  Yr <- as.matrix(Yr)

  # --- control validation ---
  if (!inherits(control, "gesearch_control")) {
    stop("'control' must be created by gesearch_control()", call. = FALSE)
  }
  
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    old_blas_threads <- blas_get_num_procs()
    if (old_blas_threads != control$blas_threads) {
      blas_set_num_threads(control$blas_threads)
      on.exit(blas_set_num_threads(old_blas_threads), add = TRUE)
    }
  } else if (Sys.info()["sysname"] == "Linux" && control$blas_threads == 1L) {
    message("Tip: Install 'RhpcBLASctl' for optimal performance on Linux:\n",
            "  install.packages('RhpcBLASctl')")
  }
  
  if (ncol(Yr) == 1L) {
    if (nrow(Xr) != length(Yr)) {
      stop("nrow(Xr) must equal length(Yr)", call. = FALSE)
    }
  } else {
    if (is.null(colnames(Yr))) {
      stop("Yr must have column names when it has multiple columns", call. = FALSE)
    }
    if (nrow(Xr) != nrow(Yr)) {
      stop("nrow(Xr) must equal nrow(Yr)", call. = FALSE)
    }
  }
  
  if (target_size > nrow(Xr)) {
    stop("'target_size' cannot be greater than nrow(Xr)", call. = FALSE)
  }
  
  if (!is.logical(verbose)) {
    stop("'verbose' must be logical", call. = FALSE)
  }
  
  if (nrow(Xu) < 3L) {
    stop("Xu must have at least 3 rows", call. = FALSE)
  }
  
  if (!is.null(Yu)) {
    Yu <- as.matrix(Yu)
    if (ncol(Yu) == 1L) {
      if (nrow(Xu) != length(Yu)) {
        stop("nrow(Xu) must equal length(Yu)", call. = FALSE)
      }
    } else {
      if (!is.null(Yu_lims)) {
        stop("'Yu_lims' only supported for single-column Yu", call. = FALSE)
      }
      if (ncol(Yr) != ncol(Yu)) {
        stop("ncol(Yr) must equal ncol(Yu)", call. = FALSE)
      }
      if (nrow(Xu) != nrow(Yu)) {
        stop("nrow(Xu) must equal nrow(Yu)", call. = FALSE)
      }
      if (is.null(colnames(Yu))) {
        stop("Yu must have column names when it has multiple columns", call. = FALSE)
      }
      if (!all(colnames(Yr) == colnames(Yu))) {
        stop("Column names of Yr and Yu must match", call. = FALSE)
      }
    }
  }
  
  # Yu required with non-missing values for response optimization
  if ("response" %in% optimization) {
    if (is.null(Yu) || !any(!is.na(Yu))) {
      stop("'Yu' with non-missing values required when optimization includes 'response'", 
           call. = FALSE)
    }
  }
  
  # No missing values allowed in Xr, Yr, or Xu
  if (anyNA(Xr)) {
    stop("'Xr' contains missing values.", call. = FALSE)
  }
  if (anyNA(Yr)) {
    stop("'Yr' contains missing values.", call. = FALSE)
  }
  if (anyNA(Xu)) {
    stop("'Xu' contains missing values.", call. = FALSE)
  }
  
  # No infinite values in any input
  if (any(is.infinite(Xr)) || any(is.infinite(Yr)) || 
      any(is.infinite(Xu)) || (!is.null(Yu) && any(is.infinite(Yu)))) {
    stop("Infinite values detected in input data", call. = FALSE)
  }
  
  # Xr and Xu must have same number of variables
  if (ncol(Xr) != ncol(Xu)) {
    stop("ncol(Xr) must equal ncol(Xu)", call. = FALSE)
  }
  
  # Validate k (gene size)
  k <- as.integer(round(k))
  if (k >= nrow(Xr)) {
    stop("'k' must be less than nrow(Xr)", call. = FALSE)
  }
  
  # Validate target_size
  target_size <- as.integer(round(target_size))
  if (target_size < k) {
    stop("'target_size' must be >= k", call. = FALSE)
  }
  
  # Validate retain proportion
  if (retain <= 0 || retain > 1) {
    stop("'retain' must be in (0, 1]", call. = FALSE)
  }
  
  # Validate optimization options
  # Note: "fresponse" is undocumented (internal/experimental)
  allowed <- c("response", "reconstruction", "similarity", "range", "fresponse")
  bad <- setdiff(optimization, allowed)
  if (length(bad) > 0L) {
    stop(
      sprintf(
        "Invalid 'optimization' option(s): %s. Allowed: %s",
        paste(shQuote(bad), collapse = ", "),
        paste(shQuote(allowed), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  # Validate fit_method
  # TODO: Currently only fit_pls() is supported; fit_wapls(), fit_gpr() may be added later
  if (!inherits(fit_method, "fit_method")) {
    stop("'fit_method' must be a fit method object (e.g., fit_pls())", call. = FALSE)
  }
  
  if (!inherits(fit_method, "fit_pls")) {
    stop("Only fit_pls() is supported for the moment", call. = FALSE)
  }
  
  # Tuning not applied during reconstruction optimization
  if ("reconstruction" %in% optimization && control$tune) {
    warning(
      "PLS components are not tuned when optimization includes 'reconstruction'; ",
      "using fixed ncomp from 'fit_method'. The 'tune' option was ignored.",
      call. = FALSE
    )
  }
  
  
  if (control$tune) {
    min_samples <- floor(min(k, dim(Xr)) * control$p) - 1
    min_cv_samples <- floor(min(k, dim(Xr)) * (1 - control$p))
    if (min_cv_samples < 3) {
      stop(paste0(
        "Local cross-validation requires at least 3 observations in ",
        "the hold-out set, the current cross-validation parameters ",
        "leave less than 3 observations for one or more resampling ",
        "iterations."
      ))
    }
  } else {
    min_samples <- floor(min(k, dim(Xr))) - 1
  }
  
  # Validate ncomp vs available samples
  # TODO: Add fit_wapls support when implemented
  if (inherits(fit_method, "fit_pls")) {
    max_ncomp <- fit_method$ncomp
  } else if (inherits(fit_method, "fit_wapls")) {
    max_ncomp <- fit_method$max_ncomp
  }
  
  if (min_samples < max_ncomp) {
    stop(
      "More PLS components than available observations in some resampling iterations. ",
      "If tuning is enabled, some observations are held out for validation.",
      call. = FALSE
    )
  }
  
  call_f <- match.call()
  
  Xr <- as.matrix(Xr)
  Xu <- as.matrix(Xu)
  
  # Proportion to remove each iteration
  r <- 1 - round(retain, .Machine$sizeof.longdouble)
  
  # Index into the spectral library (SL); we operate on indices, not copies
  sl_idx <- seq_len(nrow(Xr))
  
  # Initialize K as full SL
  k_idx <- sl_idx
  
  # Iteration counter
  outer_idx <- 0L
  
  # Sample ranking vectors (by RMSE)
  U <- V <- rep(0, length(sl_idx))
  
  # Iteration tracking
  s_to_drop <- NULL
  cuts <- NULL
  max_weakness_score <- NULL
  pp <- r
  
  # # Progress indicators
  # cat_progress <- cat_iter(c("\\", "|", "/", "-"))
  # cat_progress2 <- cat_iter(c("\\", "|", "/", "-"))
  
  # Define once outside the loop
  if (verbose) {
    # spinner <- c("⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏")
    spinner <- c(">  ", ">> ", ">>>", "   ")
    # spinner <- c("|", "/", "-", "\\")
    # spinner <- c(".", "o", "O", "o")
    # spinner <- c("[=   ]", "[ =  ]", "[  = ]", "[   =]", "[  = ]", "[ =  ]")
    # spinner <- c(".  ", ".. ", "...", "   ")
    # spinner <- c(">  ", ">> ", ">>>", "   ")
    # spinner <- c("[-]", "[\\]", "[|]", "[/]")
    
    n_spin <- length(spinner)
    prev_msg_len <- 0L
  }
  
  
  pool_sizes_log <- integer(0L)
  complete_iterations <- 0L
  selected <- list()
  # Step 2-6: Main evolutionary loop (see Step 7 for termination)
  while (length(k_idx) > (target_size * (1 + pp)) & sum(!is.na(sl_idx)) > target_size) {
    # Step 1: Update K to contain only non-dropped samples
    k_idx <- sl_idx[!is.na(sl_idx)]
    selected <- c(selected, list(k_idx))
    
    outer_idx <- outer_idx + 1L
    
    # Number of samples to remove this iteration
    cull_quantity <- round(length(k_idx) * round(r, .Machine$sizeof.longdouble))
    
    # Step 2-4: Resampling iterations (B times)
    B <- round(length(k_idx) * b / k)
    
    # Progress output
    # if (verbose) {
    #   pool_label <- if (outer_idx == 1L) "Initial set of" else "Active"
    #   
    #   msg <- paste0(
    #     "Generation ", outer_idx, ":\t",
    #     "\033[34m", cat_progress(), " ", pool_label, " genes (samples): ", length(k_idx), "\t\033[39m",
    #     "\033[34m", cat_progress2(), " Individuals to evaluate: ", B, "\033[39m"
    #   )
    #   
    #   # Pad to fixed width and return cursor to start
    #   cat(sprintf("\r%-120s", msg))
    #   flush.console()
    # }
    if (verbose) {
      spin <- spinner[(outer_idx %% n_spin) + 1L]
      pool_label <- if (outer_idx == 1L) "Initial" else "Active"
      
      # Build message without ANSI codes for length calculation
      msg_plain <- sprintf(
        "Generation %d: %s genes (samples): %d %s | Individuals: %d %s",
        outer_idx, pool_label, length(k_idx), spin, B, spin
      )
      
      # Build message with ANSI codes for display
      msg <- sprintf(
        "\rGeneration %d: \033[34m%s genes (samples): %d %s\033[39m | \033[34mIndividuals: %d %s\033[39m",
        outer_idx, pool_label, length(k_idx), spin, B, spin
      )
      
      # Pad with spaces to overwrite previous longer message
      padding <- max(0L, prev_msg_len - nchar(msg_plain))
      cat(msg, strrep(" ", padding), sep = "")
      flush.console()
      
      prev_msg_len <- nchar(msg_plain)
    }
    
    # Step 2: Sample training subsets of size k without replacement
    if (!is.null(seed)) {
      set.seed(seed + B)
    }
    
    if (outer_idx == 1L) {
      ismpl <- replicate(
        n = B,
        expr = sample(x = k_idx, size = k, replace = FALSE)
      )
    } else {
      # Replace genes dropped in last iteration with NA
      ismpl[ismpl %in% which(new_nas)] <- NA
      
      if (crossover) {
        ismpl <- replicate(
          n = B,
          expr = combine_individuals(ismpl, k, k_idx)
        )
      } else {
        ismpl <- replicate(
          n = B,
          expr = sample(x = k_idx, size = k, replace = FALSE)
        )
      }
    }
    
    # Iterator for memory-efficient parallel chunking:
    # Takes only needed slices from the matrix rather than loading
    # the full matrix into memory for each parallel worker
    if (control$allow_parallel && getDoParRegistered()) {
      leftc <- B %% pchunks
      if (leftc > 0L) {
        req_iter <- floor(B / pchunks) + leftc
      } else {
        req_iter <- B / pchunks
      }
    } else {
      req_iter <- B
      pchunks <- 1L
    }
    
    itersubs <- ithrssubsets(
      x = Xr,
      y = Yr,
      group = group,
      indx = ismpl,
      chunksize = pchunks
    )
    
    # Step 3: Calibrate PLS model using selected samples
    # Step 4: Validate on target samples
    # Step 5: Repeat steps 2-4 B times
    
    # Validate Yu_lims for range optimization
    if ("range" %in% optimization) {
      if (is.null(Yu_lims)) {
        stop("'Yu_lims' is required for 'range' optimization", call. = FALSE)
      }
      if (length(Yu_lims) != 2L) {
        stop("'Yu_lims' must be a numeric vector of length 2", call. = FALSE)
      }
    }
    
    # Set response limits for range filtering
    if (!is.null(Yu_lims) && "range" %in% optimization) {
      user_min_y <- min(Yu_lims)
      user_max_y <- max(Yu_lims)
    } else {
      user_min_y <- NULL
      user_max_y <- NULL
      optimization <- optimization[optimization != "range"]
    }
    
    results_df <- biter(
      itersubs = itersubs,
      Xu = Xu,
      Yu = Yu,
      y_lim_left = user_min_y,
      y_lim_right = user_max_y,
      iter_sequence = seq_len(req_iter),
      n = nrow(ismpl),
      optimization = optimization,
      ncomp = fit_method$ncomp,
      tune = control$tune,
      p = control$p,
      number = control$number,
      scale = fit_method$scale,
      max_iter = 1L,
      tol = 1e-6,
      seed = seed,
      algorithm = fit_method$method,
      allow_parallel = control$allow_parallel,
      pchunks = pchunks
    )
    
    # Step 6: Rank samples by weakness and identify candidates for removal
    drop_results <- apply(
      results_df$weakness_scores_subset,
      MARGIN = 2L,
      FUN = drop_indices,
      uu = U,
      vv = V,
      r = r,
      obs_indices = sl_idx,
      resampling_indices = results_df$sampleidx,
      retain_by = control$retain_by,
      cull_quantity = cull_quantity,
      percentile_type = control$percentile_type,
      max_ncomp = fit_method$ncomp
    )
    
    # Check for interruptions from drop_indices
    was_interrupted <- vapply(drop_results, function(x) x$interrupted, logical(1L))
    
    if (any(was_interrupted)) {
      # Clean up analysis names for display
      analysis_name <- gsub("rmse_|score_dis|residual_| deviation_from", "", names(drop_results))
      analysis_name <- gsub("_", " ", analysis_name)
      
      # Report which analyses were interrupted
      for (i in which(was_interrupted)) {
        message(analysis_name[i], ": ", drop_results[[i]]$message)
      }
      break
    }
    
    pool_sizes_log <- c(pool_sizes_log, length(k_idx))
    
    # Accumulate cutoff and weakness scores across iterations
    cuts <- rbind(cuts, vapply(drop_results, function(x) x$cutoff, numeric(1L)))
    max_weakness_score <- rbind(
      max_weakness_score, 
      vapply(drop_results, function(x) x$max_weakness_score, numeric(1L))
    )
    
    # Mark dropped samples as NA
    new_nas <- !complete.cases(
      vapply(drop_results, function(x) x$obs_indices, numeric(length(sl_idx)))
    )
    sl_idx[new_nas] <- NA
    
    complete_iterations <- complete_iterations + 1L
    s_to_drop <- c(s_to_drop, sum(is.na(sl_idx)))
    
    # Check for stagnation (gene pool size unchanged)
    if (length(pool_sizes_log) >= control$stagnation_limit) {
      recent_sizes <- tail(pool_sizes_log, control$stagnation_limit)
      if (all(recent_sizes == recent_sizes[1L])) {
        cat(
          "\n\033[31mEarly termination: Active genes remained at ",
          recent_sizes[1L], " for ", control$stagnation_limit,
          " consecutive iterations. Stopping.\033[39m\n",
          sep = ""
        )
        break
      }
    }
  }
  # Finalize k_idx after last iteration
  k_idx <- sl_idx[!is.na(sl_idx)]
  selected <- c(selected, list(k_idx))
  
  # Ensure final selection meets target_size
  
  msizes <- vapply(selected, length, integer(1L))
  k_idx <- selected[[max(which(msizes >= target_size))]]
  selected <- selected[msizes >= target_size]
  
  
  models_per_generation <- vector("list", length(selected) - 1)
  if (intermediate_models) {
    if (verbose) {
      cat("\nFitting intermediate models...")
    }
    for (i in 1:(length(selected) - 1)) {
      ith_selected <- selected[[i]]
      
      # Fit models for this generation
      pls_models <- list()
      validation_results <- list()
      for (j in seq_len(ncol(Yr))) {
        # Fit and validate PLS model for response j
        pls_val <- get_plsr(
          X = Xr,
          Y = Yr[, j, drop = FALSE],
          indices = ith_selected,
          fit_method = fit_method,
          number = control$number,
          group = group[ith_selected],
          retrieve = FALSE,
          tune = TRUE,
          seed = seed,
          p = control$p
        )
        validation <- list(
          val_info = list(train_resamples = pls_val$train_resamples),
          results = list(train = pls_val$cv_results)
        )
        
        # Add test set predictions if Yu available
        if (!is.null(Yu)) {
          yuhat <- predict(pls_val, newdata = Xu, ncomp = fit_method$ncomp)
          
          validation$val_info$test_predictions <- yuhat
          validation$results$test <- get_validation_metrics(yuhat, Yu)
          
          tmpx <- rbind(Xr[ith_selected, , drop = FALSE], Xu)
          tmpy <- c(Yr[ith_selected, j], Yu[, j])
          complete_idx <- complete.cases(tmpy)
          tmpx <- tmpx[complete_idx, , drop = FALSE]
          tmpy <- tmpy[complete_idx]
        } else {
          tmpx <- Xr[ith_selected, , drop = FALSE]
          tmpy <- Yr[ith_selected, j]
        }
        
        # Fit final PLS model for this generation
        intermediate_pls <- assign_pls_names(opls_get_all(
          X = tmpx,
          Y = as.matrix(tmpy),
          ncomp = fit_method$ncomp,
          scale = fit_method$scale,
          maxiter = 1L,
          tol = 1e-6,
          algorithm = if (fit_method$fit_method == "mpls") "mpls" else "pls"
        ), x_names = colnames(Xr))
        
        pls_models[[j]] <- intermediate_pls
        validation_results[[j]] <- validation
      }
      
      # Name models by response variable if available
      if (!is.null(colnames(Yu))) {
        names(pls_models) <- names(validation_results) <- colnames(Yu)
      }
      models_per_generation[[i]] <- list(
        subset_size = length(ith_selected),
        pls_models = pls_models,
        validation = validation_results
      )
    }
  }

  # Progress message for final model fitting
  if (verbose) {
    cat("\nFitting final model on", length(k_idx), "selected genes...\n")
    if (!is.null(Yu) && sum(complete.cases(Yu)) > 0L) {
      cat("and", sum(complete.cases(Yu)), "target genes with available response values...\n")
    }
  }
  
  # Fit final models on selected samples
  pls_models <- list()
  validation_results <- list()
  
  for (j in seq_len(ncol(Yr))) {
    # Fit and validate PLS model for response j
    pls_val <- get_plsr(
      X = Xr,
      Y = Yr[, j, drop = FALSE],
      indices = k_idx,
      fit_method = fit_method,
      number = control$number,
      group = group[k_idx],
      retrieve = FALSE,
      tune = TRUE,
      seed = seed,
      p = control$p
    )
    
    validation <- list(
      val_info = list(train_resamples = pls_val$train_resamples),
      results = list(train = pls_val$cv_results)
    )
    
    # Add test set predictions if Yu available
    if (!is.null(Yu)) {
      yuhat <- predict(pls_val, newdata = Xu, ncomp = fit_method$ncomp)
      
      validation$val_info$test_predictions <- yuhat
      validation$results$test <- get_validation_metrics(yuhat, Yu)
      
      tmpx <- rbind(Xr[k_idx, , drop = FALSE], Xu)
      tmpy <- c(Yr[k_idx, j], Yu[, j])
      complete_idx <- complete.cases(tmpy)
      tmpx <- tmpx[complete_idx, , drop = FALSE]
      tmpy <- tmpy[complete_idx]
    } else {
      tmpx <- Xr[k_idx, , drop = FALSE]
      tmpy <- Yr[k_idx, j]
    }

    # Fit final PLS model
    finalpls <- assign_pls_names(opls_get_all(
      X = tmpx,
      Y = as.matrix(tmpy),
      ncomp = fit_method$ncomp,
      scale = fit_method$scale,
      maxiter = 1L,
      tol = 1e-6,
      algorithm = if (fit_method$fit_method == "mpls") "mpls" else "pls"
    ), x_names = colnames(Xr))
    
    pls_models[[j]] <- finalpls
    validation_results[[j]] <- validation
  }
  
  
  
  
  
  
  # Name models by response variable if available
  if (!is.null(colnames(Yu))) {
    names(pls_models) <- colnames(Yu)
    names(validation_results) <- colnames(Yu)
  }
  
  # Compile iteration weakness statistics
  iter_weakness <- list(
    maximum = data.frame(
      iteration = seq_len(nrow(max_weakness_score)),
      max_weakness_score
    ),
    cutoff = data.frame(
      iteration = seq_len(nrow(max_weakness_score)),
      cuts
    )
  )
  
  names(selected) <- paste0("iteration_", seq_along(selected))
  
  # Assemble results
  resultsList <- list(
    x_local = Xr[k_idx, , drop = FALSE],
    y_local = Yr[k_idx, , drop = FALSE],
    indices = k_idx,
    complete_iter = complete_iterations,
    iter_weakness = iter_weakness,
    samples = selected,
    n_removed = data.frame(
      iteration = seq_along(s_to_drop),
      removed = c(s_to_drop[1L], diff(s_to_drop)),
      cumulative = s_to_drop
    ),
    control = control,
    fit_method = fit_method,
    validation_results = validation_results,
    final_models = pls_models,
    intermediate_models = models_per_generation,
    seed = seed
  )
  
  if (!intermediate_models) {
    resultsList <- resultsList[!names(resultsList) %in% "intermediate_models"]
  }
  
  attr(resultsList, "call") <- call_f
  class(resultsList) <- c("gesearch", "list")
  
  resultsList
}


#' @aliases gesearch
#' @importFrom stats quantile complete.cases diffinv na.pass model.extract 
#' model.frame model.matrix
#' @export
gesearch.formula <- function(
    formula,
    train,
    test,
    k,
    b,
    target_size,
    fit_method,
    ...,
    na_action = na.pass
) {
  if (!inherits(formula, "formula")) {
    stop("'formula' must be a formula object", call. = FALSE)
  }
  
  call_f <- match.call()
  call_env <- parent.frame() 
  
  if (missing(fit_method)) {
    stop("'fit_method' is missing", call. = FALSE)
  }
  
  mf <- match.call(expand.dots = FALSE)
  
  if (!"na_action" %in% names(mf)) {
    mf[["na_action"]] <- formals(gesearch.formula)[["na_action"]]
  }
  
  # Get the model frame
  mr <- match(x = c("formula", "train", "na_action"), table = names(mf))
  mu <- match(x = c("formula", "test"), table = names(mf))
  
  mfr <- mf[c(1L, mr)]
  mfu <- mf[c(1L, mu)]
  
  names(mfr)[names(mfr) == "na_action"] <- "na.action"
  names(mfr)[names(mfr) == "train"] <- "data"
  names(mfu)[names(mfu) == "test"] <- "data"
  
  yname <- all.vars(formula, functions = FALSE, max.names = 1L)
  
  mfr[[1L]] <- mfu[[1L]] <- as.name("model.frame")
  
  input_list <- list(...)
  
  # Handle missing response in test data
  if (!yname %in% colnames(eval(mfu$data, envir = call_env))) {
    if ("optimization" %in% names(input_list)) {
      if (input_list$optimization == "response") {
        stop("'optimization = \"response\"' requires response values in test", 
             call. = FALSE)
      }
    }
    
    test <- local({
      tmp <- make.unique(c(names(test), ".tmp_y"))[length(names(test)) + 1L]
      test[[tmp]] <- NA
      names(test)[names(test) == tmp] <- yname
      test
    })
    warning(yname, " not found in test; assigned NA.", call. = FALSE)
  }
  
  mfu <- model.frame(mfu, data = test, na.action = NULL)
  mfr <- eval(mfr, call_env)
  
  trms <- attr(mfr, "terms")
  formulaclasses <- list(attr(trms, "dataClasses"))
  
  attr(trms, "intercept") <- 0L
  xr <- model.matrix(trms, model.frame(mfr, drop.unused.levels = TRUE))
  yr <- model.extract(mfr, "response")
  
  xu <- model.matrix(trms, mfu)
  yu <- model.extract(mfu, "response")
  
  # Validate response for response optimization
  if ("optimization" %in% names(input_list)) {
    if ("response" %in% input_list$optimization) {
      if (sum(is.na(yu)) == length(yu)) {
        stop("'optimization = \"response\"' requires response values in test", 
             call. = FALSE)
      }
    }
  }
  
  rsl <- gesearch(
    Xr = xr,
    Yr = yr,
    Xu = xu,
    Yu = yu,
    k = k,
    b = b,
    target_size = target_size,
    fit_method = fit_method,
    ...
  )
  rsl$formula <- formula
  rsl$dataclasses <- formulaclasses
  
  attr(rsl, "call") <- call_f
  class(rsl) <- c("gesearch", "gesearch.formula", "list")
  
  rsl
}

## ## THIS IS TO ALLOW A LIST OF FORMULAS TO BE USED FOR MULTI-RESPONSE MODELING
## @aliases gesearch
## @importFrom stats quantile complete.cases diffinv na.pass model.extract
## model.frame model.matrix
## @export
# gesearch.list <- function(
#     formula,
#     train,
#     test,
#     k,
#     b,
#     target_size,
#     fit_method,
#     ...,
#     na_action = na.pass
# ) {
#   # Validate formula list
#   is_formula <- vapply(formula, inherits, logical(1L), what = "formula")
#   if (any(!is_formula)) {
#     stop("All elements in 'formula' must be formula objects", call. = FALSE)
#   }
#   
#   right_side <- vapply(formula, function(x) labels(terms(x)), character(1L))
#   left_side <- vapply(formula, function(x) all.vars(update(x, . ~ 1)), character(1L))
#   
#   mresponse <- paste(c("matrix_response", left_side), collapse = "_")
#   mfml <- as.formula(paste(mresponse, "~", unique(right_side)))
#   
#   if (length(unique(right_side)) > 1L) {
#     stop("Right-hand side variables must be identical across all formulas", 
#          call. = FALSE)
#   }
#   
#   call_f <- match.call()
#   
#   
#   definition <- sys.function(sys.parent())
#   mf <- match.call(expand.dots = FALSE)
#   formals <- formals(definition)
#   
#   mf[["formula"]] <- mfml
#   
#   if (!"na_action" %in% names(mf)) {
#     mf[["na_action"]] <- formals[["na_action"]]
#     match.call(definition, mf, TRUE)
#   }
#   
#   # Get the model frame
#   mr <- match(x = c("formula", "train", "na_action"), table = names(mf))
#   mu <- match(x = c("formula", "test"), table = names(mf))
#   
#   mfr <- mf[c(1L, mr)]
#   mfu <- mf[c(1L, mu)]
#   
#   names(mfr)[names(mfr) == "na_action"] <- "na.action"
#   names(mfr)[names(mfr) == "train"] <- "data"
#   names(mfu)[names(mfu) == "test"] <- "data"
#   
#   yname <- left_side
#   
#   mfr[[1L]] <- mfu[[1L]] <- as.name("model.frame")
#   
#   input_list <- list(...)
#   
#   # Handle missing response in test data
#   if (!any(yname %in% colnames(eval(mfu$data)))) {
#     if ("optimization" %in% names(input_list)) {
#       if (input_list$optimization == "response") {
#         stop("'optimization = \"response\"' requires response values in test", 
#              call. = FALSE)
#       }
#     }
#     
#     missing_vars <- setdiff(yname, names(test))
#     if (length(missing_vars) > 0L) {
#       test[missing_vars] <- NA
#       warning(
#         paste(missing_vars, collapse = ", "), " not found in test; assigned NA.",
#         call. = FALSE
#       )
#     }
#   }
#   
#   train[[mresponse]] <- as.matrix(train[, left_side])
#   test[[mresponse]] <- as.matrix(test[, left_side])
#   
#   mfu <- model.frame(mfu, data = test, na.action = NULL)
#   mfr <- model.frame(mfr, data = train, na.action = NULL)
#   
#   trms <- attr(mfr, "terms")
#   
#   # Build formula classes
#   formulaclasses <- vector("list", length(yname))
#   wformulaclasses <- attr(trms, "dataClasses")
#   wformulaclasses[[1L]] <- "nmatrix.1"
#   for (i in seq_along(yname)) {
#     formulaclasses[[i]] <- wformulaclasses
#     names(formulaclasses[[i]])[1L] <- yname[i]
#   }
#   
#   attr(trms, "intercept") <- 0L
#   xr <- model.matrix(trms, model.frame(mfr, drop.unused.levels = TRUE))
#   yr <- model.extract(mfr, "response")
#   
#   xu <- model.matrix(trms, mfu)
#   yu <- model.extract(mfu, "response")
#   
#   if (!inherits(yu, "matrix")) {
#     yu <- as.matrix(yu)
#   }
#   
#   # Validate response for response optimization
#   for (i in seq_len(ncol(yu))) {
#     if ("optimization" %in% names(input_list) && 
#         sum(is.na(yu[, i])) == length(yu[, i])) {
#       stop(
#         "'optimization = \"response\"' requires response values in test. ",
#         "Check: ", colnames(yu)[i],
#         call. = FALSE
#       )
#     }
#   }
#   
#   rsl <- gesearch(
#     Xr = xr,
#     Yr = yr,
#     Xu = xu,
#     Yu = yu,
#     k = k,
#     b = b,
#     target_size = target_size,
#     fit_method = fit_method,
#     ...
#   )
#   rsl$formula <- formula
#   rsl$dataclasses <- formulaclasses
#   
#   attr(rsl, "call") <- call_f
#   class(rsl) <- c("gesearch", "gesearch.list", "list")
#   
#   rsl
# }


#' @aliases gesearch
#' @importFrom stats .MFclass terms delete.response model.frame model.matrix
#' @export
predict.gesearch <- function(
    object, 
    newdata, 
    type = "response", 
    what = c("final", "all_generations"), 
    ...
) {
  what <- match.arg(what)
  # Validate type
  allowed_types <- c("response")
  if (!is.character(type) || length(type) != 1L || is.na(type) || 
      !(type %in% allowed_types)) {
    stop(
      "'type' must be one of: ", 
      paste(shQuote(allowed_types), collapse = ", "), 
      call. = FALSE
    )
  }
  
  if (missing(newdata)) {
    stop("'newdata' is required", call. = FALSE)
  }
  
  # Handle formula-fitted models
  if (!is.null(object$formula)) {
    dcls <- object$dataclasses[[1L]][-1L]
    
    if (!inherits(newdata, "matrix") && !inherits(newdata, "data.frame")) {
      stop("'newdata' must be a data.frame or matrix for formula-fitted models", 
           call. = FALSE)
    }
    
    if (inherits(newdata, "data.frame")) {
      missing_vars <- setdiff(names(dcls), names(newdata))
      if (length(missing_vars) > 0L) {
        stop("Missing predictor variables: ", paste(missing_vars, collapse = ", "), 
             call. = FALSE)
      }
    }
    
    # Handle non-numeric predictors
    if (any(dcls != "numeric")) {
      if (inherits(newdata, "matrix") && length(dcls) == 1L) {
        if (.MFclass(newdata) == dcls) {
          pnames <- gsub(
            names(dcls), "",
            colnames(object$final_models[[1L]]$coefficients)
          )
          if (all(pnames %in% colnames(newdata))) {
            newdata_temp <- newdata
            newdata <- data.frame(rep(NA, nrow(newdata)))
            colnames(newdata) <- names(object$dataclasses[[1L]])[1L]
            newdata[[names(dcls)]] <- newdata_temp[, pnames]
          } else {
            stop("Missing predictor variables", call. = FALSE)
          }
        }
      }
    }
    
    oterms <- terms(object$formula)
    oterms <- delete.response(oterms)
    attr(oterms, "intercept") <- 0L
    mf <- model.frame(oterms, newdata)
    newdata <- model.matrix(oterms, model.frame(mf, drop.unused.levels = TRUE))
  }
  
  if (!inherits(newdata, "matrix")) {
    stop("'newdata' must be a matrix", call. = FALSE)
  }
  
  # Generate predictions from final models
  if (type == "response") {
    preds <- vector("list", length(object$final_models))
    for (i in seq_along(object$final_models)) {
      preds[[i]] <- predict_opls(
        bo = object$final_models[[i]]$bo,
        b = t(object$final_models[[i]]$coefficients),
        ncomp = object$final_models[[i]]$ncomp,
        newdata = newdata,
        scale = object$fit_method$scale,
        Xscale = object$final_models[[i]]$transf$Xscale
      )
      rownames(preds[[i]]) <- rownames(newdata)
      colnames(preds[[i]]) <- paste0("ncomp_", seq_len(ncol(preds[[i]])))
    }
    names(preds) <- names(object$final_models)
  }
  
  # Include intermediate generation predictions if requested
  if (type == "response" && what == "all_generations") {
    if (is.null(object$intermediate_models)) {
      warning(
        "Model object does not contain intermediate models; ",
        "predictions generated only from final model.",
        call. = FALSE
      )
    } else {
      n_intermediate <- length(object$intermediate_models)
      intermediate_preds <- vector("list", n_intermediate + 1L)
      
      for (j in seq_len(n_intermediate)) {
        n_models <- length(object$intermediate_models[[j]]$pls_models)
        intermediate_preds[[j]] <- vector("list", n_models)
        
        for (i in seq_len(n_models)) {
          intermediate_preds[[j]][[i]] <- predict_opls(
            bo = object$intermediate_models[[j]]$pls_models[[i]]$bo,
            b = t(object$intermediate_models[[j]]$pls_models[[i]]$coefficients),
            ncomp = object$intermediate_models[[j]]$pls_models[[i]]$ncomp,
            newdata = newdata,
            scale = object$fit_method$scale,
            Xscale = object$intermediate_models[[j]]$pls_models[[i]]$transf$Xscale
          )
          rownames(intermediate_preds[[j]][[i]]) <- rownames(newdata)
          colnames(intermediate_preds[[j]][[i]]) <- paste0(
            "ncomp_", seq_len(ncol(intermediate_preds[[j]][[i]]))
          )
        }
        names(intermediate_preds[[j]]) <- names(
          object$intermediate_models[[j]]$pls_models
        )
      }
      
      # Add final model predictions as last generation
      intermediate_preds[[n_intermediate + 1L]] <- preds
      names(intermediate_preds) <- paste0("generation_", seq_along(intermediate_preds))
      preds <- intermediate_preds
    }
  }
  
  preds
}


#' @aliases gesearch
#' @export
plot.gesearch <- function(x, which = c("weakness", "removed"), ...) {
  which <- match.arg(which)
  
  # Save and restore graphical parameters on exit
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar), add = TRUE)
  
  if (which == "weakness") {
    nbox <- ceiling(sqrt(ncol(x$iter_weakness$maximum) - 1L))
    par(mfrow = c(nbox, nbox), mar = c(4, 5.5, 2.5, 2))
  }
  
  # Process plot arguments
  plot_dots <- list(...)
  
  if (!"col" %in% names(plot_dots)) {
    plot_dots$col <- "dodgerblue"
  }
  
  if (!"type" %in% names(plot_dots)) {
    plot_dots$type <- "b"
  }
  
  # Remove axis labels from dots (set explicitly below)
  plot_dots$xlab <- NULL
  plot_dots$ylab <- NULL
  
  # Extract or set main title
  if ("main" %in% names(plot_dots)) {
    main <- plot_dots$main
    plot_dots$main <- NULL
  } else {
    main <- "gesearch results"
  }
  
  grid_col <- rgb(0.3, 0.3, 0.3, 0.3)
  
  # Plot weakness scores
  if (which == "weakness") {
    pnames <- colnames(x$iter_weakness$maximum)[-1L]
    for (i in seq_along(pnames)) {
      ith_ylab <- gsub("_", " ", pnames[i])
      ith_ylab <- gsub("rmse ", "", ith_ylab)
      ith_ylab <- paste0("Max. ", ith_ylab)
      do.call(
        plot,
        c(
          list(
            x = x$iter_weakness$maximum[[1L]],
            y = x$iter_weakness$maximum[[1L + i]],
            xlab = "Generation",
            ylab = ith_ylab
          ),
          plot_dots
        )
      )
      grid(col = grid_col, lty = 1L)
    }
  }
  
  # Plot cumulative removals
  if (which == "removed") {
    do.call(
      plot,
      c(
        list(
          x = x$n_removed$iteration,
          y = x$n_removed$cumulative,
          xlab = "Generation",
          ylab = "Samples removed (cumulative)"
        ),
        plot_dots
      )
    )
    grid(col = grid_col, lty = 1L)
  }
  
  mtext(main, outer = TRUE, cex = 2, line = -2)
  invisible(NULL)
}
