#' @title Orthogonal projections using PCA and PLS
#'
#' @description
#' \loadmathjax
#' Performs orthogonal projections of high-dimensional data matrices using
#' principal component analysis (PCA) or partial least squares (PLS).
#'
#' @param Xr A numeric matrix of reference observations (rows) and variables
#'   (columns).
#' @param Xu An optional matrix of additional observations to project.
#' @param Yr An optional response matrix. Required for PLS methods
#'   (\code{"pls"}, \code{"mpls"}, \code{"simpls"}) and when using 
#'   \code{\link{ncomp_by_opc}()}.
#' @param ncomp Component selection method. Either:
#'   \itemize{
#'     \item A positive integer (equivalent to \code{ncomp_fixed(n)})
#'     \item An \code{ncomp_selection} object: \code{\link{ncomp_by_var}()},
#'           \code{\link{ncomp_by_cumvar}()}, \code{\link{ncomp_by_opc}()}, or
#'           \code{\link{ncomp_fixed}()}
#'   }
#'   Default is \code{ncomp_by_var(0.01)}.
#' @param method A character string specifying the projection method:
#'   \itemize{
#'     \item \code{"pca"}: PCA via singular value decomposition (default)
#'     \item \code{"pca_nipals"}: PCA via NIPALS algorithm
#'     \item \code{"pls"}: PLS via NIPALS algorithm
#'     \item \code{"mpls"}: Modified PLS via NIPALS (Shenk and Westerhaus, 1991)
#'     \item \code{"simpls"}: PLS via SIMPLS algorithm (de Jong, 1993)
#'   }
#' @param center A logical indicating whether to center the data. Default is
#'   \code{TRUE}. PLS methods always center internally regardless of this
#'   setting.
#' @param scale A logical indicating whether to scale the data to unit
#'   variance. Default is \code{FALSE}.
#' @param tol Convergence tolerance for the NIPALS algorithm. Default is
#'   \code{1e-6}. Ignored when \code{method = "simpls"}.
#' @param max_iter Maximum number of iterations for NIPALS. Default is
#'   \code{1000}. Ignored when \code{method = "simpls"}.
#' @param pc_selection `r lifecycle::badge("deprecated")` Use \code{ncomp} instead.
#' @param ... Additional arguments (currently unused).
#'
#' @details
#' ## PCA methods
#'
#' When \code{method = "pca"}, singular value decomposition factorizes the
#' data matrix \mjeqn{X}{X} as:
#'
#' \mjdeqn{X = UDV^{T}}{X = UDV^T}
#'
#' where \mjeqn{U}{U} and \mjeqn{V}{V} are orthogonal matrices (left and right
#' singular vectors), and \mjeqn{D}{D} is a diagonal matrix of singular values.
#' The score matrix is \mjeqn{UD}{UD} and the loadings are \mjeqn{V}{V}.
#'
#' When \code{method = "pca_nipals"}, the non-linear iterative partial least
#' squares (NIPALS) algorithm is used instead.
#'
#' ## PLS methods
#'
#' Three PLS variants are available:
#'
#' \itemize{
#'   \item \code{"pls"}: Standard PLS using the NIPALS algorithm with 
#'     covariance-based weights.
#'   \item \code{"mpls"}: Modified PLS using the NIPALS algorithm with 
#'     correlation-based weights, giving equal influence to all predictors 
#'     regardless of variance (Shenk and Westerhaus, 1991).
#'   \item \code{"simpls"}: SIMPLS algorithm (de Jong, 1993), which deflates 
#'     the cross-product matrix rather than X itself. Computationally faster 
#'     than NIPALS, especially for wide matrices.
#' }
#'
#' ## Component selection
#'
#' When \code{\link{ncomp_by_opc}()} is used, component selection minimizes
#' RMSD (for continuous \code{Yr}) or maximizes kappa (for categorical
#' \code{Yr}) between observations and their nearest neighbors. See
#' \code{\link{diss_evaluate}}.
#'
#' @return
#' An object of class \code{"ortho_projection"} containing:
#' \itemize{
#'   \item \code{scores}: Matrix of projected scores for \code{Xr} (and \code{Xu}).
#'   \item \code{X_loadings}: Matrix of X loadings.
#'   \item \code{Y_loadings}: Matrix of Y loadings (PLS only).
#'   \item \code{weights}: Matrix of PLS weights (PLS only).
#'   \item \code{projection_mat}: Projection matrix for new data (PLS only).
#'   \item \code{variance}: List with original and explained variance.
#'   \item \code{scores_sd}: Standard deviation of scores.
#'   \item \code{ncomp}: Number of components retained.
#'   \item \code{center}: Centering vector used.
#'   \item \code{scale}: Scaling vector used.
#'   \item \code{method}: Projection method used.
#'   \item \code{ncomp_method}: The value passed to the `ncomp` argument.
#'   \item \code{opc_evaluation}: opc optimization results (if applicable).
#' }
#'
#' @author
#' \href{https://orcid.org/0000-0002-5369-5120}{Leonardo Ramirez-Lopez}
#'
#' @references
#' de Jong, S. 1993. SIMPLS: An alternative approach to partial least squares 
#' regression. Chemometrics and Intelligent Laboratory Systems 18:251-263.
#'
#' Martens, H. 1991. Multivariate calibration. John Wiley & Sons.
#'
#' Ramirez-Lopez, L., Behrens, T., Schmidt, K., Stevens, A., Dematte, J.A.M.,
#' Scholten, T. 2013a. The spectrum-based learner: A new local approach for
#' modeling soil vis-NIR spectra of complex data sets. Geoderma 195-196:268-279.
#'
#' Ramirez-Lopez, L., Behrens, T., Schmidt, K., Viscarra Rossel, R., Dematte,
#' J.A.M., Scholten, T. 2013b. Distance and similarity-search metrics for use
#' with soil vis-NIR spectra. Geoderma 199:43-53.
#'
#' Shenk, J.S., Westerhaus, M.O. 1991. Populations structuring of near
#' infrared spectra and modified partial least squares regression. Crop
#' Science 31:1548-1555.
#'
#' @seealso
#' \code{\link{ncomp_by_var}}, \code{\link{ncomp_by_opc}},
#' \code{\link{diss_evaluate}}, \code{\link{mbl}}
#'
#' @examples
#' \donttest{
#' library(prospectr)
#' data(NIRsoil)
#'
#' # Preprocess
#' sg_det <- savitzkyGolay(
#'   detrend(NIRsoil$spc, wav = as.numeric(colnames(NIRsoil$spc))),
#'   m = 1, p = 1, w = 7
#' )
#'
#' # Split data
#' train_x <- sg_det[NIRsoil$train == 1 & !is.na(NIRsoil$CEC), ]
#' train_y <- NIRsoil$CEC[NIRsoil$train == 1 & !is.na(NIRsoil$CEC)]
#' test_x <- sg_det[NIRsoil$train == 0 & !is.na(NIRsoil$CEC), ]
#'
#' # PCA with fixed components
#' proj <- ortho_projection(train_x, ncomp = 5)
#'
#' # PCA with variance-based selection
#' proj <- ortho_projection(train_x, ncomp = ncomp_by_var(0.01))
#'
#' # PCA with OPC optimization
#' proj <- ortho_projection(train_x, Xu = test_x, Yr = train_y,
#'                          ncomp = ncomp_by_opc(40))
#'
#' # PLS projection (NIPALS)
#' proj <- ortho_projection(train_x, Xu = test_x, Yr = train_y,
#'                          method = "pls", ncomp = ncomp_by_opc(40))
#'
#' # Modified PLS
#' proj <- ortho_projection(train_x, Yr = train_y,
#'                          method = "mpls", ncomp = 10)
#'
#' # SIMPLS (faster for wide matrices)
#' proj <- ortho_projection(train_x, Yr = train_y,
#'                          method = "simpls", ncomp = 10)
#' }
#'
#' @rdname ortho_projection
#' @export
ortho_projection <- function(
    Xr, Xu = NULL, Yr = NULL,
    ncomp = ncomp_by_var(0.01),
    method = c("pca", "pca_nipals", "pls", "mpls", "simpls"),
    center = TRUE,
    scale = FALSE,
    tol = 1e-6,
    max_iter = 1000L,
    pc_selection = deprecated(),
    ...
) {
  
  # ---------------------------------------------------------------------------
  # Handle deprecated pc_selection
  # ---------------------------------------------------------------------------
  if (is_present(pc_selection)) {
    deprecate_warn(
      when = "3.0.0",
      what = "ortho_projection(pc_selection)",
      with = "ortho_projection(ncomp)",
      details = "Use ncomp_fixed(), ncomp_by_var(), ncomp_by_cumvar(), or ncomp_by_opc()."
    )
    ncomp <- .convert_pc_selection_to_ncomp(pc_selection)
  }
  
  # ---------------------------------------------------------------------------
  # Handle legacy method names
  # ---------------------------------------------------------------------------
  if (identical(method, "pca.nipals")) {
    
    deprecate_warn(
      when = "3.0.0",
      what = I('ortho_projection(method = "pca.nipals")'),
      with = I('ortho_projection(method = "pca_nipals")')
    )
    method <- "pca_nipals"
  }
  
  method <- match.arg(method)
  is_pca <- method %in% c("pca", "pca_nipals")
  is_pls <- method %in% c("pls", "mpls", "simpls")
  
  # ---------------------------------------------------------------------------
  # Coerce ncomp
  # ---------------------------------------------------------------------------
  ncomp <- .coerce_ncomp(ncomp)
  
  # ---------------------------------------------------------------------------
  # Validate inputs
  # ---------------------------------------------------------------------------
  if (!is.logical(center) || length(center) != 1L || is.na(center)) {
    stop("'center' must be TRUE or FALSE.")
  }
  if (!is.logical(scale) || length(scale) != 1L || is.na(scale)) {
    stop("'scale' must be TRUE or FALSE.")
  }
  
  if (is_pls && is.null(Yr)) {
    stop("'Yr' is required for PLS methods ('pls', 'mpls').")
  }
  
  if (inherits(ncomp, "ncomp_by_opc") && is.null(Yr)) {
    stop("'Yr' is required when using ncomp_by_opc().")
  }
  
  # ---------------------------------------------------------------------------
  # Dispatch to implementation
  # ---------------------------------------------------------------------------
  if (is_pls) {
    proj <- .pls_projection_impl(
      Xr = Xr, Xu = Xu, Yr = Yr,
      ncomp = ncomp,
      method = method,
      scale = scale,
      tol = tol,
      max_iter = max_iter
    )
  } else {
    proj <- .pca_projection_impl(
      Xr = Xr, Xu = Xu, Yr = Yr,
      ncomp = ncomp,
      method = method,
      center = center,
      scale = scale,
      tol = tol,
      max_iter = max_iter
    )
  }
  
  proj$method <- method
  proj$ncomp_method <- ncomp
  class(proj) <- c("ortho_projection", "list")
  proj
}


# =============================================================================
# Deprecated wrappers
# =============================================================================

#' @keywords internal
#' @noRd
pc_projection <- function(
    Xr, Xu = NULL, Yr = NULL,
    pc_selection = list(method = "var", value = 0.01),
    center = TRUE, scale = FALSE,
    method = "pca",
    tol = 1e-6, max_iter = 1000,
    ...
) {
  deprecate_warn(
    when = "3.0.0",
    what = "pc_projection()",
    with = "ortho_projection()"
  )
  
  ncomp <- .convert_pc_selection_to_ncomp(pc_selection)
  new_method <- if (method == "pca.nipals") "pca_nipals" else "pca"
  
  ortho_projection(
    Xr = Xr, Xu = Xu, Yr = Yr,
    ncomp = ncomp,
    method = new_method,
    center = center,
    scale = scale,
    tol = tol,
    max_iter = max_iter
  )
}


#' @keywords internal
#' @noRd
pls_projection <- function(
    Xr, Xu = NULL, Yr,
    pc_selection = list(method = "opc", value = min(dim(Xr), 40)),
    scale = FALSE,
    method = "pls",
    tol = 1e-6, max_iter = 1000,
    ...
) {
  deprecate_warn(
    when = "3.0.0",
    what = "pls_projection()",
    with = "ortho_projection()"
  )
  
  ncomp <- .convert_pc_selection_to_ncomp(pc_selection)
  
  ortho_projection(
    Xr = Xr, Xu = Xu, Yr = Yr,
    ncomp = ncomp,
    method = method,
    scale = scale,
    tol = tol,
    max_iter = max_iter
  )
}


# =============================================================================
# predict method
# =============================================================================

#' @rdname ortho_projection
#' @param object Object of class \code{"ortho_projection"}.
#' @param newdata Matrix of new observations to project.
#' @export
predict.ortho_projection <- function(object, newdata, ...) {
  if (missing(newdata)) {
    return(object$scores)
  }
  
  nms_orig <- colnames(object$X_loadings)
  nms_new <- colnames(newdata)
  
  if (any(!nms_orig %in% nms_new)) {
    stop("Variables missing in 'newdata' that are required for projection.")
  }
  
  is_pca <- grepl("pca", object$method)
  
  if (is_pca) {
    newdata <- sweep(newdata, 2L, object$center, FUN = "-")
    newdata <- sweep(newdata, 2L, object$scale, FUN = "/")
    scores <- newdata %*% t(object$X_loadings)
    colnames(scores) <- rownames(object$X_loadings)
    rownames(scores) <- rownames(newdata)
    return(scores)
  } else {
    scores <- project_opls(
      projection_mat = object$projection_mat,
      ncomp = ncol(object$projection_mat),
      newdata = newdata,
      scale = TRUE,
      Xcenter = object$center,
      Xscale = object$scale
    )
    colnames(scores) <- paste0("pls_", seq_len(ncol(scores)))
    rownames(scores) <- rownames(newdata)
    return(scores)
  }
}


# =============================================================================
# Helper: Convert legacy pc_selection to ncomp_selection
# =============================================================================

.convert_pc_selection_to_ncomp <- function(pc_selection) {
  # Handle character shorthand
  if (is.character(pc_selection) && length(pc_selection) == 1L) {
    pc_selection <- switch(pc_selection,
                           opc = list(method = "opc", value = 40L),
                           cumvar = list(method = "cumvar", value = 0.99),
                           var = list(method = "var", value = 0.01),
                           manual = list(method = "manual", value = 40L),
                           stop("Unknown pc_selection method: ", pc_selection)
    )
  }
  
  sel_method <- if (!is.null(pc_selection$method)) {
    pc_selection$method
  } else {
    pc_selection[[1]]
  }
  
  value <- if (!is.null(pc_selection$value)) {
    pc_selection$value
  } else {
    pc_selection[[2]]
  }
  
  switch(sel_method,
         manual = ncomp_fixed(as.integer(value)),
         var = ncomp_by_var(min_var = value, max_ncomp = 40L),
         cumvar = ncomp_by_cumvar(min_cumvar = value, max_ncomp = 40L),
         opc = ncomp_by_opc(max_ncomp = as.integer(value)),
         stop("Unknown pc_selection method: ", sel_method)
  )
}




# =============================================================================
# Helper: Get max_ncomp from ncomp_selection object
# =============================================================================

.get_max_ncomp <- function(ncomp) {
  if (inherits(ncomp, "ncomp_fixed")) {
    return(ncomp$ncomp)
  }
  ncomp$max_ncomp
}

# =============================================================================
# Helper: Bridge ncomp_selection to legacy pc_selection list
# =============================================================================

.ncomp_to_pc_selection <- function(ncomp_obj) {
  stopifnot(inherits(ncomp_obj, "ncomp_selection"))
  
  switch(
    class(ncomp_obj)[[1]],
    ncomp_fixed = list(method = "manual", value = ncomp_obj$ncomp),
    ncomp_by_var = list(method = "var", value = ncomp_obj$min_var),
    ncomp_by_cumvar = list(method = "cumvar", value = ncomp_obj$min_cumvar),
    ncomp_by_opc = list(method = "opc", value = ncomp_obj$max_ncomp),
    stop("Unknown ncomp_selection class: ", class(ncomp_obj)[[1]])
  )
}



# =============================================================================
# Internal PCA implementation
# =============================================================================

.pca_projection_impl <- function(
    Xr, Xu, Yr, ncomp, method, center, scale, tol, max_iter
) {
  # Bridge to legacy format
  pc_selection <- .ncomp_to_pc_selection(ncomp)
  pc_selection_method <- pc_selection$method
  max_comp <- .get_max_ncomp(ncomp)
  
  # Validation
  if (!is.null(Yr)) {
    Yr <- as.matrix(Yr)
    if (nrow(Yr) != nrow(Xr)) {
      stop("Number of rows in 'Yr' must match 'Xr'.")
    }
  }
  
  ny <- if (!is.null(Yr)) ncol(Yr) else 0L
  
  if (!is.null(Xu) && ncol(Xr) != ncol(Xu)) {
    stop("Number of columns in 'Xr' and 'Xu' must match.")
  }
  
  effective_rows_xr <- nrow(Xr)
  effective_rows_xu <- if (is.null(Xu)) 0L else nrow(Xu)
  n_cols <- ncol(Xr)
  
  # Combine Xr and Xu
  X <- rbind(Xr, Xu)
  
  # Cap max_comp at data dimensions
  max_comp <- min(max_comp, nrow(X), ncol(X))
  
  # Center
  if (center) {
    mean_vector <- colMeans(X)
    X0 <- sweep(X, 2L, mean_vector, FUN = "-")
  } else {
    mean_vector <- rep(0, n_cols)
    X0 <- X
  }
  
  # Scale
  if (scale) {
    sd_vector <- get_column_sds(X0)
    X0 <- sweep(X0, 2L, sd_vector, FUN = "/")
  } else {
    sd_vector <- rep(1, n_cols)
  }
  
  # Compute PCA
  if (method == "pca") {
    sv <- svd(X0, nu = max_comp, nv = max_comp)
    sv$d <- sv$d[1:max_comp]
    
    if (length(sv$d) == 1L) {
      pc_scores <- sv$u %*% as.matrix(sv$d)
    } else {
      pc_scores <- sv$u %*% diag(sv$d)
    }
    pc_loadings <- t(sv$v)
    
    xvariance <- overall_var(X0)
    explained_v <- (sv$d)^2 / xvariance
    cumulative_v <- cumsum(explained_v)
    variance <- rbind(
      var = (sv$d)^2,
      explained_var = explained_v,
      cumulative_explained_var = cumulative_v
    )
  } else {
    # pca_nipals
    nipals_result <- pca_nipals(
      X = X0,
      ncomp = max_comp,
      center = FALSE,
      scale = FALSE,
      maxiter = max_iter,
      tol = tol,
      pcSelmethod = pc_selection_method,
      pcSelvalue = pc_selection$value
    )
    
    pc_scores <- nipals_result$pc_scores
    pc_loadings <- nipals_result$pc_loadings
    xvariance <- nipals_result$original_x_variance
    variance <- rbind(
      var = nipals_result$pc_variance[1, ],
      explained_var = nipals_result$pc_variance[2, ],
      cumulative_explained_var = nipals_result$pc_variance[3, ]
    )
  }
  
  # Assign names
  colnames(pc_scores) <- paste0("pc_", seq_len(ncol(pc_scores)))
  rownames(pc_scores) <- c(
    paste0("Xr_", seq_len(effective_rows_xr)),
    if (effective_rows_xu > 0) paste0("Xu_", seq_len(effective_rows_xu))
  )
  colnames(pc_loadings) <- colnames(X0)
  rownames(pc_loadings) <- paste0("pc_", seq_len(nrow(pc_loadings)))
  colnames(variance) <- paste0("pc_", seq_len(ncol(variance)))
  
  # Select components
  opc_evaluation <- NULL
  
  if (pc_selection_method == "opc") {
    if (is.null(Yr) || !is.matrix(Yr)) {
      stop("'Yr' must be a matrix when using ncomp_by_opc().")
    }
    if (nrow(Yr) != effective_rows_xr) {
      stop("Number of rows in 'Yr' must match 'Xr'.")
    }
    if (!is.null(colnames(Yr)) && any(duplicated(colnames(Yr)))) {
      stop("Column names in 'Yr' must be unique.")
    }
    if (is.null(colnames(Yr))) {
      colnames(Yr) <- paste0("Yr_", seq_len(ny))
    }
    
    results <- eval_multi_pc_diss(
      pc_scores[seq_len(effective_rows_xr), 1:max_comp, drop = FALSE],
      side_info = Yr,
      method = "pc",
      check_dims = FALSE
    )
    selected_pcs <- results$best_pc
    opc_evaluation <- results$result
    
  } else if (pc_selection_method == "cumvar") {
    selected_pcs <- sum(variance["cumulative_explained_var", ] < pc_selection$value) + 1L
    selected_pcs <- min(selected_pcs, ncol(variance))
    
  } else if (pc_selection_method == "var") {
    selected_pcs <- sum(variance["explained_var", ] >= pc_selection$value)
    selected_pcs <- max(1L, selected_pcs)
    
  } else {
    # manual
    selected_pcs <- min(pc_selection$value, ncol(pc_scores))
  }
  
  # Subset to selected components
  scores_sd <- get_column_sds(pc_scores[, 1:selected_pcs, drop = FALSE])
  colnames(scores_sd) <- colnames(pc_scores)[1:selected_pcs]
  rownames(scores_sd) <- "sd"
  
  result <- list(
    scores = pc_scores[, 1:selected_pcs, drop = FALSE],
    X_loadings = pc_loadings[1:selected_pcs, , drop = FALSE],
    variance = list(
      original_x_var = xvariance,
      x_var = variance[, 1:selected_pcs, drop = FALSE]
    ),
    scores_sd = scores_sd,
    ncomp = selected_pcs,
    center = mean_vector,
    scale = sd_vector
  )
  
  if (!is.null(opc_evaluation)) {
    result$opc_evaluation <- opc_evaluation
  }
  
  result
}


# =============================================================================
# Internal PLS implementation
# =============================================================================

.pls_projection_impl <- function(
    Xr, Xu, Yr, ncomp, method, scale, tol, max_iter
) {
  # Bridge to legacy format
  pc_selection <- .ncomp_to_pc_selection(ncomp)
  pc_selection_method <- pc_selection$method
  max_comp <- .get_max_ncomp(ncomp)
  
  # Validation
  Yr <- as.matrix(Yr)
  if (!is.numeric(Yr)) {
    stop("'Yr' must be numeric.")
  }
  if (nrow(Yr) != nrow(Xr)) {
    stop("Number of rows in 'Yr' must match 'Xr'.")
  }
  
  if (!is.null(Xu) && ncol(Xr) != ncol(Xu)) {
    stop("Number of columns in 'Xr' and 'Xu' must match.")
  }
  
  ny <- ncol(Yr)
  # Handle missing values in Yr
  nas <- rowSums(is.na(Yr)) > 0
  X0 <- Xr
  Y0 <- Yr
  non_nas_yr <- seq_len(nrow(Xr))
  nas_yr <- integer(0)
  Xout <- NULL
  if (any(nas)) {
    nas_yr <- which(nas)
    non_nas_yr <- which(!nas)
    Xout <- Xr[nas_yr, , drop = FALSE]
    X0 <- Xr[non_nas_yr, , drop = FALSE]
    Y0 <- Yr[non_nas_yr, , drop = FALSE]
    
    if (pc_selection_method %in% c("opc", "manual")) {
      if (min(dim(X0)) < max_comp) {
        stop(
          "Missing values in 'Yr' reduce available observations. ",
          "The number of components exceeds available observations."
        )
      }
    }
  }
  effective_rows_xr <- nrow(X0)
  max_comp <- min(max_comp, effective_rows_xr, ncol(X0))
  
  # Set up for C++ backend
  cpp_method <- if (pc_selection_method %in% c("opc", "manual")) {
    "manual"
  } else {
    pc_selection_method
  }
  
  cpp_value <- if (pc_selection_method %in% c("opc", "manual")) {
    max_comp - 1L
  } else {
    pc_selection$value
  }
  # Run PLS
  plsp <- opls_for_projection(
    X = X0,
    Y = Y0,
    ncomp = max_comp,
    scale = scale,
    maxiter = max_iter,
    tol = tol,
    pcSelmethod = cpp_method,
    pcSelvalue = cpp_value,
    algorithm = method
  )
  
  max_comp <- plsp$ncomp
  
  # OPC selection if requested
  opc_evaluation <- NULL
  if (pc_selection_method == "opc") {
    if (is.null(colnames(Yr))) {
      colnames(Y0) <- colnames(Yr) <- paste0("Yr_", seq_len(ny))
    }
    if (any(duplicated(colnames(Yr)))) {
      stop("Column names in 'Yr' must be unique.")
    }
    
    results <- eval_multi_pc_diss(
      plsp$scores[, 1:max_comp, drop = FALSE],
      side_info = Y0,
      method = "pls",
      check_dims = FALSE
    )
    max_comp <- results$best_pc
    opc_evaluation <- results$result
  }
  
  # Extract selected components
  pls_variance <- plsp$variance$x_var[, 1:max_comp, drop = FALSE]
  weights <- plsp$weights[1:max_comp, , drop = FALSE]
  scores <- plsp$scores[, 1:max_comp, drop = FALSE]
  X_loadings <- plsp$X_loadings[1:max_comp, , drop = FALSE]
  Y_loadings <- plsp$Y_loadings[1:max_comp, , drop = FALSE]
  projection_mat <- plsp$projection_mat[, 1:max_comp, drop = FALSE]
  
  # Assign names
  comp_names <- paste0("pls_", seq_len(max_comp))
  rownames(projection_mat) <- colnames(X_loadings) <- colnames(X0)
  colnames(projection_mat) <- rownames(X_loadings) <- comp_names
  rownames(Y_loadings) <- comp_names
  colnames(Y_loadings) <- paste0("Y_pls_", seq_len(ny))
  rownames(weights) <- comp_names
  colnames(weights) <- colnames(X0)
  colnames(pls_variance) <- comp_names
  rownames(pls_variance) <- c("var", "explained_var_X", "cumulative_explained_var_X")
  
  yex <- plsp$variance$y_var[, 1:max_comp, drop = FALSE]
  colnames(yex) <- comp_names
  if (ny > 1) {
    rownames(yex) <- paste0("cumulative_explained_var_", colnames(Yr))
  } else {
    rownames(yex) <- "cumulative_explained_var_Yr"
  }
  
  # Handle observations with missing Yr
  if (length(nas_yr) > 0) {
    scores_full <- matrix(NA_real_, nrow(Xr), ncol(scores))
    scores_full[nas_yr, ] <- project_opls(
      projection_mat = projection_mat,
      ncomp = max_comp,
      newdata = Xout,
      scale = scale,
      Xcenter = plsp$transf$Xcenter,
      Xscale = plsp$transf$Xscale
    )
    scores_full[non_nas_yr, ] <- scores
    scores <- scores_full
  }
  
  rownames(scores) <- paste0("Xr_", seq_len(nrow(scores)))
  colnames(scores) <- comp_names
  
  # Project Xu if provided
  if (!is.null(Xu)) {
    if (is.vector(Xu)) Xu <- t(Xu)
    scores_Xu <- project_opls(
      projection_mat = projection_mat,
      ncomp = max_comp,
      newdata = Xu,
      scale = scale,
      Xcenter = plsp$transf$Xcenter,
      Xscale = plsp$transf$Xscale
    )
    colnames(scores_Xu) <- comp_names
    rownames(scores_Xu) <- paste0("Xu_", seq_len(nrow(Xu)))
    scores <- rbind(scores, scores_Xu)
  }
  
  # Ensure scale is a matrix
  if (!nrow(plsp$transf$Xscale)) {
    plsp$transf$Xscale <- matrix(1, 1L, length(plsp$transf$Xcenter))
  }
  
  scores_sd <- get_column_sds(scores)
  colnames(scores_sd) <- colnames(scores)
  rownames(scores_sd) <- "sd"
  
  result <- list(
    scores = scores,
    X_loadings = X_loadings,
    Y_loadings = Y_loadings,
    weights = weights,
    projection_mat = projection_mat,
    variance = list(
      original_x_var = plsp$variance$original_x_var,
      x_var = pls_variance,
      y_var = yex
    ),
    scores_sd = scores_sd,
    ncomp = max_comp,
    center = plsp$transf$Xcenter,
    scale = plsp$transf$Xscale
  )
  
  if (!is.null(opc_evaluation)) {
    result$opc_evaluation <- opc_evaluation
  }
  
  result
}