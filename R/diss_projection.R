# =============================================================================
# Orthogonal dissimilarity — redesigned API
# =============================================================================
#
# EXPORTED (public API):
#   diss_pca()            — method constructor
#   diss_pls()            — method constructor
#
# NOT EXPORTED (internal, reserved for future API):
#   diss_local_pca()      — local projection type of diss_pca
#   diss_local_pls()      — local projection type of diss_pls
#
# INTERNAL (not exported, not for direct use):
#   .ortho_diss_compute() — shared compute function
#   .ortho_diss_local()   — local projection compute helper
#   .validate_ncomp()     — ncomp argument validator
#   .validate_shared_args() — shared argument validator
#   .validate_pre_k()     — pre_k argument validator


# =============================================================================
# diss_pca — EXPORTED
# =============================================================================

#' @title PCA dissimilarity method constructor
#'
#' @description
#' Creates a configuration object for computing dissimilarity based on
#' Mahalanobis distance in a PCA-projected space. Pass the result to
#' \code{dissimilarity()} to compute the dissimilarity matrix.
#'
#' After projecting onto the PCA space, scores are standardized by their
#' standard deviation before Euclidean distance is computed — which is
#' equivalent to Mahalanobis distance in the original space.
#'
#' To compute all pairwise distances among a combined set of observations,
#' pass \code{rbind(Xr, Xu)} as the \code{Xr} argument to
#' \code{dissimilarity()} rather than using a separate \code{Xu} argument.
#'
#' @param ncomp A positive integer. Number of components to compute. When
#'   \code{ncomp_selection} is \code{NULL} (the default), exactly
#'   \code{ncomp} components are used. When an automatic selection strategy
#'   is provided via \code{ncomp_selection}, \code{ncomp} acts as the
#'   upper bound. Default \code{20}.
#' @param ncomp_selection Either \code{NULL} (use exactly \code{ncomp}
#'   components) or an \code{nncomp_selection} object created with
#'   \code{\link{ncomp_by_var}()}, \code{\link{ncomp_by_cumvar}()}, or
#'   \code{\link{ncomp_by_opc}()} for automatic selection up to \code{ncomp}.
#'   When \code{ncomp_by_opc()} is used, \code{Yr} must be passed to
#'   \code{dissimilarity()}. Default \code{NULL}.
#' @param algorithm Character. PCA algorithm. Either \code{"svd"} (default)
#'   or \code{"nipals"}.
#' @param center Logical. Center the data before projecting? Default
#'   \code{TRUE}.
#' @param scale Logical. Scale the data before projecting? Default
#'   \code{FALSE}.
#' @param return_projection Logical. Include the \code{ortho_projection}
#'   object in the output of \code{dissimilarity()}? Default \code{FALSE}.
#' @param allow_parallel Logical. Allow parallel computation? Default
#'   \code{TRUE}.
#'
#' @return An object of class \code{c("diss_pca", "diss_method")}.
#' @seealso \code{\link{diss_pls}}, \code{\link{nncomp_selection}},
#'   \code{\link{dissimilarity}}
#' @examples
#' # Fixed number of components
#' m <- diss_pca(ncomp = 10)
#'
#' # Automatic selection by OPC, evaluating up to 30 components
#' m <- diss_pca(ncomp = 30, ncomp_selection = ncomp_by_opc())
#'
#' # Automatic selection by cumulative variance, capped at 20 components
#' m <- diss_pca(ncomp = 20, ncomp_selection = ncomp_by_cumvar(min_cumvar = 0.99))
#' @export
diss_pca <- function(
    ncomp                = 20,
    ncomp_selection  = NULL,
    algorithm            = c("svd", "nipals"),
    center               = TRUE,
    scale                = FALSE,
    return_projection    = FALSE,
    allow_parallel       = TRUE
) {
  .validate_ncomp(ncomp)
  .validate_ncomp_selection(ncomp_selection)
  algorithm <- match.arg(algorithm)
  .validate_shared_args(
    center = center, scale = scale,
    return_projection = return_projection,
    allow_parallel = allow_parallel
  )
  
  .new_diss_method(
    "pca",
    center = center,
    scale  = scale,
    extra  = list(
      ncomp = as.integer(ncomp),
      ncomp_selection = ncomp_selection,
      algorithm = algorithm,
      local = FALSE,
      pre_k = NULL,
      return_projection = return_projection,
      allow_parallel = allow_parallel
    )
  )
}


# =============================================================================
# diss_pls — EXPORTED
# =============================================================================

#' @title PLS dissimilarity method constructor
#'
#' @description
#' Creates a configuration object for computing dissimilarity based on
#' Mahalanobis distance in a PLS-projected space. Pass the result to
#' \code{dissimilarity()} along with \code{Yr} (always required for PLS).
#'
#' @param ncomp A positive integer. Number of components to compute. Acts as
#'   the exact number when \code{ncomp_selection = NULL}, or as the upper
#'   bound when an automatic selection strategy is provided. Default \code{20}.
#' @param ncomp_selection Either \code{NULL} or an \code{nncomp_selection}
#'   object. See \code{\link{nncomp_selection}} for details. Default \code{NULL}.
#' @param type Character. PLS type. Either \code{"pls"} (default) or
#'   \code{"mpls"} (modified PLS, Shenk and Westerhaus, 1991).
#' @param center Logical. Center the data? Note: PLS projections are always
#'   centered internally regardless of this setting. Default \code{TRUE}.
#' @param scale Logical. Scale the data? Default \code{FALSE}.
#' @param return_projection Logical. Include the projection object in output?
#'   Default \code{FALSE}.
#' @param allow_parallel Logical. Allow parallel computation? Default
#'   \code{TRUE}.
#'
#' @return An object of class \code{c("diss_pls", "diss_method")}.
#' @seealso \code{\link{diss_pca}}, \code{\link{nncomp_selection}},
#'   \code{\link{dissimilarity}}
#' @examples
#' # Fixed number of components
#' m <- diss_pls(ncomp = 6)
#'
#' # Automatic selection by OPC, evaluating up to 30 components
#' m <- diss_pls(ncomp = 30, ncomp_selection = ncomp_by_opc())
#'
#' # Modified PLS with variance-based selection
#' m <- diss_pls(ncomp = 20, type = "mpls",
#'               ncomp_selection = ncomp_by_var(min_var = 0.01))
#' @export
diss_pls <- function(
    ncomp                = 20,
    ncomp_selection  = NULL,
    type = c("pls", "mpls"),
    center               = TRUE,
    scale                = FALSE,
    return_projection    = FALSE,
    allow_parallel       = TRUE
) {
  .validate_ncomp(ncomp)
  .validate_ncomp_selection(ncomp_selection)
  type <- match.arg(type)
  .validate_shared_args(
    center = center, scale = scale,
    return_projection = return_projection,
    allow_parallel = allow_parallel
  )
  
  .new_diss_method(
    "pls",
    center = center,
    scale  = scale,
    extra  = list(
      ncomp               = as.integer(ncomp),
      ncomp_selection = ncomp_selection,
      type = type,
      local               = FALSE,
      pre_k               = NULL,
      return_projection   = return_projection,
      allow_parallel      = allow_parallel
    )
  )
}


# =============================================================================
# diss_local_pca — NOT EXPORTED
# Reserved for future API. Accessible via resemble:::diss_local_pca()
# =============================================================================
#
# Computes dissimilarity based on local PCA projections. Unlike diss_pca(),
# each observation's dissimilarity to its neighbors is computed in a
# locally-fitted PCA space. The resulting matrix is asymmetric and cannot
# be considered a proper distance metric.
#
# @param ncomp Positive integer. Exact number or upper bound of components.
# @param ncomp_selection NULL or an nncomp_selection object.
# @param algorithm Character. "svd" or "nipals".
# @param pre_k Positive integer >= 4. Neighborhood size for local projections.
# @param center Logical.
# @param scale Logical.
# @param return_projection Logical.
# @param allow_parallel Logical.
#
# @return An object of class c("diss_local_pca", "diss_pca", "diss_method").
#
# NOTE: @export intentionally omitted — not part of public API.
diss_local_pca <- function(
    ncomp                = 20,
    ncomp_selection  = NULL,
    algorithm            = c("svd", "nipals"),
    pre_k,
    center               = TRUE,
    scale                = FALSE,
    return_projection    = FALSE,
    allow_parallel       = TRUE
) {
  .validate_ncomp(ncomp)
  .validate_ncomp_selection(ncomp_selection)
  algorithm <- match.arg(algorithm)
  
  if (missing(pre_k)) {
    stop("'pre_k' is required for diss_local_pca().")
  }
  .validate_pre_k(pre_k)
  .validate_shared_args(
    center = center, scale = scale,
    return_projection = return_projection,
    allow_parallel = allow_parallel
  )
  
  obj <- .new_diss_method(
    "pca",
    center = center,
    scale  = scale,
    extra  = list(
      ncomp               = as.integer(ncomp),
      ncomp_selection = ncomp_selection,
      algorithm           = algorithm,
      local               = TRUE,
      pre_k               = as.integer(pre_k),
      return_projection   = return_projection,
      allow_parallel      = allow_parallel
    )
  )
  class(obj) <- c("diss_local_pca", class(obj))
  obj
}


# =============================================================================
# diss_local_pls — NOT EXPORTED
# Reserved for future API. Accessible via resemble:::diss_local_pls()
# =============================================================================
#
# Computes dissimilarity based on local PLS projections. Requires Yr to be
# passed to dissimilarity(). The resulting matrix is asymmetric.
#
# @param ncomp Positive integer. Exact number or upper bound of components.
# @param ncomp_selection NULL or an nncomp_selection object.
# @param type Character. "pls" or "mpls".
# @param pre_k Positive integer >= 4.
# @param center Logical.
# @param scale Logical.
# @param return_projection Logical.
# @param allow_parallel Logical.
#
# @return An object of class c("diss_local_pls", "diss_pls", "diss_method").
#
# NOTE: @export intentionally omitted — not part of public API.
diss_local_pls <- function(
    ncomp                = 20,
    ncomp_selection  = NULL,
    type = c("pls", "mpls"),
    pre_k,
    center               = TRUE,
    scale                = FALSE,
    return_projection    = FALSE,
    allow_parallel       = TRUE
) {
  .validate_ncomp(ncomp)
  .validate_ncomp_selection(ncomp_selection)
  type <- match.arg(type)
  
  if (missing(pre_k)) {
    stop("'pre_k' is required for diss_local_pls().")
  }
  .validate_pre_k(pre_k)
  .validate_shared_args(
    center = center, scale = scale,
    return_projection = return_projection,
    allow_parallel = allow_parallel
  )
  
  obj <- .new_diss_method(
    "pls",
    center = center,
    scale  = scale,
    extra  = list(
      ncomp               = as.integer(ncomp),
      ncomp_selection = ncomp_selection,
      type = type,
      local               = TRUE,
      pre_k               = as.integer(pre_k),
      return_projection   = return_projection,
      allow_parallel      = allow_parallel
    )
  )
  class(obj) <- c("diss_local_pls", class(obj))
  obj
}


# =============================================================================
# Shared validators — NOT exported
# =============================================================================

.validate_ncomp <- function(ncomp) {
  if (!is.numeric(ncomp) || length(ncomp) != 1L || is.na(ncomp) ||
      as.integer(ncomp) < 1L) {
    stop("'ncomp' must be a single positive integer.")
  }
  invisible(NULL)
}

.validate_ncomp_selection <- function(x) {
  if (!is.null(x) && !inherits(x, "nncomp_selection")) {
    stop(
      "'ncomp_selection' must be NULL or an nncomp_selection object. ",
      "Use ncomp_by_var(), ncomp_by_cumvar(), or ncomp_by_opc()."
    )
  }
  invisible(NULL)
}

.validate_pre_k <- function(pre_k) {
  if (!is.numeric(pre_k) || length(pre_k) != 1L || as.integer(pre_k) < 4L) {
    stop("'pre_k' must be a single integer >= 4.")
  }
  invisible(NULL)
}

.validate_shared_args <- function(center, scale, return_projection, allow_parallel) {
  if (!is.logical(center) || length(center) != 1L) {
    stop("'center' must be a single logical value (TRUE or FALSE).")
  }
  if (!is.logical(scale) || length(scale) != 1L) {
    stop("'scale' must be a single logical value (TRUE or FALSE).")
  }
  if (!is.logical(return_projection) || length(return_projection) != 1L) {
    stop("'return_projection' must be a single logical value (TRUE or FALSE).")
  }
  if (!is.logical(allow_parallel) || length(allow_parallel) != 1L) {
    stop("'allow_parallel' must be a single logical value (TRUE or FALSE).")
  }
  invisible(NULL)
}


# =============================================================================
# Print methods
# =============================================================================

#' @export
print.diss_pca <- function(x, ...) {
  cat("Dissimilarity method  : pca\n")
  cat(" algorithm            :", x$algorithm, "\n")
  cat(" ncomp                :", x$ncomp, "\n")
  if (!is.null(x$ncomp_selection)) {
    cat(" ncomp_selection  : "); print(x$ncomp_selection)
  } else {
    cat(" ncomp_selection  : none (fixed)\n")
  }
  cat(" center               :", x$center, "\n")
  cat(" scale                :", x$scale, "\n")
  cat(" return_projection    :", x$return_projection, "\n")
  invisible(x)
}

#' @export
print.diss_pls <- function(x, ...) {
  cat("Dissimilarity method  : pls\n")
  cat(" type                :", x$type, "\n")
  cat(" ncomp                :", x$ncomp, "\n")
  if (!is.null(x$ncomp_selection)) {
    cat(" ncomp_selection  : "); print(x$ncomp_selection)
  } else {
    cat(" ncomp_selection  : none (fixed)\n")
  }
  cat(" center               :", x$center, "\n")
  cat(" scale                :", x$scale, "\n")
  cat(" return_projection    :", x$return_projection, "\n")
  invisible(x)
}

# Not exported — only reachable via resemble:::diss_local_pca()
print.diss_local_pca <- function(x, ...) {
  cat("Dissimilarity method  : local pca (experimental)\n")
  cat(" algorithm            :", x$algorithm, "\n")
  cat(" ncomp                :", x$ncomp, "\n")
  if (!is.null(x$ncomp_selection)) {
    cat(" ncomp_selection  : "); print(x$ncomp_selection)
  } else {
    cat(" ncomp_selection  : none (fixed)\n")
  }
  cat(" pre_k                :", x$pre_k, "\n")
  cat(" center               :", x$center, "\n")
  cat(" scale                :", x$scale, "\n")
  cat(" return_projection    :", x$return_projection, "\n")
  invisible(x)
}

# Not exported — only reachable via resemble:::diss_local_pls()
print.diss_local_pls <- function(x, ...) {
  cat("Dissimilarity method  : local pls (experimental)\n")
  cat(" type                :", x$type, "\n")
  cat(" ncomp                :", x$ncomp, "\n")
  if (!is.null(x$ncomp_selection)) {
    cat(" ncomp_selection  : "); print(x$ncomp_selection)
  } else {
    cat(" ncomp_selection  : none (fixed)\n")
  }
  cat(" pre_k                :", x$pre_k, "\n")
  cat(" center               :", x$center, "\n")
  cat(" scale                :", x$scale, "\n")
  cat(" return_projection    :", x$return_projection, "\n")
  invisible(x)
}


# =============================================================================
# Internal compute function — NOT exported
# =============================================================================
#
# Called by dissimilarity() after data-level validation is done there.
# Assumes:
#   - Xr and Xu are numeric matrices with no NAs
#   - ncol(Xu) == ncol(Xr) if Xu is provided
#   - Yr is provided when method is diss_pls/diss_local_pls or
#     ncomp_selection is ncomp_by_opc (checked in dissimilarity())
#   - pre_k < nrow(Xr) when local = TRUE (checked in dissimilarity())
#
.ortho_diss_compute <- function(Xr, Xu = NULL, Yr = NULL, method) {
  
  stopifnot(inherits(method, "diss_method"))
  stopifnot(method$method %in% c("pca", "pls"))
  
  ncomp               <- method$ncomp
  ncomp_selection <- method$ncomp_selection
  center              <- method$center
  scale               <- method$scale
  local               <- method$local
  pre_k               <- method$pre_k
  return_projection   <- method$return_projection
  allow_parallel      <- method$allow_parallel
  
  # resolve projection method string for ortho_projection() backend
  proj_method <- if (method$method == "pca") {
    if (method$algorithm == "nipals") "pca.nipals" else "pca"
  } else {
    method$type  # "pls" or "mpls"
  }
  
  # bridge to the list format ortho_projection() expects:
  #   NULL ncomp_selection → "manual" with exact ncomp
  #   nncomp_selection object   → its method + ncomp as cap
  ncomp_list <- if (is.null(ncomp_selection)) {
    list(method = "manual", value = ncomp)
  } else {
    list(method = ncomp_selection$method, value = ncomp)
  }
  
  # --- global projection -----------------------------------------------------
  projection <- ortho_projection(
    Xr = Xr,
    Yr = Yr,
    Xu = Xu,
    method = proj_method,
    pc_selection = ncomp_list,
    center = center,
    scale = scale
  )
  
  scores <- projection$scores
  
  if (method$method == "pca") {
    scores      <- sweep(scores, MARGIN = 2, STATS = projection$scores_sd, FUN = "/")
    dist_method <- diss_euclidean(center = FALSE, scale = FALSE)
  } else {
    dist_method <- diss_mahalanobis(center = FALSE, scale = FALSE)
  }
  
  ncomp <- projection$n_components
  
  # --- global distance matrix ------------------------------------------------
  if (is.null(Xu)) {
    distnc <- .f_diss_compute(Xr = scores, Xu = NULL, method = dist_method)
    dimnames(distnc) <- list(rownames(scores), rownames(scores))
  } else {
    distnc <- .f_diss_compute(
      Xr     = scores[seq_len(nrow(Xr)), , drop = FALSE],
      Xu     = scores[(nrow(Xr) + 1L):nrow(scores), , drop = FALSE],
      method = dist_method
    )
    dimnames(distnc) <- list(
      rownames(scores[seq_len(nrow(Xr)), , drop = FALSE]),
      rownames(scores[(nrow(Xr) + 1L):nrow(scores), , drop = FALSE])
    )
  }
  
  # --- local dissimilarities -------------------------------------------------
  if (local) {
    return(.ortho_diss_local(
      distnc            = distnc,
      Xr                = Xr,
      Xu                = Xu,
      Yr                = Yr,
      pre_k             = pre_k,
      proj_method       = proj_method,
      ncomp_list        = ncomp_list,
      dist_method       = dist_method,
      center            = center,
      scale             = scale,
      allow_parallel    = allow_parallel,
      ncomp             = ncomp,
      projection        = projection,
      return_projection = return_projection
    ))
  }
  
  # --- non-local output ------------------------------------------------------
  out <- list(
    ncomp                = ncomp,
    global_variance_info = projection$variance,
    dissimilarity        = distnc
  )
  if (return_projection) out$projection <- projection
  class(out) <- c("ortho_diss", "list")
  class(out$dissimilarity) <- c("ortho_diss", "matrix")
  out
}


# --- local dissimilarity helper — NOT exported -------------------------------
.ortho_diss_local <- function(
    distnc, Xr, Xu, Yr, pre_k,
    proj_method, ncomp_list, dist_method,
    center, scale, allow_parallel,
    ncomp, projection, return_projection
) {
  n_xr  <- nrow(Xr)
  pre_k <- pre_k + 1L  # include target observation in its own neighborhood indices
  
  if (!is.null(Xu)) {
    Xr <- rbind(Xr, Xu)
    Yr <- rbind(Yr, matrix(NA_real_, nrow(Xu), ncol(Yr)))
  }
  
  neighborhood_info <- k0_indices <- apply(distnc, MARGIN = 2, FUN = order)[1:pre_k, ]
  d_dimnames <- dimnames(distnc)
  rm(distnc)
  
  neighborhood_info[k0_indices <= n_xr] <- paste0("Xr_", neighborhood_info[k0_indices <= n_xr])
  neighborhood_info[k0_indices > n_xr]  <- paste0("Xu_", neighborhood_info[k0_indices > n_xr])
  
  neighborhood_info <- data.table(
    do.call("rbind", strsplit(colnames(k0_indices), "_")),
    t(neighborhood_info)
  )
  colnames(neighborhood_info) <- c(
    "Set", "Index",
    paste0("k_", seq_len(nrow(k0_indices)))
  )
  
  local_d <- local_ortho_diss(
    k_index_matrix = k0_indices,
    Xr             = Xr,
    Yr             = Yr,
    Xu             = Xu,
    diss_method    = proj_method,
    pc_selection   = ncomp_list,
    center         = center,
    scale          = scale,
    allow_parallel = allow_parallel
  )
  
  dimnames(local_d$dissimilarity_m) <- d_dimnames
  neighborhood_info <- cbind(
    neighborhood_info[, 2:1],
    local_ncomp = local_d$local_ncomp,
    neighborhood_info[, -c(1:2)]
  )
  
  out <- list(
    ncomp                = ncomp,
    global_variance_info = projection$variance,
    neighborhood_info    = neighborhood_info,
    dissimilarity        = local_d$dissimilarity_m
  )
  if (return_projection) out$projection <- projection
  class(out) <- c("ortho_diss", "list")
  class(out$dissimilarity) <- c("local_ortho_diss", "matrix")
  out
}
