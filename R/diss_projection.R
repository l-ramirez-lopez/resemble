#' @title PCA dissimilarity method constructor
#'
#' @description
#' Creates a configuration object for computing dissimilarity based on
#' Mahalanobis distance in PCA score space.
#'
#' @param ncomp Component selection. A positive integer or an
#'   \code{ncomp_selection} object. Default \code{ncomp_by_var(0.01)}.
#' @param method Character. PCA method: \code{"pca"} (default, SVD-based) or
#'   \code{"pca_nipals"} (NIPALS algorithm).
#' @param center Logical. Center data before projection? Default \code{TRUE}.
#' @param scale Logical. Scale data before projection? Default \code{FALSE}.
#' @param return_projection Logical. Return the projection object?
#'   Default \code{FALSE}.
#'
#' @return An object of class \code{c("diss_pca", "diss_method")}.
#'
#' @seealso \code{\link{diss_pls}}, \code{\link{dissimilarity}},
#'   \code{\link{ncomp_by_var}}
#'
#' @examples
#' # Default: SVD-based PCA
#' m <- diss_pca(ncomp = 10)
#'
#' # NIPALS algorithm
#' m <- diss_pca(ncomp = 10, method = "pca_nipals")
#'
#' # Automatic component selection
#' m <- diss_pca(ncomp = ncomp_by_var(0.01))
#'
#' @export
diss_pca <- function(
    ncomp = ncomp_by_var(0.01),
    method = c("pca", "pca_nipals"),
    center = TRUE,
    scale = FALSE,
    return_projection = FALSE
) {
  ncomp <- .coerce_ncomp(ncomp)
  method <- match.arg(method)
  .validate_logical_args(center, scale, return_projection)
  
  .new_diss_method(
    "pca",
    center = center,
    scale = scale,
    extra = list(
      ncomp = ncomp,
      method = method,
      local = FALSE,
      pre_k = NULL,
      return_projection = return_projection
    )
  )
}


#' @title PLS dissimilarity method constructor
#'
#' @description
#' Creates a configuration object for computing dissimilarity based on
#' Mahalanobis distance in PLS score space. Requires \code{Yr} in
#' \code{dissimilarity()}.
#'
#' @param ncomp Component selection. A positive integer or an
#'   \code{ncomp_selection} object. Default \code{ncomp_by_opc(40)}.
#' @param method Character. PLS method: \code{"pls"} (default) or
#'   \code{"mpls"} (modified PLS, Shenk & Westerhaus 1991).
#' @param center Logical. Center data? Default \code{TRUE}.
#' @param scale Logical. Scale data? Default \code{FALSE}.
#' @param return_projection Logical. Return projection object?
#'   Default \code{FALSE}.
#'
#' @return An object of class \code{c("diss_pls", "diss_method")}.
#'
#' @seealso \code{\link{diss_pca}}, \code{\link{dissimilarity}}
#'
#' @examples
#' # Standard PLS
#' m <- diss_pls(ncomp = 15)
#'
#' # Modified PLS
#' m <- diss_pls(ncomp = 10, method = "mpls")
#'
#' @export
diss_pls <- function(
    ncomp = ncomp_by_opc(40),
    method = c("pls", "mpls"),
    center = TRUE,
    scale = FALSE,
    return_projection = FALSE
) {
  ncomp <- .coerce_ncomp(ncomp)
  method <- match.arg(method)
  .validate_logical_args(center, scale, return_projection)
  
  .new_diss_method(
    "pls",
    center = center,
    scale = scale,
    extra = list(
      ncomp = ncomp,
      method = method,
      local = FALSE,
      pre_k = NULL,
      return_projection = return_projection
    )
  )
}


# =============================================================================
# Print methods
# =============================================================================

#' @export
print.diss_pca <- function(x, ...) {
  cat("Dissimilarity: pca\n")
  cat("  method             :", x$method, "\n")
  cat("  ncomp              : ")
  print(x$ncomp)
  cat("  center             :", x$center, "\n")
  cat("  scale              :", x$scale, "\n")
  cat("  return_projection  :", x$return_projection, "\n")
  invisible(x)
}

#' @export
print.diss_pls <- function(x, ...) {
  cat("Dissimilarity: pls\n")
  cat("  method             :", x$method, "\n")
  cat("  ncomp              : ")
  print(x$ncomp)
  cat("  center             :", x$center, "\n")
  cat("  scale              :", x$scale, "\n")
  cat("  return_projection  :", x$return_projection, "\n")
  invisible(x)
}

# =============================================================================
# Internal helpers
# =============================================================================

.coerce_ncomp <- function(ncomp) {
  if (is.numeric(ncomp) && length(ncomp) == 1L) {
    return(ncomp_fixed(as.integer(ncomp)))
  }
  if (!inherits(ncomp, "ncomp_selection")) {
    stop(
      "'ncomp' must be a positive integer or an ncomp_*() object.\n",
      "See ?ncomp_by_var, ?ncomp_by_cumvar, ?ncomp_by_opc, ?ncomp_fixed."
    )
  }
  ncomp
}

.validate_logical_args <- function(center, scale, return_projection) {
  if (!is.logical(center) || length(center) != 1L || is.na(center)) {
    stop("'center' must be TRUE or FALSE.")
  }
  if (!is.logical(scale) || length(scale) != 1L || is.na(scale)) {
    stop("'scale' must be TRUE or FALSE.")
  }
  if (!is.logical(return_projection) || length(return_projection) != 1L ||
    is.na(return_projection)) {
    stop("'return_projection' must be TRUE or FALSE.")
  }
  invisible(NULL)
}


# =============================================================================
# diss_local_pca — NOT EXPORTED
# =============================================================================

diss_local_pca <- function(
    ncomp = ncomp_by_var(0.01),
    method = c("pca", "pca_nipals"),
    pre_k,
    center = TRUE,
    scale = FALSE,
    return_projection = FALSE,
    allow_parallel = TRUE
) {
  ncomp <- .coerce_ncomp(ncomp)
  method <- match.arg(method)
  
  if (missing(pre_k)) {
    stop("'pre_k' is required for diss_local_pca().")
  }
  .validate_pre_k(pre_k)
  .validate_logical_args(center, scale, return_projection, allow_parallel)
  
  obj <- .new_diss_method(
    "pca",
    center = center,
    scale = scale,
    extra = list(
      ncomp = ncomp,
      method = method,
      local = TRUE,
      pre_k = as.integer(pre_k),
      return_projection = return_projection,
      allow_parallel = allow_parallel
    )
  )
  class(obj) <- c("diss_local_pca", class(obj))
  obj
}


# =============================================================================
# diss_local_pls — NOT EXPORTED
# =============================================================================

diss_local_pls <- function(
    ncomp = ncomp_by_opc(40),
    method = c("pls", "mpls"),
    pre_k,
    center = TRUE,
    scale = FALSE,
    return_projection = FALSE,
    allow_parallel = TRUE
) {
  ncomp <- .coerce_ncomp(ncomp)
  method <- match.arg(method)
  
  if (missing(pre_k)) {
    stop("'pre_k' is required for diss_local_pls().")
  }
  .validate_pre_k(pre_k)
  .validate_logical_args(center, scale, return_projection, allow_parallel)
  
  obj <- .new_diss_method(
    "pls",
    center = center,
    scale = scale,
    extra = list(
      ncomp = ncomp,
      method = method,
      local = TRUE,
      pre_k = as.integer(pre_k),
      return_projection = return_projection,
      allow_parallel = allow_parallel
    )
  )
  class(obj) <- c("diss_local_pls", class(obj))
  obj
}


# =============================================================================
# Print methods (not exported)
# =============================================================================

print.diss_local_pca <- function(x, ...) {
  cat("Dissimilarity: local pca (experimental)\n")
  cat("  method             :", x$method, "\n")
  cat("  ncomp              : ")
  print(x$ncomp)
  cat("  pre_k              :", x$pre_k, "\n")
  cat("  center             :", x$center, "\n")
  cat("  scale              :", x$scale, "\n")
  cat("  return_projection  :", x$return_projection, "\n")
  invisible(x)
}

print.diss_local_pls <- function(x, ...) {
  cat("Dissimilarity: local pls (experimental)\n")
  cat("  method             :", x$method, "\n")
  cat("  ncomp              : ")
  print(x$ncomp)
  cat("  pre_k              :", x$pre_k, "\n")
  cat("  center             :", x$center, "\n")
  cat("  scale              :", x$scale, "\n")
  cat("  return_projection  :", x$return_projection, "\n")
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
#     ncomp$method is "opc" (checked in dissimilarity())
#   - pre_k < nrow(Xr) when local = TRUE (checked in dissimilarity())
#
.ortho_diss_compute <- function(Xr, Xu = NULL, Yr = NULL, method) {
  stopifnot(inherits(method, "diss_method"))
  stopifnot(method$method %in% c("pca", "pls"))

  ncomp_obj <- method$ncomp
  center <- method$center
  scale <- method$scale
  local <- method$local
  pre_k <- method$pre_k
  return_projection <- method$return_projection
  allow_parallel <- method$allow_parallel

  # resolve projection method string for ortho_projection() backend
  proj_method <- if (method$method == "pca") {
    if (method$algorithm == "nipals") "pca.nipals" else "pca"
  } else {
    method$type
  }

  # bridge ncomp_selection object to ortho_projection()'s pc_selection format
  ncomp_list <- .ncomp_to_pc_selection(ncomp_obj)

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
    scores <- sweep(scores, MARGIN = 2, STATS = projection$scores_sd, FUN = "/")
    dist_method <- diss_euclidean(center = FALSE, scale = FALSE)
  } else {
    dist_method <- diss_mahalanobis(center = FALSE, scale = FALSE)
  }

  n_components <- projection$n_components

  # --- global distance matrix ------------------------------------------------
  if (is.null(Xu)) {
    distnc <- .f_diss_compute(Xr = scores, Xu = NULL, method = dist_method)
    dimnames(distnc) <- list(rownames(scores), rownames(scores))
  } else {
    distnc <- .f_diss_compute(
      Xr = scores[seq_len(nrow(Xr)), , drop = FALSE],
      Xu = scores[(nrow(Xr) + 1L):nrow(scores), , drop = FALSE],
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
      distnc = distnc,
      Xr = Xr,
      Xu = Xu,
      Yr = Yr,
      pre_k = pre_k,
      proj_method = proj_method,
      ncomp_list = ncomp_list,
      dist_method = dist_method,
      center = center,
      scale = scale,
      allow_parallel = allow_parallel,
      n_components = n_components,
      projection = projection,
      return_projection = return_projection
    ))
  }

  # --- non-local output ------------------------------------------------------
  out <- list(
    n_components = n_components,
    global_variance_info = projection$variance,
    dissimilarity = distnc
  )
  if (return_projection) out$projection <- projection
  class(out) <- c("ortho_diss", "list")
  class(out$dissimilarity) <- c("ortho_diss", "matrix")
  out
}


# =============================================================================
# Bridge ncomp_selection to pc_selection list format
# =============================================================================
#
# ortho_projection() expects pc_selection = list(method = "...", value = N)
# This converts the new ncomp_selection objects to that format.
#
.ncomp_to_pc_selection <- function(ncomp_obj) {
  stopifnot(inherits(ncomp_obj, "ncomp_selection"))

  switch(class(ncomp_obj)[[1]],
    ncomp_fixed = list(
      method = "manual",
      value = ncomp_obj$ncomp
    ),
    ncomp_by_var = list(
      method = "var",
      value = ncomp_obj$min_var,
      max_ncomp = ncomp_obj$max_ncomp
    ),
    ncomp_by_cumvar = list(
      method = "cumvar",
      value = ncomp_obj$min_cumvar,
      max_ncomp = ncomp_obj$max_ncomp
    ),
    ncomp_by_opc = list(
      method = "opc",
      value = ncomp_obj$max_ncomp
    ),
    stop("Unknown ncomp_selection class: ", class(ncomp_obj)[[1]])
  )
}


# =============================================================================
# Local dissimilarity helper — NOT exported
# =============================================================================

.ortho_diss_local <- function(
  distnc, Xr, Xu, Yr, pre_k,
  proj_method, ncomp_list, dist_method,
  center, scale, allow_parallel,
  n_components, projection, return_projection
) {
  n_xr <- nrow(Xr)
  pre_k <- pre_k + 1L

  if (!is.null(Xu)) {
    Xr <- rbind(Xr, Xu)
    Yr <- rbind(Yr, matrix(NA_real_, nrow(Xu), ncol(Yr)))
  }

  neighborhood_info <- k0_indices <- apply(distnc, MARGIN = 2, FUN = order)[1:pre_k, ]
  d_dimnames <- dimnames(distnc)
  rm(distnc)

  neighborhood_info[k0_indices <= n_xr] <- paste0("Xr_", neighborhood_info[k0_indices <= n_xr])
  neighborhood_info[k0_indices > n_xr] <- paste0("Xu_", neighborhood_info[k0_indices > n_xr])

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
    Xr = Xr,
    Yr = Yr,
    Xu = Xu,
    diss_method = proj_method,
    pc_selection = ncomp_list,
    center = center,
    scale = scale,
    allow_parallel = allow_parallel
  )

  dimnames(local_d$dissimilarity_m) <- d_dimnames
  neighborhood_info <- cbind(
    neighborhood_info[, 2:1],
    local_n_components = local_d$local_n_components,
    neighborhood_info[, -c(1:2)]
  )

  out <- list(
    n_components = n_components,
    global_variance_info = projection$variance,
    neighborhood_info = neighborhood_info,
    dissimilarity = local_d$dissimilarity_m
  )
  if (return_projection) out$projection <- projection
  class(out) <- c("ortho_diss", "list")
  class(out$dissimilarity) <- c("local_ortho_diss", "matrix")
  out
}
