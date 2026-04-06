#' @title PCA dissimilarity method constructor
#'
#' @description
#' Creates a configuration object for computing dissimilarity based on
#' Mahalanobis distance in PCA score space.
#'
#' @param ncomp Component selection method. Can be:
#'   \itemize{
#'     \item A positive integer for a fixed number of components
#'     \item \code{\link{ncomp_fixed}(n)}: explicit fixed selection
#'     \item \code{\link{ncomp_by_var}(min_var)}: retain components explaining
#'       at least \code{min_var} variance each (default: \code{ncomp_by_var(0.01)})
#'     \item \code{\link{ncomp_by_cumvar}(min_cumvar)}: retain components until
#'       cumulative variance reaches \code{min_cumvar}
#'     \item \code{\link{ncomp_by_opc}()}: optimize using side information
#'       (\code{Yr} required in \code{dissimilarity()})
#'   }
#' @param method Character. PCA algorithm: \code{"pca"} (default, SVD-based) or
#'   \code{"pca_nipals"} (NIPALS algorithm).
#' @param center Logical. Center data before projection? Default \code{TRUE}.
#' @param scale Logical. Scale data before projection? Default \code{FALSE}.
#' @param return_projection Logical. Return the projection object?
#'   Default \code{FALSE}.
#'
#' @return An object of class \code{c("diss_pca", "diss_method")}.
#'
#' @seealso 
#' Component selection: \code{\link{ncomp_by_var}}, \code{\link{ncomp_by_cumvar}},
#'   \code{\link{ncomp_by_opc}}, \code{\link{ncomp_fixed}}
#'
#' Other dissimilarity methods: \code{\link{diss_pls}}, 
#'   \code{\link{diss_correlation}}, \code{\link{diss_euclidean}},
#'   \code{\link{diss_cosine}}, \code{\link{diss_mahalanobis}}
#'
#' @examples
#' # Fixed number of components
#' diss_pca(ncomp = 10)
#' diss_pca(ncomp = ncomp_fixed(10))
#'
#' # Retain components explaining >= 1% variance each (default)
#' diss_pca(ncomp = ncomp_by_var(0.01))
#'
#' # Retain components until 99% cumulative variance
#' diss_pca(ncomp = ncomp_by_cumvar(0.99))
#'
#' # Optimize using side information (requires Yr)
#' diss_pca(ncomp = ncomp_by_opc(40))
#' diss_pca(ncomp = ncomp_by_opc())
#'
#' # NIPALS algorithm (useful for very large matrices)
#' diss_pca(ncomp = 10, method = "pca_nipals")
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
  
  structure(
    list(
      ncomp = ncomp,
      method = method,
      center = center,
      scale = scale,
      return_projection = return_projection
    ),
    class = c("diss_pca", "diss_method")
  )
}


#' @title PLS dissimilarity method constructor
#'
#' @description
#' Creates a configuration object for computing dissimilarity based on
#' Mahalanobis distance in PLS score space. Requires \code{Yr} in
#' \code{dissimilarity()}.
#'
#' @param ncomp Component selection method. Can be:
#'   \itemize{
#'     \item A positive integer for a fixed number of components
#'     \item \code{\link{ncomp_fixed}(n)}: explicit fixed selection
#'     \item \code{\link{ncomp_by_var}(min_var)}: retain components explaining
#'       at least \code{min_var} variance each
#'     \item \code{\link{ncomp_by_cumvar}(min_cumvar)}: retain components until
#'       cumulative variance reaches \code{min_cumvar}
#'     \item \code{\link{ncomp_by_opc}()}: optimize using side information
#'       (default; recommended for PLS since \code{Yr} is already required)
#'   }
#' @param method Character. PLS algorithm: \code{"pls"} (default) or
#'   \code{"mpls"} (modified PLS, Shenk & Westerhaus 1991).
#' @param scale Logical. Scale data? Default \code{FALSE}. Note: PLS always
#'   centers internally.
#' @param return_projection Logical. Return projection object?
#'   Default \code{FALSE}.
#'
#' @return An object of class \code{c("diss_pls", "diss_method")}.
#'
#' @author \href{https://orcid.org/0000-0002-5369-5120}{Leonardo Ramirez-Lopez}
#'
#' @seealso 
#' Component selection: \code{\link{ncomp_by_var}}, \code{\link{ncomp_by_cumvar}},
#'   \code{\link{ncomp_by_opc}}, \code{\link{ncomp_fixed}}
#'
#' Other dissimilarity methods: \code{\link{diss_pca}},
#'   \code{\link{diss_correlation}}, \code{\link{diss_euclidean}},
#'   \code{\link{diss_cosine}}, \code{\link{diss_mahalanobis}}
#'
#' @examples
#' # Default: OPC optimization (recommended)
#' diss_pls()
#'
#' # Fixed number of components
#' diss_pls(ncomp = 15)
#'
#' # Custom opc settings
#' diss_pls(ncomp = ncomp_by_opc(max_ncomp = 50))
#'
#' # Modified PLS
#' diss_pls(ncomp = 10, method = "mpls")
#'
#' @export
diss_pls <- function(
    ncomp = ncomp_by_opc(),
    method = c("pls", "mpls"),
    scale = FALSE,
    return_projection = FALSE
) {
  ncomp <- .coerce_ncomp(ncomp)
  method <- match.arg(method)
  .validate_logical_args(TRUE, scale, return_projection)
  
  structure(
    list(
      ncomp = ncomp,
      method = method,
      center = TRUE,
      scale = scale,
      return_projection = return_projection
    ),
    class = c("diss_pls", "diss_method")
  )
}


.format_ncomp <- function(x) {
  if (is.numeric(x)) {
    return(sprintf("fixed: %d", as.integer(x)))
  }
  
  ncomp_class <- class(x)[[1]]
  
  switch(
    ncomp_class,
    ncomp_by_var = sprintf("var >= %s (max: %d)", x$min_var, x$max_ncomp),
    ncomp_by_cumvar = sprintf("cumvar >= %s (max: %d)", x$min_cumvar, x$max_ncomp),
    ncomp_by_opc = sprintf("opc (%.0f%%, max: %d)", x$prop * 100, x$max_ncomp),
    ncomp_fixed = sprintf("fixed: %d", x$ncomp)
  )
}


#' @export
print.diss_pca <- function(x, ...) {
  cat("Dissimilarity: PCA\n")
  cat("  method            :", x$method, "\n")
  cat("  ncomp             :", .format_ncomp(x$ncomp), "\n")
  cat("  center            :", x$center, "\n")
  cat("  scale             :", x$scale, "\n")
  cat("  return_projection :", x$return_projection, "\n")
  invisible(x)
}


#' @export
print.diss_pls <- function(x, ...) {
  cat("Dissimilarity: PLS\n")
  cat("  method            :", x$method, "\n")
  cat("  ncomp             :", .format_ncomp(x$ncomp), "\n")
  cat("  scale             :", x$scale, "\n")
  cat("  return_projection :", x$return_projection, "\n")
  invisible(x)
}


# =============================================================================
# Internal helpers
# =============================================================================

.validate_logical_args <- function(...) {
  args <- list(...)
  nms <- as.character(substitute(list(...)))[-1]
  
  for (i in seq_along(args)) {
    val <- args[[i]]
    if (!is.logical(val) || length(val) != 1L || is.na(val)) {
      stop("'", nms[i], "' must be TRUE or FALSE.", call. = FALSE)
    }
  }
  invisible(NULL)
}


.validate_pre_k <- function(pre_k) {
  if (!is.numeric(pre_k) || length(pre_k) != 1L || pre_k < 1L) {
    stop("'pre_k' must be a positive integer.", call. = FALSE)
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
    stop("'pre_k' is required for diss_local_pca().", call. = FALSE)
  }
  .validate_pre_k(pre_k)
  .validate_logical_args(center, scale, return_projection, allow_parallel)
  
  structure(
    list(
      ncomp = ncomp,
      method = method,
      center = center,
      scale = scale,
      local = TRUE,
      pre_k = as.integer(pre_k),
      return_projection = return_projection,
      allow_parallel = allow_parallel
    ),
    class = c("diss_local_pca", "diss_pca", "diss_method")
  )
}


# =============================================================================
# diss_local_pls — NOT EXPORTED
# =============================================================================

diss_local_pls <- function(
    ncomp = ncomp_by_opc(),
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
    stop("'pre_k' is required for diss_local_pls().", call. = FALSE)
  }
  .validate_pre_k(pre_k)
  .validate_logical_args(center, scale, return_projection, allow_parallel)
  
  structure(
    list(
      ncomp = ncomp,
      method = method,
      center = TRUE,
      scale = scale,
      local = TRUE,
      pre_k = as.integer(pre_k),
      return_projection = return_projection,
      allow_parallel = allow_parallel
    ),
    class = c("diss_local_pls", "diss_pls", "diss_method")
  )
}


# =============================================================================
# Print methods for local (not exported)
# =============================================================================

#' @export
print.diss_local_pca <- function(x, ...) {
  cat("Dissimilarity: local PCA\n")
  cat("  method            :", x$method, "\n")
  cat("  ncomp             :", .format_ncomp(x$ncomp), "\n")
  cat("  pre_k             :", x$pre_k, "\n")
  cat("  center            :", x$center, "\n")
  cat("  scale             :", x$scale, "\n")
  cat("  return_projection :", x$return_projection, "\n")
  invisible(x)
}


#' @export
print.diss_local_pls <- function(x, ...) {
  cat("Dissimilarity: local PLS\n")
  cat("  method            :", x$method, "\n")
  cat("  ncomp             :", .format_ncomp(x$ncomp), "\n")
  cat("  pre_k             :", x$pre_k, "\n")
  cat("  scale             :", x$scale, "\n")
  cat("  return_projection :", x$return_projection, "\n")
  invisible(x)
}



