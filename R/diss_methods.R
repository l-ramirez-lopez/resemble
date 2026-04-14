# =============================================================================
# Euclidean, Mahalanobis and cosine dissimilarity — redesigned API
# =============================================================================
# Public:   diss_euclidean()    — method constructor
#           diss_mahalanobis()  — method constructor
#           diss_cosine()       — method constructor
# Internal: .f_diss_compute()   — shared raw computation, not exported


# -----------------------------------------------------------------------------
#' @title Euclidean dissimilarity method constructor
#'
#' @description
#' Creates a configuration object for computing Euclidean dissimilarity.
#' Pass the result to \code{dissimilarity()} to compute the dissimilarity
#' matrix.
#'
#' The scaled Euclidean dissimilarity between two observations \eqn{x_i} and
#' \eqn{x_j} is:
#'
#' \deqn{d(x_i, x_j) = \sqrt{\frac{1}{p} \sum_{k=1}^{p}(x_{i,k} - x_{j,k})^2}}
#'
#' where \eqn{p} is the number of variables. Results are equivalent to
#' \code{stats::dist()} but scaled by \eqn{1/p}.
#'
#' @param center Logical. Center the data before computing distances?
#'   Applied jointly to \code{Xr} and \code{Xu} if both are provided.
#'   Default \code{TRUE}.
#' @param scale Logical. Scale the data before computing distances?
#'   Applied jointly to \code{Xr} and \code{Xu} if both are provided.
#'   Default \code{FALSE}.
#'
#' @return An object of class \code{c("diss_euclidean", "diss_method")}.
#' @seealso \code{\link{dissimilarity}}, \code{\link{diss_mahalanobis}},
#'   \code{\link{diss_cosine}}
#' @examples
#' m <- diss_euclidean()
#' m <- diss_euclidean(center = FALSE, scale = TRUE)
#' @export
diss_euclidean <- function(center = TRUE, scale = FALSE) {
  .new_diss_method("euclidean", center = center, scale = scale)
}


# -----------------------------------------------------------------------------
#' @title Mahalanobis dissimilarity method constructor
#'
#' @description
#' Creates a configuration object for computing Mahalanobis dissimilarity.
#' Pass the result to \code{dissimilarity()} to compute the dissimilarity
#' matrix.
#'
#' The Mahalanobis distance is computed by first transforming the data into
#' Mahalanobis space via a factorization of the inverse covariance matrix
#' \eqn{M^{-1} = W^{T}W} (using SVD), then applying Euclidean distance in
#' that transformed space:
#'
#' \deqn{d(x_i, x_j) = \sqrt{\frac{1}{p}(x_i - x_j)M^{-1}(x_i - x_j)^T}}
#'
#' @section Important limitations:
#' The covariance matrix will be singular — and the distance therefore
#' uncomputable — when the number of observations is smaller than the number
#' of variables, or when variables are perfectly collinear. This is common
#' with raw spectral data; consider using \code{diss_euclidean()} on
#' PCA scores instead.
#'
#' @param center Logical. Center the data before computing distances?
#'   Applied jointly to \code{Xr} and \code{Xu} if both are provided.
#'   Default \code{TRUE}.
#' @param scale Logical. Scale the data before computing distances?
#'   Applied jointly to \code{Xr} and \code{Xu} if both are provided.
#'   Default \code{FALSE}.
#'
#' @return An object of class \code{c("diss_mahalanobis", "diss_method")}.
#' @seealso \code{\link{dissimilarity}}, \code{\link{diss_euclidean}},
#'   \code{\link{diss_cosine}}
#' @examples
#' m <- diss_mahalanobis()
#' m <- diss_mahalanobis(center = TRUE, scale = TRUE)
#' @export
diss_mahalanobis <- function(center = TRUE, scale = FALSE) {
  .new_diss_method("mahalanobis", center = center, scale = scale)
}


# -----------------------------------------------------------------------------
#' @title Cosine dissimilarity method constructor
#'
#' @description
#' Creates a configuration object for computing cosine dissimilarity
#' (also known as spectral angle mapper). Pass the result to
#' \code{dissimilarity()} to compute the dissimilarity matrix.
#'
#' The cosine dissimilarity between two observations \eqn{x_i} and
#' \eqn{x_j} is:
#'
#' \deqn{c(x_i, x_j) = \cos^{-1}
#'   \frac{\sum_{k=1}^{p} x_{i,k}\, x_{j,k}}
#'        {\sqrt{\sum_{k=1}^{p} x_{i,k}^{2}}\,
#'         \sqrt{\sum_{k=1}^{p} x_{j,k}^{2}}}}
#'
#' where \eqn{p} is the number of variables.
#'
#' @param center Logical. Center the data before computing dissimilarities?
#'   Applied jointly to \code{Xr} and \code{Xu} if both are provided.
#'   Default \code{TRUE}.
#' @param scale Logical. Scale the data before computing dissimilarities?
#'   Applied jointly to \code{Xr} and \code{Xu} if both are provided.
#'   Default \code{FALSE}.
#'
#' @return An object of class \code{c("diss_cosine", "diss_method")}.
#' @seealso \code{\link{dissimilarity}}, \code{\link{diss_euclidean}},
#'   \code{\link{diss_mahalanobis}}
#' @examples
#' m <- diss_cosine()
#' m <- diss_cosine(center = FALSE)
#' @export
diss_cosine <- function(center = TRUE, scale = FALSE) {
  .new_diss_method("cosine", center = center, scale = scale)
}




# =============================================================================
# Internal shared compute function — NOT exported
# =============================================================================
#
# Called by dissimilarity() after data-level validation is done there.
# Assumes:
#   - Xr and Xu are numeric matrices with no NAs
#   - ncol(Xu) == ncol(Xr) if Xu is provided
#   - for Mahalanobis: nrow(rbind(Xr, Xu)) > ncol(Xr) — checked in dissimilarity()
#
.f_diss_compute <- function(Xr, Xu = NULL, method) {
  
  stopifnot(inherits(method, "diss_method"))
  stopifnot(method$method %in% c("euclidean", "mahalanobis", "cosine"))
  
  center      <- method$center
  scale       <- method$scale
  method_name <- method$method
  
  # --- internal method label used for C++ dispatch --------------------------
  # Mahalanobis is transformed to Euclidean space before calling fast_diss(),
  # so the C++ backend always sees either "euclid" or "cosine"
  cpp_method <- if (method_name == "cosine") "cosine" else "euclid"
  
  # --- preprocessing --------------------------------------------------------
  # Note: Euclidean and Mahalanobis always go through rbind() because
  # Mahalanobis needs the full combined X to estimate the covariance matrix.
  # Cosine only needs it if center or scale are requested.
  needs_combined <- center || scale || method_name %in% c("euclidean", "mahalanobis")
  
  if (needs_combined) {
    X <- rbind(Xr, Xu)  # rbind(x, NULL) == x, no special casing needed
    
    if (center) {
      X <- sweep(X, MARGIN = 2, STATS = colMeans(X), FUN = "-")
    }
    if (scale) {
      X <- sweep(X, MARGIN = 2, STATS = get_col_sds(X), FUN = "/")
    }
    
    # --- Mahalanobis transform -----------------------------------------------
    # Project X into Mahalanobis space; afterwards treat as Euclidean
    if (method_name == "mahalanobis") {
      X <- try(euclid_to_mahal(X, sm_method = "svd"), silent = TRUE)
      if (!is.matrix(X)) {
        stop(
          "The covariance matrix is exactly singular and cannot be inverted. ",
          "The Mahalanobis distance cannot be computed. ",
          "Consider using diss_euclidean() on PCA scores instead."
        )
      }
    }
    
    if (!is.null(Xu)) {
      n_xu <- nrow(Xu)
      Xu   <- X[(nrow(X) - n_xu + 1L):nrow(X), , drop = FALSE]
      Xr   <- X[1L:(nrow(X) - n_xu), , drop = FALSE]
    } else {
      Xr <- X
    }
    rm(X)
  }
  
  # --- compute dissimilarity ------------------------------------------------
  if (!is.null(Xu)) {
    # abs() guard: C++ backend can return small negative values (~-1e-14)
    # due to floating point reuse in memory
    result <- abs(fast_diss(Xu, Xr, cpp_method))
    if (cpp_method == "euclid") {
      result <- sqrt(result / ncol(Xr))
    }
    rownames(result) <- paste("Xr", seq_len(nrow(Xr)), sep = "_")
    colnames(result) <- paste("Xu", seq_len(nrow(Xu)), sep = "_")
  } else {
    result <- abs(fast_diss(Xr, Xr, cpp_method))
    if (cpp_method == "euclid") {
      result <- sqrt(result / ncol(Xr))
    }
    rownames(result) <- colnames(result) <- paste("Xr", seq_len(nrow(Xr)), sep = "_")
  }
  
  # --- cosine NaN guard ------------------------------------------------------
  # The C++ backend occasionally returns NaN for identical observations
  # due to floating point precision in the arccos computation.
  # These are self-distances and are correctly zero.
  if (method_name == "cosine") {
    result[is.nan(result)] <- 0
  }
  
  result
}