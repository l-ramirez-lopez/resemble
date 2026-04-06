# =============================================================================
# Neighbor selection constructors
# =============================================================================
# Public:
#   neighbors_k()
#   neighbors_diss()
#
# Internal:
#   .new_neighbors()     - shared constructor
#   print.neighbors_k
#   print.neighbors_diss
#   print.neighbors      - fallback


#' @title Neighbor selection methods
#'
#' @description
#' These functions create configuration objects that specify how neighbors
#' are selected for memory-based learning in \code{\link{mbl}} an \code{\link{liblex}}.
#'
#' @name neighbors
#' @aliases neighbors_k neighbors_diss
#'
#' @usage
#' neighbors_k(k)
#'
#' neighbors_diss(threshold, k_min = 4L, k_max = Inf)
#'
#' @param k Integer vector. One or more neighborhood sizes to evaluate.
#'   Values will be sorted in ascending order. Minimum allowed value is 4.
#' @param threshold Numeric vector. One or more dissimilarity thresholds.
#'   Neighbors are selected if their dissimilarity to the target observation
#'   is below the threshold. Values will be sorted in ascending order.
#' @param k_min Integer. Minimum number of neighbors to retain, regardless of
#'   threshold. Default \code{4L}.
#' @param k_max Integer or \code{Inf}. Maximum number of neighbors to retain,
#'   regardless of threshold. Default \code{Inf} (no upper bound other than
#'   the size of the reference set).
#'
#' @details
#' Two strategies are available for neighbor selection:
#'
#' \strong{Fixed-k selection} (\code{neighbors_k}) 
#' 
#' A fixed number of nearest neighbors is selected for each target observation. 
#' Multiple values of \code{k} can be provided to evaluate different 
#' neighborhood sizes.
#'
#' \strong{Dissimilarity-threshold selection} (\code{neighbors_diss}) 
#' 
#' Neighbors are selected based on a dissimilarity threshold. All reference 
#' observations with dissimilarity below the threshold are included. The 
#' \code{k_min} and \code{k_max} arguments provide bounds to ensure a reasonable 
#' neighborhood size regardless of the threshold. Multiple thresholds can be 
#' provided to evaluate different settings.
#'
#' @return
#' An object of class \code{c("neighbors_k", "neighbors")} or
#' \code{c("neighbors_diss", "neighbors")}, containing the validated
#' parameters. Intended to be passed to \code{\link{mbl}}.
#'
#' @seealso \code{\link{mbl}}
#'
#' @examples
#' # Fixed neighborhood sizes
#' neighbors_k(k = 50)
#' neighbors_k(k = c(40, 60, 80, 100, 120))
#'
#' # Dissimilarity threshold with default bounds
#' neighbors_diss(threshold = 0.3)
#'
#' # Dissimilarity threshold with custom bounds
#' neighbors_diss(threshold = c(0.1, 0.2, 0.3), k_min = 10, k_max = 150)
#'
#' @export
neighbors_k <- function(k) {
  
  if (missing(k) || is.null(k)) {
    stop("'k' is required.", call. = FALSE)
  }
  
  if (!is.numeric(k)) {
    stop("'k' must be a numeric vector.", call. = FALSE)
  }
  
  if (any(k != as.integer(k))) {
    stop("'k' must contain integer values.", call. = FALSE)
  }
  
  k <- unique(sort(as.integer(k)))
  
  if (any(k < 4L)) {
    stop("All values in 'k' must be at least 4.", call. = FALSE)
  }
  
  structure(
    list(
      method = "k",
      k = k
    ),
    class = c("neighbors_k", "neighbors")
  )
}


#' @export
neighbors_diss <- function(threshold, k_min = 4L, k_max = Inf) {
  
  if (missing(threshold) || is.null(threshold)) {
    stop("'threshold' is required.", call. = FALSE)
  }
  
  if (!is.numeric(threshold)) {
    stop("'threshold' must be a numeric vector.", call. = FALSE)
  }
  
  if (any(threshold <= 0)) {
    stop("All values in 'threshold' must be positive.", call. = FALSE)
  }
  
  threshold <- unique(sort(threshold))
  
  if (!is.numeric(k_min) || length(k_min) != 1L) {
    stop("'k_min' must be a single integer.", call. = FALSE)
  }
  
  k_min <- as.integer(k_min)
  
  if (k_min < 4L) {
    stop("'k_min' must be at least 4.", call. = FALSE)
  }
  
  if (!is.numeric(k_max) || length(k_max) != 1L) {
    stop("'k_max' must be a single integer or Inf.", call. = FALSE)
  }
  
  if (is.finite(k_max)) {
    k_max <- as.integer(k_max)
  }
  
  if (k_max <= k_min) {
    stop("'k_max' must be greater than 'k_min'.", call. = FALSE)
  }
  
  structure(
    list(
      method = "diss",
      threshold = threshold,
      k_min = k_min,
      k_max = k_max
    ),
    class = c("neighbors_diss", "neighbors")
  )
}


# -----------------------------------------------------------------------------
# Print methods
# -----------------------------------------------------------------------------

#' @export
print.neighbors_k <- function(x, ...) {
  k_str <- if (length(x$k) <= 6L) {
    paste(x$k, collapse = ", ")
  } else {
    paste0(
      paste(x$k[1:3], collapse = ", "),
      ", ..., ",
      paste(x$k[(length(x$k) - 1L):length(x$k)], collapse = ", ")
    )
  }
  cat(
    "Neighbor selection: fixed k\n",
    " k :", k_str, "\n"
  )
  invisible(x)
}


#' @export
print.neighbors_diss <- function(x, ...) {
  thr_str <- if (length(x$threshold) <= 6L) {
    paste(x$threshold, collapse = ", ")
  } else {
    paste0(
      paste(x$threshold[1:3], collapse = ", "),
      ", ..., ",
      paste(x$threshold[(length(x$threshold) - 1L):length(x$threshold)], collapse = ", ")
    )
  }
  k_max_str <- if (is.infinite(x$k_max)) "Inf" else x$k_max
  cat(
    "Neighbor selection: dissimilarity threshold\n",
    "  threshold :", thr_str, "\n",
    "  k_min     :", x$k_min, "\n",
    "  k_max     :", k_max_str, "\n"
  )
  invisible(x)
}


#' @export
print.neighbors <- function(x, ...) {
  cat("Neighbor selection method:", x$method, "\n")
  for (nm in setdiff(names(x), "method")) {
    cat(" ", nm, ":", x[[nm]], "\n")
  }
  invisible(x)
}