# =============================================================================
# Correlation dissimilarity — aligned with shared constructor pattern
# =============================================================================
# Public:   diss_correlation()   — method constructor
# Internal: .cor_diss_compute()  — raw computation, not exported


#' @title Correlation dissimilarity method constructor
#'
#' @description
#' Creates a configuration object that fully specifies a correlation (or moving
#' correlation) dissimilarity method. Pass the result to \code{dissimilarity()}
#' to compute the dissimilarity matrix.
#'
#' @param ws Either \code{NULL} (default) or an odd integer greater than 2.
#'   When \code{NULL}, standard Pearson correlation dissimilarity is used.
#'   When an odd integer is provided, a moving (rolling) correlation
#'   dissimilarity is computed using a window of that size.
#' @param center Logical. Should the data be mean-centered before computing
#'   dissimilarities? Centering is applied jointly to \code{Xr} and \code{Xu}
#'   (if provided) based on their combined column means. Default is \code{TRUE}.
#' @param scale Logical. Should the data be scaled (divided by column standard
#'   deviations) before computing dissimilarities? Scaling is applied jointly
#'   to \code{Xr} and \code{Xu} (if provided). Default is \code{FALSE}.
#'
#' @return An object of class \code{c("diss_correlation", "diss_method")} — a
#'   list holding the validated method parameters. Intended to be passed to
#'   \code{dissimilarity()}, not used directly.
#'
#' @details
#' The correlation dissimilarity between two observations \eqn{x_i} and
#' \eqn{x_j} is:
#'
#' \deqn{d(x_i, x_j) = \frac{1}{2}(1 - \rho(x_i, x_j))}
#'
#' where \eqn{\rho} is the Pearson correlation coefficient. This is used when
#' \code{ws = NULL}.
#'
#' When \code{ws} is specified, the moving correlation dissimilarity is:
#'
#' \deqn{d(x_i, x_j; ws) = \frac{1}{2\,ws} \sum_{k=1}^{p - ws}
#'   \bigl(1 - \rho(x_{i,(k:k+ws)},\, x_{j,(k:k+ws)})\bigr)}
#'
#' where \eqn{ws} is the window size and \eqn{p} is the number of variables.
#'
#' @section Parallel execution:
#' The underlying C++ implementation uses OpenMP for parallel computation.
#' Thread count is controlled by the \code{OMP_NUM_THREADS} environment
#' variable. To limit threads (e.g., when calling from within a parallel
#' backend):
#' \preformatted{
#' Sys.setenv(OMP_NUM_THREADS = 1)
#' }
#'
#' @author \href{https://orcid.org/0000-0002-5369-5120}{Leonardo Ramirez-Lopez}
#'
#' @seealso \code{\link{dissimilarity}}, \code{\link{diss_euclidean}},
#'   \code{\link{diss_mahalanobis}}, \code{\link{diss_cosine}}
#'
#' @examples
#' # Standard correlation dissimilarity
#' m <- diss_correlation()
#'
#' # Moving correlation with window size 41
#' m <- diss_correlation(ws = 41)
#'
#' # Without centering
#' m <- diss_correlation(center = FALSE)
#'
#' @export
diss_correlation <- function(
    ws     = NULL,
    center = TRUE,
    scale  = FALSE
) {
  
  # --- validate ws -----------------------------------------------------------
  if (!is.null(ws)) {
    if (!is.numeric(ws) || length(ws) != 1L) {
      stop("'ws' must be a single integer value or NULL.")
    }
    ws <- as.integer(ws)
    if (ws < 3L) {
      stop("'ws' must be an odd integer greater than 2.")
    }
    if (ws %% 2L == 0L) {
      stop("'ws' must be an odd integer (e.g. 3, 5, 11, 41, ...).")
    }
  }
  
  # --- delegate center/scale validation and object construction --------------
  .new_diss_method(
    "correlation",
    center = center,
    scale  = scale,
    extra  = list(ws = ws)
  )
}


# --- print method ------------------------------------------------------------
#' @export
print.diss_correlation <- function(x, ...) {
  ws_label <- if (is.null(x$ws)) "full (no moving window)" else x$ws
  cat(
    "Dissimilarity method: correlation\n",
    "  window size (ws) :", ws_label, "\n",
    "  center           :", x$center, "\n",
    "  scale            :", x$scale, "\n"
  )
  invisible(x)
}


# =============================================================================
# Internal compute function — NOT exported
# =============================================================================
#
# Called by dissimilarity() after data-level validation has been done there.
# Assumes:
#   - Xr and Xu are numeric matrices with no NAs
#   - ncol(Xu) == ncol(Xr) if Xu is provided
#   - method$ws < ncol(Xr) if ws is not NULL (checked in dissimilarity())
#
.cor_diss_compute <- function(Xr, Xu = NULL, method) {
  
  stopifnot(inherits(method, "diss_correlation"))
  
  ws     <- method$ws
  center <- method$center
  scale  <- method$scale
  
  # --- preprocessing ---------------------------------------------------------
  if (center || scale) {
    X <- rbind(Xr, Xu)
    
    if (center) {
      X <- sweep(X, MARGIN = 2, STATS = colMeans(X), FUN = "-")
    }
    if (scale) {
      col_sds <- get_col_sds(X)
      if (any(col_sds == 0)) {
        stop("Cannot scale: one or more columns have zero standard deviation.")
      }
      X <- sweep(X, MARGIN = 2, STATS = col_sds, FUN = "/")
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
  
  # --- zero-sd guard (row-level for correlation) -----------------------------
  xr_row_sds <- get_col_sds(t(Xr))
  if (any(xr_row_sds == 0)) {
    stop(
      sprintf(
        "Correlation cannot be computed: Xr contains %d observation(s) with zero variance.",
        sum(xr_row_sds == 0)
      )
    )
  }
  
  if (!is.null(Xu)) {
    xu_row_sds <- get_col_sds(t(Xu))
    if (any(xu_row_sds == 0)) {
      stop(
        sprintf(
          "Correlation cannot be computed: Xu contains %d observation(s) with zero variance.",
          sum(xu_row_sds == 0)
        )
      )
    }
  }
  
  # --- resolve full-correlation window size ----------------------------------
  if (is.null(ws)) {
    ws <- ncol(Xr)
  }
  
  # --- dispatch to C++ backend -----------------------------------------------
  if (!is.null(Xu)) {
    result <- moving_cor_diss_xy(
      Xu, Xr, ws,
      compute_block_rows(dim(Xu)),
      compute_block_rows(dim(Xr))
    )
    rownames(result) <- paste("Xr", seq_len(nrow(Xr)), sep = "_")
    colnames(result) <- paste("Xu", seq_len(nrow(Xu)), sep = "_")
  } else {
    result <- moving_cor_diss_self_f64(Xr, ws, compute_block_rows(dim(Xr)))
    rownames(result) <- colnames(result) <- paste("Xr", seq_len(nrow(Xr)), sep = "_")
  }
  
  result
}


#' Heuristic for cache-aware tile height (block_rows)
#'
#' @description
#' Choose a tile height for tiled rolling-window correlation/distance kernels,
#' balancing per-tile cache footprint against matrix width and available
#' threads. Cross-platform (macOS, Linux, Windows).
#'
#' @details
#' Let \eqn{b} be the tile height and \eqn{T} the number of columns.
#' The working set per tile (in elements) is approximated as:
#'
#' \deqn{M(b) \approx 3 b^2 + 2 b T} if \code{include_squares = TRUE}
#' (materialised \eqn{X_i^2}, \eqn{X_j^2}), and
#' \deqn{M(b) \approx 3 b^2} otherwise.
#'
#' A per-thread budget \eqn{K} (elements) is obtained from \code{target_mb},
#' divided by the effective OpenMP thread count. The quadratic bound yields
#' \deqn{b_T = \frac{-T + \sqrt{T^2 + 3K}}{3}}
#' (or \eqn{b_T = \sqrt{K/3}} when \code{include_squares = FALSE}).
#'
#' The final choice is
#' \code{block_rows = round_to_64( min( max(64, m / min_tiles), b_T ) )},
#' clamped to \code{[64, min(max_block_cap, m)]}.
#'
#' If OpenMP is unavailable (\code{capabilities("openmp") == FALSE}),
#' the budget assumes a single thread.
#'
#' @param dimensions Integer vector of length two \code{c(m, T)}:
#'   number of rows and columns of \eqn{X}.
#' @param include_squares Logical. If \code{TRUE}, assumes the kernel
#'   materialises \code{Xi_sq} and \code{Xj_sq}; if \code{FALSE}, uses the
#'   lighter memory model.
#' @param target_mb Numeric. Target MiB of cache to devote per thread to one
#'   tile. Default is 8.
#' @param min_tiles Integer. Aim for at least this many tiles along the
#'   row dimension (default \code{12L}).
#' @param max_block_cap Integer. Hard upper bound for \code{block_rows}
#'   (default \code{1024L}).
#'
#' @return Integer scalar, the recommended \code{block_rows} (multiple of 64).
#'
#' @section Thread detection:
#' Thread count is determined by:
#' \enumerate{
#'   \item \code{OMP_NUM_THREADS} environment variable (if set)
#'   \item \code{parallel::detectCores(logical = FALSE)} (physical cores)
#'   \item Falls back to 1 if detection fails or OpenMP is unavailable
#' }
#'
#' @keywords internal
#' @noRd
compute_block_rows <- function(
    dimensions,
    include_squares = TRUE,
    target_mb       = 8,
    min_tiles       = 12L,
    max_block_cap   = 1024L
) {
  mm <- as.integer(dimensions[1])
  tt <- as.integer(dimensions[2])
  
  max_block <- min(max_block_cap, mm)
  if (max_block < 64L) return(max_block)
  
  # --- detect effective thread count -----------------------------------------
  threads <- suppressWarnings(as.integer(Sys.getenv("OMP_NUM_THREADS", NA)))
  if (is.na(threads) || threads < 1L) {
    threads <- parallel::detectCores(logical = FALSE)
    if (is.na(threads) || threads < 1L) threads <- 1L
  }
  if (!isTRUE(capabilities("openmp"))) threads <- 1L
  
  # --- per-thread element budget (double -> 8 bytes) -------------------------
  elt_size <- 8
  per_thread_mb <- max(8, as.numeric(target_mb) / threads)
  K <- (per_thread_mb * 1024^2) / elt_size
  
  # --- working-set model bound -----------------------------------------------
  b_T <- if (isTRUE(include_squares)) {
    (-tt + sqrt(tt * tt + 3 * K)) / 3
  } else {
    sqrt(K / 3)
  }
  
  # --- ensure enough tiles along rows ----------------------------------------
  b_m   <- mm / as.numeric(min_tiles)
  b_raw <- min(max(64, b_m), b_T)
  
  block <- round_to_64(b_raw)
  block <- max(64L, min(as.integer(block), max_block))
  block
}


#' Round to nearest multiple of 64
#' @param x Numeric value to round.
#' @param mult Integer multiple (default 64).
#' @return Integer rounded to nearest multiple.
#' @keywords internal
#' @noRd
round_to_64 <- function(x, mult = 64L) {
  as.integer(mult * round(x / mult))
}
