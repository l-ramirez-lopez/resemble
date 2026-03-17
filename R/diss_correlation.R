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
#' @param precision Character string, either \code{"double"} (default, 64-bit)
#'   or \code{"single"} (32-bit). Using \code{"single"} reduces memory usage
#'   and may improve performance on large datasets at the cost of numerical
#'   precision.
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
    ws        = NULL,
    center    = TRUE,
    scale     = FALSE,
    precision = c("double", "single")
) {
  
  # --- validate precision ----------------------------------------------------
  # done here before passing to .new_diss_method() since it is
  # correlation-specific (the other methods do not have this argument)
  precision <- match.arg(precision)
  
  # --- validate ws -----------------------------------------------------------
  # Note: upper bound check (ws < ncol(Xr)) is data-dependent and is
  # therefore deferred to dissimilarity()
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
  # extra fields (ws, precision) are correlation-specific and are passed
  # through to the structure via the extra argument
  .new_diss_method(
    "correlation",
    center = center,
    scale  = scale,
    extra  = list(ws = ws, precision = precision)
  )
}


# --- print method ------------------------------------------------------------
# diss_correlation has extra fields (ws, precision) that the shared
# print.diss_method does not know about, so it gets its own print method.
# This takes precedence over print.diss_method due to S3 dispatch order.
#' @export
print.diss_correlation <- function(x, ...) {
  ws_label <- if (is.null(x$ws)) "full (no moving window)" else x$ws
  cat(
    "Dissimilarity method: correlation\n",
    "  window size (ws) :", ws_label, "\n",
    "  center           :", x$center, "\n",
    "  scale            :", x$scale, "\n",
    "  precision        :", x$precision, "\n"
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
  
  ws        <- method$ws
  center    <- method$center
  scale     <- method$scale
  precision <- method$precision
  
  # --- preprocessing ---------------------------------------------------------
  if (center || scale) {
    X <- rbind(Xr, Xu)
    
    if (center) {
      X <- sweep(X, MARGIN = 2, STATS = colMeans(X), FUN = "-")
    }
    if (scale) {
      X <- sweep(X, MARGIN = 2, STATS = get_col_sds(X), FUN = "/")
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
  
  # --- zero-sd guard ---------------------------------------------------------
  # Checked after preprocessing because centering can surface zero-variance
  # observations that were non-constant before centering.
  xr_sds <- get_col_sds(t(Xr))
  if (any(xr_sds == 0)) {
    stop(
      sprintf(
        "Correlation cannot be computed: Xr contains %d observation(s) with a standard deviation of zero.",
        sum(xr_sds == 0)
      )
    )
  }
  
  if (!is.null(Xu)) {
    xu_sds <- get_col_sds(t(Xu))
    if (any(xu_sds == 0)) {
      stop(
        sprintf(
          "Correlation cannot be computed: Xu contains %d observation(s) with a standard deviation of zero.",
          sum(xu_sds == 0)
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
      compute_block_rows(dim(Xr)),
      precision = precision
    )
    rownames(result) <- paste("Xr", seq_len(nrow(Xr)), sep = "_")
    colnames(result) <- paste("Xu", seq_len(nrow(Xu)), sep = "_")
  } else {
    result <- if (precision == "double") {
      moving_cor_diss_self_f64(Xr, ws, compute_block_rows(dim(Xr)))
    } else {
      moving_cor_diss_self_f32(Xr, ws, compute_block_rows(dim(Xr)))
    }
    rownames(result) <- colnames(result) <- paste("Xr", seq_len(nrow(Xr)), sep = "_")
  }
  
  result
}