#' @title Correlation and moving correlation dissimilarity measurements  
#' @description
#' \loadmathjax
#' \ifelse{html}{\out{<a href='https://www.tidyverse.org/lifecycle/#satble'><img src='figures/lifecycle-stable.svg' alt='Stable lifecycle'></a>}}{\strong{Stable}}
#'
#' Computes correlation and moving correlation dissimilarity matrices.
#' @usage
#' cor_diss(Xr, Xu = NULL, ws = NULL,
#'          center = TRUE, scale = FALSE)
#' @param Xr a matrix.
#' @param Xu an optional matrix containing data of a second set of observations.
#' @param ws for moving correlation dissimilarity, an odd integer value which
#' specifies the window size. If \code{ws = NULL}, then the window size will be
#' equal to the number of variables (columns), i.e. instead moving correlation,
#' the normal correlation will be used. See details.
#' @param center a logical indicating if the spectral data \code{Xr} (and
#' \code{Xu} if specified) must be centered. If \code{Xu} is provided, the data
#' is scaled on the basis of \mjeqn{Xr \cup Xu}{Xr U Xu}.
#' @param scale a logical indicating if \code{Xr} (and \code{Xu} if specified)
#' must be scaled. If \code{Xu} is provided the data is scaled on the basis
#' of \mjeqn{Xr \cup Xu}{Xr U Xu}.
#' @param precision a character string indicating the numeric precision to use.
#' Possible values are \code{"double"} (default, 64-bit floating point) or
#' \code{"single"} (32-bit floating point). Using \code{"single"} reduces memory
#' usage and may improve performance for large datasets, at the cost of reduced
#' numerical precision.
#' @details
#' The correlation dissimilarity \mjeqn{d}{d} between two observations
#' \mjeqn{x_i}{x_i} and \mjeqn{x_j}{x_j} is based on the Perason's
#' correlation coefficient (\mjeqn{\rho}{\rho}) and it can be computed as
#' follows:
#'
#' \mjdeqn{d(x_i, x_j) = \frac{1}{2}((1 - \rho(x_i, x_j)))}{d(x_i, x_j) = 1/2 (1 - \rho(x_i, x_j))}
#'
#' The above formula is used when \code{ws = NULL}.
#' On the other hand (when \code{ws != NULL}) the moving correlation
#' dissimilarity between two observations \mjeqn{x_i}{x_i} and \mjeqn{x_j}{x_j}
#' is computed as follows:
#'
#' \mjdeqn{d(x_i, x_j; ws) = \frac{1}{2 ws}\sum_{k=1}^{p-ws}1 - \rho(x_{i,(k:k+ws)}, x_{j,(k:k+ws)})}{d(x_i, x_j) = 1/(2 ws)\sum_(k=1)^{p-ws}(1 - \rho(x_(i,k:k+ws), x_(j,k:k+ws)))}
#'
#' where \mjeqn{ws}{ws} represents a given window size which rolls sequentially
#' from 1 up to \mjeqn{p - ws}{p - ws} and  \mjeqn{p}{p} is the number of
#' variables of the observations.
#'
#' The function does not accept input data containing missing values.
#' @return
#' a matrix of the computed dissimilarities.
#' @author Antoine Stevens and \href{https://orcid.org/0000-0002-5369-5120}{Leonardo Ramirez-Lopez}
#' @examples
#' \donttest{
#' library(prospectr)
#' data(NIRsoil)
#'
#' Xu <- NIRsoil$spc[!as.logical(NIRsoil$train), ]
#' Xr <- NIRsoil$spc[as.logical(NIRsoil$train), ]
#'
#' cor_diss(Xr = Xr)
#'
#' cor_diss(Xr = Xr, Xu = Xu)
#'
#' cor_diss(Xr = Xr, ws = 41)
#'
#' cor_diss(Xr = Xr, Xu = Xu, ws = 41)
#' }
#' @export


######################################################################
# resemble
# Copyrigth (C) 2014 Leonardo Ramirez-Lopez and Antoine Stevens
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

## History:
## 09.03.2014 Leo     The line rslt[is.na(rslt)] <- 0 was added in order
##                    to deal with NaNs produced by the C++ code
## 21.04.2020 Leo     styler applied and Argument scaled renamed to scale
##                    the dimnames of the resulting matrix are now Xr_1... Xr_n
##                    (previusly Xr.1... Xr.n)
## 03.07.2020 Leo     FIXME: diss between the same observation in some values
##                    around 1e-15 are returned

cor_diss <- function(
    Xr, 
    Xu = NULL, 
    ws = NULL, 
    center = TRUE, 
    scale = FALSE, 
    precision = c("double", "single")
  ) {
  
  pr <- match.arg(precision, c("double", "single"))

  if (!ncol(Xr) >= 2) {
    stop("For correlation dissimilarity the number of variables must be larger than 1")
  }
  if (!is.null(Xu)) {
    if (ncol(Xu) != ncol(Xr)) {
      stop("The number of columns (variables) in Xr must be equal to the number of columns (variables) in Xu")
    }
    if (sum(is.na(Xu)) > 0) {
      stop("Input data contains missing values")
    }
  }

  if (sum(is.na(Xr)) > 0) {
    stop("Matrices with missing values are not accepted")
  }

  if (!is.logical(center)) {
    stop("'center' argument must be logical")
  }

  if (!is.logical(scale)) {
    stop("'scale' argument must be logical")
  }

  if (center | scale) {
    X <- rbind(Xr, Xu)

    if (center) {
      X <- sweep(x = X, MARGIN = 2, FUN = "-", STATS = colMeans(X))
    }

    if (scale) {
      X <- sweep(x = X, MARGIN = 2, FUN = "/", STATS = get_col_sds(X))
    }

    if (!is.null(Xu)) {
      Xu <- X[(nrow(X) - nrow(Xu) + 1):nrow(X), , drop = FALSE]
      Xr <- X[1:(nrow(X) - nrow(Xu)), ]
    } else {
      Xr <- X
    }
    rm(X)
  }
  
  xr_sds <- get_col_sds(t(Xr))
  if (any(xr_sds == 0)) {
    stop(paste0("Correlation coefficients cannot be computed. Xr contains ", sum(xr_sds == 0), " observation(s) with a standard deviation of zero."))
  }
  if (!is.null(Xu)) {
    xu_sds <- get_col_sds(t(Xu))
    if (any(xr_sds == 0)) {
      stop(paste0("Correlation coefficients cannot be computed. Xu contains ", sum(xu_sds == 0), "observation(s) with a standard deviation of zero."))
    }
  }
  
  if (is.null(ws)) {
    ws <- ncol(Xr)
  } else {
    
    if (ws < 3 | length(ws) != 1) {
      stop(paste("'ws' must be an odd value greater than 2"))
    }
    if ((ws %% 2) == 0) {
      stop("'ws' must be an odd value")
    }
    if (ws >= ncol(Xr)) {
      stop("'ws' must smaller than the number of columns (variables) in Xr")
    }
  }
  
  if (!is.null(Xu)) {
    pr <- match.arg(precision, c("double", "single"))
    rslt <- moving_cor_diss_xy(
      Xu, Xr, ws, 
      compute_block_rows(dim(Xu)),
      compute_block_rows(dim(Xr)), 
      precision = pr
    )
    colnames(rslt) <- paste("Xu", 1:nrow(Xu), sep = "_")
    rownames(rslt) <- paste("Xr", 1:nrow(Xr), sep = "_")
  } else {
    if (pr == "double") {
      rslt <- moving_cor_diss_self_f64(
        Xr, ws,
        compute_block_rows(dim(Xr))
      )
    } else  {
      rslt <- moving_cor_diss_self_f32(
        Xr, ws,
        compute_block_rows(dim(Xr))
      )
    }
    rownames(rslt) <- colnames(rslt) <- paste("Xr", 1:nrow(Xr), sep = "_")
  }

  rslt
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
#' The working set per tile (in elements of the chosen precision) is
#' approximated as:
#'
#' \deqn{M(b) \approx 3 b^2 + 2 b T} if \code{include_squares = TRUE}
#' (materialised \eqn{X_i^2}, \eqn{X_j^2}), and
#' \deqn{M(b) \approx 3 b^2} otherwise.
#'
#' A per-thread budget \eqn{K} (elements) is obtained from \code{target_mb},
#' the OS default if \code{NULL} (24 MiB on Windows, 48 MiB otherwise),
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
#' number of rows and columns of \eqn{X}.
#' @param include_squares Logical. If \code{TRUE}, assumes the kernel
#' materialises \code{Xi_sq} and \code{Xj_sq}; if \code{FALSE}, uses the
#' lighter memory model.
#' @param target_mb Numeric or \code{NULL}. Target MiB of cache to devote
#' per thread to one tile. Defaults to 24 on Windows and 48 otherwise.
#' @param min_tiles Integer. Aim for at least this many tiles along the
#' row dimension (default \code{12L}).
#' @param max_block_cap Integer. Hard upper bound for \code{block_rows}
#' (default \code{1024L}).
#'
#' @return Integer scalar, the recommended \code{block_rows} (multiple of 64).
#'
#' @section Cross-platform behaviour:
#' Uses \code{Sys.info()[["sysname"]]} when available to select a default
#' \code{target_mb}. Thread budgeting uses \code{OMP_NUM_THREADS} if set,
#' otherwise the number of physical cores; falls back to one when OpenMP
#' is not available.
#'
#' @examples
#' # Typical NIR-like shape
#' compute_block_rows(c(4096, 2150))
#'
#' # Float kernels or lighter memory model allow larger tiles
#' compute_block_rows(c(4096, 2150), precision = "single")
#' compute_block_rows(c(4096, 2150), include_squares = FALSE)
#'
#' # Override cache budget and threads explicitly
#' old <- Sys.getenv("OMP_NUM_THREADS")
#' Sys.setenv(OMP_NUM_THREADS = "16")
#' compute_block_rows(c(8192, 1800), target_mb = 32)
#' Sys.setenv(OMP_NUM_THREADS = old)
#'
#' # Small matrices respect caps
#' compute_block_rows(c(120, 2000))
#'
#' @seealso \code{\link{round_to_64}}
#' @keywords internal
#' @noRd
compute_block_rows <- function(dimensions,
                               include_squares = TRUE,
                               target_mb = 8,        # per-thread cache budget (MiB)
                               min_tiles = 12L,
                               max_block_cap = 1024L) {
  mm <- as.integer(dimensions[1])
  tt <- as.integer(dimensions[2])
  
  max_block <- min(max_block_cap, mm)
  if (max_block < 64L) return(max_block)
  
  # Effective threads for budgeting
  threads <- suppressWarnings(as.integer(Sys.getenv("OMP_NUM_THREADS", NA)))
  if (is.na(threads) || threads < 1L) {
    threads <- parallel::detectCores(logical = FALSE)
    if (is.na(threads) || threads < 1L) threads <- 1L
  }
  if (!isTRUE(capabilities("openmp"))) threads <- 1L
  
  # Per-thread element budget (double -> 8 bytes)
  elt_size <- 8
  per_thread_mb <- max(8, as.numeric(target_mb) / threads)  # keep a sane floor
  K <- (per_thread_mb * 1024^2) / elt_size                  # elements per tile
  
  # Working-set model ⇒ bound on b
  # include_squares:   3*b^2 + 2*b*T <= K  -> b_T = (-T + sqrt(T^2 + 3K))/3
  # without squares:   3*b^2          <= K  -> b_T = sqrt(K/3)
  b_T <- if (isTRUE(include_squares)) {
    (-tt + sqrt(tt * tt + 3 * K)) / 3
  } else {
    sqrt(K / 3)
  }
  
  # Also ensure enough tiles along rows
  b_m   <- mm / as.numeric(min_tiles)
  b_raw <- min(max(64, b_m), b_T)
  
  block <- round_to_64(b_raw, mult = 64L)
  block <- max(64L, min(as.integer(block), max_block))
  block
}



#' @keywords internal
#' @noRd
round_to_64 <- function(x, mult = 64L) {
  mult * round(x / mult)
}



