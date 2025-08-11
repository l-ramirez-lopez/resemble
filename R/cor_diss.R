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

  rslt[rslt < 1e-15] <- 0
  rslt
}


# R/block_rows.R

#' Choose a cache-friendly tile size (rows) from m and T (no options/env)
#'
#' Heuristic:
#' - Working-set budget K ≈ 24 MB of doubles → K = 3e6.
#' - Solve 3*b^2 + 2*T*b ≈ K for b: b_T = (-T + sqrt(T^2 + 3K)) / 3.
#' - Limit by m via ~12 tiles/axis: b_m = m / 12.
#' - Pick b* = round_to_64( min( max(64, b_m), b_T ) ), then clamp to [64, min(1024, m)].
#'
#' @param dimensions a vector of two integers representing the number of rows 
#' (samples) and the number of columns (features)
#' @return integer block_rows (multiple of 64), clamped; if m < 64, returns m
#' @examples
#' compute_block_rows(c(2475, 700))   # ~256
#' compute_block_rows(c(80000, 700))  # ~768..1024 (clamped to 1024 if needed)
#' compute_block_rows(c(50, 700))     # 50 (since m < 64)
#' @keywords internal
compute_block_rows <- function(dimensions) {
  mm <- as.integer(dimensions[1])
  tt <- as.integer(dimensions[2])
  
  max_block <- min(1024L, mm)
  if (max_block < 64L) return(max_block)
  
  K    <- 3e6                             # ~24 MB in doubles
  b_m  <- mm / 12
  b_T  <- (-tt + sqrt(tt*tt + 3*K)) / 3
  b_raw <- min(max(64, b_m), b_T)
  
  block <- round_to_64(b_raw)
  block <- max(64L, min(as.integer(block), max_block))
  block
}

#' @keywords internal
#' @noRd
round_to_64 <- function(x) {
  mult <- 64
  mult * round(x / mult)
}

