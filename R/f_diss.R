#' @title Euclidean, Mahalanobis and cosine dissimilarity measurements
#' @description
#' \loadmathjax
#' \ifelse{html}{\out{<a href='https://www.tidyverse.org/lifecycle/#satble'><img src='figures/lifecycle-stable.svg' alt='Stable lifecycle'></a>}}{\strong{Stable}}
#'
#' This function is used to compute the dissimilarity between observations
#' based on Euclidean or Mahalanobis distance measures or on cosine
#' dissimilarity measures (a.k.a spectral angle mapper).
#' @usage
#' f_diss(Xr, Xu = NULL, diss_method = "euclid",
#'        center = TRUE, scale = FALSE)
#' @param Xr a matrix containing the (reference) data.
#' @param Xu an optional matrix containing data of a second set of observations
#' (samples).
#' @param diss_method the method for computing the dissimilarity between
#' observations.
#' Options are \code{"euclid"} (Euclidean distance), \code{"mahalanobis"}
#' (Mahalanobis distance) and \code{"cosine"} (cosine distance, a.k.a spectral
#' angle mapper). See details.
#' @param center a logical indicating if the spectral data \code{Xr} (and
#' \code{Xu} if specified) must be centered. If \code{Xu} is provided, the data
#' is scaled on the basis of \mjeqn{Xr \cup Xu}{Xr U Xu}.
#' @param scale a logical indicating if \code{Xr} (and \code{Xu} if specified)
#' must be scaled. If \code{Xu} is provided the data is scaled on the basis
#' of \mjeqn{Xr \cup Xu}{Xr U Xu}.
#' @details
#' The results obtained for Euclidean dissimilarity are equivalent to those
#' returned by the [stats::dist()] function, but are scaled
#' differently. However, \code{f_diss} is considerably faster (which can be
#' advantageous when computing dissimilarities for very large matrices). The
#' final scaling of the dissimilarity scores in \code{f_diss} where
#' the number of variables is used to scale the squared dissimilarity scores. See
#' the examples section for a comparison between [stats::dist()] and
#' \code{f_diss}.
#'
#' In the case of both the Euclidean and Mahalanobis distances, the scaled
#' dissimilarity matrix \mjeqn{D}{D} between between observations in a given
#' matrix \mjeqn{X}{X} is computed as follows:
#'
#' \mjdeqn{d(x_i, x_j)^{2} = \sum (x_i - x_j)M^{-1}(x_i - x_j)^{\mathrm{T}}}{D(x_i, x_j)^{2} = \sum (x_i - x_j)M^{-1}(x_i - x_j)^T}
#' \mjdeqn{d_{scaled}(x_i, x_j) = \sqrt{\frac{1}{p}d(x_i, x_j)^{2}}}{d_scaled (x_i, x_j) = sqrt(1/p d(x_i, x_j)^2)}
#'
#' where \mjeqn{p}{p} is the number of variables in \mjeqn{X}{X}, \mjeqn{M}{M} is the identity
#' matrix in the case of the Euclidean distance and the variance-covariance
#' matrix of \mjeqn{X}{X} in the case of the Mahalanobis distance. The Mahalanobis
#' distance can also be viewed as the Euclidean distance after applying a
#' linear transformation of the original variables. Such a linear transformation
#' is done by using a factorization of the inverse covariance matrix as
#' \mjeqn{M^{-1} = W^{T}W}{M^-1 = W^TW}, where \mjeqn{M}{M} is merely the square root of
#' \mjeqn{M^{-1}}{M^{-1}} which can be found by using a singular value decomposition.
#'
#' Note that when attempting to compute the Mahalanobis distance on a data set
#' with highly correlated variables (i.e. spectral variables) the
#' variance-covariance matrix may result in a singular matrix which cannot be
#' inverted and therefore the distance cannot be computed.
#' This is also the case when the number of observations in the data set is
#' smaller than the number of variables.
#'
#' For the computation of the Mahalanobis distance, the mentioned method is
#' used.
#'
#' The cosine dissimilarity \mjeqn{c}{c} between two observations
#' \mjeqn{x_i}{x_i} and \mjeqn{x_j}{x_j} is computed as follows:
#'
#' \mjdeqn{c(x_i, x_j) = cos^{-1}{\frac{\sum_{k=1}^{p}x_{i,k} x_{j,k}}{\sqrt{\sum_{k=1}^{p} x_{i,k}^{2}} \sqrt{\sum_{k=1}^{p} x_{j,k}^{2}}}}}{c(x_i, x_j) = cos^{-1} ((sum_(k=1)^p x_(i,k) x_(j,k))/(sum_(k=1)^p x_(i,k) sum_(k=1)^p x_(j,k)))}
#'
#' where \mjeqn{p}{p} is the number of variables of the observations.
#' The function does not accept input data containing missing values.
#' NOTE: The computed distances are divided by the number of variables/columns
#' in \code{Xr}.
#'
#' @return
#' a matrix of the computed dissimilarities.
#' @author \href{https://orcid.org/0000-0002-5369-5120}{Leonardo Ramirez-Lopez} and Antoine Stevens
#' @examples
#' \donttest{
#' library(prospectr)
#' data(NIRsoil)
#'
#' Xu <- NIRsoil$spc[!as.logical(NIRsoil$train), ]
#' Xr <- NIRsoil$spc[as.logical(NIRsoil$train), ]
#'
#' # Euclidean distances between all the observations in Xr
#'
#' ed <- f_diss(Xr = Xr, diss_method = "euclid")
#'
#' # Equivalence with the dist() fucntion of R base
#' ed_dist <- (as.matrix(dist(Xr))^2 / ncol(Xr))^0.5
#' round(ed_dist - ed, 5)
#'
#' # Comparing the computational time
#' iter <- 20
#' tm <- proc.time()
#' for (i in 1:iter) {
#'   f_diss(Xr)
#' }
#' f_diss_time <- proc.time() - tm
#'
#' tm_2 <- proc.time()
#' for (i in 1:iter) {
#'   dist(Xr)
#' }
#' dist_time <- proc.time() - tm_2
#'
#' f_diss_time
#' dist_time
#'
#' # Euclidean distances between observations in Xr and observations in Xu
#' ed_xr_xu <- f_diss(Xr, Xu)
#'
#' # Mahalanobis distance computed on the first 20 spectral variables
#' md_xr_xu <- f_diss(Xr[, 1:20], Xu[, 1:20], "mahalanobis")
#'
#' # Cosine dissimilarity matrix
#' cdiss_xr_xu <- f_diss(Xr, Xu, "cosine")
#' }
#' @export

######################################################################
# resemble
# Copyright (C) 2014 Leonardo Ramirez-Lopez and Antoine Stevens
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
## 18.11.2015 Leo     Bug fixed. Code crashed When Xr was of one row
## 18.07.2019 Leo     A bug in the scaling of the euclidean distances in f_diss was detected and fixed. The distance ratios
##                    (between observations) were correctly calculated, but the final scaling of the results was not properly
##                    done. The distance between Xi and Xj were scaled by taking the squared root of the mean of the squared
##                    differences and dividing it by the number of variables i.e. sqrt(mean((Xi-Xj)^2))/ncol(Xi), however the
##                    correct calculation is done by taking the mean of the squared differences, dividing it by the number of
##                    variables and then compute the squared root i.e. sqrt(mean((Xi-Xj)^2)/ncol(Xi)).
##                    This bug had no effect on the computations of the nearest neighbors.
## 21.04.2020 Leo     styler applied and Argument scaled renamed to scale
##                    the dimnames of the resulting matrix are now Xr_1... Xr_n (previusly Xr.1... Xr.n)
## 03.07.2020 Leo     FIXME: with cosine  diss between the same observation in some
##                    cases NaN are returned and also values around 1e-08

f_diss <- function(Xr, Xu = NULL, diss_method = "euclid",
                   center = TRUE, scale = FALSE) {
  if (!is.null(Xu)) {
    if (ncol(Xu) != ncol(Xr)) {
      stop("The number of columns (variables) in Xr must be equal to \n the number of columns (variables) in Xu")
    }
    if (sum(is.na(Xu)) > 0) {
      stop("Input data contains missing values")
    }
  }
  if (sum(is.na(Xr)) > 0) {
    stop("Matrices with missing values are not accepted")
  }

  n_method <- diss_method
  if (!n_method %in% c("euclid", "mahalanobis", "cosine")) {
    stop("'diss_method' must be one of: 'euclid', 'mahalanobis' or'cosine'")

    if (length(n_method) > 1) {
    } else {
      n_method <- diss_method[[1]]
      message(paste("More than one diss_method was specified, only", n_method, "was used."))
    }
  }

  if (!is.logical(center)) {
    stop("'center' must be logical")
  }

  if (!is.logical(scale)) {
    stop("'scale' must be logical")
  }

  if (center | scale | n_method %in% c("mahalanobis", "euclid")) {
    X <- rbind(Xr, Xu)

    if (center) {
      X <- sweep(x = X, MARGIN = 2, FUN = "-", STATS = colMeans(X))
    }

    if (scale) {
      X <- sweep(x = X, MARGIN = 2, FUN = "/", STATS = get_col_sds(X))
    }

    if (n_method == "mahalanobis") {
      if (nrow(X) < ncol(X)) {
        stop("For computing the Mahalanobis distance, the total number of observations (rows) \n must be larger than the number of variables (columns).")
      }
      X <- try(euclid_to_mahal(X, sm_method = "svd"), TRUE)
      if (!is.matrix(X)) {
        stop("The covariance matrix (for the computation of the Mahalanobis distance) is exactly singular. \n Try another method.")
      }
      n_method <- "euclid"
    }

    if (!is.null(Xu)) {
      Xu <- X[(nrow(X) - nrow(Xu) + 1):nrow(X), , drop = FALSE]
      Xr <- X[1:(nrow(X) - nrow(Xu)), , drop = FALSE]
    } else {
      Xr <- X
    }
    rm(X)
  }

  if (!is.null(Xu)) {
    ## FIXME check numerical precision in Rcpp
    ## in some cases it returns 0s as -1e-14
    ## perhaps due to reuse memory?
    rslt <- abs(fast_diss(Xu, Xr, n_method))
    if (n_method == "euclid") {
      rslt <- sqrt(rslt / ncol(Xr))
    }
    colnames(rslt) <- paste("Xu", 1:nrow(Xu), sep = "_")
    rownames(rslt) <- paste("Xr", 1:nrow(Xr), sep = "_")
  } else {
    ## FIXME check numerical precision in Rcpp
    ## in some cases it returns 0s as -1e-14
    ## perhaps due to reuse memory?
    rslt <- abs(fast_diss(Xr, Xr, n_method))
    if (n_method == "euclid") {
      rslt <- sqrt(rslt / ncol(Xr))
    }
    rownames(rslt) <- paste("Xr", 1:nrow(Xr), sep = "_")
    colnames(rslt) <- rownames(rslt)
  }
  if (diss_method == "cosine") {
    rslt[is.nan(rslt)] <- 0
  }

  rslt
}


#' @title A function for transforming a matrix from its Euclidean space to
#' its Mahalanobis space
#' @description For internal use only
#' @keywords internal
#' @importFrom stats cov
euclid_to_mahal <- function(X, sm_method = c("svd", "eigen")) {
  nms <- dimnames(X)

  if (ncol(X) > nrow(X)) {
    stop("In order to project the matrix to a Mahalanobis space, the number of observations of the input matrix must larger than its number of variables")
  }

  if (length(sm_method) > 1) {
    sm_method <- sm_method[1]
  }
  if (!(sm_method %in% c("svd", "eigen"))) {
    stop("sm_method must be one of 'svd', 'eigen'")
  }

  X <- as.matrix(X)
  vcv <- cov(X)
  sq_vcv <- sqrt_sm(vcv, method = sm_method)
  sq_S <- solve(sq_vcv)
  ms_x <- X %*% sq_S
  dimnames(ms_x) <- nms

  ms_x
}


#' @title Square root of (square) symmetric matrices
#' @description For internal use only
#' @keywords internal
sqrt_sm <- function(X, method = c("svd", "eigen")) {
  if (!isSymmetric(X)) {
    stop("X must be a square symmetric matrix")
  }
  if (length(method) > 1) {
    method <- method[1]
  }
  if (!(method %in% c("svd", "eigen"))) {
    stop("method must be one of 'svd', 'eigen'")
  }

  if (method == "svd") {
    ## REPLACE BY  arma::svd(U, S, V, X, "dc")
    out <- svd(X)
    D <- diag(out$d)
    U <- out$v
    return(U %*% (D^0.5) %*% t(U))
  }

  if (method == "eigen") {
    out <- eigen(X)
    D <- diag(out$values)
    U <- out$vectors
    return(U %*% (D^0.5) %*% t(U))
  }
}
