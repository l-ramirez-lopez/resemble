#' @title A function for evaluating dissimilarity matrices (sim_eval)
#' @description
#' \loadmathjax
#'
#' \ifelse{html}{\out{<a href='https://www.tidyverse.org/lifecycle/#satble'><img src='figures/lifecycle-stable.svg' alt='Stable lifecycle'></a>}}{\strong{Stable}}
#'
#' This function searches for the most similar observation (closest neighbor) of
#' each observation in a given data set based on a dissimilarity (e.g. distance
#' matrix). The observations are compared against their corresponding closest
#' observations in terms of their side information provided. The root mean
#' square of differences and the correlation coefficient are used for continuous
#' variables and for discrete variables the kappa index is used.
#' @usage
#' sim_eval(d, side_info)
#' @param d a symmetric matrix of dissimilarity scores between observations of
#' a given data set. Alternatively, a vector of with the dissimilarity
#' scores of the lower triangle (without the diagonal values) can be used
#' (see details).
#' @param side_info a matrix containing the side information corresponding to
#' the observations in the data set from which the dissimilarity matrix was
#' computed. It can be either a numeric matrix with one or multiple
#' columns/variables or a matrix with one character variable (discrete variable).
#' If it is numeric, the root mean square of differences is used for assessing
#' the similarity between the observations and their corresponding most similar
#' observations in terms of the side information provided. If it is a character
#' variable, then the kappa index is used. See details.
#' @details
#' For the evaluation of dissimilarity matrices this function uses side
#' information (information about one variable which is available for a
#' group of observations, Ramirez-Lopez et al., 2013). It is assumed that there
#' is a (direct or indirect) correlation between this side informative variable
#' and the variables from which the dissimilarity was computed.
#' If \code{side_info} is numeric, the root mean square of differences (RMSD)
#' is used for assessing the similarity between the observations and their
#' corresponding most similar observations in terms of the side information
#' provided. It is computed as follows:
#'
#' \mjdeqn{j(i) = NN(xr_i, Xr^{\{-i\}})}{j(i) = NN(xr_i, Xr^{\{-i\}})}
#' \mjdeqn{RMSD = \sqrt{\frac{1}{m} \sum_{i=1}^n {(y_i - y_{j(i)})^2}}}{RMSD = \sqrt{1/n sum_{i=1}^m (y_i - y_{j(i)})^2}}
#'
#' where \mjeqn{NN(xr_i, Xr^{-i})}{NN(xr_i, Xr^{-i})} represents a function to
#' obtain the index  of the nearest neighbor observation found in \mjeqn{Xr}{Xr}
#' (excluding the \mjeqn{i}{i}th observation) for \mjeqn{xr_i}{xr_i},
#' \mjeqn{y_{i}}{y_i} is the value of the side variable of the \mjeqn{i}{i}th
#' observation, \mjeqn{y_{j(i)}}{y_{j(i)}} is the value of the side variable of
#' the nearest neighbor of the \mjeqn{i}{i}th observation and \mjeqn{m}{m} is
#' the total number of observations.
#'
#' If \code{side_info} is a factor the kappa index (\mjeqn{\kappa}{kappa}) is
#' used instead the RMSD. It is computed as follows:
#'
#' \mjdeqn{\kappa = \frac{p_{o}-p_{e}}{1-p_{e}}}{kappa = {p_o-p_e}/{1-p_e}}
#'
#' where both \mjeqn{p_o}{p_o} and \mjeqn{p_e}{p_e} are two different agreement
#' indices between the the side information of the observations and the side
#' information of their corresponding nearest observations (i.e. most similar
#' observations). While \mjeqn{p_o}{p_o} is the relative agreement
#' \mjeqn{p_e}{p_e} is the the agreement expected by chance.
#'
#' This functions accepts vectors to be passed to argument \code{d}, in this
#' case, the vector must represent the lower triangle of a dissimilarity matrix
#' (e.g. as returned by the [stats::dist()] function of \code{stats}).
#'
#' This function supports multi-threading based on OpenMP for retrieving the
#' closest observations.
#' @return \code{sim_eval} returns a list with the following components:
#' \itemize{
#'  \item{"\code{eval}}{ either the RMSD (and the correlation coefficient) or
#'  the kappa index}
#'  \item{\code{first_nn}}{ a matrix containing the original side
#'  informative variable in the first half of the columns, and the side
#'  informative values of the corresponding nearest neighbors in the second half
#'  of the columns.}
#'  }
#' @author \href{https://orcid.org/0000-0002-5369-5120}{Leonardo Ramirez-Lopez}
#' @references
#' Ramirez-Lopez, L., Behrens, T., Schmidt, K., Stevens, A., Dematte, J.A.M.,
#' Scholten, T. 2013a. The spectrum-based learner: A new local approach for
#' modeling soil vis-NIR spectra of complex datasets. Geoderma 195-196, 268-279.
#'
#' Ramirez-Lopez, L., Behrens, T., Schmidt, K., Viscarra Rossel, R.,
#' Dematte, J. A. M.,  Scholten, T. 2013b. Distance and similarity-search
#' metrics for use with soil vis-NIR spectra. Geoderma 199, 43-53.
#' @examples
#' \dontrun{
#' library(prospectr)
#' data(NIRsoil)
#'
#' sg <- savitzkyGolay(NIRsoil$spc, p = 3, w = 11, m = 0)
#'
#' # Replace the original spectra with the filtered ones
#' NIRsoil$spc <- sg
#'
#' Yr <- NIRsoil$Nt[as.logical(NIRsoil$train)]
#' Xr <- NIRsoil$spc[as.logical(NIRsoil$train), ]
#'
#' # Example 1
#' # Compute a principal components distance
#' pca_d <- ortho_diss(Xr, pc_selection = list("manual", 8))$dissimilarity
#'
#' # Example 1.1
#' # Evaluate the distance matrix on the baisis of the
#' # side information (Yr) associated with Xr
#' se <- sim_eval(pca_d, side_info = as.matrix(Yr))
#'
#' # The final evaluation results
#' se$eval
#'
#' # The final values of the side information (Yr) and the values of
#' # the side information corresponding to the first nearest neighbors
#' # found by using the distance matrix
#' se$first_nn
#'
#' # Example 1.2
#' # Evaluate the distance matrix on the basis of two side
#' # information (Yr and Yr2)
#' # variables associated with Xr
#' Yr_2 <- NIRsoil$CEC[as.logical(NIRsoil$train)]
#' se_2 <- sim_eval(d = pca_d, side_info = cbind(Yr, Yr_2))
#'
#' # The final evaluation results
#' se_2$eval
#'
#' # The final values of the side information variables and the values
#' # of the side information variables corresponding to the first
#' # nearest neighbors found by using the distance matrix
#' se_2$first_nn
#'
#' # Example 2
#' # Evaluate the distances produced by retaining different number of
#' # principal components (this is the same principle used in the
#' # optimized principal components approach ("opc"))
#'
#' # first project the data
#' pca_2 <- ortho_projection(Xr, pc_selection = list("manual", 30))
#'
#' results <- matrix(NA, pca_2$n_components, 3)
#' colnames(results) <- c("pcs", "rmsd", "r")
#' results[, 1] <- 1:pca_2$n_components
#' for (i in 1:pca_2$n_components) {
#'   ith_d <- f_diss(pca_2$scores[, 1:i, drop = FALSE], scale = TRUE)
#'   ith_eval <- sim_eval(ith_d, side_info = as.matrix(Yr))
#'   results[i, 2:3] <- as.vector(ith_eval$eval)
#' }
#' plot(results)
#'
#' # Example 3
#' # Example 3.1
#' # Evaluate a dissimilarity matrix computed using the correlation
#' # method
#' cd <- cor_diss(Xr)
#' eval_corr_diss <- sim_eval(cd, side_info = as.matrix(Yr))
#' eval_corr_diss$eval
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
## 09.03.2014 Leo     In the doc was specified that multi-threading is
##                    not working for mac
## 13.03.2014 Antoine The explanation of the cores argument was modified
## 18.03.2014 Antoine Add error message when input dissimilarity matrix
##                    in sim_eval is not squared
## 26.04.2020 Leo     Function was slightly refactored for tidiness
##                    cores is deprecated

sim_eval <- function(d, side_info) {
  if (!is.matrix(side_info)) {
    stop("side_inf must be a matrix")
  }
  if (is.character(side_info)) {
    if (ncol(side_info) > 1) {
      stop("For character variables in 'side_info' only one variable is allowed")
    }
    side_info <- as.factor(side_info)
    if (nlevels(side_info) < 2) {
      stop("Only one category detected in 'side_info'. At least two categories are required")
    }
    get_eval <- function(y, indices_closest) {
      get_eval_categorical(y, indices_closest)
    }
  } else {
    get_eval <- function(y, indices_closest) {
      get_eval_continuous(y, indices_closest)
    }
  }
  if (sum(colSums(!is.na(side_info)) < 4) > 0) {
    stop("At least one of the side information variables contains less than 4 observations")
  }
  if (!(is.vector(d) | is.matrix(d))) {
    stop("d must be either a matrix of distaor ")
  }
  if (is.vector(d)) {
    if (length(d) != (nrow(side_info)^2 - nrow(side_info)) / 2) {
      stop("The length of 'd' do not match with the length of a triangular matrix with the same number of observations sa in 'side_info'")
    }
    use_function <- function(X) {
      which_min_vector(X)
    }
  }
  if (is.matrix(d)) {
    if (diff(dim(d)) != 0) {
      stop("matrices passed to 'd' they must be squared")
    }
    if (nrow(d) != nrow(side_info)) {
      stop("The number of rows of the 'd' matrix does not match the number of observations in 'side_info'")
    }
    use_function <- function(X) {
      which_min(X)
    }
  }

  closest_found <- use_function(d)

  cnms <- colnames(side_info)
  if (is.null(cnms)) {
    cnms <- paste0("side_info_", 1:ncol(side_info))
  }

  rslt <- get_eval(y = side_info, indices_closest = closest_found)
  rownames(rslt) <- cnms

  if (ncol(side_info) > 1) {
    first_nn <- side_info[closest_found, , drop = FALSE]
    colnames(first_nn) <- paste0("nn_", colnames(first_nn))

    mean_vals <- data.table(mean_standardized_rmsd = mean(rslt[, "rmsd"] / get_column_sds(side_info[complete.cases(side_info), ])), mean_r = mean(rslt[, "r"]))
    final_result <- list(eval = rslt, global_eval = mean_vals, first_nn = cbind(side_info = side_info, first_nn))
  } else {
    final_result <- list(eval = rslt, first_nn = cbind(side_info = side_info, first_nn = side_info[closest_found, , drop = FALSE]))
  }

  final_result
}


#' @title get the evaluation results for categorical data
#' @description internal
#' @keywords internal
get_eval_categorical <- function(y, indices_closest) {
  tab <- as.matrix(table(y[!is.na(y)], y[indices_closest][!is.na(y)]))
  total <- sum(tab)
  tab <- tab / total
  p <- rowSums(tab) %*% t(colSums(tab))
  pra <- sum(diag(tab))
  pre <- sum(diag(p))
  kappa <- t((pra - pre) / (1 - pre))
  colnames(kappa) <- "kappa"

  kappa
}

#' @title get the evaluation results for continuous data
#' @description internal
#' @keywords internal
#' @importFrom stats cor
get_eval_continuous <- function(y, indices_closest) {
  get_ith_eval <- function(ith_y, indices_closest, ith_iter) {
    rmsd <- (mean((ith_y[indices_closest, ith_iter] - ith_y[, ith_iter])^2, na.rm = TRUE))^0.5
    r <- cor(ith_y[indices_closest, ith_iter], ith_y[, ith_iter], use = "complete.obs")
    c(rmsd = rmsd, r = r)
  }
  eval_rslt <- lapply(1:ncol(y), FUN = get_ith_eval, ith_y = y, indices_closest = indices_closest)
  eval_rslt <- do.call("rbind", eval_rslt)

  eval_rslt
}
