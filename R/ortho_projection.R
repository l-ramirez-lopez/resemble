#' @title Orthogonal projections using principal component analysis and partial
#' least squares
#' @aliases ortho_projection
#' @aliases pls_projection
#' @aliases pc_projection
#' @aliases predict.ortho_projection
#' @description
#' \loadmathjax
#' Functions to perform orthogonal projections of high dimensional data matrices
#' using principal component analysis (pca) and partial least squares (pls).
#' @usage
#' ortho_projection(Xr, Xu = NULL,
#'                  Yr = NULL,
#'                  method = "pca",
#'                  pc_selection = list(method = "cumvar", value = 0.99),
#'                  center = TRUE, scale = FALSE, ...)
#'
#' pc_projection(Xr, Xu = NULL, Yr = NULL,
#'               pc_selection = list(method = "cumvar", value = 0.99),
#'               center = TRUE, scale = FALSE,
#'               method = "pca",
#'               tol = 1e-6, max_iter = 1000, ...)
#'
#' pls_projection(Xr, Xu = NULL, Yr,
#'                pc_selection = list(method = "opc", value = min(dim(Xr), 40)),
#'                scale = FALSE,
#'                tol = 1e-6, max_iter = 1000, ...)
#'
#' \method{predict}{ortho_projection}(object, newdata, ...)
#'
#' @param Xr a matrix of observations.
#' @param Xu an optional matrix containing data of a second set of observations.
#' @param Yr if the method used in the \code{pc_selection} argument is \code{"opc"}
#' or if \code{method = "pls"}, then it must be a matrix
#' containing the side information corresponding to the spectra in \code{Xr}.
#' It is equivalent to the \code{side_info} parameter of the \code{\link{sim_eval}}
#' function. In case \code{method = "pca"}, a matrix (with one or more
#' continuous variables) can also be used as input. The root mean square of
#' differences (rmsd) is used for assessing the similarity between the observations
#' and their corresponding most similar observations in terms of the side information
#' provided. A single discrete variable of class factor can also be passed. In
#' that case, the kappa index is used. See \code{\link{sim_eval}} function for more details.
#' @param method the method for projecting the data. Options are: "pca" (principal
#' component analysis using the singular value decomposition algorithm),
#' "pca.nipals" (principal component analysis using the non-linear iterative
#' partial least squares algorithm) and "pls" (partial least squares).
#' @param pc_selection a list of length 2 which specifies the method to be used
#' for optimizing the number of components (principal components or pls factors)
#' to be retained. This list must contain two elements (in the following order):
#' \code{method} (a character indicating the method for selecting the number of
#' components) and \code{value} (a numerical value that complements the selected
#' method). The methods available are:
#' \itemize{
#'        \item{\code{"opc"}:} { optimized principal component selection based on
#'        Ramirez-Lopez et al. (2013a, 2013b). The optimal number of components
#'        of a given set of observations is the one for which its distance matrix
#'        minimizes the differences between the \code{Yr} value of each
#'        observation and the \code{Yr} value of its closest observation. In this
#'        case \code{value} must be a value (larger than 0 and
#'        below \code{min(nrow(Xr)} \code{+ nrow(Xu),} \code{ncol(Xr))} indicating 
#'        the maximum number of principal components to be tested. See details.}
#'
#'        \item{\code{"cumvar"}:}{ selection of the principal components based
#'        on a given cumulative amount of explained variance. In this case,
#'        \code{value} must be a value (larger than 0 and below or equal to 1)
#'        indicating the maximum amount of cumulative variance that the
#'        retained components should explain.}
#'
#'        \item{\code{"var"}:}{ selection of the principal components based
#'        on a given amount of explained variance. In this case,
#'        \code{value} must be a value (larger than 0 and below or equal to 1)
#'        indicating the minimum amount of variance that a component should
#'        explain in order to be retained.}
#'
#'        \item{\code{"manual"}:}{ for manually specifying a fix number of
#'        principal components. In this case, \code{value} must be a value
#'        (larger than 0 and below \code{min(nrow(Xr)} \code{+ nrow(Xu),} \code{ncol(Xr))}).
#'        indicating the minimum amount of variance that a component should
#'        explain in order to be retained.}
#'        }
#' The default list passed is \code{list(method = "cumvar", value = 0.99)}.
#' Optionally, the \code{pc_selection} argument admits \code{"opc"} or
#' \code{"cumvar"} or \code{"var"} or \code{"manual"} as a single character
#' string. In such a case the default \code{"value"} when either \code{"opc"} or
#' \code{"manual"} are used is 40. When \code{"cumvar"} is used the default
#' \code{"value"} is set to 0.99 and when \code{"var"} is used the default
#' \code{"value"} is set to 0.01.
#' @param center a logical indicating if the data \code{Xr} (and \code{Xu} if
#' specified) must be centered. If \code{Xu} is specified the data is centered
#' on the basis of \mjeqn{Xr \cup Xu}{Xr U Xu}. NOTE: This argument only applies to the
#' principal components projection. For pls projections the data is always
#' centered.
#' @param scale a logical indicating if \code{Xr} (and \code{Xu} if specified)
#' must be scaled. If \code{Xu} is specified the data is scaled on the basis of
#' \mjeqn{Xr \cup Xu}{Xr U Xu}.
#' @param tol tolerance limit for convergence of the algorithm in the nipals
#' algorithm (default is 1e-06). In the case of PLS this applies only to Yr with
#' more than one variable.
#' @param max_iter maximum number of iterations (default is 1000). In the case of
#' \code{method = "pls"} this applies only to \code{Yr} matrices with more than
#' one variable.
#' @param ... additional arguments to be passed from \code{ortho_projection}
#' to \code{pc_projection} or \code{pls_projection}.
#' @param object object of class "ortho_projection" (as returned by
#' \code{ortho_projection}, \code{pc_projection} or \code{pls_projection}).
#' @param newdata an optional data frame or matrix in which to look for variables
#' with which to predict. If omitted, the scores are used. It must contain the
#' same number of columns, to be used in the same order.
#' @details
#' In the case of \code{method = "pca"}, the algrithm used is the singular value
#' decomposition in which a given data matrix (\mjeqn{X}{X}) is factorized as follows:
#'      
#'  \mjdeqn{X = UDV^{T}}{X = UDV^{\mathrm{T}}}
#'      
#' where \mjeqn{U}{U} and \mjeqn{V}{V} are orthogonal matrices, being the left and right
#' singular vectors of \mjeqn{X}{X} respectively, \mjeqn{D}{D} is a diagonal matrix
#' containing the singular values of \mjeqn{X}{X} and \mjeqn{V}{V} is the is a matrix of
#' the right singular vectors of \mjeqn{X}{X}.
#' The matrix of principal component scores is obtained by a matrix
#' multiplication of \mjeqn{U}{U} and \mjeqn{D}{D}, and the matrix of principal component
#' loadings is equivalent to the matrix \mjeqn{V}{V}.
#'
#' When \code{method = "pca.nipals"}, the algorithm used for principal component
#' analysis is the non-linear iterative partial least squares (nipals).
#'
#' In the case of the of the partial least squares projection (a.k.a projection
#' to latent structures) the nipals regression algorithm is used. Details on the "nipals"
#' algorithm are presented in Martens (1991).
#'
#' When \code{method = "opc"}, the selection of the components is carried out by
#' using an iterative method based on the side information concept
#' (Ramirez-Lopez et al. 2013a, 2013b). First let be \mjeqn{P}{P} a sequence of
#' retained components (so that \mjeqn{P = 1, 2, ...,k }{P = 1, 2, ...,k }).
#' At each iteration, the function computes a dissimilarity matrix retaining
#' \mjeqn{p_i}{p_i} components. The values in this side information variable are
#' compared against the side information values of their most spectrally similar
#' observations (closest \code{Xr} observation).
#' The optimal number of components retrieved by the function is the one that
#' minimizes the root mean squared differences (RMSD) in the case of continuous
#' variables, or maximizes the kappa index in the case of categorical variables.
#' In this process, the \code{\link{sim_eval}} function is used.
#' Note that for the \code{"opc"} method \code{Yr} is required (i.e. the
#' side information of the observations).
#'
#' This function supports multi-threading for the computation of dissimilarities
#' via OpenMP in Rcpp.
#' @return \code{ortho_projection}, \code{pc_projection}, \code{pls_projection},
#' return a \code{list} of class \code{ortho_projection} with the following
#' components:
#' \itemize{
#'  \item{\code{scores}}{ a matrix of scores corresponding to the observations in
#'  \code{Xr} (and \code{Xu} if it was provided). The components retrieved
#'  correspond to the ones optimized or specified.}
#'  \item{\code{X_loadings}}{ a matrix of loadings corresponding to the
#'  explanatory variables. The components retrieved correspond to the ones
#'  optimized or specified.}
#'  \item{\code{Y_loadings}}{ a matrix of partial least squares loadings
#'  corresponding to \code{Yr}. The components retrieved  correspond to the
#'  ones optimized or specified.
#'  This object is only returned if the partial least squares algorithm was used.}
#'  \item{\code{weigths}}{ a matrix of partial least squares ("pls") weights.
#'  This object is only returned if the "pls" algorithm was used.}
#'  \item{\code{projection_mat}}{ a matrix that can be used to project new data
#'  onto a "pls" space. This object is only returned if the "pls" algorithm was
#'  used.}
#'  \item{\code{variance}}{ a matrix indicating the standard deviation of each
#'  component (sd), the cumulative explained variance (cumulative_explained_var) and the
#'  variance explained by each single component (explained_var). These values are
#'  computed based on the data used to create the projection matrices.
#'  For example if the "pls" method was used, then these values are computed
#'  based only on the data that contains information on \code{Yr} (i.e. the
#'  \code{Xr} data). If the principal component method is used, the this data is
#'  computed on the basis of \code{Xr} and \code{Xu} (if it applies) since both
#'  matrices are employed in the computation of the projection matrix (loadings
#'  in this case)}.
#'  \item{\code{sdv}}{ the standard deviation of the retrieved scores. This vector
#'  can be different from the "sd" in \code{variance}.}
#'  \item{\code{n_components}}{ the number of components (either principal
#'  components or partial least squares components) used for computing the
#'  global dissimilarity scores.}
#'  \item{\code{opc_evaluation}}{ a matrix containing the statistics computed
#'  for optimizing the number of principal components based on the variable(s)
#'  specified in the \code{Yr} argument. If \code{Yr} was a continuous  was a
#'  continuous vector or matrix then this object indicates the root mean square
#'  of differences (rmse) for each number of components. If \code{Yr} was a
#'  categorical variable this object indicates the kappa values for each number
#'  of components. This object is returned only if \code{"opc"} was used within
#'  the \code{pc_selection} argument. See the \code{\link{sim_eval}} function for
#'  more details.}
#'  \item{\code{method}}{ the \code{ortho_projection} method used.}
#'  }
#'  \code{predict.ortho_projection}, returns a matrix of scores proprojected for
#'  \code{newdtata}.
#' @author Leonardo Ramirez-Lopez
#' @references
#' Martens, H. (1991). Multivariate calibration. John Wiley & Sons.
#'
#' Ramirez-Lopez, L., Behrens, T., Schmidt, K., Stevens, A., Dematte, J.A.M.,
#' Scholten, T. 2013a. The spectrum-based learner: A new local approach for
#' modeling soil vis-NIR spectra of complex data sets. Geoderma 195-196, 268-279.
#'
#' Ramirez-Lopez, L., Behrens, T., Schmidt, K., Viscarra Rossel, R., Dematte,
#' J. A. M.,  Scholten, T. 2013b. Distance and similarity-search metrics for use
#' with soil vis-NIR spectra. Geoderma 199, 43-53.
#' @seealso \code{\link{ortho_diss}}, \code{\link{sim_eval}}, \code{\link{mbl}}
#' @examples
#' \dontrun{
#' library(prospectr)
#' data(NIRsoil)
#'
#' Xu <- NIRsoil$spc[!as.logical(NIRsoil$train), ]
#' Yu <- NIRsoil[!as.logical(NIRsoil$train), "CEC", drop = FALSE]
#' Yr <- NIRsoil[as.logical(NIRsoil$train), "CEC", drop = FALSE]
#' Xr <- NIRsoil$spc[as.logical(NIRsoil$train), ]
#'
#' Xu <- Xu[!is.na(Yu), ]
#' Yu <- Yu[!is.na(Yu), , drop = FALSE]
#'
#' Xr <- Xr[!is.na(Yr), ]
#' Yr <- Yr[!is.na(Yr), , drop = FALSE]
#'
#' # A principal component analysis using 5 components
#' pca_projected <- ortho_projection(Xr, pc_selection = list("manual", 5))
#' pca_projected
#'
#' # A principal components projection using the "opc" method
#' # for the selection of the optimal number of components
#' pca_projected_2 <- ortho_projection(
#'   Xr, Xu, Yr,
#'   method = "pca",
#'   pc_selection = list("opc", 40)
#' )
#' pca_projected_2
#' plot(pca_projected_2)
#'
#' # A partial least squares projection using the "opc" method
#' # for the selection of the optimal number of components
#' pls_projected <- ortho_projection(
#'   Xr, Xu, Yr,
#'   method = "pls",
#'   pc_selection = list("opc", 40)
#' )
#' pls_projected
#' plot(pls_projected)
#'
#' # A partial least squares projection using the "cumvar" method
#' # for the selection of the optimal number of components
#' pls_projected_2 <- ortho_projection(
#'   Xr = Xr, Yr = Yr, Xu = Xu,
#'   method = "pls",
#'   pc_selection = list("cumvar", 0.99)
#' )
#' }
#' @rdname ortho_projection
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
## 10.03.2014 Leo     Sanity check for the number of cases in Yr and Xr
## 13.03.2014 Antoine The explanation of the cores argument was modified
## 03.12.2015 Leo     The pls projection function is now based on Rcpp
## 04.12.2015 Leo     The center argument was removed from the pls projection
##                    function. Matrices are now always centered.
## 04.12.2015 Leo     the dist function needs to be replaced by the fDiss.
##                    dist() retrieves squared distances!
## 07.12.2015 Leo     The dist function was removed from these functions (it was
##                    implemented)
##                    for testing for a while but it was not finally released.
## 07.12.2015 Leo     The pc_projection, pls_projection and the predict.othoProjection
##                    are now visible.
## 08.09.2016 Leo     A bug in the predict.ortho_projection function was fixed.
##                    A an error was thrown when PCA preditions were requested.
##                    Furthermore, the PLS score predictions were hadnled wrognly
##                    as PCA preditions.
## 13.09.2016 Leo     Documentation generated by roxygen retuns duplicated
##                    pc_projection
##                    pls_projection docs. For the moment the workaround is to
##                    fix them in the Rd file of ortho_projection
## 24.07.2017 Leo     error in the predict function of ortho_projection the code
##                    was written as:
##                    if(length(grep("pca", object) == 0)) do PC projection...
##                    i.e. if the method is "pca" do not execute the PC projection.
##                    This was fixed by if(length(grep("pca", object) != 0))
##                    do PC projection
## 23.09.2018 Leo     add a newdata names check in predict.ortho_projection..
##                    they need to correspond to variable names in ortho_projection
##                    object
## 03.07.2019 Leo     there was a problem retriving the optimal number of pcs in
##                    pc_projection when the 'opc'method was used in combination
##                    with missing data in Yr. The neighbors retrieved were not
##                    excluding missing values generating a miscalcilation of the
##                    rmsd. This was solved by removing missing values before
##                    the calculations with the fast_diss_vector.
## 22.04.2020 Leo     Argument scaled renamed to scale
##                    Argument max_compter renamed to max_iter
##                    the call is now always part of the error message (this was
##                    not the case when the function was being called from inside
##                    other functions)
##                    Helpers were added (see ortho_helpers.R)
## 30.04.2020 Leo     pca.nipals implemented in c++
##                    cores is deprecated
ortho_projection <- function(Xr, Xu = NULL,
                             Yr = NULL,
                             method = "pca",
                             pc_selection = list(method = "cumvar", value = 0.99),
                             center = TRUE, scale = FALSE, ...) {
  method <- match.arg(method, c("pls", "pca", "pca.nipals"))

  if (method == "pls") {
    if (!is.numeric(as.matrix(Yr))) {
      stop("When pls projection is used, 'Yr' must be numeric")
    }
    proj <- pls_projection(
      Xr = Xr, Yr = Yr, Xu = Xu, pc_selection = pc_selection,
      scale = scale, ...
    )
    mthd <- "pls"
  } else {
    mthd <- if_else(method == "pca", "pca (svd)", "pca (nipals)")
    proj <- pc_projection(
      Xr = Xr, Yr = Yr, Xu = Xu, pc_selection = pc_selection,
      center = center, scale = scale, method = method, ...
    )
  }

  proj$method <- mthd
  class(proj) <- c("ortho_projection", "list")

  proj
}


#' @aliases ortho_projection
#' @aliases pls_projection
#' @aliases pc_projection
#' @aliases predict.ortho_projection
#' @importFrom stats quantile complete.cases diffinv
#' @export pc_projection
pc_projection <- function(Xr, Xu = NULL, Yr = NULL,
                          pc_selection = list(method = "cumvar", value = 0.99),
                          center = TRUE, scale = FALSE,
                          method = "pca",
                          tol = 1e-6, max_iter = 1000, ...) {
  pc_selection_method <- pc_selection[[1]]
  match.arg(pc_selection_method, c("opc", "var", "cumvar", "manual"))

  match.arg(method, c("pca", "pca.nipals"))

  if (!is.logical(center)) {
    stop("'center' must be logical")
  }

  if (!is.logical(scale)) {
    stop("'scale' must be logical")
  }

  if (!is.null(Yr)) {
    Yr <- as.matrix(Yr)
    # if (!is.matrix(Yr)) {
    #   stop("Yr must be a matrix")
    # }
    if (nrow(Yr) != nrow(Xr)) {
      stop(paste0(
        "The number of observations in 'Yr' does not match with the ",
        "number of observations in 'Xr'"
      ))
    }
  }
  
  ny <- ncol(Yr)

  if (!is.null(Xu)) {
    if (ncol(Xr) != ncol(Xu)) {
      stop("Number of columns in 'Xr' and 'Xu' do not match")
    }
  }

  effective_rows_xr <- nrow(Xr)
  n_cols_xr <- ncol(Xr)
  # here ifelse is preferred over if_else as the later returns an error
  # because it evaluates nrow(Xu) when XU == NULL
  effective_rows_xu <- ifelse(is.null(Xu), 0, nrow(Xu))
  Xr <- rbind(Xr, Xu)
  dparam <- check_pc_arguments(
    n_rows_x = effective_rows_xr + effective_rows_xu,
    n_cols_x = n_cols_xr,
    pc_selection = pc_selection,
    default_max_comp = 40,
    default_max_cumvar = 0.99,
    default_max_var = 0.01
  )
  pc_selection <- dparam$pc_selection_checked
  max_comp <- dparam$max_comp

  pc_selection_copy <- pc_selection

  # center
  if (center) {
    mean_vector <- colMeans(Xr)
    X0 <- sweep(x = Xr, MARGIN = 2, FUN = "-", STATS = mean_vector)
  } else {
    mean_vector <- rep(0, ncol(Xr))
    X0 <- Xr
  }

  if (scale) {
    sd_vector <- get_column_sds(X0)
    X0 <- sweep(x = X0, MARGIN = 2, FUN = "/", STATS = sd_vector)
  } else {
    sd_vector <- rep(1, ncol(X0))
  }

  if (method == "pca") {
    sv_decomposition <- svd(x = X0, nu = max_comp, nv = max_comp)
    sv_decomposition$d <- sv_decomposition$d[1:max_comp]
    # Loadings and scores
    pc_loadings <- t(sv_decomposition$v)
    pc_scores <- sv_decomposition$u %*% diag(sv_decomposition$d)
    # Variance of each PC variable
    sdPC <- get_column_sds(pc_scores)
    # Compute the percentage of explained variance for all the PCs
    ons <- (sv_decomposition$d)^2 / (nrow(X0) - 1)
    explained_v <- ons / sum(ons)
    cummulative_v <- diffinv(explained_v[-1], xi = explained_v[1])
    variance <- rbind(
      sd = as.vector(sdPC),
      cumulative_explained_var = cummulative_v,
      explained_var = explained_v
    )
  }

  if (method == "pca.nipals") {
    nipals_pca <- pca_nipals(
      X = X0,
      ncomp = max_comp,
      center = FALSE,
      scale = FALSE,
      maxiter = max_iter,
      tol = tol,
      pcSelmethod = pc_selection_method,
      pcSelvalue = pc_selection_copy$value
    )

    pc_scores <- nipals_pca$pc_scores
    pc_loadings <- nipals_pca$pc_loadings
    explained_v <- nipals_pca$pc_variance
    variance <- rbind(
      sd = explained_v[1, ],
      cumulative_explained_var = explained_v[2, ],
      explained_var = explained_v[3, ]
    )
    colnames(variance) <- paste0("pc_", 1:ncol(variance))
  }

  # assign names
  colnames(pc_scores) <- paste0("pc_", 1:ncol(pc_scores))
  rownames(pc_scores) <- c(
    paste0("Xr_", 1:effective_rows_xr),
    if (!is.null(Xu)) {
      paste0("Xu_", 1:effective_rows_xu)
    }
  )
  colnames(pc_loadings) <- colnames(X0)
  rownames(pc_loadings) <- paste0("pc_", 1:nrow(pc_loadings))

  if (pc_selection_method == "opc") {
    if (is.null(Yr) | !is.matrix(Yr)) {
      stop(paste0(
        "'Yr' msut be a matrix when the 'opc' method is used for ",
        "selecting the optimal number of principal components"
      ))
    }
    if (nrow(Yr) != effective_rows_xr) {
      stop("The number of rows in Xr does not match the number of cases in Yr")
    }
    if (!is.null(colnames(Yr))) {
      if (sum(duplicated(colnames(Yr))) > 0) {
        stop("column names in Yr must be different")
      }
    } else {
      colnames(Yr) <- paste0("Yr_", 1:ny)
    }
    results <- eval_multi_pc_diss(pc_scores[, 1:max_comp],
      side_info = Yr,
      method = "pc",
      check_dims = FALSE
    )
    selected_pcs <- results$best_pc
    results <- results$result
  }

  if (pc_selection_method == "cumvar") {
    selected_pcs <- variance[2, ] <= pc_selection_copy$value
    selected_pcs <- sum(selected_pcs)
  }

  if (pc_selection_method == "var") {
    selected_pcs <- variance[3, ] >= pc_selection_copy$value
    selected_pcs <- sum(selected_pcs)
  }

  if (pc_selection_method == "manual") {
    selected_pcs <- (1:ncol(pc_scores)) <= pc_selection_copy$value
    selected_pcs <- sum(selected_pcs)
  }

  if (pc_selection_method == "opc") {
    selected_pcs <- sum(selected_pcs)
  }

  scores_sd <- get_column_sds(pc_scores[, 1:selected_pcs, drop = FALSE])
  colnames(scores_sd) <- colnames(pc_scores[, 1:selected_pcs, drop = FALSE])
  rownames(scores_sd) <- "sd"

  fresults <- list(
    scores = pc_scores[, 1:selected_pcs, drop = FALSE],
    X_loadings = pc_loadings[1:selected_pcs, , drop = FALSE],
    variance = variance[, 1:selected_pcs, drop = FALSE],
    scores_sd = scores_sd,
    n_components = selected_pcs, pc_selection = pc_selection,
    center = mean_vector,
    scale = sd_vector
  )
  colnames(fresults$variance) <- rownames(fresults$X_loadings)
  fresults$method <- if_else(method == "pca", "pca (svd)", "pca (nipals)")
  if (pc_selection_method == "opc") {
    fresults$opc_evaluation <- results
  }
  class(fresults) <- c("ortho_projection", "list")

  fresults
}


#' @aliases ortho_projection
#' @aliases pls_projection
#' @aliases pc_projection
#' @aliases predict.ortho_projection
#' @export pls_projection
pls_projection <- function(Xr, Xu = NULL, Yr,
                           pc_selection = list(method = "opc", value = min(dim(Xr), 40)),
                           scale = FALSE, tol = 1e-6, max_iter = 1000, ...) {
  pc_selection_method <- pc_selection[[1]]
  match.arg(pc_selection_method, c("opc", "var", "cumvar", "manual"))

  if (!is.logical(scale)) {
    stop("'scale' argument must be logical")
  }


  if (missing(Yr)) {
    stop("'Yr' must be provided")
  } else {
    Yr <- as.matrix(Yr)
    if (!is.numeric(Yr)) {
      stop("'Yr' must be a numeric matrix")
    }
  }

  if (!is.null(Xu)) {
    if (ncol(Xr) != ncol(Xu)) {
      stop("Number of columns in 'Xr' and 'Xu' do not match")
    }
  }

  if (is.null(pc_selection$value)) {
    pc_selection_value <- pc_selection[[2]]
  } else {
    pc_selection_value <- pc_selection$value
  }

  nas <- rowSums(is.na(Yr)) > 0

  ny <- ncol(Yr)

  X0 <- Xr
  Y0 <- Yr
  non_nas_yr <- 1:nrow(Xr)
  if (sum(nas) > 0) {
    nas_yr <- (1:nrow(Xr))[nas]
    non_nas_yr <- (1:nrow(Xr))[!nas]
    Xout <- Xr[nas_yr, , drop = FALSE]
    X0 <- Xr[non_nas_yr, ]
    Y0 <- Yr[non_nas_yr, , drop = FALSE]
    if (pc_selection_method %in% c("opc", "manual")) {
      if (min(dim(X0)) < pc_selection_value) {
        stop(paste0(
          "Missing values in Yr. The number of components specified ",
          "in 'pc_selection' exceeds the number of observations with ",
          "non-missing Yr values. Try another number of components."
        ))
      }
    }
  }


  effective_rows_xr <- nrow(X0)
  n_cols_xr <- ncol(X0)

  dparam <- check_pc_arguments(
    n_rows_x = effective_rows_xr,
    n_cols_x = n_cols_xr,
    pc_selection = pc_selection,
    default_max_comp = 40,
    default_max_cumvar = 0.99,
    default_max_var = 0.01
  )
  pc_selection <- dparam$pc_selection_checked
  max_comp <- dparam$max_comp

  weights <- matrix(NA, max_comp, ncol(X0))
  scores <- matrix(NA, nrow(X0), max_comp)
  X_loadings <- matrix(NA, max_comp, ncol(X0))
  Y_loadings <- matrix(NA, max_comp, ny)
  pls_variance <- matrix(NA, 3, max_comp)

  if (pc_selection_method %in% c("opc", "manual")) {
    pc_selection$value <- pc_selection_value - 1
    cpp_method <- "manual"
  } else {
    cpp_method <- pc_selection_method
  }

  plsp <- opls_for_projection(
    X = X0,
    Y = as.matrix(Y0),
    ncomp = max_comp,
    scale = scale,
    maxiter = max_iter,
    tol = tol,
    pcSelmethod = cpp_method,
    pcSelvalue = pc_selection$value
  )
  max_comp <- plsp$ncomp

  if (pc_selection_method == "opc") {
    if (is.null(Yr) | !is.matrix(Yr)) {
      stop(paste0(
        "'Yr' msut be a matrix when the 'opc' method is used for ",
        "selecting the optimal number of principal components"
      ))
    }
    if (nrow(Yr) != nrow(Xr)) {
      stop("The number of rows in Xr does not match the number of cases in Yr")
    }
    if (!is.null(colnames(Yr))) {
      if (sum(duplicated(colnames(Yr))) > 0) {
        stop("column names in Yr must be different")
      }
    } else {
      colnames(Y0) <- colnames(Yr) <- paste0("Yr_", 1:ny)
    }
    results <- eval_multi_pc_diss(plsp$scores[, 1:max_comp],
      side_info = as.matrix(Y0),
      method = "pls",
      check_dims = FALSE
    )
    max_comp <- results$best_pc
    results <- results$result
  }


  # Select the necessary components
  pls_variance <- plsp$variance$x_var[, 1:max_comp, drop = FALSE]
  weights <- plsp$weights[1:max_comp, , drop = FALSE]
  scores <- plsp$scores[, 1:max_comp, drop = FALSE]
  X_loadings <- plsp$X_loadings[1:max_comp, , drop = FALSE]
  Y_loadings <- plsp$Y_loadings[1:max_comp, , drop = FALSE]
  plsp$projection_mat <- plsp$projection_mat[, 1:max_comp, drop = FALSE]

  # Give some names...
  colnames(X_loadings) <- colnames(X0)
  rownames(X_loadings) <- colnames(plsp$projection_mat) <- paste0("pls_", 1:max_comp)
  rownames(Y_loadings) <- rownames(X_loadings)
  colnames(Y_loadings) <- paste0("Y_pls", 1:ny)
  rownames(weights) <- rownames(X_loadings)
  colnames(weights) <- colnames(X_loadings)
  colnames(pls_variance) <- rownames(X_loadings)
  rownames(pls_variance) <- c("sd", "cumulative_explained_var_X", "explained_var_X")

  yex <- plsp$variance$y_var[, 1:max_comp, drop = FALSE]

  colnames(yex) <- rownames(Y_loadings)
  if (ny > 1) {
    rownames(yex) <- paste0("explained_var_", colnames(Yr))
  } else {
    rownames(yex) <- "explained_var_Yr"
  }

  if (sum(nas) > 0) {
    scores.a <- matrix(NA, length(c(non_nas_yr, nas_yr)), ncol(scores))
    scores.a[nas_yr, ] <- project_opls(
      projection_mat = plsp$projection_mat,
      ncomp = max_comp,
      newdata = Xout,
      scale = scale,
      Xcenter = plsp$transf$Xcenter,
      Xscale = plsp$transf$Xscale
    )
    scores.a[non_nas_yr, ] <- scores
    scores <- scores.a
  }

  rownames(scores) <- paste0("Xr_", 1:nrow(scores))
  colnames(scores) <- rownames(Y_loadings)


  if (!is.null(Xu)) {
    if (is.vector(Xu)) {
      Xu <- t(Xu)
    }
    scores_Xu <- project_opls(
      projection_mat = plsp$projection_mat,
      ncomp = max_comp,
      newdata = Xu,
      scale = scale,
      Xcenter = plsp$transf$Xcenter,
      Xscale = plsp$transf$Xscale
    )

    colnames(scores_Xu) <- rownames(Y_loadings)
    rownames(scores_Xu) <- paste0("Xu_", 1:nrow(Xu))
    scores <- rbind(scores, scores_Xu)
  }

  if (!nrow(plsp$transf$Xscale)) {
    plsp$transf$Xscale <- matrix(1, 1, length(plsp$transf$Xcenter))
  }

  scores_sd <- get_column_sds(scores)
  colnames(scores_sd) <- colnames(scores)
  rownames(scores_sd) <- "sd"
  fresults <- list(
    scores = scores,
    X_loadings = X_loadings,
    Y_loadings = Y_loadings,
    weights = weights,
    projection_mat = plsp$projection_mat[, 1:ncol(scores)],
    variance = list(x_var = pls_variance, y_var = yex),
    scores_sd = scores_sd,
    n_components = max_comp,
    pc_selection = pc_selection,
    center = plsp$transf$Xcenter,
    scale = plsp$transf$Xscale
  )
  fresults$method <- "pls"
  if (pc_selection_method == "opc") {
    fresults$opc_evaluation <- results
  }
  class(fresults) <- c("ortho_projection", "list")

  fresults
}

#' @aliases ortho_projection
#' @aliases pls_projection
#' @aliases pc_projection
#' @aliases predict.ortho_projection
#' @export
predict.ortho_projection <- function(object, newdata, ...) {
  if (missing(newdata)) {
    return(object$scores)
  }

  nms.org <- colnames(object$X_loadings)
  nms.nd <- colnames(newdata)

  if (sum(!nms.org %in% nms.nd) != 0) {
    stop("There are missing variables in new data that are required for the projection")
  }

  else {
    if (length(grep("pca", object$method)) != 0) {
      newdata <- sweep(newdata, MARGIN = 2, FUN = "-", STATS = object$center)
      newdata <- sweep(newdata, MARGIN = 2, FUN = "/", STATS = object$scale)
      return(newdata %*% t(object$X_loadings))
    } else {
      predpoj <- project_opls(
        projection_mat = object$projection_mat,
        ncomp = ncol(object$projection_mat),
        newdata = newdata,
        scale = TRUE,
        Xcenter = object$center,
        Xscale = object$scale
      )

      colnames(predpoj) <- paste0("pls", 1:ncol(predpoj))
      rownames(predpoj) <- rownames(newdata)

      return(predpoj)
    }
  }
}
