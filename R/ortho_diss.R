#' @title A function for computing dissimilarity matrices from orthogonal
#' projections (ortho_diss)
#' @description
#' \loadmathjax
#' This function computes dissimilarities (in an orthogonal space) between
#' either observations in a given set or between observations in two different
#' sets.The dissimilarities are computed based on either principal component
#' projection or partial least squares projection of the data. After projecting
#' the data, the Mahalanobis distance is applied.
#' @usage
#' ortho_diss(Xr, Xu = NULL,
#'            Yr = NULL,
#'            pc_selection = list(method = "var", value = 0.01),
#'            diss_method = "pca",
#'            .local = FALSE,
#'            pre_k,
#'            center = TRUE,
#'            scale = FALSE,
#'            compute_all = FALSE,
#'            return_projection = FALSE,
#'            allow_parallel = TRUE, ...)
#' @param Xr a matrix containing \code{n} reference observations/rows and
#' \code{p} variables/columns.
#' @param Xu an optional matrix containing data of a second set of observations
#' with \code{p} variables/columns.
#' @param Yr a matrix of \code{n} rows and one or more columns (variables) with
#' side information corresponding to the observations in \code{Xr} (e.g. response
#' variables). It can be numeric with multiple variables/columns, or character
#' with one single column. This argument is
#' required if:
#' \itemize{
#' \item{\code{diss_method == 'pls'}: \code{Yr} is required to project the variables
#' to orthogonal directions such that the covariance between the extracted pls
#' components and \code{Yr} is maximized.}
#' \item{\code{pc_selection$method == 'opc'}: \code{Yr}  is required to optimize
#' the number of components. The optimal number of projected components is the one
#' for which its distance matrix minimizes the differences between the \code{Yr}
#' value of each observation and the \code{Yr} value of its closest observation.
#' See \code{\link{sim_eval}}.}
#' }
#' @param pc_selection a list of length 2 which specifies the method to be used
#' for optimizing the number of components (principal components or pls factors)
#' to be retained. This list must contain two elements (in the following order):
#' \code{method} (a character indicating the method for selecting the number of
#' components) and \code{value} (a numerical value that complements the selected
#' method). The methods available are:
#' \itemize{
#'        \item{\code{"opc"}:} { optimized principal component selection based on
#'        Ramirez-Lopez et al. (2013a, 2013b). The optimal number of components
#'        (of a given set of observations) is the one for which its distance
#'        matrix minimizes the differences between the \code{Yr} value of each
#'        observation and the \code{Yr} value of its closest observation. In this
#'        case, \code{value} must be a value (larger than 0 and
#'        below \code{min(nrow(Xr)} \code{+ nrow(Xu),} \code{ncol(Xr))} indicating the maximum
#'        number of principal components to be tested. See the
#'        \code{\link{ortho_projection}} function for more details.}
#'
#'        \item{\code{"cumvar"}:}{ selection of the principal components based
#'        on a given cumulative amount of explained variance. In this case,
#'        \code{value} must be a value (larger than 0 and below or equal to 1)
#'        indicating the minimum amount of cumulative variance that the
#'        combination of retained components should explain.}
#'
#'        \item{\code{"var"}:}{ selection of the principal components based
#'        on a given amount of explained variance. In this case,
#'        \code{value} must be a value (larger than 0 and below or equal to 1)
#'        indicating the minimum amount of variance that a single component
#'        should explain in order to be retained.}
#'
#'        \item{\code{"manual"}:}{ for manually specifying a fix number of
#'        principal components. In this case, \code{value} must be a value
#'        (larger than 0 and below \code{min(nrow(Xr)} \code{+ nrow(Xu),} \code{ncol(Xr))}).
#'        indicating the minimum amount of variance that a component should
#'        explain in order to be retained.}
#'        }
#' The default list passed is \code{list(method = "var", value = 0.01)}.
#' Optionally, the \code{pc_selection} argument admits \code{"opc"} or
#' \code{"cumvar"} or \code{"var"} or \code{"manual"} as a single character
#' string. In such case, the default \code{"value"} when either \code{"opc"} or
#' \code{"manual"} are used is 40. When \code{"cumvar"} is used the default
#' \code{"value"} is set to 0.99 and when \code{"var"} is used, the default
#' \code{"value"} is set to 0.01.
#' @param diss_method a character value indicating the type of projection on which
#' the dissimilarities must be computed. This argument is equivalent to
#' \code{method} argument in the \code{\link{ortho_projection}} function.
#' Options are:
#' \itemize{
#' \item{\code{"pca"}}{: principal component analysis using the singular value
#' decomposition algorithm)}
#' \item{\code{"pca.nipals"}}{: principal component analysis using
#' the non-linear iterative partial least squares algorithm.}
#' \item{\code{"pls"}}{: partial least squares.}
#' }
#' See the \code{\link{ortho_projection}} function for further details on the
#' projection methods.
#' @param .local a logical indicating whether or not to compute the dissimilarities
#' locally (i.e. projecting locally the data) by using the \code{pre_k} nearest
#' neighbor observations of each target observation. Default is \code{FALSE}. See details.
#' @param pre_k if \code{.local = TRUE} a numeric integer value which indicates the
#' number of nearest neighbors to (pre-)retain for each observation to
#' compute the (local) orthogonal dissimilarities to each observation in its
#' neighborhhod.
#' @param center a logical indicating if the \code{Xr} and \code{Xu} must be
#' centered. If \code{Xu} is provided the data is centered around the mean of
#' the pooled \code{Xr} and \code{Xu} matrices (\mjeqn{Xr \cup Xu}{Xr U Xu}). For
#' dissimilarity computations based on pls, the data is always centered for
#' the projections.
#' @param scale a logical indicating if the \code{Xr} and \code{Xu} must be
#' scaled. If \code{Xu} is provided the data is scaled based on the standard
#' deviation of the the pooled \code{Xr} and \code{Xu} matrices (\mjeqn{Xr \cup Xu}{Xr U Xu}).
#' if \code{center = TRUE}, scaling is applied after centering.
#' @param compute_all a logical. In case \code{Xu} is specified it indicates
#' whether or not the distances between all the elements resulting from the
#' pooled \code{Xr} and \code{Xu} matrices (\mjeqn{Xr \cup Xu}{Xr U Xu} must be computed).
#' @param return_projection a logical. If \code{TRUE} the `ortho_projection` object
#' on which the dissimilarities are computed will be returned. Default is \code{FALSE}. Note that
#' for \code{.local = TRUE} only the initial projection is returned (i.e. local
#' projections are not).
#' @param allow_parallel a logical (default TRUE). It allows parallel computing
#' of the local distance matrices (i.e. when \code{.local = TRUE}). This is done
#' via \code{\link[foreach]{foreach}} function of the 'foreach' package.
#' @param ... additional arguments to be passed to the
#' \code{\link{ortho_projection}} function.
#' @details
#' When \code{.local = TRUE}, first a global dissimilarity matrix is computed based on
#' the parameters specified. Then, by using this matrix for each target
#' observation, a given set of nearest neighbors (\code{pre_k}) are identified.
#' These neighbors (together with the target observation) are projected
#' (from the original data space) onto a (local) orthogonal space (using the
#' same parameters specified in the function). In this projected space the
#' Mahalanobis distance between the target observation and its neighbors is
#' recomputed. A missing value is assigned to the observations that do not belong to
#' this set of neighbors (non-neighbor observations).
#' In this case the dissimilarity matrix cannot be considered as a distance
#' metric since it does not necessarily satisfies the symmetry condition for
#' distance matrices (i.e. given two observations \mjeqn{x_i}{x_i} and \mjeqn{x_j}{x_j}, the local
#' dissimilarity (\mjeqn{d}{d}) between them is relative since generally
#' \mjeqn{d(x_i, x_j) \neq d(x_j, x_i)}{d(x_i, x_j) ne d(x_j, x_i)}). On the other hand, when
#' \code{.local = FALSE}, the dissimilarity matrix obtained can be considered as
#' a distance matrix.
#'
#' In the cases where \code{"Yr"} is required to compute the dissimilarities and
#' if \code{.local = TRUE}, care must be taken as some neighborhoods might
#' not have enough observations with non-missing \code{"Yr"} values, which might retrieve
#' unreliable dissimilarity computations.
#'
#' If \code{.local = TRUE} and \code{pc_selection$method} is \code{"opc"} or
#' \code{"manual"}, the minimum number of observations with non-missing \code{"Yr"}
#' values at each neighborhood is determined by \code{pc_selection$value}
#' (i.e. the maximum number of components to compute).
#'
#'
#' @return a \code{list} of class \code{ortho_diss} with the following elements:
#' \itemize{
#'  \item{\code{n_components}}{ the number of components (either principal
#'  components or partial least squares components) used for computing the
#'  global dissimilarities.}
#'  \item{\code{global_variance_info}}{ the information about the expalined
#'  variance(s) of the projection. When \code{.local = TRUE}, the information
#'  corresponds to the global projection done prior computing the local
#'  projections.}
#'  \item{\code{local_n_components}}{ if \code{.local = TRUE}, a data.table
#'  which specifies the number of local components (either principal components
#'  or partial least squares components) used for computing the dissimilarity
#'  between each target observation and its neighbor observations.}
#'  \item{\code{dissimilarity}}{ the computed dissimilarity matrix. If
#'  \code{.local = FALSE} a distance matrix. If \code{.local = TRUE} a matrix of
#'  class \code{local_ortho_diss}. In this case, each column represent the dissimilarity
#'  between a target observation and its neighbor observations.}
#'  \item{\code{projection}}{if \code{return_projection = TRUE},
#'  an \code{ortho_projection} object.}
#'  }
#' @author \href{https://orcid.org/0000-0002-5369-5120}{Leonardo Ramirez-Lopez}
#' @references
#' Ramirez-Lopez, L., Behrens, T., Schmidt, K., Stevens, A., Dematte, J.A.M.,
#' Scholten, T. 2013a. The spectrum-based learner: A new local approach for
#' modeling soil vis-NIR spectra of complex data sets. Geoderma 195-196, 268-279.
#'
#' Ramirez-Lopez, L., Behrens, T., Schmidt, K., Viscarra Rossel, R., Dematte,
#' J. A. M.,  Scholten, T. 2013b. Distance and similarity-search metrics for use
#' with soil vis-NIR spectra. Geoderma 199, 43-53.
#' @seealso \code{\link{ortho_projection}}, \code{\link{sim_eval}}
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
#' # Computation of the orthogonal dissimilarity matrix using the
#' # default parameters
#' pca_diss <- ortho_diss(Xr, Xu)
#'
#' # Computation of a principal component dissimilarity matrix using
#' # the "opc" method for the selection of the principal components
#' pca_diss_optim <- ortho_diss(
#'   Xr, Xu, Yr,
#'   pc_selection = list("opc", 40),
#'   compute_all = TRUE
#' )
#'
#' # Computation of a partial least squares (PLS) dissimilarity
#' # matrix using the "opc" method for the selection of the PLS
#' # components
#' pls_diss_optim <- ortho_diss(
#'   Xr = Xr, Xu = Xu,
#'   Yr = Yr,
#'   pc_selection = list("opc", 40),
#'   diss_method = "pls"
#' )
#' }
#' @export

#######################################################################
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
#######################################################################

## History:
## 09.03.2014 Leo     In the doc was specified that multi-threading is
##                    not working for mac
## 13.03.2014 Antoine The explanation of the cores argument was modified
## 07.09.2014 Antoine A bug handling Yr as a matrix was fixed
## 02.12.2015 Leo     The function now outputs an the object global_variance_info
##                    which provides information on the explained variance
##                    of the projections.
## 10.08.2018 Leo     A wrong message error was found and corrceted. previously was:
##                    "Yu must be provided either when the 'opc' is used in pc_selection is used or method = 'pls'"
##                    In fact it is not Yu but Yr ( the correct argument)
## 01.05.2020 Leo     - Argument scaled renamed to scale
##                    - Argument return.all renamed to compute_all
##                    - Refined documentation
##                    - Argument cores is deprecated
##                    - When \code{.local = TRUE} a new output is produced:
##                      'neighborhood_info' which is a data.table containing
##                      relevant information about the neighborhood of
##                      each observation (e.g. neighborhood indices, number of
##                      components used at each neighborhood, etc)
##                    - Output global.variance.info has been renamed to
##                      global_variance_info

ortho_diss <- function(Xr, Xu = NULL,
                       Yr = NULL,
                       pc_selection = list(method = "var", value = 0.01),
                       diss_method = "pca",
                       .local = FALSE,
                       pre_k,
                       center = TRUE, scale = FALSE,
                       compute_all = FALSE,
                       return_projection = FALSE,
                       allow_parallel = TRUE,
                       ...) {
  if (!is.logical(.local)) {
    stop("'.local'must be logical")
  }

  if (!is.logical(compute_all)) {
    stop("'compute_all' must be logical")
  }

  if (!is.logical(center)) {
    stop("'center' must be logical")
  }

  if (!is.logical(scale)) {
    stop("'scale' must be logical")
  }

  if (.local) {
    if (pre_k < 4) {
      stop("pre_k cannot be smaller than 4")
    }
    if (pre_k > nrow(Xr)) {
      stop("pre_k cannot be larger than the numer of elements in Xr")
    }
  }

  if (!is.null(Yr)) {
    Yr <- as.matrix(Yr)
  } else {
    if (pc_selection[[1]] == "opc" | diss_method == "pls") {
      stop("Yr must be provided when the 'opc' is used in pc_selection is used or diss_method = 'pls'")
    }
  }

  projection <- ortho_projection(
    Xr = Xr,
    Yr = Yr,
    Xu = Xu,
    method = diss_method,
    pc_selection = pc_selection,
    center = center,
    scale = scale, ...
  )

  scores <- projection$scores
  if (diss_method %in% c("pca", "pca.nipals")) {
    scores <- sweep(projection$scores,
      MARGIN = 2,
      STATS = projection$scores_sd,
      FUN = "/"
    )
    use_distance_method <- "euclid"
  } else {
    use_distance_method <- "mahalanobis"
  }

  n_components <- projection$n_components

  if (is.null(Xu) | compute_all) {
    distnc <- f_diss(
      Xr = scores,
      Xu = NULL,
      diss_method = use_distance_method,
      center = FALSE, scale = FALSE
    )
    dimnames(distnc) <- list(
      rownames(scores),
      rownames(scores)
    )
  }

  if (!compute_all & !is.null(Xu)) {
    distnc <- f_diss(
      Xr = scores[1:nrow(Xr), , drop = FALSE],
      Xu = scores[(1 + nrow(Xr)):nrow(scores), , drop = FALSE],
      diss_method = use_distance_method,
      center = FALSE, scale = FALSE
    )
    dimnames(distnc) <- list(
      rownames(scores[1:nrow(Xr), , drop = FALSE]),
      rownames(scores[(1 + nrow(Xr)):nrow(scores), , drop = FALSE])
    )
  }

  if (.local) {
    n_xr <- nrow(Xr)
    if (!is.null(Xu) & compute_all) {
      pre_k <- pre_k + 1
      Xr <- rbind(Xr, Xu)
      Yr <- rbind(Yr, matrix(rep(NA, nrow(Xu))))
    }

    if (is.null(Xu)) {
      pre_k <- pre_k + 1
    }

    neighborhood_info <- k0_indices <- apply(distnc,
      MARGIN = 2,
      FUN = order
    )[1:pre_k, ]
    d_dimnames <- dimnames(distnc)
    rm(distnc)

    neighborhood_info[k0_indices <= n_xr] <- paste0(
      "Xr_",
      neighborhood_info[k0_indices <= n_xr]
    )
    neighborhood_info[k0_indices > n_xr] <- paste0(
      "Xu_",
      neighborhood_info[k0_indices > n_xr]
    )
    neighborhood_info <- data.table(
      do.call(
        "rbind",
        strsplit(colnames(k0_indices), "_")
      ),
      t(neighborhood_info)
    )
    colnames(neighborhood_info) <- c(
      "Set", "Index",
      paste0("k_", 1:nrow(k0_indices))
    )

    if ("opc" %in% pc_selection | diss_method == "pls") {
      # k0_index_matrix <- (apply(distnc[, (1 + nrow(Xr)):ncol(distnc)], 2, order))[1:pre_k, ]
      have_response <- colSums(!apply(k0_indices,
        MARGIN = 2,
        FUN = function(i, y) is.na(y[i, ]),
        y = Yr
      ))

      if (any(c("opc", "manual") %in% pc_selection)) {
        insufficient <- which(have_response < pc_selection[[2]])
        if (length(insufficient) > 0) {
          inms <- names(insufficient)

          mss <- format_xr_xu_indices(xr_xu_names = inms)

          omessage <- paste0(
            "Local disssimilarity matrix cannot be computed as ",
            "the following observation \n",
            "  indices have insufficient neighbors with Yr values:",
            mss$xr_mss,
            mss$Xu_mss
          )
          stop(omessage)
        }
      }

      if (any(have_response < 4)) {
        inms <- names(which(have_response < 3))

        mss <- format_xr_xu_indices(xr_xu_names = inms)

        omessage <- paste0(
          "Local disssimilarity matrix cannot be computed as ",
          "the following observation \n",
          "  indices have less than 3 neighbors with Yr values:",
          mss$xr_mss,
          mss$Xu_mss
        )

        stop(omessage)
      }
      cnt <- 100 * have_response / nrow(k0_indices)

      below_f1 <- sum(cnt != 100)
      if (below_f1 > 0) {
        omessage <- paste0(
          "The neighborhoods of ", below_f1,
          " observations contain missing 'Yr' values."
        )

        below_f2 <- sum(cnt < 50)
        if (below_f2 > 0) {
          omessage2 <- paste0(
            "\nFor ", below_f2, " of these observations ",
            "Yr is missing for more than 50% of their neighbors."
          )
          omessage <- paste0(omessage, omessage2)
        }
        omessage <- paste0(omessage, "\nCheck ...$neighborhood_info")
        message(omessage)
      }
      neighborhood_info <- cbind(
        missing_Yr = nrow(k0_indices) - have_response,
        neighborhood_info
      )
    }

    local_d <- local_ortho_diss(
      k_index_matrix = k0_indices,
      Xr = Xr,
      Yr = Yr,
      Xu = Xu,
      diss_method = diss_method,
      pc_selection = pc_selection,
      center = center,
      scale = scale,
      allow_parallel = allow_parallel, ...
    )

    dimnames(local_d$dissimilarity_m) <- d_dimnames

    neighborhood_info <- cbind(neighborhood_info[, 2:1],
      local_n_components = local_d$local_n_components,
      neighborhood_info[, -c(1:2)]
    )

    resultsList <- list(
      n_components = n_components,
      global_variance_info = projection$variance,
      neighborhood_info = neighborhood_info,
      dissimilarity = local_d$dissimilarity_m
    )

    if (return_projection) {
      resultsList$projection <- projection
    }

    class(resultsList) <- c("ortho_diss", "list")
    class(resultsList$dissimilarity) <- c("local_ortho_diss", "matrix")
    return(resultsList)
  } else {
    resultsList <- list(
      n_components = n_components,
      global_variance_info = projection$variance,
      dissimilarity = distnc
    )
    if (return_projection) {
      resultsList$projection <- projection
    }

    class(resultsList) <- c("ortho_diss", "list")
    class(resultsList$dissimilarity) <- c("ortho_diss", "matrix")
    return(resultsList)
  }
}
