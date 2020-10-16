#' @title A function for searching in a given reference set the neighbors of
#' another given set of observations (search_neighbors)
#' @description
#' \loadmathjax
#' This function searches in a reference set the neighbors of the observations
#' provided  in another set.
#' @usage
#' search_neighbors(Xr, Xu, diss_method = c("pca", "pca.nipals", "pls", "cor",
#'                                          "euclid", "cosine", "sid"),
#'                  Yr = NULL, k, k_diss, k_range, spike = NULL,
#'                  pc_selection = list("var", 0.01),
#'                  return_projection = FALSE, return_dissimilarity = FALSE,
#'                  ws = NULL,
#'                  center = TRUE, scale = FALSE,
#'                  documentation = character(), ...)
#'
#' @param Xr a matrix of reference (spectral) observations where the neighbors
#' of the observations in \code{Xu} are to be searched.
#' @param Xu a matrix of (spectral) observations for which its neighbors are to
#' be searched in \code{Xu}.
#' @param diss_method a character string indicating the spectral dissimilarity metric
#' to be used in the selection of the nearest neighbors of each observation.
#' \itemize{
#'        \item{\code{"pca"}:}{  Mahalanobis distance
#'        computed on the matrix of scores of a Principal Component (PC)
#'        projection of \code{Xr} and \code{Xu}. PC projection is done using the
#'        singlar value decomposition (SVD) algorithm.
#'        See \code{\link{ortho_diss}} function.}
#'
#'        \item{\code{"pca.nipals"}}{ Mahalanobis distance
#'        computed on the matrix of scores of a Principal Component (PC)
#'        projection of \code{Xr} and \code{Xu}. PC projection is done using the
#'        non-linear iterative partial least squares (niapls) algorithm.
#'        See \code{\link{ortho_diss}} function.}
#'
#'        \item{\code{"pls"}}{ Mahalanobis distance
#'        computed on the matrix of scores of a partial least squares projection
#'        of \code{Xr} and \code{Xu}. In this case, \code{Yr} is always required.
#'        See \code{\link{ortho_diss}} function.}
#'
#'        \item{\code{"cor"}}{ correlation coefficient
#'        between observations. See \code{\link{cor_diss}} function.}
#'
#'        \item{\code{"euclid"}}{ Euclidean distance
#'        between observations. See \code{\link{f_diss}} function.}
#'
#'        \item{\code{"cosine"}}{ Cosine distance
#'        between observations. See \code{\link{f_diss}} function.}
#'
#'        \item{\code{"sid"}}{ spectral information divergence between observations.
#'        See \code{\link{sid}} function.}
#'        }
#' @param Yr a numeric matrix of `n` observations used as side information of
#' \code{Xr} for the \code{\link{ortho_diss}} methods (i.e. \code{pca},
#' \code{pca.nipals} or \code{pls}). It is required when:
#' \itemize{
#'        \item{\code{diss_method = "pls"}}
#'        \item{\code{diss_method = "pca"} with \code{"opc"} used as the method
#'        in the \code{pc_selection} argument. See [ortho_diss()].}
#'        }
#' @param k an integer value indicating the k-nearest neighbors of each
#' observation in \code{Xu} that must be selected from \code{Xr}.
#' @param k_diss an integer value indicating a dissimilarity treshold.
#' For each observation in \code{Xu}, its nearest neighbors in \code{Xr}
#' are selected as those for which their dissimilarity to \code{Xu} is below
#' this \code{k_diss} threshold. This treshold depends on the corresponding
#' dissimilarity metric specified in \code{diss_method}. Either \code{k} or
#' \code{k_diss} must be specified.
#' @param k_range an integer vector of length 2 which specifies the minimum
#' (first value) and the maximum (second value) number of neighbors to be
#' retained when the \code{k_diss} is given.
#' @param spike a vector of integers indicating what observations in \code{Xr}
#' (and \code{Yr}) must be 'forced' to always be part of all the neighborhoods.
#' @param pc_selection a list of length 2 to be passed onto the
#' \code{\link{ortho_diss}} methods. It is required if the method selected in
#' \code{diss_method} is any of \code{"pca"}, \code{"pca.nipals"} or
#' \code{"pls"}. This argument is used for
#' optimizing the number of components (principal components or pls factors)
#' to be retained. This list must contain two elements in the following order:
#' \code{method} (a character indicating the method for selecting the number of
#' components) and \code{value} (a numerical value that complements the selected
#' method). The methods available are:
#' \itemize{
#'        \item{\code{"opc"}:} { optimized principal component selection based on
#'        Ramirez-Lopez et al. (2013a, 2013b). The optimal number of components
#'        (of set of observations) is the one for which its distance matrix
#'        minimizes the differences between the \code{Yr} value of each
#'        observation and the \code{Yr} value of its closest observation. In this
#'        case \code{value} must be a value  (larger than 0 and below the 
#'        minimum dimension of \code{Xr} or \code{Xr} and \code{Xu} combined) 
#'        indicating the maximum number of principal components to be tested. 
#'        See the \code{\link{ortho_projection}} function for more details.}
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
#'        (larger than 0 and below the 
#'        minimum dimension of \code{Xr} or \code{Xr} and \code{Xu} combined) 
#'        indicating the minimum amount of variance that a component should
#'        explain in order to be retained.}
#'        }
#' The default is \code{list(method = "var", value = 0.01)}.
#' 
#' Optionally, the \code{pc_selection} argument admits \code{"opc"} or
#' \code{"cumvar"} or \code{"var"} or \code{"manual"} as a single character
#' string. In such a case the default \code{"value"} when either \code{"opc"} or
#' \code{"manual"} are used is 40. When \code{"cumvar"} is used the default
#' \code{"value"} is set to 0.99 and when \code{"var"} is used, the default
#' \code{"value"} is set to 0.01.
#' @param return_projection a logical indicating if the projection(s) must be
#' returned. Projections are used if the \code{\link{ortho_diss}} methods are
#' called (i.e. \code{method = "pca"}, \code{method = "pca.nipals"} or
#' \code{method = "pls"}).
#' @param return_dissimilarity a logical indicating if the dissimilarity matrix
#' used for neighbor search must be returned.
#' @param ws  an odd integer value which specifies the window size, when
#' \code{diss_method = cor} (\code{\link{cor_diss}} method) for moving correlation
#' dissimilarity. If \code{ws = NULL} (default), then the window size will be
#' equal to the number of variables (columns), i.e. instead moving correlation,
#' the normal correlation will be used. See \code{\link{cor_diss}} function.
#' @param center a logical indicating if the \code{Xr} and \code{Xu} matrices
#'  must be centered. If \code{Xu} is provided the data is centered around the
#'  mean of the pooled \code{Xr} and \code{Xu} matrices (\mjeqn{Xr \cup Xu}{Xr U Xu}). For
#' dissimilarity computations based on \code{diss_method = pls}, the data is always
#' centered.
#' @param scale a logical indicating if the \code{Xr} and \code{Xu} matrices
#' must be scaled. If \code{Xu} is provided the data is scaled based
#' on the standard deviation of the the pooled \code{Xr} and \code{Xu} matrices
#' (\mjeqn{Xr \cup Xu}{Xr U Xu}). If \code{center = TRUE}, scaling is applied after
#' centering.
#' @param documentation an optional character string that can be used to
#' describe anything related to the \code{mbl} call (e.g. description of the
#' input data). Default: \code{character()}. NOTE: his is an experimental
#' argument.
#' @param ... further arguments to be passed to the \code{\link{dissimilarity}}
#' fucntion. See details.
#' @details
#' This function may be specially useful when the reference set (\code{Xr}) is
#' very large. In some cases the number of observations in the reference set
#' can be reduced by removing irrelevant observations (i.e. observations that are not
#' neighbors of a particular target set). For example, this fucntion can be
#' used to reduce the size of the reference set before before  running the
#' \code{\link{mbl}} function.
#'
#' This function uses the \code{\link{dissimilarity}} fucntion to compute the
#' dissimilarities between \code{Xr} and \code{Xu}. Arguments to
#' \code{\link{dissimilarity}} as well as further arguments to the functions
#' used inside \code{\link{dissimilarity}} (i.e. \code{\link{ortho_diss}}
#' \code{\link{cor_diss}} \code{\link{f_diss}} \code{\link{sid}}) can be passed to
#' those functions as additional arguments (i.e. \code{...}).
#' @return a \code{list} containing the following elements:
#' \itemize{
#'  \item{\code{neighbors_diss}}{ a matrix of the \code{Xr} dissimilarity socres
#'  corresponding to the neighbors of each observation in \code{Xu}.
#'  The neighbor dissimilarity socres are organized by columns and are sorted
#'  in ascending order.}
#'  \item{\code{neighbors}}{ a matrix of the \code{Xr} indices corresponding to
#'  the neighbors of each observation in \code{Xu}. The neighbor indices are
#'  organized by columns and are sorted in ascending order by their
#'  dissimilarity score.}
#'  \item{\code{unique_neighbors}}{ a vector of the indices in \code{Xr}
#'  identified as neighbors of any observation in \code{Xu}. This is obtained by
#'  converting the \code{neighbors} matrix into a vector and applying the
#'  \code{\link[base]{unique}} function.}
#'  \item{\code{k_diss_info}}{ a \code{data.table} that is returned only if the
#'  \code{k_diss} argument was used. It comprises three columns, the first one
#'  (\code{Xu_index}) indicates the index of the observations in \code{Xu},
#'  the second column (\code{n_k}) indicates the number of neighbors found in
#'  \code{Xr} and the third column (\code{final_n_k}) indicates the final number
#'  of neighbors selected bounded by \code{k_range}.
#'  argument.}
#'  \item{\code{dissimilarity}}{ If \code{return_dissimilarity = TRUE} the
#'  dissimilarity object used (as computed by the \code{\link{dissimilarity}}
#'  function.}
#'  \item{\code{projection}}{ an \code{ortho_projection} object. Only output if
#'        \code{return_projection = TRUE} and if \code{diss_method = "pca"},
#'        \code{diss_method = "pca.nipals"} or \code{diss_method = "pls"}.
#'        
#'        This object contains the projection used to compute
#'        the dissimilarity matrix. In case of local dissimilarity matrices,
#'        the projection corresponds to the global projection used to select the
#'        neighborhoods.  (see \code{\link{ortho_diss}} function for further
#'        details).}
#'  }
#' @author \href{https://orcid.org/0000-0002-5369-5120}{Leonardo Ramirez-Lopez}
#' @references
#' Ramirez-Lopez, L., Behrens, T., Schmidt, K., Stevens, A., Dematte, J.A.M.,
#' Scholten, T. 2013a. The spectrum-based learner: A new local approach for
#' modeling soil vis-NIR spectra of complex data sets. Geoderma 195-196, 268-279.
#'
#' Ramirez-Lopez, L., Behrens, T., Schmidt, K., Viscarra Rossel, R.,
#' Dematte, J. A. M.,  Scholten, T. 2013b. Distance and similarity-search
#' metrics for use with soil vis-NIR spectra. Geoderma 199, 43-53.
#' @seealso \code{\link{dissimilarity}} \code{\link{ortho_diss}}
#' \code{\link{cor_diss}} \code{\link{f_diss}} \code{\link{sid}}
#'  \code{\link{mbl}}
#' @examples
#' \dontrun{
#' library(prospectr)
#'
#' data(NIRsoil)
#'
#' Xu <- NIRsoil$spc[!as.logical(NIRsoil$train), ]
#' Yu <- NIRsoil$CEC[!as.logical(NIRsoil$train)]
#' Yr <- NIRsoil$CEC[as.logical(NIRsoil$train)]
#' Xr <- NIRsoil$spc[as.logical(NIRsoil$train), ]
#'
#' Xu <- Xu[!is.na(Yu), ]
#' Yu <- Yu[!is.na(Yu)]
#'
#' Xr <- Xr[!is.na(Yr), ]
#' Yr <- Yr[!is.na(Yr)]
#'
#' # Identify the neighbor observations using the correlation dissimilarity and
#' # default parameters
#' # (In this example all the observations in Xr belong at least to the
#' # first 100 neighbors of one observation in Xu)
#' ex1 <- search_neighbors(
#'   Xr = Xr, Xu = Xu,
#'   diss_method = "cor",
#'   k = 40
#' )
#'
#' # Identify the neighbor observations using principal component (PC)
#' # and partial least squares (PLS) dissimilarities, and using the "opc"
#' # approach for selecting the number of components
#' ex2 <- search_neighbors(
#'   Xr = Xr, Xu = Xu,
#'   diss_method = "pca",
#'   Yr = Yr, k = 50,
#'   pc_selection = list("opc", 40),
#'   scale = TRUE
#' )
#'
#' # Observations that do not belong to any neighborhood
#' seq(1, nrow(Xr))[!seq(1, nrow(Xr)) %in% ex2$unique_neighbors]
#'
#' ex3 <- search_neighbors(
#'   Xr = Xr, Xu = Xu,
#'   diss_method = "pls",
#'   Yr = Yr, k = 50,
#'   pc_selection = list("opc", 40),
#'   scale = TRUE
#' )
#' # Observations that do not belong to any neighborhood
#' seq(1, nrow(Xr))[!seq(1, nrow(Xr)) %in% ex3$unique_neighbors]
#'
#' # Identify the neighbor observations using local PC dissimialrities
#' # Here, 150 neighbors are used to compute a local dissimilarity matrix
#' # and then this matrix is used to select 50 neighbors
#' ex4 <- search_neighbors(
#'   Xr = Xr, Xu = Xu,
#'   diss_method = "pls",
#'   Yr = Yr, k = 50,
#'   pc_selection = list("opc", 40),
#'   scale = TRUE,
#'   .local = TRUE,
#'   pre_k = 150
#' )
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
## 25.05.2020 Leo     - function renamed from neigCleaning to search_neighbors
##                    - scaled renamed to scale
##                    - pcMethod and cores are deprecated

search_neighbors <- function(Xr, Xu, diss_method = c(
                               "pca",
                               "pca.nipals",
                               "pls",
                               "cor",
                               "euclid",
                               "cosine",
                               "sid"
                             ),
                             Yr = NULL,
                             k, k_diss, k_range,
                             spike = NULL,
                             pc_selection = list("var", 0.01),
                             return_projection = FALSE,
                             return_dissimilarity = FALSE,
                             ws = NULL,
                             center = TRUE, scale = FALSE,
                             documentation = character(), ...) {


  # Sanity checks
  match.arg(diss_method, c(
    "pca",
    "pca.nipals",
    "pls",
    "cor",
    "euclid",
    "cosine",
    "sid"
  ))

  if (missing(k)) {
    k <- NULL
  }

  if (missing(k_diss)) {
    k_diss <- NULL
  }

  if (missing(k_range)) {
    k_range <- NULL
  }

  if (!is.logical(center)) {
    stop("'center' argument must be logical")
  }

  if (!is.logical(scale)) {
    stop("'scale' argument must be logical")
  }

  if (diss_method == "cor") {
    if (!is.null(ws)) {
      if (ws < 3 | ws > (ncol(Xr) - 1) | length(ws) != 1 | (ws %% 2) == 0) {
        stop(paste(
          "ws must be an odd value between 3 and", ncol(Xr) - 1,
          "(ncol(Xr) - 1)"
        ))
      }
    }
  }

  if (!is.null(k) & !is.null(k_diss)) {
    # if k and k_diss are not called here, errors are thrown during checks
    k
    k_diss
    stop("Only one of k or k_diss can be specified")
  }

  if (is.null(k) & is.null(k_diss)) {
    stop("Either k or k_diss must be specified")
  }

  if (!is.null(k)) {
    k <- as.integer(k)
    if (k < 1) {
      stop("k must be an integer larger than 0")
    }
    if (k > nrow(Xr)) {
      stop(paste0(
        "The number of nearest neighbors cannot exceed the number ",
        "of reference observations (nrow(Xr))"
      ))
    }
    kk <- k
  }

  if (!is.null(k_diss)) {
    # if k_diss is not called here, errors are thrown during checks
    k_diss
    if (is.null(k_range)) {
      # if k_range is not called here, errors are thrown during checks
      k_range
      stop("if the k_diss argument is used, k_range must be specified")
    }
    if (length(k_range) != 2 | !is.numeric(k_range) | diff(k_range) < 0) {
      stop(paste0(
        "k_range must be a vector (of length 2) which specifies the ",
        "minimum (first value) and the maximum (second value) number ",
        "of neighbors"
      ))
    }
    k_min <- as.integer(k_range[1])
    k_max <- as.integer(k_range[2])
    if (k_min < 1) {
      stop("Minimum number of nearest neighbors must be larger than 0")
    }
    if (k_max > nrow(Xr)) {
      stop(paste0(
        "Maximum number of nearest neighbors cannot exceed the ",
        "number of reference observations"
      ))
    }
    kk <- k_max
  }
  input_dots <- list(...)
  if (".local" %in% names(input_dots)) {
    if (isTRUE(input_dots$.local)) {
      if (!"pre_k" %in% names(input_dots)) {
        stop(paste0(
          "When .local = TRUE (passed to the ortho_diss method), the ",
          "'pre_k' argument must be specified"
        ))
      }
      if (input_dots$pre_k < kk) {
        stop(paste0(
          "pre_k must be larger than ",
          if_else(is.null(k), "max(k_range)", "k")
        ))
      }
    }
  }

  if (!is.null(spike)) {
    if (!is.vector(spike)) {
      stop("spike must be a vector of integers")
    }
    if (any(!spike %% 1 == 0)) {
      stop("spike must be a vector of integers")
    }
    if (length(spike) >= nrow(Xr)) {
      stop("The lebgth of spike cannot be larger or equal to the number of rows of Xr")
    }
    if (max(spike) > nrow(Xr)) {
      stop("Argument spike contains indices subscript out of bounds of Xr")
    }
    if (!is.null(k)) {
      if (min(k) <= length(spike)) {
        stop("values for k must be larger than length(spike)")
      }
    }
    if (!is.null(k_diss)) {
      if (min(k_range) <= length(spike)) {
        stop("values for k_range must be larger than length(spike)")
      }
    }
    spike <- sort(unique(as.integer(spike)))
  }
  dsm <- dissimilarity(
    Xr = Xr,
    Xu = Xu,
    diss_method = diss_method,
    Yr = Yr,
    pc_selection = pc_selection,
    return_projection = return_projection,
    ws = ws,
    center = center,
    scale = scale,
    ...
  )

  results <- diss_to_neighbors(dsm$dissimilarity,
    k = k, k_diss = k_diss, k_range = k_range,
    spike = spike,
    return_dissimilarity = return_dissimilarity
  )

  if (return_projection & diss_method %in% c("pca", "pca.nipals", "pls")) {
    results$projection <- dsm$projection
  }
  if ("gh" %in% names(input_dots)) {
    results$gh <- dsm$gh
  }

  results
}
