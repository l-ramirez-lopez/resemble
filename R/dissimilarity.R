#' @title Dissimilarity computation between matrices
#' @description
#' \loadmathjax
#'
#' This is a wrapper to integrate the different dissimilarity functions of the
#' offered by package.It computes the dissimilarities between observations in
#' numerical matrices by using an specifed dissmilarity measure.
#' @usage
#' dissimilarity(Xr, Xu = NULL,
#'               diss_method = c("pca", "pca.nipals", "pls", "mpls",
#'                               "cor", "euclid", "cosine", "sid"),
#'               Yr = NULL, gh = FALSE, pc_selection = list("var", 0.01),
#'               return_projection = FALSE, ws = NULL,
#'               center = TRUE, scale = FALSE, documentation = character(),
#'               ...)
#' @param Xr a matrix of containing `n` observations/rows and `p`
#' variables/columns.
#' @param Xu an optional matrix containing data of a second set of observations
#' with `p` variables/columns.
#' @param diss_method a character string indicating the method to be used to
#' compute the dissimilarities between observations. Options are:
#' \itemize{
#'        \item{\code{"pca"}:}{ Mahalanobis distance
#'        computed on the matrix of scores of a Principal Component (PC)
#'        projection of \code{Xr} (and \code{Xu} if provided). PC projection is
#'        done using the singular value decomposition (SVD) algorithm.
#'        See \code{\link{ortho_diss}} function.}
#'
#'        \item{\code{"pca.nipals"}:}{ Mahalanobis distance
#'        computed on the matrix of scores of a Principal Component (PC)
#'        projection of \code{Xr} (and \code{Xu} if provided). PC projection is
#'        done using the non-linear iterative partial least squares (niapls)
#'        algorithm. See \code{\link{ortho_diss}} function.}
#'
#'        \item{\code{"pls"}}:{ Mahalanobis distance
#'        computed on the matrix of scores of a partial least squares projection
#'        of \code{Xr} (and \code{Xu} if provided). In this case, \code{Yr} is
#'        always required. See \code{\link{ortho_diss}} function.}
#'        
#'        \item{\code{"mpls"}}:{ Mahalanobis distance
#'        computed on the matrix of scores of a modified partial least squares 
#'        projection (Shenk and Westerhaus, 1991; Westerhaus, 2014)
#'        of \code{Xr} (and \code{Xu} if provided). In this case, \code{Yr} is
#'        always required. See \code{\link{ortho_diss}} function.}
#'
#'        \item{\code{"cor"}:}{ based on the correlation coefficient
#'        between observations. See \code{\link{cor_diss}} function.}
#'
#'        \item{\code{"euclid"}:}{ Euclidean distance
#'        between observations. See \code{\link{f_diss}} function.}
#'
#'        \item{\code{"cosine"}:}{ Cosine distance
#'        between observations. See \code{\link{f_diss}} function.}
#'
#'        \item{\code{"sid"}:}{ spectral information divergence between
#'        observations. See \code{\link{sid}} function.}
#'        }
#' @param Yr a numeric matrix of `n` observations used as side information of
#' \code{Xr} for the \code{\link{ortho_diss}} methods (i.e. \code{pca},
#' \code{pca.nipals} or \code{pls}). It is required when:
#' \itemize{
#'        \item{\code{diss_method = "pls"}}
#'        \item{\code{diss_method = "pca"} with \code{"opc"} used as the method
#'        in the \code{pc_selection} argument. See \code{\link{ortho_diss}.}}
#'        \item{\code{gh = TRUE}}
#'        }
#' @param gh a logical indicating if the Mahalanobis distance in the pls score
#' space between each observation and the means of the pls scores (center) must
#' be computed. This distance is also known as GH distance. Default \code{FALSE}.
#' @param pc_selection a list of length 2 to be passed onto the
#' \code{\link{ortho_diss}} methods. It is required if the method selected in
#' \code{diss_method} is any of \code{"pca"}, \code{"pca.nipals"} or
#' \code{"pls"} or if \code{gh = TRUE}. This argument is used for
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
#'        case \code{value} must be a value ((larger than 0 and
#'        below the minimum dimension of \code{Xr} or \code{Xr} and \code{Xu}
#'        combined) indicating the maximum
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
#'        (larger than 0 and
#'        below the minimum dimension of \code{Xr} or \code{Xr} and \code{Xu}
#'        combined).
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
#' called (i.e. \code{diss_method = "pca"}, \code{diss_method = "pca.nipals"} or
#' \code{diss_method = "pls"}) or when \code{gh = TRUE}.
#' In case \code{gh = TRUE} and a \code{\link{ortho_diss}} method is used (in the
#' \code{diss_method} argument), both projections are returned.
#' @param ws  an odd integer value which specifies the window size, when
#' \code{diss_method = "cor"} (\code{\link{cor_diss}} method) for moving
#' correlation dissimilarity. If \code{ws = NULL} (default), then the window
#' size will be equal to the number of variables (columns), i.e. instead moving
#' correlation, the normal correlation will be used. See \code{\link{cor_diss}}
#' function.
#' @param center a logical indicating if \code{Xr} (and \code{Xu} if provided)
#' must be centered. If \code{Xu} is provided the data is centered around the
#' mean of the pooled \code{Xr} and \code{Xu} matrices (\mjeqn{Xr \cup Xu}{Xr U Xu}). For
#' dissimilarity computations based on \code{diss_method = pls}, the data is
#' always centered.
#' @param scale a logical indicating if \code{Xr} (and \code{Xu} if
#' provided) must be  scaled. If \code{Xu} is provided the data is scaled based
#' on the standard deviation of the the pooled \code{Xr} and \code{Xu} matrices
#' (\mjeqn{Xr \cup Xu}{Xr U Xu}). If \code{center = TRUE}, scaling is applied after
#' centering.
#' @param gh a logical indicating if the Mahalanobis distance (in the pls score
#' space) between each observation and the pls centre/mean must be
#' computed.
#' @param documentation an optional character string that can be used to
#' describe anything related to the \code{mbl} call (e.g. description of the
#' input data). Default: \code{character()}. NOTE: his is an experimental
#' argument.
#' @param ... other arguments passed to the dissimilarity functions
#' (\code{\link{ortho_diss}}, \code{\link{cor_diss}}, \code{\link{f_diss}} or
#' \code{\link{sid}}).
#'
#' @details
#' This function is a wrapper for \code{\link{ortho_diss}}, \code{\link{cor_diss}},
#'  \code{\link{f_diss}}, \code{\link{sid}}. Check the documentation of these
#'  functions for further details.
#'
#' @seealso \code{\link{ortho_diss}} \code{\link{cor_diss}} \code{\link{f_diss}}
#' \code{\link{sid}}.
#'
#' @return A list with the following components:
#' \itemize{
#'        \item{\code{dissimilarity}:}{ the resulting dissimilarity matrix.}
#'
#'        \item{\code{projection}:}{ an \code{ortho_projection} object. Only output
#'        if \code{return_projection = TRUE} and if \code{diss_method = "pca"},
#'        \code{diss_method = "pca.nipals"},  \code{diss_method = "pls"} or  
#'        \code{diss_method = "mpls"}.
#'
#'        This object contains the projection used to compute
#'        the dissimilarity matrix. In case of local dissimilarity matrices,
#'        the projection corresponds to the global projection used to select the
#'        neighborhoods (see \code{\link{ortho_diss}} function for further
#'        details).}
#'
#'        \item{\code{gh}:}{ a list containing the GH distances as well as the
#'        pls projection used to compute the GH.}
#'        }
#' @references
#' Shenk, J., Westerhaus, M., and Berzaghi, P. 1997. Investigation of a LOCAL
#' calibration procedure for near infrared instruments. Journal of Near Infrared
#' Spectroscopy, 5, 223-232.
#'
#' Westerhaus, M. 2014. Eastern Analytical Symposium Award for outstanding 
#' Wachievements in near infrared spectroscopy: my contributions to 
#' Wnear infrared spectroscopy. NIR news, 25(8), 16-20.
#' 
#' @author \href{https://orcid.org/0000-0002-5369-5120}{Leonardo Ramirez-Lopez}
#' @examples
#' library(prospectr)
#' data(NIRsoil)
#'
#' # Filter the data using the first derivative with Savitzky and Golay
#' # smoothing filter and a window size of 11 spectral variables and a
#' # polynomial order of 4
#' sg <- savitzkyGolay(NIRsoil$spc, m = 1, p = 4, w = 15)
#'
#' # Replace the original spectra with the filtered ones
#' NIRsoil$spc <- sg
#'
#' Xu <- NIRsoil$spc[!as.logical(NIRsoil$train), ]
#' Yu <- NIRsoil$CEC[!as.logical(NIRsoil$train)]
#'
#' Yr <- NIRsoil$CEC[as.logical(NIRsoil$train)]
#' Xr <- NIRsoil$spc[as.logical(NIRsoil$train), ]
#'
#' Xu <- Xu[!is.na(Yu), ]
#' Xr <- Xr[!is.na(Yr), ]
#'
#' Yu <- Yu[!is.na(Yu)]
#' Yr <- Yr[!is.na(Yr)]
#'
#' dsm_pca <- dissimilarity(
#'   Xr = Xr, Xu = Xu,
#'   diss_method = c("pca"),
#'   Yr = Yr, gh = TRUE,
#'   pc_selection = list("opc", 30),
#'   return_projection = TRUE
#' )
#' @export
#'
## History:
## 22.05.2020 Leo     Hello world!


dissimilarity <- function(Xr,
                          Xu = NULL,
                          diss_method = c(
                            "pca",
                            "pca.nipals",
                            "pls",
                            "mpls",
                            "cor",
                            "euclid",
                            "cosine",
                            "sid"
                          ),
                          Yr = NULL,
                          gh = FALSE,
                          pc_selection = list("var", 0.01),
                          return_projection = FALSE,
                          ws = NULL,
                          center = TRUE,
                          scale = FALSE,
                          documentation = character(),
                          ...) {

  ## Future arguments/features?
  ## - group function to be passed to the opc methods?

  result <- list(dissimilarity = NULL)
  # Mahalanobis is excluded from this list because when used on matrices with 
  # highly correlated variables, it returns singular covariance matrices. So it
  # need to be prevented as it doe snot really make sense to compute mahalanobis 
  # on the raw spectra
  avalmethods <- c(
    "pca",
    "pca.nipals",
    "cor",
    "movcor",
    "pls",
    "mpls",
    "euclid",
    "cosine",
    "sid"
  )

  # if(!is.null(group))
  # {
  #   if(length(group) != nrow(Xr))
  #     stop("The length of 'group' must be equal to the number of observations in 'Xr'")
  # }

  if (length(diss_method) > 1) {
    warning("'diss_method' has length > 1 and only the first element will be used")
    diss_method <- diss_method[1]
  }

  if (is.null(Yr) & gh) {
    stop("for gh is necessary to supply Yr")
  }


  if (!diss_method %in% avalmethods) {
    stop(paste(
      "'diss_method' argument needs to be equal to one of the following options:\n",
      paste("'", avalmethods, "'", collapse = ", ", sep = "")
    ))
  }


  ortho_dissmethods <- c("pca", "pca.nipals", "pls", "mpls")
  f_dissmethods <- c("euclid", "cosine")
  cor_dissmethods <- c("cor")
  divergencemethods <- c("sid")

  if (diss_method %in% ortho_dissmethods) {
    if (diss_method == "mpls") {
      modified <- TRUE
    } else {
      modified <- FALSE
    }
    pcDistance <- ortho_diss(
      Xr = Xr,
      Xu = Xu,
      pc_selection = pc_selection,
      Yr = Yr,
      diss_method = diss_method,
      center = center,
      scale = scale,
      modified = modified,
      return_projection = ((gh & diss_method == "pls") | return_projection),
      ...
    )

    dmat <- pcDistance$dissimilarity

    if (return_projection) {
      result$projection <- pcDistance$projection
    }
  }

  if (gh) {
    if (diss_method == "pls") {
      pca <- pcDistance$projection
      rm(pcDistance)
    } else {
      pca <- pls_projection(
        Xr = rbind(Xr, Xu),
        Xu = NULL,
        Yr = c(
          Yr,
          rep(NA, ifelse(is.null(Xu), 0, nrow(Xu)))
        ),
        pc_selection = pc_selection,
        scale = scale,
        ...
      )
      if (!is.null(Xu)) {
        rownames(pca$scores)[-(1:nrow(Xr))] <- paste0("Xu_", 1:nrow(Xu))
      }
    }

    scores <- pca$scores

    scoresmean <- t(colMeans(scores[1:nrow(Xr), , drop = FALSE]))

    gh <- as.vector(f_diss(
      Xr = scores,
      Xu = scoresmean,
      center = FALSE,
      scale = FALSE,
      diss_method = "mahalanobis"
    ))

    result$gh <- list(gh_Xr = gh[1:nrow(Xr)])
    if (!is.null(Xu)) {
      result$gh$gh_Xu <- gh[(1 + nrow(Xr)):nrow(scores)]
    }
    result$gh$projection <- pca
    rm(scores)
  }

  if (diss_method %in% f_dissmethods) {
    dmat <- f_diss(
      Xr = Xr,
      Xu = Xu,
      center = center,
      scale = scale,
      diss_method = diss_method
    )
  }

  if (diss_method %in% cor_dissmethods) {
    dmat <- cor_diss(
      Xr = Xr,
      Xu = Xu,
      ws = ws,
      center = center,
      scale = scale
    )
    result$ws <- ws
  }

  if (diss_method %in% divergencemethods) {
    dmat <- sid(
      Xr = Xr, Xu = Xu,
      center = center,
      scale = scale,
      reg = 10^-4, ...
    )$sid
  }
  result$dissimilarity <- dmat
  result$documentation <- documentation

  result
}
#
