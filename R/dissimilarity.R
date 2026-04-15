#' @title Compute dissimilarity matrices
#'
#' @description
#' \loadmathjax
#' Computes dissimilarity matrices between observations using various methods.
#' This is the main interface for dissimilarity computation in the resemble
#' package.
#'
#' @param Xr A numeric matrix of reference observations (rows) and variables
#'   (columns).
#' @param Xu Optional matrix of additional observations with the same variables.
#' @param diss_method A dissimilarity method object created by one of:
#'   \itemize{
#'     \item \code{\link{diss_pca}()}: Mahalanobis distance in PCA space
#'     \item \code{\link{diss_pls}()}: Mahalanobis distance in PLS space
#'     \item \code{\link{diss_correlation}()}: Correlation-based dissimilarity
#'     \item \code{\link{diss_euclidean}()}: Euclidean distance
#'     \item \code{\link{diss_mahalanobis}()}: Mahalanobis distance
#'     \item \code{\link{diss_cosine}()}: Cosine dissimilarity
#'   }
#'   Default is \code{diss_pca()}.
#' @param Yr Optional response matrix. Required for PLS methods and when using
#'   \code{ncomp_by_opc()}.
#'
#' @details
#' The function dispatches to the appropriate internal computation based on the
#' class of \code{diss_method}. Each method constructor (e.g., \code{diss_pca()})
#' encapsulates all method-specific parameters including component selection,
#' centering, scaling, and whether to return projections.
#'
#' \subsection{Output dimensions}{
#' When only \code{Xr} is provided, the function computes pairwise dissimilarities
#' among all observations in \code{Xr}, returning a symmetric
#' \code{nrow(Xr)} \mjeqn{\times}{x} \code{nrow(Xr)} matrix.
#'
#' When both \code{Xr} and \code{Xu} are provided, the function computes
#' dissimilarities between each observation in \code{Xr} and each observation
#' in \code{Xu}, returning a \code{nrow(Xr)} \mjeqn{\times}{x} \code{nrow(Xu)}
#' matrix where element \mjeqn{(i, j)}{(i, j)} is the dissimilarity between the
#' \mjeqn{i}{i}-th observation in \code{Xr} and the \mjeqn{j}{j}-th observation
#' in \code{Xu}.
#' }
#'
#' \subsection{Mahalanobis distance}{
#' Note that \code{diss_mahalanobis()} computes Mahalanobis distance directly on
#' the input variables. This requires the covariance matrix to be invertible,
#' which fails when the number of variables exceeds the number of observations
#' or when variables are highly correlated (common in spectral data). For such
#' cases, use \code{diss_pca()} or \code{diss_pls()} instead.
#' }
#'
#' @return A list of class \code{"dissimilarity"} containing:
#' \describe{
#'   \item{dissimilarity}{The computed dissimilarity matrix. Dimensions are
#'     \code{nrow(Xr)} \mjeqn{\times}{x} \code{nrow(Xr)} when \code{Xu = NULL},
#'     or \code{nrow(Xr)} \mjeqn{\times}{x} \code{nrow(Xu)} otherwise.}
#'   \item{diss_method}{The \code{diss_*} constructor object used for computation.}
#'   \item{center}{Vector used to center the data.}
#'   \item{scale}{Vector used to scale the data.}
#'   \item{ncomp}{Number of components used (for projection methods).}
#'   \item{projection}{If \code{return_projection = TRUE} in the method
#'     constructor, the \code{ortho_projection} object.}
#' }
#' 
#' @author \href{https://orcid.org/0000-0002-5369-5120}{Leonardo Ramirez-Lopez} 
#' 
#' @seealso
#' \code{\link{diss_pca}}, \code{\link{diss_pls}}, 
#' \code{\link{diss_correlation}}, \code{\link{diss_euclidean}},
#' \code{\link{diss_mahalanobis}}, \code{\link{diss_cosine}}
#' 
#'
#' @references
#' Ramirez-Lopez, L., Behrens, T., Schmidt, K., Stevens, A., Dematte, J.A.M.,
#' Scholten, T. 2013a. The spectrum-based learner: A new local approach for
#' modeling soil vis-NIR spectra of complex data sets. Geoderma 195-196,
#' 268-279.
#'
#' Ramirez-Lopez, L., Behrens, T., Schmidt, K., Viscarra Rossel, R., Dematte,
#' J.A.M., Scholten, T. 2013b. Distance and similarity-search metrics for use
#' with soil vis-NIR spectra. Geoderma 199, 43-53.
#'
#'
#' @examples
#' \donttest{
#' library(prospectr)
#' data(NIRsoil)
#'
#' # Preprocess
#' sg <- savitzkyGolay(NIRsoil$spc, m = 1, p = 4, w = 15)
#'
#' Xr <- sg[as.logical(NIRsoil$train), ]
#' Xu <- sg[!as.logical(NIRsoil$train), ]
#' Yr <- NIRsoil$CEC[as.logical(NIRsoil$train)]
#' Yu <- NIRsoil$CEC[!as.logical(NIRsoil$train)]
#'
#' Xu <- Xu[!is.na(Yu), ]
#' Xr <- Xr[!is.na(Yr), ]
#' Yr <- Yr[!is.na(Yr)]
#'
#' # PCA-based dissimilarity with variance-based selection
#' d1 <- dissimilarity(Xr, Xu, diss_method = diss_pca())
#'
#' # PCA with OPC selection (requires Yr)
#' d2 <- dissimilarity(Xr, Xu,
#'   Yr = Yr,
#'   diss_method = diss_pca(
#'     ncomp = ncomp_by_opc(30),
#'     return_projection = TRUE
#'   )
#' )
#'
#' # PLS-based dissimilarity 
#' d3 <- dissimilarity(
#'   Xr, Xu,
#'   Yr = Yr,
#'   diss_method = diss_pls(
#'     ncomp = ncomp_by_opc(30)
#'   )
#' )
#'
#' # Euclidean distance
#' d4 <- dissimilarity(Xr, Xu, diss_method = diss_euclidean())
#'
#' # Correlation dissimilarity with moving window
#' d5 <- dissimilarity(Xr, Xu, diss_method = diss_correlation(ws = 41))
#'
#' # Mahalanobis distance (use only when n > p and low collinearity)
#' # d6 <- dissimilarity(Xr[, 1:20], Xu[, 1:20],
#' #                     diss_method = diss_mahalanobis())
#' }
#'
#' @export
dissimilarity <- function(
  Xr,
  Xu = NULL,
  diss_method = diss_pca(),
  Yr = NULL
) {
  # ---------------------------------------------------------------------------
  # Handle legacy character-based diss_method
  # ---------------------------------------------------------------------------
  if (is.character(diss_method)) {
    stop(
      "Character-based 'diss_method' is no longer supported.\n\n",
      "Use method constructors instead:\n\n",
      "  Old API                          -> New API\n",
      "  -------                             -------\n",
      "  diss_method = \"pca\"               -> diss_pca()\n",
      "  diss_method = \"pca.nipals\"        -> diss_pca(method = \"pca_nipals\")\n",
      "  diss_method = \"pls\"               -> diss_pls()\n",
      "  diss_method = \"mpls\"              -> diss_pls(method = \"mpls\")\n",
      "  diss_method = \"euclid\"            -> diss_euclidean()\n",
      "  diss_method = \"cosine\"            -> diss_cosine()\n",
      "  diss_method = \"cor\"               -> diss_correlation()\n\n",
      "  diss_method = \"sid\"               -> Deprecated\n",
      "  pc_selection = list(\"opc\", 30)    -> ncomp = ncomp_by_opc(30)\n",
      "  gh = TRUE                           -> Deprecated\n",
      call. = FALSE
    )
  }
  
  # ---------------------------------------------------------------------------
  # check Xr and Xu
  # ---------------------------------------------------------------------------
  if (is.numeric(Xr) && is.null(dim(Xr))) {
    Xr <- matrix(Xr, nrow = 1)
  }
  
  if (!is.null(Xu)) {
    if (is.numeric(Xu) && is.null(dim(Xu))) {
      Xu <- matrix(Xu, nrow = 1)
    }
    if (ncol(Xr) != ncol(Xu)) {
      stop("'Xr' and 'Xu' must have the same number of columns.", call. = FALSE)
    }
  } else {
    if (nrow(Xr) < 2L) {
      stop("'Xr' must have at least 2 rows when 'Xu' is not provided.", call. = FALSE)
    }
  }
  
  # ---------------------------------------------------------------------------
  # Validate diss_method
  # ---------------------------------------------------------------------------
  if (!inherits(diss_method, "diss_method")) {
    stop(
      "'diss_method' must be a dissimilarity method object.\n",
      "Use one of: diss_pca(), diss_pls(), diss_correlation(),",  
      "diss_euclidean(), diss_mahalanobis(), diss_cosine()."
    )
  }

  # ---------------------------------------------------------------------------
  # Validate Yr requirements
  # ---------------------------------------------------------------------------
  method_class <- class(diss_method)[[1]]

  requires_yr <- method_class == "diss_pls" ||
    (method_class == "diss_pca" && inherits(diss_method$ncomp, "ncomp_by_opc"))

  if (requires_yr && is.null(Yr)) {
    stop("'Yr' is required for this dissimilarity method.")
  }

  # ---------------------------------------------------------------------------
  # Dispatch to compute function
  # ---------------------------------------------------------------------------
  result <- switch(
    method_class,
    diss_pca = .diss_pca_compute(Xr, Xu, Yr, diss_method),
    diss_pls = .diss_pls_compute(Xr, Xu, Yr, diss_method),
    diss_euclidean = .diss_euclidean_compute(Xr, Xu, diss_method),
    diss_mahalanobis = .diss_mahalanobis_compute(Xr, Xu, diss_method),
    diss_cosine = .diss_cosine_compute(Xr, Xu, diss_method),
    diss_correlation = .diss_correlation_compute(Xr, Xu, diss_method),
    stop("Unknown dissimilarity method: ", method_class)
  )
  result$diss_method <- diss_method
  
  if (diss_method$center) {
    result$center <- drop(get_column_means(rbind(Xr, Xu)))
  } else {
    result$center <- NULL
  }
  
  if (diss_method$scale) {
    result$scale <- get_col_sds(rbind(Xr, Xu))
  } else {
    result$scale <- NULL
  }  

  class(result) <- c("resemble_dissimilarity", "dissimilarity", "list")
  result
}


# =============================================================================
# Internal compute functions
# =============================================================================

.diss_pca_compute <- function(Xr, Xu, Yr, method) {
  proj <- ortho_projection(
    Xr = Xr,
    Xu = Xu,
    Yr = Yr,
    ncomp = method$ncomp,
    method = method$method,
    center = method$center,
    scale = method$scale
  )

  # Standardize scores by SD before computing Euclidean distance
  scores <- proj$scores
  scores_sd <- apply(scores, 2L, sd)
  scores_sd[scores_sd < .Machine$double.eps] <- 1
  scores_std <- sweep(scores, 2L, scores_sd, "/")

  # Split scores
  n_xr <- nrow(Xr)
  scores_xr <- scores_std[seq_len(n_xr), , drop = FALSE]
  scores_xu <- if (!is.null(Xu)) {
    scores_std[(n_xr + 1L):nrow(scores_std), , drop = FALSE]
  } else {
    NULL
  }

  # Euclidean distance on standardized scores
  dmat <- f_diss(
    Xr = scores_xr,
    Xu = scores_xu,
    diss_method = "euclid",
    center = FALSE,
    scale = FALSE
  )

  result <- list(
    dissimilarity = dmat,
    ncomp = proj$ncomp
  )

  if (method$return_projection) {
    result$projection <- proj
  }

  result
}


.diss_pls_compute <- function(Xr, Xu, Yr, method) {
  if (is.null(Yr)) {
    stop("'Yr' is required for PLS dissimilarity.")
  }

  proj <- ortho_projection(
    Xr = Xr,
    Xu = Xu,
    Yr = Yr,
    ncomp = method$ncomp,
    method = method$method,
    center = TRUE,
    scale = method$scale
  )

  scores <- proj$scores
  n_xr <- nrow(Xr)

  scores_xr <- scores[seq_len(n_xr), , drop = FALSE]
  scores_xu <- if (!is.null(Xu)) {
    scores[(n_xr + 1L):nrow(scores), , drop = FALSE]
  } else {
    NULL
  }

  # Mahalanobis distance in PLS space
  dmat <- f_diss(
    Xr = scores_xr,
    Xu = scores_xu,
    diss_method = "mahalanobis",
    center = FALSE,
    scale = FALSE
  )

  result <- list(
    dissimilarity = dmat,
    ncomp = proj$ncomp
  )

  if (method$return_projection) {
    result$projection <- proj
  }
  
  result
}


.diss_euclidean_compute <- function(Xr, Xu, method) {
  dmat <- f_diss(
    Xr = Xr,
    Xu = Xu,
    diss_method = "euclid",
    center = method$center,
    scale = method$scale
  )

  list(dissimilarity = dmat)
}


.diss_mahalanobis_compute <- function(Xr, Xu, method) {
  dmat <- f_diss(
    Xr = Xr,
    Xu = Xu,
    diss_method = "mahalanobis",
    center = method$center,
    scale = method$scale
  )

  list(dissimilarity = dmat)
}


.diss_cosine_compute <- function(Xr, Xu, method) {
  dmat <- f_diss(
    Xr = Xr,
    Xu = Xu,
    diss_method = "cosine",
    center = method$center,
    scale = method$scale
  )

  list(dissimilarity = dmat)
}


.diss_correlation_compute <- function(Xr, Xu, method) {
  dmat <- .cor_diss_compute(
    Xr = Xr,
    Xu = Xu,
    method = method
  )

  result <- list(dissimilarity = dmat)

  if (!is.null(method$ws)) {
    result$ws <- method$ws
  }

  result
}

# =============================================================================
# Helper: Coerce ncomp argument
# =============================================================================

.coerce_ncomp <- function(ncomp) {
  if (is.numeric(ncomp) && length(ncomp) == 1L) {
    if (is.na(ncomp) || ncomp < 1L) {
      stop("'ncomp' must be a positive integer.")
    }
    return(ncomp_fixed(as.integer(ncomp)))
  }
  if (!inherits(ncomp, "ncomp_selection")) {
    stop(
      "'ncomp' must be a positive integer or an ncomp_*() object.\n",
      "See ?ncomp_by_var, ?ncomp_by_cumvar, ?ncomp_by_opc, ?ncomp_fixed."
    )
  }
  ncomp
}


# =============================================================================
# Print method
# =============================================================================

#' @export
print.resemble_dissimilarity <- function(x, ...) {
  dmat <- x$dissimilarity
  cat("Dissimilarity matrix\n")
  cat("  Dimensions:", nrow(dmat), "x", ncol(dmat), "\n\n")
  cat("Constructor:\n")
  print(x$diss_method)
  if (!is.null(x$ncomp)) {
    cat("\nUsed ncomp: ", x$ncomp, "\n")
  }
  invisible(x)
}
