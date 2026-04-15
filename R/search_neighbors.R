#' @title Search neighbors in a reference set
#'
#' @description
#' \loadmathjax
#' Searches for the nearest neighbors of observations in a reference set or
#' between two sets of observations.
#' 
#' @usage
#' search_neighbors(Xr, Xu = NULL,
#'                  diss_method = diss_pca(), Yr = NULL,
#'                  neighbors, spike = NULL,
#'                  return_dissimilarity = FALSE, 
#'                  k, k_diss, k_range, pc_selection, 
#'                  center, scale, documentation, ...
#'                  )
#'     
#' @param Xr A numeric matrix of reference observations (rows) and variables
#'   (columns) where the neighbor search is conducted.
#' @param Xu Optional matrix of observations for which neighbors are to be
#'   searched in \code{Xr}.
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
#' @param neighbors A neighbor selection object created by:
#'   \itemize{
#'     \item \code{\link{neighbors_k}()}: Select k nearest neighbors
#'     \item \code{\link{neighbors_diss}()}: Select neighbors by dissimilarity threshold
#'   }
#' @param spike Optional integer vector indicating observations in \code{Xr}
#'   to force into (positive indices) or exclude from (negative indices)
#'   neighborhoods.
#' @param return_dissimilarity Logical indicating whether to return the
#'   dissimilarity matrix. Default is \code{FALSE}.
#'   
#' @param k Deprecated.
#' @param k_diss Deprecated.
#' @param k_range Deprecated.
#' @param pc_selection Deprecated.
#' @param center Deprecated.
#' @param scale Deprecated.
#' @param documentation Deprecated.
#' @param ... Additional arguments (currently unused).
#'
#' @details
#' This function is useful for reducing large reference sets by identifying
#' only relevant neighbors before running \code{\link{mbl}}.
#'
#' If \code{Xu} is not provided, the function searches for neighbors within
#' \code{Xr} itself (excluding self-matches). If \code{Xu} is provided,
#' neighbors of each observation in \code{Xu} are searched in \code{Xr}.
#'
#' The \code{spike} argument allows forcing specific observations into or out
#' of all neighborhoods. Positive indices are always included; negative indices
#' are always excluded.
#'
#' @return A list containing:
#' \describe{
#'   \item{neighbors}{Matrix of \code{Xr} indices for each query observation's
#'     neighbors, sorted by dissimilarity (columns = query observations).}
#'   \item{neighbors_diss}{Matrix of dissimilarity scores corresponding to
#'     \code{neighbors}.}
#'   \item{unique_neighbors}{Vector of unique \code{Xr} indices that appear
#'     in any neighborhood.}
#'   \item{k_diss_info}{If \code{neighbors_diss()} was used, a \code{data.frame} 
#'     with columns for observation index, number of neighbors found, and final
#'     number after applying bounds.}
#'   \item{dissimilarity}{If \code{return_dissimilarity = TRUE}, the full
#'     dissimilarity object.}
#'   \item{projection}{If the dissimilarity method includes
#'     \code{return_projection = TRUE}, the projection object.}
#'   \item{gh}{If the dissimilarity method includes \code{gh = TRUE}, the
#'     GH distances.}
#' }
#'
#' @seealso \code{\link{dissimilarity}}, \code{\link{mbl}}, 
#'   \code{\link{neighbors_k}}, \code{\link{neighbors_diss}}
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
#' @author
#' \href{https://orcid.org/0000-0002-5369-5120}{Leonardo Ramirez-Lopez}
#'
#' @examples
#' \donttest{
#' library(prospectr)
#' data(NIRsoil)
#'
#' Xu <- NIRsoil$spc[!as.logical(NIRsoil$train), ]
#' Yu <- NIRsoil$CEC[!as.logical(NIRsoil$train)]
#' Yr <- NIRsoil$CEC[as.logical(NIRsoil$train)]
#' Xr <- NIRsoil$spc[as.logical(NIRsoil$train), ]
#'
#' Xu <- Xu[!is.na(Yu), ]
#' Yu <- Yu[!is.na(Yu)]
#' Xr <- Xr[!is.na(Yr), ]
#' Yr <- Yr[!is.na(Yr)]
#'
#' # Correlation-based neighbor search with k neighbors
#' ex1 <- search_neighbors(
#'   Xr = Xr, Xu = Xu,
#'   diss_method = diss_correlation(),
#'   neighbors = neighbors_k(40)
#' )
#'
#' # PCA-based with OPC selection
#' ex2 <- search_neighbors(
#'   Xr = Xr, Xu = Xu,
#'   diss_method = diss_pca(
#'     ncomp = ncomp_by_opc(40),
#'     scale = TRUE,
#'     return_projection = TRUE
#'   ),
#'   Yr = Yr,
#'   neighbors = neighbors_k(50)
#' )
#'
#' # Observations not in any neighborhood
#' setdiff(seq_len(nrow(Xr)), ex2$unique_neighbors)
#'
#' # Dissimilarity threshold-based selection
#' ex3 <- search_neighbors(
#'   Xr = Xr, Xu = Xu,
#'   diss_method = diss_pls(
#'     ncomp = ncomp_by_opc(40),
#'     scale = TRUE
#'   ),
#'   Yr = Yr,
#'   neighbors = neighbors_diss(threshold = 0.5, k_min = 10, k_max = 100)
#' )
#' }
#'
#' @export
search_neighbors <- function(
    Xr,
    Xu = NULL,
    diss_method = diss_pca(),
    Yr = NULL,
    neighbors,
    spike = NULL,
    return_dissimilarity = FALSE,
    # --- Deprecated arguments ---
    k, k_diss, k_range, pc_selection, center, scale, documentation, ...
) {
  
  
  # ---------------------------------------------------------------------------
  # Block removed arguments
  # ---------------------------------------------------------------------------
  if (!missing(k) || !missing(k_diss) || !missing(k_range)) {
    stop(
      "Arguments 'k', 'k_diss', 'k_range' have been removed.\n",
      "Use neighbors_k() or neighbors_diss() instead.\n",
      "Example: search_neighbors(..., neighbors = neighbors_k(40))",
      call. = FALSE
    )
  }
  
  if (!missing(pc_selection)) {
    stop(
      "Argument 'pc_selection' has been removed.\n",
      "Component selection is now specified in diss_*() constructors.\n",
      "Example: diss_pca(ncomp = ncomp_by_opc(40))",
      call. = FALSE
    )
  }
  
  if (!missing(center) || !missing(scale)) {
    stop(
      "Arguments 'center' and 'scale' have been removed.\n",
      "These are now set in diss_*() constructors.\n",
      "Example: diss_pca(center = TRUE, scale = FALSE)",
      call. = FALSE
    )
  }
  
  if (!missing(documentation)) {
    stop(
      "Argument 'documentation' has been removed.",
      call. = FALSE
    )
  }
  
  if (length(list(...)) > 0) {
    stop(
      "Unknown arguments passed via '...'.\n",
      "Arguments like 'ws', 'return_projection', 'gh' are now set in ",
      "diss_*() constructors.",
      call. = FALSE
    )
  }
  
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
      "  diss_method = \"cor\"               -> diss_correlation()\n",
      call. = FALSE
    )
  }
  
  # ---------------------------------------------------------------------------
  # Validate diss_method
  # ---------------------------------------------------------------------------
  if (!inherits(diss_method, "diss_method")) {
    stop(
      "'diss_method' must be a dissimilarity method object.\n",
      "Use one of: diss_pca(), diss_pls(), diss_euclidean(), ",
      "diss_mahalanobis(), diss_cosine(), diss_correlation().",
      call. = FALSE
    )
  }
  
  # ---------------------------------------------------------------------------
  # Validate neighbors
  # ---------------------------------------------------------------------------
  if (missing(neighbors)) {
    stop(
      "'neighbors' is required. Use neighbors_k() or neighbors_diss().",
      call. = FALSE
    )
  }
  
  if (!inherits(neighbors, "neighbors")) {
    stop(
      "'neighbors' must be created by neighbors_k() or neighbors_diss().",
      call. = FALSE
    )
  }
  
  # ---------------------------------------------------------------------------
  # Convert neighbors to internal format
  # ---------------------------------------------------------------------------
  n_xr <- nrow(Xr)
  
  if (inherits(neighbors, "neighbors_k")) {
    k <- max(neighbors$k)
    k_diss <- NULL
    k_range <- NULL
    
    if (k > n_xr) {
      stop("'k' (", k, ") cannot exceed nrow(Xr) (", n_xr, ").", call. = FALSE)
    }
  } else {
    # neighbors_diss
    k <- NULL
    k_diss <- max(neighbors$threshold)
    k_range <- c(neighbors$k_min, neighbors$k_max)

    if (is.infinite(k_range[2])) {
      
      if (is.null(Xu)) {
        max_neighbors  <- n_xr
        mss <- "nrow(Xr) - 1 (excluding self-matches)"
      } else {
        max_neighbors <- n_xr - 1
        mss <- "nrow(Xr)"
      }
      message(
        "setting 'k_max' (", k_range[2], ") to ", mss,  " (", n_xr, ")."
      )
      k_range[2] <- n_xr
    }
        
    if (k_range[2] > n_xr) {
      stop(
        "'k_max' (", k_range[2], ") cannot exceed nrow(Xr) (", n_xr, ").",
        call. = FALSE
      )
    }
  }
  
  # ---------------------------------------------------------------------------
  # Validate spike
  # ---------------------------------------------------------------------------
  if (!is.null(spike)) {
    if (!is.vector(spike) || !all(spike == as.integer(spike))) {
      stop("'spike' must be an integer vector.", call. = FALSE)
    }
    if (length(spike) >= n_xr) {
      stop("'spike' length must be less than nrow(Xr).", call. = FALSE)
    }
    if (max(abs(spike)) > n_xr) {
      stop("'spike' contains out-of-bounds indices.", call. = FALSE)
    }
    
    spike_hold <- spike[spike > 0L]
    spike_rm <- abs(spike[spike < 0L])
    
    if (length(intersect(spike_hold, spike_rm)) > 0L) {
      stop("'spike' contains contradictory indices.", call. = FALSE)
    }
    
    min_k <- if (!is.null(k)) min(neighbors$k) else neighbors$k_min
    if (min_k <= length(spike_hold)) {
      stop(
        "Minimum k must exceed length of positive spike indices.",
        call. = FALSE
      )
    }
    
    spike <- sort(unique(as.integer(spike)))
  }
  
  # ---------------------------------------------------------------------------
  # Compute dissimilarity
  # ---------------------------------------------------------------------------
  dsm <- dissimilarity(
    Xr = Xr,
    Xu = Xu,
    diss_method = diss_method,
    Yr = Yr
  )
  
  # ---------------------------------------------------------------------------
  # Find neighbors
  # ---------------------------------------------------------------------------
  skip_first <- is.null(Xu)
  
  results <- diss_to_neighbors(
    dsm$dissimilarity,
    k = k,
    k_diss = k_diss,
    k_range = k_range,
    spike = spike,
    return_dissimilarity = return_dissimilarity,
    skip_first = skip_first
  )
  
  # Column names
  mprefix <- if (skip_first) "Xr" else "Xu"
  colnames(results$neighbors) <- paste0(mprefix, seq_len(ncol(results$neighbors)))
  
  # ---------------------------------------------------------------------------
  # Add optional outputs from dissimilarity
  # ---------------------------------------------------------------------------
  if (return_dissimilarity) {
    results$dissimilarity <- dsm
  }
  
  if (!is.null(dsm$projection)) {
    results$projection <- dsm$projection
  }
  
  if (!is.null(dsm$gh)) {
    results$gh <- dsm$gh
  }
  
  results
}