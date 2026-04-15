#' @title Control parameters for memory-based learning
#' @description
#' \loadmathjax
#' This function controls various aspects of the memory-based learning process
#' in the \code{\link{mbl}} function.
#'
#' @usage
#' mbl_control(
#'   return_dissimilarity = FALSE,
#'   validation_type = "NNv",
#'   tune_locally = TRUE,
#'   number = 10,
#'   p = 0.75,
#'   range_prediction_limits = TRUE,
#'   allow_parallel = TRUE,
#'   blas_threads = 1L
#' )
#'
#' @param return_dissimilarity Logical indicating whether to return the
#'   dissimilarity matrix between \code{Xr} and \code{Xu}. Default is
#'   \code{FALSE}.
#' @param validation_type Character vector specifying validation method(s):
#'   \itemize{
#'     \item \code{"NNv"}: Leave-nearest-neighbor-out cross-validation (default, faster)
#'     \item \code{"local_cv"}: Local leave-group-out cross-validation
#'     \item \code{"none"}: No validation
#'   }
#'   Multiple methods can be specified (e.g., \code{c("NNv", "local_cv")}).
#'   Default is \code{"NNv"}.
#' @param tune_locally Logical indicating whether to tune PLS components
#'   locally when \code{validation_type = "local_cv"} and using
#'   \code{\link{fit_pls}} or \code{\link{fit_wapls}}. Default is \code{TRUE}.
#' @param number Integer specifying the number of sampling iterations for
#'   \code{"local_cv"} validation. Default is \code{10}.
#' @param p Numeric value between 0 and 1 indicating the proportion of
#'   observations retained at each \code{"local_cv"} iteration. Default is
#'   \code{0.75}.
#' @param range_prediction_limits Logical indicating whether predictions should
#'   be constrained to the range of response values in each neighborhood.
#'   Default is \code{TRUE}.
#' @param allow_parallel Logical indicating whether parallel execution is
#'   allowed via the \pkg{foreach} package. Default is \code{TRUE}.
#' @param blas_threads Integer specifying the number of BLAS threads to use
#'   during \code{mbl()} execution. Default is \code{1L}, which avoids thread
#'   overhead from repeated small matrix operations. Requires the
#'   \pkg{RhpcBLASctl} package to take effect. The original thread count is
#'   restored after \code{mbl()} completes. See Details.
#'
#' @details
#' ## Validation methods
#'
#' \strong{Leave-nearest-neighbor-out cross-validation (\code{"NNv"}):}
#' For each target observation, the nearest neighbor is excluded from the
#' local model, which then predicts that neighbor's value. This is faster
#' than \code{"local_cv"}. If the nearest neighbor belongs to a group
#' (specified via the \code{group} argument in \code{\link{mbl}}), all
#' group members are excluded.
#'
#' \strong{Local leave-group-out cross-validation (\code{"local_cv"}):}
#' The neighborhood is partitioned into subsets via stratified random
#' sampling. Each subset serves as validation data while the remainder
#' fits the model. This repeats \code{number} times, with \code{p}
#' controlling the training proportion. The final error is the average
#' local RMSE.
#'
#' ## BLAS threading
#'
#' On Linux systems with multi-threaded OpenBLAS, the default thread count
#' can cause significant overhead for algorithms like \code{mbl()} that
#' perform many small matrix operations. Setting \code{blas_threads = 1}
#' (the default) eliminates this overhead.
#'
#' This setting requires the \pkg{RhpcBLASctl} package. If not installed,
#' the parameter is ignored and a message is displayed. The original
#' thread count is restored when \code{mbl()} completes.
#'
#' Windows systems typically use single-threaded BLAS by default, so this
#' setting has no effect there.
#'
#' @return A list of class \code{mbl_control} with the specified control parameters.
#'
#' @author
#' \href{https://orcid.org/0000-0002-5369-5120}{Leonardo Ramirez-Lopez} and
#' Antoine Stevens
#'
#' @references
#' Ramirez-Lopez, L., Behrens, T., Schmidt, K., Stevens, A., Dematte, J.A.M.,
#' Scholten, T. 2013a. The spectrum-based learner: A new local approach for
#' modeling soil vis-NIR spectra of complex data sets. Geoderma 195-196:268-279.
#'
#' Ramirez-Lopez, L., Behrens, T., Schmidt, K., Viscarra Rossel, R., Dematte,
#' J.A.M., Scholten, T. 2013b. Distance and similarity-search metrics for use
#' with soil vis-NIR spectra. Geoderma 199:43-53.
#'
#' @seealso \code{\link{mbl}}, \code{\link{neighbors_k}},
#'   \code{\link{neighbors_diss}}
#'
#' @examples
#' # Default control parameters (NNv validation)
#' mbl_control()
#'
#' # Both validation methods
#' mbl_control(validation_type = c("NNv", "local_cv"))
#'
#' # No validation
#' mbl_control(validation_type = "none")
#'
#' # NNv validation only, no parallel
#' mbl_control(validation_type = "NNv", allow_parallel = FALSE)
#'
#' # Allow more BLAS threads (if needed for other computations)
#' mbl_control(blas_threads = 4)
#'
#' @export
mbl_control <- function(
    return_dissimilarity = FALSE,
    validation_type = "NNv",
    tune_locally = TRUE,
    number = 10,
    p = 0.75,
    range_prediction_limits = TRUE,
    allow_parallel = TRUE,
    blas_threads = 1L
) {
  
  # Validate return_dissimilarity
  if (!is.logical(return_dissimilarity) || length(return_dissimilarity) != 1L) {
    stop("'return_dissimilarity' must be TRUE or FALSE.", call. = FALSE)
  }
  
  # Validate and match validation_type
  valid_types <- c("NNv", "local_cv", "none")
  validation_type <- match.arg(validation_type, valid_types, several.ok = TRUE)
  
  if ("none" %in% validation_type && length(validation_type) > 1L) {
    stop(
      "'validation_type' cannot combine 'none' with other values.",
      call. = FALSE
    )
  }
  
  # Remove "none" if other methods specified
  if ("none" %in% validation_type && length(validation_type) > 1L) {
    validation_type <- validation_type[validation_type != "none"]
  }
  
  # Validate tune_locally
  if (!is.logical(tune_locally) || length(tune_locally) != 1L) {
    stop("'tune_locally' must be TRUE or FALSE.", call. = FALSE)
  }
  
  # Validate local_cv parameters
  if ("local_cv" %in% validation_type) {
    if (!is.numeric(number) || length(number) != 1L || number < 1L) {
      stop("'number' must be a positive integer.", call. = FALSE)
    }
    
    if (!is.numeric(p) || length(p) != 1L || p <= 0 || p >= 1) {
      stop("'p' must be a numeric value between 0 and 1 (exclusive).",
           call. = FALSE)
    }
  }
  
  # Validate range_prediction_limits
  if (!is.logical(range_prediction_limits) ||
      length(range_prediction_limits) != 1L) {
    stop("'range_prediction_limits' must be TRUE or FALSE.", call. = FALSE)
  }
  
  # Validate allow_parallel
  if (!is.logical(allow_parallel) || length(allow_parallel) != 1L) {
    stop("'allow_parallel' must be TRUE or FALSE.", call. = FALSE)
  }
  
  # Validate blas_threads
  if (!is.numeric(blas_threads) || length(blas_threads) != 1L ||
      blas_threads < 1L) {
    stop("'blas_threads' must be a positive integer.", call. = FALSE)
  }
  
  ctl <- list(
    return_dissimilarity = return_dissimilarity,
    validation_type = validation_type,
    tune_locally = tune_locally,
    number = as.integer(number),
    p = p,
    range_prediction_limits = range_prediction_limits,
    allow_parallel = allow_parallel,
    blas_threads = as.integer(blas_threads)
  )
  class(ctl) <- c("mbl_control", "list")
  ctl
}
