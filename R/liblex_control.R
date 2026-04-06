#' @title Control parameters for liblex
#'
#' @description
#' Specifies control parameters for the \code{\link{liblex}} function, including
#' output options, validation settings, tuning behavior, and parallel execution.
#'
#' @param return_dissimilarity A logical indicating whether the dissimilarity
#'   matrix should be returned in the output. Default is \code{FALSE}. Setting
#'   to \code{TRUE} can be useful for diagnostics but increases memory usage
#'   for large libraries.
#' @param mode A character string specifying the operation mode:
#'   \describe{
#'     \item{\code{"build"}}{(default) Builds the library of local models. If
#'       \code{tune = TRUE}, validation is performed first to find optimal
#'       parameters, then the library is built using those parameters. If
#'       \code{tune = FALSE}, the library is built directly using the
#'       parameters provided to \code{\link{liblex}}.}
#'     \item{\code{"validate"}}{Performs validation only without building the
#'       library. Useful for parameter exploration or when testing different
#'       configurations before committing to the full library build.}
#'   }
#' @param tune A logical indicating whether to optimize parameters via
#'   nearest-neighbor validation. Default is \code{FALSE}. When \code{TRUE},
#'   the function evaluates all combinations of \code{k} values and PLS
#'   component ranges, selecting the combination that minimizes RMSE (or
#'   maximizes R², depending on \code{metric}). See Details.
#' @param metric A character string specifying the performance metric used
#'   for parameter selection when \code{tune = TRUE}. Options are:
#'   \describe{
#'     \item{\code{"rmse"}}{(default) Root mean squared error (minimized).}
#'     \item{\code{"r2"}}{Coefficient of determination (maximized).}
#'   }
#'   Ignored when \code{tune = FALSE}.
#' @param chunk_size An integer specifying the number of local models to
#'   process per parallel task. Default is \code{1L}. Increasing this value
#'   reduces parallel overhead but may cause load imbalance. For large
#'   libraries, values between 10 and 50 often provide a good trade-off
#'   between overhead and efficiency.
#' @param allow_parallel A logical indicating whether parallel execution is
#'   permitted. Default is \code{TRUE}. Parallelization is applied to the
#'   model fitting loop using the \pkg{foreach} package. Requires a registered
#'   parallel backend (e.g., via \pkg{doParallel}).
#' @param blas_threads An integer specifying the number of threads for BLAS
#'   operations. Default is \code{1L}, which avoids thread contention when
#'   using parallel processing. Requires the \pkg{RhpcBLASctl} package to
#'   take effect. On Linux systems with multi-threaded BLAS (e.g., OpenBLAS),
#'   setting this to 1 can substantially improve performance when
#'   \code{allow_parallel = TRUE}.
#'
#' @details
#' \subsection{Nearest-neighbor validation}{
#' When \code{tune = TRUE} or \code{mode = "validate"}, the function performs
#' nearest-neighbor validation (NNv) to assess model performance. For each
#' observation in the reference set (or anchor set, if specified), the
#' procedure:
#' \enumerate{
#'   \item Identifies the k nearest neighbors of the target observation.
#'   \item Excludes the target observation (and any observations in the same
#'     group, if \code{group} is specified) from the neighbor set.
#'   \item Fits a local model using the remaining neighbors.
#'   \item Predicts the response value for the excluded target observation.
#'   \item Computes prediction errors across all observations.
#' }
#'
#' This leave-one-out style validation provides an estimate of prediction
#' performance without requiring a separate test set. When \code{tune = TRUE},
#' the parameter combination (number of neighbors and PLS component range)
#' yielding the best performance according to \code{metric} is selected for
#' building the final library.
#' }
#'
#' \subsection{Mode and tune combinations}{
#' \tabular{lll}{
#'   \code{mode} \tab \code{tune} \tab Behavior \cr
#'   \code{"build"} \tab \code{FALSE} \tab Build library using parameters as
#'     provided \cr
#'   \code{"build"} \tab \code{TRUE} \tab Validate, find optimal parameters,
#'     build library \cr
#'   \code{"validate"} \tab \code{FALSE} \tab Validate only, report
#'     performance statistics \cr
#'   \code{"validate"} \tab \code{TRUE} \tab Validate, report statistics with
#'     optimal parameters identified \cr
#' }
#' }
#'
#' \subsection{Parallel chunk size}{
#' The \code{chunk_size} parameter controls granularity of parallel work
#' distribution. When \code{allow_parallel = TRUE} and a parallel backend is
#' registered:
#' \itemize{
#'   \item \code{chunk_size = 1}: Each local model is a separate parallel task.
#'     Maximum parallelism but higher scheduling overhead.
#'   \item \code{chunk_size > 1}: Multiple models are processed sequentially
#'     within each parallel task. Reduces overhead and improves memory
#'     locality, but may cause load imbalance if the number of models is not
#'     evenly divisible.
#' }
#' When \code{allow_parallel = FALSE}, \code{chunk_size} has no effect.
#' }
#'
#' @return A list of class \code{liblex_control} containing the validated
#'   control parameters.
#'
#' @seealso \code{\link{liblex}}, \code{\link{mbl_control}}
#'
#' @author \href{https://orcid.org/0000-0002-5369-5120}{Leonardo Ramirez-Lopez}
#'
#' @examples
#' # Default settings: build library without tuning
#' liblex_control()
#'
#' # Tune parameters before building
#' liblex_control(tune = TRUE)
#'
#' # Validate only for parameter exploration
#' liblex_control(mode = "validate", tune = TRUE)
#'
#' # Include dissimilarity matrix in output
#' liblex_control(return_dissimilarity = TRUE)
#'
#' # Larger chunks for reduced parallel overhead
#' liblex_control(chunk_size = 20L)
#'
#' # Parallel settings for Linux with OpenBLAS
#' liblex_control(allow_parallel = TRUE, blas_threads = 1L)
#'
#' @export
liblex_control <- function(
    return_dissimilarity = FALSE,
    mode = c("build", "validate"),
    tune = FALSE,
    metric = c("rmse", "r2"),
    chunk_size = 1L,
    allow_parallel = TRUE,
    blas_threads = 1L
) {
  mode <- match.arg(mode)
  metric <- match.arg(metric)
  
  if (!is.logical(return_dissimilarity) || length(return_dissimilarity) != 1L ||
      is.na(return_dissimilarity)) {
    stop("'return_dissimilarity' must be TRUE or FALSE", call. = FALSE)
  }
  
  if (!is.logical(tune) || length(tune) != 1L || is.na(tune)) {
    stop("'tune' must be TRUE or FALSE", call. = FALSE)
  }
  
  if (!is.numeric(chunk_size) || length(chunk_size) != 1L ||
      is.na(chunk_size) || chunk_size < 1L) {
    stop("'chunk_size' must be a positive integer", call. = FALSE)
  }
  chunk_size <- as.integer(chunk_size)
  
  if (!is.logical(allow_parallel) || length(allow_parallel) != 1L ||
      is.na(allow_parallel)) {
    stop("'allow_parallel' must be TRUE or FALSE", call. = FALSE)
  }
  
  if (!is.numeric(blas_threads) || length(blas_threads) != 1L ||
      is.na(blas_threads) || blas_threads < 1L) {
    stop("'blas_threads' must be a positive integer", call. = FALSE)
  }
  blas_threads <- as.integer(blas_threads)
  
  structure(
    list(
      return_dissimilarity = return_dissimilarity,
      mode = mode,
      tune = tune,
      metric = metric,
      chunk_size = chunk_size,
      allow_parallel = allow_parallel,
      blas_threads = blas_threads
    ),
    class = "liblex_control"
  )
}
