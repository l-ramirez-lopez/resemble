#' @title Control parameters for gesearch
#'
#' @description
#' Creates a control object specifying algorithm parameters for
#' \code{\link{gesearch}}.
#'
#' @param retain_by A character string specifying how training observations are
#'   selected at each iteration:
#'   \describe{
#'     \item{\code{"probability"}}{(default) Retains observations with errors
#'       below a percentile estimated from a given probability. More robust to
#'       outliers.}
#'     \item{\code{"proportion"}}{Retains a fixed proportion of observations
#'       with lowest errors.}
#'   }
#' @param percentile_type An integer between 1 and 9 specifying the quantile
#'   algorithm when \code{retain_by = "probability"}. Passed to the \code{type}
#'   argument of \code{\link[stats]{quantile}}. Default is 7. Ignored when
#'   \code{retain_by = "proportion"}.
#' @param tune A logical indicating whether to tune regression parameters
#'   (e.g., number of PLS components) via cross-validation at each iteration.
#'   Increases computation time substantially. Default is \code{FALSE}.
#' @param number An integer specifying the number of groups for leave-group-out
#'   cross-validation when \code{tune = TRUE}. Default is 10.
#' @param p A numeric value in (0, 1) specifying the proportion of observations
#'   per group in leave-group-out cross-validation when \code{tune = TRUE}.
#'   Default is 0.75.
#' @param stagnation_limit An integer specifying the maximum number of
#'   consecutive iterations with no change in gene pool size before early
#'   termination. Prevents infinite loops when target size cannot be reached.
#'   Default is 5.
#' @param allow_parallel A logical indicating whether to enable parallel
#'   processing for internal resampling and calibration. The parallel backend
#'   must be registered by the user. Default is \code{TRUE}.
#' @param blas_threads An integer specifying the number of BLAS threads to use
#'   during computation. Default is 1, which avoids multi-threaded OpenBLAS
#'   overhead on Linux. Requires \pkg{RhpcBLASctl}. See Details.
#' @details
#' ## Retention strategies
#'
#' When \code{retain_by = "probability"} (default), observations with errors
#' below a percentile threshold are retained. The percentile is computed using
#' \code{\link[stats]{quantile}} with \code{probs} set to the \code{retain}
#' value from \code{\link{gesearch}}. This approach is more robust when outlier
#' observations have extreme error values.
#'
#' When \code{retain_by = "proportion"}, a fixed fraction of observations
#' (specified by the \code{retain} argument in \code{\link{gesearch}}) with the
#' lowest associated errors are kept at each iteration.
#'
#' ## Cross-validation for tuning
#'
#' When \code{tune = TRUE}, leave-group-out cross-validation is used to select
#' optimal regression parameters at each iteration. The \code{number} argument
#' controls how many CV groups are formed, and \code{p} controls the proportion
#' of observations in each group.
#' 
#' ## BLAS threading
#'
#' On Linux systems with multi-threaded OpenBLAS, the default thread count
#' can cause significant overhead for algorithms that perform many small
#' matrix operations (like the iterative PLS fits in \code{gesearch}).
#' Setting \code{blas_threads = 1} (the default) eliminates this overhead.
#'
#' This setting requires the \pkg{RhpcBLASctl} package. If not installed,
#' the parameter is ignored and a message is displayed. The original thread
#' count is restored when \code{\link{gesearch}} completes.
#' 
#' @return A list of class \code{"gesearch_control"} containing the specified
#'   parameters.
#'
#' @author
#' \href{https://orcid.org/0000-0002-5369-5120}{Leonardo Ramirez-Lopez}
#'
#' @references
#' Lobsey, C.R., Viscarra Rossel, R.A., Roudier, P., Hedley, C.B. 2017.
#' rs-local data-mines information from spectral libraries to improve local
#' calibrations. European Journal of Soil Science 68:840-852.
#'
#' @seealso \code{\link{gesearch}}, \code{\link[stats]{quantile}}
#'
#' @examples
#' # Default parameters (probability-based retention)
#' gesearch_control()
#'
#' # Proportion-based retention
#' gesearch_control(retain_by = "proportion")
#'
#' # Enable parameter tuning with custom CV settings
#' gesearch_control(tune = TRUE, number = 5, p = 0.8)
#'
#' @export

gesearch_control <- function(
    # Retention strategy
  retain_by = c("probability", "proportion"),
  percentile_type = 7L,
  # Cross-validation tuning
  tune = FALSE,
  number = 10L,
  p = 0.75,
  # Algorithm behavior
  stagnation_limit = 5L,
  allow_parallel = TRUE, 
  blas_threads = 1L
) {
  
  retain_by <- match.arg(retain_by)
  
  # Validate percentile_type
  if (retain_by == "probability") {
    if (!is.numeric(percentile_type) || length(percentile_type) != 1L ||
        !percentile_type %in% 1L:9L) {
      stop("'percentile_type' must be an integer between 1 and 9", call. = FALSE)
    }
    percentile_type <- as.integer(percentile_type)
  } else {
    percentile_type <- NULL
  }
  
  # Validate tune
  if (!is.logical(tune) || length(tune) != 1L || is.na(tune)) {
    stop("'tune' must be TRUE or FALSE", call. = FALSE)
  }
  
  # Validate number
  if (!is.numeric(number) || length(number) != 1L || number < 1L) {
    stop("'number' must be a positive integer", call. = FALSE)
  }
  number <- as.integer(number)
  
  # Validate p
  if (!is.numeric(p) || length(p) != 1L || p <= 0 || p >= 1) {
    stop("'p' must be a numeric value in (0, 1)", call. = FALSE)
  }
  
  # Validate stagnation_limit
  if (!is.numeric(stagnation_limit) || length(stagnation_limit) != 1L ||
      stagnation_limit < 1L || stagnation_limit != as.integer(stagnation_limit)) {
    stop("'stagnation_limit' must be a positive integer", call. = FALSE)
  }
  stagnation_limit <- as.integer(stagnation_limit)
  
  # Validate allow_parallel
  if (!is.logical(allow_parallel) || length(allow_parallel) != 1L ||
      is.na(allow_parallel)) {
    stop("'allow_parallel' must be TRUE or FALSE", call. = FALSE)
  }
  
  structure(
    list(
      retain_by = retain_by,
      percentile_type = percentile_type,
      tune = tune,
      number = number,
      p = p,
      stagnation_limit = stagnation_limit,
      allow_parallel = allow_parallel, 
      blas_threads = blas_threads
    ),
    class = "gesearch_control"
  )
}