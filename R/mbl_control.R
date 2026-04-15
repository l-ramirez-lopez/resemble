#' @title A function that controls some few aspects of the memory-based learning
#' process in the \code{mbl} function
#' @description
#' \strong{Maturing}
#'
#' This function is used to further control some aspects of the memory-based
#' learning process in the \code{mbl} function.
#' @usage
#' mbl_control(
#'   return_dissimilarity = FALSE,
#'   validation_type = c("NNv", "local_cv"),
#'   tune_locally = TRUE,
#'   number = 10,
#'   p = 0.75,
#'   range_prediction_limits = TRUE,
#'   progress = TRUE,
#'   allow_parallel = TRUE
#' )
#' @param return_dissimilarity A logical indicating if the dissimilarity matrix
#' between \code{Xr} and \code{Xu} must be returned.
#' @param validation_type A character vector indicating the internal validation
#' method(s) used to assess the global performance of the local models.
#' Possible options are \code{"NNv"} and \code{"local_cv"}. Alternatively,
#' \code{"none"} can be used when cross-validation is not required.
#' @param tune_locally A logical. It only applies when
#' \code{validation_type = "local_cv"} and \code{"pls"} or \code{"wapls"}
#' fitting algorithms are used. If \code{TRUE}, the parameters of the local
#' PLS-based models are tuned.
#' @param number An integer indicating the number of sampling iterations at
#' each local segment when \code{"local_cv"} is selected. Default is 10.
#' @param p A numeric value indicating the percentage of observations to be retained
#' at each sampling iteration at each local segment when \code{"local_cv"}
#' is selected. Default is 0.75 (75\%).
#' @param range_prediction_limits A logical. It indicates whether the prediction
#' limits at each local regression are determined by the range of the response
#' variable within each neighborhood. When the predicted value is outside
#' this range, it is automatically replaced with the nearest range value.
#' If \code{FALSE}, no prediction limits are imposed. Default is \code{TRUE}.
#' @param progress A logical indicating whether to print a progress bar
#' for each observation to be predicted. Default is \code{TRUE}. In case
#' parallel processing is used, these progress bars are not printed.
#' @param allow_parallel A logical indicating if parallel execution is allowed.
#' If \code{TRUE}, parallelism is applied to the loop in \code{\link{mbl}}
#' in which each iteration handles a single observation in \code{Xu}.
#' @details
#' The validation methods available for assessing predictive performance are:
#'
#' \itemize{
#'   \item \code{"NNv"}: Leave-nearest-neighbor-out cross-validation. From
#'   the group of neighbors of each observation to be predicted, the nearest
#'   observation is excluded and then a local model is fitted using the
#'   remaining neighbors. This model is then used to predict the response
#'   value of the nearest observation. These predicted values are finally
#'   cross-validated against the actual values. If the nearest sample belongs
#'   to a group of samples labeled through the \code{group} argument in
#'   \code{\link{mbl}}, then all samples in that group are excluded from the
#'   temporary local fit used for validation. This validation method is faster
#'   than \code{"local_cv"}.
#'
#'   \item \code{"local_cv"}: Local leave-group-out cross-validation. The
#'   group of neighbors of each observation to be predicted is partitioned into
#'   equal-size subsets. Each partition is selected by stratified random
#'   sampling based on response values. The selected subset is used as a local
#'   validation subset and the remaining observations are used for fitting a
#'   model. This process is repeated \eqn{m} times. In \code{\link{mbl}},
#'   \eqn{m} is controlled by the \code{number} argument and subset size is
#'   controlled by the \code{p} argument. The global error is computed as the
#'   average of the local root mean square errors.
#'
#'   \item \code{"none"}: No validation is carried out. If \code{"none"} is
#'   selected along with \code{"NNv"} or \code{"local_cv"}, it is ignored.
#' }
#' @return A \code{list} mirroring the specified parameters.
#' @author Leonardo Ramirez-Lopez and Antoine Stevens
#' @references
#' Ramirez-Lopez, L., Behrens, T., Schmidt, K., Stevens, A., Dematte, J.A.M.,
#' Scholten, T. 2013a. The spectrum-based learner: A new local approach for
#' modeling soil vis-NIR spectra of complex data sets. Geoderma 195-196, 268-279.
#'
#' Ramirez-Lopez, L., Behrens, T., Schmidt, K., Viscarra Rossel, R., Dematte,
#' J. A. M., Scholten, T. 2013b. Distance and similarity-search metrics for
#' use with soil vis-NIR spectra. Geoderma 199, 43-53.
#' @seealso \code{\link{f_diss}}, \code{\link{cor_diss}}, \code{\link{sid}},
#' \code{\link{ortho_diss}}, \code{\link{mbl}}
#' @examples
#' mbl_control()
#' @export

mbl_control <- function(
    return_dissimilarity = FALSE,
    validation_type = c("NNv", "local_cv"),
    tune_locally = TRUE,
    number = 10,
    p = 0.75,
    range_prediction_limits = TRUE,
    progress = TRUE,
    allow_parallel = TRUE
) {
  # Sanity checks
  if (!is.logical(allow_parallel)) {
    stop("allow_parallel must be a logical value")
  }
  
  if (!is.logical(return_dissimilarity)) {
    stop("'return_dissimilarity' must be logical")
  }
  
  if (sum(validation_type %in% c("NNv", "local_cv", "none")) != length(validation_type)) {
    stop("'validation_type' must be one at least one of 'NNv', 'local_cv', 'none'")
  }
  
  if ("none" %in% validation_type) {
    if (length(validation_type) > 1) {
      validation_type <- validation_type[!(validation_type == "none")]
    }
  }
  
  if ("local_cv" %in% validation_type) {
    if (!is.numeric(number) | length(number) != 1) {
      stop("The 'number' argument must be a single numeric value")
    }
    
    if (!is.numeric(p) | length(p) != 1 | p >= 1 | p <= 0) {
      stop("p must be a single numeric value larger than 0 and below than 1")
    }
  }
  
  if (!is.logical(range_prediction_limits)) {
    stop("'range_prediction_limits' must be logical")
  }
  
  if (!is.logical(progress)) {
    stop("'progress' must be logical")
  }
  cntrl <- list(
    return_dissimilarity = return_dissimilarity,
    validation_type = validation_type,
    tune_locally = tune_locally,
    number = number,
    p = p,
    range_prediction_limits = range_prediction_limits,
    progress = progress,
    allow_parallel = allow_parallel
  )
  
  cntrl
}
