#' @title A function that controls some aspects of the the \code{rslocal} function
#' @description
#' 
#' \lifecycle{maturing}
#' 
#' This function is used to further control some aspects of the \code{rslocal} function.
#' 
#' @usage
#' rs_control(retain_by = "proportion",
#'            percentile_type = if(retain_by == "probability") 7,
#'            tune = FALSE,
#'            number = 10,
#'            p = 0.75,
#'            verbose = TRUE,
#'            allow_parallel = TRUE)
#' @param retain_by a character string indicating how the training observations
#' to be kept at each iteration must be selected. Two options are available:
#' \code{'proportion'} (default) retains a fix proportion of observations at each
#' iteration and \code{'probability'} retains observations whose errors are
#' below an estimated percentile (cut point) from a given probability.
#' See details.
#' @param percentile_type if \code{retain_by = 'probability'}, an integer
#' between 1 and 9 to be passed to the the \code{type} argument of the
#' \code{\link[stats]{quantile}} function to estimate the percentile cut-off
#' value used to select the observations at each iteration. Default is 7
#' (as in \code{\link[stats]{quantile}}). If \code{retain_by != 'proportion'},
#' this argument is ignored. See details.
#' @param tune a logical indicating whether to tune the parameters of the
#' regression method (e.g. pls factors). If \code{TRUE} the parameters are tuned
#' by using cross-validation which is based on the \code{number} and \code{p}
#' arguments. Note that cross-validation greatly increases processing time.
#' Default is \code{FALSE}.
#' @param number if \code{tune = TRUE}, an integer indicating the number of groups
#' for the leave-group-out internal cross-validation used for parameter tuning
#' (e.g. pls factors) at each iteration. These groups are built from the training
#' subset selected at the respective iteration. Default is 10.
#' @param p if \code{tune = TRUE}, a value indicating the percentage of
#' observations to build each group in the the leave-group-out internal
#' cross-validation (used for parameter tuning at each iteration). These groups
#' are built from the training subset selected at the respective iteration.
#' Default is 0.75 (i.e. 75 "\%").
#' @param number an integer indicating the number of sampling subsets/iterations
#' for cross-validation. Default is 10.
#' @param verbose a logical indicating if some information about each iteration
#' must be printed.
#' @param allow_parallel set to TRUE to parallelise the internal re-sampling
#' and calibration. The parallel backend should be registered by the user.
#' @details
#' The training observations to be kept at each iteration can be selected by
#' using a fixed proportion of observations (in relation to the number of
#' observations in the training subset being used at each iteration) which
#' associated errors are the lowest, in this case the argument \code{retain_by}
#' must be set to \code{'proportion'}. This proportion (a value larger than 0
#' and below 1) is given in the in the argument \code{retain} of the
#' \code{\link{rslocal}} function.
#'
#' For selecting the training observations, it is also possible to retain the
#' ones whose associated errors are below a percentile estimated for a given
#' probability (the percentile is estimated with the
#' \code{\link[stats]{quantile}} function). In this case the probability (a
#' value larger than 0 and below 1) must be provided in the argument
#' \code{retain}, this numeric value is passed to the argument \code{probs} in
#' the \code{\link[stats]{quantile}} function. This method might be useful in
#' the presence of outlier observations with considerable deviation in their
#' associated errors.
#'
#' The internal validations for tuning the number of pls factors are based on
#' the leave-group-out cross-validation. Arguments \code{p} and \code{number}
#' are used to control the number of groups and the amount of observations
#' per group.
#' @return a \code{list} mirroring the specified parameters.
#' @author Leonardo Ramirez-Lopez
#' @references
#' Lobsey, C. R., Viscarra Rossel, R. A., Roudier, P., & Hedley, C. B. 2017.
#' rs-local data-mines information from spectral libraries to improve local
#' calibrations. European Journal of Soil Science, 68(6), 840-852.
#' @seealso \code{\link{rslocal}},  \code{\link[stats]{quantile}}
#' @examples
#' # A control list with the default parameters
#' rs_control()
#'
#' # A control list which specifies that observations must be
#' # retained based on the associated errors that are below
#' # a given probability.
#' rs_control(retain_by = "probability")
#' @export

## History:
## 2020.03.29 (Leo):  New fucntion to control rslocal
## 2020.03.29 (Leo):  New argument. Argument "verbose" was added.
## 2020.06.23 (Leo):  New argument. Argument "tune" was added. See details.


rs_control <- function(retain_by = "proportion",
                       percentile_type = if (retain_by == "probability") 7,
                       tune = FALSE,
                       number = 10,
                       p = 0.75,
                       verbose = TRUE,
                       allow_parallel = TRUE) {
  # Sanity checks
  if (!retain_by %in% c("proportion", "probability")) {
    stop("Argument 'retain_by' must be either 'proportion' or 'probability'")
  }

  percentile_type
  if (!is.null(percentile_type)) {
    if (retain_by == "proportion") {
      warning("When argument retain_by == 'proportion', the aregument percentile_type is meaningless!")
    }
    if (!percentile_type %in% c(1:9)) {
      stop("Argument 'percentile_type' must be an integer between 1 and 9")
    }
  } else {
    percentile_type <- NULL
  }


  if (!is.logical(allow_parallel)) {
    stop("allow_parallel must be a logical value")
  }

  if (!is.numeric(number) | length(number) != 1) {
    stop("The 'number' argument must be a single numeric value")
  }

  if (!is.numeric(p) | length(p) != 1 | p >= 1 | p <= 0) {
    stop("p must be a single numeric value larger than 0 and below 1")
  }

  if (!is.logical(tune)) {
    stop("tune must be a logical value")
  }

  if (!is.logical(verbose)) {
    stop("The 'progress' argument must be logical")
  }

  cntrl <- list(
    retain_by = retain_by,
    percentile_type = percentile_type,
    tune = tune,
    number = number,
    p = p,
    verbose = verbose,
    allow_parallel = allow_parallel
  )

  cntrl
}
