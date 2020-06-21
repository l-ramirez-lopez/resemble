#' @title A function that controls some aspects of the the \code{rslocal} function
#' @description
#' This function is used to further control some aspects of the \code{rslocal} function.
#' @usage 
#' rs_control(keep_by = "proportion",
#'            percentile_type = if(keep_by == "probability") 7,
#'            number = 10,
#'            p = 0.75,
#'            scaled = FALSE,
#'            pls_max_iter = 1, 
#'            pls_tol = 1e-6,
#'            verbose = TRUE,
#'            allow_parallel = TRUE)
#' @param keep_by a character string indicating how the training samples to be kept at each iteration must be selected. Two options are available: \code{'proportion'} (default) retains a fix proportion of samples at each iteration and \code{'probability'} retains samples whose errors are below an estimated percentile (cut point) from a given probability. See details.
#' @param percentile_type if \code{keep_by = 'probability'}, an integer between 1 and 9 to be passed to the the \code{type} argument of the \code{\link[stats]{quantile}} function to estimate the percentile cutoff value used to select the samples at each iteration. Default is 7 (as in \code{\link[stats]{quantile}}). If \code{keep_by != 'proportion'}, this argument is ignored. See details.
#' @param number an integer indicating the number of groups for the leave-group-out internal cross validation used for pls factor tuning (optimization) at each iteration. These groups are built from the training subset selected at the respective iteration. Default is 10.
#' @param p a value indicating the percentage of samples to build each group in the the leave-group-out internal cross validation (used for pls factor tuning (optimization) at each iteration). These groups are built from the training subset selected at the respective iteration. Default is 0.75 (i.e. 75 "\%").
#' @param scaled (BETTER DESCRIPTION REQUIRED) a logical indicating whether or not the predictor variables must be scaled at each iteration (before regression).
#' @param pls_max_iter (BETTER DESCRIPTION REQUIRED) maximum number of iterations for the partial least squares methods.
#' @param pls_tol (BETTER DESCRIPTION REQUIRED) for convergence in the partial orthogonal scores partial least squares regressions using the nipals algorithm. Default is 1e-6
#' @param number (BETTER DESCRIPTION REQUIRED) an integer indicating the number of resampling iterations for cross validation. Default is 10.
#' @param verbose a logical indicating if some information about each iteration must be printed. 
#' @param allow_parallel set to TRUE to parallelise the internal resampling and calibration. The parallel backend should be registered by the user.         
#' @details
#' The training samples to be kept at each iteration can be selected by using a fixed proportion of samples (in relation to the number of samples in the training subset being used at each iteration) whose associated errors are the lowest, in this case the argument \code{keep_by} must be set to \code{'proportion'}. This proportion (a value larger than 0 and below 1) is given in the in the argument \code{.keep} of the \code{\link{rs_local}} function. 
#' 
#' For selecting the training samples, it is also possible to retain the ones whose associated errors are below a percentile estimated for a given probability (the percentile is estimated with the \code{\link[stats]{quantile}} function). In this case the probability (a value larger than 0 and below 1) must be provided in the argument \code{.keep}, this numeric value is passed to the argument \code{probs} in the \code{\link[stats]{quantile}} function. This method might be useful in the presence of outlier samples with considerable deviation in their associated errors.
#'
#' The internal validations for tuning the number of pls factors are based on the leave-group-out cross-validation. Arguments \code{p} and \code{number} are used to control the number of groups and the amount of samples per group.
#' @return a \code{list} mirroring the specified parameters.
#' @author Leonardo Ramirez-Lopez
#' @references 
#' Lobsey C., Viscarra Rossel R.A., Pierre R., Hedley C.B. 2017. RS-LOCAL data-mines information from large spectral libraries to improve local calibrations. European Journal of Soil Science.
#' @seealso \code{\link{rslocal}},  \code{\link[stats]{quantile}}
#' @examples
#' # A control list with the default parameters
#' rs_control()
#' 
#' # A control list which specifies that samples must be 
#' # retained based on the associated errors that are below
#' # a given probability.
#' rs_control(keep_by = "probability")
#' @export

## History:
## 2020.03.29 (Leo):  New fucntion to control rslocal
## 2020.03.29 (Leo):  New argument. Argument "verbose" was added. 
## 2020.03.29 (Leo):  New argument. Argument "" was added. See details.


rs_control <- function(keep_by = "proportion",
                       percentile_type = if(keep_by == "probability") 7,
                       number = 10,
                       p = 0.75,
                       scaled = FALSE,
                       pls_max_iter = 1, 
                       pls_tol = 1e-6,
                       verbose = TRUE,
                       allow_parallel = TRUE){
  # Sanity checks
  if(!keep_by %in% c("proportion", "probability"))
    stop("Argument 'keep_by' must be either 'proportion' or 'probability'")
  
  percentile_type
  if(!is.null(percentile_type)){
    if(keep_by == "proportion")
      warning("When argument keep_by == 'proportion', the aregument percentile_type is meaningless!")
    if(!percentile_type %in% c(1:9))
      stop("Argument 'percentile_type' must be an integer between 1 and 9")
  }else{
    percentile_type <- NULL
  }
  
  
  if(!is.logical(allow_parallel))
    stop("allow_parallel must be a logical value")
  
  if(!is.logical(scaled))
    stop("'scaled' argument must be logical")
  
  if(!is.numeric(number) | length(number) != 1)
    stop("The 'number' argument must be a single numeric value")
  
  if(!is.numeric(p) | length(p) != 1 | p >= 1 | p <= 0)
    stop("p must be a single numeric value larger than 0 and below 1")
  
  if(!is.logical(verbose))
    stop("The 'progress' argument must be logical")
  
  cntrl <- list(keep_by = keep_by,
                percentile_type = percentile_type,
                number = number,
                p = p,
                scaled = scaled,
                pls_max_iter = pls_max_iter, 
                pls_tol = pls_tol,
                verbose = verbose,
                allow_parallel = allow_parallel)    
  return(cntrl)
}
