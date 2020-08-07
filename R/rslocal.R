#' @title Data-driven search and optimization in spectral libraries for
#' building site-specific calibrations (RS-LOCAL)
#' @aliases rslocal
#' @aliases rslocal.default
#' @aliases rslocal.formula
#' @aliases predict.rslocal
#' @description
#'
#' \lifecycle{maturing}
#'
#' This function implements the re-sampling local (RS-LOCAL) algorithm. The
#' algorithm selects a subset from a large calibration/reference data optimized 
#' for deriving 'local' or site-specific calibrations. It uses a data-driven 
#' approach to capture the local or site-specific relationship between the 
#' response and predictor variables often not represented by large and complex 
#' calibration data sets (e.g. spectral libraries).
#'
#' @usage
#' \method{rslocal}{formula}(formula, train, test,
#'         k, b, method,
#'         ..., na_action = na.pass)
#'
#' \method{rslocal}{default}(Xr, Yr, Xu, Yu = NULL,
#'         k, b, retain = 0.95,
#'         method = local_fit_pls(pls_c = min(dim(Xr), 10)),
#'         optimization = "reconstruction",
#'         control = rs_control(),
#'         group = NULL,
#'         scale = FALSE,
#'         documentation = character(), ...)
#'
#' \method{predict}{rslocal}(object, newdata, type = 'response', ...)
#'
#' @param formula an object of class \link[stats]{formula} which represents the
#' basic model to be use.
#' @param train a data.frame with the training data containing the variables in
#' the model.
#' @param test a data.frame with the test data (local observations) containing
#' the variables in the model.
#' @param Xr a matrix of predictor variables of the reference data
#' (observations in rows and variables in columns).
#' @param Yr a matrix of one column containing the values of the
#' response variable corresponding to the reference data.
#' @param Xu a matrix of predictor variables of the 'local' or
#' site-specific observations (observations in rows and variables in columns).
#' @param Yu a matrix of one column containing the values of the response
#' variable corresponding to the Xu data. It is only mandatory if
#' \code{optimization = "response"}. Default is \code{NULL}.
#' @param k the number of reference/training observations selected in
#' the re-sampling step of the rs-local algorithm.
#' See description.
#' @param b the number of times each observation in the training set (\code{Xr})
#' should be tested, on average, in each iteration of the rs-local
#' algorithm. See description.
#' @param retain a numeric value larger than 0 and below 1 (default 0.95).
#' This value indicates:
#'  \itemize{
#'  \item{if the \code{...$retain_by} parameter of the object passed to the
#'  \code{control} argument is \code{'proportion'}:}{ a percentage of the
#'  observations to be kept at each iteration. The proportion of observations to be
#'  removed is \code{1 - retain}. See \code{retain_by} argument of the
#'  \code{\link{rs_control}} function.}
#'  \item{if the \code{...$retain_by} parameter of the object passed to the
#'  \code{control} argument is \code{'probability'}:}{ a probability value to be
#'  used to estimate the percentile (cut point of the distribution) of
#'  associated errors. In this case, observations with associated errors below this
#'  estimated percentile value are kept and the ones equal to or above this
#'  value are removed. See \code{retain_by} argument of the
#'  \code{\link{rs_control}} function.}
#' }
#' @param method an object of class \code{\link{local_fit}} which indicates the
#' type of regression to conduct within the rs-local alrgorithm as well as
#' additional parameters affecting this regression. See \code{\link{local_fit}}
#' function. The only method allowed for the moment is pls, i.e. methods created
#' with \code{local_fit_pls()}.
#' @param optimization a character string indicating the sample search method.
#' Options are:
#' \itemize{
#'  \item{\code{'response'}: }{The observations are retained based on the root
#'  mean squared error of the prediction of the response variable in the test
#'  set (\code{Yu}). In this case, it is required to pass the response values
#'  for the test set (\code{Yu}).}
#'  \item{\code{'reconstruction'} (Default): }{The observations are retained
#'  based on the spectral reconstruction error estimated for the test set
#'  (\code{Xu}). This error is estimated by projecting the test set onto the pls
#'  space (using the pls model built in each iteration with a fixed number of
#'  pls factors) and then back-transforming (reconstruct) the projected test set
#'  to its spectral space. Finally, the root mean squared error (RMSE) of this 
#'  reconstruction is computed and used to identify the relevant observations 
#'  to be retained. This reconstruction is done with a fixed number of pls 
#'  components for all the iterations so that reconstruction errors are 
#'  comparable. In this case, it response values for the test set are not 
#'  required.}
#'  }
#' @param control a list created with the \code{\link{rs_control}} function
#' which contains additional parameters that further control some aspects of the
#' \code{rs_local} function. The default list is as returned by
#' \code{rs_control()}. See the \code{\link{rs_control}} function for more
#' details.
#' @param scale a logical indicating if the predictor variables must be scaled
#' to unit variance at each iteration before regression.
#' @param group an optional factor (or vector that can be coerced
#' to \code{\link[base]{factor}} by \code{as.factor}) that assigns to each observation in
#' the training set (i.e. \code{train} or \code{Xr}) a group/class label (e.g.
#' groups can be given by spectra collected from the same batch of measurements,
#' from the same observation, from observations with very similar origin, etc).
#' This is taken into account for internal leave-group-out cross-validation for
#' pls tuning (factor optimization) to avoid pseudo-replication. When one
#' observation is selected for cross-validation, all observations of the same
#' group are removed together and assigned to validation. The length of the
#' vector must be equal to the number of observations in the training set
#' (i.e. \code{nrow(train)} or \code{nrow(Xr)}). See details.
#' @param documentation an optional character string that can be used to
#' describe anything related to the \code{mbl} call (e.g. description of the
#' input data). Default: \code{character()}. NOTE: this is an experimental
#' argument.
#' @param ... optional parameters used in method `formula` to be passed to the
#' low level function rslocal.default. Not currently used for `predict` and
#' `default` methods of `rslocal`.
#' @param na_action  a function to specify the action to be taken if NAs are
#' found in the data passed to train. Default \code{\link[stats]{na.fail}}.
#' NOTE: \code{na_action} is not applied to data passed to \code{test}. If given,
#' this argument must be named.
#'
#' @param object an object of class `rslocal`, as that created by the function
#' \code{\link{rslocal}}.
#' @param newdata a data frame or matrix containing new data.
#' @param type a character vector indicating what to return. Options are:
#' `response` (default) and `scores`.
#'
#' @details
#' The rs-local algorithm requires a (large) reference data set (\code{Xr}, 
#' where the sample search is conducted), a subset of specific observations (\code{Xu}),
#' and three parameters, \code{k}, \code{b} and \code{retain}, which are described
#' below. 
#' The \code{Xu} observations should be representative of the entire population
#' which they originated from. They may be selected, for example, from a large
#' set of observations (with unknown response values) using sampling method such
#' Kennard-Stone (Kennard & Stone, 1969).
#'
#' The rs-local algorithm uses \code{Xu} to iteratively remove irrelevant 
#' reference observations (from \code{Yr, Xr}) until a 
#' subset of approximately \code{k} reference observations is found. Irrelevant 
#' observations are those that tend to increase the root mean square error (RMSE)
#' of either the prediction of \code{Yu} or the pls reconstruction of \code{Xu}. 
#' In rs-local, a re-sampling method is employed where the reference observations 
#' that consistently fall in the the sampling subsets with the largest errors 
#' are identified and labelled as irrelevant and the rest are retained. The 
#' argument \code{retain} of the function controls the amount of samples to be 
#' retained.
#' 
#' Since the optimizations are based on the RMSEs, the subset found by 
#' \code{rslocal} is supposed to be the subset of approximately \code{k} 
#' reference observations that maximizes the accuracy of pls mdels for samples 
#' similar to the ones provided in the test set. 
#' 
#' When choosing the values for: 
#' \itemize{
#'  \item{\code{k:}}{ The size of the randomly sampled data set used in the
#'  internal calibration and validation step of the algorithm. It is also the
#'  target number of reference/training observations returned by the algorithm
#'  i.e. when insufficient  samples remain in the reference/training set to
#'  continue re-sampling. For recommended values see
#'  Lobsey et al. 2017.}
#'  \item{\code{b:}}{ The number of times each training observation is tested in each
#'  iteration of the algorithm (on average). More consistent results are
#'  achieved with high \code{b} values, however this will increase processing
#'  time. A recommended value for \code{b} is greater than 40.}
#'  \item{\code{retain:}}{ This determines the number of observations to be kept at
#'  each iteration of the algorithm. More consistent results are achieved with
#'  large \code{retain} values, however this will increase processing time. It
#'  is recommended to use \code{retain} values larger than 0.9.}
#'  }
#' @return a \code{list} with the following elements:
#' \itemize{
#'  \item{\code{x_local}:}{ a matrix of predictor variables corresponding to the
#'  observations selected.}
#'  \item{\code{y_local}:}{ a matrix of one column with the response variable
#'  corresponding to the observations selected.}
#'  \item{\code{indices}:}{ a numeric vector with the indices of the
#'  observations selected from the original training set.}
#'  \item{\code{iter_rmse}:}{ a data.table with the maximum associated error
#'  found for the observations retained at the end of each iteration
#'  (\code{max_rmse}) and the associated error above which the observations
#'  where removed at each iteration (column \code{cut_rmse}). If the
#'  observations were retained based on a given probability (specified
#'  in the \code{retain} argument), the \code{cut_rmse} column indicates the
#'  cut-off percentile value.}
#'  \item{\code{n_removed}:}{ a data.table with the number of observations
#'  removed at the end of each iteration.}
#'  \item{\code{validation_results}:}{ a list containing some validation results.
#'  This list has two elements:
#'  \itemize{
#'  \item{\code{val_info}:}{ a list containing the indices of the observations
#'  (in the original training set) in each re-sampling iteration used for the
#'  validation of the final model. A matrix with the predictions for the test set
#'  using the selected   observations is also returned.}
#'  \item{\code{results}:}{ a list containing the results of the
#'  leave-group-out validation done for the selected training observations and
#'  if \code{Yr} was supplied, the validation results for the test set.}
#'  }
#'  }
#'  \item{\code{control}:}{ a list mirroring the one provided in the
#'  \code{control} argument.}
#'  \item{\code{scale}:}{ was scaling used?.}
#'  \item{\code{final_model}:}{ a list with the following elements of a pls
#'  regression model:
#'  #'  \itemize{
#'  \item{\code{npls}:}{ The number of pls factors used.}
#'  \item{\code{coefficients}:}{ The regression coefficients for each factor.}
#'  \item{\code{bo}:}{ The intercepts of each factor.}
#'  \item{\code{scores}:}{ The matrix of pls scores.}
#'  \item{\code{X_loadings}:}{ The matrix of pls loadings for the predictor
#'  variables.}
#'  \item{\code{Y_loadings}:}{ The matrix of pls loadings for the response
#'  variable.}
#'  \item{\code{projection_mat}:}{ The matrix for pls projections.}
#'  \item{\code{vip}:}{ The matrix of variable importance for projection of
#'  each factor.}
#'  \item{\code{selectivity_ratio}:}{ The matrix of variable importance for each
#'  factor based on the method of selectivity ratio (Rajalahti et al., 2009).}
#'  \item{\code{Y}:}{ A matrix with the response values used to fit the final
#'  pls model.}
#'  \item{\code{weights}:}{ A matrix of pls weights.}
#'  }
#'  }
#'  \item{\code{documentation}:}{ a character string mirroring the one provided
#'  in the \code{documentation} argument.}
#'  }
#' @importFrom stats quantile complete.cases diffinv na.pass
#' @author Craig Lobsey, Raphael Viscarra Rossel and Leonardo Ramirez-Lopez
#' @references
#' Lobsey, C. R., Viscarra Rossel, R. A., Roudier, P., & Hedley, C. B. 2017.
#' rs-local data-mines information from spectral libraries to improve local
#' calibrations. European Journal of Soil Science, 68(6), 840-852.
#'
#' Kennard, R.W. & Stone, L.A. 1969. Computer aided design of experiments.
#' Technometrics, 11(1), pp.137-148.
#'
#' Rajalahti, T., Arneberg, R., Berven, F. S., Myhr, K. M., Ulvik, R. J.,
#' Kvalheim, O. M. 2009. Biomarker discovery in mass spectral profiles by means
#' of selectivity ratio plot. Chemometrics and Intelligent Laboratory Systems,
#' 95(1), 35-48.
#' @examples
#' \dontrun{
#' # FIXME: GOOD EXAMPLES ARE REQUIRED
#' }
#' @export rslocal


## 2020.03.28 (Leo):    Bug fix. Optimization was wrongly set to "reconstruction"
##                      in biter(), even when "response" was selected.
## 2020.03.28 (Leo):    New output. "iter_rmse" is a data.frame containing the
##                      maximum rmse obtained at each iteration after removing
##                      the observations.
## 2020.03.29 (Leo):    Some secondary arguments moved to the new rs_control
##                      function.
## 2020.03.29 (Craig):  Bug fix. In the original code the final
##                      return(list(K.x = SL.x[k_idx,], K.y=SL.y[k_idx],
##                      k_idx=k_idx))
##                      should have been using sl_idx for indexing, not k_idx as
##                      k_idx is not updated until the begining of the next
##                      iteration.
## 2020.03.30 (Leo):    Bug fix. NAs were not properly handled to produce the
##                      final validation stats (tables with NA values were
##                      generated).
##                      Added a sanity check for missing values when
##                      optimization == "response"
## 2020.04.05 (Craig):  Input checking - changed some messages for consistency
##                      and changed check nrow(Xu)!=nrow(Yu) to be performed
##                      when optimise is response only
## 2020.04.05 (Craig):  Changed type conversion of input data variables to
##                      matrix, only convert Yu if != NULL
## 2020.04.05 (Craig):  Fixed complete.cases(Yu) to handle Yu=NULL with verbose
## 2020.06.23 (Leo):    Argument "pls_tune" was removed and passed to rs_control
##                      as tune
##                      method must be "local_fit" object


"rslocal" <-
  function(...) {
    UseMethod("rslocal")
  }


#' @aliases rslocal
#' @export
#'
## add an argument to initualize with a given amount of observations, e.g. if the library is 30.000
## you could initialize with 15.000 observations and to compensate you can increase the number of resampling iterations
## perhaps a number of initial resampling iterations can also work say we start with 10.000 for the first loop and then it goes to what b specifies
rslocal.default <- function(Xr,
                            Yr,
                            Xu,
                            Yu = NULL,
                            k,
                            b,
                            retain = 0.95,
                            method = local_fit_pls(pls_c = min(dim(Xr), 10)),
                            optimization = "reconstruction",
                            control = rs_control(),
                            group = NULL,
                            scale = FALSE,
                            documentation = character(),
                            ...) {

  # check inputs
  if (nrow(Xr) != length(Yr)) {
    stop("The number of spectra in Xr must equal to the length of Yr")
  }

  if (!is.null(Yu)) {
    if (nrow(Xu) != length(Yu)) {
      stop("The number of spectra in Xu must equal to the length of Yu")
    }
  }
  
  if (optimization == "response") {
    if (anyNA(Yu) | anyNA(Xu)) {
      stop("Missing values are not allowed in the response variables when optimization == 'response'")
    }
  }

  if (any(is.infinite(Xr)) | any(is.infinite(Yr)) | any(is.infinite(Xu)) | any(is.infinite(Yu))) {
    stop("Infinite values detected in the input data")
  }

  if (ncol(Xr) != ncol(Xu)) {
    stop("The number of variables in Xr must equal the numnber of variables in Xu")
  }

  k <- round(k)
  if (k >= nrow(Xr)) {
    stop("Argument 'k' must be an integer lower than the number of observations in the train set (Xr)")
  }


  if (retain > 1 | retain <= 0) {
    stop("Argument 'retain' must be a numerical value larger than 0 and below 1")
  }

  if (optimization == "response" & is.null(Yu)) {
    stop("When optimization = 'response', Yu values must be provided")
  }

  if (!optimization %in% c("response", "reconstruction")) {
    stop("Argument 'optimization'  must be either 'response' or 'reconstruction'")
  }

  if (!"local_fit" %in% class(method)) {
    stop("Method must be of class 'local_fit'")
  }

  if (method$method != "pls") {
    stop("The only method alowed for the moment is 'pls' generated with the 'local_fit_pls()' fucntion")
  }

  if (optimization == "reconstruction" & control$tune & method$method %in% c("pls", "wapls")) {
    warning("pls factors are not tuned when optimization = 'reconstruction', instead they are fixed to the one(s) provided in the `method` argument, therefore the `tune` option passed to control has been ignored.")
  }

  if (control$tune) {
    min_samples <- floor(min(k, dim(Xr)) * control$p) - 1
    min_cv_samples <- floor(min(k, dim(Xr)) * (1 - control$p))
    if (min_cv_samples < 3) {
      stop(paste0(
        "Local cross-validation requires at least 3 observations in ",
        "the hold-out set, the current cross-validation parameters ",
        "leave less than 3 observations for one or more resampling ",
        "iterations."
      ))
    }
  } else {
    min_samples <- floor(min(k, dim(Xr))) - 1
  }

  if (method$method %in% c("pls", "wapls")) {
    max_pls <- max(method$pls_c)
    if (any(min_samples < max_pls)) {
      stop(paste0(
        "More pls components than observations for one or more \n",
        "resampling iterations. If 'tuning' is being used, consider that some ",
        "observations \nin the resampling iterations are hold-out for ",
        "validation"
      ))
    }
  }


  call_f <- (match.call())

  Xr <- as.matrix(Xr)
  Yr <- as.matrix(Yr)
  Xu <- as.matrix(Xu)
  if (!is.null(Yu)) {
    Yu <- as.matrix(Yu)
  }

  r <- 1 - retain

  # This list tracks the current selected SL subset (K) at each iteration
  k_list <- list()

  # sl_idx is a vector representing an sample index into the SL
  # note in this implementation we operate with an index into the SL and not copies of the SL spectra
  sl_idx <- seq(1, nrow(Xr))

  #
  # Step 1 - Initialise K as a subset of the SL, initially full.
  #
  k_idx <- sl_idx

  #
  # This while loop contains step 2-6 of the algorithm (See step 7!).
  #

  # use this to keep track of the number of iterations
  outer_idx <- 0

  # initialise a vector for the sample rankings (by rmse) that corresponds with the full SL
  # also initialise a vector to track the number of times each sample was tested in the 'B' iteration
  U <- V <- rep(0, length(sl_idx))
  names(U) <- names(V) <- sl_idx
  s_to_drop <- cuts <- maxrmse <- NULL
  pp <- r
  cat_progress <- cat_iter(c("\\", "|", "/", "-"))
  cat_progress2 <- cat_iter(c("\\", "|", "/", "-"))
  

  while (length(k_idx) > (k * (1 + pp)) & sum(!is.na(sl_idx)) > k) {
    #
    # Step 1 - Initialise K as a subset of the SL, initially full. k_idx only
    # contains those SL samples still in K
    #
    # select k_idx as those in the SL not marked as dropped (zero)

    k_idx <- sl_idx[!is.na(sl_idx)]

    outer_idx <- outer_idx + 1

    # calculate the quantity of samples removed in this iteration
    cull_quantity <- round(length(k_idx) * r)

    ##  This sub-loop contains step 2-4 of the algorithm (See step 5!). This is
    ##  the B iteration.
    B <- round(length(k_idx) * b / k)
    if (control$verbose) {
      cat(paste0("Iteration ", outer_idx, ":\t"))
      if (outer_idx == 1) {
        cat(paste0("\033[34m", cat_progress(), " Initial size: ", length(k_idx), "\t\033[39m"))
      } else {
        cat(paste0("\033[34m", cat_progress(), " Current subset size: ", length(k_idx), "\t\033[39m"))
      }
      cat(paste0("\033[34m", cat_progress2(), " Resampling iterations to perform: ", B, "\t\r\033[39m"))
    }

    ## Step 2 - sample a training data set of size k from K without replacement
    ismpl <- replicate(
      n = B,
      expr = sample(x = k_idx, size = k, replace = FALSE)
    )
    ## this iterator object will avoid to use much memory in the next  foreach
    ## loop this iterator takes from the big matrix only what is needed and
    ## therefore the whole matrix is not put in the memory.
    ## This makes parallel computations more memory friendly
    itersubs <- ithrssubsets(
      x = Xr,
      y = Yr,
      group = group,
      indx = ismpl
    )

    ## Step 3 calibrate a PLS model using the selected SL samples
    ## Step 4 - validate on the 'm' site specific samples
    ## Step 5 - repeat 2 to 4 B times
    results_df <- biter(
      itersubs = itersubs,
      Xu = Xu,
      Yu = Yu,
      iter_sequence = seq(1, B),
      optimization = optimization,
      ncomp = method$pls_c,
      tune = control$tune,
      p = control$p,
      number = control$number,
      scale = scale,
      max_iter = 1,
      tol = 1e-6,
      allow_parallel = control$allow_parallel
    )

    ## Step 6 - starts here

    # iterate through all B iteration results and increment the rankings (rmse)
    # and test count for each sample
    #

    for (i in 1:nrow(results_df$sampleidx)) {
      U[results_df$sampleidx[i, ]] <- U[results_df$sampleidx[i, ]] + results_df$rmsesubset[i] # the index of samples used starts at 4
      V[results_df$sampleidx[i, ]] <- V[results_df$sampleidx[i, ]] + 1
    }
    V[V == 0] <- NA

    # normalise U using V
    U <- U / V
    u_df <- data.frame(
      idx = sl_idx,
      rmse = U
    )

    # now order / rank by rmse
    U_ordered <- u_df[order(u_df$rmse, decreasing = TRUE), ]
    U_ordered <- U_ordered[!is.na(U_ordered$rmse), ]

    if (control$retain_by == "proportion") {
      # select the poorest performing SL samples, the amount is determing by
      # cull_quantity calculated earlier
      worst_reference_set <- U_ordered[1:cull_quantity, ]
      ok_reference_set <- U_ordered[-c(1:cull_quantity), ]
      cutoff <- min(worst_reference_set$rmse)
      pp <- 1 - nrow(ok_reference_set) / nrow(U_ordered)
    } else {
      ## Instead removing a given number of samples
      ## remove samples that are avobe a given cutoff prob
      cutoff <- quantile(U, 1 - r, na.rm = TRUE, type = control$percentile_type)
      worst_reference_set <- U_ordered[U_ordered$rmse >= cutoff, ]
      ok_reference_set <- U_ordered[U_ordered$rmse < cutoff, ]
      if (nrow(ok_reference_set) < max(method$pls_c)) {
        message(paste0(
          "Iteration interrupted as the number of selected observations (",
          nrow(ok_reference_set),
          ") is lower than the number of pls factors (", max(method$pls_c), ")"
        ))
        break
      }
      pp <- 1 - nrow(ok_reference_set) / nrow(U_ordered)
    }

    s_to_drop <- c(s_to_drop, nrow(worst_reference_set))
    cuts <- c(cuts, cutoff)
    maxrmse <- c(maxrmse, max(ok_reference_set$rmse))

    # # provide some status information to the user
    # if(verbose)
    #   cat("- Train samples to drop: ", nrow(worst_reference_set), "\n")

    # if(outer_idx  == 1){
    #   rmseprobs <- quantile(U_ordered$rmse,  probs = seq(0, 1, 0.25))
    # }
    #

    # rmseprobs <- rbind(rmseprobs, quantile(ok_reference_set$rmse,  probs = seq(0, 1, 0.25), na.rm = TRUE))

    # now set the latest dropped samples in sl_idx as NA to remove them from consideration
    sl_idx[worst_reference_set$idx] <- NA

    U[] <- 0
    V[] <- 0
  }

  # Finalise k_idx following the last iteration
  k_idx <- sl_idx[!is.na(sl_idx)]

  if (control$verbose) {
    cat("\nFitting final model on", length(k_idx), "selected training samples... \n")
    if (!is.null(Yu)) {
      if (sum(complete.cases(Yu)) != 0) {
        cat("and", sum(complete.cases(Yu)), "test samples for which response values were available... \n")
      }
    }
  }

  ## TEMPORARY (THIS MUST GO WHEN A PROPER PLS FUNCTION IS IMPLEMENTED)
  plsval <- pls_cv(
    x = as.matrix(Xr[k_idx, ]),
    y = as.matrix(Yr[k_idx, ]),
    ncomp = method$pls_c,
    method = method$method,
    center = TRUE,
    scale = scale,
    min_component = 1,
    p = control$p,
    number = control$number,
    group = group[k_idx],
    retrieve = FALSE,
    tune = TRUE,
    max_iter = 1,
    tol = 1e-6
  )

  ## TEMPORARY (THIS MUST GO WHEN A PROPER PLS FUNCTION IS IMPLEMENTED)
  plsval$models <- opls_get_basics(
    X = as.matrix(Xr[k_idx, ]),
    Y = as.matrix(Yr[k_idx, ]),
    ncomp = method$pls_c,
    scale = scale,
    maxiter = 1,
    tol = 1e-6
  )

  rnms <- rownames(plsval$resamples)
  plsval$resamples <- apply(
    X = plsval$resamples,
    MARGIN = 2,
    FUN = function(r, i) i[r], i = k_idx
  )
  rownames(plsval$resamples) <- rnms

  validation <- list(
    val_info = list(train.resamples = plsval$resamples),
    results = list(train = plsval$cv_results)
  )

  if (!is.null(Yu)) {
    ## TEMPORARY (THIS MUST GO WHEN A PROPER PLS FUNCTION IS IMPLEMENTED)
    yuhat <- predict_opls(
      bo = plsval$models$bo,
      b = plsval$models$coefficients,
      newdata = Xu,
      ncomp = method$pls_c,
      scale = scale,
      Xscale = plsval$models$transf$Xscale
    )

    colnames(yuhat) <- paste("pls", 1:ncol(yuhat), sep = "")

    if (is.null(rownames(Yu))) {
      rownames(yuhat) <- 1:nrow(yuhat)
    } else {
      rownames(yuhat) <- rownames(Yu)
    }


    validation$val_info$test_predictions <- yuhat

    if (sum(is.na(Yu)) < length(Yu)) {
      testval <- data.table(
        npls = 1:method$pls_c,
        rmse = sqrt(colMeans(sweep(yuhat[complete.cases(Yu), ],
          MARGIN = 1,
          FUN = "-",
          STATS = Yu[complete.cases(Yu)]
        )^2)),
        r2 = as.vector(cor(yuhat[complete.cases(Yu), ], Yu[complete.cases(Yu)])^2),
        row.names = NULL
      )
      validation$results$test <- testval
    }

    ## TEMPORARY (THIS MUST GO WHEN A PROPER PLS FUNCTION IS IMPLEMENTED)
    tmpx <- as.matrix(rbind(Xr[k_idx, ], Xu))
    tmpy <- as.matrix(c(Yr[k_idx, ], Yu))
    tmpx <- tmpx[complete.cases(tmpy), ]
    tmpy <- tmpy[complete.cases(tmpy), , drop = FALSE]

  } else {
    tmpx <- as.matrix(Xr[k_idx, ])
    tmpy <- as.matrix(Yr[k_idx, ])
  }

  finalpls <- opls_get_all(
    X = tmpx,
    Y = tmpy,
    ncomp = method$pls_c,
    scale = scale,
    maxiter = 1,
    tol = 1e-6
  )  
  
    ncompnms <- names(finalpls) %in% "ncomp"
  if (any(ncompnms)) {
    names(finalpls)[names(finalpls) %in% "ncomp"] <- "npls"
  }

  ## TEMPORARY (THIS MUST GO WHEN A PROPER PLS FUNCTION IS IMPLEMENTED)
  finalpls$coefficients <- t(finalpls$coefficients)
  finalpls$projection_mat <- t(finalpls$projection_mat)
  finalpls$vip <- t(finalpls$vip)
  finalpls$selectivity_ratio <- t(finalpls$selectivity_ratio)

  colnames(finalpls$coefficients) <-
    colnames(finalpls$X_loadings) <-
    colnames(finalpls$projection_mat) <-
    colnames(finalpls$vip) <-
    colnames(finalpls$selectivity_ratio) <-
    colnames(finalpls$weights) <- colnames(Xr)
  rownames(finalpls$coefficients) <-
    colnames(finalpls$bo) <-
    rownames(finalpls$X_loadings) <-
    rownames(finalpls$Y_loadings) <-
    rownames(finalpls$projection_mat) <-
    rownames(finalpls$vip) <-
    rownames(finalpls$selectivity_ratio) <-
    rownames(finalpls$weights) <- 1:method$pls_c
  rownames(finalpls$bo) <- "Intercepts"
  colnames(finalpls$Y_loadings) <- "Loadings"


  colnames(finalpls$scores) <- 1:method$pls_c
  rownames(finalpls$scores) <- 1:nrow(finalpls$scores)
  
  resultsList <- list(
    x_local = Xr[k_idx, ],
    y_local = Yr[k_idx],
    indices = k_idx,
    iter_rmse = data.table(
      iteration = 1:length(maxrmse), max_rmse = maxrmse,
      cut_rmse = cuts
    ),
    n_removed = data.table(
      iteration = 1:length(s_to_drop), removed = s_to_drop,
      cummulative = diffinv(s_to_drop)[-1]
    ),
    validation_results = validation,
    control = control,
    scale = scale,
    final_model = finalpls,
    documentation = documentation
  )

  attr(resultsList, "call") <- call_f
  class(resultsList) <- c("rslocal", "list")

  resultsList
}


#' @aliases rslocal
#' @export
rslocal.formula <- function(formula,
                            train,
                            test,
                            k,
                            b,
                            method,
                            ...,
                            na_action = na.pass) {
  if (!inherits(formula, "formula")) {
    stop("'formula' is only for formula objects")
  }

  call_f <- (match.call())

  if (missing(method)) {
    stop("'method' is missing")
  }

  ## Get the model frame
  mf <- match.call(expand.dots = FALSE)

  mr <- match(x = c("formula", "train", "na_action"), table = names(mf))
  mu <- match(x = c("formula", "test"), table = names(mf))

  mfr <- mf[c(1, mr)]
  mfu <- mf[c(1, mu)]

  names(mfr)[names(mfr) %in% "na_action"] <- "na.action"
  names(mfr)[names(mfr) %in% "train"] <- "data"
  names(mfu)[names(mfu) %in% "test"] <- "data"

  yname <- all.vars(formula, functions = FALSE, max.names = 1)

  mfr[[1]] <- mfu[[1]] <- as.name("model.frame")

  input_list <- list(...)

  if (!yname %in% colnames(eval(mfu$data))) {
    if ("optimization" %in% names(input_list)) {
      if (input_list$optimization == "response") {
        stop("When optimization = 'response', response values must be provided in test")
      }
    }

    test <- test %>%
      mutate(y = NA) %>%
      rename(y = yname)
    warning(paste(
      yname, "not found in test. Missing values (NAs) were assigned to",
      yname, "in test."
    ))
  }

  mfu <- model.frame(mfu, data = test, na.action = NULL)
  mfr <- eval(mfr, parent.frame())

  trms <- attr(mfr, "terms")

  formulaclasses <- attr(trms, "dataClasses")

  attr(trms, "intercept") <- 0
  xr <- model.matrix(trms, model.frame(mfr, drop.unused.levels = T))
  yr <- model.extract(mfr, "response")

  xu <- model.matrix(trms, mfu)
  yu <- model.extract(mfu, "response")


  if ("optimization" %in% names(input_list) & sum(is.na(yu)) == length(yu)) {
    stop("When optimization = 'response', response values must be provided in test")
  }

  rsl <- rslocal(
    Xr = xr,
    Yr = yr,
    Xu = xu,
    Yu = yu,
    k = k,
    b = b,
    method = method,
    ...
  )
  rsl$formula <- formula
  rsl$dataclasses <- formulaclasses

  attr(rsl, "call") <- call_f
  class(rsl) <- c("rslocal", "rslocal.formula", "list")

  rsl
}


#' @aliases rslocal
#' @importFrom stats .MFclass terms delete.response
#' @export
predict.rslocal <- function(object, newdata, type = "response", ...) {
  
  if (missing(newdata)) {
    stop("newdata is missing")
  }

  if (!is.null(object$formula)) {
    dcls <- object$dataclasses[-1]

    if (!"matrix" %in% class(newdata) | !"data.frame" %in% class(newdata)) {
      stop(paste0(
        "When predicting from objects of class 'rslocal' fitted with ",
        "formula, the argument 'newdata' must be a 'data.frame' ",
        "or alternatively a 'matrix'"
      ))
    }
    
    if (class(newdata) == "data.frame") {
      if (!all(names(dcls) %in% names(newdata))) {
        mss <- names(dcls)[!names(dcls) %in% names(newdata)]
        stop(paste(
          "The following predictor variables are missing:",
          paste(mss, collapse = ", ")
        ))
      }
    }

    if (any(dcls != "numeric")) {
      if (class(newdata) == "matrix" & length(dcls) == 1) {
        if (.MFclass(newdata) == dcls) {
          pnames <- gsub(
            names(dcls), "",
            colnames(object$final_model$coefficients)
          )
          if (all(pnames %in% colnames(newdata))) {
            newdata_temp <- newdata
            newdata <- data.frame(rep(NA, nrow(newdata)))
            colnames(newdata) <- names(object$dataclasses)[1]
            newdata[[names(dcls)]] <- newdata_temp[, pnames]
          } else {
            stop("Missing predictor variables")
          }
        }
      }
    }

    oterms <- terms(object$formula)
    oterms <- delete.response(oterms)
    attr(oterms, "intercept") <- 0
    mf <- model.frame(oterms, newdata)
    newdata <- model.matrix(oterms, model.frame(mf, drop.unused.levels = T))
  }

  if (!"matrix" %in% class(newdata)) {
    stop("Argument 'newdata' must be a 'matrix'")
  }

  # if(all(pnames %in% colnames(newdata)))

  if (type == "response") {
    preds <- predict_opls(
      bo = object$final_model$bo,
      b = t(object$final_model$coefficients),
      ncomp = object$final_model$npls,
      newdata = as.matrix(newdata),
      scale = object$scale,
      Xscale = object$final_model$transf$Xscale
    )
  }

  if (type == "scores") {
    stop("not yet implemented") # FIXME!!
  }
  rownames(preds) <- rownames(newdata)
  colnames(preds) <- 1:ncol(preds)

  preds
}
