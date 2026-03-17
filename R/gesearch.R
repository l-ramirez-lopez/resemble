#' @title Evolutionary training sample search in spectral libraries
#' for building context-specific models 
#'
#' @aliases gesearch
#' @aliases gesearch.default
#' @aliases gesearch.formula
#' @aliases gesearch.list
#' @aliases predict.gesearch
#'
#' @description
#' This function implements an evolutionary search algorithm that selects a 
#' subset from large calibration or reference data,
#' optimized to build context-specific ('local') calibrations. It applies a
#' data-driven approach to identify relationships between the response and
#' predictor variables that are often underrepresented in large calibration
#' datasets (e.g. spectral libraries).
#'
#' @usage
#' \method{gesearch}{formula}(formula, train, test, k, b, target_size,
#'         method, ..., na_action = na.pass)
#'
#' \method{gesearch}{default}(Xr, Yr, Xu, Yu = NULL, Yu_lims = NULL,
#'         k, b, retain = 0.95, target_size = k,
#'         method = local_fit_pls(pls_c = min(dim(Xr), 10)),
#'         optimization = "reconstruction",
#'         control = search_control(), group = NULL, scale = FALSE,
#'         verbose = TRUE, seed = NULL, documentation = character(),
#'         crossover = TRUE, intermediate_models = FALSE,
#'         pchunks = 1, ...)
#'
#' \method{predict}{gesearch}(object, newdata, type = 'response', ...)
#'
#' @param formula an object of class \link[stats]{formula} defining the model. 
#' Can also be a list of formulas with identical right-hand-side terms.
#' @param train a data.frame with training data, including model variables.
#' @param test a data.frame with test (local) data, including model variables.
#' @param Xr a matrix of predictor variables for the reference data 
#' (rows = observations, columns = variables).
#' @param Yr a matrix with one or more columns of response values matching Xr.
#' @param Xu a matrix of predictor variables for 'local' observations 
#' (same structure as Xr).
#' @param Yu a matrix with the same number of columns as \code{Yr} containing 
#' response values for \code{Xu}. Only required when 
#' \code{optimization = "response"}. Default: \code{NULL}.
#' @param Yu_lims numeric range indicating expected limits of the response 
#' variable in the target population. Only used when \code{Yu} is NULL or has 
#' one column.
#' @param k integer; number of samples in each sampled subset (gene).
#' @param b integer; target average number of times each training sample is 
#' included in the population per iteration. See \code{Details}.
#' @param retain numeric in (0, 1); proportion of samples retained per iteration. 
#' Default is 0.95. See \code{search_control()} for retention strategy.
#' \itemize{
#'   \item{If \code{control$retain_by == 'proportion'}: }{a fixed proportion 
#'   \code{retain} is used to retain samples. The rest (\code{1 - retain}) 
#'   are discarded. See \code{\link{search_control}}.}
#'   \item{If \code{control$retain_by == 'probability'}: }{the value is used 
#'   as a percentile cut-off for sample errors. Observations with errors 
#'   below this percentile are retained.}
#' }
#'
#' @param target_size integer; the target number of selected training samples 
#' (gene pool size). Must be ≥ \code{k}. Defaults to \code{k}.
#'
#' @param method an object of class \code{\link{local_fit}} specifying the 
#' regression model used during the search. Only methods created with 
#' \code{local_fit_pls()} are currently supported. See \code{\link{local_fit}}.
#' @param optimization character; specifies the optimization criterion used to 
#' guide the sample search. Options include:
#' \itemize{
#'   \item{\code{"response"}: }{Retains observations based on the root mean 
#'   squared error (RMSE) of predicting the response variable in the test set 
#'   (\code{Yu}). Requires \code{Yu} to be provided.}
#'
#'   \item{\code{"reconstruction"} (default): }{Retains observations based on 
#'   the spectral reconstruction error of the test set (\code{Xu}). The test 
#'   data are projected onto the PLS space and then reconstructed using a 
#'   fixed number of PLS components. The RMSE between the original and 
#'   reconstructed spectra is used as the selection criterion. Does not require 
#'   \code{Yu}.}
#'
#'   \item{\code{"similarity"}: }{Retains observations based on the mean 
#'   spectral similarity between test samples (\code{Xu}) and the training 
#'   samples used in each iteration. Test samples are projected into PLS score 
#'   space, and Mahalanobis distances are calculated between test and training 
#'   observations. The average distance is used to determine similarity. 
#'   \code{Yu} is not required.}
#'
#'   \item{\code{"range"}: }{Removes observations that yield predictions outside 
#'   a specified response range defined by \code{Yu_lims}. Used to filter out 
#'   extreme or unrealistic predictions.}
#' }
#'
#' Multiple methods may be combined by supplying a character vector, e.g., 
#' \code{c("reconstruction", "similarity")}.
#' @param control a list created with \code{\link{search_control}} that defines 
#' additional parameters for controlling the \emph{gesearch} algorithm. Defaults to 
#' \code{search_control()}. See \code{\link{search_control}} for details.
#'
#' @param scale logical; if \code{TRUE}, predictor variables are scaled to unit 
#' variance before regression at each iteration.
#'
#' @param group an optional factor or a vector coercible to a factor (via 
#' \code{as.factor}) that assigns a group label to each training observation 
#' (in \code{Xr} or \code{train}). This is used to define cross-validation 
#' folds for PLS tuning using leave-group-out validation to avoid 
#' pseudo-replication. All observations in the same group are removed together 
#' for validation. The length of this vector must match the number of training 
#' samples (i.e., \code{nrow(Xr)} or \code{nrow(train)}).
#'
#' @param verbose logical; if \code{TRUE}, prints progress information at each 
#' iteration.
#'
#' @param seed integer; seed value for random number generation to ensure 
#' reproducibility during cross-validation and sampling. Default is 
#' \code{NULL}, meaning no seed is set.
#'
#' @param documentation character; optional string used to record metadata or 
#' comments about the \code{gesearch} call. Default is an empty string. This 
#' parameter is experimental and may be subject to change.
#'
#' @param pchunks integer; number of chunks to split each iteration into during 
#' resampling. Useful for memory-efficient parallel processing of large 
#' datasets. Default is \code{1}. Increase this value for better memory 
#' management in parallel computing environments.
#
#' @param intermediate_models store models for each intermediate generation 
#' before reaching the last one.
#' 
#' @param ... additional arguments passed to \code{gesearch.default} when 
#' calling the formula interface. Not used in \code{default} or 
#' \code{predict} methods.
#'
#' @param na_action function; defines how missing values in the training data 
#' are handled. Default is \code{\link[stats]{na.pass}}. Note: this only 
#' applies to \code{train}; no NA handling is applied to \code{test}. If used, 
#' this argument must be named.
#'
#' @param object an object of class \code{gesearch}, as returned by 
#' \code{\link{gesearch}}.
#'
#' @param newdata a data.frame or matrix containing new predictor data for 
#' prediction.
#'
#' @param type character; specifies the prediction output type. Options are 
#' \code{"response"} (default), which returns predicted values, and 
#' \code{"scores"}, which returns the PLS scores from the model.
#'
#' @details
#' The \emph{gesearch} algorithm requires a large reference dataset (\code{Xr}) where 
#' the sample search is conducted, a subset of site-specific or target 
#' observations (\code{Xu}), and three tuning parameters: \code{k}, \code{b}, 
#' and \code{retain}.
#'
#' The \code{Xu} observations should represent the target population. These may 
#' be selected, for example, via algorithms like Kennard-Stone 
#' (Kennard & Stone, 1969), especially when response values are unavailable.
#'
#' \emph{gesearch} uses \code{Xu} to iteratively remove weak or non-informative 
#' reference samples (\code{Xr}, \code{Yr}), yielding a subset of approximately 
#' \code{k} reference observations. Weak samples are identified based on one or 
#' more of the following criteria:
#' \itemize{
#'   \item{}{They increase the RMSE of prediction on \code{Yu}.}
#'   \item{}{They increase the RMSE of PLS-based reconstruction on \code{Xu}.}
#'   \item{}{They increase dissimilarity between the PLS model subset and 
#'   \code{Xu}.}
#' }
#'
#' A resampling scheme is used to identify reference samples that consistently 
#' appear in the highest-error or most-dissimilar subsets. These are labeled 
#' as weak and removed. The parameter \code{retain} controls the proportion of 
#' samples kept in each iteration.
#'
#' If multiple response variables are provided (i.e., \code{Yu} with multiple 
#' columns or a list of formulas), separate models are built per response. Weak 
#' samples identified across models are combined and removed jointly.
#'
#' Because optimization is based on prediction error and/or spectral 
#' dissimilarity, \code{gesearch} aims to find a subset of \code{k} training 
#' samples that yield the most accurate models for test-like samples 
#' (\code{Xu}).
#'
#' When setting parameters:
#' \itemize{
#'   \item{\code{k}: Number of reference observations to retain. It is also 
#'   the size of each resampling subset. For guidance, see Lobsey et al. 2017.}
#'
#'   \item{\code{b}: Average number of times each training sample is evaluated 
#'   per iteration. Higher values (e.g., >40) yield more stable results but 
#'   increase computational cost.}
#'
#'   \item{\code{retain}: Proportion of training samples retained in each 
#'   iteration. Values >0.9 are recommended for stability.}
#' }
#' @return a list with the following elements:
#' \itemize{
#'   \item{\code{x_local}: Matrix of predictor variables for selected observations.}
#'
#'   \item{\code{y_local}: Matrix of response values corresponding to \code{x_local}.}
#'
#'   \item{\code{indices}: Numeric vector of indices of selected observations from the original training set.}
#'
#'   \item{\code{complete_iter}: Number of completed resampling iterations.}
#'
#'   \item{\code{iter_weakness}: A list with iteration-level performance data:
#'     \describe{
#'       \item{\code{maximum}}{A \code{data.table} showing the highest RMSEs and dissimilarity scores across iterations. RMSEs are for response prediction (if \code{optimization} includes \code{"response"}) and PLS-based spectral reconstruction (if \code{optimization} includes \code{"reconstruction"}). Dissimilarity scores are shown only if \code{optimization} includes \code{"similarity"}.}
#'       \item{\code{cutoff}}{A \code{data.table} reporting threshold RMSEs and dissimilarity values beyond which samples were discarded. Thresholds reflect either fixed proportions or percentiles, depending on \code{control$retain_by} (see \code{\link{search_control}}).}
#'     }
#'   }
#'
#'   \item{\code{samples}: List of sample indices retained in each iteration.}
#'
#'   \item{\code{n_removed}: A \code{data.table} showing number of samples removed per iteration.}
#'
#'   \item{\code{control}: Copy of the control list passed via \code{control}.}
#'
#'   \item{\code{scale}: Logical indicating whether scaling was applied.}
#'
#'   \item{\code{validation_results}: A list with validation output:
#'     \describe{
#'       \item{\code{val_info}}{Indices used for validation in each iteration and a matrix of test-set predictions based on selected samples.}
#'       \item{\code{results}}{Internal validation outcomes from leave-group-out resampling and, if \code{Yr} was provided, prediction results for the test set.}
#'     }
#'   }
#'
#'   \item{\code{final_model}: List of models per response variable. Each model contains:
#'     \itemize{
#'       \item{\code{npls}: Number of PLS components used.}
#'       \item{\code{coefficients}: Regression coefficients for all components.}
#'       \item{\code{bo}: Model intercept(s).}
#'       \item{\code{scores}: PLS scores matrix.}
#'       \item{\code{X_loadings}: Loadings for predictor variables.}
#'       \item{\code{Y_loadings}: Loadings for response variable(s).}
#'       \item{\code{projection_mat}: Projection matrix used in PLS.}
#'       \item{\code{vip}: Variable importance in projection (VIP) matrix.}
#'       \item{\code{selectivity_ratio}: Selectivity ratio matrix per factor (Rajalahti et al., 2009).}
#'       \item{\code{Y}: Matrix of responses used to train the final PLS model.}
#'       \item{\code{weights}: Matrix of PLS weights.}
#'     }
#'   }
#'
#'   \item{\code{seed}: Value of RNG seed used (if provided).}
#'
#'   \item{\code{documentation}: Optional string copied from input \code{documentation} argument.}
#' }
#' @author Leonardo Ramirez-Lopez, Claudio Orellano, Craig Lobsey and Raphael Viscarra Rossel 
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
#' # NOTE: These examples are based on a small dataset where the function may  
#' # not demonstrate substantial benefits. They are intended solely to  
#' # illustrate how the function works.
#' 
#' library(prospectr)
#' data(NIRsoil)
#'
#' # Preprocess the data using detrend + first derivative + Savitzky-Golay
#' sg_det <- savitzkyGolay(
#'   detrend(
#'     NIRsoil$spc,
#'     wav = as.numeric(colnames(NIRsoil$spc))
#'   ),
#'   m = 1, p = 1, w = 7
#' )
#' NIRsoil$spc_pr <- sg_det
#'
#' # Split into training, validation, and testing sets
#' test_x <- NIRsoil$spc_pr[NIRsoil$train == 0 & !is.na(NIRsoil$Ciso), ]
#' test_y <- NIRsoil$Ciso[NIRsoil$train == 0 & !is.na(NIRsoil$Ciso)]
#'
#' set.seed(24)
#' sel <- sample(nrow(test_x), 30)
#'
#' val_x <- test_x[-sel, ]
#' val_y <- test_x[-sel]
#'
#' test_x <- test_x[sel, ]
#' test_y <- test_y[sel]
#'
#' train_y <- NIRsoil$Ciso[NIRsoil$train == 1 & !is.na(NIRsoil$Ciso)]
#' train_x <- NIRsoil$spc_pr[NIRsoil$train == 1 & !is.na(NIRsoil$Ciso), ]
#'
#' # Search based on test set with reference values
#' my_gs <- gesearch(
#'   Xr = train_x, Yr = train_y,
#'   Xu = test_x, Yu = test_y,
#'   k = 50, b = 100, retain = 0.97,
#'   target_size = 350,
#'   method = local_fit_pls(15, modified = TRUE),
#'   optimization = c("reconstruction", "similarity", "response"),
#'   control = search_control(retain_by = "probability"),
#'   scale = FALSE,
#'   verbose = TRUE,
#'   seed = 24
#' )
#' my_preds_gs <- predict(my_gs, val_x)
#'
#' # Search without reference values
#' my_gs_no_response <- gesearch(
#'   Xr = train_x, Yr = train_y,
#'   Xu = test_x, Yu = NULL,
#'   k = 50, b = 100, retain = 0.97,
#'   target_size = 350,
#'   method = local_fit_pls(15, modified = TRUE),
#'   optimization = c("reconstruction", "similarity"),
#'   control = search_control(retain_by = "probability"),
#'   scale = FALSE,
#'   verbose = TRUE,
#'   seed = 24
#' )
#' my_preds_gs_no_response <- predict(my_gs_no_response, val_x)
#'
#' # Search using a more comprehensive set (no reference values)
#' my_gs_no_response2 <- gesearch(
#'   Xr = train_x, Yr = train_y,
#'   Xu = val_x, Yu = NULL,
#'   k = 50, b = 100, retain = 0.97,
#'   target_size = 350,
#'   method = local_fit_pls(15, modified = TRUE),
#'   optimization = c("reconstruction", "similarity"),
#'   control = search_control(retain_by = "probability"),
#'   scale = FALSE,
#'   verbose = TRUE,
#'   seed = 24
#' )
#' my_preds_gs_no_response2 <- predict(my_gs_no_response2, val_x)
#' 
#' 
#' # Using parallel processors
#' n_cores <- 2
#' 
#' if (parallel::detectCores() < 2) {
#'   n_cores <- 1
#' }
#' 
#' # Alternatively:
#' # n_cores <- parallel::detectCores() - 1
#' # if (n_cores == 0) {
#' #  n_cores <- 1
#' # }
#' 
#' library(doParallel)
#' clust <- makeCluster(n_cores)
#' registerDoParallel(clust)
#' 
#' # Alernatively:
#' # library(doSNOW)
#' # clust <- makeCluster(n_cores, type = "SOCK")
#' # registerDoSNOW(clust)
#' # getDoParWorkers()
#' 
#' # search based on a small test set with reference values
#' my_gs <- gesearch(
#'   Xr = train_x,
#'   Yr = train_y,
#'   Xu = test_x,
#'   Yu = test_y,
#'   k = 50, b = 100, retain = 0.97,
#'   target_size = 350,
#'   method = local_fit_pls(15, modified = TRUE),
#'   optimization = c("reconstruction", "similarity", "response"),
#'   control = search_control(retain_by = "probability"),
#'   scale = FALSE,
#'   verbose = TRUE,
#'   seed = 24,
#'   pchunks = 3
#' )
#' }
#' 
#' @importFrom stats as.formula update
#' @export gesearch


## 2020.03.28 (Leo):    Bug fix. Optimization was wrongly set to "reconstruction"
##                      in biter(), even when "response" was selected.
## 2020.03.28 (Leo):    New output. "iter_rmse" is a data.frame containing the
##                      maximum rmse obtained at each iteration after removing
##                      the observations.
## 2020.03.29 (Leo):    Some secondary arguments moved to the new search_control
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
## 2020.06.23 (Leo):    Argument "pls_tune" was removed and passed to search_control
##                      as tune
##                      method must be "local_fit" object
## 2020.10.05 (Leo):    Bug fix:  when optimization = "reconstruction" the local
##                      spectra were not scaled before the reconstruction
## 2020.11.05 (Leo):    Sample search can be conducted using both
##                      "reconstruction" and "response" in argument optimization
## 2020.11.05 (Leo):    verbose argument moved from search_control to gesearch
## 2020.12.03 (Leo):    Multiresponse modeling added

"gesearch" <-
  function(...) {
    UseMethod("gesearch")
  }


#' @aliases gesearch
#' @export
#'
## add an argument to initualize with a given amount of observations, e.g. if 
## the library is 30.000 you could initialize with 15.000 observations and to 
## compensate you can increase the number of resampling iterations
## perhaps a number of initial resampling iterations can also work say we 
## start with 10.000 for the first loop and then it goes to what b specifies
gesearch.default <- function(
    Xr,
    Yr,
    Xu,
    Yu = NULL,
    Yu_lims = NULL,
    k,
    b,
    retain = 0.95,
    target_size = k,
    method = local_fit_pls(pls_c = min(dim(Xr), 10)),
    optimization = "reconstruction",
    control = search_control(),
    group = NULL,
    scale = FALSE,
    verbose = TRUE,
    seed = NULL,
    documentation = character(),
    pchunks = 1,
    intermediate_models = FALSE,
    ...
) {
  crossover <- TRUE # This was previously a parameter
  Yr <- as.matrix(Yr)
  # check inputs
  
  if (ncol(Yr) == 1) {
    if (nrow(Xr) != length(Yr)) {
      stop("The number of spectra in Xr must equal to the length of Yr")
    }
  } else {
    if (is.null(colnames(Yr))) {
      stop("missing column names in Yr")
    }
    if (nrow(Xr) != nrow(Yr)) {
      stop("The number of spectra in Xr must equal to the number of rows in Yr")
    }
  }
  
  if (target_size > nrow(Xr)) {
    stop("'target_size' cannot ne greater than nrow(Xr)")
  }
  
  if (!is.logical(verbose)) {
    stop("verbose argument must be logical")
  }
  
  if (nrow(Xu) < 3) {
    stop(paste0("Too few genes (", nrow(Xu), ") in the target set. Minimum required is 3."))
  }
  
  if (!is.null(Yu)) {
    Yu <- as.matrix(Yu)
    if (ncol(Yu) == 1) {
      if (nrow(Xu) != length(Yu)) {
        stop("The number of spectra in Xu must equal to the length of Yu")
      }
    } else {
      if (!is.null(Yu_lims)) {
        stop("Yu_lims can only be used for one single response variable for the moment.")
      }
      if (ncol(Yr) != ncol(Yu)) {
        stop("Different number of columns in Yr and Yu")
      }
      
      if (nrow(Xu) != nrow(Yu)) {
        stop("The number of spectra in Xu must equal to the number of rows in Yu")
      }
      
      if (is.null(colnames(Yu))) {
        stop("Missing column names in Yu")
      }
      if (!all(colnames(Yr) == colnames(Yu))) {
        stop("Column names in Yr and Yu do not match")
      }
    }
  }
  
  if ("response" %in% optimization) {
    if (!any(!is.na(Yu))) {
      stop("No response values, these are required when optimization == 'response'")
    }
  }
  
  if (anyNA(Xu)) {
    stop("Missing values detected in the predictor variables of test/Xu")
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
  
  target_size <- round(target_size)
  if (target_size < k) {
    stop("Argument 'taget_size' must be equal or larger than 'k'")
  }
  
  if (retain > 1 | retain <= 0) {
    stop("Argument 'retain' must be a numerical value larger than 0 and below 1")
  }
  
  if ("response" %in% optimization & is.null(Yu)) {
    stop("When optimization = 'response', Yu values must be provided")
  }
  
  allowed <- c("response", "reconstruction", "similarity", "range", "fresponse")
  
  bad <- setdiff(optimization, allowed)
  if (length(bad) > 0L) {
    stop(
      sprintf(
        "Invalid 'optimization' option(s): %s. Allowed: %s",
        paste(shQuote(bad), collapse = ", "),
        paste(shQuote(allowed), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  if (!"local_fit" %in% class(method)) {
    stop("Method must be of class 'local_fit'")
  }
  
  if (method$method != "pls") {
    stop(
      "The only method alowed for the moment is 'pls'", 
      " generated with the 'local_fit_pls()' fucntion"
    )
  }
  
  if ("reconstruction" %in% optimization & control$tune & method$method %in% c("pls", "wapls")) {
    warning(
      "pls factors are not tuned when optimization = 'reconstruction'",
      " instead they are fixed to the one(s) provided in the `method` argument,",
      " therefore the `tune` option passed to control has been ignored."
    )
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
  Xu <- as.matrix(Xu)
  
  r <- 1 - round(retain, .Machine$sizeof.longdouble)
  
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
  # names(U) <- names(V) <- sl_idx
  s_to_drop <- cuts <- max_weakness_score <- NULL
  pp <- r
  cat_progress <- cat_iter(c("\\", "|", "/", "-"))
  cat_progress2 <- cat_iter(c("\\", "|", "/", "-"))
  
  pool_sizes_log <- c()
  complete_iterations <- 0
  selected <- list()
  while (length(k_idx) > (target_size * (1 + pp)) & sum(!is.na(sl_idx)) > target_size) {
    # Step 1 - Initialise K as a subset of the SL, initially full. k_idx only
    # contains those SL samples still in K
    #
    # select k_idx as those in the SL not marked as dropped (zero)
    k_idx <- sl_idx[!is.na(sl_idx)]
    selected <- c(selected, list(k_idx))
    
    outer_idx <- outer_idx + 1
    
    # calculate the quantity of samples removed in this iteration
    cull_quantity <- round(length(k_idx) * round(r, .Machine$sizeof.longdouble))
    
    ##  This sub-loop contains step 2-4 of the algorithm (See step 5!). This is
    ##  the B iteration.
    B <- round(length(k_idx) * b / k)
    if (verbose) {
      # Clear the full current line
      cat("\r\033[K")  # \033[K erases to end of line
      
      cat(paste0("Generation ", outer_idx, ":\t"))
      
      if (outer_idx == 1) {
        cat(paste0(
          "\033[34m", cat_progress(),
          " Initial gene pool size: ", length(k_idx), "\t\033[39m"
        ))
      } else {
        cat(paste0(
          "\033[34m", cat_progress(),
          " Current gene pool size: ", length(k_idx), "\t\033[39m"
        ))
      }
      
      cat(paste0(
        "\033[34m", cat_progress2(),
        " Individuals to evaluate: ", B, "\033[39m"
      ))
      cat("\r")
    }
    
    
    ## Step 2 - sample a training data set of size k from K without replacement
    if (!is.null(seed)) {
      set.seed(seed + B)
    }
    
    if (outer_idx == 1) {
      ismpl <- replicate(
        n = B,
        expr = sample(x = k_idx, size = k, replace = FALSE)
      )
    } else {
      # Replace genes that were dropped in last iteration by NA
      ismpl[ismpl %in% which(new_nas)] <- NA
      
      if (crossover) {
        ismpl <- replicate(
          n = B,
          expr = combine_individuals(ismpl, k, k_idx)
        )
      } else {
        ismpl <- replicate(
          n = B,
          expr = sample(x = k_idx, size = k, replace = FALSE)
        )
      }
    }
    ## this iterator object will avoid to use much memory in the next foreach
    ## loop this iterator takes from the big matrix only what is needed and
    ## therefore the whole matrix is not put in the memory.
    ## This makes parallel computations more memory friendly
    if (control$allow_parallel & getDoParRegistered()) {
      leftc <- B %% pchunks
      if (leftc > 0) {
        req_iter <- floor(B / pchunks) + leftc
      } else {
        req_iter <- B / pchunks
      }
    } else {
      req_iter <- B
      pchunks <- 1
    }
    
    itersubs <- ithrssubsets(
      x = Xr,
      y = Yr,
      group = group,
      indx = ismpl,
      chunksize = pchunks
    )
    
    ## Step 3 calibrate a PLS model using the selected SL samples
    ## Step 4 - validate on the 'm' site specific samples
    ## Step 5 - repeat 2 to 4 B times
    
    if ("range" %in% optimization) {
      if (is.null(Yu_lims)) {
        stop("Yu_lims is missing for 'range' optimization")
      }
      if (!length(Yu_lims) == 2) {
        stop("Yu_lims must be a numerical vector of length 2")
      }
    }
    
    
    if (!is.null(Yu_lims) & "range" %in% optimization) {
      user_min_y <- min(Yu_lims)
      user_max_y <- max(Yu_lims)
    } else {
      user_min_y <- user_max_y <- NULL
      optimization <- optimization[!optimization %in% "range"]
    }
    
    results_df <- biter(
      itersubs = itersubs,
      Xu = Xu,
      Yu = Yu,
      y_lim_left = user_min_y,
      y_lim_right = user_max_y,
      iter_sequence = seq(1, req_iter),
      n = nrow(ismpl),
      optimization = optimization,
      ncomp = method$pls_c,
      tune = control$tune,
      p = control$p,
      number = control$number,
      scale = scale,
      max_iter = 1,
      tol = 1e-6,
      seed = seed,
      modified = method$modified,
      allow_parallel = control$allow_parallel,
      pchunks = pchunks
    )
    
    ## Step 6 - starts here
    # iterate through all B iteration results and increment the rankings (rmse)
    # and test count for each sample
    drop_results <- apply(
      results_df$weakness_scores_subset,
      MARGIN = 2,
      FUN = drop_indices,
      uu = U,
      vv = V,
      r = r,
      obs_indices = sl_idx,
      resampling_indices = results_df$sampleidx,
      retain_by = control$retain_by,
      cull_quantity = cull_quantity,
      percentile_type = control$percentile_type,
      max_pls_c = max(method$pls_c)
    )
    
    was_interrupted <- sapply(drop_results, FUN = function(x) x$interrupted)
    
    if (any(was_interrupted)) {
      analysis_name <- gsub(
        "rmse_|score_dis|residual_| deviation_from", 
        "", names(drop_results)
      )
      analysis_name <- gsub("_", " ", analysis_name)
      messages <- sapply(
        which(was_interrupted),
        FUN = function(x, name, i) {
          message(paste0(name[i], ": ", x[[i]]$message))
        },
        x = drop_results,
        name = analysis_name
      )
      break
    }
    pool_sizes_log <- c(pool_sizes_log, length(k_idx))
    
    cuts <- rbind(cuts, sapply(drop_results, FUN = function(x) x$cutoff))
    max_weakness_score <- rbind(
      max_weakness_score, 
      sapply(drop_results, FUN = function(x) x$max_weakness_score)
    )
    
    new_nas <- !complete.cases(
      sapply(drop_results, FUN = function(x) x$obs_indices)
    )
    sl_idx[new_nas] <- NA
    complete_iterations <- complete_iterations + 1
    s_to_drop <- c(s_to_drop, sum(is.na(sl_idx)))
    
    # Inside the while loop (after updating pool_sizes_log)
    if (length(pool_sizes_log) >= control$stagnation_limit) {
      recent_sizes <- tail(pool_sizes_log, control$stagnation_limit)
      if (all(recent_sizes == recent_sizes[1])) {
        cat(paste0(
          "\n\033[31mEarly termination: Gene pool size remained at ",
          recent_sizes[1], " for ", control$stagnation_limit,
          " consecutive iterations. Stopping.\033[39m\n"
        ))
        break
      }
    }
  }
  
  # Finalise k_idx following the last iteration
  k_idx <- sl_idx[!is.na(sl_idx)]
  # make sure the final k_idx is at least of target_size
  selected <- c(selected, list(k_idx))
  msizes <- sapply(selected, FUN = length)
  k_idx <- selected[[max(which(msizes >= target_size))]]
  selected <- selected[which(msizes >= target_size)]
  
  models_per_generation <- vector("list", length(selected) - 1)
  if (intermediate_models) {
    if (verbose) {
      cat("\nFitting intermediate models...")
    }
    for (i in 1:(length(selected) - 1)) {
      ith_selected <- selected[[i]]
      
      # Prepare the final models
      pls_models <- NULL
      validation_results <- NULL
      for (j in 1:ncol(Yr)) {
        
        pls_val <- get_plsr(
          X = as.matrix(Xr),
          Y = as.matrix(Yr[, j, drop = FALSE]),
          indices = ith_selected,
          method = method, # this is an object of class 'local_fit'
          center = TRUE,
          scale = scale,
          number = control$number,
          group = group[ith_selected],
          retrieve = FALSE,
          tune = TRUE,
          max_iter = 1,
          tol = 1e-6,
          seed = seed, 
          p = control$p
        )
        
        validation <- list(
          val_info = list(train_resamples = pls_val$train_resamples),
          results = list(train = pls_val$cv_results)
        )
        
        if (!is.null(Yu)) {
          yuhat <- predict(pls_val, newdata = Xu, ncomp = method$pls_c)
          
          validation$val_info$test_predictions <- yuhat
          validation$results$test <- get_validation_metrics(yuhat, Yu)
          
          tmpx <- as.matrix(rbind(Xr[ith_selected, ], Xu))
          tmpy <- as.matrix(c(Yr[ith_selected, j, drop = FALSE], Yu[, j, drop = FALSE]))
          tmpx <- tmpx[complete.cases(tmpy), ]
          tmpy <- tmpy[complete.cases(tmpy), , drop = FALSE]
        } else {
          tmpx <- as.matrix(Xr[ith_selected, ])
          tmpy <- as.matrix(Yr[ith_selected, j, drop = FALSE])
        }
        
        intermediate_pls <- opls_get_all(
          X = tmpx,
          Y = tmpy,
          ncomp = method$pls_c,
          scale = scale,
          maxiter = 1,
          tol = 1e-6,
          algorithm = ifelse(method$modified, "mpls", "pls")
        ) |> assign_pls_names(x_names = colnames(Xr))
        
        pls_models[[j]] <- intermediate_pls
        validation_results[[j]] <- validation
      }
      if (!is.null(colnames(Yu))) {
        names(pls_models) <- names(validation_results) <- colnames(Yu)
      }
      models_per_generation[[i]] <- list(
        subset_size = length(ith_selected),
        pls_models = pls_models,
        validation = validation_results
      )
    }
  }
  
  
  if (verbose) {
    cat("\nFitting final model on", length(k_idx), "selected genes... \n")
    if (!is.null(Yu)) {
      if (sum(complete.cases(Yu)) != 0) {
        cat("and", sum(complete.cases(Yu)), "target genes for which response values were available... \n")
      }
    }
  }
  
  
  # Prepare the final models
  pls_models <- NULL
  validation_results <- NULL
  for (j in 1:ncol(Yr)) {
    
    pls_val <- get_plsr(
      X = as.matrix(Xr),
      Y = as.matrix(Yr[, j, drop = FALSE]),
      indices = k_idx,
      method = method, # this is an object of class 'local_fit'
      center = TRUE,
      scale = scale,
      number = control$number,
      group = group[k_idx],
      retrieve = FALSE,
      tune = TRUE,
      max_iter = 1,
      tol = 1e-6,
      seed = seed, 
      p = control$p
    )
    
    validation <- list(
      val_info = list(train_resamples = pls_val$train_resamples),
      results = list(train = pls_val$cv_results)
    )
    
    if (!is.null(Yu)) {
      yuhat <- predict(pls_val, newdata = Xu, ncomp = method$pls_c)
      
      validation$val_info$test_predictions <- yuhat
      validation$results$test <- get_validation_metrics(yuhat, Yu)
      
      tmpx <- as.matrix(rbind(Xr[k_idx, ], Xu))
      tmpy <- as.matrix(c(Yr[k_idx, j, drop = FALSE], Yu[, j, drop = FALSE]))
      tmpx <- tmpx[complete.cases(tmpy), ]
      tmpy <- tmpy[complete.cases(tmpy), , drop = FALSE]
    } else {
      tmpx <- as.matrix(Xr[k_idx, ])
      tmpy <- as.matrix(Yr[k_idx, j, drop = FALSE])
    }
    
    finalpls <- opls_get_all(
      X = tmpx,
      Y = tmpy,
      ncomp = method$pls_c,
      scale = scale,
      maxiter = 1,
      tol = 1e-6,
      algorithm = ifelse(method$modified, "mpls", "pls")
    ) |> assign_pls_names(x_names = colnames(Xr))
    
    pls_models[[j]] <- finalpls
    validation_results[[j]] <- validation
  }
  
  if (!is.null(colnames(Yu))) {
    names(pls_models) <- names(validation_results) <- colnames(Yu)
  }
  
  
  iter_weakness <- NULL
  iter_weakness$maximum <- data.table(
    iteration = 1:nrow(max_weakness_score),
    max_weakness_score
  )
  
  iter_weakness$cutoff <- data.table(
    iteration = 1:nrow(max_weakness_score),
    cuts
  )
  
  names(selected) <- paste0("iteration_", 1:length(selected))
  
  resultsList <- list(
    x_local = Xr[k_idx, ],
    y_local = Yr[k_idx, , drop = FALSE],
    indices = k_idx,
    complete_iter = complete_iterations,
    iter_weakness = iter_weakness,
    samples = selected,
    n_removed = data.table(
      iteration = 1:length(s_to_drop),
      removed = c(s_to_drop[1], diff(s_to_drop)),
      cummulative = s_to_drop
    ),
    control = control,
    scale = scale,
    validation_results = validation_results,
    final_models = pls_models,
    intermediate_models = models_per_generation, 
    seed = seed,
    documentation = documentation
  )
  
  if (!intermediate_models) {
    resultsList <- resultsList[!names(resultsList) %in% "intermediate_models"]
  }
  
  attr(resultsList, "call") <- call_f
  class(resultsList) <- c("gesearch", "list")
  
  resultsList
}


#' @aliases gesearch
#' @importFrom stats quantile complete.cases diffinv na.pass model.extract 
#' model.frame model.matrix
#' @export
gesearch.formula <- function(
    formula,
    train,
    test,
    k,
    b,
    target_size,
    method,
    ...,
    na_action = na.pass
) {
  if (!inherits(formula, "formula")) {
    stop("'formula' is only for formula objects")
  }
  
  call_f <- match.call()
  
  if (missing(method)) {
    stop("'method' is missing")
  }
  definition <- sys.function(sys.parent())
  mf <- match.call(expand.dots = FALSE)
  formals <- formals(definition)
  
  if (!"na_action" %in% names(mf)) {
    mf[["na_action"]] <- formals[["na_action"]]
    match.call(definition, mf, TRUE)
  }
  
  ## Get the model frame
  
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

    test <- local({
      tmp <- make.unique(c(names(test), ".tmp_y"))[length(names(test)) + 1L]
      test[[tmp]] <- NA
      names(test)[names(test) == tmp] <- yname
      test
    })
    warning(paste(
      yname, "not found in test. Missing values (NAs) were assigned to",
      yname, "in test."
    ))
  }
  
  mfu <- model.frame(mfu, data = test, na.action = NULL)
  mfr <- eval(mfr, parent.frame())
  
  trms <- attr(mfr, "terms")
  
  formulaclasses <- list(attr(trms, "dataClasses"))
  
  attr(trms, "intercept") <- 0
  xr <- model.matrix(trms, model.frame(mfr, drop.unused.levels = T))
  yr <- model.extract(mfr, "response")
  
  xu <- model.matrix(trms, mfu)
  yu <- model.extract(mfu, "response")
  
  
  if ("optimization" %in% names(input_list)) {
    if ("response" %in% input_list$optimization) {
      if (sum(is.na(yu)) == length(yu)) {
        stop("When optimization = 'response', response values must be provided in test")
      }
    }
  }
  
  rsl <- gesearch(
    Xr = xr,
    Yr = yr,
    Xu = xu,
    Yu = yu,
    k = k,
    b = b,
    target_size = target_size,
    method = method,
    ...
  )
  rsl$formula <- formula
  rsl$dataclasses <- formulaclasses
  
  attr(rsl, "call") <- call_f
  class(rsl) <- c("gesearch", "gesearch.formula", "list")
  
  rsl
}

#' @aliases gesearch
#' @importFrom stats quantile complete.cases diffinv na.pass model.extract model.frame model.matrix
#' @export
gesearch.list <- function(formula,
                          train,
                          test,
                          k,
                          b,
                          method,
                          ...,
                          na_action = na.pass) {
  is_formula <- sapply(formula, FUN = function(x) inherits(x, "formula"))
  if (any(!is_formula)) {
    stop("non-forumla objects in formula")
  }
  
  right_side <- sapply(formula, FUN = function(x) labels(terms(x)))
  left_side <- sapply(formula, FUN = function(x) all.vars(update(x, . ~ 1)))
  
  mresponse <- paste(c("matrix_response", left_side), collapse = "_")
  mfml <- as.formula(paste(mresponse, "~", unique(right_side)))
  
  if (length(unique(right_side)) > 1) {
    stop("The right hand side variables must be the same across all the formulas in the list")
  }
  
  call_f <- match.call()
  
  
  if (missing(method)) {
    stop("'method' is missing")
  }
  definition <- sys.function(sys.parent())
  mf <- match.call(expand.dots = FALSE)
  formals <- formals(definition)
  
  mf[["formula"]] <- mfml
  
  if (!"na_action" %in% names(mf)) {
    mf[["na_action"]] <- formals[["na_action"]]
    match.call(definition, mf, TRUE)
  }
  
  ## Get the model frame
  mr <- match(x = c("formula", "train", "na_action"), table = names(mf))
  mu <- match(x = c("formula", "test"), table = names(mf))
  
  mfr <- mf[c(1, mr)]
  mfu <- mf[c(1, mu)]
  
  names(mfr)[names(mfr) %in% "na_action"] <- "na.action"
  names(mfr)[names(mfr) %in% "train"] <- "data"
  names(mfu)[names(mfu) %in% "test"] <- "data"
  
  yname <- left_side
  
  mfr[[1]] <- mfu[[1]] <- as.name("model.frame")
  
  input_list <- list(...)
  
  if (!any(yname %in% colnames(eval(mfu$data)))) {
    if ("optimization" %in% names(input_list)) {
      if (input_list$optimization == "response") {
        stop("When optimization = 'response', response values must be provided in test")
      }
    }
    
    # for (i in 1:length(yname)) {
    #   test <- cbind(test, yyyyyyyy801124 = NA)
    #   names(test)[names(test) == "yyyyyyyy801124"] <- yname[i]
    #   warning(paste(
    #     yname[i], "not found in test. Missing values (NAs) were assigned to",
    #     yname[i], "in test."
    #   ))
    # }
    missing <- setdiff(yname, names(test))
    
    if (length(missing) > 0L) {
      test[missing] <- NA
      warning(sprintf(
        "%s not found in test. Assigned NA.",
        paste(missing, collapse = ", ")
      ))
    }
  }
  
  train[[mresponse]] <- as.matrix(train[, left_side])
  test[[mresponse]] <- as.matrix(test[, left_side])
  
  mfu <- model.frame(mfu, data = test, na.action = NULL)
  mfr <- model.frame(mfr, data = train, na.action = NULL)
  
  trms <- attr(mfr, "terms")
  
  formulaclasses <- NULL
  wformulaclasses <- attr(trms, "dataClasses")
  wformulaclasses[[1]] <- "nmatrix.1"
  for (i in 1:length(yname)) {
    formulaclasses[[i]] <- wformulaclasses
    names(formulaclasses[[i]])[1] <- yname[i]
  }
  formulaclasses
  
  attr(trms, "intercept") <- 0
  xr <- model.matrix(trms, model.frame(mfr, drop.unused.levels = T))
  yr <- model.extract(mfr, "response")
  
  xu <- model.matrix(trms, mfu)
  yu <- model.extract(mfu, "response")
  
  if (!class(yu) %in% "matrix") {
    yu <- as.matrix(yu)
  }
  
  for (i in 1:ncol(yu)) {
    if ("optimization" %in% names(input_list) & sum(is.na(yu[, i])) == length(yu[, i])) {
      mss <- paste(
        "When optimization = 'response', response values must be provided in test.",
        "Check", colnames(yu)[i]
      )
      stop(mss)
    }
  }
  
  rsl <- gesearch(
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
  class(rsl) <- c("gesearch", "gesearch.list", "list")
  
  rsl
}


#' Predict from a \code{gesearch} object
#'
#' @description
#' Generates predictions from a fitted \code{gesearch} object using the stored
#' final PLS models, or (optionally) returns predictions for all intermediate
#' generations as well. If the model was fitted with a formula, \code{newdata}
#' may be a \code{data.frame} (which will be processed via the stored design
#' matrix terms) or a \code{matrix}. Otherwise, \code{newdata} must be a
#' \code{matrix}.
#'
#' @param object A fitted object of class \code{gesearch}. It is expected to
#'   contain elements such as \code{$final_models}, \code{$intermediate_models}
#'   (optional), \code{$formula} (optional), \code{$dataclasses}, and
#'   \code{$scale}.
#' @param newdata New data on which to obtain predictions. If \code{object}
#'   was created with a formula, \code{newdata} may be a \code{data.frame}
#'   (containing all required predictor variables) or a \code{matrix}. If no
#'   formula was used, \code{newdata} must be a \code{matrix}.
#' @param type Character; prediction type. Currently only \code{"response"} is
#'   supported. 
#' @param what Character; which models to use for prediction. Either
#'   \code{"final"} (default; return predictions from the final models) or
#'   \code{"all_generations"} (return a list of predictions for each
#'   intermediate generation plus the final generation).
#' @param ... Ignored (reserved for future use).
#'
#' @details
#' If the model was fitted with a formula, the function validates and transforms
#' \code{newdata} to the appropriate model matrix using the stored terms
#' (no intercept). Factor levels not present in \code{newdata} are dropped via
#' \code{drop.unused.levels = TRUE}. For non-numeric predictors, the function
#' supports the case where \code{newdata} is a matrix with appropriately named
#' columns (matching the fitted design matrix), otherwise a \code{data.frame}
#' is required.
#'
#' For \code{type = "response"}, predictions are generated by delegating to
#' \code{predict_opls()} for each stored PLS model. Column names of the output
#' are of the form \code{"pls_1"}, \code{"pls_2"}, \dots. Row names (if present
#' in \code{newdata}) are propagated to the prediction matrices.
#'
#' When \code{what = "all_generations"}, the return value is a named list with
#' one element per generation (e.g., \code{"generation_1"}, \code{"generation_2"},
#' \dots), where each element is itself a list of prediction matrices (one per
#' response), and the last element corresponds to the final models.
#'
#' @return
#' \describe{
#'   \item{\code{type = "response"} and \code{what = "final"}}{
#'     A named list of numeric matrices (one per response), each with
#'     \code{nrow(newdata)} rows and columns corresponding to the number of PLS
#'     components used in the respective model.
#'   }
#'   \item{\code{type = "response"} and \code{what = "all_generations"}}{
#'     A named list of generations; each generation is a named list of
#'     prediction matrices as above.
#'   }
#'   \item{\code{type = "scores"}}{
#'     An error is raised (not yet implemented).
#'   }
#' }
#' For \code{type = "scores"}: an error is raised (not yet implemented).
#'
#' @section Input validation:
#' The function checks that all required predictor variables are present in
#' \code{newdata} when a formula was used. It also enforces that \code{newdata}
#' is a \code{matrix} when no formula is available.
#'
#' @seealso
#' \code{\link{predict_opls}}, \code{\link[stats]{terms}},
#' \code{\link[stats]{model.matrix}}, \code{\link[stats]{delete.response}}
#'
#' @aliases gesearch
#' @importFrom stats .MFclass terms delete.response
#' @export
predict.gesearch <- function(
    object, 
    newdata, 
    type = "response", 
    what = c("final", "all_generations"), 
    ...
  ) {
  
  what <- match.arg(what)
  
  allowed_types <- c("response")
  if (!is.character(type) || length(type) != 1L || is.na(type) || !(type %in% allowed_types)) {
    stop(
      "`type` must be one of: ", 
      paste(shQuote(allowed_types), collapse = ", "), 
      call. = FALSE
    )
  }
  
  if (missing(newdata)) {
    stop("newdata is missing")
  }
  
  if (!is.null(object$formula)) {
    dcls <- object$dataclasses[[1]][-1]
    
    if (!("matrix" %in% class(newdata) | "data.frame" %in% class(newdata))) {
      stop(paste0(
        "When predicting from objects of class 'gesearch' fitted with ",
        "formula, the argument 'newdata' must be a 'data.frame' ",
        "or alternatively a 'matrix'"
      ))
    }
    
    if ("data.frame" %in% class(newdata)) {
      if (!all(names(dcls) %in% names(newdata))) {
        mss <- names(dcls)[!names(dcls) %in% names(newdata)]
        stop(paste(
          "The following predictor variables are missing:",
          paste(mss, collapse = ", ")
        ))
      }
    }
    
    if (any(dcls != "numeric")) {
      if ("matrix" %in% class(newdata) & length(dcls) == 1) {
        if (.MFclass(newdata) == dcls) {
          pnames <- gsub(
            names(dcls), "",
            colnames(object$final_models[[1]]$coefficients)
          )
          if (all(pnames %in% colnames(newdata))) {
            newdata_temp <- newdata
            newdata <- data.frame(rep(NA, nrow(newdata)))
            colnames(newdata) <- names(object$dataclasses[[1]])[1]
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
    newdata <- model.matrix(oterms, model.frame(mf, drop.unused.levels = TRUE))
  }
  
  if (!"matrix" %in% class(newdata)) {
    stop("Argument 'newdata' must be a 'matrix'")
  }
  
  # if(all(pnames %in% colnames(newdata)))
  if (type == "response") {
    preds <- NULL
    for (i in 1:length(object$final_models)) {
      preds[[i]] <- predict_opls(
        bo = object$final_models[[i]]$bo,
        b = t(object$final_models[[i]]$coefficients),
        ncomp = object$final_models[[i]]$npls,
        newdata = as.matrix(newdata),
        scale = object$scale,
        Xscale = object$final_models[[i]]$transf$Xscale
      )
      rownames(preds[[i]]) <- rownames(newdata)
      colnames(preds[[i]]) <- paste0("pls_", 1:ncol(preds[[i]]))
    }
    names(preds) <- names(object$final_models)
  }
  
  
  
  if (type == "response" & what == "all_generations") {
    if (is.null(object$intermediate_models)) {
      warning(
        "The model object does not contain intermediate models.", 
        "Predictions were generated only from the final model."
      )
    }
    intermediate_preds <- NULL
    for (j in 1:length(object$intermediate_models)) {
      intermediate_preds[[j]] <- list()
      for (i in 1:length(object$intermediate_models[[j]]$pls_models)) {
        intermediate_preds[[j]][[i]] <- predict_opls(
          bo = object$intermediate_models[[j]]$pls_models[[i]]$bo,
          b = t(object$intermediate_models[[j]]$pls_models[[i]]$coefficients),
          ncomp = object$intermediate_models[[j]]$pls_models[[i]]$npls,
          newdata = as.matrix(newdata),
          scale = object$scale,
          Xscale = object$intermediate_models[[j]]$pls_models[[i]]$transf$Xscale
        )
        rownames(intermediate_preds[[j]][[i]]) <- rownames(newdata)
        colnames(intermediate_preds[[j]][[i]]) <- paste0(
          "pls_", 1:ncol(intermediate_preds[[j]][[i]])
        )
      }
      names(intermediate_preds[[j]]) <- names(
        object$intermediate_models[[j]]$pls_models
      )
    }
    intermediate_preds[[j + 1]] <- preds
    names(intermediate_preds) <- paste0(
      "generation_", 1:length(intermediate_preds)
    )
    preds <- intermediate_preds
  }
  
  # if (type == "scores") {
  #   stop("not yet implemented") # FIXME!!
  # }
  preds
}
