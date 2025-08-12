#' @title Build predictive/fucntional libraries
#' @description
#' This function builds a library of predictive functions based on memory-based learning
#' @usage
#' fit_library(Xr, Yr, sm, k, k_diss, k_range, ws = NULL, pls_c, 
#'             pc_selection = list("opc", min(dim(Xr), 40)),
#'             group, build_around, gh = TRUE, center = TRUE,
#'             scale = FALSE, diss_predictors = FALSE,
#'             metric = "rmse", return_best = TRUE,
#'             pls_max_iter = 1, pls_tol = 1e-6,
#'             documentation = character(), ...)
#' @param Xr a numeric matrix (or data.frame) of predictor variables of dimentions `n * p` corresponding to the reference data (observations in rows and variables in columns). 
#' @param Yr a numeric vector of length `n` containing the values of the response variable corresponding to the reference data. Missing values might be allowed (see details).
#' @param diss_method a character string indicating the spectral dissimilarity metric to be used in the selection of the nearest neighbours of each observation. 
#'        Options are: 
#'        \itemize{
#'        \item{`"pc"`: Principal components dissimilarity: Mahalanobis dissimilarity computed on the principal components space.}
#'        \item{`"cor"`: Correlation dissimilarity.}
#'        \item{`"pls"`: Partial least squares dissimilarity: Mahalanobis dissimilarity computed on the partial least squares space.}
#'        \item{`"euclid"`: Euclidean dissimilarity.}
#'        \item{`"cosine"`: Cosine dissimilarity.}
#'        \item{`"sidF"`: Spectral information divergence computed on the spectral variables.}
#'        \item{`"sidD"`: Spectral information divergence computed on the density distributions of the spectra.}
#'        \item{`"loc.pc"`: Dissimilarity estimation based on local principal components.}
#'        \item{`"loc.pls"` Dissimilarity estimation based on local partial least squares.}
#'        }
#'        The `"pc"` spectral dissimilarity metric is the default. If the `"sidD"` is chosen, the default parameters of the `sid` function are used however they can be specified by passing them them as additional arguments `...`.  
#'        
#'        ======= NOT YET IMPLEMENTED!!! ======= --> --> --> Alternatively a squared dissimilarity matrix of size `n * n` can also be passed. In this case, each row of the matrix must correspond to the observations. The nearest neighbors search is conducted column-wise.
#' @param k a vector of integers specifying the sequence of k nearest neighbours to be tested. Either `k` or `k_diss` must be specified. This `vector` will be automatically sorted into ascending order. The values of vectors of class numeric, will be coerced to the next upper inetegers.
#' @param k_diss NOT YET IMPLEMENTED a numeric `vector` specifying the sequence of dissimilarity thresholds to be tested for the selection of the nearest neighbors around each observation. For a given observation, its neighbors are those that exhibit a dissimilarity value equal to or below a given threshold. These thresholds depend on the corresponding dissimilarity measure specified in `sm`. Either `k` or `k_diss` must be specified. 
#' @param k_range NOT YET IMPLEMENTED a vector of two integers specifying the minimum (first value) and the maximum (second value) number of neighbours allowed when the `k_diss` argument is used. 
#' @param ws an optional odd integer value which specifies a window size when the correlation dissimilarity is used (i.e `sm = "cor"`). If not specified, the standard correlation dissimilarity is computed (i.e no moving window). Default is `NULL`.
#' @param pls_c a vector of two integers indicating the minimum (first value) and maximum (second value) number of pls components to be used in the weighted average pls (wapls) regressions. 
#' @param pc_selection a list specifying the details of the method to be used for identifying the number of components to be retained for computing the dissimilarities between samples when the method passed to `sm` is one of the following options: `"pc"`, `"loc.pc"`, ` "pls"` or ` "loc.pls"`. This list must contain two elements in the following order: \itemize{
#'        \item{`method`:}{the method for selecting the number of components. Possible options are:  `"opc"` (optimized pc selection based on Ramirez-Lopez et al. (2013a, 2013b). See the \code{\link{orthoProjection}} function for more details;  `"cumvar"` (for selecting the number of principal components based on a given cumulative amount of explained variance); `"var"` (for selecting the number of principal components based on a given amount of explained variance); and  `"manual"` (for specifying manually the desired number of principal components)}
#'        \item{`value`:}{a numerical value that complements the selected method. If `"opc"` is chosen, it must be a value indicating the maximal number of principal components to be tested (see Ramirez-Lopez et al., 2013a, 2013b). If `"cumvar"` is chosen, it must be a value (larger than 0 and below 1) indicating the maximum amount of cumulative variance that the retained components should explain. If `"var"` is chosen, it must be a value (larger than 0 and below 1) indicating that components that explain (individually) a variance lower than this threshold must be excluded. If `"manual"` is chosen, it must be a value specifying the desired number of principal components to retain.
#'        }}
#' The default method for the `pc_selection` argument is `"opc"` and the maximum number of principal components to be tested is set to `min(dim(Xr), 40)`, where code{Xr} is the matrix of reference predictors in the \code{\link{mbl}} function.
#' Optionally, the `pc_selection` argument admits `"opc"` or `"cumvar"` or `"var"` or `"manual"` as a single character string. In such a case the default for `"value"` when either `"opc"` or `"manual"` are used is `min(dim(Xr), 40)`. When `"cumvar"` is used the default `"value"` is set to 0.99 and when `"var"` is used the default `"value"` is set to 0.01.
#' @param group (PUT THE DETAILS OF THE 1_NN VALIDATION IN DETAILS) an optional `factor` (or `vector` that can be coerced to a `factor` by `as.factor`) that assigns to each observation in `Xr` a group/class label (e.g. groups can be given by spectra collected from the same batch of measurements, from the same sample, from samples with very similar origin, etc). This is taken into account for internal validation of pls models (and factor optimization) to avoid pseudo-replication. For validation, the model of a given target observation is first fitted based on its neighbors and excluding that observation (and all the samples belonging to its group), then the model is used to predict the response variable of the target observation. See details.
#' @param build_around NOT YET IMPLEMENTED! a vector of integers of length `m` (`m < n`) indicating the indices of the samples in the reference set for which the local models must be fitted. Default is `NULL` which means that local models are fitted for all the `n` observations in the reference set. See details.
#' @param gh a logical indicating whether or not to compute and return the Mahalanobis distance (in the pls space) between each element in `Xr` and the center of `Xr`.
#' @param center a logical indicating whether or not the predictor variables must be centered at each local segment (before regression).
#' @param scale a logical indicating whether or not the predictor variables must be scaled to unit variance at each local segment (before regression).
#' @param diss_predictors a logical indicating if the local (square symmetric) dissimilarity matrix corresponding the selected neighbourhood shall be used as source of additional predictors (i.e the columns of this local matrix are treated as predictor variables). In some cases this may result in an improvement  of the prediction performance (Ramirez-Lopez et al., 2013a). 
#' @param metric a character value indicating what model perfoamance statistics shall be used to select the best set of parameters (i.e. number of neighbors or threshold distances and pls factors)  to fit the final model. Options are `"rmse"` (default) or `"r2"` (coefficient of determination).
#' @param return_best a logical indicating if the final library of functions using the optimal parameters found shall be returned.
#' @param pls_max_iter (BETTER DESCRIPTION REQUIRED) maximum number of iterations for the partial least squares methods.
#' @param pls_tol (BETTER DESCRIPTION REQUIRED) for convergence in the partial orthogonal scores partial least squares regressions using the nipals algorithm. Default is 1e-6
#' @param documentation (BETTER DESCRIPTION REQUIRED) an optional character string for documentating the call to this function.
#' @param ...  arguments passed to the \code{\link{dissimilarity}} function.
#' @details
#' By default, this function fits a local model for each of the `n` observations in the reference set, i.e. a library of `n` fucntions is built. Each local model is fitted with the query observation and its nearest neighbors.
#' Alternatively, the user may specificy the samples in the reference set around which the library of functions must be built. This is done by indicating in the `build_around` argument the indices of the samples in the reference set. In tis case, each function is fitted by using the nearest neighbors found in the entire set of `n` reference observations. The `build_around` argument may be useful in cases where the number of observations in the reference set is extremely large.
#' Missing values in `Yr` are allowed. If they exceed 25% of the total number of observations a warning message will be printed. The local model of an observation with missing `Yr` value is fited only with its nearest neighbors (i.e without including this observation).
#' @references 
#' Rajalahti, T., Arneberg, R., Berven, F. S., Myhr, K. M., Ulvik, R. J., & Kvalheim, O. M. (2009). Biomarker discovery in mass spectral profiles by means of selectivity ratio plot. Chemometrics and Intelligent Laboratory Systems, 95(1), 35-48.
#'
#' @return DESCRIOTION REQUIRED
#' @author Leonardo Ramirez-Lopez
#' @examples
#' \dontrun{
#' #' ## GOOD EXAMPLES ARE REQUIRED
#' }
#' @export fit_library

fit_library <- function(
    Xr,
    Yr,
    diss_method,
    k,
    k_diss,
    k_range,
    ws = NULL,
    pls_c,
    pc_selection = list(method = "var", value = 0.01),
    group,
    build_around,
    gh = TRUE,
    center = TRUE,
    scale = FALSE,
    diss_predictors = FALSE,
    metric = "rmse",
    return_best = TRUE,
    pls_max_iter = 1,
    pls_tol = 1e-6,
    documentation = character(),
    ...
) {
  
  call.f <-(match.call())    
  
  
  ## Insert sanity check HERE
  
  if (is.null(colnames(Xr))) {
    stop("column names are mandatory for Xr")
  }
  
  if (missing(group)) {
    group <- as.factor(paste("g", 1:nrow(Xr), sep = ""))
  }
  ## Compute the dissimilarity matrix for Xr
  
  if (class(Yr) != "matrix") {
    Yr <- matrix(Yr, nrow = nrow(Xr))
  }
  
  
  if (missing(k) & missing(k_diss)) 
    stop("Either k or k_diss must be specified")
  
  if (!missing(k)) {
    if (!missing(k_diss))
      stop("Only one of k or k_diss can be specified")  
    if (!is.numeric(k)) {
      stop("k must be a vector of integers")
    }else{
      k <- unique(sort(ceiling(k)))
    }
    k <- sort(k)
  }
  
  
  if (!missing(k_diss)) {
    dtc <- k.diss
    if (missing(k_range))
      stop("If the k_diss argument is used, k_range must be specified")
    if (length(k_range) != 2 | !is.numeric(k_range) | diff(k_range) < 0)
      stop("k_range must be a vector of two 2 integer values specifying the minimum (first value) and the maximum (second value) number of neighbours") 
    k.min <- as.integer(k.range[1])
    k.max <- as.integer(k.range[2])
    if (k.min < 10)
      stop("Minimum number of nearest neighbours allowed is 10")
    if (k.max > nrow(Xr))
      stop("Maximum number of nearest neighbours cannot exceed the number of reference observations")
  } else{
    dtc <- NULL
  }
  
  
  if (pcSel %in% c("cumvar", "var"))
  {
    if (length(pc_selection) == 1) {
      if (pcSel == "cumvar") {
        pc_selection <- list("cumvar", 0.99)
        message(paste("Since the value of the 'pc_selection' argument is missing, the amount of cumulative variance that the components to be retained should explain was automatically set to 0.99 (99%)"))
      }else{
        pc_selection <- list("var", 0.01)
        message(paste("Since the value of the 'pc_selection' argument is missing, the amount of variance that the last component to be retained should explain was automatically set to 0.01 (1%)"))
      }
      names(pc_selection) <- c("method", "value")
    } else {
      if (!is.list(pc_selection))
        stop("The 'pc_selection' argument must be a list in which the first object indicates the selection method and the second object indicates the parameter value of the method. Optionally, instead a list, a character string specifying only the method can be used, in this case the parameter value is set automatically")     
      pc_selection <- list(pc_selection[[1]], pc_selection[[2]])
      names(pc_selection) <- c("method", "value")
      if (!is.numeric(pc_selection$value)) 
        stop("The second object in 'pc_selection' must be a numeric value")
      if (pc_selection$value > 1 | pc_selection$value <= 0) 
        stop(paste("When the method for 'pc_selection' is either 'var' or 'cumvar' the value in 'pc_selection' must be a number larger than 0 and below or equal to 1"))
    }
  }
  
  
  if (pcSel == "manual")
  {
    if (is.list(pc_selection)) {
      if (pc_selection$value < 2) 
        stop(paste("the number of principal components must be an integer value larger than 1")) 
    } else{
      pc_selection <- list(method = "manual", value = 3)
      message(paste("the number of principal components to be retained was set to 3.", "\n", "Note: An user-defined number of principal components to be retained can be specified in the pc_selection argument with a list in which the first object indicates the method 'manual' and the second object indicates the number of principal components"))
    }
  }
  
  
  cat("Computing dissimilarities... \n")
  dsm <- dissimilarity(
    Xr = Xr,
    diss_method = diss_method,
    Yr = Yr,
    center = center,
    scale = scale,
    gh = gh,
    pc_selection = pc_selection,
    return_projection = TRUE,
    ws = ws
  )
  
  dsm <- dsm[!names(dsm) %in% "documentation"]
  sml <- list(diss_method = diss_method)
  dsm <- append(sml, dsm)
  names(dsm)
  
  
  if (any(is.na(Yr))) {
    
    orderf <- function(x, mk, skip = which(is.na(Yr))) {
      ord <- order(x, na.last = TRUE)
      ord <- ord[c(TRUE, !ord[-1] %in% skip)][1:mk]
      ord
    }
    
    if (sum(!is.na(Yr)) < max(k)) {
      stop(paste("the amount of maximum neighbors selected (", max(k), ") is larger than the number of non-NA observations in Yr. Try with less neighbors.", sep = ""))
    }
  }else{
    orderf <- function(x, mk) {
      ord <- order(x, na.last = TRUE)[1:mk]
      ord
    }
  }
  
  ## find the indices of the nearest neighbors
  kidxmat <- apply(dsm$dissimilarity, 
                   MARGIN = 2, 
                   FUN = orderf,
                   mk = max(k))
  
  kdissmat <- sapply(1:ncol(dsm$dissimilarity), 
                     FUN = function(x, y, ..i..) {
                       x[,..i..][y[,..i..]]
                     },
                     x = dsm$dissimilarity,
                     y = kidxmat)
  
  #browser()
  ## for each sample in Xu show what of its
  ## nearest neighbors samples belong to its group
  kidxgrop <- sapply(1:ncol(kidxmat),
                     FUN = function(..i.., kidxmat, group) {
                       gmat <- which(group[..i..] == group)
                       return(!kidxmat[,..i..] %in% gmat)
                     }, 
                     group = group,
                     kidxmat = kidxmat)
  
  nnstats <- lapply(1:length(k),
                    FUN = i_nn_stats,
                    kidxgrop = kidxgrop,
                    kidxmat = kidxmat,
                    Yr = Yr,
                    k = k)
  
  names(nnstats) <- paste("k.", k, sep = "")
  
  ## RMSE/(q"75%" - q"25%")
  ## Similar to the ratio of performance to 
  ## inter-quartile distance (RPIQ)
  # Bellon-Maurel, V., Fernandez-Ahumada, E., 
  # Palagos, B., Roger, J.M., McBratney, A., 2010. 
  # Critical review of chemometric indicators commonly 
  # used for assessing the quality of the prediction of 
  # soil attributes by NIR spectroscopy. TrAC Trends 
  # in Analytical Chemistry, 29(9), pp.1073-1081.
  itq <- abs(sapply(nnstats,
                    FUN = function(x) {
                      x[,"75%"] - x[,"25%"]
                    }))
  
  if (diss_predictors) {
    dssm <- dsm$dissimilarity
  }else{
    dssm <- NULL
    npredictors <- ncol(Xr)
  }
  
  addit <- ifelse(return_best, 1, 0)
  
  pb <- txtProgressBar(min = 0, max = length(k) + addit, char = "-")
  
  minF <- min(pls_c)
  maxF <- max(pls_c)
  sgrid <- expand.grid(minpls = minF:maxF, maxpls = minF:maxF)
  sgrid <- sgrid[sgrid$minpls <= sgrid$maxpls,]
  row.names(sgrid) <- 1:nrow(sgrid)
  
  emgrid <- t(sapply(1:nrow(sgrid),
                     FUN = function(q, x, wv) {
                       wv[x[q, 1]:x[q, 2]] <- 1
                       wv
                     },
                     x = sgrid,
                     wv = rep(0, maxF)))
  emgrid <- emgrid[,minF:maxF]
  
  
  cat("Fitting models... \n")
  
  ## perform the nearest neghbor predictions  
  nnpreds <- sapply(1:length(k), 
                    FUN = .get_all_fits,
                    Xr = Xr, 
                    Yr = Yr, 
                    k = k,
                    min_component = minF, 
                    max_component = maxF, 
                    emgrid = emgrid,
                    scale = scale, 
                    maxiter = pls_max_iter, 
                    tol = pls_tol, 
                    regression = TRUE, 
                    pc_selection = pc_selection,
                    kidxmat = kidxmat,
                    kidxgrop = kidxgrop,
                    dissimilarity_mat = dssm,
                    pb = pb)
  
  ## Organize the results (in nnpreds)
  pparam <- matrix(NA, nrow(emgrid), 4)
  sstats <- function(y, yhat, itqk, pparam) {
    me <- sweep(-yhat, MARGIN = 1, STATS = y, FUN = "+", check.margin = FALSE)
    pparam[,1] <- cor(y, yhat, use = "complete.obs")^2
    pparam[,2] <- colMeans(me^2, na.rm = TRUE)^0.5
    pparam[,3] <- colMeans(me, na.rm = TRUE)
    pparam[,4] <- colMeans(sweep(me, 
                                 MARGIN = 1, 
                                 STATS = itqk, 
                                 FUN = "/", 
                                 check.margin = FALSE)^2, 
                           na.rm = TRUE)^0.5
    return(pparam)
  }
  
  isubset3Row <- function(x1, x2, x3, x4) {
    sargs <- names(match.call())[-1]
    nextEl <- function(..ii..) {
      sapply(sargs, FUN = function(x) get(x)[,..ii..], simplify = FALSE)
    }
    obj <- list(nextElem = nextEl)
  }
  
  itr <- isubset3Row(x1 = nnpreds, x2 = itq)
  
  kpredstats <- function(..k.., 
                         itr,
                         pparam,
                         y) {
    
    ne <- itr$nextElem(..k..)
    statsresults <- sstats(y = y, 
                           yhat = t(matrix(ne$x1, nrow(pparam))), 
                           itqk = ne$x2, 
                           pparam = pparam)  
    statsresults
  }
  
  predperformance <- lapply(1:ncol(nnpreds), 
                            FUN = kpredstats,
                            itr = itr,
                            pparam = pparam,
                            y = Yr)
  
  predperformance <- data.frame(do.call("rbind", predperformance))
  colnames(predperformance) <- c("r2", "rmse", "me", "st.rmse")
  predperformance <- data.frame(minpls = rep(sgrid$minpls, times = length(k)),
                                maxpls = rep(sgrid$maxpls, times = length(k)),
                                k = rep(k, each = nrow(pparam)),
                                predperformance)
  
  
  # ## store results
  # kresults <- data.frame(k = k, 
  #                        r2 = cor(nnpreds , Yr[,], use = "complete.obs")^2,
  #                        rmse = colMeans(me^2, na.rm = TRUE)^0.5,
  #                        rmse.st = colMeans((me/itq)^2, na.rm = TRUE)^0.5,
  #                        check.names = FALSE)
  # plot(kresults$k, kresults$rmse)
  
  # find optinmal parameters
  
  if (metric == c("rmse")) {
    bestp <- predperformance[which.min(predperformance[,metric]),][1,]
  }
  
  if (metric == c("r2")) {
    bestp <- predperformance[which.max(predperformance[,metric]),][1,]
  }
  
  optimalk <- bestp$k
  optimalminpls <- bestp$minpls
  optimalmaxpls <- bestp$maxpls
  
  
  ## Extract the vector of predictions corresponding to the best predictions
  ## (optimal k, optimal pls range)
  plsitemn <- which(sgrid$minpls == optimalminpls & sgrid$maxpls == optimalmaxpls)
  plsitemn <- seq(plsitemn, by = nrow(pparam), length.out = nrow(Xr))
  bestpreds <- itr$nextElem(which(k == optimalk))$x1[plsitemn]
  bestpredsresiduals <- Yr - bestpreds   
  
  # ithbarrio <- ith_subsets(x = Xr,
  #                          y = Yr,
  #                          kindx = kidxmat[1:optimalk,],
  #                          D = dssm)
  
  
  if (return_best) {
    ## build the pls librarby
    # plslib <- sapply(1:nrow(Xr), 
    #                  FUN = final_fits,
    #                  ithsubset = ithbarrio,
    #                  minF = optimalminpls, 
    #                  maxF = optimalmaxpls, 
    #                  scale = group, 
    #                  maxiter = pls_max_iter, 
    #                  tol = pls_tol, 
    #                  pc_selection = pc_selection)
    
    plslib <- foreach(i = 1:nrow(Xr), 
                      .export = c("i_nn_stats",
                                  "ith_pred_subsets",
                                  "ith_subsets",
                                  "ith_subsets_by_group",
                                  ".get_all_fits",
                                  "ith_local_fit"),
                      ithbarrio = ith_subsets(x = Xr, y = Yr, kindx = kidxmat[1:optimalk,], D = dssm)) %dopar%{ 
                        
                        iplslib <- final_fits(
                          ilocalsubset = ithbarrio,
                          min_component = optimalminpls, 
                          max_component = optimalmaxpls, 
                          scale = group, 
                          maxiter = pls_max_iter, 
                          tol = pls_tol
                        )
                        
                        # print(format(object.size(iplslib), "b"))
                        iplslib
                      }
    
    plslib <- do.call("cbind", plslib)
    
    setTxtProgressBar(pb, length(k) + 1)
    
    
    plslib <- t(plslib)
    
    namesk <- NULL
    if (diss_predictors) {
      npredictors <- optimalk + ncol(Xr)
      namesk <- paste("k", 1:optimalk, sep = "")
    }
    
    if (group) {
      xscale <- plslib[ ,-c(1:((3 * npredictors) + 1))]
      plslib <- plslib[ ,c(1:((3 * npredictors) + 1))]
    }else{
      xscale <- matrix(1, nrow(plslib), npredictors)
    }
    plsvips <- plslib[ ,-c(1:(npredictors + 1), (ncol(plslib) - npredictors + 1):ncol(plslib))]
    plssratios <- plslib[ ,c((ncol(plslib) - npredictors + 1):ncol(plslib))]
    plslib <- plslib[ ,c(1:(npredictors + 1))]
    
    colnames(plslib) <- c("b0", namesk, colnames(Xr))
    colnames(xscale) <- c(namesk, colnames(Xr))
    colnames(plssratios) <- colnames(plsvips) <- c(namesk, colnames(Xr))
    
    if (diss_predictors) {
      bs <- list(B0 = plslib[,1],
                 B = plslib[,colnames(plslib) %in% colnames(Xr)],
                 Bk = plslib[,colnames(plslib) %in% namesk])
      
    }else{
      bs <- list(B0 = plslib[,1],
                 B = plslib[,colnames(plslib) %in% colnames(Xr)])
    }
    
    if (center) {
      center <- colMeans(Xr)
    }else{
      center <- rep(0, ncol(Xr))
    }
    
    if (group) {
      gscale <- colSds(Xr)
    }else{
      gscale <- rep(1, ncol(Xr))
    }
    
    iscale <- list(centre = center,
                   scale = gscale,
                   local.x.scale = xscale[ ,colnames(xscale) %in% colnames(Xr)])
    if (diss_predictors) {
      iscale$local.diss.scale <- xscale[ ,colnames(xscale) %in% namesk]
    }
    
    if (diss_method %in% c("pca", "pls")) {
      fresults <- list(dissimilatity = dsm[!names(dsm) %in% "gh"],
                       gh = dsm$gh, 
                       results = predperformance,
                       best = bestp,
                       sel.param = list(optimalk = optimalk,
                                        optimal.factors = c(min.pls = bestp$minpls, 
                                                            max.pls = bestp$maxpls)),
                       residuals = bestpredsresiduals,
                       functionlibrary = bs,
                       functionvips = plsvips,
                       functionselectivityrs = plssratios,
                       scale =  iscale,
                       yu.nnstats = nnstats,
                       documentation = documentation)
    }else{
      fresults <- list(dissimilatity = dsm[!names(dsm) %in% "gh"],
                       gh = dsm$gh, 
                       results = predperformance,
                       best = bestp,
                       sel.param = list(optimalk = optimalk,
                                        optimal.factors = c(min.pls = bestp$minpls, 
                                                            max.pls = bestp$maxpls)),
                       residuals = bestpredsresiduals,
                       functionlibrary = bs,
                       functionvips = plsvips,
                       functionselectivityrs = plssratios,
                       scale = iscale,
                       Xr = Xr,
                       yu.nnstats = nnstats,
                       documentation = documentation)
    }
    
    if (!is.null(rownames(Xr))) {
      namesrows <- rownames(Xr)
    }else{
      namesrows <- 1:nrow(Xr)
    }
    
    fresults$yu.nnstats <- lapply(fresults$yu.nnstats, 
                                  FUN = function(x, nms) {
                                    rownames(x) <- nms
                                    x}, 
                                  nms = namesrows)
    
    rownames(fresults$functionvips) <- 
      rownames(fresults$functionselectivityrs) <- 
      rownames(fresults$residuals) <- 
      rownames(fresults$functionlibrary$B) <- 
      names(fresults$functionlibrary$B0) <- namesrows
    if (!is.null(fresults$functionlibrary$Bk))
      rownames(fresults$functionlibrary$Bk) <- namesrows
    class(fresults) <- c("funlib","list")
  }else{
    
    fresults <- list(dissimilatity = dsm[!names(dsm) %in% "gh"],
                     gh = dsm$gh, 
                     results = predperformance,
                     best = bestp,
                     residuals = bestpredsresiduals,
                     yu.nnstats = nnstats,
                     documentation = documentation)
    
    if (!is.null(rownames(Xr))) {
      namesrows <- rownames(Xr)
    }else{
      namesrows <- 1:nrow(Xr)
    }
    
    fresults$yu.nnstats <- lapply(fresults$yu.nnstats, 
                                  FUN = function(x, nms) {
                                    rownames(x) <- nms
                                    x}, 
                                  nms = namesrows)
    
    rownames(fresults$residuals) <- namesrows
    class(fresults) <- c("validationfunlib","list")
  }
  
  cat("\nDone!")
  attr(fresults, "call") <- call.f
  return(fresults)
}
