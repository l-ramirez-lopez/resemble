#' @title A function for memory-based learning (mbl)
#' @description
#' \loadmathjax
#' This function is implemented for memory-based learning (a.k.a.
#' instance-based learning or local regression) which is a non-linear lazy
#' learning approach for predicting a given response variable from a set of
#' predictor variables. For each observation in a prediction set, a specific
#' local regression is carried out based on a subset of similar observations
#' (nearest neighbors) selected from a reference set. The local model is
#' then used to predict the response value of the target (prediction)
#' observation. Therefore this function does not yield a global
#' regression model.
#' @usage
#' mbl(Xr, Yr, Xu, Yu = NULL, k, k_diss, k_range, spike = NULL,
#'     method = local_fit_wapls(min_pls_c = 3, max_pls_c = min(dim(Xr), 15)),
#'     diss_method = "pca", diss_usage = "predictors", gh = TRUE,
#'     pc_selection = list(method = "opc", value = min(dim(Xr), 40)),
#'     control = mbl_control(), group = NULL, center = TRUE, scale = FALSE,
#'     verbose = TRUE, documentation = character(), seed = NULL, ...)
#'
#' @param Xr a matrix of predictor variables of the reference data
#' (observations in rows and variables in columns).
#' @param Yr a numeric matrix of one column containing the values of the
#' response variable corresponding to the reference data.
#' @param Xu a matrix of predictor variables of the data to be predicted
#' (observations in rows and variables in columns).
#' @param Yu an optional matrix of one column containing the values of the
#' response variable corresponding to the data to be predicted. Default is
#' \code{NULL}.
#' @param k a vector of integers specifying the sequence of k-nearest
#' neighbors to be tested. Either \code{k} or \code{k_diss} must be specified.
#' This vector will be automatically sorted into ascending order. If
#' non-integer numbers are passed, they will be coerced to the next upper
#' integers.
#' @param k_diss a numeric vector specifying the sequence of dissimilarity
#' thresholds to be tested for the selection of the nearest neighbors found in
#' \code{Xr} around each observation in \code{Xu}. These thresholds depend on
#' the corresponding dissimilarity measure specified in the object passed to
#' \code{control}. Either \code{k} or \code{k_diss} must be specified.
#' @param k_range an integer vector of length 2 which specifies the minimum
#' (first value) and the maximum (second value) number of neighbors to be
#' retained when the \code{k_diss} is given.
#' @param spike an integer vector (with positive and/or negative values) indicating
#' the indices of observations in \code{Xr} that must be either be forced into
#' or avoided in the neighborhoods of every \code{Xu} observation. Default is
#' \code{NULL} (i.e. no observations are forced or avoided). Note
#' that this argument is not intended for increasing or reducing the neighborhood
#'  size which is only controlled by \code{k} or \code{k_diss} and \code{k_range}.
#' By forcing observations into the neighborhood, some of the farthest
#' observations may be forced out of the neighborhood. In contrast, by avoiding
#' observations in the neighborhood,  some of farthest
#' observations may be included into the neighborhood. See details.
#' @param method an object of class \code{\link{local_fit}} which indicates the
#' type of regression to conduct at each local segment as well as additional
#' parameters affecting this regression. See \code{\link{local_fit}} function.
#' @param diss_method a character string indicating the spectral dissimilarity
#' metric to be used in the selection of the nearest neighbors of each
#' observation. Options are:
#' \itemize{
#'        \item{\code{"pca"} (Default): Mahalanobis distance
#'        computed on the matrix of scores of a Principal Component (PC)
#'        projection of \code{Xr} and \code{Xu}. PC projection is done using the
#'        singular value decomposition (SVD) algorithm.
#'        See \code{\link{ortho_diss}} function.}

#'        \item{\code{"pca.nipals"}: Mahalanobis distance
#'        computed on the matrix of scores of a Principal Component (PC)
#'        projection of \code{Xr} and \code{Xu}. PC projection is done using the
#'        non-linear iterative partial least squares (nipals) algorithm.
#'        See \code{\link{ortho_diss}} function.}

#'        \item{\code{"pls"}: Mahalanobis distance
#'        computed on the matrix of scores of a partial least squares projection
#'        of \code{Xr} and \code{Xu}. In this case, \code{Yr} is always
#'        required. See \code{\link{ortho_diss}} function.}

#'        \item{\code{"cor"}: correlation coefficient
#'        between observations. See \code{\link{cor_diss}} function.}

#'        \item{\code{"euclid"}: Euclidean distance
#'        between observations. See \code{\link{f_diss}} function.}

#'        \item{\code{"cosine"}: Cosine distance
#'        between observations. See \code{\link{f_diss}} function.}

#'        \item{\code{"sid"}: spectral information divergence between
#'        observations. See \code{\link{sid}} function.}
#'        }
#' Alternatively, a matrix of dissimilarities can also be passed to this
#' argument. This matrix is supposed to be a user-defined matrix
#' representing the dissimilarities between observations in \code{Xr} and
#' \code{Xu}. When \code{diss_usage = "predictors"}, this matrix must be squared
#' (derived from a matrix of the form \code{rbind(Xr, Xu)}) for which the
#' diagonal values are zeros (since the dissimilarity between an object and
#' itself must be 0). On the other hand, if \code{diss_usage} is set to either
#' \code{"weights"} or \code{"none"}, it must be a matrix representing the
#' dissimilarity of each observation in \code{Xu} to each observation in
#' \code{Xr}. The number of columns of the input matrix must be equal to the
#' number of rows in \code{Xu} and the number of rows equal to the number of
#' rows in \code{Xr}.
#' @param diss_usage a character string specifying how the dissimilarity
#' information shall be used. The possible options are: \code{"predictors"},
#' \code{"weights"} and \code{"none"} (see details below).
#' Default is \code{"predictors"}.
#' @param control a list created with the \code{\link{mbl_control}} function
#' which contains additional parameters that control some few aspects of the
#' \code{mbl} function (cross-validation, parameter tuning, etc).
#' The default list is as returned by \code{mbl_control()}.
#' See the \code{\link{mbl_control}} function for more details.
#' @param gh a logical indicating if the global Mahalanobis distance (in the pls
#' score space) between each observation and the pls mean (centre) must be
#' computed. This metric is known as the GH distance in the literature. Note
#' that this computation is based on the number of pls components determined by
#' using the \code{pc_selection} argument. See details.
#' @param pc_selection a list of length 2 used for the computation of GH (if
#' \code{gh = TRUE}) as well as in the computation of the dissimilarity methods
#' based on \code{\link{ortho_diss}} (i.e. when \code{diss_method} is one of:
#' \code{"pca"}, \code{"pca.nipals"} or \code{"pls"}) or when \code{gh = TRUE}.
#' This argument is used for optimizing the number of components (principal
#' components or pls factors) to be retained for dissimilarity/distance
#' computation purposes only (i.e not for regression).
#' This list must contain two elements in the following order:
#' \code{method} (a character indicating the method for selecting the number of
#' components) and \code{value} (a numerical value that complements the selected
#' method). The methods available are:
#' \itemize{
#'        \item{\code{"opc"}: optimized principal component selection based
#'        on Ramirez-Lopez et al. (2013a, 2013b). The optimal number of
#'        components (of set of observations) is the one for which its distance
#'        matrix minimizes the differences between the \code{Yr} value of each
#'        observation and the \code{Yr} value of its closest observation. In
#'        this case \code{value} must be a value (larger than 0 and
#'        below the minimum dimension of \code{Xr} or \code{Xr} and \code{Xu}
#'        combined) indicating the maximum
#'        number of principal components to be tested. See the
#'        \code{\link{ortho_projection}} function for more details.}

#'        \item{\code{"cumvar"}: selection of the principal components based
#'        on a given cumulative amount of explained variance. In this case,
#'        \code{value} must be a value (larger than 0 and below or equal to 1)
#'        indicating the minimum amount of cumulative variance that the
#'        combination of retained components should explain.}

#'        \item{\code{"var"}: selection of the principal components based
#'        on a given amount of explained variance. In this case,
#'        \code{value} must be a value (larger than 0 and below or equal to 1)
#'        indicating the minimum amount of variance that a single component
#'        should explain in order to be retained.}

#'        \item{\code{"manual"}: for manually specifying a fix number of
#'        principal components. In this case, \code{value} must be a value
#'        (larger than 0 and below the minimum dimension of \code{Xr} or
#'        \code{Xr} and \code{Xu} combined).
#'        indicating the minimum amount of variance that a component should
#'        explain in order to be retained.}
#'        }
#' The list
#' \code{list(method = "opc", value = min(dim(Xr), 40))} is the default.
#' Optionally, the \code{pc_selection} argument admits \code{"opc"} or
#' \code{"cumvar"} or \code{"var"} or \code{"manual"} as a single character
#' string. In such a case the default \code{"value"} when either \code{"opc"} or
#' \code{"manual"} are used is 40. When \code{"cumvar"} is used the default
#' \code{"value"} is set to 0.99 and when \code{"var"} is used, the default
#' \code{"value"} is set to 0.01.
#' @param group an optional factor (or character vector vector
#' that can be coerced to \code{\link[base]{factor}} by \code{as.factor}) that
#' assigns a group/class label to each observation in \code{Xr}
#' (e.g. groups can be given by spectra collected from the same batch of
#' measurements, from the same observation, from observations with very similar
#' origin, etc). This is taken into account for internal leave-group-out cross
#' validation for pls tuning (factor optimization) to avoid pseudo-replication.
#' When one observation is selected for cross-validation, all observations of
#' the same group are removed together and assigned to validation. The length
#' of the vector must be equal to the number of observations in the
#' reference/training set (i.e. \code{nrow(Xr)}). See details.
#' @param center a logical if the predictor variables must be centred at each
#' local segment (before regression). In addition, if \code{TRUE}, \code{Xr}
#' and \code{Xu} will be centred for  dissimilarity computations.
#' @param scale a logical indicating if the predictor variables must be scaled
#' to unit variance at each local segment (before regression). In addition, if
#' \code{TRUE}, \code{Xr} and \code{Xu} will be scaled for  dissimilarity
#' computations.
#' @param verbose a logical indicating whether or not to print a progress bar
#' for each observation to be predicted. Default is \code{TRUE}. Note: In case
#' parallel processing is used, these progress bars will not be printed.
#' @param seed an integer value containing the random number generator (RNG)
#' state for random number generation. This argument can be used for
#' reproducibility purposes (for random sampling) in the cross-validation
#' results. Default is \code{NULL}, i.e. no RNG is applied.
#' @param documentation an optional character string that can be used to
#' describe anything related to the \code{mbl} call (e.g. description of the
#' input data). Default: \code{character()}. NOTE: his is an experimental
#' argument.
#' @param ... further arguments to be passed to the \code{\link{dissimilarity}}
#' function. See details.
#'
#' @details
#' The argument \code{spike} can be used to indicate what reference observations
#' in \code{Xr} must be kept in the neighborhood of every single \code{Xu}
#' observation. If a vector of length \mjeqn{m}{m} is passed to this argument,
#' this means that the \mjeqn{m}{m} original neighbors with the largest
#' dissimilarities to the target observations will be forced out of the
#' neighborhood. Spiking might be useful in cases where
#' some reference observations are known to be somehow related to the ones in
#' \code{Xu} and therefore might be relevant for fitting the local models. See
#' Guerrero et al. (2010) for an example on the benefits of spiking.
#'
#' The \code{mbl} function uses the \code{\link{dissimilarity}} function to
#' compute the dissimilarities between \code{Xr} and \code{Xu}. The dissimilarity
#' method to be used is specified in the \code{diss_method} argument.
#' Arguments to \code{\link{dissimilarity}} as well as further arguments to the
#' functions used inside \code{\link{dissimilarity}}
#' (i.e. \code{\link{ortho_diss}} \code{\link{cor_diss}} \code{\link{f_diss}}
#' \code{\link{sid}}) can be passed to those functions by using \code{...}.
#'
#' The \code{diss_usage} argument is used to specify whether the dissimilarity
#' information must be used within the local regressions and, if so, how.
#' When \code{diss_usage = "predictors"} the local (square symmetric)
#' dissimilarity matrix corresponding the selected neighborhood is used as
#' source of additional predictors (i.e the columns of this local matrix are
#' treated as predictor variables). In some cases this results in an improvement
#' of the prediction performance (Ramirez-Lopez et al., 2013a).
#' If \code{diss_usage = "weights"}, the neighbors of the query point
#' (\mjeqn{xu_{j}}{xu_j}) are weighted according to their dissimilarity to
#' \mjeqn{xu_{j}}{xu_j} before carrying out each local regression. The following
#' tricubic function (Cleveland and Delvin, 1988; Naes et al., 1990) is used for
#' computing the final weights based on the measured dissimilarities:
#'
#' \mjdeqn{W_{j}  =  (1 - v^{3})^{3}}{W_j  =  (1 - v^3)^3}
#'
#' where if \mjeqn{{xr_{i} \in }}{xr_i in} neighbors of \mjeqn{xu_{j}}{xu_j}:
#'
#' \mjdeqn{v_{j}(xu_{j})  =  d(xr_{i}, xu_{j})}{v_j(xu_j)  =  d(xr_i, xu_j)}
#'
#' otherwise:
#'
#' \mjdeqn{v_{j}(xu_{j})  =  0}{v_j(xu_j)  =  0}
#'
#' In the above formulas \mjeqn{d(xr_{i}, xu_{j})}{d(xr_i, xu_j)} represents the
#' dissimilarity between the query point and each object in \mjeqn{Xr}{Xr}.
#' When \code{diss_usage = "none"} is chosen the dissimilarity information is
#' not used.
#'
#' The global Mahalanobis distance (a.k.a GH) is computed based on the scores
#' of a pls projection. A pls projection model is built with for \code{\{Yr\}, \{Xr\}}
#' and this model is used to obtain the pls scores of the \code{Xu}
#' observations. The Mahalanobis distance between each \code{Xu} observation in
#' (the pls space) and the centre of \code{Xr} is then computed. The number of
#' pls components is optimized based on the parameters passed to the
#' \code{pc_selection} argument. In addition, the \code{mbl} function also
#' reports the GH distance for the observations in \code{Xr}.
#'
#' Some aspects of the mbl process, such as the type of internal validation,
#' parameter tuning, what extra objects to return, permission for parallel
#' execution, prediction limits, etc, can be specified by using the
#' \code{\link{mbl_control}} function.
#'
#' By using the \code{group} argument one can specify groups of observations
#' that have something in common (e.g. observations with very similar origin).
#' The purpose of \code{group} is to avoid biased cross-validation results due
#' to pseudo-replication. This argument allows to select calibration points
#' that are independent from the validation ones. In this regard, when
#' \code{validation_type = "local_cv"} (used in \code{\link{mbl_control}}
#' function), then the \code{p} argument refers to the percentage of groups of
#' observations (rather than single observations) to be retained in each
#' sampling iteration at each local segment.
#'
#' @return a \code{list} of class \code{mbl} with the following components
#' (sorted either by \code{k} or \code{k_diss}):
#'
#' \itemize{
#'  \item{\code{call}: the call to mbl.}
#'  \item{\code{cntrl_param}: the list with the control parameters passed to
#'  control.}
#'  \item{\code{Xu_neighbors}: a list containing two elements: a matrix of
#'  \code{Xr} indices corresponding to the neighbors of \code{Xu} and a matrix
#'  of dissimilarities between each \code{Xu} observation and its corresponding
#'  neighbor in \code{Xr}.}
#'  \item{\code{dissimilarities}: a list with the method used to obtain the
#'  dissimilarity matrices and the dissimilarity matrix corresponding to
#'  \mjeqn{D(Xr, Xu)}{D(Xr, Xu)}. This object is returned only if the
#'  \code{return_dissimilarity} argument in the \code{control} list was set
#'  to \code{TRUE}.}
#'  \item{\code{n_predictions}: the total number of observations predicted.}
#'  \item{\code{gh}: if \code{gh = TRUE}, a list containing the global
#'  Mahalanobis distance values for the observations in \code{Xr} and \code{Xu}
#'  as well as the results of the global pls projection object used to obtain
#'  the GH values.}
#'  \item{\code{validation_results}: a list of validation results for
#'  "local cross validation" (returned if the \code{validation_type} in
#'  \code{control} list was set to \code{"local_cv"}),
#'  "nearest neighbor validation" (returned if the \code{validation_type}
#'  in \code{control} list was set to \code{"NNv"}) and
#'  "Yu prediction statistics" (returned  if \code{Yu} was supplied).}``
#'  \item{\code{results}: a list of data tables containing the results of the
#'  predictions for each either \code{k} or \code{k_diss}. Each data table
#'  contains the following columns:}
#'  \itemize{
#'    \item{\code{o_index}: The index of the predicted observation.}
#'    \item{\code{k_diss}: This column is only output if the \code{k_diss}
#'    argument is used. It indicates the corresponding dissimilarity threshold
#'    for selecting the neighbors.}
#'    \item{\code{k_original}: This column is only output if the \code{k_diss}
#'    argument is used. It indicates the number of neighbors that were originally
#'    found when the given dissimilarity threshold is used.}
#'    \item{\code{k}: This column indicates the final number of neighbors
#'    used.}
#'    \item{\code{npls}: This column is only output if the \code{pls}
#'    regression method was used. It indicates the final number of pls
#'    components used.}
#'    \item{\code{min_pls}: This column is only output if \code{wapls}
#'    regression method was used. It indicates the final number of minimum pls
#'    components used. If no optimization was set, it retrieves the original
#'    minimum pls components passed to the \code{method} argument.}
#'    \item{\code{max_pls}: This column is only output if the \code{wapls}
#'    regression method was used. It indicates the final number of maximum pls
#'    components used. If no optimization was set, it retrieves the original
#'    maximum pls components passed to the \code{method} argument.}
#'    \item{\code{yu_obs}: The input values given in \code{Yu} (the response
#'    variable corresponding to the data to be predicted). If \code{Yu = NULL},
#'    then \code{NA}s are retrieved.}
#'    \item{\code{pred}: The predicted Yu values.}
#'    \item{\code{yr_min_obs}: The minimum reference value (of the response
#'    variable) in the neighborhood.}
#'    \item{\code{yr_max_obs}: The maximum reference value (of the response
#'    variable) in the neighborhood.}
#'    \item{\code{index_nearest_in_Xr}: The index of the nearest neighbor found
#'    in \code{Xr}.}
#'    \item{\code{index_farthest_in_Xr}: The index of the farthest neighbor
#'    found in \code{Xr}.}
#'    \item{\code{y_nearest}: The reference value (\code{Yr}) corresponding to
#'    the nearest neighbor found in \code{Xr}.}
#'    \item{\code{y_nearest_pred}: This column is only output if the
#'    validation method in the object passed to \code{control} was set to
#'    \code{"NNv"}. It represents the predicted value of the nearest neighbor
#'    observation found in \code{Xr}. This prediction come from model fitted
#'    with the remaining observations in the neighborhood of the target
#'    observation in \code{Xu}.}
#'    \item{\code{loc_rmse_cv}: This column is only output if the validation
#'    method in the object passed to \code{control} was set to
#'    \code{'local_cv'}. It represents the RMSE of the cross-validation
#'    computed for the neighborhood of the target observation in \code{Xu}.}
#'    \item{\code{loc_st_rmse_cv}: This column is only output if the
#'    validation method in the object passed to \code{control} was set to
#'    \code{'local_cv'}. It represents the standardized RMSE of the
#'    cross-validation computed for the neighborhood of the target observation
#'    in \code{Xu}.}
#'    \item{\code{dist_nearest}: The distance to the nearest neighbor.}
#'    \item{\code{dist_farthest}: The distance to the farthest neighbor.}
#'    \item{\code{loc_n_components}: This column is only output if the
#'    dissimilarity method used is one of \code{"pca"}, \code{"pca.nipals"} or
#'    \code{"pls"} and in addition the dissimilarities are requested to be
#'    computed locally by passing \code{.local = TRUE} to the \code{mbl}
#'    function.
#'    See \code{.local} argument in the \code{\link{ortho_diss}} function.}
#'    }
#'  \item{\code{seed}: a value mirroring the one passed to seed.}
#'  \item{\code{documentation}: a character string mirroring the one provided
#'  in the \code{documentation} argument.}
#'  }
#' When the \code{k_diss} argument is used, the printed results show a table
#' with a column named '\code{p_bounded}. It represents the percentage of
#' observations for which the neighbors selected by the given dissimilarity
#' threshold were outside the boundaries specified in the \code{k_range}
#' argument.
#' @author \href{https://orcid.org/0000-0002-5369-5120}{Leonardo Ramirez-Lopez}
#' and Antoine Stevens
#' @references
#' Cleveland, W. S., and Devlin, S. J. 1988. Locally weighted regression: an
#' approach to regression analysis by local fitting. Journal of the American
#' Statistical Association, 83, 596-610.
#'
#' Guerrero, C., Zornoza, R., Gómez, I., Mataix-Beneyto, J. 2010. Spiking of
#' NIR regional models using observations from target sites: Effect of model
#' size on prediction accuracy. Geoderma, 158(1-2), 66-77.
#'
#' Naes, T., Isaksson, T., Kowalski, B. 1990. Locally weighted regression and
#' scatter correction for near-infrared reflectance data. Analytical Chemistry
#' 62, 664-673.
#'
#' Ramirez-Lopez, L., Behrens, T., Schmidt, K., Stevens, A., Dematte, J.A.M.,
#' Scholten, T. 2013a. The spectrum-based learner: A new local approach for
#' modeling soil vis-NIR spectra of complex data sets. Geoderma 195-196,
#' 268-279.
#'
#' Ramirez-Lopez, L., Behrens, T., Schmidt, K., Viscarra Rossel, R., Dematte,
#' J. A. M.,  Scholten, T. 2013b. Distance and similarity-search metrics for
#' use with soil vis-NIR spectra. Geoderma 199, 43-53.
#'
#' Rasmussen, C.E., Williams, C.K. Gaussian Processes for Machine Learning.
#' Massachusetts Institute of Technology: MIT-Press, 2006.
#'
#' Shenk, J., Westerhaus, M., and Berzaghi, P. 1997. Investigation of a LOCAL
#' calibration procedure for near infrared instruments. Journal of Near
#' Infrared Spectroscopy, 5, 223-232.
#'
#' @seealso \code{\link{mbl_control}}, \code{\link{f_diss}},
#' \code{\link{cor_diss}}, \code{\link{sid}}, \code{\link{ortho_diss}},
#' \code{\link{search_neighbors}},  \code{\link{local_fit}}
#' @examples
#' \donttest{
#' library(prospectr)
#' data(NIRsoil)
#'
#' # Proprocess the data using detrend plus first derivative with Savitzky and
#' # Golay smoothing filter
#' sg_det <- savitzkyGolay(
#'   detrend(NIRsoil$spc,
#'     wav = as.numeric(colnames(NIRsoil$spc))
#'   ),
#'   m = 1,
#'   p = 1,
#'   w = 7
#' )
#'
#' NIRsoil$spc_pr <- sg_det
#'
#' # split into training and testing sets
#' test_x <- NIRsoil$spc_pr[NIRsoil$train == 0 & !is.na(NIRsoil$CEC), ]
#' test_y <- NIRsoil$CEC[NIRsoil$train == 0 & !is.na(NIRsoil$CEC)]
#'
#' train_y <- NIRsoil$CEC[NIRsoil$train == 1 & !is.na(NIRsoil$CEC)]
#' train_x <- NIRsoil$spc_pr[NIRsoil$train == 1 & !is.na(NIRsoil$CEC), ]
#'
#' # Example 1
#' # A mbl implemented in Ramirez-Lopez et al. (2013,
#' # the spectrum-based learner)
#' # Example 1.1
#' # An exmaple where Yu is supposed to be unknown, but the Xu
#' # (spectral variables) are known
#' my_control <- mbl_control(validation_type = "NNv")
#'
#' ## The neighborhood sizes to test
#' ks <- seq(40, 140, by = 20)
#'
#' sbl <- mbl(
#'   Xr = train_x,
#'   Yr = train_y,
#'   Xu = test_x,
#'   k = ks,
#'   method = local_fit_gpr(),
#'   control = my_control,
#'   scale = TRUE
#' )
#' sbl
#' plot(sbl)
#' get_predictions(sbl)
#'
#' # Example 1.2
#' # If Yu is actually known...
#' sbl_2 <- mbl(
#'   Xr = train_x,
#'   Yr = train_y,
#'   Xu = test_x,
#'   Yu = test_y,
#'   k = ks,
#'   method = local_fit_gpr(),
#'   control = my_control
#' )
#' sbl_2
#' plot(sbl_2)
#'
#' # Example 2
#' # the LOCAL algorithm (Shenk et al., 1997)
#' local_algorithm <- mbl(
#'   Xr = train_x,
#'   Yr = train_y,
#'   Xu = test_x,
#'   Yu = test_y,
#'   k = ks,
#'   method = local_fit_wapls(min_pls_c = 3, max_pls_c = 15),
#'   diss_method = "cor",
#'   diss_usage = "none",
#'   control = my_control
#' )
#' local_algorithm
#' plot(local_algorithm)
#'
#' # Example 3
#' # A variation of the LOCAL algorithm (using the optimized pc
#' # dissmilarity matrix) and dissimilarity matrix as source of
#' # additional preditors
#' local_algorithm_2 <- mbl(
#'   Xr = train_x,
#'   Yr = train_y,
#'   Xu = test_x,
#'   Yu = test_y,
#'   k = ks,
#'   method = local_fit_wapls(min_pls_c = 3, max_pls_c = 15),
#'   diss_method = "pca",
#'   diss_usage = "predictors",
#'   control = my_control
#' )
#' local_algorithm_2
#' plot(local_algorithm_2)
#'
#' # Example 4
#' # Running the mbl function in parallel with example 2
#'
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
#' local_algorithm_par <- mbl(
#'   Xr = train_x,
#'   Yr = train_y,
#'   Xu = test_x,
#'   Yu = test_y,
#'   k = ks,
#'   method = local_fit_wapls(min_pls_c = 3, max_pls_c = 15),
#'   diss_method = "cor",
#'   diss_usage = "none",
#'   control = my_control
#' )
#' local_algorithm_par
#'
#' registerDoSEQ()
#' try(stopCluster(clust))
#'
#' # Example 5
#' # Using local pls distances
#' with_local_diss <- mbl(
#'   Xr = train_x,
#'   Yr = train_y,
#'   Xu = test_x,
#'   Yu = test_y,
#'   k = ks,
#'   method = local_fit_wapls(min_pls_c = 3, max_pls_c = 15),
#'   diss_method = "pls",
#'   diss_usage = "predictors",
#'   control = my_control,
#'   .local = TRUE,
#'   pre_k = 150,
#' )
#' with_local_diss
#' plot(with_local_diss)
#' }
#' @export

######################################################################
# resemble
# Copyright (C) 2014 Leonardo Ramirez-Lopez and Antoine Stevens
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
######################################################################

## History:
## 09.03.2014 Leo     Doc examples  were formated with a max. line width
## 13.03.2014 Antoine The explanation of the cores argument was modified
## 23.04.2014 Leo     Added default variable names when they are missing
## 08.09.2014 Leo     A bug related with the computations of the weights
##                    for wapls2 was fixed
## 23.09.2014 Leo     A bug that prevented the mbl function of using
##                    the 'dissimilarityM' argument was fixed
## 03.10.2014 Antoine Fix bug when scale = T and add allow_parallel argument
## 12.10.2014 Leo     noise_variance was missing in the locFit function used for
##                    the nearest neighbor validation
## 16.11.2015 Leo     Now the scale argument for gaussian process is a
##                    indicates that both x and y variables must be scaled
##                    to zero mean and unit variance. Before it only the x
##                    variables were scaled to unit variance
## 18.11.2015 Leo     The mbl examples were corrected (unnecessary arguments
##                    were deleted)
## 01.12.2015 Leo     The wapls2 was removed from the options of methods
## 10.12.2015 Leo     The locFit function has been renamed to locFitnpred
##                    and now it always performs a prediction.
## 10.12.2015 Leo     Several redundant/repaeated sanity checks (ocurring
##                    when combining functions) were deactivated. This does not
##                    impact the finaly sanity checks of the overall mbl
##                    function.
## 11.12.2015 Leo     A bug when the k_diss argument was used was corrected.
##                    The results about the percentage of observations that were
##                    bounded by k_range was not not correct.
## 11.12.2015 Leo     The explanation of the output variables in the results
##                    element of the mbl objects was extended. The rep variable
##                    is not output anymore in the results element.
## 03.01.2016 Leo     Now it is possible to optimize the max and min pls
##                    components of wapls1
## 04.02.2016 Leo     An extrange bug was fixed. The object pred_obs
##                    (in the parallel loop) had a variable named pls_c
##                    (pred_obs$pls_c). When when method = "gpr" was used,
##                    and mbl was runing in parallel it retrieved and error
##                    saying that pls_c was missing!!! This was perhaps due to
##                    the fact that the pls_c  was variable (in pred_obs) and
##                    an argument.
## 16.02.2016 Leo     Bug fixed. It caused the mbl function to return an error
##                    (sometimes) when the group argument was used together
##                    with local cross-validation. The errors occurred when
##                    groups containing very few observations (e.g. 1 or 2) were used.
## 09.03.2018 Leo     A new output (XuneighborList) has been added. It was
##                    requested by Eva Ampe and Miriam de Winter.
## 16.05.2018 Leo     A parameter called documentation has been added.
## 21.06.2020 Leo     - pls.max.iter, pls.tol and noise.v were moved to mbl from
##                      mbl_control()
##                    - Argument scaled (from mbl_control()) renamed to .scale
##                      and moved to mbl
##                    - new arguments: gh and spike
##                    - order of the Yr, Xr, Yu and Xu arguments has changed to
##                      Xr, Yr, Xu and Yu
##                    - input type for the argument method has changed.
##                      Previously it received a character string  indicating
##                      the type of local regression (i.e. "pls",
##                      "wapls1" or "gpr"). Now it receives an object of class
##                      local_fit which is output by the new local_fit functions.
##                    - dissimilarityM has been deprecated. It was used to pass
##                      a dissimilarity matrix computed outside the mbl
##                      function. This can be done now with the new argument
##                      diss_method of mbl which was previously named "sm" and
##                      it was in mbl_control()
##                    - the warning message coming from the foreach loop about no
##                      parallel backend registered is now avoided by checking
##                      first if there is any parallel backend registered
## 22.06.2020 Leo     - Updated examples

mbl <- function(Xr, Yr, Xu = NULL, Yu = NULL,
                k, k_diss, k_range,
                spike = NULL,
                method = local_fit_wapls(
                  min_pls_c = 3,
                  max_pls_c = min(dim(Xr), 15)
                ),
                diss_method = "pca",
                diss_usage = "predictors",
                gh = TRUE,
                pc_selection = list(
                  method = "opc",
                  value = min(dim(Xr), 40)
                ),
                control = mbl_control(),
                group = NULL,
                center = TRUE,
                scale = FALSE,
                verbose = TRUE,
                documentation = character(),
                seed = NULL,
                ...) {
  f_call <- match.call()

  "%mydo%" <- get("%do%")
  if (control$allow_parallel & getDoParRegistered()) {
    "%mydo%" <- get("%dopar%")
  }
  if (!is.logical(verbose)) {
    stop("'verbose' must be logical")
  }

  if (missing(k)) {
    k <- NULL
  }

  if (missing(k_diss)) {
    k_diss <- NULL
  }

  if (missing(k_range)) {
    k_range <- NULL
  }

  input_dots <- list(...)
  ini_cntrl <- control
  ortho_diss_methods <- c("pca", "pca.nipals", "pls")

  if (".local" %in% names(input_dots)) {
    if (isTRUE(input_dots$.local)) {
      if (!"pre_k" %in% names(input_dots)) {
        stop("When '.local = TRUE', argument 'pre_k' needs to be provided. See ortho_diss documentation")
      }
      if (!is.null(k)) {
        if (input_dots$pre_k < max(k)) {
          stop("'k' cannot be larger than 'pre_k'")
        }
      }
    }
  }

  # Sanity checks
  if (!is.logical(center)) {
    stop("'center' argument must be logical")
  }

  if (!is.logical(scale)) {
    stop("'scale' argument must be logical")
  }

  
  if (missing(Xu)) {
    stop("Xu is missing")
  }
  
  if (!is.null(Xu)) {
    if (ncol(Xr) != ncol(Xu)) {
      stop("The number of predictor variables in Xr must be equal to the number of variables in Xu")
    }
  }

  if (ncol(Xr) < 4) {
    stop("This function works only with matrices containing more than 3 predictor variables")
  }

  if (length(Yr) != nrow(Xr)) {
    stop("length(Yr) must be equal to nrow(Xr)")
  }

  if (!is.null(Xu)) {
    if (any(is.na(Yr))) {
      stop("The current version of the mbl function does not handle NAs in the response variable of the reference observations (Yr)")
    }
  }


  Xr <- as.matrix(Xr)
  Yr <- as.matrix(Yr)

  rownames(Xr) <- 1:nrow(Xr)

  if (is.null(colnames(Xr))) {
    colnames(Xr) <- 1:ncol(Xr)
  }

  validation_type <- control$validation_type
  is_local_cv <- "local_cv" %in% validation_type
  is_nnv_val <- "NNv" %in% validation_type

  if (all(c("local_cv", "NNv") %in% control$validation_type)) {
    validation_type <- "both"
  }

  if (!is.null(Xu)) {
    pre_nms_ng <- "Xu_"
    n_xu <- ln <- nrow(Xu)

    Xu <- as.matrix(Xu)
    rownames(Xu) <- 1:nrow(Xu)

    if (is.null(colnames(Xu))) {
      colnames(Xu) <- 1:ncol(Xu)
    }

    if (sum(!colnames(Xu) == colnames(Xr)) != 0) {
      stop("Variable names in Xr do not match those in Xu")
    }

    if (validation_type %in% c("NNv", "both") & nrow(Xu) < 3) {
      stop(paste0(
        "For nearest neighbor validation (control$validation_type == 'NNv')",
        " Xu must contain at least 3 observations"
      ))
    }

    if (!is.null(Yu)) {
      Yu <- as.matrix(Yu)
      if (length(Yu) != nrow(Xu)) {
        stop("Number of observations in Yu and Xu differ")
      }
    }
    constellation <- FALSE
    first_nn <- 1
    observed <- Yu
    y_output_name <- "yu_obs"
    y_hat_output_name <- "pred"
    val_summary_name <- "Yu_prediction_statistics"
  } else {
    pre_nms_ng <- "Xr_"
    ln <- nrow(Xr)

    n_xu <- 0
    constellation <- TRUE
    first_nn <- 2
    observed <- Yr
    y_output_name <- "yr_obs"
    y_hat_output_name <- "fitted"
    val_summary_name <- "Yr_fitted_statistics"
  }



  n_xr <- nrow(Xr)
  n_total <- n_xr + n_xu


  diss_methods <- c(
    "pca", "pca.nipals", "pls", "cor",
    "euclid", "cosine", "sid"
  )

  if (!is.character(diss_method) & !is.matrix(diss_method)) {
    mtds <- paste(diss_methods, collapse = ", ")
    stop(paste0(
      "'diss_method' must be one of: ",
      mtds,
      " or a matrix"
    ))
  }

  if (!is.null(group)) {
    if (length(group) != nrow(Xr)) {
      stop(paste0(
        "The length of 'group' must be equal to the number of ",
        "observations in 'Xr'"
      ))
    }
  }

  if (length(pc_selection) != 2 | !is.list(pc_selection)) {
    stop("'pc_selection' must be a list of length 2")
  }

  if (!all(names(pc_selection) %in% c("method", "value")) | is.null(names(pc_selection))) {
    names(pc_selection)[sapply(pc_selection, FUN = is.character)] <- "method"
    names(pc_selection)[sapply(pc_selection, FUN = is.numeric)] <- "value"
  }

  pc_sel_method <- match.arg(pc_selection$method, c(
    "opc",
    "var",
    "cumvar",
    "manual"
  ))
  pc_threshold <- pc_selection$value

  if (pc_sel_method %in% c("opc", "manual") & pc_selection$value > min(n_total, ncol(Xr))) {
    pc_threshold <- min(n_total, ncol(Xr), n_total)

    if (!is.null(Xu)) {
      message_pc <- paste0(
        "When pc_selection$method is 'opc' or 'manual', the value ",
        "specified in \npc_selection$value cannot be larger than ",
        "min(nrow(Xr) + nrow(Xu), ncol(Xr)) \n(i.e ",
        pc_threshold,
        "). Therefore the value was reset to ",
        pc_threshold
      )
    } else {
      message_pc <- paste0(
        "When pc_selection$method is 'opc' or 'manual', the value ",
        "specified in \npc_selection$value cannot be larger than ",
        "min(dim(Xr)) \n(i.e ",
        pc_threshold,
        "). Therefore the value was reset to ",
        pc_threshold
      )
    }
    warning(message_pc)
  }


  match.arg(diss_usage, c("predictors", "weights", "none"))

  if (is.null(k) & is.null(k_diss)) {
    stop("Either k or k_diss must be specified")
  }

  k_max <- NULL
  if (!is.null(k)) {
    if (!is.null(k_diss)) {
      stop("Only one of k or k_diss can be specified")
    }
    if (!is.numeric(k)) {
      stop("k must be a vector of integers")
    } else {
      k <- unique(sort(ceiling(k)))
    }
    k <- sort(k)
    k_max <- max(k)
  }

  k_diss_max <- NULL
  if (!is.null(k_diss)) {
    k_diss <- unique(sort(k_diss))
    if (is.null(k_range)) {
      stop("If the k_diss argument is used, k_range must be specified")
    }
    if (length(k_range) != 2 | !is.numeric(k_range) | any(k_range < 1)) {
      stop("k_range must be a vector of length 2 which specifies the minimum (first value, larger than 0) and the maximum (second value) number of neighbors")
    }
    k_range <- sort(k_range)
    k_min_range <- as.integer(k_range[1])
    k_max_range <- as.integer(k_range[2])
    if (k_min_range < 4) {
      stop("Minimum number of nearest neighbors allowed is 4")
    }
    if (k_max_range > nrow(Xr)) {
      stop("Maximum number of nearest neighbors cannot exceed the number of reference observations")
    }
    k_diss_max <- max(k_diss)
  }

  if (".local" %in% names(input_dots)) {
    if (isTRUE(input_dots$local)) {
      if (!"pre_k" %in% names(input_dots)) {
        stop(paste0(
          "When .local = TRUE (passed to the ortho_diss method), the ",
          "'pre_k' argument must be specified"
        ))
      }
      if (input_dots$pre_k < k_max) {
        stop(paste0(
          "pre_k must be larger than ",
          ifelse(is.null(k), "max(k_range)", "max(k)")
        ))
      }
    }
  }

  if (!"local_fit" %in% class(method)) {
    stop("Object passed to method must be of class local_fit")
  }

  if (!is.null(k)) {
    k <- as.integer(k)
    if (min(k) < 4) {
      stop("Minimum number of nearest neighbors allowed is 3")
    }
    if (max(k) > nrow(Xr)) {
      stop(paste0(
        "The number of nearest neighbors cannot exceed the number ",
        "of observations in Xr"
      ))
    }
  }

  has_projection <- FALSE


  if (!is.matrix(diss_method)) {
    # when .local = TRUE, k_max is replaced with k_pre inside get_neighbor_info()

    if (is.null(Xu)) {
      rdiss <- FALSE
    } else {
      rdiss <- control$return_dissimilarity
    }

    spike <- c(spike, -which(is.na(Yr)))
    if (length(spike) == 0) {
      spike <- NULL
    }

    neighborhoods <- get_neighbor_info(
      Xr = Xr, Xu = Xu,
      diss_method = diss_method, Yr = Yr,
      k = k_max, k_diss = k_diss_max,
      k_range = k_range,
      spike = spike,
      pc_selection = pc_selection,
      return_dissimilarity = rdiss,
      center = center, scale = scale,
      gh = gh, diss_usage = diss_usage,
      allow_parallel = control$allow_parallel,
      ...
    )

    if (is.null(Xu)) {
      neighborhoods$neighbors <- rbind(
        k_0 = 1:ncol(neighborhoods$neighbors),
        neighborhoods$neighbors
      )
      neighborhoods$neighbors_diss <- rbind(k_0 = 0, neighborhoods$neighbors_diss)
      k <- k + 1 # each target observation does not count as a neighbor to itself
      diss_xr_xtarget <- neighborhoods$diss_xr_xr
    } else {
      diss_xr_xtarget <- neighborhoods$dissimilarity
    }

    if (!is.null(neighborhoods$projection)) {
      diss_xr_xtarget_projection <- neighborhoods$projection
      has_projection <- TRUE
    }
  } else {
    diss_xr_xr <- NULL

    dim_diss <- dim(diss_method)
    if (diss_usage == "predictors") {
      if (diff(dim_diss) != 0 | dim_diss[1] != n_total | any(diag(diss_method) != 0)) {
        stop(paste0(
          "If a matrix is passed to 'diss_method' ",
          "and diss_usage = 'predictors', this matrix must be ",
          "squared symmetric zeroes in its diagonal"
        ))
      }
      diss_xr_xr <- diss_method[1:nrow(Xr), 1:nrow(Xr)]

      if (!is.null(Xu)) {
        diss_method <- diss_method[1:nrow(Xr), (1 + nrow(Xr)):ncol(diss_method)]
      } else {
        diss_method <- diss_xr_xr
      }
    }

    if (!is.null(Xu)) {
      if (diss_usage %in% c("weights", "none")) {
        if (dim_diss[1] != n_xr & dim_diss[2] != n_xu) {
          stop(paste0(
            "If a matrix is passed to 'diss_method' ",
            "and 'diss_usage' argument is set to either 'weights' or  ",
            "'none', the number of rows and columns of this matrix ",
            "must be equal to the number of rows of 'Xr' and the ",
            "number of rows of 'Xu' respectively"
          ))
        }
      }
    }

    diss_xr_xtarget <- diss_method
    diss_method <- "external_matrix"

    neighborhoods <-
      diss_to_neighbors(
        diss_xr_xtarget,
        k = k_max, k_diss = k_diss_max,
        k_range = k_range,
        spike = spike,
        return_dissimilarity = control$return_dissimilarity,
        skip_first = ifelse(is.null(Xu), TRUE, FALSE)
      )

    if (is.null(Xu)) {
      neighborhoods$neighbors <- rbind(k_0 = 0, neighborhoods$neighbors)
      neighborhoods$neighbors_diss <- rbind(k_0 = 0, neighborhoods$neighbors_diss)
    }

    if (gh) {
      neighborhoods$gh$projection <- pls_projection(
        Xr = Xr, Xu = Xu,
        Yr = Yr,
        pc_selection = pc_selection,
        scale = scale, ...
      )
      neighborhoods$gh$gh_Xr <- f_diss(neighborhoods$gh$projection$scores,
        Xu = t(colMeans(neighborhoods$gh$projection$scores[1:nrow(Xr), ])),
        diss_method = "mahalanobis",
        center = FALSE, scale = FALSE
      )
      if (!is.null(Xu)) {
        neighborhoods$gh$gh_Xu <- neighborhoods$gh$gh_Xr[-c(1:nrow(Xr))]
      } else {
        neighborhoods$gh$gh_Xu <- NULL
      }

      neighborhoods$gh$gh_Xr <- neighborhoods$gh$gh_Xr[c(1:nrow(Xr))]
      neighborhoods$gh <- neighborhoods$gh[c("gh_Xr", "gh_Xu", "projection")]
    }

    neighborhoods$diss_xr_xr <- diss_xr_xr
    rm(diss_xr_xr)
    gc()
  }

  if (!is.null(k)) {
    smallest_neighborhood <- neighborhoods$neighbors[1:min(k), , drop = FALSE]
    smallest_n_neighbors <- colSums(!is.na(smallest_neighborhood))
  }

  if (!is.null(k_diss)) {
    min_diss <- neighborhoods$neighbors_diss <= min(k_diss)
    if (!is.null(spike)) {
      min_diss[1:length(spike), ] <- TRUE
    }
    smallest_neighborhood <- neighborhoods$neighbors
    smallest_neighborhood[!min_diss] <- NA
    smallest_n_neighbors <- colSums(!is.na(smallest_neighborhood))
    smallest_n_neighbors[smallest_n_neighbors < min(k_range)] <- min(k_range)
    smallest_n_neighbors[smallest_n_neighbors > max(k_range)] <- max(k_range)
  }


  if (is_local_cv) {
    min_n_samples <- floor(min(smallest_n_neighbors) * control$p) - 1
    min_cv_samples <- floor(min(k, k_range) * (1 - control$p))
    if (min_cv_samples < 3) {
      stop(paste0(
        "Local cross-validation requires at least 3 observations in ",
        "the hold-out set, the current cross-validation parameters ",
        "leave less than 3 observations in some neighborhoods."
      ))
    }
  } else {
    min_n_samples <- smallest_n_neighbors - 1
  }

  if (method$method %in% c("pls", "wapls")) {
    max_pls <- max(method$pls_c)
    if (any(min_n_samples < max_pls)) {
      stop(paste0(
        "More pls components than observations in some neighborhoods.\n",
        "If 'local_cv' is being used, consider that some ",
        "observations \nin the neighborhoods are hold-out for local ",
        "validation"
      ))
    }
  }


  if (!".local" %in% names(input_dots)) {
    iter_neighborhoods <- ith_mbl_neighbor(
      Xr = Xr, Xu = Xu, Yr = Yr, Yu = Yu,
      diss_usage = diss_usage,
      neighbor_indices = neighborhoods$neighbors,
      neighbor_diss = neighborhoods$neighbors_diss,
      diss_xr_xr = neighborhoods$diss_xr_xr,
      group = group
    )
  } else {
    iter_neighborhoods <- ith_mbl_neighbor(
      Xr = Xr, Xu = Xu, Yr = Yr, Yu = Yu,
      diss_usage = "none",
      neighbor_indices = neighborhoods$neighbors,
      neighbor_diss = neighborhoods$neighbors_diss,
      group = group
    )
  }

  r_fields <- c(
    "o_index", "k_diss", "k_original", "k", "npls", "min_pls", "max_pls",
    y_output_name, y_hat_output_name, "yr_min_obs", "yr_max_obs",
    "index_nearest_in_Xr", "index_farthest_in_Xr",
    "y_nearest", "y_nearest_pred",
    "y_farthest", "diss_nearest", "diss_farthest",
    "loc_rmse_cv", "loc_st_rmse_cv", "loc_n_components", "rep"
  )

  n_ith_result <- ifelse(is.null(k_diss), length(k), length(k_diss))

  template_pred_results <- data.table(
    matrix(
      NA, n_ith_result, length(r_fields),
      dimnames = list(NULL, r_fields)
    )
  )

  template_pred_results$rep[1] <- 0

  if (!is.null(k_diss)) {
    template_pred_results$k_diss <- k_diss
  } else {
    template_pred_results$k <- k
  }
  pg_bar_width <- 10
  # to_erase <- getOption("width") - pg_bar_width - (2 * nchar(nrow(Xu))) - 2
  if (!is.null(Xu)) {
    n_characters <- nchar(n_xu)
    n_iter <- n_xu
  } else {
    n_characters <- nchar(n_xr)
    n_iter <- n_xr
  }

  to_erase <- pg_bar_width + (2 * n_characters) + 8
  to_erase <- paste(rep(" ", to_erase), collapse = "")

  if (verbose) {
    cat("\033[32m\033[3mPredicting...\n\033[23m\033[39m")
  }
  
  pred_obs <- foreach(
    i = 1:n_iter,
    ith_observation = iter_neighborhoods,
    .inorder = FALSE,
    .export = c(
      "ortho_diss", "fit_and_predict", "pls_cv",
      "get_col_sds", "get_wapls_weights"
    ),
    .noexport = c("Xr", "Xu")
  ) %mydo% {
  ################
    # it <- iter_neighborhoods
    # pred_obs <- vector("list", n_iter)  # %do% returns a list; mimic that
    # for (i in seq_len(n_iter)) {
    #   ith_observation <- tryCatch(
    #     nextElem(it),
    #     error = function(e) {
    #       if (inherits(e, "StopIteration")) return(NULL)
    #       stop(e)
    #     }
    #   )
    #   if (is.null(ith_observation)) {
    #     pred_obs <- pred_obs[seq_len(i - 1L)]
    #     break
    #   }
    # #   ################
    
    
    ith_pred_results <- template_pred_results
    additional_results <- NULL
    ith_pred_results$o_index[] <- i

    if (".local" %in% names(input_dots) & diss_method %in% ortho_diss_methods) {
      ith_observation <- get_ith_local_neighbors(
        ith_xr = ith_observation$ith_xr,
        ith_xu = ith_observation$ith_xu,
        ith_yr = ith_observation$ith_yr,
        ith_yu = ith_observation$ith_yu,
        diss_usage = diss_usage,
        ith_neig_indices = ith_observation$ith_neig_indices,
        k = k_max, k_diss = k_diss_max,
        k_range = k_range,
        spike = spike,
        diss_method = diss_method,
        pc_selection = pc_selection,
        center = center, scale = scale,
        ith_group = ith_observation$ith_group,
        ...
      )

      ith_pred_results$loc_n_components[] <- ith_observation$ith_components
      additional_results$ith_neig_indices <- ith_observation$ith_neig_indices
      additional_results$ith_neigh_diss <- ith_observation$ith_neigh_diss
    }

    if (verbose) {
      cat(paste0("\033[34m\033[3m", i, "/", n_iter, "\033[23m\033[39m"))
      pb <- txtProgressBar(width = pg_bar_width, char = "\033[34m_\033[39m")
    }

    if (!is.null(k_diss)) {
      ith_diss <- ith_observation$ith_neigh_diss
      if (!is.null(spike)) {
        ith_diss[1:length(spike)] <- 0
      }
      ith_pred_results$k_original <- sapply(k_diss, FUN = function(x, d) sum(d < x), d = ith_diss)
      ith_pred_results$k <- ith_pred_results$k_original
      ith_pred_results$k[ith_pred_results$k_original < min(k_range)] <- min(k_range)
      ith_pred_results$k[ith_pred_results$k_original > max(k_range)] <- max(k_range)
    } else {
      ith_pred_results$k <- k
    }

    for (kk in 1:nrow(ith_pred_results)) {
      if (verbose) {
        setTxtProgressBar(pb, kk / nrow(ith_pred_results))
      }

      # If the sample has not been predicted before,
      # then create a model and predict it (useful only when k_diss is used)
      current_k <- ith_pred_results$k[kk]
      if (current_k != ifelse(kk == 1, 0, ith_pred_results$k[kk - 1])) {
        if (diss_usage == "predictors") {
          keep_cols <- c(
            1:current_k,
            (1 + ith_observation$n_k):ncol(ith_observation$ith_xr)
          )
          i_k_xr <- ith_observation$ith_xr[1:current_k, keep_cols]
          i_k_xu <- ith_observation$ith_xu[, keep_cols, drop = FALSE]
        } else {
          i_k_xr <- ith_observation$ith_xr[1:current_k, ]
          i_k_xu <- ith_observation$ith_xu
        }

        # for extracting some basic stats
        i_k_yr <- ith_observation$ith_yr[first_nn:current_k, , drop = FALSE]

        i_k_yu <- ith_observation$ith_yu
        kth_diss <- ith_observation$ith_neigh_diss[first_nn:current_k]
        i_idx <- ith_observation$ith_neig_indices[first_nn:current_k]


        ith_pred_results$rep[kk] <- 0
        ith_yr_range <- range(i_k_yr)
        ith_pred_results$yr_min_obs[kk] <- ith_yr_range[1]
        ith_pred_results$yr_max_obs[kk] <- ith_yr_range[2]
        ith_pred_results$diss_farthest[kk] <- max(kth_diss)
        ith_pred_results$diss_nearest[kk] <- min(kth_diss)
        ith_pred_results$y_farthest[kk] <- i_k_yr[which.max(kth_diss)]
        ith_pred_results$y_nearest[kk] <- i_k_yr[which.min(kth_diss)]
        ith_pred_results$index_nearest_in_Xr[kk] <- i_idx[which.min(kth_diss)]
        ith_pred_results$index_farthest_in_Xr[kk] <- i_idx[which.max(kth_diss)]

        # use the final y vector (in case Xu was not provided)
        i_k_yr <- ith_observation$ith_yr[1:current_k, , drop = FALSE]

        if (!is.null(group)) {
          i_k_group <- factor(ith_observation$ith_group[1:current_k])
        } else {
          i_k_group <- NULL
        }

        if (diss_usage == "weights") {
          if (is.null(Xu)) {
            stop("'weights' are not yet supported for diss_usage")
          }
          # Weights are defined according to a tricubic function
          # as in Cleveland and Devlin (1988) and Naes and Isaksson (1990).
          std_kth_diss <- kth_diss / max(kth_diss)
          kth_weights <- (1 - (std_kth_diss^3))^3
          kth_weights[which(kth_weights == 0)] <- 1e-04
        } else {
          kth_weights <- rep(1, current_k)
        }



        # FIXME: WHEN TO HAVE NNv AND HOW... HOW TO GET THE FINAL MODELS --> HOW
        # TO NOT CONFUSE A FINAL MODEL (WITHOUT VALIDATION) WITH A SINGLE MODEL
        # DOES REQUIRE VALIDATION?
        ith_observation$ith_neig_indices
        # local fit
        i_k_pred <- fit_and_predict(
          x = i_k_xr,
          y = i_k_yr,
          pred_method = method$method,
          scale = scale,
          pls_c = method$pls_c,
          weights = kth_weights,
          newdata = i_k_xu,
          CV = is_local_cv,
          tune = control$tune_locally,
          group = i_k_group,
          p = control$p,
          number = control$number,
          noise_variance = method$noise_variance,
          range_prediction_limits = control$range_prediction_limits,
          pls_max_iter = 1,
          pls_tol = 1e-6,
          seed = seed,
          modified = ifelse(is.null(method$modified), FALSE, method$modified) ## applies to pls only
        )

        # this first one will rutun the maximum number of components
        # of the one that was optimized

        i_k_pred$validation$models$coefficients

        # the second one is to rertrieve only NN validation stats

        # for the the first one, there must be one argument that indicates if
        # model is to be retrieved --> warnings when multiple k are run must be thrown.


        ith_pred_results[[y_hat_output_name]][kk] <- i_k_pred$prediction

        selected_pls <- NULL
        if (is_local_cv) {
          if (control$tune_locally) {
            best_row <- which.min(i_k_pred$validation$cv_results$rmse_cv)
          } else {
            best_row <- ifelse(method$method == "pls", method$pls_c, 1)
          }

          if (method$method == "pls") {
            ith_pred_results$npls[kk] <- i_k_pred$validation$cv_results$npls[best_row]
            selected_pls <- ith_pred_results$npls[kk]
          }
          if (method$method == "wapls") {
            ith_pred_results$min_pls[kk] <- i_k_pred$validation$cv_results$min_component[best_row]
            ith_pred_results$max_pls[kk] <- i_k_pred$validation$cv_results$max_component[best_row]
            selected_pls <- i_k_pred$validation$cv_results[best_row, 1:2]
          }

          ith_pred_results$loc_rmse_cv[kk] <- i_k_pred$validation$cv_results$rmse_cv[best_row]
          ith_pred_results$loc_st_rmse_cv[kk] <- i_k_pred$validation$cv_results$st_rmse_cv[best_row]
        } else {
          if (method$method == "pls") {
            ith_pred_results$npls[kk] <- method$pls_c
            selected_pls <- ith_pred_results$npls[kk]
          }
          if (method$method == "wapls") {
            ith_pred_results$min_pls[kk] <- method$pls_c[[1]]
            ith_pred_results$max_pls[kk] <- method$pls_c[[2]]
            selected_pls <- method$pls_c
          }
        }

        if (is_nnv_val) {
          if (!is.null(group)) {
            out_group <- which(i_k_group == i_k_group[[ith_observation$local_index_nearest]])
          } else {
            out_group <- ith_observation$local_index_nearest
          }

          nearest_pred <- fit_and_predict(
            x = i_k_xr[-out_group, ],
            y = i_k_yr[-out_group, , drop = FALSE],
            pred_method = method$method,
            scale = scale,
            pls_c = selected_pls,
            noise_variance = method$noise_variance,
            newdata = i_k_xr[ith_observation$local_index_nearest, , drop = FALSE],
            CV = FALSE,
            tune = FALSE,
            range_prediction_limits = control$range_prediction_limits,
            pls_max_iter = 1,
            pls_tol = 1e-6,
            seed = seed,
            modified = ifelse(is.null(method$modified), FALSE, method$modified) ## applies to pls only
          )$prediction

          ith_pred_results$y_nearest_pred[kk] <- nearest_pred / kth_weights[1]
        }
      } else {
        ith_k_diss <- ith_pred_results$k_diss[kk]
        ith_pred_results[kk, ] <- ith_pred_results[kk - 1, ]
        ith_pred_results$rep[kk] <- 1
        ith_pred_results$k_diss[kk] <- ith_k_diss
      }
    }

    if (verbose) {
      if (kk == nrow(ith_pred_results) & i != n_iter) {
        cat("\r", to_erase, "\r")
      }

      if (i == n_iter) {
        cat("\n")
      }
      # do not use close() (it prints a new line)
      ## close(pb)
    }
    list(
      results = ith_pred_results,
      additional_results = additional_results
    )
  }

  iteration_order <- sapply(
    pred_obs,
    FUN = function(x) x$results$o_index[1]
  )

  pred_obs <- pred_obs[order(iteration_order, decreasing = FALSE)]

  results_table <- do.call(
    "rbind", lapply(pred_obs, FUN = function(x) x$results)
  )

  if (is.null(Xu) & !is.null(k)) {
    results_table$k <- results_table$k - 1
    fix_k <- 1
  } else {
    fix_k <- 0
  }


  if (".local" %in% names(input_dots) & diss_method %in% ortho_diss_methods) {
    diss_xr_xtarget <- do.call(
      "cbind",
      lapply(iteration_order,
        FUN = function(x, m, ii) {
          idc <- x[[ii]]$additional_results$ith_neig_indices
          d <- x[[ii]]$additional_results$ith_neigh_diss
          m[idc] <- d
          m
        },
        x = pred_obs,
        m = matrix(NA, nrow(Xr), 1)
      )
    )
    class(diss_xr_xtarget) <- c("local_ortho_diss", "matrix")

    dimnames(diss_xr_xtarget) <- list(
      paste0("Xr_", 1:nrow(diss_xr_xtarget)),
      paste0(pre_nms_ng, 1:ncol(diss_xr_xtarget))
    )

    neighborhoods$neighbors <- do.call(
      "cbind", lapply(iteration_order,
        FUN = function(x, m, ii) {
          idc <- x[[ii]]$additional_results$ith_neig_indices
          m[1:length(idc)] <- idc
          m
        },
        x = pred_obs,
        m = matrix(NA, max(results_table$k), 1)
      )
    )
  }


  out <- c(
    if (is.null(Yu) & !is.null(Xu)) {
      "yu_obs"
    },
    if (all(is.na(results_table$k_original))) {
      "k_original"
    },
    if (!(validation_type %in% c("NNv", "both"))) {
      "y_nearest_pred"
    },
    if (method$method != "wapls") {
      c("min_pls", "max_pls")
    },
    if (method$method != "pls") {
      "npls"
    },
    if (!(validation_type %in% c("local_cv", "both"))) {
      c("loc_rmse_cv", "loc_st_rmse_cv")
    },
    "rep"
  )

  results_table[, (out) := NULL]
  if (is.null(Xu)) {
    names(results_table)[names(results_table) %in% "yu_obs"] <- "yr_obs"
    names(results_table)[names(results_table) %in% "pred"] <- "fitted"
  }


  if (!is.null(k_diss)) {
    param <- "k_diss"
    results_table <- lapply(
      get(param),
      FUN = function(x, sel, i) x[x[[sel]] == i, ],
      x = results_table,
      sel = param
    )
    names(results_table) <- paste0("k_diss_", k_diss)
    p_bounded <- sapply(
      results_table,
      FUN = function(x, k_range) {
        sum(x$k_original <= k_range[1] | x$k_original >= k_range[2])
      },
      k_range = k_range
    )
    col_ks <- data.table(
      k_diss = k_diss,
      p_bounded = paste0(round(100 * p_bounded / nrow(Xu), 3), "%")
    )
  } else {
    param <- "k"
    results_table <- lapply(
      get(param) - fix_k,
      FUN = function(x, sel, i) x[x[[sel]] == i, ],
      x = results_table,
      sel = param
    )
    names(results_table) <- paste0("k_", k - fix_k)
    col_ks <- data.table(k = k - fix_k)
  }

  if (validation_type %in% c("NNv", "both")) {
    nn_stats <- function(x) {
      nn_rmse <- (mean((x$y_nearest - x$y_nearest_pred)^2))^0.5
      nn_st_rmse <- nn_rmse / diff(range(x$y_nearest))
      nn_rsq <- (cor(x$y_nearest, x$y_nearest_pred))^2
      c(nn_rmse = nn_rmse, nn_st_rmse = nn_st_rmse, nn_rsq = nn_rsq)
    }

    loc_nn_res <- do.call("rbind", lapply(results_table, FUN = nn_stats))
    loc_nn_res <- cbind(
      col_ks,
      rmse = loc_nn_res[, "nn_rmse"],
      st_rmse = loc_nn_res[, "nn_st_rmse"],
      r2 = loc_nn_res[, "nn_rsq"]
    )
  } else {
    loc_nn_res <- NULL
  }

  if (validation_type %in% c("local_cv", "both")) {
    mean_loc_res <- function(x) {
      mean_loc_rmse <- mean(x$loc_rmse_cv)
      mean_loc_st_rmse <- mean(x$loc_st_rmse_cv)
      c(loc_rmse = mean_loc_rmse, loc_st_rmse = mean_loc_st_rmse)
    }
    loc_res <- do.call("rbind", lapply(results_table, mean_loc_res))
    loc_res <- cbind(
      col_ks,
      rmse = loc_res[, "loc_rmse"],
      st_rmse = loc_res[, "loc_st_rmse"]
    )
  } else {
    loc_res <- NULL
  }

  if (!is.null(observed)) {
    for (i in 1:length(results_table)) {
      results_table[[i]][[y_output_name]] <- observed
    }
    yu_stats <- function(x, y_hat, y) {
      y_rmse <- mean((x[[y_hat]] - x[[y]])^2, na.rm = TRUE)^0.5
      y_st_rmse <- y_rmse / diff(range(x[[y_hat]]), na.rm = TRUE)
      y_rsq <- cor(x[[y_hat]], x[[y]], use = "complete.obs")^2
      c(rmse = y_rmse, st_rmse = y_st_rmse, rsq = y_rsq)
    }

    pred_res <- do.call(
      "rbind",
      lapply(results_table, yu_stats, y_hat = y_hat_output_name, y = y_output_name)
    )

    pred_res <- cbind(
      col_ks,
      rmse = pred_res[, "rmse"],
      st_rmse = pred_res[, "st_rmse"],
      r2 = pred_res[, "rsq"]
    )
  } else {
    pred_res <- NULL
  }

  if ("local_ortho_diss" %in% class(diss_xr_xtarget)) {
    diss_method <- paste0(diss_method, " (locally computed)")
  }

  if (control$return_dissimilarity) {
    diss_list <- list(
      diss_method = diss_method,
      diss_xr_xu = diss_xr_xtarget
    )
    if (has_projection) {
      diss_list$global_projection <- diss_xr_xtarget_projection
    }
  } else {
    diss_list <- NULL
  }

  colnames(neighborhoods$neighbors) <- paste0(pre_nms_ng, 1:ln)
  rownames(neighborhoods$neighbors) <- paste0("k_", 1:nrow(neighborhoods$neighbors))


  val_list <- structure(
    list(
      loc_res,
      loc_nn_res,
      pred_res
    ),
    names = c(
      "local_cross_validation",
      "nearest_neighbor_validation",
      val_summary_name
    )
  )

  results_list <- list(
    call = f_call,
    cntrl_param = control,
    dissimilarities = diss_list,
    Xu_neighbors = list(
      neighbors = neighborhoods$neighbors,
      neighbors_diss = neighborhoods$neighbors_diss
    ),
    n_predictions = nrow(Xu),
    gh = neighborhoods$gh,
    validation_results = val_list,
    results = results_table,
    documentation = documentation,
    seed = seed
  )

  if (is.null(Xu)) {
    names(results_list)[names(results_list) %in% "Xu_neighbors"] <- "Xr_neighbors"
  }

  attr(results_list, "call") <- f_call
  class(results_list) <- c("mbl", "list")

  results_list
}
