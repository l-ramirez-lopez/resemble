% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mbl.R
\name{mbl}
\alias{mbl}
\title{A function for memory-based learning (mbl)}
\usage{
mbl(Xr, Yr, Xu, Yu = NULL, k, k_diss, k_range, spike = NULL,
    method = local_fit_wapls(min_pls_c = 3, max_pls_c = min(dim(Xr), 15)),
    diss_method = "pca", diss_usage = "predictors", gh = TRUE,
    pc_selection = list(method = "opc", value = min(dim(Xr), 40)),
    control = mbl_control(), group = NULL, center = TRUE, scale = FALSE,
    verbose = TRUE, documentation = character(), seed = NULL, ...)
}
\arguments{
\item{Xr}{a matrix of predictor variables of the reference data
(observations in rows and variables in columns).}

\item{Yr}{a numeric matrix of one column containing the values of the
response variable corresponding to the reference data.}

\item{Xu}{a matrix of predictor variables of the data to be predicted
(observations in rows and variables in columns).}

\item{Yu}{an optional matrix of one column containing the values of the
response variable corresponding to the data to be predicted. Default is
\code{NULL}.}

\item{k}{a vector of integers specifying the sequence of k-nearest
neighbors to be tested. Either \code{k} or \code{k_diss} must be specified.
This vector will be automatically sorted into ascending order. If
non-integer numbers are passed, they will be coerced to the next upper
integers.}

\item{k_diss}{a numeric vector specifying the sequence of dissimilarity
thresholds to be tested for the selection of the nearest neighbors found in
\code{Xr} around each observation in \code{Xu}. These thresholds depend on
the corresponding dissimilarity measure specified in the object passed to
\code{control}. Either \code{k} or \code{k_diss} must be specified.}

\item{k_range}{an integer vector of length 2 which specifies the minimum
(first value) and the maximum (second value) number of neighbors to be
retained when the \code{k_diss} is given.}

\item{spike}{an integer vector (with positive and/or negative values) indicating
the indices of observations in \code{Xr} that must be either be forced into
or avoided in the neighborhoods of every \code{Xu} observation. Default is
\code{NULL} (i.e. no observations are forced or avoided). Note
that this argument is not intended for increasing or reducing the neighborhood
size which is only controlled by \code{k} or \code{k_diss} and \code{k_range}.
By forcing observations into the neighborhood, some of the farthest
observations may be forced out of the neighborhood. In contrast, by avoiding
observations in the neighborhood,  some of farthest
observations may be included into the neighborhood. See details.}

\item{method}{an object of class \code{\link{local_fit}} which indicates the
type of regression to conduct at each local segment as well as additional
parameters affecting this regression. See \code{\link{local_fit}} function.}

\item{diss_method}{a character string indicating the spectral dissimilarity
metric to be used in the selection of the nearest neighbors of each
observation. Options are:
\itemize{
\item{\code{"pca"} (Default):}{ Mahalanobis distance
computed on the matrix of scores of a Principal Component (PC)
projection of \code{Xr} and \code{Xu}. PC projection is done using the
singular value decomposition (SVD) algorithm.
See \code{\link{ortho_diss}} function.}
\item{\code{"pca.nipals"}}{ Mahalanobis distance
computed on the matrix of scores of a Principal Component (PC)
projection of \code{Xr} and \code{Xu}. PC projection is done using the
non-linear iterative partial least squares (nipals) algorithm.
See \code{\link{ortho_diss}} function.}
\item{\code{"pls"}}{ Mahalanobis distance
computed on the matrix of scores of a partial least squares projection
of \code{Xr} and \code{Xu}. In this case, \code{Yr} is always
required. See \code{\link{ortho_diss}} function.}
\item{\code{"cor"}}{ correlation coefficient
between observations. See \code{\link{cor_diss}} function.}
\item{\code{"euclid"}}{ Euclidean distance
between observations. See \code{\link{f_diss}} function.}
\item{\code{"cosine"}}{ Cosine distance
between observations. See \code{\link{f_diss}} function.}
\item{\code{"sid"}}{ spectral information divergence between
observations. See \code{\link{sid}} function.}
}
Alternatively, a matrix of dissimilarities can also be passed to this
argument. This matrix is supposed to be a user-defined matrix
representing the dissimilarities between observations in \code{Xr} and
\code{Xu}. When \code{diss_usage = "predictors"}, this matrix must be squared
(derived from a matrix of the form \code{rbind(Xr, Xu)}) for which the
diagonal values are zeros (since the dissimilarity between an object and
itself must be 0). On the other hand, if \code{diss_usage} is set to either
\code{"weights"} or \code{"none"}, it must be a matrix representing the
dissimilarity of each observation in \code{Xu} to each observation in
\code{Xr}. The number of columns of the input matrix must be equal to the
number of rows in \code{Xu} and the number of rows equal to the number of
rows in \code{Xr}.}

\item{diss_usage}{a character string specifying how the dissimilarity
information shall be used. The possible options are: \code{"predictors"},
\code{"weights"} and \code{"none"} (see details below).
Default is \code{"predictors"}.}

\item{gh}{a logical indicating if the global Mahalanobis distance (in the pls
score space) between each observation and the pls mean (centre) must be
computed. This metric is known as the GH distance in the literature. Note
that this computation is based on the number of pls components determined by
using the \code{pc_selection} argument. See details.}

\item{pc_selection}{a list of length 2 used for the computation of GH (if
\code{gh = TRUE}) as well as in the computation of the dissimilarity methods
based on \code{\link{ortho_diss}} (i.e. when \code{diss_method} is one of:
\code{"pca"}, \code{"pca.nipals"} or \code{"pls"}) or when \code{gh = TRUE}.
This argument is used for optimizing the number of components (principal
components or pls factors) to be retained for dissimilarity/distance
computation purposes only (i.e not for regression).
This list must contain two elements in the following order:
\code{method} (a character indicating the method for selecting the number of
components) and \code{value} (a numerical value that complements the selected
method). The methods available are:
\itemize{
\item{\code{"opc"}:} { optimized principal component selection based
on Ramirez-Lopez et al. (2013a, 2013b). The optimal number of
components (of set of observations) is the one for which its distance
matrix minimizes the differences between the \code{Yr} value of each
observation and the \code{Yr} value of its closest observation. In
this case \code{value} must be a value (larger than 0 and
below the minimum dimension of \code{Xr} or \code{Xr} and \code{Xu}
combined) indicating the maximum
number of principal components to be tested. See the
\code{\link{ortho_projection}} function for more details.}
\item{\code{"cumvar"}:}{ selection of the principal components based
on a given cumulative amount of explained variance. In this case,
\code{value} must be a value (larger than 0 and below or equal to 1)
indicating the minimum amount of cumulative variance that the
combination of retained components should explain.}
\item{\code{"var"}:}{ selection of the principal components based
on a given amount of explained variance. In this case,
\code{value} must be a value (larger than 0 and below or equal to 1)
indicating the minimum amount of variance that a single component
should explain in order to be retained.}
\item{\code{"manual"}:}{ for manually specifying a fix number of
principal components. In this case, \code{value} must be a value
(larger than 0 and below the minimum dimension of \code{Xr} or
\code{Xr} and \code{Xu} combined).
indicating the minimum amount of variance that a component should
explain in order to be retained.}
}
The list
\code{list(method = "opc", value = min(dim(Xr), 40))} is the default.
Optionally, the \code{pc_selection} argument admits \code{"opc"} or
\code{"cumvar"} or \code{"var"} or \code{"manual"} as a single character
string. In such a case the default \code{"value"} when either \code{"opc"} or
\code{"manual"} are used is 40. When \code{"cumvar"} is used the default
\code{"value"} is set to 0.99 and when \code{"var"} is used, the default
\code{"value"} is set to 0.01.}

\item{control}{a list created with the \code{\link{mbl_control}} function
which contains additional parameters that control some few aspects of the
\code{mbl} function (cross-validation, parameter tuning, etc).
The default list is as returned by \code{mbl_control()}.
See the \code{\link{mbl_control}} function for more details.}

\item{group}{an optional factor (or character vector vector
that can be coerced to \code{\link[base]{factor}} by \code{as.factor}) that
assigns a group/class label to each observation in \code{Xr}
(e.g. groups can be given by spectra collected from the same batch of
measurements, from the same observation, from observations with very similar
origin, etc). This is taken into account for internal leave-group-out cross
validation for pls tuning (factor optimization) to avoid pseudo-replication.
When one observation is selected for cross-validation, all observations of
the same group are removed together and assigned to validation. The length
of the vector must be equal to the number of observations in the
reference/training set (i.e. \code{nrow(Xr)}). See details.}

\item{center}{a logical if the predictor variables must be centred at each
local segment (before regression). In addition, if \code{TRUE}, \code{Xr}
and \code{Xu} will be centred for  dissimilarity computations.}

\item{scale}{a logical indicating if the predictor variables must be scaled
to unit variance at each local segment (before regression). In addition, if
\code{TRUE}, \code{Xr} and \code{Xu} will be scaled for  dissimilarity
computations.}

\item{verbose}{a logical indicating whether or not to print a progress bar
for each observation to be predicted. Default is \code{TRUE}. Note: In case
parallel processing is used, these progress bars will not be printed.}

\item{documentation}{an optional character string that can be used to
describe anything related to the \code{mbl} call (e.g. description of the
input data). Default: \code{character()}. NOTE: his is an experimental
argument.}

\item{seed}{an integer value containing the random number generator (RNG)
state for random number generation. This argument can be used for
reproducibility purposes (for random sampling) in the cross-validation
results. Default is \code{NULL}, i.e. no RNG is applied.}

\item{...}{further arguments to be passed to the \code{\link{dissimilarity}}
function. See details.}
}
\value{
a \code{list} of class \code{mbl} with the following components
(sorted either by \code{k} or \code{k_diss}):

\itemize{
\item{\code{call}:}{ the call to mbl.}
\item{\code{cntrl_param}:}{ the list with the control parameters passed to
control.}
\item{\code{Xu_neighbors}:}{ a list containing two elements: a matrix of
\code{Xr} indices corresponding to the neighbors of \code{Xu} and a matrix
of dissimilarities between each \code{Xu} observation and its corresponding
neighbor in \code{Xr}.}
\item{\code{dissimilarities}:}{ a list with the method used to obtain the
dissimilarity matrices and the dissimilarity matrix corresponding to
\mjeqn{D(Xr, Xu)}{D(Xr, Xu)}. This object is returned only if the
\code{return_dissimilarity} argument in the \code{control} list was set
to \code{TRUE}.}
\item{\code{n_predictions}}{ the total number of observations predicted.}
\item{\code{gh}:}{ if \code{gh = TRUE}, a list containing the global
Mahalanobis distance values for the observations in \code{Xr} and \code{Xu}
as well as the results of the global pls projection object used to obtain
the GH values.}
\item{\code{validation_results}:}{ a list of validation results for
"local cross validation" (returned if the \code{validation_type} in
\code{control} list was set to \code{"local_cv"}),
"nearest neighbor validation" (returned if the \code{validation_type}
in \code{control} list was set to \code{"NNv"}) and
"Yu prediction statistics" (returned  if \code{Yu} was supplied).}``
\item{\code{results}:}{ a list of data tables containing the results of the
predictions for each either \code{k} or \code{k_diss}. Each data table
contains the following columns:}
\itemize{
\item{\code{o_index}:}{ The index of the predicted observation.}
\item{\code{k_diss}:}{ This column is only output if the \code{k_diss}
argument is used. It indicates the corresponding dissimilarity threshold
for selecting the neighbors.}
\item{\code{k_original}:}{ This column is only output if the \code{k_diss}
argument is used. It indicates the number of neighbors that were originally
found when the given dissimilarity threshold is used.}
\item{\code{k}:}{ This column indicates the final number of neighbors
used.}
\item{\code{npls}:}{ This column is only output if the \code{pls}
regression method was used. It indicates the final number of pls
components used.}
\item{\code{min_pls}:}{ This column is only output if \code{wapls}
regression method was used. It indicates the final number of minimum pls
components used. If no optimization was set, it retrieves the original
minimum pls components passed to the \code{method} argument.}
\item{\code{max_pls}:}{ This column is only output if the \code{wapls}
regression method was used. It indicates the final number of maximum pls
components used. If no optimization was set, it retrieves the original
maximum pls components passed to the \code{method} argument.}
\item{\code{yu_obs}:}{ The input values given in \code{Yu} (the response
variable corresponding to the data to be predicted). If \code{Yu = NULL},
then \code{NA}s are retrieved.}
\item{\code{pred}:}{ The predicted Yu values.}
\item{\code{yr_min_obs}:}{ The minimum reference value (of the response
variable) in the neighborhood.}
\item{\code{yr_max_obs}:}{ The maximum reference value (of the response
variable) in the neighborhood.}
\item{\code{index_nearest_in_Xr}}{ The index of the nearest neighbor found
in \code{Xr}.}
\item{\code{index_farthest_in_Xr}}{ The index of the farthest neighbor
found in \code{Xr}.}
\item{\code{y_nearest}:}{ The reference value (\code{Yr}) corresponding to
the nearest neighbor found in \code{Xr}.}
\item{\code{y_nearest_pred}:}{ This column is only output if the
validation method in the object passed to \code{control} was set to
\code{"NNv"}. It represents the predicted value of the nearest neighbor
observation found in \code{Xr}. This prediction come from model fitted
with the remaining observations in the neighborhood of the target
observation in \code{Xu}.}
\item{\code{loc_rmse_cv}:}{ This column is only output if the validation
method in the object passed to \code{control} was set to
\code{'local_cv'}. It represents the RMSE of the cross-validation
computed for the neighborhood of the target observation in \code{Xu}.}
\item{\code{loc_st_rmse_cv}:}{ This column is only output if the
validation method in the object passed to \code{control} was set to
\code{'local_cv'}. It represents the standardized RMSE of the
cross-validation computed for the neighborhood of the target observation
in \code{Xu}.}
\item{\code{dist_nearest}:}{ The distance to the nearest neighbor.}
\item{\code{dist_farthest}:}{ The distance to the farthest neighbor.}
\item{\code{loc_n_components}:}{ This column is only output if the
dissimilarity method used is one of \code{"pca"}, \code{"pca.nipals"} or
\code{"pls"} and in addition the dissimilarities are requested to be
computed locally by passing \code{.local = TRUE} to the \code{mbl}
function.
See \code{.local} argument in the \code{\link{ortho_diss}} function.}
}
\item{\code{seed}}{ a value mirroing the one passed to seed.}
\item{\code{documentation}:}{ a character string mirroring the one provided
in the \code{documentation} argument.}
}
When the \code{k_diss} argument is used, the printed results show a table
with a column named '\code{p_bounded}. It represents the percentage of
observations for which the neighbors selected by the given dissimilarity
threshold were outside the boundaries specified in the \code{k_range}
argument.
}
\description{
\loadmathjax
This function is implemented for memory-based learning (a.k.a.
instance-based learning or local regression) which is a non-linear lazy
learning approach for predicting a given response variable from a set of
predictor variables. For each observation in a prediction set, a specific
local regression is carried out based on a subset of similar observations
(nearest neighbors) selected from a reference set. The local model is
then used to predict the response value of the target (prediction)
observation. Therefore this function does not yield a global
regression model.
}
\details{
The argument \code{spike} can be used to indicate what reference observations
in \code{Xr} must be kept in the neighborhood of every single \code{Xu}
observation. If a vector of length \mjeqn{m}{m} is passed to this argument,
this means that the \mjeqn{m}{m} original neighbors with the largest
dissimilarities to the target observations will be forced out of the
neighborhood. Spiking might be useful in cases where
some reference observations are known to be somehow related to the ones in
\code{Xu} and therefore might be relevant for fitting the local models. See
Guerrero et al. (2010) for an example on the benefits of spiking.

The \code{mbl} function uses the \code{\link{dissimilarity}} function to
compute the dissimilarities between \code{Xr} and \code{Xu}. The dissimilarity
method to be used is specified in the \code{diss_method} argument.
Arguments to \code{\link{dissimilarity}} as well as further arguments to the
functions used inside \code{\link{dissimilarity}}
(i.e. \code{\link{ortho_diss}} \code{\link{cor_diss}} \code{\link{f_diss}}
\code{\link{sid}}) can be passed to those functions by using \code{...}.

The \code{diss_usage} argument is used to specify whether the dissimilarity
information must be used within the local regressions and, if so, how.
When \code{diss_usage = "predictors"} the local (square symmetric)
dissimilarity matrix corresponding the selected neighborhood is used as
source of additional predictors (i.e the columns of this local matrix are
treated as predictor variables). In some cases this results in an improvement
of the prediction performance (Ramirez-Lopez et al., 2013a).
If \code{diss_usage = "weights"}, the neighbors of the query point
(\mjeqn{xu_{j}}{xu_j}) are weighted according to their dissimilarity to
\mjeqn{xu_{j}}{xu_j} before carrying out each local regression. The following
tricubic function (Cleveland and Delvin, 1988; Naes et al., 1990) is used for
computing the final weights based on the measured dissimilarities:

\mjdeqn{W_{j}  =  (1 - v^{3})^{3}}{W_j  =  (1 - v^3)^3}

where if \mjeqn{{xr_{i} \in }}{xr_i in} neighbors of \mjeqn{xu_{j}}{xu_j}:

\mjdeqn{v_{j}(xu_{j})  =  d(xr_{i}, xu_{j})}{v_j(xu_j)  =  d(xr_i, xu_j)}

otherwise:

\mjdeqn{v_{j}(xu_{j})  =  0}{v_j(xu_j)  =  0}

In the above formulas \mjeqn{d(xr_{i}, xu_{j})}{d(xr_i, xu_j)} represents the
dissimilarity between the query point and each object in \mjeqn{Xr}{Xr}.
When \code{diss_usage = "none"} is chosen the dissimilarity information is
not used.

The global Mahalanobis distance (a.k.a GH) is computed based on the scores
of a pls projection. A pls projection model is built with for \code{\{Yr\}, \{Xr\}}
and this model is used to obtain the pls scores of the \code{Xu}
observations. The Mahalanobis distance between each \code{Xu} observation in
(the pls space) and the centre of \code{Xr} is then computed. The number of
pls components is optimized based on the parameters passed to the
\code{pc_selection} argument. In addition, the \code{mbl} function also
reports the GH distance for the observations in \code{Xr}.

Some aspects of the mbl process, such as the type of internal validation,
parameter tuning, what extra objects to return, permission for parallel
execution, prediction limits, etc, can be specified by using the
\code{\link{mbl_control}} function.

By using the \code{group} argument one can specify groups of observations
that have something in common (e.g. observations with very similar origin).
The purpose of \code{group} is to avoid biased cross-validation results due
to pseudo-replication. This argument allows to select calibration points
that are independent from the validation ones. In this regard, when
\code{validation_type = "local_cv"} (used in \code{\link{mbl_control}}
function), then the \code{p} argument refers to the percentage of groups of
observations (rather than single observations) to be retained in each
sampling iteration at each local segment.
}
\examples{
\donttest{
library(prospectr)
data(NIRsoil)

# Proprocess the data using detrend plus first derivative with Savitzky and
# Golay smoothing filter
sg_det <- savitzkyGolay(
  detrend(NIRsoil$spc,
    wav = as.numeric(colnames(NIRsoil$spc))
  ),
  m = 1,
  p = 1,
  w = 7
)

NIRsoil$spc_pr <- sg_det

# split into training and testing sets
test_x <- NIRsoil$spc_pr[NIRsoil$train == 0 & !is.na(NIRsoil$CEC), ]
test_y <- NIRsoil$CEC[NIRsoil$train == 0 & !is.na(NIRsoil$CEC)]

train_y <- NIRsoil$CEC[NIRsoil$train == 1 & !is.na(NIRsoil$CEC)]
train_x <- NIRsoil$spc_pr[NIRsoil$train == 1 & !is.na(NIRsoil$CEC), ]

# Example 1
# A mbl implemented in Ramirez-Lopez et al. (2013,
# the spectrum-based learner)
# Example 1.1
# An exmaple where Yu is supposed to be unknown, but the Xu
# (spectral variables) are known
my_control <- mbl_control(validation_type = "NNv")

## The neighborhood sizes to test
ks <- seq(40, 140, by = 20)

sbl <- mbl(
  Xr = train_x,
  Yr = train_y,
  Xu = test_x,
  k = ks,
  method = local_fit_gpr(),
  control = my_control,
  scale = TRUE
)
sbl
plot(sbl)
get_predictions(sbl)

# Example 1.2
# If Yu is actually known...
sbl_2 <- mbl(
  Xr = train_x,
  Yr = train_y,
  Xu = test_x,
  Yu = test_y,
  k = ks,
  method = local_fit_gpr(),
  control = my_control
)
sbl_2
plot(sbl_2)

# Example 2
# the LOCAL algorithm (Shenk et al., 1997)
local_algorithm <- mbl(
  Xr = train_x,
  Yr = train_y,
  Xu = test_x,
  Yu = test_y,
  k = ks,
  method = local_fit_wapls(min_pls_c = 3, max_pls_c = 15),
  diss_method = "cor",
  diss_usage = "none",
  control = my_control
)
local_algorithm
plot(local_algorithm)

# Example 3
# A variation of the LOCAL algorithm (using the optimized pc
# dissmilarity matrix) and dissimilarity matrix as source of
# additional preditors
local_algorithm_2 <- mbl(
  Xr = train_x,
  Yr = train_y,
  Xu = test_x,
  Yu = test_y,
  k = ks,
  method = local_fit_wapls(min_pls_c = 3, max_pls_c = 15),
  diss_method = "pca",
  diss_usage = "predictors",
  control = my_control
)
local_algorithm_2
plot(local_algorithm_2)

# Example 4
# Running the mbl function in parallel with example 2

n_cores <- 2

if (parallel::detectCores() < 2) {
  n_cores <- 1
}

# Alternatively:
# n_cores <- parallel::detectCores() - 1
# if (n_cores == 0) {
#  n_cores <- 1
# }

library(doParallel)
clust <- makeCluster(n_cores)
registerDoParallel(clust)

# Alernatively:
# library(doSNOW)
# clust <- makeCluster(n_cores, type = "SOCK")
# registerDoSNOW(clust)
# getDoParWorkers()

local_algorithm_par <- mbl(
  Xr = train_x,
  Yr = train_y,
  Xu = test_x,
  Yu = test_y,
  k = ks,
  method = local_fit_wapls(min_pls_c = 3, max_pls_c = 15),
  diss_method = "cor",
  diss_usage = "none",
  control = my_control
)
local_algorithm_par

registerDoSEQ()
try(stopCluster(clust))

# Example 5
# Using local pls distances
with_local_diss <- mbl(
  Xr = train_x,
  Yr = train_y,
  Xu = test_x,
  Yu = test_y,
  k = ks,
  method = local_fit_wapls(min_pls_c = 3, max_pls_c = 15),
  diss_method = "pls",
  diss_usage = "predictors",
  control = my_control,
  .local = TRUE,
  pre_k = 150,
)
with_local_diss
plot(with_local_diss)
}
}
\references{
Cleveland, W. S., and Devlin, S. J. 1988. Locally weighted regression: an
approach to regression analysis by local fitting. Journal of the American
Statistical Association, 83, 596-610.

Guerrero, C., Zornoza, R., Gómez, I., Mataix-Beneyto, J. 2010. Spiking of
NIR regional models using observations from target sites: Effect of model
size on prediction accuracy. Geoderma, 158(1-2), 66-77.

Naes, T., Isaksson, T., Kowalski, B. 1990. Locally weighted regression and
scatter correction for near-infrared reflectance data. Analytical Chemistry
62, 664-673.

Ramirez-Lopez, L., Behrens, T., Schmidt, K., Stevens, A., Dematte, J.A.M.,
Scholten, T. 2013a. The spectrum-based learner: A new local approach for
modeling soil vis-NIR spectra of complex data sets. Geoderma 195-196,
268-279.

Ramirez-Lopez, L., Behrens, T., Schmidt, K., Viscarra Rossel, R., Dematte,
J. A. M.,  Scholten, T. 2013b. Distance and similarity-search metrics for
use with soil vis-NIR spectra. Geoderma 199, 43-53.

Rasmussen, C.E., Williams, C.K. Gaussian Processes for Machine Learning.
Massachusetts Institute of Technology: MIT-Press, 2006.

Shenk, J., Westerhaus, M., and Berzaghi, P. 1997. Investigation of a LOCAL
calibration procedure for near infrared instruments. Journal of Near
Infrared Spectroscopy, 5, 223-232.
}
\seealso{
\code{\link{mbl_control}}, \code{\link{f_diss}},
\code{\link{cor_diss}}, \code{\link{sid}}, \code{\link{ortho_diss}},
\code{\link{search_neighbors}},  \code{\link{local_fit}}
}
\author{
\href{https://orcid.org/0000-0002-5369-5120}{Leonardo Ramirez-Lopez}
and Antoine Stevens
}
