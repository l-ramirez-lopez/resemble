\name{mbl}
\alias{mbl}
\title{A function for memory-based learning (mbl)}
\usage{
mbl(Yr, Xr, Yu = NULL, Xu,
    mblCtrl = mblControl(),
    dissimilarityM = NULL,
    dissUsage = "predictors",
    k, k.diss, k.range,
    method,
    pls.c,
    noise.v = 0.001,
    ...)
}
\arguments{
  \item{Xr}{input \code{matrix} (or \code{data.frame}) of
  predictor variables of the reference data (observations
  in rows and variables in columns).}

  \item{Yr}{a numeric \code{vector} containing the values
  of the response variable corresponding to the reference
  data}

  \item{Xu}{input \code{matrix} (or \code{data.frame}) of
  predictor variables of the data to be predicted
  (observations in rows and variables in columns).}

  \item{Yu}{an optional numeric \code{vector} containing
  the values of the response variable corresponding to the
  data to be predicted}

  \item{mblCtrl}{a list (created with the
  \code{\link{mblControl}} function) which contains some
  parameters that control the some aspects of the
  \code{mbl} function. See the \code{\link{mblControl}}
  function for more details.}

  \item{dissimilarityM}{(optional) a dissimilarity matrix.
  This argument can be used in case a user-defined
  dissimilarity matrix is preferred over the automatic
  dissimilarity matrix computation specified in the
  \code{sm} argument of the \code{\link{mblControl}}
  function. When \code{dissUsage = "predictors"},
  \code{dissimilarityM} must be a square symmetric
  dissimilarity matrix (derived from a matrix of the form
  \code{rbind(Xr, Xu)}) for which the diagonal values are
  zeros (since the dissimilarity between an object and
  itself must be 0).  On the other hand if \code{dissUsage}
  is set to either \code{"weights"} or \code{"none"},
  \code{dissimilarityM} must be a \code{matrix}
  representing the dissimilarity of each element in
  \code{Xu} to each element in \code{Xr}. The number of
  columns of the object correspondent to
  \code{dissimilarityM} must be equal to the number of rows
  in \code{Xu} and the number of rows equal to the number
  of rows in \code{Xr}.  If both \code{dissimilarityM} and
  \code{sm} are specified, only the \code{dissimilarityM}
  argument will be taken into account.}

  \item{dissUsage}{specifies how the dissimilarity
  information shall be used. The possible options are:
  \code{"predictors"}, \code{"weights"} and \code{"none"}
  (see details below).  Default is "predictors".}

  \item{k}{a numeric (integer) \code{vector} containing the
  sequence of k nearest neighbours to be tested. Either
  \code{k} or \code{k.diss} must be specified. Numbers with
  decimal values will be coerced to their next higher
  integer values. This vector will be automatically sorted
  into ascending order.}

  \item{k.diss}{a \code{vector} containing the sequence of
  dissimilarity thresholds to be tested. When the
  dissimilarity between a sample in \code{Xr} and a sample
  \code{Xu} is below the given threshold, the sample in
  sample in \code{Xr} is treated as a neighbour of the
  sample in \code{Xu}, otherwise it is ignored.  These
  thresholds depend on the corresponding dissimilarity
  measure specified in \code{sm}. Either \code{k} or
  \code{k.diss} must be specified.}

  \item{k.range}{a vector of length 2 which specifies the
  minimum (first value) and the maximum (second value)
  number of neighbours allowed when the \code{k.diss}
  argument is used.}

  \item{method}{a character indicating the method to be
  used at each local multivariate regression. Options are:
  \code{"pls"}, \code{"wapls1"}, \code{"wapls2"} and
  \code{"gpr"} (see details below).}

  \item{pls.c}{the number of pls components to be used in
  the regressions if one of the following methods is used:
  \code{"pls"}, \code{"wapls1"} or \code{"wapls2"}.  When
  \code{"pls"} is used, this argument must be a single
  numerical value. When either \code{"wapls1"} or
  \code{"wapls2"} are used this argument must be a vector
  of length 2 indicating the minimum (first value) and the
  maximum (second value) number of pls components used for
  the regressions (see details below).}

  \item{noise.v}{a value indicating the variance of the
  noise for Gaussian process regression. Default is 0.001.}

  \item{...}{additional arguments to be passed to other
  functions.}
}
\value{
a \code{list} of class \code{mbl} with the following
components (sorted by either \code{k} or \code{k.diss}
according to the case): \itemize{ \item{\code{call}}{ the
call used.} \item{\code{cntrlParam}}{ the list with the
control parameters used. If one or more control parameters
were reset automatically, then a list containing a list
with the initial control parameters specified and a list
with the parameters which were finally used.}
\item{\code{dissimilarities}}{ a list with the method used
to obtain the dissimilarity matrices and the dissimilarity
matrices corresponding to \eqn{D(Xr, Xu)} and
\eqn{D(Xr,Xr)} if \code{dissUsage = "predictors"}. This
object is returned only if the \code{returnDiss} argument
in the \code{mblCtrl} list was set to \code{TRUE} in the
the call used.} \item{\code{totalSamplesPredicted}}{ the
total number of samples predicted.}
\item{\code{pcAnalysis}}{ a list containing the results of
the principal component analysis. The first two objects
(\code{scores_Xr} and \code{scores_Xu}) are the scores of
the \code{Xr} and \code{Xu} matrices. It also contains the
number of principal components used
(\code{n.componentsUsed}) and another object which is a
\code{vector} containing the standardized Mahalanobis
dissimilarities (also called GH, Global H distance) between
each sample in \code{Xu} and the centre of \code{Xr}.}
\item{\code{components}}{ a list containing either the
number of principal components or partial least squares
components used for the computation of the orthogonal
dissimilarities. This object is only returned if the
dissimilarity meausre specified in \code{mblCtrl} is any of
the following options: \code{'pc'}, \code{'loc.pc'},
\code{"pls"}, \code{'loc.pls'}. If any of the local
orthogonal dissimilarities was used (\code{'loc.pc'} or
\code{"pls"}) a \code{data.frame} is also returned in his
list. This object is equivalent to the
\code{loc.n.components} object returned by the
\code{\link{orthoDiss}} function. It specifies the number
of local components (either principal components or partial
least squares components) used for computing the
dissimilarity between each query sample and its neighbour
samples, as returned by the \code{\link{orthoDiss}}
function.  } \item{\code{nnValStats}}{ a data frame
containing the statistics of the nearest neighbour
cross-validation for each either \code{k} or \code{k.diss}
depending on the arguments specified in the call. It is
returned only if \code{'NNv'} or \code{'both'} were
selected as validation method}
\item{\code{localCrossValStats}}{ a data frame containing
the statistics of the local leave-group-out cross
validation for each either \code{k} or \code{k.diss}
depending on the arguments specified in the call. It is
returned only if \code{'local_crossval'} or \code{'both'}
were selected as validation method}
\item{\code{YuPredictionStats}}{ a data frame containing
the statistics of the cross-validation of the prediction of
\code{Yu} for each either \code{k} or \code{k.diss}
depending on the arguments specified in the call. It is
returned only if \code{Yu} was provided.}
\item{\code{results}}{ a list of data frames which contains
the results of the predictions for each either \code{k} or
\code{k.diss}.} }
}
\description{
This function is implemented for memory-based learning
(a.k.a. instance-based learning or local regression) which
is a non-linear lazy learning approach for predicting a
given response variable from a set of (spectral) predictor
variables. For each sample in an prediction set a specific
local regression is carried out based on a subset of
similar samples (nearest neighbours) selected from a
reference set. The local model is then used to predict the
response value of the target (prediction) sample. Therefore
this function does not yield a global regression model.
}
\details{
\code{dissUsage} is used to specifiy wheter or not and how
to use dissimilarity information for local regressions.
When \code{dissUsage = "predictors"} the local (square
symmetric) dissimilarity matrix corresponding the selected
neighbourhood is used as source of additional predictors
(i.e the columns of this local matrix are treated as
predictor variables). In some cases this may result in an
improvement of the prediction performance (Ramirez-Lopez et
al., 2013a). If \code{dissUsage = "weights"}, the
neighbours of the query point (\eqn{xu_{j}}) are weighted
according to their dissimilarity (e.g. distance) to
\eqn{xu_{j}} prior carrying out each local regression. The
following tricubic function (Cleveland and Delvin, 1988;
Naes et al., 1990) is used for computing the final weights
based on the measured dissimilarities: \deqn{W_{j} = (1 -
v^{3})^{3}} where if \eqn{xr_{i} \in} neighbours of
\eqn{xu_{j}}: \deqn{v_{j}(xu_{j}) = d(xr_{i}, xu_{j})}
otherwise: \deqn{v_{j}(xu_{j}) = 0} In the above formulas
\eqn{d(xr_{i}, xu_{j})} represents the dissimilarity
between the query point and each object in \eqn{Xr}. When
\code{dissUsage = "none"} is chosen the dissimilarity
information is not used. The possible options for
performing regressions at each local segment implemented in
the \code{mbl} function are described as follows: \itemize{
\item{Partial least squares (\code{"pls"}):}{ It uses the
\code{'oscorespls'} algorithm (equivalent to the non-linear
iterative partial least sqaures algorithm) option
implemented in the \code{\link[pls]{plsr}} function of the
\code{pls} package. The only parameter which needs to be
optimized is the number of pls components. This can be done
by cross-validation at each local segment.} \item{Weighted
average pls 1 (\code{"wapls1"}):}{ It uses multiple models
generated by multiple pls components (i.e. between a
minimum and a maximum number of pls components). At each
local partition the final predicted value is a weighted
average of all the predicted values generated by the
multiple pls models. The weight for each component is
calculated as follows: \deqn{ w_{j} =
\frac{1}{s_{1:j}\times g_{j}} } where \eqn{s_{1:j}} is the
root mean square of the spectral residuals of the unknown
(or target) sample when a total of \eqn{j} pls components
are used and \eqn{g_{j}} is the root mean square of the
regression coefficients corresponding to the \eqn{j}th pls
component (see Shenk et al., 1997 and Zhang et al., 2004
for more details). \code{"wapls1"} is not compatible with
\code{valMethod = "loc_crossval"} since the weights are
computed based on the sample to be predicted at each local
iteration.} \item{Weighted average pls 2
(\code{"wapls2"}):}{ It uses multiple models generated by
multiple pls components (i.e. between a minimum and a
maximum number of pls components). At each local partition
the final predicted value is a weighted average of all the
predicted values generated by the multiple pls models. The
weights are calculated according to the root mean square
error (RMSE) of internal cross-validations calculated for
each pls component as described in Zhang et al. (2004).
The equation for computing the weights is as follows:
\deqn{ w_{j} = \frac{1}{RMSE_{j}} } where \eqn{j} is the
\eqn{j}th pls component.} \item{Gaussian process with dot
product covariance (\code{"gpr"}):}{ Gaussian process
regression is a probabilistic and non-parametric Bayesian
approach. It is commonly described as a collection of
random variables which have a joint Gaussian distribution
and it is characterized by both a mean and a covariance
function (Williams and Rasmussen, 1996). The covariance
function used in the implemented method is the dot product,
which inplies that there are no parameters to be optimized
for the computation of the covariance.  Here, the process
for predicting the response variable of a new sample
(\eqn{y_{new}}) from its predictor variables
(\eqn{x_{new}}) is carried out first by computing a
prediction vector (\eqn{A}). It is derived from a set of
reference spectra (\eqn{X}) and their respective response
vector (\eqn{Y}) as follows: \deqn{ A = (X X^\textup{T} +
\sigma^2 I)^{-1} Y } where \eqn{\sigma^2} denotes the
variance of the noise and \eqn{I} the identity matrix (with
dimensions equal to the number of observations in \eqn{X}).
The prediction of \eqn{y_{new}} is then carried out by:
\deqn{ y_{new} = (x_{new}x_{new}^\textup{T}) A } } } The
loop used to iterate over the \code{Xu} samples in
\code{mbl} uses the \code{\%dopar\%} operator of the
\code{\link{foreach}} package, which can be used to
parallelize this internal loop. The last example given in
the \code{\link{mbl}} function ilustrates how to
parallelize the \code{\link{mbl}} function.
}
\examples{
\dontrun{
require(prospectr)

data(NIRsoil)

# Filter the data using the Savitzky and Golay smoothing filter with
# a window size of 11 spectral variables and a polynomial order of 3
# (no differentiation).
sg <- savitzkyGolay(NIRsoil$spc, p = 3, w = 11, m = 0)

# Replace the original spectra with the filtered ones
NIRsoil$spc <- sg

Xu <- NIRsoil$spc[!as.logical(NIRsoil$train),]
Yu <- NIRsoil$CEC[!as.logical(NIRsoil$train)]

Yr <- NIRsoil$CEC[as.logical(NIRsoil$train)]
Xr <- NIRsoil$spc[as.logical(NIRsoil$train),]

Xu <- Xu[!is.na(Yu),]
Xr <- Xr[!is.na(Yr),]

Yu <- Yu[!is.na(Yu)]
Yr <- Yr[!is.na(Yr)]

# Example 1
# A mbl approach (the spectrum-based learner) as implemented
# in Ramirez-Lopez et al. (2013)
# Example 1.1
# An exmaple where Yu is supposed to be unknown, but the Xu
# (spectral variables) are known
ctrl1 <- mblControl(sm = "pc", pcSelection = list("opc", 40),
                    valMethod = "NNv",
                    scaled = TRUE, center = TRUE)

sbl.u <- mbl(Yr = Yr, Xr = Xr, Yu = NULL, Xu = Xu,
             mblCtrl = ctrl1,
             dissUsage = "predictors",
             k = seq(40, 150, by = 10),
             method = "gpr")
sbl.u


# Example 1.2
# If Yu is actually known...
sbl.u2 <- mbl(Yr = Yr, Xr = Xr, Yu = Yu, Xu = Xu,
              mblCtrl = ctrl1,
              dissUsage = "predictors",
              k = seq(40, 150, by = 10),
              method = "gpr")
sbl.u2

# Example 1.3
# A variation of the spectrum-based learner implemented in
# Ramirez-Lopez et al. (2013)where the dissimilarity matrices are
# recomputed based on partial least squares scores
ctrl_1.3 <- mblControl(sm = "pls", pcSelection = list("opc", 40),
                       valMethod = "NNv",
                       scaled = TRUE, center = TRUE)

sbl_1.3 <- mbl(Yr = Yr, Xr = Xr, Yu = Yu, Xu = Xu,
               mblCtrl = ctrl_1.3,
               dissUsage = "predictors",
               k = seq(40, 150, by = 10),
               method = "gpr",
               valMethod = "NNv")
sbl_1.3

# Example 2
# A mbl approach similar to the ones implemented in
# Ramirez-Lopez et al. (2013)
# and Fernandez Pierna and Dardenne (2008)
ctrl.mbl <- mblControl(sm = "cor",
                       pcSelection = list("cumvar", 0.999),
                       valMethod = "NNv",
                       scaled = TRUE, center = TRUE)

local.mbl <- mbl(Yr = Yr, Xr = Xr, Yu = Yu, Xu = Xu,
                 mblCtrl = ctrl.mbl,
                 dissUsage = "none",
                 k = seq(40, 150, by = 10),
                 pls.c = c(7, 20),
                 method = "wapls1",
                 valMethod = "NNv")
local.mbl

# Example 3
# A WA-LOCAL approach as implemented in Zhang et al. (2004)
ctrl.wa <- mblControl(sm = "cor",
                      pcSelection = list("cumvar", 0.999),
                      valMethod = c("NNv", "loc_crossval"),
                      resampling = 10, p = 0.75,
                      scaled = TRUE, center = TRUE)

wa.local <- mbl(Yr = Yr, Xr = Xr, Yu = Yu, Xu = Xu,
                mblCtrl = ctrl.wa,
                dissUsage = "none",
                k = seq(40, 150, by = 10),
                pls.c = c(7, 20),
                method = "wapls2")
wa.local

# Example 4
# Using the function with user-defined dissimilarities
# Examples 4.1 - 4.2: Compute a square symetric matrix of
# dissimilarities between
# all the elements in Xr and Xu (dissimilarities will be used as
# additional predictor variables later in the mbl function)
# Examples 4.3 - 4.4: Derive a dissimilarity value of each element
# in Xu to each element in Xr (in this case dissimilarities will
# not be used as additional predictor variables later in the
# mbl function)
# Example 4.1
# the manhattan distance
manhattanD <- dist(rbind(Xr, Xu), method = "manhattan")
manhattanD <- as.matrix(manhattanD)

ctrl.udd <- mblControl(sm = "none",
                       pcSelection = list("cumvar", 0.999),
                       valMethod = c("NNv", "loc_crossval"),
                       resampling = 10, p = 0.75,
                       scaled = TRUE, center = TRUE)

mbl.udd1 <- mbl(Yr = Yr, Xr = Xr, Yu = Yu, Xu = Xu,
                mblCtrl = ctrl.udd,
                dissimilarityM = manhattanD,
                dissUsage = "predictors",
                k = seq(40, 150, by = 10),
                method = "gpr")
mbl.udd1

#Example 4.2
# first derivative spectra
der.sp <- t(diff(t(rbind(Xr, Xu)), lag = 1, differences = 1))

# The euclidean dissimilarity on the derivative spectra
# (a.k.a spectral dissimilarity)
spc.dist <- fDiss(Xr = der.sp, method = "euclid",
                  center = FALSE, scale = FALSE)

mbl.udd2 <- mbl(Yr = Yr, Xr = Xr, Yu = Yu, Xu = Xu,
                mblCtrl = ctrl.udd,
                dissimilarityM = spc.dist,
                dissUsage = "predictors",
                k = seq(40, 150, by = 10),
                method = "gpr")

#Example 4.3
# first derivative spectra
der.Xr <- t(diff(t(Xr), lag = 1, differences = 1))
der.Xu <- t(diff(t(Xu), lag = 1, differences = 1))
# the sid on the derivative spectra
der.sid <- sid(Xr = der.Xr, X2 = der.Xu, mode = "density",
               center = TRUE, scaled = TRUE)
der.sid <- der.sid$sid

mbl.udd3 <- mbl(Yr = Yr, Xr = Xr, Yu = Yu, Xu = Xu,
                mblCtrl = ctrl.udd,
                dissimilarityM = der.sid,
                dissUsage = "none",
                k = seq(40, 150, by = 10),
                method = "gpr")
mbl.udd3

# Example 5
# For running the sbl function in parallel
n.cores <- 2   # two cores

# Set the number of cores according to the OS
if (.Platform[["OS.type"]] == "windows") {
  library(doParallel)
  cl <- makeCluster(n.cores)
  registerDoParallel(cl)
}else{
  library(doSNOW)
  cluster <- makeCluster(n.cores, type = "SOCK")
  registerDoSNOW(cluster)
  ncores <- getDoParWorkers()
}

ctrl <- mblControl(sm = "pc", pcSelection = list("opc", 40),
                   valMethod = "NNv",
                   scaled = TRUE, center = TRUE)

mbl.p <- mbl(Yr = Yr, Xr = Xr, Yu = Yu, Xu = Xu,
             mblCtrl = ctrl,
             dissUsage = "none",
             k = seq(40, 150, by = 10),
             method = "gpr")
registerDoSEQ()
try(stopCluster(cl))
mbl.p
}
}
\author{
Leonardo Ramirez-Lopez and Antoine Stevens
}
\references{
Cleveland, W. S., and Devlin, S. J. 1988. Locally weighted
regression: an approach to regression analysis by local
fitting. Journal of the American Statistical Association,
83, 596-610.

Fernandez Pierna, J.A., Dardenne, P. 2008. Soil parameter
quantification by NIRS as a Chemometric challenge at
"Chimiomitrie 2006". Chemometrics and Intelligent
Laboratory Systems 91, 94-98

Naes, T., Isaksson, T., Kowalski, B. 1990. Locally weighted
regression and scatter correction for near-infrared
reflectance data. Analytical Chemistry 62, 664-673.

Ramirez-Lopez, L., Behrens, T., Schmidt, K., Stevens, A.,
Dematte, J.A.M., Scholten, T. 2013a. The spectrum-based
learner: A new local approach for modeling soil vis-NIR
spectra of complex datasets. Geoderma 195-196, 268-279.

Ramirez-Lopez, L., Behrens, T., Schmidt, K., Viscarra
Rossel, R., Dematte, J. A. M., Scholten, T. 2013b. Distance
and similarity-search metrics for use with soil vis-NIR
spectra. Geoderma 199, 43-53.

Rasmussen, C.E., Williams, C.K. Gaussian Processes for
Machine Learning. Massachusetts Institute of Technology:
MIT-Press, 2006.

Shenk, J., Westerhaus, M., and Berzaghi, P. 1997.
Investigation of a LOCAL calibration procedure for near
infrared instruments. Journal of Near Infrared
Spectroscopy, 5, 223-232.

Zhang, M.H., Xu, Q.S., Massart, D.L. 2004. Averaged and
weighted average partial least squares. Analytica Chimica
Acta 504, 279-289.
}
\seealso{
\code{\link{fDiss}}, \code{\link{corDiss}},
\code{\link{sid}}, \code{\link{orthoDiss}},
\code{\link{neigCleaning}}
}
