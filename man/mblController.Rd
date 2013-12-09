\name{mblController}
\alias{mblController}
\title{A function that controls some aspects of the memory-based learning process in the \code{mbl} function}
\usage{
  mblController(sm = 'pc', pcSelection = list('opc', 40),
    ws = if(sm == 'movcor') 41, k0, returnDiss = FALSE,
    center = TRUE, scaled = TRUE, valMethod = c('NNv',
    'loc_crossval'), resampling = 10, p = 0.75, progress =
    TRUE)
}
\arguments{
  \item{sm}{a character string indicating the spectral
  similarity/dissimilarity metric to be used in the
  selection of the nearest neighbours of each observation
  in \code{'Xu'} (see the \code{\link{mbl}} function).
  Options are: \code{'euclid'} (Euclidean distance),
  \code{'cosine'} (cosine distance), \code{'sidF'}
  (spectral information divergence computed on the spectral
  variables), \code{'sidD'} (spectral information
  divergence computed on the density distributions of the
  spectra) \code{'cor'} (correlation
  similarity/dissimilarity), \code{'movcor'} (moving window
  correlation similarity/dissimilarity), \code{'pc'}
  (principal components distance: Mahalanobis distance
  computed on the principal components space) and
  \code{'loc.pc'} (similarity/dissimilarity estimation
  based on local principal components), \code{'pls'}
  (partial least squares distance: Mahalanobis distance
  computed on the partial least squares space) and
  \code{'loc.pls'} (similarity/dissimilarity estimation
  based on local partial least squares). The \code{'pc'}
  spectral similarity/dissimilarity metric is the default.
  If the \code{'sidD'} is chosen, the default parameters of
  the \code{sid} function are used. This argument can also
  be set to \code{NULL}, in such a case, a
  similarity/dissimilarity matrix must be specified in the
  \code{distanceM} argument of the \code{\link{mbl}}
  function.}

  \item{pcSelection}{a list which specifies the method to
  be used for identifying the number of principal
  components to be retained for computing the Mahalanobis
  distance of each sample in \code{sm = 'Xu'} to the centre
  of \code{sm = 'Xr'}. It also specifies the number of
  components in any of the following cases: \code{sm =
  'pc'}, \code{sm = 'loc.pc'}, \code{sm = 'pls'} and
  \code{sm = 'loc.pls'}. This list must contain two objects
  in the following order: \itemize{
  \item{\code{method}:}{the method for selecting the number
  of components. Possible options are: \code{'opc'}
  (optimized pc selection based on Ramirez-Lopez et al.
  (2013a, 2013b). See the \code{\link{orthoProjection}}
  function for more details; \code{'cumvar'} (for selecting
  the number of principal components based on a given
  cumulative amount of explained variance); \code{'var'}
  (for selecting the number of principal components based
  on a given amount of explained variance); and
  \code{'manual'} (for specifying manually the desired
  number of principal components)} \item{\code{value}:}{a
  numerical value that complements the selected method. If
  \code{'opc'} is chosen, it must be a value indicating the
  maximal number of principal components to be tested (see
  Ramirez-Lopez et al., 2013a, 2013b). If \code{'cumvar'}
  is chosen, it must be a value (higher than 0 and lower
  than 1) indicating the maximum amount of cumulative
  variance that the retained components should explain. If
  \code{'var'} is chosen, it must be a value (higher than 0
  and lower than 1) indicating that components that explain
  (individually) a variance lower than this threshold must
  be excluded. If \code{'manual'} is chosen, it must be a
  value specifying the desired number of principal
  components to retain.  }} The default method for the
  \code{pcSelection} argument is \code{'opc'} and the
  maximal number of principal components to be tested is
  set to 40.  Optionally, the \code{pcSelection} argument
  admits \code{'opc'} or \code{'cumvar'} or \code{'var'} or
  \code{'manual'} as a single character string. In such a
  case the default for \code{'value'} when either
  \code{'opc'} or \code{'manual'} are used is 40. When
  \code{'cumvar'} is used the default \code{'value'} is set
  to 0.99 and when \code{'var'} is used the default
  \code{'value'} is set to 0.01.}

  \item{ws}{an odd integer value which specifies the window
  size when the moving window correlation
  similarity/dissimilarity is used (i.e \code{sm =
  'movcor'}). The default is 41.}

  \item{k0}{if any of the local similarity/dissimilarity
  methods is used (i.e. either \code{sm = 'loc.pc'} or
  \code{sm = 'loc.pls'}) a numeric integer value. This
  argument controls the number of initial
  neighbours(\eqn{k0}) to retain in order to compute the
  local principal components (at each neighbourhood).}

  \item{returnDiss}{a logical indicating if the distance
  matrices must be returned.}

  \item{center}{a logical indicating whether or not the
  predictor variables must be centered at each local
  segment (before regression).}

  \item{scaled}{a logical indicating whether or not the
  predictor variables must be scaled at each local segment
  (before regression).}

  \item{valMethod}{a character vector which indicates the
  (internal) validation method(s) to be used for assessing
  the global performance of the local models. Possible
  options are: \code{'NNv'} and \code{'loc_crossval'}. For
  no validation \code{'none'} is used (see details below).}

  \item{resampling}{a value indicating the number of
  resampling iterations at each local segment when
  \code{'loc_crossval'} is selected in the \code{valMethod}
  argument. Default is 10.}

  \item{p}{a value indicating the percentage of samples to
  be retained in each resampling iteration at each local
  segment when \code{'loc_crossval'} is selected in the
  \code{valMethod} argument. Default is 0.75 (i.e. 75
  "\%")}

  \item{progress}{a logical indicating whether or not to
  print a progress bar for each sample to be predicted.
  Default is \code{TRUE}. Note: In case multicore
  processing is used, this progress bar will not be
  printed.}
}
\value{
  \code{mblController} returns a \code{list} of class
  \code{mbl} with the specified parameters
}
\description{
  This function is used to specify various aspects in the
  memory-based learning process in the \code{mbl} function
}
\details{
  The validation methods avaliable for assessing the
  predictive performance of the memory-based learning
  method used are described as follows: \itemize{
  \item{Leave-nearest-neighbour-out cross validation
  (\code{'NNv'}):}{ From the group of neighbours of each
  sample to be predicted, the nearest sample (i.e. the most
  similar sample) is excluded and then a local model is
  fitted using the remaining neighbours. This model is then
  used to predict the value of the target response variable
  of the nearest sample. These predicted values are finally
  cross validated with the actual values (See Ramirez-Lopez
  et al. (2013a) for additional details). This method is
  faster than \code{"loc_crossval"}} \item{Local
  leave-group-out cross validation
  (\code{'loc_crossval'}):}{ The group of neighbours of
  each sample to be predicted is partitioned into different
  equal size subsets. Each partition is selected based on a
  stratified random sampling which takes into account the
  values of the response variable of the corresponding set
  of neighbours. The selected local subset is used as local
  validation subset and the remaining samples are used for
  fitting a model. This model is used to predict the target
  response variable values of the local validation subset
  and the local root mean square error is computed. This
  process is repeated \eqn{m} times and the final local
  error is computed as the average of the local root mean
  square error of all the \eqn{m} iterations. In the
  \code{mbl} function \eqn{m} is controlled by the
  \code{resampling} argument and the size of the subsets is
  controlled by the \code{p} argument which indicates the
  percentage of samples to be selected from the subset of
  nearest neighbours. The global error of the predictions
  is computed as the average of the local root mean square
  errors.} \item{No validation (\code{'none'}):}{ No
  validation is carried out. If \code{'none'} is seleceted
  along with \code{'NNv'} and/or \code{'loc_crossval'},
  then it will be ignored and the respective validation(s)
  will be carried out.} }
}
\examples{
#A control list with the default parameters
mblController()

#A control list which specifies the moving/rolling correlation similarity/dissimilarity metric with a moving window of 30
mblController(sm = "movcor", ws = 31)
}
\author{
  Leonardo Ramirez-Lopez and Antoine Stevens
}
\references{
  Ramirez-Lopez, L., Behrens, T., Schmidt, K., Stevens, A.,
  Dematte, J.A.M., Scholten, T. 2013a. The spectrum-based
  learner: A new local approach for modeling soil vis-NIR
  spectra of complex datasets. Geoderma 195-196, 268-279.

  Ramirez-Lopez, L., Behrens, T., Schmidt, K., Viscarra
  Rossel, R., Dematte, J. A. M., Scholten, T. 2013b.
  Distance and similarity-search metrics for use with soil
  vis-NIR spectra. Geoderma 199, 43-53.
}
\seealso{
  \code{\link{fDiss}}, \code{\link{corDiss}},
  \code{\link{sid}}, \code{\link{orthoDiss}},
  \code{\link{mbl}}
}

