\name{orthoDiss}
\alias{orthoDiss}
\title{A function for computing orthogonal dissimilarity matrices (orthoDiss)}
\usage{
orthoDiss(Xr, X2 = NULL,
          Yr = NULL,
          pcSelection = list("cumvar", 0.99),
          method = "pca",
          local = FALSE,
          k0,
          center = TRUE, scaled = TRUE,
          return.all = FALSE, cores = 1, ...)
}
\arguments{
  \item{Xr}{a \code{matrix} (or \code{data.frame})
  containing the (reference) data.}

  \item{X2}{an optional \code{matrix} (or
  \code{data.frame}) containing data of a second set of
  observations(samples).}

  \item{Yr}{either if the method used in the
  \code{pcSelection} argument is \code{"opc"} or if the
  \code{sm} argument is either \code{"pls"} or
  \code{"loc.pls"}, then it must be a \code{vector}
  containing the side information corresponding to the
  spectra in \code{Xr}. It is equivalent to the
  \code{sideInf} parameter of the \code{\link{simEval}}
  function. It can be a numeric \code{vector} or
  \code{matrix} (regarding one or more continuous
  variables). The root mean square of differences (rmsd) is
  used for assessing the similarity between the samples and
  their corresponding most similar samples in terms of the
  side information provided. When \code{sm = "pc"}, this
  parameter can also be a single discrete variable of class
  \code{factor}. In such a case the kappa index is used.
  See \code{\link{simEval}} function for more details.}

  \item{pcSelection}{a list which specifies the method to
  be used for identifying the number of principal
  components to be retained for computing the Mahalanobis
  distance of each sample in \code{sm = "Xu"} to the centre
  of \code{sm = "Xr"}. It also specifies the number of
  components in any of the following cases: \code{sm =
  "pc"}, \code{sm = "loc.pc"}, \code{sm = "pls"} and
  \code{sm = "loc.pls"}. This list must contain two objects
  in the following order: \itemize{
  \item{\code{method}:}{the method for selecting the number
  of components. Possible options are: \code{"opc"}
  (optimized pc selection based on Ramirez-Lopez et al.
  (2013a, 2013b). See the \code{\link{orthoProjection}}
  function for more details; \code{"cumvar"} (for selecting
  the number of principal components based on a given
  cumulative amount of explained variance); \code{"var"}
  (for selecting the number of principal components based
  on a given amount of explained variance); and
  \code{"manual"} (for specifying manually the desired
  number of principal components)} \item{\code{value}:}{a
  numerical value that complements the selected method. If
  \code{"opc"} is chosen, it must be a value indicating the
  maximal number of principal components to be tested (see
  Ramirez-Lopez et al., 2013a, 2013b). If \code{"cumvar"}
  is chosen, it must be a value (higher than 0 and lower
  than 1) indicating the maximum amount of cumulative
  variance that the retained components should explain. If
  \code{"var"} is chosen, it must be a value (higher than 0
  and lower than 1) indicating that components that explain
  (individually) a variance lower than this threshold must
  be excluded. If \code{"manual"} is chosen, it must be a
  value specifying the desired number of principal
  components to retain.  }} The default method for the
  \code{pcSelection} argument is \code{"opc"} and the
  maximal number of principal components to be tested is
  set to 40.  Optionally, the \code{pcSelection} argument
  admits \code{"opc"} or \code{"cumvar"} or \code{"var"} or
  \code{"manual"} as a single character string. In such a
  case the default for \code{"value"} when either
  \code{"opc"} or \code{"manual"} are used is 40. When
  \code{"cumvar"} is used the default \code{"value"} is set
  to 0.99 and when \code{"var"} is used the default
  \code{"value"} is set to 0.01.}

  \item{method}{the method for projecting the data. Options
  are: "pca" (principal component analysis using the
  singular value decomposition algorithm), "pca.nipals"
  (principal component analysis using the non-linear
  iterative partial least squares algorithm) and "pls"
  (partial least squares). See the
  \code{\link{orthoProjection}} function for further
  details on the projection methods.}

  \item{local}{a logical indicating whether or not to
  compute the distances locally (i.e. projecting locally
  the data) by using the \eqn{k0} nearest neighbour samples
  of each sample. Default is \code{FALSE}. See details.}

  \item{k0}{if \code{local = TRUE} a numeric integer value
  which indicates the number of nearest
  neighbours(\eqn{k0}) to retain in order to recompute the
  local orthogonal distances.}

  \item{center}{a logical indicating if the spectral data
  \code{Xr} (and \code{X2} if specified) must be centered.
  If \code{X2} is specified the data is scaled on the basis
  of \eqn{Xr \cup Xu}.}

  \item{scaled}{a logical indicating if \code{Xr} (and
  \code{X2} if specified) must be scaled. If \code{X2} is
  specified the data is scaled on the basis of \eqn{Xr \cup
  Xu}.}

  \item{return.all}{a logical. In case \code{X2} is
  specified it indicates whether or not the distances
  between all the elements resulting from \eqn{Xr \cup Xu}
  must be computed.}

  \item{cores}{number of cores used when \code{method} in
  \code{pcSelection} is \code{"opc"} (which can be
  computationally intensive) and \code{local = FALSE}
  (default = 1). Dee details.}

  \item{...}{additional arguments to be passed to the
  \code{\link{orthoProjection}} function.}
}
\value{
a \code{list} of class \code{orthoDiss} with the following
components: \itemize{ \item{\code{n.components}}{ the
number of components (either principal components or
partial least squares components) used for computing the
global distances.} \item{\code{loc.n.components}}{ if
\code{local = TRUE}, a \code{data.frame} which specifies
the number of local components (either principal components
or partial least squares components) used for computing the
dissimilarity between each target sample and its neighbour
samples.} \item{\code{dissimilarity}}{ the computed
dissimilarity matrix. If \code{local = FALSE} a distance
\code{matrix}. If \code{local = TRUE} a \code{matrix} of
class \code{orthoDiss}. In this case each column represent
the dissimilarity between a target sample and its
neighbourhood.} } Multi-core processing works currently
only on windows and linux (at least for the C++ code used
in the package that use multi-threading).
}
\description{
This function computes orthogonal dissimilarities between
either observations in a given set or between observations
in two different sets. The dissimilarities are computed
based on either principal component projection or partial
least squares projection of the data. After projecting the
data, the Mahalanobis distance is applied.
}
\details{
When \code{local = TRUE}, first a global distance matrix is
computed based on the parameters specified. Then, by using
this matrix for each target observation, a given set of
nearest neighbours (\eqn{k0}) are identified. These
neighbours (together with the target observation) are
projected (from the original data space) onto a (local)
orthogonal space (using the same parameters specified in
the function). In this projected space the Mahalanobis
distance between the target sample and the neighbours is
recomputed. A missing value is assigned to the samples that
do not belong to this set of neighbours (non-neighbour
samples). In this case the dissimilarity matrix cannot be
considered as a distance metric since it does not
necessarily satisfies the symmetry condition for distance
matrices (i.e. given two samples \eqn{x_i} and \eqn{x_j},
the local dissimilarity (\eqn{d}) between them is relative
since generally \eqn{d(x_i, x_j) \neq d(x_j, x_i)}). On the
other hand, when \code{local = FALSE}, the dissimilarity
matrix obtained can be considered as a distance matrix.
}
\examples{
\dontrun{
require(prospectr)

data(NIRsoil)

Xu <- NIRsoil$spc[!as.logical(NIRsoil$train),]
Yu <- NIRsoil$CEC[!as.logical(NIRsoil$train)]
Yr <- NIRsoil$CEC[as.logical(NIRsoil$train)]
Xr <- NIRsoil$spc[as.logical(NIRsoil$train),]

Xu <- Xu[!is.na(Yu),]
Yu <- Yu[!is.na(Yu)]

Xr <- Xr[!is.na(Yr),]
Yr <- Yr[!is.na(Yr)]

# Computation of the orthogonal dissimilarity matrix using the
# default parameters
ex1 <- orthoDiss(Xr = Xr, X2 = Xu)

# Computation of a principal component dissimilarity matrix using
# the "opc" method for the selection of the principal components
ex2 <- orthoDiss(Xr = Xr, X2 = Xu,
                 Yr = Yr,
                 pcSelection = list("opc", 40),
                 method = "pca",
                 return.all = TRUE)

# Computation of a partial least squares (PLS) dissimilarity
# matrix using the "opc" method for the selection of the PLS
# components
ex3 <- orthoDiss(Xr = Xr, X2 = Xu,
                 Yr = Yr,
                 pcSelection = list("opc", 40),
                 method = "pls")

# Computation of a partial least squares (PLS) local dissimilarity
# matrix using the "opc" method for the selection of the PLS
# components
ex4 <- orthoDiss(Xr = Xr, X2 = Xu,
                 Yr = Yr,
                 pcSelection = list("opc", 40),
                 method = "pls",
                 local = TRUE,
                 k0 = 200)
}
}
\author{
Leonardo Ramirez-Lopez
}
\references{
Ramirez-Lopez, L., Behrens, T., Schmidt, K., Stevens, A.,
Dematte, J.A.M., Scholten, T. 2013a. The spectrum-based
learner: A new local approach for modeling soil vis-NIR
spectra of complex datasets. Geoderma 195-196, 268-279.

Ramirez-Lopez, L., Behrens, T., Schmidt, K., Viscarra
Rossel, R., Dematte, J. A. M., Scholten, T. 2013b. Distance
and similarity-search metrics for use with soil vis-NIR
spectra. Geoderma 199, 43-53.
}
\seealso{
\code{\link{orthoProjection}}, \code{\link{simEval}}
}

