\name{fDiss}
\alias{fDiss}
\title{Euclidean, Mahalanobis and cosine dissimilarity measurements}
\usage{
  fDiss(Xr, X2 = NULL, method = 'euclid', center = TRUE,
    scaled = TRUE)
}
\arguments{
  \item{Xr}{a \code{matrix} (or \code{data.frame})
  containing the (reference) data.}

  \item{X2}{an optional \code{matrix} (or
  \code{data.frame}) containing data of a second set of
  observations(samples).}

  \item{method}{the method for computing the dissimilarity
  matrix. Options are \code{'euclid'} (Euclidean distance),
  \code{'mahalanobis'} (Mahalanobis distance) and
  \code{'cosine'} (cosine distance, a.k.a spectral angle
  mapper).}

  \item{center}{a logical indicating if the spectral data
  \code{Xr} (and \code{X2} if specified) must be centered.
  If X2 is specified the data is scaled on the basis of
  \eqn{Xr \cup Xu}.}

  \item{scaled}{a logical indicating if \code{Xr} (and
  \code{X2} if specified) must be scaled. If X2 is
  specified the data is scaled on the basis of \eqn{Xr \cup
  Xu}.}
}
\value{
  a \code{matrix} of the computed dissimilarities.
}
\description{
  This function is used to compute the dissimilarity
  between observations based on Euclidean or Mahalanobis
  distance measures or on cosine dissimilarity measures
  (a.k.a spectral angle mapper).
}
\details{
  In the case of both the Euclidean and Mahalanobis
  distances, the dissimilarity matrix \eqn{D} between
  between samples in a given matrix \eqn{X} is computed as
  follows: \deqn{ D(x_i, x_j) = \sqrt{(x_i - x_j)M^{-1}(x_i
  - x_j)^{\mathrm{T}}} } where \eqn{M} is the identity
  matrix in the case of the Euclidean distance and the
  variance-covariance matrix of \eqn{M} in the case of the
  Mahalanobis distance. The Mahalanobis distance can also
  be viewed as the Euclidean distance after applying a
  linear transformation of the original variables. Such a
  linear transformation is carried by using a factorization
  of the inverse covariance matrix as \eqn{M^{-1} =
  W^{\mathrm{T}}W}, where \eqn{M} is merely the square root
  of \eqn{M^{-1}} which can be found by using a singular
  value decomposition. Note that when attempting to compute
  the Mahalanobis distance on a dataset with highly
  correlated variables (i.e. spectral variables) the
  variance-covariance matrix may result in a singular
  matrix which cannot be inverted and therefore the
  distance cannot be computed. This is also the case when
  the number of samples in the dataset is smaller than the
  number of variables. For the computation of the
  Mahalanobis distance, the mentioned method is used. On
  the other hand the cosine dissimilarity \eqn{S} between
  two obsvervations \eqn{x_i} and \eqn{x_j} is computed as
  follows: \deqn{ S(x_i, x_j) =
  cos^{-1}{\frac{\sum_{k=1}^{p}x_{i,k}
  x_{j,k}}{\sqrt{\sum_{k=1}^{p} x_{i,k}^{2}}
  \sqrt{\sum_{k=1}^{p} x_{j,k}^{2}}}} } where \eqn{p} is
  the number of variables of the observations.
}
\examples{
require(prospectr)
data(NIRsoil)

Xu <- NIRsoil$spc[!as.logical(NIRsoil$train),]
Xr <- NIRsoil$spc[as.logical(NIRsoil$train),]

# Euclidean distances between all the samples in Xr
ed <- fDiss(Xr = Xr, method = "euclid", center = TRUE, scaled = TRUE)

# Euclidean distances between samples in Xr and samples in Xu
ed.xr.xu <- fDiss(Xr = Xr, X2 = Xu, method = "euclid", center = TRUE, scaled = TRUE)

# Mahalanobis distance computed on the first 20 spectral variables
md.xr.xu <- fDiss(Xr = Xr[,1:20], X2 = Xu[,1:20], method = "mahalanobis", center = TRUE, scaled = TRUE)

# Cosine dissimilarity matrix
cdiss.xr.xu <- fDiss(Xr = Xr, X2 = Xu, method = "cosine", center = TRUE, scaled = TRUE)
}
\author{
  Leonardo Ramirez-Lopez and Antine Stevens
}

