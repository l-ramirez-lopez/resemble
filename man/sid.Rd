\name{sid}
\alias{sid}
\title{A function for computing the spectral information divergence between spectra (sid)}
\usage{
sid(Xr, X2 = NULL,
    mode = "density",
    center = FALSE, scaled = TRUE,
    kernel = "gaussian",
    n = if(mode == "density") round(0.5 * ncol(Xr)),
    bw = "nrd0",
    reg = 1e-04,
    ...)
}
\arguments{
  \item{Xr}{a \code{matrix} (or \code{data.frame})
  containing the spectral (reference) data.}

  \item{X2}{an optional \code{matrix} (or
  \code{data.frame}) containing the spectral data of a
  second set of samples.}

  \item{mode}{the method to be used for computing the
  spectral information divergence. Options are
  \code{"density"} (default) for computing the divergence
  values on the density distributions of the spectral
  observations, and \code{"feature"} for computing the
  divergence vales on the spectral features. See details.}

  \item{center}{a logical indicating if the computations
  must be carried out on the centered \code{X} and
  \code{X2} (if specified) matrices. If \code{mode =
  "feature"} centering is not carried out since this option
  does not accept negative values which are generated after
  centering the matrices. Default is FALSE. See details.}

  \item{scaled}{a logical indicating if the computations
  must be carried out on the variance scaled \code{X} and
  \code{X2} (if specified) matrices. Default is TRUE.}

  \item{kernel}{if \code{mode = "density"} a character
  string indicating the smoothing kernel to be used. It
  must be one of \code{"gaussian"} (default),
  \code{"rectangular"}, \code{"triangular"},
  \code{"epanechnikov"}, \code{"biweight"}, \code{"cosine"}
  or \code{"optcosine"}. See the
  \code{\link[stats]{density}} function of the \code{stats}
  package.}

  \item{n}{if \code{mode = "density"} a numerical value
  indicating the number of equally spaced points at which
  the density is to be estimated. See the
  \code{\link[stats]{density}} function of the \code{stats}
  package for further details. Default is \code{round(0.5 *
  ncol(X))}.}

  \item{bw}{if \code{mode = "density"} a numerical value
  indicating the smoothing kernel bandwidth to be used.
  Optionally the character string \code{"nrd0"} can be
  used, it computes the bandwidth using the
  \code{\link[stats]{bw.nrd0}} function of the \code{stats}
  package. See the \code{\link[stats]{density}} and the
  \code{\link[stats]{bw.nrd0}} functions for more details.
  By default \code{"nrd0"} is used, in this case the
  bandwidth is computed as \code{bw.nrd0(as.vector(X))}, if
  \code{X2} is specified the bandwidth is computed as
  \code{bw.nrd0(as.vector(rbind(X, X2)))}.}

  \item{reg}{a numerical value higher than 0 which
  indicates a regularization parameter. Values
  (probabilities) below this threshold are replaced by this
  value for numerical stability. Default is 1e-4.}

  \item{...}{additional arguments to be passed to the
  \code{\link[stats]{density}} function of the base
  package.}
}
\value{
a \code{list} with the following components: \itemize{
\item{\code{sid}}{ if only \code{"X"} is specified (i.e.
\code{X2 = NULL}), a square symmetric matrix of SID
distances between all the components in \code{"X"}. If both
\code{"X"} and \code{"X2"} are specified, a matrix of SID
distances between the components in \code{"X"} and the
components in \code{"X2"}) where the rows represent the
objects in \code{"X"} and the columns represent the objects
in \code{"X2"}} \item{\code{Xr}}{ the (centered and/or
scaled if specified) spectral \code{X} matrix}
\item{\code{X2}}{ the (centered and/or scaled if specified)
spectral \code{X2} matrix} \item{\code{densityDisXr}}{ if
\code{mode = "density"}, the computed density distributions
of \code{Xr}} \item{\code{densityDisX2}}{ if \code{mode =
"density"}, the computed density distributions of
\code{X2}} }
}
\description{
This function computes the spectral information divergence
(distance) between spectra based on the kullback-leibler
divergence algorithm (see details).
}
\details{
This function computes the spectral information divergence
(distance) between spectra. When \code{mode = "density"},
the function first computes the probability distribution of
each spectrum which result in a matrix of density
distribution estimates. The density distributions of all
the samples in the datasets are compared based on the
kullback-leibler divergence algorithm. When \code{mode =
"feature"}, the kullback-leibler divergence between all the
samples is computed directly on the spectral variables. The
spectral information divergence (SID) algorithm (Chang,
2000) uses the Kullback-Leibler divergence (\eqn{KL}) or
relative entropy (Kullback and Leibler, 1951) to account
for the vis-NIR information provided by each spectrum. The
SID between two spectra (\eqn{x_{i}} and \eqn{x_{j}}) is
computed as follows: \deqn{ SID(x_{i},x_{j}) = KL(x_{i}
\left |\right | x_{j}) + KL(x_{j} \left |\right | x_{i}) }
\deqn{ SID(x_{i},x_{j}) = \sum_{l=1}^{k} p_l \
log(\frac{p_l}{q_l}) + \sum_{l=1}^{k} q_l \
log(\frac{q_l}{p_l}) } where \eqn{k} represents the number
of variables or spectral features, \eqn{p} and \eqn{q} are
the probability vectors of \eqn{x_{i}} and \eqn{x_{j}}
respectively which are calculated as: \deqn{ p =
\frac{x_i}{\sum_{l=1}^{k} x_{i,l}} } \deqn{ q =
\frac{x_j}{\sum_{l=1}^{k} x_{j,l}} } From the above
equations it can be seen that the original SID algorithm
assumes that all the components in the data matrices are
nonnegative. Therefore centering cannot be applied when
\code{mode = "feature"}. If a data matrix with negative
values is provided and \code{mode = "feature"}, the
\code{sid} function automatically scales the matrix as
follows: \deqn{ X_s = \frac{X-min(X)}{max(X)-min(X)} } or
\deqn{ X_{s} = \frac{X-min(X, X2)}{max(X, X2)-min(X, X2)} }
\deqn{ X2_{s} = \frac{X2-min(X, X2)}{max(X, X2)-min(X, X2)}
} if \code{X2} is specified. The 0 values are replaced by a
regularization parameter (\code{reg} argument) for
numerical stability. The default of the \code{sid} function
is to compute the SID based on the density distributions of
the spectra (\code{mode = "density"}). For each spectrum in
\code{X} the density distribution is computed using the
\code{\link[stats]{density}} function of the \code{stats}
package. The 0 values of the estimated density
distributions of the spectra are replaced by a
regularization parameter (\code{"reg"} argument) for
numerical stability. Finally the divergence between the
computed spectral histogramas is computed using the SID
algorithm. Note that if \code{mode = "density"}, the
\code{sid} function will accept negative values and matrix
centering will be possible.
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
Xr <- Xr[!is.na(Yr),]

# Example 1
# Compute the SID distance between all the samples in Xr
xr.sid <- sid(Xr = Xr)
xr.sid

# Example 2
# Compute the SID distance between the samples in Xr and the samples in Xu
xru.sid <- sid(Xr = Xr, X2 = Xu)
xru.sid

# Example 3
# Compute the SID distance between the samples in Xr and the samples in Xu using the histograms
xru.sid.hist <- sid(Xr = Xr, X2 = Xu, mode = "feature")
xru.sid.hist
}
}
\author{
Leonardo Ramirez-Lopez
}
\references{
Chang, C.I. 2000. An information theoretic-based approach
to spectral variability, similarity and discriminability
for hyperspectral image analysis. IEEE Transactions on
Information Theory 46, 1927-1932.
}
\seealso{
\code{\link[stats]{density}}
}

