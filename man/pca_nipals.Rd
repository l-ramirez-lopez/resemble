% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{pca_nipals}
\alias{pca_nipals}
\title{Principal components based on  the non-linear iterative partial least squares (nipals) algorithm}
\usage{
pca_nipals(X, ncomp, center, scale,
           maxiter, tol,
           pcSelmethod = "var",
           pcSelvalue = 0.01)
}
\arguments{
\item{X}{a matrix of predictor variables.}

\item{ncomp}{the number of pls components.}

\item{scale}{logical indicating whether \code{X} must be scaled.}

\item{maxiter}{maximum number of iterations.}

\item{tol}{limit for convergence of the algorithm in the nipals algorithm.}

\item{pcSelmethod}{the method for selecting the number of components.
Options are: \code{'cumvar'} (for selecting the number of principal components based on a given
cumulative amount of explained variance) and \code{"var"} (for selecting the number of principal
components based on a given amount of explained variance). Default is \code{'var'}}

\item{pcSelvalue}{a numerical value that complements the selected method (\code{pcSelmethod}).
If \code{"cumvar"} is chosen, it must be a value (larger than 0 and below 1) indicating the maximum
amount of cumulative variance that the retained components should explain. If \code{"var"} is chosen,
it must be a value (larger than 0 and below 1) indicating that components that explain (individually)
a variance lower than this threshold must be excluded. If \code{"manual"} is chosen, it must be a value
specifying the desired number of principal components to retain. Default is 0.01.}

\item{Y}{a matrix of either a single or multiple response variables.}
}
\value{
a list containing the following elements:
\itemize{
\item{\code{pc_scores}}{ a matrix of principal component scores.}
\item{\code{pc_loadings}}{ a matrix of of principal component loadings.}
\item{\code{variance}}{ a matrix of the variance of the principal components.}
\item{\code{scale}}{ a \code{list} conating two objects: \code{center} and \code{scale}, which correspond to the vectors used to center and scale the input matrix.}
}
}
\description{
Computes orthogonal socres partial least squares (opls) regressions with the NIPALS algorithm. It allows multiple response variables.
For internal use only!
}
\author{
Leonardo Ramirez-Lopez
}
\keyword{internal}
