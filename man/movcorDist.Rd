\name{movcorDist}
\alias{movcorDist}
\title{Moving/rolling correlation distance of two matrices}
\usage{
  movcorDist(X,Y,w)
}
\arguments{
  \item{X}{a \code{matrix}}

  \item{Y}{a \code{matrix}}

  \item{w}{window size (must be odd)}
}
\value{
  a \code{matrix} of correlation distance
}
\description{
  Computes a moving window correlation distance between two
  data matrices
}
\author{
  Leonardo Ramirez-Lopez and Antoine Stevens
}

