\name{which_minV}
\alias{which_minV}
\title{A function to compute indices of minimum values of a distance vector}
\usage{
which_minV(X,cores)
}
\arguments{
  \item{X}{a \code{vector} of distance (as computed in
  \code{resemble:::fastDistVV} or \code{base::dist})}

  \item{cores}{number of cores used to run the computation}
}
\value{
a \code{vector} of the indices of the nearest neighbours
}
\description{
A function to compute indices of minimum values of a
distance vector
}
\details{
Used internally to find the nearest neighbours
}
\author{
Antoine Stevens
}
\keyword{internal}

