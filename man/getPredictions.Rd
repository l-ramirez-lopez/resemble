\name{getPredictions}
\alias{getPredictions}
\title{Extract predictions from an object of class \code{mbl}}
\usage{
  getPredictions(object)
}
\arguments{
  \item{object}{an object of class \code{mbl} as returned
  by \code{mbl}}
}
\value{
  a \code{data.frame} of predicted values according to
  either \code{k} or \code{k.dist}
}
\description{
  Extract predictions from an object of class \code{mbl}
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

ctrl <- mblController(sm = "pls", ws = 51, pcSelection = list("opc", 40), valMethod = c("NNv"), scaled = TRUE, center = TRUE)

ex1 <- mbl(Yr = Yr, Xr = Xr, Yu = NULL, Xu = Xu,
           mblCtrl = ctrl,
           distUsage = "predictors",
           k = seq(30, 150, 15),
           method = "wapls1",
           plsF = c(7, 20))

getPredictions(ex1)
}
}
\author{
  Leonardo Ramirez-Lopez and Antoine Stevens
}
\seealso{
  \code{\link{mbl}}
}

