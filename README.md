# resemble: Regression and similarity evaluation for memory-based learning in spectral chemometrics
_Leonardo Ramirez-Lopez & Antoine Stevens_

Visit the [`resemble` site here](http://l-ramirez-lopez.github.io/resemble/)

You can install the `resemble` package directly from github using [`devtools`](http://cran.r-project.org/web/packages/devtools/index.html) (with a proper installed version of [Rtools](http://cran.r-project.org/bin/windows/Rtools/)):

```
require("devtools")
install_github("resemble","l-ramirez-lopez")
```

You can also download the [binary (.zip) file from here](https://github.com/l-ramirez-lopez/resemble/blob/master/Installers/resemble_1.0.zip?raw=true) or the [source file (.tar.gz) from here](https://github.com/l-ramirez-lopez/resemble/blob/master/Installers/resemble_1.0.tar.gz?raw=true). Remeber you should have R>3.0.0. Supose you downloaded the binary file to 'C:/MyFolder/', then you should be able to install the package as follows:

First, if you do not have the following packages installed, you should install them first
```
install.packages('Rcpp')
install.packages('RcppArmadillo')
install.packages('pls')
install.packages('foreach')
install.packages('iterators')
```
Then, install `resemble`

```
install.packages('C:/MyFolder/resemble_1.0.zip')
````
or

```
install.packages('C:/MyFolder/resemble_1.0.tar.gz', type = 'source')
```

After installing `resemble` you should be also able to run the following lines:

```
require(resemble)

help(mbl)

#install.packages('prospectr')
require(prospectr)

data(NIRsoil)

Xu <- NIRsoil$spc[!as.logical(NIRsoil$train),]
Yu <- NIRsoil$CEC[!as.logical(NIRsoil$train)]
Yr <- NIRsoil$CEC[as.logical(NIRsoil$train)]
Xr <- NIRsoil$spc[as.logical(NIRsoil$train),]

Xu <- Xu[!is.na(Yu),]
Xr <- Xr[!is.na(Yr),]

Yu <- Yu[!is.na(Yu)]
Yr <- Yr[!is.na(Yr)]

# Example of the mbl function
# A mbl approach (the spectrum-based learner) as implemented in Ramirez-Lopez et al. (2013)
# An exmaple where Yu is supposed to be unknown, but the Xu (spectral variables) are known
ctrl1 <- mblController(sm = 'pc', pcSelection = list('opc', 40),
                       valMethod = 'NNv',
                       scaled = TRUE, center = TRUE)

sbl.u <- mbl(Yr = Yr, Xr = Xr, Yu = NULL, Xu = Xu,
             mblCtrl = ctrl1,
             dissUsage = 'predictors',
             k = seq(40, 150, by = 10),
             method = 'gpr')

getPredictions(sbl.u)
````

[`resemble`](http://l-ramirez-lopez.github.io/resemble/) implements a function dedicated to non-linear modelling of complex visible and infrared spectral data based on memory-based learning (MBL, _a.k.a_ instance-based learning or local modelling in the chemometrics literature). The package also includes functions for: computing and evaluate spectral similarity/dissimilarity matrices; projecting the spectra onto low dimensional orthogonal variables; removing irrelevant spectra from a reference set; etc. 


The functions for computing and evaluate spectral similarity/dissimilarity matrices can be summarized as follows:

`fDiss`:                  Euclidean and Mahalanobis distances as well as the cosine dissimilarity (_a.k.a_ spectral angle mapper)              
`corDiss`:                correlation and moving window correlation dissimilarity                                                 
`sid`:                    spectral information divergence between spectra or between the probability distributions of spectra      
`orthoDiss`:              principal components and partial least squares dissimilarity (including several options)                  
`simEval`:                evaluates a given similarity/dissimilarity matrix based on the concept of side information                 

The functions for projecting the spectra onto low dimensional orthogonal variables are:

`pcProjection`:            projects the spectra onto a principal component space                                                                              
`plsProjection`:           projects the spectra onto a partial least squares component space  (_a.k.a_ projection to latent structures)                                       
`orthoProjection`:         reproduces either the `pcProjection` or the `plsProjection` functions                                          

The projection functions also offer different options for optimizing/selecting the number of components involved in the projection.

The functions modelling the spectra using memory-based learning are:

`mblControl`:              controls some modelling aspects of the `mbl` function                         
`mbl`:                     models the spectra by memory-based learning                                                    

Some additional miscellaneous functions are:

`print.mbl`:               prints a summary of the results obtained by the `mbl` function                              
`plot.mbl`:                plots a summary of the results obtained by the `mbl` function                 
`print.localOrthoDiss`:    prints local distance matrices generated with the `orthoDiss` function 

In order to expand a little bit more the explanation on the `mbl` function, let's define first the basic input datasets:

* __Reference (training) set__: Dataset with *n* reference samples (e.g. spectral library) to be used in the calibration of spectral models. Xr represents the matrix of samples (containing the spectral predictor variables) and Yr represents a given response variable corresponding to Xr.

* __Prediction set__ : Data set with _m_ samples where the response variable (Yu) is unknown. However it can be predicted by applying a spectral model (calibrated by using Xr and Yr) on the spectra of these samples (Xu). 

In order to predict each value in Yu, the `mbl` function takes each sample in Xu and searches in Xr for its _k_-nearest neighbours (most spectrally similar samples). Then a (local) model is calibrated with these (reference) neighbours and it immediately predicts the correspondent value in Yu from Xu. In the function, the _k_-nearest neighbour search is performed by computing spectral similarity/dissimilarity matrices between samples. The `mbl` function offers the following regression options for calibrating the (local) models:
                          
`'gpr'`:                                   Gaussian process with linear kernel        
`'pls'`:                                   Partial least squares                      
`'wapls1'`:                                Weighted average partial least squares 1   
`'wapls2'`:                                Weighted average partial least squares 2   

## Keywords
* _Infrared spectroscopy_
* _Chemometrics_
* _Local modelling_
* _Spectral library_
* _Lazy learning_
* _Soil spectroscopy_

## News

2014-03: The package was released on CRAN!

## Other R'elated stuff
* Check our other project called [`prospectr`](http://antoinestevens.github.io/prospectr/)

## Bug report and development version

You can send an email to the package maintainer (<leonardo.ramirez@usys.ethz.ch>; <leonardo.ramirez@wsl.ch>) or create an [issue](https://github.com/l-ramirez-lopez/resemble/issues) on github. 
