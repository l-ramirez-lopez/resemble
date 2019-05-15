---
output:
  pdf_document: default
  html_document: default
---
[![Travis-CI Build Status](https://travis-ci.org/l-ramirez-lopez/resemble.svg?branch=master)](https://travis-ci.org/l-ramirez-lopez/resemble/)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/resemble)](https://CRAN.R-project.org/package=resemble)
[![Total_Downloads](https://cranlogs.r-pkg.org/badges/grand-total/resemble)](https://CRAN.R-project.org/package=resemble)
![](https://raw.githubusercontent.com/l-ramirez-lopez/resemble/master/resemble_logo.png)

# Regression and Similarity Evaluation for Memory-Based Learning in Spectral Chemometrics
_Leo Ramirez-Lopez & Antoine Stevens_

_Last update: 08.03.2018_

Visit the [`resemble` site here](https://l-ramirez-lopez.github.io/resemble/)

Installing the package is very simple:
  ```
install.packages('resemble')
```
If you do not have the following packages installed, in some cases it is better to install them first
```
install.packages('Rcpp')
install.packages('RcppArmadillo')
install.packages('foreach')
install.packages('iterators')
```
__Note__: Apart from these packages we stronly recommend to download and install Rtools ([directly from here](https://cran.r-project.org/bin/windows/Rtools). 
This is important for obtaining the proper C++ toolchain that you might need for using `resemble`.

Then, install `resemble`

For the development version, download the [binary (.zip) file from here](https://github.com/l-ramirez-lopez/resemble/archive/1.2.2.zip) or the [source file (.tar.gz) from here](https://github.com/l-ramirez-lopez/resemble/archive/1.2.2.tar.gz). Remember you should have [R>=3.2.2](https://cran.r-project.org/). Suppose you downloaded the binary file to 'C:/MyFolder/', then you should be able to install the package as follows:
  
  ```
install.packages('C:/MyFolder/resemble-1.2.2.zip', repos = NULL)
````
or

```
install.packages('C:/MyFolder/resemble-1.2.2.tar.gz', type = 'source', repos = NULL)
```

You can also install the `resemble` package directly from github using [`devtools`](https://CRAN.R-project.org/package=devtools) (along with a proper installed version of [Rtools](https://cran.r-project.org/bin/windows/Rtools/)):
  
  ```
require("devtools")
install_github("resemble","l-ramirez-lopez")
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
ctrl <- mblControl(sm = 'pc', pcSelection = list('opc', 40),
                   valMethod = 'NNv', center = TRUE)

sbl.u <- mbl(Yr = Yr, Xr = Xr, Yu = NULL, Xu = Xu,
             mblCtrl = ctrl,
             dissUsage = 'predictors',
             k = seq(40, 150, by = 10),
             method = 'gpr')

getPredictions(sbl.u)
````

[`resemble`](https://l-ramirez-lopez.github.io/resemble/) implements a function dedicated to non-linear modelling of complex visible and infrared spectral data based on memory-based learning (MBL, _a.k.a_ instance-based learning or local modelling in the chemometrics literature). The package also includes functions for: computing and evaluate spectral similarity/dissimilarity matrices; projecting the spectra onto low dimensional orthogonal variables; removing irrelevant spectra from a reference set; etc. 


The functions for computing and evaluate spectral similarity/dissimilarity matrices can be summarized as follows:
  
  __`fDiss`__:                  Euclidean and Mahalanobis distances as well as the cosine dissimilarity (_a.k.a_ spectral angle mapper)              
__`corDiss`__:                correlation and moving window correlation dissimilarity                                                 
__`sid`__:                    spectral information divergence between spectra or between the probability distributions of spectra      
__`orthoDiss`__:              principal components and partial least squares dissimilarity (including several options)                  
__`simEval`__:                evaluates a given similarity/dissimilarity matrix based on the concept of side information                 

The functions for projecting the spectra onto low dimensional orthogonal variables are:
  
  __`pcProjection`__:            projects the spectra onto a principal component space                                                                              
__`plsProjection`__:           projects the spectra onto a partial least squares component space  (_a.k.a_ projection to latent structures)                                       
__`orthoProjection`__:         reproduces either the `pcProjection` or the `plsProjection` functions                                          

The projection functions also offer different options for optimizing/selecting the number of components involved in the projection.

The functions modelling the spectra using memory-based learning are:
  
  __`mblControl`__:              controls some modelling aspects of the `mbl` function                         
__`mbl`__:                     models the spectra by memory-based learning                                                    

Some additional miscellaneous functions are:
  
  __`print.mbl`__:               prints a summary of the results obtained by the `mbl` function                              
__`plot.mbl`__:                plots a summary of the results obtained by the `mbl` function                 
__`print.localOrthoDiss`__:    prints local distance matrices generated with the `orthoDiss` function 

In order to expand a little bit more the explanation on the `mbl` function, let's define first the basic input datasets:

* __Reference (training) set__: Dataset with *n* reference samples (e.g. spectral library) to be used in the calibration of spectral models. Xr represents the matrix of samples (containing the spectral predictor variables) and Yr represents a given response variable corresponding to Xr.

* __Prediction set__ : Data set with _m_ samples where the response variable (Yu) is unknown. However it can be predicted by applying a spectral model (calibrated by using Xr and Yr) on the spectra of these samples (Xu). 

In order to predict each value in Yu, the `mbl` function takes each sample in Xu and searches in Xr for its _k_-nearest neighbours (most spectrally similar samples). Then a (local) model is calibrated with these (reference) neighbours and it immediately predicts the correspondent value in Yu from Xu. In the function, the _k_-nearest neighbour search is performed by computing spectral similarity/dissimilarity matrices between samples. The `mbl` function offers the following regression options for calibrating the (local) models:

__`'gpr'`__:                                   Gaussian process with linear kernel        
__`'pls'`__:                                   Partial least squares                      
__`'wapls1'`__:                                Weighted average partial least squares 1   
__`'wapls2'`__:                                Weighted average partial least squares 2 (no longer supported)

## Keywords
* _Infrared spectroscopy_
* _Chemometrics_
* _Local modelling_
* _Spectral library_
* _Lazy learning_
* _Soil spectroscopy_

## `resemble` in science
* 2018-03: In [this scientific paper, Dotto et al. (2018),] (https://www.sciencedirect.com/science/article/pii/S0016706117308571) sucesfully used `resemble` to predict soil organic carbon content at a regional scale in Brazil. 
* 2017-03: [Kopf et al. (2017),] (http://publica.fraunhofer.de/eprints/urn_nbn_de_0011-n-4390395.pdf) sucessfully used `resemble` for developing multi-product brix NIR predictions.
* 2017-11: In [this scientific paper, Jaconi et al. (2017),] (http://onlinelibrary.wiley.com/doi/10.1111/ejss.12485/full) sucessfully used `resemble` to predict soil organic carbon content at the country scale in Germany. 
* 2017-11: In [this scientific paper, Lobsey et al. (2017),] (https://onlinelibrary.wiley.com/doi/10.1111/ejss.12490/full) used `resemble` to combine a local regression method with a spiking technic to predict soil organic carbon content in two farms in Australia using a national spectral library. 
* 2016-08: In [this scientific paper, Clairotte, et al (2016)] (https://www.sciencedirect.com/science/article/pii/S001670611630180X) sucesfully used `resemble` to predict soil organic carbon content at national scale in France. 
* 2016-04: In [this scientific paper, Gholizadeh et al. (2016),] (https://www.mdpi.com/2072-4292/8/4/341) used `resemble` to predict soil texture at regional scale in Czech Republic.
* 2016-04: In [this blog] (https://nir-quimiometria.blogspot.com/), the author shows some exmaples on the use `resemble`
* 2014-04: A short note on the resemble and prospectr packages was published in [this newsletter](http://www.pedometrics.org/Pedometron/Pedometron34.pdf). There we provide some examples on representative subset selection and on how to reproduce the LOCAL and spectrum-based learner algorithms. In those examples the dataset of the Chemometric challenge of 'Chimiométrie 2006' (included in the prospectr package) is used.



## News
*2017-06 A crash course entitled "Using R for infrared spectroscopy research" was delivered during the pedometrics conference 2018 in Netherlands. The course included some lessos on `resemble`. The course  materials can be found [here](https://github.com/AlexandreWadoux/pedometrics2017_workshop)
* 2016-04: [This paper] (https://www.mdpi.com/2072-4292/8/4/341) shows some interesting results on appliying memory-based learning to predict soil properties.
* 2016-02: `resemble 1.2 (alma-de-coco)` is available on [CRAN](https://CRAN.R-project.org/package=resemble).
* 2016-01: The version 1.2 (alma-de-coco) was submitted to CRAN
* 2015-11: A pre-release of the version 1.2.0 (1.2.0.9000 alma-de-coco) is now available! `resemble` is now faster! Some critical functions (e.g. pls and gaussian process regressions were re-written in C++ using `Rcpp`. This time the new version will be available at CRAN very soon!.
* 2015-11 Well, the version 1.1.3 was never released on CRAN since we decided to carry out major improvements in terms of computational performance. 
* 2014-10: A pre-release of the version 1.1.3 of the package is already available at this website. We hope it will be available at CRAN very soon!
* 2014-06: Check  [this video](https://www.youtube.com/watch?v=7sCIEeNehgE&feature=youtu.be) where a renowned NIR scientist talks about local calibrations.
<iframe width="640" height="360" src="https://www.youtube.com/embed/7sCIEeNehgE" frameborder="0" allowfullscreen></iframe>
* 2014-03: The package was released on CRAN!

## Other R'elated stuff
* [Check our other project called `prospectr`.](https://antoinestevens.github.io/prospectr/)
* [Check this presentation in which we used the resemble package to predict soil attributes from large scale soil spectral libraries.](http://www.fao.org/fileadmin/user_upload/GSP/docs/Spectroscopy_dec13/SSW2013_f.pdf)

## Bug report and development version

You can send an e-mail to the package maintainer (<ramirez.lopez.leo@gmail.com>) or create an [issue](https://github.com/l-ramirez-lopez/resemble/issues) on github. 
