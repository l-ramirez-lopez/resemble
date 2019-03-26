[![Travis-CI Build Status](https://travis-ci.org/l-ramirez-lopez/resemble.svg?branch=master)](https://travis-ci.org/l-ramirez-lopez/resemble/)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/resemble)](http://cran.r-project.org/web/packages/resemble)
[![Total_Downloads](http://cranlogs.r-pkg.org/badges/grand-total/resemble)](http://cran.r-project.org/web/packages/resemble)
![alt text](https://raw.githubusercontent.com/l-ramirez-lopez/resemble/master/resemble_logo.png)

# Regression and Similarity Evaluation for Memory-Based Learning in Spectral Chemometrics
_Leo Ramirez-Lopez & Antoine Stevens_

_Last update: 17.05.2019 :::: 10:09 GMT+2_

Visit the [`resemble` site here](http://l-ramirez-lopez.github.io/resemble/)

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
__Note__: Apart from these packages we stronly recommend to download and install Rtools ([directly from here](http://cran.r-project.org/web/packages/devtools/index.html) or from CRAN [https://cran.r-project.org/bin/windows/Rtools/](https://cran.r-project.org/bin/windows/Rtools/)). 
This is important for obtaining the proper C++ toolchain that you might need for using `resemble`.

Then, install `resemble`

For the development version, download the [binary (.zip) file from here](https://github.com/l-ramirez-lopez/resemble/archive/1.2.2.zip) or the [source file (.tar.gz) from here](https://github.com/l-ramirez-lopez/resemble/archive/1.2.2.tar.gz). Remember you should have [R>=3.2.2](http://cran.r-project.org/). Suppose you downloaded the binary file to 'C:/MyFolder/', then you should be able to install the package as follows:
  
  ```
install.packages('C:/MyFolder/resemble-1.2.2.zip', repos = NULL)
````
or

```
install.packages('C:/MyFolder/resemble-1.2.2.tar.gz', type = 'source', repos = NULL)
```

You can also install the `resemble` package directly from github using [`devtools`](http://cran.r-project.org/web/packages/devtools/index.html) (along with a proper installed version of [Rtools](http://cran.r-project.org/bin/windows/Rtools/)):
  
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

[`resemble`](http://l-ramirez-lopez.github.io/resemble/) implements a function dedicated to non-linear modelling of complex visible and infrared spectral data based on memory-based learning (MBL, _a.k.a_ instance-based learning or local modelling in the chemometrics literature). The package also includes functions for: computing and evaluate spectral similarity/dissimilarity matrices; projecting the spectra onto low dimensional orthogonal variables; removing irrelevant spectra from a reference set; etc. 


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

## News
* 2019-03: Another paper using `resemble`... I just published a [scientific paper](https://onlinelibrary.wiley.com/doi/10.1111/ejss.12752) were we used memory-based learning (MBL) for digital soil mapping. Here we use MBL to remove local calibration outliers rather than using this approach to overcome the typical complexity of large spectral datasets. (Ramirez‐Lopez, L., Wadoux, A. C., Franceschini, M. H. D., Terra, F. S., Marques, K. P. P., Sayão, V. M., & Demattê, J. A. M. (2019). Robust soil mapping at the farm scale with vis–NIR spectroscopy. European Journal of Soil Science. 70, 378–393).
* 2019-03:[Jaconi et al. (2019)](https://www.sciencedirect.com/science/article/pii/S0016706118311340?via%3Dihub) implemented a memory-based learning algorithm (using `resemble`) to conduct accurate NIR predictions of soil texture at National scale in Germany. (Jaconi, A., Vos, C. and Don, A., 2019. Near infrared spectroscopy as an easy and precise method to estimate soil texture. Geoderma, 337, pp.906-913).
*2018-12: [Chen, et al. (2018)](https://www.sciencedirect.com/science/article/pii/S0341816218304867) implemented a memory-based learning algorithm (using `resemble`) to improve the accuracy of NIR predictions of soil organic matter in China. (Hong, Y., Chen, S., Liu, Y., Zhang, Y., Yu, L., Chen, Y., Liu, Y., Cheng, H. and Liu, Y., 2019. Combination of fractional order derivative and memory-based learning algorithm to improve the estimation accuracy of soil organic matter by visible and near-infrared spectroscopy. Catena, 174, pp.104-116).
* 2018-11: In [this recent scientific paper](https://www.frontiersin.org/articles/10.3389/fpls.2018.01642/full) the authors used `resemble` to predict the chemoical composition of Common Beans in Spain. (Rivera, A., Plans, M., Sabaté, J., Casañas, F., Casals, J., Rull, A., & Simó, J. (2018). The Spanish core collection of common beans (Phaseolus vulgaris L.): an important source of variability for breeding chemical composition. Frontiers in Plant Science, 9).
* 2018-07: Another use-case of `resemble` is presented [Gholizadeh et al.(2018)](https://www.mdpi.com/2072-4292/10/8/1172/htm). (Gholizadeh, A., Saberioon, M., Carmon, N., Boruvka, L. and Ben-Dor, E., 2018. Examining the Performance of PARACUDA-II Data-Mining Engine versus Selected Techniques to Model Soil Carbon from Reflectance Spectra. Remote Sensing, 10(8), p.1172.).
* 2018-01: [Dotto, et al. (2018)](http://ufgrunwald.com/wp-content/uploads/2018/02/Dotto-et-al-2018.-Vis-NIR-Geoderma.pdf) have implemented memory-based learning with `resemble` to accurately predict soil organic Carbon at a region in Brazil. (Dotto, A. C., Dalmolin, R. S. D., ten Caten, A., & Grunwald, S. (2018). A systematic study on the application of scatter-corrective and spectral-derivative preprocessing for multivariate prediction of soil organic carbon by Vis-NIR spectra. Geoderma, 314, 262-274).
* 2017-11: [Here](https://books.google.ch/books?hl=en&lr=&id=xlKODgAAQBAJ&oi=fnd&pg=PA129&dq=Evaluation+and+comparison+of+different+approaches+to+multi-product+brix+calibration+in+near-infrared+spectroscopy&ots=TGWzsRpP7X&sig=3u-pfuwkMZoxmhbfcx40eLiMOsY#v=onepage&q=Evaluation%20and%20comparison%20of%20different%20approaches%20to%20multi-product%20brix%20calibration%20in%20near-infrared%20spectroscopy&f=false) the authors predicted brix values in differet food products using memory-based learning implemented with `resemble`. (Kopf, M., Gruna, R., Längle, T. and Beyerer, J., 2017, March. Evaluation and comparison of different approaches to multi-product brix calibration in near-infrared spectroscopy. In OCM 2017-Optical Characterization of Materials-conference proceedings (p. 129). KIT Scientific Publishing).
* 2016-05: In [this recent scientific paper](http://www.sciencedirect.com/science/article/pii/S001670611630180X) the authors sucesfully used `resemble` to predict soil organic carbon content at national scale in France. (Clairotte, M., Grinand, C., Kouakoua, E., Thébault, A., Saby, N. P., Bernoux, M., & Barthès, B. G. (2016). National calibration of soil organic carbon concentration using diffuse infrared reflectance spectroscopy. Geoderma, 276, 41-52)
* 2016-04: [This paper](http://www.mdpi.com/2072-4292/8/4/341) shows some interesting results on applying memory-based learning to predict soil properties.
* 2016-04: In some recent entries of [this blog](http://nir-quimiometria.blogspot.com/), the author shows some exmaples on the use `resemble`
* 2016-02: As promised, `resemble 1.2 (alma-de-coco)` is now available on [CRAN](https://cran.r-project.org/web/packages/resemble/index.html).
* 2016-01: The version 1.2 (alma-de-coco) has been submitted to CRAN and is available from the github repository!
* 2015-11: A pre-release of the version 1.2.0 (1.2.0.9000 alma-de-coco) is now available! `resemble` is now faster! Some critical functions (e.g. pls and gaussian process regressions were re-written in C++ using `Rcpp`. This time the new version will be available at CRAN very soon!.
* 2015-11 Well, the version 1.1.3 was never released on CRAN since we decided to carry out major improvements in terms of computational performance. 
* 2014-10: A pre-release of the version 1.1.3 of the package is already available at this website. We hope it will be available at CRAN very soon!
* 2014-06: Check  [this video](https://www.youtube.com/watch?v=7sCIEeNehgE&feature=youtu.be) where a renowned NIR scientist talks about local calibrations.
<iframe width="640" height="360" src="https://www.youtube.com/embed/7sCIEeNehgE" frameborder="0" allowfullscreen></iframe>
* 2014-04: A short note on the resemble and prospectr packages was published in [this newsletter](http://www.pedometrics.org/Pedometron/Pedometron34.pdf). There we provide some examples on representative subset selection and on how to reproduce the LOCAL and spectrum-based learner algorithms. In those examples the dataset of the Chemometric challenge of 'Chimiométrie 2006' (included in the prospectr package) is used.
* 2014-03: The package was released on CRAN!

## Other R'elated stuff
* [Check our other project called `prospectr`.](http://antoinestevens.github.io/prospectr/)
* [Check this presentation in which we used the resemble package to predict soil attributes from large scale soil spectral libraries.](http://www.fao.org/fileadmin/user_upload/GSP/docs/Spectroscopy_dec13/SSW2013_f.pdf)

## Bug report and development version

You can send an e-mail to the package maintainer (<ramirez.lopez.leo@gmail.com>) or create an [issue](https://github.com/l-ramirez-lopez/resemble/issues) on github. 
