# resemble: Regression and similarity evaluation for memory-based learning of spectral data

Click [here](http://l-ramirez-lopez.github.io/resemble_v0.1/) to visit the page where you can download the [`resemble`](http://l-ramirez-lopez.github.io/resemble_v0.1/) package!

[`resemble`](http://l-ramirez-lopez.github.io/resemble_v0.1/) implements a function dedicated to non-linear modelling of complex visible and infrared spectral data based on memory-based learning (MBL, _a.k.a_ instance-based learning or local modelling in the chemometrics literature). The package also includes functions for: computing and evaluate spectral similarity/dissimilarity matrices; projecting the spectra onto low dimensional orthogonal variables; removing irrelevant spectra from a reference set; etc. 


The functions for computing and evaluate spectral similarity/dissimilarity matrices can be summarized as follows:

| Function                 | Description/computes...                                                                                  |
| -----------------------  | -------------------------------------------------------------------------------------------------------  |
| `fDiss`                  | Euclidean and Mahalanobis distances as well as the cosine dissimilarity (_a.k.a_ spectral angle mapper)  |            
| `corDiss`                | correlation and moving window correlation dissimilarity                                                  |
| `sid`                    | spectral information divergence between spectra or between the probability distributions of spectra      |
| `orthoDiss`              | principal components and partial least squares dissimilarity (including several options)                 | 
| `simEval`                | evaluates a given similarity/dissimilarity matrix based on the concept of side information               |  

The functions for projecting the spectra onto low dimensional orthogonal variables are:

| Function                 | Description                                                                                                  |
| -----------------------  | ------------------------------------------------------------------------------------------------------------ |
| `pcProjection`           | projects the spectra onto a principal component space                                                        |                      
| `plsProjection`          | projects the spectra onto a partial least squares component space  (_a.k.a_ projection to latent structures) |                                      
| `orthoProjection`        | reproduces either the `pcProjection` or the `plsProjection` functions                                        |  

The projection functions also offer different options for optimizing/selecting the number of components involved in the projection.

The functions modelling the spectra using memory-based learning are:

| Function                 | Description                                              |
| -----------------------  | -------------------------------------------------------  |
| `mblController`          | controls some modelling aspects of the `mbl` function    |                     
| `mbl`                    | models the spectra by memory-based learning              |                                      

Some additional miscellaneous functions are:

| Function                 | Description                                                            |
| -----------------------  | ---------------------------------------------------------------------  |
| `print.mbl`              | prints a summary of the results obtained by the `mbl` function         |                     
| `plot.mbl`               | plots a summary of the results obtained by the `mbl` function          |       
| `print.localOrthoDiss`   | prints local distance matrices generated with the `orthoDiss` function |

In order to expand a little bit more the explanation on the `mbl` function, let's define first the basic input datasets:

* __Reference (training) set__: Dataset with *n* reference samples (e.g. spectral library) to be used in the calibration of a spectral models. Xr represents the matrix of samples (containing the spectral predictor variables) and Yr represents a given response variable corresponding to Xr.

* __Prediction set__ : Data set with _m_ samples where the response variable (Yu) is unknown. However it can be predicted by applying a spectral model (calibrated by using Xr and Yr) on the spectra of these samples (Xu). 

In order to predict each value in Yu, the `mbl` function takes each sample in Xu and searches in Xr for its _k_-nearest neighbours (most spectrally similar samples). Then a (local) model is calibrated with these (reference) neighbours and it immediately predicts the correspondent value in Yu from Xu. In the function, the _k_-nearest neighbour search is performed by computing spectral similarity/dissimilarity matrices between samples. The `mbl` function offers the following regression options for calibrating the (local) models:
                          
| Regression methods in the `mbl` function | Description                                                                             
| ---------------------------------------  | ----------------------------------------- | 
| `'gpr'`                                  | Gaussian process with linear kernel       | 
| `'pls'`                                  | Partial least squares                     | 
| `'wapls1'`                               | Weighted average partial least squares 1  | 
| `'wapls2'`                               | Weighted average partial least squares 2  | 

## Keywords
* _Infrared spectroscopy_
* _Chemometrics_
* _Local modelling_
* _Spectral library_
* _Lazy learning_
* _Soil spectroscopy_

## Bug report and development version

You can send an email to the package maintainer (<leonardo.ramirez@usys.ethz.ch>; <leonardo.ramirez@wsl.ch>) or create an [issue](https://github.com/l-ramirez-lopez/resemble_v0.1/issues) on github. 