# Global spectral calibration model

Fits a global calibration model for spectral data using partial least
squares (PLS) or Gaussian process regression (GPR). Unlike
[`mbl`](https://l-ramirez-lopez.github.io/resemble/reference/mbl.md),
which builds local models for each prediction, `model()` fits a single
global model to the reference data.

## Usage

``` r
model(Xr, Yr, fit_method, control = model_control(), verbose = TRUE)

# S3 method for class 'resemble_model'
predict(object, newdata, ncomp = NULL, ...)
```

## Arguments

- Xr:

  A numeric matrix of predictor variables, with observations in rows and
  variables in columns. Typically spectral data.

- Yr:

  A numeric matrix with one column containing the response variable
  values corresponding to the observations in `Xr`.

- fit_method:

  An object created by
  [`fit_pls`](https://l-ramirez-lopez.github.io/resemble/reference/fit_methods.md)
  or
  [`fit_gpr`](https://l-ramirez-lopez.github.io/resemble/reference/fit_methods.md)
  specifying the regression method and its parameters, including
  centring and scaling options.
  [`fit_wapls()`](https://l-ramirez-lopez.github.io/resemble/reference/fit_methods.md)
  is not supported for global models because it requires target
  observations to compute weights.

- control:

  A list created by
  [`model_control()`](https://l-ramirez-lopez.github.io/resemble/reference/model_control.md)
  specifying the cross-validation settings. The default is
  [`model_control()`](https://l-ramirez-lopez.github.io/resemble/reference/model_control.md).

- verbose:

  Logical indicating whether progress information should be printed.
  Default is `TRUE`.

- object:

  A `resemble_model` object returned by `model`.

- newdata:

  A numeric matrix of new observations with the same number of columns
  as the training data.

- ncomp:

  For PLS models, the number of components to use for prediction. The
  default is the number used during fitting. Ignored for GPR models. For
  prediction, this can be a vector of integers representing the
  predictions coming from models with the requested components.

- ...:

  Additional arguments, currently unused.

## Value

For `model()`, an object of class `resemble_model` containing:

- `fit_method`: Fitting method object used.

- `control`: Model control object used.

- `model`: Fitted model object.

- `cv_results`: Cross-validation results, if `validation_type = "lgo"`.

- `n_obs`: Number of observations used.

- `n_vars`: Number of predictor variables.

For `predict.resemble_model()`, a numeric matrix of predictions. For PLS
models, columns correspond to different numbers of components.

## Details

`model()` provides a straightforward interface for fitting global
calibration models, in contrast to the local modelling approach
implemented in
[`mbl`](https://l-ramirez-lopez.github.io/resemble/reference/mbl.md).
This is useful when:

- the relationship between spectra and the response is approximately
  linear across the full dataset

- a single portable model is needed for deployment

- computational efficiency is prioritised over predictive performance in
  heterogeneous datasets

### Cross-validation

When `validation_type = "lgo"`, stratified random sampling is used to
create training-validation splits that preserve the distribution of the
response variable. Cross-validation results include RMSE, R², and
standardised RMSE for each component in PLS models, or overall in GPR
models.

### PLS models

When
[`fit_pls()`](https://l-ramirez-lopez.github.io/resemble/reference/fit_methods.md)
is used, the function fits a PLS model with the specified number of
components. If cross-validation is enabled, the optimal number of
components is selected based on the minimum RMSE.

### GPR models

When
[`fit_gpr()`](https://l-ramirez-lopez.github.io/resemble/reference/fit_methods.md)
is used, the function fits a Gaussian process regression model with a
dot-product covariance function. The noise variance parameter controls
regularisation.

## References

de Jong, S. (1993). SIMPLS: An alternative approach to partial least
squares regression. *Chemometrics and Intelligent Laboratory Systems*,
18(3), 251–263.

Rasmussen, C. E., & Williams, C. K. I. (2006). *Gaussian Processes for
Machine Learning*. MIT Press.

Shenk, J. S., & Westerhaus, M. O. (1991). Population structuring of near
infrared spectra and modified partial least squares regression. *Crop
Science*, 31(6), 1548–1555.

## See also

[`fit_pls`](https://l-ramirez-lopez.github.io/resemble/reference/fit_methods.md)
and
[`fit_gpr`](https://l-ramirez-lopez.github.io/resemble/reference/fit_methods.md)
for specifying the fitting method;
[`model_control`](https://l-ramirez-lopez.github.io/resemble/reference/model_control.md)
for controlling aspects of the modelling process;
[`mbl`](https://l-ramirez-lopez.github.io/resemble/reference/mbl.md) for
memory-based (local) learning.

## Author

[Leonardo Ramirez-Lopez](https://orcid.org/0000-0002-5369-5120)

## Examples

``` r
if (FALSE) { # \dontrun{
library(prospectr)
data(NIRsoil)

# Preprocess spectra
Xr <- savitzkyGolay(NIRsoil$spc, m = 1, p = 2, w = 11)
Yr <- NIRsoil$CEC

# Remove missing values
ok <- !is.na(Yr)
Xr <- Xr[ok, ]
Yr <- as.matrix(Yr[ok])

# Fit a PLS model with 10 components and cross-validation
# Scaling is controlled via fit_pls()
pls_mod <- model(
  Xr = Xr,
  Yr = Yr,
  fit_method = fit_pls(ncomp = 10, scale = FALSE),
  control = model_control(validation_type = "lgo", number = 10)
)

# View cross-validation results
pls_mod$cv_results

# Fit a GPR model (centring/scaling controlled via fit_gpr())
gpr_mod <- model(
  Xr = Xr,
  Yr = Yr,
  fit_method = fit_gpr(noise_variance = 0.001, scale = TRUE),
  control = model_control(validation_type = "lgo")
)
} # }
```
