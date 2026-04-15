# Local fitting method constructors

These functions create configuration objects that specify how local
regression models are fitted within the
[`mbl`](https://l-ramirez-lopez.github.io/resemble/reference/mbl.md)
function.

## Usage

``` r
fit_pls(ncomp, method = c("pls", "mpls", "simpls"), 
        scale = FALSE, max_iter = 100L, tol = 1e-6)

fit_wapls(min_ncomp, max_ncomp, method = c("mpls", "pls", "simpls"),
          scale = FALSE, max_iter = 100L, tol = 1e-6)

fit_gpr(noise_variance = 0.001, center = TRUE, scale = TRUE)
```

## Arguments

- ncomp:

  an integer indicating the number of PLS components to use in local
  regressions when `fit_pls` is used.

- min_ncomp:

  an integer indicating the minimum number of PLS components to use in
  local regressions when `fit_wapls` is used. See details.

- max_ncomp:

  an integer indicating the maximum number of PLS components to use in
  local regressions when `fit_wapls` is used. See details.

- method:

  a character string indicating the PLS algorithm to use. Options are:

  - `'pls'`: standard PLS using covariance between X and Y for weight
    computation (NIPALS algorithm).

  - `'mpls'`: modified PLS using correlation between X and Y for weight
    computation (NIPALS algorithm). See Shenk and Westerhaus (1991).

  - `'simpls'`: SIMPLS algorithm (de Jong, 1993). Computationally faster
    as it avoids iterative X deflation. Parameters `max_iter` and `tol`
    are ignored when this method is used.

  Default is `'pls'` for `fit_pls` and `'mpls'` for `fit_wapls`.

- scale:

  logical indicating whether predictors must be scaled. Default is
  `FALSE` for PLS methods and `TRUE` for GPR.

- max_iter:

  an integer indicating the maximum number of iterations for convergence
  in the NIPALS algorithm. Only used when `method = 'pls'` or
  `method = 'mpls'`. Default is 100.

- tol:

  a numeric value indicating the convergence tolerance for calculating
  scores in the NIPALS algorithm. Only used when `method = 'pls'` or
  `method = 'mpls'`. Default is 1e-6.

- noise_variance:

  a numeric value indicating the variance of the noise for Gaussian
  process local regressions (`fit_gpr`). Default is 0.001.

- center:

  logical indicating whether predictors should be centered before
  fitting. Only used for `fit_gpr`. Default is `TRUE`.

## Value

An object of class `c("fit_<method>", "fit_method")` containing the
specified parameters. This object is passed to
[`mbl`](https://l-ramirez-lopez.github.io/resemble/reference/mbl.md) to
configure local model fitting.

## Details

These functions create configuration objects that are passed to
[`mbl`](https://l-ramirez-lopez.github.io/resemble/reference/mbl.md) to
specify how local regression models are fitted.

There are three fitting methods available:

### Partial least squares (`fit_pls`)

Uses orthogonal scores partial least squares regression. Three algorithm
variants are available:

- **Standard PLS** (`method = 'pls'`): Uses the NIPALS algorithm with
  covariance-based weights.

- **Modified PLS** (`method = 'mpls'`): Uses the NIPALS algorithm with
  correlation-based weights. Proposed by Shenk and Westerhaus (1991),
  this approach gives equal influence to all predictors regardless of
  their variance scale.

- **SIMPLS** (`method = 'simpls'`): Uses the SIMPLS algorithm (de Jong,
  1993), which deflates the cross-product matrix rather than X itself.
  This is computationally faster, especially for wide matrices, and
  produces identical predictions to standard PLS.

The only parameter to optimise is the number of PLS components
(`ncomp`).

### Weighted average PLS (`fit_wapls`)

This method was developed by Shenk et al. (1997) and is used as the
regression method in the LOCAL algorithm. It fits multiple PLS models
using different numbers of components (from `min_ncomp` to `max_ncomp`).
The final prediction is a weighted average of predictions from all
models, where the weight for component \\j\\ is:

\\w\_{j} = \frac{1}{s\_{1:j} \times g\_{j}}\\

where \\s\_{1:j}\\ is the root mean square of the spectral
reconstruction error of the target observation(s) when \\j\\ PLS
components are used, and \\g\_{j}\\ is the root mean square of the
squared regression coefficients for the \\j\\th component.

The same algorithm variants (`'pls'`, `'mpls'`, `'simpls'`) are
available. The default is `'mpls'` following the original LOCAL
implementation.

### Gaussian process regression (`fit_gpr`)

Gaussian process regression is a non-parametric Bayesian method
characterised by a mean and covariance function. This implementation
uses a dot product covariance.

The prediction vector \\A\\ is computed from training data (\\X\\,
\\Y\\) as:

\\A = (X X^{T} + \sigma^2 I)^{-1} Y\\

where \\\sigma^2\\ is the noise variance and \\I\\ is the identity
matrix. Prediction for a new observation \\x\_{u}\\ is:

\\\hat{y}\_{u} = x\_{u} X^{T} A\\

The only parameter is the noise variance (`noise_variance`).

## References

de Jong, S. (1993). SIMPLS: An alternative approach to partial least
squares regression. Chemometrics and Intelligent Laboratory Systems,
18(3), 251-263.

Rasmussen, C.E., Williams, C.K. (2006). Gaussian Processes for Machine
Learning. MIT Press.

Shenk, J.S., & Westerhaus, M.O. (1991). Populations structuring of near
infrared spectra and modified partial least squares regression. Crop
Science, 31(6), 1548-1555.

Shenk, J., Westerhaus, M., & Berzaghi, P. (1997). Investigation of a
LOCAL calibration procedure for near infrared instruments. Journal of
Near Infrared Spectroscopy, 5, 223-232.

Westerhaus, M. (2014). Eastern Analytical Symposium Award for
outstanding achievements in near infrared spectroscopy: my contributions
to near infrared spectroscopy. NIR news, 25(8), 16-20.

## See also

[`mbl`](https://l-ramirez-lopez.github.io/resemble/reference/mbl.md)

## Author

[Leonardo Ramirez-Lopez](https://orcid.org/0000-0002-5369-5120)

## Examples

``` r
# PLS with 10 components using standard algorithm
fit_pls(ncomp = 10)
#> Fitting method: pls
#>    ncomp    : 10 
#>    method   : pls 
#>    scale    : FALSE 
#>    max_iter : 100 
#>    tol      : 1e-06 

# PLS with modified algorithm (correlation-based weights)
fit_pls(ncomp = 10, method = "mpls")
#> Fitting method: pls
#>    ncomp    : 10 
#>    method   : mpls 
#>    scale    : FALSE 
#>    max_iter : 100 
#>    tol      : 1e-06 

# PLS with SIMPLS (faster, no iteration)
fit_pls(ncomp = 10, method = "simpls")
#> Fitting method: pls
#>    ncomp    : 10 
#>    method   : simpls 
#>    scale    : FALSE 
#>    max_iter : 100 
#>    tol      : 1e-06 

# Weighted average PLS (LOCAL-style)
fit_wapls(min_ncomp = 3, max_ncomp = 12)
#> Fitting method: wapls
#>    min_ncomp : 3 
#>    max_ncomp : 12 
#>    method    : mpls 
#>    scale     : FALSE 
#>    max_iter  : 100 
#>    tol       : 1e-06 

# Weighted average PLS with SIMPLS
fit_wapls(min_ncomp = 3, max_ncomp = 15, method = "simpls")
#> Fitting method: wapls
#>    min_ncomp : 3 
#>    max_ncomp : 15 
#>    method    : simpls 
#>    scale     : FALSE 
#>    max_iter  : 100 
#>    tol       : 1e-06 

# Gaussian process regression
fit_gpr()
#> Fitting method: gpr
#>    noise_variance : 0.001 
#>    center         : TRUE 
#>    scale          : TRUE 
fit_gpr(noise_variance = 0.01)
#> Fitting method: gpr
#>    noise_variance : 0.01 
#>    center         : TRUE 
#>    scale          : TRUE 
```
