# Evaluate dissimilarity matrices

**\[stable\]**

Evaluates a dissimilarity matrix by comparing each observation to its
nearest neighbor based on side information. For continuous variables,
RMSD and correlation are computed; for categorical variables, the kappa
index is used.

## Usage

``` r
diss_evaluate(diss, side_info)

sim_eval(d, side_info)
```

## Arguments

- diss:

  A symmetric dissimilarity matrix. Alternatively, a vector containing
  the lower triangle values (as returned by
  [`dist`](https://rdrr.io/r/stats/dist.html)).

- side_info:

  A matrix of side information corresponding to the observations. Can be
  numeric (one or more columns) or character (single column for
  categorical data).

- d:

  Deprecated. Use `diss` in `diss_evaluate()` instead.

## Value

A list with the following components:

- eval:

  For numeric side information: a matrix with columns `rmsd` and `r`
  (correlation). For categorical: a matrix with column `kappa`.

- global_eval:

  If multiple numeric side information variables are provided, summary
  statistics across variables.

- first_nn:

  A matrix with the original side information and the side information
  of each observation's nearest neighbor.

## Details

This function assesses whether a dissimilarity matrix captures
meaningful structure by examining the side information of nearest
neighbor pairs (Ramirez-Lopez et al., 2013). If observations that are
similar in the dissimilarity space also have similar side information
values, the dissimilarity is considered effective.

For numeric `side_info`, the root mean square of differences (RMSD)
between each observation and its nearest neighbor is computed:

\\j(i) = NN(x_i, X^{{-i}})\\ \\RMSD = \sqrt{\frac{1}{m} \sum\_{i=1}^{m}
(y_i - y\_{j(i)})^2}\\

where \\NN(x_i, X^{-i})\\ returns the index of the nearest neighbor of
observation \\i\\ (excluding itself), \\y_i\\ is the side information
value for observation \\i\\, and \\m\\ is the number of observations.

For categorical `side_info`, the kappa index is computed:

\\\kappa = \frac{p_o - p_e}{1 - p_e}\\

where \\p_o\\ is the observed agreement and \\p_e\\ is the agreement
expected by chance.

## References

Ramirez-Lopez, L., Behrens, T., Schmidt, K., Stevens, A., Dematte,
J.A.M., Scholten, T. 2013a. The spectrum-based learner: A new local
approach for modeling soil vis-NIR spectra of complex datasets. Geoderma
195-196, 268-279.

Ramirez-Lopez, L., Behrens, T., Schmidt, K., Viscarra Rossel, R.,
Dematte, J.A.M., Scholten, T. 2013b. Distance and similarity-search
metrics for use with soil vis-NIR spectra. Geoderma 199, 43-53.

## See also

[`dissimilarity`](https://l-ramirez-lopez.github.io/resemble/reference/dissimilarity.md),
[`ncomp_by_opc`](https://l-ramirez-lopez.github.io/resemble/reference/ncomp_selection.md)

## Author

[Leonardo Ramirez-Lopez](https://orcid.org/0000-0002-5369-5120)

## Examples

``` r
# \donttest{
library(prospectr)
#> prospectr version 0.2.8 -- galo
#> check the package repository at: https://github.com/l-ramirez-lopez/prospectr
data(NIRsoil)

sg <- savitzkyGolay(NIRsoil$spc, p = 3, w = 11, m = 0)
NIRsoil$spc <- sg

Yr <- NIRsoil$Nt[as.logical(NIRsoil$train)]
Xr <- NIRsoil$spc[as.logical(NIRsoil$train), ]

# Compute PCA-based dissimilarity
d <- dissimilarity(Xr, diss_method = diss_pca(ncomp = 8))

# Evaluate using side information
ev <- diss_evaluate(d$dissimilarity, side_info = as.matrix(Yr))
ev$eval
#>                 rmsd         r
#> side_info_1 0.644268 0.8404257

# Evaluate with multiple side information variables
Yr_2 <- NIRsoil$CEC[as.logical(NIRsoil$train)]
ev_2 <- diss_evaluate(d$dissimilarity, side_info = cbind(Yr, Yr_2))
ev_2$eval
#>          rmsd         r
#> Yr   0.644268 0.8404257
#> Yr_2 4.718934 0.6586788
ev_2$global_eval
#>   mean_standardized_rmsd    mean_r
#> 1              0.5814115 0.7495522
# }
```
