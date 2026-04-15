# Compute dissimilarity matrices

Computes dissimilarity matrices between observations using various
methods. This is the main interface for dissimilarity computation in the
resemble package.

## Usage

``` r
dissimilarity(Xr, Xu = NULL, diss_method = diss_pca(), Yr = NULL)
```

## Arguments

- Xr:

  A numeric matrix of reference observations (rows) and variables
  (columns).

- Xu:

  Optional matrix of additional observations with the same variables.

- diss_method:

  A dissimilarity method object created by one of:

  - [`diss_pca()`](https://l-ramirez-lopez.github.io/resemble/reference/diss_pca.md):
    Mahalanobis distance in PCA space

  - [`diss_pls()`](https://l-ramirez-lopez.github.io/resemble/reference/diss_pls.md):
    Mahalanobis distance in PLS space

  - [`diss_correlation()`](https://l-ramirez-lopez.github.io/resemble/reference/diss_correlation.md):
    Correlation-based dissimilarity

  - [`diss_euclidean()`](https://l-ramirez-lopez.github.io/resemble/reference/diss_euclidean.md):
    Euclidean distance

  - [`diss_mahalanobis()`](https://l-ramirez-lopez.github.io/resemble/reference/diss_mahalanobis.md):
    Mahalanobis distance

  - [`diss_cosine()`](https://l-ramirez-lopez.github.io/resemble/reference/diss_cosine.md):
    Cosine dissimilarity

  Default is
  [`diss_pca()`](https://l-ramirez-lopez.github.io/resemble/reference/diss_pca.md).

- Yr:

  Optional response matrix. Required for PLS methods and when using
  [`ncomp_by_opc()`](https://l-ramirez-lopez.github.io/resemble/reference/ncomp_selection.md).

## Value

A list of class `"dissimilarity"` containing:

- dissimilarity:

  The computed dissimilarity matrix. Dimensions are `nrow(Xr)`
  \\\times\\ `nrow(Xr)` when `Xu = NULL`, or `nrow(Xr)` \\\times\\
  `nrow(Xu)` otherwise.

- diss_method:

  The `diss_*` constructor object used for computation.

- center:

  Vector used to center the data.

- scale:

  Vector used to scale the data.

- ncomp:

  Number of components used (for projection methods).

- projection:

  If `return_projection = TRUE` in the method constructor, the
  `ortho_projection` object.

## Details

The function dispatches to the appropriate internal computation based on
the class of `diss_method`. Each method constructor (e.g.,
[`diss_pca()`](https://l-ramirez-lopez.github.io/resemble/reference/diss_pca.md))
encapsulates all method-specific parameters including component
selection, centering, scaling, and whether to return projections.

### Output dimensions

When only `Xr` is provided, the function computes pairwise
dissimilarities among all observations in `Xr`, returning a symmetric
`nrow(Xr)` \\\times\\ `nrow(Xr)` matrix.

When both `Xr` and `Xu` are provided, the function computes
dissimilarities between each observation in `Xr` and each observation in
`Xu`, returning a `nrow(Xr)` \\\times\\ `nrow(Xu)` matrix where element
\\(i, j)\\ is the dissimilarity between the \\i\\-th observation in `Xr`
and the \\j\\-th observation in `Xu`.

### Mahalanobis distance

Note that
[`diss_mahalanobis()`](https://l-ramirez-lopez.github.io/resemble/reference/diss_mahalanobis.md)
computes Mahalanobis distance directly on the input variables. This
requires the covariance matrix to be invertible, which fails when the
number of variables exceeds the number of observations or when variables
are highly correlated (common in spectral data). For such cases, use
[`diss_pca()`](https://l-ramirez-lopez.github.io/resemble/reference/diss_pca.md)
or
[`diss_pls()`](https://l-ramirez-lopez.github.io/resemble/reference/diss_pls.md)
instead.

## References

Ramirez-Lopez, L., Behrens, T., Schmidt, K., Stevens, A., Dematte,
J.A.M., Scholten, T. 2013a. The spectrum-based learner: A new local
approach for modeling soil vis-NIR spectra of complex data sets.
Geoderma 195-196, 268-279.

Ramirez-Lopez, L., Behrens, T., Schmidt, K., Viscarra Rossel, R.,
Dematte, J.A.M., Scholten, T. 2013b. Distance and similarity-search
metrics for use with soil vis-NIR spectra. Geoderma 199, 43-53.

## See also

[`diss_pca`](https://l-ramirez-lopez.github.io/resemble/reference/diss_pca.md),
[`diss_pls`](https://l-ramirez-lopez.github.io/resemble/reference/diss_pls.md),
[`diss_correlation`](https://l-ramirez-lopez.github.io/resemble/reference/diss_correlation.md),
[`diss_euclidean`](https://l-ramirez-lopez.github.io/resemble/reference/diss_euclidean.md),
[`diss_mahalanobis`](https://l-ramirez-lopez.github.io/resemble/reference/diss_mahalanobis.md),
[`diss_cosine`](https://l-ramirez-lopez.github.io/resemble/reference/diss_cosine.md)

## Author

[Leonardo Ramirez-Lopez](https://orcid.org/0000-0002-5369-5120)

## Examples

``` r
# \donttest{
library(prospectr)
data(NIRsoil)

# Preprocess
sg <- savitzkyGolay(NIRsoil$spc, m = 1, p = 4, w = 15)

Xr <- sg[as.logical(NIRsoil$train), ]
Xu <- sg[!as.logical(NIRsoil$train), ]
Yr <- NIRsoil$CEC[as.logical(NIRsoil$train)]
Yu <- NIRsoil$CEC[!as.logical(NIRsoil$train)]

Xu <- Xu[!is.na(Yu), ]
Xr <- Xr[!is.na(Yr), ]
Yr <- Yr[!is.na(Yr)]

# PCA-based dissimilarity with variance-based selection
d1 <- dissimilarity(Xr, Xu, diss_method = diss_pca())

# PCA with OPC selection (requires Yr)
d2 <- dissimilarity(Xr, Xu,
  Yr = Yr,
  diss_method = diss_pca(
    ncomp = ncomp_by_opc(30),
    return_projection = TRUE
  )
)

# PLS-based dissimilarity 
d3 <- dissimilarity(
  Xr, Xu,
  Yr = Yr,
  diss_method = diss_pls(
    ncomp = ncomp_by_opc(30)
  )
)

# Euclidean distance
d4 <- dissimilarity(Xr, Xu, diss_method = diss_euclidean())

# Correlation dissimilarity with moving window
d5 <- dissimilarity(Xr, Xu, diss_method = diss_correlation(ws = 41))

# Mahalanobis distance (use only when n > p and low collinearity)
# d6 <- dissimilarity(Xr[, 1:20], Xu[, 1:20],
#                     diss_method = diss_mahalanobis())
# }
```
