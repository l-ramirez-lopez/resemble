# Search neighbors in a reference set

Searches for the nearest neighbors of observations in a reference set or
between two sets of observations.

## Usage

``` r
search_neighbors(Xr, Xu = NULL,
                 diss_method = diss_pca(), Yr = NULL,
                 neighbors, spike = NULL,
                 return_dissimilarity = FALSE, 
                 k, k_diss, k_range, pc_selection, 
                 center, scale, documentation, ...
                 )
```

## Arguments

- Xr:

  A numeric matrix of reference observations (rows) and variables
  (columns) where the neighbor search is conducted.

- Xu:

  Optional matrix of observations for which neighbors are to be searched
  in `Xr`.

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

- neighbors:

  A neighbor selection object created by:

  - [`neighbors_k()`](https://l-ramirez-lopez.github.io/resemble/reference/neighbors.md):
    Select k nearest neighbors

  - [`neighbors_diss()`](https://l-ramirez-lopez.github.io/resemble/reference/neighbors.md):
    Select neighbors by dissimilarity threshold

- spike:

  Optional integer vector indicating observations in `Xr` to force into
  (positive indices) or exclude from (negative indices) neighborhoods.

- return_dissimilarity:

  Logical indicating whether to return the dissimilarity matrix. Default
  is `FALSE`.

- k:

  Deprecated.

- k_diss:

  Deprecated.

- k_range:

  Deprecated.

- pc_selection:

  Deprecated.

- center:

  Deprecated.

- scale:

  Deprecated.

- documentation:

  Deprecated.

- ...:

  Additional arguments (currently unused).

## Value

A list containing:

- neighbors:

  Matrix of `Xr` indices for each query observation's neighbors, sorted
  by dissimilarity (columns = query observations).

- neighbors_diss:

  Matrix of dissimilarity scores corresponding to `neighbors`.

- unique_neighbors:

  Vector of unique `Xr` indices that appear in any neighborhood.

- k_diss_info:

  If
  [`neighbors_diss()`](https://l-ramirez-lopez.github.io/resemble/reference/neighbors.md)
  was used, a `data.frame` with columns for observation index, number of
  neighbors found, and final number after applying bounds.

- dissimilarity:

  If `return_dissimilarity = TRUE`, the full dissimilarity object.

- projection:

  If the dissimilarity method includes `return_projection = TRUE`, the
  projection object.

- gh:

  If the dissimilarity method includes `gh = TRUE`, the GH distances.

## Details

This function is useful for reducing large reference sets by identifying
only relevant neighbors before running
[`mbl`](https://l-ramirez-lopez.github.io/resemble/reference/mbl.md).

If `Xu` is not provided, the function searches for neighbors within `Xr`
itself (excluding self-matches). If `Xu` is provided, neighbors of each
observation in `Xu` are searched in `Xr`.

The `spike` argument allows forcing specific observations into or out of
all neighborhoods. Positive indices are always included; negative
indices are always excluded.

## References

Ramirez-Lopez, L., Behrens, T., Schmidt, K., Stevens, A., Dematte,
J.A.M., Scholten, T. 2013a. The spectrum-based learner: A new local
approach for modeling soil vis-NIR spectra of complex data sets.
Geoderma 195-196, 268-279.

Ramirez-Lopez, L., Behrens, T., Schmidt, K., Viscarra Rossel, R.,
Dematte, J.A.M., Scholten, T. 2013b. Distance and similarity-search
metrics for use with soil vis-NIR spectra. Geoderma 199, 43-53.

## See also

[`dissimilarity`](https://l-ramirez-lopez.github.io/resemble/reference/dissimilarity.md),
[`mbl`](https://l-ramirez-lopez.github.io/resemble/reference/mbl.md),
[`neighbors_k`](https://l-ramirez-lopez.github.io/resemble/reference/neighbors.md),
[`neighbors_diss`](https://l-ramirez-lopez.github.io/resemble/reference/neighbors.md)

## Author

[Leonardo Ramirez-Lopez](https://orcid.org/0000-0002-5369-5120)

## Examples

``` r
# \donttest{
library(prospectr)
data(NIRsoil)

Xu <- NIRsoil$spc[!as.logical(NIRsoil$train), ]
Yu <- NIRsoil$CEC[!as.logical(NIRsoil$train)]
Yr <- NIRsoil$CEC[as.logical(NIRsoil$train)]
Xr <- NIRsoil$spc[as.logical(NIRsoil$train), ]

Xu <- Xu[!is.na(Yu), ]
Yu <- Yu[!is.na(Yu)]
Xr <- Xr[!is.na(Yr), ]
Yr <- Yr[!is.na(Yr)]

# Correlation-based neighbor search with k neighbors
ex1 <- search_neighbors(
  Xr = Xr, Xu = Xu,
  diss_method = diss_correlation(),
  neighbors = neighbors_k(40)
)

# PCA-based with OPC selection
ex2 <- search_neighbors(
  Xr = Xr, Xu = Xu,
  diss_method = diss_pca(
    ncomp = ncomp_by_opc(40),
    scale = TRUE,
    return_projection = TRUE
  ),
  Yr = Yr,
  neighbors = neighbors_k(50)
)

# Observations not in any neighborhood
setdiff(seq_len(nrow(Xr)), ex2$unique_neighbors)
#>  [1]   3  16  20  23  26  27  42  44  59  61  86  96 103 111 123 149 231 267 279
#> [20] 298 310 311 326 328 330

# Dissimilarity threshold-based selection
ex3 <- search_neighbors(
  Xr = Xr, Xu = Xu,
  diss_method = diss_pls(
    ncomp = ncomp_by_opc(40),
    scale = TRUE
  ),
  Yr = Yr,
  neighbors = neighbors_diss(threshold = 0.5, k_min = 10, k_max = 100)
)
# }
```
