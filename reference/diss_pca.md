# PCA dissimilarity method constructor

Creates a configuration object for computing dissimilarity based on
Mahalanobis distance in PCA score space.

## Usage

``` r
diss_pca(
  ncomp = ncomp_by_var(0.01),
  method = c("pca", "pca_nipals"),
  center = TRUE,
  scale = FALSE,
  return_projection = FALSE
)
```

## Arguments

- ncomp:

  Component selection method. Can be:

  - A positive integer for a fixed number of components

  - [`ncomp_fixed`](https://l-ramirez-lopez.github.io/resemble/reference/ncomp_selection.md)`(n)`:
    explicit fixed selection

  - [`ncomp_by_var`](https://l-ramirez-lopez.github.io/resemble/reference/ncomp_selection.md)`(min_var)`:
    retain components explaining at least `min_var` variance each
    (default: `ncomp_by_var(0.01)`)

  - [`ncomp_by_cumvar`](https://l-ramirez-lopez.github.io/resemble/reference/ncomp_selection.md)`(min_cumvar)`:
    retain components until cumulative variance reaches `min_cumvar`

  - [`ncomp_by_opc()`](https://l-ramirez-lopez.github.io/resemble/reference/ncomp_selection.md):
    optimize using side information (`Yr` required in
    [`dissimilarity()`](https://l-ramirez-lopez.github.io/resemble/reference/dissimilarity.md))

- method:

  Character. PCA algorithm: `"pca"` (default, SVD-based) or
  `"pca_nipals"` (NIPALS algorithm).

- center:

  Logical. Center data before projection? Default `TRUE`.

- scale:

  Logical. Scale data before projection? Default `FALSE`.

- return_projection:

  Logical. Return the projection object? Default `FALSE`.

## Value

An object of class `c("diss_pca", "diss_method")`.

## See also

Component selection:
[`ncomp_by_var`](https://l-ramirez-lopez.github.io/resemble/reference/ncomp_selection.md),
[`ncomp_by_cumvar`](https://l-ramirez-lopez.github.io/resemble/reference/ncomp_selection.md),
[`ncomp_by_opc`](https://l-ramirez-lopez.github.io/resemble/reference/ncomp_selection.md),
[`ncomp_fixed`](https://l-ramirez-lopez.github.io/resemble/reference/ncomp_selection.md)

Other dissimilarity methods:
[`diss_pls`](https://l-ramirez-lopez.github.io/resemble/reference/diss_pls.md),
[`diss_correlation`](https://l-ramirez-lopez.github.io/resemble/reference/diss_correlation.md),
[`diss_euclidean`](https://l-ramirez-lopez.github.io/resemble/reference/diss_euclidean.md),
[`diss_cosine`](https://l-ramirez-lopez.github.io/resemble/reference/diss_cosine.md),
[`diss_mahalanobis`](https://l-ramirez-lopez.github.io/resemble/reference/diss_mahalanobis.md)

## Examples

``` r
# Fixed number of components
diss_pca(ncomp = 10)
#> Dissimilarity: PCA
#>   method            : pca 
#>   ncomp             : fixed: 10 
#>   center            : TRUE 
#>   scale             : FALSE 
#>   return_projection : FALSE 
diss_pca(ncomp = ncomp_fixed(10))
#> Dissimilarity: PCA
#>   method            : pca 
#>   ncomp             : fixed: 10 
#>   center            : TRUE 
#>   scale             : FALSE 
#>   return_projection : FALSE 

# Retain components explaining >= 1% variance each (default)
diss_pca(ncomp = ncomp_by_var(0.01))
#> Dissimilarity: PCA
#>   method            : pca 
#>   ncomp             : var >= 0.01 (max: 40) 
#>   center            : TRUE 
#>   scale             : FALSE 
#>   return_projection : FALSE 

# Retain components until 99% cumulative variance
diss_pca(ncomp = ncomp_by_cumvar(0.99))
#> Dissimilarity: PCA
#>   method            : pca 
#>   ncomp             : cumvar >= 0.99 (max: 40) 
#>   center            : TRUE 
#>   scale             : FALSE 
#>   return_projection : FALSE 

# Optimize using side information (requires Yr)
diss_pca(ncomp = ncomp_by_opc(40))
#> Dissimilarity: PCA
#>   method            : pca 
#>   ncomp             :  
#>   center            : TRUE 
#>   scale             : FALSE 
#>   return_projection : FALSE 
diss_pca(ncomp = ncomp_by_opc())
#> Dissimilarity: PCA
#>   method            : pca 
#>   ncomp             :  
#>   center            : TRUE 
#>   scale             : FALSE 
#>   return_projection : FALSE 

# NIPALS algorithm (useful for very large matrices)
diss_pca(ncomp = 10, method = "pca_nipals")
#> Dissimilarity: PCA
#>   method            : pca_nipals 
#>   ncomp             : fixed: 10 
#>   center            : TRUE 
#>   scale             : FALSE 
#>   return_projection : FALSE 
```
