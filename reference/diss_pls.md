# PLS dissimilarity method constructor

Creates a configuration object for computing dissimilarity based on
Mahalanobis distance in PLS score space. Requires `Yr` in
[`dissimilarity()`](https://l-ramirez-lopez.github.io/resemble/reference/dissimilarity.md).

## Usage

``` r
diss_pls(
  ncomp = ncomp_by_opc(),
  method = c("pls", "mpls"),
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

  - [`ncomp_by_cumvar`](https://l-ramirez-lopez.github.io/resemble/reference/ncomp_selection.md)`(min_cumvar)`:
    retain components until cumulative variance reaches `min_cumvar`

  - [`ncomp_by_opc()`](https://l-ramirez-lopez.github.io/resemble/reference/ncomp_selection.md):
    optimize using side information (default; recommended for PLS since
    `Yr` is already required)

- method:

  Character. PLS algorithm: `"pls"` (default) or `"mpls"` (modified PLS,
  Shenk & Westerhaus 1991).

- scale:

  Logical. Scale data? Default `FALSE`. Note: PLS always centers
  internally.

- return_projection:

  Logical. Return projection object? Default `FALSE`.

## Value

An object of class `c("diss_pls", "diss_method")`.

## See also

Component selection:
[`ncomp_by_var`](https://l-ramirez-lopez.github.io/resemble/reference/ncomp_selection.md),
[`ncomp_by_cumvar`](https://l-ramirez-lopez.github.io/resemble/reference/ncomp_selection.md),
[`ncomp_by_opc`](https://l-ramirez-lopez.github.io/resemble/reference/ncomp_selection.md),
[`ncomp_fixed`](https://l-ramirez-lopez.github.io/resemble/reference/ncomp_selection.md)

Other dissimilarity methods:
[`diss_pca`](https://l-ramirez-lopez.github.io/resemble/reference/diss_pca.md),
[`diss_correlation`](https://l-ramirez-lopez.github.io/resemble/reference/diss_correlation.md),
[`diss_euclidean`](https://l-ramirez-lopez.github.io/resemble/reference/diss_euclidean.md),
[`diss_cosine`](https://l-ramirez-lopez.github.io/resemble/reference/diss_cosine.md),
[`diss_mahalanobis`](https://l-ramirez-lopez.github.io/resemble/reference/diss_mahalanobis.md)

## Author

[Leonardo Ramirez-Lopez](https://orcid.org/0000-0002-5369-5120)

## Examples

``` r
# Default: OPC optimization (recommended)
diss_pls()
#> Dissimilarity: PLS
#>   method            : pls 
#>   ncomp             :  
#>   scale             : FALSE 
#>   return_projection : FALSE 

# Fixed number of components
diss_pls(ncomp = 15)
#> Dissimilarity: PLS
#>   method            : pls 
#>   ncomp             : fixed: 15 
#>   scale             : FALSE 
#>   return_projection : FALSE 

# Custom opc settings
diss_pls(ncomp = ncomp_by_opc(max_ncomp = 50))
#> Dissimilarity: PLS
#>   method            : pls 
#>   ncomp             :  
#>   scale             : FALSE 
#>   return_projection : FALSE 

# Modified PLS
diss_pls(ncomp = 10, method = "mpls")
#> Dissimilarity: PLS
#>   method            : mpls 
#>   ncomp             : fixed: 10 
#>   scale             : FALSE 
#>   return_projection : FALSE 
```
