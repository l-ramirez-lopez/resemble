# Component selection methods

Constructor functions for specifying how to select the number of
components in projection-based dissimilarity methods
([`diss_pca()`](https://l-ramirez-lopez.github.io/resemble/reference/diss_pca.md),
[`diss_pls()`](https://l-ramirez-lopez.github.io/resemble/reference/diss_pls.md)).

## Usage

``` r
ncomp_by_var(min_var = 0.01, max_ncomp = 40L)

ncomp_by_cumvar(min_cumvar = 0.99, max_ncomp = 40L)

ncomp_by_opc(max_ncomp = 40L)

ncomp_fixed(ncomp)
```

## Arguments

- min_var:

  Numeric in (0, 1\]. Minimum variance a single component must explain
  to be retained.

- max_ncomp:

  Positive integer. Maximum number of components to compute or evaluate.

- min_cumvar:

  Numeric in (0, 1\]. Minimum cumulative variance that the retained
  components must explain.

- ncomp:

  Positive integer. Exact number of components to use.

## Value

An object of class `"ncomp_selection"` with a subclass indicating the
method:

- `ncomp_by_var`: class `c("ncomp_by_var", "ncomp_selection")`

- `ncomp_by_cumvar`: class `c("ncomp_by_cumvar", "ncomp_selection")`

- `ncomp_by_opc`: class `c("ncomp_by_opc", "ncomp_selection")`

- `ncomp_fixed`: class `c("ncomp_fixed", "ncomp_selection")`

## Details

Four selection methods are available:

- `ncomp_by_var()`:

  Retains components that individually explain at least `min_var`
  proportion of variance.

- `ncomp_by_cumvar()`:

  Retains the minimum number of components whose combined explained
  variance reaches `min_cumvar`.

- `ncomp_by_opc()`:

  Optimized principal component selection based on side information
  (Ramirez-Lopez et al., 2013). The optimal number of components
  minimizes the RMSD between each observation's response and its nearest
  neighbor's response in the projected space. Requires `Yr`.

- `ncomp_fixed()`:

  Uses exactly `ncomp` components with no automatic selection.
  Equivalent to passing an integer directly.

At runtime, `max_ncomp` is capped at `min(max_ncomp, nrow(X), ncol(X))`.

## References

Ramirez-Lopez, L., Behrens, T., Schmidt, K., Stevens, A., Dematte,
J.A.M., Scholten, T. 2013. The spectrum-based learner: A new local
approach for modeling soil vis-NIR spectra of complex data sets.
Geoderma 195-196, 268-279.

## See also

[`diss_pca()`](https://l-ramirez-lopez.github.io/resemble/reference/diss_pca.md),
[`diss_pls()`](https://l-ramirez-lopez.github.io/resemble/reference/diss_pls.md),
[`dissimilarity()`](https://l-ramirez-lopez.github.io/resemble/reference/dissimilarity.md)

## Author

[Leonardo Ramirez-Lopez](https://orcid.org/0000-0002-5369-5120)

## Examples

``` r
# Retain components explaining >= 1% variance each
ncomp_by_var(0.01)
#> Component selection: by per-component variance >= 0.01 (max: 40) 

# Retain enough components for 99% cumulative variance
ncomp_by_cumvar(0.99)
#> Component selection: by cumulative variance >= 0.99 (max: 40) 

# Optimize using side information (requires Yr)
ncomp_by_opc(max_ncomp = 40)
#> Component selection: by OPC (max: 40) 

# Fix at exactly 10 components
ncomp_fixed(10)
#> Component selection: fixed: 10 

# Usage in dissimilarity constructors
diss_pca(ncomp = ncomp_by_var(0.01))
#> Dissimilarity: PCA
#>   method            : pca 
#>   ncomp             : var >= 0.01 (max: 40) 
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
diss_pca(ncomp = 10)  
#> Dissimilarity: PCA
#>   method            : pca 
#>   ncomp             : fixed: 10 
#>   center            : TRUE 
#>   scale             : FALSE 
#>   return_projection : FALSE 
```
