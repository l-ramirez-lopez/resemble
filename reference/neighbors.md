# Neighbor selection methods

These functions create configuration objects that specify how neighbors
are selected for memory-based learning in
[`mbl`](https://l-ramirez-lopez.github.io/resemble/reference/mbl.md) an
[`liblex`](https://l-ramirez-lopez.github.io/resemble/reference/liblex.md).

## Usage

``` r
neighbors_k(k)

neighbors_diss(threshold, k_min = 4L, k_max = Inf)
```

## Arguments

- k:

  Integer vector. One or more neighborhood sizes to evaluate. Values
  will be sorted in ascending order. Minimum allowed value is 4.

- threshold:

  Numeric vector. One or more dissimilarity thresholds. Neighbors are
  selected if their dissimilarity to the target observation is below the
  threshold. Values will be sorted in ascending order.

- k_min:

  Integer. Minimum number of neighbors to retain, regardless of
  threshold. Default `4L`.

- k_max:

  Integer or `Inf`. Maximum number of neighbors to retain, regardless of
  threshold. Default `Inf` (no upper bound other than the size of the
  reference set).

## Value

An object of class `c("neighbors_k", "neighbors")` or
`c("neighbors_diss", "neighbors")`, containing the validated parameters.
Intended to be passed to
[`mbl`](https://l-ramirez-lopez.github.io/resemble/reference/mbl.md).

## Details

Two strategies are available for neighbor selection:

**Fixed-k selection** (`neighbors_k`)

A fixed number of nearest neighbors is selected for each target
observation. Multiple values of `k` can be provided to evaluate
different neighborhood sizes.

**Dissimilarity-threshold selection** (`neighbors_diss`)

Neighbors are selected based on a dissimilarity threshold. All reference
observations with dissimilarity below the threshold are included. The
`k_min` and `k_max` arguments provide bounds to ensure a reasonable
neighborhood size regardless of the threshold. Multiple thresholds can
be provided to evaluate different settings.

## See also

[`mbl`](https://l-ramirez-lopez.github.io/resemble/reference/mbl.md)

## Examples

``` r
# Fixed neighborhood sizes
neighbors_k(k = 50)
#> Neighbor selection: fixed k
#>   k : 50 
neighbors_k(k = c(40, 60, 80, 100, 120))
#> Neighbor selection: fixed k
#>   k : 40, 60, 80, 100, 120 

# Dissimilarity threshold with default bounds
neighbors_diss(threshold = 0.3)
#> Neighbor selection: dissimilarity threshold
#>    threshold : 0.3 
#>    k_min     : 4 
#>    k_max     : Inf 

# Dissimilarity threshold with custom bounds
neighbors_diss(threshold = c(0.1, 0.2, 0.3), k_min = 10, k_max = 150)
#> Neighbor selection: dissimilarity threshold
#>    threshold : 0.1, 0.2, 0.3 
#>    k_min     : 10 
#>    k_max     : 150 
```
