# A function that controls some aspects of the the `gesearch` function

This function is used to further control some aspects of the `gesearch`
function.

## Usage

``` r
search_control(
  retain_by = "proportion",
  percentile_type = if(retain_by == "probability") 7,
  tune = FALSE,
  number = 10,
  p = 0.75,
  stagnation_limit = 5,
  allow_parallel = TRUE
)
```

## Arguments

- retain_by:

  a character string indicating how the training observations to be kept
  at each iteration must be selected. Two options are available:
  `'proportion'` (default) retains a fix proportion of observations at
  each iteration and `'probability'` retains observations whose errors
  are below an estimated percentile (cut point) from a given
  probability. See details.

- percentile_type:

  if `retain_by = 'probability'`, an integer between 1 and 9 to be
  passed to the the `type` argument of the
  [`quantile`](https://rdrr.io/r/stats/quantile.html) function to
  estimate the percentile cut-off value used to select the observations
  at each iteration. Default is 7 (as in
  [`quantile`](https://rdrr.io/r/stats/quantile.html)). If
  `retain_by != 'proportion'`, this argument is ignored. See details.

- tune:

  a logical indicating whether to tune the parameters of the regression
  method (e.g. pls factors). If `TRUE` the parameters are tuned by using
  cross-validation which is based on the `number` and `p` arguments.
  Note that cross-validation greatly increases processing time. Default
  is `FALSE`.

- number:

  if `tune = TRUE`, an integer indicating the number of groups for the
  leave-group-out internal cross-validation used for parameter tuning
  (e.g. pls factors) at each iteration. These groups are built from the
  training subset selected at the respective iteration. Default is 10.

- p:

  if `tune = TRUE`, a value indicating the proportion of observations to
  build each group in the the leave-group-out internal cross-validation
  (used for parameter tuning at each iteration). These groups are built
  from the training subset selected at the respective iteration. Default
  is 0.75 (i.e. 75 percent).

- stagnation_limit:

  Integer. Maximum number of consecutive iterations/generations during
  which the gene pool size is unchanged before early termination is
  triggered. This prevents infinite loops in cases where the target size
  cannot be reached due to no further removable samples. Default is 5.

- allow_parallel:

  set to TRUE to parallelise the internal re-sampling and calibration.
  The parallel backend should be registered by the user.

## Value

a `list` mirroring the specified parameters.

## Details

The training observations to be kept at each iteration can be selected
by using a fixed proportion of observations (in relation to the number
of observations in the training subset being used at each iteration)
which associated errors are the lowest, in this case the argument
`retain_by` must be set to `'proportion'`. This proportion (a value
larger than 0 and below 1) is given in the in the argument `retain` of
the
[`gesearch`](https://l-ramirez-lopez.github.io/resemble/reference/gesearch.md)
function.

For selecting the training observations, it is also possible to retain
the ones whose associated errors are below a percentile estimated for a
given probability (the percentile is estimated with the
[`quantile`](https://rdrr.io/r/stats/quantile.html) function). In this
case the probability (a value larger than 0 and below 1) must be
provided in the argument `retain`, this numeric value is passed to the
argument `probs` in the
[`quantile`](https://rdrr.io/r/stats/quantile.html) function. This
method might be useful in the presence of outlier observations with
considerable deviation in their associated errors.

The internal validations for tuning the number of pls factors are based
on the leave-group-out cross-validation. Arguments `p` and `number` are
used to control the number of groups and the amount of observations per
group.

## References

Lobsey, C. R., Viscarra Rossel, R. A., Roudier, P., & Hedley, C. B.
2017. rs-local data-mines information from spectral libraries to improve
local calibrations. European Journal of Soil Science, 68(6), 840-852.

## See also

[`gesearch`](https://l-ramirez-lopez.github.io/resemble/reference/gesearch.md),
[`quantile`](https://rdrr.io/r/stats/quantile.html)

## Author

Leonardo Ramirez-Lopez

## Examples

``` r
if (FALSE) { # \dontrun{
# A control list with the default parameters
search_control()

# A control list which specifies that observations must be
# retained based on the associated errors that are below
# a given probability.
search_control(retain_by = "probability")
} # }
```
