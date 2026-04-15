# Control parameters for global model fitting

Specifies cross-validation settings for the
[`model`](https://l-ramirez-lopez.github.io/resemble/reference/model.md)
function.

## Usage

``` r
model_control(validation_type = c("lgo", "none"), number = 10L, p = 0.75)
```

## Arguments

- validation_type:

  a character string specifying the validation method:

  - `"lgo"`: Leave-group-out cross-validation. At each iteration, a
    proportion `p` of observations is retained for training and the
    remainder is used for validation. This is repeated `number` times.

  - `"none"`: No cross-validation is performed.

- number:

  an integer indicating the number of cross-validation iterations. Only
  used when `validation_type = "lgo"`. Default is 10.

- p:

  a numeric value between 0 and 1 indicating the proportion of
  observations to retain for training at each cross-validation
  iteration. Only used when `validation_type = "lgo"`. Default is 0.75.

## Value

A list of class `"model_control"` containing the specified parameters.

## See also

[`model`](https://l-ramirez-lopez.github.io/resemble/reference/model.md)

## Author

[Leonardo Ramirez-Lopez](https://orcid.org/0000-0002-5369-5120)

## Examples

``` r
# Default settings (leave-group-out CV with 10 iterations)
model_control()
#> Model control parameters:
#>   validation_type : lgo 
#>   number          : 10 
#>   p               : 0.75 

# No cross-validation
model_control(validation_type = "none")
#> Model control parameters:
#>   validation_type : none 

# Custom CV settings
model_control(validation_type = "lgo", number = 20, p = 0.80)
#> Model control parameters:
#>   validation_type : lgo 
#>   number          : 20 
#>   p               : 0.8 
```
