# Simple global models

> *It is vain to do with more what can be done with less* – William of
> Ockham

![](logo.png)

## 1 Introduction

The
[`model()`](https://l-ramirez-lopez.github.io/resemble/reference/model.md)
function provides a simple interface for building global calibration
models using partial least squares (PLS) or Gaussian process regression
(GPR). Unlike
[`mbl()`](https://l-ramirez-lopez.github.io/resemble/reference/mbl.md),
which builds local models for each prediction,
[`model()`](https://l-ramirez-lopez.github.io/resemble/reference/model.md)
fits a single model to the entire reference set.

## 2 Fitting methods

Regression methods are specified via constructor functions. The
available methods are:

| Constructor                                                                        | Method               | Description                                                                                                                                                                                           | Reference                                                   |
|------------------------------------------------------------------------------------|----------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------|
| `fit_pls(method = "pls")`                                                          | PLS                  | Non-linear iterative partial least squares (NIPALS)                                                                                                                                                   | Wold ([1975](#ref-wold1975soft))                            |
| `fit_pls(method = "simpls")`                                                       | SIMPLS               | Straightforward implementation of PLS                                                                                                                                                                 | De Jong ([1993](#ref-de1993simpls))                         |
| `fit_pls(method = "mpls")`                                                         | Modified PLS         | Uses correlation instead of covariance for weights                                                                                                                                                    | Shenk and Westerhaus ([1991](#ref-shenk1991populations))    |
| `fit_wapls(method = "mpls")`                                                       | Weighted average PLS | Predictions are a weighted average of multiple models from multiple PLS factors. **NOT YET AVAIABLE FOR THE** [`model()`](https://l-ramirez-lopez.github.io/resemble/reference/model.md) **FUNCTION** | Shenk et al. ([1997](#ref-shenk1997investigation))          |
| [`fit_gpr()`](https://l-ramirez-lopez.github.io/resemble/reference/fit_methods.md) | GPR                  | Gaussian process regression with dot product kernel                                                                                                                                                   | Rasmussen and Williams ([2006](#ref-rasmussen2006gaussian)) |

``` r
# PLS with 15 components
fit_pls(ncomp = 15)
```

``` ansi
Fitting method: pls
   ncomp    : 15
   method   : pls
   scale    : FALSE
   max_iter : 100
   tol      : 1e-06 
```

``` r
# SIMPLS with scaling
fit_pls(ncomp = 15, method = "simpls", scale = TRUE)
```

``` ansi
Fitting method: pls
   ncomp    : 15
   method   : simpls
   scale    : TRUE
   max_iter : 100
   tol      : 1e-06 
```

``` r
# mPLS with scaling
fit_pls(ncomp = 15, method = "mpls")
```

``` ansi
Fitting method: pls
   ncomp    : 15
   method   : mpls
   scale    : FALSE
   max_iter : 100
   tol      : 1e-06 
```

``` r
# GPR with default noise variance
fit_gpr()
```

``` ansi
Fitting method: gpr
   noise_variance : 0.001
   center         : TRUE
   scale          : TRUE 
```

## 3 Example

### 3.1 Data preparation

``` r
library(prospectr)

wavs <- as.numeric(colnames(NIRsoil$spc))

# Preprocess: detrend + first derivative
NIRsoil$spc_pr <- savitzkyGolay(
  detrend(NIRsoil$spc, wav = wavs),
  m = 1, p = 1, w = 7
)

# Split data
train_x <- NIRsoil$spc_pr[NIRsoil$train == 1, ]
train_y <- NIRsoil$Ciso[NIRsoil$train == 1]
test_x  <- NIRsoil$spc_pr[NIRsoil$train == 0, ]
test_y  <- NIRsoil$Ciso[NIRsoil$train == 0]

# Remove missing values
ok_train <- !is.na(train_y)
ok_test  <- !is.na(test_y)
train_x <- train_x[ok_train, ]
train_y <- train_y[ok_train]
test_x  <- test_x[ok_test, ]
test_y  <- test_y[ok_test]
```

### 3.2 Fitting a PLS model

``` r
set.seed(1124) # guarantee same CV splits for all methods
pls_mod <- model(
 Xr = train_x,
 Yr = train_y,
 fit_method = fit_pls(ncomp = 12, method = "mpls", scale = TRUE),
 control = model_control(validation_type = "lgo", number = 10)
)
```

    Running cross-validation...

    Fitting model...

``` r
pls_mod
```

``` ansi
--- Global resemble model ---
Method:       pls
Observations: 548
Variables:    694
_______________________________________________________
Fit method
Fitting method: pls
   ncomp    : 12
   method   : mpls
   scale    : TRUE
   max_iter : 100
   tol      : 1e-06
_______________________________________________________
Cross-validation
  Best ncomp: 10

 ncomp rmse rmse_sd st_rmse    r2
     1 1.34   0.160  0.1166 0.447
     2 1.21   0.162  0.1048 0.558
     3 1.14   0.169  0.0991 0.606
     4 1.14   0.163  0.0989 0.610
     5 1.10   0.174  0.0947 0.641
     6 1.09   0.171  0.0944 0.644
     7 1.06   0.185  0.0911 0.667
     8 1.06   0.172  0.0917 0.665
     9 1.07   0.156  0.0925 0.663
    10 1.06   0.163  0.0912 0.672
    11 1.07   0.158  0.0923 0.668
    12 1.07   0.141  0.0926 0.670
_______________________________________________________ 
```

### 3.3 Fitting a GPR model

``` r
set.seed(1124) # guarantee same CV splits for all methods
gpr_mod <- model(
 Xr = train_x,
 Yr = train_y,
 fit_method = fit_gpr(noise_variance = 0.5),
 control = model_control(validation_type = "lgo", number = 10)
)
```

    Running cross-validation...

    Fitting model...

``` r
gpr_mod
```

``` ansi
--- Global resemble model ---
Method:       gpr
Observations: 548
Variables:    694
_______________________________________________________
Fit method
Fitting method: gpr
   noise_variance : 0.5
   center         : TRUE
   scale          : TRUE
_______________________________________________________
Cross-validation
 rmse rmse_sd st_rmse    r2
 1.15   0.184   0.104 0.654
_______________________________________________________ 
```

### 3.4 Prediction

``` r
# Predict using optimal number of components (from CV)
pls_pred <- predict(pls_mod, newdata = test_x)

# Predict with GPR
gpr_pred <- predict(gpr_mod, newdata = test_x)

# Compare predictions
data.frame(
  observed = test_y,
  pls = pls_pred[, which(pls_mod$cv_results$optimal)],
  gpr = as.vector(gpr_pred)
) |> head(10)
```

``` ansi
    observed       pls        gpr
619     0.15 2.2893489 1.56360594
620     0.39 1.6013062 0.08271392
623     0.59 0.6469817 0.85006225
625     0.70 1.8391486 1.01129163
629     0.64 0.2144320 0.31474586
634     0.85 0.2687261 0.70627228
635     0.89 2.5970168 1.62297992
636     0.91 0.6575313 0.31100348
637     0.90 1.2757025 1.71785327
638     0.90 1.6560021 1.07326142
```

### 3.5 Validation statistics

``` r
# Function to compute stats
eval_pred <- function(obs, pred) {
  data.frame(
    rmse = sqrt(mean((obs - pred)^2)),
    r2 = cor(obs, pred)^2
  )
}
```

``` r
## PLS evaluation
eval_pred(test_y, pls_pred[, which(pls_mod$cv_results$optimal)])
```

``` ansi
       rmse        r2
1 0.8955579 0.6682104
```

``` r
## GPR evaluation
eval_pred(test_y, gpr_pred)
```

``` ansi
       rmse predicted
1 0.8215557 0.7225666
```

## References

De Jong, S., 1993. SIMPLS: An alternative approach to partial least
squares regression. Chemometrics and intelligent laboratory systems 18,
251–263.

Rasmussen, C.E., Williams, C.K.I., 2006. Gaussian processes for machine
learning. MIT Press, Cambridge, MA.

Shenk, J.S., Westerhaus, M.O., 1991. Populations structuring of near
infrared spectra and modified partial least squares regression. Crop
science 31, 1548–1555.

Shenk, J.S., Westerhaus, M.O., Berzaghi, P., 1997. Investigation of a
LOCAL calibration procedure for near infrared instruments. Journal of
Near Infrared Spectroscopy 5, 223–232.

Wold, H., 1975. Soft modelling by latent variables: The non-linear
iterative partial least squares (NIPALS) approach. Journal of Applied
Probability 12, 117–142.
