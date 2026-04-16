# Evolutionary subset search and modeling using gesearch

> *Particular individuals do not recur, but their building blocks do* –
> ([Holland, 1995](#ref-holland1995hidden))

![](logo.png)

## 1 Introduction

Large reference datasets, such as spectral libraries, can support model
development when only a limited number of target-domain samples are
available. The
[`gesearch()`](https://l-ramirez-lopez.github.io/resemble/reference/gesearch.md)
function implements the *gesearch* algorithm described in Ramirez-Lopez
et al. ([2026](#ref-ramirezlopez2026a)), which builds on the
sample-selection idea introduced by Lobsey et al.
([2017](#ref-lobsey2017rs)). The algorithm uses an evolutionary search
to identify a subset of library samples that is most relevant to a
target domain and then uses the selected subset to fit a
context-specific model.

In its current implementation,
[`gesearch()`](https://l-ramirez-lopez.github.io/resemble/reference/gesearch.md)
uses partial least squares (PLS) models to evaluate candidate subsets.
Sample selection is guided by how well a subset predicts the target
response, reconstructs the target spectra, or aligns with the target
domain in latent space. The method can therefore be used for calibration
transfer or domain adaptation, with or without response values from the
target domain.

It is important to note that the *gesearch* algorithm is designed to
identify a single subset of samples that is most relevant to the target
domain. Accordingly, some data from the target domain must be available.
In addition, **an important assumption underlying** ***gesearch*** **is
that the relationships between spectra and response values in the target
domain are expected to be predominantly linear** ([Ramirez-Lopez et al.,
2026](#ref-ramirezlopez2026a)), even though the search process itself is
non-linear.

This vignette introduces the main workflow for `gesearch`, including how
to define the search, control the selection process, inspect the
selected subsets, and generate predictions from the final model. It also
highlights key tuning arguments such as `k`, `b`, `retain`, and the
control settings available through
[`gesearch_control()`](https://l-ramirez-lopez.github.io/resemble/reference/gesearch_control.md).

## 2 Glossary and conventions

As *gesearch* is an evolutionary search algorithm, it draws on
terminology from evolutionary optimization and sample selection. This
glossary summarizes the main terms used throughout the vignette.

| Term                                 | Meaning                                                                                                                                                                                              |
|:-------------------------------------|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Target set ($\mathcal{D}_{u}$)       | The population of interest (`Xu`, `Yu`), typically small ($m$ samples). Response values may be partially or entirely unavailable. The subscript $u$ denotes the target domain (partially “unknown”). |
| Spectral library ($\mathcal{D}_{r}$) | A large, heterogeneous collection of candidate samples for model training (`Xr`, `Yr`), with $n \gg m$. The subscript $r$ denotes the reference library.                                             |
| Gene                                 | A sample from the spectral library.                                                                                                                                                                  |
| Gene pool                            | The set of genes currently eligible for selection.                                                                                                                                                   |
| Individual                           | A subset of `k` genes used to fit a candidate PLS model.                                                                                                                                             |
| Population                           | The collection of all individuals in a given generation.                                                                                                                                             |
| Silenced gene                        | A gene permanently excluded from all subsequent generations.                                                                                                                                         |
| Weakness score                       | A measure of how poorly a gene contributes to model performance for the target set.                                                                                                                  |
| Generation                           | One iteration of the evolutionary cycle.                                                                                                                                                             |
| Incidence count                      | The number of individuals in which a gene appears.                                                                                                                                                   |

## 3 The *gesearch* algorithm: How it works

### 3.1 Evolutionary search overview

The algorithm treats sample selection as an evolutionary process:

1.  **Initialization**: Create a population of individuals, each
    containing `k` genes randomly selected from the library. A gene, and
    therefore its corresponding sample, can appear in multiple
    individuals. The representation factor `b` controls how often each
    gene is expected to appear across the population, so that genes are
    represented approximately equally.

2.  **Evaluation**: Fit a PLS model for each individual. One or more
    figures of merit are then computed for that model using the
    available target samples. The weakness score of a gene is defined as
    the average performance of all individuals in which that gene is
    active. This procedure is repeated for all genes in the active gene
    pool.

3.  **Silencing**: Permanently mark genes with high weakness scores as
    inactive, preventing them from being selected in subsequent
    generations. In practical terms, the corresponding samples are
    permanently excluded from further selection.

4.  **Reproduction**: Generate new individuals through crossover
    (combining genes from two parents) and mutation (adding genes from
    the remaining pool).

5.  **Termination**: Stop when the gene pool shrinks to `target_size`.

The surviving genes form the final calibration subset, which is expected
to deliver near-optimal performance for the target population.

### 3.2 Weakness score computation

The key idea behind `gesearch` is that gene quality is assessed
indirectly: each individual (subset) receives a performance score, and a
gene’s weakness is defined as the average performance across all
individuals containing that gene.

#### 3.2.1 Step 1: Evaluate each individual

For each individual $S_{i}^{(\gamma)}$ in generation $\gamma$, a PLS
model is fitted using its `k` active genes. This model is then evaluated
on the target spectra $\mathbf{X}_{u}$, producing an individual-level
score, for example, prediction error, reconstruction error, or distance
to the target set.

#### 3.2.2 Step 2: Aggregate to gene-level weakness

Each gene appears in multiple individuals, with its frequency controlled
by the representation factor `b`. The weakness of gene $g_{j}$ is the
mean of the individual-level scores across all individuals containing
that gene:

$$w^{(\gamma)}(g_{j}) = \frac{1}{h_{j}^{(\gamma)}}\sum\limits_{i:\, g_{j} \in S_{i}^{(\gamma)}}\text{score}(S_{i}^{(\gamma)})$$

where $h_{j}^{(\gamma)}$ is the incidence count, that is, the number of
individuals containing gene $g_{j}$.

#### 3.2.3 Interpretation

A gene that consistently appears in poorly performing individuals
accumulates a high weakness score. Conversely, a gene that contributes
to well-performing individuals has a low weakness score. This averaging
mechanism helps isolate the marginal contribution of each gene despite
the combinatorial nature of subset selection.

### 3.3 Weakness functions

The `optimization` argument controls which weakness criteria are
applied. Multiple criteria can be combined, for example
`c("reconstruction", "similarity")`. At least one criterion must be
specified:

| Criterion                    | Individual-level score                                               | Requires `Yu`?   |
|:-----------------------------|:---------------------------------------------------------------------|:-----------------|
| `"reconstruction"` (default) | Error of reconstructing $\mathbf{X}_{u}$ from the PLS latent space   | No               |
| `"similarity"` (default)     | Mahalanobis distance between the target and training score centroids | No               |
| `"response"`                 | RMSE of predicting $\mathbf{y}_{u}$                                  | Yes              |
| `"range"`                    | Penalty for predictions outside `Yu_lims`                            | No (bounds only) |

### 3.4 Gene silencing

After weakness scores are computed, genes are retained or silenced based
on the `retain` parameter and the retention strategy specified in
[`gesearch_control()`](https://l-ramirez-lopez.github.io/resemble/reference/gesearch_control.md):

- **Probability-based** (`retain_by = "probability"`, default): Genes
  with weakness scores below the `retain` quantile are kept. This
  approach is more robust to outliers.
- **Proportion-based** (`retain_by = "proportion"`): A fixed proportion
  of genes with the lowest weakness scores is kept.

Silenced genes are permanently excluded:

$$G^{(\gamma + 1)} = G^{(\gamma)}\backslash U^{(\gamma)}$$

### 3.5 Final model

Evolution terminates when the gene pool falls to `target_size` or when
no additional genes can be silenced (stagnation). A final PLS model is
then fitted using the surviving genes, with the optimal number of
components selected by cross-validation independently of the fixed
`ncomp` used during evolution.

## 4 Key parameters

| Parameter     | Description                                                                                               |
|:--------------|:----------------------------------------------------------------------------------------------------------|
| `k`           | Number of active genes per individual                                                                     |
| `b`           | Gene-representation factor, that is, the target average number of individuals in which each gene appears  |
| `retain`      | Retention threshold, that is, the proportion of genes kept per generation (values \> 0.9 are recommended) |
| `target_size` | Minimum gene pool size at which evolution terminates                                                      |
| `ncomp`       | Number of PLS components fixed during evolution for comparability                                         |

The population size at generation $\gamma$ is approximately:

$$p^{(\gamma)} \approx \frac{b \cdot n^{(\gamma)}}{k}$$

where $n^{(\gamma)}$ is the number of active genes remaining.

## 5 Example workflow

The `NIRSoil` dataset is also used here to illustrate how to perform an
evolutionary subset search with
[`gesearch()`](https://l-ramirez-lopez.github.io/resemble/reference/gesearch.md).
The workflow includes defining the search, controlling the selection
process, inspecting the selected subsets, and generating predictions
from the final model.

It is important to note that the `NIRSoil` dataset may not be an ideal
test case for
[`gesearch()`](https://l-ramirez-lopez.github.io/resemble/reference/gesearch.md),
because of two reasons: *i*) the target set is a subset of the library
and the sample size is relatively small. *ii*) the test set, which could
be used as the target set, appears to be dominated by non-linear
relationships between spectra and response values, which may limit the
performance of the search. This is mainly because the dataset originates
from samples from a large geographical area, which may introduce complex
interactions between soil properties and spectral features.

### 5.1 The dataset, the target set, and its practical constraints

To better illustrate the main workflow and key parameters of the method,
the target set was extracted from the test set using a simple filtering
strategy intended to reduce potential non-linearities between spectra
and response values. Spectra that were substantially dissimilar to the
median spectrum of the test set were removed. Although this does not
guarantee the elimination of non-linear relationships, it provides a
straightforward way to define a more spectrally homogeneous target set
for demostration purposes.

``` r
library(prospectr)
data(NIRsoil)

# Preprocess spectra
NIRsoil$spc_pr <- savitzkyGolay(
  detrend(NIRsoil$spc, wav = as.numeric(colnames(NIRsoil$spc))),
  m = 1, p = 1, w = 7
)

# Missing values in the response are allowed
train_x <- NIRsoil$spc_pr[NIRsoil$train == 1, ]
train_y <- NIRsoil$Ciso[NIRsoil$train == 1]
test_x  <- NIRsoil$spc_pr[NIRsoil$train == 0, ]
test_y  <- NIRsoil$Ciso[NIRsoil$train == 0]

md <- dissimilarity(
  test_x,
  apply(test_x, 2, median),
  diss_method = diss_correlation(ws = 101, center = T, scale = T)
)

threshold <- 0.4
```

``` r
op <- par(mfrow = c(1, 2), mar = c(4.5, 4.5, 2, 1))

plot(
  seq_along(md$dissimilarity),
  md$dissimilarity, 
  ylim = c(0, 1), 
  col = "steelblue",
  pch = 16,
  cex = 1.5,
  xlab = "Sample index", 
  ylab = "Dissimilarity to median spectrum"
)
points(
  which(md$dissimilarity >= threshold),
  md$dissimilarity[md$dissimilarity >= threshold],
  col = "tomato",
  pch = 16,
  cex = 1.5
)
abline(h = threshold, col = "tomato", lty = 2, lwd = 2)
grid(lty = 1)

matplot(
  as.numeric(colnames(train_x)),
  t(test_x[md$dissimilarity >= threshold, ]),
  col = "tomato",
  ylab = "Preprocessed spectra",
  xlab = "Wavelengths, nm",
  lty = 1,
  type = "l"
)
grid(lty = 1)
matlines(
  as.numeric(colnames(train_x)),
  t(test_x[md$dissimilarity < threshold, ]),
  col = "steelblue",
  lty = 1,
  type = "l"
)

par(op)
```

![](evolutionary-subset-search_files/figure-html/fig-diss-threshold-1.png)

Figure 1: Left: Dissimilarity of the test spectra to the median test
spectrum. Samples above the threshold (dashed line and red) were
excluded from the target set as a simple attempt to reduce spectral
heterogeneity and potential non-linearities in the spectra-response
relationship. Right: Spectra of the samples above (red, to be excluded)
and below (blue, to be kept) the dissimilarity threshold.

To construct the final target set, samples whose dissimilarity to the
median spectrum falls below the correlation-dissimilarity threshold
(0.4) are retained. This set is then split into two subsets: a very
small subset, used to guide the search within
[`gesearch()`](https://l-ramirez-lopez.github.io/resemble/reference/gesearch.md),
and a second subset, used to validate the models built from the
reference samples selected by
[`gesearch()`](https://l-ramirez-lopez.github.io/resemble/reference/gesearch.md).
The use of only a very small subset to drive the search reflects a
typical practical constraint in target-domain applications, where only a
limited number of samples with laboratory reference measurements are
available within a short time frame. The “target” subsets are then
obtained as follows:

``` r
keep <- md$dissimilarity < threshold

cat("Number of samples retained:", sum(keep))
```

    Number of samples retained: 119

``` r
test_x_test <- test_x[keep, ]
test_y_test <- test_y[keep]

# sample a very small subset to guide the search, 
# and use the rest for testing the final model
set.seed(1124)
ref_ind <- sample(sum(keep), 8)
test_x_ref <- test_x_test[ref_ind, ]
test_y_ref <- test_y_test[ref_ind]

test_x_test <- test_x_test[-ref_ind, ]
test_y_test <- test_y_test[-ref_ind]
```

### 5.2 A simple linear model for predicting in the target set (baseline)

Here, the full training set is combined with the small subset of target
samples to fit a simple PLS model. This serves as a baseline for
comparison with the models generated by
[`gesearch()`](https://l-ramirez-lopez.github.io/resemble/reference/gesearch.md).
The optimal number of components is selected by cross-validation, and
the resulting model is evaluated on the test subset of the target set.

``` r
set.seed(1124)

pls_model <- model(
  Xr = rbind(
    train_x[!is.na(train_y), ],
    test_x_ref[!is.na(test_y_ref), ]
  ), 
  Yr = c(
    train_y[!is.na(train_y)],
    test_y_ref[!is.na(test_y_ref)]
  ),
  fit_method = fit_pls(ncomp = 15, method = "mpls", scale = TRUE),
  control = model_control(validation_type = "lgo", number = 10)
)
```

    Running cross-validation...

    Fitting model...

``` r
best_ncomp <- which(pls_model$cv_results$optimal)

pred <- predict(pls_model, test_x_test, ncomp = best_ncomp)
```

A quick function for computing the $R^{2}$ and root mean square error
($RMSE$) of the predictions is defined here for convenience:

``` r
reg_metrics <- function(obs, pred, na.rm = TRUE) {
  if (na.rm) {
    ok <- complete.cases(obs, pred)
    obs <- obs[ok]
    pred <- pred[ok]
  }
  c(RMSE = sqrt(mean((obs - pred)^2)), R2 = cor(obs, pred)^2)
}
```

The $RMSE$ and $R^{2}$ of the globla PLS model:

``` r
reg_metrics(test_y_test, pred)
```

         RMSE        R2
    0.4570417 0.7675839 

### 5.3 The model built with the samples found by *gesearch*

The
[`gesearch()`](https://l-ramirez-lopez.github.io/resemble/reference/gesearch.md)
function is then used to identify a subset of samples from the library
that is most relevant to the target set. In the following example, the
search is guided by three optimisation criteria: reconstruction error,
similarity in latent space, and response prediction error. The
parameters used are `k = 50`, meaning that each individual contains 50
active genes or samples; `b = 100`, meaning that each gene is expected
to appear in approximately 100 individuals; and `retain = 0.98`, meaning
that, for each weakness metric, samples with weakness scores above the
98th quantile are silenced. The search terminates when the gene pool
shrinks to `target_size = 80` samples. During evolution, a PLS model
with a fixed number of components (`ncomp = 15`) is used for
comparability across individuals, but the optimal number of components
for the final model is selected by cross-validation independently of
this fixed value.

Using the
[`gesearch_control()`](https://l-ramirez-lopez.github.io/resemble/reference/gesearch_control.md)
additional aspects of the search can be controlled, such as the strategy
for retaining genes based on their weakness scores. In this example, a
probability-based retention strategy is used, where genes with weakness
scores below the `retain` quantile are kept. This approach is more
robust to outliers compared to a fixed proportion-based strategy. Here,
tuning the number of factors for each individual is also disbaled for
comparability across individuals, but it can be enabled by setting
`tune = TRUE` in the control settings. The way the final models (and
also the intermediate models during evolution in case `tune = TRUE`) is
controlled by the parameters `number` and `p`, which specify the number
of groups and the proportion of observations per group for
leave-group-out cross-validation, respectively.

**Note** that `tune = FALSE` is recommended for most applications (and
it is the default), as it enables a more direct comparison of the
performance of different subsets during evolution. This is especially
advisable when the search is guided by the `"reconstruction"` and/or
`"similarity"` criteria. The reason is that tuning targets only the
minimization of the RMSE of the response variable. As a result, models
with different numbers of components may be selected during evolution.
Models with too few components may reconstruct the target set poorly,
which can in turn bias the weakness scores and, consequently, the
selection of samples. This may also artificially increase the apparent
similarity of some subsets to the target set in latent space. In
additon, when `tune = FALSE`, the search process is considerable faster
than if that is set to `TRUE`.

``` r
my_control <- gesearch_control(
  retain_by = "probability",
  tune = FALSE,
  number = 10L,
  p = 0.75
)
```

``` r
# Basic search with reconstruction optimization
gs <- gesearch(
  Xr = train_x[!is.na(train_y), ], # the spectra of the reference "set/library"
  Yr = train_y[!is.na(train_y)],   # the response of the reference "set/library"
  Xu = test_x_ref, # the available target domain spectra
  Yu = test_y_ref, # the available target domain response (accepts NAs)
  k = 50, b = 100, retain = 0.98,
  target_size = 80,
  fit_method = fit_pls(ncomp = 20, scale = TRUE),
  optimization = c("reconstruction", "similarity", "response"),
  control = my_control,
  seed = 1124
)
```

[Figure 2](#fig-evolution-gs) presents the evolution of the weakness
scores through the generations of the search.

``` r
plot(gs)
```

![](evolutionary-subset-search_files/figure-html/fig-evolution-gs-1.png)

Figure 2: Evolution of the three weakness metrics (response,
reconstruction and dissimilarity) used in the search.

The final model is validated in two ways. First, it is internally
validated on the selected samples using leave-group-out
cross-validation. Second, a model fitted only on the selected samples is
applied to the target samples that drove the search and for which
response values are available (`Yu`).

The validation results can be accessed from the final object. For the
internal leave-group-out cross-validation (CV) carried out on the
training set, the following summary metrics are provided: the mean RMSE
across CV groups (`rmse_cv`), the standardised RMSE (`st_rmse_cv`),
computed as the RMSE divided by the range of the response values
(maximum minus minimum), the standard deviation of the RMSE across CV
groups (`rmse_sd_cv`), and the mean (R^2) across CV groups (`r2_cv`).

``` r
# the internal leave-group out CV results for the data found
gs$validation_results[[1]]$results$train
```

       ncomp   rmse_cv st_rmse_cv rmse_sd_cv     r2_cv
    1      1 0.7221250 0.12680817 0.40804619 0.6503387
    2      2 0.7163662 0.12826806 0.38718606 0.7121203
    3      3 0.6886790 0.12312154 0.37755168 0.7352152
    4      4 0.6157296 0.10955431 0.34994044 0.7684410
    5      5 0.5364958 0.09483074 0.30778869 0.8204205
    6      6 0.4943293 0.08839862 0.26269455 0.8483245
    7      7 0.4634445 0.08488671 0.21472527 0.8613798
    8      8 0.4128501 0.07728087 0.16143597 0.8837915
    9      9 0.3750232 0.07237421 0.11441093 0.8954057
    10    10 0.3327559 0.06658169 0.07281027 0.9070918
    11    11 0.3160122 0.06347594 0.07008877 0.9137224
    12    12 0.2938045 0.05974065 0.06885816 0.9219629
    13    13 0.3161099 0.06444232 0.06818090 0.9111418
    14    14 0.3229156 0.06598235 0.06264075 0.9060152
    15    15 0.3287992 0.06688175 0.06974288 0.9046427
    16    16 0.3487916 0.06975661 0.09623761 0.9003171
    17    17 0.3584635 0.07119997 0.10093481 0.9003744
    18    18 0.3783390 0.07457618 0.11573189 0.8958179
    19    19 0.3917115 0.07705482 0.11449587 0.8887094
    20    20 0.3963061 0.07710436 0.12276786 0.8890305

``` r
# the CV on the target samples that drove the search (note that this is not a proper test set, as these samples were used to guide the search, but it can provide some insight into the performance of the final model on the target domain)
gs$validation_results[[1]]$results$test
```

       ncomp        r2       rmse         me
    1      1 0.7894760 0.28058875 0.13604327
    2      2 0.8816770 0.24483076 0.15165217
    3      3 0.7339895 0.28523730 0.08009277
    4      4 0.8393261 0.23616530 0.03298316
    5      5 0.8246959 0.23559046 0.08508896
    6      6 0.8597384 0.23692129 0.11163843
    7      7 0.9315969 0.16365480 0.08461275
    8      8 0.8662769 0.20784713 0.09851803
    9      9 0.8972302 0.18693274 0.06764719
    10    10 0.8827052 0.18761602 0.05044255
    11    11 0.9052043 0.17175485 0.02959933
    12    12 0.9149548 0.17115238 0.04752007
    13    13 0.8953919 0.17617524 0.04193763
    14    14 0.9200482 0.15702861 0.04448261
    15    15 0.9317147 0.13709674 0.03730481
    16    16 0.9552837 0.11041106 0.01228111
    17    17 0.9648178 0.10832365 0.02386677
    18    18 0.9656547 0.10175148 0.03941513
    19    19 0.9721626 0.10535369 0.03908109
    20    20 0.9769381 0.09696647 0.03315576

Assume that the final model is chosen as the one with the lowest CV RMSE
obtained from the samples found:

``` r
best_ncomp <- which.min(gs$validation_results[[1]]$results$train$rmse_cv)
```

The predictions on the samples that represent the target set and that
did not participate in the search process:

``` r
pred_gs <- predict(gs, test_x_test)[[1]]

reg_metrics(test_y_test, pred_gs[, best_ncomp])
```

         RMSE        R2
    0.2572124 0.8190477 

[Figure 3](#fig-pred-comparison) compares the observed versus predicted
values for the model fitted using all available training samples and the
model fitted using the subset of samples selected by
[`gesearch()`](https://l-ramirez-lopez.github.io/resemble/reference/gesearch.md).

``` r
op <- par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))

rng <- range(pred, pred_gs, test_y_test, na.rm = TRUE)

plot(
  pred, test_y_test,
  xlim = rng,
  ylim = rng,
  xlab = "Predicted Total Carbon, %",
  ylab = "Observed Total Carbon, %",
  main = "Model with all training \nsamples",
  cex = 1.5,
  pch = 16,
  col = rgb(0.5, 0.5, 0.5, 0.6)
)
grid(lty = 1)
abline(0, 1, col = "red")

plot(
  pred_gs[, best_ncomp], test_y_test,
  xlim = rng,
  ylim = rng,
  xlab = "Predicted Total Carbon, %",
  ylab = "Observed Total Carbon, %",
  main = "Model with gesearch-selected \nsamples",
  cex = 1.5,
  pch = 16,
  col = rgb(0.5, 0.5, 0.5, 0.6)
)
grid(lty = 1)
abline(0, 1, col = "red")

par(op)
```

![](evolutionary-subset-search_files/figure-html/fig-pred-comparison-1.png)

Figure 3: Comparison of observed versus predicted total carbon values.
Left: model fitted using all available training samples. Right: model
fitted using the subset of samples selected by gesearch().

### 5.4 Further examples of *gesearch*

#### 5.4.1 No response values available for the target set

The search process can be conducted when no information on the response
variable is available for the target set:

``` r
# Basic search with reconstruction optimization
gs_nr <- gesearch(
  Xr = train_x[!is.na(train_y), ], # the spectra of the reference "set/library"
  Yr = train_y[!is.na(train_y)],   # the response of the reference "set/library"
  Xu = test_x_ref, # the available target domain spectra
  Yu = NULL, # NO available target domain response (accepts NAs)
  k = 50, b = 100, retain = 0.98,
  target_size = 120,
  fit_method = fit_pls(ncomp = 20, scale = TRUE),
  optimization = c("reconstruction", "similarity"),
  control = my_control,
  seed = 1124
)
```

``` r
best_ncomp_nr <- which.min(gs_nr$validation_results[[1]]$results$train$rmse_cv)
best_ncomp_nr
```

    [1] 7

``` r
pred_gs_nr <- predict(gs_nr, test_x_test)[[1]]

reg_metrics(test_y_test, pred_gs_nr[, best_ncomp_nr])
```

         RMSE        R2
    0.3926205 0.7136982 

#### 5.4.2 No response values in the target set but with information about the response range

In some cases, the response values for the target set may not be
available, but the range of possible response values may be known. This
information can be used to guide the search by penalizing predictions
that fall outside the specified range. The `Yu_lims` argument can be
used for this purpose, and the `"range"` criterion must be included in
the `optimization` argument. The `"range"` option penalizes samples that
tend to generate models that predict values for the target spectra that
fall outside the specified limits.

``` r
# Basic search with reconstruction optimization
gs_nr2 <- gesearch(
  Xr = train_x[!is.na(train_y), ], # the spectra of the reference "set/library"
  Yr = train_y[!is.na(train_y)],   # the response of the reference "set/library"
  Xu = test_x_ref, # the available target domain spectra
  Yu = NULL, # NO available target domain response (accepts NAs)
  Yu_lims = c(0, 5), # the response range in the target domain
  k = 50, b = 100, retain = 0.98,
  target_size = 120,
  fit_method = fit_pls(ncomp = 20, scale = TRUE),
  optimization = c("reconstruction", "similarity", "range"),
  control = my_control,
  seed = 1124
)
```

``` r
best_ncomp_nr2 <- which.min(gs_nr2$validation_results[[1]]$results$train$rmse_cv)
best_ncomp_nr2
```

    [1] 8

``` r
pred_gs_nr2 <- predict(gs_nr2, test_x_test)[[1]]

reg_metrics(test_y_test, pred_gs_nr2[, best_ncomp_nr2])
```

         RMSE        R2
    0.2828443 0.7968024 

#### 5.4.3 Using SIMPLS for speeding up the search

The
[`fit_pls()`](https://l-ramirez-lopez.github.io/resemble/reference/fit_methods.md)
function allows the use of the SIMPLS algorithm for fitting PLS models
during the search. SIMPLS is a computationally efficient algorithm for
fitting PLS models, which can speed up the search process, especially
when dealing with very large datasets or when a large number of
individuals need to be evaluated per generation. The simpls algorithm
works by directly computing the PLS components without the need for
iterative deflation, which can significantly reduce computation time
while still providing accurate estimates of the PLS components. The
results obtained with simpls will differ from those obatained with the
default method (modified PLS, `"mpls"`) as the PLS weights are computed
differently. To use SIMPLS, simply set the `method` argument to
`"simpls"` in the
[`fit_pls()`](https://l-ramirez-lopez.github.io/resemble/reference/fit_methods.md)
function:

``` r
gs_simpls <- gesearch(
  Xr = train_x[!is.na(train_y), ], # the spectra of the reference "set/library"
  Yr = train_y[!is.na(train_y)],   # the response of the reference "set/library"
  Xu = test_x_ref, # the available target domain spectra
  Yu = test_y_ref, # the available target domain response (accepts NAs)
  k = 50, b = 100, retain = 0.98,
  target_size = 80,
  fit_method = fit_pls(ncomp = 20, method = "simpls", scale = TRUE),
  optimization = c("reconstruction", "similarity", "response"),
  control = my_control,
  seed = 1124
)
```

``` r
best_ncomp_simpls <- which.min(gs_simpls$validation_results[[1]]$results$train$rmse_cv)
best_ncomp_simpls
```

    [1] 15

``` r
pred_gs_simpls <- predict(gs_simpls, test_x_test)[[1]]

reg_metrics(test_y_test, pred_gs_simpls[, best_ncomp_simpls])
```

         RMSE        R2
    0.2660343 0.8172405 

## 6 Supported parallel processing

The
[`gesearch()`](https://l-ramirez-lopez.github.io/resemble/reference/gesearch.md)
function supports parallel execution for evaluating individuals during
the evolutionary search. To improve memory efficiency, the data are not
sent in full to each worker. Instead, only the subsets of samples
required for the individuals being evaluated are extracted and passed to
the parallel workers. The `pchunks` argument controls how many
individuals are processed together at each iterator step during parallel
execution. Its effect on speed and memory usage depends on the size of
the subsets and the available computing resources. By default,
`pchunks = 1L`, meaning that individuals are processed one at a time.
Increasing `pchunks` can reduce parallelisation overhead, but it may
also increase memory usage because more data may need to be sent to each
worker at once. The optimal value depends on the dataset and
computational environment, and may require some experimentation.

In the example below, the [doParallel
package](https://CRAN.R-project.org/package=doParallel) is used to
register the cores for parallel execution. The [doSNOW
package](https://CRAN.R-project.org/package=doSNOW) can also be used.
The example illustrates how to run
[`gesearch()`](https://l-ramirez-lopez.github.io/resemble/reference/gesearch.md)
using multiple cores.

``` r
# Running gesearch() using multiple cores

# Execute with two cores, if available, ...
n_cores <- 2

# ... if not then go with 1 core
if (parallel::detectCores() < 2) {
  n_cores <- 1
}

# Set the number of cores
library(doParallel)
clust <- makeCluster(n_cores)
registerDoParallel(clust)

# Alternatively:
# library(doSNOW)
# clust <- makeCluster(n_cores, type = "SOCK")
# registerDoSNOW(clust)
# getDoParWorkers()

my_control <- gesearch_control(
  retain_by = "probability",
  tune = FALSE,
  number = 10L,
  p = 0.75
)

gs_p <- gesearch(
  Xr = train_x[!is.na(train_y), ],
  Yr = train_y[!is.na(train_y)],
  Xu = test_x_ref,
  Yu = test_y_ref,
  k = 50,
  b = 100,
  retain = 0.98,
  target_size = 80,
  fit_method = fit_pls(ncomp = 20, scale = TRUE),
  optimization = c("reconstruction", "similarity", "response"),
  control = my_control,
  seed = 1124,
  pchunks = 2L
)

# Go back to sequential processing
registerDoSEQ()
try(stopCluster(clust))

gs_p
```

Parallel processing is particularly useful when the search requires
fitting a large number of intermediate models per generation. However,
the computational gain depends on the number of available cores, the
size of the population, and the cost of fitting each individual model.

## References

Holland, J.H., 1995. Hidden order: How adaptation builds complexity.
Addison-Wesley, Reading, MA.

Lobsey, C., Viscarra-Rossel, R., Roudier, P., Hedley, C., 2017. Rs-local
data-mines information from spectral libraries to improve local
calibrations. European Journal of Soil Science 68, 840–852.

Ramirez-Lopez, L., Viscarra Rossel, R., Behrens, T., Orellano, C.,
Perez-Fernandez, E., Kooijman, L., Wadoux, A.M.J.-C., Breure, T.,
Summerauer, L., Safanelli, J.L., Plans, M., 2026. When spectral
libraries are too complex to search: Evolutionary subset selection for
domain-adaptive calibration. Analytica Chimica Acta.
