# Package index

## Package

- [`resemble-package`](https://l-ramirez-lopez.github.io/resemble/reference/resemble-package.md)
  [`resemble`](https://l-ramirez-lopez.github.io/resemble/reference/resemble-package.md)
  : Overview of the functions in the resemble package

## Dissimilarity Computation

Compute dissimilarity matrices using various methods

- [`dissimilarity()`](https://l-ramirez-lopez.github.io/resemble/reference/dissimilarity.md)
  : Compute dissimilarity matrices
- [`diss_pca()`](https://l-ramirez-lopez.github.io/resemble/reference/diss_pca.md)
  : PCA dissimilarity method constructor
- [`diss_pls()`](https://l-ramirez-lopez.github.io/resemble/reference/diss_pls.md)
  : PLS dissimilarity method constructor
- [`diss_correlation()`](https://l-ramirez-lopez.github.io/resemble/reference/diss_correlation.md)
  : Correlation dissimilarity method constructor
- [`diss_euclidean()`](https://l-ramirez-lopez.github.io/resemble/reference/diss_euclidean.md)
  : Euclidean dissimilarity method constructor
- [`diss_mahalanobis()`](https://l-ramirez-lopez.github.io/resemble/reference/diss_mahalanobis.md)
  : Mahalanobis dissimilarity method constructor
- [`diss_cosine()`](https://l-ramirez-lopez.github.io/resemble/reference/diss_cosine.md)
  : Cosine dissimilarity method constructor
- [`diss_evaluate()`](https://l-ramirez-lopez.github.io/resemble/reference/diss_evaluate.md)
  [`sim_eval()`](https://l-ramirez-lopez.github.io/resemble/reference/diss_evaluate.md)
  **\[stable\]** : Evaluate dissimilarity matrices
- [`sid()`](https://l-ramirez-lopez.github.io/resemble/reference/sid.md)
  : A function for computing the spectral information divergence between
  spectra (sid)

## Neighbour Search

Find nearest neighbours in reference sets

- [`search_neighbors()`](https://l-ramirez-lopez.github.io/resemble/reference/search_neighbors.md)
  : Search neighbors in a reference set

- [`search_control()`](https://l-ramirez-lopez.github.io/resemble/reference/search_control.md)
  :

  A function that controls some aspects of the the `gesearch` function

- [`neighbors_k()`](https://l-ramirez-lopez.github.io/resemble/reference/neighbors.md)
  [`neighbors_diss()`](https://l-ramirez-lopez.github.io/resemble/reference/neighbors.md)
  : Neighbor selection methods

## Memory-Based Learning

Local modelling with per-query model fitting

- [`mbl()`](https://l-ramirez-lopez.github.io/resemble/reference/mbl.md)
  [`plot(`*`<mbl>`*`)`](https://l-ramirez-lopez.github.io/resemble/reference/mbl.md)
  [`get_predictions()`](https://l-ramirez-lopez.github.io/resemble/reference/mbl.md)
  : Memory-based learning (mbl)
- [`mbl_control()`](https://l-ramirez-lopez.github.io/resemble/reference/mbl_control.md)
  : Control parameters for memory-based learning

## Library of Local Experts

Pre-computed local models for retrieval-gated prediction

- [`liblex()`](https://l-ramirez-lopez.github.io/resemble/reference/liblex.md)
  [`predict(`*`<liblex>`*`)`](https://l-ramirez-lopez.github.io/resemble/reference/liblex.md)
  [`plot(`*`<liblex>`*`)`](https://l-ramirez-lopez.github.io/resemble/reference/liblex.md)
  : Build a precomputed library of localised experts using memory-based
  learning
- [`liblex_control()`](https://l-ramirez-lopez.github.io/resemble/reference/liblex_control.md)
  : Control parameters for liblex

## Evolutionary Subset Selection

Domain-adaptive sample selection from spectral libraries

- [`gesearch(`*`<default>`*`)`](https://l-ramirez-lopez.github.io/resemble/reference/gesearch.md)
  [`gesearch(`*`<formula>`*`)`](https://l-ramirez-lopez.github.io/resemble/reference/gesearch.md)
  [`predict(`*`<gesearch>`*`)`](https://l-ramirez-lopez.github.io/resemble/reference/gesearch.md)
  [`plot(`*`<gesearch>`*`)`](https://l-ramirez-lopez.github.io/resemble/reference/gesearch.md)
  : Evolutionary sample search for context-specific calibrations
- [`gesearch_control()`](https://l-ramirez-lopez.github.io/resemble/reference/gesearch_control.md)
  : Control parameters for gesearch

## Global Modelling

Global calibration models with cross-validation

- [`model()`](https://l-ramirez-lopez.github.io/resemble/reference/model.md)
  [`predict(`*`<resemble_model>`*`)`](https://l-ramirez-lopez.github.io/resemble/reference/model.md)
  : Global spectral calibration model
- [`model_control()`](https://l-ramirez-lopez.github.io/resemble/reference/model_control.md)
  : Control parameters for global model fitting

## Projection Methods

Dimensionality reduction via PCA and PLS

- [`ortho_projection()`](https://l-ramirez-lopez.github.io/resemble/reference/ortho_projection.md)
  [`predict(`*`<ortho_projection>`*`)`](https://l-ramirez-lopez.github.io/resemble/reference/ortho_projection.md)
  [`plot(`*`<ortho_projection>`*`)`](https://l-ramirez-lopez.github.io/resemble/reference/ortho_projection.md)
  : Orthogonal projections

## Fitting Methods

Constructors for local regression methods

- [`fit_pls()`](https://l-ramirez-lopez.github.io/resemble/reference/fit_methods.md)
  [`fit_wapls()`](https://l-ramirez-lopez.github.io/resemble/reference/fit_methods.md)
  [`fit_gpr()`](https://l-ramirez-lopez.github.io/resemble/reference/fit_methods.md)
  : Local fitting method constructors

## Component Selection

Methods for selecting the number of components

- [`ncomp_by_var()`](https://l-ramirez-lopez.github.io/resemble/reference/ncomp_selection.md)
  [`ncomp_by_cumvar()`](https://l-ramirez-lopez.github.io/resemble/reference/ncomp_selection.md)
  [`ncomp_by_opc()`](https://l-ramirez-lopez.github.io/resemble/reference/ncomp_selection.md)
  [`ncomp_fixed()`](https://l-ramirez-lopez.github.io/resemble/reference/ncomp_selection.md)
  : Component selection methods
