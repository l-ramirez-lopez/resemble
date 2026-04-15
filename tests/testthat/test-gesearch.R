context("test-gesearch")
library(foreach)
library(foreach)
library(RhpcBLASctl)
registerDoSEQ()
# =============================================================================
# Setup helper
# =============================================================================

.setup_gesearch_data <- function() {
  skip_if_not_installed("prospectr")
  
  data("NIRsoil", package = "prospectr")
  
  # Preprocess spectra
  sg_det <- prospectr::savitzkyGolay(
    prospectr::detrend(NIRsoil$spc, wav = as.numeric(colnames(NIRsoil$spc))),
    m = 1, p = 1, w = 7
  )
  NIRsoil$spc_pr <- sg_det
  
  # Split data
  train_idx <- NIRsoil$train == 1 & !is.na(NIRsoil$Ciso)
  test_idx <- NIRsoil$train == 0 & !is.na(NIRsoil$Ciso)
  
  list(
    train_x = NIRsoil$spc_pr[train_idx, ],
    train_y = NIRsoil$Ciso[train_idx],
    test_x = NIRsoil$spc_pr[test_idx, ],
    test_y = NIRsoil$Ciso[test_idx]
  )
}


# =============================================================================
# gesearch input validation tests
# =============================================================================

test_that("gesearch validates input dimensions", {
  skip_if_not_installed("prospectr")
  
  d <- .setup_gesearch_data()
  
  # nrow(Xr) must equal length(Yr)
  expect_error(
    gesearch(
      Xr = d$train_x,
      Yr = d$train_y[1:10],
      Xu = d$test_x,
      k = 30, b = 10, target_size = 350,
      fit_method = fit_pls(ncomp = 5)
    ),
    "nrow.*Xr.*must equal"
  )
  
  # ncol(Xr) must equal ncol(Xu)
  expect_error(
    gesearch(
      Xr = d$train_x,
      Yr = d$train_y,
      Xu = d$test_x[, 1:10],
      k = 30, b = 10, target_size = 350,
      fit_method = fit_pls(ncomp = 5)
    ),
    "ncol.*Xr.*must equal.*ncol.*Xu"
  )
})


test_that("gesearch validates k parameter", {
  skip_if_not_installed("prospectr")
  
  d <- .setup_gesearch_data()
  
  expect_error(
    gesearch(
      Xr = d$train_x,
      Yr = d$train_y,
      Xu = d$test_x,
      k = nrow(d$train_x) + 10,
      b = 10, target_size = 350,
      fit_method = fit_pls(ncomp = 5)
    ),
    "k.*must be less than"
  )
})


test_that("gesearch validates target_size parameter", {
  skip_if_not_installed("prospectr")
  
  d <- .setup_gesearch_data()
  
  # target_size must be >= k
  expect_error(
    gesearch(
      Xr = d$train_x,
      Yr = d$train_y,
      Xu = d$test_x,
      k = 50, b = 10, target_size = 30,
      fit_method = fit_pls(ncomp = 5)
    ),
    "target_size.*must be.*>= k"
  )
  
  # target_size cannot exceed nrow(Xr)
  expect_error(
    gesearch(
      Xr = d$train_x,
      Yr = d$train_y,
      Xu = d$test_x,
      k = 30, b = 10, target_size = nrow(d$train_x) + 100,
      fit_method = fit_pls(ncomp = 5)
    ),
    "target_size.*cannot be greater"
  )
})


test_that("gesearch validates retain parameter", {
  skip_if_not_installed("prospectr")
  
  d <- .setup_gesearch_data()
  
  expect_error(
    gesearch(
      Xr = d$train_x,
      Yr = d$train_y,
      Xu = d$test_x,
      k = 30, b = 10, retain = 0, target_size = 350,
      fit_method = fit_pls(ncomp = 5)
    ),
    "retain.*must be in"
  )
  
  expect_error(
    gesearch(
      Xr = d$train_x,
      Yr = d$train_y,
      Xu = d$test_x,
      k = 30, b = 10, retain = 1.5, target_size = 350,
      fit_method = fit_pls(ncomp = 5)
    ),
    "retain.*must be in"
  )
})


test_that("gesearch validates optimization parameter", {
  skip_if_not_installed("prospectr")
  
  d <- .setup_gesearch_data()
  
  expect_error(
    gesearch(
      Xr = d$train_x,
      Yr = d$train_y,
      Xu = d$test_x,
      k = 30, b = 10, target_size = 350,
      fit_method = fit_pls(ncomp = 5),
      optimization = "invalid_option"
    ),
    "Invalid.*optimization"
  )
})


test_that("gesearch requires Yu for response optimization", {
  skip_if_not_installed("prospectr")
  
  d <- .setup_gesearch_data()
  
  expect_error(
    gesearch(
      Xr = d$train_x,
      Yr = d$train_y,
      Xu = d$test_x,
      Yu = NULL,
      k = 30, b = 10, target_size = 350,
      fit_method = fit_pls(ncomp = 5),
      optimization = "response"
    ),
    "Yu.*required.*response"
  )
})


test_that("gesearch validates fit_method", {
  skip_if_not_installed("prospectr")
  
  d <- .setup_gesearch_data()
  
  expect_error(
    gesearch(
      Xr = d$train_x,
      Yr = d$train_y,
      Xu = d$test_x,
      k = 30, b = 10, target_size = 350,
      fit_method = "pls"
    ),
    "fit_method.*must be"
  )
})


test_that("gesearch rejects missing values in Xu", {
  skip_if_not_installed("prospectr")
  
  d <- .setup_gesearch_data()
  
  test_x_na <- d$test_x
  test_x_na[1, 1] <- NA
  
  expect_error(
    gesearch(
      Xr = d$train_x,
      Yr = d$train_y,
      Xu = test_x_na,
      k = 30, b = 10, target_size = 350,
      fit_method = fit_pls(ncomp = 5)
    ),
    "Xu.*missing values"
  )
})


test_that("gesearch rejects infinite values", {
  skip_if_not_installed("prospectr")
  
  d <- .setup_gesearch_data()
  
  test_x_inf <- d$test_x
  test_x_inf[1, 1] <- Inf
  
  expect_error(
    gesearch(
      Xr = d$train_x,
      Yr = d$train_y,
      Xu = test_x_inf,
      k = 30, b = 10, target_size = 350,
      fit_method = fit_pls(ncomp = 5)
    ),
    "Infinite values"
  )
})


test_that("gesearch validates Xu has at least 3 rows", {
  skip_if_not_installed("prospectr")
  
  d <- .setup_gesearch_data()
  
  expect_error(
    gesearch(
      Xr = d$train_x,
      Yr = d$train_y,
      Xu = d$test_x[1:2, ],
      k = 30, b = 10, target_size = 350,
      fit_method = fit_pls(ncomp = 5)
    ),
    "Xu must have at least 3 rows"
  )
})


test_that("gesearch validates verbose parameter", {
  skip_if_not_installed("prospectr")
  
  d <- .setup_gesearch_data()
  
  expect_error(
    gesearch(
      Xr = d$train_x,
      Yr = d$train_y,
      Xu = d$test_x,
      k = 30, b = 10, target_size = 350,
      fit_method = fit_pls(ncomp = 5),
      verbose = "yes"
    ),
    "verbose.*must be logical"
  )
})


# =============================================================================
# Basic functionality tests (skip on CRAN)
# =============================================================================

test_that("gesearch runs with reconstruction optimization", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_gesearch_data()
  
  gs <- gesearch(
    Xr = d$train_x,
    Yr = d$train_y,
    Xu = d$test_x,
    k = 30, b = 20, retain = 0.90,
    target_size = 350,
    fit_method = fit_pls(ncomp = 10),
    optimization = "reconstruction",
    verbose = FALSE,
    seed = 42
  )
  
  expect_s3_class(gs, "gesearch")
  expect_true("x_local" %in% names(gs))
  expect_true("y_local" %in% names(gs))
  expect_true("indices" %in% names(gs))
  expect_true("final_models" %in% names(gs))
  expect_true("validation_results" %in% names(gs))
  expect_gte(length(gs$indices), 50)
})


test_that("gesearch runs with different PLS methods", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_gesearch_data()
  
  # mpls
  gs_mpls <- gesearch(
    Xr = d$train_x,
    Yr = d$train_y,
    Xu = d$test_x,
    k = 30, b = 20, retain = 0.90,
    target_size = 350,
    fit_method = fit_pls(ncomp = 10, method = "mpls"),
    optimization = "reconstruction",
    verbose = FALSE,
    seed = 42
  )
  expect_s3_class(gs_mpls, "gesearch")
  
  # simpls
  gs_simpls <- gesearch(
    Xr = d$train_x,
    Yr = d$train_y,
    Xu = d$test_x,
    k = 30, b = 20, retain = 0.90,
    target_size = 350,
    fit_method = fit_pls(ncomp = 10, method = "simpls"),
    optimization = "reconstruction",
    verbose = FALSE,
    seed = 42
  )
  expect_s3_class(gs_simpls, "gesearch")
})


test_that("gesearch runs with response optimization", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_gesearch_data()
  
  gs <- gesearch(
    Xr = d$train_x,
    Yr = d$train_y,
    Xu = d$test_x,
    Yu = d$test_y,
    k = 30, b = 20, retain = 0.90,
    target_size = 350,
    fit_method = fit_pls(ncomp = 10),
    optimization = "response",
    verbose = FALSE,
    seed = 42
  )
  
  expect_s3_class(gs, "gesearch")
  expect_gte(length(gs$indices), 50)
})


test_that("gesearch runs with similarity optimization", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_gesearch_data()
  
  gs <- gesearch(
    Xr = d$train_x,
    Yr = d$train_y,
    Xu = d$test_x,
    k = 30, b = 20, retain = 0.90,
    target_size = 350,
    fit_method = fit_pls(ncomp = 10),
    optimization = "similarity",
    verbose = FALSE,
    seed = 42
  )
  expect_s3_class(gs, "gesearch")
})


test_that("gesearch runs with combined optimization criteria", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_gesearch_data()
  
  gs <- gesearch(
    Xr = d$train_x,
    Yr = d$train_y,
    Xu = d$test_x,
    Yu = d$test_y,
    k = 30, b = 20, retain = 0.90,
    target_size = 350,
    fit_method = fit_pls(ncomp = 10),
    optimization = c("reconstruction", "similarity"),
    verbose = FALSE,
    seed = 42
  )
  
  expect_s3_class(gs, "gesearch")
})


test_that("gesearch with scale = TRUE works and predict works too", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_gesearch_data()
  
  gs <- gesearch(
    Xr = d$train_x,
    Yr = d$train_y,
    Xu = d$test_x,
    k = 30, b = 20, retain = 0.90,
    target_size = 350,
    fit_method = fit_pls(ncomp = 10, scale = TRUE),
    optimization = "reconstruction",
    verbose = FALSE,
    seed = 42
  )
  
  expect_s3_class(gs, "gesearch")
})


# =============================================================================
# predict.gesearch tests
# =============================================================================



test_that("predict.gesearch validations", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_gesearch_data()
  
  gs <- gesearch(
    Xr = d$train_x,
    Yr = d$train_y,
    Xu = d$test_x,
    k = 30, b = 20, retain = 0.90,
    target_size = 350,
    fit_method = fit_pls(ncomp = 10),
    optimization = "reconstruction",
    verbose = FALSE,
    seed = 42
  )
  
  expect_error(
    predict(gs, newdata = d$test_x, type = "invalid"),
    "type.*must be"
  )
  
  preds <- predict(gs, newdata = d$test_x)
  
  expect_type(preds, "list")
  expect_equal(nrow(preds[[1]]), nrow(d$test_x))
  expect_equal(ncol(preds[[1]]), 10)  # ncomp columns
  expect_error(predict(gs), "newdata.*required")
  
  
})


test_that("predict.gesearch with what = 'all_generations' works", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_gesearch_data()
  
  gs <- gesearch(
    Xr = d$train_x,
    Yr = d$train_y,
    Xu = d$test_x,
    Yu = d$test_y,
    k = 30, b = 20, retain = 0.90,
    target_size = 350,
    fit_method = fit_pls(ncomp = 10),
    optimization = "reconstruction",
    intermediate_models = TRUE,
    verbose = FALSE,
    seed = 42
  )
  
  preds <- predict(gs, newdata = d$test_x, what = "all_generations")
  
  expect_type(preds, "list")
  expect_true(length(preds) == gs$complete_iter)
  expect_true(all(grepl("^generation_", names(preds))))
})


test_that("predict.gesearch warns without intermediate_models", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_gesearch_data()
  
  gs <- gesearch(
    Xr = d$train_x,
    Yr = d$train_y,
    Xu = d$test_x,
    k = 30, b = 20, retain = 0.90,
    target_size = 350,
    fit_method = fit_pls(ncomp = 10),
    optimization = "reconstruction",
    intermediate_models = FALSE,
    verbose = FALSE,
    seed = 42
  )
  
  expect_warning(
    predict(gs, newdata = d$test_x, what = "all_generations"),
    "intermediate models"
  )
})


# =============================================================================
# plot.gesearch tests
# =============================================================================

test_that("plot.gesearch works", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_gesearch_data()
  
  gs <- gesearch(
    Xr = d$train_x,
    Yr = d$train_y,
    Xu = d$test_x,
    k = 30, b = 20, retain = 0.90,
    target_size = 350,
    fit_method = fit_pls(ncomp = 10),
    optimization = "reconstruction",
    verbose = FALSE,
    seed = 42
  )
  
  expect_silent(plot(gs, which = "weakness"))
  expect_silent(plot(gs, which = "removed"))
})


# =============================================================================
# Reproducibility tests
# =============================================================================

test_that("gesearch is reproducible with seed", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_gesearch_data()
  
  gs1 <- gesearch(
    Xr = d$train_x,
    Yr = d$train_y,
    Xu = d$test_x,
    k = 30, b = 20, retain = 0.90,
    target_size = 350,
    fit_method = fit_pls(ncomp = 10),
    optimization = "reconstruction",
    verbose = FALSE,
    seed = 42
  )
  
  gs2 <- gesearch(
    Xr = d$train_x,
    Yr = d$train_y,
    Xu = d$test_x,
    k = 30, b = 20, retain = 0.90,
    target_size = 350,
    fit_method = fit_pls(ncomp = 10),
    optimization = "reconstruction",
    verbose = FALSE,
    seed = 42
  )
  
  expect_equal(gs1$indices, gs2$indices)
  expect_equal(gs1$complete_iter, gs2$complete_iter)
})


# =============================================================================
# Output structure tests
# =============================================================================

test_that("gesearch output contains expected elements", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_gesearch_data()
  
  gs <- gesearch(
    Xr = d$train_x,
    Yr = d$train_y,
    Xu = d$test_x,
    Yu = d$test_y,
    k = 30, b = 20, retain = 0.90,
    target_size = 350,
    fit_method = fit_pls(ncomp = 10),
    optimization = "reconstruction",
    verbose = FALSE,
    seed = 42
  )
  
  expected_elements <- c(
    "x_local", "y_local", "indices", "complete_iter",
    "iter_weakness", "samples", "n_removed", "control",
    "fit_method", "validation_results", "final_models", "seed"
  )
  
  for (elem in expected_elements) {
    expect_true(elem %in% names(gs), info = paste("Missing element:", elem))
  }
  
  # Check dimensions
  expect_equal(nrow(gs$x_local), length(gs$indices))
  expect_equal(nrow(gs$y_local), length(gs$indices))
  
  # Check validation results structure
  expect_true("results" %in% names(gs$validation_results[[1]]))
  expect_true("train" %in% names(gs$validation_results[[1]]$results))
  expect_true("test" %in% names(gs$validation_results[[1]]$results))
})


test_that("gesearch n_removed is a data.frame with correct structure", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_gesearch_data()
  
  gs <- gesearch(
    Xr = d$train_x,
    Yr = d$train_y,
    Xu = d$test_x,
    k = 30, b = 20, retain = 0.90,
    target_size = 350,
    fit_method = fit_pls(ncomp = 10),
    optimization = "reconstruction",
    verbose = FALSE,
    seed = 42
  )
  
  expect_s3_class(gs$n_removed, "data.frame")
  expect_true(all(c("iteration", "removed", "cumulative") %in% names(gs$n_removed)))
  expect_equal(nrow(gs$n_removed), gs$complete_iter)
  
  # Cumulative should be non-decreasing
  expect_true(all(diff(gs$n_removed$cumulative) >= 0))
})


test_that("gesearch with intermediate_models stores generations", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_gesearch_data()
  
  gs <- gesearch(
    Xr = d$train_x,
    Yr = d$train_y,
    Xu = d$test_x,
    Yu = d$test_y,
    k = 30, b = 20, retain = 0.90,
    target_size = 350,
    fit_method = fit_pls(ncomp = 10),
    optimization = "reconstruction",
    intermediate_models = TRUE,
    verbose = FALSE,
    seed = 42
  )
  
  expect_true("intermediate_models" %in% names(gs))
  expect_type(gs$intermediate_models, "list")
  expect_gt(length(gs$intermediate_models), 0)
  
  # Check structure of intermediate models
  for (gen in gs$intermediate_models) {
    expect_true("subset_size" %in% names(gen))
    expect_true("pls_models" %in% names(gen))
    expect_true("validation" %in% names(gen))
  }
})















test_that("gesearch produces similar results with simpls and pls for optimization = 'response'", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_gesearch_data()
  
  gs_simpls <- gesearch(
    Xr = d$train_x,
    Yr = d$train_y,
    Xu = d$test_x,
    Yu = d$test_y,
    k = 30, b = 20, retain = 0.90,
    target_size = 350,
    fit_method = fit_pls(ncomp = 10, method = "simpls"),
    optimization = "response",
    intermediate_models = TRUE,
    verbose = FALSE,
    seed = 42
  )
  
  gs_pls <- gesearch(
    Xr = d$train_x,
    Yr = d$train_y,
    Xu = d$test_x,
    Yu = d$test_y,
    k = 30, b = 20, retain = 0.90,
    target_size = 350,
    fit_method = fit_pls(ncomp = 10, method = "pls"),
    optimization = "response",
    intermediate_models = TRUE,
    verbose = FALSE,
    seed = 42
  )
  
  expect_true(identical(gs_simpls$indices, gs_pls$indices))
})



test_that("gesearch produces different results with pls and mpls for optimization = 'response'", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_gesearch_data()
  
  gs_pls <- gesearch(
    Xr = d$train_x,
    Yr = d$train_y,
    Xu = d$test_x,
    Yu = d$test_y,
    k = 30, b = 20, retain = 0.90,
    target_size = 350,
    fit_method = fit_pls(ncomp = 10, method = "pls"),
    optimization = "response",
    intermediate_models = TRUE,
    verbose = FALSE,
    seed = 42
  )
  
  gs_mpls <- gesearch(
    Xr = d$train_x,
    Yr = d$train_y,
    Xu = d$test_x,
    Yu = d$test_y,
    k = 30, b = 20, retain = 0.90,
    target_size = 350,
    fit_method = fit_pls(ncomp = 10, method = "mpls"),
    optimization = "response",
    intermediate_models = TRUE,
    verbose = FALSE,
    seed = 42
  )
  
  expect_false(identical(gs_mpls$indices, gs_pls$indices))
})
