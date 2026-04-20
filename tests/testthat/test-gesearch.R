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















































































































# =============================================================================
# gesearch_control validation tests
# =============================================================================

test_that("gesearch validates control class", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_gesearch_data()
  
  expect_error(
    gesearch(
      Xr = d$train_x,
      Yr = d$train_y,
      Xu = d$test_x,
      k = 30, b = 10, target_size = 350,
      fit_method = fit_pls(ncomp = 5),
      control = list(tune = FALSE)
    ),
    "'control' must be created by gesearch_control\\(\\)"
  )
})

# =============================================================================
# Multi-column Yr validation tests
# =============================================================================

test_that("gesearch validates multi-column Yr has column names", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_gesearch_data()
  
  # Multi-column Yr without column names
  Yr_multi <- cbind(d$train_y, d$train_y * 2)
  
  expect_error(
    gesearch(
      Xr = d$train_x,
      Yr = Yr_multi,
      Xu = d$test_x,
      k = 30, b = 10, target_size = 350,
      fit_method = fit_pls(ncomp = 5),
      verbose = FALSE
    ),
    "Yr must have column names when it has multiple columns"
  )
})

test_that("gesearch validates multi-column Yr row count", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_gesearch_data()
  
  # Multi-column Yr with wrong row count
  Yr_multi <- cbind(y1 = d$train_y[1:10], y2 = d$train_y[1:10])
  
  expect_error(
    gesearch(
      Xr = d$train_x,
      Yr = Yr_multi,
      Xu = d$test_x,
      k = 30, b = 10, target_size = 350,
      fit_method = fit_pls(ncomp = 5),
      verbose = FALSE
    ),
    "nrow\\(Xr\\) must equal nrow\\(Yr\\)"
  )
})

# =============================================================================
# Multi-column Yu validation tests
# =============================================================================

test_that("gesearch validates Yu length matches Xu rows", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_gesearch_data()
  
  expect_error(
    gesearch(
      Xr = d$train_x,
      Yr = d$train_y,
      Xu = d$test_x,
      Yu = d$test_y[1:5],
      k = 30, b = 10, target_size = 350,
      fit_method = fit_pls(ncomp = 5),
      verbose = FALSE
    ),
    "nrow\\(Xu\\) must equal length\\(Yu\\)"
  )
})

test_that("gesearch rejects Yu_lims with multi-column Yu", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_gesearch_data()
  
  Yu_multi <- cbind(y1 = d$test_y, y2 = d$test_y * 2)
  Yr_multi <- cbind(y1 = d$train_y, y2 = d$train_y * 2)
  
  expect_error(
    gesearch(
      Xr = d$train_x,
      Yr = Yr_multi,
      Xu = d$test_x,
      Yu = Yu_multi,
      Yu_lims = c(0, 10),
      k = 30, b = 10, target_size = 350,
      fit_method = fit_pls(ncomp = 5),
      verbose = FALSE
    ),
    "'Yu_lims' only supported for single-column Yu"
  )
})

test_that("gesearch validates multi-column Yu ncol matches Yr", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_gesearch_data()
  
  Yr_multi <- cbind(y1 = d$train_y, y2 = d$train_y * 2)
  Yu_wrong <- cbind(y1 = d$test_y, y2 = d$test_y * 2, y3 = d$test_y * 3)
  
  expect_error(
    gesearch(
      Xr = d$train_x,
      Yr = Yr_multi,
      Xu = d$test_x,
      Yu = Yu_wrong,
      k = 30, b = 10, target_size = 350,
      fit_method = fit_pls(ncomp = 5),
      verbose = FALSE
    ),
    "ncol\\(Yr\\) must equal ncol\\(Yu\\)"
  )
})

test_that("gesearch validates multi-column Yu nrow matches Xu", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_gesearch_data()
  
  Yr_multi <- cbind(y1 = d$train_y, y2 = d$train_y * 2)
  Yu_wrong <- cbind(y1 = d$test_y[1:5], y2 = d$test_y[1:5] * 2)
  
  expect_error(
    gesearch(
      Xr = d$train_x,
      Yr = Yr_multi,
      Xu = d$test_x,
      Yu = Yu_wrong,
      k = 30, b = 10, target_size = 350,
      fit_method = fit_pls(ncomp = 5),
      verbose = FALSE
    ),
    "nrow\\(Xu\\) must equal nrow\\(Yu\\)"
  )
})

test_that("gesearch validates multi-column Yu has column names", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_gesearch_data()
  
  Yr_multi <- cbind(y1 = d$train_y, y2 = d$train_y * 2)
  Yu_no_names <- cbind(d$test_y, d$test_y * 2)
  
  expect_error(
    gesearch(
      Xr = d$train_x,
      Yr = Yr_multi,
      Xu = d$test_x,
      Yu = Yu_no_names,
      k = 30, b = 10, target_size = 350,
      fit_method = fit_pls(ncomp = 5),
      verbose = FALSE
    ),
    "Yu must have column names when it has multiple columns"
  )
})

test_that("gesearch validates Yu and Yr column names match", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_gesearch_data()
  
  Yr_multi <- cbind(y1 = d$train_y, y2 = d$train_y * 2)
  Yu_wrong_names <- cbind(a = d$test_y, b = d$test_y * 2)
  
  expect_error(
    gesearch(
      Xr = d$train_x,
      Yr = Yr_multi,
      Xu = d$test_x,
      Yu = Yu_wrong_names,
      k = 30, b = 10, target_size = 350,
      fit_method = fit_pls(ncomp = 5),
      verbose = FALSE
    ),
    "Column names of Yr and Yu must match"
  )
})

# =============================================================================
# fit_method validation tests
# =============================================================================

test_that("gesearch rejects non-fit_pls fit_method", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_gesearch_data()
  
  expect_error(
    gesearch(
      Xr = d$train_x,
      Yr = d$train_y,
      Xu = d$test_x,
      k = 30, b = 10, target_size = 350,
      fit_method = fit_gpr(),
      verbose = FALSE
    ),
    "Only fit_pls\\(\\) is supported"
  )
})

test_that("gesearch warns when tune is used with reconstruction", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_gesearch_data()
  
  expect_warning(
    gesearch(
      Xr = d$train_x,
      Yr = d$train_y,
      Xu = d$test_x,
      k = 30, b = 20, retain = 0.90,
      target_size = 350,
      fit_method = fit_pls(ncomp = 10),
      optimization = "reconstruction",
      control = gesearch_control(tune = TRUE),
      verbose = FALSE,
      seed = 42
    ),
    "PLS components are not tuned when optimization includes 'reconstruction'"
  )
})

# =============================================================================
# Local CV validation tests
# =============================================================================

test_that("gesearch errors when local CV leaves fewer than 3 observations", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_gesearch_data()
  
  expect_error(
    gesearch(
      Xr = d$train_x,
      Yr = d$train_y,
      Xu = d$test_x,
      Yu = d$test_y,
      k = 10,  # small k
      b = 10,
      target_size = 350,
      fit_method = fit_pls(ncomp = 5),
      optimization = "response",
      control = gesearch_control(tune = TRUE, p = 0.95),  # leaves ~5% for hold-out
      verbose = FALSE
    ),
    "Local cross-validation requires at least 3 observations"
  )
})

test_that("gesearch errors when ncomp exceeds available samples", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_gesearch_data()
  
  expect_error(
    gesearch(
      Xr = d$train_x,
      Yr = d$train_y,
      Xu = d$test_x,
      k = 10,  # small k
      b = 10,
      target_size = 350,
      fit_method = fit_pls(ncomp = 15),  # more components than k
      control = gesearch_control(tune = FALSE),
      verbose = FALSE
    ),
    "More PLS components than available observations"
  )
})

# =============================================================================
# Range optimization validation tests
# =============================================================================

test_that("gesearch requires Yu_lims for range optimization", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_gesearch_data()
  
  expect_error(
    gesearch(
      Xr = d$train_x,
      Yr = d$train_y,
      Xu = d$test_x,
      k = 30, b = 10, target_size = 350,
      fit_method = fit_pls(ncomp = 5),
      optimization = "range",
      verbose = FALSE
    ),
    "'Yu_lims' is required for 'range' optimization"
  )
})

test_that("gesearch validates Yu_lims length", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_gesearch_data()
  
  expect_error(
    gesearch(
      Xr = d$train_x,
      Yr = d$train_y,
      Xu = d$test_x,
      k = 30, b = 10, target_size = 350,
      fit_method = fit_pls(ncomp = 5),
      optimization = "range",
      Yu_lims = c(0, 5, 10),
      verbose = FALSE
    ),
    "'Yu_lims' must be a numeric vector of length 2"
  )
})

test_that("gesearch with range optimization works", {
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
    optimization = c("reconstruction", "range"),
    Yu_lims = range(d$train_y),
    verbose = FALSE,
    seed = 42
  )
  
  expect_s3_class(gs, "gesearch")
})

# =============================================================================
# gesearch.formula tests
# =============================================================================

test_that("gesearch.formula validates formula class", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_gesearch_data()
  
  train_df <- data.frame(y = d$train_y, spc = I(d$train_x))
  test_df <- data.frame(y = d$test_y, spc = I(d$test_x))
  
  expect_error(
    gesearch(
      formula = "y ~ spc",  # character, not formula
      train = train_df,
      test = test_df,
      k = 30, b = 10, target_size = 350,
      fit_method = fit_pls(ncomp = 5)
    ),
    "'formula' argument must be a formula object"
  )
})

test_that("gesearch.formula requires fit_method", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_gesearch_data()
  
  train_df <- data.frame(y = d$train_y, spc = I(d$train_x))
  test_df <- data.frame(y = d$test_y, spc = I(d$test_x))
  
  expect_error(
    gesearch(
      formula = y ~ spc,
      train = train_df,
      test = test_df,
      k = 30, b = 10, target_size = 350
    ),
    "'fit_method' is missing"
  )
})

test_that("gesearch.formula warns when response missing in test", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_gesearch_data()
  
  train_df <- data.frame(y = d$train_y, spc = I(d$train_x))
  test_df <- data.frame(spc = I(d$test_x))  # No y column
  
  expect_warning(
    gs <- resemble::gesearch(
      formula = y ~ spc,
      train = train_df,
      test = test_df,
      k = 30, b = 20, retain = 0.90,
      target_size = 350,
      fit_method = fit_pls(ncomp = 10),
      optimization = "reconstruction",
      verbose = FALSE,
      seed = 42
    ),
    "y not found in test; assigned NA."
  )
  
  expect_s3_class(gs, "gesearch")
})

test_that("gesearch.formula errors when response optimization without response in test", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_gesearch_data()
  
  train_df <- data.frame(y = d$train_y, spc = I(d$train_x))
  test_df <- data.frame(spc = I(d$test_x))  # No y column
  
  expect_error(
    gesearch(
      formula = y ~ spc,
      train = train_df,
      test = test_df,
      k = 30, b = 10, target_size = 350,
      fit_method = fit_pls(ncomp = 5),
      optimization = "response"
    ),
    "'optimization = \"response\"' requires response values in test"
  )
})

test_that("gesearch.formula works with data.frame input", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_gesearch_data()
  
  train_df <- data.frame(y = d$train_y, spc = I(d$train_x))
  test_df <- data.frame(y = d$test_y, spc = I(d$test_x))
  
  gs <- gesearch(
    formula = y ~ spc,
    train = train_df,
    test = test_df,
    k = 30, b = 20, retain = 0.90,
    target_size = 350,
    fit_method = fit_pls(ncomp = 10),
    optimization = "reconstruction",
    verbose = FALSE,
    seed = 42
  )
  
  expect_s3_class(gs, "gesearch")
  expect_s3_class(gs, "gesearch.formula")
  expect_true("formula" %in% names(gs))
  expect_true("dataclasses" %in% names(gs))
})

# =============================================================================
# predict.gesearch formula model tests
# =============================================================================

test_that("predict.gesearch validates newdata for formula models", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_gesearch_data()
  
  train_df <- data.frame(y = d$train_y, spc = I(d$train_x))
  test_df <- data.frame(y = d$test_y, spc = I(d$test_x))
  
  gs <- gesearch(
    formula = y ~ spc,
    train = train_df,
    test = test_df,
    k = 30, b = 20, retain = 0.90,
    target_size = 350,
    fit_method = fit_pls(ncomp = 10),
    optimization = "reconstruction",
    verbose = FALSE,
    seed = 42
  )
  
  # Wrong type for newdata
  expect_error(
    predict(gs, newdata = as.list(d$test_x)),
    "'newdata' must be a data.frame or matrix"
  )
})

test_that("predict.gesearch validates missing variables for formula models", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_gesearch_data()
  
  train_df <- data.frame(y = d$train_y, spc = I(d$train_x))
  test_df <- data.frame(y = d$test_y, spc = I(d$test_x))
  
  gs <- gesearch(
    formula = y ~ spc,
    train = train_df,
    test = test_df,
    k = 30, b = 20, retain = 0.90,
    target_size = 350,
    fit_method = fit_pls(ncomp = 10),
    optimization = "reconstruction",
    verbose = FALSE,
    seed = 42
  )
  
  # Missing predictor variables
  wrong_df <- data.frame(wrong_name = I(d$test_x))
  
  expect_error(
    predict(gs, newdata = wrong_df),
    "Missing predictor variables"
  )
})

test_that("predict.gesearch works with formula model and data.frame", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_gesearch_data()
  
  train_df <- data.frame(y = d$train_y, spc = I(d$train_x))
  test_df <- data.frame(y = d$test_y, spc = I(d$test_x))
  
  gs <- gesearch(
    formula = y ~ spc,
    train = train_df,
    test = test_df,
    k = 30, b = 20, retain = 0.90,
    target_size = 350,
    fit_method = fit_pls(ncomp = 10),
    optimization = "reconstruction",
    verbose = FALSE,
    seed = 42
  )
  
  preds <- predict(gs, newdata = test_df)
  
  expect_type(preds, "list")
  expect_equal(nrow(preds[[1]]), nrow(d$test_x))
})

test_that("predict.gesearch works with formula model and matrix", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_gesearch_data()
  
  train_df <- data.frame(y = d$train_y, spc = I(d$train_x))
  test_df <- data.frame(y = d$test_y, spc = I(d$test_x))
  
  gs <- gesearch(
    formula = y ~ spc,
    train = train_df,
    test = test_df,
    k = 30, b = 20, retain = 0.90,
    target_size = 350,
    fit_method = fit_pls(ncomp = 10),
    optimization = "reconstruction",
    verbose = FALSE,
    seed = 42
  )
  
  # Predict with matrix
  preds <- predict(gs, newdata = d$test_x)
  
  expect_type(preds, "list")
  expect_equal(nrow(preds[[1]]), nrow(d$test_x))
})

test_that("predict.gesearch rejects non-matrix newdata for non-formula models", {
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
    predict(gs, newdata = as.data.frame(d$test_x)),
    "'newdata' must be a matrix"
  )
})

# =============================================================================
# Verbose output tests
# =============================================================================

test_that("gesearch verbose output works", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_gesearch_data()
  
  expect_output(
    gs <- gesearch(
      Xr = d$train_x,
      Yr = d$train_y,
      Xu = d$test_x,
      k = 30, b = 20, retain = 0.90,
      target_size = 350,
      fit_method = fit_pls(ncomp = 10),
      optimization = "reconstruction",
      verbose = TRUE,
      seed = 42
    ),
    "Generation"
  )
})

# =============================================================================
# Missing Xr/Yr validation tests
# =============================================================================

test_that("gesearch rejects missing values in Xr", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_gesearch_data()
  
  train_x_na <- d$train_x
  train_x_na[1, 1] <- NA
  
  expect_error(
    gesearch(
      Xr = train_x_na,
      Yr = d$train_y,
      Xu = d$test_x,
      k = 30, b = 10, target_size = 350,
      fit_method = fit_pls(ncomp = 5),
      verbose = FALSE
    ),
    "'Xr' contains missing values"
  )
})

test_that("gesearch rejects missing values in Yr", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_gesearch_data()
  
  train_y_na <- d$train_y
  train_y_na[1] <- NA
  
  expect_error(
    gesearch(
      Xr = d$train_x,
      Yr = train_y_na,
      Xu = d$test_x,
      k = 30, b = 10, target_size = 350,
      fit_method = fit_pls(ncomp = 5),
      verbose = FALSE
    ),
    "'Yr' contains missing values"
  )
})
