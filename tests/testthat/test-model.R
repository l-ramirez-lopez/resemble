# test-model.R
# Unit tests for model(), model_control(), and predict.resemble_model()

context("test-model")

# =============================================================================
# Setup helper
# =============================================================================

.setup_model_data <- function(n_xr = 80, n_xu = 20, preprocess = TRUE) {
  data("NIRsoil", package = "prospectr")
  
  if (preprocess) {
    NIRsoil$spc <- prospectr::savitzkyGolay(NIRsoil$spc, p = 3, w = 11, m = 0)
  }
  
  Xu <- NIRsoil$spc[!as.logical(NIRsoil$train), ]
  Yu <- NIRsoil$CEC[!as.logical(NIRsoil$train)]
  Yr <- NIRsoil$CEC[as.logical(NIRsoil$train)]
  Xr <- NIRsoil$spc[as.logical(NIRsoil$train), ]
  
  Xu <- Xu[!is.na(Yu), ][seq_len(n_xu), ]
  Xr <- Xr[!is.na(Yr), ][seq_len(n_xr), ]
  Yu <- Yu[!is.na(Yu)][seq_len(n_xu)]
  Yr <- Yr[!is.na(Yr)][seq_len(n_xr)]
  
  list(Xr = Xr, Xu = Xu, Yr = Yr, Yu = Yu)
}


# =============================================================================
# model_control tests
# =============================================================================

test_that("model_control constructor works", {
  ctrl <- model_control()
  
  expect_type(ctrl, "list")
  expect_s3_class(ctrl, "model_control")
  expect_equal(ctrl$validation_type, "lgo")
  expect_equal(ctrl$number, 10L)
  expect_equal(ctrl$p, 0.75)
})


test_that("model_control accepts valid arguments", {
  ctrl_none <- model_control(validation_type = "none")
  expect_equal(ctrl_none$validation_type, "none")
  
  ctrl_custom <- model_control(validation_type = "lgo", number = 5, p = 0.8)
  expect_equal(ctrl_custom$number, 5L)
  expect_equal(ctrl_custom$p, 0.8)
})


test_that("model_control validates inputs", {
  expect_error(model_control(validation_type = "invalid"), "arg")
  expect_error(model_control(number = 0), "positive integer")
  expect_error(model_control(number = -5), "positive integer")
  expect_error(model_control(number = "ten"), "positive integer")
  expect_error(model_control(p = 0), "between 0 and 1")
  expect_error(model_control(p = 1), "between 0 and 1")
  expect_error(model_control(p = 1.5), "between 0 and 1")
  expect_error(model_control(p = -0.5), "between 0 and 1")
})


test_that("print.model_control works", {
  ctrl <- model_control()
  expect_output(print(ctrl), "validation_type")
  expect_output(print(ctrl), "number")
  expect_output(print(ctrl), "p")
  
  ctrl_none <- model_control(validation_type = "none")
  output <- capture.output(print(ctrl_none))
  expect_true(any(grepl("none", output)))
})


# =============================================================================
# model() input validation tests
# =============================================================================

test_that("model validates Xr and Yr dimensions", {
  skip_on_cran()
  Xr <- matrix(rnorm(100), nrow = 10)
  Yr <- matrix(rnorm(5), nrow = 5)
  
  expect_error(
    model(Xr, Yr, fit_method = fit_pls(ncomp = 3), verbose = FALSE),
    "same number of rows"
  )
})


test_that("model rejects Yr with multiple columns", {
  Xr <- matrix(rnorm(100), nrow = 10)
  Yr <- matrix(rnorm(20), nrow = 10, ncol = 2)
  
  expect_error(
    model(Xr, Yr, fit_method = fit_pls(ncomp = 3), verbose = FALSE),
    "exactly one column"
  )
})


test_that("model rejects NA values", {
  Xr <- matrix(rnorm(100), nrow = 10)
  Yr <- matrix(rnorm(10), nrow = 10)
  
  Xr_na <- Xr
  Xr_na[1, 1] <- NA
  expect_error(
    model(Xr_na, Yr, fit_method = fit_pls(ncomp = 3), verbose = FALSE),
    "missing values"
  )
  
  Yr_na <- Yr
  Yr_na[1] <- NA
  expect_error(
    model(Xr, Yr_na, fit_method = fit_pls(ncomp = 3), verbose = FALSE),
    "missing values"
  )
})


test_that("model rejects invalid fit_method", {
  Xr <- matrix(rnorm(100), nrow = 10)
  Yr <- matrix(rnorm(10), nrow = 10)
  
  expect_error(
    model(Xr, Yr, fit_method = list(method = "pls"), verbose = FALSE),
    "fit_method.*fit_\\*\\(\\)"
  )
})


test_that("model rejects fit_wapls", {
  Xr <- matrix(rnorm(100), nrow = 10)
  Yr <- matrix(rnorm(10), nrow = 10)
  
  expect_error(
    model(Xr, Yr, fit_method = fit_wapls(3, 10), verbose = FALSE),
    "fit_wapls.*not supported"
  )
})


test_that("model rejects invalid control", {
  Xr <- matrix(rnorm(100), nrow = 10)
  Yr <- matrix(rnorm(10), nrow = 10)
  
  expect_error(
    model(Xr, Yr, fit_method = fit_pls(ncomp = 3), 
          control = list(validation_type = "lgo"), verbose = FALSE),
    "model_control"
  )
})


# =============================================================================
# Basic model functionality tests - PLS
# =============================================================================

test_that("model with fit_pls without CV works", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_model_data()
  
  mod <- model(
    Xr = d$Xr,
    Yr = d$Yr,
    fit_method = fit_pls(ncomp = 6, scale = FALSE),
    control = model_control(validation_type = "none"),
    verbose = FALSE
  )
  
  expect_s3_class(mod, "resemble_model")
  expect_equal(mod$n_obs, nrow(d$Xr))
  expect_equal(mod$n_vars, ncol(d$Xr))
  expect_null(mod$cv_results)
  expect_true(!is.null(mod$model))
})


test_that("model with fit_pls with CV works", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_model_data()
  
  ctrl <- model_control(validation_type = "lgo", number = 5, p = 0.75)
  
  mod <- model(
    Xr = d$Xr,
    Yr = d$Yr,
    fit_method = fit_pls(ncomp = 8, scale = TRUE),
    control = ctrl,
    verbose = FALSE
  )
  
  expect_s3_class(mod, "resemble_model")
  expect_true(!is.null(mod$cv_results))
  expect_equal(nrow(mod$cv_results), 8)
  expect_true("rmse" %in% names(mod$cv_results))
  expect_true("r2" %in% names(mod$cv_results))
  expect_true("optimal" %in% names(mod$cv_results))
  expect_equal(sum(mod$cv_results$optimal), 1)
  expect_equal(which(mod$cv_results$optimal), which.min(mod$cv_results$rmse))
})


test_that("model with fit_pls (mpls method) works", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_model_data()
  
  mod <- model(
    Xr = d$Xr,
    Yr = d$Yr,
    fit_method = fit_pls(ncomp = 6, method = "mpls"),
    control = model_control(validation_type = "none"),
    verbose = FALSE
  )
  
  expect_s3_class(mod, "resemble_model")
})


test_that("model with fit_pls (simpls method) works", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_model_data()
  
  mod <- model(
    Xr = d$Xr,
    Yr = d$Yr,
    fit_method = fit_pls(ncomp = 6, method = "simpls"),
    control = model_control(validation_type = "none"),
    verbose = FALSE
  )
  
  expect_s3_class(mod, "resemble_model")
})


# =============================================================================
# Basic model functionality tests - GPR
# =============================================================================

test_that("model with fit_gpr without CV works", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_model_data(n_xr = 60)
  
  mod <- model(
    Xr = d$Xr,
    Yr = d$Yr,
    fit_method = fit_gpr(noise_variance = 0.01, scale = TRUE),
    control = model_control(validation_type = "none"),
    verbose = FALSE
  )
  
  expect_s3_class(mod, "resemble_model")
  expect_equal(mod$n_obs, nrow(d$Xr))
  expect_null(mod$cv_results)
  expect_true(!is.null(mod$model))
})


test_that("model with fit_gpr with CV works", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_model_data(n_xr = 60)
  
  ctrl <- model_control(validation_type = "lgo", number = 5, p = 0.7)
  
  mod <- model(
    Xr = d$Xr,
    Yr = d$Yr,
    fit_method = fit_gpr(noise_variance = 0.001, scale = TRUE),
    control = ctrl,
    verbose = FALSE
  )
  
  expect_s3_class(mod, "resemble_model")
  expect_true(!is.null(mod$cv_results))
  
  # GPR CV results: single row (no components)
  expect_equal(nrow(mod$cv_results), 1)
  expect_true("rmse" %in% names(mod$cv_results))
  expect_true("r2" %in% names(mod$cv_results))
  expect_true(mod$cv_results$rmse > 0)
})


# =============================================================================
# predict.resemble_model tests
# =============================================================================

test_that("predict works for PLS model", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_model_data()
  
  mod <- model(
    Xr = d$Xr,
    Yr = d$Yr,
    fit_method = fit_pls(ncomp = 6),
    control = model_control(validation_type = "none"),
    verbose = FALSE
  )
  
  
  preds <- predict(mod, d$Xu)
  
  expect_true(is.matrix(preds))
  expect_equal(nrow(preds), nrow(d$Xu))
  expect_equal(ncol(preds), 6)
  expect_equal(colnames(preds), paste0("ncomp", 1:6))
})


test_that("predict respects ncomp argument for PLS", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_model_data()
  
  mod <- model(
    Xr = d$Xr,
    Yr = d$Yr,
    fit_method = fit_pls(ncomp = 10),
    control = model_control(validation_type = "none"),
    verbose = FALSE
  )
  
  # Request specific number of components
  preds_5 <- predict(mod, d$Xu, ncomp = 5)
  expect_equal(ncol(preds_5), 1)
  
  # Vector of components
  preds_vec <- predict(mod, d$Xu, ncomp = c(2, 5, 8))
  expect_equal(ncol(preds_vec), 3)
  expect_equal(colnames(preds_vec), c("ncomp2", "ncomp5", "ncomp8"))
})


test_that("predict rejects ncomp exceeding model", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_model_data()
  
  mod <- model(
    Xr = d$Xr,
    Yr = d$Yr,
    fit_method = fit_pls(ncomp = 5),
    control = model_control(validation_type = "none"),
    verbose = FALSE
  )
  
  expect_error(predict(mod, d$Xu, ncomp = 10), "exceeds")
})


test_that("predict works for GPR model", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_model_data(n_xr = 60)
  
  mod <- model(
    Xr = d$Xr,
    Yr = d$Yr,
    fit_method = fit_gpr(noise_variance = 0.01, scale = TRUE),
    control = model_control(validation_type = "none"),
    verbose = FALSE
  )
  
  preds <- predict(mod, d$Xu)
  
  expect_true(is.matrix(preds))
  expect_equal(nrow(preds), nrow(d$Xu))
  expect_equal(ncol(preds), 1)
  expect_equal(colnames(preds), "predicted")
})


test_that("predict rejects wrong number of columns", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_model_data()
  
  mod <- model(
    Xr = d$Xr,
    Yr = d$Yr,
    fit_method = fit_pls(ncomp = 5),
    control = model_control(validation_type = "none"),
    verbose = FALSE
  )
  
  Xu_wrong <- matrix(rnorm(100), nrow = 10, ncol = 10)
  expect_error(predict(mod, Xu_wrong), "columns")
})


test_that("predict preserves rownames from newdata", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_model_data(n_xu = 5)
  rownames(d$Xu) <- paste0("sample_", 1:5)
  
  mod <- model(
    Xr = d$Xr,
    Yr = d$Yr,
    fit_method = fit_pls(ncomp = 5),
    control = model_control(validation_type = "none"),
    verbose = FALSE
  )
  
  preds <- predict(mod, d$Xu)
  expect_equal(rownames(preds), paste0("sample_", 1:5))
})


# =============================================================================
# Expected results tests (skipped on CRAN)
# =============================================================================
## Sanity checks ensuring results stay within plausible bounds rather 
## than testing for exact values. This catches regressions where something 
## breaks catastrophically.

test_that("model delivers expected results", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  data("NIRsoil", package = "prospectr")
  NIRsoil$spc <- prospectr::savitzkyGolay(NIRsoil$spc, p = 3, w = 11, m = 0)
  
  Xu <- NIRsoil$spc[!as.logical(NIRsoil$train), ]
  Yu <- NIRsoil$CEC[!as.logical(NIRsoil$train)]
  Yr <- NIRsoil$CEC[as.logical(NIRsoil$train)]
  Xr <- NIRsoil$spc[as.logical(NIRsoil$train), ]
  
  Xu <- Xu[!is.na(Yu), ]
  Xr <- Xr[!is.na(Yr), ]
  Yu <- Yu[!is.na(Yu)]
  Yr <- Yr[!is.na(Yr)]
  
  ctrl <- model_control(validation_type = "lgo", number = 10, p = 0.75)
  
  # PLS
  pls_mod <- model(
    Xr = Xr,
    Yr = Yr,
    fit_method = fit_pls(ncomp = 15, scale = TRUE),
    control = ctrl,
    verbose = FALSE
  )
  
  # GPR
  gpr_mod <- model(
    Xr = Xr,
    Yr = Yr,
    fit_method = fit_gpr(noise_variance = 0.001, scale = TRUE),
    control = ctrl,
    verbose = FALSE
  )
  
  # Check PLS CV RMSE bounds
  expect_true(all(pls_mod$cv_results$rmse < 6))
  expect_true(min(pls_mod$cv_results$rmse) < 4.5)
  
  # Check PLS CV R2 bounds
  expect_true(max(pls_mod$cv_results$r2) > 0.5)
  
  # Check GPR CV RMSE bounds
  expect_true(gpr_mod$cv_results$rmse < 5)
  expect_true(gpr_mod$cv_results$rmse > 1)
  
  # Check GPR CV R2 bounds
  expect_true(gpr_mod$cv_results$r2 > 0.4)
  
  # Predictions should be in reasonable range
  pls_preds <- predict(pls_mod, Xu)
  gpr_preds <- predict(gpr_mod, Xu)
  
  # Check prediction correlation with Yu
  pls_cor <- cor(pls_preds[, which.min(pls_mod$cv_results$rmse)], Yu)
  gpr_cor <- cor(gpr_preds[, 1], Yu)
  
  
  expect_true(pls_cor > 0.6)
  expect_true(gpr_cor > 0.5)
})


test_that("model with different PLS methods delivers consistent results", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_model_data()
  
  ctrl <- model_control(validation_type = "none")
  
  # NIPALS (default)
  mod_nipals <- model(
    Xr = d$Xr,
    Yr = d$Yr,
    fit_method = fit_pls(ncomp = 6, method = "pls"),
    control = ctrl,
    verbose = FALSE
  )
  
  # SIMPLS
  mod_simpls <- model(
    Xr = d$Xr,
    Yr = d$Yr,
    fit_method = fit_pls(ncomp = 6, method = "simpls"),
    control = ctrl,
    verbose = FALSE
  )
  
  # MPLS
  mod_mpls <- model(
    Xr = d$Xr,
    Yr = d$Yr,
    fit_method = fit_pls(ncomp = 6, method = "mpls"),
    control = ctrl,
    verbose = FALSE
  )
  
  # Predictions should be similar (not identical due to different algorithms)
  preds_nipals <- predict(mod_nipals, d$Xu)
  preds_simpls <- predict(mod_simpls, d$Xu)
  preds_mpls <- predict(mod_mpls, d$Xu)
  
  # Correlation between methods should be very high
  expect_true(cor(preds_nipals[, 6], preds_simpls[, 6]) > 0.99)
  expect_true(cor(preds_nipals[, 6], preds_mpls[, 6]) > 0.99)
})


test_that("model predictions are numerically consistent across calls", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_model_data()
  
  mod <- model(
    Xr = d$Xr,
    Yr = d$Yr,
    fit_method = fit_pls(ncomp = 5),
    control = model_control(validation_type = "none"),
    verbose = FALSE
  )
  
  preds1 <- predict(mod, d$Xu)
  preds2 <- predict(mod, d$Xu)
  
  expect_equal(preds1, preds2)
})


test_that("model stores fit_method and control in output", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_model_data()
  
  fm <- fit_pls(ncomp = 6, scale = TRUE, method = "simpls")
  ctrl <- model_control(validation_type = "none")
  
  mod <- model(d$Xr, d$Yr, fit_method = fm, control = ctrl, verbose = FALSE)
  
  expect_equal(mod$fit_method$ncomp, 6)
  expect_equal(mod$fit_method$scale, TRUE)
  expect_equal(mod$fit_method$method, "simpls")
  expect_equal(mod$control$validation_type, "none")
})


test_that("model stores call attribute", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_model_data()
  
  mod <- model(
    Xr = d$Xr,
    Yr = d$Yr,
    fit_method = fit_pls(ncomp = 5),
    control = model_control(validation_type = "none"),
    verbose = FALSE
  )
  
  expect_true(!is.null(attr(mod, "call")))
  expect_true(inherits(attr(mod, "call"), "call"))
})
