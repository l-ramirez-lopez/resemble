context("test-fit-constructors")

# =============================================================================
# fit_pls constructor tests
# =============================================================================

test_that("fit_pls constructor works", {
  f1 <- fit_pls(ncomp = 10)
  expect_s3_class(f1, "fit_pls")
  expect_s3_class(f1, "fit_method")
  expect_equal(f1$ncomp, 10L)
  expect_equal(f1$method, "pls")
  expect_false(f1$scale)
  expect_equal(f1$max_iter, 100L)
  expect_equal(f1$tol, 1e-6)
  
  f2 <- fit_pls(ncomp = 15, method = "mpls", scale = TRUE)
  expect_equal(f2$ncomp, 15L)
  expect_equal(f2$method, "mpls")
  expect_true(f2$scale)
  
  f3 <- fit_pls(ncomp = 10, method = "simpls")
  expect_equal(f3$method, "simpls")
})


test_that("fit_pls requires ncomp", {
  expect_error(fit_pls(), "ncomp.*required")
})


test_that("fit_pls validates ncomp", {
  expect_error(fit_pls(ncomp = 0), "positive integer")
  expect_error(fit_pls(ncomp = -5), "positive integer")
  expect_error(fit_pls(ncomp = "ten"), "positive integer")
  expect_error(fit_pls(ncomp = c(5, 10)), "positive integer")
})


test_that("fit_pls validates scale", {
  expect_error(fit_pls(ncomp = 10, scale = "yes"), "TRUE or FALSE")
  expect_error(fit_pls(ncomp = 10, scale = c(TRUE, FALSE)), "TRUE or FALSE")
})


test_that("fit_pls validates method", {
  expect_error(fit_pls(ncomp = 10, method = "invalid"))
})


test_that("print.fit_pls works", {
  f <- fit_pls(ncomp = 10, method = "mpls")
  expect_output(print(f), "Fitting method: pls")
  expect_output(print(f), "ncomp.*: 10")
  expect_output(print(f), "method.*: mpls")
})


# =============================================================================
# fit_wapls constructor tests
# =============================================================================

test_that("fit_wapls constructor works", {
  f1 <- fit_wapls(min_ncomp = 3, max_ncomp = 15)
  expect_s3_class(f1, "fit_wapls")
  expect_s3_class(f1, "fit_method")
  expect_equal(f1$min_ncomp, 3L)
  expect_equal(f1$max_ncomp, 15L)
  expect_equal(f1$method, "mpls")  # default for wapls
  expect_false(f1$scale)
  
  f2 <- fit_wapls(min_ncomp = 5, max_ncomp = 20, method = "simpls", scale = TRUE)
  expect_equal(f2$method, "simpls")
  expect_true(f2$scale)
})


test_that("fit_wapls requires both min_ncomp and max_ncomp", {
  expect_error(fit_wapls(min_ncomp = 3), "required")
  expect_error(fit_wapls(max_ncomp = 15), "required")
})


test_that("fit_wapls validates min_ncomp < max_ncomp", {
  expect_error(fit_wapls(min_ncomp = 15, max_ncomp = 10), "less than")
  expect_error(fit_wapls(min_ncomp = 10, max_ncomp = 10), "less than")
})


test_that("fit_wapls validates ncomp values", {
  expect_error(fit_wapls(min_ncomp = 0, max_ncomp = 10), "positive integer")
  expect_error(fit_wapls(min_ncomp = 3, max_ncomp = -5), "positive integer")
})


test_that("print.fit_wapls works", {
  f <- fit_wapls(min_ncomp = 3, max_ncomp = 15)
  expect_output(print(f), "Fitting method: wapls")
  expect_output(print(f), "min_ncomp.*: 3")
  expect_output(print(f), "max_ncomp.*: 15")
})


# =============================================================================
# fit_gpr constructor tests
# =============================================================================

test_that("fit_gpr constructor works", {
  f1 <- fit_gpr()
  expect_s3_class(f1, "fit_gpr")
  expect_s3_class(f1, "fit_method")
  expect_equal(f1$noise_variance, 0.001)
  expect_true(f1$center)
  expect_true(f1$scale)
  
  f2 <- fit_gpr(noise_variance = 0.01, center = FALSE, scale = FALSE)
  expect_equal(f2$noise_variance, 0.01)
  expect_false(f2$center)
  expect_false(f2$scale)
})


test_that("fit_gpr validates noise_variance", {
  expect_error(fit_gpr(noise_variance = 0), "positive")
  expect_error(fit_gpr(noise_variance = -0.01), "positive")
  expect_error(fit_gpr(noise_variance = "small"), "numeric")
  expect_error(fit_gpr(noise_variance = c(0.001, 0.01)), "single numeric")
})


test_that("fit_gpr validates center and scale", {
  expect_error(fit_gpr(center = "yes"), "TRUE or FALSE")
  expect_error(fit_gpr(scale = "yes"), "TRUE or FALSE")
})


test_that("print.fit_gpr works", {
  f <- fit_gpr(noise_variance = 0.01)
  expect_output(print(f), "Fitting method: gpr")
  expect_output(print(f), "noise_variance.*: 0.01")
})
