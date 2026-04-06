context("fit-and-k-constructors")

# =============================================================================
# neighbors_* constructor tests
# =============================================================================

test_that("neighbors_k constructor works", {
  ng <- neighbors_k(50)
  expect_s3_class(ng, c("neighbors_k", "neighbors"))
  expect_equal(ng$k, 50L)
  expect_equal(ng$method, "k")
  
  ng2 <- neighbors_k(c(100, 50, 80))
  expect_equal(ng2$k, c(50L, 80L, 100L))  # sorted
})


test_that("neighbors_k validates inputs", {
  expect_error(neighbors_k(), "required")
  expect_error(neighbors_k("a"), "numeric")
  expect_error(neighbors_k(3), "at least 4")
  expect_error(neighbors_k(c(10, 2)), "at least 4")
  expect_error(neighbors_k(10.5), "integer values")
})


test_that("neighbors_diss constructor works", {
  ng <- neighbors_diss(0.3)
  expect_s3_class(ng, c("neighbors_diss", "neighbors"))
  expect_equal(ng$threshold, 0.3)
  expect_equal(ng$k_min, 4L)
  expect_equal(ng$k_max, Inf)
  expect_equal(ng$method, "diss")
  
  ng2 <- neighbors_diss(c(0.3, 0.1, 0.2), k_min = 10, k_max = 150)
  expect_equal(ng2$threshold, c(0.1, 0.2, 0.3))  # sorted
  expect_equal(ng2$k_min, 10L)
  expect_equal(ng2$k_max, 150L)
})


test_that("neighbors_diss validates inputs", {
  expect_error(neighbors_diss(), "required")
  expect_error(neighbors_diss("a"), "numeric")
  expect_error(neighbors_diss(-0.1), "positive")
  expect_error(neighbors_diss(0.1, k_min = 2), "at least 4")
  expect_error(neighbors_diss(0.1, k_min = 10, k_max = 5), "greater than")
})


# =============================================================================
# fit_* constructor tests
# =============================================================================

test_that("fit_pls constructor works", {
  m <- fit_pls(ncomp = 10)
  expect_s3_class(m, c("fit_pls", "fit_method"))
  expect_equal(m$ncomp, 10L)
  expect_equal(m$method, "pls")
  expect_false(m$scale)
  expect_equal(m$max_iter, 100L)
  expect_equal(m$tol, 1e-6)
  
  m2 <- fit_pls(ncomp = 5, method = "mpls", scale = TRUE)
  expect_equal(m2$method, "mpls")
  expect_true(m2$scale)
  
  m3 <- fit_pls(ncomp = 5, method = "simpls")
  expect_equal(m3$method, "simpls")
})


test_that("fit_pls validates inputs", {
  expect_error(fit_pls(), "required")
  expect_error(fit_pls(0), "positive integer")
  expect_error(fit_pls("a"), "positive integer")
  expect_error(fit_pls(10, scale = "yes"), "TRUE or FALSE")
})


test_that("fit_wapls constructor works", {
  m <- fit_wapls(min_ncomp = 3, max_ncomp = 15)
  expect_s3_class(m, c("fit_wapls", "fit_method"))
  expect_equal(m$min_ncomp, 3L)
  expect_equal(m$max_ncomp, 15L)
  expect_equal(m$method, "mpls")  # default for wapls
  expect_false(m$scale)
  
  m2 <- fit_wapls(min_ncomp = 5, max_ncomp = 20, method = "simpls", scale = TRUE)
  expect_equal(m2$method, "simpls")
  expect_true(m2$scale)
})


test_that("fit_wapls validates inputs", {
  expect_error(fit_wapls(), "required")
  expect_error(fit_wapls(min_ncomp = 5), "required")
  expect_error(fit_wapls(min_ncomp = 0, max_ncomp = 10), "positive integer")
  expect_error(fit_wapls(min_ncomp = 10, max_ncomp = 5), "less than")
  expect_error(fit_wapls(min_ncomp = 5, max_ncomp = 5), "less than")
})


test_that("fit_gpr constructor works", {
  m <- fit_gpr()
  expect_s3_class(m, c("fit_gpr", "fit_method"))
  expect_equal(m$noise_variance, 0.001)
  expect_true(m$center)
  expect_true(m$scale)
  
  m2 <- fit_gpr(noise_variance = 0.01, center = FALSE, scale = FALSE)
  expect_equal(m2$noise_variance, 0.01)
  expect_false(m2$center)
  expect_false(m2$scale)
})


test_that("fit_gpr validates inputs", {
  expect_error(fit_gpr(noise_variance = 0), "positive")
  expect_error(fit_gpr(noise_variance = -1), "positive")
  expect_error(fit_gpr(noise_variance = "a"), "single numeric")
  expect_error(fit_gpr(center = "yes"), "TRUE or FALSE")
  expect_error(fit_gpr(scale = "yes"), "TRUE or FALSE")
})

