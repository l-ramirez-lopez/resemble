context("test-ncomp_by-constructors")

# =============================================================================
# ncomp_* constructor tests
# =============================================================================

# test-constructors.R  
test_that("ncomp_by_var returns correct class", {
  x <- ncomp_by_var(0.01)
  expect_s3_class(x, "ncomp_by_var")
})

test_that("ncomp_by_var constructor works", {
  nc <- ncomp_by_var(0.01)
  expect_s3_class(nc, c("ncomp_by_var", "ncomp_selection"))
  expect_equal(nc$min_var, 0.01)
  expect_equal(nc$max_ncomp, 40L)
  
  nc2 <- ncomp_by_var(0.05, max_ncomp = 20)
  expect_equal(nc2$min_var, 0.05)
  expect_equal(nc2$max_ncomp, 20L)
})


test_that("ncomp_by_var validates inputs", {
  expect_error(ncomp_by_var(0), "\\(0, 1\\]")
  expect_error(ncomp_by_var(1.5), "\\(0, 1\\]")
  expect_error(ncomp_by_var("a"), "\\(0, 1\\]")
  expect_error(ncomp_by_var(0.01, max_ncomp = 0), "positive integer")
})


test_that("ncomp_by_cumvar constructor works", {
  nc <- ncomp_by_cumvar(0.99)
  expect_s3_class(nc, c("ncomp_by_cumvar", "ncomp_selection"))
  expect_equal(nc$min_cumvar, 0.99)
  expect_equal(nc$max_ncomp, 40L)
  
  nc2 <- ncomp_by_cumvar(0.95, max_ncomp = 50)
  expect_equal(nc2$min_cumvar, 0.95)
  expect_equal(nc2$max_ncomp, 50L)
})


test_that("ncomp_by_cumvar validates inputs", {
  expect_error(ncomp_by_cumvar(0), "\\(0, 1\\]")
  expect_error(ncomp_by_cumvar(1.5), "\\(0, 1\\]")
  expect_error(ncomp_by_cumvar(0.99, max_ncomp = -1), "positive integer")
})


test_that("ncomp_by_opc constructor works", {
  nc <- ncomp_by_opc()
  expect_s3_class(nc, c("ncomp_by_opc", "ncomp_selection"))
  expect_equal(nc$max_ncomp, 40L)
  
  nc2 <- ncomp_by_opc(max_ncomp = 30)
  expect_equal(nc2$max_ncomp, 30L)
})


test_that("ncomp_by_opc validates inputs", {
  expect_error(ncomp_by_opc(max_ncomp = 0), "positive integer")
  expect_error(ncomp_by_opc(max_ncomp = "a"), "positive integer")
})


test_that("ncomp_fixed constructor works", {
  nc <- ncomp_fixed(10)
  expect_s3_class(nc, c("ncomp_fixed", "ncomp_selection"))
  expect_equal(nc$ncomp, 10L)
})


test_that("ncomp_fixed validates inputs", {
  expect_error(ncomp_fixed(), "required")
  expect_error(ncomp_fixed(0), "positive integer")
  expect_error(ncomp_fixed(-5), "positive integer")
  expect_error(ncomp_fixed("a"), "positive integer")
})
