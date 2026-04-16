context("test-neighbor-constructors")

test_that("neighbors_k validates input", {
  expect_error(neighbors_k(-1))
  expect_s3_class(neighbors_k(10), "neighbors_k")
})


test_that("neighbors_diss validates input", {
  expect_error(neighbors_diss(-1, 50, 100))
  expect_s3_class(neighbors_diss(10), "neighbors_diss")
})
