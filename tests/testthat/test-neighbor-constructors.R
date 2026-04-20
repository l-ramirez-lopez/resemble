context("test-neighbor-constructors")


test_that("neighbors_k validates input", {
  expect_error(neighbors_k(-1))
  expect_s3_class(neighbors_k(10), "neighbors_k")
})


test_that("neighbors_diss validates input", {
  expect_error(neighbors_diss(-1, 50, 100))
  expect_s3_class(neighbors_diss(10), "neighbors_diss")
})


# tests/testthat/test-neighbors-print.R
test_that("print.neighbors_k prints compact output for short k", {
  x <- neighbors_k(c(4, 6, 8))
  
  out <- capture.output(ret <- print(x))
  
  expect_identical(
    out,
    c(
      "Neighbor selection: fixed k", "  k : 4, 6, 8 " 
    )
  )
  expect_identical(ret, x)
})

test_that("print.neighbors_k truncates long k vectors", {
  x <- neighbors_k(c(4, 5, 6, 7, 8, 9, 10))
  
  out <- capture.output(ret <- print(x))
  
  expect_identical(
    out,
    c(
      "Neighbor selection: fixed k", 
      "  k : 4, 5, 6, ..., 9, 10 " 
    )
  )
  expect_identical(ret, x)
})

test_that("print.neighbors_diss prints finite k_max", {
  x <- neighbors_diss(
    threshold = c(0.1, 0.2, 0.3),
    k_min = 4L,
    k_max = 12L
  )
  
  out <- capture.output(ret <- print(x))
  
  expect_identical(
    out,
    c(
      "Neighbor selection: dissimilarity threshold", 
      "   threshold : 0.1, 0.2, 0.3 ",               
      "   k_min     : 4 ", 
      "   k_max     : 12 "        
    )
  )
  expect_identical(ret, x)
})

test_that("print.neighbors_diss truncates long threshold vectors", {
  x <- neighbors_diss(
    threshold = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7),
    k_min = 4L,
    k_max = 10L
  )
  
  out <- capture.output(ret <- print(x))
  
  expect_identical(
    out,
    c(
      "Neighbor selection: dissimilarity threshold", 
      "   threshold : 0.1, 0.2, 0.3, ..., 0.6, 0.7 ",
      "   k_min     : 4 ", 
      "   k_max     : 10 "  
    )
  )
  expect_identical(ret, x)
})

test_that("print.neighbors_diss prints Inf k_max", {
  x <- neighbors_diss(
    threshold = c(0.1, 0.2),
    k_min = 4L,
    k_max = Inf
  )
  
  out <- capture.output(ret <- print(x))
  
  expect_identical(
    out,
    c(
      "Neighbor selection: dissimilarity threshold", 
      "   threshold : 0.1, 0.2 ",
      "   k_min     : 4 ", 
      "   k_max     : Inf "
    )
  )
  expect_identical(ret, x)
})

test_that("print.neighbors prints all fields for generic neighbors object", {
  x <- structure(
    list(
      method = "custom",
      alpha = 0.25,
      beta = 7L
    ),
    class = "neighbors"
  )
  
  out <- capture.output(ret <- print(x))
  
  expect_identical(
    out,
    c(
      "Neighbor selection method: custom ",
      "  alpha : 0.25 ",
      "  beta : 7 "
    )
  )
  expect_identical(ret, x)
})
test_that("print.neighbors handles vector fields in generic neighbors objects", {
  x <- structure(
    list(
      method = "custom",
      ids = c(4L, 5L, 6L)
    ),
    class = "neighbors"
  )
  
  out <- capture.output(ret <- print(x))
  
  expect_identical(
    out,
    c(
      "Neighbor selection method: custom ",
      "  ids : 4 5 6 "
    )
  )
  expect_identical(ret, x)
})
