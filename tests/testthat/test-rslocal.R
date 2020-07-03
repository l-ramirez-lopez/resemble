context("test-rslocal")

test_that("rslocal works", {
  nirdata <- data("NIRsoil", package = "prospectr")
  
  Xu <- NIRsoil$spc[!as.logical(NIRsoil$train), ]
  Yu <- NIRsoil$CEC[!as.logical(NIRsoil$train)]
  
  Yr <- NIRsoil$CEC[as.logical(NIRsoil$train)]
  Xr <- NIRsoil$spc[as.logical(NIRsoil$train), ]
  
  Xu <- Xu[!is.na(Yu), ][1:20, ]
  Xr <- Xr[!is.na(Yr), ][1:40, ]
  
  Yu <- Yu[!is.na(Yu)][1:20]
  Yr <- Yr[!is.na(Yr)][1:40]
  
  ctrl_1 <- rs_control(retain_by = "probability",
                       percentile_type = 7,
                       tune = F,
                       number = 10,
                       p = 0.75,
                       verbose = FALSE,
                       allow_parallel = FALSE)
  
  rslocal_1 <- rslocal(Xr, Yr, Xu, Yu = Yu,
                       k = 14, b = 10, retain = 0.95,
                       method = local_fit_pls(pls_c = 5),
                       optimization = "response",
                       control = ctrl_1,
                       scale = FALSE)
  
  ctrl_2 <- rs_control(retain_by = "proportion",
                       tune = F,
                       number = 10,
                       p = 0.75,
                       verbose = FALSE,
                       allow_parallel = FALSE)
  
  rslocal_2 <- rslocal(Xr, Yr, Xu, Yu = Yu,
                       k = 14, b = 10, retain = 0.95,
                       method = local_fit_pls(pls_c = 5),
                       optimization = "response",
                       control = ctrl_2,
                       scale = FALSE)
  
  names(rslocal_2)
  
  output_names <- names(rslocal_1)
  expected_names <- c("y_local", "x_local", "indices", "iter_rmse", 
                      "n_removed", "validation_results", "control", 
                      "final_model", "documentation")
  
  expect_is(rslocal_1, "list")
  expect_is(rslocal_2, "list")
  expect_true(all(expected_names %in% output_names))
})
