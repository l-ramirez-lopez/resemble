context("test-mbl")

test_that("mbl works", {
  nirdata <- data("NIRsoil", package = "prospectr")

  Xu <- NIRsoil$spc[!as.logical(NIRsoil$train), ]
  Yu <- NIRsoil$CEC[!as.logical(NIRsoil$train)]

  Yr <- NIRsoil$CEC[as.logical(NIRsoil$train)]
  Xr <- NIRsoil$spc[as.logical(NIRsoil$train), ]

  Xu <- Xu[!is.na(Yu), ][1:20, ]
  Xr <- Xr[!is.na(Yr), ][1:40, ]

  Yu <- Yu[!is.na(Yu)][1:20]
  Yr <- Yr[!is.na(Yr)][1:40]

  k_test <- seq(25, 35, by = 10)
  k_diss_test <- 0.1
  k_range_test <- c(15, 30)

  ctrl_1 <- mbl_control(
    validation_type = c("NNv", "local_cv"),
    number = 4, p = 0.5,
    allow_parallel = FALSE
  )

  gpr <- mbl(
    Xr = Xr, Yr = Yr, Xu = Xu, Yu = Yu,
    k = k_test,
    method = local_fit_gpr(),
    control = ctrl_1,
    verbose = FALSE
  )

  pls <- mbl(
    Xr = Xr, Yr = Yr, Xu = Xu, Yu = Yu,
    k = k_test,
    method = local_fit_pls(5),
    control = ctrl_1,
    verbose = FALSE
  )
  
  mpls <- mbl(
    Xr = Xr, Yr = Yr, Xu = Xu, Yu = Yu,
    k = k_test,
    method = local_fit_pls(5, modified = TRUE),
    control = ctrl_1,
    verbose = FALSE
  )
  
  wapls <- mbl(
    Xr = Xr, Yr = Yr, Xu = Xu, Yu = Yu,
    k = k_test,
    method = local_fit_wapls(3, 5),
    control = ctrl_1,
    verbose = FALSE
  )

  wampls <- mbl(
    Xr = Xr, Yr = Yr, Xu = Xu, Yu = Yu,
    k = k_test,
    method = local_fit_wapls(3, 5, modified = TRUE),
    control = ctrl_1,
    verbose = FALSE
  )
  
  gpr_k_diss <- mbl(
    Xr = Xr, Yr = Yr, Xu = Xu, Yu = Yu,
    k_diss = k_diss_test, k_range = k_range_test,
    method = local_fit_gpr(),
    control = ctrl_1,
    verbose = FALSE
  )

  pls_k_diss <- mbl(
    Xr = Xr, Yr = Yr, Xu = Xu, Yu = Yu,
    k_diss = k_diss_test, k_range = k_range_test,
    method = local_fit_pls(5),
    control = ctrl_1,
    verbose = FALSE
  )

  wapls_k_diss <- mbl(
    Xr = Xr, Yr = Yr, Xu = Xu, Yu = Yu,
    k_diss = k_diss_test, k_range = k_range_test,
    method = local_fit_wapls(3, 5),
    control = ctrl_1,
    verbose = FALSE
  )


  group_test <- mbl(
    Xr = Xr, Yr = Yr, Xu = Xu, Yu = Yu,
    k = k_test,
    method = local_fit_pls(5),
    control = ctrl_1,
    group = rep(c(1, 2), length.out = nrow(Xr)),
    verbose = FALSE
  )

  output_names <- names(gpr)
  expected_names <- c(
    "call", "cntrl_param", "dissimilarities", "Xu_neighbors",
    "n_predictions", "gh", "validation_results", "results",
    "documentation"
  )

  expect_is(gpr, "list")
  expect_is(pls, "list")
  expect_is(wapls, "list")
  expect_is(mpls, "list")
  expect_is(wampls, "list")
  expect_is(gpr_k_diss, "list")
  expect_is(pls_k_diss, "list")
  expect_is(wapls_k_diss, "list")
  expect_is(group_test, "list")
  expect_true(all(expected_names %in% output_names))
})

test_that("mbl delivers expeted results", {
  skip_on_cran()
  skip_on_travis()
  require(prospectr)
  nirdata <- data("NIRsoil", package = "prospectr")
  NIRsoil$spc <- prospectr::savitzkyGolay(NIRsoil$spc, p = 3, w = 11, m = 0)

  Xu <- NIRsoil$spc[!as.logical(NIRsoil$train), ]
  Yu <- NIRsoil$CEC[!as.logical(NIRsoil$train)]

  Yr <- NIRsoil$CEC[as.logical(NIRsoil$train)]
  Xr <- NIRsoil$spc[as.logical(NIRsoil$train), ]

  Xu <- Xu[!is.na(Yu), ]
  Xr <- Xr[!is.na(Yr), ]

  Yu <- Yu[!is.na(Yu)]
  Yr <- Yr[!is.na(Yr)]

  # tune locally
  ctrl_1 <- mbl_control(
    validation_type = c("NNv", "local_cv"),
    number = 4, p = 0.8,
    tune_locally = TRUE,
    allow_parallel = FALSE
  )

  k_test <- c(40, 150)
  k_diss_test <- 0.1
  k_range_test <- c(20, 100)
  pls_wapls <- c(3, 15)
  pls_pls <- c(10)
  grpnoisevar <- 0.0001

  tseed <- 141020

  set.seed(tseed)
  gpr <- mbl(
    Xr = Xr, Yr = Yr, Xu = Xu, Yu = Yu,
    k = k_test,
    method = local_fit_gpr(grpnoisevar),
    control = ctrl_1,
    verbose = FALSE
  )

  set.seed(tseed)
  pls <- mbl(
    Xr = Xr, Yr = Yr, Xu = Xu, Yu = Yu,
    k = k_test,
    method = local_fit_pls(pls_pls),
    control = ctrl_1,
    verbose = FALSE
  )

  set.seed(tseed)
  wapls <- mbl(
    Xr = Xr, Yr = Yr, Xu = Xu, Yu = Yu,
    k = k_test,
    method = local_fit_wapls(pls_wapls[1], pls_wapls[2]),
    control = ctrl_1,
    verbose = FALSE
  )
  
  set.seed(tseed)
  mpls <- mbl(
    Xr = Xr, Yr = Yr, Xu = Xu, Yu = Yu,
    k = k_test,
    method = local_fit_pls(pls_pls, modified = TRUE),
    control = ctrl_1,
    verbose = FALSE
  )
  
  set.seed(tseed)
  wampls <- mbl(
    Xr = Xr, Yr = Yr, Xu = Xu, Yu = Yu,
    k = k_test,
    method = local_fit_wapls(pls_wapls[1], pls_wapls[2], modified = TRUE),
    control = ctrl_1,
    verbose = FALSE
  )
  

  set.seed(tseed)
  gpr_k_diss <- mbl(
    Xr = Xr, Yr = Yr, Xu = Xu, Yu = Yu,
    k_diss = k_diss_test, k_range = k_range_test,
    method = local_fit_gpr(),
    control = ctrl_1,
    verbose = FALSE
  )

  set.seed(tseed)
  pls_k_diss <- mbl(
    Xr = Xr, Yr = Yr, Xu = Xu, Yu = Yu,
    k_diss = k_diss_test, k_range = k_range_test,
    method = local_fit_pls(pls_pls),
    control = ctrl_1,
    verbose = FALSE
  )

  set.seed(tseed)
  wapls_k_diss <- mbl(
    Xr = Xr, Yr = Yr, Xu = Xu, Yu = Yu,
    k_diss = k_diss_test, k_range = k_range_test,
    method = local_fit_wapls(pls_wapls[1], pls_wapls[2]),
    control = ctrl_1,
    verbose = FALSE
  )

  set.seed(tseed)
  xgroup <- rep((1:(floor(nrow(Xr) / 2))), each = 2)
  pls_group <- mbl(
    Xr = Xr, Yr = Yr, Xu = Xu, Yu = Yu,
    k = k_test,
    method = local_fit_pls(pls_pls),
    control = ctrl_1, group = xgroup,
    verbose = FALSE
  )

  set.seed(tseed)
  pls_group_local <- mbl(
    Xr = Xr, Yr = Yr, Xu = Xu, Yu = Yu,
    k = k_test,
    method = local_fit_pls(pls_pls),
    control = ctrl_1, group = xgroup,
    .local = TRUE, pre_k = 200,
    verbose = FALSE
  )


  cv_gpr <- c(
    gpr$validation_results$local_cross_validation$rmse < 1.8,
    gpr$validation_results$local_cross_validation$rmse > 1.5
  )

  cv_pls <- c(
    pls$validation_results$local_cross_validation$rmse < 1.8,
    pls$validation_results$local_cross_validation$rmse > 1.4
  )

  cv_wapls <- c(
    wapls$validation_results$local_cross_validation$rmse < 1.8,
    wapls$validation_results$local_cross_validation$rmse > 1.4
  )
  
  
  cv_mpls <- c(
    mpls$validation_results$local_cross_validation$rmse < 1.8,
    mpls$validation_results$local_cross_validation$rmse > 1.5
  )
  
  cv_wampls <- c(
    wampls$validation_results$local_cross_validation$rmse < 1.68,
    wampls$validation_results$local_cross_validation$rmse > 1.45
  )

  cv_gpr_k_diss <- c(
    gpr_k_diss$validation_results$local_cross_validation$rmse < 1.9,
    gpr_k_diss$validation_results$local_cross_validation$rmse > 1.6
  )

  cv_pls_k_diss <- c(
    pls_k_diss$validation_results$local_cross_validation$rmse < 1.7,
    pls_k_diss$validation_results$local_cross_validation$rmse > 1.4
  )

  cv_wapls_k_diss <- c(
    wapls_k_diss$validation_results$local_cross_validation$rmse < 2,
    wapls_k_diss$validation_results$local_cross_validation$rmse > 1.4
  )

  cv_group <- c(
    pls_group$validation_results$local_cross_validation$rmse < 2,
    pls_group$validation_results$local_cross_validation$rmse > 1.4
  )

  cv_group_local <- c(
    pls_group_local$validation_results$local_cross_validation$rmse < 2,
    pls_group_local$validation_results$local_cross_validation$rmse > 1
  )

  nnv_group_local <- pls_group_local$validation_results$nearest_neighbor_validation$r2 > 0.4

  nnv_group <- pls_group$validation_results$nearest_neighbor_validation$r2 > 0.7

  nnv_gpr <- gpr$validation_results$nearest_neighbor_validation$r2 > 0.81

  nnv_pls <- pls$validation_results$nearest_neighbor_validation$r2 > 0.74

  nnv_wapls <- wapls$validation_results$nearest_neighbor_validation$r2 > 0.80

  nnv_gpr_k_diss <- gpr_k_diss$validation_results$nearest_neighbor_validation$r2 > 0.81

  nnv_pls_k_diss <- pls_k_diss$validation_results$nearest_neighbor_validation$r2 > 0.81

  nnv_wapls_k_diss <- wapls_k_diss$validation_results$nearest_neighbor_validation$r2 > 0.81


  yuv_gpr <- gpr$validation_results$Yu_prediction_statistics$r2 > 0.72

  yuv_pls <- pls$validation_results$Yu_prediction_statistics$r2 > 0.67

  yuv_wapls <- wapls$validation_results$Yu_prediction_statistics$r2 > 0.69

  yuv_gpr_k_diss <- gpr_k_diss$validation_results$Yu_prediction_statistics$r2 > 0.72

  yuv_pls_k_diss <- pls_k_diss$validation_results$Yu_prediction_statistics$r2 > 0.60

  yuv_wapls_k_diss <- wapls_k_diss$validation_results$Yu_prediction_statistics$r2 > 0.65

  expect_true(all(cv_gpr))
  expect_true(all(cv_pls))
  expect_true(all(cv_wapls))
  expect_true(all(cv_mpls))
  expect_true(all(cv_wampls))
  expect_true(all(cv_gpr_k_diss))
  expect_true(all(cv_pls_k_diss))
  expect_true(all(cv_wapls_k_diss))
  expect_true(all(cv_group))
  expect_true(all(cv_group_local))

  expect_true(all(nnv_gpr))
  expect_true(all(nnv_pls))
  expect_true(all(nnv_wapls))
  expect_true(all(nnv_gpr_k_diss))
  expect_true(all(nnv_pls_k_diss))
  expect_true(all(nnv_wapls_k_diss))
  expect_true(all(nnv_group))
  expect_true(all(nnv_group_local))

  expect_true(all(yuv_gpr))
  expect_true(all(yuv_pls))
  expect_true(all(yuv_wapls))
  expect_true(all(yuv_gpr_k_diss))
  expect_true(all(yuv_pls_k_diss))
  expect_true(all(yuv_wapls_k_diss))
})


test_that("mbl with external disstances works", {
  tol <- 1e-10
  nirdata <- data("NIRsoil", package = "prospectr")
  
  # Proprocess the data using detrend plus first derivative with Savitzky and
  # Golay smoothing filter
  sg_det <- savitzkyGolay(
    detrend(NIRsoil$spc,
            wav = as.numeric(colnames(NIRsoil$spc))
    ),
    m = 1,
    p = 1,
    w = 7
  )
  
  NIRsoil$spc_pr <- sg_det
  
  # split into training and testing sets
  test_x <- NIRsoil$spc_pr[NIRsoil$train == 0 & !is.na(NIRsoil$CEC), ]
  test_y <- NIRsoil$CEC[NIRsoil$train == 0 & !is.na(NIRsoil$CEC)]
  
  train_y <- NIRsoil$CEC[NIRsoil$train == 1 & !is.na(NIRsoil$CEC)]
  train_x <- NIRsoil$spc_pr[NIRsoil$train == 1 & !is.na(NIRsoil$CEC), ]
  
  my_control <- mbl_control(validation_type = "NNv")
  
  ## The neighborhood sizes to test
  ks <- seq(40, 140, by = 20)
  
  ext_d <- dissimilarity(
    rbind(train_x, test_x),  Xu = rbind(train_x, test_x),
    diss_method = "cor",
    center = FALSE, scale = FALSE
  )$dissimilarity
  
  dim(ext_d)
  diag(ext_d) <- 0
  
  
  sbl_external_diss <- mbl(
    Xr = train_x,
    Yr = train_y,
    Xu = test_x,
    k = ks,
    spike = 1:5,
    method = local_fit_gpr(),
    diss_method = ext_d,
    diss_usage = "predictors",
    control = my_control,
    scale = FALSE, 
    center = FALSE
  )
  
  sbl_internal_diss <- mbl(
    Xr = train_x,
    Yr = train_y,
    Xu = test_x,
    k = ks,
    spike = 1:5,
    method = local_fit_gpr(),
    diss_method = "cor",
    diss_usage = "predictors",
    control = my_control,
    scale = FALSE, 
    center = FALSE
  )
  
  r_ext <- sbl_internal_diss$validation_results$nearest_neighbor_validation 
  r_int <- sbl_external_diss$validation_results$nearest_neighbor_validation
  
  expect_true(sum(abs(r_ext - r_int)) < tol)
})
