context("test-search_neighbors")

test_that("search_neighbors works", {
  nirdata <- data("NIRsoil", package = "prospectr")
  
  Xu <- NIRsoil$spc[!as.logical(NIRsoil$train), ]
  Yu <- NIRsoil$CEC[!as.logical(NIRsoil$train)]
  
  Yr <- NIRsoil$CEC[as.logical(NIRsoil$train)]
  Xr <- NIRsoil$spc[as.logical(NIRsoil$train), ]
  
  Xu <- Xu[!is.na(Yu), ]
  Xr <- Xr[!is.na(Yr), ]
  
  Yu <- Yu[!is.na(Yu)]
  Yr <- Yr[!is.na(Yr)]
  
  spikes <- c(5, 10, 15)
  
  nn <- search_neighbors(Xr = Xr, Xu = Xu, 
                         diss_method = c("pca.nipals"), 
                         Yr = Yr, k = 50,
                         spike = spikes,
                         pc_selection = list("manual", 10),
                         center = TRUE, scale = FALSE)
  
  nn_k_diss <- search_neighbors(Xr = Xr, Xu = Xu, 
                         diss_method = c("pca.nipals"), 
                         Yr = Yr, k_diss = 0.1, k_range = c(10, 50),
                         spike = spikes,
                         pc_selection = list("manual", 10),
                         center = TRUE, scale = FALSE)
  
  first_nns_pca <- search_neighbors(Xr = Xr, Xu = Xu, 
                                diss_method = c("pca"), 
                                Yr = Yr, k = 10,
                                pc_selection = list("manual", 10),
                                center = TRUE, scale = FALSE)
  expected_sum_of_indices_pca <- 17853
  
  first_nns_pls <- search_neighbors(Xr = Xr, Xu = Xu, 
                                    diss_method = c("pls"), 
                                    Yr = Yr, k = 10,
                                    pc_selection = list("manual", 10),
                                    center = TRUE, scale = FALSE)
  expected_sum_of_indices_pls <- 18410
  
  first_nns_cor <- search_neighbors(Xr = Xr, Xu = Xu, 
                                    diss_method = c("cor"), 
                                    Yr = Yr, k_diss = 0.1, k_range = c(10, 50),
                                    pc_selection = list("manual", 10),
                                    center = TRUE, scale = FALSE)
  expected_sum_of_indices_cor <- 18036
  
  
  first_nns_euclid <- search_neighbors(Xr = Xr, Xu = Xu, 
                                    diss_method = c("euclid"), 
                                    Yr = Yr, k_diss = 0.1, k_range = c(10, 50),
                                    pc_selection = list("manual", 10),
                                    center = TRUE, scale = FALSE)
  expected_sum_of_indices_euclid <- 19586
  
  first_nns_cosine <- search_neighbors(Xr = Xr, Xu = Xu, 
                                       diss_method = c("cosine"), 
                                       Yr = Yr, k_diss = 0.1, k_range = c(10, 50),
                                       pc_selection = list("manual", 10),
                                       center = TRUE, scale = FALSE)
  expected_sum_of_indices_cosine <- 18910
  
  expect_is(nn, "list")
  expect_is(nn_k_diss, "list")
  expect_true(all(rowMeans(nn$neighbors[1:length(spikes),]) == spikes))
  expect_true(all(rowMeans(nn_k_diss$neighbors[1:length(spikes),]) == spikes))
  expect_true(sum(first_nns_pca$neighbors[1,]) == expected_sum_of_indices_pca)
  expect_true(sum(first_nns_pls$neighbors[1,]) == expected_sum_of_indices_pls)
  expect_true(sum(first_nns_cor$neighbors[1,]) == expected_sum_of_indices_cor)
  expect_true(sum(first_nns_euclid$neighbors[1,]) == expected_sum_of_indices_euclid)
  expect_true(sum(first_nns_cosine$neighbors[1,]) == expected_sum_of_indices_cosine)
})
