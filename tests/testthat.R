Sys.setenv(OMP_NUM_THREADS = 2)
library(testthat)
library(prospectr)
library(resemble)

test_check("resemble")
