# to comply with CRAN's 2-core limit
Sys.setenv(OMP_NUM_THREADS = "2")
Sys.setenv(OMP_THREAD_LIMIT = "2")
options(mc.cores = 2)
