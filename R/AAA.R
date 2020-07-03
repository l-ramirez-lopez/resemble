# .RESEMBLE_CACHE <- new.env(FALSE, parent = globalenv())

.onAttach <- function(lib, pkg) {
  # assign("gpclib", FALSE, envir=.RESEMBLE_CACHE)
  resemble_v <- read.dcf(file = system.file("DESCRIPTION", package = pkg), fields = "Version")
  packageStartupMessage(paste0("\033[34m", pkg, " version ", resemble_v, " -- 'gordillo'\033[39m"))
  packageStartupMessage("\033[34mcheck the package repository at https://github.com/l-ramirez-lopez/resemble/\033[39m")
}

# .onUnload <- function(libpath) {
#     rm(.RESEMBLE_CACHE)
# }
