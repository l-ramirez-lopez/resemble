# .RESEMBLE_CACHE <- new.env(FALSE, parent = globalenv())

.onAttach <- function(lib, pkg) {
  # assign("gpclib", FALSE, envir=.RESEMBLE_CACHE)
  resemble_v <- read.dcf(
    file = system.file("DESCRIPTION", package = pkg),
    fields = "Version"
  )
  mss <- paste0(
    "\033[34m",
    pkg, " version ",
    resemble_v,
    " -- 'piapia'\033[39m"
  )
  mss2 <- paste0(
    "\033[34mcheck the package repository at: ",
    "https://github.com/l-ramirez-lopez/resemble/\033[39m"
  )
  packageStartupMessage(mss)
  packageStartupMessage(mss2)
}

# .onUnload <- function(libpath) {
#     rm(.RESEMBLE_CACHE)
# }
