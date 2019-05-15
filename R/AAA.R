# .RESEMBLE_CACHE <- new.env(FALSE, parent = globalenv())

.onAttach <- function(lib, pkg) {
    #assign("gpclib", FALSE, envir=.RESEMBLE_CACHE)
    resemble.v <- read.dcf(file = system.file("DESCRIPTION", package=pkg), fields="Version")
    packageStartupMessage(paste(pkg, "version", resemble.v, "-- 'milagros'"))
    packageStartupMessage("check the package website at http://l-ramirez-lopez.github.io/resemble/")
}

# .onUnload <- function(libpath) {
#     rm(.RESEMBLE_CACHE)
# }
