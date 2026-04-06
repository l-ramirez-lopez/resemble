.onAttach <- function(lib, pkg) {
  pkg_v <- pkg_info()
  
  mss <- paste0(
    "\033[34m",
    pkg, " version ",
    paste(pkg_v[1:2], collapse = " \U002D\U002D "),
    "\033[39m"
  )
  
  mss2 <- paste0(
    "\033[34mRepository: ",
    pkg_v[, "URL"],
    "\033[39m"
  )
  
  mss3 <- paste0(
    "\033[33m",
    "Note: Version 3.0 has deprecated several functions and arguments.\n",
    "For previous functionality, install version 2.2.5:\n",
    "  remotes::install_github(\"l-ramirez-lopez/resemble@2.2.5\")",
    "\033[39m"
  )
  
  packageStartupMessage(mss)
  packageStartupMessage(mss2)
  packageStartupMessage(mss3)
}
