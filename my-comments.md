# resemble
 
# Rhub checks for release of `resemble 2.0.0` (`gordillo`)

29.10.2020
As requested by CRAN:
-The length of the title is now below
  65 characters
- A <doi:...> has been added in the description field of DESCRIPTION
- \donttest{} is now used (instead of \dontrun{}) for those examples
  taking more than 5 seconds 
- verobse argument has been added to the functions to easily suppress any
  message different from error warnings or messages.
- on.exit() is now called properly to reset to user
  parameters when the functions are exited
- User's options() are reset in the examples in the vignette
  that require changes in those options()
- Examples are now explicitly using maximum two cores (if
  available).

The package has been checked in multiple platforms using rhub::check(), as well 
as in the win-builders provided by CRAN.



21.10.2020

The checks were conducted in the following platforms trhough rhub:

- "debian-clang-devel"

- "debian-gcc-devel"

- "fedora-gcc-devel"

- "debian-gcc-devel-nold"

- "debian-gcc-patched"

- "debian-gcc-release"

- "linux-x86_64-rocker-gcc-san" 

- "macos-highsierra-release-cran" 

- "solaris-x86-patched-ods" 

- "ubuntu-gcc-release"

- "windows-x86_64-devel"

For example, for checks with "fedora-gcc-devel"", the following code was used::
```
rhub::check("./resemble_2.0.0.tar.gz", 
            platform = c("fedora-gcc-devel"), 
            email = "ramirez.lopez.leo@gmail.com")
```

Some of the unit tests for  `pls_projection` and `pc_projection` were failing
in three platfroms apparently due to numeric accuracy and the use of OPENMP in 
Rcpp. These platforms were:
 
 - "debian-gcc-devel"
 
 - "debian-gcc-devel-nold"
 
 - "linux-x86_64-rocker-gcc-san" 
 
The remaining platforms were passing all the tests sucessfully. For the above 
three platforms, the solution was to disable OPENMP.

All platforms pass the checks successfully for the release. 

## Size of the package
To reduce the size of the package, Makevars was modified and Makevars.win was 
added.

## NOTE for compiled code
An strange not was thrown when the check was done locally on windows with R `4.0.3`
It is apparently a problem in R core and not related to the package nor Rcpp. 
The issue was reported here:
https://stackoverflow.com/questions/64402688/information-on-o-files-for-x64-is-not-available-note-on-r-package-checks-using/64419033#64419033


