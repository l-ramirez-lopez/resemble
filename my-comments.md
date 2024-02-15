# resemble

# version 2.2.3

# submission message:
Dear CRAN maintainers,
I am submitting my package "resemble" to CRAN. This version introduces a new 
feature in one of the functions and one bug fix.
In previous submissions, the date in DESCRIPTION was not in the correct format and there were some broken URL links. These issues have been fixed.
Prior to this submission, this tarball has been checked with in the winbuilder service. Apart from that it has been also submitted to extensive tests in rhub.
A first submission of this version failed (for "r-devel-linux-x86_64-debian-gcc"), 
therefore following platforms were tested for a second submission using Rhub: 
- Debian Linux, R-devel, GCC ASAN/UBSAN
- Debian Linux, R-devel, GCC, no long double
- Debian Linux, R-devel, clang, ISO-8859-15 locale
- Debian Linux, R-devel, GCC
For this second submission the package passed all the tests in the above platforms. 
Reverse dependencies have also been checked. 
Best regards,
Leonardo



## Package was built using: 
```
devtools::build(
  pkg = ".",
  path = NULL,
  binary = FALSE,
  vignettes = TRUE,
  manual = TRUE,
  args = NULL,
  quiet = FALSE
)
```

# R win builder checks for release of `resemble 2.2.3` (`embryo`) 19.04.2023 
passed all the checks without notes.

# Rhub checks for release of `resemble 2.2.2` (`Sky`) 19.04.2023
The checks were conducted in the following platforms through rhub:

```
rhub::check(paste0(gsub("/resemble$", "/", getwd()), "resemble_2.2.2.tar.gz"), 
            platform = c("fedora-gcc-devel"), 
            email = "ramirez.lopez.leo@gmail.com")
```
- "linux-x86_64-rocker-gcc-san" Test cannot be executed as the docker failed 
the error was as follows: ‘RcppArmadillo’ is not available for package ‘resemble’

- "fedora-gcc-devel" NOTE
* checking installed package size ... NOTE
  installed size is  9.8Mb
  sub-directories of 1Mb or more:
    doc    1.6Mb
    libs   7.5Mb

- "windows-x86_64-devel" OK

- "macos-highsierra-release-cran" OK

- "windows-x86_64-release" OK

- "ubuntu-gcc-release" NOTE
* checking installed package size ... NOTE
  installed size is 14.4Mb
  sub-directories of 1Mb or more:
    doc    1.6Mb
    libs  12.2Mb

- "solaris-x86-patched-ods" 
* checking package dependencies ... ERROR
Packages suggested but not available: 'testthat', 'rmarkdown', 'bookdown'

The suggested packages are required for a complete check.
Checking can be attempted without them by setting the environment
variable _R_CHECK_FORCE_SUGGESTS_ to a false value.

See section ‘The DESCRIPTION file’ in the ‘Writing R Extensions’
manual.





# version 2.2.2

# submission message:
Dear CRAN maintainers,
I am submitting my package "resemble" to CRAN. This version introduces a new 
feature in one of the functions and one bug fix.
In previous submissions, the date in DESCRIPTION was not in the correct format and there were some broken URL links. These issues have been fixed.
Prior to this submission, this tarball has been checked with in the winbuilder service. Apart from that it has been also submitted to extensive tests in rhub.
A first submission of this version failed (for "r-devel-linux-x86_64-debian-gcc"), 
therefore following platforms were tested for a second submission using Rhub: 
- Debian Linux, R-devel, GCC ASAN/UBSAN
- Debian Linux, R-devel, GCC, no long double
- Debian Linux, R-devel, clang, ISO-8859-15 locale
- Debian Linux, R-devel, GCC
For this second submission the package passed all the tests in the above platforms. 
Reverse dependencies have also been checked. 
Best regards,
Leonardo



## Package was built using: 
```
devtools::build(
  pkg = ".",
  path = NULL,
  binary = FALSE,
  vignettes = TRUE,
  manual = TRUE,
  args = NULL,
  quiet = FALSE
)
```

# R win builder checks for release of `resemble 2.2.2` (`Sky`) 19.04.2023 
passed all the checks without notes.

# Rhub checks for release of `resemble 2.2.2` (`Sky`) 19.04.2023
The checks were conducted in the following platforms through rhub:

```
rhub::check(paste0(gsub("/resemble$", "/", getwd()), "resemble_2.2.2.tar.gz"), 
            platform = c("fedora-gcc-devel"), 
            email = "ramirez.lopez.leo@gmail.com")
```
- "linux-x86_64-rocker-gcc-san" Test cannot be executed as the docker failed 
the error was as follows: ‘RcppArmadillo’ is not available for package ‘resemble’

- "fedora-gcc-devel" NOTE
* checking installed package size ... NOTE
  installed size is  9.8Mb
  sub-directories of 1Mb or more:
    doc    1.6Mb
    libs   7.5Mb

- "windows-x86_64-devel" OK

- "macos-highsierra-release-cran" OK

- "windows-x86_64-release" OK

- "ubuntu-gcc-release" NOTE
* checking installed package size ... NOTE
  installed size is 14.4Mb
  sub-directories of 1Mb or more:
    doc    1.6Mb
    libs  12.2Mb

- "solaris-x86-patched-ods" 
* checking package dependencies ... ERROR
Packages suggested but not available: 'testthat', 'rmarkdown', 'bookdown'

The suggested packages are required for a complete check.
Checking can be attempted without them by setting the environment
variable _R_CHECK_FORCE_SUGGESTS_ to a false value.

See section ‘The DESCRIPTION file’ in the ‘Writing R Extensions’
manual.



# version 2.2.1

# submission message:
Dear CRAN maintainers,
I am submitting my package "resemble" to CRAN. This version accounts for 
problems found in Rd files auto-generated with roxygen2 7.1.2 (not compatible 
with HTML5). The new Rd files are now compatible with HTML5 (as Rd files 
are generated with roxygen2_7.2.0 ). 
Prior to this submission, this tarball has been checked with in the winbuilder service. Apart from that it has been also submitted to extensive tests in rhub.
A first submission of this version failed (for "r-devel-linux-x86_64-debian-gcc"), 
therefore following platforms were tested for a second submission using Rhub: 
- Debian Linux, R-devel, GCC ASAN/UBSAN
- Debian Linux, R-devel, GCC, no long double
- Debian Linux, R-devel, clang, ISO-8859-15 locale
- Debian Linux, R-devel, GCC
For this second submission the package passed all the tests in the above platforms. 
Reverse dependencies have also been checked. 
Best regards,
Leonardo



## Package was built using: 
```
devtools::build(
  pkg = ".",
  path = NULL,
  binary = FALSE,
  vignettes = TRUE,
  manual = TRUE,
  args = NULL,
  quiet = FALSE
)
```

# R win builder checks for release of `resemble 2.2.1` (`Fix-Hodges`) 30.08.2022 
passed all the checks without notes.

# Rhub checks for release of `resemble 2.2.1` (`Fix-Hodges`) 30.08.2022
The checks were conducted in the following platforms through rhub:

```
rhub::check(paste0(gsub("/resemble$", "/", getwd()), "resemble_2.2.1.tar.gz"), 
            platform = c("fedora-gcc-devel"), 
            email = "ramirez.lopez.leo@gmail.com")
```
- "linux-x86_64-rocker-gcc-san" OK

- "fedora-gcc-devel" NOTE
  installed size is 11.7Mb
  sub-directories of 1Mb or more:
    doc    1.6Mb
    libs   9.5Mb

- "windows-x86_64-devel" OK

- "macos-highsierra-release-cran" OK

- "windows-x86_64-release" OK

- "ubuntu-gcc-release" NOTE
     installed size is 13.5Mb
     sub-directories of 1Mb or more:
       doc    1.6Mb
       libs  11.3Mb


- "solaris-x86-patched-ods" Package suggested but not available: ‘testthat’
   
   The suggested packages are required for a complete check.
   Checking can be attempted without them by setting the environment
   variable _R_CHECK_FORCE_SUGGESTS_ to a false value.
   
   See section ‘The DESCRIPTION file’ in the ‘Writing R Extensions’
   manual.



# version 2.1.1

## Package was built using: 
```
devtools::build(
  pkg = ".",
  path = NULL,
  binary = FALSE,
  vignettes = TRUE,
  manual = TRUE,
  args = NULL,
  quiet = FALSE
)
```

26.11.2021
## Size of the package
According to Brian Ripley the package was violating CRAN policies
The CRAN policy contains

- Packages should not attempt to disable compiler diagnostics, nor to
remove other diagnostic information such as symbols in shared objects.

The package was stripping some symbols for Rcpp functions in Makevars in order 
to reduce the installation size of the package. Now these lines have been 
commented to comply with CRAN policies:
#strippedLib: $(SHLIB)
#		if test -e "/usr/bin/strip" & test -e "/bin/uname" & [[ `uname` == "Linux" ]]; then /usr/bin/strip --strip-debug $(SHLIB); fi
#.phony: strippedLib

# Rhub checks for release of `resemble 2.1.1` (`piapia`)

The checks were conducted in the following platforms trhough rhub:

- "windows-x86_64-devel"

- "macos-highsierra-release-cran"

- "linux-x86_64-rocker-gcc-san" ## Not checked Rhub throwed a PREPERROR

- "debian-clang-devel" OK

- "fedora-gcc-devel"
* checking installed package size ... NOTE
installed size is 12.2Mb
sub-directories of 1Mb or more:
  doc    2.0Mb
libs   9.6Mb

- "solaris-x86-patched-ods" OK

For example, for checks with "fedora-gcc-devel"", the following code was used::
```
rhub::check("/home/rl_leonardo/github/resemble_2.1.1.tar.gz", 
            platform = c("fedora-gcc-devel"), 
            email = "ramirez.lopez.leo@gmail.com")
```



# version 2.0 
# Rhub checks for release of `resemble 2.0.0` (`gordillo`)

29.10.2020
As requested by CRAN:
-The length of the title is now below
  65 characters
- A <doi:...> has been added in the description field of DESCRIPTION
- \donttest{} is now used (instead of \dontrun{}) for those examples
  taking more than 5 seconds 
- verbose argument has been added to the functions to easily suppress any
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


