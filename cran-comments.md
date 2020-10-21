# resemble
 
# Rhub checks for release of `resemble 2.0.0` (`gordillo`)

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
rhub::check("C:/Users/raml/Documents/GitHub/prospectr_0.2.1.tar.gz", 
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
