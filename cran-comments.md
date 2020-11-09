# resemble
 
# checks for release of `resemble 2.0.0` (`gordillo`)

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