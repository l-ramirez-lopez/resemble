## resemble 3.0.0

## R CMD check results

0 errors | 0 warnings | 0 notes

## Test environments

- Local Ubuntu 24.04, R 4.5.3
- win-builder R-devel (Windows Server 2022, R 4.6.0 beta)
- Debian Linux R-devel (clang 21.1.8)
- GitHub Actions (ubuntu, windows, macos)

## Notes

Checktime on Windows is ~14 min due to rebuilding 8 Quarto vignettes. Vignettes are pre-built and included in the package. Added `skip_on_cran()` to 119 computationally intensive tests.

## Downstream dependencies

No reverse dependencies.