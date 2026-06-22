# Project Instructions

## Git and Branch Rules

- **Never create a new branch without explicit instruction from the user.**
- **Never use "claude" anywhere** in branch names, commit messages, PR titles, code comments, or any other artifact pushed to the repository.
- Commit messages must read as normal developer work — no AI attribution, no session references.

## Known CRAN Check Notes (r-devel-linux-x86_64-debian-gcc)

These notes are known and tracked:

- **Check: tests** — `Running 'testthat.R'` had CPU time 2.6 times elapsed time
- **Check: re-building of vignette outputs** — Re-building vignettes had CPU time 4 times elapsed time

These are parallelism-related timing notes on CRAN's Debian GCC build. Keep in mind when working on tests or vignettes.
