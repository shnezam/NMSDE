# NMSDE 0.1.0

* First CRAN submission.
* Documented all exported functions and the `dat` and `datPar` data sets;
  added a vignette and tests.
* Removed the dependency on the archived package 'dse'. The Shiny app now
  simulates VAR series with base R.
* `initial()`: `R` must now be supplied (the previous default referred to an
  object of the Shiny app); the difference penalty now uses order `a`
  instead of always 2; the derivative penalty now uses the same knots as the
  basis and the derivative order `a`; true spectra are now computed from the
  `paras` and `sigma` components of `tru` (as stored in `datPar`); input
  checks were added; the result gains `elbow`, `pen` and `a`.
* `iterate()`: new argument `verbose`; `a` now defaults to the value used in
  `initial()`.
* `NMSDEBoot()`: fixed an error (undefined `trueF`) with
  `method = "separate"`; new argument `verbose`.
* The step-halving loop is now capped at 30 halvings.
* The C++ code now compiles with C++17/C++20 compilers (ambiguous `span`).
* The C++ code now finds the R function `loglike()` in the package namespace,
  so the common-bases model works without attaching the package.
* `plot_sds()` restores the graphical parameters on exit.
* Fixed the package unloading hook, and the Shiny app no longer writes to
  the global environment.
