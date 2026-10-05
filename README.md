# NMSDE

<!-- badges: start -->
<!-- badges: end -->

**NMSDE** provides nonparametric estimation of spectral density matrices for
multivariate stationary time series. It implements:

- **NMSDE**, which estimates the spectral matrix of a single multivariate
  series (Nezampour et al., 2024,
  [doi:10.1007/s00180-023-01451-4](https://doi.org/10.1007/s00180-023-01451-4)).
- **NMCSDE**, which estimates the spectral matrices of a whole collection of
  multivariate series together, borrowing strength across series and
  clustering them in the frequency domain. It extends the collective
  approach of Maadooliat et al. (2018,
  [doi:10.1002/sim.7972](https://doi.org/10.1002/sim.7972)) to multivariate
  series.

The inverse spectral matrix of each series is parameterized by its complex
Cholesky factor. The Cholesky components are expanded in a B-spline basis
and lie in a common low-dimensional subspace. They are estimated by
minimizing a penalized Whittle likelihood with Newton–Raphson iterations
written in C++ (RcppArmadillo). The resulting estimates are always Hermitian
and positive definite. The smoothing parameter can be updated within the
iterations, which avoids a grid search, and spectral bootstrap confidence
intervals are available.

## Installation

Install the development version from GitHub with:

``` r
# install.packages("remotes")
remotes::install_github("shnezam/NMSDE")
```

## Example

``` r
library(NMSDE)

X <- dat$X2                       # 15 simulated bivariate series of length 200
init <- initial(X, nbasis = 10, R = 3, tru = datPar$modelPar2)
plot(init$elbow, type = "b")      # elbow method for the rank R

fit <- iterate(init, lambda = 1)  # lambda is updated within the iterations

# spectral estimates of all series: P^2 x K x M array
G <- tA2G(fit$thetas[[fit$r.iter]], fit$As[[fit$r.iter]], init)
spec <- array(apply(G, 3, inChol2spec, init = init),
              dim = c(init$P^2, init$K, init$m))

plot_sds(fit, init)               # estimates (black) and true spectra (blue)
```

See `vignette("NMSDE")` for a complete walk-through, which covers the
common- and different-bases models, noncollective estimates, clustering and
the bootstrap. To use the interactive Shiny application, run:

``` r
launch_app()
```
