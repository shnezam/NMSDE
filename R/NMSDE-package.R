#' NMSDE: Nonparametric Multivariate Spectral Density Estimation
#'
#' Nonparametric estimation of the spectral density matrices of a collection
#' of \eqn{M} stationary \eqn{P}-variate time series (NMCSDE, nonparametric
#' multivariate collective spectral density estimation).
#'
#' For each series \eqn{i}, the inverse spectral matrix is factored as
#' \eqn{f_i^{-1}(\omega) = G_i(\omega) G_i^*(\omega)}, where \eqn{G_i} is a
#' complex lower-triangular Cholesky factor. Its \eqn{P^2} real and imaginary
#' components \eqn{g_j^{(i)}} are evaluated at the Fourier frequencies
#' \eqn{\omega_k = k/T}, \eqn{k = 1, \dots, K}, and modelled as
#' \deqn{\Gamma_j = B \Theta_j A_j^\top, \quad j = 1, \dots, P^2,}{Gamma_j = B Theta_j t(A_j),  j = 1, ..., P^2,}
#' where \eqn{B} (\eqn{K \times L}{K x L}) holds a rich B-spline basis,
#' \eqn{\Theta_j} (\eqn{L \times R}{L x R}) defines \eqn{R \ll L} data-driven
#' basis functions and \eqn{A_j} (\eqn{M \times R}{M x R}) holds the scores of
#' the \eqn{M} series. Two variants are available:
#' \describe{
#'   \item{different bases}{a separate \eqn{\Theta_j} for every Cholesky
#'     component (\code{theta_P2 = TRUE} in \code{\link{initial}});}
#'   \item{common bases}{one \eqn{\Theta} shared by all components
#'     (\code{theta_P2 = FALSE}).}
#' }
#' The parameters minimize the penalized Whittle negative log-likelihood
#' \deqn{\sum_{i=1}^M \sum_{k=1}^K \left\{-\log|G_i G_i^*| + Y_{i,k}^* G_i G_i^* Y_{i,k}\right\} + \lambda\, \mathrm{tr}(\Theta^\top D \Theta),}{sum_i sum_k { -log|G_i G_i^*| + Y_ik^* G_i G_i^* Y_ik } + lambda tr(t(Theta) D Theta),}
#' where \eqn{Y_{i,k}} is the discrete Fourier transform of series \eqn{i} at
#' \eqn{\omega_k} and \eqn{D} is a roughness penalty matrix. Estimation uses
#' Newton-Raphson (Fisher scoring) iterations with step halving, and the
#' smoothing parameter \eqn{\lambda} may be updated within the iterations.
#' Because the estimate is rebuilt from a Cholesky factor it is always
#' Hermitian and positive definite. With \eqn{M = 1}, or by fitting each series
#' separately, the method reduces to the noncollective NMSDE estimator of
#' Nezampour et al. (2024).
#'
#' @section Typical workflow:
#' \enumerate{
#'   \item \code{\link{initial}} computes the DFTs, periodograms, B-spline
#'     basis, penalty matrix and starting values.
#'   \item \code{\link{iterate}} runs the penalized Whittle likelihood
#'     minimization.
#'   \item \code{\link{tA2G}} and \code{\link{inChol2spec}} turn the fitted
#'     parameters into spectral density estimates; \code{\link{plot_sds}}
#'     plots them.
#'   \item \code{\link{NMSDEBoot}} gives spectral bootstrap replicates for
#'     pointwise confidence intervals.
#' }
#' \code{\link{launch_app}} starts an interactive 'shiny' application with the
#' same steps, plus clustering of the series.
#'
#' @section Vectorized spectral and Cholesky matrices:
#' Throughout the package a \eqn{P \times P}{P x P} Hermitian spectral matrix,
#' or a lower-triangular Cholesky factor, at one frequency is stored as a
#' real vector of length \eqn{P^2}: first the real parts of the
#' \eqn{P(P+1)/2} lower-triangular entries (column by column), then the
#' imaginary parts of the \eqn{P(P-1)/2} strictly lower-triangular entries.
#' The ordering is given by \code{\link{Cholindexer}}. Over \eqn{K}
#' frequencies these vectors form \eqn{P^2 \times K}{P^2 x K} matrices.
#'
#' @references
#' Nezampour, S., Nematollahi, A., Krafty, R. T. and Maadooliat, M. (2024).
#' A new approach to nonparametric estimation of multivariate spectral density
#' function using basis expansion. \emph{Computational Statistics}, 39,
#' 3625--3641. \doi{10.1007/s00180-023-01451-4}
#'
#' Maadooliat, M., Sun, Y. and Chen, T. (2018). Nonparametric collective
#' spectral density estimation with an application to clustering the brain
#' signals. \emph{Statistics in Medicine}, 37(30), 4789--4806.
#' \doi{10.1002/sim.7972}
#'
#' Nezampour, S., Krafty, R. T. and Maadooliat, M. Nonparametric collective
#' spectral density estimation of multiple multivariate time series.
#' Manuscript.
#'
#' @keywords internal
#' @import shiny
#' @importFrom graphics axis lines pairs par plot points title
#' @importFrom stats as.dendrogram cmdscale cutree dist fft hclust plot.ts rnorm ts.plot
#' @importFrom utils data head
#' @importFrom Rcpp sourceCpp
#' @useDynLib NMSDE, .registration = TRUE
"_PACKAGE"
