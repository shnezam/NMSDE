# initial -----------------------------------------------------------------

#' Initialize the NMCSDE algorithm
#'
#' Sets up everything needed to estimate the spectral density matrices of a
#' collection of multivariate time series: the discrete Fourier transforms
#' and periodograms, the B-spline basis and roughness penalty, the
#' Nadaraya-Watson starting estimate, and the starting values of the
#' low-rank parameters \eqn{(\Theta, A)}.
#'
#' The Nadaraya-Watson estimate of each spectral matrix (\code{\link{NWspec}})
#' is inverted and factored by a complex Cholesky decomposition
#' (\code{\link{spec2inChol}}). Its components are projected onto the B-spline
#' basis, and a truncated singular value decomposition of the projected
#' coefficients gives the starting values of \eqn{\Theta} and \eqn{A}. The
#' columns of \eqn{\Theta} are orthonormal and those of \eqn{A} are ordered by
#' decreasing singular value, which makes the factorization identifiable.
#'
#' @param X A list of \eqn{M} numeric matrices, one per time series. Each
#'   matrix is \eqn{P \times T}{P x T}: rows are the \eqn{P \ge 2} component
#'   series and columns are time points. All series must have the same
#'   dimensions. The series are assumed to have mean zero.
#' @param nbasis Number of B-spline basis functions \eqn{L}.
#' @param norder Order of the B-splines (degree + 1); \code{4} gives cubic
#'   splines.
#' @param R Number of data-driven basis functions (the rank of the collective
#'   model). It is closely related to the number of clusters among the series;
#'   \code{init$elbow} can help to choose it. Must not exceed \eqn{M} when
#'   \code{theta_P2 = TRUE}.
#' @param pen Type of roughness penalty matrix \eqn{D}: \code{"diff"} for the
#'   difference penalty \eqn{D = L_a^\top L_a}{D = t(L_a) L_a} built from
#'   \eqn{a}-th order differences of adjacent coefficients, or any other value
#'   (e.g. \code{"deriv"}) for the integrated squared \eqn{a}-th derivative of
#'   the B-splines.
#' @param a Order of the differences (or derivative) used in the penalty.
#'   With \code{pen = "diff"}, \code{a = 0} gives a ridge penalty
#'   \eqn{D = I}. The value is stored in the result and used by
#'   \code{\link{iterate}} when updating \eqn{\lambda}.
#' @param theta_P2 Logical. \code{TRUE} (default) fits the \emph{different
#'   bases} model, with one coefficient matrix \eqn{\Theta_j} per Cholesky
#'   component; \code{FALSE} fits the \emph{common bases} model, with one
#'   \eqn{\Theta} shared by all \eqn{P^2} components.
#' @param tru Optional list describing the true vector autoregressive models
#'   of simulated data, used only to compute the true spectra for comparison.
#'   It has components \code{paras} (a list with one element per model, each a
#'   list of \eqn{P \times P}{P x P} lag-polynomial matrices
#'   \eqn{I, -\Phi_1, \dots, -\Phi_d}), \code{sigma} (innovation covariance
#'   matrix) and \code{rnd} (the model index of each series). See
#'   \code{\link{datPar}} for an example.
#' @param f.r Optional integer vector \code{c(first, last)} restricting the
#'   analysis to Fourier frequencies \eqn{\omega_k = k/T} with
#'   \code{first <= k <= last}. By default all \eqn{K = \lfloor (T-1)/2 \rfloor}{K = floor((T-1)/2)}
#'   frequencies are used.
#'
#' @return A list (referred to as \code{init} elsewhere in this package) with
#'   components
#'   \item{P, m, K, R}{the dimension of the series, the number of series
#'     \eqn{M}, the number of frequencies and the rank.}
#'   \item{omegas}{the \eqn{K} frequencies, in cycles per unit time.}
#'   \item{f.r}{the frequency index range used.}
#'   \item{B, knots}{the \eqn{K \times L}{K x L} B-spline basis matrix and its
#'     knots.}
#'   \item{D, pen, a}{the \eqn{L \times L}{L x L} penalty matrix and the
#'     arguments used to build it.}
#'   \item{indexer, Nreal}{output of \code{\link{Cholindexer}} and the number
#'     \eqn{P(P+1)/2} of real components.}
#'   \item{drivG}{derivatives used by the Newton-Raphson steps (see
#'     \code{\link{drivG}}).}
#'   \item{periods, tildeX}{lists with the vectorized periodograms
#'     (\eqn{P^2 \times K}{P^2 x K}) and DFTs (\eqn{P \times K}{P x K}) of each
#'     series (see \code{\link{time2freq}}).}
#'   \item{NW}{list of Nadaraya-Watson estimates (see \code{\link{NWspec}}).}
#'   \item{thetas, As}{starting values: \code{thetas} is an
#'     \eqn{L \times R \times P^2}{L x R x P^2} array (different bases) or an
#'     \eqn{L \times R}{L x R} matrix (common bases); \code{As} is an
#'     \eqn{M \times R \times P^2}{M x R x P^2} array.}
#'   \item{theta.tA}{\eqn{L \times M \times P^2}{L x M x P^2} array of basis
#'     coefficients of the projected starting estimate.}
#'   \item{elbow}{total within-cluster sum of squares from hierarchical
#'     clustering (Ward's method) of the starting estimates into
#'     \eqn{1, 2, \dots} clusters; plot it to choose \code{R} by the elbow
#'     method.}
#'   \item{truspec, rnd}{only when \code{tru} is supplied: a
#'     \eqn{P^2 \times K \times}{P^2 x K x} (number of models) array of true
#'     spectra (see \code{\link{psd}}) and the model index of each series.}
#'
#' @references Nezampour, S., Nematollahi, A., Krafty, R. T. and Maadooliat,
#'   M. (2024). A new approach to nonparametric estimation of multivariate
#'   spectral density function using basis expansion. \emph{Computational
#'   Statistics}, 39, 3625--3641. \doi{10.1007/s00180-023-01451-4}
#' @seealso \code{\link{iterate}}
#' @examples
#' # six bivariate series of length 200 (see ?dat)
#' X <- dat$X1
#' init <- initial(X, nbasis = 15, R = 3, tru = datPar$modelPar1)
#' dim(init$thetas)
#' plot(init$elbow, type = "b", xlab = "Number of clusters",
#'      ylab = "Within-cluster sum of squares")
#' @export initial
initial <- function(X, nbasis = 10, norder = 4, R,
                    pen = "diff", a = 2, theta_P2 = TRUE, tru = NULL, f.r = NULL) {
  # set up some basic values
  P <- nrow(X[[1]])
  if (is.null(P) || P < 2) stop("each element of 'X' must be a P x T matrix with P >= 2")
  P2 <- P^2
  Nreal <- P * (P + 1) / 2
  len <- ncol(X[[1]])
  m <- length(X)
  if (R > nbasis || (theta_P2 && R > m)) {
    stop("'R' must not exceed 'nbasis' or, when 'theta_P2' is TRUE, the number of series")
  }
  K <- floor((len - 1) / 2)
  omegas <- (1:K) / len
  indexer <- Cholindexer(P)
  if (is.null(f.r)) {
    f.r <- c(1, length(omegas))
  } else {
    omegas <- omegas[f.r[1]:f.r[2]]
    K <- length(omegas)
  }

  # Bases(B splines)[K*L] and Penalty matrix [L*L]
  knots <- seq(omegas[1] - (1 / len), omegas[length(omegas)] + (1 / len),
              length.out = nbasis - norder + 2)
  B <- fda::bsplineS(omegas, breaks = knots, norder = norder)

  # Penalty functions
  if (pen == "diff") {
    if (a != 0) {
      L <- diff(diag(nbasis), differences = a)
    } else {
      L <- diag(nbasis)
    }
    D <- t(L) %*% L
  } else {
    D <- fda::bsplinepen(fda::create.bspline.basis(rangeval = range(knots), breaks = knots,
                                                   norder = norder), Lfdobj = a)
  }

  init <- list(P = P, m = m, K = K, R= R, omegas = omegas, f.r = f.r, B = B, D = D,
               pen = pen, a = a, indexer = indexer, Nreal = Nreal, knots = knots)
  
  # some required derivatives for Fisher-scoring algorithm
  init$drivG <- drivG(init)
  
  # go to frequency domain from time domain (DFTs and Periodograms)
  freqdom <- lapply(X, time2freq, init)
  init$periods <- lapply(freqdom, function(element) element$perds)
  init$tildeX <- lapply(freqdom, function(element) element$tildeX)
  
  # Nadaraya-Watson kernel estimate of spectral matrix
  NWF <- lapply(init$periods, NWspec, init)
  init$NW <- NWF
  # Cholesky components of inverse of Nadaraya-Watson estimate
  NWG_m <- sapply(NWF, spec2inChol, init, simplify = "array")
  
  # get initial values by projecting the Cholesky components (NW) 
  #onto the linear subspace (of the basis)
  init <- append(init, G2tA(NWG_m, theta_P2, init))

  # within-cluster sum of squares of the initial estimates (elbow method for R)
  init$elbow <- elbow_wss(init)

  # get true spectrum
  if (!is.null(tru)) {
    if (is.null(tru$parasmdl)) {
      tru$parasmdl <- lapply(tru$paras, function(Phi) list(Phi, list(diag(P)), tru$sigma))
    }
    init$truspec <- sapply(tru$parasmdl, psd, init, simplify = "array")

    init$rnd <- tru$rnd
  }

  return(init)
}

# elbow_wss ---------------------------------------------------------------

# Total within-cluster sum of squares of Ward clusterings of the initial
# (projected Nadaraya-Watson) Cholesky components, for 1, ..., max.k clusters
elbow_wss <- function(init, max.k = 10) {
  feat <- do.call(cbind, lapply(seq_len(init$P^2), function(j) t(init$B %*% init$theta.tA[, , j])))
  if (nrow(feat) < 2) return(0)
  max.k <- min(max.k, nrow(feat))
  clust <- hclust(dist(feat), method = "ward.D2")
  sapply(seq_len(max.k), function(k) {
    cl <- cutree(clust, k = k)
    sum(sapply(unique(cl), function(g) {
      x <- feat[cl == g, , drop = FALSE]
      sum(sweep(x, 2, colMeans(x))^2)
    }))
  })
}

# time2freq ---------------------------------------------------------------

#' Discrete Fourier transform and periodogram of a multivariate time series
#'
#' Computes the scaled discrete Fourier transform
#' \eqn{Y_k = T^{-1/2} \sum_{t=1}^T X_t e^{-2\pi i \omega_k t}}{Y_k = T^(-1/2) sum_t X_t exp(-2 pi i omega_k t)}
#' and the periodogram matrices \eqn{S(\omega_k) = Y_k Y_k^*} at the
#' frequencies selected in \code{init}.
#'
#' @param x A \eqn{P \times T}{P x T} numeric matrix (rows are components).
#' @param init List returned by \code{\link{initial}} (only \code{P},
#'   \code{K}, \code{f.r}, \code{Nreal} and \code{indexer} are used).
#'
#' @return A list with
#'   \item{tildeX}{\eqn{P \times K}{P x K} complex matrix of DFTs.}
#'   \item{perds}{\eqn{P^2 \times K}{P^2 x K} matrix of vectorized
#'     periodograms (see \code{\link{Cholindexer}} for the ordering).}
#' @seealso \code{\link{NWspec}}
#' @examples
#' init <- initial(dat$X1, nbasis = 10, R = 2)
#' fd <- time2freq(dat$X1[[1]], init)
#' matplot(init$omegas, t(fd$perds), type = "l", lty = 1,
#'         xlab = "Frequency", ylab = "Periodogram")
#' @export time2freq
time2freq <- function(x, init) {
  
  # for tildeX
  tildex <- matrix(0, nrow = init$P, ncol = ncol(x))
  for (p in 1:init$P) tildex[p, ] <- (fft(x[p, ])) / sqrt(ncol(x))
  tildex <- tildex[, (init$f.r[1] + 1):(init$f.r[2] + 1)]

  # for perds
  perds <- matrix(0, nrow = init$P^2, ncol = init$K)
  for (k in 1:init$K) {
    P_mat <- tildex[, k] %*% (Conj(t(tildex[, k])))
    for (p in 1:init$Nreal) perds[p, k] <- Re(P_mat[init$indexer[p, 3], 
                                                    init$indexer[p, 4]])
    for (p in (init$Nreal + 1):init$P^2) perds[p, k] <- Im(P_mat[init$indexer[p, 3], 
                                                                 init$indexer[p, 4]])
  }
  
  return(list(tildeX = tildex, perds = perds))
}


# NWspec ------------------------------------------------------------------

#Nadaraya-Watson kernel regression estimate

#' Nadaraya-Watson kernel estimate of a spectral matrix
#'
#' Smooths the vectorized periodogram across frequencies with a Gaussian
#' kernel, using bandwidth \eqn{h = (20/K)^2}. The result is used to start the
#' Newton-Raphson iterations and as a simple competing estimator.
#'
#' @param perds \eqn{P^2 \times K}{P^2 x K} matrix of vectorized periodograms,
#'   as returned in \code{time2freq()$perds}.
#' @param init List returned by \code{\link{initial}} (only \code{omegas} is
#'   used).
#'
#' @return A \eqn{P^2 \times K}{P^2 x K} matrix with the vectorized smoothed
#'   spectral matrices.
#' @seealso \code{\link{time2freq}}, \code{\link{spec2inChol}}
#' @examples
#' init <- initial(dat$X1, nbasis = 10, R = 2)
#' nw <- NWspec(init$periods[[1]], init)
#' all.equal(nw, init$NW[[1]])
#' @export NWspec
NWspec <- function(perds,init) {
  P2 <- nrow(perds)
  K <- ncol(perds)
  logh <- 2 * log(20 / K)
  h <- exp(logh)
  hatF <- matrix(0, nrow = P2, ncol = K)
  # smooth the periodogram
  for (k in 1:K) {
    z <- kerf((init$omegas[k] - init$omegas) / h)
    hatF[, k] <- perds %*% z / sum(z)
  }
  
  return(hatF)
}


# Gaussian kernel-------------
kerf <- function(z) {
  kerval <- exp(-(z^2) / 2) / sqrt(2 * pi)
  return(kerval)
}


# Chol --------------------------------------------------------------------

# Cholesky decomposition (for a complex matrix)

#' Complex Cholesky decomposition
#'
#' Computes the lower-triangular Cholesky factor \eqn{G} of a Hermitian
#' positive-definite (possibly complex) matrix, so that \eqn{f = G G^*}.
#'
#' @param f A \eqn{P \times P}{P x P} Hermitian positive-definite matrix with
#'   \eqn{P \ge 2}.
#'
#' @return The \eqn{P \times P}{P x P} lower-triangular factor \eqn{G}, with a
#'   real positive diagonal.
#' @seealso \code{\link{spec2inChol}}
#' @examples
#' f <- matrix(c(2, 1 - 1i, 1 + 1i, 3), 2, 2)
#' G <- Chol(f)
#' all.equal(G %*% Conj(t(G)), f)
#' @export Chol
Chol <- function(f) {
  Gw <- sqrt(f[1, 1])
  
  for (k in 2:nrow(f)) {
    fk <- f[1:k - 1, k]
    fkk <- f[k, k]
    gk <- solve(Gw, fk)
    gkk <- sqrt(fkk - ((Conj(t(gk))) %*% gk))
    Gw <- rbind(cbind(Gw, matrix(0, nrow = nrow(as.matrix(Gw)), ncol = 1)), cbind(Conj(t(gk)), gkk))
  }
  
  return(Gw)
}


# spec2inChol -------------------------------------------------------------

# Takes spectral components and gives the inverse Cholesky components

#' Convert between spectral matrices and inverse Cholesky components
#'
#' \code{spec2inChol} maps vectorized spectral matrices \eqn{f(\omega_k)} to
#' the vectorized components of \eqn{G(\omega_k)}, the Cholesky factor of
#' \eqn{f^{-1}(\omega_k) = G G^*}. \code{inChol2spec} is the inverse map, and
#' is used to rebuild spectral estimates from fitted Cholesky components (see
#' \code{\link{tA2G}}). Any real input to \code{inChol2spec} with a nonzero
#' diagonal gives a Hermitian positive-definite spectral matrix.
#'
#' @param spec \eqn{P^2 \times K}{P^2 x K} matrix of vectorized spectral
#'   matrices (see \code{\link{Cholindexer}} for the ordering).
#' @param inChol \eqn{P^2 \times K}{P^2 x K} matrix of vectorized Cholesky
#'   components.
#' @param init List returned by \code{\link{initial}} (only \code{P},
#'   \code{K}, \code{Nreal} and \code{indexer} are used).
#'
#' @return A \eqn{P^2 \times K}{P^2 x K} matrix: the vectorized inverse
#'   Cholesky components (\code{spec2inChol}) or spectral matrices
#'   (\code{inChol2spec}).
#' @examples
#' init <- initial(dat$X1, nbasis = 10, R = 2)
#' G <- spec2inChol(init$NW[[1]], init)
#' all.equal(inChol2spec(G, init), init$NW[[1]])
#' @export spec2inChol
spec2inChol <- function(spec, init) {
  
  P2 <- init$P^2
  inChol <- matrix(0, nrow = P2, ncol = init$K)
  
  for (k in 1:init$K) {
    # spectral matrix at freq k
    F_k <- matrix(0, nrow = init$P, ncol = init$P)
    
    for (p in 1:init$Nreal) {
      F_k[init$indexer[p, 3], init$indexer[p, 4]] <- spec[p, k]
      F_k[init$indexer[p, 4], init$indexer[p, 3]] <- spec[p, k]
    }
    
    for (p in (init$Nreal + 1):P2) {
      F_k[init$indexer[p, 3], init$indexer[p, 4]] <- F_k[init$indexer[p, 3], init$indexer[p, 4]] + 1i * spec[p, k]
      F_k[init$indexer[p, 4], init$indexer[p, 3]] <- F_k[init$indexer[p, 4], init$indexer[p, 3]] - 1i * spec[p, k]
    }
    
    # inv cholesky at freq k
    cholMat <- (Chol(solve(F_k)))
    for (p in 1:init$Nreal) inChol[p, k] <- Re(cholMat[init$indexer[p, 3], init$indexer[p, 4]])
    for (p in (init$Nreal + 1):P2) inChol[p, k] <- Im(cholMat[init$indexer[p, 3], init$indexer[p, 4]])
  }
  
  return(inChol)
}


# inChol2spec -------------------------------------------------------------

# Takes the inverse Cholesky components and returns the components of the 
# spectral matrix

#' @rdname spec2inChol
#' @export inChol2spec
inChol2spec <- function(inChol, init) {
  P2 <- init$P^2
  spec <- matrix(0, nrow = P2, ncol = init$K)

  for (k in 1:init$K) {

    # inv cholesky at freq k
    cholMat <- matrix(0, nrow = init$P, ncol = init$P)

    for (p in 1:init$Nreal) {
      cholMat[init$indexer[p, 3], init$indexer[p, 4]] <- inChol[p, k]
    }

    for (p in (init$Nreal + 1):P2) {
      cholMat[init$indexer[p, 3], init$indexer[p, 4]] <- cholMat[init$indexer[p, 3], init$indexer[p, 4]] + 1i * inChol[p, k]
    }

    # spec matrix at freq k
    F_k <- solve(cholMat %*% Conj(t(cholMat)))

    for (p in 1:init$Nreal) {
      spec[p, k] <- Re(F_k[init$indexer[p, 3], init$indexer[p, 4]])
    }

    for (p in (init$Nreal + 1):P2) {
      spec[p, k] <- Im(F_k[init$indexer[p, 3], init$indexer[p, 4]])
    }
  }

  return(spec)
}


# psd ---------------------------------------------------------------------

# True power spectral density(PSD) for VARMA(p,q)

#' True spectral density matrix of a VARMA process
#'
#' Evaluates the spectral density matrix
#' \eqn{f(\omega) = \Phi(\omega)^{-1} \Theta(\omega) \Sigma \Theta(\omega)^* \Phi(\omega)^{-*}}{f(w) = Phi(w)^-1 Theta(w) Sigma Theta(w)^* Phi(w)^-*}
#' of the VARMA process \eqn{\Phi(L) X_t = \Theta(L) \epsilon_t}{Phi(L) X_t = Theta(L) e_t}
#' at the frequencies in \code{init}, where
#' \eqn{\Phi(\omega) = \sum_{j} \Phi_j e^{-2\pi i \omega j}}{Phi(w) = sum_j Phi_j exp(-2 pi i w j)}.
#'
#' @param paras A list of three elements: a list of the AR lag-polynomial
#'   matrices \eqn{I, -\Phi_1, \dots, -\Phi_d}, a list of the MA lag-polynomial
#'   matrices \eqn{I, \Theta_1, \dots, \Theta_q}, and the innovation covariance
#'   matrix \eqn{\Sigma}.
#' @param init List returned by \code{\link{initial}} (only \code{P},
#'   \code{K}, \code{omegas}, \code{Nreal} and \code{indexer} are used).
#'
#' @return A \eqn{P^2 \times K}{P^2 x K} matrix of vectorized spectral
#'   matrices.
#' @examples
#' init <- initial(dat$X1, nbasis = 10, R = 2)
#' # bivariate VAR(1) with correlated innovations
#' Phi1 <- matrix(c(0, -0.6, 0.5, -0.2), 2, 2)
#' Sigma <- matrix(c(1, 0.9, 0.9, 1), 2, 2)
#' f <- psd(list(list(diag(2), -Phi1), list(diag(2)), Sigma), init)
#' matplot(init$omegas, t(f), type = "l", lty = 1,
#'         xlab = "Frequency", ylab = "Spectral components")
#' @export psd
psd <- function(paras, init) {
  Phi <- paras[[1]]
  Theta <- paras[[2]]
  sigma <- paras[[3]]
  P2 <- init$P^2

  # True spectral density function for VARMA(p,q)
  F <- function(w) {
    A <- B <- matrix(0, nrow = init$P, ncol = init$P)
    for (j in 1:(length(Phi))) {
      A <- A + Phi[[j]] * exp(-2 * pi * 1i * w * (j - 1))
    }
    for (j in 1:(length(Theta))) {
      B <- B + Theta[[j]] * exp(-2 * pi * 1i * w * (j - 1))
    }

    invA <- solve(A)
    f <- invA %*% B %*% sigma %*% Conj(t(B)) %*% Conj(t(invA))
    return(f)
  }

  ## spectral density at each frequency
  Flist <- lapply(init$omegas, F)
  F
  # vectorize P*P spectral matrix at all frequencies(Flist) and get a P2*K matrix(spec)
  spec <- matrix(0, nrow = P2, ncol = init$K)

  for (k in 1:init$K) {
    for (p in 1:init$Nreal) spec[p, k] <- Re(Flist[[k]][init$indexer[p, 3], init$indexer[p, 4]])
    for (p in (init$Nreal + 1):P2) spec[p, k] <- Im(Flist[[k]][init$indexer[p, 3], init$indexer[p, 4]])
  }

  return(spec)
}

# loglike -----------------------------------------------------------------

# Whittle negative loglikelihood

#' Penalized Whittle negative log-likelihood
#'
#' Evaluates the Whittle negative log-likelihood
#' \deqn{\sum_{i=1}^M \sum_{k=1}^K \left\{-\log|G_i G_i^*| + Y_{i,k}^* G_i G_i^* Y_{i,k}\right\}}{sum_i sum_k { -log|G_i G_i^*| + Y_ik^* G_i G_i^* Y_ik }}
#' of the model \eqn{\Gamma_j = B \Theta_j A_j^\top}{Gamma_j = B Theta_j t(A_j)}
#' and adds the roughness penalty
#' \eqn{\lambda\,\mathrm{tr}(\Theta_j^\top D \Theta_j)}{lambda tr(t(Theta_j) D Theta_j)}.
#' This is the reference R implementation; \code{\link{loglike_c_cube}} is a
#' faster C++ version for the different bases model.
#'
#' @param thetas Basis coefficients: an \eqn{L \times R \times P^2}{L x R x P^2}
#'   array (different bases) or an \eqn{L \times R}{L x R} matrix (common
#'   bases).
#' @param As \eqn{M \times R \times P^2}{M x R x P^2} array of scores.
#' @param lambda Smoothing parameter. For different bases: a scalar, a vector
#'   of length \eqn{P^2} or a \eqn{P^2 \times R}{P^2 x R} matrix. For common
#'   bases: a scalar or a vector of length \eqn{R}.
#' @param init List returned by \code{\link{initial}}.
#' @param tildeX List of \eqn{M} complex \eqn{P \times K}{P x K} matrices of
#'   DFTs, usually \code{init$tildeX}.
#'
#' @return A list with components \code{ploglike} (penalized) and
#'   \code{loglike} (unpenalized) negative log-likelihood.
#' @examples
#' init <- initial(dat$X1, nbasis = 10, R = 2)
#' loglike(init$thetas, init$As, lambda = 1, init, init$tildeX)
#' @export loglike
loglike <- function(thetas, As, lambda, init, tildeX) {
  hatG_m <- tA2G(thetas, As, init)
  l <- matrix(0, nrow = init$K, ncol = nrow(As))
  for (i in 1:nrow(As)) {
    for (k in 1:init$K) {
      # get the inverse cholesky matrix
      Cholmat_i <- matrix(0, nrow = init$P, ncol = init$P)
      for (p in 1:init$Nreal) Cholmat_i[init$indexer[p, 3], init$indexer[p, 4]] <- hatG_m[p, k, i]
      for (p in (init$Nreal + 1):init$P^2) Cholmat_i[init$indexer[p, 3], init$indexer[p, 4]] <- Cholmat_i[init$indexer[p, 3], init$indexer[p, 4]] + 1i * hatG_m[p, k, i]
      ## print(Cholmat_i)

      # get the loglikelihood
      l[k, i] <- -log(prod(eigen(Cholmat_i %*% Conj(t(Cholmat_i)), only.values = TRUE)$values)) + Conj(t(tildeX[[i]][, k])) %*% Cholmat_i %*% Conj(t(Cholmat_i)) %*% tildeX[[i]][, k]
    }
  }
  loglike <- Re(sum(l))

  # Penalized loglikelihood
  if (!is.matrix(thetas)) {
    lam.pen <- rep(0, 0)
    for (j in 1:init$P^2) {
      if (is.vector(lambda) && length(lambda) == 1) {
        lam.pen[j] <- sum(lambda * diag(t(thetas[, , j]) %*% init$D %*% thetas[, , j])) # lambda: scalar
      } else if (is.vector(lambda)) {
        lam.pen[j] <- sum(lambda[j] * diag(t(thetas[, , j]) %*% init$D %*% thetas[, , j])) # lambda_P2: vector
      } else if (is.matrix(lambda)) {
        lam.pen[j] <- sum(diag(lambda[j, ], ncol(thetas)) %*% diag(t(thetas[, , j]) %*% init$D %*% thetas[, , j])) # lambda_P2*R: matrix
      }
    }
    ploglike <- Re(loglike + sum(lam.pen))
  } else {
    if (is.vector(lambda) && length(lambda) == 1) lambda <- rep(lambda, ncol(thetas))
    ploglike <- Re(loglike + sum(diag(lambda, ncol(thetas)) %*% diag(t(thetas) %*% init$D %*% thetas)))
  }

  return(list(ploglike = ploglike, loglike = loglike))
}

# Cholindexer -------------------------------------------------------------

# column of indexer:
# [,1]: index of the vector (1,...,P^2); [,2]: index of real or not {0,1}; [,3]: row number(1,...,P)
# [,4]: column number (1,...,(row number)); [,5]: index of diagonal elements {1,0}

#' Index of the vectorized Cholesky and spectral components
#'
#' Describes how a \eqn{P \times P}{P x P} Hermitian matrix, or a
#' lower-triangular complex matrix, is stored as a real vector of length
#' \eqn{P^2}: the real parts of the \eqn{P(P+1)/2} lower-triangular entries
#' come first (column by column), followed by the imaginary parts of the
#' \eqn{P(P-1)/2} strictly lower-triangular entries.
#'
#' @param P Dimension of the time series (\eqn{P \ge 2}).
#'
#' @return A \eqn{P^2 \times 5}{P^2 x 5} integer-valued matrix whose columns
#'   give, for each element of the vector: its position; \code{1} for a real
#'   part and \code{0} for an imaginary part; the row of the matrix entry; the
#'   column of the matrix entry; and \code{1} if the entry is on the diagonal.
#' @examples
#' Cholindexer(2)
#' Cholindexer(3)
#' @export Cholindexer
Cholindexer <- function(P) {
  Nreal <- P * (P + 1) / 2
  indexer <- matrix(0, nrow = P^2, ncol = 5)
  indexer[, 1] <- 1:P^2
  indexer[1:Nreal, 2] <- 1
  bigguy <- matrix(c(kronecker(rep(1, P), 1:P), kronecker(1:P, rep(1, P))), nrow = P^2)
  indexer[1:Nreal, 3:4] <- bigguy[bigguy[, 1] >= bigguy[, 2], ]
  indexer[(Nreal + 1):P^2, 3:4] <- bigguy[bigguy[, 1] > bigguy[, 2], ]
  indexer[indexer[, 3] == indexer[, 4], 5] <- 1

  return(indexer)
}


# drivG -------------------------------------------------------------------

# Required derivatives for Fisher-scoring algorithm
#' Derivatives of the Cholesky factor
#'
#' Computes the derivatives needed by the Newton-Raphson (Fisher scoring)
#' steps: the derivative \eqn{\mathcal{D}G_j}{DG_j} of the Cholesky factor
#' \eqn{G} with respect to each of its \eqn{P^2} vectorized components, and
#' the second derivatives of \eqn{G G^*}. Called by \code{\link{initial}}.
#'
#' @param init A list with at least the components \code{P} and
#'   \code{indexer} (see \code{\link{Cholindexer}}).
#'
#' @return A list with
#'   \item{dG}{list of \eqn{P^2} matrices \eqn{\partial G / \partial g_j}{dG/dg_j}.}
#'   \item{HGG}{\eqn{P \times P \times P^2 \times P^2}{P x P x P^2 x P^2} array
#'     of \eqn{\partial^2 (G G^*) / \partial g_{j_1} \partial g_{j_2}}{d^2(GG^*)/dg_j1 dg_j2}.}
#'   \item{HGGvec_c}{the same second derivatives with each \eqn{P \times P}{P x P}
#'     matrix vectorized.}
#' @examples
#' d <- drivG(list(P = 2, indexer = Cholindexer(2)))
#' d$dG
#' @export drivG
drivG <- function(init) {
  
  # Derivative of Cholesky components(G) w.r.t. its vectorized elements(vecG)
  dG <- list()
  for (p in 1:init$P^2) {
    dGp <- matrix(0, nrow = init$P, ncol = init$P)
    dGp[init$indexer[p, 3], init$indexer[p, 4]] <- init$indexer[p, 2] + (1 - init$indexer[p, 2]) * 1i
    dG[[p]] <- dGp
  }
  
  # 2nd derivative of (GG*) w.r.t. the vectorized elements(vecG)
  HGG <- array(0, dim = c(init$P, init$P, init$P^2, init$P^2))
  for (j1 in 1:init$P^2) {
    for (j2 in 1:init$P^2) {
      HGG[, , j1, j2] <- dG[[j1]] %*% Conj(t(dG[[j2]])) + dG[[j2]] %*% Conj(t(dG[[j1]]))
    }
  }
  HGGvec_c <- apply(HGG, c(3, 4), as.vector)
  
  return(list(dG = dG, HGG = HGG, HGGvec_c = HGGvec_c))
}

# G2tA --------------------------------------------------------------------

# Projects G onto linear subspace (of the bases)
G2tA <- function(G_m, theta_P2, init) {
  As <- array(0, dim = c(init$m, init$R, init$P^2))
  smooth <- solve(t(init$B) %*% init$B) %*% t(init$B)
  theta.tA <- array(0, dim = c(ncol(init$B), init$m, init$P^2))
  
  for (j in 1:init$P^2) theta.tA[, , j] <- smooth %*% G_m[j, , ]
  if (isTRUE(theta_P2)) {
    svd.theta.tA <- list()
    thetas <- array(0, dim = c(ncol(init$B), init$R, init$P^2))
    
    for (j in 1:init$P^2) {
      svd.theta.tA[[j]] <- svd(theta.tA[, , j])
      thetas[, , j] <- as.matrix(svd.theta.tA[[j]]$u[, 1:init$R])
      As[, , j] <- as.matrix(svd.theta.tA[[j]]$v %*% diag(svd.theta.tA[[j]]$d, ncol(svd.theta.tA[[j]]$v)))[, 1:init$R]
    }
    
  } else {
    mergG <- matrix(0, init$K, init$m * init$P^2)
    
    for (j in 1:init$P^2) mergG[, (((j - 1) * (init$m)) + 1):(j * init$m)] <- G_m[j, , ]
   
    theta.tmergA <- smooth %*% mergG
    svd.theta.tmergA <- svd(theta.tmergA)
    thetas <- as.matrix(svd.theta.tmergA$u[, 1:init$R])
    mergA <- as.matrix(svd.theta.tmergA$v %*% diag(svd.theta.tmergA$d, ncol(svd.theta.tmergA$v)))
    mergA <- as.matrix(mergA[, 1:init$R])
    
    for (j in 1:init$P^2) As[, , j] <- mergA[((j - 1) * init$m + 1):(j * init$m), ]
  }
  
  return(list(thetas = thetas, As = As, theta.tA = theta.tA))
}


# tA2G --------------------------------------------------------------------

# Takes thetas and As and reconstructs G
#' Reconstruct the Cholesky components from the model parameters
#'
#' Evaluates \eqn{\Gamma_j = B \Theta_j A_j^\top}{Gamma_j = B Theta_j t(A_j)},
#' \eqn{j = 1, \dots, P^2}, the vectorized inverse Cholesky components of all
#' series at all frequencies. Apply \code{\link{inChol2spec}} to each slice to
#' obtain the spectral density estimates.
#'
#' @param thetas Basis coefficients: an \eqn{L \times R \times P^2}{L x R x P^2}
#'   array (different bases) or an \eqn{L \times R}{L x R} matrix (common
#'   bases).
#' @param As \eqn{M \times R \times P^2}{M x R x P^2} array of scores.
#' @param init List returned by \code{\link{initial}}.
#'
#' @return A \eqn{P^2 \times K \times M}{P^2 x K x M} array; slice \code{[, , i]}
#'   holds the vectorized Cholesky components of series \eqn{i}.
#' @seealso \code{\link{iterate}}
#' @examples
#' init <- initial(dat$X1, nbasis = 10, R = 2)
#' G <- tA2G(init$thetas, init$As, init)
#' # spectral estimates of all series: P^2 x K x M
#' spec <- array(apply(G, 3, inChol2spec, init = init),
#'               dim = c(init$P^2, init$K, init$m))
#' dim(spec)
#' @export tA2G
tA2G <- function(thetas, As, init) {
  hatG_m <- array(0, dim = c(init$P^2, init$K, nrow(As)))
  theta.tA <- array(0, dim = c(ncol(init$B), nrow(As), init$P^2))
  
  for (j in 1:init$P^2) {
    
    if (!is.matrix(thetas)) {
      theta.tA[, , j] <- thetas[, , j] %*% t(As[, , j])
    } else {
      theta.tA[, , j] <- thetas %*% t(As[, , j])
    }
    
    hatG_m[j, , ] <- init$B %*% theta.tA[, , j]
  }
  
  return(hatG_m)
}


# init.theta.A_i ----------------------------------------------------------

init.theta.A_i <- function(init, i = 1) {
  mtilde <- length(i)
  R <- min(mtilde, ncol(init$thetas))
  As <- array(0, dim = c(mtilde, R, init$P^2))
  
  if (!is.matrix(init$thetas)) {
    thetas <- array(0, dim = c(ncol(init$B), R, init$P^2))
     
    for (j in 1:init$P^2) {
      #for (j in 1:init$P^2) theta.tA[, , j] <- smooth %*% G_m[j, , ]
      svd.theta.tA <- svd(init$theta.tA[, i, j])
      thetas[, , j] <- as.matrix(svd.theta.tA$u[, 1:R])
      As[, , j] <- as.matrix(svd.theta.tA$v %*% diag(svd.theta.tA$d, ncol(svd.theta.tA$v)))[, 1:R]
    }
    
  } else {
    merg.theta.tA <- matrix(0, nrow(init$thetas), mtilde * init$P^2)
    
    for (j in 1:init$P^2) merg.theta.tA[, (((j - 1) * (mtilde)) + 1):(j * mtilde)] <- init$theta.tA[, i, j]
    
    svd.merg.theta.tA <- svd(merg.theta.tA)
    thetas <- as.matrix(svd.merg.theta.tA$u[, 1:R])
    mergA <- as.matrix(svd.merg.theta.tA$v %*% diag(svd.merg.theta.tA$d, ncol(svd.merg.theta.tA$v))[, 1:R])
    
    for (j in 1:init$P^2) As[, , j] <- mergA[((j - 1) * mtilde + 1):(j * mtilde), ]
  }
  
  return(list(thetas = thetas, As = As))
}


#----------------------------------------------------------------------------------------------
##Workingfish
Workingfish <- function(thetas, As, lambda, init, tildeX) {
  nbasis <- ncol(init$B)
  norder <- nbasis - length(init$knots) + 2
  P2 <- init$P^2
  R <- ncol(init$thetas)
  if (!is.matrix(thetas)) {
    U.theta <- matrix(0, nrow = nbasis * P2, ncol = R)
    H.theta <- array(0, dim = c(nbasis * P2, nbasis * P2, R))
    F.theta <- array(0, dim = c(nbasis * P2, nbasis * P2, R))
    tmp1i <- F.theta
    theta.chng <- matrix(0, nrow = nbasis * P2, ncol = R)
    curr.thetas <- array(0, dim = c(nbasis, R, P2))
    theta_L.P2_R <- apply(thetas, 2, as.vector) # L*R*P2 array (thetas) into L.P2*R matrix
  } else {
    U.theta <- matrix(0, nrow = nbasis, ncol = R)
    H.theta <- array(0, dim = c(nbasis, nbasis, R))
    F.theta <- array(0, dim = c(nbasis, nbasis, R))
    tmp1i <- F.theta
    theta.chng <- matrix(0, nrow = nbasis, ncol = R)
    curr.thetas <- matrix(0, nrow = nbasis, ncol = R)
  }

  A.chng <- matrix(0, nrow = R * P2, ncol = nrow(As))
  curr.As <- array(0, dim = c(nrow(As), R, P2))
  A_R.P2_m <- apply(As, 1, as.vector) # m*R*P2 array (As) into R.P2*m matrix
  hatG_m <- tA2G(thetas, As, init)
  hatG_m_list <- lapply(seq(dim(hatG_m)[3]), function(x) hatG_m[, , x])
  spec <- sapply(hatG_m_list, inChol2spec, init, simplify = "array")

  for (i in 1:nrow(As)) {
    # get the score,Hessian and negative Fisher matrices for A
    U.A_i <- matrix(0, nrow = R * P2, ncol = 1)
    FA_k_i <- array(0, dim = c(R, P2, R, P2, init$K))
    HA_k_i <- array(0, dim = c(R, P2, R, P2, init$K))
    G_i <- matrix(0, nrow = init$P, ncol = init$P)
    F_i <- matrix(0, nrow = init$P, ncol = init$P)
    temp1_i <- matrix(0, nrow = P2, ncol = init$K) # components of the first derivative of loglikelihood w.r.t Cholesky components(vecG)
    temp4_i <- array(0, dim = c(nbasis, P2, init$K, R)) # components of the first derivative of loglikelihood w.r.t thetas (chain rule)
    temp6_i <- array(0, dim = c(nbasis, P2, nbasis, P2, init$K, R)) # components of the second derivative of loglikelihood w.r.t thetas (chain rule)
    temp9_i <- array(0, dim = c(nbasis, P2, nbasis, P2, init$K, R)) # components of negative Fisher matrix for thetas (expectation of Hessian)
    for (k in 1:init$K) {
      tildeX_k_i <- tildeX[[i]][, k]
      CjtildeX_k_i <- Conj(t(tildeX_k_i))

      # get the G and F matrices for ith time series
      for (p in 1:init$Nreal) {
        G_i[init$indexer[p, 3], init$indexer[p, 4]] <- hatG_m[p, k, i]
        F_i[init$indexer[p, 3], init$indexer[p, 4]] <- spec[p, k, i]
        F_i[init$indexer[p, 4], init$indexer[p, 3]] <- spec[p, k, i]
      }
      for (p in (init$Nreal + 1):P2) {
        G_i[init$indexer[p, 3], init$indexer[p, 4]] <- G_i[init$indexer[p, 3], init$indexer[p, 4]] + 1i * hatG_m[p, k, i]
        F_i[init$indexer[p, 3], init$indexer[p, 4]] <- F_i[init$indexer[p, 3], init$indexer[p, 4]] + 1i * spec[p, k, i]
        F_i[init$indexer[p, 4], init$indexer[p, 3]] <- F_i[init$indexer[p, 4], init$indexer[p, 3]] - 1i * spec[p, k, i]
      }

      UA_k_i <- matrix(0, nrow = R, ncol = P2)
      temp2_k_i <- matrix(0, nrow = P2, ncol = P2)
      temp3_k_i <- temp2_k_i
      temp8_k_i <- temp2_k_i

      b <- max(which(init$knots <= init$omegas[k])) # select nonzero basis at freq k
      if (!is.matrix(thetas)) {
        B.theta_k <- matrix(0, nrow = R, ncol = P2)
        for (j in 1:P2) for (r in 1:R) B.theta_k[r, j] <- sum(init$B[k, b:(norder + b - 1)] * thetas[b:(norder + b - 1), r, j])
      } else {
        B.theta_k <- rep(0, 0)
        for (r in 1:R) B.theta_k[r] <- sum(init$B[k, b:(norder + b - 1)] * thetas[b:(norder + b - 1), r])
      }


      # calculate temp-s
      for (j1 in 1:P2) {
        temp1_i[j1, k] <- Re(CjtildeX_k_i %*% (init$drivG$dG[[j1]] %*% Conj(t(G_i)) + G_i %*% Conj(t(init$drivG$dG[[j1]]))) %*% tildeX_k_i + init$indexer[j1, 5] * ((-2) / hatG_m[j1, k, i]))
        diag(temp2_k_i)[j1] <- init$indexer[j1, 5] * (2 / (hatG_m[j1, k, i]^2))
        for (l in which(init$B[k, ] != 0)) {
          for (r in 1:R) {
            temp4_i[l, j1, k, r] <- temp1_i[j1, k] * As[i, r, j1] * init$B[k, l]
          }
        }
        for (j2 in 1:P2) {
          temp3_k_i[j1, j2] <- Re(CjtildeX_k_i %*% init$drivG$HGG[, , j1, j2] %*% tildeX_k_i) + temp2_k_i[j1, j2] # imaginary part is zero
          temp8_k_i[j1, j2] <- Re(sum(diag(((init$drivG$HGG[, , j1, j2]) %*% F_i))) + temp2_k_i[j1, j2]) # imaginary part is zero
          if (!is.matrix(thetas)) {
            for (r1 in 1:R) {
              UA_k_i[r1, j1] <- (temp1_i[j1, k]) * B.theta_k[r1, j1]
              for (r2 in 1:R) {
                HA_k_i[r1, j1, r2, j2, k] <- temp3_k_i[j1, j2] * B.theta_k[r1, j1] * B.theta_k[r2, j2]
                FA_k_i[r1, j1, r2, j2, k] <- temp8_k_i[j1, j2] * B.theta_k[r1, j1] * B.theta_k[r2, j2]
              }
            }
          } else {
            for (r1 in 1:R) {
              UA_k_i[r1, j1] <- (temp1_i[j1, k]) * B.theta_k[r1]
              for (r2 in 1:R) {
                HA_k_i[r1, j1, r2, j2, k] <- temp3_k_i[j1, j2] * B.theta_k[r1] * B.theta_k[r2]
                FA_k_i[r1, j1, r2, j2, k] <- temp8_k_i[j1, j2] * B.theta_k[r1] * B.theta_k[r2]
              }
            }
          }

          for (l1 in which(init$B[k, ] != 0)) {
            for (l2 in which(init$B[k, ] != 0)) {
              for (r in 1:R) {
                temp6_i[l1, j1, l2, j2, k, r] <- temp3_k_i[j1, j2] * init$B[k, l1] * init$B[k, l2] * As[i, r, j1] * As[i, r, j2]
                temp9_i[l1, j1, l2, j2, k, r] <- temp8_k_i[j1, j2] * init$B[k, l1] * init$B[k, l2] * As[i, r, j1] * As[i, r, j2]
              }
            }
          }
        }
      }

      U.A_i <- U.A_i + as.vector(UA_k_i)
    } # k-loop is closed

    H.A_i <- apply(apply(apply(HA_k_i, c(1, 2, 3, 4), sum), c(1, 2), as.vector), 1, as.vector)
    F.A_i <- apply(apply(apply(FA_k_i, c(1, 2, 3, 4), sum), c(1, 2), as.vector), 1, as.vector) # R*P2*R*P2*K array into R.P2*R.P2 matrix

    if (!is.matrix(thetas)) {
      temp5_i <- apply(apply(temp4_i, c(1, 2, 4), sum), 3, as.vector) # L*P2*K*R array into L.P2*R matrix (scores for thetas)
      temp10_i <- apply(apply(apply(temp9_i, c(1, 2, 3, 4, 6), sum), c(1, 2, 5), as.vector), c(1, 4), as.vector) # L*P2*L*P2*K*R array into L.P2*L.P2*R array
      temp7_i <- apply(apply(apply(temp6_i, c(1, 2, 3, 4, 6), sum), c(1, 2, 5), as.vector), c(1, 4), as.vector)
    } else {
      temp5_i <- apply(temp4_i, c(1, 4), sum) # L*P2*K*R array into L*R matrix (scores for thetas)
      temp10_i <- apply(temp9_i, c(1, 3, 6), sum) # L*P2*L*P2*K*R array into L*L*R array
      temp7_i <- apply(temp6_i, c(1, 3, 6), sum)
    }
    U.theta <- U.theta + temp5_i
    F.theta <- F.theta + temp10_i
    H.theta <- H.theta + temp7_i
    A.chng[, i] <- solve(F.A_i) %*% U.A_i

    if (i == nrow(As)) {
      if (!is.matrix(thetas)) {
        DF <- matrix(0, nrow = P2, ncol = R)
        if (is.vector(lambda)) lambda <- matrix(lambda, nrow = P2, ncol = R)
        for (r in 1:R) {
          tmp1i[, , r] <- solve(F.theta[, , r] + kronecker(diag(lambda[, r], P2), init$D))
          theta.chng[, r] <- tmp1i[, , r] %*% (U.theta[, r] + kronecker(diag(lambda[, r], P2), init$D) %*% theta_L.P2_R[, r])
          tmp1i.F <- tmp1i[, , r] %*% F.theta[, , r]
          for (j in 1:P2) DF[j, r] <- sum(diag(tmp1i.F[((j - 1) * nbasis + 1):(j * nbasis), ((j - 1) * nbasis + 1):(j * nbasis)]))
        }
      } else {
        DF <- rep(0, 0)
        if (length(lambda) == 1) lambda <- rep(lambda, R)
        for (r in 1:R) {
          tmp1i[, , r] <- solve(F.theta[, , r] + lambda[r] * init$D)
          theta.chng[, r] <- tmp1i[, , r] %*% (U.theta[, r] + lambda[r] * init$D %*% thetas[, r])
          DF[r] <- sum(diag(tmp1i[, , r] %*% F.theta[, , r]))
        }
      }
    }
  } # i-loop is closed

  # Step-halving algorithm
  repit <- TRUE
  counter <- 0
  while (repit && counter < 30) { # cap the number of step halvings
    tau <- (0.5)^counter
    for (i in 1:nrow(As)) curr.As[i, , ] <- matrix(A_R.P2_m[, i] - tau * A.chng[, i], nrow = R, ncol = P2)

    if (!is.matrix(thetas)) {
      for (r in 1:R) curr.thetas[, r, ] <- matrix(theta_L.P2_R[, r] - tau * theta.chng[, r], nrow = nbasis, ncol = P2)
    } else {
      for (r in 1:R) curr.thetas[, r] <- thetas[, r] - tau * theta.chng[, r]
    }
    if (loglike_c_cube(curr.thetas, curr.As, lambda, init, tildeX)$ploglike <= loglike_c_cube(thetas, As, lambda, init, tildeX)$ploglike) {
      repit <- FALSE
    } else {
      counter <- counter + 1
    }
  }
  return(list(thetas = curr.thetas, As = curr.As, DF = Re(DF), stepH = counter))
}

# iterate -----------------------------------------------------------------

# Iterartive algorithm
#' Fit the NMCSDE model by penalized Whittle likelihood
#'
#' Minimizes the penalized Whittle negative log-likelihood with Newton-Raphson
#' (Fisher scoring) iterations and step halving. Each iteration updates the
#' scores \eqn{A} and the basis coefficients \eqn{\Theta}, and optionally the
#' smoothing parameter.
#'
#' When \code{update.lam = TRUE}, the smoothing parameter is updated within
#' the iterations by
#' \deqn{\lambda_{new}^{-1} = \frac{\mathrm{tr}(\hat\Theta_{old}^\top D \hat\Theta_{old})}{\mathrm{df}(\hat\lambda_{old}) - (a-1)},}{1/lambda_new = tr(t(Theta_old) D Theta_old) / (df(lambda_old) - (a - 1)),}
#' where \eqn{\mathrm{df}(\lambda)}{df(lambda)} are the effective degrees of
#' freedom (Schall, 1991). This avoids a grid search. With
#' \code{update.lam = FALSE} the value of \code{lambda} is kept fixed and can
#' be chosen by minimizing the returned \code{AIC} over a grid.
#'
#' Iterations stop when the decrease of the penalized negative
#' log-likelihood is at most \eqn{0.2\sqrt{K}}{0.2 sqrt(K)} or after
#' \code{n.iter} iterations.
#'
#' @param init List returned by \code{\link{initial}}.
#' @param lambda Starting (or fixed) smoothing parameter. For the different
#'   bases model: a scalar or a vector of length \eqn{P^2} (one value per
#'   Cholesky component). For the common bases model: a scalar or a vector of
#'   length \eqn{R}.
#' @param a Order of the penalty, used in the update of \eqn{\lambda}. Defaults
#'   to the value used in \code{\link{initial}}.
#' @param n.iter Maximum number of iterations.
#' @param which.spec Optional integer vector of series indices. If given, only
#'   these series are fitted (for example \code{which.spec = i} gives the
#'   noncollective NMSDE estimate of series \eqn{i}); starting values are
#'   recomputed from \code{init} for this subset.
#' @param update.lam Logical; update \eqn{\lambda} within the iterations?
#' @param tildeX Optional list of DFT matrices to use instead of
#'   \code{init$tildeX} (used by \code{\link{NMSDEBoot}}).
#' @param verbose Logical; print the progress of each iteration?
#'
#' @return A list with
#'   \item{thetas, As}{lists with the parameter values at each iteration;
#'     the last element is the final estimate.}
#'   \item{ploglikes}{penalized negative log-likelihood at each iteration.}
#'   \item{Lambda}{list of smoothing parameter values at each iteration.}
#'   \item{AIC}{penalized negative log-likelihood plus twice the effective
#'     degrees of freedom, at each iteration.}
#'   \item{r.iter}{number of parameter sets stored (iterations + 1).}
#'   \item{which.spec}{the \code{which.spec} argument.}
#'
#' @references Schall, R. (1991). Estimation in generalized linear models
#'   with random effects. \emph{Biometrika}, 78(4), 719--727.
#'   \doi{10.1093/biomet/78.4.719}
#' @seealso \code{\link{initial}}, \code{\link{tA2G}}, \code{\link{plot_sds}}
#' @examples
#' X <- dat$X1
#' init <- initial(X, nbasis = 10, R = 3)
#' fit <- iterate(init, lambda = 1, n.iter = 5)
#' fit$ploglikes
#'
#' # spectral estimates of all series (P^2 x K x M)
#' last <- fit$r.iter
#' G <- tA2G(fit$thetas[[last]], fit$As[[last]], init)
#' spec <- array(apply(G, 3, inChol2spec, init = init),
#'               dim = c(init$P^2, init$K, init$m))
#' @export iterate
iterate <- function(init, lambda = 0, a = if (is.null(init$a)) 2 else init$a, n.iter = 40,
                    which.spec = NULL, update.lam = TRUE, tildeX = NULL, verbose = TRUE) {

  eps <- 0.2*sqrt(init$K)
  #eps <- 2
  thetas.lst <- list()
  As.lst <- list()
  ploglikes <- rep(0, 0)
  AIC <- rep(0, 0)
  Lambda <- list()
  l <- 1
  if (!is.null(which.spec)) {
    init_i <- init.theta.A_i(init, i = which.spec)
    init$thetas <- init_i$thetas
    init$As <- init_i$As
    init$m <- length(which.spec)
    #init$freqdom <- init$freqdom[which.spec]
    init$perds <- init$perds[which.spec]
    init$tildeX <- init$tildeX[which.spec]
    if (is.null(tildeX)) {
      tildeX <- init$tildeX
    }else{
      tildeX <- tildeX[which.spec]
    }
  }
  if (is.null(tildeX)) tildeX <- init$tildeX
  
  thetas.lst[[l]] <- Re(init$thetas)
  As.lst[[l]] <- Re(init$As)
  
  Lambda[[l]] <- lambda

  if (is.matrix(init$thetas)) {
    ploglikes[l] <- loglike_c_mat(thetas.lst[[l]], As.lst[[l]], Lambda[[l]], init, tildeX)$ploglike
  } else{
    ploglikes[l] <- loglike_c_cube(thetas.lst[[l]], As.lst[[l]], Lambda[[l]], init, tildeX)$ploglike
  }

  AIC[l] <- ploglikes[l]
  diff <- eps + ncol(init$thetas)

  while (abs(diff) > eps && l < n.iter) {
    if (is.matrix(init$thetas)) {
      workvals <- Workingfish_mat(thetas.lst[[l]], As.lst[[l]], Lambda[[l]], init, tildeX)
    } else {
      workvals <- Workingfish_cube(thetas.lst[[l]], As.lst[[l]], Lambda[[l]], init, tildeX)
      #workvals <- Workingfish(thetas.lst[[l]], As.lst[[l]], Lambda[[l]], init, tildeX)
    }
   
    thetas.lst[[l + 1]] <- Re(workvals$thetas)
    As.lst[[l + 1]] <- Re(workvals$As)

    if (update.lam) {
      if (is.matrix(init$thetas) && length(lambda) == 1) {
        Lambda[[l + 1]] <- rep(1 / (sum(diag(t(thetas.lst[[(l)]]) %*% init$D %*% thetas.lst[[(l)]])) / (sum(workvals$DF) - (a - 1))), ncol(init$thetas)) # lambda(scalar)
      } else if (is.matrix(init$thetas)) {
        Lambda[[l + 1]] <- 1 / (diag(t(thetas.lst[[(l)]]) %*% init$D %*% thetas.lst[[(l)]]) / (as.vector(workvals$DF) - (a - 1))) # lambda_R(vector)
      } else {
        diagtemp <- matrix(0, nrow = init$P^2, ncol = ncol(init$thetas))
        for (r in 1:ncol(init$thetas)) {
          diagtemp[, r] <- diag(t(thetas.lst[[(l)]][, r, ]) %*% init$D %*% thetas.lst[[(l)]][, r, ])
        }
        if (length(lambda) == 1) {
          Lambda[[l + 1]] <- 1 / (sum(diagtemp) / (sum(workvals$DF) - (a - 1))) # lambda(scalar)
        } else if (is.vector(lambda)) {
          Lambda[[l + 1]] <- 1 / (apply(diagtemp, 1, sum) / (apply(workvals$DF, 1, sum) - (a - 1))) # lambda_P2(vector)
        } else {
          Lambda[[l + 1]] <- matrix(0, nrow = init$P^2, ncol = ncol(init$thetas))
          for (r in 1:ncol(init$thetas)) Lambda[[l + 1]][, r] <- 1 / (diagtemp[, r] / (workvals$DF[, r] - (a - 1))) # lambda_P2*R(matrix)
        }
      }
    } else {
      Lambda[[l + 1]] <- as.vector(Lambda[[l]])
    }
    
    if (is.matrix(init$thetas)) {
      ploglikes[l + 1] <- loglike_c_mat(thetas.lst[[l + 1]], As.lst[[l + 1]], Lambda[[l + 1]], init, tildeX)$ploglike
    } else{
      ploglikes[l + 1] <- loglike_c_cube(thetas.lst[[l + 1]], As.lst[[l + 1]], Lambda[[l + 1]], init, tildeX)$ploglike
    }
    
    AIC[(l + 1)] <- Re(ploglikes[l + 1] + 2 * sum(workvals$DF))
    diff <- (ploglikes[l] - ploglikes[l + 1])
    l <- l + 1; #
    
    if (verbose) cat(paste("\n l=", (l), " pLogL=", round(ploglikes[(l)]), " AIC=", round(AIC[(l)]), " lambda=", round(Lambda[[(l)]], 2), " DF=", round(workvals$DF), " stepH=", workvals$stepH, " diff=", round(diff, 2), sep = ""), "\n")
  }
  
  return(list(thetas = thetas.lst, As = As.lst, ploglikes = ploglikes, Lambda = Lambda, AIC = AIC, r.iter = l, which.spec = which.spec))
}


# NMSDEBoot ---------------------------------------------------------------
#' Spectral bootstrap for NMSDE and NMCSDE estimates
#'
#' Multivariate spectral bootstrap (Berkowitz and Diebold, 1998) for pointwise
#' confidence intervals of the spectral density matrices. For large \eqn{T},
#' \eqn{Y_k \approx \hat f(\omega_k)^{1/2} \epsilon_k}{Y_k ~ f(w_k)^(1/2) e_k}
#' with independent standard complex normal \eqn{\epsilon_k}{e_k}. Each
#' bootstrap sample draws new \eqn{\epsilon_k}{e_k}, builds bootstrap DFTs
#' from the Cholesky square root of \code{hatF}, and re-estimates the spectra
#' with \code{\link{iterate}}.
#'
#' @param hatF \eqn{P^2 \times K \times M}{P^2 x K x M} array with the
#'   vectorized spectral estimates of the \eqn{M} series (for example from an
#'   NMCSDE fit; see the examples).
#' @param nBoot Number of bootstrap samples.
#' @param init List returned by \code{\link{initial}}.
#' @param lambda Smoothing parameter passed to \code{\link{iterate}}. If
#'   \code{update_lam = FALSE}, a vector of candidate values; the one with the
#'   smallest final AIC is used for each bootstrap sample.
#' @param theta_P2 Currently unused: the model (common or different bases) is
#'   the one set up in \code{init}.
#' @param update_lam Logical; update \eqn{\lambda} within the iterations (see
#'   \code{\link{iterate}})?
#' @param method \code{"collective"} (default) re-estimates all series jointly
#'   (NMCSDE); \code{"separate"} re-estimates each series on its own (NMSDE).
#' @param verbose Logical; print the progress of the iterations?
#'
#' @return A list of \eqn{P^2} arrays, one per vectorized spectral component,
#'   each of dimension \eqn{nBoot \times K \times M}{nBoot x K x M}. Pointwise
#'   percentile intervals are obtained with \code{apply(x, c(2, 3), quantile,
#'   probs = c(0.025, 0.975))}.
#' @references Berkowitz, J. and Diebold, F. X. (1998). Bootstrapping
#'   multivariate spectra. \emph{Review of Economics and Statistics}, 80(4),
#'   664--666. \doi{10.1162/003465398557753}
#' @examples
#' \donttest{
#' init <- initial(dat$X1, nbasis = 10, R = 3)
#' fit <- iterate(init, lambda = 1, n.iter = 5, verbose = FALSE)
#' G <- tA2G(fit$thetas[[fit$r.iter]], fit$As[[fit$r.iter]], init)
#' hatF <- array(apply(G, 3, inChol2spec, init = init),
#'               dim = c(init$P^2, init$K, init$m))
#' boot <- NMSDEBoot(hatF, nBoot = 2, init = init, lambda = 1,
#'                   update_lam = TRUE)
#' # 95% pointwise interval for the first spectral component of series 1
#' ci <- apply(boot[[1]][, , 1, drop = FALSE], 2, quantile, probs = c(0.025, 0.975))
#' }
#' @export NMSDEBoot
NMSDEBoot <- function(hatF, nBoot, init, lambda, theta_P2, update_lam, method = "collective",
                      verbose = FALSE) {
  omegas <- init$omegas
  P2 <- nrow(hatF)

  bootcell <- list(P2)
  for (j in 1:P2) {
    bootcell[[j]] <- array(0, dim = c(nBoot, init$K, init$m))
  }

  Fsqrt <- list(init$m)
  for (i in 1:init$m) {
    # Compute sqrt matrices
    Fsqrt[[i]] <- list(init$K)
    for (k in 1:init$K) {
      ## SDM at freq k
      f <- matrix(0, nrow = init$P, ncol = init$P)
      for (p in 1:init$Nreal) {
        f[init$indexer[p, 3], init$indexer[p, 4]] <- hatF[p, k, i]
        f[init$indexer[p, 4], init$indexer[p, 3]] <- hatF[p, k, i]
      }
      for (p in (init$Nreal + 1):P2) {
        f[init$indexer[p, 3], init$indexer[p, 4]] <- f[init$indexer[p, 3], init$indexer[p, 4]] +
          1i * hatF[p, k, i]
        f[init$indexer[p, 4], init$indexer[p, 3]] <- f[init$indexer[p, 4], init$indexer[p, 3]] -
          1i * hatF[p, k, i]
      }
      ## Cholesky of SDM at freq k
      Fsqrt[[i]][[k]] <- Chol(f)
    }
  }

  #################################################
  # Run for the different bootstraps
  for (nb in 1:nBoot) {
    tildeXBoot <- list(init$m)

    for (i in 1:init$m) {
      # pull the DFT residuals
      tildeXBoot[[i]] <- matrix(0, nrow = init$P, ncol = init$K)
      rand1 <- rnorm(n = init$P * init$K, mean = 0, sd = sqrt(1 / 2))
      rand1_mat <- matrix(rand1, nrow = init$P, ncol = init$K)
      rand2 <- rnorm(n = init$P * init$K, mean = 0, sd = sqrt(1 / 2))
      rand2_mat <- matrix(rand2, nrow = init$P, ncol = init$K)

      # Reshape the vector into a PxK matrix
      resids <- rand1_mat + 1i * rand2_mat

      # get the fourier transform and periodograms
      for (k in 1:init$K) {
        tildeXBoot[[i]][, k] <- Fsqrt[[i]][[k]] %*% resids[, k]
      }
    }

    if (method == "separate") {
      res <- list()
      for (t in 1:init$m) {

        # finding the best lambda
        if (!update_lam) { # (min AIC)
          iter <- list()
          AIC <- rep(0, 0)
          for (l in 1:length(lambda)) {
            iter[[l]] <- iterate(init, lambda = lambda[l], which.spec = t, update.lam = update_lam, tildeX = tildeXBoot, verbose = verbose)
            AIC.l <- unlist(iter[[l]])[substr(names(unlist(iter[[l]])), 1, 3) == "AIC"]
            AIC[l] <- AIC.l[length(AIC.l)]
          }
          res[[t]] <- iter[[which.min(AIC)]]
        } else { # update_lambda
          iter <- iterate(init, lambda = lambda, which.spec = t, update.lam = update_lam, tildeX = tildeXBoot, verbose = verbose)
          res[[t]] <- iter
        }

        # get hatF and save
        hatFBoot <- array(apply(tA2G(res[[t]]$thetas[[length(res[[t]]$thetas)]], res[[t]]$As[[length(res[[t]]$thetas)]], init), 3, inChol2spec, init), dim = c(init$P^2, init$K, 1))
        for (j in 1:P2) {
          bootcell[[j]][nb, , t] <- hatFBoot[j, , 1]
        }
      }
    } else {
      if (!update_lam) {
        iter <- list()
        AIC <- rep(0, 0)
        for (i in 1:length(lambda)) {
          iter[[i]] <- iterate(init, lambda = lambda[i], update.lam = update_lam, tildeX = tildeXBoot, verbose = verbose)
          AIC.i <- unlist(iter[[i]])[substr(names(unlist(iter[[i]])), 1, 3) == "AIC"]
          AIC[i] <- AIC.i[length(AIC.i)]
        }
        res <- iter[[which.min(AIC)]]
      } else {
        iter <- iterate(init, lambda = lambda, update.lam = update_lam, tildeX = tildeXBoot, verbose = verbose)
        res <- iter
      }

      # get hatF and save
      hatFBoot <- array(apply(tA2G(res$thetas[[length(res$thetas)]], res$As[[length(res$thetas)]], init), 3, inChol2spec, init), dim = c(init$P^2, init$K, init$m))
      for (j in 1:P2) {
        bootcell[[j]][nb, , ] <- hatFBoot[j, , ]
      }
    }
  }
  return(bootcell)
}

# plot_sds ----------------------------------------------------------------

# Plots estimate of the components spectral matrix and Trues
#' Plot estimated spectral density matrices
#'
#' For each series, draws the \eqn{P^2} vectorized components of the final
#' spectral estimate from \code{\link{iterate}} on a \eqn{P \times P}{P x P}
#' grid of panels (real parts of the lower triangle, then imaginary parts;
#' see \code{\link{Cholindexer}}). If \code{init} contains true spectra
#' (argument \code{tru} of \code{\link{initial}}), they are added in blue.
#' One page is drawn per series.
#'
#' @param res List returned by \code{\link{iterate}}.
#' @param init List returned by \code{\link{initial}} and used to fit
#'   \code{res}.
#'
#' @return Called for its side effect (plots); returns \code{NULL}
#'   invisibly.
#' @examples
#' init <- initial(dat$X1, nbasis = 10, R = 3, tru = datPar$modelPar1)
#' fit <- iterate(init, lambda = 1, n.iter = 3, verbose = FALSE)
#' plot_sds(fit, init)
#' @export plot_sds
plot_sds <- function(res, init) {
  oldpar <- par(mfrow = c(init$P, init$P), oma = c(0.5, 0.5, 2.5, 0.5))
  on.exit(par(oldpar))
  if (!is.null(res$which.spec)) {
    init$m <- length(res$which.spec)
    if (!is.null(init$rnd)) init$rnd <- init$rnd[res$which.spec]
  }

  
  for (i in 1:init$m) {
    for (j in 1:init$P^2) {
      
      spec_j <- inChol2spec(tA2G(res$thetas[[length(res$thetas)]], res$As[[length(res$As)]], init)[, , i], init)[j, ]
      ylim <- range(spec_j, if (!is.null(init$rnd)) init$truspec[j, , init$rnd[i]])
      plot(init$omegas, spec_j, "l", xlab = "", ylab = "", main = j, ylim = ylim)
      if (!is.null(init$rnd)) lines(init$omegas, init$truspec[j, , init$rnd[i]], type = "l", lty = 1, col = "blue")
      title(ifelse(is.null(res$which.spec), i, res$which.spec[i]), outer = TRUE)
    }
  }
}

.onUnload <- function(libpath) {
  library.dynam.unload("NMSDE", libpath)
}

