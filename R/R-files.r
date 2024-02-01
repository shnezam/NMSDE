#' @importFrom graphics lines par plot title
#' @importFrom stats fft
#' @importFrom Rcpp sourceCpp
#' @useDynLib NMSDE, .registration = TRUE


# initial -----------------------------------------------------------------

#' @export initial
initial <- function(X, nbasis = 10, norder = 4, R = length(input$model), 
                    pen = "diff", a = 2, theta_P2 = T, tru = NULL, f.r = NULL) {
  # set up some basic values
  P <- nrow(X[[1]])
  P2 <- P^2
  Nreal <- P * (P + 1) / 2
  len <- ncol(X[[1]])
  m <- length(X)
  K <- floor((len - 1) / 2)
  omegas <- (1:K) / len
  
  if (is.null(f.r)) {
    f.r <- c(1, length(omegas))
  } else {
    omegas <- omegas[f.r[1]:f.r[2]]
    K <- length(omegas)
  }
  
  indexer <- Cholindexer(P)
  
  # Bases(B splines)[K*L] and Penalty matrix [L*L]
  knots <- seq(omegas[1] - (1 / len), omegas[length(omegas)] + (1 / len),
               length.out = nbasis - norder + 2)
  B <- fda::bsplineS(omegas, breaks = knots, norder = norder)
  #B1 <- pbs(omegas, knots = knots, degree = 3, intercept = FALSE,
  #Boundary.knots = range(omegas), periodic = TRUE)
  
  # Define parameters for the Fourier basis
  
  # Create a Fourier basis
  # knots <- seq(omegas[1] - (1 / len), omegas[length(omegas)] + (1 / len), 
  #              length.out = nbasis)
  # fourierbasis <- fda::create.fourier.basis(rangeval = c(min(knots), max(knots)), nbasis = length(knots))
  # B <- eval.basis(omegas, fourierbasis)
  # nbasis <- ncol(B)
  
  # Penalty functions
  if (pen == "diff") {
    if (a != 0) {
      L <- diff(diag(nbasis), differences = 2)
    } else {
      L <- diag(nbasis)
    }
    D <- t(L) %*% L
  } else {
    D <- fda::bsplinepen(fda::create.bspline.basis(nbasis = nbasis, norder = norder))
  }
  
  init <- list(P = P, m = m, K = K, omegas = omegas, f.r = f.r, B = B, D = D, 
               indexer = indexer, Nreal = Nreal, knots = knots)
  
  # some required derivatives for Fisher-scoring algorithm
  init$drivG <- drivG(init)
  
  # go to frequency domain from time domain (DFTs and Periodograms)
  init$freqdom <- lapply(X, time2freq, init)
  
  # Nadaraya-Watson kernel estimate of spectral matrix
  NWF <- lapply(X, NWspec, init)
  init$NW <- NWF
  # Cholesky components of inverse of Nadaraya-Watson estimate
  NWG_m <- sapply(NWF, spec2inChol, init, simplify = "array")
  
  # get initial values by projecting the Cholesky components (NW) 
  #onto the linear subspace (of the basis)
  init <- append(init, G2tA(NWG_m, theta_P2, init, R))
  
  # elbow
  merg.NWF <- matrix(0, nrow = init$m, ncol = init$P^2 * init$K)
  for (i in 1:init$m) merg.NWF[i, ] <- as.vector(t(NWF[[i]]))
  #init$elbow<-fviz_nbclust(merg.NWF,FUNcluster=hcut,method="wss",k.max=min(12,init$m-1))+geom_vline(xintercept=R,linetype=2)
  
  
  # get true spectrum
  if (!is.null(tru)) {
    init$truspec <- sapply(tru$parasmdl, psd, init, simplify = "array")
    
    init$rnd <- tru$rnd
  }
  
  return(init)
}

# time2freq ---------------------------------------------------------------

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

#' @export NWspec
NWspec <- function(x, init) {
  logh <- 2 * log(20 / init$K)
  h <- exp(logh)
  hatF <- matrix(0, nrow = init$P^2, ncol = init$K)
  
  # smooth the periodogram
  for (k in 1:init$K) {
    z <- kerf((init$omegas[k] - init$omegas) / h)
    hatF[, k] <- time2freq(x, init)$perds %*% z / sum(z)
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

#' @export inChol2spec
inChol2spec <- function(inChol, init) {
  
  P2 <- init$P^2
  spec <- matrix(0, nrow = P2, ncol = init$K)
  
  for (k in 1:init$K) {
    
    # inv cholesky at freq k
    cholMat <- matrix(0, nrow = init$P, ncol = init$P)
    
    for (p in 1:init$Nreal){
      cholMat[init$indexer[p, 3], init$indexer[p, 4]] <- inChol[p, k]
    }
    
    for (p in (init$Nreal + 1):P2){
      cholMat[init$indexer[p, 3], init$indexer[p, 4]] <- cholMat[init$indexer[p, 3], init$indexer[p, 4]] + 1i * inChol[p, k]
    }
    
    # spec matrix at freq k
    F_k <- solve(cholMat %*% Conj(t(cholMat)))
    
    for (p in 1:init$Nreal){ 
      spec[p, k] <- Re(F_k[init$indexer[p, 3], init$indexer[p, 4]])
    }
    
    for (p in (init$Nreal + 1):P2){ 
      spec[p, k] <- Im(F_k[init$indexer[p, 3], init$indexer[p, 4]])
    }
    
  }
  
  return(spec)
}


# psd ---------------------------------------------------------------------

# True power spectral density(PSD) for VARMA(p,q)

#' @export psd
psd <- function(paras, init) {
  Phi <- paras[[1]]
  Theta <- paras[[2]]
  sigma <- paras[[3]]
  P2 <- init$P^2
  
  # True spectral density function for VARMA(p,q) 
  F <- function(w) {
    A <- B <- matrix(0, nrow = init$P, ncol = init$P)
    for (j in 1:(length(Phi))){
      A <- A + Phi[[j]] * exp(-2 * pi * 1i * w * (j-1))
    }
    for (j in 1:(length(Theta))){
      B <- B + Theta[[j]] * exp(-2 * pi * 1i * w * (j-1))
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

#' @export loglike
#'
loglike <- function(thetas, As, lambda, init) {
  
  hatG_m <- tA2G(thetas, As, init)
  l <- matrix(0, nrow = init$K, ncol = nrow(As))
  for (i in 1:nrow(As)) {
    for (k in 1:init$K) {
      # get the inverse cholesky matrix
      Cholmat_i <- matrix(0, nrow = init$P, ncol = init$P)
      for (p in 1:init$Nreal) Cholmat_i[init$indexer[p, 3], init$indexer[p, 4]] <- hatG_m[p, k, i]
      for (p in (init$Nreal + 1):init$P^2) Cholmat_i[init$indexer[p, 3], init$indexer[p, 4]] <- Cholmat_i[init$indexer[p, 3], init$indexer[p, 4]] + 1i * hatG_m[p, k, i]
      
      # get the loglikelihood
      l[k, i] <- -log(prod(eigen(Cholmat_i %*% Conj(t(Cholmat_i)), only.values = T)$values)) + Conj(t(init$freqdom[[i]]$tildeX[, k])) %*% Cholmat_i %*% Conj(t(Cholmat_i)) %*% init$freqdom[[i]]$tildeX[, k]
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
G2tA <- function(G_m, theta_P2, init, R) {
  As <- array(0, dim = c(init$m, R, init$P^2))
  smooth <- solve(t(init$B) %*% init$B) %*% t(init$B)
  theta.tA <- array(0, dim = c(ncol(init$B), init$m, init$P^2))
  
  for (j in 1:init$P^2) theta.tA[, , j] <- smooth %*% G_m[j, , ]
  if (theta_P2 == T) {
    svd.theta.tA <- list()
    thetas <- array(0, dim = c(ncol(init$B), R, init$P^2))
    
    for (j in 1:init$P^2) {
      svd.theta.tA[[j]] <- svd(theta.tA[, , j])
      thetas[, , j] <- as.matrix(svd.theta.tA[[j]]$u[, 1:R])
      As[, , j] <- as.matrix(svd.theta.tA[[j]]$v %*% diag(svd.theta.tA[[j]]$d, ncol(svd.theta.tA[[j]]$v)))[, 1:R]
    }
    
  } else {
    mergG <- matrix(0, init$K, init$m * init$P^2)
    
    for (j in 1:init$P^2) mergG[, (((j - 1) * (init$m)) + 1):(j * init$m)] <- G_m[j, , ]
   
    theta.tmergA <- smooth %*% mergG
    svd.theta.tmergA <- svd(theta.tmergA)
    thetas <- as.matrix(svd.theta.tmergA$u[, 1:R])
    mergA <- as.matrix(svd.theta.tmergA$v %*% diag(svd.theta.tmergA$d, ncol(svd.theta.tmergA$v)))
    mergA <- as.matrix(mergA[, 1:R])
    
    for (j in 1:init$P^2) As[, , j] <- mergA[((j - 1) * init$m + 1):(j * init$m), ]
  }
  
  return(list(thetas = thetas, As = As, theta.tA = theta.tA))
}


# tA2G --------------------------------------------------------------------

# Takes thetas and As and reconstructs G
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

# Workingfish -------------------------------------------------------------

Workingfish <- function(thetas = ini$thetas, As = ini$As, lambda, init) {
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
      tildeX_k_i <- init$freqdom[[i]]$tildeX[, k]
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
  
  while (repit) {
    tau <- (0.5)^counter
    for (i in 1:nrow(As)) curr.As[i, , ] <- matrix(A_R.P2_m[, i] - tau * A.chng[, i], nrow = R, ncol = P2)

    if (!is.matrix(thetas)) {
      for (r in 1:R) curr.thetas[, r, ] <- matrix(theta_L.P2_R[, r] - tau * theta.chng[, r], nrow = nbasis, ncol = P2)
    } else {
      for (r in 1:R) curr.thetas[, r] <- thetas[, r] - tau * theta.chng[, r]
    }
    
    if (loglike(curr.thetas, curr.As, lambda, init)$ploglike <= loglike(thetas, As, lambda, init)$ploglike) {
      repit <- FALSE
    } else {
      counter <- counter + 1
    }
    
  }
  
  return(list(thetas = curr.thetas, As = curr.As, DF = Re(DF), stepH = counter))
}



# iterate -----------------------------------------------------------------

# Iterartive algorithm
#' @export iterate
iterate <- function(init, lambda = 0, a = 2, eps = 0.5, n.iter = 100, which.spec = NULL, update.lam = T) {
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
    init$freqdom <- init$freqdom[which.spec]
  }
  
  thetas.lst[[l]] <- Re(init$thetas)
  As.lst[[l]] <- Re(init$As)
  
  Lambda[[l]] <- lambda
  
  if (is.matrix(init$thetas)) {
    ploglikes[l] <- loglike(thetas.lst[[l]], As.lst[[l]], Lambda[[l]], init)$ploglike
  } else{
    ploglikes[l] <- loglike_c_cube(thetas.lst[[l]], As.lst[[l]], Lambda[[l]], init)$ploglike
  }
  
  AIC[l] <- ploglikes[l]
  diff <- eps + ncol(init$thetas)
  
  while (abs(diff) > eps && l < n.iter) {
    if (is.matrix(init$thetas)) {
      workvals <- Workingfish_mat(thetas.lst[[l]], As.lst[[l]], Lambda[[l]], init)
    } else {
      workvals <- Workingfish_cube(thetas.lst[[l]], As.lst[[l]], Lambda[[l]], init)
      print(workvals$stepH)
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
      #Lambda[[l + 1]] <- Lambda[[l]]
      Lambda[[l + 1]] <- as.vector(Lambda[[l]])
    }
    #ploglikes[l + 1] <- loglike(thetas.lst[[l + 1]], As.lst[[l + 1]], Lambda[[l + 1]], init)$ploglike;
    if (is.matrix(init$thetas)) {
      ploglikes[l + 1] <- loglike(thetas.lst[[l + 1]], As.lst[[l + 1]], Lambda[[l + 1]], init)$ploglike
    } else{
      ploglikes[l + 1] <- loglike_c_cube(thetas.lst[[l + 1]], As.lst[[l + 1]], Lambda[[l + 1]], init)$ploglike
    }
    
    AIC[(l + 1)] <- Re(ploglikes[l + 1] + 2 * sum(workvals$DF))
    diff <- ploglikes[l] - ploglikes[l + 1]
    l <- l + 1; #
    #incProgress(1/n.iter, detail = paste("iteration", l))
    
    cat(paste("\n l=", (l), " pLogL=", round(ploglikes[(l)]), " AIC=", round(AIC[(l)]), " lambda=", round(Lambda[[(l)]], 2), " DF=", round(workvals$DF), " stepH=", workvals$stepH, " diff=", round(diff, 2), sep = ""), "\n")
  }
  
  return(list(thetas = thetas.lst, As = As.lst, ploglikes = ploglikes, Lambda = Lambda, AIC = AIC, r.iter = l, which.spec = which.spec))
}


# plot.sds ----------------------------------------------------------------

# Plots estimate of the components spectral matrix and Trues
#' @export plot_sds
plot_sds <- function(res, init) {
  par(mfrow = c(init$P, init$P), oma = c(0.5, 0.5, 2.5, 0.5))
  if (!is.null(res$which.spec)) {
    init$m <- length(res$which.spec)
    if (!is.null(init$rnd)) init$rnd <- init$rnd[res$which.spec]
  }
  
  for (i in 1:init$m) {
    for (j in 1:init$P^2) {
      
      spec_j <- inChol2spec(tA2G(res$thetas[[length(res$thetas)]], res$As[[length(res$As)]], init)[, , i], init)[j, ]
      plot(init$omegas, spec_j, "l", xlab = "", ylab = "", main = j)
      if (!is.null(init$rnd)) lines(init$omegas, init$truspec[j, , init$rnd[i]], type = "l", lty = 1, col = "blue")
      title(ifelse(is.null(res$which.spec), i, res$which.spec[i]), outer = TRUE)
      ylim = c(0, max(init$truspec[j, , init$rnd[i]],spec_j))
    }
  }
}

.onUnload <- function(libpath) {
  library.dynam.unload("combinIT", libpath)
}