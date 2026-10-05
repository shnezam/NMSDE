X <- dat$X1[1:4]
init <- initial(X, nbasis = 8, R = 2, tru = datPar$modelPar1)

test_that("Chol gives a lower-triangular factor of a Hermitian matrix", {
  f <- matrix(c(2, 1 - 1i, 1 + 1i, 3), 2, 2)
  G <- Chol(f)
  expect_equal(G %*% Conj(t(G)), f)
  expect_equal(unname(G[1, 2]), 0 + 0i)
})

test_that("spec2inChol and inChol2spec are inverse maps", {
  expect_equal(inChol2spec(spec2inChol(init$NW[[1]], init), init), init$NW[[1]])
})

test_that("initial returns objects of the documented dimensions", {
  L <- 8; R <- 2; P2 <- 4
  expect_equal(dim(init$thetas), c(L, R, P2))
  expect_equal(dim(init$As), c(length(X), R, P2))
  expect_equal(dim(init$B), c(init$K, L))
  expect_equal(dim(init$truspec), c(P2, init$K, 3))
  expect_length(init$elbow, length(X))
  expect_error(initial(X, nbasis = 8, R = 5), "R")
  expect_error(initial(list(matrix(rnorm(50), 1)), R = 1), "P >= 2")
})

test_that("the penalty order 'a' is used", {
  D3 <- initial(X, nbasis = 8, R = 2, a = 3)$D
  expect_equal(qr(D3)$rank, 8 - 3)
})

test_that("C++ and R likelihoods agree", {
  r <- loglike(init$thetas, init$As, 2, init, init$tildeX)
  c <- loglike_c_cube(init$thetas, init$As, 2, init, init$tildeX)
  expect_equal(r, c)
})

test_that("iterate decreases the penalized likelihood with fixed lambda", {
  for (tp in c(TRUE, FALSE)) {
    ini <- initial(X, nbasis = 8, R = 2, theta_P2 = tp)
    fit <- iterate(ini, lambda = 1, n.iter = 4, update.lam = FALSE, verbose = FALSE)
    expect_true(all(diff(fit$ploglikes) <= 1e-8))
  }
})

test_that("estimates are Hermitian positive definite", {
  fit <- iterate(init, lambda = 1, n.iter = 3, verbose = FALSE)
  G <- tA2G(fit$thetas[[fit$r.iter]], fit$As[[fit$r.iter]], init)
  s <- inChol2spec(G[, , 1], init)
  f <- matrix(c(s[1, 1], s[2, 1] + 1i * s[4, 1], s[2, 1] - 1i * s[4, 1], s[3, 1]), 2, 2)
  expect_true(all(eigen(f, only.values = TRUE)$values > 0))
})
