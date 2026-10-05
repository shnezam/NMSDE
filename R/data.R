#' Simulated collections of bivariate time series
#'
#' Example data used by \code{\link{launch_app}} and in the examples. Each
#' collection contains bivariate series simulated from one of three vector
#' autoregressive (VAR) models with standard normal innovations of
#' correlation 0.9; the true model parameters are in \code{\link{datPar}}.
#'
#' @format \code{dat} is a list of four collections. Each collection is a
#'   list of \eqn{2 \times T}{2 x T} numeric matrices (rows are the two
#'   component series, columns are time points), in the form expected by
#'   \code{\link{initial}}:
#'   \describe{
#'     \item{X1}{6 series of length 200.}
#'     \item{X2}{15 series of length 200.}
#'     \item{X3}{15 series of length 200 (true parameters not stored).}
#'     \item{X4}{30 series of length 400.}
#'   }
#' @source Simulated by the package authors.
#' @examples
#' str(dat, max.level = 1)
#' matplot(t(dat$X1[[1]]), type = "l", lty = 1, ylab = "")
"dat"

#' True VAR parameters of the simulated data
#'
#' The models used to simulate \code{\link{dat}}, in the form expected by
#' argument \code{tru} of \code{\link{initial}}.
#'
#' @format A list of four elements, \code{modelPar1}, \dots,
#'   \code{modelPar4}, matching \code{dat$X1}, \dots, \code{dat$X4}. Each
#'   element (except \code{modelPar3}, which is \code{list(NULL)}) is a list
#'   with components
#'   \describe{
#'     \item{paras}{list of the three VAR models; each model is a list of
#'       \eqn{2 \times 2}{2 x 2} lag-polynomial matrices
#'       \eqn{I, -\Phi_1, \dots, -\Phi_d}.}
#'     \item{sigma}{\eqn{2 \times 2}{2 x 2} innovation covariance matrix.}
#'     \item{rnd}{index of the model that generated each series.}
#'   }
#' @source Simulated by the package authors.
#' @examples
#' datPar$modelPar1$paras[[2]]
#' table(datPar$modelPar1$rnd)
"datPar"
