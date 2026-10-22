#' Pseudo-Observation AFT Regression under Sequential Truncation
#'
#' Estimate the coefficient of a single numeric covariate in an accelerated
#' failure time model for an event time subject to sequential truncation.
#'
#' @param data A data frame or numeric matrix with at least three rows. The first
#'   four columns contain finite numeric values of \eqn{(L, X, R, Z)}, in that
#'   order. The response \eqn{L} must be strictly positive for the log transform,
#'   \eqn{Z} is a single numeric covariate, and every row must satisfy
#'   \eqn{L < X < R}. Column names are optional; extra columns are ignored.
#' @return A named list with \code{Coefficient estimate}, the coefficient of
#'   \eqn{Z} in the model for \eqn{\log(L)}, and \code{Standard Error}, its
#'   approximate jackknife standard error from \code{geepack::geese}.
#' @details
#' The sequential-truncation model has \eqn{P(L < X) = 1}, and only observations
#' satisfying \eqn{L < X < R} are sampled. See the reference for the independence
#' and identifiability assumptions and the conditions for each method.
#'
#' \code{seqTrun.POreg.aft} constructs leave-one-out pseudo-observations for the
#' estimated mean of \eqn{\log(L)}, using trapezoidal integration of the estimated
#' survival function on the observed log-time grid. The simple method can be
#' biased when truncation depends on the covariate.
#'
#' \code{seqTrun.modPOreg.aft} regresses \eqn{\log(L)} on \eqn{Z} using inverse
#' estimated truncation-survival weights. Estimated probabilities are bounded
#' below by \eqn{10^{-4}} when constructing these weights. Both functions use
#' Gaussian GEE with an identity mean link, working independence, and one cluster
#' per subject. Truncation survival probabilities use the linear interpolation
#' in the supplied regression implementation.
#' @references
#' Qian, J., Parner, E. T., Overgaard, M. and Betensky, R. A. (2026).
#' Pseudo-observation regression for sequentially truncated data.
#' \emph{Biometrics}, 82(2), ujag084. \doi{10.1093/biomtc/ujag084}.
#' @seealso \code{\link{seqTrun.POreg.cox}}, \code{\link{seqTrunNPMLE}}
#' @examples
#' # See ?seqTrun.POreg.cox for a simulated example using all four methods.
#' @export
seqTrun.POreg.aft <- function(data) {
  temp.dat <- .seqTrun_po_data(data, log.time = TRUE)
  n <- nrow(temp.dat)
  theta <- .seqTrun_po_mean_log(temp.dat$L, temp.dat$X, temp.dat$R)
  po.logT <- numeric(n)
  for (i in seq_len(n)) {
    dat2 <- temp.dat[-i, , drop = FALSE]
    theta.i <- .seqTrun_po_mean_log(dat2$L, dat2$X, dat2$R)
    po.logT[i] <- n * theta - (n - 1) * theta.i
  }
  dat.long <- data.frame(temp.dat, po.logT = po.logT, id = seq_len(n))
  fit <- geepack::geese(
    stats::as.formula("po.logT ~ Z"),
    id = dat.long$id, data = dat.long, scale.fix = TRUE,
    family = stats::gaussian(), jack = TRUE,
    mean.link = "identity", corstr = "independence"
  )
  .seqTrun_po_result(fit)
}

#' @rdname seqTrun.POreg.aft
#' @export
seqTrun.modPOreg.aft <- function(data) {
  temp.dat <- .seqTrun_po_data(data, log.time = TRUE)
  po.wt <- 1 / pmax(.seqTrun_po_survival(temp.dat$X, temp.dat$R), 1e-4)
  dat.long <- data.frame(temp.dat, logL = log(temp.dat$L),
                         po.wt = po.wt, id = seq_len(nrow(temp.dat)))
  fit <- geepack::geese(
    stats::as.formula("logL ~ Z"),
    id = dat.long$id, data = dat.long, scale.fix = TRUE,
    family = stats::gaussian(), jack = TRUE, weights = dat.long$po.wt,
    mean.link = "identity", corstr = "independence"
  )
  .seqTrun_po_result(fit)
}
