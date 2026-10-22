#' Pseudo-Observation Cox Regression under Sequential Truncation
#'
#' Estimate the coefficient of a single numeric covariate in a Cox proportional
#' hazards model for an event time subject to sequential truncation.
#'
#' @param data A data frame or numeric matrix with at least three rows. The first
#'   four columns contain finite numeric values of \eqn{(L, X, R, Z)}, in that
#'   order, where \eqn{L} is the response event time and \eqn{Z} is a single
#'   covariate. Every observed row must satisfy \eqn{L < X < R}. Column names are
#'   optional; additional columns are ignored.
#' @return A named list with \code{Coefficient estimate}, the coefficient of
#'   \eqn{Z}, and \code{Standard Error}, its approximate jackknife standard error
#'   from \code{geepack::geese}.
#' @details
#' The sequential-truncation model has \eqn{P(L < X) = 1}, and only observations
#' satisfying \eqn{L < X < R} are sampled. See the reference for the independence
#' and identifiability assumptions and the conditions for each regression method.
#'
#' \code{seqTrun.POreg.cox} constructs leave-one-out survival pseudo-observations
#' at five empirical quantiles of \eqn{L}, with probabilities equally spaced
#' from 0.1 to 0.8. Their complements are fitted using Gaussian GEE with a
#' complementary log-log mean link. The simple method can be biased when
#' truncation depends on the covariate.
#'
#' \code{seqTrun.modPOreg.cox} fits event indicators using binomial GEE and
#' inverse estimated truncation-survival weights. The estimated probabilities
#' are bounded below by \eqn{10^{-4}} when constructing these weights.
#' Both methods use a separate intercept at each distinct evaluation time,
#' working independence, and subject identifiers to group repeated observations.
#' Truncation survival probabilities use the linear interpolation in the supplied
#' regression implementation.
#' @references
#' Qian, J., Parner, E. T., Overgaard, M. and Betensky, R. A. (2026).
#' Pseudo-observation regression for sequentially truncated data.
#' \emph{Biometrics}, 82(2), ujag084. \doi{10.1093/biomtc/ujag084}.
#' @seealso \code{\link{seqTrun.POreg.aft}}, \code{\link{seqTrunNPMLE}}
#' @example inst/examples/ex_seqTrunPOReg.R
#' @export
seqTrun.POreg.cox <- function(data) {
  temp.dat <- .seqTrun_po_data(data)
  n <- nrow(temp.dat)
  t0 <- .seqTrun_po_times(temp.dat$L)
  surv <- 1 - .seqTrun_po_cdf(temp.dat$L, temp.dat$X, temp.dat$R, t0)
  rho <- matrix(NA_real_, n, length(t0))
  for (i in seq_len(n)) {
    dat2 <- temp.dat[-i, , drop = FALSE]
    surv2 <- 1 - .seqTrun_po_cdf(dat2$L, dat2$X, dat2$R, t0)
    rho[i, ] <- n * surv - (n - 1) * surv2
  }
  dat.long <- do.call(rbind, lapply(seq_along(t0), function(j) {
    data.frame(temp.dat, ipseudo = 1 - rho[, j],
               tpseudo = t0[j], id = seq_len(n))
  }))
  dat.long <- dat.long[order(dat.long$id), ]
  fit <- geepack::geese(
    stats::as.formula("ipseudo ~ as.factor(tpseudo) + Z"),
    id = dat.long$id, data = dat.long, scale.fix = TRUE,
    family = stats::gaussian(), jack = TRUE,
    mean.link = "cloglog", corstr = "independence"
  )
  .seqTrun_po_result(fit)
}

#' @rdname seqTrun.POreg.cox
#' @export
seqTrun.modPOreg.cox <- function(data) {
  temp.dat <- .seqTrun_po_data(data)
  n <- nrow(temp.dat)
  t0 <- .seqTrun_po_times(temp.dat$L)
  po.wt <- 1 / pmax(.seqTrun_po_survival(temp.dat$X, temp.dat$R), 1e-4)
  dat.long <- do.call(rbind, lapply(seq_along(t0), function(j) {
    data.frame(temp.dat, irho = as.numeric(temp.dat$L <= t0[j]),
               trho = t0[j], po.wt = po.wt, id = seq_len(n))
  }))
  dat.long <- dat.long[order(dat.long$id), ]
  fit <- geepack::geese(
    stats::as.formula("irho ~ as.factor(trho) + Z"),
    id = dat.long$id, data = dat.long, scale.fix = TRUE,
    family = stats::binomial(), jack = TRUE, weights = dat.long$po.wt,
    mean.link = "cloglog", corstr = "independence"
  )
  .seqTrun_po_result(fit)
}
