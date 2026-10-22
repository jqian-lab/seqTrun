# Internal helpers for pseudo-observation regression. Keep these separate from
# the NPMLE/SPMLE helpers so that regression code cannot replace them.
.seqTrun_po_data <- function(data, log.time = FALSE) {
  if ((!is.data.frame(data) && !is.matrix(data)) ||
      ncol(data) < 4L || nrow(data) < 3L) {
    stop("data must contain at least three rows and four columns (L, X, R, Z).",
         call. = FALSE)
  }
  data <- as.data.frame(data[, seq_len(4L), drop = FALSE])
  names(data) <- c("L", "X", "R", "Z")
  if (!all(vapply(data, is.numeric, logical(1))) ||
      any(!is.finite(as.matrix(data)))) {
    stop("The first four columns (L, X, R, Z) must be finite numeric values.",
         call. = FALSE)
  }
  if (any(data$L >= data$X | data$X >= data$R)) {
    stop("Every observed row must satisfy L < X < R.", call. = FALSE)
  }
  if (log.time && any(data$L <= 0)) {
    stop("L must be positive for AFT regression of log(L).", call. = FALSE)
  }
  if (length(unique(data$Z)) < 2L) {
    stop("Z must contain at least two distinct values.", call. = FALSE)
  }
  data
}

.seqTrun_po_survival <- function(x, r) {
  fit <- survival::survfit(
    survival::Surv(x, r, event = rep(1, length(x))) ~ 1,
    timefix = FALSE
  )
  if (length(unique(fit$time)) < 2L) {
    stop("At least two distinct R times are needed for survival interpolation.",
         call. = FALSE)
  }
  # Retain the linear interpolation used in the supplied regression scripts.
  stats::approx(fit$time, fit$surv, xout = x,
                method = "linear", rule = 2)$y
}

.seqTrun_po_times <- function(l) {
  t0 <- as.numeric(stats::quantile(l, seq(0.1, 0.8, length.out = 5)))
  if (length(unique(t0)) < 2L) {
    stop("Cox regression requires at least two distinct evaluation quantiles of L.",
         call. = FALSE)
  }
  t0
}

.seqTrun_po_cdf <- function(l, x, r, a) {
  s.r.x <- .seqTrun_po_survival(x, r)
  if (any(!is.finite(s.r.x) | s.r.x <= 0)) {
    stop("Estimated truncation survival probabilities must be positive.",
         call. = FALSE)
  }
  weights <- 1 / s.r.x
  vapply(a, function(t) sum((l <= t) * weights) / sum(weights), numeric(1))
}

.seqTrun_po_mean_log <- function(l, x, r) {
  tt <- sort(l)
  surv <- 1 - .seqTrun_po_cdf(l, x, r, tt)
  # Trapezoidal integration on the log-time scale, as in the supplied AFT code.
  log(tt[1L]) + sum((surv[-length(tt)] + surv[-1L]) * diff(log(tt)) / 2)
}

.seqTrun_po_result <- function(fit) {
  if (fit$error != 0L) {
    warning("The GEE fit did not converge; interpret the returned estimates with care.",
            call. = FALSE)
  }
  coef.table <- summary(fit)$mean
  list("Coefficient estimate" = unname(coef.table["Z", "estimate"]),
       "Standard Error" = unname(coef.table["Z", "ajs.se"]))
}
