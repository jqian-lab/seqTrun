library(seqTrun)

# Results were checked against the supplied Cox/AFT scripts after correcting
# their global-data references and removing the trailing placeholder functions.
set.seed(2026)
n <- 200
Z <- rnorm(n)
L <- exp(1 + 0.3 * Z + rnorm(n, sd = 0.4))
X <- L + runif(n, 0.5, 1.5)
R <- runif(n, 3, 12)
test.data <- data.frame(L, X, R, Z)
test.data <- test.data[test.data$X < test.data$R, ]

methods <- c("seqTrun.POreg.cox", "seqTrun.modPOreg.cox",
             "seqTrun.POreg.aft", "seqTrun.modPOreg.aft")
expected <- list(
  c(-0.832573598703963, 0.136599158671311),
  c(-0.802623819022748, 0.112663430520492),
  c(0.303528275927412, 0.0450479516432239),
  c(0.279208774824217, 0.030882410357695)
)

renamed <- as.matrix(test.data)
colnames(renamed) <- c("onset", "entry", "exit", "covariate")
set.seed(12)
permuted <- test.data[sample(nrow(test.data)), ]

for (j in seq_along(methods)) {
  fun <- get(methods[j])
  actual <- fun(test.data)
  stopifnot(identical(names(actual), c("Coefficient estimate", "Standard Error")),
            all(is.finite(unlist(actual))), actual[[2L]] > 0,
            isTRUE(all.equal(unname(unlist(actual)), expected[[j]],
                             tolerance = 1e-6)),
            isTRUE(all.equal(actual, fun(renamed), tolerance = 1e-7)),
            isTRUE(all.equal(actual, fun(permuted), tolerance = 1e-7)))
}

# The modified AFT coefficient also solves an ordinary weighted least-squares
# problem under the identity link and working independence.
surv.fit <- survival::survfit(
  survival::Surv(X, R, rep(1, nrow(test.data))) ~ 1, data = test.data,
  timefix = FALSE
)
prob <- stats::approx(surv.fit$time, surv.fit$surv,
                      xout = test.data$X, method = "linear", rule = 2)$y
weighted.fit <- stats::lm(log(L) ~ Z, data = test.data,
                          weights = 1 / pmax(prob, 1e-4))
stopifnot(isTRUE(all.equal(
  seqTrun.modPOreg.aft(test.data)[[1L]],
  unname(stats::coef(weighted.fit)["Z"]), tolerance = 1e-8
)))

expect_error <- function(expr, pattern) {
  err <- tryCatch({force(expr); NULL}, error = identity)
  stopifnot(inherits(err, "error"), grepl(pattern, conditionMessage(err)))
}
expect_error(seqTrun.POreg.cox(test.data[1:2, ]), "at least three rows")
bad <- test.data
bad$Z[1] <- NA_real_
expect_error(seqTrun.modPOreg.cox(bad), "finite numeric")
bad <- test.data
bad$X[1] <- bad$L[1]
expect_error(seqTrun.POreg.aft(bad), "L < X < R")
bad <- test.data
bad$L[1] <- 0
expect_error(seqTrun.modPOreg.aft(bad), "L must be positive")
bad <- test.data
bad$Z <- 1
expect_error(seqTrun.POreg.cox(bad), "two distinct values")
bad <- test.data
bad$L <- min(test.data$L)
expect_error(seqTrun.POreg.cox(bad), "two distinct evaluation quantiles")
