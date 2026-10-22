# Generate a cohort with L < X by construction and independent R; retain X < R.
set.seed(2026)
n <- 200
Z <- rnorm(n)
L <- exp(1 + 0.3 * Z + rnorm(n, sd = 0.4))
X <- L + runif(n, 0.5, 1.5)
R <- runif(n, 3, 12)
example.data <- data.frame(L, X, R, Z)
example.data <- example.data[example.data$X < example.data$R, ]

# Cox model: simple and modified pseudo-observation regression
seqTrun.POreg.cox(example.data)
seqTrun.modPOreg.cox(example.data)

# AFT model: simple and modified pseudo-observation regression
seqTrun.POreg.aft(example.data)
seqTrun.modPOreg.aft(example.data)
