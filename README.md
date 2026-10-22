# seqTrun

## Estimation and regression with sequentially truncated survival data

`seqTrun` implements nonparametric and semiparametric maximum likelihood
estimation for sequentially truncated survival data, together with simple and
modified pseudo-observation regression for Cox and accelerated failure time
(AFT) models.

## Installation

```r
install.packages("remotes") # if needed
remotes::install_github("jqian-lab/seqTrun")
library(seqTrun)
```

The package imports `survival` and `geepack`; these dependencies are installed
automatically.

## Available functions

| Function | Purpose |
| --- | --- |
| `seqTrunNPMLE()` | Nonparametric estimation of the distribution of L |
| `seqTrunSPMLE()` | Semiparametric estimation of the distribution of L |
| `seqTrun.POreg.cox()` | Simple pseudo-observation Cox regression |
| `seqTrun.modPOreg.cox()` | Modified pseudo-observation Cox regression |
| `seqTrun.POreg.aft()` | Simple pseudo-observation AFT regression |
| `seqTrun.modPOreg.aft()` | Modified pseudo-observation AFT regression |

## Pseudo-observation regression

The regression functions accept a data frame or numeric matrix whose first
four columns are `(L, X, R, Z)`, in that order. Each observed row must satisfy
`L < X < R`. `L` is the response event time and `Z` is a single numeric covariate;
the AFT functions require `L > 0`. Column names are optional.

```r
set.seed(2026)
n <- 200
Z <- rnorm(n)
L <- exp(1 + 0.3 * Z + rnorm(n, sd = 0.4))
X <- L + runif(n, 0.5, 1.5)
R <- runif(n, 3, 12)
example.data <- data.frame(L, X, R, Z)
example.data <- example.data[example.data$X < example.data$R, ]

cox.simple   <- seqTrun.POreg.cox(example.data)
cox.modified <- seqTrun.modPOreg.cox(example.data)
aft.simple   <- seqTrun.POreg.aft(example.data)
aft.modified <- seqTrun.modPOreg.aft(example.data)

cox.modified[["Coefficient estimate"]]
cox.modified[["Standard Error"]]
```

Each regression function returns the coefficient estimate for `Z` and its
approximate jackknife standard error from `geepack::geese()`. The Cox functions
use five empirical evaluation quantiles between the 10th and 80th percentiles
of `L`, with repeated observations grouped by subject. The modified functions
use inverse estimated truncation-survival weights with a probability floor
of `1e-4`.

The simple approach can be biased when truncation depends on the covariate.
See Qian et al. (2026) and the package help pages for the assumptions required
by each method.

## References

* Betensky, R. A., Qian, J. and Hou, J. (2023). Nonparametric and semiparametric
  estimation with sequentially truncated survival data. *Biometrics*,
  79(2), 1000–1013.
* Qian, J., Parner, E. T., Overgaard, M. and Betensky, R. A. (2026).
  Pseudo-observation regression for sequentially truncated data. *Biometrics*,
  82(2), ujag084. [doi:10.1093/biomtc/ujag084](https://doi.org/10.1093/biomtc/ujag084).
