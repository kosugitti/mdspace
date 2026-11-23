# mdspace <img src="man/figures/logo.png" align="right" height="139" alt="mdspace logo" />

Multidimensional scaling (MDS) in R, with a focus on visualizing asymmetric
proximities. The package is under development.

## Installation

```r
# install.packages("remotes")
remotes::install_github("kosugitti/mdspace")
```

## Models

| Function | Model | Reference |
|---|---|---|
| `circle_mds()` | Circle model (radius-distance model); one-mode two-way and two-mode three-way | Okada & Imaizumi (1987, 1997) |
| `ellipse_mds()` | Ellipse model | Okada (1990) |
| `vonmises_mds()` | Von Mises model (AMISESCAL) | Shojima (2011, 2012) |
| `hfm_mds()` | Hermitian form model (HFM) | Chino & Shiraiwa (1993) |
| `gstatis_mds()` | GSTATIS | Chino, Grorud & Yoshino (1996) |

Helpers: `as_proximity()`, `asym_decomp()`, `missing_pattern()`.

## Example

```r
library(mdspace)

# x: a square matrix of asymmetric similarities
fit <- circle_mds(x, ndim = 2)
fit
plot(fit)
```

## License

MIT
