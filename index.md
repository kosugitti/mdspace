# mdspace

Multidimensional scaling (MDS) in R, with a focus on visualizing
asymmetric proximities. The package is under development.

## Installation

``` r

# install.packages("remotes")
remotes::install_github("kosugitti/mdspace")
```

## Models

| Function | Model | Reference |
|----|----|----|
| [`circle_mds()`](https://kosugitti.github.io/mdspace/reference/circle_mds.md) | Circle model (radius-distance model); one-mode two-way and two-mode three-way | Okada & Imaizumi (1987, 1997) |
| [`ellipse_mds()`](https://kosugitti.github.io/mdspace/reference/ellipse_mds.md) | Ellipse model | Okada (1990) |
| [`vonmises_mds()`](https://kosugitti.github.io/mdspace/reference/vonmises_mds.md) | Von Mises model (AMISESCAL) | Shojima (2011, 2012) |
| [`hfm_mds()`](https://kosugitti.github.io/mdspace/reference/hfm_mds.md) | Hermitian form model (HFM) | Chino & Shiraiwa (1993) |
| [`gstatis_mds()`](https://kosugitti.github.io/mdspace/reference/gstatis_mds.md) | GSTATIS | Chino, Grorud & Yoshino (1996) |

Helpers:
[`as_proximity()`](https://kosugitti.github.io/mdspace/reference/as_proximity.md),
[`asym_decomp()`](https://kosugitti.github.io/mdspace/reference/asym_decomp.md),
[`missing_pattern()`](https://kosugitti.github.io/mdspace/reference/missing_pattern.md).

## Example

``` r

library(mdspace)

# x: a square matrix of asymmetric similarities
fit <- circle_mds(x, ndim = 2)
fit
plot(fit)
```

## License

MIT
