# Nonmetric asymmetric MDS by the radius-distance (circle) model

Each object `j` is a point `x_j` in a `p`-dimensional Euclidean space
with a circle (sphere, hypersphere) of radius `r_j` centered at it. The
asymmetric proximity from `j` to `k` is represented by
`m_jk = d_jk - r_j + r_k`, the distance from the near edge of the circle
of `j` to the far edge of the circle of `k`, which is assumed to be
monotonically decreasing in similarities (increasing in
dissimilarities). The fit is measured by
`S = sqrt(sum((m_jk - mhat_jk)^2) / sum((m_jk - mbar)^2))` (Kruskal's
stress formula 2 applied to `m_jk`), where `mhat_jk` are the disparities
obtained by monotone regression (primary approach to ties). `S` is
minimized over the coordinates and radii.

## Usage

``` r
circle_mds(
  x,
  ndim = 2,
  type = c("similarity", "dissimilarity"),
  ndim_max = ndim,
  nrandom = 10,
  maxit = 1000,
  reltol = 1e-10,
  seed = 1
)
```

## Arguments

- x:

  A square matrix of asymmetric proximities. The diagonal is ignored;
  other missing cells are left out of the fit.

- ndim:

  Number of dimensions.

- type:

  `"similarity"` (default) or `"dissimilarity"`. If `x` has a `type`
  attribute set by
  [`as_proximity()`](https://kosugitti.github.io/mdspace/reference/as_proximity.md),
  that is used.

- ndim_max:

  Dimensionality at which the rational initial configuration is built;
  the solutions are reduced from there down to `ndim`.

- nrandom:

  Number of random starts tried in addition to the rational start (`0`
  uses the rational start only).

- maxit, reltol:

  Passed to [`stats::optim()`](https://rdrr.io/r/stats/optim.html).

- seed:

  Random seed for the random starts. The state of the random number
  generator is restored afterwards. `NULL` uses (and advances) the
  current state of the generator.

## Value

An object of class `c("mdspace_circle", "mdspace")` with `conf` (n x
ndim), `radius`, `stress` (`S`), `niter` (gradient evaluations of the
chosen solution), `mjk` and `disparity` (n x n matrices, `NA` on the
diagonal and missing cells), `path` (`S` at each dimensionality from
`ndim_max` to `ndim` for the rational start), `starts` (`S` of each
start, the rational one first), `best_start` (which start was chosen),
and `converged`.

## Details

The initial configuration is "rational" as in Okada and Imaizumi (1987):
Torgerson's method applied to the symmetrized dissimilarities, with the
additive constant that makes them Euclidean (Cailliez, 1983, as in
[`stats::cmdscale()`](https://rdrr.io/r/stats/cmdscale.html) with
`add = TRUE`; Okada & Imaizumi, 1997, add Torgerson's constant), gives
the configuration in `ndim_max` dimensions, and the least squares
solution of `r_k - r_j = a_jk` (`a_jk` the skew-symmetric part of the
dissimilarities) gives the radii. The solution in `ndim_max` dimensions
is then reduced one dimension at a time by its principal axes down to
`ndim`. Because the radius-distance model has local minima, `nrandom`
random starts (10 by default) are also tried in `ndim` dimensions and
the solution with the smallest `S` is returned (Okada & Imaizumi, 2024,
Section 4.1.1).

With few objects the solution may not be unique: a simulation with
`n = 5` (20 off-diagonal cells) reached `S = 0` even with error while
the coordinates differed from the true ones by up to 0.8 (normalized
scale), and one case with `n = 10` in three dimensions did the same. The
configuration should then be interpreted with care. See
`develop/20261117_circle_recovery.R` in the source repository.

The solution is normalized so that the configuration is centered with
the sum of squared coordinates equal to `n`, and the smallest radius is
0 (the radii are determined only up to an additive constant), as in the
nASYMS program.

The minimization uses the analytic gradient of `S` with the disparities
held fixed, which is the gradient of `S` as a function of the
coordinates and radii (Kruskal, 1964b), and
[`stats::optim()`](https://rdrr.io/r/stats/optim.html) (BFGS) instead of
the steepest descent of the original program.

## References

Okada, A., & Imaizumi, T. (1987). Nonmetric multidimensional scaling of
asymmetric proximities. *Behaviormetrika*, 14(21), 81-96.

Okada, A., & Imaizumi, T. (2024). *Applied multidimensional scaling of
asymmetric relationships*. Springer.

Cailliez, F. (1983). The analytical solution of the additive constant
problem. *Psychometrika*, 48(2), 305-308.

Okada, A., & Imaizumi, T. (1997). Asymmetric multidimensional scaling of
two-mode three-way proximities. *Journal of Classification*, 14,
195-224.

Kruskal, J. B. (1964b). Nonmetric multidimensional scaling: A numerical
method. *Psychometrika*, 29(2), 115-129.

## Examples

``` r
x <- matrix(c(NA, 5, 2, 1,
              3, NA, 4, 2,
              1, 2, NA, 6,
              2, 1, 4, NA), 4, byrow = TRUE)
fit <- circle_mds(x, ndim = 1)
fit
#> mdspace: circle 
#> Call: circle_mds(x = x, ndim = 1)
#> Objects: 4  Dimensions: 1 
#> Stress: 0.1915 
#> Iterations: 16 
#> Starts: 11 (rational + 10 random); best = random 8, S range 0.1915 to 0.6396
```
