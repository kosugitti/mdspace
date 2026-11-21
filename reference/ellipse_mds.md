# Nonmetric asymmetric MDS by the ellipse-distance (Model II) and radius-distance (Model III) models

Okada and Imaizumi (1984; see Chino & Okada, 1996, Section 3.5) wrote
the asymmetric proximity from object `j` to object `k` as
`m_jk = d_jk + alpha * c(j, k) + beta * c(k, j)`, where `d_jk` is the
Euclidean distance between the points `x_j` and `x_k`, and named three
families by how the asymmetric term `c` is restricted:

## Usage

``` r
ellipse_mds(
  x,
  ndim = 2,
  model = c("II", "III", "I"),
  alpha = -1,
  beta = 1,
  type = c("similarity", "dissimilarity"),
  ndim_max = ndim,
  nrandom = 10,
  maxit = 2000,
  reltol = 1e-10,
  seed = 1,
  kappa = Inf,
  axis_max = Inf,
  check_orientation = TRUE,
  angles = c(10, 15, 20)
)
```

## Arguments

- x:

  A square matrix of asymmetric proximities. The diagonal is ignored;
  other missing cells are left out of the fit.

- ndim:

  Number of dimensions.

- model:

  `"II"` (ellipse-distance model, default) or `"III"` (radius-distance
  model). `"I"` is reserved and gives an error.

- alpha, beta:

  Signs of the asymmetric terms of the object itself and of the other
  object, each `-1`, `0` or `1` (not both 0).

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

- kappa, axis_max:

  Model II only: optional upper bounds of the aspect ratio
  (`max_t r_jt / min_t r_jt`) and of the size (geometric mean of the
  semiaxes, in units of the root mean square distance of the points from
  the centroid) of each ellipse. The default `Inf` leaves them unbounded
  as in the original method. `S` hardly determines the shape of the
  ellipses: without bounds some ellipses can become very long and thin
  or very large while `S` decreases only slightly, and with bounds the
  solution usually sits on the bound. A warning is given when an
  unbounded solution has an aspect ratio above 10 or a size above 3, and
  when a bound is reached.

- check_orientation:

  Logical; for Model II with `ndim >= 2`, compute `S` after rotating the
  semiaxis configuration (see Details).

- angles:

  Rotation angles in degrees for `check_orientation`.

## Value

An object of class `c("mdspace_ellipse", "mdspace")` with `conf` (n x
ndim), `semiaxis` (n x ndim; Model II) or `radius` (Model III), `stress`
(`S`), `niter`, `mjk` and `disparity` (n x n matrices, `NA` on the
diagonal and missing cells), `path` (`S` at each dimensionality from
`ndim_max` to `ndim` for the rational start), `starts` (`S` of each
start, the rational one first), `best_start`, `converged`, `model`,
`alpha`, `beta`, `kappa`, `axis_max`, `at_bound` (whether some ellipse
reached `kappa` or `axis_max`), and `orientation` (a data frame with the
dimensions, angle, and `S` after counterclockwise and clockwise
rotation; `NULL` when not computed).

## Details

- **Model II** (ellipse-distance model; Okada, 1990): each object `j`
  has an ellipse (ellipsoid, hyperellipsoid) centered at `x_j` with
  semiaxes `r_j1, ..., r_jp` along the dimensions, and `c(j, k)` is the
  radius of the ellipse of `j` in the direction of `k`,
  `rho_jk = d_jk v_jk` with
  `v_jk = 1 / sqrt(sum_t ((x_jt - x_kt) / r_jt)^2)` (Okada & Imaizumi,
  2024, Equations 3.10 to 3.13). When all semiaxes of `j` are equal to
  `r_j`, `rho_jk = r_j` and Model II reduces to Model III.

- **Model III** (radius-distance model; Okada & Imaizumi, 1987):
  `c(j, k) = r_j`, the radius of a circle (sphere, hypersphere) centered
  at `x_j`.

- **Model I**: the most general model, in which `c` depends on the pair
  of objects and the dimension. It is not implemented yet, because how
  it is restricted to be estimable is given only in Okada and Imaizumi
  (1984).

The signs `alpha` and `beta` give the submodels. The defaults
`alpha = -1, beta = 1` are the models used by Okada and Imaizumi (1987,
2024; Model IIIB): `m_jk` is the distance from the near edge of the
ellipse (circle) of `j` to the far edge of that of `k`, and it is
assumed to be monotonically decreasing in similarities (increasing in
dissimilarities). With `alpha = beta` the model is symmetric (Saito's
cluster model for Model III).

The fit is `S = sqrt(sum((m_jk - mhat_jk)^2) / sum((m_jk - mbar)^2))`
with disparities `mhat_jk` from monotone regression (primary approach to
ties), as in
[`circle_mds()`](https://kosugitti.github.io/mdspace/reference/circle_mds.md).
`S` is minimized over the coordinates and the semiaxes (radii) by
[`stats::optim()`](https://rdrr.io/r/stats/optim.html) (BFGS) with the
analytic gradient, with the semiaxes parameterized by their logarithms
so that they stay positive.

**Initial values** (Okada & Imaizumi, 2024, Section 5.4). The rational
start is the Torgerson configuration of the symmetrized dissimilarities
(with the additive constant) in `ndim_max` dimensions and the radii from
the skew-symmetric part, which serve as initial semiaxes for all
dimensions (a small constant is added so that no semiaxis is zero). The
solution is reduced one dimension at a time down to `ndim`: for Model II
by dropping the dimension with the smallest
`sum_j x_jt^2 + sum_j r_jt^2` (Equation 5.23), because the orientation
of the dimensions is determined by the model; for Model III by
projecting on the principal axes. `nrandom` random starts are also tried
in `ndim` dimensions, because the models have local minima; a typical
one for Model II puts most points near the origin with almost concentric
ellipses (Section 5.3).

**Normalization.** The configuration is centered with the sum of squared
coordinates equal to `n`, and the semiaxes (radii) are rescaled with it.
For Model III the smallest radius is set to 0 (the radii are determined
only up to an additive constant); the semiaxes of Model II have no such
indeterminacy and are not shifted. The dimensions of Model II are not
rotated, because rotating them changes `S`; they are ordered by
decreasing variance of the coordinates. Model III is rotated to
principal axes.

**Orientation of dimensions (Model II).** If the semiaxes of two
dimensions are nearly equal for all objects (the points of the semiaxis
configuration lie near the diagonal), the ellipses are nearly circles
and the orientation of the two dimensions is not determined. As in Okada
and Imaizumi (2024, Sections 3.1.2 and 4.1.2), for each pair of
dimensions the semiaxis configuration is rotated by `angles` degrees
counterclockwise and clockwise, the coordinates are refitted with the
rotated semiaxes fixed (starting from the solution), and the resulting
`S` is reported. If it is clearly larger than the `S` of the solution,
the orientation is uniquely determined up to permutations and
reflections of the dimensions.

## References

Okada, A. (1990). A generalization of asymmetric multidimensional
scaling. In M. Schader & W. Gaul (Eds.), *Knowledge, data and
computer-assisted decisions* (pp. 127-138). Springer.

Okada, A., & Imaizumi, T. (1987). Nonmetric multidimensional scaling of
asymmetric proximities. *Behaviormetrika*, 14(21), 81-96.

Okada, A., & Imaizumi, T. (2024). *Applied multidimensional scaling of
asymmetric relationships*. Springer.

千野直仁・岡太彬訓 (1996). 非対称多次元尺度構成法とその周辺. 行動計量学,
23(2), 130-152.

## See also

[`circle_mds()`](https://kosugitti.github.io/mdspace/reference/circle_mds.md)
for Model IIIB with the two-mode three-way extension.

## Examples

``` r
x <- matrix(c(NA, 5, 2, 1, 3,
              3, NA, 4, 2, 2,
              1, 2, NA, 6, 3,
              2, 1, 4, NA, 5,
              4, 2, 2, 3, NA), 5, byrow = TRUE)
fit <- ellipse_mds(x, ndim = 2, nrandom = 2)
fit
#> mdspace: ellipse 
#> Call: ellipse_mds(x = x, ndim = 2, nrandom = 2)
#> Objects: 5  Dimensions: 2 
#> Stress: 0.0000 
#> Iterations: 8 
#> Model II (ellipse-distance), alpha = -1, beta = 1
#> Starts: 3 (rational + 2 random); best = rational, S range 0.0000 to 0.0000
#> S after rotating the semiaxis configuration (solution S = 0.0000):
#>  dims angle counterclockwise clockwise
#>   1-2    10           0.0046    0.0110
#>   1-2    15           0.0526    0.0484
#>   1-2    20           0.0612    0.0576
```
