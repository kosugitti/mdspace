# Asymmetric MDS by the von Mises model (AMISESCAL)

Each object `i` is a point `x_i` in two dimensions with a von Mises
distribution of directions, with mean direction `mu_i` and concentration
`kappa_i`, which expresses toward which direction the object's one-sided
relations point (Shojima, 2011, 2012). With `theta_ij` the direction
from `x_i` to `x_j` and `delta_ij` their Euclidean distance, the model
for the dissimilarity from `i` to `j` is `xi_ij = (1 - pi_ij) delta_ij`,
where `pi_ij = c_i exp{kappa_i (cos(theta_ij - mu_i) - 1)}`. `pi_ij` is
the von Mises density
`omega_i exp{kappa_i cos(theta - mu_i)} / I0(kappa_i)` of Shojima (2012)
rewritten with its peak `c_i = omega_i exp(kappa_i) / I0(kappa_i)`;
`c_i < 1` keeps `pi_ij < 1`. An object `j` lying in the direction `mu_i`
looks closer from `i` than `i` looks from `j`.

## Usage

``` r
vonmises_mds(
  x,
  type = c("dissimilarity", "similarity"),
  omega = c("free", "fixed"),
  prior = list(),
  nrandom = 10,
  maxit = 2000,
  reltol = 1e-10,
  hessian = FALSE,
  seed = 1
)
```

## Arguments

- x:

  A square matrix of asymmetric proximities. Rows evaluate columns
  (`x[i, j]` is from `i` to `j`). The diagonal is ignored; other missing
  cells are left out of the fit.

- type:

  `"dissimilarity"` (default) or `"similarity"`. Similarities are
  reversed as `max + min - x`. If `x` has a `type` attribute set by
  [`as_proximity()`](https://kosugitti.github.io/mdspace/reference/as_proximity.md),
  that is used.

- omega:

  `"free"` (default; Shojima, 2012) or `"fixed"` (`1 / (2 pi)`; Shojima,
  2011). See Details.

- prior:

  A list with `x_sd`, `kappa_shape`, `kappa_rate`, `c_a` and `c_b`. See
  Details. Missing elements take the defaults
  `list(x_sd = 1, kappa_shape = 2, kappa_rate = 0.5, c_a = 2, c_b = 4)`.

- nrandom:

  Number of random starts in addition to the rational start.

- maxit, reltol:

  Passed to [`stats::optim()`](https://rdrr.io/r/stats/optim.html).

- hessian:

  Logical; compute approximate standard errors.

- seed:

  Random seed for the random starts. The state of the random number
  generator is restored afterwards. `NULL` uses (and advances) the
  current state.

## Value

An object of class `c("mdspace_vonmises", "mdspace")` with `conf` (n x
2, in units of the scaled dissimilarities), `mu` (radians), `kappa`, `c`
(peaks of `pi`), `sigma`, `stress` (`sqrt(sum((g - xi)^2) / sum(g^2))`
on the scaled dissimilarities), `logpost` (log posterior at the mode, up
to a constant), `niter`, `pi` and `xi` (n x n; `xi` in the original
units), `scale` (the root mean square used for scaling), `anchor`
(labels of the three anchors), `starts` (`logpost` of each start, the
rational one first), `best_start`, `converged`, and with
`hessian = TRUE`, `se` (a data frame of standard errors of `mu`, `kappa`
and `c`).

## Details

With `omega = "fixed"` the density is the ordinary von Mises density
(`omega_i = 1 / (2 pi)`, Shojima, 2011), so that `c_i` is determined by
`kappa_i`, and `kappa_i < 6.5` keeps `c_i < 1`. With `omega = "free"`
(Shojima, 2012) `c_i` is a parameter of its own. Leaving it free without
further constraint lets the solution degenerate (points move apart while
all `pi_ij` approach 1), which the original program prevents by an
additional term in the stress function. Here the estimates are instead
the mode of a posterior distribution (maximum a posteriori, MAP) with
weak priors, the same as the Bayesian version written in Stan
(`develop/vonmises_mds.stan` in the source repository):
`g_ij ~ N(xi_ij, sigma^2)`, deviations of the coordinates from their
centroid `~ N(0, prior$x_sd^2)`,
`kappa_i ~ Gamma(prior$kappa_shape, prior$kappa_rate)` truncated at the
upper bound, `c_i ~ Beta(prior$c_a, prior$c_b)` and
`sigma ~ half-t(4, 0, 1)`. The dissimilarities are divided by their root
mean square before fitting, so that the priors are on a common scale.

The configuration is identified as in the Bayesian version by a Cholesky
factor (lower triangular with positive diagonal; Stan Development Team,
Reference Manual, "Cholesky factors of covariance matrices"): one anchor
object is placed at the origin, the second on the positive side of the
first axis and the third on the positive side of the second axis. The
anchors are chosen from the Torgerson solution of the symmetric part so
that they are far apart and not on a line. This anchoring is only a
device for computation: the likelihood and the prior (placed on the
deviations from the centroid, which is isotropic) do not change under
translation, rotation and reflection, so the solution does not depend on
which objects are the anchors. The mean directions rotate with the
configuration, so they are identified too. For presentation the returned
configuration is centered and rotated to its principal axes, and the
mean directions are rotated with it.

The posterior mode has many local optima, so the rational start (the
Torgerson solution with each mean direction pointing at the objects
toward which the object has one-sided relations) and `nrandom` random
starts are tried and the best is returned. The minimization uses the
analytic gradient and
[`stats::optim()`](https://rdrr.io/r/stats/optim.html) (BFGS). With
`hessian = TRUE` approximate standard errors of `mu`, `kappa` and `c`
are obtained from the Hessian of the negative log posterior (Laplace
approximation; delta method for the transformed parameters).

## References

Shojima, K. (2011). Asymmetric von Mises scaling. *Proceedings of the
39th Annual Meeting of the Behaviormetric Society of Japan*, 261-262.

Shojima, K. (2012). On the stress function of asymmetric von Mises
scaling. Paper presented at the Joint Meeting of the German
Classification Society and the Japanese Classification Society.

Mardia, K. V., & Jupp, P. E. (2000). *Directional statistics*. Wiley.

Stan Development Team. *Stan reference manual*: Cholesky factors of
covariance matrices. <https://mc-stan.org/docs/reference-manual/>

## Examples

``` r
x <- matrix(c(NA, 5, 4, 5, 3,
              6, NA, 4, 2, 1,
              4, 4, NA, 3, 3,
              4, 1, 2, NA, 1,
              7, 1, 2, 1, NA), 5, byrow = TRUE,
            dimnames = list(LETTERS[1:5], LETTERS[1:5]))
fit <- vonmises_mds(x, nrandom = 2)
fit
#> mdspace: vonmises 
#> Call: vonmises_mds(x = x, nrandom = 2)
#> Objects: 5  Dimensions: 2 
#> Stress: 0.0000 
#> Iterations: 611 
#> omega: free (free; Shojima, 2012) 
#> sigma: 0.0000  log posterior: 371.87 
#> Starts: 3 (rational + 2 random); best = rational
```
