# Hermitian form model for asymmetric proximities

The square matrix `X` is split into its symmetric part `S` and its
skew-symmetric part `K`, and the Hermitian matrix `H = S + iK` is
eigen-decomposed: `H = U diag(lambda) U*`. The eigenvalues are real and
are ordered by absolute value. Each dimension is a complex plane; object
`j` is the point `z_j = u_j * sqrt(|lambda|)`. For a dimension with a
positive eigenvalue, `Re(z_j Conj(z_k))` reproduces the symmetric part
and `Im(z_j Conj(z_k))` the skew-symmetric part, so the angle from `k`
to `j` (counterclockwise) shows the direction of asymmetry.

## Usage

``` r
hfm_mds(
  x,
  ndim = 2,
  rotate = c("std", "varimax", "none"),
  missing = c("none", "impute"),
  impute_ndim = ndim,
  tol = 1e-04,
  maxit = 1000
)
```

## Arguments

- x:

  A square similarity matrix (rows = from, columns = to). The diagonal
  is used as is.

- ndim:

  Number of complex dimensions to keep.

- rotate:

  `"std"`, `"varimax"`, or `"none"`.

- missing:

  `"none"` (default; missing values are an error) or `"impute"`. With
  `"impute"`, an off-diagonal cell whose transposed cell is observed is
  filled with that value (so the pair has no skew-symmetric part), and
  the diagonal cells and pairs missing on both sides are estimated by
  the MINERS method: starting from the mean of the observed cells in the
  same row and column, they are replaced by the rank-`p` HFM
  reconstruction `Re(H) + Im(H)` until the largest change is below
  `tol`. This is the refactoring method of Chino and Harshman (1996), a
  generalization of Thomson's refactoring method, and the MINERS option
  of the 2001 Visual Basic program. An object whose whole row and column
  are missing is an error. See
  [`missing_pattern()`](https://kosugitti.github.io/mdspace/reference/missing_pattern.md).

- impute_ndim:

  Number of HFM dimensions `p` used in the MINERS reconstruction;
  defaults to `ndim`.

- tol, maxit:

  Convergence criterion and maximum number of iterations of the MINERS
  method.

## Value

An object of class `c("mdspace_hfm", "mdspace")` with `coord` (n x ndim
complex matrix), `conf` (n x 2\*ndim real matrix: real and imaginary
parts of each dimension), `eigen` (all eigenvalues), `gof`
(`lambda^2 / sum(lambda^2)` for each dimension), `unexplained` (share of
the sum of squares of `H` not reproduced by the kept dimensions, i.e.
1 - cumulative GOF), `hermitian` (`H`), and, with `missing = "impute"`,
`imputed` (data frame of the filled cells with their type and value) and
`impute_niter`. HFM is not a distance model, so `stress` is `NA`.

## Details

An eigenvector is determined only up to a factor `exp(i * theta)`, i.e.
a rotation of its complex plane. `rotate = "std"` (default) turns each
plane so that the principal axis of the points lies on the real axis
(the sum of squared real parts is maximal). `rotate = "varimax"` turns
each plane so that the variance of `Re(z^2)` over objects is maximal,
which pushes the points toward the real and imaginary axes. Both follow
the 2001 Visual Basic program (`rotateSTD` and `rotateMAX`). In both
cases the plane is then flipped, if needed, so that the object with the
largest absolute real part lies on the positive side.

## References

Chino, N., & Harshman, R. A. (1996). A refactoring method for fitting
HFM that allows self-similarity estimation and missing values.
*Proceedings of the 24th Annual Meeting of the Behaviormetric Society of
Japan*, 206-207.

Chino, N., & Shiraiwa, K. (1993). Geometrical structures of some
non-distance models for asymmetric MDS. *Behaviormetrika*, 20(1), 35-47.

## Examples

``` r
x <- matrix(c(0, 1, 1, 7, 1, 0, 1, 7, 7, 7, 0, 1, 1, 1, 7, 0), 4, byrow = TRUE)
fit <- hfm_mds(x, ndim = 2)
fit
#> mdspace: hfm 
#> Call: hfm_mds(x = x, ndim = 2)
#> Objects: 4  Dimensions: 2 
#> Unexplained SS (1 - cumulative GOF): 0.0144 
#> 
#>                  D1      D2
#> eigenvalue -11.4501 10.8290
#> gof          0.5203  0.4653
#> cumulative   0.5203  0.9856
#> Note: negative eigenvalues reverse the reading of angles in those planes.
plot(fit, dim = 1)
```
