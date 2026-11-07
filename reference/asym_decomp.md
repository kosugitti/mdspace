# Decompose a square matrix into symmetric and skew-symmetric parts

`X = S + K` with `S = (X + t(X)) / 2` and `K = (X - t(X)) / 2`. The
share of asymmetry is `sum(K^2) / sum(X^2)` (Gower, 1977). Missing
values are ignored in the sums.

## Usage

``` r
asym_decomp(x)
```

## Arguments

- x:

  A square proximity matrix (see
  [`as_proximity()`](https://kosugitti.github.io/mdspace/reference/as_proximity.md)).

## Value

An object of class `mdspace_asym` with elements `S`, `K`, `ss` (sums of
squares of `X`, `S`, `K`) and `asymmetry`.

## References

Gower, J. C. (1977). The analysis of asymmetry and orthogonality. In J.
R. Barra et al. (Eds.), *Recent developments in statistics* (pp.
109-123). North-Holland.

## Examples

``` r
x <- matrix(c(0, 1, 1, 7, 1, 0, 1, 7, 7, 7, 0, 1, 1, 1, 7, 0), 4, byrow = TRUE)
asym_decomp(x)
#> Symmetric / skew-symmetric decomposition
#> Sum of squares: total = 252  symmetric = 162  skew = 90 
#> Share of asymmetry: 0.3571 
```
