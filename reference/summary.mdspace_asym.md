# Summarize a symmetric / skew-symmetric decomposition

The sums of squares with their shares, and for each object the row sum
of the skew-symmetric part `K` (positive when the object's row exceeds
its column on the whole) with the row sum of `|K|`.

## Usage

``` r
# S3 method for class 'mdspace_asym'
summary(object, ...)
```

## Arguments

- object:

  An object of class `mdspace_asym`.

- ...:

  Not used.

## Value

An object of class `summary.mdspace_asym` with elements `ss` (sums of
squares and shares), `asymmetry`, and `objects` (row sums of `K` and of
`|K|`).

## References

Gower, J. C. (1977). The analysis of asymmetry and orthogonality. In J.
R. Barra et al. (Eds.), *Recent developments in statistics* (pp.
109-123). North-Holland.
