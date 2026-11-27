# Summarize an HFM solution

The configuration (real and imaginary parts of each complex dimension)
with the eigenvalues, the goodness of fit of the retained dimensions,
and the share of the sum of squares left unexplained. `stress` is `NA`
because the HFM is not a distance model.

## Usage

``` r
# S3 method for class 'mdspace_hfm'
summary(object, ...)
```

## Arguments

- object:

  An object of class `mdspace_hfm`.

- ...:

  Not used.

## Value

An object of class `summary.mdspace`; the model-specific tables are in
`extra`.

## References

Chino, N., & Shiraiwa, K. (1993). Geometrical structures of some
non-distance models for asymmetric MDS. *Behaviormetrika*, 20(1), 35-47.
