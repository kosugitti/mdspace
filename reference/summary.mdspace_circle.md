# Summarize a radius-distance (circle) model solution

The configuration with the radii. For the two-mode three-way model the
weights of the sources (`w`: symmetry, `u`: asymmetry), the aspect ratio
of their ellipses and `S` by source are added.

## Usage

``` r
# S3 method for class 'mdspace_circle'
summary(object, ...)
```

## Arguments

- object:

  An object of class `mdspace_circle`.

- ...:

  Not used.

## Value

An object of class `summary.mdspace`; for the two-mode three-way model
the tables of the sources are in `extra`.

## References

Okada, A., & Imaizumi, T. (1987). Nonmetric multidimensional scaling of
asymmetric proximities. *Behaviormetrika*, 14(21), 81-96.

Okada, A., & Imaizumi, T. (1997). Asymmetric multidimensional scaling of
two-mode three-way proximities. *Journal of Classification*, 14,
195-224.
