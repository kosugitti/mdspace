# Plot one complex plane of an HFM solution

Each object is drawn as an arrow from the origin. For a positive
eigenvalue, the counterclockwise angle from object k to object j is
proportional to the skew-symmetric part `K[j, k]`.

## Usage

``` r
# S3 method for class 'mdspace_hfm'
plot(x, dim = 1, ...)
```

## Arguments

- x:

  An object of class `mdspace_hfm`.

- dim:

  Which complex dimension to plot.

- ...:

  Passed to
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

`x`, invisibly.
