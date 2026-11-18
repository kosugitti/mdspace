# Plot a von Mises model solution

Each object is drawn as a point with a small circle, and the curve
outside the circle shows `pi_i(theta)`, the strength of its one-sided
relations in each direction (as in the figures of Shojima, 2011, 2012).
A line from the point shows the mean direction `mu_i`.

## Usage

``` r
# S3 method for class 'mdspace_vonmises'
plot(x, labels = TRUE, rho = NULL, ...)
```

## Arguments

- x:

  An object of class `mdspace_vonmises`.

- labels:

  Logical; draw object labels.

- rho:

  Radius of the base circles. `NULL` uses 6% of the range of the
  configuration.

- ...:

  Passed to
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

`x`, invisibly.
