# Plot an ellipse-distance or radius-distance model solution

`which = "conf"` draws each object as a point with its ellipse (Model
II) or circle (Model III) in two dimensions. `which = "semiaxis"` draws
the semiaxis configuration of Model II with the diagonal; when the
points lie near the diagonal the ellipses are nearly circles and the
orientation of the two dimensions is not determined (Okada & Imaizumi,
2024, Fig. 3.4).

## Usage

``` r
# S3 method for class 'mdspace_ellipse'
plot(x, which = c("conf", "semiaxis"), dims = c(1, 2), labels = TRUE, ...)
```

## Arguments

- x:

  An object of class `mdspace_ellipse`.

- which:

  `"conf"` or `"semiaxis"`.

- dims:

  Two dimensions to plot.

- labels:

  Logical; draw object labels.

- ...:

  Passed to
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

`x`, invisibly.
