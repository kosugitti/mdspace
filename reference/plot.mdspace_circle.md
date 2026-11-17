# Plot a radius-distance (circle) model solution

Each object is drawn as a point with a circle of its radius. A larger
circle means that the object is reached more easily from the others than
it reaches them (for similarities).

## Usage

``` r
# S3 method for class 'mdspace_circle'
plot(x, dims = c(1, 2), labels = TRUE, ...)
```

## Arguments

- x:

  An object of class `mdspace_circle`.

- dims:

  Two dimensions to plot.

- labels:

  Logical; draw object labels.

- ...:

  Passed to
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

`x`, invisibly.
