# Plot a GSTATIS solution

`type = "interstructure"` plots the sources on the first two
eigenvectors of the GRV matrix. `type = "compromise"` plots one complex
plane of a compromise like
[`plot.mdspace_hfm()`](https://kosugitti.github.io/mdspace/reference/plot.mdspace_hfm.md).
`type = "trajectory"` draws, for each object, its images over the
sources joined by lines (step 3); the filled point is the image of the
compromise row, i.e. the weighted mean of the images with the GRV
weights.

## Usage

``` r
# S3 method for class 'mdspace_gstatis'
plot(
  x,
  type = c("trajectory", "compromise", "interstructure"),
  dim = 1,
  view = 1,
  objects = NULL,
  scale = c("none", "max"),
  ...
)
```

## Arguments

- x:

  An object of class `mdspace_gstatis`.

- type:

  `"trajectory"`, `"compromise"`, or `"interstructure"`.

- dim:

  Which complex dimension to plot.

- view:

  Which compromise (point of view).

- objects:

  Objects to draw in `type = "trajectory"` (names or indices); all by
  default.

- scale:

  `"none"` or `"max"`: with `"max"` each object's trajectory is divided
  by its largest norm over sources (over all kept dimensions), as
  printed by Chino's program.

- ...:

  Passed to
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

`x`, invisibly.
