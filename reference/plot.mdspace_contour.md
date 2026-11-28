# Plot a contour map model

Draws the valences of one person over a plane of the configuration, as a
contour map (`type = "contour"`) with the contour of valence 0 drawn
thick, or as a valence surface (`type = "persp"`) with each object drawn
at its rating, filled when the rating is above the surface and open when
it is below (as in Kosugi & Fujihara, 2004). For three or more
dimensions the plane passes through the centroid of the configuration.

## Usage

``` r
# S3 method for class 'mdspace_contour'
plot(
  x,
  person = 1,
  type = c("contour", "persp"),
  dims = c(1, 2),
  grid = NULL,
  lim = NULL,
  nlevels = 10,
  labels = TRUE,
  theta = 30,
  phi = 30,
  ...
)
```

## Arguments

- x:

  An object of class `mdspace_contour`.

- person:

  Row number or name of the person to draw.

- type:

  `"contour"` or `"persp"`.

- dims:

  Two dimensions spanning the plane.

- grid:

  Number of grid points on each axis. `NULL` uses 101 for
  `type = "contour"` and 41 for `type = "persp"`.

- lim:

  Range of both axes. `NULL` extends the range of the configuration by
  30%.

- nlevels:

  Approximate number of contour levels.

- labels:

  Logical; draw object labels.

- theta, phi:

  Viewing angles for `type = "persp"`.

- ...:

  Passed to
  [`graphics::contour()`](https://rdrr.io/r/graphics/contour.html) or
  [`graphics::persp()`](https://rdrr.io/r/graphics/persp.html).

## Value

The grid (`x`, `y`, `z`), invisibly.
