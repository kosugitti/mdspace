# Plot a radius-distance (circle) model solution

Each object is drawn as a point with a circle of its radius. A larger
circle means that the object is reached more easily from the others than
it reaches them (for similarities).

## Usage

``` r
# S3 method for class 'mdspace_circle'
plot(
  x,
  dims = c(1, 2),
  labels = TRUE,
  source = NULL,
  which = c("conf", "weights"),
  ...
)
```

## Arguments

- x:

  An object of class `mdspace_circle`.

- dims:

  Two dimensions to plot.

- labels:

  Logical; draw object (or source) labels.

- source:

  Two-mode three-way model only: `NULL` (common circles), the index or
  name of a source, or `"all"`. Used with `which = "conf"`.

- which:

  `"conf"` (default) or `"weights"` (two-mode three-way model only).

- ...:

  Passed to
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

`x`, invisibly.

## Details

For the two-mode three-way model, `source` chooses what is drawn around
the points: the common circles of radius `r_j` (`NULL`), the ellipses of
one source with semiaxes `u_it r_j` (an index or a name), or the
ellipses of all sources overlaid in different colours (`"all"`). The
axis limits are the same for every choice of `source`, so that plots of
different sources can be compared side by side (set `par(mfrow = )` and
call the function once per source). In one dimension the ellipse of a
source is the circle of radius `u_i1 r_j`.

`which = "weights"` (two-mode three-way model, two or more dimensions)
draws the sources as points in the plane of their asymmetry weights
`u_it`, with the diagonal. A source on the diagonal keeps the circles;
the farther from the diagonal, the more elongated its ellipses are along
the dimension with the larger weight. This is the counterpart of the
weight space of INDSCAL and of `which = "semiaxis"` in
[`plot.mdspace_ellipse()`](https://kosugitti.github.io/mdspace/reference/plot.mdspace_ellipse.md).

## References

Okada, A., & Imaizumi, T. (1987). Nonmetric multidimensional scaling of
asymmetric proximities. *Behaviormetrika*, 14(21), 81-96.

Okada, A., & Imaizumi, T. (1997). Asymmetric multidimensional scaling of
two-mode three-way proximities. *Journal of Classification*, 14,
195-224.
