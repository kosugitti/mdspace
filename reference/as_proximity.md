# Coerce data to a square proximity matrix

Accepts a square matrix, a data frame, or a `dist` object and returns a
square numeric matrix with row and column labels. Rows are "from"
objects and columns are "to" objects, as in asymmetric MDS.

## Usage

``` r
as_proximity(x, type = c("dissimilarity", "similarity"), diag = NULL)
```

## Arguments

- x:

  A square matrix, data frame, or `dist` object.

- type:

  Either `"dissimilarity"` or `"similarity"`. Stored as an attribute;
  the values are not transformed.

- diag:

  Value to put on the diagonal. `NULL` keeps the diagonal as is.

## Value

A square numeric matrix with attribute `type`.

## Examples

``` r
d <- as_proximity(eurodist)
attr(d, "type")
#> [1] "dissimilarity"
```
