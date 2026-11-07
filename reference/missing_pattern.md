# Missing-value patterns of asymmetric proximity matrices

Classifies the missing cells of a square matrix (or each matrix of a
list) into four types:

## Usage

``` r
missing_pattern(x)
```

## Arguments

- x:

  A square matrix, or a list of square matrices.

## Value

A data frame with one row per matrix: `n`, `one_side`, `both_sides`,
`diagonal`, `absent` (number of absent objects), `absent_objects` (their
labels), and `rate` (share of missing cells).

## Details

- `one_side`: off-diagonal cell `(j, k)` is missing but `(k, j)` is
  observed (this includes an object that did not respond, whose row is
  missing while its column is observed);

- `both_sides`: both `(j, k)` and `(k, j)` are missing (counted once per
  cell, so a pair counts 2);

- `diagonal`: missing diagonal cells;

- `absent`: objects whose whole row and column (off the diagonal) are
  missing, i.e. objects absent at that occasion. Their cells are not
  counted in the other types.

## Examples

``` r
x <- matrix(c(NA, 1, 2, NA, NA, 3, 4, 5, 6), 3, byrow = TRUE)
missing_pattern(x)
#>   matrix n one_side both_sides diagonal absent absent_objects      rate
#> 1     T1 3        1          0        2      0                0.3333333
```
