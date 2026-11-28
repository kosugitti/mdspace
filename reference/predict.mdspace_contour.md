# Valences of a contour map model at given points

Valences of a contour map model at given points

## Usage

``` r
# S3 method for class 'mdspace_contour'
predict(object, newdata = NULL, ...)
```

## Arguments

- object:

  An object of class `mdspace_contour`.

- newdata:

  A matrix of points with as many columns as the configuration. `NULL`
  returns the heights at the objects.

- ...:

  Not used.

## Value

A matrix with one row per point and one column per person.
