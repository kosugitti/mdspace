# Summarize a GSTATIS solution

The configuration of the first compromise with the eigenvalues of the RV
matrix (inter-structure), the weights of the sources in each compromise,
the eigenvalues and goodness of fit of the first compromise, and the
distance of each source from each compromise. `stress` is `NA` because
GSTATIS is not a distance model.

## Usage

``` r
# S3 method for class 'mdspace_gstatis'
summary(object, ...)
```

## Arguments

- object:

  An object of class `mdspace_gstatis`.

- ...:

  Not used.

## Value

An object of class `summary.mdspace`; the model-specific tables are in
`extra`.

## References

Chino, N., Grorud, A., & Yoshino, R. (1996). A complex analysis for
two-mode three-way asymmetric relational data. *Abstracts of the talks
at the 5th Conference of the International Federation of Classification
Societies*, Vol. 2, 83-86.
