# Check if INLA is installed and meets the minimum version requirement

Checks whether the INLA package is installed and whether the installed
version is \>= 23.4.24. This is an internal utility used by functions
that depend on INLA.

## Usage

``` r
check_inla(error = FALSE)
```

## Arguments

- error:

  Logical. If `TRUE`, throws an error when INLA is not installed or the
  version is too old. If `FALSE` (default), returns `FALSE` silently.

## Value

`TRUE` if INLA is installed and meets the version requirement, `FALSE`
otherwise (when `error = FALSE`).
