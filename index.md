# bdlnm

The goal of `bdlnm` is to provide tools for fitting **Bayesian
distributed lag non-linear models (B-DLNM)** using the **INLA**
framework. The package allows users to estimate exposure-lag-response
associations, generate posterior predictions, compute attributable risk
measures, and estimate optimal exposure values within a fully Bayesian
framework.

## Key features

The package provides functions to:

- Fit Bayesian distributed lag linear and non-linear models using `INLA`
- Generate predictions of exposure-lag-response associations
- Estimate attributable fractions and attributable numbers
- Estimate optimal exposure values (e.g., minimum risk exposure)
- Visualize exposure–lag–response associations

More detailed examples and tutorials are available in the package
vignettes.

## Installation

You can install the development version of bdlnm like so:

``` r
devtools::install_github("pasahe/bdlnm")
```

> At least the stable version of INLA 23.4.24 (or newest) must be
> installed beforehand. You can install the newest stable INLA version
> by:

``` r
install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
```

## Citation

If you use `bdlnm` in your research, please cite:

``` r
citation("bdlnm")
```

## Acknowledgements

This package extends the distributed lag modelling framework implemented
in the `dlnm` package developed by [Gasparrini
(2011)](https://doi.org/10.18637/jss.v043.i08) to a Bayesian setting.
