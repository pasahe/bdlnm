# Calculate the exposure value that minimizes or maximizes the overall cumulative effect of a Bayesian distributed lag non-linear model (B-DLNM)

Find exposure values that optimize the overall effect for each posterior
sample drawn from a Bayesian distributed lag non-linear model
([`bdlnm()`](https://pasahe.github.io/bdlnm/reference/bdlnm.md)). The
function returns the exposure value that minimizes or maximizes the
overall cumulative effect (summed across lags) for each posterior
sample, together with summary statistics (mean, sd, credible-interval
quantiles and mode). When used to find the minimum effect in
temperature–mortality analyses this optimal exposure value is commonly
called the Minimum Mortality Temperature (MMT).

## Usage

``` r
optimal_exposure(
  object,
  basis = NULL,
  exp_at = NULL,
  lag_at = NULL,
  which = "min",
  local_optimal = FALSE,
  ci.level = 0.95
)
```

## Arguments

- object:

  A fitted `"bdlnm"` object returned by
  [`bdlnm()`](https://pasahe.github.io/bdlnm/reference/bdlnm.md).

- basis:

  If the `bdlnm` model has more than one basis, the name of the basis to
  use to compute predictions. It must be one of `names(object$basis)`.
  If the model contains only one basis it is selected automatically.

- exp_at:

  Numeric vector of exposure values at which to evaluate predictions. If
  `NULL`, the exposure range is extracted from the attributes of the
  specified `basis` and a grid of 50 values is constructed using
  [`pretty()`](https://rdrr.io/r/base/pretty.html).

- lag_at:

  Numeric vector of integer lag values ver which to compute the overall
  cumulative effect. Only used when `basis` is a `crossbasis`. If
  `NULL`, the overall effect is computed by summing over the full lag
  range stored in `basis` (with step size `1`).

- which:

  Selection criterion to calculate the optimal exposure: `"min"`
  (default) chooses the exposure with the minimum overall cumulative
  effect, `"max"` chooses the exposure with maximum overall cumulative
  effect.

- local_optimal:

  Logical (default `FALSE`). When `TRUE` find a local optimal (minimum
  or maximum) point with the optimal effect instead of the absolute
  optimal point. If a local optimal point doesn't exist it will fall
  back to finding the absolute optimal point.

- ci.level:

  Numeric in `(0,1)` giving the credible-interval level (default
  `0.95`). Credible interval quantiles are computed from the posterior
  samples.

## Value

An S3 object of class `"optimal_exposure"` containing:

- `est`: numeric vector with the optimal exposure value for each
  posterior sample (named sample1, sample2, ...).

- `summary`: a one-row data frame with summary statistics for the
  optimal values across all samples (mean, sd, quantiles, mode).

## Details

The function internally calls
[`bcrosspred()`](https://pasahe.github.io/bdlnm/reference/bcrosspred.md)
to compute the posterior distribution of the overall cumulative exposure
effect for the grid specified by `exp_at`. For each posterior sample the
function calculates the exposure value that optimizes (minimizes or
maximizes) the overall cumulative effect and then summarizes these
optimal values across samples using mean, sd, credible-interval
quantiles and the mode (most frequent observed value).

The overall cumulative effect is computed by summing for each exposure
the lag-specific effects over the lags specified in `lag_at`. If
`lag_at` is `NULL`, the cumulative effect is computed for each exposure
by summing over the full lag range stored in `basis` (with step size 1).
If `basis` is a `onebasis`, the function optimizes the exposure-response
association stored in `$matfit`, and `lag_at` is ignored.

The function searches for the absolute optimal value (minimum or
maximum) of each sample in the posterior distribution, by default. If
`local_optimal` is set to `TRUE`, the function searches for a local
optimal point instead. If more than one optimal point is found, the
function will return the one with the optimal effect. If a posterior
sample has no local optimal values, the function returns the absolute
optimal value.

This optimal exposure value can be used as the reference exposure value
to estimate effects passing it to the
[`bcrosspred()`](https://pasahe.github.io/bdlnm/reference/bcrosspred.md)
and
[`attributable()`](https://pasahe.github.io/bdlnm/reference/attributable.md)
functions as the center exposure. In temperature-mortality studies, for
example, the minimum exposure value is typically used as the optimal
exposure value to center the effects and it's called Minimum Mortality
Temperature (MMT). However, note that in the Bayesian framework, this
reference temperature is characterized by a full posterior distribution
(in contrast to the frequentist approach, where the association is
centered on a single point estimate). This distribution may be
asymmetric and non-unimodal, so reporting a single summary statistic
(e.g., the median) as the reference value can be misleading in such
cases. Therefore, before selecting an optimal exposure value as the
center, it is recommended that you visualize the distribution of the
optimal exposure values using
[`plot.optimal_exposure()`](https://pasahe.github.io/bdlnm/reference/plot.optimal_exposure.md).

This function cannot be used when the specified basis function is one of
`thr`, `strata`, `integer`, or `lin`. The exposure-response relationship
is discrete, piecewise, or strictly linear in these situations, so
searching for an optimum is not meaningful.

## References

Quijal-Zamorano M, Martinez-Beneito MA, Ballester J, Marí-Dell’Olmo M.
Spatial Bayesian distributed lag non-linear models (SB-DLNM) for
small-area exposure-lag-response epidemiological modelling.
International Journal of Epidemiology. 2024;53(3):dyae061.

Gasparrini A (2011). Distributed lag linear and non-linear models in R:
the package dlnm. Journal of Statistical Software, 43(8), 1–20.

Armstrong B. Models for the relationship between ambient temperature and
daily mortality. Epidemiology. 2006;17(6):624-31.

## See also

[`plot.optimal_exposure()`](https://pasahe.github.io/bdlnm/reference/plot.optimal_exposure.md)
to plot the optimal exposure values stored in a `"optimal_exposure"`
object.

[`bcrosspred()`](https://pasahe.github.io/bdlnm/reference/bcrosspred.md)
to predict exposure–lag–response associations for a `"bdlnm"` object,

[`bdlnm()`](https://pasahe.github.io/bdlnm/reference/bdlnm.md) to fit a
Bayesian distributed lag non-linear model (`"bdlnm"`).

[`attributable()`](https://pasahe.github.io/bdlnm/reference/attributable.md)
to calculate attributable fractions and numbers for a `"bdlnm"` object.

## Author

Pau Satorra, Marcos Quijal-Zamorano.

## Examples

``` r
# Set exposure-response and lag-response spline parameters
 dlnm_var <- list(
   var_prc = c(10, 75, 90),
   var_fun = "ns",
   lag_fun = "ns",
   max_lag = 21,
   lagnk = 3
 )


# Set cross-basis parameters
 argvar <- list(fun = dlnm_var$var_fun,
                knots = stats::quantile(london$tmean,
                                 dlnm_var$var_prc/100, na.rm = TRUE),
                Bound = range(london$tmean, na.rm = TRUE))

 arglag <- list(fun = dlnm_var$lag_fun,
                knots = dlnm::logknots(dlnm_var$max_lag, nk = dlnm_var$lagnk))

 # Create crossbasis
 cb <- dlnm::crossbasis(london$tmean, lag = dlnm_var$max_lag, argvar, arglag)

 # Seasonality of mortality time series
 seas <- splines::ns(london$date, df = round(8 * length(london$date) / 365.25))

 # Prediction values (equidistant points)
 temp <- round(seq(min(london$tmean), max(london$tmean), by = 0.1), 1)

if (bdlnm:::check_inla()) {
 # Fit the model
 mod <- bdlnm(mort_75plus ~ cb + factor(dow) + seas, data = london, family = "poisson",
             sample.arg = list(seed = 432, seed = 1L))

 # Find minimum risk exposure value
 mmt <- optimal_exposure(mod, "cb", exp_at = temp)
}
#> Warning: Since 'seed!=0', parallel model is disabled and serial model is selected, num.threads='1:1'
```
