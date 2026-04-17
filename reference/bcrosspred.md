# Predict exposure-lag-response effects from a Bayesian distributed-lag model (B-DLNM).

Calculate predictions from a fitted Bayesian distributed lag non-linear
model ([`bdlnm()`](https://pasahe.github.io/bdlnm/reference/bdlnm.md)).
Predicted associations are computed on a grid of values of the exposure
and lags, relative to a reference exposure center value. The function
gives posterior samples of exposure–lag-specific associations, overall
cumulative associations (summed across lags) and optionally incremental
cumulative associations, together with summary statistics (mean, sd,
credible interval quantiles and mode).

## Usage

``` r
bcrosspred(
  object,
  basis = NULL,
  exp_at = NULL,
  lag_at = NULL,
  cen = NULL,
  model.link = NULL,
  ci.level = 0.95,
  cumul = FALSE
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
  `NULL`, the exposure range is extracted from the attributes of `basis`
  and a grid of 50 values is constructed using
  [`pretty()`](https://rdrr.io/r/base/pretty.html).

- lag_at:

  Numeric vector of integer lag values at which to evaluate predictions.
  Only used when `basis` is a `crossbasis`. If `NULL`, the full lag
  range stored in `basis` is used with a step of `1`.

- cen:

  Centering exposure value for predictions. If `NULL` the centering
  value depends on the exposure basis function or set to a mid-range
  value (see Details).

- model.link:

  Optional character specifying the model link (if `NULL` it is inferred
  from the fitted model).

- ci.level:

  Numeric in `(0,1)` giving the credible interval level (default
  `0.95`). Credible interval quantiles are computed from the posterior
  samples.

- cumul:

  Logical; if `TRUE` compute incremental cumulative effects along lags
  (default `FALSE`). It will give the cumulative effect from lag 0 up to
  each subsequent lag (e.g., lag `0`, lag `0–1`, lag `0–2`, etc.).

## Value

An object of class `"bcrosspred"` (a list) with elements including:

- `exp_at`: the exposure grid used for prediction as supplied via the
  `exp_at` argument (vector).

- `lag_at`: if basis is `crossbasis`, the lag grid used for prediction
  as supplied via the `lag_at` argument (vector).

- `cen`: the exposure centering value used for prediction (number).

- `coefficients`: matrix of posterior coefficient draws (columns =
  samples).

- `coefficients.summary`: matrix of coefficient summaries (mean, sd,
  quantiles, mode).

- `matfit`: 3D array of sampled lag-specific effects (exp × lag ×
  samples).

- `matfit.summary`: 3D array of summaries for `matfit` (exp × lag ×
  summary-statistics).

- `allfit`: matrix of sampled overall cumulative (summed across lags)
  effects (exp × samples).

- `allfit.summary`: matrix of summaries for `allfit` (exp x
  summary-statistics).

- `cumfit`: (optional) 3D array of sampled incremental cumulative
  effects (exp × lag × samples).

- `cumfit.summary`: (optional) 3D array of summaries for `cumfit` (exp ×
  lag × summary-statistics).

- `matRRfit`, `allRRfit`, `matRRfit.summary`, `allRRfit.summary`,
  `cumRRfit.summary` (optional), `cumRRfit.summary` (optional:
  relative-risk versions (only when link is `log` or `logit`).

- `ci.level`, `model.class`, `model.link`.

## Details

The function computes predictions for specific combinations of exposure
and lag values specified in `exp_at` and `lag_at.` All values must lie
within the range defined by the specified `basis.` Note that if the
specified `basis` is a `onebasis`, `lag_at` is ignored. If either
argument is `NULL`, the grid is derived from the attributes of the
specified `basis`: for exposure values, a grid of approximately 50
values is constructed using
[`pretty()`](https://rdrr.io/r/base/pretty.html) within the exposure
range; for lag values, an integer sequence covering the full lag range
with a step of 1 is used.

Predictions are computed relative to a centering/reference value
(`cen`). If `NULL`, the default `cen` depends on the exposure-response
basis function: for `strata`, `thr` and `integer` the reference
corresponds to the reference category, and for `lin` the reference is
set to 0. For other choices, such as `ns`, `bs`, `poly` or other
existing or user-defined functions, the default centering value is set
to an approximate mid-range value. For non-linear exposure-response
associations is sometimes recommended to manually set the centering
value to a data-driven center such as the optimal exposure value (see
[`optimal_exposure()`](https://pasahe.github.io/bdlnm/reference/optimal_exposure.md)).

Posterior sample of the predicted associations are stored as matrices
for the overall cumulative effect and 3D arrays for the
exposure-lag-specific predictions. Summaries across these samples are
computed using the mean, sd, credible-interval quantiles (the mid and
the lower/upper tails according to `ci.level`) and an approximate mode
obtained from a kernel density estimate. Relative risks versions of
these associations (exponentiated predictions) are also included if
`model.link` is equal to `"log"` or `"logit"`. The `model.link` can be
manually specified or, if `NULL`, it is tried to be inferred from the
`model` type in `object`.

Be aware of memory usage: exposure-lag-specific predictions are stored
in 3D arrays of dimension `length(predvar)` × `length(predlag)` ×
`n_sim`. For dense grids and many posterior samples this can be
computationally intensive.

This function can also be used to compute predictions for models with
simple uni-dimensional basis functions not including lags, if the
specified basis is `"onebasis"` instead of `"crossbasis"`. In this case,
only unlagged predicted associations are returned.

## Note

This function is inspired by
[`dlnm::crosspred()`](https://rdrr.io/pkg/dlnm/man/crosspred.html)
developed by Gasparrini (2011) <doi:10.18637/jss.v043.i08>. It has been
adapted to work in a Bayesian framework within the bdlnm package.

## References

Gasparrini A. (2011). Distributed lag linear and non-linear models in R:
the package dlnm. *Journal of Statistical Software*, 43(8), 1-20.
<doi:10.18637/jss.v043.i08>.

Quijal-Zamorano M., Martinez-Beneito M.A., Ballester J., Marí-Dell'Olmo
M. (2024). Spatial Bayesian distributed lag non-linear models (SB-DLNM)
for small-area exposure-lag-response epidemiological modelling.
*International Journal of Epidemiology*, 53(3), dyae061.
<doi:10.1093/ije/dyae061>.

## See also

[`plot.bcrosspred()`](https://pasahe.github.io/bdlnm/reference/plot.bcrosspred.md)
to plot the predicted associations stored in a `"bcrosspred"` object,

[`bdlnm()`](https://pasahe.github.io/bdlnm/reference/bdlnm.md) to fit a
Bayesian distributed lag non-linear model (`"bdlnm"`).

[`attributable()`](https://pasahe.github.io/bdlnm/reference/attributable.md)
to calculate attributable fractions and numbers for a `"bdlnm"` object,

[`optimal_exposure()`](https://pasahe.github.io/bdlnm/reference/optimal_exposure.md)
to estimate exposure values that optimize the predicted effect for a
`"bdlnm"` object.

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
 # Ensure it falls inside the range of temperatures after rounding:
 temp <- temp[temp >= min(london$tmean) & temp <= max(london$tmean)]

if (check_inla()) {
 # Fit the model
 mod <- bdlnm(mort_75plus ~ cb + factor(dow) + seas, data = london, family = "poisson",
             sample.arg = list(seed = 432, seed = 1L))

 # Prediction
 cpred <- bcrosspred(mod, exp_at = temp)
}
#> Warning: Since 'seed!=0', parallel model is disabled and serial model is selected, num.threads='1:1'
#> Warning: Centering value unspecified (`cen`). Automatically set to: 12.85.
```
