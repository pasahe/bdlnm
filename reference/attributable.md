# Calculate attributable number and fractions from a Bayesian distributed-lag model (B-DLNM).

Compute attributable numbers (AN) and attributable fractions (AF) from a
fitted Bayesian distributed lag non-linear model
([`bdlnm()`](https://pasahe.github.io/bdlnm/reference/bdlnm.md)). The
function uses posterior predicted relative risks from
[`bcrosspred()`](https://pasahe.github.io/bdlnm/reference/bcrosspred.md)
and applies a forward or backward lag algorithm to compute per-time and
total attributable measures, optionally filtering the calculation to a
subset of time points.

## Usage

``` r
attributable(
  object,
  data,
  name_date = NULL,
  name_exposure,
  name_cases = NULL,
  name_filter = NULL,
  dir = "back",
  basis = NULL,
  cen,
  range = NULL,
  lag_average = TRUE
)
```

## Arguments

- object:

  A fitted `"bdlnm"` class object returned by
  [`bdlnm()`](https://pasahe.github.io/bdlnm/reference/bdlnm.md).

- data:

  A data frame containing the temporal series needed to calculate
  attributable measures. It can include a date column (optional, see
  `name_date`), the exposure values (mandatory, see `name_exposure`),
  the number of cases (optional, see `name_cases`) and a binary (`0/1`)
  filter column (optional, see `name_filter`).

- name_date:

  Optional single string with the name of the column in `data`
  containing the date. The column must be of class `Date` or `POSIXt`.
  When provided the function checks that the series is regularly spaced
  (see `Details`).

- name_exposure:

  Single string with the name of the exposure column in `data`.

- name_cases:

  Optional single string with the name of the cases column in `data`. If
  not provided, the function returns only attributable fractions (AF)
  per time point.

- name_filter:

  Optional single string with the name of a binary (`0/1`) column in
  `data`. Only rows with value `1` are used to compute total AF and AN
  and per-time results are returned only for the filtered rows.

- dir:

  Character; direction of the algorithm to calculate attributable
  measures. `"back"` (default) calculates AF and AN attributing
  current-time outcome to past exposures; `"forw"` calculates AF and AN
  attributing current-time exposure to future outcomes.

- basis:

  If the `bdlnm` model has more than one basis, the name of the basis to
  use to compute predictions. It must be one of `names(object$basis)`.
  If the model contains only one basis it is selected automatically.

- cen:

  Numeric scalar; centering exposure value used to compute predictions.
  If missing the function will attempt to read it from the attributes of
  the specified `basis`. If no centering is available the function
  aborts.

- range:

  Optional numeric vector of length 2 with the exposure range for which
  attributable measures will be calculated. Values outside `range` are
  coerced to `cen` before prediction.

- lag_average:

  Logical (default `TRUE`). When `TRUE` use lag-averaged contributions
  to compute AN in the forward algorithm; when `FALSE` use the full
  lag-structured contributions instead.

## Value

A list with components:

- `af`: matrix (rows = time points, columns = posterior samples) with
  attributable fractions per time point.

- `an`: matrix (rows = time points, columns = posterior samples) with
  attributable numbers per time point.

- `aftotal`: numeric vector (length = number of posterior samples) with
  posterior samples of total AF across all the selected period.

- `antotal`: numeric vector (length = number of posterior samples) with
  posterior samples of total AN across all the selected period.

- `af.summary`: data.frame with summary statistics (mean, sd, quantiles,
  mode) for AF per time point.

- `an.summary`: data.frame with summary statistics (mean, sd, quantiles,
  mode) for AN per time point.

- `aftotal.summary`: data.frame with summary statistics (mean, sd,
  quantiles, mode) for total AF.

- `antotal.summary`: data.frame with summary statistics (mean, sd,
  quantiles, mode) for total AN.

## Details

The function first obtains posterior predicted effects at the observed
exposure values by calling
[`bcrosspred()`](https://pasahe.github.io/bdlnm/reference/bcrosspred.md)
with `exp_at = data[[name_exposure]]`. Predictions must include
relative-risk scale predictions so the previously fitted `"bdlnm"` model
must have a `log` or `logit` link; otherwise the function aborts. These
predictions require a defined centering value (`cen`), as attributable
measures are always computed with respect to a reference exposure. This
reference exposure is usually an optimal exposure computed with the
[`optimal_exposure()`](https://pasahe.github.io/bdlnm/reference/optimal_exposure.md)
function, such as the Minimum Mortality Temperature (MMT).

Two different algorithms can be chosen to calculate attributable
measures:

- Backward (`dir = "back"`): for each time point, contributions from
  past exposures (over the lag window) are combined to calculate the
  daily AF/AN.

- Forward (`dir = "forw"`): for each time point, the contribution of
  that time point exposure to future outcomes (over the lag window) is
  combined to calculate the daily AF/AN.

Both algorithms are fully described by Gasparrini and Leone (2014)
<doi:10.1186/1471-2288-14-55>.

Required columns to calculate `AF` and `AN` are `name_exposure` and
`name_cases` columns. If `name_cases` is not supplied only AF per time
can be computed and the output will only contain two elements: `$af` and
`$af.summary`. If `name_date` is provided the function checks that dates
are equispaced (checks seconds, minutes, hours, days, weeks, months or
years). Time series have to be equispaced because the algorithms used to
calculate attributable measures rely on consecutive time points over the
lag window,. For example, if you only have seasonal observations (e.g.,
summers) expand the data to the full sequence and insert `NA` for
missing exposures/cases and use `name_filter` to compute measures only
for the seasonal subset.

Only `"bdlnm"` objects fitted with a cross-basis are supported; models
fitted with a one-basis (no lag) are not suitable for attributable
calculations.

## Note

This function is inspired by `attrdl()` developed by Gasparrini and
Leone (2014) <doi:10.1186/1471-2288-14-55>. It has been adapted to work
in a Bayesian framework within the bdlnm package.

## References

Gasparrini A., Leone M. (2014). Attributable risk from distributed lag
models. *BMC Medical Research Methodology*, 14, 55.
<doi:10.1186/1471-2288-14-55>.

Quijal-Zamorano M., Martinez-Beneito M.A., Ballester J., Marí-Dell'Olmo
M. (2024). Spatial Bayesian distributed lag non-linear models (SB-DLNM)
for small-area exposure-lag-response epidemiological modelling.
*International Journal of Epidemiology*, 53(3), dyae061.
<doi:10.1093/ije/dyae061>.

## See also

[`bcrosspred()`](https://pasahe.github.io/bdlnm/reference/bcrosspred.md)
to predict exposure–lag–response associations for a `"bdlnm"` object,

[`bdlnm()`](https://pasahe.github.io/bdlnm/reference/bdlnm.md) to fit a
Bayesian distributed lag non-linear model (`"bdlnm"`),

[`optimal_exposure()`](https://pasahe.github.io/bdlnm/reference/optimal_exposure.md)
to estimate exposure values that optimize the predicted effect for a
`"bdlnm"` object.

## Author

Pau Satorra, Marcos Quijal-Zamorano.

## Examples

``` r
# Filter the dataset to reduce computational time:

# Exposure-response and lag-response spline parameters
dlnm_var <- list(
  var_prc = c(10, 75, 90),
  var_fun = "ns",
  lag_fun = "ns",
  max_lag = 21,
  lagnk = 3
)

# Cross-basis parameters
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

# Model

if (check_inla()) {
mod <- bdlnm(mort_75plus ~ cb + factor(dow) + seas,
             data = london,
             family = "poisson",
             sample.arg = list(seed = 432, seed = 1L))

# Predict
cpred <- bcrosspred(mod, exp_at = temp)

# compute centering (MMT) using optimal_exposure
mmt <- optimal_exposure(mod, exp_at = temp)
cen <- mmt$summary[["0.5quant"]]

# Attributable numbers and fractions (using the backwards algorithm):
attr <- attributable(mod, london, name_date = "date",
name_exposure = "tmean", name_cases = "mort_75plus", cen = cen, dir = "back")
}
#> Warning: Since 'seed!=0', parallel model is disabled and serial model is selected, num.threads='1:1'
#> Warning: Centering value unspecified (`cen`). Automatically set to: 12.85.
#> Registered S3 method overwritten by 'crs':
#>   method    from
#>   print.crs sf  

```
