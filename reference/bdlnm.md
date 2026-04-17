# Fit a Bayesian distributed lag non-linear model (B-DLNM)

Fit a distributed lag non-linear model (DLNM) using a Bayesian
framework. The function calls
[`INLA::inla()`](https://rdrr.io/pkg/INLA/man/inla.html) to fit the
model and then draws posterior samples of the model fixed effects with
[`INLA::inla.posterior.sample()`](https://rdrr.io/pkg/INLA/man/posterior.sample.html).
See the package vignette for worked examples and recommended workflows.

## Usage

``` r
bdlnm(
  formula,
  data,
  family = "gaussian",
  sample.arg = list(n = 1000, seed = 0L),
  ci.level = 0.95,
  na.action = getOption("na.action"),
  ...
)
```

## Arguments

- formula:

  A model formula (as for
  [`INLA::inla()`](https://rdrr.io/pkg/INLA/man/inla.html)). The model
  must be a distributed lag linear or non-linear model (DLNM), so a
  `"crossbasis"`
  ([`dlnm::crossbasis()`](https://rdrr.io/pkg/dlnm/man/crossbasis.html))
  or a `"onebasis"`
  ([`dlnm::onebasis()`](https://rdrr.io/pkg/dlnm/man/onebasis.html))
  must be included.

- data:

  an optional data frame, list or environment (or object coercible by
  [as.data.frame](https://rdrr.io/r/base/as.data.frame.html) to a data
  frame) containing the variables in the model. If not found in data,
  the variables are taken from `environment(formula)`, typically the
  environment from which `bdlnm` is called.

- family:

  Character. Family name passed to
  [`INLA::inla()`](https://rdrr.io/pkg/INLA/man/inla.html) (default
  `"gaussian"`).

- sample.arg:

  List of arguments passed to
  [`INLA::inla.posterior.sample()`](https://rdrr.io/pkg/INLA/man/posterior.sample.html).
  Defaults to `list(n = 1000, seed = 0L)` (draws `1000` posterior
  samples; seed at random). For reproducible sampling set a non-zero
  numeric `seed`.

- ci.level:

  Numeric in `(0,1)` giving the credible interval level (default
  `0.95`). Credible interval quantiles are computed to summarize
  coefficients from the posterior samples.

- na.action:

  A function specifying how to handle NA values when constructing the
  model frame. The default is taken from the na.action setting of
  [options](https://rdrr.io/r/base/options.html), which is by default
  [na.omit](https://rdrr.io/r/stats/na.fail.html) (drops rows with any
  `NA` among the variables referenced in formula). In the presence of a
  random effect term `f()` in the formula, this argument is ignored and
  `NA`s are not discarded. When rows containing missing values are
  retained, `INLA` handles them internally (see `Details` below).

- ...:

  Additional arguments passed to
  [`INLA::inla()`](https://rdrr.io/pkg/INLA/man/inla.html).

## Value

An S3 object of class `"bdlnm"` with the following components:

- `model`: the fitted `INLA` model returned by
  [`INLA::inla()`](https://rdrr.io/pkg/INLA/man/inla.html).

- `basis`: a named list containing all basis of class `crossbasis` or
  `onebasis` included in the model `formula`.

- `coefficients`: a matrix whose columns are posterior sample draws
  returned by
  [`INLA::inla.posterior.sample()`](https://rdrr.io/pkg/INLA/man/posterior.sample.html)
  (named `sample1`, `sample2`, ...) and whose rows are all model
  coefficients.

- `coefficients.summary`: a matrix of summary statistics for all the
  posterior samples stored in `coefficients` (mean, sd, quantiles,
  mode).

## Distributed lag non-linear model

The fitted model must be a distributed lag linear or non-linear model
(DLNM). DLNMs describe potentially non-linear and delayed (lagged)
associations between an exposure and an outcome, commonly referred to as
exposure–lag–response relationships. This modelling framework is based
on the definition of a cross-basis (a bi-dimensional space of functions)
constructed with
[`dlnm::crossbasis()`](https://rdrr.io/pkg/dlnm/man/crossbasis.html),
which defines the exposure–response and lag–response functions
simultaneously. The cross-basis object must be created beforehand and
supplied as an object in the calling environment (not as a column inside
data) and explicitly included in the model formula (e.g., y ~ cb + ...).
A basis object constructed with
[`dlnm::onebasis()`](https://rdrr.io/pkg/dlnm/man/onebasis.html) can be
used instead when the model is restricted to a uni-dimensional
exposure–response relationship (i.e., without lagged effects). All basis
objects included in the model `formula` are stored and returned as a
named list in the `basis` component. Any of these basis objects can
later be supplied to
[`bcrosspred()`](https://pasahe.github.io/bdlnm/reference/bcrosspred.md)
to extract predictions for the corresponding exposure–lag–response (or
exposure-response, if created with `onebasis`) association.

## INLA

Models are fit using Integrated Nested Laplace approximation (INLA) via
[`INLA::inla()`](https://rdrr.io/pkg/INLA/man/inla.html). INLA is a
method for approximate Bayesian inference. In the last years it has
established itself as an alternative to other methods such as Markov
chain Monte Carlo because of its speed and ease of use via the R-INLA
package ([What is INLA?](https://www.r-inla.org/what-is-inla)).

Additional arguments supplied via `...` are forwarded to
[`INLA::inla()`](https://rdrr.io/pkg/INLA/man/inla.html) (see
documentation for all available arguments). Internally, the function
ensures that `control.compute = list(config = TRUE)` in order to enable
posterior sample drawing with
[`INLA::inla.posterior.sample()`](https://rdrr.io/pkg/INLA/man/posterior.sample.html).

In the presence of missing values in variables referenced in formula,
the `na.action` argument controls how the model frame is constructed. If
`na.action` is set to [na.omit](https://rdrr.io/r/stats/na.fail.html)
(the default set by [options](https://rdrr.io/r/base/options.html)), a
complete-case analysis is performed and any row with a missing value in
a variable appearing in formula is removed. If `na.action` is set to
[na.pass](https://rdrr.io/r/stats/na.fail.html) instead, missing values
are retained in the model frame and passed to INLA.

Note that when the model formula includes a random effect term,
specified via `f(k, model = ...)`, the `na.action` is ignored and rows
with missing values are not dropped.

When missing values are present in the data supplied to
[`INLA::inla()`](https://rdrr.io/pkg/INLA/man/inla.html) (either because
a random effect term is included or because `na.action = na.pass`),
`INLA` handles them internally as follows:

- If `NA` values occur in the response, the corresponding observation
  contributes nothing to the likelihood (the response is treated as
  unobserved for that observation).

- If `NA` values occur in fixed-effect covariates,
  [`INLA::inla()`](https://rdrr.io/pkg/INLA/man/inla.html) replaces them
  internally with zero so that the covariate does not contribute to the
  linear predictor for that observation.

- If `NA` values occur in a fixed-effect covariate that is a factor,
  this is not allowed unless `NA` is explicitly included as a level, or
  `control.fixed = list(expand.factor.strategy = "inla")` is specified.
  With this option, `NA` is interpreted similarly as in the fixed-effect
  case, producing no contribution from that covariate to the linear
  predictor.

- If `NA` values occur in a random effect, the random effect does not
  contribute to the linear predictor for the corresponding observation.

## Posterior samples

After fitting the model, the function draw samples from the approximate
posterior distribution of the latent field via
[`INLA::inla.posterior.sample()`](https://rdrr.io/pkg/INLA/man/posterior.sample.html).
These samples are collected into a matrix and summarized across samples
(mean, sd, quantiles and mode). For a `"crossbasis"` built from an
exposure basis with C parameters and a lag basis with L parameters,
there will be C × L cross-basis associated coefficients (named e.g.
v1.l1, ..., vC.lL). For a `"onebasis"` object the coefficients follow
the simpler form b1, ... bC.

Additional arguments supplied via `sample.arg` are forwarded to
[`INLA::inla.posterior.sample()`](https://rdrr.io/pkg/INLA/man/posterior.sample.html)
(see documentation for all available arguments). By default, the number
of samples is `1000`. Be aware of the computation and memory cost when
increasing the number of samples drawn. By default, the seed is set at
random. For reproducible samplings, you need to set a non-zero numeric
`seed` in `sample.arg`.

Posterior sample estimations are then summarized across samples using
mean, sd, credible-interval quantiles (the mid and the lower/upper tails
according to `ci.level`) and an approximate mode obtained from a kernel
density estimate.

## Requirements

The INLA package must be installed from the R-INLA repository ([R-INLA
Project](https://www.r-inla.org/)); if not available the function aborts
with a short instruction on how to install it.

## References

Quijal-Zamorano M., Martinez-Beneito M.A., Ballester J., Marí-Dell'Olmo
M. (2024). Spatial Bayesian distributed lag non-linear models (SB-DLNM)
for small-area exposure-lag-response epidemiological modelling.
*International Journal of Epidemiology*, 53(3), dyae061.
<doi:10.1093/ije/dyae061>.

Quijal-Zamorano M., Martinez-Beneito M.A., Ballester J., Marí-Dell'Olmo
M. (2025). Spatial Bayesian distributed lag non-linear models with
R-INLA. *International Journal of Epidemiology*, 54(4), dyaf120.
<doi:10.1093/ije/dyaf120>.

Gasparrini A. (2011). Distributed lag linear and non-linear models in R:
the package dlnm. *Journal of Statistical Software*, 43(8), 1-20.
<doi:10.18637/jss.v043.i08>.

Rue H., Martino S., Chopin N. (2009). Approximate Bayesian inference for
latent Gaussian models by using integrated nested Laplace
approximations. *Journal of the Royal Statistical Society: Series B*,
71(2), 319-392. <doi:10.1111/j.1467-9868.2008.00700.x>.

## See also

[`bcrosspred()`](https://pasahe.github.io/bdlnm/reference/bcrosspred.md)
to predict exposure–lag–response associations for a `bdlnm` object,

[`attributable()`](https://pasahe.github.io/bdlnm/reference/attributable.md)
to calculate attributable fractions and numbers for a `bdlnm` object,

[`optimal_exposure()`](https://pasahe.github.io/bdlnm/reference/optimal_exposure.md)
to estimate exposure values that optimise the predicted effect for a
`bdlnm` object.

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
}
#> Warning: Since 'seed!=0', parallel model is disabled and serial model is selected, num.threads='1:1'


```
