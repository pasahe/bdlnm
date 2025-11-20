# Use data internal for tests:

# filter london dataset to run faster
slondon <- london[london$year >= 2012,]

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
               knots = stats::quantile(slondon$tmean,
                                dlnm_var$var_prc/100, na.rm = TRUE),
               Bound = range(slondon$tmean, na.rm = TRUE))

arglag <- list(fun = dlnm_var$lag_fun,
               knots = dlnm::logknots(dlnm_var$max_lag, nk = dlnm_var$lagnk))

# Create crossbasis
cb <- dlnm::crossbasis(slondon$tmean, lag = dlnm_var$max_lag, argvar, arglag)

# Seasonality of mortality time series
seas <- splines::ns(slondon$date, df = round(8 * length(slondon$date) / 365.25))

# Prediction values (equidistant points)
temp <- round(seq(min(slondon$tmean), max(slondon$tmean), by = 0.1), 1)

# Model

# Fit bdlnm with a small number of posterior samples to keep tests fast
n_sim <- 10

mod <- bdlnm(mort_75plus ~ cb + factor(dow) + seas,
           basis = cb,
           data = slondon,
           family = "poisson",
           sample.arg = list(n = n_sim, seed = 1L))

# Predict
cpred <- bcrosspred(mod, cb, at = temp)

# compute centering (MMT) using minimum_effect
mmt <- minimum_effect(mod, cb, at = temp)
cen <- mmt$min.summary[["0.5quant"]]
