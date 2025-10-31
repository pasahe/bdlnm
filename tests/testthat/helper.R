# Use data internal for tests:


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
temp <- seq(round(min(london$tmean), 1), round(max(london$tmean), 1), by = 0.1)
