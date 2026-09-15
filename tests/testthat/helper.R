# Use data internal for tests:

# filter london dataset to run faster
slondon <- london[london$year >= 2011, ]

# Exposure-response and lag-response spline parameters
dlnm_var <- list(
  var_prc = c(10, 75, 90),
  var_fun = "ns",
  lag_fun = "ns",
  max_lag = 21,
  lagnk = 3
)


# Cross-basis parameters
argvar <- list(
  fun = dlnm_var$var_fun,
  knots = stats::quantile(slondon$tmean, dlnm_var$var_prc / 100, na.rm = TRUE),
  Bound = range(slondon$tmean, na.rm = TRUE)
)

arglag <- list(
  fun = dlnm_var$lag_fun,
  knots = dlnm::logknots(dlnm_var$max_lag, nk = dlnm_var$lagnk)
)

# Create crossbasis
cb <- dlnm::crossbasis(slondon$tmean, lag = dlnm_var$max_lag, argvar, arglag)

# Seasonality of mortality time series
seas <- splines::ns(slondon$date, df = round(8 * length(slondon$date) / 365.25))

# Prediction values (equidistant points)
temp <- round(seq(min(slondon$tmean), max(slondon$tmean), by = 0.1), 1)
# Ensure it falls inside the range of temperatures after rounding:
temp <- temp[temp >= min(slondon$tmean) & temp <= max(slondon$tmean)]

# Model

# Fit bdlnm with a small number of posterior samples to keep tests fast
n_sim <- 10

if (check_inla()) {
  mod <- bdlnm(
    mort_75plus ~ cb + factor(dow) + seas,
    data = slondon,
    family = "poisson",
    sample.arg = list(n = n_sim, seed = 1L)
  )

  # Predict
  cpred <- bcrosspred(mod, exp_at = temp)

  # Compute centering (MMT) using optimal_exposure
  mmt <- optimal_exposure(cpred)
  cen <- mmt$summary[["0.5quant"]]

  # Shared fixed-interaction fixtures
  interaction_data <- slondon
  interaction_data$group <- factor(
    rep(c("a", "b"), length.out = nrow(interaction_data))
  )
  interaction_data$id_group <- rep(
    seq_len(5), length.out = nrow(interaction_data)
  )

  mod_interaction <- suppressWarnings(
    bdlnm(
      mort_75plus ~ cb * group + factor(dow) + seas +
        f(id_group, model = "iid"),
      data = interaction_data,
      family = "poisson",
      sample.arg = list(n = 2L, seed = 1L)
    )
  )

  cpred_fixed_interaction <- suppressWarnings(
    bcrosspred(
      mod_interaction,
      combine_effects = ~ cb * group,
      exp_at = interaction_data$tmean,
      cen = cen
    )
  )

  # Shared random-interaction fixtures
  random_data <- slondon
  random_data$id_group <- rep(seq_len(5), length.out = nrow(random_data))

  mod_random <- suppressWarnings(
    bdlnm(
      mort_75plus ~ cb + factor(dow) + seas +
        f(id_group, cb, model = "iid"),
      data = random_data,
      family = "poisson",
      sample.arg = list(n = 2L, seed = 1L)
    )
  )

  cpred_random_interaction <- suppressWarnings(
    bcrosspred(
      mod_random,
      combine_effects = ~ cb + f(id_group, cb, model = "iid"),
      exp_at = random_data$tmean,
      cen = cen
    )
  )

  # Centered prediction used by attributable tests
  cpred_centered <- bcrosspred(mod, exp_at = unique(slondon$tmean), cen = cen)
}
