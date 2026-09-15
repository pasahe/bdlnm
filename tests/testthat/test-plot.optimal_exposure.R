test_that("does a plot of the optimal effect exposure value (crossbasis)", {
  skip_on_cran()
  skip_if_not(check_inla(), "INLA not available")

  expect_no_error(plot(
    mmt,
    xlab = "Temperature (ºC)",
    main = paste0("MMT (Median = ", round(mmt$summary[["0.5quant"]], 1), "ºC)")
  ))
})

test_that("selects a fixed-interaction optimal exposure to plot", {
  skip_on_cran()
  skip_if_not(check_inla(), "INLA not available")

  optimal_fixed <- optimal_exposure(cpred_fixed_interaction)

  expect_no_error(
    plot(optimal_fixed, select = "group[a]")
  )
  expect_no_error(
    plot(optimal_fixed, select = "group[b]")
  )
  expect_error(
    plot(optimal_fixed),
    "select.*Available optimal-exposure results"
  )
  expect_error(
    plot(optimal_fixed, select = "group[c]"),
    "Optimal-exposure result.*was not found"
  )
})

test_that("selects a random-interaction optimal exposure to plot", {
  skip_on_cran()
  skip_if_not(check_inla(), "INLA not available")

  optimal_random <- optimal_exposure(cpred_random_interaction)

  expect_no_error(plot(optimal_random, select = "global"))
  expect_no_error(plot(optimal_random, select = "id_group[1]"))
})

test_that("does a plot of the optimal effect exposure value (onebasis)", {
  skip_on_cran()
  skip_if_not(check_inla(), "INLA not available")

  ob <- dlnm::onebasis(
    slondon$tmean,
    fun = dlnm_var$var_fun,
    knots = stats::quantile(
      slondon$tmean,
      dlnm_var$var_prc / 100,
      na.rm = TRUE
    ),
    Bound = range(slondon$tmean, na.rm = TRUE)
  )
  mod_2 <- suppressWarnings(bdlnm(
    mort_75plus ~ ob + factor(dow) + seas,
    data = slondon,
    family = "poisson",
    sample.arg = list(n = n_sim, seed = 1L)
  ))
  expect_warning(cpred_2 <- bcrosspred(mod_2, exp_at = temp))
  mmt_2 <- optimal_exposure(cpred_2)
  expect_no_error(plot(
    mmt_2,
    xlab = "Temperature (ºC)",
    main = paste0("MMT (Median = ", round(mmt$summary[["0.5quant"]], 1), "ºC)")
  ))
})
