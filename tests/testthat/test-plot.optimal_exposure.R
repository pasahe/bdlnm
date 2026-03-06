test_that("does a plot of the optimal effect exposure value (crossbasis)", {
  skip_on_cran()
  skip_if_not(check_inla(), "INLA not available")

  expect_no_error(plot(
    mmt,
    xlab = "Temperature (ºC)",
    main = paste0("MMT (Median = ", round(mmt$summary[["0.5quant"]], 1), "ºC)")
  ))
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
  mod_2 <- bdlnm(
    mort_75plus ~ ob + factor(dow) + seas,
    data = slondon,
    family = "poisson",
    sample.arg = list(n = n_sim)
  )
  expect_warning(cpred_2 <- bcrosspred(mod_2, exp_at = temp))
  mmt_2 <- optimal_exposure(mod_2, exp_at = temp)
  expect_no_error(plot(
    mmt_2,
    xlab = "Temperature (ºC)",
    main = paste0("MMT (Median = ", round(mmt$summary[["0.5quant"]], 1), "ºC)")
  ))
})
