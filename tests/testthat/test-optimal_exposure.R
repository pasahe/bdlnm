test_that("optimal_exposure returns object of class optimal_exposure with expected components (crossbasis)", {
  skip_on_cran()
  skip_if_not(check_inla(), "INLA not available")

  # Basic structure
  expect_s3_class(mmt, "optimal_exposure")
  expect_true(is.list(mmt))
  expect_equal(names(mmt), c("est", "summary"))

  # samples length: one optimum per posterior sample
  n_sim <- attr(mod, "n_sim")
  expect_equal(length(mmt$est), n_sim)

  # names of mmt vector follow sample convention
  expect_equal(names(mmt$est), paste0("sample", seq_len(n_sim)))

  # optimal values belong to the prediction grid (exp)
  expect_equal(attr(mmt, "exp_at"), temp)
  expect_true(all(mmt$est %in% as.numeric(attr(mmt, "exp"))))

  expect_equal(
    attr(mmt, "lag_at"),
    seq(attr(cb, "lag")[1], attr(cb, "lag")[2], by = 1)
  )

  # summary format: must contain mean, sd, mode and at least one quantile
  ms_cols <- names(mmt$summary)
  expect_true(all(c("mean", "sd", "mode") %in% ms_cols))
  expect_true(any(grepl("quant$", ms_cols)))

})

test_that("optimal_exposure returns object of class optimal_exposure with expected components (onebasis)", {
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

  # Basic structure
  expect_s3_class(mmt_2, "optimal_exposure")
  expect_true(is.list(mmt_2))
  expect_equal(names(mmt_2), c("est", "summary"))

  # samples length: one optimum per posterior sample
  n_sim <- attr(mod, "n_sim")
  expect_equal(length(mmt_2$est), n_sim)

  # names of mmt vector follow sample convention
  expect_equal(names(mmt_2$est), paste0("sample", seq_len(n_sim)))

  # optimal values belong to the prediction grid (exp)
  expect_equal(attr(mmt_2, "exp_at"), temp)
  expect_true(all(mmt_2$est %in% as.numeric(attr(mmt_2, "exp_at"))))

  expect_equal(attr(mmt_2, "lag_at"), NULL)

  # summary format: must contain mean, sd, mode and at least one quantile
  ms_cols <- names(mmt_2$summary)
  expect_true(all(c("mean", "sd", "mode") %in% ms_cols))
  expect_true(any(grepl("quant$", ms_cols)))
})

test_that("optimal_exposure uses the grid stored in bcrosspred", {
  skip_on_cran()
  skip_if_not(check_inla(), "INLA not available")

  cpred_default <- suppressWarnings(bcrosspred(mod))
  expect_silent(mmt2 <- optimal_exposure(cpred_default))

  # Must return class and exp attribute
  expect_s3_class(mmt2, "optimal_exposure")

  expect_equal(length(attr(mmt2, "exp_at")), 48)
})

test_that("optimal_exposure calculates fixed-interaction results", {
  skip_on_cran()
  skip_if_not(check_inla(), "INLA not available")

  optimal_fixed <- optimal_exposure(cpred_fixed_interaction)

  expect_s3_class(optimal_fixed, "optimal_exposure")
  expect_identical(
    names(optimal_fixed),
    c("group[a]", "group[b]")
  )

  for (result_name in names(optimal_fixed)) {
    expect_s3_class(optimal_fixed[[result_name]], "optimal_exposure")
    expected <- cpred_fixed_interaction[[result_name]]$exp_at[apply(
      cpred_fixed_interaction[[result_name]]$allfit,
      2L,
      which.min
    )]
    names(expected) <- paste0("sample", seq_along(expected))
    expect_identical(optimal_fixed[[result_name]]$est, expected)
  }
})

test_that("optimal_exposure calculates random-interaction results", {
  skip_on_cran()
  skip_if_not(check_inla(), "INLA not available")

  optimal_random <- optimal_exposure(cpred_random_interaction)

  expected_names <- c("global", paste0("id_group[", seq_len(5), "]"))
  expect_s3_class(optimal_random, "optimal_exposure")
  expect_identical(names(optimal_random), expected_names)

  for (result_name in expected_names) {
    expect_s3_class(optimal_random[[result_name]], "optimal_exposure")
  }
})

test_that("optimal_exposure validates its object and prediction arguments", {
  skip_on_cran()
  skip_if_not(check_inla(), "INLA not available")

  # NULL x should error
  expect_error(
    optimal_exposure(exp_at = temp),
    "argument .*object.* is missing"
  )

  # Direct model use is retained only as a deprecated compatibility path
  expect_warning(
    optimal_exposure(mod, exp_at = temp),
    "bdlnm.*deprecated"
  )

  expect_error(
    optimal_exposure(mod_interaction, exp_at = temp),
    "containing basis interactions cannot be used directly"
  )

  expect_error(
    optimal_exposure(cpred, exp_at = temp),
    "Prediction arguments cannot be supplied"
  )
})
