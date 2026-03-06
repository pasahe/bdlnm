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
  mod_2 <- bdlnm(
    mort_75plus ~ ob + factor(dow) + seas,
    data = slondon,
    family = "poisson",
    sample.arg = list(n = n_sim)
  )
  expect_warning(cpred_2 <- bcrosspred(mod_2, exp_at = temp))
  mmt_2 <- optimal_exposure(mod_2, exp_at = temp)

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

test_that("optimal_exposure works with default at (reconstructs grid from basis)", {
  skip_on_cran()
  skip_if_not(check_inla(), "INLA not available")

  # Should not error when 'at' is not provided
  expect_silent(mmt2 <- optimal_exposure(mod))

  # Must return class and exp attribute
  expect_s3_class(mmt2, "optimal_exposure")

  expect_equal(length(attr(mmt2, "exp_at")), 48)
})

test_that("optimal_exposure errors when object is not a bdlnm output or basis wrong", {
  skip_on_cran()
  skip_if_not(check_inla(), "INLA not available")

  # NULL x should error
  expect_snapshot_error(optimal_exposure(exp_at = temp))

  # NULL basis should not error
  expect_silent(optimal_exposure(mod, exp_at = temp))

  # wrong basis type should error
  expect_snapshot_error(optimal_exposure(mod, basis = "ob", exp_at = temp))
})
