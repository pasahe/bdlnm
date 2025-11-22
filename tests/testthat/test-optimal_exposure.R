
test_that("optimal_exposure returns object of class optimal_exposure with expected components", {

  # Basic structure
  expect_s3_class(mmt, "optimal_exposure")
  expect_true(is.list(mmt))
  expect_equal(names(mmt), c("est", "summary"))

  # samples length: one optimum per posterior sample
  n_sim <- attr(mod, "n_sim")
  expect_equal(length(mmt$est), n_sim)

  # names of mmt vector follow sample convention
  expect_equal(names(mmt$est), paste0("sample", seq_len(n_sim)))

  # optimal values belong to the prediction grid (xvar)
  expect_equal(attr(mmt, "xvar"), temp)
  expect_true(all(mmt$est %in% as.numeric(attr(mmt, "xvar"))))

  # summary format: must contain mean, sd, mode and at least one quantile
  ms_cols <- names(mmt$summary)
  expect_true(all(c("mean", "sd", "mode") %in% ms_cols))
  expect_true(any(grepl("quant$", ms_cols)))
})

test_that("optimal_exposure works with default at (reconstructs grid from basis)", {

  # Should not error when 'at' is not provided
  expect_silent(mmt2 <- optimal_exposure(mod, cb))

  # Must return class and xvar attribute
  expect_s3_class(mmt2, "optimal_exposure")

  expect_equal(length(attr(mmt2, "xvar")), 48)

})

test_that("optimal_exposure errors when object is not a bdlnm output or basis wrong", {

  # NULL x should error
  expect_snapshot_error(optimal_exposure(cb, at = temp))

  # NULL basis should error
  expect_snapshot_error(optimal_exposure(mod, at = temp))

  # wrong basis type should error
  expect_snapshot_error(optimal_exposure(mod, basis = 1:5, at = temp))

})
