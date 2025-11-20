
test_that("minimum_effect returns object of class min.risk with expected components", {

  # Basic structure
  expect_s3_class(mmt, "min.risk")
  expect_true(is.list(mmt))
  expect_equal(names(mmt), c("min", "min.summary"))

  # samples length: one minimum per posterior sample
  n_sim <- attr(mod, "n_sim")
  expect_equal(length(mmt$min), n_sim)

  # names of min vector follow sample convention
  expect_equal(names(mmt$min), paste0("sample", seq_len(n_sim)))

  # min values belong to the prediction grid (xvar)
  expect_equal(attr(mmt, "xvar"), temp)
  expect_true(all(mmt$min %in% as.numeric(attr(mmt, "xvar"))))

  # min.summary format: must contain mean, sd, mode and at least one quantile
  ms_cols <- colnames(mmt$min.summary)
  expect_true(all(c("mean", "sd", "mode") %in% ms_cols))
  expect_true(any(grepl("quant$", ms_cols)))
})

test_that("minimum_effect works with default at (reconstructs grid from basis)", {

  # Should not error when 'at' is not provided
  expect_silent(mmt2 <- minimum_effect(mod, cb))

  # Must return class and xvar attribute
  expect_s3_class(mmt2, "min.risk")

  expect_equal(length(attr(mmt2, "xvar")), 54)

})

test_that("minimum_effect errors when x is not a bdlnm output or basis wrong", {

  # NULL x should error
  expect_snapshot_error(minimum_effect(cb, at = temp))

  # NULL basis should error
  expect_snapshot_error(minimum_effect(mod, at = temp))

  # wrong basis type should error
  expect_snapshot_error(minimum_effect(mod, basis = 1:5, at = temp))

})
