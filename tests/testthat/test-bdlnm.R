test_that("ensure bdlnm returned expected structure", {

  expect_type(mod, "list")
  expect_equal(length(mod), 3L)
  expect_equal(attr(mod, "n_sim"), n_sim)
  expect_equal(class(mod$model), "inla")
  expect_equal(class(mod$coefficients), c("matrix", "array"))
  expect_equal(dim(mod$coefficients), c(32, n_sim))
  expect_equal(class(mod$coefficients.summary), c("matrix", "array"))
  expect_equal(dim(mod$coefficients.summary), c(32, 6L))

})

test_that("bdlnm errors when required arguments are missing or inappropiate", {

  expect_snapshot_error(
    bdlnm(mort_75plus ~ cb + factor(dow) + seas, basis = cb)
  )

  expect_snapshot_error(
    bdlnm(mort_75plus ~ cb + factor(dow) + seas, data = slondon)
  )

  expect_snapshot_error(
    bdlnm(mort_75plus ~ cb + factor(dow) + seas, basis = cb, data = slondon, sample.arg = 5)
  )

  expect_snapshot_error(
    bdlnm(mort_75plus ~ cb + factor(dow) + seas, basis = cb, data = slondon, sample.arg = list(not_an_argument = 5))
  )

  expect_snapshot_error(
    bdlnm(mort_75plus ~ cb + factor(dow) + seas, basis = cb, data = slondon, not_an_argument = 5)
  )

})

test_that("bdlnm honors sample.arg", {

  mod2 <- bdlnm(
    mort_75plus ~ cb + factor(dow) + seas,
    basis = cb,
    data = slondon,
    family = "poisson",
    sample.arg = list(n = 5L)
  )

  expect_equal(attr(mod2, "n_sim"), 5L)
  expect_equal(ncol(mod2$coefficients), 5L)

})

test_that("bdlnm throws informative error if control.compute config = FALSE is passed", {

  expect_snapshot_error(
    bdlnm(
      mort_75plus ~ cb + factor(dow) + seas,
      basis = cb,
      data = slondon,
      family = "poisson",
      control.compute = list(config = FALSE)
    )
  )

})

test_that("bdlnm returns sensible summaries (means, sd, quantiles, mode)", {

  coefsum <- mod$coefficients.summary

  # check summary columns include mean, sd and at least one quantile and mode
  expect_true(all(c("mean", "sd", "mode") %in% colnames(coefsum)))
  expect_true(any(grepl("quant$", colnames(coefsum))))

})

