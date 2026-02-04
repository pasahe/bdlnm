test_that("ensure bdlnm returned expected structure (crossbasis example)", {

  # mod build in helper.R
  expect_type(mod, "list")
  expect_equal(length(mod), 3L)
  expect_equal(attr(mod, "n_sim"), n_sim)
  expect_equal(class(mod$model), "inla")
  expect_equal(class(mod$coefficients), c("matrix", "array"))
  expect_equal(dim(mod$coefficients), c(35, n_sim))
  expect_equal(class(mod$coefficients.summary), c("matrix", "array"))
  expect_equal(dim(mod$coefficients.summary), c(35, 6L))

})

test_that("ensure bdlnm returned expected structure (onebasis example)", {

  ob <- dlnm::onebasis(slondon$tmean, "strata", breaks = c(5, 10, 20))
  expect_warning(
    mod_2 <- bdlnm(
      mort_75plus ~ ob + factor(dow) + seas,
      basis = ob,
      data = slondon,
      family = "poisson",
      sample.arg = list(n = n_sim, seed = 1L)
    ))
  expect_type(mod_2, "list")
  expect_equal(length(mod_2), 3L)
  expect_equal(attr(mod_2, "n_sim"), n_sim)
  expect_equal(class(mod_2$model), "inla")
  expect_equal(class(mod_2$coefficients), c("matrix", "array"))
  expect_equal(dim(mod_2$coefficients), c(18, n_sim))
  expect_equal(class(mod_2$coefficients.summary), c("matrix", "array"))
  expect_equal(dim(mod_2$coefficients.summary), c(18, 6L))

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
