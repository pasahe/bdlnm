test_that("ensure bdlnm returned expected structure (crossbasis example)", {
  # mod build in helper.R
  expect_type(mod, "list")
  expect_equal(length(mod), 4L)
  expect_equal(attr(mod, "n_sim"), n_sim)
  expect_equal(class(mod$model), "inla")
  expect_equal(class(mod$coefficients), c("matrix", "array"))
  expect_equal(dim(mod$coefficients), c(35, n_sim))
  expect_equal(class(mod$coefficients.summary), c("matrix", "array"))
  expect_equal(dim(mod$coefficients.summary), c(35, 6L))
})

test_that("ensure bdlnm returned expected structure (onebasis example)", {
  ob <- dlnm::onebasis(slondon$tmean, "strata", breaks = c(5, 10, 20))
  mod_2 <- bdlnm(
    mort_75plus ~ ob + factor(dow) + seas,
    data = slondon,
    family = "poisson",
    sample.arg = list(n = n_sim)
  )
  expect_type(mod_2, "list")
  expect_equal(length(mod_2), 4L)
  expect_equal(attr(mod_2, "n_sim"), n_sim)
  expect_equal(class(mod_2$model), "inla")
  expect_equal(class(mod_2$coefficients), c("matrix", "array"))
  expect_equal(dim(mod_2$coefficients), c(18, n_sim))
  expect_equal(class(mod_2$coefficients.summary), c("matrix", "array"))
  expect_equal(dim(mod_2$coefficients.summary), c(18, 6L))
})


test_that("works with two different basis", {
  ob <- dlnm::onebasis(slondon$tmean, "strata", breaks = c(5, 10, 20))
  mod_2 <- bdlnm(
    mort_75plus ~ cb + ob + factor(dow) + seas,
    data = slondon,
    family = "poisson",
    sample.arg = list(n = n_sim)
  )
  expect_type(mod_2, "list")
  expect_length(mod_2$basis, 2)
  expect_equal(names(mod_2$basis), c("cb", "ob"))
  expect_equal(mod_2$basis[["cb"]], cb)
  expect_equal(mod_2$basis[["ob"]], ob)
})

test_that("bdlnm errors when required arguments are missing or inappropiate", {
  expect_snapshot_error(
    bdlnm(mort_75plus ~ cb + factor(dow) + seas, data = slondon, sample.arg = 5)
  )
})

test_that("bdlnm honors sample.arg", {
  mod2 <- bdlnm(
    mort_75plus ~ cb + factor(dow) + seas,
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

test_that("na.action = na.omit drops NA rows and returns consistent dimensions, and gives same coefficients as complete data model", {
  expect_equal(nrow(mod$model$model.matrix), nrow(na.omit(cb)))

  na_rows <- which(!complete.cases(cb))
  cb_clean <- cb[-na_rows, ]
  class(cb_clean) <- "crossbasis"
  seas_clean <- seas[-na_rows, ]
  slondon_clean <- slondon[-na_rows, ]

  mod_clean <- suppressWarnings(
    bdlnm(
      mort_75plus ~ cb_clean + factor(dow) + seas_clean,
      data = slondon_clean,
      family = "poisson",
      sample.arg = list(n = n_sim, seed = 1L)
    )
  )

  expect_identical(
    round(as.numeric(mod_clean$coefficients), 2),
    round(as.numeric(mod$coefficients), 2)
  )
})

test_that("na.action = na.pass keeps NA rows and returns consistent dimensions", {
  mod_na <- bdlnm(
    mort_75plus ~ cb + factor(dow) + seas,
    data = slondon,
    family = "poisson",
    na.action = na.pass,
    sample.arg = list(n = n_sim)
  )

  expect_equal(nrow(mod_na$model$model.matrix), nrow(slondon))
})

test_that("na.action = na.fail throws an error when NAs are present", {
  expect_error(
    bdlnm(
      mort_75plus ~ cb + factor(dow) + seas,
      data = slondon,
      family = "poisson",
      na.action = na.fail,
      sample.arg = list(n = n_sim)
    )
  )
})

test_that("na.action is ignored when a random effect is included", {
  slondon$id <- seq_len(nrow(slondon))
  expect_message(
    mod_rt <- bdlnm(
      mort_75plus ~ cb + factor(dow) + seas + f(id),
      data = slondon,
      family = "poisson",
      sample.arg = list(n = n_sim)
    )
  )
  expect_equal(nrow(mod_rt$model$model.matrix), nrow(slondon))
})
