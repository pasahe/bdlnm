library(splines)
library(dlnm)
library(INLA)

dlnm_var <- list(
  var_prc = c(10, 75, 90),
  var_fun = "ns",
  lag_fun = "ns",
  max_lag = 21,
  lagnk = 3
)

seas <- ns(london$date, df = round(8 * length(london$date) / 365.25))

# Cross-basis parameters
argvar <- list(fun = dlnm_var$var_fun,
               knots = quantile(london$tmean,
                                dlnm_var$var_prc/100, na.rm = TRUE),
               Bound = range(london$tmean, na.rm = TRUE))

arglag <- list(fun = dlnm_var$lag_fun,
               knots = logknots(dlnm_var$max_lag, nk = dlnm_var$lagnk))

# Prediction values
temp <-
  c(quantile(london$tmean, seq(0, 0.9, by = 0.1)/100),
    quantile(london$tmean, seq(1, 99, by = 1)/100),
    quantile(london$tmean, seq(99.1, 100, by = 0.1)/100))

lags <- 0:dlnm_var$max_lag

# Create crossbasis
cb <- crossbasis(london$tmean, lag = dlnm_var$max_lag, argvar, arglag)

# Omit NA's
attr_names <- names(attributes(cb))
attr_names <- attr_names[! attr_names %in% c("dim", "dimnames")]
attr <- attributes(cb)[attr_names]

cb <- cb[-c(1:max(dlnm_var$max_lag)),]

class(cb) <- c("crossbasis","matrix")
attributes(cb)[attr_names] <- attr

london_2 <- london[-c(1:dlnm_var$max_lag),]
seas <- seas[-c(1:dlnm_var$max_lag),]

# Fit the model (75+ years)
mod <- inla(mort_75plus ~ cb + factor(dow) + seas, data = london_2, family = "poisson", control.compute = list(config = TRUE))

names_sel <- grep("^v", rownames(mod$summary.fixed), value = TRUE)
list_sel <- lapply(1:length(names_sel), function(x) 1)
names(list_sel) <- names_sel

coef <- inla.posterior.sample(n = 100, mod, selection = list_sel)

test_that("works with inla model", {

  expect_no_error(cpred <- bcrosspred(cb, mod, at = temp))
  expect_s3_class(cpred, c("bcrosspred", "crosspred"))
  expect_equal(cpred$predvar, as.numeric(temp))
  expect_equal(cpred$cen, 10)
  expect_equal(cpred$lag, c(0, 21))
  expect_equal(cpred$bylag, 1)
  #Check matrix dimensions
  expect_equal(dim(cpred$coefficients), c(20, 100))
  expect_equal(dim(cpred$matfit), c(2618, 100))
  expect_equal(dim(cpred$allfit), c(119, 100))
  expect_equal(dim(cpred$matRRfit), c(2618, 100))
  expect_equal(dim(cpred$allRRfit), c(119, 100))
  expect_equal(cpred$model.class, "inla")
  expect_equal(cpred$model.link, "log")
})

test_that("works with coefficient", {
  expect_no_error(cpred <- bcrosspred(cb, coef = coef, model.link = "log", at = temp))
  expect_s3_class(cpred, c("bcrosspred", "crosspred"))
  expect_equal(cpred$predvar, as.numeric(temp))
  expect_equal(cpred$cen, 10)
  expect_equal(cpred$lag, c(0, 21))
  expect_equal(cpred$bylag, 1)
  #Check matrix dimensions
  expect_equal(dim(cpred$coefficients), c(20, 100))
  expect_equal(dim(cpred$matfit), c(2618, 100))
  expect_equal(dim(cpred$allfit), c(119, 100))
  expect_equal(dim(cpred$matRRfit), c(2618, 100))
  expect_equal(dim(cpred$allRRfit), c(119, 100))
  expect_equal(cpred$model.class, "inla")
  expect_equal(cpred$model.link, "log")
})

test_that("error if no models nor coefficients are supplied", {
  expect_snapshot(bcrosspred(cb, at = temp), error = TRUE)
})

test_that("error if both model and coefficients are supplied", {
  expect_snapshot(bcrosspred(cb, mod = mod, coef = coef, at = temp), error = TRUE)
})

test_that("error if another kind of model is supplied", {
  mod_2 <- glm(mort_75plus ~ cb + factor(dow) + seas, data = london_2, family = "poisson")
  expect_snapshot(bcrosspred(cb, mod = mod_2, at = temp), error = TRUE)
})

test_that("error if inla model doesn't have `control.config=TRUE`", {
  mod_2 <- inla(mort_75plus ~ cb + factor(dow) + seas, data = london_2, family = "poisson")
  expect_snapshot(bcrosspred(cb, mod = mod_2, at = temp), error = TRUE)
})

test_that("error if 'coef' is provided without a ' model.link'", {
  expect_snapshot(bcrosspred(cb, coef = coef, at = temp), error = TRUE)
})

test_that("no error if the prediction points are not specified", {
  expect_no_error(bcrosspred(cb, mod))
})
