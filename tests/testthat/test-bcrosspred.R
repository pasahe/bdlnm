test_that("works with inla model", {

  expect_no_error(cpred <- bcrosspred(bdlnm:::cb_london, bdlnm:::mod_london, at = temp))
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
  expect_no_error(cpred <- bcrosspred(bdlnm:::cb_london, coef = bdlnm:::coef_london, model.link = "log", at = temp))
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
  expect_snapshot(bcrosspred(bdlnm:::cb_london, at = temp), error = TRUE)
})

test_that("error if both model and coefficients are supplied", {
  expect_snapshot(bcrosspred(bdlnm:::cb_london, mod = bdlnm:::mod_london, coef = bdlnm:::coef_london, at = temp), error = TRUE)
})

test_that("error if another kind of model is supplied", {

  seas <- splines::ns(london$date, df = round(8 * length(london$date) / 365.25))
  london_2 <- london[-c(1:21),]
  seas_2 <- seas[-c(1:21),]
  mod_2 <- glm(mort_75plus ~ bdlnm:::cb_london + factor(dow) + seas_2, data = london_2, family = "poisson")
  expect_snapshot(bcrosspred(bdlnm:::cb_london, mod = mod_2, at = temp), error = TRUE)

})

test_that("error if 'coef' is provided without a ' model.link'", {
  expect_snapshot(bcrosspred(bdlnm:::cb_london, coef = bdlnm:::coef_london, at = temp), error = TRUE)
})

test_that("no error if the prediction points are not specified", {
  expect_no_error(bcrosspred(bdlnm:::cb_london, bdlnm:::mod_london))
})
