
test_that("bcrosspred basic structure (crossbasis)", {

  #cpred build in helper.R

  expect_s3_class(cpred, "bcrosspred")
  expect_equal(length(cpred), 17)
  expect_equal(names(cpred), c("predvar", "cen", "lag", "bylag", "coefficients", "matfit", "allfit", "matRRfit", "allRRfit", "coefficients.summary", "matfit.summary", "allfit.summary", "matRRfit.summary", "allRRfit.summary", "ci.level", "model.class", "model.link"))

  expect_equal(cpred$predvar, temp)
  expect_equal(cpred$cen, 12.5)
  expect_equal(cpred$lag, attr(cb, "lag"))
  expect_equal(cpred$bylag, 1)
  expect_equal(dim(cpred$coefficients), c(20, n_sim))
  expect_equal(rownames(cpred$coefficients), colnames(cb))
  expect_equal(dim(cpred$matfit), c(length(temp), 22, n_sim))
  expect_equal(rownames(cpred$matfit), as.character(temp))
  expect_equal(dim(cpred$allfit), c(length(temp), n_sim))
  expect_equal(rownames(cpred$allfit), as.character(temp))
  expect_true(all(cpred$allRRfit > 0))
  expect_equal(dim(cpred$matfit.summary), c(length(temp), 22, 6))
  expect_equal(cpred$ci.level, 0.95)
  expect_equal(cpred$model.class, "inla")
  expect_equal(cpred$model.link, "log")

})

test_that("bcrosspred basic structure (onebasis)", {

  ob <- dlnm::onebasis(slondon$tmean, "strata", breaks = c(5, 10, 20))
  expect_warning(
    mod_2 <- bdlnm(
      mort_75plus ~ ob + factor(dow) + seas,
      basis = ob,
      data = slondon,
      family = "poisson",
      sample.arg = list(n = n_sim, seed = 1L)
    )
  )
  cpred_2 <- bcrosspred(mod_2, ob, at = temp)

  expect_s3_class(cpred_2, "bcrosspred")
  expect_equal(length(cpred_2), 16)
  expect_equal(names(cpred_2), c("predvar", "lag", "bylag", "coefficients", "matfit", "allfit", "matRRfit", "allRRfit", "coefficients.summary", "matfit.summary", "allfit.summary", "matRRfit.summary", "allRRfit.summary", "ci.level", "model.class", "model.link"))

  expect_equal(cpred_2$predvar, temp)
  expect_equal(cpred_2$lag, c(0, 0))
  expect_equal(cpred_2$bylag, 1)
  expect_equal(dim(cpred_2$coefficients), c(3, n_sim))
  expect_equal(rownames(cpred_2$coefficients), colnames(ob))
  expect_equal(dim(cpred_2$matfit), c(length(temp), 1, n_sim))
  expect_equal(rownames(cpred_2$matfit), as.character(temp))
  expect_equal(dim(cpred_2$allfit), c(length(temp), n_sim))
  expect_equal(rownames(cpred_2$allfit), as.character(temp))
  expect_true(all(cpred_2$allRRfit > 0))
  expect_equal(dim(cpred_2$matfit.summary), c(length(temp), 1, 6))
  expect_equal(cpred_2$ci.level, 0.95)
  expect_equal(cpred_2$model.class, "inla")
  expect_equal(cpred_2$model.link, "log")

})

test_that("bcrosspred errors when some argument is missing or invalid", {

  expect_snapshot_error(bcrosspred(NULL, cb, at = temp))

  expect_snapshot_error(bcrosspred(list(a = 1), cb, at = temp))

  expect_snapshot_error(bcrosspred(mod, NULL, at = temp))

  expect_snapshot_error(bcrosspred(mod, basis = 1:5, at = temp))

  expect_warning(cpred2 <- bcrosspred(mod, cb))

  expect_equal(length(cpred2$predvar), 48)

})
