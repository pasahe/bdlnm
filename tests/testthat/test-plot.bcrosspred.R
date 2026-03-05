#review
test_that("does an overall plot (crossbasis)", {

  expect_no_error(plot(cpred, "overall", xlab = "Temperature (ºC)", ylab = "Relative Risk", col = 4, main="Overall", log = "y"))
  expect_no_error(plot(cpred, "overall", xlab = "Temperature (ºC)", ylab = "Relative Risk", col = 4, main="Overall", log = "y", ci.level = 0.99))
  expect_no_error(plot(cpred, "overall", xlab = "Temperature (ºC)", ylab = "Relative Risk", col = 4, main="Overall", log = "y", ci = "sampling"))

})


test_that("does a 3d plot", {

  expect_no_error(plot(cpred, "3d", zlab = "Relative risk", col = 4, lphi = 60, cex.axis = 0.6, xlab = "Temperature (ºC)", main = "3D graph of temperature effect"))
})

test_that("does a contour plot", {

  expect_no_error(plot(cpred, "contour", xlab = "Temperature (ºC)", ylab = "Lag", main="Contour plot"))

})

test_that("does a slice plot", {

  htemp <- 23
  expect_no_error(plot(cpred , "slices", exp_at = htemp, col=3, ylab="RR",
       main=paste0("Association for a high temperature (", htemp, "ºC)")))

  expect_no_error(plot(cpred , "slices", lag_at = 0, col=4, ylab="RR",
       main=paste0("Association at Lag 0")))

})

test_that("onebasis plots", {

  ob <- dlnm::onebasis(slondon$tmean, "strata", breaks = c(5, 10, 20))
  expect_warning(
    mod_2 <- bdlnm(
      mort_75plus ~ ob + factor(dow) + seas,
      data = slondon,
      family = "poisson",
      sample.arg = list(n = n_sim, seed = 1L)
    )
  )
  cpred_2 <- bcrosspred(mod_2, "ob", exp_at = temp)

  # overall
  expect_no_error(plot(cpred_2, "overall", xlab = "Temperature (ºC)", ylab = "Relative Risk", col = 4, log = "y"))
  expect_no_error(plot(cpred_2, "overall", xlab = "Temperature (ºC)", ylab = "Relative Risk", col = 4, log = "y", ci.level = 0.99))
  expect_no_error(plot(cpred_2, "overall", xlab = "Temperature (ºC)", ylab = "Relative Risk", col = 4, log = "y", ci = "sampling"))

  # 3d plot
  expect_error(plot(cpred_2, "3d", zlab = "Relative risk", col = 4, lphi = 60, cex.axis = 0.6, xlab = "Temperature (ºC)", main = "3D graph of temperature effect"))

  # contour plot
  expect_error(plot(cpred_2, "contour", xlab = "Temperature (ºC)", ylab = "Lag", main="Contour plot"))

  # slices
  htemp <- 23
  expect_error(plot(cpred_2, "slices", exp_at = htemp, col=3, ylab="RR", main=paste0("Association for a high temperature (", htemp, "ºC)")))
  expect_error(plot(cpred_2, "slices", lag_at = 0, col=4, ylab="RR", main=paste0("Association at Lag 0")))

})
