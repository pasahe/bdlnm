#review
test_that("does an overall plot", {

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
  expect_no_error(plot(cpred , "slices", var = htemp, col=3, ylab="RR",
       main=paste0("Association for a high temperature (", htemp, "ºC)")))

  expect_no_error(plot(cpred , "slices", lag = 0, col=4, ylab="RR",
       main=paste0("Association at Lag 0")))

})

