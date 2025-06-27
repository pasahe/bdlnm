# Prediction values
temp <-
  c(quantile(london$tmean, seq(0, 0.9, by = 0.1)/100),
    quantile(london$tmean, seq(1, 99, by = 1)/100),
    quantile(london$tmean, seq(99.1, 100, by = 0.1)/100))

cpred <- bcrosspred(bdlnm:::cb_london, bdlnm:::mod_london, at = temp)

test_that("does an overall plot", {

  expect_no_error(plot(cpred, "overall", xlab = "Temperature (ºC)", ylab = "Relative Risk", col = 4, main="Overall", log = "y")
)
})

test_that("does a 3d plot", {

  expect_no_error(plot(cpred, zlab = "Relative risk", col = 4, lphi = 60, cex.axis = 0.6, xlab = "Temperature (ºC)")
)
})

test_that("does a contour plot", {

  expect_no_error(plot(cpred, "contour", xlab = "Temperature (ºC)", ylab = "Lag", main="Contour plot")
)
})

test_that("does a slice plot", {

  expect_no_error(plot(cpred , "slices", var = temp["99%"], col=3, ylab="RR", ci.arg=list(density=15,lwd=2),
                       main=paste0("Association for a high temperature (", round(temp["99%"], 0), "ºC)"))
  )
})
