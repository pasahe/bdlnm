#review
test_that("does a plot of the optimal effect exposure value", {

  expect_no_error(plot(mmt, xlab = "Temperature (ºC)", main = paste0("MMT (Median = ", round(mmt$summary[["0.5quant"]], 1), "ºC)")))

})
