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

cpred <- bcrosspred(cb, mod, at = temp)

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
