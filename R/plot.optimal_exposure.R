#' Plot posterior distribution of optimal effect exposure values
#'
#' Plot the posterior distribution of the optimal effect exposure values returned by `optimal_exposure()`.
#'
#'  Plot an histogram showing the posterior distribution of the optimal effect exposure values returned by `optimal_exposure()`. The original prediction grid will be used as axis breaks. A dashed red vertical line marks the posterior median of these values.
#' @param x x of class `optimal_exposure` as returned by [optimal_exposure].
#' @param line.arg List of graphical arguments for the plotting of the median vertical line passed to [graphics::abline].
#' @param ... Additional graphical parameters passed to [hist].
#'
#' @export
#'
#' @examples
#'
#' # Set exposure-response and lag-response spline parameters
#'  dlnm_var <- list(
#'    var_prc = c(10, 75, 90),
#'    var_fun = "ns",
#'    lag_fun = "ns",
#'    max_lag = 21,
#'    lagnk = 3
#'  )
#'
#'
#' # Set cross-basis parameters
#'  argvar <- list(fun = dlnm_var$var_fun,
#'                 knots = stats::quantile(london$tmean,
#'                                  dlnm_var$var_prc/100, na.rm = TRUE),
#'                 Bound = range(london$tmean, na.rm = TRUE))
#'
#'  arglag <- list(fun = dlnm_var$lag_fun,
#'                 knots = dlnm::logknots(dlnm_var$max_lag, nk = dlnm_var$lagnk))
#'
#'  # Create crossbasis
#'  cb <- dlnm::crossbasis(london$tmean, lag = dlnm_var$max_lag, argvar, arglag)
#'
#'  # Seasonality of mortality time series
#'  seas <- splines::ns(london$date, df = round(8 * length(london$date) / 365.25))
#'
#'  # Prediction values (equidistant points)
#'  temp <- round(seq(min(london$tmean), max(london$tmean), by = 0.1), 1)
#'
#'  # Fit the model
#'  mod <- bdlnm(mort_75plus ~ cb + factor(dow) + seas, basis = cb, data = london,
#'  family = "poisson")
#'
#'  # Find minimum risk exposure value
#'  mmt <- optimal_exposure(mod, cb, at = temp)
#'
#'  # Plot
#'  plot(mmt, xlab = "Temperature (ºC)",
#'  main = paste0("MMT (Median = ", round(mmt$summary[["0.5quant"]], 1), "ºC)"))
#'
#'
plot.optimal_exposure <- function(x, line.arg = NULL, ...) {

  ## ---------------------------
  ## Basic checks
  ## ---------------------------

  if (missing(x) || !inherits(x, "optimal_exposure")) {
    cli::cli_abort("{.arg x} must be an x of class {.cls optimal_exposure} as returned by {.fn optimal_exposure}.")
  }

  ## ---------------------------
  ## Plot
  ## ---------------------------

  # by default set the following arguments
  which <- attr(x, "which")

  plot.arg <- list(
    x = x$est,
    breaks = attr(x, "xvar"),
    xaxt = "n",
    xlab = "Exposure",
    main = paste0(which, "effect values (Median = ", x$summary[["0.5quant"]], ")"),
    xlim = c(floor(min(x$est)), ceiling(max(x$est)))
  )

  # merge with user arguments
  plot.arg <- utils::modifyList(plot.arg, list(...))

  # plot histogram
  do.call("hist", plot.arg)

  # x-axis ticks at the supplied grid points
  graphics::axis(side = 1, at = attr(x, "xvar"))

  # by default set the following arguments
  if(is.null(line.arg)) line.arg <- list()

  line.arg.def <- list(
    col = "red",
    lwd = 3,
    lty = "dashed"
  )

  # merge with user arguments
  line.arg <- utils::modifyList(line.arg.def, line.arg)

  line.arg <- c(list(v = x$summary[["0.5quant"]]), line.arg)

  # draw median line
  do.call(graphics::abline, line.arg)

}
