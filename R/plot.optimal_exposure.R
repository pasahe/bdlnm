#' Plot posterior distribution of optimal effect exposure values
#'
#' It plots a histogram of the posterior distribution of the optimal effect exposure values returned by [optimal_exposure()].
#'
#' @param x An object of class `"optimal_exposure"` returned by [optimal_exposure()].
#' @param show_median Logical. If `TRUE` (default) the function draws a vertical line for the median of all the posterior samples. If `FALSE` it doesn't draw any additional line.
#' @param vline.arg Optional list of graphical arguments passed to [graphics::abline()] when drawing the median vertical line.
#' @param ... Optional graphical parameters passed to [graphics::hist()].
#'
#' @details
#'
#' The histogram uses the original prediction grid in `attr(object, "xvar")` as the x-axis values, ensuring that the bars align with prediction exposure values. The function plots the posterior distribution of the optimal exposure values (stored in `x$est`) and highlights the posterior median across samples (stored in `x$summary[["0.5quant"]]`) with a vertical line if `show_median = TRUE`. Use `vline.arg` to change the appearance of that line passed to [graphics::abline()] and `...` to change the graphical parameters of the histogram passed to [graphics::hist()] (to control axis labels, title, colours, etc.). See the original functions for a complete list of the arguments. Some arguments, if not specified, are set to different default values than the original functions.
#'
#' @author Pau Satorra, Marcos Quijal-Zamorano.
#'
#' @references
#'
#' Quijal-Zamorano M, Martinez-Beneito MA, Ballester J, Marí-Dell’Olmo M. Spatial Bayesian distributed lag non-linear models (SB-DLNM) for small-area exposure-lag-response epidemiological modelling. International Journal of Epidemiology. 2024;53(3):dyae061.
#'
#' Gasparrini A (2011). Distributed lag linear and non-linear models in R: the package dlnm. Journal of Statistical Software, 43(8), 1–20.
#'
#' Armstrong B. Models for the relationship between ambient temperature and daily mortality. Epidemiology. 2006;17(6):624-31.
#'
#' @seealso [optimal_exposure()] to estimate exposure values that optimize the predicted effect for a `"bdlnm"` object.
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
plot.optimal_exposure <- function(x, show_median = TRUE, vline.arg = NULL, ...) {

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

  # if vertical line with median has to be drawn
  if(show_median) {

    # by default set the following arguments
    if(is.null(vline.arg)) vline.arg <- list()

    vline.arg.def <- list(
      col = "red",
      lwd = 3,
      lty = "dashed"
    )

    # merge with user arguments
    vline.arg <- utils::modifyList(vline.arg.def, vline.arg)

    vline.arg <- c(list(v = x$summary[["0.5quant"]]), vline.arg)

    # draw median line
    do.call(graphics::abline, vline.arg)

  }

}
