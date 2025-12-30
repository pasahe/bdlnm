#' Calculate the exposure value that minimizes or maximizes the effect of a Bayesian distributed lag non-linear model (B-DLNM)
#'
#' Find exposure values that optimize the overall effect for each posterior sample drawn from a Bayesian distributed lag non-linear model ([bdlnm()]). The function returns the exposure value that minimizes or maximizes the overall cumulative effect (summed across lags) for each posterior sample, together with summary statistics (mean, sd, credible-interval quantiles and mode). When used to find the minimum effect in temperature–mortality analyses this optimal exposure value is commonly called the Minimum Mortality Temperature (MMT).
#'
#' @param object A fitted `"bdlnm"` object returned by [bdlnm()].
#' @param basis A DLNM basis object produced by `dlnm`. It must be of class `"crossbasis"` ([dlnm::crossbasis()]) or `"onebasis"` ([dlnm::onebasis()]).
#' @param at Numeric vector (or matrix) of exposure values at which to compute predictions. If `NULL` the function reconstructs a grid using `from`, `to`, `by` together with the `basis` attributes.
#' @param from,to,by Optional numeric used to construct `at` when not provided.
#' @param which Selection criterion to calculate the optimal exposure: `"min"` (default) chooses the exposure with the minimum overall cumulative effect, `"max"` chooses the exposure with maximum overall cumulative effect.
#' @param ci.level Numeric in `(0,1)` giving the credible-interval level (default `0.95`). Credible interval quantiles are computed from the posterior samples.

#' @details
#'
#' The function internally calls [bcrosspred()] to compute the posterior distribution of the overall exposure effect for the grid specified by `at` (or reconstructed using `from`, `to`, `by` and the attributes of `basis`). For each posterior sample the function calculates the exposure value that optimizes (minimizes or maximizes) the overall cumulative effect and then summarizes these optimal values across samples using mean, sd, credible-interval quantiles and the mode (most frequent observed value).
#'
#' If `basis` is a `crossbasis`, the function works on the overall cumulative effect of each exposure summed across all the lags, stored by [bcrosspred()] in `$allfit`. If `basis` is a `onebasis` instead, then the function uses the exposure effect stored in `$matfit`.
#'
#' In the presence of a non-linear association between exposure and response, this optimal exposure value can be used as the reference exposure value for estimating effects. Therefore, it can be passed to the [bcrosspred()] and [attributable()] functions as the centre exposure value. However, be aware that, in the presence of uncertainty, the optimal exposure range across all samples can be wide, so providing one summary statistic (e.g. the median) as the center reference value can be misleading. It is recommended to visualize the distribution of these optimal exposure values using [plot.optimal_exposure()] before using an optimal exposure value as the center.
#'
#' @return An S3 object of class `"optimal_exposure"` containing:
#'  - `est`: numeric vector with the optimal exposure value for each posterior sample (named sample1, sample2, ...).
#'  - `summary`: a one-row data frame with summary statistics for the optimal values across all samples (mean, sd, quantiles, mode).
#'
#' @author Pau Satorra, Marcos Quijal.
#'
#' @references
#'
#' Gasparrini A (2011). Distributed lag linear and non-linear models in R: the package dlnm. Journal of Statistical Software, 43(8), 1–20.
#'
#' Armstrong B. Models for the relationship between ambient temperature and daily mortality. Epidemiology. 2006;17(6):624-31.
#'
#' @seealso [plot.optimal_exposure()] to plot the optimal exposure values stored in a `"optimal_exposure"` object.
#' @seealso [bcrosspred()] to predict exposure–lag–response associations for a `"bdlnm"` object,
#' @seealso [bdlnm()] to fit a Bayesian distributed lag non-linear model (`"bdlnm"`).
#' @seealso [attributable()] to calculate attributable fractions and numbers for a `"bdlnm"` object.
#'
#' @export
#'
#' @examples
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
#'  mod <- bdlnm(mort_75plus ~ cb + factor(dow) + seas, basis = cb, data = london, family = "poisson")
#'
#'  # Find minimum risk exposure value
#'  mmt <- optimal_exposure(mod, cb, at = temp)
#'
#'
optimal_exposure <- function(object, basis, at = NULL, from = NULL, to = NULL, by = NULL, which = "min", ci.level = 0.95) {

  ## ---------------------------
  ## Basic checks
  ## ---------------------------

  # check object
  check_bdlnm(object)

  # check basis
  if(missing(basis) || !inherits(basis, "crossbasis")) {
    cli::cli_abort("The {.arg basis} argument must be an object of class {.cls crossbasis}.")
  }

  # extract basis attributes
  attr <- attributes(basis)
  range <- attr(basis,"range")
  lag <- attr(basis,"lag")

  # determine number of posterior samples
  n_sample <- attr(object, "n_sim")

  # Set at if not provided
  if (is.null(at)) {
    if (is.null(from))
      from <- range[1]
    if (is.null(to))
      to <- range[2]
    nobs <- ifelse(is.null(by), 50, max(1, diff(range) / by))
    pretty <- pretty(c(from, to), n = nobs)
    pretty <- pretty[pretty >= from & pretty <= to]
    at <- if (is.null(by))
      pretty
    else
      seq(from = min(pretty),
          to = to,
          by = by)
  } else {
    if (!is.numeric(at)) {
      cli::cli_abort("{.arg at} must be an integer vector or matrix.")
    } else {
      if (is.matrix(at)) {
        if (dim(at)[2] != diff(lag) + 1L)
          cli::cli_abort("matrix in {.arg at} must have {.val {diff(lag)+1}} columns.")
        if (is.null(rownames(at)))
          rownames(at) <- seq(nrow(at))
      } else {
        at <- sort(unique(at))
      }
    }
  }

  if(! which %in% c("min", "max")){
    cli::cli_abort("{.arg which} has to be either {.val min} or {.val max}.")
  }

  # 0 < ci.level < 1
  if (!is.numeric(ci.level) || length(ci.level) != 1 || ci.level <= 0 || ci.level >= 1) {
    cli::cli_abort("{.arg ci.level} must be a single numeric value strictly between 0 and 1.")
  }

  ## ---------------------------
  ## Predict using bcrosspred()
  ## ---------------------------

  # define the matrix of temperatures and lags in which predictions will be made
  predvar <- if (is.matrix(at)) rownames(at) else at

  predlag <- seq(from = lag[1], to = lag[2], by = 1)

  # prediction
  cpred <- tryCatch({
    suppressWarnings(bcrosspred(object, basis, at = at, ci.level = ci.level))
  }, error = function(e) {
    cli::cli_abort("Failed to compute predictions via bcrosspred: {conditionMessage(e)}")
  })

  ## ---------------------------
  ## Find the optimal
  ## ---------------------------

  # each column of allfit corresponds to a posterior sample; find index of optimal
  if(which == "min") {
    which.fun <- which.min
  } else {
    which.fun <- which.max
  }

  opt_index <- apply(cpred$allfit, 2, which.fun)
  opt_values <- predvar[opt_index]

  names(opt_values) <- paste0("sample", seq_len(n_sample))

  ## -----------------------
  ## summaries
  ## -----------------------

  optsum <- numeric(ncol(cpred$allfit.summary))
  names(optsum) <- colnames(cpred$allfit.summary)

  optsum["mean"] <- mean(opt_values)
  optsum["sd"] <- stats::sd(opt_values)

  quant_cols <- grep("quant$", colnames(cpred$allfit.summary), value = TRUE)
  quant_val <- as.numeric(gsub("quant$", "", quant_cols))

  for(i in seq_along(quant_cols)) {
    optsum[quant_cols[i]] <- stats::quantile(opt_values, quant_val[i])
  }

  # approximate mode as most frequent observed value among samples
  optsum["mode"] <- unique(opt_values)[which.max(tabulate(match(opt_values, unique(opt_values))))]

  #
  res <- list(est = opt_values, summary = optsum)

  attr(res, "xvar") <- predvar

  attr(res, "which") <- which

  #Define a class in order to do a plot method
  class(res) <- "optimal_exposure"

  return(res)
}
