#' Optimal effect exposure value estimation
#'
#' The function finds the exposure value that optimizes the overall effect for each posterior sample, together with a small summary of these values. The user can decide to choose the optimal exposure as the value that has either the minimum or the maximum estimated effect. In case of finding the minimum exposure effect in the context of the relationship between temperature and mortality, this value is called Minimum Mortality Temperature (MMT).
#'
#' @param object A fitted `"bdlnm"` class object returned by [bdlnm].
#' @param basis A DLNM basis object produced by `dlnm`. It can be one of [dlnm::crossbasis] or [dlnm::onebasis].
#' @param at Values (or matrix) of the predictor at which to predict; can be `NULL` and reconstructed from `from`, `to`, `by` and the basis attributes.
#' @param from,to,by Optional numeric used to construct `at` when not provided.
#' @param which Character string specifying if the optimal exposure value is chosen as the one that minimizes the effect or as the one that maximizes it. By default is `"min"`.
#'
#' @returns A list of class `optimal_exposure` containing:
#'  - `est`: numeric vector with the optimal exposure value for each posterior sample
#'  - `summary`: data frame with summary statistics for these optimal values
#'
#'  The returned object has attribute `xvar` containing the exposure grid used.
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
optimal_exposure <- function(object, basis, at = NULL, from = NULL, to = NULL, by = NULL, which = "min") {

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

  ## ---------------------------
  ## Predict using bcrosspred()
  ## ---------------------------

  # define the matrix of temperatures and lags in which predictions will be made
  predvar <- if (is.matrix(at)) rownames(at) else at

  predlag <- seq(from = lag[1], to = lag[2], by = 1)

  # prediction
  cpred <- tryCatch({
    suppressWarnings(bcrosspred(object, basis, at = at))
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
