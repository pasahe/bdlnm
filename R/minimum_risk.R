#' Minimum-risk exposure value estimation
#'
#' The function finds the exposure value that minimises the overall effect for each posterior sample, together with a small summary of these minima. This exposure value is called Minimum Mortality Temperature (MMT) in case of the relationship between temperature and mortality.
#'
#' @param x A fitted object returned by [bdlnm] (list with components `model` and `coef`).
#' @param basis A DLNM basis object produced by `dlnm`. It can be one of [dlnm::crossbasis] or [dlnm::onebasis].
#' @param at Values (or matrix) of the predictor at which to predict; can be `NULL` and reconstructed from `from`, `to`, `by` and the basis attributes.
#' @param from,to,by Optional numeric used to construct `at` when not provided.
#'
#' @returns A list of class `min.risk` containing:
#'  - `min`: numeric vector with the minimum-exposure value for each posterior sample
#'  - `min.summary`: data frame with summary statistics for the minima
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
#'  temp <- seq(round(min(london$tmean), 1), round(max(london$tmean), 1), by = 0.1)
#'
#'  # Fit the model
#'  mod <- bdlnm(mort_75plus ~ cb + factor(dow) + seas, basis = cb, data = london, family = "poisson")
#'
#'  # Find minimum risk exposure value
#'  mmt <- minimum_risk(x, cb, at = temp)
#'
minimum_risk <- function(x, basis, at = NULL, from = NULL, to = NULL, by = NULL) {

  ## ---------------------------
  ## Basic checks
  ## ---------------------------

  # check object
  check_bdlnm(x)

  # check basis
  if(missing(basis) || !inherits(basis, "crossbasis")) {
    cli::cli_abort("The {.arg basis} argument must be an object of class {.cls crossbasis}.")
  }

  # extract basis attributes
  attr <- attributes(basis)
  range <- attr(basis,"range")
  lag <- attr(basis,"lag")

  # determine number of posterior samples
  n_sample <- ncol(x$coef)

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

  ## ---------------------------
  ## Predict using bcrosspred()
  ## ---------------------------

  # define the matrix of temperatures and lags in which predictions will be made
  predvar <- if (is.matrix(at))
    rownames(at)
  else
    at

  predlag <- seq(from = lag[1], to = lag[2], by = 1)

  # prediction
  cpred <- tryCatch({
    suppressWarnings(bcrosspred(x, basis, at = at))
  }, error = function(e) {
    cli::cli_abort("Failed to compute predictions via bcrosspred: {conditionMessage(e)}")
  })

  ## ---------------------------
  ## Find the minimum
  ## ---------------------------

  # each column of allfit corresponds to a posterior sample; find index of min
  min_index <- apply(cpred$allfit, 2, which.min)
  min_values <- predvar[min_index]

  names(min_values) <- paste0("sample", seq_len(n_sample))

  ## -----------------------
  ## summaries
  ## -----------------------

  minsum <- matrix(nrow = 1, ncol = ncol(cpred$allfit.summary))
  colnames(minsum) <- colnames(cpred$allfit.summary)

  minsum[,"mean"] <- mean(min_values)
  minsum[,"sd"] <- stats::sd(min_values)

  quant_cols <- grep("quant$", colnames(cpred$allfit.summary), value = TRUE)
  quant_val <- as.numeric(gsub("quant$", "", quant_cols))

  for(i in seq_along(quant_cols)) {
    minsum[, quant_cols[i]] <- stats::quantile(min_values, quant_val[i])
  }

  # approximate mode as most frequent observed value among samples
  minsum[, "mode"] <- unique(min_values)[which.max(tabulate(match(min_values, unique(min_values))))]

  #
  res <- list(min = min_values, min.summary = as.data.frame(minsum))

  attr(res, "xvar") <- predvar

  #Define a class in order to do a plot
  class(res) <- "min.risk"

  return(res)
}
