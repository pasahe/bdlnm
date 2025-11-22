#' Fit a bayesian distributed lag non-linear model (B-DLNM)
#'
#' Fit a distributed lag non-linear model (DLNM) that includes a basis object from `dlnm`. It uses `INLA` to fit the bayesian model and to extract posterior samples of the basis coefficients.
#'
#' @param formula A model formula (as for [INLA::inla]).
#' @param basis A DLNM basis object produced by `dlnm`. It can be one of [dlnm::crossbasis] or [dlnm::onebasis].
#' @param family Family name passed to [INLA::inla] (default `"gaussian"`).
#' @param data A data frame used for the model fitted with [INLA::inla].
#' @param sample.arg List containing the arguments passed to [INLA::inla.posterior.sample]. Default sets `n = 1000` as the number of samples and `seed = 0L` (set the seed at 'random').
#' @param ci.level Numeric in `(0,1)` giving the credible interval level (default `0.95`).
#' @param ... Additional arguments passed to [INLA::inla].
#'
#' @returns An object of class `"bdlnm"` with components:
#' - `model`: The fitted `INLA` model.
#' - `coefficients`: A matrix with posterior samples for the basis coefficients
#'
#' @details
#' If the basis indicates some lag, the first `max(lags)` rows of `data` are set to `NA` to omit the basis first missing rows for model fitting.
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
#'  mod <- bdlnm(mort_75plus ~ cb + factor(dow) + seas, basis = cb, data = london, family = "poisson")
#'
#'
#' @export
#'
bdlnm <- function(formula,
                  basis,
                  data,
                  family = "gaussian",
                  sample.arg = list(n = 1000, seed = 0L),
                  ci.level = 0.95,
                  ...) {

  # ----------------------------
  # Basic checks
  # ----------------------------
  if (missing(data) || is.null(data)) {
    cli::cli_abort("A {.arg data} data.frame must be provided.")
  }

  if (missing(basis)) {
    cli::cli_abort(
      "A basis of class {.cls 'crossbasis'} or {.cls 'onebasis'} must be provided to {.arg basis}."
    )
  }

  # Determine basis type
  type <- if (inherits(basis, "crossbasis")) {
    "cb"
  } else if (inherits(basis, "onebasis")) {
    "one"
  } else {
    cli::cli_abort(
      "Unsupported class for {.arg basis}. Expected {.cls 'crossbasis'} or {.cls 'onebasis'}."
    )
  }

  if (!is.list(sample.arg)) {
    cli::cli_abort(
      "{.arg sample.arg} must be a {.val list}."
    )
  }

  # ----------------------------
  # Remove first NA lags
  # ----------------------------

  # Extract lag attribute
  lag <- switch(type, cb = attr(basis, "lag"), one = c(0, 0))

  # Compute prediction lags: (revisar el by)
  predlag <- seq(from = lag[1], to = lag[2], by = 1)

  # Insert NA's at the data first max(predlag) rows (revisar. Crec que si els lags són negatius s'haurà de treure les files del final...)
  if (any(predlag > 0)) {
    na_rows <- predlag[predlag > 0]
    data[na_rows, ] <- NA
  }

  # ----------------------------
  # Check for INLA availability
  # ----------------------------
  if (!requireNamespace("INLA", quietly = TRUE)) {
    cli::cli_abort(
      c(
        "Package {.pkg INLA} is required to fit the model but is not installed.",
        "i" = "Install from the R-INLA repository (https://www.r-inla.org/) and restart R."
      )
    )
  }

  # ----------------------------
  # Fit INLA model
  # ----------------------------

  # Collect INLA arguments from '...'
  inla_options <- list(...)

  # Ensure control.compute$config = TRUE
  cc <- inla_options[["control.compute"]]

  if (!is.null(cc)) {
    if (!cc$config) {
      cli::cli_abort(
        "{.arg control.compute(config = FALSE)} cannot be provided: {.arg config} must be {.val TRUE} to use {.fn inla.posterior.sample}."
      )
    }
  } else {
    cc <- list(config = TRUE)
  }

  inla_options[["control.compute"]] <- cc

  inla_options <- c(list(
    formula = formula,
    family = family,
    data = data
  ),
  inla_options)

  #Fit model
  model <- tryCatch(
    do.call(INLA::inla, inla_options),
    error = function(e) {
      cli::cli_abort(
        c("Failed to fit INLA model: {.emph {e$message}}", "i" = "Check the model formula, family, data, and additional options passed via `...`.")
      )
    }
  )

  # ----------------------------
  # Draw posterior samples
  # ----------------------------

  # Validate presence of fixed-effect summaries
  if (is.null(model$summary.fixed) || nrow(model$summary.fixed) == 0L) {
    cli::cli_abort(
      "Fitted {.pkg INLA} model contains no fixed-effect summaries in {.code model$summary.fixed}; cannot extract basis coefficients."
    )
  }

  names_sel <- rownames(model$summary.fixed)

  if (is.null(names_sel)) {
    cli::cli_abort(
      "Could not determine fixed-effect names from the fitted {.pkg INLA} model: {.code rownames(model$summary.fixed)} is {.val NULL}."
    )
  }

  # Create selection list for posterior sampling
  list_sel <- as.list(rep(1L, length(names_sel)))
  names(list_sel) <- names_sel

  # set defaults of sampling arguments
  sample.arg.def <- list(
    n = 1000,
    seed = 0L
  )

  # merge with user arguments
  sample.arg <- utils::modifyList(sample.arg.def, sample.arg)

  sample.arg <- c(list(result = model, selection = list_sel), sample.arg)

  posterior <- tryCatch({
    do.call(INLA::inla.posterior.sample, sample.arg)
  }, error = function(e) {
    cli::cli_abort(
      c(
        "Failed to draw posterior samples via {.fn INLA::inla.posterior.sample}: {.emph {e$message}}",
        "i" = "Check the additional options passed via {.arg sample.arg}."
      )
    )
  })

  # Combine posterior latent vectors into a matrix (columns = samples)
  coef <- do.call(cbind, lapply(posterior, function (x) x$latent))

  if (is.null(ncol(coef))) {
    # If only one sample, ensure matrix form
    coef <- matrix(coef, ncol = 1L)
  }
  rownames(coef) <- gsub("\\:1$", "", rownames(coef))
  colnames(coef) <- paste0("sample", seq_len(ncol(coef)))

  ## -----------------------
  ## summaries
  ## -----------------------

  quantiles <- c((1 - ci.level) / 2, 0.5, 1 - (1 - ci.level) / 2)
  sumcols <- c("mean", "sd", paste0(quantiles, "quant"), "mode")

  # sampled coefficients summary
  coefsum <- matrix(nrow = nrow(coef), ncol = length(sumcols))
  rownames(coefsum) <- rownames(coef)
  colnames(coefsum) <- sumcols

  coefsum[, "mean"] <- apply(coef, 1, mean)
  coefsum[, "sd"] <- apply(coef, 1, stats::sd)

  for (q in quantiles) {
    coefsum[, paste0(q, "quant")] <- apply(coef, 1, stats::quantile, probs = q)
  }

  #calculate mode (using default kernel density estimate, revisar...)
  coefsum[, "mode"] <- apply(coef, 1, function(v) {
    dv <- stats::density(v)
    with(dv, x[which.max(y)])
  })

  res <- list(model = model, coefficients = coef, coefficients.summary = coefsum)

  attr(res, "n_sim") <- sample.arg$n

  class(res) <- c("bdlnm")

  return(res)

}
