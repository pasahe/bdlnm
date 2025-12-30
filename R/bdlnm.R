#' Fit a Bayesian distributed lag non-linear model (B-DLNM)
#'
#' Fit a distributed lag non-linear model (DLNM) using a Bayesian framework. The function calls [INLA::inla()] to fit the model and then draws posterior samples of the model latent field with [INLA::inla.posterior.sample()]. See the package vignette for worked examples and recommended workflows.
#'
#' @param formula A model formula (as for [INLA::inla()]).
#' @param basis A DLNM basis object produced by `dlnm`. It must be of class `"crossbasis"` ([dlnm::crossbasis()]) or `"onebasis"` ([dlnm::onebasis()]).
#' @param data A data frame containing the variables referenced in `formula` used by [INLA::inla()].
#' @param family Character. Family name passed to [INLA::inla()] (default `"gaussian"`).
#' @param sample.arg List of arguments passed to [INLA::inla.posterior.sample()]. Defaults to `list(n = 1000, seed = 0L)` (draws `1000` posterior samples; seed at random). For reproducible sampling set a non-zero numeric `seed`.
#' @param ci.level Numeric in `(0,1)` giving the credible interval level (default `0.95`). Credible interval quantiles are computed from the posterior samples.
#' @param ... Additional arguments passed to [INLA::inla()].
#'
#' @section Distributed lag non-linear model:
#' Distributed lag non-linear models (DLNMs) describe simultaneous non-linear and delayed (lagged) dependencies, commonly called exposure-lag-response associations. This modelling framework is based on the definition of a cross-basis (a bi-dimensional space of functions) built with [dlnm::crossbasis()] that encondes the dependency along the space of the predictor (exposure-response) and along lags (lag-response). This cross-basis matrix has to be supplied in the `basis` argument and included in `formula`. A basis object built with [dlnm::onebasis()] can be used instead, when we want to simplify the modelling of a uni-dimensional exposure-response relationship.
#'
#' @section INLA:
#' Models are fit using Integrated Nested Laplace approximation (INLA) via [INLA::inla()]. INLA is a method for approximate Bayesian inference. In the last years it has established itself as an alternative to other methods such as Markov chain Monte Carlo because of its speed and ease of use via the R-INLA package (\href{https://www.r-inla.org/what-is-inla}{What is INLA?}).
#'
#' Additional arguments supplied via `...` are forwarded to [INLA::inla()] (see documentation for all available arguments). Internally, the function ensures that `control.compute = list(config = TRUE)` in order to enable posterior sample drawing with [INLA::inla.posterior.sample()].
#'
#' @section Posterior samples:
#'
#' After fitting the model, the function draw samples from the approximate posterior distribution of the latent field via [INLA::inla.posterior.sample()]. These samples are collected into a matrix and summarized across samples (mean, sd, quantiles and mode). For a cross-basis built from an exposure basis with C parameters and a lag basis with L parameters, there will be C × L cross-basis associated coefficients (named e.g. v1.l1, ..., vC.lL). For a `"onebasis"` object the coefficients follow the simpler form b1, ... bC.
#'
#' Additional arguments supplied via `sample.arg` are forwarded to [INLA::inla.posterior.sample()] (see documentation for all available arguments). By default, the number of samples is `1000`, be aware of the computation and memory cost when increasing the number of samples drawn. By default, the seed is set at random. For reproducible samplings, you need to set a non-zero numeric `seed` in `sample.arg`.
#'
#' Posterior sample estimations are then summarized across samples using mean, sd, credible-interval quantiles (the mid and the lower/upper tails according to `ci.level`) and an approximate mode obtained from a kernel density estimate.
#'
#' @section Requirements:
#' The INLA package must be installed from the R-INLA repository (\href{https://www.r-inla.org/}{R-INLA Project}); if not available the function aborts with a short instruction on how to install it.
#'
#' @return An S3 object of class `"bdlnm"` with the following components:
#' - `model`: the fitted `INLA` model returned by [INLA::inla()].
#' - `coefficients`: a matrix whose columns are posterior sample draws returned by [INLA::inla.posterior.sample()] (named `sample1`, `sample2`, ...) and whose rows are all model coefficients.
#' - `coefficients.summary`: a matrix of summary statistics for all the posterior samples stored in `coefficients` (mean, sd, quantiles, mode).
#'
#' @author Pau Satorra, Marcos Quijal.
#'
#'
#' @references
#'
#' Gasparrini A. Distributed lag linear and non-linear models in R: the package dlnm. Journal of Statistical Software. 2011; 43(8):1-20.
#'
#' Havard Rue, Sara Martino, and Nicholas Chopin (2009), Approximate Bayesian Inference for Latent Gaussian Models Using Integrated Nested Laplace Approximations (with discussion), Journal of the Royal Statistical Society B, 71, 319-392.
#'
#' @seealso [bcrosspred()] to predict exposure–lag–response associations for a `bdlnm` object,
#' @seealso [attributable()] to calculate attributable fractions and numbers for a `bdlnm` object,
#' @seealso [optimal_exposure()] to estimate exposure values that optimise the predicted effect for a `bdlnm` object.
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
#'  mod <- bdlnm(mort_75plus ~ cb + factor(dow) + seas, basis = cb, data = london, family = "poisson")
#'
#'
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
  # Remove observations when onebasis/crossbasis is NA
  # ----------------------------

  # Insert NA's at the data first max(predlag) rows
  if (type == "cb") {
    na_rows <- which(apply(basis, 1, function(x) all(is.na(x))))
  } else if (type == "one") {
    na_rows <- which(is.na(basis))
  }
  data[na_rows,] <- NA

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
