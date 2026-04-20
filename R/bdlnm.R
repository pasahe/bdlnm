#' Fit a Bayesian distributed lag non-linear model (B-DLNM)
#'
#' Fit a distributed lag non-linear model (DLNM) using a Bayesian framework. The function calls [INLA::inla()] to fit the model and then draws posterior samples of the model fixed effects with [INLA::inla.posterior.sample()]. See the package vignette for worked examples and recommended workflows.
#'
#' @param formula A model formula (as for [INLA::inla()]). The model must be a distributed lag linear or non-linear model (DLNM), so a `"crossbasis"` ([dlnm::crossbasis()]) or a `"onebasis"` ([dlnm::onebasis()]) must be included.
#' @param data an optional data frame, list or environment (or object coercible by [as.data.frame] to a data frame) containing the variables in the model. If not found in data, the variables are taken from `environment(formula)`, typically the environment from which `bdlnm` is called.
#' @param family Character. Family name passed to [INLA::inla()] (default `"gaussian"`).
#' @param sample.arg List of arguments passed to [INLA::inla.posterior.sample()]. Defaults to `list(n = 1000, seed = 0L)` (draws `1000` posterior samples; seed at random). For reproducible sampling set a non-zero numeric `seed`.
#' @param ci.level Numeric in `(0,1)` giving the credible interval level (default `0.95`). Credible interval quantiles are computed to summarize coefficients from the posterior samples.
#' @param na.action A function specifying how to handle NA values when constructing the model frame. The default is taken from the na.action setting of [options], which is by default [na.omit] (drops rows with any `NA` among the variables referenced in formula). In the presence of a random effect term `f()` in the formula, this argument is ignored and `NA`s are not discarded. When rows containing missing values are retained, `INLA` handles them internally (see `Details` below).
#' @param ... Additional arguments passed to [INLA::inla()].
#'
#' @section Distributed lag non-linear model:
#' The fitted model must be a distributed lag linear or non-linear model (DLNM). DLNMs describe potentially non-linear and delayed (lagged) associations between an exposure and an outcome, commonly referred to as exposure–lag–response relationships. This modelling framework is based on the definition of a cross-basis (a bi-dimensional space of functions) constructed with [dlnm::crossbasis()], which defines the exposure–response and lag–response functions simultaneously. The cross-basis object must be created beforehand and supplied as an object in the calling environment (not as a column inside data) and explicitly included in the model formula (e.g., y ~ cb + ...). A basis object constructed with [dlnm::onebasis()] can be used instead when the model is restricted to a uni-dimensional exposure–response relationship (i.e., without lagged effects). All basis objects included in the model `formula` are stored and returned as a named list in the `basis` component. Any of these basis objects can later be supplied to [bcrosspred()] to extract predictions for the corresponding exposure–lag–response (or exposure-response, if created with `onebasis`) association.
#'
#' @section INLA:
#' Models are fit using Integrated Nested Laplace approximation (INLA) via [INLA::inla()]. INLA is a method for approximate Bayesian inference. In the last years it has established itself as an alternative to other methods such as Markov chain Monte Carlo because of its speed and ease of use via the R-INLA package (\href{https://www.r-inla.org/what-is-inla}{What is INLA?}).
#'
#' Additional arguments supplied via `...` are forwarded to [INLA::inla()] (see documentation for all available arguments). Internally, the function ensures that `control.compute = list(config = TRUE)` in order to enable posterior sample drawing with [INLA::inla.posterior.sample()].
#'
#' In the presence of missing values in variables referenced in formula, the `na.action` argument controls how the model frame is constructed. If `na.action` is set to [na.omit] (the default set by [options]), a complete-case analysis is performed and any row with a missing value in a variable appearing in formula is removed. If `na.action` is set to [na.pass] instead, missing values are retained in the model frame and passed to INLA.
#'
#' Note that when the model formula includes a random effect term, specified via `f(k, model = ...)`, the `na.action` is ignored and rows with missing values are not dropped.
#'
#' When missing values are present in the data supplied to [INLA::inla()] (either because a random effect term is included or because `na.action = na.pass`), `INLA` handles them internally as follows:
#'
#' - If `NA` values occur in the response, the corresponding observation contributes nothing to the likelihood (the response is treated as unobserved for that observation).
#' - If `NA` values occur in fixed-effect covariates, [INLA::inla()] replaces them internally with zero so that the covariate does not contribute to the linear predictor for that observation.
#' - If `NA` values occur in a fixed-effect covariate that is a factor, this is not allowed unless `NA` is explicitly included as a level, or `control.fixed = list(expand.factor.strategy = "inla")` is specified. With this option, `NA` is interpreted similarly as in the fixed-effect case, producing no contribution from that covariate to the linear predictor.
#' - If `NA` values occur in a random effect, the random effect does not contribute to the linear predictor for the corresponding observation.
#'
#' @section Posterior samples:
#'
#' After fitting the model, the function draw samples from the approximate posterior distribution of the latent field via [INLA::inla.posterior.sample()]. These samples are collected into a matrix and summarized across samples (mean, sd, quantiles and mode). For a `"crossbasis"` built from an exposure basis with C parameters and a lag basis with L parameters, there will be C × L cross-basis associated coefficients (named e.g. v1.l1, ..., vC.lL). For a `"onebasis"` object the coefficients follow the simpler form b1, ... bC.
#'
#' Additional arguments supplied via `sample.arg` are forwarded to [INLA::inla.posterior.sample()] (see documentation for all available arguments). By default, the number of samples is `1000`. Be aware of the computation and memory cost when increasing the number of samples drawn. By default, the seed is set at random. For reproducible samplings, you need to set a non-zero numeric `seed` in `sample.arg`.
#'
#' Posterior sample estimations are then summarized across samples using mean, sd, credible-interval quantiles (the mid and the lower/upper tails according to `ci.level`) and an approximate mode obtained from a kernel density estimate.
#'
#' @section Requirements:
#' The INLA package must be installed from the R-INLA repository (\href{https://www.r-inla.org/}{R-INLA Project}); if not available the function aborts with a short instruction on how to install it.
#'
#' @return An S3 object of class `"bdlnm"` with the following components:
#' - `model`: the fitted `INLA` model returned by [INLA::inla()].
#' - `basis`: a named list containing all basis of class `crossbasis` or `onebasis` included in the model `formula`.
#' - `coefficients`: a matrix whose columns are posterior sample draws returned by [INLA::inla.posterior.sample()] (named `sample1`, `sample2`, ...) and whose rows are all model coefficients.
#' - `coefficients.summary`: a matrix of summary statistics for all the posterior samples stored in `coefficients` (mean, sd, quantiles, mode).
#'
#' @author Pau Satorra, Marcos Quijal-Zamorano.
#'
#'
#' @references
#'
#' Quijal-Zamorano M., Martinez-Beneito M.A., Ballester J., Marí-Dell'Olmo M. (2024). Spatial Bayesian distributed lag non-linear models (SB-DLNM) for small-area exposure-lag-response epidemiological modelling. _International Journal of Epidemiology_, 53(3), dyae061. <doi:10.1093/ije/dyae061>.
#'
#' Quijal-Zamorano M., Martinez-Beneito M.A., Ballester J., Marí-Dell'Olmo M. (2025). Spatial Bayesian distributed lag non-linear models with R-INLA. _International Journal of Epidemiology_, 54(4), dyaf120. <doi:10.1093/ije/dyaf120>.
#'
#' Gasparrini A. (2011). Distributed lag linear and non-linear models in R: the package dlnm. _Journal of Statistical Software_, 43(8), 1-20. <doi:10.18637/jss.v043.i08>.
#'
#' Rue H., Martino S., Chopin N. (2009). Approximate Bayesian inference for latent Gaussian models by using integrated nested Laplace approximations. _Journal of the Royal Statistical Society: Series B_, 71(2), 319-392. <doi:10.1111/j.1467-9868.2008.00700.x>.
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
#'  # Ensure it falls inside the range of temperatures after rounding:
#'  temp <- temp[temp >= min(london$tmean) & temp <= max(london$tmean)]
#'
#' if (check_inla()) {
#'  # Fit the model
#'  mod <- bdlnm(mort_75plus ~ cb + factor(dow) + seas, data = london, family = "poisson",
#'              sample.arg = list(seed = 432, seed = 1L))
#' }
#'
#'
#'
bdlnm <- function(
  formula,
  data,
  family = "gaussian",
  sample.arg = list(n = 1000, seed = 0L),
  ci.level = 0.95,
  na.action = getOption("na.action"),
  ...
) {
  # ----------------------------
  # Basic checks
  # ----------------------------
  if (is.null(formula)) {
    cli::cli_abort(
      "A {.arg formula} must be provided."
    )
  } else {
    if (!inherits(formula, "formula")) {
      cli::cli_abort(
        "An object of class {.cls formula} must be provided in {.arg formula}."
      )
    }
  }

  if (!is.list(sample.arg)) {
    cli::cli_abort(
      "{.arg sample.arg} must be a {.val list}."
    )
  }

  # detect that formula must contain a crossbasis/onebasis and store them in a list to be returned
  terms <- all.vars(formula)
  envir <- parent.frame()
  basis <- list()
  for (t in terms) {
    obj <- NULL

    # Search first in the dataframe and then in the environment
    if (!is.null(data) && t %in% names(as.data.frame(data))) {
      obj <- data[[t]]
    } else if (exists(t, envir = envir, inherits = TRUE)) {
      obj <- get(t, envir = envir, inherits = TRUE)
    }

    if (!is.null(obj) && inherits(obj, c("crossbasis", "onebasis"))) {
      basis[[t]] <- obj
    }
  }

  if (length(basis) == 0) {
    cli::cli_abort(
      "A basis of class {.cls crossbasis} and/or {.cls onebasis} must be included in {.arg formula}."
    )
  }

  # ----------------------------
  # Build model frame, model design matrix and response
  # ----------------------------

  # In case there is no random effect, we can fit the model using model frame:

  fterms <- terms(formula, specials = "f")

  if (is.null(attr(fterms, "specials")$f)) {
    has_random_terms <- FALSE
    .bdlnm_mframe <- stats::model.frame(
      formula,
      data = data,
      na.action = na.action
    )
    .bdlnm_mmatrix <- stats::model.matrix(formula, data = .bdlnm_mframe)
    .bdlnm_mresponse <- stats::model.response(.bdlnm_mframe)
  } else {
    has_random_terms <- TRUE
    cli::cli_inform(
      c(
        "A random term has been detected in {.arg formula}, so {.arg na.action} will be ignored.",
        "i" = "Missing values will be treated as documented in the help page."
      )
    )
    # add the name of the basis to its coefficient names, recreating model frame
    # new environment in which the formula will be evaluated to reassign the renamed basis
    eval_env <- new.env(parent = envir)
    for (t in names(basis)) {
      base <- basis[[t]]
      colnames(base) <- paste0(t, colnames(base))
      assign(t, base, envir = eval_env)
    }
    environment(formula) <- eval_env
  }

  # ----------------------------
  # Check for INLA availability
  # ----------------------------
  check_inla(error = TRUE)

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

  if (!has_random_terms) {
    inla_options <- c(
      list(
        formula = .bdlnm_mresponse ~ -1 + .bdlnm_mmatrix,
        family = family,
        data = list(
          .bdlnm_mresponse = .bdlnm_mresponse,
          .bdlnm_mmatrix = .bdlnm_mmatrix
        )
      ),
      inla_options
    )
  } else {
    inla_options <- c(
      list(
        formula = formula,
        family = family,
        data = data
      ),
      inla_options
    )
  }

  #Fit model
  model <- tryCatch(
    do.call(INLA::inla, inla_options),
    error = function(e) {
      cli::cli_abort(
        c(
          "Failed to fit INLA model: {.emph {e$message}}",
          "i" = "Check the model formula, family, data, and additional options passed via `...`."
        )
      )
    }
  )

  # ----------------------------
  # Draw posterior samples
  # ----------------------------

  # Validate presence of fixed-effect summaries
  if (is.null(model$summary.fixed) || nrow(model$summary.fixed) == 0L) {
    cli::cli_abort(
      "Fitted {.pkg INLA} model contains no fixed-effect summaries in {.code model$summary.fixed}."
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

  posterior <- tryCatch(
    {
      do.call(INLA::inla.posterior.sample, sample.arg)
    },
    error = function(e) {
      cli::cli_abort(
        c(
          "Failed to draw posterior samples via {.fn INLA::inla.posterior.sample}: {.emph {e$message}}",
          "i" = "Check the additional options passed via {.arg sample.arg}."
        )
      )
    }
  )

  # Combine posterior latent vectors into a matrix (columns = samples)
  coef <- do.call(cbind, lapply(posterior, function(x) x$latent))

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

  res <- list(
    model = model,
    basis = basis,
    coefficients = coef,
    coefficients.summary = coefsum
  )

  attr(res, "n_sim") <- sample.arg$n

  class(res) <- c("bdlnm")

  return(res)
}
