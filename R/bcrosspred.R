#' Predict exposure-lag-response effects from a Bayesian distributed-lag model (B-DLNM).
#'
#' Calculate predictions from a fitted Bayesian distributed lag non-linear model ([bdlnm()]). Predicted associations are computed on a grid of values of the exposure and lags, relative to a reference exposure center value. The function gives posterior samples of exposure–lag-specific associations, overall cumulative associations (summed across lags) and optionally incremental cumulative associations, together with summary statistics (mean, sd, credible interval quantiles and mode).
#'
#' @param object A fitted `"bdlnm"` object returned by [bdlnm()].
#' @param basis If the `bdlnm` model has more than one basis, the name of the basis to use to compute predictions. It must be one of `names(object$basis)`. If the model contains only one basis it is selected automatically.
#' @param exp_at Numeric vector of exposure values at which to evaluate predictions. If `NULL`, the exposure range is extracted from the attributes of `basis` and a grid of 50 values is constructed using [pretty()].
#' @param lag_at Numeric vector of integer lag values at which to evaluate predictions. Only used when `basis` is a `crossbasis`. If `NULL`, the full lag range stored in `basis` is used with a step of `1`.
#' @param cen Centering exposure value for predictions. If `NULL` the centering value depends on the exposure basis function or set to a mid-range value (see Details).
#' @param model.link Optional character specifying the model link (if `NULL` it is inferred from the fitted model).
#' @param ci.level Numeric in `(0,1)` giving the credible interval level (default `0.95`). Credible interval quantiles are computed from the posterior samples.
#' @param cumul Logical; if `TRUE` compute incremental cumulative effects along lags (default `FALSE`). It will give the cumulative effect from lag 0 up to each subsequent lag (e.g., lag `0`, lag `0–1`, lag `0–2`, etc.).
#'
#' @details
#'
#' The function computes predictions for specific combinations of exposure and lag values specified in `exp_at` and `lag_at.` All values must lie within the range defined by the specified `basis.` Note that if the specified `basis` is a `onebasis`, `lag_at` is ignored. If either argument is `NULL`, the grid is derived from the attributes of the specified `basis`: for exposure values, a grid of approximately 50 values is constructed using [pretty()] within the exposure range; for lag values, an integer sequence covering the full lag range with a step of 1 is used.
#'
#' Predictions are computed relative to a centering/reference value (`cen`). If `NULL`, the default `cen` depends on the exposure-response basis function: for `strata`, `thr` and `integer` the reference corresponds to the reference category, and for `lin` the reference is set to 0. For other choices, such as `ns`, `bs`, `poly` or other existing or user-defined functions, the default centering value is set to an approximate mid-range value. For non-linear exposure-response associations is sometimes recommended to manually set the centering value to a data-driven center such as the optimal exposure value (see [optimal_exposure()]).
#'
#' Posterior sample of the predicted associations are stored as matrices for the overall cumulative effect and 3D arrays for the exposure-lag-specific predictions. Summaries across these samples are computed using the mean, sd, credible-interval quantiles (the mid and the lower/upper tails according to `ci.level`) and an approximate mode obtained from a kernel density estimate. Relative risks versions of these associations (exponentiated predictions) are also included if `model.link` is equal to `"log"` or `"logit"`. The `model.link` can be manually specified or, if `NULL`, it is tried to be inferred from the `model` type in `object`.
#'
#' Be aware of memory usage: exposure-lag-specific predictions are stored in 3D arrays of dimension `length(predvar)` × `length(predlag)` × `n_sim`. For dense grids and many posterior samples this can be computationally intensive.
#'
#' This function can also be used to compute predictions for models with simple uni-dimensional basis functions not including lags, if the specified basis is `"onebasis"` instead of `"crossbasis"`. In this case, only unlagged predicted associations are returned.
#'
#'
#' @return An object of class `"bcrosspred"` (a list) with elements including:
#' - `exp_at`: the exposure grid used for prediction as supplied via the `exp_at` argument (vector).
#' - `lag_at`: if basis is `crossbasis`, the lag grid used for prediction as supplied via the `lag_at` argument (vector).
#' - `cen`: the exposure centering value used for prediction (number).
#' - `coefficients`: matrix of posterior coefficient draws (columns = samples).
#' - `coefficients.summary`: matrix of coefficient summaries (mean, sd, quantiles, mode).
#' - `matfit`: 3D array of sampled lag-specific effects (exp × lag × samples).
#' - `matfit.summary`: 3D array of summaries for `matfit` (exp × lag × summary-statistics).
#' - `allfit`: matrix of sampled overall cumulative (summed across lags) effects (exp × samples).
#' - `allfit.summary`: matrix of summaries for `allfit` (exp x summary-statistics).
#' - `cumfit`: (optional) 3D array of sampled incremental cumulative effects (exp × lag × samples).
#' - `cumfit.summary`: (optional) 3D array of summaries for `cumfit` (exp × lag × summary-statistics).
#' - `matRRfit`, `allRRfit`, `matRRfit.summary`, `allRRfit.summary`, `cumRRfit.summary` (optional), `cumRRfit.summary` (optional: relative-risk versions (only when link is `log` or `logit`).
#' - `ci.level`, `model.class`, `model.link`.
#'
#' @author Pau Satorra, Marcos Quijal-Zamorano.
#'
#' @note This function is inspired by [dlnm::crosspred()] developed by Gasparrini (2011) <doi:10.18637/jss.v043.i08>. It has been adapted to work in a Bayesian framework within the \pkg{bdlnm} package.
#'
#' @references
#'
#' Gasparrini A. (2011). Distributed lag linear and non-linear models in R: the package dlnm. _Journal of Statistical Software_, 43(8), 1-20. <doi:10.18637/jss.v043.i08>.
#'
#' Quijal-Zamorano M., Martinez-Beneito M.A., Ballester J., Marí-Dell'Olmo M. (2024). Spatial Bayesian distributed lag non-linear models (SB-DLNM) for small-area exposure-lag-response epidemiological modelling. _International Journal of Epidemiology_, 53(3), dyae061. <doi:10.1093/ije/dyae061>.
#'
#' @seealso [plot.bcrosspred()]  to plot the predicted associations stored in a `"bcrosspred"` object,
#' @seealso [bdlnm()] to fit a Bayesian distributed lag non-linear model (`"bdlnm"`).
#' @seealso [attributable()] to calculate attributable fractions and numbers for a `"bdlnm"` object,
#' @seealso [optimal_exposure()] to estimate exposure values that optimize the predicted effect for a `"bdlnm"` object.
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
#'  # Ensure it falls inside the range of temperatures after rounding:
#'  temp <- temp[temp >= min(london$tmean) & temp <= max(london$tmean)]
#'
#' if (check_inla()) {
#'  # Fit the model
#'  mod <- bdlnm(mort_75plus ~ cb + factor(dow) + seas, data = london, family = "poisson",
#'              sample.arg = list(seed = 432, seed = 1L))
#'
#'  # Prediction
#'  cpred <- bcrosspred(mod, exp_at = temp)
#' }
#'

bcrosspred <- function(
  object,
  basis = NULL,
  exp_at = NULL,
  lag_at = NULL,
  cen = NULL,
  model.link = NULL,
  ci.level = 0.95,
  cumul = FALSE
) {
  ## -----------------------
  ## Basic checks
  ## -----------------------

  # check object
  check_bdlnm(object)

  #Get model basis and coefficients
  model <- object$model

  obj_basis <- object$basis

  if (!is.null(basis)) {
    if (is.character(basis) && length(basis) == 1) {
      if (!basis %in% names(obj_basis)) {
        cli::cli_abort(
          "The name of the basis in {.arg basis} must match the name of the basis stored in {.arg object$basis}: {.or {.val {names(obj_basis)}}}."
        )
      } else {
        basis <- obj_basis[[basis]]
      }
    } else {
      cli::cli_abort(
        "The {.arg basis} argument must be a character element specifying the name of the basis stored in {.arg object$basis}."
      )
    }
  } else {
    if (length(obj_basis) > 1) {
      cli::cli_abort(
        "{.arg basis} must be provided containing the name of the basis to select for predictions (stored in {.arg object$basis}): {.or {.val {names(obj_basis)}}}."
      )
    } else {
      basis <- obj_basis[[1]]
    }
  }

  # Basis type:
  if (inherits(basis, "crossbasis")) {
    type <- "cb"
  } else if (inherits(basis, "onebasis")) {
    type <- "one"
  } else {
    cli::cli_abort(
      "Unsupported {.arg basis} class. Expected {.cls crossbasis} or {.cls onebasis}."
    )
  }

  #Get only CB coefficients
  coef <- extract_coef(object, basis)

  #Get model link
  if (is.null(model.link)) {
    model.link <- get_link(model)
  }

  #Get number of posterior samples
  n_sample <- attr(object, "n_sim")

  # 0 < ci.level < 1
  if (
    !is.numeric(ci.level) ||
      length(ci.level) != 1 ||
      ci.level <= 0 ||
      ci.level >= 1
  ) {
    cli::cli_abort(
      "{.arg ci.level} must be a single numeric value strictly between 0 and 1."
    )
  }

  ## -----------------------
  ## Construct prediction grid & centering
  ## -----------------------

  #Get range of the exposure stored in the attributes of the crossbasis
  exp_range <- attr(basis, "range")
  exp_from <- exp_range[1]
  exp_to <- exp_range[2]

  # Set exposure at if not provided
  if (is.null(exp_at)) {
    exp_pretty <- pretty(c(exp_from, exp_to), n = 50)
    predvar <- exp_pretty[exp_pretty >= exp_from & exp_pretty <= exp_to]
  } else {
    if (!is.numeric(exp_at)) {
      cli::cli_abort("{.arg exp_at} must be a numeric vector.")
    } else {
      # check if exposure is contained in range
      if (any(exp_at > exp_to | exp_at < exp_from)) {
        cli::cli_abort(
          "{.arg exp_at} is outside the range of exposure values defined by the specified {.arg basis}."
        )
      } else {
        # if values are duplicated or not sorted
        predvar <- sort(unique(exp_at))
      }
    }
  }

  lag_range <- switch(type, cb = attr(basis, "lag"), one = c(0, 0))
  lag_from <- lag_range[1]
  lag_to <- lag_range[2]

  if (type == "one" && !is.null(lag_at)) {
    cli::cli_warn(
      "{.arg lag_at} is ignored when {.arg basis} is a {.cls onebasis}."
    )
    lag_at <- NULL
  }

  if (is.null(lag_at)) {
    predlag <- seq(lag_from, lag_to, by = 1)
  } else {
    if (!is.numeric(lag_at) || !all(lag_at == floor(lag_at))) {
      cli::cli_abort(
        "{.arg lag_at} must be a numeric vector containing integer values."
      )
    } else {
      # check if lag values are contained in range
      if (any(lag_at > lag_to | lag_at < lag_from)) {
        cli::cli_abort(
          "{.arg lag_at} is outside the range of lag values defined by the specified {.arg basis}."
        )
      } else {
        # if values are duplicated or not sorted
        predlag <- sort(unique(lag_at))
      }
    }
  }

  #Define centering value
  nocen <- FALSE
  # If NULL try to extract it from basis
  if (is.null(cen)) {
    nocen <- TRUE
    cen <- switch(
      type,
      cb = attributes(basis)$argvar$cen,
      one = attributes(basis)$cen
    )
  }

  fun <- switch(
    type,
    cb = attributes(basis)$argvar$fun,
    one = attributes(basis)$fun
  )

  if (!is.null(fun) && fun %in% c("thr", "strata", "integer", "lin")) {
    if (is.logical(cen)) cen <- NULL
  } else {
    # If we cannot infer it from the function basis, set to median of exposure values
    if (is.null(cen) || (is.logical(cen) && cen)) {
      cen <- stats::median(predvar)
    }
    # If FALSE set to NULL
    if (is.logical(cen) && !cen) {
      cen <- NULL
    }
  }

  # If intercept is present, set to NULL
  int <- switch(
    type,
    cb = attributes(basis)$argvar$intercept,
    one = attributes(basis)$intercept
  )

  if (is.logical(int) && int) {
    cen <- NULL
  }

  if (!is.null(cen) && (!is.numeric(cen) || length(cen) != 1L)) {
    cli::cli_abort("{.arg cen} must be a numeric scalar.")
  }

  if (nocen && !is.null(cen)) {
    cli::cli_warn(
      "Centering value unspecified ({.arg cen}). Automatically set to: {cen}."
    )
  }

  ## -----------------------
  ## Prediction of lag-specific effects
  ## -----------------------

  #Create the matrix of transformed centered exposure-lags values

  varvec <- rep(predvar, length(predlag))
  lagvec <- rep(predlag, each = length(predvar))
  # Reconstruct basis with the specified grid predictor values in exp_at and lag_at
  if (type == "cb") {
    # Marginal basis with onebasis (DLNM) and tensor product
    basisvar <- do.call(
      dlnm::onebasis,
      c(list(x = varvec), attr(basis, "argvar"))
    )
    basislag <- do.call(
      dlnm::onebasis,
      c(list(x = lagvec), attr(basis, "arglag"))
    )
    if (!is.null(cen)) {
      basiscen <- do.call(
        dlnm::onebasis,
        c(list(x = cen), attr(basis, "argvar"))
      )
      basisvar <- scale(basisvar, center = basiscen, scale = FALSE)
    }
    Xpred <- crs::tensor.prod.model.matrix(list(basisvar, basislag))
  } else if (type == "one") {
    #Get the argument names needed for each function
    arg_fun <- switch(
      fun,
      ns = names(formals(splines::ns)),
      bs = names(formals(splines::bs)),
      poly = names(formals(stats::poly)),
      ps = c("x", "df", "knots", "degree", "intercept", "fx", "S", "diff"),
      cr = c("x", "df", "knots", "intercept", "fx", "S"),
      strata = c("x", "df", "breaks", "ref", "intercept"),
      thr = c("x", "thr.value", "side", "intercept"),
      integer = c("x", "values", "intercept"),
      lin = c("x", "intercept"),
      stop("Unsupported function")
    )
    # Call onebasis function (DLNM)
    ind <- match(c("fun", arg_fun), names(attributes(basis)), nomatch = 0)
    basisvar <- do.call(
      dlnm::onebasis,
      c(list(x = varvec), attributes(basis)[ind])
    )
    if (!is.null(cen)) {
      basiscen <- do.call(
        dlnm::onebasis,
        c(list(x = cen), attributes(basis)[ind])
      )
      basisvar <- scale(basisvar, center = basiscen, scale = FALSE)
    }
    Xpred <- basisvar
  }

  # Create the estimated lag-specific effects (each sample in each column)
  matfit <- array(
    Xpred %*% coef,
    dim = c(length(predvar), length(predlag), n_sample),
    dimnames = list(
      predvar,
      paste0("lag", predlag),
      paste0("sample", seq_len(n_sample))
    )
  )

  ## ----------------------
  ## Overall & cumulative
  ## ----------------------

  # Xpred arranged as blocks for each lag: sum over lags for overall effect

  Xpredall <- 0

  cumfit <- matrix(0, nrow(Xpred), n_sample)

  for (i in seq(length(predlag))) {
    ind <- seq(length(predvar)) + length(predvar) * (i - 1)
    Xpredall <- Xpredall + Xpred[ind, , drop = FALSE]

    if (cumul) {
      cumfit[ind, ] <- Xpredall %*% coef
    }
  }

  # convert to array
  cumfit <- array(
    cumfit,
    dim = c(length(predvar), length(predlag), n_sample),
    dimnames = list(
      predvar,
      paste0("lag", predlag),
      paste0("sample", seq_len(n_sample))
    )
  )

  allfit <- Xpredall %*% coef
  rownames(allfit) <- predvar
  colnames(allfit) <- paste0("sample", seq_len(n_sample))

  ## ----------------------
  ## Build output object
  ## ----------------------

  # Initial list then add components
  res <- list(exp_at = predvar)

  if (type == "cb") {
    res$lag_at <- predlag
  }

  if (!is.null(cen)) {
    res$cen <- cen
  }

  colnames(coef) <- paste0("sample", seq_len(n_sample))

  res <- c(
    res,
    list(
      coefficients = coef,
      matfit = matfit,
      allfit = allfit
    )
  )

  if (cumul) {
    res$cumfit <- cumfit
  }

  ## inverse link
  if (!is.na(model.link) && model.link %in% c("log", "logit")) {
    res$matRRfit <- exp(matfit)
    res$allRRfit <- exp(allfit)
    if (cumul) res$cumRRfit <- exp(cumfit)
  }

  ## -----------------------
  ## summaries
  ## -----------------------
  coefsum <- extract_coef(object, basis, "coefficients.summary")

  # recalculate quantiles if another ci.level is provided
  quantiles <- c((1 - ci.level) / 2, 0.5, 1 - (1 - ci.level) / 2)
  quantiles_orig <- as.numeric(gsub(
    "[^0-9.]",
    "",
    grep("quant", colnames(coefsum), value = TRUE)
  ))
  ci_compute <- !identical(quantiles, quantiles_orig)

  if (ci_compute) {
    for (i in seq_along(quantiles)) {
      colnames(coefsum)[
        colnames(coefsum) == paste0(quantiles_orig[i], "quant")
      ] <- paste0(quantiles[i], "quant")
      coefsum[, paste0(quantiles[i], "quant")] <- apply(
        coef,
        1,
        stats::quantile,
        probs = quantiles[i]
      )
    }
  }

  sumcols <- colnames(coefsum)

  res$coefficients.summary <- coefsum

  #matfit summary
  matfitsum <- array(
    dim = c(length(predvar), length(predlag), length(sumcols)),
    dimnames = list(predvar, paste0("lag", predlag), sumcols)
  )

  matfitsum[,, "mean"] <- apply(res$matfit, c(1, 2), mean)
  matfitsum[,, "sd"] <- apply(res$matfit, c(1, 2), stats::sd)

  for (q in quantiles) {
    matfitsum[,, paste0(q, "quant")] <- apply(
      res$matfit,
      c(1, 2),
      stats::quantile,
      probs = q
    )
  }

  #calculate mode (using default kernel density estimate, revisar...)
  matfitsum[,, "mode"] <- apply(res$matfit, c(1, 2), function(v) {
    dv <- stats::density(v)
    with(dv, x[which.max(y)])
  })

  res$matfit.summary <- matfitsum

  if (cumul) {
    #cumfit summary
    cumfitsum <- array(
      dim = c(length(predvar), length(predlag), length(sumcols)),
      dimnames = list(predvar, paste0("lag", predlag), sumcols)
    )

    cumfitsum[,, "mean"] <- apply(res$cumfit, c(1, 2), mean)
    cumfitsum[,, "sd"] <- apply(res$cumfit, c(1, 2), stats::sd)

    for (q in quantiles) {
      cumfitsum[,, paste0(q, "quant")] <- apply(
        res$cumfit,
        c(1, 2),
        stats::quantile,
        probs = q
      )
    }

    #calculate mode (using default kernel density estimate, revisar...)
    cumfitsum[,, "mode"] <- apply(res$cumfit, c(1, 2), function(v) {
      dv <- stats::density(v)
      with(dv, x[which.max(y)])
    })

    res$cumfit.summary <- cumfitsum
  }

  #allfit summary
  allfitsum <- matrix(nrow = nrow(res$allfit), ncol = length(sumcols))
  rownames(allfitsum) <- rownames(res$allfit)
  colnames(allfitsum) <- sumcols

  allfitsum[, "mean"] <- apply(res$allfit, 1, mean)
  allfitsum[, "sd"] <- apply(res$allfit, 1, stats::sd)

  for (q in quantiles) {
    allfitsum[, paste0(q, "quant")] <- apply(
      res$allfit,
      1,
      stats::quantile,
      probs = q
    )
  }

  #calculate mode (using default kernel density estimate, revisar...)
  allfitsum[, "mode"] <- apply(res$allfit, 1, function(v) {
    dv <- stats::density(v)
    with(dv, x[which.max(y)])
  })

  res$allfit.summary <- allfitsum

  # Relative-risk summaries
  if (!is.na(model.link) && model.link %in% c("log", "logit")) {
    res$matRRfit.summary <- exp(res$matfit.summary)
    if (cumul) {
      res$cumRRfit.summary <- exp(res$cumfit.summary)
    }
    res$allRRfit.summary <- exp(res$allfit.summary)
  }

  res$ci.level <- ci.level
  res$model.class <- "inla"
  res$model.link <- model.link

  class(res) <- "bcrosspred"

  return(res)
}
