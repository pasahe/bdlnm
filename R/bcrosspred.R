#' Predict exposure-lag-response effects from a Bayesian distributed-lag models (B-DLNM).
#'
#' Calculate predictions from a fitted Bayesian distributed lag non-linear model ([bdlnm()]). Predicted associations are computed on a grid of values of the exposure and lags, relative to a reference exposure center value. The function gives posterior samples of exposure–lag-specific associations, overall cumulative associations (summed across lags) and optionally incremental cumulative associations, together with summary statistics (mean, sd, credible interval quantiles and mode).
#'
#' @param object A fitted `"bdlnm"` object returned by [bdlnm()].
#' @param basis A DLNM basis object produced by `dlnm`. It must be of class `"crossbasis"` ([dlnm::crossbasis()]) or `"onebasis"` ([dlnm::onebasis()]).
#' @param model.link Optional character specifying the model link (if `NULL` it is inferred from the fitted model).
#' @param at Numeric vector (or matrix) of exposure values at which to compute predictions. If `NULL` the function reconstructs a grid using `from`, `to`, `by` together with the `basis` attributes.
#' @param from,to,by Optional numeric used to construct `at` when not provided.
#' @param lag Optional lag specification; when `NULL` the original lag from `basis` is used. If supplied it must be of length 1 or 2 and represent the lag interval.
#' @param bylag Integer lag step (default `1L`).
#' @param cen Centering exposure value for predictions. If `NULL` the centering value depends on the exposure basis function or set to a mid-range value (see Details).
#' @param ci.level Numeric in `(0,1)` giving the credible interval level (default `0.95`). Credible interval quantiles are computed from the posterior samples.
#' @param cumul Logical; if `TRUE` compute incremental cumulative predictions along lags (default `FALSE`).
#'
#' @details
#'
#' The function computes predictions for specific combinations of the requested exposure and lag grids. The values in `at` can be provided as a vector; in this case, they are replicated for each lag. Alternatively, `at` can be provided as a matrix of complete exposure histories over the same lag period used for estimation to compute the association with a specific exposure pattern.
#'
#' Predictions are computed relative to a centering/reference value (`cen`). If `NULL`, the default `cen` depends on the exposure-response basis function: for `strata`, `thr` and `integer` the reference corresponds to the reference region, and for `lin` the reference is set to 0. For other choices, such as `ns`, `bs`, `poly` or other existing or user-defined functions, the default centering value is set to an approximate mid-range value. For non-linear exposure-response associations is sometimes recommended to manually set the centering value to a data-driven center such as the optimal exposure value (see [optimal_exposure()]).
#'
#' Posterior sample of the predicted associations are stored as matrices for the overall cumulative effect and 3D arrays for the exposure-lag-specific predictions. Summaries across these samples are computed using the mean, sd, credible-interval quantiles (the mid and the lower/upper tails according to `ci.level`) and an approximate mode obtained from a kernel density estimate. Relative risks versions of these associations (exponentiated predictions) are also included if `model.link` is equal to `"log"` or `"logit"`. The `model.link` can be manually specified or, if `NULL`, it is tried to be inferred from the `model` type in `object`.
#'
#' Be aware of memory usage: exposure-lag-specific predictions are stored in 3D arrays of dimension `length(predvar)` × `length(predlag)` × `n_sim`. For dense grids and many posterior samples this can be computationally intensive.
#'
#' This function can also be used to compute predictions for models with simple uni-dimensional basis functions not including lags, if the basis supplied is `"onebasis"` instead of `"crossbasis"`. In this case, only unlagged predicted associations are returned.
#'
#'
#' @return An object of class `"bcrosspred"` (a list) with elements including:
#' - `predvar`: the exposure grid used for prediction (vector).
#' - `lag`: numeric vector of length 2 with the lag interval used.
#' - `bylag`: the lag step used.
#' - `coefficients`: matrix of posterior coefficient draws (columns = samples).
#' - `coefficients.summary`: matrix of coefficient summaries (mean, sd, quantiles, mode).
#' - `matfit`: 3D array of sampled lag-specific effects (predvar × predlag × samples).
#' - `matfit.summary`: 3D array of summaries for `matfit` (predvar × predlag × summary-statistics).
#' - `allfit`: matrix of sampled overall cumulative (summed across lags) effects (predvar × samples).
#' - `allfit.summary`: matrix of summaries for `allfit` (predvar x summary-statistics).
#' - `cumfit`: (optional) 3D array of sampled incremental cumulative effects (predvar × predlag × samples).
#' - `cumfit.summary`: (optional) 3D array of summaries for `cumfit` (predvar × predlag × summary-statistics).
#' - `matRRfit`, `allRRfit`, `matRRfit.summary`, `allRRfit.summary`, `cumRRfit.summary` (optional), `cumRRfit.summary` (optional: relative-risk versions (only when link is `log` or `logit`).
#' - `cen`: centering value used.
#' - `ci.level`, `model.class`, `model.link`.
#'
#' @author Pau Satorra, Marcos Quijal.
#'
#' @note This function is inspired by [dlnm::crosspred()] (Gasparrini 2011). It has been adapted to work in a Bayesian framework within the \pkg{bdlnm} package.
#'
#' @references
#'
#' Gasparrini A. Distributed lag linear and non-linear models in R: the package dlnm. Journal of Statistical Software. 2011; 43(8):1-20.
#'
#' @seealso [plot.bcrosspred()]  to plot the predicted associations stored in a `"bcrosspred"` object,
#' @seealso [bdlnm()] to fit a Bayesian distributed lag non-linear model.
#' @seealso [attributable()] to calculate attributable fractions and numbers for a `bdlnm` object,
#' @seealso [optimal_exposure()] to estimate exposure values that optimise the predicted effect for a `bdlnm` object.
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
#'  # Prediction
#'  cpred <- bcrosspred(mod, cb, at = temp)
#'

bcrosspred <- function(object, basis, model.link = NULL, at = NULL, from = NULL, to = NULL, by = NULL, lag = NULL, bylag = 1L, cen = NULL, ci.level = 0.95, cumul = FALSE) {

  ## -----------------------
  ## Basic checks
  ## -----------------------

  check_bdlnm(object)

  #Get model and coefficients
  model <- object$model

  if (missing(basis)) {
    cli::cli_abort(
      "A basis of class {.cls 'crossbasis'} or {.cls 'onebasis'} must be provided to {.arg basis}."
    )
  }

  # Basis type:
  if (inherits(basis, "crossbasis")) {
    type <- "cb"
  } else if (inherits(basis, "onebasis")) {
    type <- "one"
  } else {
    cli::cli_abort(
      "Unsupported {.arg basis} class. Expected {.cls 'crossbasis'} or {.cls 'onebasis'}."
    )
  }

  #Get only CB coefficients
  coef <- extract_coef(object$coefficients, basis)

  #Get model link
  if (is.null(model.link)) model.link <- get_link(model)

  #Get number of posterior samples
  n_sample <- attr(object, "n_sim")


  #Get default lag
  origlag <- switch(type, cb  = attr(basis, "lag"), one = c(0, 0))

  # if missing we use original lag; if provided, attempt to convert
  if (is.null(lag)) {
    lag <- origlag
  } else {
    #  Lag must be a postive integer vector
    if (any(!is.numeric(lag)) || length(lag) > 2) {
      cli::cli_abort("{.arg lag} must be a integer vector of length {.val 2} or {.val 1}")
    }
    if (length(lag) == 1L) {
      lag <- if (lag < 0L) c(lag, 0L)
    } else {
      c(0L, lag)
    }

    if (diff(lag) < 0L) cli::cli_abort("{.code lag[1]} must be lower or equal than {.code lag[2]}")

    lag <- round(lag[1L:2L])
  }

  # cumul not allowed if a subperiod of lags is specified
  if (!identical(origlag, lag) && cumul) {
    cli::cli_abort(
      "{.arg cumul = TRUE} is not allowed when {.arg lag} differs from the original lag specified in the basis attribute."
    )
  }

  lagfun <- switch(type, cb = attr(basis, "arglag")$fun, one = NULL)
  # check on lagfun 'integer'
  if (bylag != 1L && !is.null(lagfun) && identical(lagfun, "integer")) {
    cli::cli_abort(
      "{.arg bylag} != 1 is not allowed when the lag function for {.arg basis} is {.val integer}."
    )
  }

  # 0 < ci.level < 1
  if (!is.numeric(ci.level) || length(ci.level) != 1 || ci.level <= 0 || ci.level >= 1) {
    cli::cli_abort("{.arg ci.level} must be a single numeric value strictly between 0 and 1.")
  }


  ## -----------------------
  ## Construct prediction grid (at / predvar / predlag) & centering
  ## -----------------------

  # RANGE
  #Get range of the exposure stored in the attributes of the crossbasis
  range <- attr(basis, "range")

  # Set at if not provided
  if (is.null(at)) {
    if (is.null(from)) from <- range[1]
    if (is.null(to)) to <- range[2]
    nobs <- ifelse(is.null(by), 50, max(1, diff(range) / by))
    pretty <- pretty(c(from, to), n = nobs)
    pretty <- pretty[pretty >= from & pretty <= to]
    at <- if (is.null(by)) pretty
    else {
      seq(from = min(pretty),
          to = to,
          by = by)
    }
  } else {
    if (!is.numeric(at)) {
      cli::cli_abort("{.arg at} must be an integer vector or matrix.")
    } else {
      if (is.matrix(at)) {
        if (dim(at)[2] != diff(lag) + 1L) cli::cli_abort("matrix in {.arg at} must have {.val {diff(lag)+1}} columns.")
        if (bylag != 1) cli::cli_abort("{.arg bylag} cannot be different from {val 1}, if {.arg at} is in matrix form")
        if (is.null(rownames(at))) rownames(at) <- seq(nrow(at))
      } else {
        at <- sort(unique(at))
      }
    }
  }

  # define the matrix of temperatures and lags in which predictions will be made
  predvar <- if (is.matrix(at)) rownames(at) else at

  predlag <- seq(from = lag[1], to = lag[2], by = 1)

  #Define centering value
  nocen <- FALSE
  # If NULL try to extract it from basis
  if (is.null(cen)) {
    nocen <- TRUE
    cen <- switch(type,
                  cb = attributes(basis)$argvar$cen,
                  one = attributes(basis)$cen)

  }

  fun <- switch(type,
                cb = attributes(basis)$argvar$fun,
                one = attributes(basis)$fun)

  if (!is.null(fun) && fun %in% c("thr", "strata", "integer", "lin")) {
    if (is.logical(cen)) cen <- NULL
  } else {
    # If we cannot infer it from the function basis, set to mid-range (approximately)
    if (is.null(cen) || (is.logical(cen) && cen)) {
      cen <- stats::median(pretty(range))
    }
    # If FALSE set to NULL
    if (is.logical(cen) && !cen)
      cen <- NULL
  }

  # If intercept is present, set to NULL
  int <- switch(type,
                cb = attributes(basis)$argvar$intercept,
                one = attributes(basis)$intercept)

  if (is.logical(int) && int) cen <- NULL

  if (!is.numeric(cen) || length(cen) != 1L) {
    cli::cli_abort("{.arg cen} must be a numeric scalar.")
  }

  if (nocen && !is.null(cen))
    cli::cli_warn("Centering value unspecified ({.arg cen}). Automatically set to: {cen}.")

  ## -----------------------
  ## Prediction of lag-specific effects
  ## -----------------------

  #Create the matrix of transformed centered exposure-lags values

  varvec <- if (is.matrix(at)) as.numeric(at) else rep(at, length(predlag))
  lagvec <- rep(predlag, each = length(predvar))
  #
  if (type == "cb") {
    # Marginal basis with onebasis (DLNM) and tensor product
    basisvar <- do.call(dlnm::onebasis, c(list(x = varvec), attr(basis, "argvar")))
    basislag <- do.call(dlnm::onebasis, c(list(x = lagvec), attr(basis, "arglag")))
    if (!is.null(cen)) {
      basiscen <- do.call(dlnm::onebasis, c(list(x = cen), attr(basis, "argvar")))
      basisvar <- scale(basisvar, center = basiscen, scale = FALSE)
    }
    Xpred <- crs::tensor.prod.model.matrix(list(basisvar, basislag))
  } else if (type == "one") {
    # Call onebasis function (DLNM)
    ind <- match(c("fun", names(formals(
      attr(basis, "fun")
    ))), names(attributes(basis)), nomatch = 0)
    basisvar <- do.call(dlnm::onebasis, c(list(x = varvec), attributes(basis)[ind]))
    if (!is.null(cen)) {
      basiscen <- do.call(dlnm::onebasis, c(list(x = cen), attributes(basis)[ind]))
      basisvar <- scale(basisvar, center = basiscen, scale = FALSE)
    }
    Xpred <- basisvar
  }

  # Create the estimated lag-specific effects (each sample in each column)
  matfit <- array(Xpred %*% coef,
                  dim = c(length(predvar), length(predlag), n_sample),
                  dimnames = list(predvar,
                                  paste0("lag", predlag),
                                  paste0("sample", seq_len(n_sample))))

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
  cumfit <- array(cumfit,
                  dim = c(length(predvar), length(predlag), n_sample),
                  dimnames = list(predvar,
                                  paste0("lag", predlag),
                                  paste0("sample", seq_len(n_sample))))

  allfit <- Xpredall %*% coef
  rownames(allfit) <- predvar
  colnames(allfit) <- paste0("sample", seq_len(n_sample))

  ## ----------------------
  ## Build output object
  ## ----------------------

  # Initial list then add components
  res <- list(predvar = predvar)
  if (!is.null(cen)) res$cen <- cen

  colnames(coef) <- paste0("sample", seq_len(n_sample))

  res <- c(res,
           list(
             lag = lag,
             bylag = bylag,
             coefficients = coef,
             matfit = matfit,
             allfit = allfit
           ))

  if (cumul) res$cumfit <- cumfit

  ## inverse link
  if (!is.na(model.link) && model.link %in% c("log", "logit")) {
    res$matRRfit <- exp(matfit)
    res$allRRfit <- exp(allfit)
    if(cumul) res$cumRRfit <- exp(cumfit)
  }


  ## -----------------------
  ## summaries
  ## -----------------------
  coefsum <- extract_coef(object$coefficients.summary, basis)

  # recalculate quantiles if another ci.level is provided
  quantiles <- c((1 - ci.level) / 2, 0.5, 1 - (1 - ci.level) / 2)
  quantiles_orig <- as.numeric(gsub("[^0-9.]", "", grep("quant", colnames(coefsum), value = TRUE)))
  ci_compute <- !identical(quantiles, quantiles_orig)

  if(ci_compute) {
    for (i in seq_along(quantiles)) {
      colnames(coefsum)[colnames(coefsum) == paste0(quantiles_orig[i], "quant")] <- paste0(quantiles[i], "quant")
      coefsum[, paste0(quantiles[i], "quant")] <- apply(coef, 1, stats::quantile, probs = quantiles[i])
    }
  }

  sumcols <- colnames(coefsum)

  res$coefficients.summary <- coefsum

  #matfit summary
  matfitsum <- array(dim = c(length(predvar), length(predlag), length(sumcols)), dimnames = list(predvar, paste0("lag", predlag), sumcols))

  matfitsum[, , "mean"] <- apply(res$matfit, c(1, 2), mean)
  matfitsum[, , "sd"] <- apply(res$matfit, c(1, 2), stats::sd)

  for (q in quantiles) {
    matfitsum[, , paste0(q, "quant")] <- apply(res$matfit, c(1, 2), stats::quantile, probs = q)
  }

  #calculate mode (using default kernel density estimate, revisar...)
  matfitsum[, , "mode"] <- apply(res$matfit, c(1, 2), function(v) {
    dv <- stats::density(v)
    with(dv, x[which.max(y)])
  })

  res$matfit.summary <- matfitsum

  if(cumul) {
    #cumfit summary
    cumfitsum <- array(dim = c(length(predvar), length(predlag), length(sumcols)), dimnames = list(predvar, paste0("lag", predlag), sumcols))

    cumfitsum[, , "mean"] <- apply(res$cumfit, c(1, 2), mean)
    cumfitsum[, , "sd"] <- apply(res$cumfit, c(1, 2), stats::sd)

    for (q in quantiles) {
      cumfitsum[, , paste0(q, "quant")] <- apply(res$cumfit, c(1, 2), stats::quantile, probs = q)
    }

    #calculate mode (using default kernel density estimate, revisar...)
    cumfitsum[, , "mode"] <- apply(res$cumfit, c(1, 2), function(v) {
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
    allfitsum[, paste0(q, "quant")] <- apply(res$allfit, 1, stats::quantile, probs = q)
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
    if(cumul) res$cumRRfit.summary <- exp(res$cumfit.summary)
    res$allRRfit.summary <- exp(res$allfit.summary)
  }

  res$ci.level <- ci.level
  res$model.class <- "inla"
  res$model.link <- model.link

  class(res) <- "bcrosspred"

  return(res)

}
