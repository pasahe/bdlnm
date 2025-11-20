#' Generate predicted effects from Bayesian distributed-lag models (B-DLNM).
#'
#' It produces lag-specific, overall and (optionally) cumulative predictions from an object produced by `bdlnm()`.
#'
#' @param x A fitted object returned by [bdlnm] (list with components `model` and `coef`).
#' @param basis A DLNM basis object produced by `dlnm`. It can be one of [dlnm::crossbasis] or [dlnm::onebasis].
#' @param model.link Optional character; model link (if `NULL` it is
#'   inferred from the fitted model).
#' @param at Values (or matrix) of the predictor at which to predict; can be `NULL` and reconstructed from `from`, `to`, `by` and the basis attributes.
#' @param from,to,by Optional numeric used to construct `at` when
#'   not provided.
#' @param lag Optional lag specification; when `NULL` the original lag
#'   from `basis` is used.
#' @param bylag Integer step for lags (default `1`).
#' @param cen Centering value for predictions; if `NULL` a centering
#'   value is constructed from `basis` attributes.
#' @param ci.level Numeric in `(0,1)` giving the credible interval level (default `0.95`).
#' @param cumul Logical; if `TRUE` compute cumulative predictions (default `FALSE`).
#'
#' @returns An object of class `"bcrosspred"` containing
#'   coefficients, fitted matrices, summaries and relative-risks (if link is
#'   log/logit).
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

bcrosspred <- function(x,
                       basis,
                       model.link = NULL,
                       at = NULL,
                       from = NULL,
                       to = NULL,
                       by = NULL,
                       lag = NULL,
                       bylag = 1L,
                       cen = NULL,
                       ci.level = 0.95,
                       cumul = FALSE) {

  ## -----------------------
  ## Basic checks
  ## -----------------------

  check_bdlnm(x)

  #Get model and coefficients
  model <- x$model

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
  coef <- extract_coef(x$coefficients, basis)

  #Get model link
  if (is.null(model.link)) model.link <- get_link(model)

  #Get number of posterior samples
  n_sample <- attr(x, "n_sim")


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
    # If we cannot infer it from the basis, set to mid-range (approximately)
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

  #Revisar aquest tros (amb l'exemple, per entendre-ho bé)
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
  coefsum <- extract_coef(x$coefficients.summary, basis)

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
