#' Calculate attributable risks from Bayesian distributed-lag models (B-DLNM).
#'
#' It computes the attributable number (AN) and fraction (AF) from a fitted B-DLNM returned by `bdlnm()`.
#'
#' @param x A fitted object returned by [bdlnm] (list with components `model` and `coef`).
#' @param basis A DLNM crossbasis object produced by `dlnm`. It has to be of class [dlnm::crossbasis].
#' @param exp Numeric vector of exposure values (one per observed time point) at which to compute attributable risk.
#' @param cases Numeric vector of counts of the modelled outcome (one per observed time point). Has to be of the same length of `exp`.
#' @param tot Logical; if TRUE (default) returns total attributable number / fraction across the time series. If FALSE returns values for each time point.
#' @param dir Character; "back" (default) or "forw" direction in the lag dimension for calculating attributable risks.
#' @param average Logical. When TRUE the function uses average (lag-averaged) contributions to calculate attributable risks. For forward perspective, average defaults to FALSE.
#' @param cen Numeric scalar; centering value for predictions. If missing the function will attempt to read it from attr(basis, "argvar")$cen.
#' @param range Optional numeric vector of range 2. It gives the exposure value range for which attributable risks will be calculated.
#'
#' @return A list with elements:
#'   - `af`: attributable fraction (AF)
#'   - `an`: attributable number (AN)
#'   - `af.summary`: data frame with summary statistics for AF
#'   - `an.summary`: data frame with summary statistics for AN
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
#'  # Center at Minimum Mortality Temperature (MMT)
#'  mmt <- minimum_risk(mod, cb, at = temp)
#'  cen <- mmt$min.summary[["0.5quant"]]
#'
#'  # Attributables risk (using the backwards algorithm):
#'  ar <- attributable_risk(mod, cb, london$tmean, london$mort_75plus, cen = cen, dir = "back")
#'
#'
#'
attributable_risk <- function(x, basis, exp, cases, tot = TRUE, dir = "back", average = FALSE, cen, range = NULL) {
# TODO: I WOULD RETHINK THE NAME OF THE FUNCTION, I THINK attributable / attributable_measures or attributable_impacts
# TODO: I don't like "x" for the model, i will change it to "model" or "coef" if only the coefficients are necessary

  # TODO 2: We should rewrite the function to add the option type = "af" / "an"
  # (many time you are only interested in the af and in that case for examples is not necessary to provide the cases)

  # is more suitable
  ## -----------------------
  ## Basic checks
  ## -----------------------

  # x
  check_bdlnm(x)

  # basis
  if (missing(basis) || !inherits(basis, "crossbasis")) {
    cli::cli_abort("{.arg basis} must be an object of class {.cls 'crossbasis'}.")
  }

  # exp
  if (missing(exp) || !is.numeric(exp)) {
    cli::cli_abort("{.arg exp} must be a numeric vector of exposure values (one per observed time point).")
  }

  # cases
  if (missing(cases) || !is.numeric(cases)) {
    cli::cli_abort("{.arg cases} must be a numeric vector of counts (one per observed time point).")
  }
  if (length(cases) != length(exp)) {
    cli::cli_abort("Length mismatch: {.arg exp} and {.arg cases} must have the same length.")
  }

  # direction
  if (! dir %in% c("back", "forw")) {
    cli::cli_abort("{.arg cases} must be one of: {.val 'back'}, {.val 'forw'}.")
  }

  # average
  if (!is.logical(average)) {
    cli::cli_abort("{.arg average} must be logical ({.val TRUE}/{.val FALSE}).")
  }

  # define centering value

  # get range of the exposure stored in the attributes of the crossbasis
  range_val <- attr(basis, "range")

  # if NULL try to extract it from basis
  if (is.null(cen)) {
    cen <- attributes(basis)$argvar$cen
    if(is.null(cen)) cli::cli_abort("Centering value must be specified ({.arg cen}) or specified as an attribute of the basis.")
  }

  if (!is.numeric(cen) || length(cen) != 1L) {
    cli::cli_abort("{.arg cen} must be a numeric scalar.")
  }


  ## -----------------------
  ## Prepare
  ## -----------------------

  # get number of samples
  n_sample <- ncol(x$coef)

  # get lags
  lag <- attr(basis,"lag")
  predlag <- seq(from = lag[1], to = lag[2], by = 1)

  # obtain centered predictions at the exposure values
  cpred <- tryCatch(
    suppressMessages(bcrosspred(x, basis, at = exp, cen = cen)),
    error = function(e) cli::cli_abort("Failed computing predictions with {.fn bcrosspred}: {conditionMessage(e)}")
  )

  cp <- cpred$matfit

  if (!is.null(cpred$matRRfit)) {
    cp_rr <- cpred$matRRfit
  } else {
    # if there is no RR available in the model we can not calculate attributable numbers
    cli::cli_abort("An attributable risk cannot be computed because predicted effects are not on the relative-risk scale. Ensure the link of the {.arg x} model is {.val 'log'} or {.val 'logit'}.")

  }

  ## -----------------------
  ## Calculate AF/AN
  ## -----------------------

  # initialize matrices
  M_an <- M_af <- matrix(nrow = length(exp), ncol = n_sample)

  if(!tot) {
    an <- af <- M_an
  } else {
    an <- af <- matrix(nrow = 1L, ncol = n_sample)
  }

  # forward perspective: contributions from the current day to future days
  if(dir == "forw") {

    #Compute the Lagged matrix of daily cases for that day and the next max_lag days
    # TODO: WARMING, WE ARE ASSUMING THAT THE TIME SERIES OF TEMPEARTURES IS
    # CONTINUOS, IF IS NOT THE USER NEED TO PROVIDE THE LAG CASES
    # TODO: maybe we should think of an idea to deal with the issue when the
    # time series is not only continuos (e.g. only summers)
    lagged_cases <- tsModel::Lag(cases, seq(-lag[1], -lag[2]))

    # if(!average) {

      for(i in seq_len(n_sample)) {

        # filter for i-th sample
        rr_sample <- cp_rr[match(exp, cpred$predvar),,i]
        af_sample <- (rr_sample - 1) / rr_sample

        # multiply element-wise by lagged cases
        an_sample <- af_sample * lagged_cases

        # sum across lags to get AN per time point (ignoring NA to include rows in the end)
        # TODO: for me makes more sense to leave an NA in the case there are some NAs in the lagged_cases
        M_an[,i] <- rowSums(an_sample)

        # total if requested
        if(tot) {
          isna <- is.na(M_an[, i])
          an[1L, i] <- sum(M_an[!isna, i])
          af[1L, i] <- if(sum(cases[!isna]) > 0) an[1L, i]/sum(rowMeans(lagged_cases)[!isna]) else NA
          # an <- af*den # revisar: no entenc això del denominador ajustat
        } else {
          an[, i] <- M_an[,i]
          af[, i] <- exp(rowSums(log(rr_sample)))
        }
      }

    # } else {
    #
    #   for(i in seq_len(n_sample)) {
    #
    #     # filter for i-th sample
    #     cp_sample <- cp[, , i]
    #     # sum across lags then exponentiate
    #     cp_rr_all <- exp(rowSums(cp_sample))
    #     M_af[,i] <- (cp_rr_all - 1) / cp_rr_all
    #     # AN per time point using average of lagged cases (mean across lag columns)
    #     M_an[,i] <- M_af[,i] * rowMeans(lagged_cases, na.rm = TRUE)
    #
    #     # total if requested
    #     if(tot) {
    #       isna <- is.na(M_an[,i])
    #       an[1L, i] <- sum(M_an[!isna, i], na.rm = TRUE)
    #       af[1L, i] <- if(sum(cases[!isna]) > 0) an[,i]/sum(cases[!isna]) else NA
    #       # an <- af*den # revisar: no entenc això del denominador ajustat
    #     } else {
    #       an <- M_an
    #       af <- M_af
    #     }
    #   }
    #
    # }

  # backward perspective: contributions from past exposures to current day
  } else {

    # TODO: WARMING, WE ARE ASSUMING THAT THE TIME SERIES OF TEMPEARTURES IS
    # CONTINUOS, IF IS NOT THE USER NEED TO PROVIDE THE BACK_LAGGED TEMPERATURES
    # TODO: ADD back_lagged_temp as a parameter of the functions (default = NULL)
    # calculate the matrix of backward temperatures
    back_lagged_temp <- NULL
    if(is.null(back_lagged_temp)) {
      back_lagged_temp <- tsModel::Lag(exp, seq(lag[1], lag[2]))
      colnames(back_lagged_temp) <- paste0("lag", seq(lag[1], lag[2]))
    }

    for(i in seq_len(n_sample)) {

      # calculate the matrix of backward rr
      back_lagged_rr <- sapply(seq(lag[1], lag[2]), function(i_lag) {
        col_lag <- paste0("lag", i_lag)
        rr <- cpred$matRRfit[
          match(back_lagged_temp[,col_lag], cpred$predvar), col_lag, i]
        return(rr)
      })
      back_rr <- exp(rowSums(log(back_lagged_rr)))

      # filter for the i-th sample
      af_sample <- (back_rr - 1) / back_rr

      # easy to visualize if we do it in a 3x3 example. It's clear that we don't want those row-column elements that sum more than the number of rows
      M_af[,i] <- af_sample
      M_an[,i] <- M_af[,i] * cases

      # total if requested
      if(tot) {
        isna <- is.na(M_an[, i])
        an[1L, i] <- sum(M_an[!isna, i], na.rm = TRUE)
        af[1L, i] <- if(sum(cases[!isna]) > 0) an[,i]/sum(cases[!isna]) else NA
        # an <- af*den #No entenc això del denominador ajustat
      } else {
        an <- M_an
        af <- M_af
      }

    }

  }

  ## -----------------------
  ## summaries
  ## -----------------------

  # TODO: Rethink if we want to give the summary statistics of the attributable
  # measures. I don't see it that relevant (any user can do it easily) and for
  # me is more clean the code and output without that. But maybe i'm bias
  # and people will find it useful.

  ansum <- afsum <-  matrix(nrow = nrow(an), ncol = ncol(cpred$allfit.summary))
  colnames(ansum) <- colnames(afsum) <- colnames(cpred$allfit.summary)

  ansum[,"mean"] <- apply(an, 1, mean)
  ansum[,"sd"] <- apply(an, 1, stats::sd)

  quant_cols <- grep("quant$", colnames(cpred$allfit.summary), value = TRUE)
  quant_val <- as.numeric(gsub("quant$", "", quant_cols))

  for(i in seq_along(quant_cols)) {
    ansum[, quant_cols[i]] <- apply(an, 1, function (x) stats::quantile(x, quant_val[i]))
  }

  #calculate mode (using default kernel density estimate, revisar...)
  ansum[, "mode"] <- apply(an, 1, function(v) {
    dv <- stats::density(v)
    with(dv, x[which.max(y)])
  })
  afsum[, "mode"] <- apply(af, 1, function(v) {
    dv <- stats::density(v)
    with(dv, x[which.max(y)])
  })

  afsum[, "mean"] <- apply(af, 1, mean)
  afsum[, "sd"] <- apply(af, 1, stats::sd)

  for(i in seq_along(quant_cols)) {
    afsum[, quant_cols[i]] <- apply(af, 1, function (x) stats::quantile(x, quant_val[i]))
  }

  # an and af has to be a matrix (tot=FALSE) or row-vector (tot=TRUE)
  if(tot) {
    an <- an[1L,]
    af <- af[1L,]
    names(an) <- names(af) <- paste0("sample", seq_len(n_sample))
  } else {
    rownames(an) <- rownames(af) <- rownames(ansum) <- rownames(afsum) <- paste0("day", seq_along(exp))
    colnames(an) <- colnames(af) <- paste0("sample", seq_len(n_sample))
  }

  # Return results as data.frames for summaries
  res <- list(
    af = af,
    an = an,
    af.summary = as.data.frame(afsum, stringsAsFactors = FALSE),
    an.summary = as.data.frame(ansum, stringsAsFactors = FALSE)
  )

  return(res)

}
