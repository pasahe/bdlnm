#' Calculate attributable number and fractions from Bayesian distributed-lag models (B-DLNM).
#'
#' It computes the attributable number (AN) and fraction (AF) from a fitted B-DLNM returned by `bdlnm()`.
#'
#' @param object A fitted `"bdlnm"` class object returned by [bdlnm].
#' @param basis A DLNM crossbasis object produced by `dlnm`. It has to be of class [dlnm::crossbasis].
#' @param data A dataframe containing the time series of exposure values (see `name_exposure`) and number of cases (see `name_cases`). Ensure that the time series are provided on a regular basis for attributable measures to be calculated properly. It can also include a column with the measured time point (see `name_date`) or a column with a `0/1` indicator that filters the calculation of attributable measures for specific time periods (see `name_filter`).
#' @param name_date A character with the name of the column with the date (class `Date` or `POSICXt`) of the time series measurement (optional).
#' @param name_exposure A character with the name of the column with the time series exposure values.
#' @param name_cases A character with the name of the column with the time series case values. If not provided, only attributable fractions per time point can be calculated.
#' @param name_filter A character with the name of the column with the indicator that can filter the time points in which to calculate attributable measures (optional).
#' @param tot Logical; if TRUE (default) returns total attributable number / fraction across the time series. If FALSE returns values for each time point.
#' @param dir Character; "back" (default) or "forw" direction in the lag dimension for calculating attributable numbers and fractions.
#' @param cen Numeric scalar; centering value for predictions. If missing the function will attempt to read it from attr(basis, "argvar")$cen.
#' @param range Optional numeric vector of range 2. It gives the exposure value range for which attributable numbers and fractions will be calculated.
#'
#' @return A list with elements:
#'   - `af`: attributable fraction (AF)
#'   - `an`: attributable number (AN)
#'   - `af.summary`: data frame with summary statistics for AF
#'   - `an.summary`: data frame with summary statistics for AN
#'
#'
#' @export
#'
#' @examples
#'
#' # Filter the dataset to reduce computational time:
#'
#' # Exposure-response and lag-response spline parameters
#' dlnm_var <- list(
#'   var_prc = c(10, 75, 90),
#'   var_fun = "ns",
#'   lag_fun = "ns",
#'   max_lag = 21,
#'   lagnk = 3
#' )
#'
#' # Cross-basis parameters
#' argvar <- list(fun = dlnm_var$var_fun,
#'                knots = stats::quantile(london$tmean,
#'                                 dlnm_var$var_prc/100, na.rm = TRUE),
#'                Bound = range(london$tmean, na.rm = TRUE))
#'
#' arglag <- list(fun = dlnm_var$lag_fun,
#'                knots = dlnm::logknots(dlnm_var$max_lag, nk = dlnm_var$lagnk))
#'
#' # Create crossbasis
#' cb <- dlnm::crossbasis(london$tmean, lag = dlnm_var$max_lag, argvar, arglag)
#'
#' # Seasonality of mortality time series
#' seas <- splines::ns(london$date, df = round(8 * length(london$date) / 365.25))
#'
#' # Prediction values (equidistant points)
#' temp <- round(seq(min(london$tmean), max(london$tmean), by = 0.1), 1)
#'
#' # Model
#'
#' mod <- bdlnm(mort_75plus ~ cb + factor(dow) + seas,
#'             basis = cb,
#'              data = london,
#'              family = "poisson",
#'              sample.arg = list(seed = 432))
#'
#' # Predict
#' cpred <- bcrosspred(mod, cb, at = temp)
#'
#' # compute centering (MMT) using optimal_exposure
#' mmt <- optimal_exposure(mod, cb, at = temp)
#' cen <- mmt$summary[["0.5quant"]]
#'
#' # Attributable numbers and fractions (using the backwards algorithm):
#' ar <- attributable(mod, cb, london, name_date = "date",
#' name_exposure = "tmean", name_cases = "mort_75plus", cen = cen, dir = "back")
#'
#'
#'
attributable <- function(object, basis, data, name_date = NULL, name_exposure, name_cases = NULL, name_filter = NULL, tot = TRUE, dir = "back", cen, range = NULL) {

  ## -----------------------
  ## Basic checks
  ## -----------------------

  # object
  check_bdlnm(object)

  # basis
  if (missing(basis) || !inherits(basis, "crossbasis")) {
    cli::cli_abort("{.arg basis} must be an object of class {.cls crossbasis}.")
  }

  # data
  if (missing(data)) {
    cli::cli_abort("{.arg data} must be provided: a data.frame containing the time index and temporal series of exposures and (optionally) cases.")
  }

  if (!inherits(data, "data.frame")) {
    cli::cli_abort("{.arg data} must be a {.cls data.frame}.")
  }

  # name_exposure (required)
  if (missing(name_exposure) || !is.character(name_exposure) || length(name_exposure) != 1L) {
    cli::cli_abort("{.arg name_exposure} must be a single string with the name of the exposure column in {.arg data}.")
  }

  if (!name_exposure %in% colnames(data)) {
    cli::cli_abort("Exposure column {.val {name_exposure}} not found in {.arg data}.")
  }

  exp <- data[[name_exposure]]

  # name_date (optional)
  if (!is.null(name_date)) {
    if (!is.character(name_date) || length(name_date) != 1L) {
      cli::cli_abort("{.arg name_date} must be a single string naming the date column in {.arg data}.")
    }
    if (!name_date %in% colnames(data)) {
      cli::cli_abort("Date column {.val {name_date}} not found in {.arg data}.")
    }

    date <- data[[name_date]]

  } else {

    cli::cli_warn(c("Ensure that {.arg data} contains time series measured on a regular basis (e.g., every day) to properly calculate attributable measures.",
                    "i" = "If you have only seasonal observations (e.g., summers only), expand your data to the full sequence inserting NA for missing exposures/cases."))

    date <- NULL
  }

  # name_cases (optional)
  if (!is.null(name_cases)) {
    if (!is.character(name_cases) || length(name_cases) != 1L) {
      cli::cli_abort("{.arg name_cases} must be a single string naming the cases column in {.arg data}.")
    }
    if (!name_cases %in% colnames(data)) {
      cli::cli_abort("Cases column {.val {name_cases}} not found in {.arg data}.")
    }
    cases <- data[[name_cases]]
  } else {
    cli::cli_warn("Only attributable fractions per time point can be calculated as {.arg name_cases} is not provided.")
    cases <- NULL
    tot <- FALSE
  }

  # name_filter (optional)
  if (!is.null(name_filter)) {
    if (!is.character(name_filter) || length(name_filter) != 1L) {
      cli::cli_abort("{.arg name_filter} must be a single string naming the filter column in {.arg data}.")
    }
    if (!name_filter %in% colnames(data)) {
      cli::cli_abort("Filter column {.val {name_filter}} not found in {.arg data}.")
    }
    filter <- data[[name_filter]]
    cli::cli_warn("Attributable fractions and numbers will only be calculated for time points filtered by {.val {name_filter}}")
  } else {
    filter <- NULL
  }

  # validate exposure and cases types
  if (!is.numeric(exp)) {
    cli::cli_abort("Exposure column {.val {name_exposure}} must be numeric.")
  }
  if (!is.null(cases) && !is.numeric(cases)) {
    cli::cli_abort("Cases column {.val {name_cases}} must be numeric.")
  }

  # validate date type
  if (!is.null(date) && !inherits(date, "Date") && !inherits(date, "POSIXt")) {
    cli::cli_abort("Date column {.val {name_date}} must be of class {.cls Date} or {.cls POSIXt}.")
  }

  # validate filter
  if(!is.null(filter) && (!is.numeric(filter) | !all(filter %in% c(0, 1)))) {
    cli::cli_abort("Filter column {.val {name_filter}} must contain only {.val 0} and {.val 1}.")
  }

  if(!is.null(filter) && all(filter == 0)) {
    cli::cli_abort("Filter column  {.val {name_filter}} must contain some row equal to {.val 1}, corresponding to the filtered time points.")
  }

  # Here we require differences to be constant in time. Time can be on a seconds, minutes, hourly, daily, weekly, monthly, yearly basis.
  if (!is.null(date)) {

    # Check if it's sorted
    if(any(diff(date) < 1)) {
      cli::cli_abort("The provided date ({.val {name_date}}) is not ordered in time.")
    }

    #Let's check regularity with different units
    is_regular <- length(unique(diff(date, units = "secs"))) == 1 ||
      length(unique(diff(date, units = "mins"))) == 1 ||
      length(unique(diff(date, units = "hours"))) == 1 ||
      length(unique(diff(date, units = "days"))) == 1 ||
      length(unique(diff(date, units = "weeks"))) == 1

    #Let's check if it's on a monthly or yearly basis:
    month_basis <- seq(date[1L], date[length(date)], by = "month")
    is_month_basis <- length(month_basis) == length(date) && all(month_basis == date)

    year_basis <- seq(date[1L], date[length(date)], by = "year")
    is_year_basis <- length(year_basis) == length(date) && all(year_basis == date)

    if(!is_regular && !is_month_basis && !is_year_basis) {
      cli::cli_abort(c(
        "The provided date ({.val {name_date}}) is not provided on a regular basis: time steps between dates are not equal.",
        "i" = "If you have only seasonal observations (e.g., summers only), expand your data to the full sequence inserting NA for missing exposures/cases."
      ))
    }

  }

  # direction
  if (! dir %in% c("back", "forw")) {
    cli::cli_abort("{.arg dir} must be one of: {.val 'back'}, {.val 'forw'}.")
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
  n_sample <- attr(object, "n_sim")

  # get lags
  lag <- attr(basis,"lag")
  predlag <- seq(from = lag[1], to = lag[2], by = 1)

  # obtain centered predictions at the exposure values
  cpred <- tryCatch(
    suppressMessages(bcrosspred(object, basis, at = exp, cen = cen)),
    error = function(e) cli::cli_abort("Failed computing predictions with {.fn bcrosspred}: {conditionMessage(e)}")
  )

  cp <- cpred$matfit

  if (!is.null(cpred$matRRfit)) {
    cp_rr <- cpred$matRRfit
  } else {
    # if there is no RR available in the model we can not calculate attributable numbers
    cli::cli_abort("An attributable fraction cannot be computed because predicted effects are not on the relative-risk scale. Ensure the link of the {.arg object} model is {.val 'log'} or {.val 'logit'}.")

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

    if(!is.null(cases)) {
      #Compute the Lagged matrix of daily cases for that day and the next max_lag days
      lagged_cases <- tsModel::Lag(cases, seq(-lag[1], -lag[2]))
    }

    for(i in seq_len(n_sample)) {

      # filter for i-th sample
      rr_sample <- cp_rr[match(exp, cpred$predvar),,i]
      af_sample <- (rr_sample - 1) / rr_sample

      if(!is.null(cases)) {
        # multiply element-wise by lagged cases
        an_sample <- af_sample * lagged_cases
        # sum across lags to get AN per time point
        M_an[,i] <- rowSums(an_sample)
      }

      # total if requested
      if(tot) {
        ind <- if(is.null(filter)) !is.na(M_an[, i]) else !is.na(M_an[, i]) & (filter == 1)

        if(sum(ind) == length(cases)) cli::cli_abort(c(
          "All time points are excluded by the filter or have a missing attributable number.",
          "i" = "Remember that if the number of cases is missing in some of the lags of a time point, an attributable number for that time point will be missing.")
        )

        sum_cases <- sum(rowMeans(lagged_cases)[ind])

        if(sum_cases == 0) cli::cli_abort(c(
          "The number of cases aggregated in the filtered time points with a non-missing attributable number is zero.",
          "i" = "Remember that if the number of cases is missing in some of the lags of a time point, an attributable number for that time point will be missing.")
        )


        af[1L, i] <- sum(M_an[ind, i])/sum_cases
        an[1L, i] <- af[1L, i] * sum(cases, na.rm = TRUE)
      } else {
        an[, i] <- M_an[,i]
        af[, i] <- exp(rowSums(log(rr_sample)))
        # insert missings if time point is not selected
        if(!is.null(filter)) an[filter == 0, ] <- NA
      }
    }

  # backward perspective: contributions from past exposures to current day
  } else {

    #Compute the Lagged matrix of daily exposures for that day and the previous max_lag days
    back_lagged_temp <- tsModel::Lag(exp, seq(lag[1], lag[2]))
    colnames(back_lagged_temp) <- paste0("lag", seq(lag[1], lag[2]))

    for(i in seq_len(n_sample)) {

      # calculate the matrix of backward rr: get the rr associated to each temperature and lag specified in the column of back_lagged_temp
      back_lagged_rr <- sapply(colnames(back_lagged_temp), function(col_lag) {
        rr <- cp_rr[match(back_lagged_temp[,col_lag], cpred$predvar), col_lag, i]
        return(rr)
      })
      back_rr <- exp(rowSums(log(back_lagged_rr)))

      # filter for the i-th sample
      af_sample <- (back_rr - 1) / back_rr

      M_af[,i] <- af_sample

      if(!is.null(cases)) M_an[,i] <- M_af[,i] * cases

      # total if requested
      if(tot) {
        ind <- if(is.null(filter)) !is.na(M_an[, i]) else !is.na(M_an[, i]) & (filter == 1)

        if(sum(ind) == length(cases)) cli::cli_abort(c(
          "All time points are excluded by the filter or have a missing attributable number.",
          "i" = "Remember that if the number of cases is missing in some of the lags of a time point, an attributable number for that time point will be missing.")
        )

        sum_cases <- sum(cases[ind])

        if(sum_cases == 0) cli::cli_abort(c(
          "The number of cases aggregated in the filtered time points with a non-missing attributable number is zero.",
          "i" = "Remember that if the number of cases is missing in some of the lags of a time point, an attributable number for that time point will be missing.")
        )

        af[1L, i] <- sum(M_an[ind, i])/sum_cases
        an[1L, i] <- af[1L, i] * sum(cases, na.rm = TRUE)
      } else {
        an[, i] <- M_an[,i]
        af[, i] <- M_af[,i]
        # insert missings if time point is not selected
        if(!is.null(filter)) an[filter == 0, ] <- NA
      }

    }

  }

  #Remove rows with NA (will be in the beginning or in the end depending on the algorithm type, and also if a filter is specified)
  if(!tot) {
    if(!is.null(date)) {
      rownames(an) <- rownames(af) <- as.character(date)
    } else {
      rownames(an) <- rownames(af) <- paste0("time", seq_along(exp))
    }
    colnames(an) <- colnames(af) <- paste0("sample", seq_len(n_sample))

    if(!is.null(cases)) {
      ind <- !rowSums(is.na(an))
    } else {
      ind <- !rowSums(is.na(af))
    }

    an <- an[ind,]
    af <- af[ind,]

  }

  ## -----------------------
  ## summaries
  ## -----------------------

  ansum <- afsum <-  matrix(nrow = nrow(af), ncol = ncol(cpred$allfit.summary))
  rownames(ansum) <- rownames(afsum) <- rownames(an)
  colnames(ansum) <- colnames(afsum) <- colnames(cpred$allfit.summary)

  quant_cols <- grep("quant$", colnames(cpred$allfit.summary), value = TRUE)
  quant_val <- as.numeric(gsub("quant$", "", quant_cols))

  if(!is.null(cases)) {
    ansum[,"mean"] <- apply(an, 1, mean)
    ansum[,"sd"] <- apply(an, 1, stats::sd)

    for(i in seq_along(quant_cols)) {
      ansum[, quant_cols[i]] <- apply(an, 1, function (x) stats::quantile(x, quant_val[i]))
    }

    #calculate mode (using default kernel density estimate, revisar...)
    ansum[, "mode"] <- apply(an, 1, function(v) {
      dv <- stats::density(v)
      with(dv, x[which.max(y)])
    })
  }

  afsum[, "mode"] <- apply(af, 1, function(v) {
    dv <- stats::density(v)
    with(dv, x[which.max(y)])
  })

  afsum[, "mean"] <- apply(af, 1, mean)
  afsum[, "sd"] <- apply(af, 1, stats::sd)

  for(i in seq_along(quant_cols)) {
    afsum[, quant_cols[i]] <- apply(af, 1, function (x) stats::quantile(x, quant_val[i]))
  }

  # an and af has to be a row-vector (tot=TRUE)
  if(tot) {
    an <- an[1L,]
    af <- af[1L,]
    names(an) <- names(af) <- paste0("sample", seq_len(n_sample))
  }

  # Return results as data.frames for summaries
  if(!is.null(cases)) {
    res <- list(
      af = af,
      an = an,
      af.summary = as.data.frame(afsum, stringsAsFactors = FALSE),
      an.summary = as.data.frame(ansum, stringsAsFactors = FALSE)
    )
  } else {
    res <- list(
      af = af,
      af.summary = as.data.frame(afsum, stringsAsFactors = FALSE)
    )
  }


  return(res)

}
