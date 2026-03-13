#' Calculate attributable number and fractions from a Bayesian distributed-lag model (B-DLNM).
#'
#' Compute attributable numbers (AN) and attributable fractions (AF) from a fitted Bayesian distributed lag non-linear model ([bdlnm()]). The function uses posterior predicted relative risks from [bcrosspred()] and applies a forward or backward lag algorithm to compute per-time and total attributable measures, optionally filtering the calculation to a subset of time points.
#'
#' @param object A fitted `"bdlnm"` class object returned by [bdlnm()].
#' @param data A data frame containing the temporal series needed to calculate attributable measures. It can include a date column (optional, see `name_date`), the exposure values (mandatory, see `name_exposure`), the number of cases (optional, see `name_cases`) and a binary (`0/1`) filter column (optional, see `name_filter`).
#' @param name_date Optional single string with the name of the column in `data` containing the date. The column must be of class `Date` or `POSIXt`. When provided the function checks that the series is regularly spaced (see `Details`).
#' @param name_exposure Single string with the name of the exposure column in `data`.
#' @param name_cases Optional single string with the name of the cases column in `data`. If not provided, the function returns only attributable fractions (AF) per time point.
#' @param name_filter Optional single string with the name of a binary (`0/1`) column in `data`. Only rows with value `1` are used to compute total AF and AN and per-time results are returned only for the filtered rows.
#' @param dir Character; direction of the algorithm to calculate attributable measures. `"back"` (default) calculates AF and AN attributing current-time outcome to past exposures; `"forw"` calculates AF and AN attributing current-time exposure to future outcomes.
#' @param basis If the `bdlnm` model has more than one basis, the name of the basis to use to compute predictions. It must be one of `names(object$basis)`. If the model contains only one basis it is selected automatically.
#' @param cen Numeric scalar; centering exposure value used to compute predictions. If missing the function will attempt to read it from the attributes of the specified `basis`. If no centering is available the function aborts.
#' @param range Optional numeric vector of length 2 with the exposure range for which attributable measures will be calculated. Values outside `range` are coerced to `cen` before prediction.
#' @param lag_average Logical (default `TRUE`). When `TRUE` use lag-averaged contributions to compute AN in the forward algorithm; when `FALSE` use the full lag-structured contributions instead.
#'
#' @details
#'
#' The function first obtains posterior predicted effects at the observed exposure values by calling [bcrosspred()] with `exp_at = data[[name_exposure]]`. Predictions must include relative-risk scale predictions so the previously fitted `"bdlnm"` model must have a `log` or `logit` link; otherwise the function aborts. These predictions require a defined centering value (`cen`), as attributable measures are always computed with respect to a reference exposure. This reference exposure is usually an optimal exposure computed with the [optimal_exposure()] function, such as the Minimum Mortality Temperature (MMT).
#'
#' Two different algorithms can be chosen to calculate attributable measures:
#'
#'  - Backward (`dir = "back"`): for each time point, contributions from past exposures (over the lag window) are combined to calculate the daily AF/AN.
#'  - Forward (`dir = "forw"`): for each time point, the contribution of that time point exposure to future outcomes (over the lag window) is combined to calculate the daily AF/AN.
#'
#' Both algorithms are fully described by Gasparrini and Leone (2014) <doi:10.1186/1471-2288-14-55>.
#'
#' Required columns to calculate `AF` and `AN` are `name_exposure` and `name_cases` columns. If `name_cases` is not supplied only AF per time can be computed and the output will only contain two elements: `$af` and `$af.summary`. If `name_date` is provided the function checks that dates are equispaced (checks seconds, minutes, hours, days, weeks, months or years). Time series have to be equispaced because the algorithms used to calculate attributable measures rely on consecutive time points over the lag window,. For example, if you only have seasonal observations (e.g., summers) expand the data to the full sequence and insert `NA` for missing exposures/cases and use `name_filter` to compute measures only for the seasonal subset.
#'
#' Only `"bdlnm"` objects fitted with a cross-basis are supported; models fitted with a one-basis (no lag) are not suitable for attributable calculations.
#'
#'
#' @return A list with components:
#'   - `af`: matrix (rows = time points, columns = posterior samples) with attributable fractions per time point.
#'   - `an`: matrix (rows = time points, columns = posterior samples) with attributable numbers per time point.
#'   - `aftotal`: numeric vector (length = number of posterior samples) with posterior samples of total AF across all the selected period.
#'   - `antotal`: numeric vector (length = number of posterior samples) with posterior samples of total AN across all the selected period.
#'   - `af.summary`: data.frame with summary statistics (mean, sd, quantiles, mode) for AF per time point.
#'   - `an.summary`: data.frame with summary statistics (mean, sd, quantiles, mode) for AN per time point.
#'   - `aftotal.summary`: data.frame with summary statistics (mean, sd, quantiles, mode) for total AF.
#'   - `antotal.summary`: data.frame with summary statistics (mean, sd, quantiles, mode) for total AN.
#'
#' @author Pau Satorra, Marcos Quijal-Zamorano.
#'
#' @note This function is inspired by `attrdl()` developed by Gasparrini and Leone (2014) <doi:10.1186/1471-2288-14-55>. It has been adapted to work in a Bayesian framework within the \pkg{bdlnm} package.
#'
#' @references
#'
#' Gasparrini A., Leone M. (2014). Attributable risk from distributed lag models. _BMC Medical Research Methodology_, 14, 55. <doi:10.1186/1471-2288-14-55>.
#' 
#' Quijal-Zamorano M., Martinez-Beneito M.A., Ballester J., Marí-Dell'Olmo M. (2024). Spatial Bayesian distributed lag non-linear models (SB-DLNM) for small-area exposure-lag-response epidemiological modelling. _International Journal of Epidemiology_, 53(3), dyae061. <doi:10.1093/ije/dyae061>.
#'
#' @seealso [bcrosspred()] to predict exposure–lag–response associations for a `"bdlnm"` object,
#' @seealso [bdlnm()] to fit a Bayesian distributed lag non-linear model (`"bdlnm"`),
#' @seealso [optimal_exposure()] to estimate exposure values that optimize the predicted effect for a `"bdlnm"` object.
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
#' if (check_inla()) {
#' mod <- bdlnm(mort_75plus ~ cb + factor(dow) + seas,
#'              data = london,
#'              family = "poisson",
#'              sample.arg = list(seed = 432, seed = 1L))
#'
#' # Predict
#' cpred <- bcrosspred(mod, exp_at = temp)
#'
#' # compute centering (MMT) using optimal_exposure
#' mmt <- optimal_exposure(mod, exp_at = temp)
#' cen <- mmt$summary[["0.5quant"]]
#'
#' # Attributable numbers and fractions (using the backwards algorithm):
#' attr <- attributable(mod, london, name_date = "date",
#' name_exposure = "tmean", name_cases = "mort_75plus", cen = cen, dir = "back")
#' }
#'
#'
attributable <- function(
  object,
  data,
  name_date = NULL,
  name_exposure,
  name_cases = NULL,
  name_filter = NULL,
  dir = "back",
  basis = NULL,
  cen,
  range = NULL,
  lag_average = TRUE
) {
  ## -----------------------
  ## Basic checks
  ## -----------------------

  # check object
  check_bdlnm(object)

  # basis
  name_basis <- basis
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

  # data
  if (missing(data)) {
    cli::cli_abort(
      "{.arg data} must be provided: a data.frame containing the time index and temporal series of exposures and (optionally) cases."
    )
  }

  if (!inherits(data, "data.frame")) {
    cli::cli_abort("{.arg data} must be a {.cls data.frame}.")
  }

  # name_exposure (required)
  if (
    missing(name_exposure) ||
      !is.character(name_exposure) ||
      length(name_exposure) != 1L
  ) {
    cli::cli_abort(
      "{.arg name_exposure} must be a single string with the name of the exposure column in {.arg data}."
    )
  }

  if (!name_exposure %in% colnames(data)) {
    cli::cli_abort(
      "Exposure column {.val {name_exposure}} not found in {.arg data}."
    )
  }

  exp <- data[[name_exposure]]

  # name_date (optional)
  if (!is.null(name_date)) {
    if (!is.character(name_date) || length(name_date) != 1L) {
      cli::cli_abort(
        "{.arg name_date} must be a single string naming the date column in {.arg data}."
      )
    }
    if (!name_date %in% colnames(data)) {
      cli::cli_abort("Date column {.val {name_date}} not found in {.arg data}.")
    }

    date <- data[[name_date]]
  } else {
    cli::cli_warn(c(
      "Ensure that {.arg data} contains time series measured on a regular basis (e.g., every day) to properly calculate attributable measures.",
      "i" = "If you have only seasonal observations (e.g., summers only), expand your data to the full sequence inserting NA for missing exposures/cases."
    ))

    date <- NULL
  }

  # name_cases (optional)
  if (!is.null(name_cases)) {
    if (!is.character(name_cases) || length(name_cases) != 1L) {
      cli::cli_abort(
        "{.arg name_cases} must be a single string naming the cases column in {.arg data}."
      )
    }
    if (!name_cases %in% colnames(data)) {
      cli::cli_abort(
        "Cases column {.val {name_cases}} not found in {.arg data}."
      )
    }
    cases <- data[[name_cases]]
    tot <- TRUE
  } else {
    cli::cli_warn(
      "Only attributable fractions per time point can be calculated as {.arg name_cases} is not provided."
    )
    cases <- NULL
    tot <- FALSE
  }

  # name_filter (optional)
  if (!is.null(name_filter)) {
    if (!is.character(name_filter) || length(name_filter) != 1L) {
      cli::cli_abort(
        "{.arg name_filter} must be a single string naming the filter column in {.arg data}."
      )
    }
    if (!name_filter %in% colnames(data)) {
      cli::cli_abort(
        "Filter column {.val {name_filter}} not found in {.arg data}."
      )
    }
    filter <- data[[name_filter]]
    cli::cli_warn(
      "Attributable fractions and numbers will only be calculated for time points filtered by {.val {name_filter}}"
    )
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
    cli::cli_abort(
      "Date column {.val {name_date}} must be of class {.cls Date} or {.cls POSIXt}."
    )
  }

  # validate filter
  if (!is.null(filter) && (!is.numeric(filter) | !all(filter %in% c(0, 1)))) {
    cli::cli_abort(
      "Filter column {.val {name_filter}} must contain only {.val 0} and {.val 1}."
    )
  }

  if (!is.null(filter) && all(filter == 0)) {
    cli::cli_abort(
      "Filter column  {.val {name_filter}} must contain some row equal to {.val 1}, corresponding to the filtered time points."
    )
  }

  # Here we require differences to be constant in time. Time can be on a seconds, minutes, hourly, daily, weekly, monthly, yearly basis.
  if (!is.null(date)) {
    # Check if it's sorted
    if (any(diff(date) < 1)) {
      cli::cli_abort(
        "The provided date ({.val {name_date}}) is not ordered in time."
      )
    }

    #Let's check regularity with different units
    is_regular <- length(unique(diff(date, units = "secs"))) == 1 ||
      length(unique(diff(date, units = "mins"))) == 1 ||
      length(unique(diff(date, units = "hours"))) == 1 ||
      length(unique(diff(date, units = "days"))) == 1 ||
      length(unique(diff(date, units = "weeks"))) == 1

    #Let's check if it's on a monthly or yearly basis:
    month_basis <- seq(date[1L], date[length(date)], by = "month")
    is_month_basis <- length(month_basis) == length(date) &&
      all(month_basis == date)

    year_basis <- seq(date[1L], date[length(date)], by = "year")
    is_year_basis <- length(year_basis) == length(date) &&
      all(year_basis == date)

    if (!is_regular && !is_month_basis && !is_year_basis) {
      cli::cli_abort(c(
        "The provided date ({.val {name_date}}) is not provided on a regular basis: time steps between dates are not equal.",
        "i" = "If you have only seasonal observations (e.g., summers only), expand your data to the full sequence inserting NA for missing exposures/cases."
      ))
    }
  }

  # direction
  if (!dir %in% c("back", "forw")) {
    cli::cli_abort("{.arg dir} must be one of: {.val 'back'}, {.val 'forw'}.")
  }

  # define centering value

  # get range of the exposure stored in the attributes of the crossbasis
  range_val <- attr(basis, "range")

  # if NULL try to extract it from basis
  if (is.null(cen)) {
    cen <- attributes(basis)$argvar$cen
    if (is.null(cen)) {
      cli::cli_abort(
        "Centering value must be specified ({.arg cen}) or specified as an attribute of the basis."
      )
    }
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
  lag <- attr(basis, "lag")
  predlag <- seq(from = lag[1], to = lag[2], by = 1)

  # obtain centered predictions at the exposure values
  cpred <- tryCatch(
    suppressMessages(bcrosspred(
      object,
      basis = name_basis,
      exp_at = exp,
      cen = cen
    )),
    error = function(e) {
      cli::cli_abort(
        "Failed computing predictions with {.fn bcrosspred}: {conditionMessage(e)}"
      )
    }
  )

  cp <- cpred$matfit

  if (!is.null(cpred$matRRfit)) {
    cp_rr <- cpred$matRRfit
  } else {
    # if there is no RR available in the model we can not calculate attributable numbers
    cli::cli_abort(
      "An attributable fraction cannot be computed because predicted effects are not on the relative-risk scale. Ensure the link of the {.arg object} model is {.val 'log'} or {.val 'logit'}."
    )
  }

  ## -----------------------
  ## Calculate AF/AN
  ## -----------------------

  # initialize matrices
  M_an <- M_af <- matrix(nrow = length(exp), ncol = n_sample) # per time
  an <- af <- numeric(n_sample) # total
  names(an) <- names(af) <- paste0("sample", seq_len(n_sample))

  # forward perspective: contributions from the current time point to future time points
  if (dir == "forw") {
    if (!is.null(cases)) {
      #Compute the Lagged matrix of daily cases for that time point and the next max_lag time points
      lagged_cases <- tsModel::Lag(cases, seq(-lag[1], -lag[2]))
    }

    for (i in seq_len(n_sample)) {
      # filter for i-th sample
      rr_sample <- cp_rr[match(exp, cpred$exp_at), , i]
      af_sample <- (rr_sample - 1) / rr_sample

      all_sample <- exp(rowSums(log(rr_sample)))
      M_af[, i] <- (all_sample - 1) / all_sample

      if (!is.null(cases)) {
        if (lag_average) {
          # average of lagged cases
          M_an[, i] <- M_af[, i] * rowMeans(lagged_cases)
        } else {
          # multiply element-wise by lagged cases and sum by row
          M_an[, i] <- rowSums(af_sample * lagged_cases)
        }
      }

      # total if allowed
      if (tot) {
        ind <- if (is.null(filter)) {
          !is.na(M_an[, i])
        } else {
          !is.na(M_an[, i]) & (filter == 1)
        }

        if (sum(ind) == length(cases)) {
          cli::cli_abort(c(
            "All time points are excluded by the filter or have a missing attributable number.",
            "i" = "Remember that if the number of cases is missing in some of the lags of a time point, an attributable number for that time point will be missing."
          ))
        }

        sum_cases <- sum(rowMeans(lagged_cases)[ind])

        if (sum_cases == 0) {
          cli::cli_abort(c(
            "The number of cases aggregated in the filtered time points with a non-missing attributable number is zero.",
            "i" = "Remember that if the number of cases is missing in some of the lags of a time point, an attributable number for that time point will be missing."
          ))
        }

        af[i] <- sum(M_an[ind, i]) / sum_cases
        an[i] <- af[i] * sum(cases, na.rm = TRUE)
      }
    }

    # backward perspective: contributions from past exposures to current time point
  } else {
    #Compute the Lagged matrix of daily exposures for that time point and the previous max_lag time points
    back_lagged_temp <- tsModel::Lag(exp, seq(lag[1], lag[2]))
    colnames(back_lagged_temp) <- paste0("lag", seq(lag[1], lag[2]))

    for (i in seq_len(n_sample)) {
      # calculate the matrix of backward rr: get the rr associated to each temperature and lag specified in the column of back_lagged_temp
      back_lagged_rr <- sapply(colnames(back_lagged_temp), function(col_lag) {
        rr <- cp_rr[
          match(back_lagged_temp[, col_lag], cpred$exp_at),
          col_lag,
          i
        ]
        return(rr)
      })
      back_rr <- exp(rowSums(log(back_lagged_rr)))

      # filter for the i-th sample
      af_sample <- (back_rr - 1) / back_rr

      M_af[, i] <- af_sample

      if (!is.null(cases)) {
        M_an[, i] <- M_af[, i] * cases
      }

      # total if allowed
      if (tot) {
        ind <- if (is.null(filter)) {
          !is.na(M_an[, i])
        } else {
          !is.na(M_an[, i]) & (filter == 1)
        }

        if (sum(ind) == length(cases)) {
          cli::cli_abort(c(
            "All time points are excluded by the filter or have a missing attributable number.",
            "i" = "Remember that if the number of cases is missing in some of the lags of a time point, an attributable number for that time point will be missing."
          ))
        }

        sum_cases <- sum(cases[ind])

        if (sum_cases == 0) {
          cli::cli_abort(c(
            "The number of cases aggregated in the filtered time points with a non-missing attributable number is zero.",
            "i" = "Remember that if the number of cases is missing in some of the lags of a time point, an attributable number for that time point will be missing."
          ))
        }

        af[i] <- sum(M_an[ind, i]) / sum_cases
        an[i] <- af[i] * sum(cases, na.rm = TRUE)
      }
    }
  }

  if (!is.null(date)) {
    rownames(M_an) <- rownames(M_af) <- as.character(date)
  } else {
    rownames(M_an) <- rownames(M_af) <- paste0("time", seq_along(exp))
  }
  colnames(M_an) <- colnames(M_af) <- paste0("sample", seq_len(n_sample))

  # Remove rows if filter is specified
  if (!is.null(filter)) {
    M_an <- M_an[filter == 1, ]
    M_af <- M_af[filter == 1, ]
  }

  ## -----------------------
  ## summaries
  ## -----------------------

  # only calculate summaries in points without missing values

  #per time
  M_ansum <- M_afsum <- matrix(
    nrow = nrow(M_af),
    ncol = ncol(cpred$allfit.summary)
  )
  rownames(M_ansum) <- rownames(M_afsum) <- rownames(M_an)
  colnames(M_ansum) <- colnames(M_afsum) <- colnames(cpred$allfit.summary)

  quant_cols <- grep("quant$", colnames(cpred$allfit.summary), value = TRUE)
  quant_val <- as.numeric(gsub("quant$", "", quant_cols))

  # AF don't have the same missings as AN (e.g, in forward perspective we have non-missing AF in the last time points whereas we have missing AN)
  ind_af <- !rowSums(is.na(M_af))

  if (!is.null(cases)) {
    ind_an <- !rowSums(is.na(M_an))

    M_ansum[ind_an, "mean"] <- apply(M_an[ind_an, ], 1, mean)
    M_ansum[ind_an, "sd"] <- apply(M_an[ind_an, ], 1, stats::sd)

    for (i in seq_along(quant_cols)) {
      M_ansum[ind_an, quant_cols[i]] <- apply(M_an[ind_an, ], 1, function(x) {
        stats::quantile(x, quant_val[i])
      })
    }

    #calculate mode (using default kernel density estimate, revisar...)
    M_ansum[ind_an, "mode"] <- apply(M_an[ind_an, ], 1, function(v) {
      dv <- stats::density(v)
      with(dv, x[which.max(y)])
    })
  }

  M_afsum[ind_af, "mode"] <- apply(M_af[ind_af, ], 1, function(v) {
    dv <- stats::density(v)
    with(dv, x[which.max(y)])
  })

  M_afsum[ind_af, "mean"] <- apply(M_af[ind_af, ], 1, mean)
  M_afsum[ind_af, "sd"] <- apply(M_af[ind_af, ], 1, stats::sd)

  for (i in seq_along(quant_cols)) {
    M_afsum[ind_af, quant_cols[i]] <- apply(M_af[ind_af, ], 1, function(x) {
      stats::quantile(x, quant_val[i])
    })
  }

  #total (if allowed)
  if (tot) {
    ansum <- afsum <- numeric(ncol(cpred$allfit.summary))
    names(ansum) <- names(afsum) <- colnames(cpred$allfit.summary)

    quant_cols <- grep("quant$", colnames(cpred$allfit.summary), value = TRUE)
    quant_val <- as.numeric(gsub("quant$", "", quant_cols))

    if (!is.null(cases)) {
      ansum["mean"] <- mean(an)
      ansum["sd"] <- stats::sd(an)

      for (i in seq_along(quant_cols)) {
        ansum[quant_cols[i]] <- stats::quantile(an, quant_val[i])
      }

      #calculate mode (using default kernel density estimate, revisar...)
      dv <- stats::density(an)
      ansum["mode"] <- with(dv, x[which.max(y)])
    }

    #calculate mode (using default kernel density estimate, revisar...)
    dv <- stats::density(af)
    afsum["mode"] <- with(dv, x[which.max(y)])

    afsum["mean"] <- mean(af)
    afsum["sd"] <- stats::sd(af)

    for (i in seq_along(quant_cols)) {
      afsum[quant_cols[i]] <- stats::quantile(af, quant_val[i])
    }
  }

  # Return results as data.frames for summaries
  if (!is.null(cases)) {
    if (tot) {
      res <- list(
        af = M_af,
        an = M_an,
        aftotal = af,
        antotal = an,
        af.summary = M_afsum,
        an.summary = M_ansum,
        aftotal.summary = afsum,
        antotal.summary = ansum
      )
    } else {
      res <- list(
        af = M_af,
        an = M_an,
        af.summary = M_afsum,
        an.summary = M_ansum
      )
    }
  } else {
    res <- list(
      af = M_af,
      af.summary = as.data.frame(M_afsum, stringsAsFactors = FALSE)
    )
  }

  return(res)
}
