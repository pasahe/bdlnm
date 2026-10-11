#' Calculate attributable numbers and fractions from posterior predictions
#'
#' Compute attributable numbers (AN) and attributable fractions (AF) from posterior predicted relative risks returned by [bcrosspred()]. The function applies a forward or backward lag algorithm to compute per-time and total attributable measures, optionally filtering the calculation to a subset of time points.
#'
#' @param object A `"bcrosspred"` object returned by [bcrosspred()]. Passing a fitted `"bdlnm"` object is deprecated.
#' @param data A data frame containing the temporal series needed to calculate attributable measures.
#' @param name_date,name_exposure,name_cases,name_filter Deprecated column-name arguments. Use `columns` instead.
#' @param columns Named list mapping the roles required by `attributable()` to columns in `data`. The `exp` element is required; `cases`, `date`, `filter`, and `int` are optional. See `Details` for the definition and requirements of each element.
#' @param dir Character; direction of the algorithm to calculate attributable measures. `"back"` (default) calculates AF and AN attributing current-time outcome to past exposures; `"forw"` calculates AF and AN attributing current-time exposure to future outcomes.
#' @param basis Deprecated. Only used when `object` is a `"bdlnm"` model. If the model has more than one basis, the name of the basis to use.
#' @param cen Deprecated. Only used when `object` is a `"bdlnm"` model. Numeric scalar specifying the centering exposure value.
#' @param lag_average Logical (default `TRUE`). When `TRUE` use lag-averaged contributions to compute AN in the forward algorithm; when `FALSE` use the full lag-structured contributions instead.
#'
#' @details
#'
#' The supported elements of `columns` are:
#'
#' - `exp`: required single string naming the numeric exposure column.
#' - `cases`: optional single string naming the numeric cases column. If omitted, only attributable fractions per time point are returned.
#' - `date`: optional single string naming a `Date` or `POSIXt` column. When supplied, the function checks that observations are ordered and regularly spaced.
#' - `filter`: optional single string naming a numeric binary (`0`/`1`) column. Only observations equal to `1` are included in total attributable measures and returned in the per-time results.
#' - `int`: optional single string naming the interaction column that assigns each observation to the corresponding interaction-specific prediction. It is required when `object` contains interaction-specific predictions, and its name must match the interaction variable recorded in `object`.
#'
#' `data` must represent one regularly spaced time series; multiple concurrent series or panel data are not supported. If only seasonal observations are available, expand the data to the full sequence, insert `NA` for missing exposures or cases, and use `filter` to select the relevant observations.
#'
#' The prediction grid in `object` must include every non-missing exposure value observed in the data's exposure column. Values are matched directly to the `exp_at` grid of each relevant prediction. Therefore, before calculating attributable measures, construct centered relative-risk predictions with an exposure grid that includes the unique non-missing values in `data[[columns$exp]]` (ex. `bcrosspred(cpred, exp_at = unique(data[[columns$exp]]))`). 
#'
#' Two different algorithms can be chosen to calculate attributable measures:
#'
#'  - Backward (`dir = "back"`): for each time point, contributions from past exposures (over the lag window) are combined to calculate the daily AF/AN.
#'  - Forward (`dir = "forw"`): for each time point, the contribution of that time point exposure to future outcomes (over the lag window) is combined to calculate the daily AF/AN.
#'
#' Both algorithms are fully described by Gasparrini and Leone (2014) <doi:10.1186/1471-2288-14-55>.
#'
#' For interaction-specific predictions, the backward algorithm uses the interaction value at the outcome time for all lagged exposures contributing to that outcome. For each exposure time, the forward algorithm uses the interaction value at the future outcome time corresponding to each lag.
#'
#' Only predictions from a cross-basis are supported; one-basis predictions with no lag dimension are not suitable for attributable calculations.
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
#'  # Prediction values (equidistant points)
#'  temp <- round(seq(min(london$tmean), max(london$tmean), by = 0.1), 1)
#'  # Ensure it falls inside the range of temperatures after rounding:
#'  temp <- temp[temp >= min(london$tmean) & temp <= max(london$tmean)]
#'
#' # Model
#'
#' \dontrun{
#' if (check_inla()) {
#' mod <- bdlnm(mort_75plus ~ cb + factor(dow) + seas,
#'              data = london,
#'              family = "poisson",
#'              sample.arg = list(n = 1000, seed = 1L))
#'
#' # compute centering (MMT) using optimal_exposure
#' cpred <- bcrosspred(mod, exp_at = temp)
#' mmt <- optimal_exposure(cpred)
#' cen <- mmt$summary[["0.5quant"]]
#'
#' # Centered predictions covering the observed exposures
#' cpred <- bcrosspred(mod, exp_at = london$tmean, cen = cen)
#'
#' # Attributable numbers and fractions (using the backwards algorithm):
#' attr <- attributable(
#'   cpred,
#'   london,
#'   columns = list(
#'     date = "date",
#'     exp = "tmean",
#'     cases = "mort_75plus"
#'   ),
#'   dir = "back"
#' )
#' }
#' }
#'
#'
attributable <- function(
  object,
  data,
  name_date = NULL,
  name_exposure = NULL,
  name_cases = NULL,
  name_filter = NULL,
  dir = "back",
  basis = NULL,
  cen,
  lag_average = TRUE,
  columns = NULL
) {
  ## -----------------------
  ## Basic checks
  ## -----------------------

  cen_missing <- missing(cen)
  is_prediction <- inherits(object, "bcrosspred")

  if (is_prediction) {
    if (!is.list(object)) {
      cli::cli_abort(
        "{.arg object} must be a valid {.cls bcrosspred} object returned by {.fn bcrosspred}."
      )
    }

    if (!is.null(basis) || !cen_missing) {
      cli::cli_abort(c(
        "Prediction arguments cannot be supplied when {.arg object} is a {.cls bcrosspred} object.",
        "i" = "Select the basis and centering value when calling {.fn bcrosspred}."
      ))
    }

    required_components <- c("exp_at", "lag_at", "matfit", "allfit.summary")
    interaction_info <- attr(object, "basis_interaction")
    is_interaction_prediction <- !is.null(interaction_info)

    if (is_interaction_prediction) {
      prediction_names <- named_result_names(
        object,
        required_components = required_components
      )

      if (length(prediction_names) == 0L) {
        cli::cli_abort(
          "{.arg object} must be a valid cross-basis {.cls bcrosspred} object returned by {.fn bcrosspred}."
        )
      }

      predictions <- object[prediction_names]
    } else {
      if (!all(required_components %in% names(object))) {
        cli::cli_abort(
          "{.arg object} must be a valid cross-basis {.cls bcrosspred} object returned by {.fn bcrosspred}."
        )
      }
      prediction_names <- character(0)
      predictions <- list(global = object)
    }

    cpred <- predictions[[1L]]
    cen <- cpred$cen
    if (is.null(cen)) {
      cli::cli_abort(c(
        "{.arg object} must contain predictions centered at a reference exposure.",
        "i" = "Rerun {.fn bcrosspred} with a numeric {.arg cen} value."
      ))
    }
  } else {
    is_interaction_prediction <- FALSE
    check_bdlnm(object)

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
    } else if (length(obj_basis) > 1) {
      cli::cli_abort(
        "{.arg basis} must be provided containing the name of the basis to select for predictions (stored in {.arg object$basis}): {.or {.val {names(obj_basis)}}}."
      )
    } else {
      basis <- obj_basis[[1]]
    }

    selected_basis_name <- if (is.null(name_basis)) {
      names(obj_basis)[[1L]]
    } else {
      name_basis
    }
    basis_interactions <- attr(object, "basis_interaction")[[
      selected_basis_name
    ]]
    if (!is.null(basis_interactions)) {
      cli::cli_abort(c(
        "A {.cls bdlnm} model containing basis interactions cannot be used directly in {.fn attributable}.",
        "i" = "Call {.fn bcrosspred} with the appropriate {.arg combine_effects}, then pass the resulting prediction object to {.fn attributable}."
      ))
    }

    if (cen_missing || is.null(cen)) {
      cen <- attributes(basis)$argvar$cen
      if (is.null(cen)) {
        cli::cli_abort(
          "Centering value must be specified ({.arg cen}) or specified as an attribute of the basis."
        )
      }
    }

    cli::cli_warn(c(
      "Using a {.cls bdlnm} object directly in {.fn attributable} is deprecated.",
      "i" = "Call {.fn bcrosspred} first with an exposure grid covering the observed exposures and a centering value, then pass its result to {.fn attributable}.",
      "i" = "For interaction models, select the required effects with {.arg combine_effects} when calling {.fn bcrosspred}."
    ))
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

  legacy_columns <- c(
    date = !missing(name_date),
    exp = !missing(name_exposure),
    cases = !missing(name_cases),
    filter = !missing(name_filter)
  )

  if (!is.null(columns)) {
    if (!is.list(columns) || is.null(names(columns)) || any(names(columns) == "")) {
      cli::cli_abort(
        "{.arg columns} must be a named list."
      )
    }

    allowed_columns <- c("exp", "cases", "date", "filter", "int")
    unknown_columns <- setdiff(names(columns), allowed_columns)
    if (length(unknown_columns) > 0L) {
      cli::cli_abort(c(
        "Unknown element{?s} in {.arg columns}: {.val {unknown_columns}}.",
        "i" = "Allowed elements are {.or {.val {allowed_columns}}}."
      ))
    }

    if (anyDuplicated(names(columns))) {
      cli::cli_abort("Elements in {.arg columns} must have unique names.")
    }

    if (is.null(columns$exp)) {
      cli::cli_abort("{.arg columns} must contain an {.val exp} element.")
    }

    overlapping_columns <- intersect(names(columns), names(legacy_columns)[legacy_columns])
    if (length(overlapping_columns) > 0L) {
      cli::cli_abort(c(
        "Column roles cannot be supplied in both {.arg columns} and deprecated {.arg name_*} arguments.",
        "x" = "Supplied twice: {.val {overlapping_columns}}"
      ))
    }

    name_exposure <- columns$exp
    if (!is.null(columns$date)) name_date <- columns$date
    if (!is.null(columns$cases)) name_cases <- columns$cases
    if (!is.null(columns$filter)) name_filter <- columns$filter
    name_int <- columns$int
  } else {
    name_int <- NULL
  }

  if (any(legacy_columns)) {
    cli::cli_warn(c(
      "The {.arg name_date}, {.arg name_exposure}, {.arg name_cases}, and {.arg name_filter} arguments are deprecated.",
      "i" = "Supply column names through {.arg columns}."
    ))
  }

  column_values <- columns[!vapply(columns, is.null, logical(1))]
  if (length(column_values) > 0L) {
    invalid_columns <- names(column_values)[!vapply(
      column_values,
      function(x) is.character(x) && length(x) == 1L,
      logical(1)
    )]
    if (length(invalid_columns) > 0L) {
      cli::cli_abort(
        "Each element of {.arg columns} must be a single string naming a column in {.arg data}."
      )
    }

    missing_columns <- unlist(column_values, use.names = FALSE)
    missing_columns <- missing_columns[!missing_columns %in% names(data)]
    if (length(missing_columns) > 0L) {
      cli::cli_abort(
        "Column{?s} {.val {missing_columns}} from {.arg columns} {?is/are} not present in {.arg data}."
      )
    }
  }

  # name_exposure (required)
  if (
    is.null(name_exposure) ||
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

  if (!is.null(name_int)) {
    int <- as.character(data[[name_int]])
  } else {
    int <- NULL
  }

  if (is_interaction_prediction) {
    if (is.null(interaction_info) || is.null(interaction_info$variable)) {
      cli::cli_abort(
        "{.arg object} does not contain the interaction metadata required by {.fn attributable}."
      )
    }

    if (is.null(name_int)) {
      cli::cli_abort(c(
        "{.arg columns} must contain an {.val int} element for interaction predictions.",
        "i" = "Use the interaction variable {.val {interaction_info$variable}}."
      ))
    }

    if (!identical(name_int, interaction_info$variable)) {
      cli::cli_abort(c(
        "The {.val int} element of {.arg columns} does not match the interaction in {.arg object}.",
        "x" = "Expected {.val {interaction_info$variable}}, not {.val {name_int}}."
      ))
    }

    observation_prediction <- paste0(name_int, "[", int, "]")
    missing_predictions <- unique(
      observation_prediction[!is.na(int) & !observation_prediction %in% prediction_names]
    )
    if (length(missing_predictions) > 0L) {
      cli::cli_abort(c(
        "{.arg object} does not contain predictions for every interaction value in {.arg data}.",
        "x" = "{length(missing_predictions)} interaction value{?s} {?is/are} missing."
      ))
    }
    observation_prediction[is.na(int)] <- NA_character_
  } else {
    if (!is.null(name_int)) {
      cli::cli_abort(
        "The {.val int} element of {.arg columns} can only be used with interaction predictions."
      )
    }
    observation_prediction <- rep("global", nrow(data))
  }

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

  if (!is.numeric(cen) || length(cen) != 1L) {
    cli::cli_abort("{.arg cen} must be a numeric scalar.")
  }

  ## -----------------------
  ## Prepare
  ## -----------------------

  if (!is_prediction) {
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
    predictions <- list(global = cpred)
  }

  observed_exp <- unique(exp[!is.na(exp)])
  missing_exp <- observed_exp[!observed_exp %in% cpred$exp_at]
  if (length(missing_exp) > 0L) {
    cli::cli_abort(c(
      "{.arg object} does not contain predictions for every observed exposure in {.arg data}.",
      "x" = "{length(missing_exp)} unique exposure value{?s} {?is/are} missing.",
      "i" = "Call {.fn bcrosspred} with an {.arg exp_at} grid equal to {.code unique(data[[name_exposure]])}."
    ))
  }

  lag <- range(cpred$lag_at)

  has_relative_risk <- vapply(
    predictions,
    function(prediction) !is.null(prediction$matRRfit),
    logical(1L)
  )
  if (!all(has_relative_risk)) {
    # if there is no RR available in the model we can not calculate attributable numbers
    cli::cli_abort(
      "An attributable fraction cannot be computed because predictions are not on the relative-risk scale. Ensure {.fn bcrosspred} uses a model with a {.val 'log'} or {.val 'logit'} link."
    )
  }

  n_sample <- dim(predictions[[1L]]$matRRfit)[3L]
  lag_names <- paste0("lag", seq(lag[1], lag[2]))

  get_relative_risk <- function(exposure, prediction_id, lag_name, sample) {
    relative_risk <- rep(NA_real_, length(exposure))
    valid <- !is.na(exposure) & !is.na(prediction_id)

    for (id in unique(prediction_id[valid])) {
      index <- valid & !is.na(prediction_id) & prediction_id == id
      prediction <- predictions[[id]]
      relative_risk[index] <- prediction$matRRfit[
        match(exposure[index], prediction$exp_at),
        lag_name,
        sample
      ]
    }

    relative_risk
  }

  ## -----------------------
  ## Calculate AF/AN
  ## -----------------------

  # initialize matrices
  M_an <- M_af <- matrix(nrow = length(exp), ncol = n_sample) # per time
  an <- af <- numeric(n_sample) # total
  names(an) <- names(af) <- paste0("sample", seq_len(n_sample))

  if (tot) {
    case_ind <- if (is.null(filter)) {
      rep(TRUE, length(cases))
    } else {
      filter == 1
    }
    total_cases <- sum(cases[case_ind], na.rm = TRUE)
  }

  # forward perspective: contributions from the current time point to future time points
  if (dir == "forw") {
    if (!is.null(cases)) {
      #Compute the Lagged matrix of daily cases for that time point and the next max_lag time points
      lagged_cases <- tsModel::Lag(cases, seq(-lag[1], -lag[2]))
    }

    if (is_interaction_prediction) {
      future_prediction <- matrix(
        NA_character_,
        nrow = length(exp),
        ncol = length(lag_names),
        dimnames = list(NULL, lag_names)
      )
      time_index <- seq_along(exp)
      for (j in seq_along(lag_names)) {
        outcome_index <- time_index + lag[[1L]] + j - 1L
        valid_index <- outcome_index >= 1L & outcome_index <= length(exp)
        future_prediction[valid_index, j] <- observation_prediction[
          outcome_index[valid_index]
        ]
      }
    } else {
      future_prediction <- matrix(
        "global",
        nrow = length(exp),
        ncol = length(lag_names),
        dimnames = list(NULL, lag_names)
      )
    }

    for (i in seq_len(n_sample)) {
      # filter for i-th sample
      rr_sample <- vapply(seq_along(lag_names), function(j) {
        get_relative_risk(
          exposure = exp,
          prediction_id = future_prediction[, j],
          lag_name = lag_names[[j]],
          sample = i
        )
      }, numeric(length(exp)))
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

        if (!any(ind)) {
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
        an[i] <- af[i] * total_cases
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
        get_relative_risk(
          exposure = back_lagged_temp[, col_lag],
          prediction_id = observation_prediction,
          lag_name = col_lag,
          sample = i
        )
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

        if (!any(ind)) {
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
        an[i] <- af[i] * total_cases
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

  if (is_interaction_prediction) {
    attr(res, "basis_interaction") <- interaction_info
  }

  return(res)
}
