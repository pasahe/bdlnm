### Util functions for the R package bdlnm

# Function to check if INLA is installed
check_inla <- function(error = FALSE) {
  if (!requireNamespace("INLA", quietly = TRUE)) {
    if (error) {
      cli::cli_abort(c(
        "Package {.pkg INLA} is required but is not installed.",
        "i" = "Install from the R-INLA repository (https://www.r-inla.org/) and restart R."
      ))
    }
    return(FALSE)
  }

  version <- tryCatch(
    utils::packageVersion("INLA"),
    error = function(e) NA
  )

  if (is.na(version) || version < "23.4.24") {
    if (error) {
      cli::cli_abort(c(
        "Installed {.pkg INLA} version ({version}) is too old. Version >= 23.4.24 is required.",
        "i" = "Install a newest version from the R-INLA repository (https://www.r-inla.org/) and restart R."
      ))
    }
    return(FALSE)
  }

  TRUE
}

# Function to select only the coefficients from the basis from all the posterior samples
extract_coef <- function(object, basis, type = "coefficients") {
  # find if basis is present in object$basis
  which_basis <- NULL
  for (b in names(object$basis)) {
    if (identical(object$basis[[b]], basis)) {
      which_basis <- b
      break
    }
  }

  if (is.null(which_basis)) {
    cli::cli_abort(
      c(
        "The provided {.arg basis} was not found in the previous fitted model.",
        "Make sure the provided {.arg basis} is the same of the basis used previously in the {.fun bdlnm()} fitted model."
      )
    )
  }

  coef <- object[[type]]

  # extract only selected basis coefficients
  coef_names <- paste0(which_basis, colnames(basis))
  missing <- setdiff(coef_names, rownames(coef))
  if (length(missing) > 0) {
    cli::cli_abort(
      c(
        "The following basis coefficients don't match between the basis included in the previous fitted model and the provided {.arg basis}",
        "x" = "{.val {missing}}"
      )
    )
  }
  coef_cb <- coef[coef_names, , drop = FALSE]

  return(coef_cb)
}

# Verify that an object is the expected result from bdlnm().
check_bdlnm <- function(object) {
  if (is.null(object)) {
    cli::cli_abort(
      "The object returned by {.fn bdlnm} must be provided as the {.arg object} argument."
    )
  }

  if (
    !inherits(object, "bdlnm") ||
      is.null(object$model) ||
      is.null(object$basis) ||
      is.null(object$coefficients) ||
      is.null(object$coefficients.summary)
  ) {
    cli::cli_abort(
      "{.arg object} must be an object of class {.cls bdlnm} returned by {.fn bdlnm}."
    )
  }

  model <- object$model
  basis <- object$basis
  coef <- object$coefficients

  # check model
  if (!inherits(model, "inla")) {
    cli::cli_abort(
      "{.arg object$model} must be an {.cls inla} model (object returned by {.fn INLA::inla})."
    )
  } else {
    if (!is.null(model$.args) && !is.null(model$.args$control.compute)) {
      if (!model$.args$control.compute$config) {
        cli::cli_abort(
          "{.arg object$model} must be fitted with {.arg control.compute = list(config = TRUE)} so that {.fn INLA::inla.posterior.sample} can be run afterwards."
        )
      }
    } else {
      cli::cli_abort(
        "{.code $.args$control.compute} cannot be extracted from {.arg object$model}."
      )
    }
  }

  # check basis
  if (is.list(basis)) {
    is_basis <- vapply(
      basis,
      inherits,
      logical(1),
      what = c("crossbasis", "onebasis")
    )
    if (!any(is_basis)) {
      cli::cli_abort(
        "{.arg object$basis} must be a list containing basis of class {.cls crossbasis} or {.cls onebasis}."
      )
    }
  } else {
    cli::cli_abort(
      "{.arg object$basis} must be a list containing basis of class {.cls crossbasis} or {.cls onebasis}."
    )
  }

  # check coefficients
  if (!is.matrix(coef) || !is.numeric(coef)) {
    cli::cli_abort(
      "{.arg object$coefficients} must be a numeric matrix of posterior samples as returned by {.fn INLA::inla.posterior.sample} (columns = samples)."
    )
  }
}

# Extract the model link for an INLA model object.
get_link <- function(model) {
  # get model link from model. INLA can store the family/link in different places
  link <- NA

  # 1) explicit control.family structure
  if (
    !is.null(model$.args) &&
      !is.null(model$.args$control.family) &&
      !is.null(model$.args$control.family[[1]]$link) &&
      model$.args$control.family[[1]]$link != "default"
  ) {
    link <- model$.args$control.family[[1]]$link
  }

  # 2) get from model$.args$family
  if (!is.null(model$.args$family)) {
    fam <- tolower(model$.args$family)

    link <- if (fam %in% c("poisson", "coxph", "exponential")) {
      "log"
    } else if (fam %in% c("binomial")) {
      "logit"
    } else if (fam %in% c("gaussian")) {
      "identity"
    }
  }

  if (is.na(link)) {
    cli::cli_warn(
      "Could not unambiguously determine the model link from the fitted {.pkg INLA} object. If you know the link, specify it in {.arg model.link}."
    )
  }

  return(link)
}


# Draw credible-interval panels in plots.
# All credits for this function go to the dlnm package. We have only included the option for 'sampling'.
fci <- function(ci, x, y, high, low, ci.arg, plot.arg, noeff = NULL) {
  if (ci == "area") {
    polygon.arg <- utils::modifyList(
      list(col = grDevices::grey(0.9), border = NA),
      ci.arg
    )
    polygon.arg <- utils::modifyList(
      polygon.arg,
      list(x = c(x, rev(x)), y = c(high, rev(low)))
    )
    do.call(graphics::polygon, polygon.arg)
  } else if (ci == "bars") {
    range <- diff(range(x)) / 300
    segments.arg <- utils::modifyList(
      ci.arg,
      list(x0 = x, y0 = high, x1 = x, y1 = low)
    )
    do.call(graphics::segments, segments.arg)
    segments.arg <- utils::modifyList(
      segments.arg,
      list(x0 = x - range, y0 = high, x1 = x + range, y1 = high)
    )
    do.call(graphics::segments, segments.arg)
    segments.arg <- utils::modifyList(
      segments.arg,
      list(x0 = x - range, y0 = low, x1 = x + range, y1 = low)
    )
    do.call(graphics::segments, segments.arg)
  } else if (ci == "lines") {
    lines.arg <- list(lty = 2)
    if (!is.null(plot.arg$col)) {
      lines.arg$col <- plot.arg$col
    }
    lines.arg <- utils::modifyList(lines.arg, ci.arg)
    lines.arg <- utils::modifyList(lines.arg, list(x = x, y = high))
    do.call(graphics::lines, lines.arg)
    lines.arg <- utils::modifyList(lines.arg, list(x = x, y = low))
    do.call(graphics::lines, lines.arg)
  } else if (ci == "sampling") {
    lines.arg <- list(col = grDevices::grey(0.9), lty = 1)
    lines.arg <- utils::modifyList(lines.arg, ci.arg)
    for (i in 1:ncol(y)) {
      lines.arg <- utils::modifyList(lines.arg, list(x = x, y = y[, i]))
      do.call(graphics::lines, lines.arg)
    }
  }
  if (!is.null(noeff)) graphics::abline(h = noeff)
}
