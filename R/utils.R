### Util functions for the R package bdlnm

# Function to select only the coefficients from the basis from all the posterior samples
extract_coef <- function (coef, basis) {

  # get estimated coefficients name
  names_sel <- rownames(coef)

  #Revisar per onebasis si tmb té colnames
  basis_cols <- colnames(basis)

  # intersection between model fixed names and basis column names
  names_sel <- names_sel[names_sel %in% basis_cols]

  if (length(names_sel) == 0L) {
    cli::cli_abort(
      "No matching coefficient names found between the fitted {.fun bdlnm()} model fixed effects and the provided basis columns. Make sure the {.arg basis} columns are included in the {.arg formula} of the previously {.fun bdlnm()} fitted model."
    )
  }

  if(length(names_sel) != length(basis_cols)) {
    cli::cli_abort(
      "Not all parameters of the provided basis are in the estimated coefficients from the fitted {.fun bdlnm()} model. Make sure the provided {.arg basis} is the same of the basis provided previously in the {.fun bdlnm()} fitted model."
    )
  }

  # extract only selected basis coefficients
  coef_cb <- coef[names_sel,]

  return(coef_cb)

}

# Verify that an object is the expected result from bdlnm().
check_bdlnm <- function(object) {

  if(is.null(object)) {
    cli::cli_abort("The object returned by {.fn bdlnm} must be provided as the {.arg object} argument.")
  }

  if (!inherits(object, "bdlnm") || is.null(object$model) || is.null(object$coefficients) || is.null(object$coefficients.summary)) {
    cli::cli_abort("{.arg object} must be an object of class {.cls bdlnm} returned by {.fn bdlnm}.")
  }

  model <- object$model
  coef <- object$coefficients

  # check model
  if (!inherits(model, "inla")) {
    cli::cli_abort("{.arg object$model} must be an {.cls inla} model (object returned by {.fn INLA::inla}).")
  } else {
    if (!is.null(model$.args) && !is.null(model$.args$control.compute)) {
      if(!model$.args$control.compute$config) cli::cli_abort("{.arg object$model} must be fitted with {.arg control.compute = list(config = TRUE)} so that {.fn INLA::inla.posterior.sample} can be run afterwards.")
    } else {
      cli::cli_abort("{.code $.args$control.compute} cannot be extracted from {.arg object$model}.")
    }

  }

  # check coefficients
  if (!is.matrix(coef) || !is.numeric(coef)) {
    cli::cli_abort("{.arg object$coefficients} must be a numeric matrix of posterior samples as returned by {.fn INLA::inla.posterior.sample} (columns = samples).")
  }

}

# Extract the model link for an INLA model object.
get_link <- function(model) {

  # get model link from model. INLA can store the family/link in different places
  link <- NA

  # 1) explicit control.family structure
  if(!is.null(model$.args) && !is.null(model$.args$control.family) && !is.null(model$.args$control.family[[1]]$link) && model$.args$control.family[[1]]$link != "default") {
    link <- model$.args$control.family[[1]]$link
  }

  # 2) get from model$.args$family
  if(!is.null(model$.args$family)) {

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
    cli::cli_warn("Could not unambiguously determine the model link from the fitted {.pkg INLA} object. If you know the link, specify it in {.arg model.link}.")
  }

  return(link)

}


# Draw credible-interval panels in plots.
# All credits for this function go to the dlnm package. We have only included the option for 'sampling'.
fci <- function(ci, x, y, high, low, ci.arg, plot.arg, noeff = NULL) {

    if(ci=="area") {
      polygon.arg <- utils::modifyList(list(col=grDevices::grey(0.9),border=NA),ci.arg)
      polygon.arg <- utils::modifyList(polygon.arg,
                                list(x=c(x,rev(x)),y=c(high,rev(low))))
      do.call(graphics::polygon,polygon.arg)
    } else if(ci=="bars") {
      range <- diff(range(x))/300
      segments.arg <- utils::modifyList(ci.arg,list(x0=x,y0=high,x1=x,y1=low))
      do.call(graphics::segments,segments.arg)
      segments.arg <- utils::modifyList(segments.arg,list(x0=x-range,y0=high,
                                                   x1=x+range,y1=high))
      do.call(graphics::segments,segments.arg)
      segments.arg <- utils::modifyList(segments.arg,list(x0=x-range,y0=low,
                                                   x1=x+range,y1=low))
      do.call(graphics::segments,segments.arg)
    } else if(ci=="lines") {
      lines.arg <- list(lty=2)
      if(!is.null(plot.arg$col)) lines.arg$col <- plot.arg$col
      lines.arg <- utils::modifyList(lines.arg,ci.arg)
      lines.arg <- utils::modifyList(lines.arg,list(x=x,y=high))
      do.call(graphics::lines,lines.arg)
      lines.arg <- utils::modifyList(lines.arg,list(x=x,y=low))
      do.call(graphics::lines,lines.arg)
    } else if(ci == "sampling") {
      lines.arg <- list(col=grDevices::grey(0.9), lty = 1)
      lines.arg <- utils::modifyList(lines.arg,ci.arg)
      for(i in 1:ncol(y)) {
        lines.arg <- utils::modifyList(lines.arg,list(x=x,y=y[,i]))
        do.call(graphics::lines, lines.arg)
      }
    }
    if(!is.null(noeff)) graphics::abline(h=noeff)
  }

