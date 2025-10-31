### Util functions for the R package bdlnm

# Verify that an object is the expected result from bdlnm().
check_bdlnm <- function(x) {

  if(is.null(x)) {
    cli::cli_abort("The object returned by {.fn bdlnm} must be provided as the {.arg x} argument.")
  }

  if (!is.list(x) || is.null(x$model) || is.null(x$coef)) {
    cli::cli_abort("{.arg x} must be the list returned by {.fn bdlnm}: it should contain the fitted INLA model {.val $model} and the posterior samples matrix {.val $coef}.")
  }

  model <- x$model
  coef <- x$coef

  # check model
  if (!inherits(model, "inla")) {
    cli::cli_abort("{.arg x$model} must be an {.cls inla} model (object returned by {.fn INLA::inla}).")
  } else {
    if (!is.null(model$.args) && !is.null(model$.args$control.compute)) {
      if(!model$.args$control.compute$config) cli::cli_abort("{.arg x$model} must be fitted with {.arg control.compute = list(config = TRUE)} so that {.fn INLA::inla.posterior.sample} can be run afterwards.")
    } else {
      cli::cli_abort("{.code $.args$control.compute} cannot be extracted from {.arg x$model}.")
    }

  }

  # check coefficients
  if (!is.matrix(coef) || !is.numeric(coef)) {
    cli::cli_abort("{.arg x$coef} must be a numeric matrix of posterior samples as returned by {.fn INLA::inla.posterior.sample} (columns = samples).")
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

