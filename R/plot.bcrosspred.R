
#' plot.bcrosspred
#' This function generates plots for bayesian distributed lag non-linear model's predictions stored in a bcrosspred object.
#' @param x
#' @param ptype
#' @param var
#' @param lag
#' @param ci
#' @param ci.arg
#' @param ci.level
#' @param cumul
#' @param exp
#' @param ...
#'
#' @returns
#' @export
#'
#' @examples
plot.bcrosspred <- function(x, ptype, var=NULL, lag=NULL, ci="area", ci.arg, ci.level = 0.95, cumul=FALSE, exp=NULL, ...) {

  ######################
  #CHECK

  #Check if it matches the possible choices for ci
  ci <- match.arg(ci, c("area","bars","lines","n"))

  # SETTING DEFAULT FOR ptype: OVERALL FOR NO LAG, SLICES FOR VAR/LAG, OTHERWISE 3D
  if(missing(ptype)) {
    if(!is.null(var)||!is.null(lag)) {
      ptype <- "slices"
    } else if(diff(x$lag)==0) {
      ptype <- "overall"
    } else ptype <- "3d"
  }

  #Check if matches the possible choices for ptype
  ptype <- match.arg(ptype, c("slices","3d","contour","overall"))

  if(is.null(var)&&is.null(lag)&&(ptype=="slices")) {
    stop("at least 'var' or 'lag' must be provided when ptype='slices'")
  }
  if(!is.null(var)&&!is.numeric(var)&&length(var)>4&&ptype=="slices") {
    stop("'var' and 'lag' must be numeric and of length <=4")
  }
  if(!is.null(lag)&&!is.numeric(lag)&&length(lag)>4&&ptype=="slices") {
    stop("'var' and 'lag' must be numeric and of length <=4")
  }
  if((!is.null(var)&!is.null(lag))&&length(var)!=length(lag)&&ptype=="slices") {
    stop("if both are provided, length of 'var' and 'lag' must be the same")
  }
  if(!is.null(var)&&sum(var%in%x$predvar)!=length(var)&&(ptype=="slices")) {
    stop("'var' must match values used for prediction")
  }
  if(!is.null(lag)&&sum(lag%in%seqlag(x$lag,x$bylag))!=length(lag)&&(ptype=="slices")) {
    stop("'lag' must match values used for prediction")
  }
  if(missing(ci.arg)) {
    ci.arg <- list()
  } else if(!is.list(ci.arg)) stop("'ci.arg' must be a list")
  if(!is.numeric(ci.level)||ci.level>=1||ci.level<=0) {
    stop("'ci.level' must be numeric and between 0 and 1")
  }
  if(cumul==TRUE) {
    # SET THE LAG STEP EQUAL TO 1
    x$bylag <- 1
    if(is.null(x$cumfit)) {
      stop("Cumulative outcomes can be plotted if predicted in the 'crosspred'
  object. Set the argument 'cumul=TRUE' in the function crosspred()")
    }
  }
  if(!is.null(exp)&&!is.logical(exp)) stop("'exp' must be logical")
  #

  ###########################
  # COMPUTE OUTCOMES
  #
  # CUMULATIVE IF CUMUL==T
  if(cumul==TRUE) {
    x$matfit <- x$cumfit
  }

  noeff <- 0
  #
  # EXPONENTIAL
  if((is.null(exp)&&!is.null(x$model.link)&&x$model.link%in%c("log","logit"))||
     (!is.null(exp)&&exp==TRUE)) {

    x$matfit <- x$matRRfit
    x$allfit <- x$allRRfit

    noeff <- 1

  }

  ###################
  # GRAPHS

  if(ptype == "overall") {
    # SET DEFAULT VALUES IF NOT INCLUDED BY THE USER
    min_all <- min(sapply(x$allfit, min))
    max_all <- max(sapply(x$allfit, max))
    median_all <- do.call(cbind, lapply(x$allfit, as.data.frame))
    median_all <- apply(median_all, 1, median)

    plot.arg <- list(type="l", ylim=c(min_all, max_all),
                     xlab="Var", ylab="Outcome", bty="l")

    plot.arg <- modifyList(plot.arg,c(list(x=x$predvar,y=x$allfit[[1]], col = "grey80")))
    do.call("plot", plot.arg)

    for(i in 2:length(x$allfit)) {
      plot.arg <- modifyList(plot.arg,c(list(x=x$predvar,y=x$allfit[[i]], col = "grey80")))
      # PLOT
      do.call("lines", plot.arg)
    }

    plot.arg <- modifyList(plot.arg, c(list(x = x$predvar, y = median_all, col = "black")))

    do.call("lines", plot.arg)
  }


}
