
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
plot.bcrosspred <- function(x, ptype, var=NULL, lag=NULL, ci="area", ci.arg, ci.level = x$ci.level, cumul=FALSE, exp=NULL, ...) {

  ######################
  #CHECK

  #Check if it matches the possible choices for ci
  ci <- match.arg(ci, c("area","bars","lines","n", "sampling"))

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
  if(!is.null(lag)&&sum(lag%in%dlnm:::seqlag(x$lag,x$bylag))!=length(lag)&&(ptype=="slices")) {
    stop("'lag' must match values used for prediction")
  }
  if(missing(ci.arg)) {
    ci.arg <- list()
  } else if(!is.list(ci.arg)) stop("'ci.arg' must be a list")
  if(!is.numeric(ci.level)||ci.level>=1||ci.level<=0) {
    stop("'ci.level' must be numeric and between 0 and 1")
  }

  # CHECK IF LOW AND HIGH CREDIBLE INTERVALS ARE ALREADY CALCULATED IN SUMMARY ELEMENT
  if(ci.level != x$ci.level) {
    ci_compute <- TRUE
  } else {
    ci_compute <- FALSE
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
    x$matfit.summary <- x$matRRfit.summary
    x$allfit.summary <- x$allRRfit.summary

    noeff <- 1

  }

  ###################
  # GRAPHS

  # SLICES
  if(ptype=="slices") {

    # SET FRAME AND GREY SCALE
    mar.old <- par()$mar
    mfrow.old <- par()$mfrow
    mgp.old <- par()$mgp
    grey <- grey(0.9)
    if(length(lag)+length(var)>1) {
      layout(matrix(1:(length(var)+length(lag)),
                    ncol=sum(!is.null(var),!is.null(lag))))
      grey <- grey(0.8)
      par(mgp=c(2,0.7,0),mar=c(4.1,4.1,2.1,1.1))
    }

    predlag <- dlnm:::seqlag(x$lag, x$bylag)

    #Get the median coefficients
    median <- matrix(x$matfit.summary[, "0.5quant"], length(x$predvar), length(predlag))

    if(!ci_compute) {
      quantiles <- grep("quant$", colnames(x$matfit.summary))
      low <- matrix(x$matfit.summary[, quantiles[1]], length(x$predvar), length(predlag))
      high <- matrix(x$matfit.summary[, quantiles[3]], length(x$predvar), length(predlag))
    } else {
      quantiles <- c((1 - ci.level)/2, 1 - (1 - ci.level)/2)
      low <- matrix(sapply(1:nrow(x$matfit), function (i) quantile(x$matfit[i,], quantiles[1])), length(x$predvar), length(predlag))
      high <- matrix(sapply(1:nrow(x$matfit), function (i) quantile(x$matfit[i,], quantiles[2])), length(x$predvar), length(predlag))
    }

    # LAG
    if(!is.null(lag)) {
      # START LOOP FOR LAG
      for(i in lag) {

        #Find effect for the corresponding lags
        ind <- seq(length(x$predvar)*(which(predlag == i) - 1) + 1, length(x$predvar)*(which(predlag == i)), by = 1)

        # SET DEFAULT VALUES IF NOT INCLUDED BY THE USER
        plot.arg <- list(type="l",xlab="Var",ylab="Outcome",
                         ylim=c(min(x$matfit[ind,]), max(x$matfit[ind,])),bty="l")

        if(length(lag)+length(var)>1)  plot.arg$cex.axis <- 0.7

        #Get plot arguments supplied to the function
        plot.arg <- modifyList(plot.arg,list(...))

        # SET CONFIDENCE INTERVALS
        smatfit <- matrix(x$matfit[ind,], length(x$predvar), ncol(x$matfit))

        ci.list <- list(panel.first=call("fci",ci=ci,x=x$predvar, y=smatfit,
                                         high=high[,which(predlag == i)],low=low[,which(predlag == i)],ci.arg,plot.arg,
                                         noeff=noeff))

        plot.arg <- modifyList(plot.arg,c(ci.list,
                                          list(x = x$predvar, y = median[,which(predlag == i)])))

        col <- plot.arg$col

        if(length(lag)+length(var)>1) {
          plot.arg$main <- ""
          plot.arg$xlab <- "Var"
        }

        # PLOT
        do.call("plot", plot.arg)

        if(length(lag)>1) mtext(paste("Lag =",i), cex=0.8)

      }
    }

    # VAR
    if(!is.null(var)) {
      # START LOOP FOR VAR
      for(i in var) {

        #Find effect for the corresponding temperature values
        ind <- which(x$predvar == i) + 119*(0:(length(predlag) - 1))

        # SET DEFAULT VALUES IF NOT INCLUDED BY THE USER
        plot.arg <- list(type="l",xlab="Lag",ylab="Outcome",
                         ylim=c(min(x$matfit[ind,]),max(x$matfit[ind,])),bty="l")

        if(length(lag)+length(var)>1)  plot.arg$cex.axis <- 0.7

        #Get plot arguments supplied to the function
        plot.arg <- modifyList(plot.arg, list(...))

        # SET CONFIDENCE INTERVALS
        smatfit <- matrix(x$matfit[ind,], length(predlag), ncol(x$matfit))

        ci.list <- list(panel.first=call("fci",ci=ci,x=predlag, y=smatfit,
                                         high=high[which(x$predvar == i),],low=low[which(x$predvar == i),],ci.arg,plot.arg,
                                         noeff=noeff))

        plot.arg <- modifyList(plot.arg,c(ci.list,
                                          list(x = predlag, y = median[which(x$predvar == i),])))


        if(length(lag)+length(var)>1) {
          plot.arg$main <- ""
          plot.arg$xlab <- "Lag"
        }

        # PLOT
        do.call("plot",plot.arg)

        if(length(lag)>1) mtext(paste("Var =",var[i]),cex=0.8)

      }
    }
    if(length(lag)+length(var)>1) {
      par(mar=mar.old,mfrow=mfrow.old,mgp=mgp.old)
    }
  }

  # OVERALL
  if(ptype == "overall") {
    # SET DEFAULT VALUES IF NOT INCLUDED BY THE USER
    min_all <- min(x$allfit)
    max_all <- max(x$allfit)

    plot.arg <- list(type="l", ylim=c(min_all, max_all),
                     xlab="Var", ylab="Outcome", bty="l")

    #Get plot arguments supplied to the function
    plot.arg <- modifyList(plot.arg,list(...))

    col <- plot.arg$col

    median <- x$allfit.summary[, "0.5quant"]

    if(!ci_compute) {
      quantiles <- grep("quant$", colnames(x$allfit.summary))
      low <- x$allfit.summary[, quantiles[1]]
      high <- x$allfit.summary[, quantiles[3]]
    } else {
      quantiles <- c((1 - ci.level)/2, 1 - (1 - ci.level)/2)
      low <- sapply(1:nrow(x$allfit), function (i) quantile(x$allfit[i,], quantiles[1]))
      high <- sapply(1:nrow(x$allfit), function (i) quantile(x$allfit[i,], quantiles[2]))
    }

    # SET CONFIDENCE INTERVALS
    ci.list <- list(panel.first=call("fci",ci=ci,x=x$predvar,y=x$allfit,
                                     high=high,low=low,ci.arg,plot.arg,noeff=noeff))
    plot.arg <- modifyList(plot.arg,c(ci.list,
                                      list(x=x$predvar, y=median)))

    # PLOT
    do.call("plot",plot.arg)

  }

  #3D
  if(ptype=="3d") {

    if(diff(x$lag)==0) stop("3D plot not conceivable for unlagged associations")

    predlag <- dlnm:::seqlag(x$lag, x$bylag)

    #Get the median coefficients
    summaryfit <- matrix(x$matfit.summary[, "0.5quant"], length(x$predvar), length(predlag))
    rownames(summaryfit) <- x$predvar
    colnames(summaryfit) <- outer("lag", predlag, paste, sep="")

    # SET DEFAULT VALUES IF NOT INCLUDED BY THE USER
    plot.arg <- list(ticktype="detailed",theta=210,phi=30,xlab="Var",
                     ylab="Lag",zlab="Outcome",col="lightskyblue",
                     zlim=c(min(summaryfit),max(summaryfit)),ltheta=290,shade=0.75,r=sqrt(3),d=5)

    #Get plot arguments supplied to the function
    plot.arg <- modifyList(plot.arg,list(...))

    plot.arg <- modifyList(plot.arg,list(x=x$predvar, y=predlag,
                                         z=summaryfit))

    # PLOT
    do.call("persp",plot.arg)
  }

  # CONTOURPLOT
  if(ptype=="contour") {
    if(x$lag[2]==0) stop("contour plot not conceivable for unlagged associations")

    predlag <- dlnm:::seqlag(x$lag, x$bylag)

    #Get the median coefficients
    summaryfit <- matrix(x$matfit.summary[, "0.5quant"], length(x$predvar), length(predlag))
    rownames(summaryfit) <- x$predvar
    colnames(summaryfit) <- outer("lag", predlag, paste, sep="")

    # SET DEFAULT VALUES (NOT TO BE SPECIFIED BY THE USER)
    levels <- pretty(summaryfit, 20)
    col1 <- colorRampPalette(c("blue","white"))
    col2 <- colorRampPalette(c("white","red"))
    col <- c(col1(sum(levels<noeff)),col2(sum(levels>noeff)))
    filled.contour(x=x$predvar,y=dlnm:::seqlag(x$lag,x$bylag),z=summaryfit,col=col,
                   levels=levels,...)
  }

}
