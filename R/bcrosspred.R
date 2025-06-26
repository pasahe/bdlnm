
#' bcrosspred
#' This function generates predictions from bayesian distributed lag linear (DLMs) and distributed lag non-linear models (DLNM).
#' @param basis
#' @param model
#' @param coef
#' @param model.link
#' @param at
#' @param from
#' @param to
#' @param by
#' @param lag
#' @param bylag
#' @param cen
#' @param cumul
#' @param n_sample
#'
#' @returns
#' @export
#'
#' @examples

bcrosspred <- function(basis, model = NULL, coef = NULL, model.link = NULL, at = NULL, from = NULL, to = NULL, by = NULL, lag, bylag = 1, cen = NULL, cumul = FALSE, n_sample = NULL) {

  #######################################
  #Determine the type of model and checks

  # TYPE OF PREDICTION: CROSSBASIS or ONEBASIS
  type <- if(any(class(basis)%in%"crossbasis")) "cb" else
    if(any(class(basis)%in%"onebasis")) "one"

  # SET name
  name <- deparse(substitute(basis))

  #  EXTRACT ORIGINAL lag (DEPENDENT ON TYPE)
  origlag <- switch(type,
                    cb = attr(basis,"lag"),
                    one = c(0,0)
  )
  lag <- if(missing(lag)) origlag else dlnm:::mklag(lag)

  # CHECKS ON lag AND bylag
  if(!all(lag==origlag) && cumul)
    stop("cumulative prediction not allowed for lag sub-period")
  lagfun <- switch(type, cb=attr(basis,"arglag")$fun, one=NULL)
  if(bylag!=1L && !is.null(lagfun) && lagfun=="integer")
    stop("prediction for non-integer lags not allowed for type 'integer'")
  #
  # OTHER COHERENCE CHECKS
  if(is.null(model) && (is.null(coef)))
    stop("At least 'model' or 'coef' must be provided")

  if(!is.null(model)) {
    if(any(class(model) != "inla")) stop("'mod' must be an inla model.") else
      if(!model$.args$control.compute$config) stop("'mod' must be an inla model computed with option `control.compute=list(config = TRUE)`.")
    if(!is.null(coef)) stop("'mod' and 'coef' cannot be provided at the same time.")
  } else {

    if(!"latent" %in% names(coef[[1]])) stop("'coef' must be an object created by `inla.posterior.sample()`.")

    if(is.null(model.link)) stop("'model.link' has to be provided if a model is not supplied.")

  }

  ###############################
  # SET COEF AND LINK FOR EVERY TYPE OF MODELS
  if(!is.null(model)) {

    if(is.null(n_sample)) n_sample <- 100

    # WRITE CONDITIONS (DEPENDENT ON TYPE AND IF MATRIX/VECTOR)
    #Regex type condition to find the estimates of the crossbasis coefficients in the model
    cond <- if(ncol(basis)==1L) name else "v[0-9]{1,2}\\.l[0-9]{1,2}"

    #Sample from the posterior distribution to get the models
    names_sel <- grep("^v", rownames(model$summary.fixed), value = TRUE)
    list_sel <- lapply(1:length(names_sel), function(x) 1)
    names(list_sel) <- names_sel

    inla_res <- INLA::inla.posterior.sample(n = n_sample, model, selection = list_sel)
    coef <- do.call(cbind, lapply(inla_res, function (x) x$latent))

    # Define for only some specific classes
    if(is.null(model.link)) {
      model.link <- if(model$.args$family %in% c("poisson", "coxph", "exponential")) "log" else if(model$.args$family %in% c("binomial"))
        "logit" else if(model$.args$family == "gaussian")
          "identity" else
            NA
    }

    model.class <- "inla"

  } else {

    #Get coef as a matrix
    coef <- do.call(cbind, lapply(coef, function (x) x$latent))

    model.class <- "inla"

    if(!is.null(n_sample)) {
      if(n_sample != ncol(coef)) stop("'n_sample' must match the number of samples provided in 'coef'")
    } else {
      n_sample <- ncol(coef)
    }

  }

  #
  # CHECK
  #Number of parameters of the crossbasis
  npar <- ncol(basis)

  if(nrow(coef) != npar || any(is.na(coef))) {
      #It the number of parameters of the original crossbasis and the number of parameters estimated by the model don't match
      stop("coef/vcov not consistent with basis matrix. See help(crosspred)")
  }

  #####################
  # AT, PREDVAR, PREDLAG AND CENTERING

  # RANGE
  #Get range of temperature stored in the attributes of the crossbasis
  range <- attr(basis,"range")

  # SET at, predvar AND predlag
  #If at is not given, it is reconstructed from the other arguments. If at is given it doesn't change anything
  at <- dlnm:::mkat(at,from,to,by,range,lag,bylag)
  #Define the matrix of temperatures and lags in which predictions will be made
  predvar <- if(is.matrix(at)) rownames(at) else at
  predlag <- dlnm:::seqlag(lag,bylag)

  # DEFINE CENTERING VALUE (NULL IF UNCENTERED), AND REMOVE INFO FROM BASIS
  #It is centered in the value: median(pretty(range)) (in my example, 10)
  cen <- dlnm:::mkcen(cen, type, basis, range)
  if(type=="one") attributes(basis)$cen <- NULL
  if(type=="cb") attributes(basis)$argvar$cen <- NULL

  #################################
  # PREDICTION OF LAG-SPECIFIC EFFECTS

  # CREATE THE MATRIX OF TRANSFORMED CENTRED VARIABLES (DEPENDENT ON TYPE)
  #They center the effects previous to estimating the effects (we do it later)
  #This matrix is created with a tensor product and it will be (length(predvar) x length(predlag)) x npar (we do it manually with a loop)
  Xpred <- dlnm:::mkXpred(type,basis,at,predvar,predlag,cen)

  matfit <- Xpred %*% coef
  colnames(matfit) <- outer("sample", seq(1, n_sample), paste, sep="")

  #################################
  # PREDICTION OF OVERALL+CUMULATIVE EFFECTS

  # RE-CREATE LAGGED VALUES (NB: ONLY LAG INTEGERS)
  predlag <- dlnm:::seqlag(lag)

  # CREATE THE MATRIX OF TRANSFORMED VARIABLES (DEPENDENT ON TYPE)
  Xpred <- dlnm:::mkXpred(type,basis,at,predvar,predlag,cen)

  # CREATE OVERALL AND (OPTIONAL) CUMULATIVE EFFECTS AND SE
  Xpredall <- 0

  cumfit <- matrix(0, nrow(matfit), ncol(matfit))

  for(i in seq(length(predlag))) {

    ind <- seq(length(predvar)) + length(predvar)*(i-1)
    Xpredall <- Xpredall + Xpred[ind,,drop=FALSE]

    if(cumul) {
      cumfit[ind,] <- Xpredall %*% coef
      colnames(cumfit) <- outer("sample", seq(1, n_sample), paste, sep="")
    }
  }

  allfit <- Xpredall %*% coef
  rownames(allfit) <- predvar
  colnames(allfit) <- outer("sample", seq(1, n_sample), paste, sep="")

  ########################
  # CREATE THE OBJECT

  # INITIAL LIST, THEN ADD COMPONENTS
  list <- list(predvar=predvar)
  if(!is.null(cen)) list$cen <- cen

  list <- c(list, list(lag=lag, bylag=bylag, coefficients=coef, matfit=matfit, allfit=allfit))

  if(cumul) list <- c(list, list(cumfit=cumfit))

  link.inv <- if(!is.null(model.link) && model.link %in% c("log","logit")) exp else identity

  list$matRRfit <- link.inv(matfit)
  list$allRRfit <- link.inv(allfit)

  if(cumul) list$cumRRfit <- link.inv(cumfit)

  list$model.class <- model.class
  list$model.link <- model.link

  class(list) <- c("bcrosspred", "crosspred")

  return(list)

}
