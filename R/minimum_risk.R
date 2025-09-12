#' minimum_risk: function to obtain the minimum exposure value in which there is the minimum overall risk effect.
#'
#' @param basis
#' @param model
#' @param at
#' @param from
#' @param to
#' @param by
#' @param n_sample
#' @param seed
#'
#' @returns
#' @export
#'
#' @examples
minimum_risk <- function(basis, model = NULL, at = NULL, from = NULL, to = NULL, by = NULL, n_sample = 1000, seed = 0L) {

  ################################################################################
  # CREATE THE BASIS AND EXTRACT COEF-VCOV
  #
  # CHECK AND DEFINE BASIS
  if(!any(class(basis) == "crossbasis"))
    stop("the first argument must be an object of class 'crossbasis'")
  #
  # INFO
  attr <- attributes(basis)
  range <- attr(basis,"range")
  if(is.null(by)) by <- 0.1
  lag <- attr(basis,"lag")
  if(is.null(model)) stop("'model' must be provided")
  name <- deparse(substitute(basis))

  predvar <- if(is.matrix(at)) rownames(at) else at
  predlag <- dlnm:::seqlag(lag,by=1)

  #
  # PREDICT
  cpred <- suppressMessages(bcrosspred(basis, model, at = at, n_sample = n_sample, seed = seed))

  ################################################################################
  # FIND THE MINIMUM
  #
  ind <- apply(cpred$allfit, 2, which.min)
  min <- predvar[ind]
  names(min) <- NULL

  # SUMMARY.FIXED
  min.summary <- matrix(nrow = 1, ncol = ncol(cpred$allfit.summary))
  colnames(min.summary) <- colnames(cpred$allfit.summary)
  min.summary[,"mean"] <- mean(min)
  min.summary[,"sd"] <- sd(min)
  quantiles <- as.numeric(gsub("quant", "", colnames(cpred$allfit.summary)[grep("quant", colnames(cpred$allfit.summary))]))
  for(q in quantiles) {
    min.summary[, paste0(q, "quant")] <- quantile(min, q)
  }
  min.summary[, "mode"] <- unique(min)[which.max(tabulate(match(min, unique(min))))]

  #
  res <- list(res = min, res.summary = min.summary |> as.data.frame())
  return(res)
}
