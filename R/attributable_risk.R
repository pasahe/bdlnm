#' attributable_risk function
#'
#' @param x
#' @param basis
#' @param cases
#' @param model
#' @param dir
#' @param average
#' @param tot
#' @param cen
#' @param range
#' @param n_sample
#' @param seed
#'
#' @returns
#' @export
#'
#' @examples
attributable_risk <- function(x, basis, cases, model, dir = "back", average = NULL, tot = TRUE, cen, range = NULL, n_sample = 1000, seed = 0L) {

  ################################################################################
  # CHECK

  if(is.null(model)) stop("'model' must be provided")
  # EXTRACT NAME AND CHECK type AND dir
  name <- deparse(substitute(basis))
  dir <- match.arg(dir,c("back","forw"))
  #
  # DEFINE CENTERING
  if(missing(cen) && is.null(cen <- attr(basis,"argvar")$cen))
    stop("'cen' must be provided")
  if(!is.numeric(cen) && length(cen)>1L) stop("'cen' must be a numeric scalar")
  attributes(basis)$argvar$cen <- NULL
  #
  # SELECT RANGE (FORCE TO CENTERING VALUE OTHERWISE, MEANING NULL RISK)
  if(!is.null(range)) x[x<range[1]|x>range[2]] <- cen
  #
  if(length(cases)!=NROW(x)) stop("'x' and 'cases' not consistent")

  if(!is.null(average)) {
    if(!is.logical(average)) stop("'average' must be TRUE/FALSE")
  } else {
    if(dir == "forw") average <- FALSE
  }

  lag <- attr(basis,"lag")
  predlag <- dlnm:::seqlag(lag)

  #
  ################################################################################
  #Calculate the centered RR at each of the values of the daily cases

  cpred <- bcrosspred(basis, model, at = x, cen = cen, n_sample = n_sample, seed = seed)

  # cp <- matrix(cpred$matfit.summary[,"0.5quant"], length(x), length(predlag))
  # cp_rr <- matrix(cpred$matRRfit.summary[,"0.5quant"], length(x), length(predlag))

  cp <- cpred$matfit
  cp_rr <- cpred$matRRfit

  ################################################################################
  #
  # COMPUTE AF AND AN

  #Initialize matrices:
  M_an <- M_af <- matrix(nrow = length(x), ncol = ncol(cp))

  if(!tot) {
    an <- af <- M_an
  } else {
    an <- af <- matrix(nrow = 1, ncol = ncol(cp))
  }

  #Forward perspective
  if(dir == "forw") {

    #Compute the Lagged matrix of daily cases for that day and the next max_lag days
    lagged_cases <- tsModel::Lag(cases, seq(lag[1], -lag[2]))

    if(!average) {

      af_cp <- (cp_rr - 1) / cp_rr

      for(i in 1:ncol(cp)) {
        an_sample <- matrix(af_cp[,i], length(x), length(predlag)) * lagged_cases
        #Hem d'ometre els missings de lagged_cases perquè sino les últimes files donaran missing i no contribuiran al AN
        M_an[,i] <- rowSums(an_sample, na.rm = TRUE)
        #Comprovar que està bé amb Marcos (per una sample i)
        M_af[,i] <- M_an[,i]/rowSums(lagged_cases, na.rm = TRUE)

        #
        # TOTAL
        #   - SELECT NON-MISSING OBS CONTRIBUTING TO COMPUTATION
        #   - DERIVE TOTAL AF
        #   - COMPUTE TOTAL AN WITH ADJUSTED DENOMINATOR (OBSERVED TOTAL NUMBER)
        if(tot) {
          isna <- is.na(M_an[,i])
          an[,i] <- sum(M_an[!isna, i])
          af[,i] <- an[,i]/sum(cases[!isna])
          # an <- af*den #No entenc això del denominador ajustat
        } else {
          an <- M_an
          af <- M_af
        }
      }

    } else {

      for(i in 1:ncol(cp)) {
        cp_sample <- matrix(cp[,i], length(x), length(predlag))
        cp_rr_all <- exp(rowSums(cp_sample))
        M_af[,i] <- (cp_rr_all - 1) / cp_rr_all
        M_an[,i] <- M_af[,i] * rowMeans(lagged_cases, na.rm = TRUE)

        if(tot) {
          isna <- is.na(M_an[,i])
          an[,i] <- sum(M_af[,i] * rowMeans(lagged_cases, na.rm = TRUE), na.rm = TRUE)
          af[,i] <- an[,i]/sum(cases[!isna])
          # an <- af*den #No entenc això del denominador ajustat
        } else {
          an <- M_an
          af <- M_af
        }
      }

    }

  } else if(dir == "back") {

      af_cp <- (cp_rr - 1) / cp_rr

      for(i in 1:ncol(cp)) {

        af_sample <- matrix(af_cp[,i], length(x), length(predlag))

        #Sum of anti-diagonals of AF
        #The key idea is that all anti-diagonal element share the same value of row_index + col_index so we can use tapply
        af_sample <- tapply(af_sample, row(af_sample) + col(af_sample), sum)
        #If we do it in a 3x3 example it's clear that we don't want those row-column elements that sum more than the number of rows
        M_af[,i] <- af_sample[1:length(cases)]

        M_an[,i] <- M_af[,i] * cases

        #
        # TOTAL
        #   - SELECT NON-MISSING OBS CONTRIBUTING TO COMPUTATION
        #   - DERIVE TOTAL AF
        #   - COMPUTE TOTAL AN WITH ADJUSTED DENOMINATOR (OBSERVED TOTAL NUMBER)
        if(tot) {
          isna <- is.na(M_an[,i])
          an[,i] <- sum(M_an[!isna, i])
          af[,i] <- an[,i]/sum(cases[!isna])
          # an <- af*den #No entenc això del denominador ajustat
        } else {
          an <- M_an
          af <- M_af
        }

      }

  }

  an_res <- af_res <-  matrix(nrow = nrow(an), ncol = ncol(cpred$allfit.summary))
  colnames(an_res) <- colnames(af_res) <- colnames(cpred$allfit.summary)

  if(!tot) {
    rownames(an_res) <- rownames(af_res) <- x
  }

  an_res[,"mean"] <- apply(an, 1, mean)
  an_res[,"sd"] <- apply(an, 1, sd)
  quantiles <- as.numeric(gsub("quant", "", colnames(cpred$allfit.summary)[grep("quant", colnames(cpred$allfit.summary))]))
  for(q in quantiles) {
    an_res[, paste0(q, "quant")] <- apply(an, 1, function (i) quantile(i, q))
  }
  an_res[, "mode"] <- apply(an, 1, function (i) unique(i)[which.max(tabulate(match(i, unique(i))))])

  af_res[,"mean"] <- apply(af, 1, mean)
  af_res[,"sd"] <- apply(af, 1, sd)
  quantiles <- as.numeric(gsub("quant", "", colnames(cpred$allfit.summary)[grep("quant", colnames(cpred$allfit.summary))]))
  for(q in quantiles) {
    af_res[, paste0(q, "quant")] <- apply(af, 1, function (i) quantile(i, q))
  }
  af_res[, "mode"] <- apply(af, 1, function (i) unique(i)[which.max(tabulate(match(i, unique(i))))])

  res <- list("af" = af_res |> as.data.frame(), "an" = an_res |> as.data.frame())
  return(res)

}
