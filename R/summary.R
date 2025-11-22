#' summary methods for bdlnm objects
#'
#' @name summary_bdlnm
#' @keywords internal
#' @param object An object created using bdlnm functions
#' @param ... Not used
NULL

#' @rdname summary_bdlnm
#' @export
summary.bdlnm <- function(object, ...) {

  cat("Fitted INLA model:\n\n")

  print(object$model)

  cat("\n")

  cat("Posterior sampled coefficients:\n\n")

  if(ncol(object$coefficients) > 5) {

    print(object$coefficients[, seq_len(5)], row.names = TRUE)
    cat("... (showing first 5 of", ncol(object$coefficients), "posterior samples)\n")

  } else {

    print(object$coefficients, row.names = TRUE)

  }

  cat("\n")

  cat("Summary of posterior sampled coefficients:\n\n")

  print(object$coefficients.summary, row.names = TRUE)

  cat("\n")

  invisible(object)
}


#' @rdname summary_bdlnm
#' @export
summary.bcrosspred <- function(object, ...) {

  cat("Summary of prediction characteristics:\n\n")
  cat(" Number of basis coefficients:", nrow(object$coefficients), "\n")
  cat(" Number of posterior samples:", ncol(object$coefficients), "\n")
  cat(" Exposure range:", paste(range(object$predvar), collapse = ", "), paste0("(",length(object$predvar), " values)\n"))
  cat(" Centered at:", object$cen, "\n")
  cat(" Lag range:", paste(object$lag, collapse = ", "), paste0("(by = ", object$bylag, ")\n"))

  invisible(object)
}

#' @rdname summary_bdlnm
#' @export
summary.optimal_exposure <- function(object, ...) {

  which_text <- switch(attr(object, "which"),
                       "min" = "Minimum",
                       "max" = "Maximum"
  )

  cat("Optimal criterion:", which_text, "\n\n")

  cat("Posterior sampling optimal exposures:\n\n")

  if(length(object$est) > 5) {

    n_sim <- ifelse(length(object$est) > 5, 5, length(object$est))

    print(object$est[seq_len(n_sim)])
    cat("... (showing first 5 of", length(object$est), "posterior samples)")

  } else {

    print(object$est)

  }

  cat("\n\n")

  cat("Summary of posterior sampling optimal exposures:\n\n")

  print(object$summary)

  cat("\n")

  invisible(object)

}


