### Util functions for the R package bdlnm

#' Check if INLA is installed and meets the minimum version requirement
#'
#' Checks whether the \pkg{INLA} package is installed and whether the installed version is >= 23.4.24. This is an internal utility used by functions that
#' depend on \pkg{INLA}.
#'
#' @param error Logical. If `TRUE`, throws an error when \pkg{INLA} is not installed or the version is too old. If `FALSE` (default), returns `FALSE` silently.
#'
#' @return `TRUE` if \pkg{INLA} is installed and meets the version requirement, `FALSE` otherwise (when `error = FALSE`).
#'
#' @keywords internal
#' @export
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
extract_coef <- function(object, basis, combine_effects) {
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

	coef <- object$coefficients

	cb_cols   <- paste0(which_basis, colnames(basis))
	n_cols    <- length(cb_cols)
	n_samples <- ncol(coef)

	# if combine_effects is specified look for an interaction term

	if (!is.null(combine_effects)) {
	
		# parse combine_effects to identify fixed and random cb terms
		fe_terms <- terms(combine_effects, specials = "f")
		fe_vars  <- attr(fe_terms, "variables")
		fe_f_idx <- attr(fe_terms, "specials")$f

		# identify fixed cb term: is which_basis present as a non-f() term?
		fe_labels  <- attr(fe_terms, "term.labels")
		fe_f_labels <- if (!is.null(fe_f_idx)) {
			vapply(fe_f_idx, function(i) deparse(fe_vars[[i + 1]], width.cutoff = 500), character(1))
		} else character(0)
		has_global <- which_basis %in% setdiff(fe_labels, fe_f_labels) && all(cb_cols %in% rownames(coef))

		# identify random effect terms: f(id_var, cb, ...) where cb matches which_basis
		re_id_vars <- character(0)
		if (!is.null(fe_f_idx)) {
			for (i in fe_f_idx) {
				f_call         <- fe_vars[[i + 1]]
				id_var_name    <- deparse(f_call[[2]])
				cb_arg_name    <- deparse(f_call[[3]])
				if (cb_arg_name == which_basis) re_id_vars <- c(re_id_vars, id_var_name)
			}
		}

		# for each random effect id var, find unit ids from coef rownames
		re_prefix <- paste0(re_id_vars, "_", cb_cols[1])   # use first cb col to find units
		re_rows   <- rownames(coef)[grepl(paste(paste0("^", re_prefix, ":"), collapse = "|"), rownames(coef))]
		if (length(re_rows) == 0) {
			cli::cli_abort(c(
				"No random effect coefficients found for {.val {re_id_vars}} and the provided basis.",
				"i" = "Expected rownames like {.val {paste0(re_prefix[1], ':<id>')}} in the model coefficients."
			))
		}
		unit_ids <- unique(gsub(".*:", "", re_rows))
		n_units  <- length(unit_ids)

		# global/mean slice
		if (has_global) {
			global_coef <- coef[cb_cols, , drop = FALSE]
			# build array
			all_ids <- c("global", paste0(re_id_vars[1], "[", unit_ids, "]"))
			coef_cb <- array(
				NA_real_,
				dim      = c(n_cols, n_samples, n_units + 1L),
				dimnames = list(cb_cols, colnames(coef), all_ids)
			)
			coef_cb[, , "global"] <- global_coef
		}  else {
			all_ids <- paste0(re_id_vars[1], "[", unit_ids, "]")
			coef_cb <- array(
				NA_real_,
				dim      = c(n_cols, n_samples, n_units),
				dimnames = list(cb_cols, colnames(coef), all_ids)
			)
		}


		replace_call <- function(expr, target, replacement) {
  		if (identical(expr, target)) return(replacement)
  		if (is.call(expr)) {
    		for (i in seq_along(expr)) expr[[i]] <- replace_call(expr[[i]], target, replacement)
  		}
  		expr
		}

		for (k in unit_ids) {

			k_named  <- paste0(re_id_vars[1], "[", k, "]")
			eval_env <- new.env(parent = parent.frame())
			expr <- combine_effects[[2]]

			# assign fixed cb term -> its coefficient matrix
			if (has_global) {
				assign(which_basis, global_coef, envir = eval_env)
			}

			# assign each f() term -> its coefficient matrix for unit k
			for (i in fe_f_idx) {
				f_call      <- fe_vars[[i + 1]]
				id_var_name <- deparse(f_call[[2]])
				cb_arg_name <- deparse(f_call[[3]])
				if (cb_arg_name != which_basis) next

				re_row_names_k <- paste0(id_var_name, "_", cb_cols, ":", k)
				missing_re     <- setdiff(re_row_names_k, rownames(coef))
				if (length(missing_re) > 0) {
					cli::cli_abort(c(
						"Random effect coefficients not found for {.val {id_var_name}} = {.val {k}}.",
						"x" = "{.val {missing_re}}"
					))
				}
				re_coef_k           <- coef[re_row_names_k, , drop = FALSE]
				rownames(re_coef_k) <- cb_cols

				safe_name <- make.names(deparse(f_call, width.cutoff = 500))
    		assign(safe_name, re_coef_k, envir = eval_env)
				expr_str <- gsub(deparse(f_call, width.cutoff = 500), safe_name, 
                 deparse(expr, width.cutoff = 500), fixed = TRUE)
				expr     <- str2lang(expr_str)
			}

			# evaluate the RHS of combine_effects in that environment
			coef_cb[, , k_named] <- with(as.list(eval_env), eval(expr))
		}
		
	} else {
		# extract only the raw basis coefficients
		coef_names <- paste0(which_basis, colnames(basis))
		missing <- setdiff(coef_names, rownames(coef))
		if (length(missing) > 0) {
			cli::cli_abort(
				c(
					"No basis coefficients found in the model for the provided {.arg basis}",
					"x" = "{.val {missing}}"
				)
			)
		}

		coef_cb <- array(
			NA_real_,
			dim      = c(n_cols, n_samples, 1L),
			dimnames = list(cb_cols, colnames(coef), 1)
		)

		coef_cb[, , 1] <- coef[coef_names, , drop = FALSE]

	}

	return(coef_cb) # array

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
