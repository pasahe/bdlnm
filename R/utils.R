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

# Identify named results stored inside a collection object.
named_result_names <- function(x, required_components) {
	names(x)[vapply(
		x,
		function(element) {
			is.list(element) &&
				all(required_components %in% names(element))
		},
		logical(1L)
	)]
}

# Normalize the orientation of fixed basis-interaction coefficient names.
normalize_basis_interaction_names <- function(coefficient_names, basis) {
	basis_coefficient_names <- unlist(
		Map(
			function(basis_name, basis_object) {
				paste0(basis_name, colnames(basis_object))
			},
			names(basis),
			basis
		),
		use.names = FALSE
	)

	vapply(
		coefficient_names,
		function(coefficient_name) {
			components <- strsplit(
				coefficient_name,
				":",
				fixed = TRUE
			)[[1L]]
			basis_position <- which(components %in% basis_coefficient_names)

			if (
				length(components) >= 2L &&
					length(basis_position) == 1L &&
					basis_position != 1L
			) {
				components <- c(
					components[[basis_position]],
					components[-basis_position]
				)
			}

			paste(components, collapse = ":")
		},
		character(1L),
		USE.NAMES = FALSE
	)
}

# Select one result from a named collection while preserving standalone objects.
select_named_result <- function(
	x,
	select = NULL,
	required_components,
	metadata_names = character(0),
	result_class,
	result_label = "result"
) {
	is_interaction_result <- !is.null(attr(x, "basis_interaction"))
	result_label_sentence <- paste0(
		toupper(substr(result_label, 1L, 1L)),
		substr(result_label, 2L, nchar(result_label))
	)

	if (!is_interaction_result) {
		if (!is.null(select)) {
			cli::cli_abort(
				"{.arg select} can only be used when {.arg x} contains multiple {result_label}s."
			)
		}
		return(x)
	}

	result_names <- named_result_names(x, required_components)
	if (length(result_names) == 0L) {
		cli::cli_abort(
			"Interaction-specific {result_label}s expected in {.arg x} could not be found."
		)
	}

	if (is.null(select)) {
		if (length(result_names) > 1L) {
			cli::cli_abort(c(
				"{.arg select} must identify the {result_label} to use.",
				"i" = "Available {result_label}s: {.and {result_names}}."
			))
		}
		select <- result_names[[1L]]
	} else {
		if (
			!is.character(select) ||
				length(select) != 1L ||
				is.na(select) ||
				!nzchar(select)
		) {
			cli::cli_abort(
				"{.arg select} must be a single non-empty character string."
			)
		}

		if (!select %in% result_names) {
			cli::cli_abort(c(
				"{result_label_sentence} {.val {select}} was not found in {.arg x}.",
				"i" = "Available {result_label}s: {.and {result_names}}."
			))
		}
	}

	selected_result <- x[[select]]
	for (metadata_name in metadata_names) {
		if (is.null(selected_result[[metadata_name]])) {
			selected_result[[metadata_name]] <- x[[metadata_name]]
		}
	}
	class(selected_result) <- result_class
	selected_result
}

#' Combine formula expressions with addition
#'
#' Combines a list of language objects into a single expression joined by `+`.
#' This internal utility is used when expanding basis terms in model formulas.
#'
#' @param expressions A non-empty list of language objects.
#'
#' @return A language object containing the combined expression.
#'
#' @keywords internal
#' @export
sum_formula_expressions <- function(expressions) {
	if (length(expressions) == 1L) {
		return(expressions[[1L]])
	}

	Reduce(
		function(left, right) call("+", left, right),
		expressions
	)
}

#' Expand basis terms in a model formula
#'
#' Expands basis symbols and basis-specific random-effect terms by manipulating
#' the formula language objects. Generated grouping columns required by expanded
#' random-effect terms are added to `data`.
#'
#' @param formula A model formula.
#' @param basis A named list of `crossbasis` or `onebasis` objects.
#' @param data A data frame containing the variables referenced in `formula`.
#'
#' @return A list containing the expanded `formula` and updated `data`.
#'
#' @keywords internal
#' @export
expand_basis_formula <- function(formula, basis, data) {
	formula_terms <- stats::terms(formula, specials = "f", keep.order = TRUE)
	term_labels <- attr(formula_terms, "term.labels")
	expanded_terms <- list()

	for (term_label in term_labels) {
		term_expression <- str2lang(term_label)

		# Expand basis-specific random-effect terms: f(id, basis, ...).
		if (
			is.call(term_expression) &&
			identical(term_expression[[1L]], as.name("f")) &&
			length(term_expression) >= 3L &&
			is.symbol(term_expression[[3L]]) &&
			as.character(term_expression[[3L]]) %in% names(basis)
		) {
			id_expression <- term_expression[[2L]]
			basis_name <- as.character(term_expression[[3L]])

			if (!is.symbol(id_expression)) {
				cli::cli_abort(
					"The first argument of an {.fn f} term involving a basis must be a column name."
				)
			}

			id_name <- as.character(id_expression)

			if (!id_name %in% names(data)) {
				cli::cli_abort(c(
					"Column {.var {id_name}} not found in {.arg data}.",
					"i" = "Found in formula term {.code {deparse1(term_expression)}}.",
					"i" = "Available columns: {.and {names(data)}}."
				))
			}

			basis_columns <- paste0(basis_name, colnames(basis[[basis_name]]))

			for (basis_column in basis_columns) {
				generated_id <- paste0(basis_column, "_", id_name)

				if (generated_id %in% names(data)) {
					cli::cli_abort(
						"Cannot generate column {.var {generated_id}} because it already exists in {.arg data}."
					)
				}

				data[[generated_id]] <- data[[id_name]]

				random_term <- term_expression
				random_term[[2L]] <- as.name(generated_id)
				random_term[[3L]] <- as.name(basis_column)
				expanded_terms[[length(expanded_terms) + 1L]] <- random_term
			}

			next
		}

		# Iteratively substitute every basis occurring in a fixed-effect term.
		term_basis <- intersect(all.vars(term_expression), names(basis))
		term_expansions <- list(term_expression)

		for (basis_name in term_basis) {
			basis_columns <- paste0(basis_name, colnames(basis[[basis_name]]))
			new_expansions <- list()

			for (term_expansion in term_expansions) {
				for (basis_column in basis_columns) {
					replacement <- stats::setNames(list(as.name(basis_column)), basis_name)
					new_expansions[[length(new_expansions) + 1L]] <- do.call(
						substitute,
						list(term_expansion, replacement)
					)
				}
			}

			term_expansions <- new_expansions
		}

		expanded_terms <- c(expanded_terms, term_expansions)
	}

	expanded_rhs <- sum_formula_expressions(expanded_terms)

	if (attr(formula_terms, "intercept") == 0L) {
		expanded_rhs <- call("-", expanded_rhs, 1L)
	}

	expanded_formula <- call("~", formula[[2L]], expanded_rhs)
	class(expanded_formula) <- "formula"
	environment(expanded_formula) <- environment(formula)

	list(
		formula = expanded_formula,
		data = data
	)
}

#' Describe interactions involving a basis
#'
#' Extracts compact metadata for fixed and random interactions involving a
#' cross-basis or one-basis in a model formula.
#'
#' @param formula A model formula.
#' @param basis A named list of basis objects.
#' @param data A data frame containing the variables referenced in `formula`.
#'
#' @return A nested list indexed by basis and interaction-variable names.
#'
#' @keywords internal
basis_interaction_info <- function(formula, basis, data) {
	formula_terms <- stats::terms(formula, specials = "f", keep.order = TRUE)
	term_expressions <- lapply(attr(formula_terms, "term.labels"), str2lang)
	interaction_info <- list()

	variable_values <- function(variable, fixed = FALSE) {
		if (!variable %in% names(data)) {
			cli::cli_abort(
				"Interaction variable {.var {variable}} was not found in {.arg data}."
			)
		}

		x <- data[[variable]]
		if (fixed && (is.factor(x) || is.character(x) || is.logical(x))) {
			x <- if (is.factor(x)) x else factor(x)
			values <- levels(x)
			list(values = values, reference = values[[1L]])
		} else if (fixed) {
			# Continuous fixed interactions do not have categorical values or a
			# reference level.
			list(values = NULL, reference = NULL)
		} else {
			values <- if (is.factor(x)) levels(x) else unique(x[!is.na(x)])
			list(values = values)
		}
	}

	for (term in term_expressions) {
		if (
			is.call(term) &&
			identical(term[[1L]], as.name("f")) &&
			length(term) >= 3L &&
			is.symbol(term[[2L]]) &&
			is.symbol(term[[3L]]) &&
			as.character(term[[3L]]) %in% names(basis)
		) {
			variable <- as.character(term[[2L]])
			basis_name <- as.character(term[[3L]])
			value_info <- variable_values(variable)

			interaction_info[[basis_name]][[variable]] <- list(
				type = "random",
				values = value_info$values
			)
			next
		}

		if (!is.call(term) || !identical(term[[1L]], as.name(":"))) next

		term_variables <- all.vars(term)
		basis_names <- intersect(term_variables, names(basis))
		if (length(basis_names) != 1L) next

		interaction_variables <- setdiff(term_variables, basis_names)
		if (length(interaction_variables) != 1L) next

		basis_name <- basis_names[[1L]]
		variable <- interaction_variables[[1L]]
		value_info <- variable_values(variable, fixed = TRUE)

		interaction_info[[basis_name]][[variable]] <- list(
			type = "fixed",
			values = value_info$values,
			reference = value_info$reference
		)
	}

	interaction_info
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
	
		combine_terms <- stats::terms(
			combine_effects,
			specials = "f",
			keep.order = TRUE
		)
		term_labels <- attr(combine_terms, "term.labels")
		term_expressions <- lapply(term_labels, str2lang)

		global_terms <- vapply(
			term_expressions,
			function(term) is.symbol(term) && identical(as.character(term), which_basis),
			logical(1)
		)
		has_standalone_basis <- any(global_terms)

		random_terms <- logical(length(term_expressions))
		re_id_vars <- character(0)

		for (i in seq_along(term_expressions)) {
			term <- term_expressions[[i]]

			if (
				is.call(term) &&
				identical(term[[1L]], as.name("f")) &&
				length(term) >= 3L
			) {
				if (
					is.symbol(term[[3L]]) &&
					identical(as.character(term[[3L]]), which_basis)
				) {
					if (!is.symbol(term[[2L]])) {
						cli::cli_abort(
							"The first argument of an {.fn f} term in {.arg combine_effects} must be a grouping-variable name."
						)
					}

					random_terms[[i]] <- TRUE
					re_id_vars <- c(re_id_vars, as.character(term[[2L]]))
				} else {
					cli::cli_abort(
						"All {.fn f} terms in {.arg combine_effects} must involve the selected basis {.val {which_basis}}."
					)
				}
			}
		}

		fixed_interaction_terms <- logical(length(term_expressions))
		interaction_expression <- NULL
		interaction_variable <- NULL

		for (i in seq_along(term_expressions)) {
			term <- term_expressions[[i]]

			if (is.call(term) && identical(term[[1L]], as.name(":")) && length(term) == 3L) {
				left_is_basis <- is.symbol(term[[2L]]) &&
					identical(as.character(term[[2L]]), which_basis)
				right_is_basis <- is.symbol(term[[3L]]) &&
					identical(as.character(term[[3L]]), which_basis)

				if (left_is_basis || right_is_basis) {
					fixed_interaction_terms[[i]] <- TRUE
					interaction_expression <- if (left_is_basis) term[[3L]] else term[[2L]]
					interaction_variables <- all.vars(interaction_expression)

					if (length(interaction_variables) != 1L) {
						cli::cli_abort(
							"A fixed interaction in {.arg combine_effects} must involve the selected basis and one variable."
						)
					}

					interaction_variable <- interaction_variables[[1L]]
				}
			}
		}

		if (sum(random_terms) > 1L) {
			cli::cli_abort(
				"Only one basis-specific {.fn f} term is currently supported in {.arg combine_effects}."
			)
		}

		if (sum(fixed_interaction_terms) > 1L) {
			cli::cli_abort(
				"Only one fixed interaction is currently supported in {.arg combine_effects}."
			)
		}

		has_random <- any(random_terms)
		has_fixed_interaction <- any(fixed_interaction_terms)

		if (has_random && has_fixed_interaction) {
			cli::cli_abort(
				"Fixed and random basis interactions cannot currently be combined in {.arg combine_effects}."
			)
		}

		allowed_terms <- global_terms | random_terms | fixed_interaction_terms

		if (has_fixed_interaction) {
			interaction_main_term <- vapply(
				term_expressions,
				function(term) identical(term, interaction_expression),
				logical(1)
			)
			allowed_terms <- allowed_terms | interaction_main_term
		}

		if (any(!allowed_terms)) {
			cli::cli_abort(c(
				"The {.arg combine_effects} formula contains unsupported terms.",
				"x" = "Unsupported: {.val {term_labels[!allowed_terms]}}"
			))
		}

		if (!has_standalone_basis && !has_random && !has_fixed_interaction) {
			cli::cli_abort(
				"The selected basis {.val {which_basis}} was not found in {.arg combine_effects}."
			)
		}

		if (has_standalone_basis) {
			missing_global <- setdiff(cb_cols, rownames(coef))

			if (length(missing_global) > 0L) {
				cli::cli_abort(c(
          "Global basis coefficients were requested in {.arg combine_effects}, but were not found for {.val {which_basis}}.",
          "x" = "Missing: {.val {missing_global}}",
          "i" = "Remove the global {.val {which_basis}} term from {.arg combine_effects} if the fitted model only contains basis-specific effects."
        ))
			}
		}
		

		if (has_random) {

			# for each random effect id var, find unit ids from coef rownames
			re_prefix <- paste0(cb_cols[1], "_", re_id_vars[[1L]], ":")
			re_rows <- rownames(coef)[startsWith(rownames(coef), re_prefix)]
			if (length(re_rows) == 0) {
				cli::cli_abort(c(
					"No random effect coefficients found for {.val {re_id_vars}} and the provided basis.",
					"i" = "Expected rownames starting with {.val {re_prefix}} in the model coefficients."
				))
			}
			unit_ids <- unique(substring(re_rows, nchar(re_prefix) + 1L))
			n_units  <- length(unit_ids)

			# global/mean slice
			if (has_standalone_basis) {
				global_coef <- coef[cb_cols, , drop = FALSE]
				# build array
				all_ids <- c("global", paste0(re_id_vars[[1L]], "[", unit_ids, "]"))
				coef_cb <- array(
					NA_real_,
					dim      = c(n_cols, n_samples, n_units + 1L),
					dimnames = list(cb_cols, colnames(coef), all_ids)
				)
				coef_cb[, , "global"] <- global_coef
			}  else {
				all_ids <- paste0(re_id_vars[[1L]], "[", unit_ids, "]")
				coef_cb <- array(
					NA_real_,
					dim      = c(n_cols, n_samples, n_units),
					dimnames = list(cb_cols, colnames(coef), all_ids)
				)
			}

			for (k in unit_ids) {
				k_named <- paste0(re_id_vars[[1L]], "[", k, "]")
				re_row_names_k <- paste0(cb_cols, "_", re_id_vars[[1L]], ":", k)
				missing_re <- setdiff(re_row_names_k, rownames(coef))

				if (length(missing_re) > 0) {
					cli::cli_abort(c(
						"Random effect coefficients not found for {.val {re_id_vars[[1L]]}} = {.val {k}}.",
						"x" = "{.val {missing_re}}"
					))
				}

				re_coef_k <- coef[re_row_names_k, , drop = FALSE]
				rownames(re_coef_k) <- cb_cols

				if (has_standalone_basis) {
					re_coef_k <- global_coef + re_coef_k
				}

				coef_cb[, , k_named] <- re_coef_k
			}
		} else if (has_fixed_interaction) {
			interaction_prefix <- deparse1(interaction_expression)
			int_levels <- character(0)
			has_continuous_interaction <- FALSE

			for (coefficient_name in rownames(coef)) {
				name_parts <- strsplit(coefficient_name, ":", fixed = TRUE)[[1L]]

				if (length(name_parts) != 2L) next

				if (identical(name_parts[[1L]], cb_cols[[1L]])) {
					interaction_component <- name_parts[[2L]]
				} else if (identical(name_parts[[2L]], cb_cols[[1L]])) {
					interaction_component <- name_parts[[1L]]
				} else {
					next
				}

				if (startsWith(interaction_component, interaction_prefix)) {
					level <- substring(
						interaction_component,
						nchar(interaction_prefix) + 1L
					)

					if (nzchar(level)) {
						int_levels <- c(int_levels, level)
					} else {
						has_continuous_interaction <- TRUE
					}
				}
			}

			int_levels <- unique(int_levels)

			if (length(int_levels) == 0L && has_continuous_interaction) {
				cli::cli_abort(c(
					"A continuous interaction was found for {.val {which_basis}} and {.val {interaction_variable}}.",
					"i" = "Continuous interactions are not currently supported by {.arg combine_effects} because no value of {.val {interaction_variable}} is specified for prediction.",
					"i" = "If {.val {interaction_variable}} represents categories, include it as a factor in the fitted model."
				))
			}

			if (length(int_levels) == 0L) {
				cli::cli_abort(c(
					"No interaction coefficients were found for {.val {which_basis}} and {.val {interaction_variable}}.",
					"i" = "Make sure the fitted model contains an interaction between these terms, or remove the interaction from {.arg combine_effects}."
				))
			}
			
			if (has_standalone_basis) {
				interaction_metadata <- attr(object, "basis_interaction")[[which_basis]][[interaction_variable]]
				reference_level <- interaction_metadata$reference
				if (is.null(reference_level)) {
					reference_level <- "reference"
				}
				reference_id <- paste0(interaction_variable, "[", reference_level, "]")

				all_ids <- c(
					reference_id,
					paste0(interaction_variable, "[", int_levels, "]")
				)
			} else {
				all_ids <- paste0(interaction_variable, "[", int_levels, "]")
			}

			n_levels <- length(all_ids)

			coef_cb <- array(
				NA_real_,
				dim      = c(n_cols, n_samples, n_levels),
				dimnames = list(cb_cols, colnames(coef), all_ids)
			)

			# reference level: just global cb coefficients
			if (has_standalone_basis) {
				global_coef <- coef[cb_cols, , drop = FALSE]
				coef_cb[, , reference_id] <- global_coef
			}

			# other levels: global cb + interaction coefficients
			for (lv in int_levels) {
				interaction_component <- paste0(interaction_prefix, lv)
				int_row_names <- character(length(cb_cols))

				for (i in seq_along(cb_cols)) {
					basis_first <- paste0(cb_cols[[i]], ":", interaction_component)
					basis_second <- paste0(interaction_component, ":", cb_cols[[i]])

					if (basis_first %in% rownames(coef)) {
						int_row_names[[i]] <- basis_first
					} else if (basis_second %in% rownames(coef)) {
						int_row_names[[i]] <- basis_second
					} else {
						cli::cli_abort(
							"Interaction coefficient not found for {.val {cb_cols[[i]]}} and level {.val {lv}}."
						)
					}
				}

				int_coef <- coef[int_row_names, , drop = FALSE]
				rownames(int_coef) <- cb_cols

				if (has_standalone_basis) {
					int_coef <- global_coef + int_coef
				}

				coef_cb[, , paste0(interaction_variable, "[", lv, "]")] <- int_coef
			}

		} else {

			#simple case: just extract global cb coefficients
			if (!identical(term_labels, which_basis)) {
				cli::cli_abort(c(
					"The {.arg combine_effects} formula contains terms that are not recognized.",
					"i" = "Expected only {.val {which_basis}} or interaction terms with it."
				))
			}

			coef_cb <- array(
				coef[cb_cols, , drop = FALSE],
				dim = c(n_cols, n_samples, 1L),
				dimnames = list(cb_cols, colnames(coef), "global")
			)

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
