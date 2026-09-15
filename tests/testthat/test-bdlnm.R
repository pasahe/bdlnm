test_that("ensure bdlnm returned expected structure (crossbasis example)", {
  skip_on_cran()
  skip_if_not(check_inla(), "INLA not available")

  # mod build in helper.R
  expect_type(mod, "list")
  expect_equal(length(mod), 4L)
  expect_equal(attr(mod, "n_sim"), n_sim)
  expect_equal(class(mod$model), "inla")
  expect_equal(class(mod$coefficients), c("matrix", "array"))
  expect_equal(dim(mod$coefficients), c(35, n_sim))
  expect_equal(class(mod$coefficients.summary), c("matrix", "array"))
  expect_equal(dim(mod$coefficients.summary), c(35, 6L))
})

test_that("ensure bdlnm returned expected structure (onebasis example)", {
  skip_on_cran()
  skip_if_not(check_inla(), "INLA not available")

  ob <- dlnm::onebasis(slondon$tmean, "strata", breaks = c(5, 10, 20))
  mod_2 <- suppressWarnings(
    bdlnm(
      mort_75plus ~ ob + factor(dow) + seas,
      data = slondon,
      family = "poisson",
      sample.arg = list(n = n_sim, seed = 1L)
    )
  )
  expect_type(mod_2, "list")
  expect_equal(length(mod_2), 4L)
  expect_equal(attr(mod_2, "n_sim"), n_sim)
  expect_equal(class(mod_2$model), "inla")
  expect_equal(class(mod_2$coefficients), c("matrix", "array"))
  expect_equal(dim(mod_2$coefficients), c(18, n_sim))
  expect_equal(class(mod_2$coefficients.summary), c("matrix", "array"))
  expect_equal(dim(mod_2$coefficients.summary), c(18, 6L))
})


test_that("works with two different basis", {
  skip_on_cran()
  skip_if_not(check_inla(), "INLA not available")

  ob <- dlnm::onebasis(slondon$tmean, "strata", breaks = c(5, 10, 20))
  mod_2 <- suppressWarnings(
    bdlnm(
      mort_75plus ~ cb + ob + factor(dow) + seas,
      data = slondon,
      family = "poisson",
      sample.arg = list(n = n_sim, seed = 1L)
    )
  )
  expect_type(mod_2, "list")
  expect_length(mod_2$basis, 2)
  expect_equal(names(mod_2$basis), c("cb", "ob"))
  expect_equal(mod_2$basis[["cb"]], cb)
  expect_equal(mod_2$basis[["ob"]], ob)
})

test_that("bdlnm errors when required arguments are missing or inappropiate", {
  skip_on_cran()
  skip_if_not(check_inla(), "INLA not available")

  expect_snapshot_error(
    bdlnm(mort_75plus ~ cb + factor(dow) + seas, data = slondon, sample.arg = 5)
  )
})

test_that("bdlnm honors sample.arg", {
  skip_on_cran()
  skip_if_not(check_inla(), "INLA not available")

  mod2 <- suppressWarnings(
    bdlnm(
      mort_75plus ~ cb + factor(dow) + seas,
      data = slondon,
      family = "poisson",
      sample.arg = list(n = 5L, seed = 1L)
    )
  )

  expect_equal(attr(mod2, "n_sim"), 5L)
  expect_equal(ncol(mod2$coefficients), 5L)
})

test_that("bdlnm throws informative error if control.compute config = FALSE is passed", {
  skip_on_cran()
  skip_if_not(check_inla(), "INLA not available")

  expect_snapshot_error(
    bdlnm(
      mort_75plus ~ cb + factor(dow) + seas,
      data = slondon,
      family = "poisson",
      control.compute = list(config = FALSE)
    )
  )
})

test_that("bdlnm returns sensible summaries (means, sd, quantiles, mode)", {
  skip_on_cran()
  skip_if_not(check_inla(), "INLA not available")

  coefsum <- mod$coefficients.summary

  # check summary columns include mean, sd and at least one quantile and mode
  expect_true(all(c("mean", "sd", "mode") %in% colnames(coefsum)))
  expect_true(any(grepl("quant$", colnames(coefsum))))
})

test_that("na.action = na.omit drops NA rows and returns consistent dimensions, and gives same coefficients as complete data model", {
  skip_on_cran()
  skip_if_not(check_inla(), "INLA not available")

  expect_equal(nrow(mod$model$model.matrix), nrow(na.omit(cb)))

  na_rows <- which(!complete.cases(cb))
  cb_clean <- cb[-na_rows, ]
  class(cb_clean) <- "crossbasis"
  seas_clean <- seas[-na_rows, ]
  slondon_clean <- slondon[-na_rows, ]

  mod_clean <- suppressWarnings(
    bdlnm(
      mort_75plus ~ cb_clean + factor(dow) + seas_clean,
      data = slondon_clean,
      family = "poisson",
      sample.arg = list(n = n_sim, seed = 1L)
    )
  )

  expect_identical(
    round(as.numeric(mod_clean$coefficients), 2),
    round(as.numeric(mod$coefficients), 2)
  )
})

test_that("na.action = na.pass keeps NA rows and returns consistent dimensions", {
  skip_on_cran()
  skip_if_not(check_inla(), "INLA not available")

  mod_na <- suppressWarnings(
    bdlnm(
      mort_75plus ~ cb + factor(dow) + seas,
      data = slondon,
      family = "poisson",
      na.action = na.pass,
      sample.arg = list(n = n_sim, seed = 1L)
    )
  )

  expect_equal(nrow(mod_na$model$model.matrix), nrow(slondon))
})

test_that("na.action = na.fail throws an error when NAs are present", {
  skip_on_cran()
  skip_if_not(check_inla(), "INLA not available")

  expect_error(
    bdlnm(
      mort_75plus ~ cb + factor(dow) + seas,
      data = slondon,
      family = "poisson",
      na.action = na.fail,
      sample.arg = list(n = n_sim)
    )
  )
})

test_that("na.action is ignored when a random effect is included", {
  skip_on_cran()
  skip_if_not(check_inla(), "INLA not available")

  slondon$id <- seq_len(nrow(slondon))
  expect_message(
    suppressWarnings(
      mod_rt <- bdlnm(
        mort_75plus ~ cb + factor(dow) + seas + f(id, model = "iid"),
        data = slondon,
        family = "poisson",
        sample.arg = list(n = n_sim, seed = 1L)
      )
    )
  )
  expect_equal(nrow(mod_rt$model$model.matrix), nrow(slondon))
})

test_that("bdlnm expands a fixed basis interaction when random effects are present", {
  skip_on_cran()
  skip_if_not(check_inla(), "INLA not available")

  interaction_data <- slondon
  interaction_data$group <- factor(rep(c("a", "b"), length.out = nrow(interaction_data)))
  interaction_data$id_group <- rep(seq_len(5), length.out = nrow(interaction_data))

  mod_interaction <- suppressWarnings(
    bdlnm(
      mort_75plus ~ factor(group) * cb + factor(dow) + seas +
        f(id_group, model = "iid"),
      data = interaction_data,
      family = "poisson",
      sample.arg = list(n = 2L, seed = 1L)
    )
  )

  coefficient_names <- rownames(mod_interaction$coefficients)
  basis_names <- paste0("cb", colnames(cb))

  fitted_terms <- terms(
    mod_interaction$model$.args$formula,
    specials = "f",
    keep.order = TRUE
  )
  expected_terms <- c(
    "factor(group)",
    basis_names,
    paste0("factor(group):", basis_names),
    "factor(dow)",
    "seas",
    'f(id_group, model = "iid")'
  )

  expect_identical(deparse1(attr(fitted_terms, "variables")[[2L]]), "mort_75plus")
  expect_identical(attr(fitted_terms, "intercept"), 1L)
  expect_identical(attr(fitted_terms, "term.labels"), expected_terms)

  expect_true(all(basis_names %in% coefficient_names))
  expect_true(all(
    paste0(basis_names, ":factor(group)b") %in% coefficient_names
  ))
  expect_false(any(
    paste0("factor(group)b:", basis_names) %in% coefficient_names
  ))
  expect_true(any(grepl("^id_group:", coefficient_names)))
})

test_that("bdlnm expands basis-specific random effects", {
  skip_on_cran()
  skip_if_not(check_inla(), "INLA not available")

  random_data <- slondon
  random_data$id_group <- rep(seq_len(5), length.out = nrow(random_data))

  mod_random <- suppressWarnings(
    bdlnm(
      mort_75plus ~ cb + factor(dow) + seas +
        f(id_group, cb, model = "iid"),
      data = random_data,
      family = "poisson",
      sample.arg = list(n = 2L, seed = 1L)
    )
  )

  coefficient_names <- rownames(mod_random$coefficients)
  basis_names <- paste0("cb", colnames(cb))

  fitted_terms <- terms(
    mod_random$model$.args$formula,
    specials = "f",
    keep.order = TRUE
  )
  expected_terms <- c(
    basis_names,
    "factor(dow)",
    "seas",
    paste0(
      "f(", basis_names, "_id_group, ", basis_names,
      ', model = "iid")'
    )
  )

  expect_identical(deparse1(attr(fitted_terms, "variables")[[2L]]), "mort_75plus")
  expect_identical(attr(fitted_terms, "intercept"), 1L)
  expect_identical(attr(fitted_terms, "term.labels"), expected_terms)

  expect_true(all(basis_names %in% coefficient_names))

  for (basis_name in basis_names) {
    expect_true(any(grepl(
      paste0("^", basis_name, "_id_group:"),
      coefficient_names
    )))
  }
})

test_that("bdlnm stores metadata for basis interactions", {
  skip_on_cran()
  skip_if_not(check_inla(), "INLA not available")

  fixed_info <- attr(mod_interaction, "basis_interaction")$cb$group
  expect_identical(fixed_info$type, "fixed")
  expect_identical(fixed_info$values, c("a", "b"))
  expect_identical(fixed_info$reference, "a")

  random_info <- attr(mod_random, "basis_interaction")$cb$id_group
  expect_identical(random_info$type, "random")
  expect_identical(random_info$values, seq_len(5))
  expect_false("reference" %in% names(random_info))

  expect_null(attr(mod, "basis_interaction"))
})

test_that("bdlnm expands crossbasis and onebasis random effects together", {
  skip_on_cran()
  skip_if_not(check_inla(), "INLA not available")

  two_basis_data <- slondon
  two_basis_data$id_group <- rep(
    seq_len(5),
    length.out = nrow(two_basis_data)
  )
  ob <- dlnm::onebasis(
    two_basis_data$tmean,
    "strata",
    breaks = c(5, 10, 20)
  )

  mod_two_basis <- suppressWarnings(
    bdlnm(
      mort_75plus ~ cb + ob +
        f(id_group, cb, model = "iid") +
        f(id_group, ob, model = "iid") +
        factor(dow) + seas,
      data = two_basis_data,
      family = "poisson",
      sample.arg = list(n = 2L, seed = 1L)
    )
  )

  expect_identical(names(mod_two_basis$basis), c("cb", "ob"))

  cb_names <- paste0("cb", colnames(cb))
  ob_names <- paste0("ob", colnames(ob))
  expected_terms <- c(
    cb_names,
    ob_names,
    paste0(
      "f(", cb_names, "_id_group, ", cb_names,
      ', model = "iid")'
    ),
    paste0(
      "f(", ob_names, "_id_group, ", ob_names,
      ', model = "iid")'
    ),
    "factor(dow)",
    "seas"
  )
  fitted_terms <- terms(
    mod_two_basis$model$.args$formula,
    specials = "f",
    keep.order = TRUE
  )

  expect_identical(deparse1(attr(fitted_terms, "variables")[[2L]]), "mort_75plus")
  expect_identical(attr(fitted_terms, "intercept"), 1L)
  expect_identical(attr(fitted_terms, "term.labels"), expected_terms)

  coefficient_names <- rownames(mod_two_basis$coefficients)

  expect_true(all(cb_names %in% coefficient_names))
  expect_true(all(ob_names %in% coefficient_names))

  for (basis_name in c(cb_names, ob_names)) {
    expect_true(any(grepl(
      paste0("^", basis_name, "_id_group:"),
      coefficient_names
    )))
  }
})
