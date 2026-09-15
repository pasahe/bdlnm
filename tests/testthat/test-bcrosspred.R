test_that("bcrosspred basic structure (crossbasis)", {
  skip_on_cran()
  skip_if_not(check_inla(), "INLA not available")

  #cpred build in helper.R

  expect_s3_class(cpred, "bcrosspred")
  expect_equal(length(cpred), 16)
  expect_equal(
    names(cpred),
    c(
      "exp_at",
      "lag_at",
      "cen",
      "coefficients",
      "matfit",
      "allfit",
      "matRRfit",
      "allRRfit",
      "coefficients.summary",
      "matfit.summary",
      "allfit.summary",
      "matRRfit.summary",
      "allRRfit.summary",
      "ci.level",
      "model.class",
      "model.link"
    )
  )

  expect_equal(cpred$exp_at, temp)
  expect_equal(cpred$cen, 12.85)
  expect_equal(
    cpred$lag_at,
    seq(attr(cb, "lag")[1], attr(cb, "lag")[2], by = 1)
  )
  expect_equal(dim(cpred$coefficients), c(20, n_sim))
  expect_equal(rownames(cpred$coefficients), paste0("cb", colnames(cb)))
  expect_equal(dim(cpred$matfit), c(length(temp), 22, n_sim))
  expect_equal(rownames(cpred$matfit), as.character(temp))
  expect_equal(dim(cpred$allfit), c(length(temp), n_sim))
  expect_equal(rownames(cpred$allfit), as.character(temp))
  expect_true(all(cpred$allRRfit > 0))
  expect_equal(dim(cpred$matfit.summary), c(length(temp), 22, 6))
  expect_equal(cpred$ci.level, 0.95)
  expect_equal(cpred$model.class, "inla")
  expect_equal(cpred$model.link, "log")
})

test_that("bcrosspred basic structure (onebasis)", {
  skip_on_cran()
  skip_if_not(check_inla(), "INLA not available")

  ob <- dlnm::onebasis(slondon$tmean, "strata", breaks = c(5, 10, 20))
  mod_2 <- suppressWarnings(bdlnm(
    mort_75plus ~ ob + factor(dow) + seas,
    data = slondon,
    family = "poisson",
    sample.arg = list(n = n_sim, seed = 1L)
  ))
  cpred_2 <- bcrosspred(mod_2, exp_at = temp)

  expect_s3_class(cpred_2, "bcrosspred")
  expect_equal(length(cpred_2), 14)
  expect_equal(
    names(cpred_2),
    c(
      "exp_at",
      "coefficients",
      "matfit",
      "allfit",
      "matRRfit",
      "allRRfit",
      "coefficients.summary",
      "matfit.summary",
      "allfit.summary",
      "matRRfit.summary",
      "allRRfit.summary",
      "ci.level",
      "model.class",
      "model.link"
    )
  )

  expect_equal(cpred_2$exp_at, temp)
  expect_equal(cpred_2$lag_at, NULL)
  expect_equal(dim(cpred_2$coefficients), c(3, n_sim))
  expect_equal(rownames(cpred_2$coefficients), paste0("ob", colnames(ob)))
  expect_equal(dim(cpred_2$matfit), c(length(temp), 1, n_sim))
  expect_equal(rownames(cpred_2$matfit), as.character(temp))
  expect_equal(dim(cpred_2$allfit), c(length(temp), n_sim))
  expect_equal(rownames(cpred_2$allfit), as.character(temp))
  expect_true(all(cpred_2$allRRfit > 0))
  expect_equal(dim(cpred_2$matfit.summary), c(length(temp), 1, 6))
  expect_equal(cpred_2$ci.level, 0.95)
  expect_equal(cpred_2$model.class, "inla")
  expect_equal(cpred_2$model.link, "log")
})

test_that("works with two different basis", {
  skip_on_cran()
  skip_if_not(check_inla(), "INLA not available")

  ob <- dlnm::onebasis(slondon$tmean, "strata", breaks = c(5, 10, 20))
  mod_2 <- suppressWarnings(bdlnm(
    mort_75plus ~ cb + ob + factor(dow) + seas,
    data = slondon,
    family = "poisson",
    sample.arg = list(n = n_sim, seed = 1L)
  ))

  expect_error(bcrosspred(mod_2, exp_at = temp))

  cpred_2 <- bcrosspred(mod_2, "ob", exp_at = temp)

  expect_s3_class(cpred_2, "bcrosspred")
  expect_equal(length(cpred_2), 14)
  expect_equal(
    names(cpred_2),
    c(
      "exp_at",
      "coefficients",
      "matfit",
      "allfit",
      "matRRfit",
      "allRRfit",
      "coefficients.summary",
      "matfit.summary",
      "allfit.summary",
      "matRRfit.summary",
      "allRRfit.summary",
      "ci.level",
      "model.class",
      "model.link"
    )
  )

  expect_equal(cpred_2$exp_at, temp)
  expect_equal(cpred_2$lag_at, NULL)
  expect_equal(dim(cpred_2$coefficients), c(3, n_sim))
  expect_equal(rownames(cpred_2$coefficients), paste0("ob", colnames(ob)))
  expect_equal(dim(cpred_2$matfit), c(length(temp), 1, n_sim))
  expect_equal(rownames(cpred_2$matfit), as.character(temp))
  expect_equal(dim(cpred_2$allfit), c(length(temp), n_sim))
  expect_equal(rownames(cpred_2$allfit), as.character(temp))
  expect_true(all(cpred_2$allRRfit > 0))
  expect_equal(dim(cpred_2$matfit.summary), c(length(temp), 1, 6))
  expect_equal(cpred_2$ci.level, 0.95)
  expect_equal(cpred_2$model.class, "inla")
  expect_equal(cpred_2$model.link, "log")
})

test_that("bcrosspred errors when some argument is missing or invalid", {
  skip_on_cran()
  skip_if_not(check_inla(), "INLA not available")

  expect_snapshot_error(bcrosspred(NULL, exp_at = temp))

  expect_snapshot_error(bcrosspred(list(a = 1), "cb", exp_at = temp))

  expect_snapshot_error(bcrosspred(mod, basis = "cb1", exp_at = temp))

  expect_warning(cpred2 <- bcrosspred(mod))

  expect_equal(length(cpred2$exp_at), 48)
})

test_that("bcrosspred works when a random effect is included", {
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
  expect_warning(cpred_rt <- bcrosspred(mod_rt, exp_at = temp))
  expect_s3_class(cpred_rt, "bcrosspred")
  expect_equal(length(cpred_rt), 16)
  expect_equal(
    names(cpred_rt),
    c(
      "exp_at",
      "lag_at",
      "cen",
      "coefficients",
      "matfit",
      "allfit",
      "matRRfit",
      "allRRfit",
      "coefficients.summary",
      "matfit.summary",
      "allfit.summary",
      "matRRfit.summary",
      "allRRfit.summary",
      "ci.level",
      "model.class",
      "model.link"
    )
  )
  expect_equal(cpred_rt$exp_at, temp)
  expect_equal(cpred_rt$cen, 12.85)
  expect_equal(
    cpred_rt$lag_at,
    seq(attr(cb, "lag")[1], attr(cb, "lag")[2], by = 1)
  )
  expect_equal(dim(cpred_rt$coefficients), c(20, n_sim))
  expect_equal(rownames(cpred_rt$coefficients), paste0("cb", colnames(cb)))
  expect_equal(dim(cpred_rt$matfit), c(length(temp), 22, n_sim))
  expect_equal(rownames(cpred_rt$matfit), as.character(temp))
  expect_equal(dim(cpred_rt$allfit), c(length(temp), n_sim))
  expect_equal(rownames(cpred_rt$allfit), as.character(temp))
  expect_true(all(cpred_rt$allRRfit > 0))
  expect_equal(dim(cpred_rt$matfit.summary), c(length(temp), 22, 6))
  expect_equal(cpred_rt$ci.level, 0.95)
  expect_equal(cpred_rt$model.class, "inla")
  expect_equal(cpred_rt$model.link, "log")

})

test_that("bcrosspred combines fixed basis interactions", {
  skip_on_cran()
  skip_if_not(check_inla(), "INLA not available")

  cpred_basis_second <- suppressWarnings(
    bcrosspred(
      mod_interaction,
      combine_effects = ~ group * cb,
      exp_at = temp
    )
  )
  cpred_basis_extended <- suppressWarnings(
    bcrosspred(
      mod_interaction,
      combine_effects = ~ cb + cb:group,
      exp_at = temp
    )
  )
  cpred_interaction_only <- suppressWarnings(
    bcrosspred(
      mod_interaction,
      combine_effects = ~ cb:group,
      exp_at = temp
    )
  )

  expect_identical(
    names(cpred_fixed_interaction)[seq_len(2L)],
    c("group[a]", "group[b]")
  )
  expect_identical(
    attr(cpred_fixed_interaction, "basis_interaction")$variable,
    "group"
  )
  expect_identical(
    cpred_fixed_interaction[["group[a]"]]$coefficients,
    mod_interaction$coefficients[paste0("cb", colnames(cb)), , drop = FALSE]
  )
  expect_identical(
    cpred_fixed_interaction[["group[a]"]]$coefficients,
    cpred_basis_second[["group[a]"]]$coefficients,
    cpred_basis_extended[["group[a]"]]$coefficients
  )
  expect_identical(
    cpred_fixed_interaction[["group[b]"]]$coefficients,
    cpred_basis_second[["group[b]"]]$coefficients,
    cpred_basis_extended[["group[b]"]]$coefficients
  )

  expect_warning(
    cpred_global <- bcrosspred(
      mod_interaction,
      combine_effects = ~ cb,
      exp_at = temp,
      cen = 12.85
    ),
    "contains interaction effects"
  )
  expect_identical(
    cpred_global$coefficients,
    mod_interaction$coefficients[paste0("cb", colnames(cb)), , drop = FALSE]
  )

  expect_identical(
    names(cpred_interaction_only)[[1L]],
    "group[b]"
  )

  interaction_rows <- paste0("cb", colnames(cb), ":groupb")
  expected_interaction_coef <- mod_interaction$coefficients[
    interaction_rows,
    ,
    drop = FALSE
  ]
  rownames(expected_interaction_coef) <- paste0("cb", colnames(cb))

  expect_identical(
    cpred_interaction_only[["group[b]"]]$coefficients,
    expected_interaction_coef
  )

})

test_that("bcrosspred combines a fixed interaction without a standalone basis", {
  skip_on_cran()
  skip_if_not(check_inla(), "INLA not available")

  interaction_data <- slondon
  interaction_data$group <- factor(
    rep(c("a", "b"), length.out = nrow(interaction_data))
  )
  interaction_data$id_group <- rep(
    seq_len(5), length.out = nrow(interaction_data)
  )

  mod_pure_interaction <- suppressWarnings(
    bdlnm(
      mort_75plus ~ cb:group + factor(dow) + seas +
        f(id_group, model = "iid"),
      data = interaction_data,
      family = "poisson",
      sample.arg = list(n = 2L, seed = 1L)
    )
  )

  cpred_pure_interaction <- suppressWarnings(
    bcrosspred(
      mod_pure_interaction,
      combine_effects = ~ cb:group,
      exp_at = temp
    )
  )

  expect_identical(
    names(cpred_pure_interaction)[seq_len(2L)],
    c("group[a]", "group[b]")
  )

  for (group_level in c("a", "b")) {
    interaction_rows <- paste0("cb", colnames(cb), ":group", group_level)
    expected_coef <- mod_pure_interaction$coefficients[
      interaction_rows,
      ,
      drop = FALSE
    ]
    rownames(expected_coef) <- paste0("cb", colnames(cb))

    expect_identical(
      cpred_pure_interaction[[paste0("group[", group_level, "]")]]$coefficients,
      expected_coef
    )
  }

  expect_error(
    suppressWarnings(
      bcrosspred(
        mod_pure_interaction,
        combine_effects = ~ cb,
        exp_at = temp
      )
    ),
    "Global basis coefficients were requested"
  )
})

test_that("bcrosspred combines basis-specific random effects", {
  skip_on_cran()
  skip_if_not(check_inla(), "INLA not available")

  cpred_random_only <- suppressWarnings(
    bcrosspred(
      mod_random,
      combine_effects = ~ f(id_group, cb, model = "iid"),
      exp_at = temp
    )
  )

  expected_names <- c("global", paste0("id_group[", seq_len(5), "]"))
  expect_identical(
    names(cpred_random_interaction)[seq_along(expected_names)],
    expected_names
  )
  expect_identical(
    attr(cpred_random_interaction, "basis_interaction")$variable,
    "id_group"
  )
  expect_identical(
    names(cpred_random_only)[seq_len(5L)],
    expected_names[-1L]
  )

  global_coef <- mod_random$coefficients[paste0("cb", colnames(cb)), , drop = FALSE]
  random_rows <- paste0("cb", colnames(cb), "_id_group:1")
  random_coef <- mod_random$coefficients[random_rows, , drop = FALSE]
  rownames(random_coef) <- rownames(global_coef)

  expect_equal(
    cpred_random_interaction[["id_group[1]"]]$coefficients,
    global_coef + random_coef
  )
  expect_identical(
    cpred_random_only[["id_group[1]"]]$coefficients,
    random_coef
  )

})

test_that("bcrosspred combines a random interaction without a standalone basis", {
  skip_on_cran()
  skip_if_not(check_inla(), "INLA not available")

  random_data <- slondon
  random_data$id_group <- rep(seq_len(5), length.out = nrow(random_data))

  mod_without_standalone <- suppressWarnings(
    bdlnm(
      mort_75plus ~ factor(dow) + seas +
        f(id_group, cb, model = "iid"),
      data = random_data,
      family = "poisson",
      sample.arg = list(n = 2L, seed = 1L)
    )
  )
  cpred_without_standalone <- suppressWarnings(
    bcrosspred(
      mod_without_standalone,
      combine_effects = ~ f(id_group, cb, model = "iid"),
      exp_at = temp
    )
  )

  expect_identical(
    names(cpred_without_standalone)[seq_len(5L)],
    paste0("id_group[", seq_len(5), "]")
  )

  random_only_rows <- paste0("cb", colnames(cb), "_id_group:1")
  expected_random_only <- mod_without_standalone$coefficients[
    random_only_rows,
    ,
    drop = FALSE
  ]
  rownames(expected_random_only) <- paste0("cb", colnames(cb))

  expect_identical(
    cpred_without_standalone[["id_group[1]"]]$coefficients,
    expected_random_only
  )

  expect_error(
    suppressWarnings(
      bcrosspred(
        mod_without_standalone,
        combine_effects = ~ cb,
        exp_at = temp
      )
    ),
    "Global basis coefficients were requested"
  )

})

test_that("bcrosspred rejects multiple basis-specific random effects", {
  skip_on_cran()
  skip_if_not(check_inla(), "INLA not available")

  multiple_random_data <- slondon
  multiple_random_data$id_group <- rep(
    seq_len(5), length.out = nrow(multiple_random_data)
  )
  multiple_random_data$id_other <- rep(
    seq_len(4), length.out = nrow(multiple_random_data)
  )

  mod_multiple_random <- suppressWarnings(
    bdlnm(
      mort_75plus ~ cb + factor(dow) + seas +
        f(id_group, cb, model = "iid") +
        f(id_other, cb, model = "iid"),
      data = multiple_random_data,
      family = "poisson",
      sample.arg = list(n = 2L, seed = 1L)
    )
  )

  expect_error(
    bcrosspred(
      mod_multiple_random,
      combine_effects = ~ cb +
        f(id_group, cb, model = "iid") +
        f(id_other, cb, model = "iid"),
      exp_at = temp
    ),
    "Only one basis-specific"
  )
})

test_that("bcrosspred rejects simultaneous fixed and random basis interactions", {
  skip_on_cran()
  skip_if_not(check_inla(), "INLA not available")

  fixed_random_data <- slondon
  fixed_random_data$group <- factor(
    rep(c("a", "b"), length.out = nrow(fixed_random_data))
  )
  fixed_random_data$id_group <- rep(
    seq_len(5), length.out = nrow(fixed_random_data)
  )

  mod_fixed_random <- suppressWarnings(
    bdlnm(
      mort_75plus ~ cb * group + factor(dow) + seas +
        f(id_group, cb, model = "iid"),
      data = fixed_random_data,
      family = "poisson",
      sample.arg = list(n = 2L, seed = 1L)
    )
  )

  expect_error(
    bcrosspred(
      mod_fixed_random,
      combine_effects = ~ cb * group +
        f(id_group, cb, model = "iid"),
      exp_at = temp
    ),
    "Fixed and random basis interactions"
  )
})
