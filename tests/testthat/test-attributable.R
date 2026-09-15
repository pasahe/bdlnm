test_that("attributable returns expected structure and summaries (backwards)", {
  skip_on_cran()
  skip_if_not(check_inla(), "INLA not available")

  expect_silent(
    attr <- attributable(
      cpred_centered,
      slondon,
      columns = list(date = "date", exp = "tmean", cases = "mort_75plus")
    )
  )

  # basic structure
  expect_type(attr, "list")
  expect_equal(
    names(attr),
    c(
      "af",
      "an",
      "aftotal",
      "antotal",
      "af.summary",
      "an.summary",
      "aftotal.summary",
      "antotal.summary"
    )
  )

  # sample count and names
  expect_equal(dim(attr$an), c(nrow(slondon), n_sim))
  expect_equal(dim(attr$af), c(nrow(slondon), n_sim))
  expect_equal(colnames(attr$an), paste0("sample", seq_len(n_sim)))
  expect_equal(colnames(attr$af), paste0("sample", seq_len(n_sim)))

  # summaries present and contain mean, sd and at least one quantile + mode
  expect_true(all(c("mean", "sd", "mode") %in% colnames(attr$an.summary)))
  expect_true(any(grepl("quant$", colnames(attr$an.summary))))
  expect_true(all(c("mean", "sd", "mode") %in% colnames(attr$af.summary)))
  expect_true(any(grepl("quant$", colnames(attr$af.summary))))

  #NA at the beginning
  expect_equal(
    names(which(is.na(attr$af.summary[, "0.5quant"]))),
    as.character(slondon$date[1:21])
  )
  expect_equal(
    names(which(is.na(attr$an.summary[, "0.5quant"]))),
    as.character(slondon$date[1:21])
  )
})

test_that("attributable returns expected structure and summaries (forward)", {
  skip_on_cran()
  skip_if_not(check_inla(), "INLA not available")

  expect_silent(
    attr <- attributable(
      cpred_centered,
      slondon,
      columns = list(date = "date", exp = "tmean", cases = "mort_75plus"),
      dir = "forw"
    )
  )

  # basic structure
  expect_type(attr, "list")
  expect_equal(
    names(attr),
    c(
      "af",
      "an",
      "aftotal",
      "antotal",
      "af.summary",
      "an.summary",
      "aftotal.summary",
      "antotal.summary"
    )
  )

  # sample count and names
  expect_equal(dim(attr$an), c(nrow(slondon), n_sim))
  expect_equal(dim(attr$af), c(nrow(slondon), n_sim))
  expect_equal(colnames(attr$an), paste0("sample", seq_len(n_sim)))
  expect_equal(colnames(attr$af), paste0("sample", seq_len(n_sim)))

  # summaries present and contain mean, sd and at least one quantile + mode
  expect_true(all(c("mean", "sd", "mode") %in% colnames(attr$an.summary)))
  expect_true(any(grepl("quant$", colnames(attr$an.summary))))
  expect_true(all(c("mean", "sd", "mode") %in% colnames(attr$af.summary)))
  expect_true(any(grepl("quant$", colnames(attr$af.summary))))

  #NA at the end
  expect_equal(sum(is.na(attr$af.summary[, "0.5quant"])), 0)
  expect_equal(
    names(which(is.na(attr$an.summary[, "0.5quant"]))),
    as.character(sort(rev(slondon$date)[1:21]))
  )
})


test_that("attributable returns filtered time-series matrices when filter is specified", {
  skip_on_cran()
  skip_if_not(check_inla(), "INLA not available")

  # filter only for summer
  summer_dates <- slondon$date[
    slondon$date >= as.Date("2011-06-01") &
      slondon$date <= as.Date("2011-09-30")
  ]
  slondon$summer <- ifelse(slondon$date %in% summer_dates, 1, 0)

  expect_warning(
    attr2 <- attributable(
      cpred_centered,
      slondon,
      columns = list(
        date = "date",
        exp = "tmean",
        cases = "mort_75plus",
        filter = "summer"
      )
    )
  )

  expect_equal(dim(attr2$af), c(length(summer_dates), n_sim))
  expect_equal(dim(attr2$an), c(length(summer_dates), n_sim))

  expect_equal(colnames(attr2$an), paste0("sample", seq_len(n_sim)))
  expect_equal(colnames(attr2$af), paste0("sample", seq_len(n_sim)))

  expect_equal(rownames(attr2$an), as.character(summer_dates))
  expect_equal(rownames(attr2$af), as.character(summer_dates))
  expect_equal(rownames(attr2$an.summary), as.character(summer_dates))
  expect_equal(rownames(attr2$af.summary), as.character(summer_dates))

  # no missings
  expect_equal(sum(is.na(attr2$af.summary[, "0.5quant"])), 0)
  expect_equal(sum(is.na(attr2$an.summary[, "0.5quant"])), 0)
})

test_that("attributable works when cases = NULL (only AF returned)", {
  skip_on_cran()
  skip_if_not(check_inla(), "INLA not available")

  # Expect a warning (informing that only AF will be calculated)
  expect_warning(
    attr3 <- attributable(
      cpred_centered,
      slondon,
      columns = list(date = "date", exp = "tmean")
    )
  )

  expect_type(attr3, "list")
  expect_equal(names(attr3), c("af", "af.summary"))
  expect_equal(dim(attr3$af), c(nrow(slondon), n_sim))
  expect_equal(dim(attr3$af.summary), c(nrow(slondon), 6))
})


test_that("attributable gives error when time series are not ordered or not provided on a regular basis", {
  skip_on_cran()
  skip_if_not(check_inla(), "INLA not available")

  slondon2 <- slondon[order(slondon$tmean), ]
  expect_snapshot_error(attributable(
    cpred_centered,
    slondon2,
    columns = list(date = "date", exp = "tmean", cases = "mort_75plus")
  ))

  slondon2 <- slondon[-2, ]
  expect_snapshot_error(attributable(
    cpred_centered,
    slondon2,
    columns = list(date = "date", exp = "tmean", cases = "mort_75plus")
  ))

  # Works if date is provided on a weekly basis
  slondon2 <- slondon
  slondon2$date <- seq(
    as.Date("1900-01-01"),
    as.Date("2010-01-01"),
    by = "week"
  )[seq_len(nrow(slondon2))]
  expect_silent(attributable(
    cpred_centered,
    slondon2,
    columns = list(date = "date", exp = "tmean", cases = "mort_75plus")
  ))

  # Works if date is provided on a monthly basis
  slondon2 <- slondon
  slondon2$date <- seq(
    as.Date("1900-01-01"),
    as.Date("2010-01-01"),
    by = "month"
  )[seq_len(nrow(slondon2))]
  expect_silent(attributable(
    cpred_centered,
    slondon2,
    columns = list(date = "date", exp = "tmean", cases = "mort_75plus")
  ))

  # Works if date is provided on a yearly basis
  slondon2 <- slondon
  slondon2$date <- seq(
    as.Date("1500-01-01"),
    as.Date("2010-01-01"),
    by = "year"
  )[seq_len(nrow(slondon2))]
  expect_silent(attributable(
    cpred_centered,
    slondon2,
    columns = list(date = "date", exp = "tmean", cases = "mort_75plus")
  ))
})

test_that("attributable requires predictions when int is supplied", {
  skip_on_cran()
  skip_if_not(check_inla(), "INLA not available")

  expect_error(
    attributable(
      mod_interaction,
      interaction_data,
      columns = list(
        date = "date",
        exp = "tmean",
        cases = "mort_75plus",
        int = "group"
      ),
      cen = cen
    ),
    "containing basis interactions cannot be used directly"
  )
})

test_that("attributable must be calculated from centered predictions", {
  skip_on_cran()
  skip_if_not(check_inla(), "INLA not available")

  cpred_without_center <- cpred_centered
  cpred_without_center$cen <- NULL

  expect_error(
    attributable(
      cpred_without_center,
      slondon,
      "date",
      "tmean",
      "mort_75plus"
    ),
    "centered at a reference exposure"
  )

})

test_that("attributable uses fixed interaction predictions", {
  skip_on_cran()
  skip_if_not(check_inla(), "INLA not available")

  fixed_attributable <- attributable(
    cpred_fixed_interaction,
    interaction_data,
    columns = list(
      date = "date",
      exp = "tmean",
      cases = "mort_75plus",
      int = "group"
    )
  )

  expect_equal(dim(fixed_attributable$af), c(nrow(interaction_data), 2L))
  expect_equal(dim(fixed_attributable$an), c(nrow(interaction_data), 2L))
})

test_that("attributable uses random interaction predictions", {
  skip_on_cran()
  skip_if_not(check_inla(), "INLA not available")

  random_attributable <- attributable(
    cpred_random_interaction,
    random_data,
    columns = list(
      date = "date",
      exp = "tmean",
      cases = "mort_75plus",
      int = "id_group"
    ),
    dir = "forw"
  )

  expect_equal(dim(random_attributable$af), c(nrow(random_data), 2L))
  expect_equal(dim(random_attributable$an), c(nrow(random_data), 2L))
})

test_that("attributable validates interaction columns", {
  skip_on_cran()
  skip_if_not(check_inla(), "INLA not available")

  expect_error(
    attributable(
      cpred_fixed_interaction,
      interaction_data,
      columns = list(
        date = "date",
        exp = "tmean",
        cases = "mort_75plus"
      )
    ),
    "must contain an.*int.*element"
  )

  expect_error(
    attributable(
      cpred_fixed_interaction,
      interaction_data,
      columns = list(
        date = "date",
        exp = "tmean",
        cases = "mort_75plus",
        int = "id_group"
      )
    ),
    "does not match the interaction"
  )

  interaction_data_missing_level <- interaction_data
  interaction_data_missing_level$group <- as.character(
    interaction_data_missing_level$group
  )
  interaction_data_missing_level$group[1L] <- "c"
  expect_error(
    attributable(
      cpred_fixed_interaction,
      interaction_data_missing_level,
      columns = list(
        date = "date",
        exp = "tmean",
        cases = "mort_75plus",
        int = "group"
      )
    ),
    "does not contain predictions for every interaction value"
  )
})
