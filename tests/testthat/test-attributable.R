test_that("attributable returns expected structure and summaries (backwards)", {
  skip_on_cran()
  skip_if_not(check_inla(), "INLA not available")

  expect_silent(
    attr <- attributable(
      mod,
      slondon,
      "date",
      "tmean",
      "mort_75plus",
      cen = cen
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
      mod,
      slondon,
      "date",
      "tmean",
      "mort_75plus",
      cen = cen,
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
      mod,
      slondon,
      "date",
      "tmean",
      "mort_75plus",
      "summer",
      cen = cen
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
    attr3 <- attributable(mod, slondon, "date", "tmean", cen = cen)
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
    mod,
    slondon2,
    "date",
    "tmean",
    "mort_75plus",
    cen = cen
  ))

  slondon2 <- slondon[-2, ]
  expect_snapshot_error(attributable(
    mod,
    slondon2,
    "date",
    "tmean",
    "mort_75plus",
    cen = cen
  ))

  # Works if date is provided on a weekly basis
  slondon2 <- slondon
  slondon2$date <- seq(
    as.Date("1900-01-01"),
    as.Date("2010-01-01"),
    by = "week"
  )[seq_len(nrow(slondon2))]
  expect_silent(attributable(
    mod,
    slondon2,
    "date",
    "tmean",
    "mort_75plus",
    cen = cen
  ))

  # Works if date is provided on a monthly basis
  slondon2 <- slondon
  slondon2$date <- seq(
    as.Date("1900-01-01"),
    as.Date("2010-01-01"),
    by = "month"
  )[seq_len(nrow(slondon2))]
  expect_silent(attributable(
    mod,
    slondon2,
    "date",
    "tmean",
    "mort_75plus",
    cen = cen
  ))

  # Works if date is provided on a yearly basis
  slondon2 <- slondon
  slondon2$date <- seq(
    as.Date("1500-01-01"),
    as.Date("2010-01-01"),
    by = "year"
  )[seq_len(nrow(slondon2))]
  expect_silent(attributable(
    mod,
    slondon2,
    "date",
    "tmean",
    "mort_75plus",
    cen = cen
  ))
})
