test_that("attributable returns expected structure and summaries (backwards, tot = TRUE)", {

  expect_silent(ar <- attributable(mod, cb, slondon, "date", "tmean", "mort_75plus", cen = cen))

  # basic structure
  expect_type(ar, "list")
  expect_equal(names(ar), c("af", "an", "af.summary", "an.summary"))

  # sample count and names
  expect_equal(length(ar$an), n_sim)
  expect_equal(length(ar$af), n_sim)
  expect_equal(names(ar$an), paste0("sample", seq_len(n_sim)))
  expect_equal(names(ar$af), paste0("sample", seq_len(n_sim)))

  # summaries present and contain mean, sd and at least one quantile + mode
  expect_true(all(c("mean", "sd", "mode") %in% colnames(ar$an.summary)))
  expect_true(any(grepl("quant$", colnames(ar$an.summary))))
  expect_true(all(c("mean", "sd", "mode") %in% colnames(ar$af.summary)))
  expect_true(any(grepl("quant$", colnames(ar$af.summary))))
})

test_that("attributable returns expected structure and summaries (forward, tot = TRUE)", {

  expect_silent(ar <- attributable(mod, cb, slondon, "date", "tmean", "mort_75plus", cen = cen, dir = "forw"))

  # basic structure
  expect_type(ar, "list")
  expect_equal(names(ar), c("af", "an", "af.summary", "an.summary"))

  # sample count and names
  expect_equal(length(ar$an), n_sim)
  expect_equal(length(ar$af), n_sim)
  expect_equal(names(ar$an), paste0("sample", seq_len(n_sim)))
  expect_equal(names(ar$af), paste0("sample", seq_len(n_sim)))

  # summaries present and contain mean, sd and at least one quantile + mode
  expect_true(all(c("mean", "sd", "mode") %in% colnames(ar$an.summary)))
  expect_true(any(grepl("quant$", colnames(ar$an.summary))))
  expect_true(all(c("mean", "sd", "mode") %in% colnames(ar$af.summary)))
  expect_true(any(grepl("quant$", colnames(ar$af.summary))))
})

test_that("attributable returns time-series matrices when tot = FALSE (backward)", {

  expect_silent(ar2 <- attributable(mod, cb, slondon, "date", "tmean", "mort_75plus", cen = cen, tot = FALSE))

  # should return matrices (or vectors if trimmed) for af/an
  n_row <- nrow(slondon) - dlnm_var$max_lag
  expect_equal(dim(ar2$af), c(n_row , n_sim))
  expect_equal(dim(ar2$an), c(n_row, n_sim))

  expect_equal(colnames(ar2$an), paste0("sample", seq_len(n_sim)))
  expect_equal(colnames(ar2$af), paste0("sample", seq_len(n_sim)))

  row_dates <- as.character(slondon$date[(dlnm_var$max_lag + 1):nrow(slondon)])
  expect_equal(rownames(ar2$an), row_dates)
  expect_equal(rownames(ar2$af), row_dates)
  expect_equal(rownames(ar2$an.summary), row_dates)
  expect_equal(rownames(ar2$af.summary), row_dates)

})

test_that("attributable returns time-series matrices when tot = FALSE (forward)", {

  expect_silent(ar2 <- attributable(mod, cb, slondon, "date", "tmean", "mort_75plus", cen = cen, tot = FALSE, dir = "forw"))

  # should return matrices (or vectors if trimmed) for af/an
  n_row <- nrow(slondon) - dlnm_var$max_lag
  expect_equal(dim(ar2$af), c(n_row , n_sim))
  expect_equal(dim(ar2$an), c(n_row, n_sim))

  expect_equal(colnames(ar2$an), paste0("sample", seq_len(n_sim)))
  expect_equal(colnames(ar2$af), paste0("sample", seq_len(n_sim)))

  row_dates <- as.character(slondon$date[1:(nrow(slondon) - dlnm_var$max_lag)])
  expect_equal(rownames(ar2$an), row_dates)
  expect_equal(rownames(ar2$af), row_dates)
  expect_equal(rownames(ar2$an.summary), row_dates)
  expect_equal(rownames(ar2$af.summary), row_dates)

})


test_that("attributable works when cases = NULL (only AF returned)", {

  # Expect a warning (informing that only AF will be calculated)
  expect_warning(
    ar3 <- attributable(mod, cb, slondon, "date", "tmean", cen = cen, tot = TRUE)
  )

  n_row <- nrow(slondon) - dlnm_var$max_lag

  expect_type(ar3, "list")
  expect_equal(names(ar3), c("af", "af.summary"))
  expect_equal(dim(ar3$af), c(n_row, n_sim))
  expect_equal(dim(ar3$af.summary), c(n_row, 6))

})


test_that("attributable gives error when time series are not ordered or not provided on a regular basis", {

  slondon2 <- slondon[order(slondon$tmean),]
  expect_snapshot_error(attributable(mod, cb, slondon2, "date", "tmean", "mort_75plus", cen = cen))

  slondon2 <- slondon[-2,]
  expect_snapshot_error(attributable(mod, cb, slondon2, "date", "tmean", "mort_75plus", cen = cen))

  # Works if date is provided on a weekly basis
  slondon2 <- slondon
  slondon2$date <- seq(as.Date("1900-01-01"), as.Date("2010-01-01"), by = "week")[seq_len(nrow(slondon2))]
  expect_silent(attributable(mod, cb, slondon2, "date", "tmean", "mort_75plus", cen = cen))

  # Works if date is provided on a monthly basis
  slondon2 <- slondon
  slondon2$date <- seq(as.Date("1900-01-01"), as.Date("2010-01-01"), by = "month")[seq_len(nrow(slondon2))]
  expect_silent(attributable(mod, cb, slondon2, "date", "tmean", "mort_75plus", cen = cen))

  # Works if date is provided on a yearly basis
  slondon2 <- slondon
  slondon2$date <- seq(as.Date("1500-01-01"), as.Date("2010-01-01"), by = "year")[seq_len(nrow(slondon2))]
  expect_silent(attributable(mod, cb, slondon2, "date", "tmean", "mort_75plus", cen = cen))

})
