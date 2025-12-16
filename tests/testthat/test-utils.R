test_that("df2ldl produces expected numerical outout", {
  data(england, package = 'tsgc')
  lng <- df2ldl(england$cum_cases[1:10])
  x <- as.numeric(lng[2])
  expect_equal(x,log((7915-6877)/6877))
})

test_that("df2ldl returns xts object", {
  data(england, package = 'tsgc')
  lng <- df2ldl(england$cum_cases[1:10])
  expect_s3_class(lng,"xts")
})

test_that("get_timeframe works correctly", {
  data(england, package = 'tsgc')
  x <- get_timeframe(england, '2020-03-20', '2020-03-29')
  y <- head(england,10)
  expect_equal(x,y)
})

test_that("reinitialise_dataframe works correctly", {
  data(england, package = 'tsgc')
  reinit <- reinitialise_dataframe(england$cum_cases,as.Date('2022-03-14'))
  x <- as.numeric(reinit[1,1])
  expect_equal(x,79242)
})

test_that("argmax returns expected output", {
  x <- xts::xts(1:10, order.by = seq(as.Date("2021-01-01"),
                                     length.out = 10, by = 1))
  x[5,] <- 20
  expect_identical(str(zoo::index(argmax(x))), str(zoo::index(x)[5]))
  expect_identical(
    str(zoo::index(argmax(x, decreasing = FALSE))),
    str(zoo::index(x)[1])
  )
})

test_that("estimate_r0 returns correct data frame structure", {
  data(gauteng, package = "tsgc")
  cumulative_cases <- gauteng[, 1]
  gen_int <- 4
  ndays <- 7
  
  # Estimate model
  model_q <- SSModelDynamicGompertz$new(
    Y = cumulative_cases,
    q = 0.005,
    start.date = as.Date("2021-02-01"),
    end.date = as.Date("2021-04-19")
  )
  res_q <- estimate(model_q)
  # Call the function for data frame output
  r_t_df <- estimate_r0(res_q, gen_int, ndays, show_plot = FALSE)
  
  # Check the class and dimensions
  expect_s3_class(r_t_df, "data.frame")
  expect_equal(nrow(r_t_df), ndays)
  expect_equal(ncol(r_t_df), 4)
  
  # Check column names
  expected_names <- c("Date", "Rt", "lower", "upper")
  expect_equal(names(r_t_df), expected_names)
  
  # Check column types
  expect_s3_class(r_t_df$Date, "Date")
  expect_type(r_t_df$Rt, "double")
  expect_type(r_t_df$lower, "double")
  expect_type(r_t_df$upper, "double")
  
  # Check logical constraint: lower <= Rt <= upper
  expect_true(all(r_t_df$lower <= r_t_df$Rt))
  expect_true(all(r_t_df$Rt <= r_t_df$upper))
})

test_that("estimate_r0 returns a ggplot object when show_plot is TRUE", {
  data(gauteng, package = "tsgc")
  cumulative_cases <- gauteng[, 1]
  gen_int <- 4
  ndays <- 7
  
  # Estimate model
  model_q <- SSModelDynamicGompertz$new(
    Y = cumulative_cases,
    q = 0.005,
    start.date = as.Date("2021-02-01"),
    end.date = as.Date("2021-04-19")
  )
  res_q <- estimate(model_q)
  
  # Call the function for plot output
  r_t_plot <- estimate_r0(res_q, gen_int, ndays, show_plot = TRUE)
  
  # Check the class of the returned object
  expect_s3_class(r_t_plot, "ggplot")
})

test_that("get_time_resolution handles invalid inputs correctly", {
  # Not a vector of dates (length > 1)
  expect_error(get_time_resolution(1:5), "Input is not a vector of dates.")
  # Single date
  expect_error(get_time_resolution(Sys.Date()), "Input is not a vector of dates.")
  # Dates with unequal time separation
  dates_unequal <- as.Date(c("2024-01-01", "2024-01-03", "2024-01-07"))
  expect_error(get_time_resolution(dates_unequal), "The dates are not separated by the same time resolution.")
})

test_that("get_time_resolution identifies daily resolution", {
  daily_dates <- seq(as.Date("2024-01-01"), by = "day", length.out = 5)
  expect_equal(get_time_resolution(daily_dates), "daily")
})

test_that("get_time_resolution identifies monthly resolution", {
  monthly_dates <- seq(as.yearmon("2024-01"), by = 1/12, length.out = 4)
  expect_equal(get_time_resolution(monthly_dates), "monthly")
})

test_that("get_time_resolution identifies quarterly resolution", {
  quarterly_dates <- seq(as.yearqtr("2024-01"), by = 0.25, length.out = 4)
  expect_equal(get_time_resolution(quarterly_dates), "quarterly")
})

test_that("get_time_resolution identifies yearly resolution (yearmon)", {
  yearly_mon_dates <- seq(as.yearmon("2024-01"), by = 1, length.out = 3)
  expect_equal(get_time_resolution(yearly_mon_dates), "yearly")
})

test_that("get_time_resolution identifies yearly resolution (yearqtr)", {
  yearly_qtr_dates <- seq(as.yearqtr("2024-01"), by = 1, length.out = 3)
  expect_equal(get_time_resolution(yearly_qtr_dates), "yearly")
})

test_that("qtr2date correctly converts yearqtr to Date", {
  # Input: Single yearqtr object
  qtr_date <- zoo::yearqtr(2024.25) # Q2 2024
  expected_date <- as.Date("2024-04-01") # First day of Q2
  converted_date <- qtr2date(qtr_date)
  
  expect_s3_class(converted_date, "Date")
  expect_equal(converted_date, expected_date)
  
  # Input: Multiple yearqtr objects
  multi_qtr_dates <- c(zoo::yearqtr(2020.0), zoo::yearqtr(2020.5)) # Q1, Q3 2020
  expected_dates <- as.Date(c("2020-01-01", "2020-07-01"))
  expect_equal(qtr2date(multi_qtr_dates), expected_dates)
})

test_that("qtr2date correctly converts yearmon to Date", {
  # Input: Multiple yearmon objects
  mon_dates <- zoo::yearmon(c(2024 + 5/12, 2024 + 8/12)) # June, September 2024
  expected_dates <- as.Date(c("2024-06-01", "2024-09-01"))
  converted_dates <- qtr2date(mon_dates)
  
  expect_s3_class(converted_dates, "Date")
  expect_equal(converted_dates, expected_dates)
})

test_that("seq_dates handles input constraints correctly", {
  from_date <- as.Date("2025-01-01")
  
  # Constraint 1: Both 'length.out' and 'to' cannot be empty
  expect_error(seq_dates(from_date, "daily"), "Both length.out and to inputs cannot be empty.")
  
  # Constraint 2: Only one of 'length.out' or 'to' can be supplied
  to_date <- as.Date("2025-01-05")
  expect_error(seq_dates(from_date, "daily", to = to_date, length.out = 5), "Please supply only one of length.out or to.")
})

test_that("seq_dates generates daily sequences correctly", {
  from_date <- as.Date("2025-01-01")
  
  # Scenario 1: Using length.out
  expected_len <- seq(from_date, by = 'day', length.out = 5)
  expect_equal(seq_dates(from_date, "daily", length.out = 5), expected_len)
  
  # Scenario 2: Using 'to'
  to_date <- as.Date("2025-01-05")
  expected_to <- seq(from_date, to_date, by = 'day')
  expect_equal(seq_dates(from_date, "daily", to = to_date), expected_to)
})

test_that("seq_dates generates quarterly sequences correctly", {
  from_qtr <- zoo::as.yearqtr("2024-01-01") # Q1 2024
  
  # Scenario 1: Using length.out
  expected_len <- zoo::as.yearqtr(seq(as.numeric(from_qtr), by = 0.25, length.out = 3))
  expect_equal(seq_dates(from_qtr, "quarterly", length.out = 3), expected_len)
  
  # Scenario 2: Using 'to'
  to_qtr <- zoo::as.yearqtr(as.Date("2024-07-01")) # Q3 2024
  expected_to <- zoo::as.yearqtr(seq(as.numeric(from_qtr), as.numeric(to_qtr), by = 0.25))
  expect_equal(seq_dates(from_qtr, "quarterly", to = to_qtr), expected_to)
})

test_that("seq_dates generates monthly sequences correctly", {
  from_mon <- zoo::as.yearmon("2025-03")
  
  # Scenario 1: Using length.out
  expected_len <- zoo::as.yearmon(seq(as.numeric(from_mon), by = 1/12, length.out = 4))
  expect_equal(seq_dates(from_mon, "monthly", length.out = 4), expected_len)
  
  # Scenario 2: Using 'to'
  to_mon <- zoo::as.yearmon("2025-06")
  expected_to <- zoo::as.yearmon(seq(as.numeric(from_mon), as.numeric(to_mon), by = 1/12))
  expect_equal(seq_dates(from_mon, "monthly", to = to_mon), expected_to)
})

test_that("seq_dates generates yearly sequences correctly (from yearmon)", {
  from_mon <- zoo::as.yearmon("2024-01")
  
  # Using length.out
  expected_len <- zoo::as.yearmon(seq(as.numeric(from_mon), by = 1, length.out = 3))
  expect_equal(seq_dates(from_mon, "yearly", length.out = 3), expected_len)
})



test_that("cross_val reports the correct criterion (e.g., 'rmse')", {
  data(ukitaly, package = "tsgc")
  Yuk <- tsgc::ukitaly[, "UK"]
  est.start <- as.Date("2020-02-25")
  est.end <- as.Date("2020-04-01")
  
  # Create a list of models for comparison
  cv_models <- list()
  cv_models[["Vanilla_q"]] <- SSModelDynamicGompertz$new(Y = Yuk, q = 0.005, start.date = est.start, end.date = est.end)
  cv_models[["Vanilla_ar1"]] <- SSModelDynamicGompertz$new(Y = Yuk, start.date = est.start, end.date = est.end, ar1 = TRUE)
  cv_models[["Lag7"]] <- SSModelLeadingIndicator$new(Y = ukitaly, start.date = est.start, end.date = est.end, n.lag = 7)
  
  n.ahead <- 7
  n.estimate <- 3
  gap <- 2
  
  # Run validation using RMSE
  cv_result_rmse <- cross_val(
    Y = ukitaly,
    model_list = cv_models,
    est.end.date = est.end,
    n.ahead = n.ahead,
    n.estimate = 1, # Use 1 estimate for a simple check
    gap = gap,
    criterion = "rmse"
  )
  
  # Check if the result is numeric and positive
  expect_type(cv_result_rmse$`2020-04-01`, "double")
  expect_true(all(cv_result_rmse$`2020-04-01` > 0))
  
  # Spot-check one value (e.g., Vanilla_q RMSE for the first date)
  # Calculate the expected value manually or use a previously validated constant
  model_q_test <- SSModelDynamicGompertz$new(Y = Yuk, q = 0.005, start.date = est.start, end.date = est.end)
  res_q_test <- estimate(model_q_test)
  expected_rmse <- round(mapes(res_q_test, n.ahead, Yuk)[["rmse"]], 2)
  
  # Expect the cross_val result to be close to the manual result
  expect_equal(cv_result_rmse[cv_result_rmse$Model == "Vanilla_q", "2020-04-01"][[1]], expected_rmse)
})