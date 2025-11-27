library(testthat)
library(tsgc)
library(KFAS)
library(xts)
library(zoo)  

test_that("FilterResultsLI::predict_level matches manual construction (daily, seasonal, no xpred)", {
  data(england, package = "tsgc")
  
  ## Leading indicator setup from the example in the class docs
  Y <- england[, 1:2]  # cumulative cases (leading) and cumulative hospitalisations (target)
  
  est.start <- as.Date("2021-04-30")
  est.end   <- as.Date("2021-07-24")
  n.lag     <- 4
  nf        <- 7
  
  mod <- SSModelLeadingIndicator(
    Y          = Y,
    n.lag      = n.lag,
    q          = NULL,
    LeadIndCol = 1,
    sea.period = 7,
    start.date = est.start,
    end.date   = est.end
  )
  
  res <- estimate(mod)   # this should return a FilterResultsLI object
  
  ## ------------------------------------------------------------------
  ## KFAS-based reference: get forecasts of log growth of the target
  ## ------------------------------------------------------------------
  # predict_all() here is used only to get the KFAS forecast object y.hat.kfas
  kfas_fc <- res$predict_all(
    n.ahead          = nf,
    sea.on           = TRUE,
    return.all       = FALSE,
    confidence.level = 0.68
  )$y.hat.kfas
  
  # Columns of y.hat.kfas$LDLtarg are: fit, lwr, upr
  delta_fit <- as.vector(kfas_fc$LDLtarg[, 1])
  delta_lwr <- as.vector(kfas_fc$LDLtarg[, 2])
  delta_upr <- as.vector(kfas_fc$LDLtarg[, 3])
  
  ## ------------------------------------------------------------------
  ## Recover last cumulative target value from data_xts
  ## ------------------------------------------------------------------
  last_row <- get_timeframe(res$data_xts, res$end.date)[1, ]
  YT       <- as.numeric(last_row$cTarg)  # last cumulative hospitalisations
  
  ## ------------------------------------------------------------------
  ## Manually replicate FilterResultsLI$predict_level() logic
  ## (Harvey/Kattuman multiplicative formula)
  ## ------------------------------------------------------------------
  # point forecasts
  cp_fit   <- cumprod(1 + exp(delta_fit[1:(nf - 1)]))
  mult_fit <- c(1, cp_fit)
  forc_fit <- YT * exp(delta_fit) * mult_fit
  
  # lower bound
  cp_lwr   <- cumprod(1 + exp(delta_lwr[1:(nf - 1)]))
  mult_lwr <- c(1, cp_lwr)
  forc_lwr <- YT * exp(delta_lwr) * mult_lwr
  
  # upper bound
  cp_upr   <- cumprod(1 + exp(delta_upr[1:(nf - 1)]))
  mult_upr <- c(1, cp_upr)
  forc_upr <- YT * exp(delta_upr) * mult_upr
  
  ## ------------------------------------------------------------------
  ## Compare with res$predict_level()
  ## ------------------------------------------------------------------
  forc_tsgc <- res$predict_level(
    n.ahead          = nf,
    confidence.level = 0.68,
    sea.on           = TRUE
  )
  
  # FilterResultsLI$predict_level returns an xts with columns forc, lwr, upr
  # It rounds forecasts to 2 decimal places, so we round the reference as well.
  expect_equal(as.numeric(forc_tsgc$forc), round(forc_fit,  2))
  expect_equal(as.numeric(forc_tsgc$lwr),  round(forc_lwr,  2))
  expect_equal(as.numeric(forc_tsgc$upr),  round(forc_upr,  2))
})

test_that("FilterResultsLI::get_growth_y toggles smoothed correctly", {
  data(england, package = "tsgc")
  
  Y <- england[, 1:2]
  
  est.start <- as.Date("2021-04-30")
  est.end   <- as.Date("2021-07-24")
  n.lag     <- 4
  
  mod <- SSModelLeadingIndicator(
    Y          = Y,
    n.lag      = n.lag,
    q          = NULL,
    LeadIndCol = 1,
    sea.period = 7,
    start.date = est.start,
    end.date   = est.end
  )
  res <- estimate(mod)
  
  smooth    <- res$get_growth_y(smoothed = TRUE,  return.components = TRUE)
  filt_only <- res$get_growth_y(smoothed = FALSE, return.components = FALSE)
  filt_all  <- res$get_growth_y(smoothed = FALSE, return.components = TRUE)
  
  expect_equal(length(smooth), 3)
  expect_false(is.list(filt_only))
  expect_equal(names(smooth[[1]]), "smoothed gy.t")
  expect_equal(names(smooth[[2]]), "smoothed g.t")
  expect_equal(names(smooth[[3]]), "smoothed gamma.t")
  expect_equal(names(filt_all[[1]]), "filtered gy.t")
  expect_equal(names(filt_all[[2]]), "filtered g.t")
  expect_equal(names(filt_all[[3]]), "filtered gamma.t")
  expect_false(isTRUE(all.equal(smooth[[1]], filt_only)))
})

test_that("FilterResultsLI::get_gy_ci toggles smoothed correctly", {
  data(england, package = "tsgc")
  
  Y <- england[, 1:2]
  
  est.start <- as.Date("2021-04-30")
  est.end   <- as.Date("2021-07-24")
  n.lag     <- 4
  
  mod <- SSModelLeadingIndicator(
    Y          = Y,
    n.lag      = n.lag,
    q          = NULL,
    LeadIndCol = 1,
    sea.period = 7,
    start.date = est.start,
    end.date   = est.end
  )
  res <- estimate(mod)
  
  smooth <- res$get_gy_ci(smoothed = TRUE)
  filt   <- res$get_gy_ci(smoothed = FALSE)
  
  expect_equal(colnames(smooth), c("fit", "lower", "upper"))
  expect_false(isTRUE(all.equal(smooth, filt)))
})

test_that("FilterResultsLI print/summary do not error", {
  data(england, package = "tsgc")
  Y <- england[, 1:2]
  
  est.start <- as.Date("2021-04-30")
  est.end   <- as.Date("2021-07-24")
  n.lag     <- 4
  
  mod <- SSModelLeadingIndicator(
    Y          = Y,
    n.lag      = n.lag,
    q          = NULL,
    LeadIndCol = 1,
    sea.period = 7,
    start.date = est.start,
    end.date   = est.end
  )
  res <- estimate(mod)
  
  expect_no_error(res$print_estimation_results())
  expect_no_error(print(res))
  expect_no_error(summary(res))
})

test_that("FilterResultsLI plotting methods do not error", {
  data(england, package = "tsgc")
  Y <- england[, 1:2]
  
  est.start <- as.Date("2021-04-30")
  est.end   <- as.Date("2021-07-24")
  n.lag     <- 4
  
  mod <- SSModelLeadingIndicator(
    Y          = Y,
    n.lag      = n.lag,
    q          = NULL,
    LeadIndCol = 1,
    sea.period = 7,
    start.date = est.start,
    end.date   = est.end
  )
  res <- estimate(mod)
  
  expect_no_error(res$plot_forecast())
  expect_no_error(res$plot_log_forecast(Y = Y))
  expect_no_error(res$plot_gy_ci())
  expect_no_error(res$plot_gy_components())
  expect_no_error(res$plot_holdout(Y = Y))
})

test_that("FilterResultsLI::mapes works", {
  data(england, package = "tsgc")
  Y <- england[, 1:2]
  
  est.start <- as.Date("2021-04-30")
  est.end   <- as.Date("2021-07-24")
  n.lag     <- 4
  
  mod <- SSModelLeadingIndicator(
    Y          = Y,
    n.lag      = n.lag,
    q          = NULL,
    LeadIndCol = 1,
    sea.period = 7,
    start.date = est.start,
    end.date   = est.end
  )
  res <- estimate(mod)
  
  errs <- res$mapes(n.ahead = 7, Y = Y)
  expect_equal(length(errs), 5)
})
