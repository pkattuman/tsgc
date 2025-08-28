library(KFAS)
library(xts)

test_that("tsgc produces same LI output as KFAS", {
  data(ukitaly, package = 'tsgc')
  
  est.start <- as.Date("2020-02-25") 
  est.end <- as.Date("2020-04-01")
  
  n.lag <- 14
  
  tsgc_mod <- SSModelLeadingIndicator(Y = ukitaly, n.lag = n.lag, q = NULL,
                                      sea.period = 0, start.date = est.start, 
                                      end.date = est.end, LeadIndCol = 1)
  tsgc_est <- tsgc_mod$estimate()
  
  idx <- (zoo::index(ukitaly) >= est.start) & (zoo::index(ukitaly) <= est.end)
  ukitaly_ldl <- df2ldl(ukitaly)[idx,]
  ukitaly_ldl$Italy <- lag(as.vector(ukitaly_ldl$Italy), n.lag)
  ukitaly_ldl <- na.omit(ukitaly_ldl)
  kfas_mod <- SSModel(as.matrix(ukitaly_ldl) ~ SSMtrend(degree = 2, 
                                        Q = matrix(c(0,0,0,NA),2,2),
                                        type = 'common') +
                   SSMtrend(degree = 1, Q = matrix(NA),index=1),
                 H = matrix(c(NA,0,0,NA),2,2))
  npar <- sum(is.na(kfas_mod$Q)) + sum(is.na(kfas_mod$H))
  kfas_fit <- fitSSM(kfas_mod, rep(0,npar))
  kfas_est <- KFS(kfas_fit$model)
  
  expect_equal(unname(as.matrix(tsgc_est$output$alphahat)), 
               unname(as.matrix(kfas_est$alphahat)))
})

test_that("tsgc produces same LI output as KFAS with fixed q", {
  data(ukitaly, package = 'tsgc')
  
  est.start <- as.Date("2020-02-25") 
  est.end <- as.Date("2020-04-01")
  
  n.lag <- 14
  
  tsgc_mod <- SSModelLeadingIndicator(Y = ukitaly, n.lag = n.lag, q = 0.005,
                                      sea.period = 0, start.date = est.start, 
                                      end.date = est.end, LeadIndCol = 1)
  tsgc_est <- tsgc_mod$estimate()
  
  updatesn=function(pars, model, snr, order, index){
    if(any(is.na(model$Q))){
      Q <- as.matrix(model$Q[,,1])
      naQd  <- which(is.na(diag(Q)))
      naQnd <- which(upper.tri(Q[naQd,naQd]) & is.na(Q[naQd,naQd]))
      Q[naQd,naQd][lower.tri(Q[naQd,naQd])] <- 0
      diag(Q)[naQd] <- exp(0.5 * pars[1:length(naQd)])
      Q[naQd,naQd][naQnd] <- pars[length(naQd)+1:length(naQnd)]
      model$Q[naQd,naQd,1] <- crossprod(Q[naQd,naQd])
    }
    if(!identical(model$H,'Omitted') && any(is.na(model$H))){
      H<-as.matrix(model$H[,,1])
      naHd  <- which(is.na(diag(H)))
      naHnd <- which(upper.tri(H[naHd,naHd]) & is.na(H[naHd,naHd]))
      H[naHd,naHd][lower.tri(H[naHd,naHd])] <- 0
      diag(H)[naHd] <-
        exp(0.5 * pars[length(naQd)+length(naQnd)+1:length(naHd)])
      H[naHd,naHd][naHnd] <-
        pars[length(naQd)+length(naQnd)+length(naHd)+1:length(naHnd)]
      model$H[naHd,naHd,1] <- crossprod(H[naHd,naHd])
      model$Q[order,order,1] <- snr*crossprod(H[index,index])
    }
    model
  }
  updateli = updatesn %>% purrr::partial(snr=0.005,order=2,index=2)
  
  idx <- (zoo::index(ukitaly) >= est.start) & (zoo::index(ukitaly) <= est.end)
  ukitaly_ldl <- df2ldl(ukitaly)[idx,]
  ukitaly_ldl$Italy <- lag(as.vector(ukitaly_ldl$Italy), n.lag)
  ukitaly_ldl <- na.omit(ukitaly_ldl)
  kfas_mod <- SSModel(as.matrix(ukitaly_ldl) ~ SSMtrend(degree = 2, 
                                                    Q = matrix(c(0,0,0,NA),2,2),
                                                    type = 'common') +
                        SSMtrend(degree = 1, Q = matrix(NA),index=1),
                      H = matrix(c(NA,0,0,NA),2,2))
  npar <- sum(is.na(kfas_mod$Q)) + sum(is.na(kfas_mod$H))
  kfas_fit <- fitSSM(kfas_mod, rep(0,npar), updatefn = updateli)
  kfas_est <- KFS(kfas_fit$model)
  
  tsgc_snr <- tsgc_est$output$model$Q[2,2,1]/tsgc_est$output$model$H[2,2,1]
  kfas_snr <- kfas_est$model$Q[2,2,1]/kfas_est$model$H[2,2,1]
  
  expect_equal(tsgc_snr, 0.005)
  expect_equal(kfas_snr, 0.005)
  expect_equal(unname(as.matrix(tsgc_est$output$alphahat)), 
               unname(as.matrix(kfas_est$alphahat)))
})

test_that("Summary method works", {
  data(ukitaly, package = 'tsgc')
  
  est.start <- as.Date("2020-02-25") 
  est.end <- as.Date("2020-04-01")
  
  n.lag <- 14
  
  tsgc_mod <- SSModelLeadingIndicator(Y = ukitaly, n.lag = n.lag, q = NULL,
                                      sea.period = 0, start.date = est.start, 
                                      end.date = est.end, LeadIndCol = 1)

  expect_no_error(expect_no_warning(tsgc_mod$summary()))
})

test_that("Print method works", {
  data(ukitaly, package = 'tsgc')
  
  est.start <- as.Date("2020-02-25") 
  est.end <- as.Date("2020-04-01")
  
  n.lag <- 14
  
  tsgc_mod <- SSModelLeadingIndicator(Y = ukitaly, n.lag = n.lag, q = NULL,
                                      sea.period = 0, start.date = est.start, 
                                      end.date = est.end, LeadIndCol = 1)
  
  expect_no_error(expect_no_warning(tsgc_mod$print()))
})

test_that("Plot method works", {
  data(ukitaly, package = 'tsgc')
  
  est.start <- as.Date("2020-02-25") 
  est.end <- as.Date("2020-04-01")
  
  n.lag <- 14
  
  tsgc_mod <- SSModelLeadingIndicator(Y = ukitaly, n.lag = n.lag, q = NULL,
                                      sea.period = 0, start.date = est.start, 
                                      end.date = est.end, LeadIndCol = 1)
  
  expect_no_error(suppressWarnings(tsgc_mod$plot()))
})

test_that("Leading indicator model + seasonal has correct number of elements", {
  eng <- tsgc::england[, 1:2]
  
  mod <- SSModelLeadingIndicator$new(eng, n.lag = 5)
})