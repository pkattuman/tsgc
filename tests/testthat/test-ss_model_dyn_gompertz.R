library(KFAS)

test_that("tsgc gives same output as KFAS", {
  data(england, package = 'tsgc')
  dat <- england$cum_cases[1:100]
  dat.ldl <- na.omit(df2ldl(dat))
  
  tsgc.model <- SSModelDynamicGompertz$new(Y = dat,
                                       q = NULL,
                                       sea.period = 0)
  tsgc.est <- tsgc.model$estimate()
  final.states.tsgc <- tsgc.est$output$alphahat[100,]
  
  kfas.model <- SSModel(as.matrix(dat.ldl) ~
                                SSMtrend(degree = 2, 
                                               Q = list(matrix(0), matrix(NA))), 
                              H = matrix(NA))
  kfas.fit <- fitSSM(kfas.model,inits=c(0,0))
  kfas.out <- KFS(kfas.fit$model)
  final.states.kfas <- kfas.out$alphahat[99,]
  
  expect_equal(final.states.tsgc,final.states.kfas,tolerance=0.01)
})

test_that({"tsgc enforces signal-to-noise restrictions correctly"}, {
  data(england, package = 'tsgc')
  dat <- england$cum_cases[1:100]
  
  q.choose <- 0.005
  
  tsgc.model <- SSModelDynamicGompertz$new(Y = dat,
                                           q = q.choose,
                                           sea.period = 0)
  tsgc.est <- tsgc.model$estimate()
  
  q.out <- as.vector(tsgc.est$output$model$Q[2,2,1]/tsgc.est$output$model$H)
  
  expect_equal(q.out,q.choose)
})

test_that({"tsgc summary method functions"}, {
  data(england, package = 'tsgc')
  dat <- england$cum_cases[1:100]
  
  q.choose <- 0.005
  
  tsgc.model <- SSModelDynamicGompertz$new(Y = dat,
                                           q = q.choose,
                                           sea.period = 0)
  tsgc.model$summary()
})

test_that({"tsgc print method functions"}, {
  data(england, package = 'tsgc')
  dat <- england$cum_cases[1:100]
  
  q.choose <- 0.005
  
  tsgc.model <- SSModelDynamicGompertz$new(Y = dat,
                                           q = q.choose,
                                           sea.period = 0)
  tsgc.model$print()
})

test_that({"tsgc plot method functions"}, {
  data(england, package = 'tsgc')
  dat <- england$cum_cases[1:100]
  
  q.choose <- 0.005
  
  tsgc.model <- SSModelDynamicGompertz$new(Y = dat,
                                           q = q.choose,
                                           sea.period = 0)
  suppressWarnings(tsgc.model$plot())
})

test_that({"tsgc enforces signal-to-noise restrictions correctly"}, {
  data(england, package = 'tsgc')
  dat <- england$cum_cases[1:100]
  
  q.choose <- 0.005
  
  tsgc.model <- SSModelDynamicGompertz$new(Y = dat,
                                           q = q.choose,
                                           sea.period = 0)
  tsgc.est <- tsgc.model$estimate()
  
  q.out <- as.vector(tsgc.est$output$model$Q[2,2,1]/tsgc.est$output$model$H)
  
  expect_equal(q.out,q.choose)
})


test_that({"Model with seasonal has expected number of seasonal components"},{
  data(england, package = 'tsgc')
  dat <- england$cum_cases[1:100]
  
  sea.choose <- 7
  
  tsgc.model <- SSModelDynamicGompertz$new(Y = dat,
                                           q = 0.005,
                                           sea.period = sea.choose)
  tsgc.est <- tsgc.model$estimate()
  sea.est <- length(tsgc.est$output$model["a1","seasonal"])
  
  expect_equal(sea.est,2*floor(sea.choose/2))
})
