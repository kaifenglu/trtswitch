library(dplyr, warn.conflicts = FALSE)
library(survival)

test_that("liferegr: eligible distributions", {
  for (dist in c("exponential", "weibull", "lognormal", "loglogistic")) {
    fit1 <- liferegr(ovarian, time="futime", event="fustat",
                     covariates=c("ecog.ps", "rx"), dist=dist)
    
    fit2 <- survreg(Surv(futime, fustat) ~ ecog.ps + rx, data=ovarian, 
                    dist=dist)
    
    expect_equal(fit1$beta[1:(fit1$nvar+1)], fit2$coefficients)
    if (dist != "exponential") {
      expect_equal(exp(fit1$parest$beta[(fit1$nvar+2)]), fit2$scale)
    }
  }
})


test_that("liferegr: left-censored data", {
  tobin <- tobin %>% mutate(time = ifelse(durable > 0, durable, NA))
  
  for (dist in c("gaussian", "logistic")) {
    fit1 <- liferegr(tobin, time="time", time2="durable",
                     covariates=c("age", "quant"), dist=dist)
    
    fit2 <- survreg(Surv(durable, durable>0, type='left') ~ age + quant, 
                    data=tobin, dist=dist)
    
    expect_equal(fit1$beta[1:(fit1$nvar+1)], fit2$coefficients)
    expect_equal(exp(fit1$parest$beta[(fit1$nvar+2)]), fit2$scale)
  }
})


test_that("liferegr: stratification", {
  lung1 <- lung %>% mutate(event = 1*(status == 2))
  
  fit1 <- liferegr(lung1, stratum="sex", time="time", event="event",
                   covariates=c("ph.ecog", "age"))
  
  fit2 <- survreg(Surv(time, status) ~ ph.ecog + age + strata(sex), 
                  data=lung)
  
  expect_equal(fit1$beta[1:(fit1$nvar+1)], fit2$coefficients)
  expect_equal(exp(fit1$parest$beta[(fit1$nvar+2):(fit1$p)]), 
               as.vector(fit2$scale))
})


test_that("liferegr: robust variance", {
  diabetic <- diabetic %>% mutate(juvenile = 1*(age < 20))
  
  fit1 <- liferegr(diabetic, time="time", event="status", 
                   covariates=c("trt", "juvenile"), id="id", 
                   robust=TRUE)
  fit2 <- survreg(Surv(time, status) ~ trt + juvenile, cluster=id, 
                  data=diabetic, robust=TRUE)
  
  expect_equal(fit1$beta[1:(fit1$nvar+1)], fit2$coefficients)
  expect_equal(exp(fit1$parest$beta[(fit1$nvar+2)]), fit2$scale)
  expect_equal(fit1$vbeta, fit2$var)
})
