library(dplyr, warn.conflicts = FALSE)
library(survival)

test_that("rpsftm: control to active switch", {
  data1 <- immdef %>%
    mutate(rx = 1-xoyrs/progyrs) %>%
    arrange(imm)
  
  fit1 <- rpsftm(
    data1, time = "progyrs", event = "prog", treat = "imm", 
    rx = "rx", censor_time = "censyrs", boot = FALSE)
  
  # log-rank for ITT
  fit_lr <- survdiff(Surv(progyrs, prog) ~ imm, data = data1)
  z_lr = (fit_lr$obs - fit_lr$exp)[2]/sqrt(fit_lr$var[2,2])
  
  f <- function(psi) {
    data1 %>%
      mutate(u_star = progyrs*((1-rx) + rx*exp(psi)),
             c_star = ifelse(imm==0, pmin(censyrs, censyrs*exp(psi)), 1e10),
             time_ts = pmin(u_star, c_star),
             event_ts = prog*(u_star <= c_star))
  }
  
  g <- function(psi) {
    data2 <- f(psi)
    fit_lr <- survdiff(Surv(time_ts, event_ts) ~ imm, data = data2)
    (fit_lr$obs - fit_lr$exp)[2]/sqrt(fit_lr$var[2,2])
  }
  
  # psi based on log-rank test
  psi <- uniroot(g, c(-3,3), tol = 1e-6)$root
  
  # observed on treated and counterfactual on control
  data2 <- data1 %>%
    filter(imm == 0) %>%
    mutate(u_star = xoyrs + (progyrs - xoyrs)*exp(psi),
           c_star = pmin(censyrs, censyrs*exp(psi)),
           time_ts = pmin(u_star, c_star),
           event_ts = prog*(u_star <= c_star)) %>%
    select(-c("u_star", "c_star")) %>%
    bind_rows(data1 %>%
                filter(imm == 1) %>%
                mutate(time_ts = progyrs, 
                       event_ts = prog))
  
  fit <- coxph(Surv(time_ts, event_ts) ~ imm, data = data2)
  beta = as.numeric(fit$coefficients[1])
  se = beta/z_lr
  
  zcrit = qnorm(0.975)
  hr1 <- exp(c(beta, beta - zcrit*se, beta + zcrit*se))
  expect_equal(hr1, c(fit1$hr, fit1$hr_CI))
})
