library(dplyr, warn.conflicts = FALSE)
library(survival)

testthat::test_that("rpsftm: control to active switch", {
  data1 <- immdef %>% mutate(rx = 1-xoyrs/progyrs)
  
  fit1 <- rpsftm(
    data1, time = "progyrs", event = "prog", treat = "imm", 
    rx = "rx", censor_time = "censyrs", boot = FALSE)
  
  # log-rank for ITT
  fit_lr <- survdiff(Surv(progyrs, prog) ~ imm, data = data1)
  z_lr <- sqrt(fit_lr$chisq)
  if (fit_lr$obs[2] < fit_lr$exp[2]) z_lr <- -z_lr
  
  f <- function(psi) {
    data1 %>%
      mutate(u_star = progyrs*((1-rx) + rx*exp(psi)),
             c_star = ifelse(imm==0, pmin(censyrs, censyrs*exp(psi)), 1e10),
             t_star = pmin(u_star, c_star),
             d_star = prog*(u_star <= c_star))
  }
  
  g <- function(psi) {
    data2 <- f(psi)
    fit_lr <- survdiff(Surv(t_star, d_star) ~ imm, data = data2)
    z_lr <- sqrt(fit_lr$chisq)
    if (fit_lr$obs[2] < fit_lr$exp[2]) z_lr <- -z_lr
    z_lr
  }
  
  # psi based on log-rank test
  psi <- uniroot(g, c(-1,1), tol = 1e-6)$root
  
  # observed on treated and counterfactual on control
  data2 <- data1 %>%
    filter(imm == 0) %>%
    mutate(u_star = xoyrs + (progyrs - xoyrs)*exp(psi),
           c_star = pmin(censyrs, censyrs*exp(psi)),
           t_star = pmin(u_star, c_star),
           d_star = prog*(u_star <= c_star)) %>%
    select(-c("u_star", "c_star")) %>%
    bind_rows(data1 %>%
                filter(imm == 1) %>%
                mutate(t_star = progyrs, 
                       d_star = prog))
  
  fit <- phregr(data2, time = "t_star", event = "d_star", covariates = "imm")
  beta = as.numeric(fit$beta[1])
  se = beta/z_lr
  
  zcrit = qnorm(0.975)
  hr1 <- exp(c(beta, beta - zcrit*se, beta + zcrit*se))
  testthat::expect_equal(hr1, c(fit1$hr, fit1$hr_CI))
})
