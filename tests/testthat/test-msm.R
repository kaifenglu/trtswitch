library(dplyr, warn.conflicts = FALSE)
library(tidyr)
library(splines)
library(survival)

testthat::test_that("msm: pooled logistic regression switching model", {
  sim1 <- tssim(
    tdxo = 1, coxo = 1, allocation1 = 1, allocation2 = 1,
    p_X_1 = 0.3, p_X_0 = 0.3, 
    rate_T = 0.002, beta1 = -0.5, beta2 = 0.3, 
    gamma0 = 0.3, gamma1 = -0.9, gamma2 = 0.7, gamma3 = 1.1, gamma4 = -0.8,
    zeta0 = -3.5, zeta1 = 0.5, zeta2 = 0.2, zeta3 = -0.4, 
    alpha0 = 0.5, alpha1 = 0.5, alpha2 = 0.4, 
    theta1_1 = -0.4, theta1_0 = -0.4, theta2 = 0.2,
    rate_C = 0.0000855,  accrualIntensity = 20/30, 
    fixedFollowup = FALSE, plannedTime = 1350, days = 30,
    n = 500, NSim = 100, seed = 314159)
  
  data1 <- sim1[[1]] %>% 
    mutate(os = died, ostime = timeOS, event = Y)
  
  fit1 <- msm(
    data1, id = "id", tstart = "tstart", 
    tstop = "tstop", event = "Y", treat = "trtrand", 
    swtrt = "xo", swtrt_time = "xotime", 
    base_cov = "bprog", numerator = "bprog", 
    denominator = c("bprog", "L"), 
    ns_df = 3, swtrt_control_only = TRUE, boot = FALSE)
  
  # exclude observations after treatment switch
  data2 <- data1 %>%
    filter(ifelse(xo == 1, tstart < xotime, tstop < ostime)) %>%
    group_by(id) %>%
    mutate(condition = row_number() == n() & xo == 1 & tstop >= xotime,
           cross = ifelse(condition, 1, 0)) %>%
    ungroup()
  
  # fit pooled logistic regression switching models
  data3 <- data2 %>% filter(trtrand == 0)
  
  s1 <- ns(data3$tstop[data3$cross == 1], df = 3)
  s2 <- ns(data3$tstop, knots = attr(s1, "knots"), 
           Boundary.knots = attr(s1, "Boundary.knots"))
  switch1 <- glm(cross ~ bprog + L + s2, 
                 family = binomial, data = data3)
  phat1 <- as.numeric(predict(switch1, newdata = data3, type = "response"))
  
  switch2 <- glm(cross ~ bprog + s2, family = binomial, data = data3)
  phat2 <- as.numeric(predict(switch2, newdata = data3, type = "response"))
  
  # stabilized weights and time-dependent covariates
  data4 <- data1 %>% 
    filter(trtrand == 0) %>% 
    select(id, trtrand, bprog, tstart, tstop, event, xo, xotime) %>%
    left_join(data3 %>% 
                mutate(tstart = tstop,
                       o_den = ifelse(cross, phat1, 1-phat1), 
                       o_num = ifelse(cross, phat2, 1-phat2)) %>%
                select(id, tstart, o_den, o_num), 
              by = c("id", "tstart")) %>%
    mutate(o_den = ifelse(is.na(o_den), 1, o_den),
           o_num = ifelse(is.na(o_num), 1, o_num)) %>%
    group_by(id) %>%
    mutate(p_den = cumprod(o_den), p_num = cumprod(o_num),
           stabilized_weight = p_num/p_den) %>%
    select(id, tstart, tstop, event, stabilized_weight, bprog, 
           trtrand, xo, xotime) %>%
    ungroup() %>%
    bind_rows(data1 %>% 
                filter(trtrand == 1) %>%
                mutate(stabilized_weight = 1) %>%
                select(id, tstart, tstop, event, stabilized_weight, bprog, 
                       trtrand, xo, xotime)) %>%
    mutate(cross = ifelse(xo == 1, tstart >= xotime, 0))
  
  fit <- coxph(Surv(tstart, tstop, event) ~ trtrand + bprog + cross, 
               data = data4, weight = stabilized_weight,
                id = id, ties = "efron", robust = TRUE)
  
  hr1 <- as.numeric(exp(cbind(fit$coefficients, confint(fit)))["trtrand",])
  
  testthat::expect_equal(data4$stabilized_weight, 
                         fit1$data_outcome$stabilized_weight)
  
  testthat::expect_equal(hr1, c(fit1$hr, fit1$hr_CI))
})

