library(dplyr, warn.conflicts = FALSE)
library(tidyr)
library(splines)
library(survival)

testthat::test_that("msm: pooled logistic regression switching model", {
  sim1 <- tssim(
    tdxo = 0, coxo = 0, p_R = 0.5, p_X_1 = 0.3, p_X_0 = 0.3, 
    rate_T = 0.002, beta1 = -0.5, beta2 = 0.3, 
    gamma0 = 0.3, gamma1 = -0.9, gamma2 = 0.7, gamma3 = 1.1, gamma4 = -0.8,
    zeta0 = -3.5, zeta1 = 0.5, zeta2 = 0.2, zeta3 = -0.4, 
    alpha0 = 0.5, alpha1 = 0.5, alpha2 = 0.4, 
    theta1_1 = -0.4, theta1_0 = -0.4, theta2 = 0.2,
    rate_C = 0.0000855, followup = 20, days = 30,
    n = 500, NSim = 100, seed = 314159)
  
  fit1 <- msm(
    sim1[[1]], id = "id", tstart = "tstart", 
    tstop = "tstop", event = "Y", treat = "trtrand", 
    swtrt = "xo", swtrt_time = "xotime", base_cov = "bprog", 
    numerator = "bprog", denominator = c("bprog", "L"), 
    ns_df = 3, swtrt_control_only = TRUE, boot = FALSE)
  
  # exclude observations after treatment switch
  data1 <- sim1[[1]] %>%
    filter(xo == 0 | tstart < xotime) %>%
    group_by(id) %>%
    mutate(condition = row_number() == n() & xo == 1 & tstop >= xotime,
           cross = ifelse(condition, 1, 0),
           tstop = ifelse(condition, xotime, tstop))
  
  # fit pooled logistic regression switching models
  data2 <- data1 %>% filter(trtrand == 0)
  
  ns1 <- ns(data2$tstop[data2$cross == 1], df = 3)
  ns2 <- ns(data2$tstop, knots = attr(ns1, "knots"), 
            Boundary.knots = attr(ns1, "Boundary.knots"))
  switch1 <- glm(cross ~ bprog + L + ns2, family = binomial, 
                 data = data2)
  h_den <- as.numeric(predict(switch1, newdata = data2, type = "response"))
  
  switch2 <- glm(cross ~ bprog + ns2, family = binomial, data = data2)
  h_num <- as.numeric(predict(switch2, newdata = data2, type = "response"))
  
  # stablized weights
  data3 <- data2 %>% 
    ungroup() %>%
    mutate(pstay_den = 1-h_den, pstay_num = 1-h_num) %>%
    mutate(p0_den = ifelse(cross == 1, h_den, pstay_den),
           p0_num = ifelse(cross == 1, h_num, pstay_num)) %>%
    group_by(id) %>%
    mutate(p_den = cumprod(p0_den), p_num = cumprod(p0_num),
           stabilized_weight = p_num/p_den) %>%
    select(id, tstart, tstop, Y, stabilized_weight, bprog, trtrand)
  
  data4 <- data1 %>%
    left_join(data3 %>% 
                select(id, tstop, stabilized_weight), 
              by = c("id", "tstop")) %>%
    mutate(sw = ifelse(is.na(stabilized_weight), 1, stabilized_weight))
  
  
  # set up time-dependent treatment switching indicators
  dt1 <- sim1[[1]] %>%
    mutate(cross = ifelse(xo == 1 & tstop >= xotime, 1, 0))
  
  # weights do not change after treatment switch
  dt2 <- dt1 %>% 
    filter(trtrand == 0) %>%
    left_join(data4 %>% select(id, tstop, sw), 
              by = c("id", "tstop")) %>%
    group_by(id) %>%
    arrange(tstop, .by_group = TRUE) %>%
    fill(sw, .direction = "down") %>%   # LOCF
    bind_rows(dt1 %>% 
                filter(trtrand == 1) %>%
                mutate(sw = 1) %>%
                select(id, tstart, tstop, Y, sw, bprog, trtrand, cross))
  
  fit <- coxph(Surv(tstart, tstop, Y) ~ trtrand + bprog + cross, 
               data = dt2, weight = sw,
                id = id, ties = "efron", robust = TRUE)
  
  hr1 <- as.numeric(exp(cbind(fit$coefficients, confint(fit)))["trtrand",])
  
  testthat::expect_equal(dt2$sw, fit1$data_outcome$stabilized_weight)
  
  testthat::expect_equal(hr1, c(fit1$hr, fit1$hr_CI))
})

