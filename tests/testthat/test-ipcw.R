library(dplyr, warn.conflicts = FALSE)
library(splines)
library(survival)

testthat::test_that("ipcw: pooled logistic regression switching model", {
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
    mutate(ostime = timeOS)
  
  fit1 <- ipcw(
    data1, id = "id", tstart = "tstart", 
    tstop = "tstop", event = "event", treat = "trtrand", 
    swtrt = "xo", swtrt_time = "xotime", 
    base_cov = "bprog", numerator = "bprog", 
    denominator = c("bprog", "L"),
    logistic_switching_model = TRUE, ns_df = 3,
    swtrt_control_only = TRUE, boot = FALSE)
  
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
  
  # stablized weights
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
    filter(ifelse(xo == 1, tstart < xotime, TRUE)) %>%
    select(id, tstart, tstop, event, stabilized_weight, bprog, trtrand) %>%
    ungroup() %>%
    bind_rows(data1 %>% 
                filter(trtrand == 1) %>%
                mutate(stabilized_weight = 1) %>%
                select(id, tstart, tstop, event, 
                       stabilized_weight, bprog, trtrand))
  
  fit <- coxph(Surv(tstart, tstop, event) ~ trtrand + bprog, 
               data = data4, weight = stabilized_weight,
                id = id, ties = "efron", robust = TRUE)
  
  hr1 <- as.numeric(exp(cbind(fit$coefficients, confint(fit)))["trtrand",])
  
  testthat::expect_equal(data4$stabilized_weight, 
                         fit1$data_outcome$stabilized_weight)
  
  testthat::expect_equal(hr1, c(fit1$hr, fit1$hr_CI))
})


testthat::test_that("ipcw: time-dependent covariates Cox switching model", {
  fit2 <- ipcw(
    shilong, id = "id", tstart = "tstart", tstop = "tstop", 
    event = "event", treat = "bras.f", swtrt = "co", 
    swtrt_time = "dco", 
    base_cov = c("agerand", "sex.f", "tt_Lnum", "rmh_alea.c", "pathway.f"),
    numerator = c("agerand", "sex.f", "tt_Lnum", "rmh_alea.c", "pathway.f"),
    denominator = c("agerand", "sex.f", "tt_Lnum", "rmh_alea.c", "pathway.f",
                    "ps", "ttc", "tran"),
    swtrt_control_only = FALSE, boot = FALSE)
  
  # exclude observations after treatment switch
  data1 <- shilong %>%
    mutate(treated = ifelse(bras.f == "CT", 0, 1)) %>%
    filter(ifelse(co, tstart < dco, TRUE)) %>%
    arrange(treated, id) %>%
    group_by(treated, id) %>%
    mutate(condition = row_number() == n() & co & tstop >= dco,
           cross = ifelse(condition, 1, 0),
           event = ifelse(condition, 0, event),
           tstop = ifelse(condition, dco, tstop)) %>%
    ungroup()
  
  # replicate event times within each subject
  cut <- sort(unique(data1$tstop[data1$event == 1]))
  a1 <- survsplit(data1$tstart, data1$tstop, cut)
  data2 <- data1[a1$row+1,]
  data2$tstart = a1$start
  data2$tstop = a1$end
  data2$event[a1$censor] = 0
  data2$cross[a1$censor] = 0
  
  tablist <- lapply(0:1, function(h) {
    df1 <- data2 %>% 
      filter(treated == h)
    fit_den <- coxph(Surv(tstart, tstop, cross) ~ agerand + sex.f + 
                       tt_Lnum + rmh_alea.c + pathway.f + ps + ttc + tran, 
                     data = df1, id = id, ties = "efron", robust = TRUE)
    fit_num <- coxph(Surv(tstart, tstop, cross) ~ agerand + sex.f + 
                       tt_Lnum + rmh_alea.c + pathway.f, 
                     data = df1, id = id, ties = "efron", robust = TRUE)
    
    km_den <- survfit(fit_den, newdata = df1, id = id)
    km_num <- survfit(fit_num, newdata = df1, id = id)
    
    id <- rep(as.numeric(names(km_den$strata)), km_den$strata)
    
    tibble(id = id, tstop = km_den$time, 
           surv_den = km_den$surv, surv_num = km_num$surv)
  })
  
  data3 <- do.call(rbind, tablist)
  
  data4 <- data2 %>% 
    left_join(data3, by = c("id", "tstop")) %>%
    mutate(stabilized_weight = surv_num/surv_den)
  
  fit <- coxph(Surv(tstart, tstop, event) ~ treated + agerand + 
                 sex.f + tt_Lnum + rmh_alea.c + pathway.f, 
               data = data4, weight = stabilized_weight, 
               id = id, ties = "efron", robust = TRUE)
  
  hr1 <- as.numeric(exp(cbind(fit$coefficients, confint(fit)))["treated",])
  
  testthat::expect_equal(data4$stabilized_weight, 
                         fit2$data_outcome$stabilized_weight)
  
  testthat::expect_equal(hr1, c(fit2$hr, fit2$hr_CI))
})
