library(dplyr, warn.conflicts = FALSE)
library(survival)

testthat::test_that("tsesimp: weibull aft", {
  # modify pd and dpd based on co and dco
  shilong <- shilong %>%
    mutate(dpd = ifelse(co & !pd, dco, dpd),
           pd = ifelse(co & !pd, 1, pd)) %>%
    mutate(dpd = ifelse(pd & co & dco < dpd, dco, dpd))
  
  # the eventual survival time
  shilong1 <- shilong %>%
    arrange(bras.f, id, tstop) %>%
    group_by(bras.f, id) %>%
    slice(n()) %>%
    select(-c("ps", "ttc", "tran"))
  
  # the last value of time-dependent covariates before pd
  shilong2 <- shilong %>%
    filter(pd == 0 | tstart <= dpd) %>%
    arrange(bras.f, id, tstop) %>%
    group_by(bras.f, id) %>%
    slice(n()) %>%
    select(bras.f, id, ps, ttc, tran)
  
  # combine baseline and time-dependent covariates
  shilong3 <- shilong1 %>%
    left_join(shilong2, by = c("bras.f", "id"))
  
  # apply the two-stage method
  fit1 <- tsesimp(
    data = shilong3, id = "id", time = "tstop", event = "event",
    treat = "bras.f", censor_time = "dcut", pd = "pd",
    pd_time = "dpd", swtrt = "co", swtrt_time = "dco",
    base_cov = c("agerand", "sex.f", "tt_Lnum", "rmh_alea.c",
                 "pathway.f"),
    base2_cov = c("agerand", "sex.f", "tt_Lnum", "rmh_alea.c",
                  "pathway.f", "ps", "ttc", "tran"),
    aft_dist = "weibull", alpha = 0.05,
    recensor = TRUE, swtrt_control_only = FALSE, offset = 1,
    boot = FALSE)

  # numeric code of treatment and apply administrative censoring
  data1 <- shilong3 %>% 
    mutate(treated = 1*(bras.f == "MTA"),
           swtrt = 1*co)
  
  tablist <- lapply(0:1, function(h) {
    df1 <- data1 %>% 
      filter(treated == h & pd == 1) %>%
      mutate(time = tstop - dpd + 1)
    
    fit_aft <- survreg(Surv(time, event) ~ swtrt + agerand + sex.f + 
                         tt_Lnum + rmh_alea.c + pathway.f + 
                         ps + ttc + tran, data = df1)
    
    psi = -fit_aft$coefficients[2]

    data1 %>% 
      filter(treated == h) %>%
      mutate(u_star = ifelse(swtrt==1, dpd-1 + (tstop-dpd+1)*exp(psi), tstop),
             c_star = pmin(dcut, dcut*exp(psi)),
             t_star = pmin(u_star, c_star),
             d_star = ifelse(c_star < u_star, 0, event))
  })
  
  data2 <- do.call(rbind, tablist)
  
  fit <- coxph(Surv(t_star, d_star) ~ treated + agerand + sex.f + 
                 tt_Lnum + rmh_alea.c + pathway.f, data = data2)
    
  hr1 <- as.numeric(exp(cbind(fit$coefficients, confint(fit)))["treated",])
  testthat::expect_equal(hr1, c(fit1$hr, fit1$hr_CI))
})


testthat::test_that("tsesimp: boot", {
  # modify pd and dpd based on co and dco
  shilong <- shilong %>%
    mutate(dpd = ifelse(co & !pd, dco, dpd),
           pd = ifelse(co & !pd, 1, pd)) %>%
    mutate(dpd = ifelse(pd & co & dco < dpd, dco, dpd))
  
  # the eventual survival time
  shilong1 <- shilong %>%
    arrange(bras.f, id, tstop) %>%
    group_by(bras.f, id) %>%
    slice(n()) %>%
    select(-c("ps", "ttc", "tran"))
  
  # the last value of time-dependent covariates before pd
  shilong2 <- shilong %>%
    filter(pd == 0 | tstart <= dpd) %>%
    arrange(bras.f, id, tstop) %>%
    group_by(bras.f, id) %>%
    slice(n()) %>%
    select(bras.f, id, ps, ttc, tran)
  
  # combine baseline and time-dependent covariates
  shilong3 <- shilong1 %>%
    left_join(shilong2, by = c("bras.f", "id"))
  
  fit2 <- tsesimp(
    data = shilong3, id = "id", time = "tstop", event = "event",
    treat = "bras.f", censor_time = "dcut", pd = "pd",
    pd_time = "dpd", swtrt = "co", swtrt_time = "dco",
    base_cov = c("agerand", "sex.f", "tt_Lnum", "rmh_alea.c",
                 "pathway.f"),
    base2_cov = c("agerand", "sex.f", "tt_Lnum", "rmh_alea.c",
                  "pathway.f", "ps", "ttc", "tran"),
    aft_dist = "weibull", alpha = 0.05,
    recensor = TRUE, swtrt_control_only = FALSE, offset = 1,
    boot = TRUE, n_boot = 100, seed = 12345)
  
  hr2 <- c(0.9125816, 0.5485876, 1.5180896)
  testthat::expect_equal(hr2, c(fit2$hr, fit2$hr_CI))
})
