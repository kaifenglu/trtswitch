library(dplyr, warn.conflicts = FALSE)

testthat::test_that("tsesimp: weibull aft", {
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
    data = shilong3, time = "tstop", event = "event",
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
    
    fit_aft <- liferegr(df1, time = "time", event = "event", 
                        covariates = c("swtrt", "agerand", "sex.f", 
                                       "tt_Lnum", "rmh_alea.c", "pathway.f",
                                       "ps", "ttc", "tran"), 
                        dist = "weibull")
    
    psi = -fit_aft$beta[2]

    data1 %>% 
      filter(treated == h) %>%
      mutate(u_star = ifelse(swtrt == 1, 
                             ifelse(pd == 1, dpd-1 + (tstop-dpd+1)*exp(psi), 
                                    dco-1 + (tstop-dco+1)*exp(psi)), tstop),
             c_star = pmin(dcut, dcut*exp(psi)),
             t_star = pmin(u_star, c_star),
             d_star = event*(u_star <= c_star))
  })
  
  data2 <- do.call(rbind, tablist)
  
  fit <- phregr(data2, time = "t_star", event = "d_star", 
                covariates = c("treated", "agerand", "sex.f", 
                               "tt_Lnum", "rmh_alea.c", "pathway.f"))
  
  hr1 <- exp(as.numeric(c(fit$parest$beta[1], fit$parest$lower[1], 
                          fit$parest$upper[1])))
  testthat::expect_equal(hr1, c(fit1$hr, fit1$hr_CI))
})
