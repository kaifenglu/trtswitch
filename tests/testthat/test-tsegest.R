library(dplyr, warn.conflicts = FALSE)
library(geepack)

testthat::test_that("tsegest: logistic g-estimation", {
  sim1 <- tsegestsim(
    n = 500, allocation1 = 2, allocation2 = 1, pbprog = 0.5, 
    trtlghr = -0.5, bprogsl = 0.3, shape1 = 1.8, 
    scale1 = 0.000025, shape2 = 1.7, scale2 = 0.000015, 
    pmix = 0.5, admin = 5000, pcatnotrtbprog = 0.5, 
    pcattrtbprog = 0.25, pcatnotrt = 0.2, pcattrt = 0.1, 
    catmult = 0.5, tdxo = 1, ppoor = 0.1, pgood = 0.04, 
    ppoormet = 0.4, pgoodmet = 0.2, xomult = 1.4188308, 
    milestone = 546, swtrt_control_only = TRUE,
    outputRawDataset = 1, seed = 2000)
  
  fit1 <- tsegest(
    data = sim1$paneldata, id = "id", 
    tstart = "tstart", tstop = "tstop", event = "died", 
    treat = "trtrand", censor_time = "censor_time", 
    pd = "progressed", pd_time = "timePFSobs", swtrt = "xo", 
    swtrt_time = "xotime", swtrt_time_upper = "xotime_upper",
    base_cov = "bprog", conf_cov = "bprog*catlag", 
    low_psi = -3, hi_psi = 3, strata_main_effect_only = TRUE,
    recensor = TRUE, admin_recensor_only = TRUE, 
    swtrt_control_only = TRUE, alpha = 0.05, ties = "efron", 
    tol = 1.0e-6, boot = FALSE)
  
  # last observation within each subject
  data1 <- sim1$paneldata %>%
    group_by(id) %>%
    slice(n()) %>%
    mutate(os = died, ostime = tstop) %>% 
    select(id, trtrand, progressed, timePFSobs, os, ostime, 
           censor_time, xo, xotime, bprog)
  
  # post progression data up to switching 
  data2 <- sim1$paneldata %>%
    filter(trtrand == 0 & progressed == 1 & tstop >= timePFSobs & 
             ((xo == 0 & tstop <= xotime_upper) | 
                (xo == 1 & tstop <= xotime)))
  
  # re-baseline
  data3a <- data2 %>% 
    left_join(data1 %>% 
                select(id, os, ostime), 
              by = "id") %>%
    group_by(id) %>%
    slice(1) %>%
    mutate(ostime = ostime - timePFSobs + 1,
           censor_time = censor_time - timePFSobs + 1,
           xotime = xotime - timePFSobs + 1) %>%
    ungroup()
  
  # setup treatment switching indicators
  data3b <- data2 %>%
    left_join(data1 %>% 
                select(id, os, ostime), 
              by = "id") %>%
    mutate(y = ifelse(xo == 1 & tstop == xotime, 1, 0))
  
  f <- function(psi, target) {
    data4a <- data3a %>%
      mutate(u_star = ifelse(xo == 1, xotime - 1 + 
                               (ostime - xotime + 1)*exp(psi), 
                             ostime), 
             c_star = pmin(censor_time, censor_time*exp(psi)),
             t_star = pmin(u_star, c_star),
             d_star = os*(u_star <= c_star))
    
    fit_cox <- phregr(data4a, time = "t_star", event = "d_star", 
                      est_resid = TRUE)
    resid <- fit_cox$residuals
    
    data4b <- data3b %>%
      left_join(data4a %>% 
                  mutate(resid = resid) %>%
                  select(id, resid), by = "id")
    
    fit_lgs <- geeglm(y ~ resid + bprog*catlag,
                      family = binomial, data = data4b, id = id)
    
    z_lgs <- coefficients(fit_lgs)/sqrt(diag(vcov(fit_lgs)))
    
    as.numeric(z_lgs["resid"]) - target
  }
  
  psi <- uniroot(f, c(-3, 3), 0, tol = 1.0e-6)$root
  
  data4 <- data1 %>%
    filter(trtrand == 0) %>%
    mutate(u_star = ifelse(xo == 1, xotime - 1 + 
                             (ostime - xotime + 1)*exp(psi), 
                           ostime), 
           c_star = pmin(censor_time, censor_time*exp(psi)),
           t_star = pmin(u_star, c_star),
           d_star = os*(u_star <= c_star)) %>%
    select(id, t_star, d_star, trtrand, bprog) %>%
    bind_rows(data1 %>% 
                filter(trtrand == 1) %>%
                mutate(t_star = ostime, d_star = os) %>%
                select(id, t_star, d_star, trtrand, bprog))
  
  fit <- phregr(data4, time = "t_star", event = "d_star", 
                covariates = c("trtrand", "bprog"), ties = "efron")
  
  hr1 <- exp(as.numeric(c(fit$parest$beta[1], fit$parest$lower[1], 
                          fit$parest$upper[1])))
  testthat::expect_equal(hr1, c(fit1$hr, fit1$hr_CI))
})
