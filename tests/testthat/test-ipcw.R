library(dplyr, warn.conflicts = FALSE)
library(splines)
library(survival)

testthat::test_that("ipcw: pooled logistic regression switching model", {
  sim1 <- tsegestsim(
    n = 500, allocation1 = 2, allocation2 = 1, pbprog = 0.5, 
    trtlghr = -0.5, bprogsl = 0.3, shape1 = 1.8, 
    scale1 = 360, shape2 = 1.7, scale2 = 688, 
    pmix = 0.5, admin = 5000, pcatnotrtbprog = 0.5, 
    pcattrtbprog = 0.25, pcatnotrt = 0.2, pcattrt = 0.1, 
    catmult = 0.5, tdxo = 1, ppoor = 0.1, pgood = 0.04, 
    ppoormet = 0.4, pgoodmet = 0.2, xomult = 1.4188308, 
    milestone = 546, outputRawDataset = 1, seed = 2000)
  
  fit1 <- ipcw(
    sim1$paneldata, id = "id", tstart = "tstart", 
    tstop = "tstop", event = "event", treat = "trtrand", 
    swtrt = "xo", swtrt_time = "xotime", base_cov = "bprog", 
    numerator = "bprog", denominator = "bprog*catlag", 
    logistic_switching_model = TRUE, ns_df = 3,
    swtrt_control_only = TRUE, boot = FALSE)
  
  # exclude observations after treatment switch
  data1 <- sim1$paneldata %>%
    filter(xo == 0 | tstart < xotime) %>%
    group_by(id) %>%
    mutate(condition = row_number() == n() & xo == 1 & tstop >= xotime,
           cross = ifelse(condition, 1, 0),
           event = ifelse(condition, 0, event),
           tstop = ifelse(condition, xotime, tstop))
  
  # fit pooled logistic regression switching models
  data2 <- data1 %>% filter(trtrand == 0)
  
  ns1 <- ns(data2$tstop[data2$cross == 1], df = 3)
  ns2 <- ns(data2$tstop, knots = attr(ns1, "knots"), 
            Boundary.knots = attr(ns1, "Boundary.knots"))
  switch1 <- glm(cross ~ bprog*catlag + ns2, family = binomial, data = data2)
  phat1 <- as.numeric(predict(switch1, newdata = data2, type = "response"))
  
  switch2 <- glm(cross ~ bprog + ns2, family = binomial, data = data2)
  phat2 <- as.numeric(predict(switch2, newdata = data2, type = "response"))
  
  # stablized weights
  data3 <- data1 %>% 
    filter(trtrand == 0) %>% 
    select(id, trtrand, bprog, tstart, tstop, event) %>%
    left_join(data2 %>% 
                ungroup() %>%
                mutate(pstay_den = 1-phat1, pstay_num = 1-phat2) %>%
                select(id, tstop, pstay_den, pstay_num), 
              by = c("id", "tstop")) %>%
    mutate(pstay_den = ifelse(is.na(pstay_den), 1, pstay_den),
           pstay_num = ifelse(is.na(pstay_num), 1, pstay_num)) %>%
    group_by(id) %>%
    mutate(surv_den = cumprod(pstay_den), surv_num = cumprod(pstay_num),
           stabilized_weight = surv_num/surv_den) %>%
    select(id, tstart, tstop, event, stabilized_weight, bprog, trtrand) %>%
    bind_rows(data1 %>% 
                filter(trtrand == 1) %>%
                mutate(stabilized_weight = 1) %>%
                select(id, tstart, tstop, event, 
                       stabilized_weight, bprog, trtrand))
  
  fit <- coxph(Surv(tstart, tstop, event) ~ trtrand + bprog, 
               data = data3, weight = stabilized_weight,
                id = id, ties = "efron", robust = TRUE)
  
  hr1 <- as.numeric(exp(cbind(fit$coefficients, confint(fit)))["trtrand",])
  
  testthat::expect_equal(data3$stabilized_weight, 
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
    filter(!co | tstart < dco) %>%
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
  data2$event[a1$censor] = 0;
  data2$cross[a1$censor] = 0;
  
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
