# Simple Two-Stage Estimation

``` r

library(trtswitch)
library(dplyr, warn.conflicts = FALSE)
```

## Introduction

The two-stage estimation (TSE) approach first estimates the effect of
treatment switching on survival, then uses this estimate to derive
counterfactual survival times if switching had not occurred. A key
requirement of the TSE approach is that treatment switching must occur
at or after a disease-related “secondary baseline” time point, typically
defined by disease progression.

The TSE method assumes the no unmeasured confounding condition, meaning
treatment switching must be independent of potential outcomes, provided
all relevant patient characteristics at the secondary baseline are
measured and included in the model. The simple TSE method further
assumes that no time-dependent confounding exists between the secondary
baseline and the time of switching.

## Estimation of \\\psi\\

The simple TSE method involves applying an accelerated failure time
(AFT) model that compares post-progression survival in control group
switchers with post-progression survival in control group non-switchers.
Prognostic variables measured at the secondary baseline are included to
account for differences between the groups. The treatment effect of
switching is estimated as a time ratio and used to adjust the survival
times of switchers.

## Estimation of Hazard Ratio

Once \\\psi\\ has been estimated, we can derive an adjusted data set and
fit a (potentially stratified) Cox proportional hazards model to the
adjusted data set to obtain an estimate of the hazard ratio. The
confidence interval for the hazard ratio can be derived by bootstrapping
the entire adjustment and subsequent model-fitting process.

## Example

We start by preparing the data.

``` r

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
```

Next we apply the simple TSE method.

``` r

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
```

We can examine the Weibull AFT model fits and the corresponding value of
\\\hat{\psi}\\.

``` r

# control group
fit1$fit_aft[[1]]$fit$parest[, c("param", "beta", "sebeta", "z")]
#>                     param         beta      sebeta          z
#> 1             (Intercept)  5.572896483 0.758612068  7.3461743
#> 2                   swtrt  1.067653070 0.237997891  4.4859770
#> 3                 agerand  0.001927251 0.007493015  0.2572064
#> 4             sex.fFemale  0.316148630 0.240086216  1.3168129
#> 5                 tt_Lnum -0.023248007 0.040491209 -0.5741495
#> 6              rmh_alea.c -0.894680233 0.193398728 -4.6260916
#> 7             pathway.fHR -0.484778956 0.358671753 -1.3515950
#> 8  pathway.fPI3K.AKT.mTOR -0.604986176 0.382244440 -1.5827207
#> 9                      ps -0.548628111 0.136464658 -4.0202945
#> 10                    ttc  0.050010291 0.269239313  0.1857466
#> 11                   tran  0.437778730 0.385005847  1.1370703
#> 12             Log(scale) -0.428099622 0.105960290 -4.0401892

fit1$psi
#> [1] -1.067653

# experimental group
fit1$fit_aft[[2]]$fit$parest[, c("param", "beta", "sebeta", "z")]
#>                     param         beta      sebeta          z
#> 1             (Intercept)  4.819297083 0.760270308  6.3389258
#> 2                   swtrt  0.983747428 0.271723182  3.6204030
#> 3                 agerand  0.001253863 0.009513484  0.1317985
#> 4             sex.fFemale  0.395790503 0.237599761  1.6657866
#> 5                 tt_Lnum  0.066813220 0.051717213  1.2918952
#> 6              rmh_alea.c -0.714216404 0.209567497 -3.4080495
#> 7             pathway.fHR -0.301948752 0.296989948 -1.0166969
#> 8  pathway.fPI3K.AKT.mTOR -0.076210121 0.281075448 -0.2711376
#> 9                      ps -0.401333183 0.114040814 -3.5192066
#> 10                    ttc  0.414584864 0.377881586  1.0971290
#> 11                   tran -0.615854802 0.390436760 -1.5773484
#> 12             Log(scale) -0.414571271 0.116718016 -3.5519047

fit1$psi_trt
#> [1] -0.9837474
```

Now we fit the outcome Cox model and compare the treatment hazard ratio
estimate with the reported.

``` r

fit1$fit_outcome$parest[, c("param", "beta", "sebeta", "z")]
#>                    param          beta      sebeta           z
#> 1                treated -0.0914777971 0.184469500 -0.49589660
#> 2                agerand  0.0045435260 0.008026488  0.56606649
#> 3            sex.fFemale -0.4181004068 0.191327703 -2.18525806
#> 4                tt_Lnum -0.0008255737 0.036590479 -0.02256253
#> 5             rmh_alea.c  0.6881577040 0.191026441  3.60242121
#> 6            pathway.fHR  0.1580245496 0.297687989  0.53083952
#> 7 pathway.fPI3K.AKT.mTOR  0.1873258331 0.303906799  0.61639237

c(fit1$hr, fit1$hr_CI)
#> [1] 0.9125816 0.6356982 1.3100637
```

Finally, to ensure the uncertainty is accurately represented, the entire
adjustment process and subsequent survival modeling can be bootstrapped.

``` r

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
  boot = TRUE, n_boot = 1000, seed = 12345)

c(fit2$hr, fit2$hr_CI)
#> [1] 0.9125816 0.5570529 1.4950198
```
