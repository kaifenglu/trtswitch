# Simple Two-Stage Estimation (TSEsimp) for Treatment Switching

Estimates the causal parameter by fitting an accelerated failure time
(AFT) model comparing post-progression survival between switchers and
non-switchers, and derives the adjusted hazard ratio from the Cox model
using counterfactual *unswitched* survival times based on the estimated
causal parameter.

## Usage

``` r
tsesimp(
  data,
  id = "id",
  stratum = "",
  time = "time",
  event = "event",
  treat = "treat",
  censor_time = "censor_time",
  pd = "pd",
  pd_time = "pd_time",
  swtrt = "swtrt",
  swtrt_time = "swtrt_time",
  base_cov = "",
  base2_cov = "",
  aft_dist = "weibull",
  strata_main_effect_only = TRUE,
  recensor = TRUE,
  admin_recensor_only = TRUE,
  swtrt_control_only = TRUE,
  alpha = 0.05,
  ties = "efron",
  offset = 1,
  boot = TRUE,
  n_boot = 1000,
  seed = 0,
  nthreads = 0
)
```

## Arguments

- data:

  The input data frame that contains the following variables:

  - `id`: The subject id.

  - `stratum`: The stratum.

  - `time`: The survival time for right censored data.

  - `event`: The event indicator, 1=event, 0=no event.

  - `treat`: The randomized treatment indicator, 1=treatment, 0=control.

  - `censor_time`: The administrative censoring time. It should be
    provided for all subjects including those who had events.

  - `pd`: The disease progression indicator, 1=PD, 0=no PD.

  - `pd_time`: The time from randomization to disease progression.

  - `swtrt`: The treatment switch indicator, 1=switch, 0=no switch.

  - `swtrt_time`: The time from randomization to treatment switch.

  - `base_cov`: The baseline covariates (excluding treat).

  - `base2_cov`: The baseline and secondary baseline covariates
    (excluding treat).

- id:

  The name of the id variable in the input data.

- stratum:

  The name(s) of the stratum variable(s) in the input data.

- time:

  The name of the time variable in the input data.

- event:

  The name of the event variable in the input data.

- treat:

  The name of the treatment variable in the input data.

- censor_time:

  The name of the censor_time variable in the input data.

- pd:

  The name of the pd variable in the input data.

- pd_time:

  The name of the pd_time variable in the input data.

- swtrt:

  The name of the swtrt variable in the input data.

- swtrt_time:

  The name of the swtrt_time variable in the input data.

- base_cov:

  The names of baseline covariates (excluding treat) in the input data
  for the outcome Cox model.

- base2_cov:

  The names of baseline and secondary baseline covariates (excluding
  treat) in the input data for the AFT model for post-progression
  survival.

- aft_dist:

  The assumed distribution for time to event for the AFT model. Options
  include "exponential", "weibull" (default), "loglogistic", and
  "lognormal".

- strata_main_effect_only:

  Whether to only include the strata main effects in the AFT model.
  Defaults to `TRUE`, otherwise all possible strata combinations will be
  considered in the AFT model.

- recensor:

  Whether to apply recensoring to counterfactual survival times.
  Defaults to `TRUE`.

- admin_recensor_only:

  Whether to apply recensoring to administrative censoring times only.
  Defaults to `TRUE`. If `FALSE`, recensoring will be applied to the
  actual censoring times for dropouts.

- swtrt_control_only:

  Whether treatment switching occurred only in the control group. The
  default is `TRUE`.

- alpha:

  The significance level to calculate confidence intervals.

- ties:

  The method for handling ties in the Cox model, either "breslow" or
  "efron" (default).

- offset:

  The offset to calculate the time disease progression to death or
  censoring. We can set `offset` equal to 0 (no offset), and 1
  (default), 1/30.4375, or 1/365.25 if the time unit is day, month, or
  year, respectively.

- boot:

  Whether to use bootstrap to obtain the confidence interval for hazard
  ratio. Defaults to `TRUE`.

- n_boot:

  The number of bootstrap samples.

- seed:

  The seed to reproduce the bootstrap results.

- nthreads:

  The number of threads to use in bootstrapping (0 means the default
  RcppParallel behavior)

## Value

A list with the following components:

- `psi`: The estimated causal parameter for the control group.

- `psi_CI`: The confidence interval for `psi`.

- `psi_CI_type`: The type of confidence interval for `psi`, i.e., "AFT
  model" or "bootstrap".

- `pvalue`: The two-sided p-value.

- `pvalue_type`: The type of two-sided p-value for treatment effect,
  i.e., "Cox model" or "bootstrap".

- `hr`: The estimated hazard ratio from the Cox model.

- `hr_CI`: The confidence interval for hazard ratio.

- `hr_CI_type`: The type of confidence interval for hazard ratio, either
  "Cox model" or "bootstrap".

- `event_summary`: A data frame containing the count and percentage of
  deaths, disease progressions, and switches by treatment arm.

- `data_aft`: A list of input data for the AFT model by treatment group.
  The variables include `id`, `stratum`, `"pps"`, `"event"`, `"swtrt"`,
  `base2_cov`, `pd_time`, `swtrt_time`, and `time`.

- `fit_aft`: A list of fitted AFT models by treatment group.

- `res_aft`: A list of deviance residuals from the AFT models by
  treatment group.

- `data_outcome`: The input data for the outcome Cox model of
  counterfactual unswitched survival times. The variables include `id`,
  `stratum`, `"t_star"`, `"d_star"`, `"treated"`, `base_cov`, and
  `treat`.

- `km_outcome`: The Kaplan-Meier estimates of the survival functions for
  the treatment and control groups based on the counterfactual
  unswitched survival times.

- `lr_outcome`: The log-rank test results for the treatment effect based
  on the counterfactual unswitched survival times.

- `fit_outcome`: The fitted outcome Cox model.

- `fail`: Whether a model fails to converge.

- `psimissing`: Whether the `psi` parameter cannot be estimated.

- `settings`: A list containing the input parameter values.

- `psi_trt`: The estimated causal parameter for the experimental group
  if `swtrt_control_only` is `FALSE`.

- `psi_trt_CI`: The confidence interval for `psi_trt` if
  `swtrt_control_only` is `FALSE`.

- `fail_boots`: The indicators for failed bootstrap samples if `boot` is
  `TRUE`.

- `fail_boots_data`: The data for failed bootstrap samples if `boot` is
  `TRUE`.

- `hr_boots`: The bootstrap hazard ratio estimates if `boot` is `TRUE`.

- `psi_boots`: The bootstrap `psi` estimates if `boot` is `TRUE`.

- `psi_trt_boots`: The bootstrap `psi_trt` estimates if `boot` is `TRUE`
  and `swtrt_control_only` is `FALSE`.

## Details

Assuming one-way switching from control to treatment, the hazard ratio
and confidence interval under a no-switching scenario are obtained as
follows:

- Estimate the causal parameter \\\psi\\ by fitting an AFT model
  comparing post-progression survival between switchers and
  non-switchers in the control group who experienced disease
  progression.

- Compute counterfactual survival times for control patients using the
  estimated \\\psi\\.

- Fit a Cox model to the observed survival times for the treatment group
  and the counterfactual survival times for the control group to
  estimate the hazard ratio.

- When bootstrapping is used, derive the confidence interval and p-value
  for the hazard ratio from a t-distribution with `n_boot - 1` degrees
  of freedom.

If treatment switching occurs before or in the absence of recorded
disease progression, the patient is considered to have progressed at the
time of treatment switching.

## References

Nicholas R Latimer, KR Abrams, PC Lambert, MK Crowther, AJ Wailoo, JP
Morden, RL Akehurst, and MJ Campbell. Adjusting for treatment switching
in randomised controlled trials - A simulation study and a simplified
two-stage method. Statistical Methods in Medical Research.
2017;26(2):724-751.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r

library(dplyr)

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

fit1
#>             n event_n event_pct pd_n pd_pct switch_n switch_pct event_out_n
#> Control    93      63      67.7   85   91.4       68       73.1          58
#> Treatment 100      67      67.0   83   83.0       25       25.0          63
#>           event_out_pct bras.f
#> Control            62.4     CT
#> Treatment          63.0    MTA
#> 
#>                                              Estimate Lower 95% Upper 95%
#> Causal parameter psi for control arm         -1.068   -1.534    -0.601   
#> Causal survival time ratio for control arm   2.909    1.824     4.637    
#> Causal parameter psi for treatment arm       -0.984   -1.516    -0.451   
#> Causal survival time ratio for treatment arm 2.674    1.570     4.555    
#> Hazard ratio (HR)                            0.913    0.636     1.310    
#> P-value (Cox model)                          0.6200                      
```
