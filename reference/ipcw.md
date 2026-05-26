# Inverse Probability of Censoring Weights (IPCW) for Treatment Switching

Excludes data after treatment switching and fits a switching model to
estimate the probability of not switching. The inverse of these
probabilities (inverse probability of censoring weights) are then used
as weights in a weighted Cox model to obtain the adjusted hazard ratio.

## Usage

``` r
ipcw(
  data,
  id = "id",
  stratum = "",
  tstart = "tstart",
  tstop = "tstop",
  event = "event",
  treat = "treat",
  swtrt = "swtrt",
  swtrt_time = "swtrt_time",
  base_cov = "",
  numerator = "",
  denominator = "",
  logistic_switching_model = FALSE,
  strata_main_effect_only = TRUE,
  ns_df = 3,
  firth = FALSE,
  flic = FALSE,
  stabilized_weights = TRUE,
  trunc = 0,
  trunc_upper_only = TRUE,
  swtrt_control_only = TRUE,
  alpha = 0.05,
  ties = "efron",
  boot = FALSE,
  n_boot = 1000,
  seed = 0,
  nthreads = 0
)
```

## Arguments

- data:

  The input data frame that contains the following variables:

  - `id`: The id to identify observations belonging to the same subject
    for counting process data with time-dependent covariates.

  - `stratum`: The stratum.

  - `tstart`: The starting time of the time interval for
    counting-process data with time-dependent covariates.

  - `tstop`: The stopping time of the time interval for counting-process
    data with time-dependent covariates.

  - `event`: The event indicator, 1=event, 0=no event.

  - `treat`: The randomized treatment indicator, 1=treatment, 0=control.

  - `swtrt`: The treatment switch indicator, 1=switch, 0=no switch.

  - `swtrt_time`: The time from randomization to treatment switch.

  - `base_cov`: The baseline covariates (excluding treat) used in the
    outcome model.

  - `numerator`: The baseline covariates (excluding treat) used in the
    numerator switching model for stabilized weights.

  - `denominator`: The baseline (excluding treat) and time-dependent
    covariates used in the denominator switching model.

- id:

  The name of the id variable in the input data.

- stratum:

  The name(s) of the stratum variable(s) in the input data.

- tstart:

  The name of the tstart variable in the input data.

- tstop:

  The name of the tstop variable in the input data.

- event:

  The name of the event variable in the input data.

- treat:

  The name of the treatment variable in the input data.

- swtrt:

  The name of the swtrt variable in the input data.

- swtrt_time:

  The name of the swtrt_time variable in the input data.

- base_cov:

  The names of baseline covariates (excluding treat) in the input data
  for the Cox model.

- numerator:

  The names of baseline covariates (excluding treat) in the input data
  for the numerator switching model for stabilized weights.

- denominator:

  The names of baseline (excluding treat) and time-dependent covariates
  in the input data for the denominator switching model.

- logistic_switching_model:

  Whether a pooled logistic regression switching model is used.

- strata_main_effect_only:

  Whether to only include the strata main effects in the logistic
  regression switching model. Defaults to `TRUE`, otherwise all possible
  strata combinations will be considered in the switching model.

- ns_df:

  Degrees of freedom for the natural cubic spline for visit-specific
  intercepts of the pooled logistic regression model. Defaults to 3 for
  two internal knots at the 33 and 67 percentiles of the treatment
  switching times.

- firth:

  Whether the Firth's bias reducing penalized likelihood should be used.

- flic:

  Whether to apply intercept correction to obtain more accurate
  predicted probabilities.

- stabilized_weights:

  Whether to use the stabilized weights. The default is `TRUE`.

- trunc:

  The truncation fraction of the weight distribution. Defaults to 0 for
  no truncation in weights.

- trunc_upper_only:

  Whether to truncate the weights from the upper end of the weight
  distribution only. Defaults to `TRUE`, otherwise the weights will be
  truncated from both the lower and upper ends of the distribution.

- swtrt_control_only:

  Whether treatment switching occurred only in the control group. The
  default is `TRUE`.

- alpha:

  The significance level to calculate confidence intervals.

- ties:

  The method for handling ties in the Cox model, either "breslow" or
  "efron" (default).

- boot:

  Whether to use bootstrap to obtain the confidence interval for hazard
  ratio. Defaults to `FALSE`.

- n_boot:

  The number of bootstrap samples.

- seed:

  The seed to reproduce the bootstrap results.

- nthreads:

  The number of threads to use in bootstrapping (0 means the default
  RcppParallel behavior).

## Value

A list with the following components:

- `pvalue`: The two-sided p-value.

- `pvalue_type`: The type of two-sided p-value for treatment effect,
  i.e., "Cox model" or "bootstrap".

- `hr`: The estimated hazard ratio from the Cox model.

- `hr_CI`: The confidence interval for hazard ratio.

- `hr_CI_type`: The type of confidence interval for hazard ratio, either
  "Cox model" or "bootstrap".

- `event_summary`: A data frame containing the count and percentage of
  deaths and switches by treatment arm.

- `data_switch`: A list of input data for the switching models by
  treatment group. The variables include `id`, `stratum`, `"tstart"`,
  `"tstop"`, `"cross"`, `denominator`, `swtrt`, and `swtrt_time`. For
  logistic switching models, `stratum` variables are converted to dummy
  variables, and natural cubic spline basis variables are created for
  the visit-specific intercepts.

- `fit_switch`: A list of fitted switching models for the denominator
  and numerator by treatment group.

- `data_outcome`: The input data for the outcome Cox model including the
  inverse probability of censoring weights. The variables include `id`,
  `stratum`, `"tstart"`, `"tstop"`, `"event"`, `"treated"`,
  `"unstablized_weight"`, `"stabilized_weight"`, `base_cov`, and
  `treat`.

- `weight_summary`: A data frame summarizing the weights by treatment
  arm.

- `km_outcome`: The Kaplan-Meier estimates of the survival functions for
  the treatment and control groups based on the weighted outcome data.

- `lr_outcome`: The log-rank test results for the treatment effect based
  on the weighted outcome data.

- `fit_outcome`: The fitted outcome Cox model.

- `fail`: Whether a model fails to converge.

- `settings`: A list containing the input parameter values.

- `fail_boots`: The indicators for failed bootstrap samples if `boot` is
  `TRUE`.

- `fail_boots_data`: The data for failed bootstrap samples if `boot` is
  `TRUE`.

- `hr_boots`: The bootstrap hazard ratio estimates if `boot` is `TRUE`.

## Details

The hazard ratio and confidence interval under a no-switching scenario
are obtained as follows:

- Exclude all observations after treatment switch.

- Define the crossover and event indicators for the last time interval
  of each subject.

- For time-dependent Cox switching models, replicate unique event times
  across treatment arms within each subject.

- Fit the denominator switching model (and numerator model for
  stabilized weights) to estimate inverse probability of censoring
  weights. Either a Cox model with time-dependent covariates or a pooled
  logistic regression model can be used.

  - For the pooled logistic regression model, the probability of
    remaining uncensored (i.e., not switching) is calculated as \\1 -
    \hat{p}\_{\text{switch}}\\ and accumulated over time up to the start
    of each interval.

  - For the time-dependent Cox model, the probability of remaining
    unswitched is derived from the estimated baseline hazard and
    predicted risk score up to the end of each interval.

- Fit a weighted Cox model to the outcome survival times (excluding data
  after switching) to estimate the hazard ratio.

- Construct the p-value and confidence interval for the hazard ratio
  using either robust sandwich variance or bootstrapping. When
  bootstrapping is used, the confidence interval and p-value are based
  on a t-distribution with `n_boot - 1` degrees of freedom.

## References

James M. Robins and Dianne M. Finkelstein. Correcting for noncompliance
and dependent censoring in an AIDS clinical trial with inverse
probability of censoring weighted (IPCW) log-rank tests. Biometrics.
2000;56(3):779-788.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
# Example 1: pooled logistic regression switching model

sim1 <- tssim(
  tdxo = TRUE, coxo = TRUE, allocation1 = 1, allocation2 = 1,
  p_X_1 = 0.3, p_X_0 = 0.3, 
  rate_T = 0.002, beta1 = -0.5, beta2 = 0.3, 
  gamma0 = 0.3, gamma1 = -0.9, gamma2 = 0.7, gamma3 = 1.1, gamma4 = -0.8,
  zeta0 = -3.5, zeta1 = 0.5, zeta2 = 0.2, zeta3 = -0.4, 
  alpha0 = 0.5, alpha1 = 0.5, alpha2 = 0.4, 
  theta1_1 = -0.4, theta1_0 = -0.4, theta2 = 0.2,
  rate_C = 0.0000855, accrualIntensity = 20/30,
  fixedFollowup = FALSE, plannedTime = 1350, days = 30,
  n = 500, NSim = 100, seed = 314159)
  
fit1 <- ipcw(
  sim1[[1]], id = "id", tstart = "tstart", 
  tstop = "tstop", event = "event", treat = "trtrand", 
  swtrt = "xo", swtrt_time = "xotime", 
  base_cov = "bprog", numerator = "bprog", 
  denominator = c("bprog", "L"),
  logistic_switching_model = TRUE, ns_df = 3,
  swtrt_control_only = TRUE, boot = FALSE)
  
fit1 
#>             n event_n event_pct switch_n switch_pct event_out_n event_out_pct
#> Control   250     208      83.2       77       30.8         161          64.4
#> Treatment 250     185      74.0        0        0.0         185          74.0
#>           trtrand
#> Control         0
#> Treatment       1
#> 
#>                                 Weight summary
#>            N    Min     Q1 Median   Mean     Q3    Max trtrand
#> Control 2093 0.9358 0.9959 1.0005 1.0005 1.0069 1.0499       0
#> 
#>                     Estimate Lower 95% Upper 95%
#> Hazard ratio (HR)   0.553    0.444     0.689    
#> P-value (Cox model) <.0001                      

# Example 2: time-dependent covariates Cox switching model

fit2 <- ipcw(
  shilong, id = "id", tstart = "tstart", tstop = "tstop", 
  event = "event", treat = "bras.f", swtrt = "co", 
  swtrt_time = "dco", 
  base_cov = c("agerand", "sex.f", "tt_Lnum", "rmh_alea.c", 
               "pathway.f"),
  numerator = c("agerand", "sex.f", "tt_Lnum", "rmh_alea.c", 
                "pathway.f"),
  denominator = c("agerand", "sex.f", "tt_Lnum", "rmh_alea.c",
                  "pathway.f", "ps", "ttc", "tran"),
  swtrt_control_only = FALSE, boot = FALSE)

fit2
#>             n event_n event_pct switch_n switch_pct event_out_n event_out_pct
#> Control    93      63      67.7       68       73.1          23          24.7
#> Treatment 100      67      67.0       25       25.0          53          53.0
#>           bras.f
#> Control       CT
#> Treatment    MTA
#> 
#>                                   Weight summary
#>              N    Min     Q1 Median   Mean     Q3    Max bras.f
#> Control   3213 0.7223 0.9654 0.9944 0.9971 1.0079 1.8117     CT
#> Treatment 4301 0.6950 0.9917 0.9991 0.9988 1.0004 1.9053    MTA
#> 
#>                     Estimate Lower 95% Upper 95%
#> Hazard ratio (HR)   1.428    0.866     2.355    
#> P-value (Cox model) 0.1627                      
```
