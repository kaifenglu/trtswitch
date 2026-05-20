# Marginal Structural Model (MSM) for Treatment Switching

Excludes data after treatment switching when fitting the switching model
to estimate the probabilities of not switching and then switching. The
inverse of these probabilities (inverse probability of treatment
weights) are then used as weights in a Cox model including data after
switching to estimate the adjusted hazard ratio.

## Usage

``` r
msm(
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
  strata_main_effect_only = TRUE,
  ns_df = 3,
  firth = FALSE,
  flic = FALSE,
  stabilized_weights = TRUE,
  trunc = 0,
  trunc_upper_only = TRUE,
  swtrt_control_only = TRUE,
  treat_alt_interaction = TRUE,
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

- treat_alt_interaction:

  Whether to include an interaction between randomized and alternative
  treatments in the outcome model when both randomized arms can switch
  to alternative treatment.

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
  RcppParallel behavior)

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
  `"tstop"`, `"cross"`, `denominator`, `swtrt`, and `swtrt_time`. In
  addition, `stratum` variables are converted to dummy variables, and
  natural cubic spline basis variables are created for the
  visit-specific intercepts.

- `fit_switch`: A list of fitted switching models for the denominator
  and numerator by treatment group.

- `data_outcome`: The input data for the outcome Cox model including the
  inverse probability of censoring weights. The variables include `id`,
  `stratum`, `"tstart"`, `"tstop"`, `"event"`, `"treated"`, `"crossed"`,
  `"unstablized_weight"`, `"stabilized_weight"`, `base_cov`, and
  `treat`. If `treat_alt_interaction` is `TRUE`, the data set also
  includes the `"treated_crossed"` variable.

- `weight_summary`: A data frame summarizing the weights by treatment
  arm.

- `km_outcome`: The Kaplan-Meier estimates of the survival functions for
  the treatment and control groups based on the weighted outcome data
  truncated at time of treatment switching.

- `lr_outcome`: The log-rank test results for the treatment effect based
  on the weighted outcome data truncated at time of treatment switching.

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

- Exclude observations after treatment switch when fitting the switching
  model.

- Define crossover indicators for the last time interval of each
  subject.

- Fit the denominator switching model (and numerator model for
  stabilized weights) using a pooled logistic regression model to
  estimate the inverse probability of treatment weights (IPTWs).

  - The probability of remaining unswitched is calculated as \\1 -
    \hat{p}\_{\text{switch}}\\ and multiplied over time before treatment
    switch.

  - At the time of switching, this product is multiplied by the
    predicted probability of switching.

  - After treatment switch, the IPTW remains constant.

  - The inverse of the probability at the start of each interval is used
    as the interval weight.

- Fit a weighted Cox model to the outcome survival times, including data
  after treatment switch, to estimate the hazard ratio.

- Construct the p-value and confidence interval for the hazard ratio
  using either robust sandwich variance or bootstrapping. When
  bootstrapping is used, the confidence interval and p-value are based
  on a t-distribution with `n_boot - 1` degrees of freedom.

## References

James M. Robins, Miguel Angel Hernan, and Babette Brumback. Marginal
structural models and causal inference in epidemiology. Epidemiology.
2000;11(5):550-560.

Miguel Angel Hernan, Babette Brumback, and James M. Robins. Marginal
structural modesl to estimate the causual effect of zidovudine on the
survival of HIV-positive men. Epidemiology. 2000;11(5):561-570.

Jing Xu, Guohui Liu, and Bingxia Wang. Bias and Type I error control in
correcting treatment effect for treatment switching using marginal
structural models in Phase III oncology trials. Journal of
Biopharmaceutical Statistics. 2022;32(6):897-914.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
sim1 <- tssim(
  tdxo = 1, coxo = 1, allocation1 = 1, allocation2 = 1,
  p_X_1 = 0.3, p_X_0 = 0.3, 
  rate_T = 0.002, beta1 = -0.5, beta2 = 0.3, 
  gamma0 = 0.3, gamma1 = -0.9, gamma2 = 0.7, gamma3 = 1.1, gamma4 = -0.8,
  zeta0 = -3.5, zeta1 = 0.5, zeta2 = 0.2, zeta3 = -0.4, 
  alpha0 = 0.5, alpha1 = 0.5, alpha2 = 0.4, 
  theta1_1 = -0.4, theta1_0 = -0.4, theta2 = 0.2,
  rate_C = 0.0000855, accrualIntensity = 20/30,
  fixedFollowup = FALSE, plannedTime = 1350, days = 30,
  n = 500, NSim = 100, seed = 314159)

fit1 <- msm(
  sim1[[1]], id = "id", tstart = "tstart", 
  tstop = "tstop", event = "event", treat = "trtrand", 
  swtrt = "xo", swtrt_time = "xotime", 
  base_cov = "bprog", numerator = "bprog", 
  denominator = c("bprog", "L"), 
  ns_df = 3, swtrt_control_only = TRUE, boot = FALSE)
  
fit1
#>             n event_n event_pct switch_n switch_pct event_out_n event_out_pct
#> Control   250     208      83.2       77       30.8         208          83.2
#> Treatment 250     185      74.0        0        0.0         185          74.0
#>           trtrand
#> Control         0
#> Treatment       1
#> 
#>                                 Weight summary
#>            N    Min     Q1 Median   Mean     Q3    Max trtrand
#> Control 3148 0.8813 0.9686 0.9976 0.9958 1.0046 1.3166       0
#> 
#>                     Estimate Lower 95% Upper 95%
#> Hazard ratio (HR)   0.555    0.447     0.690    
#> P-value (Cox model) <.0001                      
```
