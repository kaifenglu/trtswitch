# Two-Stage Estimation with g-Estimation (TSEgest) for Treatment Switching

Estimates the causal parameter using g-estimation by fitting a pooled
logistic regression switching model that includes counterfactual
*unswitched* survival times and time-dependent confounders as
covariates. The adjusted hazard ratio is then obtained from the Cox
model using counterfactual *unswitched* survival times based on the
estimated causal parameter.

## Usage

``` r
tsegest(
  data,
  id = "id",
  stratum = "",
  tstart = "tstart",
  tstop = "tstop",
  event = "event",
  treat = "treat",
  censor_time = "censor_time",
  pd = "pd",
  pd_time = "pd_time",
  swtrt = "swtrt",
  swtrt_time = "swtrt_time",
  base_cov = "",
  conf_cov = "",
  strata_main_effect_only = TRUE,
  ns_df = 3,
  firth = FALSE,
  flic = FALSE,
  low_psi = -2,
  hi_psi = 2,
  n_eval_z = 101,
  recensor = TRUE,
  admin_recensor_only = TRUE,
  swtrt_control_only = TRUE,
  gridsearch = TRUE,
  root_finding = "brent",
  alpha = 0.05,
  ties = "efron",
  tol = 1e-06,
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

  - `id`: The id to identify observations belonging to the same subject
    for counting process data with time-dependent covariates.

  - `stratum`: The stratum.

  - `tstart`: The starting time of the time interval for
    counting-process data with time-dependent covariates.

  - `tstop`: The stopping time of the time interval for counting-process
    data with time-dependent covariates.

  - `event`: The event indicator, 1=event, 0=no event.

  - `treat`: The randomized treatment indicator, 1=treatment, 0=control.

  - `censor_time`: The administrative censoring time. It should be
    provided for all subjects including those who had events.

  - `pd`: The disease progression indicator, 1=PD, 0=no PD.

  - `pd_time`: The time from randomization to disease progression.

  - `swtrt`: The treatment switch indicator, 1=switch, 0=no switch.

  - `swtrt_time`: The time from randomization to treatment switch.

  - `base_cov`: The baseline covariates (excluding treat).

  - `conf_cov`: The confounding variables (excluding treat) for
    predicting treatment switching.

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
  for the Cox model.

- conf_cov:

  The names of confounding variables (excluding treat) in the input data
  for the logistic regression switching model.

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

- low_psi:

  The lower limit of the causal parameter.

- hi_psi:

  The upper limit of the causal parameter.

- n_eval_z:

  The number of points between `low_psi` and `hi_psi` (inclusive) at
  which to evaluate the Wald statistics for the coefficient of the
  counterfactual in the logistic regression switching model.

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

- gridsearch:

  Whether to use grid search to estimate the causal parameter `psi`.
  Defaults to `TRUE`, otherwise, a root finding algorithm will be used.

- root_finding:

  Character string specifying the univariate root-finding algorithm to
  use. Options are `"brent"` (default) for Brent's method, or
  `"bisection"` for the bisection method.

- alpha:

  The significance level to calculate confidence intervals.

- ties:

  The method for handling ties in the Cox model, either "breslow" or
  "efron" (default).

- tol:

  The desired accuracy (convergence tolerance) for `psi` for the root
  finding algorithm.

- offset:

  The offset to calculate the time from disease progression to death or
  censoring, the time from disease progression to treatment switch, and
  the time from treatment switch to death or censoring. We can set
  `offset` equal to 0 (no offset), and 1 (default), 1/30.4375, or
  1/365.25 if the time unit is day, month, or year, respectively.

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

- `psi_roots`: Vector of `psi` values for the control group at which the
  Z-statistic is zero, identified using grid search and linear
  interpolation.

- `psi_CI`: The confidence interval for `psi`.

- `psi_CI_type`: The type of confidence interval for `psi`, i.e., "grid
  search", "root finding", or "bootstrap".

- `logrank_pvalue`: The two-sided p-value of the log-rank test for the
  ITT analysis.

- `cox_pvalue`: The two-sided p-value for treatment effect based on the
  Cox model applied to counterfactual unswitched survival times. If
  `boot` is `TRUE`, this value represents the bootstrap p-value.

- `hr`: The estimated hazard ratio from the Cox model.

- `hr_CI`: The confidence interval for hazard ratio.

- `hr_CI_type`: The type of confidence interval for hazard ratio, either
  "Cox model" or "bootstrap".

- `event_summary`: A data frame containing the count and percentage of
  deaths, disease progressions, and switches by treatment arm.

- `data_switch`: The list of input data for the time from disease
  progression to switch by treatment group. The variables include `id`,
  `stratum`, `"swtrt"`, and `"swtrt_time"`. If `swtrt == 0`, then
  `swtrt_time` is censored at the time from disease progression to death
  or censoring.

- `km_switch`: The list of Kaplan-Meier plot data for the time from
  disease progression to switch by treatment group.

- `eval_z`: The list of data by treatment group containing the Wald
  statistics for the coefficient of the counterfactual in the logistic
  regression switching model, evaluated at a sequence of `psi` values.
  Used to plot and check if the range of `psi` values to search for the
  solution and limits of confidence interval of `psi` need be modified.

- `data_nullcox`: The list of input data for counterfactual survival
  times for the null Cox model by treatment group. The variables include
  `id`, `stratum`, `"t_star"` and `"d_star"`.

- `fit_nullcox`: The list of fitted null Cox models for counterfactual
  survival times by treatment group, which contains the martingale
  residuals.

- `data_logis`: The list of input data for pooled logistic regression
  models for treatment switching using g-estimation. The variables
  include `id`, `stratum`, `"tstart"`, `"tstop"`, `"cross"`,
  `"counterfactual"`, `conf_cov`, `ns`, `pd_time`, `swtrt`, and
  `swtrt_time`.

- `fit_logis`: The list of fitted pooled logistic regression models for
  treatment switching using g-estimation.

- `data_outcome`: The input data for the outcome Cox model of
  counterfactual unswitched survival times. The variables include `id`,
  `stratum`, `"t_star"`, `"d_star"`, `"treated"`, `base_cov` and
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

- `psi_trt_roots`: Vector of `psi_trt` values for the experimental group
  at which the Z-statistic is zero, identified using grid search and
  linear interpolation, if `swtrt_control_only` is `FALSE`.

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

- Fit a pooled logistic regression switching model among control-arm
  patients who experienced disease progression:
  \$\$\text{logit}(p(E\_{ik})) = \alpha U\_{i,\psi} + \sum\_{j} \beta_j
  x\_{ijk}\$\$ where \\E\_{ik}\\ is the switch indicator for subject
  \\i\\ at observation \\k\\, \$\$U\_{i,\psi} = T\_{C_i} +
  e^{\psi}T\_{E_i}\$\$ is the counterfactual survival time given a
  specific \\\psi\\, and \\x\_{ijk}\\ represents the time-dependent
  confounders. Natural cubic splines of time can be included to model
  time-varying baseline hazards. \\U\_{i,\psi}\\ is defined relative to
  the secondary baseline at disease progression and represents
  post-progression counterfactual survival, where \\T\_{C_i}\\ and
  \\T\_{E_i}\\ correspond to time spent after progression on control and
  experimental treatments, respectively. Martingale residuals may be
  used in place of counterfactual survival times to account for
  censoring.

- Identify the value of \\\psi\\ for which the Z-statistic of \\\alpha\\
  is approximately zero. This value is the causal parameter estimate.

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

If grid search is used to estimate \\\psi\\, the estimated \\\psi\\ is
the one with the smallest absolute value among those at which the
Z-statistic is zero based on linear interpolation. If root finding is
used, the estimated \\\psi\\ is the solution to the equation where the
Z-statistic is zero.

## References

NR Latimer, IR White, K Tilling, and U Siebert. Improved two-stage
estimation to adjust for treatment switching in randomised trials:
g-estimation to address time-dependent confounding. Statistical Methods
in Medical Research. 2020;29(10):2900-2918.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r

library(dplyr)

sim1 <- tsegestsim(
  n = 500, allocation1 = 2, allocation2 = 1, pbprog = 0.5, 
  trtlghr = -0.5, bprogsl = 0.3, shape1 = 1.8, 
  scale1 = 360, shape2 = 1.7, scale2 = 688, 
  pmix = 0.5, admin = 5000, pcatnotrtbprog = 0.5, 
  pcattrtbprog = 0.25, pcatnotrt = 0.2, pcattrt = 0.1, 
  catmult = 0.5, tdxo = 1, ppoor = 0.1, pgood = 0.04, 
  ppoormet = 0.4, pgoodmet = 0.2, xomult = 1.4188308, 
  milestone = 546, seed = 2000)
  
data1 <- sim1$paneldata %>%
  mutate(visit7on = ifelse(progressed, tstop > timePFSobs + 105, 0))

fit1 <- tsegest(
  data = data1, id = "id", 
  tstart = "tstart", tstop = "tstop", event = "event", 
  treat = "trtrand", censor_time = "censor_time", 
  pd = "progressed", pd_time = "timePFSobs", 
  swtrt = "xo", swtrt_time = "xotime", 
  base_cov = "bprog", 
  conf_cov = c("bprog*cattdc", "timePFSobs", "visit7on"), 
  ns_df = 3, low_psi = -1, hi_psi = 1, n_eval_z = 101,
  recensor = TRUE, admin_recensor_only = TRUE, 
  swtrt_control_only = TRUE, alpha = 0.05, ties = "efron", 
  tol = 1.0e-6, offset = 0, boot = FALSE)
  
fit1
#>             n event_n event_pct pd_n pd_pct switch_n switch_pct event_out_n
#> Control   167     167     100.0  165   98.8       86       51.5         167
#> Treatment 333     333     100.0  333  100.0        0        0.0         333
#>           event_out_pct trtrand
#> Control           100.0       0
#> Treatment         100.0       1
#> 
#>                            Estimate Lower 95% Upper 95%
#> Causal parameter psi       -0.414   -0.654    -0.188   
#> Causal survival time ratio 1.513    1.207     1.923    
#> Hazard ratio (HR)          0.548    0.452     0.664    
#> P-value (Cox model)        <.0001                      
```
