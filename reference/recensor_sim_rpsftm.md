# Simulation Study to Evaluate Recensoring Rules in RPSFTM

Simulates datasets to evaluate the performance of various recensoring
strategies under the Rank Preserving Structural Failure Time Model
(RPSFTM) for handling treatment switching in survival analysis.

## Usage

``` r
recensor_sim_rpsftm(
  nsim = 100L,
  n = 400L,
  shape = 1.5,
  scale = 553.9,
  gamma = 0.001,
  tfmin = 407.5,
  tfmax = 407.5,
  psi = -0.4621,
  omega = 0,
  pswitch = 0.7,
  a = 2,
  b = 4,
  low_psi = -5,
  hi_psi = 5,
  treat_modifier = 1,
  recensor_type = 1L,
  admin_recensor_only = TRUE,
  autoswitch = TRUE,
  alpha = 0.05,
  ties = "efron",
  tol = 1e-06,
  boot = TRUE,
  n_boot = 100L,
  seed = 0L
)
```

## Arguments

- nsim:

  Number of simulated datasets.

- n:

  Number of subjects per simulation.

- shape:

  Shape parameter of the Weibull distribution for time to death.

- scale:

  Scale parameter of the Weibull distribution for time to death in the
  control group.

- gamma:

  Rate parameter of the exponential distribution for random dropouts in
  the control group.

- tfmin:

  Minimum planned follow-up time (in days).

- tfmax:

  Maximum planned follow-up time (in days).

- psi:

  Log time ratio of death time for control vs experimental treatment.

- omega:

  Log time ratio of dropout time for control vs experimental treatment.

- pswitch:

  Probability of treatment switching at disease progression.

- a:

  Shape parameter 1 of the Beta distribution for time to disease
  progression as a fraction of time to death.

- b:

  Shape parameter 2 of the Beta distribution for time to disease
  progression.

- low_psi:

  Lower bound for the search interval of the causal parameter \\\psi\\.

- hi_psi:

  Upper bound for the search interval of the causal parameter \\\psi\\.

- treat_modifier:

  Sensitivity parameter modifying the constant treatment effect
  assumption.

- recensor_type:

  Type of recensoring to apply:

  - 0: No recensoring

  - 1: Recensor all control-arm subjects

  - 2: Recensor only switchers in the control arm

  - 3: Recensor only control-arm switchers whose counterfactual survival
    exceeds the planned follow-up time

- admin_recensor_only:

  Logical. If `TRUE`, recensoring is applied only to administrative
  censoring times. If `FALSE`, it is also applied to dropout times.

- autoswitch:

  Logical. If `TRUE`, disables recensoring in arms without any treatment
  switching.

- alpha:

  Significance level for confidence interval calculation (default is
  0.05).

- ties:

  Method for handling tied event times in the Cox model. Options are
  `"efron"` (default) or `"breslow"`.

- tol:

  Convergence tolerance for root-finding in estimation of \\\psi\\.

- boot:

  Logical. If `TRUE`, bootstrap is used to estimate the confidence
  interval for the hazard ratio. If `FALSE`, the confidence interval is
  matched to the log-rank p-value.

- n_boot:

  Number of bootstrap samples, used only if `boot = TRUE`.

- seed:

  Optional. Random seed for reproducibility.

## Value

A data frame summarizing the simulation results, including:

- `recensor_type`, `admin_recensor_only`: Settings used in the
  simulation.

- Event rates: `p_event_1`, `p_dropout_1`, `p_admin_censor_1`,
  `p_event_0`, `p_dropout_0`, `p_admin_censor_0`.

- Progression and switching: `p_pd_0`, `p_swtrt_0`, `p_recensored_0`.

- Causal parameter (\\\psi\\) estimates: `psi`, `psi_est`, `psi_bias`,
  `psi_se`, `psi_mse`.

- Log hazard ratio estimates: `loghr`, `loghr_est`, `loghr_se`,
  `loghr_mse`.

- Hazard ratio metrics: `hr`, `hr_est` (geometric mean), `hr_pctbias`
  (percent bias).

- Standard errors of log hazard ratio: `loghr_se_cox`, `loghr_se_lr`,
  `loghr_se_boot`.

- Coverage probabilities: `hr_ci_cover_cox`, `hr_ci_cover_lr`,
  `hr_ci_cover_boot`.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
# \donttest{
result <- recensor_sim_rpsftm(
  nsim = 10, n = 400, shape = 1.5, scale = exp(6.3169),
  gamma = 0.001, tfmin = 407.5, tfmax = 407.5,
  psi = log(0.5) / 1.5, omega = log(1), pswitch = 0.7,
  a = 2, b = 4, low_psi = -5, hi_psi = 5,
  treat_modifier = 1, recensor_type = 1,
  admin_recensor_only = TRUE, autoswitch = TRUE,
  alpha = 0.05, tol = 1e-6, boot = TRUE,
  n_boot = 10, seed = 314159)
# }
```
