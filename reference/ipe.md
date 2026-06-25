# Iterative Parameter Estimation (IPE) for Treatment Switching

Estimates the causal parameter by iteratively fitting an accelerated
failure time (AFT) model to counterfactual *unswitched* survival times,
and derives the adjusted hazard ratio from the Cox model using
counterfactual *unswitched* survival times based on the estimated causal
parameter.

## Usage

``` r
ipe(
  data,
  id = "id",
  stratum = "",
  time = "time",
  event = "event",
  treat = "treat",
  rx = "rx",
  censor_time = "censor_time",
  base_cov = "",
  aft_dist = "weibull",
  strata_main_effect_only = TRUE,
  low_psi = -2,
  hi_psi = 2,
  treat_modifier = 1,
  recensor = TRUE,
  admin_recensor_only = TRUE,
  autoswitch = TRUE,
  root_finding = "brent",
  alpha = 0.05,
  ties = "efron",
  tol = 1e-06,
  boot = FALSE,
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

  - `rx`: The proportion of time on active treatment.

  - `censor_time`: The administrative censoring time. It should be
    provided for all subjects including those who had events.

  - `base_cov`: The baseline covariates (excluding treat).

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

- rx:

  The name of the rx variable in the input data.

- censor_time:

  The name of the censor_time variable in the input data.

- base_cov:

  The names of baseline covariates (excluding treat) in the input data
  for the causal AFT model and the outcome Cox model.

- aft_dist:

  The assumed distribution for time to event for the AFT model. Options
  include "exponential", "weibull" (default), "loglogistic", and
  "lognormal".

- strata_main_effect_only:

  Whether to only include the strata main effects in the AFT model.
  Defaults to `TRUE`, otherwise all possible strata combinations will be
  considered in the AFT model.

- low_psi:

  The lower limit of the causal parameter.

- hi_psi:

  The upper limit of the causal parameter.

- treat_modifier:

  The optional sensitivity parameter for the constant treatment effect
  assumption.

- recensor:

  Whether to apply recensoring to counterfactual survival times.
  Defaults to `TRUE`.

- admin_recensor_only:

  Whether to apply recensoring to administrative censoring times only.
  Defaults to `TRUE`. If `FALSE`, recensoring will be applied to the
  actual censoring times for dropouts.

- autoswitch:

  Whether to exclude recensoring for treatment arms with no switching.
  Defaults to `TRUE`.

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

- boot:

  Whether to use bootstrap to obtain the confidence interval for hazard
  ratio. Defaults to `FALSE`, in which case, the confidence interval
  will be constructed to match the log-rank test p-value.

- n_boot:

  The number of bootstrap samples.

- seed:

  The seed to reproduce the bootstrap results.

- nthreads:

  The number of threads to use in bootstrapping (0 means the default
  RcppParallel behavior).

## Value

A list with the following components:

- `psi`: The estimated causal parameter.

- `psi_CI`: The confidence interval for `psi`.

- `psi_CI_type`: The type of confidence interval for `psi`, i.e.,
  "log-rank p-value" or "bootstrap".

- `pvalue`: The two-sided p-value.

- `pvalue_type`: The type of two-sided p-value for treatment effect,
  i.e., "log-rank" or "bootstrap".

- `hr`: The estimated hazard ratio from the Cox model.

- `hr_CI`: The confidence interval for hazard ratio.

- `hr_CI_type`: The type of confidence interval for hazard ratio, either
  "log-rank p-value" or "bootstrap".

- `event_summary`: A data frame containing the count and percentage of
  deaths and switches by treatment arm.

- `Sstar`: A data frame containing the counterfactual untreated survival
  times and event indicators for each treatment group. The variables
  include `id`, `stratum`, `"t_star"`, `"d_star"`, `"treated"`,
  `base_cov`, and `treat`.

- `kmstar`: A data frame containing the Kaplan-Meier estimates based on
  the counterfactual untreated survival times by treatment arm.

- `data_aft`: The input data for the AFT model for estimating `psi`. The
  variables include `id`, `stratum`, `"t_star"`, `"d_star"`,
  `"treated"`, `base_cov`, and `treat`.

- `fit_aft`: The fitted AFT model for estimating `psi`.

- `res_aft`: The deviance residuals from the fitted AFT model.

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

- `fail_boots`: The indicators for failed bootstrap samples if `boot` is
  `TRUE`.

- `fail_boots_data`: The data for failed bootstrap samples if `boot` is
  `TRUE`.

- `hr_boots`: The bootstrap hazard ratio estimates if `boot` is `TRUE`.

- `psi_boots`: The bootstrap `psi` estimates if `boot` is `TRUE`.

## Details

Assuming one-way switching from control to treatment, the hazard ratio
and confidence interval under a no-switching scenario are obtained as
follows:

- Estimate the causal parameter \\\psi\\ by iteratively fitting an AFT
  model to the observed survival times for the treatment arm and the
  counterfactual survival times for the control arm: \$\$U\_{i,\psi} =
  T\_{C_i} + e^{\psi}T\_{E_i}\$\$

- Compute counterfactual survival times for control patients using the
  estimated \\\psi\\.

- Fit a Cox model to the observed survival times for the treatment group
  and the counterfactual survival times for the control group to
  estimate the hazard ratio.

- Obtain the confidence interval for the hazard ratio using either the
  ITT log-rank test p-value or bootstrap. When bootstrapping, the
  interval and p-value are derived from a t-distribution with
  `n_boot - 1` degrees of freedom.

The causal parameter \\\psi\\ is estimated by solving the fixed-point
equation \\\psi\_{new}(\psi) = \psi\\, where \\\psi\_{new}\\ is the
treatment coefficient from the AFT model fitted at a given \\\psi\\.
Brent's method (or bisection) is applied to the bracketed interval
\\\[\psi\_{lo}, \psi\_{hi}\]\\, where \\\psi\_{new}(\psi) - \psi\\ is
positive at \\\psi\_{lo}\\ and negative at \\\psi\_{hi}\\. When multiple
fixed points exist in the bracket, the algorithm converges to one of
them without any selection preference for the root closest to zero.

## References

Michael Branson and John Whitehead. Estimating a treatment effect in
survival studies in which patients switch treatment. Statistics in
Medicine. 2002;21(17):2449-2463.

Ian R White. Letter to the Editor: Estimating treatment effects in
randomized trials with treatment switching. Statistics in Medicine.
2006;25(9):1619-1622.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
library(dplyr)
#> 
#> Attaching package: ‘dplyr’
#> The following objects are masked from ‘package:stats’:
#> 
#>     filter, lag
#> The following objects are masked from ‘package:base’:
#> 
#>     intersect, setdiff, setequal, union

# Example 1: one-way treatment switching (control to active)

data <- immdef %>% mutate(rx = 1-xoyrs/progyrs)

fit1 <- ipe(
  data, id = "id", time = "progyrs", event = "prog", treat = "imm", 
  rx = "rx", censor_time = "censyrs", aft_dist = "weibull",
  boot = FALSE)

fit1
#>             n event_n event_pct switch_n switch_pct event_out_n event_out_pct
#> Control   500     169      33.8      189       37.8         142          28.4
#> Treatment 500     143      28.6        0        0.0         143          28.6
#>           imm
#> Control     0
#> Treatment   1
#> 
#>                            Estimate Lower 95% Upper 95%
#> Causal parameter psi       -0.183   -0.370    0.004    
#> Causal survival time ratio 1.201    0.996     1.448    
#> Hazard ratio (HR)          0.766    0.583     1.006    
#> P-value (log-rank)         0.0556                      

# Example 2: two-way treatment switching (illustration only)

# the eventual survival time
shilong1 <- shilong %>%
  arrange(bras.f, id, tstop) %>%
  group_by(bras.f, id) %>%
  slice(n()) %>%
  select(-c("ps", "ttc", "tran"))

shilong2 <- shilong1 %>%
  mutate(rx = ifelse(co, ifelse(bras.f == "MTA", dco/ady, 
                                1 - dco/ady),
                     ifelse(bras.f == "MTA", 1, 0)))

fit2 <- ipe(
  shilong2, id = "id", time = "tstop", event = "event",
  treat = "bras.f", rx = "rx", censor_time = "dcut",
  base_cov = c("agerand", "sex.f", "tt_Lnum", "rmh_alea.c",
               "pathway.f"),
  aft_dist = "weibull", boot = FALSE)

fit2
#>             n event_n event_pct switch_n switch_pct event_out_n event_out_pct
#> Control    93      63      67.7       68       73.1          58          62.4
#> Treatment 100      67      67.0       25       25.0          63          63.0
#>           bras.f
#> Control       CT
#> Treatment    MTA
#> 
#>                            Estimate Lower 95% Upper 95%
#> Causal parameter psi       0.953    -0.457    2.363    
#> Causal survival time ratio 0.386    0.094     1.579    
#> Hazard ratio (HR)          2.716    0.620     11.900   
#> P-value (log-rank)         0.1851                      
```
