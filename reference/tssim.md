# Simulate Data for Treatment Switching

Simulates data for studies involving treatment switching, incorporating
time-dependent confounding. The generated data can be used to evaluate
methods for handling treatment switching in survival analysis.

## Usage

``` r
tssim(
  tdxo = FALSE,
  coxo = TRUE,
  allocation1 = 1L,
  allocation2 = 1L,
  p_X_1 = 0.3,
  p_X_0 = 0.3,
  rate_T = 0.002,
  beta1 = -0.5,
  beta2 = 0.3,
  gamma0 = 0.3,
  gamma1 = -0.9,
  gamma2 = 0.7,
  gamma3 = 1.1,
  gamma4 = -0.8,
  zeta0 = -3.5,
  zeta1 = 0.5,
  zeta2 = 0.2,
  zeta3 = -0.4,
  alpha0 = 0.5,
  alpha1 = 0.5,
  alpha2 = 0.4,
  theta1_1 = -0.4,
  theta1_0 = -0.4,
  theta2 = 0.2,
  rate_C = 8.55e-05,
  accrualTime = 0L,
  accrualIntensity = NA_real_,
  followupTime = NA_real_,
  fixedFollowup = FALSE,
  plannedTime = 1350,
  days = 30,
  n = 500L,
  NSim = 1000L,
  seed = 0L
)
```

## Arguments

- tdxo:

  Logical indicator for timing of treatment switching:

  - 1: Treatment switching can occur at or after disease progression.

  - 0: Treatment switching is restricted to the time of disease
    progression.

- coxo:

  Logical indicator for arm-specific treatment switching:

  - 1: Treatment switching occurs only in the control arm.

  - 0: Treatment switching is allowed in both arms.

- allocation1:

  Number of subjects in the active treatment group in a randomization
  block. Defaults to 1 for equal randomization.

- allocation2:

  Number of subjects in the control group in a randomization block.
  Defaults to 1 for equal randomization.

- p_X_1:

  Probability of poor baseline prognosis in the experimental arm.

- p_X_0:

  Probability of poor baseline prognosis in the control arm.

- rate_T:

  Baseline hazard rate for time to death.

- beta1:

  Log hazard ratio for randomized treatment (`R`).

- beta2:

  Log hazard ratio for baseline covariate (`X`).

- gamma0:

  Intercept for the time-dependent covariate model (`L`).

- gamma1:

  Coefficient for lagged treatment switching (`Alag`) in the `L` model.

- gamma2:

  Coefficient for lagged `L` (`Llag`) in the `L` model.

- gamma3:

  Coefficient for baseline covariate (`X`) in the `L` model.

- gamma4:

  Coefficient for randomized treatment (`R`) in the `L` model.

- zeta0:

  Intercept for the disease progression model (`Z`).

- zeta1:

  Coefficient for `L` in the `Z` model.

- zeta2:

  Coefficient for baseline covariate (`X`) in the `Z` model.

- zeta3:

  Coefficient for randomized treatment (`R`) in the `Z` model.

- alpha0:

  Intercept for the treatment switching model (`A`).

- alpha1:

  Coefficient for `L` in the `A` model.

- alpha2:

  Coefficient for baseline covariate (`X`) in the `A` model.

- theta1_1:

  Negative log time ratio for `A` for the experimental arm.

- theta1_0:

  Negative log time ratio for `A` for the control arm.

- theta2:

  Negative log time ratio for `L`.

- rate_C:

  Hazard rate for random (dropout) censoring.

- accrualTime:

  A vector that specifies the starting time of piecewise Poisson
  enrollment time intervals. Must start with 0, e.g., `c(0, 3)` breaks
  the time axis into 2 accrual intervals: \[0, 3) and \[3, Inf).

- accrualIntensity:

  A vector of accrual intensities. One for each accrual time interval.

- followupTime:

  Follow-up time for a fixed follow-up design.

- fixedFollowup:

  Whether a fixed follow-up design is used. Defaults to 0 for variable
  follow-up.

- plannedTime:

  The calendar time for the analysis.

- days:

  Number of days in each treatment cycle.

- n:

  Number of subjects per simulation.

- NSim:

  Number of simulated datasets.

- seed:

  Random seed for reproducibility.

## Value

A list of data frames, each containing simulated longitudinal covariate,
pd status, alternative therapy status, and event history data with the
following variables:

- `id`: Subject identifier.

- `arrival_time`: The enrollment time for the subject.

- `trtrand`: Randomized treatment assignment (0 = control, 1 =
  experimental)

- `bprog`: Baseline prognosis (0 = good, 1 = poor).

- `tpoint`: Treatment cycle index.

- `tstart`: Start day of the treatment cycle.

- `tstop`: End day of the treatment cycle.

- `L`: Time-dependent covariate at `tstart` predicting survival and
  switching; affected by treatment switching.

- `Llag`: Lagged value of `L`.

- `Z`: Disease progression status at `tstart`.

- `A`: Treatment switching status at `tstart`.

- `Alag`: Lagged value of `A`.

- `event`: Death indicator at `tstop`.

- `timeOS`: Observed time to death or censoring.

- `died`: Indicator of death by end of follow-up.

- `progressed`: Indicator of disease progression by end of follow-up.

- `timePD`: Observed time to progression or censoring.

- `xo`: Indicator for whether treatment switching occurred.

- `xotime`: Time of treatment switching (if applicable).

- `censor_time`: Administrative censoring time.

## Details

The simulation algorithm is adapted from Xu et al. (2022) and simulates
disease progression status while incorporating the multiplicative
effects of both baseline and time-dependent covariates on survival time.
The design options `tdxo` and `coxo` specify the timing of treatment
switching and the study arm eligibility for switching, respectively.
Time is measured in days.

In a fixed follow-up design, all subjects share the same follow-up
duration. In contrast, under a variable follow-up design, follow-up time
also depends on each subject's enrollment date. The number of treatment
cycles for a subject is determined by dividing the follow-up time by the
number of days in each cycle.

1.  At randomization, subjects are assigned to treatment arms using
    block randomization, with `allocation1` patients assigned to active
    treatment and `allocation2` to control within each randomized block.
    A baseline covariate is also generated for each subject: \$\$X_i
    \sim \mbox{Bernoulli}(p_1 R_i + p_0 (1-R_i))\$\$

2.  The initial survival time is drawn from an exponential distribution
    with hazard: \$\$\lambda_T \exp(\beta_1 R_i + \beta_2 X_i)\$\$ We
    define the event indicator at cycle \\j\\ as \$\$Y\_{i,j} = I(T_i
    \leq j\times days)\$\$

3.  The initial states are set to \\L\_{i,0} = 0\\, \\Z\_{i,0} = 0\\,
    \\A\_{i,0} = 0\\, \\Y\_{i,0} = 0\\. For each treatment cycle
    \\j=1,\ldots,J\\, we set \\tstart = (j-1) \times days\\.

4.  Generate time-dependent covariates: \$\$\mbox{logit}
    P(L\_{i,j}=1\|\mbox{history}) = \gamma_0 + \gamma_1 A\_{i,j-1} +
    \gamma_2 L\_{i,j-1} + \gamma_3 X_i + \gamma_4 R_i\$\$

5.  If \\T_i \leq j \times days\\, set \\tstop = T_i\\ and \\Y\_{i,j} =
    1\\, which completes data generation for subject \\i\\.

6.  If \\T_i \> j \times days\\, set \\tstop = j\times days\\,
    \\Y\_{i,j} = 0\\, and perform the following before proceeding to the
    next cycle for the subject.

7.  Generate disease progression status: If \\Z\_{i,j-1} = 0\\,
    \$\$\mbox{logit} P(Z\_{i,j}=1 \| \mbox{history}) = \zeta_0 + \zeta_1
    L\_{i,j} + \zeta_2 X_i + \zeta_3 R_i\$\$ Otherwise, set \\Z\_{i,j} =
    1\\.

8.  Generate alternative therapy status: If \\A\_{i,j-1} = 0\\,
    determine switching eligibility based on design options. If
    switching is allowed: \$\$\mbox{logit} P(A\_{i,j} = 1 \|
    \mbox{history}) = \alpha_0 + \alpha_1 L\_{i,j} + \alpha_2 X_i\$\$ If
    switching is now allowed, set \\A\_{i,j} = 0\\. If \\A\_{i,j-1} =
    1\\, set \\A\_{i,j} = 1\\ (once switched to alternative therapy,
    remain on alternative therapy).

9.  Update survival time based on changes in alternative therapy status
    and time-dependent covariates: \$\$T_i = j\times days + (T_i -
    j\times days) \exp\\ -(\theta\_{1,1}R_i +
    \theta\_{1,0}(1-R_i))(A\_{i,j} - A\_{i,j-1}) -\theta_2 (L\_{i,j} -
    L\_{i,j-1})\\\$\$

Additional random censoring times are generated from an exponential
distribution with hazard rate \\\lambda_C\\.

An extra record is generated when the minimum of the latent survival
time, the random censoring time, and the administrative censoring time
is greater than the number of regular treatment cycles times days per
cycle.

Finally we apply the lag function so that \\Z\_{i,j}\\ and \\A\_{i,j}\\
represent the PD status and alternative therapy status at the start of
cycle \\j\\ (and thus remain appplicable for the entire cycle \\j\\) for
subject \\i\\.

To estimate the true treatment effect in a Cox marginal structural
model, one can set \\\alpha_0 = -\infty\\, which effectively forces
\\A\_{i,j} = 0\\ (disabling treatment switching). The coefficient for
the randomized treatment can then be estimated using a Cox proportional
hazards model.

## References

Jessica G. Young, and Eric J. Tchetgen Tchetgen. Simulation from a known
Cox MSM using standard parametric models for the g-formula. Statistics
in Medicine. 2014;33(6):1001-1014.

NR Latimer, IR White, K Tilling, and U Siebert. Improved two-stage
estimation to adjust for treatment switching in randomised trials:
g-estimation to address time-dependent confounding. Statistical Methods
in Medical Research. 2020;29(10):2900-2918.

Jing Xu, Guohui Liu, and Bingxia Wang. Bias and type I error control in
correcting treatment effect for treatment switching using marginal
structural models in Phse III oncology trials. Journal of
Biopharmaceutical Statistics. 2022;32(6):897-914.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r

library(dplyr)

simulated.data <- tssim(
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
  
simulated.data[[1]] %>% filter(id == 1)
#>   id arrival_time trtrand bprog tpoint tstart tstop L Llag Z A Alag event
#> 1  1            3       1     1      1      0    30 0    0 0 0    0     0
#> 2  1            3       1     1      2     30    60 1    0 0 0    0     0
#> 3  1            3       1     1      3     60    90 0    1 0 0    0     0
#> 4  1            3       1     1      4     90   120 1    0 0 0    0     0
#> 5  1            3       1     1      5    120   150 1    1 0 0    0     0
#> 6  1            3       1     1      6    150   180 1    1 0 0    0     0
#> 7  1            3       1     1      7    180   207 1    1 0 0    0     1
#>   timeOS died progressed timePD xo xotime censor_time
#> 1    207    1          0    NaN  0    NaN        1347
#> 2    207    1          0    NaN  0    NaN        1347
#> 3    207    1          0    NaN  0    NaN        1347
#> 4    207    1          0    NaN  0    NaN        1347
#> 5    207    1          0    NaN  0    NaN        1347
#> 6    207    1          0    NaN  0    NaN        1347
#> 7    207    1          0    NaN  0    NaN        1347
```
