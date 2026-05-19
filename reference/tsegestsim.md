# Simulate Survival Data for Two-Stage Estimation with g-estimation

Obtains the simulated data for baseline prognosis, disease progression,
treatment switching, death, and time-dependent covariates.

## Usage

``` r
tsegestsim(
  n = 500L,
  allocation1 = 2L,
  allocation2 = 1L,
  pbprog = 0.5,
  trtlghr = -0.5,
  bprogsl = 0.3,
  shape1 = 1.8,
  scale1 = 360,
  shape2 = 1.7,
  scale2 = 688,
  pmix = 0.5,
  admin = 5000,
  pcatnotrtbprog = 0.5,
  pcattrtbprog = 0.25,
  pcatnotrt = 0.2,
  pcattrt = 0.1,
  catmult = 0.5,
  tdxo = 1,
  ppoor = 0.1,
  pgood = 0.04,
  ppoormet = 0.4,
  pgoodmet = 0.2,
  xomult = 1.4188308,
  milestone = 546,
  seed = 0L
)
```

## Arguments

- n:

  The total sample size for two treatment arms combined.

- allocation1:

  The number of subjects in the active treatment group in a
  randomization block.

- allocation2:

  The number of subjects in the control group in a randomization block.

- pbprog:

  The probability of having poor prognosis at baseline.

- trtlghr:

  The treatment effect in terms of log hazard ratio.

- bprogsl:

  The poor prognosis effect in terms of log hazard ratio.

- shape1:

  The shape parameter for the Weibull event distribution for the first
  component.

- scale1:

  The scale parameter for the Weibull event distribution for the first
  component.

- shape2:

  The shape parameter for the Weibull event distribution for the second
  component.

- scale2:

  The scale parameter for the Weibull event distribution for the second
  component.

- pmix:

  The mixing probability of the first component Weibull distribution.

- admin:

  The administrative censoring time.

- pcatnotrtbprog:

  The probability of developing metastatic disease on control treatment
  with poor baseline prognosis.

- pcattrtbprog:

  The probability of developing metastatic disease on active treatment
  with poor baseline prognosis.

- pcatnotrt:

  The probability of developing metastatic disease on control treatment
  with good baseline prognosis.

- pcattrt:

  The probability of developing metastatic disease on active treatment
  with good baseline prognosis.

- catmult:

  The impact of metastatic disease on shortening remaining survival
  time.

- tdxo:

  Whether treatment crossover depends on time-dependent covariates
  between disease progression and treatment switching.

- ppoor:

  The probability of switching for poor baseline prognosis with no
  metastatic disease.

- pgood:

  The probability of switching for good baseline prognosis with no
  metastatic disease.

- ppoormet:

  The probability of switching for poor baseline prognosis after
  developing metastatic disease.

- pgoodmet:

  The probability of switching for good baseline prognosis after
  developing metastatic disease.

- xomult:

  The direct effect of crossover on extending remaining survival time.

- milestone:

  The milestone to calculate restricted mean survival time.

- seed:

  The seed to reproduce the simulation results.

## Value

A list with two data frames.

- `sumdata`: A summary data frame with the following variables:

  - `simtrueconstmean`: The true control group restricted mean survival
    time (RMST).

  - `simtrueconstlb`: The lower bound for control group RMST.

  - `simtrueconstub`: The upper bound for control group RMST.

  - `simtrueconstse`: The standard error for control group RMST.

  - `simtrueexpstmean`: The true experimental group restricted mean
    survival time (RMST).

  - `simtrueexpstlb`: The lower bound for experimental group RMST.

  - `simtrueexpstub`: The upper bound for experimental group RMST.

  - `simtrueexpstse`: The standard error for experimental group RMST.

  - `simtrue_coxwbprog_hr`: The treatment hazard ratio from the Cox
    model adjusting for baseline prognosis.

  - `simtrue_cox_hr`: The treatment hazard ratio from the Cox model
    without adjusting for baseline prognosis.

  - `simtrue_aftwbprog_af`: The average acceleration factor from the
    Weibull AFT model adjusting for baseline prognosis.

  - `simtrue_aft_af`: The average acceleration factor from the Weibull
    AFT model without adjusting for baseline prognosis.

- `paneldata`: A counting process style subject-level data frame with
  the following variables:

  - `id`: The subject ID.

  - `trtrand`: The randomized treatment arm.

  - `bprog`: Whether the patient had poor baseline prognosis.

  - `tstart`: The left end of time interval.

  - `tstop`: The right end of time interval.

  - `event`: Whether the patient died at the end of the interval.

  - `timeOS`: The observed survival time.

  - `died`: Whether the patient died during the study.

  - `progressed`: Whether the patient had disease progression.

  - `timePFSobs`: The observed time of disease progression at regular
    scheduled visits.

  - `progtdc`: The time-dependent covariate for progression.

  - `catevent`: Whether the patient developed metastatic disease.

  - `cattime`: When the patient developed metastatic disease.

  - `cattdc`: The time-dependent covariate for cat event.

  - `xo`: Whether the patient switched treatment.

  - `xotime`: When the patient switched treatment.

  - `xotdc`: The time-dependent covariate for treatment switching.

  - `censor_time`: The administrative censoring time.

## References

NR Latimer, IR White, K Tilling, and U Siebert. Improved two-stage
estimation to adjust for treatment switching in randomised trials:
g-estimation to address time-dependent confounding. Statistical Methods
in Medical Research. 2020;29(10):2900-2918.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r

sim1 <- tsegestsim(
  n = 500, allocation1 = 2, allocation2 = 1, pbprog = 0.5, 
  trtlghr = -0.5, bprogsl = 0.3, shape1 = 1.8, 
  scale1 = 360, shape2 = 1.7, scale2 = 688, 
  pmix = 0.5, admin = 5000, pcatnotrtbprog = 0.5, 
  pcattrtbprog = 0.25, pcatnotrt = 0.2, pcattrt = 0.1, 
  catmult = 0.5, tdxo = 1, ppoor = 0.1, pgood = 0.04, 
  ppoormet = 0.4, pgoodmet = 0.2, xomult = 1.4188308, 
  milestone = 546, seed = 2000)
```
