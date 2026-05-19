# Prepare Survival Data With Time-Dependent Covariates

This function prepares a counting-process style survival dataset for
analyses with time-dependent covariates. It merges baseline and
longitudinal data, fills in missing covariate values using
last-observation-carried-forward (LOCF), restricts to time points where
covariates change (optional), and constructs `tstart`, `tstop`, and
`event` variables suitable for use in survival models.

## Usage

``` r
preptdc(
  adsl,
  adtdc,
  id = "SUBJID",
  randdt = "RANDDT",
  trtsdt = "TRTSDT",
  pddt = "PDDT",
  xodt = "XODT",
  osdt = "OSDT",
  died = "DIED",
  dcutdt = "DCUTDT",
  paramcd = "PARAMCD",
  adt = "ADT",
  aval = "AVAL",
  nodup = TRUE,
  offset = TRUE
)
```

## Arguments

- adsl:

  A data set containing baseline subject-level information. It should
  include, at a minimum, subject ID (`id`), randomization date
  (`randdt`), treatment start date (`trtsdt`), progression date
  (`pddt`), treatment switch date (`xodt`), survival outcome (`osdt`,
  `died`), and data cut-off date (`dcutdt`).

- adtdc:

  A data set containing longitudinal time-dependent covariate data, with
  subject ID (`id`), parameter code (`paramcd`), analysis date (`adt`),
  and covariate value (`aval`).

- id:

  Character string specifying the column name for subject ID.

- randdt:

  Character string specifying the column name for randomization date.

- trtsdt:

  Character string specifying the column name for treatment start date.

- pddt:

  Character string specifying the column name for progression date.

- xodt:

  Character string specifying the column name for treatment
  crossover/switch date.

- osdt:

  Character string specifying the column name for overall survival date
  (death date or last known alive date).

- died:

  Character string specifying the column name for death indicator (0 =
  alive/censored, 1 = died).

- dcutdt:

  Character string specifying the column name for data cut-off date.

- paramcd:

  Character string specifying the column name for parameter code
  (identifying different covariates).

- adt:

  Character string specifying the column name for analysis date in the
  time-dependent covariate dataset.

- aval:

  Character string specifying the column name for analysis value
  (covariate values).

- nodup:

  Logical; if `TRUE` (default), only rows where at least one covariate
  changes compared to the previous row (within each subject) are
  retained, along with the first row per subject (baseline).

- offset:

  Logical; if `TRUE` (default), add 1-day offset when computing analysis
  day variables (`ady`, `osdy`, etc.).

## Value

A data set with one row per subject and time interval, including:

- `tstart`, `tstop` — interval start and stop times (days from
  randomization).

- `event` — event indicator (0/1).

- Covariates expanded to wide format.

- Auxiliary variables such as progression indicator (`pd`), treatment
  switch indicator (`swtrt`), and administrative censoring time.

## Details

The function performs the following steps:

1.  Merge `adsl` and `adtdc` to obtain randomization date and treatment
    start date.

2.  Define `adt2` as `adt` if `adt > trtsdt`, and `randdt` if
    `adt <= trtsdt` (i.e., baseline time point). This ensures that the
    baseline covariate value is the last non-missing value at or before
    the treatment start date. Post-baseline covariate values are
    anchored at their actual analysis dates. The first record per
    subject corresponds to survival time zero at randomization and
    ensures availability of baseline covariates at randomization.

3.  Keep the last record per subject, `adt2`, and `paramcd`.

4.  Construct a complete skeleton so all covariates are present for each
    subject and time point.

5.  Fill missing covariate values using LOCF.

6.  Pivot to wide format with one row per subject and time point.

7.  Optionally drop rows without covariate changes (`nodup = TRUE`).

8.  Merge survival outcomes from `adsl`.

9.  Compute time-to-event variables (`ady`, `osdy`, etc.), as well as
    counting-process style variables `tstart`, `tstop`, and `event`.

## Examples

``` r

surv_data <- preptdc(adsl, adtdc, nodup = TRUE)
head(surv_data)
#>       SUBJID     RANDDT       adt2 ECOG101 LDH SEX STRAT1V     STRAT2V TRT01P
#> 1 000028-240 2020-07-27 2020-07-27       1 202   F       1 Carboplatin Active
#> 2 000028-240 2020-07-27 2020-08-18       1 278   F       1 Carboplatin Active
#> 3 000028-240 2020-07-27 2020-09-08       1 269   F       1 Carboplatin Active
#> 4 000028-240 2020-07-27 2020-09-29       1 223   F       1 Carboplatin Active
#> 5 000028-240 2020-07-27 2020-10-08       1 267   F       1 Carboplatin Active
#> 6 000028-240 2020-07-27 2020-10-27       1 164   F       1 Carboplatin Active
#>       TRTSDT       PDDT       XODT       OSDT DIED     DCUTDT ady osdy tstart
#> 1 2020-07-27 2021-05-23 2021-06-17 2021-09-15    1 2023-04-20   1  416      1
#> 2 2020-07-27 2021-05-23 2021-06-17 2021-09-15    1 2023-04-20  23  416     23
#> 3 2020-07-27 2021-05-23 2021-06-17 2021-09-15    1 2023-04-20  44  416     44
#> 4 2020-07-27 2021-05-23 2021-06-17 2021-09-15    1 2023-04-20  65  416     65
#> 5 2020-07-27 2021-05-23 2021-06-17 2021-09-15    1 2023-04-20  74  416     74
#> 6 2020-07-27 2021-05-23 2021-06-17 2021-09-15    1 2023-04-20  93  416     93
#>   tstop event   pd pd_time swtrt swtrt_time censor_time
#> 1    23     0 TRUE     301  TRUE        326         998
#> 2    44     0 TRUE     301  TRUE        326         998
#> 3    65     0 TRUE     301  TRUE        326         998
#> 4    74     0 TRUE     301  TRUE        326         998
#> 5    93     0 TRUE     301  TRUE        326         998
#> 6   114     0 TRUE     301  TRUE        326         998
```
