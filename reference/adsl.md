# Baseline subject-level data

This data set contains baseline subject-level data. Of note, `PDDT` can
be derived from the `ADT` variable of the `ADTTE` data set by selecting
`PARAMCD == "INPFS" & CNSR == 0 & EVNTDESC == "PROGRESSIVE DISEASE"`.
Additionally, `OSDT` and `DIED` can be derived from the `ADT` and `CNSR`
variables of the `ADTTE` data set by selecting `PARAMCD == "OS"`.

## Usage

``` r
adsl
```

## Format

An object of class `tbl_df` (inherits from `tbl`, `data.frame`) with 412
rows and 12 columns.

## Details

- `SUBJID`:

  subject ID

- `SEX`:

  sex: "M" or "F"

- `STRAT1V`:

  stratification factor 1: ECOG PS

- `STRAT2V`:

  stratification factor 2: inv. chosen chemotherapy

- `RANDDT`:

  randomization date

- `TRT01P`:

  planned treatment: Active or Placebo

- `TRTSDT`:

  treatment start date

- `PDDT`:

  date of disease progression

- `XODT`:

  date of treatment crossover

- `OSDT`:

  date of death or censoring

- `DIED`:

  whether the patient died

- `DCUTDT`:

  date of data cut
