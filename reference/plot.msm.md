# Plot method for msm objects

Generate histogram of weights and weighted Kaplan-Meier (KM) plot for
censored outcomes of an msm object.

## Usage

``` r
# S3 method for class 'msm'
plot(x, time_unit = "day", show_hr = TRUE, show_risk = TRUE, ...)
```

## Arguments

- x:

  An object of class `msm`.

- time_unit:

  The time unit used in the input data. Options are "day" (default),
  "week", "month", or "year".

- show_hr:

  Logical; whether to show hazard ratio on the KM plot. Default is TRUE.

- show_risk:

  Logical; whether to show number at risk table below the KM plot.
  Default is TRUE.

- ...:

  Ensures that all arguments starting from "..." are named.

## Value

A list of two ggplot2 objects, one for histogram and the other for KM
plot.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>
