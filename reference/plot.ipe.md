# Plot method for ipe objects

Generate box plot for AFT model deviance residuals and Kaplan-Meier (KM)
plot for potential outcomes of an ipe object.

## Usage

``` r
# S3 method for class 'ipe'
plot(x, time_unit = "day", show_hr = TRUE, show_risk = TRUE, ...)
```

## Arguments

- x:

  An object of class `ipe`.

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

A list of two ggplot2 objects, one for box plot and the other for KM
plot.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>
