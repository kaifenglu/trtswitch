# Split a survival data set at specified cut points

For a given survival dataset and specified cut times, each record is
split into multiple subrecords at each cut time. The resulting dataset
is in counting process format, with each subrecord containing a start
time, stop time, and event status. This is adapted from the survsplit.c
function from the survival package.

## Usage

``` r
survsplit(tstart, tstop, cut)
```

## Arguments

- tstart:

  The starting time of the time interval for counting-process data.

- tstop:

  The stopping time of the time interval for counting-process data.

- cut:

  The vector of cut points.

## Value

A data frame with the following variables:

- `row`: The row number of the observation in the input data (starting
  from 0).

- `start`: The starting time of the resulting subrecord.

- `end`: The ending time of the resulting subrecord.

- `censor`: Whether the subrecord lies strictly within a record in the
  input data (1 for all but the last interval and 0 for the last
  interval).

- `interval`: The interval number derived from cut (starting from 0 if
  the interval lies to the left of the first cutpoint).

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
survsplit(15, 60, c(10, 30, 40))
#>   row start end censor interval
#> 1   0    15  30      1        1
#> 2   0    30  40      1        2
#> 3   0    40  60      0        3
```
