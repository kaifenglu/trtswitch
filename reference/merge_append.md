# Append columns from B to A by keyed join (fast, data.table-based)

A and B may be data.frame or data.table.

- If A is a data.table it will be modified by reference.

- If A is a data.frame a modified data.frame is returned.

## Usage

``` r
merge_append(A, B, by_vars, new_vars, overwrite = TRUE, first_match = TRUE)
```

## Arguments

- A:

  left table (data.frame or data.table) — rows preserved/order preserved

- B:

  right table (data.frame or data.table) — must contain by_vars and
  new_vars

- by_vars:

  character vector of join columns (e.g. c("id","tstart","tstop"))

- new_vars:

  character vector of column names from B to append to A

- overwrite:

  logical: if TRUE, overwrite any existing columns in A with same avar
  names (default TRUE)

- first_match:

  logical: if TRUE use first match when B has duplicate keys; if FALSE,
  error on duplicates

## Value

If A was a data.frame, returns modified data.frame; if A was data.table
returns invisible(NULL).
