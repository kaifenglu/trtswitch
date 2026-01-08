#' helper: process a base covariate specification and ensure design columns exist in df
#' @keywords internal
process_cov <- function(base_cov, df) {
  miss <- length(base_cov) == 0 || (length(base_cov) == 1 && base_cov[1] == "")
  if (miss) {
    return(list(varnames = "", vnames = character(0), df = df))
  }
  
  fml <- as.formula(paste("~", paste(base_cov, collapse = "+")))
  vnames <- rownames(attr(terms(fml), "factors"))
  
  # QUICK PATH: if all covariates present in df and are numeric, avoid model.matrix
  cov_present <- base_cov %in% names(df)
  all_numeric <- FALSE
  if (all(cov_present)) {
    all_numeric <- all(vapply(df[base_cov], is.numeric, logical(1)))
  }
  
  if (all_numeric) {
    varnames <- base_cov
    return(list(varnames = varnames, vnames = vnames, df = df))
  }
  
  # FALLBACK: robust behavior via model.frame + model.matrix
  mf <- model.frame(fml, data = df, na.action = na.pass)
  mm <- model.matrix(fml, mf)
  colnames(mm) <- make.names(colnames(mm))
  varnames <- if (ncol(mm) > 1) colnames(mm)[-1] else character(0)
  
  missing_cols <- setdiff(varnames, names(df))
  if (length(missing_cols) > 0) {
    for (vn in missing_cols) df[[vn]] <- mm[, vn, drop = TRUE]
  }
  
  list(varnames = varnames, vnames = vnames, df = df)
}


#' Append columns from B to A by keyed join (fast, data.table-based)
#'
#' A and B may be data.frame or data.table.
#' - If A is a data.table it will be modified by reference.
#' - If A is a data.frame a modified data.frame is returned.
#'
#' @param A left table (data.frame or data.table) — rows preserved/order preserved
#' @param B right table (data.frame or data.table) — must contain by_vars and new_vars
#' @param by_vars character vector of join columns (e.g. c("id","tstart","tstop"))
#' @param new_vars character vector of column names from B to append to A
#' @param overwrite logical: if TRUE, overwrite any existing columns in A 
#'   with same avar names (default TRUE)
#' @param first_match logical: if TRUE use first match when B has duplicate keys; 
#'   if FALSE, error on duplicates
#' @return If A was a data.frame, returns modified data.frame; if A was 
#'   data.table returns invisible(NULL).
merge_append <- function(A, B, by_vars, new_vars,
                         overwrite = TRUE, first_match = TRUE) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("data.table is required")
  }
  dt <- data.table::as.data.table
  
  # prepare B: select only required columns to reduce memory
  if (!all(by_vars %in% colnames(B))) stop("Some by_vars not in B")
  if (!all(new_vars %in% colnames(B))) stop("Some new_vars not in B")
  Bdt <- dt(B)[, c(by_vars, new_vars), with = FALSE]
  
  # set key on B for extremely fast lookup (by reference)
  data.table::setkeyv(Bdt, by_vars)
  
  # prepare A (remember whether original was data.table)
  A_is_dt <- data.table::is.data.table(A)
  Adt <- if (A_is_dt) A else dt(A)
  
  # optional duplicate check on B: if first_match=FALSE error when duplicates in B
  if (!first_match) {
    dup_keys <- Bdt[duplicated(Bdt, by = by_vars), .N, by = by_vars]
    if (nrow(dup_keys) > 0) stop(paste("Duplicate keys found in B for join;",
                                       "set first_match=TRUE to allow first-match"))
  }
  
  # perform keyed join: Bdt[Adt] returns rows aligned to Adt (keeps Adt order)
  # For speed we only need the new_vars from the join result.
  joined <- Bdt[Adt, on = by_vars, nomatch = NA, mult = "first"]
  
  # If A already contains any of new_vars and overwrite == FALSE, skip those
  if (!overwrite) {
    new_vars_to_assign <- setdiff(new_vars, intersect(new_vars, colnames(Adt)))
  } else {
    new_vars_to_assign <- new_vars
  }
  
  # Assign columns into Adt by reference (fast)
  for (col in new_vars_to_assign) {
    # if the column name already exists in Adt and overwrite==TRUE, we replace it
    # joined[[col]] has the values aligned to Adt rows
    Adt[, (col) := joined[[col]]]
  }
  
  # If input A was data.frame, return modified data.frame; 
  # if data.table, we modified by-ref
  if (!A_is_dt) {
    # convert back to data.frame
    return(as.data.frame(Adt))
  } else {
    invisible(NULL)
  }
}

#' @title Split a survival data set at specified cut points
#' @description For a given survival dataset and specified cut times, 
#' each record is split into multiple subrecords at each cut time. 
#' The resulting dataset is in counting process format, with each 
#' subrecord containing a start time, stop time, and event status.
#' This is adapted from the survsplit.c function from the survival package.
#'
#' @param tstart The starting time of the time interval for 
#'   counting-process data.
#' @param tstop The stopping time of the time interval for 
#'   counting-process data.
#' @param cut The vector of cut points.
#'
#' @return A data frame with the following variables:
#'
#' * \code{row}: The row number of the observation in the input data 
#'   (starting from 0).
#'
#' * \code{start}: The starting time of the resulting subrecord.
#'
#' * \code{end}: The ending time of the resulting subrecord.
#'
#' * \code{censor}: Whether the subrecord lies strictly within a record
#'   in the input data (1 for all but the last interval and 0 for the 
#'   last interval).
#'
#' * \code{interval}: The interval number derived from cut (starting 
#'   from 0 if the interval lies to the left of the first cutpoint).
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @keywords internal
#'
#' @examples
#'
#' survsplit(15, 60, c(10, 30, 40))
#'
#' @export
survsplit <- function(tstart, tstop, cut) {
  survsplitRcpp(tstart, tstop, cut)
}

# LOCF within each id and paramcd
locf_safe <- function(x) {
  v <- !is.na(x)
  # Find indices of non-NA values
  idx <- which(v)
  
  # If no non-NA values exist, or the first value is already the only value
  if (length(idx) == 0) return(x)
  
  # Identify where the first non-NA value starts
  first_val_idx <- idx[1]
  
  # For indices before the first non-NA, keep them as NA
  # For indices after, use the cumulative sum to find the latest non-NA
  res <- x
  if (length(idx) > 0) {
    # Generate mapping: 0 for leading NAs, then 1, 1, 2, 3...
    vv <- cumsum(v)
    # Subset only the portion from the first non-NA onwards
    res[first_val_idx:length(x)] <- x[idx[vv[first_val_idx:length(x)]]]
  }
  return(res)
}
