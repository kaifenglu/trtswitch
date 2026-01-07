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
