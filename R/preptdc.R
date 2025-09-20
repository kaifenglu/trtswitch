#' @title Prepare Survival Data With Time-Dependent Covariates
#'
#' @description
#' This function prepares a counting-process style survival dataset
#' for analyses with time-dependent covariates. It merges baseline
#' and longitudinal data, fills in missing covariate values using
#' last-observation-carried-forward (LOCF), restricts to time points
#' where covariates change (optional), and constructs `tstart`, `tstop`,
#' and `event` variables suitable for use in survival models.
#'
#' @param adsl A data set containing baseline subject-level data, 
#'   including subject ID (`id`), randomization date (`randdt`),
#'   survival outcome (`osdt`, `died`), progression date (`pddt`), 
#'   treatment switch date (`xodt`), and data cut-off date (`dcutdt`).
#' @param adtdc A data set containing longitudinal
#'   time-dependent covariate data, with subject ID (`id`),
#'   analysis date (`adt`), parameter code (`paramcd`), and covariate
#'   value (`aval`).   
#' @param id Character string specifying the column name for subject ID.
#' @param randdt Character string specifying the column name for 
#'   randomization date.
#' @param osdt Character string specifying the column name for overall
#'   survival date (death date or last known alive date).
#' @param died Character string specifying the column name for death
#'   indicator (0 = alive/censored, 1 = died).
#' @param pddt Character string specifying the column name for progression
#'   date.
#' @param xodt Character string specifying the column name for treatment
#'   crossover/switch date.
#' @param dcutdt Character string specifying the column name for data
#'   cut-off date.
#' @param adt Character string specifying the column name for analysis
#'   date in the time-dependent covariate dataset.
#' @param paramcd Character string specifying the column name for parameter
#'   code (identifying different covariates).
#' @param aval Character string specifying the column name for analysis
#'   value (covariate values).
#' @param nodup Logical; if `TRUE` (default), only rows where at least
#'   one covariate changes compared to the previous row (within each subject)
#'   are retained, along with the first row per subject (baseline).
#' @param offset Logical; if `TRUE` (default), add 1-day offset when 
#'   computing analysis day variables (`ady`, `osdy`, etc.).
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Merge `adsl` and `adtdc` to obtain randomization date.
#'   \item Define `adt2` as the maximum of `adt` and `randdt`.
#'   \item Keep the last record per subject, `adt2`, and `paramcd`.
#'   \item Construct a complete skeleton so all covariates are present
#'         for each subject and time point.
#'   \item Fill missing covariate values using LOCF.
#'   \item Pivot to wide format with one row per subject and time point.
#'   \item Optionally drop rows without covariate changes (`nodup = TRUE`).
#'   \item Merge survival outcomes from `adsl`.
#'   \item Compute time-to-event variables (`ady`, `osdy`, etc.), as well as
#'         counting-process style variables `tstart`, `tstop`, and `event`.
#' }
#'
#' @return A data set with one row per subject and time interval, including:
#' \itemize{
#'   \item `tstart`, `tstop` — interval start and stop times 
#'         (days from randomization).
#'   \item `event` — event indicator (0/1).
#'   \item Covariates expanded to wide format.
#'   \item Auxiliary variables such as progression indicator (`pd`),
#'         treatment switch indicator (`swtrt`), and administrative 
#'         censoring time.
#' }
#'
#' @examples
#' 
#' surv_data <- preptdc(adsl, adtdc, nodup = TRUE)
#' head(surv_data)
#' 
#' @export
preptdc <- function(adsl, adtdc, id = "SUBJID", randdt = "RANDDT", 
                    osdt = "OSDT", died = "DIED", pddt = "PDDT", 
                    xodt = "XODT", dcutdt = "DCUTDT", adt = "ADT", 
                    paramcd = "PARAMCD", aval = "AVAL",
                    nodup = TRUE, offset = TRUE) {
  
  # convert input data into data.table
  data.table::setDT(adsl)
  data.table::setDT(adtdc)
  
  # verify whether the input data have required columns
  cols <- colnames(adsl)
  req_cols <- c(id, randdt, osdt, died, pddt, xodt, dcutdt)
  if (!all(req_cols %in% cols)) {
    stop(paste("The following columns are missing from adsl:", 
               paste(req_cols[!(req_cols %in% cols)], collapse = ", ")))
  }
  
  cols <- colnames(adtdc)
  req_cols <- c(id, adt, paramcd, aval)
  if (!all(req_cols %in% cols)) {
    stop(paste("The following columns are missing from adtdc:", 
               paste(req_cols[!(req_cols %in% cols)], collapse = ", ")))
  }
  
  # obtain randdt
  cols <- c(id, randdt)
  data1 <- merge(adtdc, adsl[, mget(cols)], by = id)
  
  # define adt2 to differentiate baseline and post-baseline time points
  data1[, `:=`(adt2 = pmax(get(adt), get(randdt)))]
  
  # keep last record with each subject, adt2, and paramcd
  data.table::setorderv(data1, c(id, "adt2", paramcd, adt))
  data1 <- data1[, .SD[.N], by = c(id, "adt2", paramcd)]
  
  # build skeleton
  cols <- c(id, randdt, "adt2")
  pars <- unique(data1[[paramcd]])
  data2 <- merge(
    unique(data1[, mget(cols)])[, `:=`(dummy = 1)],
    data.table::data.table(temp_col = pars, dummy = 1),
    by = "dummy", allow.cartesian = TRUE)
  data.table::setnames(data2, "temp_col", paramcd)
  data2[, `:=`(dummy = NULL)]
  
  # right join
  data3 <- merge(data1, data2, 
                 by = c(id, randdt, "adt2", paramcd), 
                 all.y = TRUE)
  
  # LOCF
  data.table::setorderv(data3, c(id, paramcd, "adt2"))
  data3[, `:=`(temp_col = data.table::nafill(get(aval), type = "locf")), 
        by = c(id, paramcd)]
  data3[[aval]] <- NULL
  data.table::setnames(data3, "temp_col", aval)
  
  # wide format
  fml <- paste(paste(c(id, randdt, "adt2"), collapse = " + "), "~", paramcd)
  data4 <- data.table::dcast(data3, fml, value.var = aval)
  
  # de-dup
  if (nodup) {
    # only keep rows where at least one paramcd changes value
    data4[, `:=`(change = rowSums(
      .SD != data.table::shift(.SD, type = "lag"), na.rm = TRUE) > 0),
      by = id, .SDcols = pars]
    
    # keep only change rows (plus the very first row per subject 
    # to establish baseline)
    data4 <- data4[get("change") | !duplicated(get(id))]
    data4[, `:=`(change = NULL)]   # drop helper column if not needed
  }
  
  # merge survival info
  data5 <- merge(data4, adsl, by = c(id, randdt))
  data5[, `:=`(ady = as.integer(difftime(get("adt2"), get(randdt), 
                                         units = "days") + offset), 
               osdy = as.integer(difftime(get(osdt), get(randdt), 
                                          units = "days")) + offset)]
  
  # create tstart and tstop relative to randdt
  # time-dependent covariate values at tstart
  # last interval ends with osdt     
  data5[, `:=`(tstart = get("ady"))]
  data5[, `:=`(tstop = data.table::shift(get("ady"), type = "lead")), 
        by = id]
  data5[, `:=`(tstop = data.table::fifelse(is.na(get("tstop")), 
                                           get("osdy"), get("tstop")))]
  
  # create event
  data5[, `:=`(event = data.table::fifelse(
    seq_len(.N) == .N, get(died), 0L)), by = id]
  
  # set up pd, swtrt, and censoring time
  data5[, `:=`(pd = !is.na(get(pddt)), 
               pd_time = as.integer(difftime(
                 get(pddt), get(randdt), units = "days")) + offset,
               swtrt = !is.na(get(xodt)),
               swtrt_time = as.integer(difftime(
                 get(xodt), get(randdt), units = "days")) + offset, 
               censor_time = as.integer(difftime(
                 get(dcutdt), get(randdt), units = "days")) + offset)]
  
  data5
}
