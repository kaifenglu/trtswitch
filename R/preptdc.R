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
#' @param adsl A data set containing baseline subject-level information. 
#'   It should include, at a minimum, subject ID (`id`), 
#'   randomization date (`randdt`), treatment start date (`trtsdt`),
#'   progression date (`pddt`), treatment switch date (`xodt`), 
#'   survival outcome (`osdt`, `died`), and data cut-off date (`dcutdt`).
#' @param adtdc A data set containing longitudinal
#'   time-dependent covariate data, with subject ID (`id`),
#'   parameter code (`paramcd`), analysis date (`adt`), and covariate
#'   value (`aval`).
#' @param id Character string specifying the column name for subject ID.
#' @param randdt Character string specifying the column name for 
#'   randomization date.
#' @param trtsdt Character string specifying the column name for 
#'   treatment start date.
#' @param pddt Character string specifying the column name for progression
#'   date.
#' @param xodt Character string specifying the column name for treatment
#'   crossover/switch date.
#' @param osdt Character string specifying the column name for overall
#'   survival date (death date or last known alive date).
#' @param died Character string specifying the column name for death
#'   indicator (0 = alive/censored, 1 = died).
#' @param dcutdt Character string specifying the column name for data
#'   cut-off date.
#' @param paramcd Character string specifying the column name for parameter
#'   code (identifying different covariates).
#' @param adt Character string specifying the column name for analysis
#'   date in the time-dependent covariate dataset.
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
#'   \item Merge `adsl` and `adtdc` to obtain randomization date and 
#'         treatment start date.
#'   \item Define `adt2` as `adt` if `adt > trtsdt`, 
#'         and `randdt` if `adt <= trtsdt` (i.e., baseline time point).
#'         This ensures that the baseline covariate value is the last 
#'         non-missing value at or before the treatment start date.
#'         Post-baseline covariate values are anchored at their actual
#'         analysis dates. The first record per subject corresponds to 
#'         survival time zero at randomization and ensures availability 
#'         of baseline covariates at randomization.
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
                    trtsdt = "TRTSDT", pddt = "PDDT", xodt = "XODT", 
                    osdt = "OSDT", died = "DIED", dcutdt = "DCUTDT", 
                    paramcd = "PARAMCD", adt = "ADT", aval = "AVAL", 
                    nodup = TRUE, offset = TRUE) {
  
  # convert input data into data.table
  data.table::setDT(adsl)
  data.table::setDT(adtdc)
  
  # verify whether the input data have required columns
  cols <- colnames(adsl)
  req_cols <- c(id, randdt, trtsdt, pddt, xodt, osdt, died, dcutdt)
  if (!all(req_cols %in% cols)) {
    stop(paste("The following columns are missing from adsl:", 
               paste(req_cols[!(req_cols %in% cols)], collapse = ", ")))
  }
  
  cols <- colnames(adtdc)
  req_cols <- c(id, paramcd, adt, aval)
  if (!all(req_cols %in% cols)) {
    stop(paste("The following columns are missing from adtdc:", 
               paste(req_cols[!(req_cols %in% cols)], collapse = ", ")))
  }
  
  # obtain randdt and trtsdt for each record in adtdc
  cols <- c(id, randdt, trtsdt)
  data1 <- merge(adtdc, adsl[, mget(cols)], by = id)
  
  # define adt2 to differentiate baseline and post-baseline time points
  data1[, `:=`(adt2 = data.table::fifelse(get(adt) <= get(trtsdt), 
                                          get(randdt), get(adt)))]
  
  # keep last record with each subject, adt2, and paramcd
  data.table::setorderv(data1, c(id, "adt2", paramcd, adt))
  data1 <- data1[, .SD[.N], by = c(id, "adt2", paramcd)]
  
  # build skeleton for id, randdt, trtsdt, adt2, and paramcd
  cols <- c(id, randdt, trtsdt, "adt2")
  pars <- unique(data1[[paramcd]])
  data2 <- merge(
    unique(data1[, mget(cols)])[, `:=`(dummy = 1)],
    data.table::data.table(temp_col = pars, dummy = 1),
    by = "dummy", allow.cartesian = TRUE)
  data.table::setnames(data2, "temp_col", paramcd)
  data2[, `:=`(dummy = NULL)]
  
  # right join
  data3 <- merge(data1, data2, 
                 by = c(id, randdt, trtsdt, "adt2", paramcd), 
                 all.y = TRUE)
  
  data.table::setorderv(data3, c(id, paramcd, "adt2"))
  
  
  data_list <- split(data3, by = c(id, paramcd))
  data_list <- lapply(data_list, function(sub) {
    # Perform assignment in standard R context
    sub[[aval]] = data.table::nafill(sub[[aval]], type = "locf")
    return(sub)
  })
  data3 <- data.table::rbindlist(data_list)
  
  # wide format
  fml <- paste(paste(c(id, randdt, "adt2"), collapse = " + "), "~", paramcd)
  data4 <- data.table::dcast(data3, fml, value.var = aval)
  
  # de-dup
  if (nodup) {
    # keep only change rows (plus the very first row per subject for baseline)
    idx <- data.table::rleidv(data4, c(id, pars))
    idxprev <- data.table::shift(idx, type = "lag", fill = 0L)
    data4 <- data4[which(idx != idxprev), ]
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
  
  data_list <- split(data5, by = id)
  data_list <- lapply(data_list, function(sub) {
    # Perform assignment in standard R context
    tstop <- data.table::shift(sub$ady, type = "lead")
    sub$tstop <- data.table::fifelse(is.na(tstop), sub$osdy, tstop)
    return(sub)
  })
  data5 <- data.table::rbindlist(data_list)
  
  # create event
  last_rows <- data5[, .I[.N], by = id]$V1
  data5[, event := 0L]
  data5[last_rows, event := as.integer(get(died))]  
  
  # set up pd, swtrt, and censoring time
  data5[, `:=`(pd = !is.na(get(pddt)), 
               pd_time = as.integer(difftime(
                 get(pddt), get(randdt), units = "days")) + offset,
               swtrt = !is.na(get(xodt)),
               swtrt_time = as.integer(difftime(
                 get(xodt), get(randdt), units = "days")) + offset, 
               censor_time = as.integer(difftime(
                 get(dcutdt), get(randdt), units = "days")) + offset)]
  
  as.data.frame(data5)
}
