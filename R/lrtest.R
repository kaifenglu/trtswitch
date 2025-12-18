#' @title Log-Rank Test of Survival Curve Difference
#' @description Obtains the log-rank test using the Fleming-Harrington
#' family of weights.
#'
#' @param data The input data frame or list of data frames that contains 
#' the following variables:
#'
#'   * \code{stratum}: The stratum.
#'
#'   * \code{treat}: The treatment.
#'
#'   * \code{time}: The follow-up time for right censored data, or
#'     the left end of each interval for counting process data.
#'
#'   * \code{time2}: The right end of each interval for counting process
#'     data. Intervals are assumed to be open on the left
#'     and closed on the right, and event indicates whether an event
#'     occurred at the right end of each interval.
#'
#'   * \code{event}: The event indicator, 1=event, 0=no event.
#'
#'   * \code{weight}: The weight for each observation.
#'
#' @param stratum The name(s) of the stratum variable(s) in the input data.
#' @param treat The name of the treatment variable in the input data.
#' @param time The name of the time variable or the left end of each
#'   interval for counting process data in the input data.
#' @param time2 The name of the right end of each interval for counting
#'   process data in the input data.
#' @param event The name of the event variable in the input data.
#' @param weight The name of the weight variable in the input data.
#' @param weight_readj Whether the weight variable at each event time
#'   will be readjusted to be proportional to the number at risk by
#'   treatment group. Defaults to `FALSE`.
#' @param rho1 The first parameter of the Fleming-Harrington family of
#'   weighted log-rank test. Defaults to 0 for conventional log-rank test.
#' @param rho2 The second parameter of the Fleming-Harrington family of
#'   weighted log-rank test. Defaults to 0 for conventional log-rank test.
#' @param nthreads The number of threads to use in the computation (0 means 
#'   the default RcppParallel behavior)
#'   
#' @return A data frame (or list of data frames if the input is a list 
#' of data frames) with the following variables:
#'
#' * \code{uscore}: The numerator of the log-rank test statistic.
#'
#' * \code{vscore}: The variance of the log-rank score test statistic.
#'
#' * \code{logRankZ}: The Z-statistic value.
#'
#' * \code{logRankPValue}: The two-sided p-value.
#'
#' * \code{weight_readj}: Whether the weight variable will be readjusted.
#'
#' * \code{rho1}: The first parameter of the Fleming-Harrington weights.
#'
#' * \code{rho2}: The second parameter of the Fleming-Harrington weights.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' library(data.table)
#' dt <- as.data.table(rawdata)
#' data <- split(dt, by = "iterationNumber", keep.by = TRUE, flatten = FALSE)
#'
#' df <- lrtest(data, stratum = "stratum", treat = "treatmentGroup", 
#'              time = "timeUnderObservation", event = "event", 
#'              rho1 = 0.5, rho2 = 0)
#'              
#' df_all <- rbindlist(df)
#' df_all
#'
#' @export
lrtest <- function(data, stratum = "", treat = "", time = "time", 
                   time2 = "", event = "event", weight = "",  
                   weight_readj = FALSE, rho1 = 0, rho2 = 0, nthreads = 0) {
  
  # Validate data: must be a data.frame or a list of data.frames
  if (!inherits(data, "data.frame") && !is.list(data)) {
    stop("`data` must be either a data.frame or a list of data.frames.")
  }
  if (is.list(data) && !inherits(data, "data.frame")) {
    # ensure every element in the list is a data.frame (supports tibbles etc.)
    ok <- vapply(data, function(x) inherits(x, "data.frame"), logical(1))
    if (!all(ok)) {
      stop("When 'data' is a list, every element must be a data.frame.")
    }
  }
  
  # Respect user-requested number of threads (best effort)
  if (nthreads > 0) {
    RcppParallel::setThreadOptions(min(nthreads, parallel::detectCores(logical = FALSE)))
  }
  
  lrtestRcpp(data, stratum = stratum, treat = treat, time = time, 
             time2 = time2, event = event, weight = weight, 
             weight_readj = weight_readj, rho1 = rho1, rho2 = rho2)
}
