#' @title Estimate of Milestone Survival Difference
#' @description Obtains the estimate of milestone survival difference
#' between two treatment groups.
#'
#' @param data The input data frame that contains the following variables:
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
#' @param stratum The name of the stratum variable in the input data.
#' @param treat The name of the treatment variable in the input data.
#' @param time The name of the time variable or the left end of each
#'   interval for counting process data in the input data.
#' @param time2 The name of the right end of each interval for counting
#'   process data in the input data.
#' @param event The name of the event variable in the input data.
#' @param weight The name of the weight variable in the input data.
#' @param milestone The milestone time at which to calculate the
#'   survival probability.
#' @param survDiffH0 The difference in milestone survival probabilities
#'   under the null hypothesis. Defaults to 0 for superiority test.
#' @param conflev The level of the two-sided confidence interval for
#'   the difference in milestone survival probabilities. Defaults to 0.95.
#' @param nthreads The number of threads to use in the computation (0 means 
#'   the default RcppParallel behavior)   
#'
#' @return A data frame (or list of data frames if the input is a list 
#' of data frames) with the following variables:
#'
#' * \code{milestone}: The milestone time relative to randomization.
#'
#' * \code{survDiffH0}: The difference in milestone survival probabilities
#'   under the null hypothesis.
#'
#' * \code{surv1}: The estimated milestone survival probability for
#'   the treatment group.
#'
#' * \code{surv2}: The estimated milestone survival probability for
#'   the control group.
#'
#' * \code{survDiff}: The estimated difference in milestone survival
#'   probabilities.
#'
#' * \code{vsurv1}: The variance for surv1.
#'
#' * \code{vsurv2}: The variance for surv2.
#'
#' * \code{sesurvDiff}: The standard error for survDiff.
#'
#' * \code{survDiffZ}: The Z-statistic value.
#'
#' * \code{survDiffPValue}: The two-sided p-value.
#'
#' * \code{lower}: The lower bound of confidence interval.
#'
#' * \code{upper}: The upper bound of confidence interval.
#'
#' * \code{conflev}: The level of confidence interval.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' df <- kmdiff(data = rawdata, stratum = "stratum", 
#'              treat = "treatmentGroup", time = "timeUnderObservation", 
#'              event = "event", milestone = 12)
#' head(df)
#'
#' @export
kmdiff <- function(data, stratum = "", treat = "", time = "time", 
                   time2 = "", event = "event", weight = "",  
                   milestone, survDiffH0 = 0, conflev = 0.95, nthreads = 0) {
  
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
  
  kmdiffRcpp(data, stratum = stratum, treat = treat, time = time, 
             time2 = time2, event = event, weight = weight, 
             milestone = milestone, survDiffH0 = survDiffH0,
             conflev = conflev)
}
