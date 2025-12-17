#' @title Kaplan-Meier Estimates of Survival Curve
#' @description Obtains the Kaplan-Meier estimates of the survival curve.
#'
#' @param data The input data frame or list of data frames that contains 
#' the following variables:
#'
#'   * \code{stratum}: The stratum.
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
#' @param time The name of the time variable or the left end of each
#'   interval for counting process data in the input data.
#' @param time2 The name of the right end of each interval for counting
#'   process data in the input data.
#' @param event The name of the event variable in the input data.
#' @param weight The name of the weight variable in the input data.
#' @param conftype The type of the confidence interval. One of "none",
#'   "plain", "log", "log-log" (the default), or "arcsin".
#'   The arcsin option bases the intervals on asin(sqrt(survival)).
#' @param conflev The level of the two-sided confidence interval for
#'   the survival probabilities. Defaults to 0.95.
#' @param keep_censor Whether to retain the censoring time in the output
#'   data frame.
#' @param nthreads The number of threads to use in the computation (0 means 
#'   the default RcppParallel behavior)
#'
#' @return A data frame (or list of data frames if the input is a list 
#' of data frames) with the following variables:
#'
#' * \code{size}: The number of subjects in the stratum.
#'
#' * \code{time}: The event time.
#'
#' * \code{nrisk}: The number of subjects at risk.
#'
#' * \code{nevent}: The number of subjects having the event.
#'
#' * \code{ncensor}: The number of censored subjects.
#'
#' * \code{surv}: The Kaplan-Meier estimate of the survival probability.
#'
#' * \code{sesurv}: The standard error of the estimated survival
#'   probability based on the Greendwood formula.
#'
#' * \code{lower}: The lower bound of confidence interval if requested.
#'
#' * \code{upper}: The upper bound of confidence interval if requested.
#'
#' * \code{conflev}: The level of confidence interval if requested.
#'
#' * \code{conftype}: The type of confidence interval if requested.
#'
#' * \code{stratum}: The stratum.
#'
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' kmest(data = aml, stratum = "x", time = "time", event = "status")
#'
#' @export
kmest <- function(data, stratum = "", time = "time", time2 = "", 
                  event = "event", weight = "", conftype = "log-log",
                  conflev = 0.95, keep_censor = FALSE, nthreads = 0) {
  
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
  
  kmestRcpp(data, stratum = stratum, time = time, time2 = time2,
            event = event, weight = weight, conftype = conftype,
            conflev = conflev, keep_censor = keep_censor)
}
