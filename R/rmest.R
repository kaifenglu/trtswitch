#' @title Estimate of Restricted Mean Survival Time
#' @description Obtains the estimate of restricted means survival time
#' for each stratum.
#'
#' @param data The input data frame or list of data frames that contains 
#' the following variables:
#'
#'   * \code{stratum}: The stratum.
#'
#'   * \code{time}: The possibly right-censored survival time.
#'
#'   * \code{event}: The event indicator.
#'
#' @param stratum The name of the stratum variable in the input data.
#' @param time The name of the time variable in the input data.
#' @param event The name of the event variable in the input data.
#' @param milestone The milestone time at which to calculate the
#'   restricted mean survival time.
#' @param conflev The level of the two-sided confidence interval for
#'   the survival probabilities. Defaults to 0.95.
#' @param biascorrection Whether to apply bias correction for the
#'   variance estimate. Defaults to no bias correction.
#' @param nthreads The number of threads to use in the computation (0 means 
#'   the default RcppParallel behavior)
#'
#' @return A data frame (or list of data frames if the input is a list 
#' of data frames) with the following variables:
#'
#' * \code{stratum}: The stratum variable.
#'
#' * \code{size}: The number of subjects in the stratum.
#'
#' * \code{milestone}: The milestone time relative to randomization.
#'
#' * \code{rmst}: The estimate of restricted mean survival time.
#'
#' * \code{stderr}: The standard error of the estimated rmst.
#'
#' * \code{lower}: The lower bound of confidence interval if requested.
#'
#' * \code{upper}: The upper bound of confidence interval if requested.
#'
#' * \code{conflev}: The level of confidence interval if requested.
#'
#' * \code{biascorrection}: Whether to apply bias correction for the
#'   variance estimate.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' rmest(data = aml, stratum = "x",
#'       time = "time", event = "status", milestone = 24)
#'
#' @export
rmest <- function(data, stratum = "", time = "time", event = "event", 
                  milestone, conflev = 0.95, biascorrection = FALSE, 
                  nthreads = 0) {
  
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
  
  rmestRcpp(data, stratum = stratum, time = time, event = event, 
            milestone = milestone, conflev = conflev, 
            biascorrection = biascorrection)
}
