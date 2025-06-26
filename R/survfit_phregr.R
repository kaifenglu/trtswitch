#' @title Survival Curve for Proportional Hazards Regression Models
#' @description Obtains the predicted survivor function for a proportional
#' hazards regression model.
#'
#' @param object The output from the \code{phregr} call.
#' @param newdata A data frame with the same variable names as those that
#'   appear in the \code{phregr} call. For right-censored data, one curve is
#'   produced per row to represent a cohort whose covariates correspond to
#'   the values in \code{newdata}. For counting-process data, one curve is
#'   produced per \code{id} in \code{newdata} to present the survival curve
#'   along the path of time-dependent covariates at the observed event
#'   times in the data used to fit \code{phregr}.
#' @param sefit Whether to compute the standard error of the survival
#'   estimates.
#' @param conftype The type of the confidence interval. One of \code{"none"},
#'   \code{"plain"}, \code{"log"}, \code{"log-log"} (the default), or
#'   \code{"arcsin"}. The \code{arcsin} option bases the intervals on
#'   \code{asin(sqrt(surv))}.
#' @param conflev The level of the two-sided confidence interval for
#'   the survival probabilities. Defaults to 0.95.
#'
#' @details
#' If \code{newdata} is not provided and there is no covariate, survival
#' curves based on the \code{basehaz} data frame will be produced.
#'
#' @return A data frame with the following variables:
#'
#' * \code{id}: The id of the subject for counting-process data with
#'   time-dependent covariates.
#'
#' * \code{time}: The observed times in the data used to fit
#'   \code{phregr}.
#'
#' * \code{nrisk}: The number of patients at risk at the time point in the
#'   data used to fit \code{phregr}.
#'
#' * \code{nevent}: The number of patients having event at the time point
#'   in the data used to fit \code{phregr}.
#'
#' * \code{cumhaz}: The cumulative hazard at the time point.
#'
#' * \code{surv}: The estimated survival probability at the time point.
#'
#' * \code{sesurv}: The standard error of the estimated survival probability.
#'
#' * \code{lower}: The lower confidence limit for survival probability.
#'
#' * \code{upper}: The upper confidence limit for survival probability.
#'
#' * \code{conflev}: The level of the two-sided confidence interval.
#'
#' * \code{conftype}: The type of the confidence interval.
#'
#' * \code{covariates}: The values of covariates based on \code{newdata}.
#'
#' * \code{stratum}: The stratum of the subject.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @references
#' Terry M. Therneau and Patricia M. Grambsch.
#' Modeling Survival Data: Extending the Cox Model.
#' Springer-Verlag, 2000.
#'
#' @examples
#'
#' library(dplyr)
#'
#' # Example 1 with right-censored data
#' fit1 <- phregr(data = rawdata %>% filter(iterationNumber == 1) %>%
#'                  mutate(treat = 1*(treatmentGroup == 1)),
#'                stratum = "stratum",
#'                time = "timeUnderObservation", event = "event",
#'                covariates = "treat")
#'
#' surv1 <- survfit_phregr(fit1,
#'                         newdata = data.frame(
#'                           stratum = as.integer(c(1,1,2,2)),
#'                           treat = c(1,0,1,0)))
#'
#' # Example 2 with counting process data and robust variance estimate
#' fit2 <- phregr(data = heart %>% mutate(rx = as.numeric(transplant) - 1),
#'                time = "start", time2 = "stop", event = "event",
#'                covariates = c("rx", "age"), id = "id", robust = TRUE)
#'
#' surv2 <- survfit_phregr(fit2,
#'                         newdata = data.frame(
#'                           id = c(4,4,11,11),
#'                           age = c(-7.737,-7.737,-0.019,-0.019),
#'                           start = c(0,36,0,26),
#'                           stop = c(36,39,26,153),
#'                           rx = c(0,1,0,1)))
#'
#' @export
survfit_phregr <- function(object, newdata, sefit = TRUE,
                           conftype = "log-log", conflev = 0.95) {

  p = object$p
  if (p == 0) {
    beta = 0
    vbeta = 0
  } else {
    beta = object$beta
    vbeta = object$vbeta
  }

  basehaz = object$basehaz

  covariates = object$covariates
  stratum = object$stratum
  offset = object$offset
  id = object$id

  if (id != "") {
    tstart = object$time
    tstop = object$time2
  } else {
    tstart = ""
    tstop = ""
  }

  nvar = length(covariates)
  if (missing(covariates) || is.null(covariates) || (nvar == 1 && (
    covariates[1] == "" || tolower(covariates[1]) == "none"))) {
    p3 = 0
  } else {
    fml1 = formula(paste("~", paste(covariates, collapse = "+")))
    p3 = length(rownames(attr(terms(fml1), "factors")))
  }

  if (p >= 1 && p3 >= 1 && !(missing(newdata) || is.null(newdata))) {
    df = newdata
    mf1 <- model.frame(fml1, data = df, na.action = na.pass, 
                       xlev = object$xlevels)
    mm <- model.matrix(fml1, mf1)
    colnames(mm) = make.names(colnames(mm))
    varnames = colnames(mm)[-1]
    for (i in 1:length(varnames)) {
      if (!(varnames[i] %in% names(df))) {
        df[,varnames[i]] = mm[,varnames[i]]
      }
    }
  } else {
    beta = NA
    vbeta = NA
    varnames = ""
  }

  if (!is.matrix(vbeta)) vbeta = as.matrix(vbeta)

  if (missing(basehaz) || is.null(basehaz)) {
    stop("basehaz must be provided")
  }

  if (missing(newdata) || is.null(newdata)) {
    if (p > 0) {
      stop("newdata must be provided for Cox models with covariates")
    } else {
      p_stratum = length(stratum);
      if (p_stratum == 1 && (stratum[1] == "" || stratum[1] == "none")) {
        df = data.frame(dummy_x_ = 0)
      } else {
        df = unique(basehaz[, stratum, drop = FALSE])
      }
    }
  }

  survfit_phregcpp(p = p, beta = beta, vbeta = vbeta, basehaz = basehaz,
                   newdata = df, covariates = varnames,
                   stratum = stratum, offset = offset, id = id,
                   tstart = tstart, tstop = tstop, sefit = sefit,
                   conftype = conftype, conflev = conflev)
}
