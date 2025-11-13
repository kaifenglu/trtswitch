#' @title Proportional Hazards Regression Models
#' @description Obtains the hazard ratio estimates from the proportional
#' hazards regression model with right censored or counting process data.
#'
#' @param data The input data frame that contains the following variables:
#'
#'   * \code{rep}: The replication for by-group processing.
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
#'   * \code{covariates}: The values of baseline covariates (and
#'     time-dependent covariates in each interval for counting
#'     process data).
#'
#'   * \code{weight}: The weight for each observation.
#'
#'   * \code{offset}: The offset for each observation.
#'
#'   * \code{id}: The optional subject ID for counting process data
#'     with time-dependent covariates.
#'
#' @param rep The name(s) of the replication variable(s) in the input data.
#' @param stratum The name(s) of the stratum variable(s) in the input data.
#' @param time The name of the time variable or the left end of each
#'   interval for counting process data in the input data.
#' @param time2 The name of the right end of each interval for counting
#'   process data in the input data.
#' @param event The name of the event variable in the input data.
#' @param covariates The vector of names of baseline and time-dependent
#'   covariates in the input data.
#' @param weight The name of the weight variable in the input data.
#' @param offset The name of the offset variable in the input data.
#' @param id The name of the id variable in the input data.
#' @param ties The method for handling ties, either "breslow" or 
#'   "efron" (default).
#' @param init The vector of initial values. Defaults to zero for all 
#'   variables.
#' @param robust Whether a robust sandwich variance estimate should be
#'   computed. In the presence of the id variable, the score residuals
#'   will be aggregated for each id when computing the robust sandwich
#'   variance estimate.
#' @param est_basehaz Whether to estimate the baseline hazards.
#'   Defaults to \code{TRUE}.
#' @param est_resid Whether to estimate the martingale residuals.
#'   Defaults to \code{TRUE}.
#' @param firth Whether to use Firth’s penalized likelihood method.
#'   Defaults to \code{FALSE}.
#' @param plci Whether to obtain profile likelihood confidence interval.
#' @param alpha The two-sided significance level.
#' @param maxiter The maximum number of iterations.
#' @param eps The tolerance to declare convergence. 
#'
#' @return A list with the following components:
#'
#' * \code{sumstat}: The data frame of summary statistics of model fit
#'   with the following variables:
#'
#'     - \code{n}: The number of observations.
#'
#'     - \code{nevents}: The number of events.
#'
#'     - \code{loglik0}: The (penalized) log-likelihood under null.
#'
#'     - \code{loglik1}: The maximum (penalized) log-likelihood.
#'
#'     - \code{scoretest}: The score test statistic.
#'
#'     - \code{niter}: The number of Newton-Raphson iterations.
#'
#'     - \code{ties}: The method for handling ties, either "breslow" or
#'       "efron".
#'
#'     - \code{p}: The number of columns of the Cox model design matrix.
#'
#'     - \code{robust}: Whether to use the robust variance estimate.
#'
#'     - \code{firth}: Whether to use Firth's penalized likelihood method.
#'     
#'     - \code{fail}: Whether the model fails to converge.
#'
#'     - \code{loglik0_unpenalized}: The unpenalized log-likelihood under null.
#'
#'     - \code{loglik1_unpenalized}: The maximum unpenalized log-likelihood.
#'
#'     - \code{rep}: The replication.
#'
#' * \code{parest}: The data frame of parameter estimates with the
#'   following variables:
#'
#'     - \code{param}: The name of the covariate for the parameter estimate.
#'
#'     - \code{beta}: The log hazard ratio estimate.
#'
#'     - \code{sebeta}: The standard error of log hazard ratio estimate.
#'
#'     - \code{z}: The Wald test statistic for log hazard ratio.
#'
#'     - \code{expbeta}: The hazard ratio estimate.
#'
#'     - \code{vbeta}: The covariance matrix for parameter estimates.
#'
#'     - \code{lower}: The lower limit of confidence interval.
#'
#'     - \code{upper}: The upper limit of confidence interval.
#'
#'     - \code{p}: The p-value from the chi-square test.
#'
#'     - \code{method}: The method to compute the confidence interval and
#'       p-value.
#'
#'     - \code{sebeta_naive}: The naive standard error of log hazard ratio
#'       estimate if robust variance is requested.
#'
#'     - \code{vbeta_naive}: The naive covariance matrix for parameter
#'       estimates if robust variance is requested.
#'
#'     - \code{rep}: The replication.
#'
#' * \code{basehaz}: The data frame of baseline hazards with the following
#'   variables (if est_basehaz is TRUE):
#'
#'     - \code{time}: The observed event time.
#'
#'     - \code{nrisk}: The number of patients at risk at the time point.
#'
#'     - \code{nevent}: The number of events at the time point.
#'
#'     - \code{haz}: The baseline hazard at the time point.
#'
#'     - \code{varhaz}: The variance of the baseline hazard at the time point
#'       assuming the parameter beta is known.
#'
#'     - \code{gradhaz}: The gradient of the baseline hazard with respect to
#'       beta at the time point (in the presence of covariates).
#'
#'     - \code{stratum}: The stratum.
#'
#'     - \code{rep}: The replication.
#'
#' * \code{residuals}: The martingale residuals.
#'
#' * \code{p}: The number of parameters.
#'
#' * \code{param}: The parameter names.
#'
#' * \code{beta}: The parameter estimate.
#'
#' * \code{vbeta}: The covariance matrix for parameter estimates.
#'
#' * \code{vbeta_naive}: The naive covariance matrix for parameter estimates.
#'
#' * \code{terms}: The terms object.
#'
#' * \code{xlevels}: A record of the levels of the factors used in fitting.
#' 
#' * \code{settings}: A list containing the input parameter values.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @references
#' Per K. Anderson and Richard D. Gill.
#' Cox's regression model for counting processes, a large sample study.
#' Annals of Statistics 1982; 10:1100-1120.
#'
#' Terry M. Therneau and Patricia M. Grambsch.
#' Modeling Survival Data: Extending the Cox Model.
#' Springer-Verlag, 2000.
#'
#' @examples
#'
#' library(dplyr)
#'
#' # Example 1 with right-censored data
#' (fit1 <- phregr(
#'   data = rawdata %>% mutate(treat = 1*(treatmentGroup == 1)),
#'   rep = "iterationNumber", stratum = "stratum",
#'   time = "timeUnderObservation", event = "event",
#'   covariates = "treat", est_basehaz = FALSE, est_resid = FALSE))
#'
#' # Example 2 with counting process data and robust variance estimate
#' (fit2 <- phregr(
#'   data = heart %>% mutate(rx = as.numeric(transplant) - 1),
#'   time = "start", time2 = "stop", event = "event",
#'   covariates = c("rx", "age"), id = "id",
#'   robust = TRUE, est_basehaz = TRUE, est_resid = TRUE))
#'
#' @export
phregr <- function(data, rep = "", stratum = "",
                   time = "time", time2 = "", event = "event",
                   covariates = "", weight = "", offset = "",
                   id = "", ties = "efron", 
                   init = NA_real_,  robust = FALSE,
                   est_basehaz = TRUE, est_resid = TRUE, firth = FALSE,
                   plci = FALSE, alpha = 0.05, 
                   maxiter = 50, eps = 1.0e-9) {
  
  rownames(data) = NULL
  
  elements = c(rep, stratum, time, event, covariates, weight, offset, id)
  elements = unique(elements[elements != "" & elements != "none"])
  fml = formula(paste("~", paste(elements, collapse = "+")))
  mf = model.frame(fml, data = data, na.action = na.omit)
  
  rownum = as.integer(rownames(mf))
  df = data[rownum,]
  
  nvar = length(covariates)
  if (missing(covariates) || is.null(covariates) || (nvar == 1 && (
    covariates[1] == "" || tolower(covariates[1]) == "none"))) {
    p = 0
    t1 = terms(formula("~1"))
  } else {
    fml1 = formula(paste("~", paste(covariates, collapse = "+")))
    p = length(rownames(attr(terms(fml1), "factors")))
    t1 = terms(fml1)
  }
  
  if (p >= 1) {
    mf1 <- model.frame(fml1, data = df, na.action = na.pass)
    mm <- model.matrix(fml1, mf1)
    xlevels = mf1$xlev
    param = colnames(mm)
    colnames(mm) = make.names(colnames(mm))
    varnames = colnames(mm)[-1]
    for (i in 1:length(varnames)) {
      if (!(varnames[i] %in% names(df))) {
        df[,varnames[i]] = mm[,varnames[i]]
      }
    }
  } else {
    xlevels = NULL
    param = "(Intercept)"
    varnames = ""
  }
  
  fit <- phregcpp(data = df, rep = rep, stratum = stratum, time = time,
                  time2 = time2, event = event, covariates = varnames,
                  weight = weight, offset = offset, id = id, 
                  ties = ties, init = init, 
                  robust = robust, est_basehaz = est_basehaz,
                  est_resid = est_resid, firth = firth,
                  plci = plci, alpha = alpha, 
                  maxiter = maxiter, eps = eps)
  
  fit$p = fit$sumstat$p[1]
  
  if (fit$p > 0) {
    fit$param = param[-1]
    fit$beta = fit$parest$beta
    names(fit$beta) = rep(fit$param, length(fit$beta)/fit$p)
    
    if (fit$p > 1) {
      fit$vbeta = as.matrix(fit$parest[, paste0("vbeta.", 1:fit$p)])
      if (robust) {
        fit$vbeta_naive = as.matrix(fit$parest[, paste0("vbeta_naive.",
                                                        1:fit$p)])
      }
    } else {
      fit$vbeta = as.matrix(fit$parest[, "vbeta"])
      if (robust) {
        fit$vbeta_naive = as.matrix(fit$parest[, "vbeta_naive"])
      }
    }
    
    dimnames(fit$vbeta) = list(names(fit$beta), fit$param)
    if (robust) {
      dimnames(fit$vbeta_naive) = list(names(fit$beta), fit$param)
    }
  }
  
  fit$terms = t1
  if (fit$p > 0) fit$xlevels = xlevels
  
  fit$settings <- list(
    data = data, rep = rep, stratum = stratum, time = time, 
    time2 = time2, event = event, covariates = covariates,
    weight = weight, offset = offset, id = id, ties = ties,
    iniy = init, robust = robust, est_basehaz = est_basehaz,
    est_resid = est_resid, firth = firth, plci = plci,
    alpha = alpha, maxiter = maxiter, eps = eps
  )
  
  class(fit) = "phregr"
  fit
}
