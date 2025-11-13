#' @title Parametric Regression Models for Failure Time Data
#' @description Obtains the parameter estimates from parametric
#' regression models with uncensored, right censored, left censored, or
#' interval censored data.
#'
#' @param data The input data frame that contains the following variables:
#'
#'   * \code{rep}: The replication for by-group processing.
#'
#'   * \code{stratum}: The stratum.
#'
#'   * \code{time}: The follow-up time for right censored data, or
#'     the left end of each interval for interval censored data.
#'
#'   * \code{time2}: The right end of each interval for interval
#'     censored data.
#'
#'   * \code{event}: The event indicator, 1=event, 0=no event.
#'
#'   * \code{covariates}: The values of baseline covariates.
#'
#'   * \code{weight}: The weight for each observation.
#'
#'   * \code{offset}: The offset for each observation.
#'
#'   * \code{id}: The optional subject ID to group the score residuals
#'     in computing the robust sandwich variance.
#'
#' @param rep The name(s) of the replication variable(s) in the input data.
#' @param stratum The name(s) of the stratum variable(s) in the input data.
#' @param time The name of the time variable or the left end of each
#'   interval for interval censored data in the input data.
#' @param time2 The name of the right end of each interval for
#'   interval censored data in the input data.
#' @param event The name of the event variable in the input data
#'   for right censored data.
#' @param covariates The vector of names of baseline covariates
#'   in the input data.
#' @param weight The name of the weight variable in the input data.
#' @param offset The name of the offset variable in the input data.
#' @param id The name of the id variable in the input data.
#' @param dist The assumed distribution for time to event. Options include
#'   "exponential", "weibull", "lognormal", and "loglogistic" to be
#'   modeled on the log-scale, and "normal" and "logistic" to be modeled
#'   on the original scale.
#' @param init A vector of initial values for the model parameters, 
#'   including regression coefficients and the log scale parameter. 
#'   By default, initial values are derived from an intercept-only model. 
#'   If this approach fails, ordinary least squares (OLS) estimates, 
#'   ignoring censoring, are used instead.
#' @param robust Whether a robust sandwich variance estimate should be
#'   computed. In the presence of the id variable, the score residuals
#'   will be aggregated for each id when computing the robust sandwich
#'   variance estimate.
#' @param plci Whether to obtain profile likelihood confidence interval.
#' @param alpha The two-sided significance level.
#' @param maxiter The maximum number of iterations.
#' @param eps The tolerance to declare convergence. 
#'
#' @details There are two ways to specify the model, one for right censored
#' data through the time and event variables, and the other for interval
#' censored data through the time (lower) and time2 (upper) variables.
#' For the second form, we follow the convention used in SAS PROC LIFEREG:
#'
#' * If lower is not missing, upper is not missing, and lower is equal
#'   to upper, then there is no censoring and the event occurred at
#'   time lower.
#'
#' * If lower is not missing, upper is not missing, and lower < upper,
#'   then the event time is censored within the interval (lower, upper).
#'
#' * If lower is missing, but upper is not missing, then upper will be
#'   used as the left censoring value.
#'
#' * If lower is not missing, but upper is missing, then lower will be
#'   used as the right censoring value.
#'
#' * If lower is not missing, upper is not missing, but lower > upper,
#'   or if both lower and upper are missing, then the observation will
#'   not be used.
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
#'     - \code{loglik0}: The log-likelihood under null.
#'
#'     - \code{loglik1}: The maximum log-likelihood.
#'
#'     - \code{niter}: The number of Newton-Raphson iterations.
#'
#'     - \code{dist}: The assumed distribution.
#'
#'     - \code{p}: The number of parameters, including the intercept,
#'       regression coefficients associated with the covariates, and
#'       the log scale parameters for the strata.
#'
#'     - \code{nvar}: The number of regression coefficients associated
#'       with the covariates (excluding the intercept).
#'
#'     - \code{robust}: Whether the robust sandwich variance estimate
#'       is requested.
#'       
#'     - \code{fail}: Whether the model fails to converge.
#'
#'     - \code{rep}: The replication.
#'
#' * \code{parest}: The data frame of parameter estimates with the
#'   following variables:
#'
#'     - \code{param}: The name of the covariate for the parameter estimate.
#'
#'     - \code{beta}: The parameter estimate.
#'
#'     - \code{sebeta}: The standard error of parameter estimate.
#'
#'     - \code{z}: The Wald test statistic for the parameter.
#'
#'     - \code{expbeta}: The exponentiated parameter estimate.
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
#'     - \code{sebeta_naive}: The naive standard error of parameter estimate
#'       if robust variance is requested.
#'
#'     - \code{vbeta_naive}: The naive covariance matrix for parameter
#'       estimates if robust variance is requested.
#'
#'     - \code{rep}: The replication.
#'
#' * \code{p}: The number of parameters.
#'
#' * \code{nvar}: The number of columns of the design matrix excluding
#'   the intercept.
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
#' John D. Kalbfleisch and Ross L. Prentice.
#' The Statistical Analysis of Failure Time Data.
#' Wiley: New York, 1980.
#'
#' @examples
#'
#' library(dplyr)
#'
#' # right censored data
#' (fit1 <- liferegr(
#'   data = rawdata %>% mutate(treat = 1*(treatmentGroup == 1)),
#'   rep = "iterationNumber", stratum = "stratum",
#'   time = "timeUnderObservation", event = "event",
#'   covariates = "treat", dist = "weibull"))
#'
#' # tobit regression for left censored data
#' (fit2 <- liferegr(
#'   data = tobin %>% mutate(time = ifelse(durable>0, durable, NA)),
#'   time = "time", time2 = "durable",
#'   covariates = c("age", "quant"), dist = "normal"))
#'
#' @export
liferegr <- function(data, rep = "", stratum = "",
                     time = "time", time2 = "", event = "event",
                     covariates = "", weight = "", offset = "",
                     id = "", dist = "weibull", init = NA_real_, 
                     robust = FALSE, plci = FALSE, alpha = 0.05, 
                     maxiter = 50, eps = 1.0e-9) {
  rownames(data) = NULL
  
  elements = c(rep, stratum, covariates, weight, offset, id)
  elements = unique(elements[elements != "" & elements != "none"])
  if (!(length(elements) == 0)) {
    fml = formula(paste("~", paste(elements, collapse = "+")))
    mf = model.frame(fml, data = data, na.action = na.omit)
  } else {
    mf = model.frame(formula("~1"), data = data)
  }
  
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
  
  fit <- liferegcpp(data = df, rep = rep, stratum = stratum, time = time,
                    time2 = time2, event = event, covariates = varnames,
                    weight = weight, offset = offset, id = id, dist = dist,
                    init = init, robust = robust, plci = plci, 
                    alpha = alpha, maxiter = maxiter, eps = eps)
  
  fit$p <- fit$sumstat$p[1]
  fit$nvar <- fit$sumstat$nvar[1]
  
  if (fit$p > 0) {
    par = fit$parest$param[1:fit$p]
    if (length(par) > length(param)) {
      fit$param = c(param, par[(1+length(param)):length(par)])
    } else {
      fit$param = param
    }
    
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
    data = data, rep = rep, stratum = stratum, time = time, time2 = time2,
    event = event, covariates = covariates, weight = weight, offset = offset,
    id = id, dist = dist, init = init, robust = robust, plci = plci, 
    alpha = alpha, maxiter = maxiter, eps = eps
  )
  
  class(fit) = "liferegr"
  fit
}
