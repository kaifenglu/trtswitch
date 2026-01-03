#' @title Parametric Regression Models for Failure Time Data
#' @description Obtains the parameter estimates from parametric
#' regression models with uncensored, right censored, left censored, or
#' interval censored data.
#'
#' @param data The input data frame that contains the following variables:
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
#' * \code{linear_predictors}: The vector of linear predictors.
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
#'   data = rawdata %>% filter(iterationNumber == 1) %>% 
#'          mutate(treat = (treatmentGroup == 1)),
#'   stratum = "stratum",
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
liferegr <- function(data, stratum = "", time = "time", time2 = "", 
                     event = "event", covariates = "", weight = "", 
                     offset = "", id = "", dist = "weibull", 
                     init = NA_real_, robust = FALSE, plci = FALSE, 
                     alpha = 0.05, maxiter = 50, eps = 1.0e-9) {
  
  # validate input
  if (!inherits(data, "data.frame")) {
    stop("Input 'data' must be a data frame");
  }
  
  if (inherits(data, "data.table") || inherits(data, "tbl") || 
      inherits(data, "tbl_df")) {
    df <- as.data.frame(data)
  } else {
    df <- data
  }
  
  for (nm in c(time, time2, event, weight, offset, id)) {
    if (!is.character(nm) || length(nm) != 1) {
      stop(paste(nm, "must be a single character string."));
    }
  }
  
  # select complete cases for the relevant variables
  elements <- unique(c(stratum, covariates, weight, offset, id))
  elements <- elements[elements != ""]
  fml_all <- formula(paste("~", paste(elements, collapse = "+")))
  var_all <- all.vars(fml_all)
  
  # check if the input data contains the required columns
  missing_cols <- setdiff(var_all, names(df))
  if (length(missing_cols) > 0) {
    stop(paste0("The following required columns are missing in the input data: ",
                paste(missing_cols, collapse = ", ")))
  }
  
  # use complete.cases on the subset of columns we care about
  rows_ok <- which(complete.cases(df[, var_all, drop = FALSE]))
  if (length(rows_ok) == 0) stop("No complete cases found for the specified variables.")
  df <- df[rows_ok, , drop = FALSE]

  # Determine if covariates were provided (empty string or NULL means no covariates)
  misscovariates <- length(covariates) == 0 || 
    (length(covariates) == 1 && (covariates[1] == ""))
  
  # build design matrix and extract variable names
  if (misscovariates) {
    t1 <- terms(formula("~1"))
    param <- "(Intercept)"
    varnames <- ""
    xlevels <- NULL
  } else {
    fml_cov <- as.formula(paste("~", paste(covariates, collapse = "+")))
               
    # QUICK PATH: if all covariates present in df and are numeric, avoid model.matrix
    cov_present <- covariates %in% names(df)
    all_numeric <- FALSE
    if (all(cov_present)) {
      all_numeric <- all(vapply(df[ covariates ], is.numeric, logical(1)))
    }
    
    if (all_numeric) {
      # Build design columns directly from numeric covariates (intercept + columns)
      # This avoids model.matrix and is valid when covariates are simple numeric columns.
      param <- c("(Intercept)", covariates)
      varnames <- covariates
      t1 <- terms(fml_cov)
      xlevels <- NULL      
    } else {
      # FALLBACK (existing robust behavior): use model.frame + model.matrix on df
      mf <- model.frame(fml_cov, data = df, na.action = na.pass)
      mm <- model.matrix(fml_cov, mf)
      param <- colnames(mm)
      colnames(mm) <- make.names(colnames(mm))
      varnames <- colnames(mm)[-1]
      t1 <- terms(fml_cov)
      xlevels <- mf$xlev
      # copy model-matrix columns into df only if they are missing
      missing_cols <- setdiff(varnames, names(df))
      if (length(missing_cols) > 0) {
        for (vn in missing_cols) df[[vn]] <- mm[, vn, drop = TRUE]
      }
    }
  }
  
  # call the core fitting function
  fit <- liferegRcpp(
    data = df, 
    stratum = stratum, 
    time = time,
    time2 = time2, 
    event = event, 
    covariates = varnames,
    weight = weight, 
    offset = offset, 
    id = id, 
    dist = dist,
    init = init, 
    robust = robust, 
    plci = plci, 
    alpha = alpha, 
    maxiter = maxiter, 
    eps = eps)
  
  # post-process the output
  fit$p <- fit$sumstat$p[1]
  fit$nvar <- fit$sumstat$nvar[1]
  
  if (fit$p > 0) {
    par <- fit$parest$param[1:fit$p]
    if (length(par) > length(param)) {
      fit$param <- c(param, par[(1+length(param)):length(par)])
    } else {
      fit$param <- param
    }
    
    fit$beta <- fit$parest$beta
    names(fit$beta) <- fit$param
    
    if (fit$p > 1) {
      fit$vbeta <- as.matrix(fit$parest[, paste0("vbeta.", seq_len(fit$p))])
      if (robust) {
        fit$vbeta_naive <- as.matrix(fit$parest[, paste0("vbeta_naive.", seq_len(fit$p))])
      }
    } else {
      fit$vbeta <- as.matrix(fit$parest[, "vbeta", drop = FALSE])
      if (robust) {
        fit$vbeta_naive <- as.matrix(fit$parest[, "vbeta_naive", drop = FALSE])
      }
    }
    
    dimnames(fit$vbeta) <- list(fit$param, fit$param)
    if (robust) {
      dimnames(fit$vbeta_naive) <- list(fit$param, fit$param)
    }
  }
  
  fit$terms <- t1
  if (fit$p > 0) fit$xlevels <- xlevels
  
  fit$settings <- list(
    data = data, 
    stratum = stratum, 
    time = time, 
    time2 = time2,
    event = event, 
    covariates = covariates, 
    weight = weight, 
    offset = offset,
    id = id, 
    dist = dist, 
    init = init, 
    robust = robust, 
    plci = plci, 
    alpha = alpha, 
    maxiter = maxiter, 
    eps = eps
  )
  
  class(fit) <- "liferegr"
  fit
}
