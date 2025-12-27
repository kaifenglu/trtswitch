#' @title Proportional Hazards Regression Models
#' @description Obtains the hazard ratio estimates from the proportional
#' hazards regression model with right censored or counting process data.
#'
#' @param data The input data frame that contains the following variables:
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
#' * \code{residuals}: The martingale residuals.
#' 
#' * \code{linear_predictors}: The vector of linear predictors.
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
#'   data = rawdata %>% filter(iterationNumber == 1) %>% 
#'     mutate(treat = 1*(treatmentGroup == 1)),
#'   stratum = "stratum",
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
phregr <- function(data, stratum = "", time = "time", time2 = "", 
                   event = "event", covariates = "", weight = "", 
                   offset = "", id = "", ties = "efron", 
                   init = NA_real_,  robust = FALSE,
                   est_basehaz = TRUE, est_resid = TRUE, 
                   firth = FALSE, plci = FALSE, alpha = 0.05, 
                   maxiter = 50, eps = 1.0e-9) {
  
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
  
  if (length(time) > 1) {
    stop("Only one time variable can be specified.")
  }
  if (length(time2) > 1) {
    stop("Only one time2 variable can be specified.")
  }
  
  if (length(event) == 0 || (length(event) == 1 && event[1] == "")) {
    stop("The event variable must be specified.")
  } else if (length(event) > 1) {
    stop("Only one event variable can be specified.")
  }
  if (length(weight) > 1) {
    stop("Only one weight variable can be specified.")
  }
  if (length(offset) > 1) {
    stop("Only one offset variable can be specified.")
  }
  if (length(id) > 1) {
    stop("Only one id variable can be specified.")
  }
  
  # select complete cases for the relevant variables
  elements <- c(stratum, time, time2, event, covariates, weight, offset, id)
  elements <- unique(elements[elements != "" & elements != "none"])
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
      # Use model.matrix to handle factors and interactions
      mf <- model.frame(fml_cov, data = df, na.action = na.pass)
      mm <- model.matrix(fml_cov, mf)
      param <- colnames(mm)
      colnames(mm) <- make.names(colnames(mm))
      varnames <- colnames(mm)[-1]
      t1 <- terms(fml_cov)
      xlevels <- mf$xlev
      # Add derived columns to df if not already present
      missing_cols <- setdiff(varnames, names(df))
      if (length(missing_cols) > 0) {
        for (vn in missing_cols) df[[vn]] <- mm[, vn, drop = TRUE]
      }
    }
  }
  
  # call the core fitting function
  fit <- phregRcpp(
    data = df, 
    stratum = stratum, 
    time = time,
    time2 = time2, 
    event = event, 
    covariates = varnames,
    weight = weight, 
    offset = offset, 
    id = id, 
    ties = ties, 
    init = init, 
    robust = robust, 
    est_basehaz = est_basehaz,
    est_resid = est_resid, 
    firth = firth,
    plci = plci, 
    alpha = alpha, 
    maxiter = maxiter, 
    eps = eps)
  
  # post-process the output
  fit$p <- fit$sumstat$p[1]
  
  if (fit$p > 0) {
    fit$param <- param[-1]
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
    ties = ties,
    iniy = init, 
    robust = robust, 
    est_basehaz = est_basehaz,
    est_resid = est_resid, 
    firth = firth, 
    plci = plci,
    alpha = alpha, 
    maxiter = maxiter, 
    eps = eps
  )
  
  class(fit) <- "phregr"
  fit
}
