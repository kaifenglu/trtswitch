#' @title Logistic Regression Models for Binary Data
#' @description Obtains the parameter estimates from logistic regression
#' models with binary data.
#'
#' @param data The input data frame that contains the following variables:
#'
#'   * \code{rep}: The replication for by-group processing.
#'
#'   * \code{event}: The event indicator, 1=event, 0=no event.
#'
#'   * \code{covariates}: The values of baseline covariates.
#'
#'   * \code{freq}: The frequency for each observation.
#'
#'   * \code{weight}: The weight for each observation.
#'
#'   * \code{offset}: The offset for each observation.
#'
#'   * \code{id}: The optional subject ID to group the score residuals
#'     in computing the robust sandwich variance.
#'
#' @param rep The name(s) of the replication variable(s) in the input data.
#' @param event The name of the event variable in the input data.
#' @param covariates The vector of names of baseline covariates
#'   in the input data.
#' @param freq The name of the frequency variable in the input data.
#'   The frequencies must be the same for all observations within each
#'   cluster as indicated by the id. Thus freq is the cluster frequency.
#' @param weight The name of the weight variable in the input data.
#' @param offset The name of the offset variable in the input data.
#' @param id The name of the id variable in the input data.
#' @param link The link function linking the response probabilities to the
#'   linear predictors. Options include "logit" (default), "probit", and
#'   "cloglog" (complementary log-log).
#' @param init A vector of initial values for the model parameters. 
#'   By default, initial values are derived from an intercept-only model. 
#' @param robust Whether a robust sandwich variance estimate should be
#'   computed. In the presence of the id variable, the score residuals
#'   will be aggregated for each id when computing the robust sandwich
#'   variance estimate.
#' @param firth Whether the firth's bias reducing penalized likelihood
#'   should be used. The default is \code{FALSE}.
#' @param flic Whether to apply intercept correction to obtain more
#'   accurate predicted probabilities. The default is \code{FALSE}.
#' @param plci Whether to obtain profile likelihood confidence interval.
#' @param alpha The two-sided significance level.
#' @param maxiter The maximum number of iterations.
#' @param eps The tolerance to declare convergence. 
#'
#' @details
#' Fitting a logistic regression model using Firth's bias reduction method
#' is equivalent to penalization of the log-likelihood by the Jeffreys prior.
#' Firth's penalized log-likelihood is given by
#' \deqn{l(\beta) + \frac{1}{2} \log(\mbox{det}(I(\beta)))}
#' and the components of the gradient \eqn{g(\beta)} are computed as
#' \deqn{g(\beta_j) + \frac{1}{2} \mbox{trace}\left(I(\beta)^{-1}
#' \frac{\partial I(\beta)}{\partial \beta_j}\right)}
#' The Hessian matrix is not modified by this penalty.
#'
#' Firth's method reduces bias in maximum likelihood estimates of
#' coefficients, but it introduces a bias toward one-half in the
#' predicted probabilities.
#'
#' A straightforward modification to Firth’s logistic regression to
#' achieve unbiased average predicted probabilities involves a post hoc
#' adjustment of the intercept. This approach, known as Firth’s logistic
#' regression with intercept correction (FLIC), preserves the
#' bias-corrected effect estimates. By excluding the intercept from
#' penalization, it ensures that we don't sacrifice the accuracy of
#' effect estimates to improve the predictions.
#'
#' @return A list with the following components:
#'
#' * \code{sumstat}: The data frame of summary statistics of model fit
#'   with the following variables:
#'
#'     - \code{n}: The number of subjects.
#'
#'     - \code{nevents}: The number of events.
#'
#'     - \code{loglik0}: The (penalized) log-likelihood under null.
#'
#'     - \code{loglik1}: The maximum (penalized) log-likelihood.
#'
#'     - \code{niter}: The number of Newton-Raphson iterations.
#'
#'     - \code{p}: The number of parameters, including the intercept,
#'       and regression coefficients associated with the covariates.
#'
#'     - \code{link}: The link function.
#'
#'     - \code{robust}: Whether a robust sandwich variance estimate should
#'       be computed.
#'
#'     - \code{firth}: Whether the firth's penalized likelihood is used.
#'
#'     - \code{flic}: Whether to apply intercept correction.
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
#'     - \code{sebeta_naive}: The naive standard error of parameter estimate.
#'
#'     - \code{vbeta_naive}: The naive covariance matrix of parameter
#'       estimates.
#'
#'     - \code{rep}: The replication.
#'
#' * \code{fitted}: The data frame with the following variables:
#'
#'     - \code{linear_predictors}: The linear fit on the link function scale.
#'
#'     - \code{fitted_values}: The fitted probabilities of having an event,
#'       obtained by transforming the linear predictors by the inverse of
#'       the link function.
#'
#'     - \code{rep}: The replication.
#'
#' * \code{p}: The number of parameters.
#'
#' * \code{link}: The link function.
#'
#' * \code{param}: The parameter names.
#'
#' * \code{beta}: The parameter estimate.
#'
#' * \code{vbeta}: The covariance matrix for parameter estimates.
#'
#' * \code{vbeta_naive}: The naive covariance matrix for parameter estimates.
#'
#' * \code{linear_predictors}: The linear fit on the link function scale.
#'
#' * \code{fitted_values}: The fitted probabilities of having an event.
#'
#' * \code{terms}: The terms object.
#'
#' * \code{xlevels}: A record of the levels of the factors used in fitting.
#'
#' * \code{data}: The input data.
#'
#' * \code{rep}: The name(s) of the replication variable(s).
#'
#' * \code{event}: The name of the event variable.
#'
#' * \code{covariates}: The names of baseline covariates.
#'
#' * \code{freq}: The name of the freq variable.
#'
#' * \code{weight}: The name of the weight variable.
#'
#' * \code{offset}: The name of the offset variable.
#'
#' * \code{id}: The name of the id variable.
#'
#' * \code{robust}: Whether a robust sandwich variance estimate should be
#'   computed.
#'
#' * \code{firth}: Whether to use the firth's bias reducing penalized
#'   likelihood.
#'
#' * \code{flic}: Whether to apply intercept correction.
#'
#' * \code{plci}: Whether to obtain profile likelihood confidence interval.
#'
#' * \code{alpha}: The two-sided significance level.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @references
#' David Firth.
#' Bias Reduction of Maximum Likelihood Estimates.
#' Biometrika 1993; 80:27–38.
#'
#' Georg Heinze and Michael Schemper.
#' A solution to the problem of separation in logistic regression.
#' Statistics in Medicine 2002;21:2409–2419.
#'
#' Rainer Puhr, Georg Heinze, Mariana Nold, Lara Lusa, and
#' Angelika Geroldinger.
#' Firth's logistic regression with rare events: accurate effect
#' estimates and predictions?
#' Statistics in Medicine 2017; 36:2302-2317.
#'
#' @examples
#'
#' (fit1 <- logisregr(
#'   ingots, event = "NotReady", covariates = "Heat*Soak", freq = "Freq"))
#'
#' @export
logisregr <- function(data, rep = "", event = "event", covariates = "",
                      freq = "", weight = "", offset = "", id = "",
                      link = "logit", init = NA_real_, 
                      robust = FALSE, firth = FALSE,
                      flic = FALSE, plci = FALSE, alpha = 0.05, 
                      maxiter = 50, eps = 1.0e-9) {
  
  rownames(data) = NULL
  
  elements = c(rep, event, covariates, freq, weight, offset)
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
  
  fit <- logisregcpp(data = df, rep = rep, event = event,
                     covariates = varnames, freq = freq, weight = weight,
                     offset = offset, id = id, link = link, 
                     init = init, robust = robust,
                     firth = firth, flic = flic, plci = plci,
                     alpha = alpha, maxiter = maxiter, eps = eps)
  
  fit$p <- fit$sumstat$p[1]
  fit$link <- fit$sumstat$link[1]
  
  if (fit$p > 0) {
    fit$param = param
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
  
  fit$linear_predictors = fit$fitted$linear_predictors
  fit$fitted_values = fit$fitted$fitted_values
  
  fit$terms = t1
  if (fit$p > 0) fit$xlevels = xlevels
  
  fit$data = data
  fit$rep = rep
  fit$event = event
  fit$covariates = covariates
  fit$freq = freq
  fit$weight = weight
  fit$offset = offset
  fit$id = id
  fit$robust = robust
  fit$firth = firth
  fit$flic = flic
  fit$plci = plci
  fit$alpha = alpha
  
  class(fit) = "logisregr"
  fit
}
