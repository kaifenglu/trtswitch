#' @title Logistic Regression Models for Binary Data
#' @description Obtains the parameter estimates from logistic regression
#' models with binary data.
#'
#' @param data The input data frame that contains the following variables:
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
#' * \code{fitted}: The data frame with the following variables:
#'
#'     - \code{linear_predictors}: The linear fit on the link function scale.
#'
#'     - \code{fitted_values}: The fitted probabilities of having an event,
#'       obtained by transforming the linear predictors by the inverse of
#'       the link function.
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
#' * \code{settings}: A list containing the input parameter values.
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
logisregr <- function(data, event = "event", covariates = "",
                      freq = "", weight = "", offset = "", id = "",
                      link = "logit", init = NA_real_,
                      robust = FALSE, firth = FALSE,
                      flic = FALSE, plci = FALSE, alpha = 0.05,
                      maxiter = 50, eps = 1.0e-9) {
  
  # Determine if covariates were provided (empty string or NULL means no covariates)
  misscovariates <- length(covariates) == 0 || 
    (length(covariates) == 1 && (covariates[1] == ""))
  
  # Precompute the formulas used repeatedly by prepare_df
  elements <- unique(c(event, covariates, freq, weight, offset, id))
  elements <- elements[elements != ""]
  fml_all <- formula(paste("~", paste(elements, collapse = "+")))
  
  if (!misscovariates) {
    fml_cov <- formula(paste("~", paste(covariates, collapse = "+")))
  } else {
    fml_cov <- NULL
  }
  
  # Helper: prepare a single data.frame (subset rows with complete cases for
  # the variables of interest and add model-matrix columns if needed).
  prepare_df <- function(df_in) {
    # normalize rownames
    if (is.data.frame(df_in) && (inherits(df_in, "data.table") || 
                                 inherits(df_in, "tbl") ||
                                 inherits(df_in, "tbl_df"))) {
      df <- as.data.frame(df_in)
    } else {
      df <- df_in
    }
    rownames(df) <- NULL
    
    # FAST ROW SELECTION: use complete.cases on the relevant columns
    # (much cheaper than model.frame)
    # elements vector was computed earlier in the wrapper scope 
    # (event, covariates, freq, weight, offset)
    sel_cols <- intersect(elements, names(df))        # only existing columns
    if (length(sel_cols) > 0) {
      # use base::complete.cases on the subset of columns we care about
      rows_ok <- which(complete.cases(df[, sel_cols, drop = FALSE]))
      if (length(rows_ok) == 0) {
        # no complete rows -> return zero-row frame with same columns
        df <- df[integer(0), , drop = FALSE]
      } else {
        df <- df[rows_ok, , drop = FALSE]
      }
    }
    
    if (misscovariates) {
      # Intercept-only model: minimal metadata
      list(df = df,
           t1 = terms(formula("~1")),
           param = "(Intercept)",
           varnames = character(0),
           xlevels = NULL)
    } else {
      # QUICK PATH: if all covariates present in df and are numeric, avoid model.matrix
      cov_present <- covariates %in% names(df)
      all_numeric <- FALSE
      if (all(cov_present)) {
        all_numeric <- all(vapply(df[ covariates ], is.numeric, logical(1)))
      }
      
      if (all_numeric) {
        # Build design columns directly from numeric covariates (intercept + columns)
        # This avoids model.matrix and is valid when covariates are simple numeric columns.
        param <- c("(Intercept)", make.names(covariates))
        varnames <- if (length(covariates) > 0) make.names(covariates) else character(0)
        list(df = df,
             t1 = terms(fml_cov),
             param = param,
             varnames = varnames,
             xlevels = NULL)
      } else {
        # FALLBACK (existing robust behavior): use model.frame + model.matrix on df
        mf1 <- model.frame(fml_cov, data = df, na.action = na.pass)
        mm <- model.matrix(fml_cov, mf1)
        param <- colnames(mm)
        colnames(mm) <- make.names(colnames(mm))
        varnames <- if (ncol(mm) > 1) colnames(mm)[-1] else character(0)
        # copy model-matrix columns into df only if they are missing
        if (length(varnames) > 0) {
          missing_cols <- setdiff(varnames, names(df))
          if (length(missing_cols) > 0) {
            for (vn in missing_cols) df[[vn]] <- mm[, vn, drop = TRUE]
          }
        }
        list(df = df,
             t1 = terms(fml_cov),
             param = param,
             varnames = varnames,
             xlevels = mf1$xlev)
      }
    }
  }
  
  # Helper: post-process a raw C++ result (ListCpp wrapped into an R list)
  # into the user-friendly flat fit object produced previously.
  postprocess_one <- function(raw_fit, meta) {
    fit <- raw_fit
    fit$p <- fit$sumstat$p[1]
    fit$link <- fit$sumstat$link[1]
    
    if (fit$p > 0) {
      fit$param <- meta$param
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
    
    fit$linear_predictors <- fit$fitted$linear_predictors
    fit$fitted_values <- fit$fitted$fitted_values
    fit$terms <- meta$t1
    if (fit$p > 0) fit$xlevels <- meta$xlevels
    
    fit$settings <- list(
      data = data,
      event = event,
      covariates = covariates,
      freq = freq,
      weight = weight,
      offset = offset,
      id = id,
      link = link,
      init = init,
      robust = robust,
      firth = firth,
      flic = flic,
      plci = plci,
      alpha = alpha,
      maxiter = maxiter,
      eps = eps
    )
    
    class(fit) <- "logisregr"
    fit
  }
  
  if (inherits(data, "data.frame")) {
    meta <- prepare_df(data)
    raw <- logisregRcpp(
      data = meta$df,
      event = event,
      covariates = meta$varnames,
      freq = freq,
      weight = weight,
      offset = offset,
      id = id,
      link = link,
      init = init,
      robust = robust,
      firth = firth,
      flic = flic,
      plci = plci,
      alpha = alpha,
      maxiter = maxiter,
      eps = eps
    )
    fit <- postprocess_one(raw, meta)
    return(fit)
  } else {
    stop("Input 'data' must be a data frame");
  }
}
