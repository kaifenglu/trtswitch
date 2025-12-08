#' @title Logistic Regression Models for Binary Data
#' @description Obtains the parameter estimates from logistic regression
#' models with binary data.
#'
#' @param data The input data frame or list of data frames that contains 
#' the following variables:
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
#' @param nthreads The number of threads to use in the computation.
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
#' @return A list (or list of lists if the input is a list of data frames) 
#' with the following components:
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
                      maxiter = 50, eps = 1.0e-9,
                      nthreads = 0) {
  
  if (nthreads > 0) {
    RcppParallel::setThreadOptions(min(nthreads, parallel::detectCores(logical = FALSE)))
  }
  
  misscovariates = missing(covariates) || is.null(covariates) ||
    (length(covariates) == 1 && covariates[1] == "");
  
  # Helper: prepare a single data.frame (subset rows with complete cases for
  # the variables of interest and add model-matrix columns if needed).
  prepare_df <- function(df) {
    rownames(df) <- NULL
    elements <- c(event, covariates, freq, weight, offset, id)
    elements <- unique(elements[elements != ""])
    fml <- formula(paste("~", paste(elements, collapse = "+")))
    mf <- model.frame(fml, data = df, na.action = na.omit)
    rownum <- as.integer(rownames(mf))
    df2 <- df[rownum, , drop = FALSE]
    
    # determine covariate expansion (model matrix) and produce the columns used
    if (misscovariates) {
      p <- 0
      t1 <- terms(formula("~1"))
      xlevels <- NULL
      param <- "(Intercept)"
      varnames <- ""
    } else {
      fml1 <- formula(paste("~", paste(covariates, collapse = "+")))
      mf1 <- model.frame(fml1, data = df2, na.action = na.pass)
      mm <- model.matrix(fml1, mf1)
      xlevels <- mf1$xlev
      param <- colnames(mm)
      colnames(mm) <- make.names(colnames(mm))
      varnames <- colnames(mm)[-1]   # exclude intercept name
      # copy model-matrix columns into df2 if needed (so C++ can find them)
      for (vn in varnames) {
        if (!(vn %in% names(df2))) {
          df2[[vn]] <- mm[, vn, drop = TRUE]
        }
      }
      p <- length(varnames)
      t1 <- terms(fml1)
    }
    
    list(df = df2, p = p, t1 = t1, param = param, varnames = varnames, xlevels = xlevels)
  }
  
  # Helper: post-process a raw C++ result (ListCpp wrapped into an R list)
  # into the user-friendly flat fit object produced previously.
  postprocess_one <- function(raw_fit, meta) {
    # raw_fit: list returned from C++ (sumstat, parest, fitted)
    # meta: list containing p, t1, param, varnames, xlevels, robust
    fit <- raw_fit
    fit$p <- fit$sumstat$p[1]
    fit$link <- fit$sumstat$link[1]
    
    if (fit$p > 0) {
      fit$param <- meta$param
      fit$beta <- fit$parest$beta
      names(fit$beta) <- rep(fit$param, length(fit$beta) / fit$p)
      
      if (fit$p > 1) {
        # parest contains columns vbeta.1 ... vbeta.p in the C output style
        fit$vbeta <- as.matrix(fit$parest[, paste0("vbeta.", 1:fit$p)])
        if (meta$robust) {
          fit$vbeta_naive <- as.matrix(fit$parest[, paste0("vbeta_naive.", 1:fit$p)])
        }
      } else {
        fit$vbeta <- as.matrix(fit$parest[, "vbeta"])
        if (meta$robust) {
          fit$vbeta_naive <- as.matrix(fit$parest[, "vbeta_naive"])
        }
      }
      
      dimnames(fit$vbeta) <- list(names(fit$beta), fit$param)
      if (meta$robust) {
        dimnames(fit$vbeta_naive) <- list(names(fit$beta), fit$param)
      }
    }
    
    fit$linear_predictors <- fit$fitted$linear_predictors
    fit$fitted_values <- fit$fitted$fitted_values
    fit$terms <- meta$t1
    if (fit$p > 0) fit$xlevels <- meta$xlevels
    
    fit$settings <- list(
      data = meta$orig_data,
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
  
  # Decide whether data is a single data.frame or a list of data.frames
  if (inherits(data, "data.frame")) {
    meta <- prepare_df(data)
    # Call the parallel-capable C++ function but with a single data.frame (it will run on main thread)
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
    # raw is the C output as a flat list
    fit <- postprocess_one(raw, c(meta, list(orig_data = data, robust = robust)))
    return(fit)
  }
  
  # If it's a list, check elements are data.frames
  if (is.list(data)) {
    if (length(data) == 0) return(list())
    # Prepare each dataset, accumulate prepared dfs and metadata
    prepared <- vector("list", length(data))
    metas <- vector("list", length(data))
    for (i in seq_along(data)) {
      if (!inherits(data[[i]], "data.frame")) {
        stop("When 'data' is a list, every element must be a data.frame (or inherit 'data.frame').")
      }
      prepared_i <- prepare_df(data[[i]])
      prepared[[i]] <- prepared_i$df
      metas[[i]] <- c(prepared_i, list(orig_data = data[[i]], robust = robust))
    }
    
    # Check whether all prepared varnames are identical across datasets.
    varnames_list <- lapply(metas, function(m) m$varnames)
    same_varnames <- TRUE
    if (length(varnames_list) > 0) {
      ref <- varnames_list[[1]]
      same_varnames <- all(vapply(varnames_list, function(x) identical(x, ref), logical(1)))
    } else {
      same_varnames <- TRUE
    }
    
    if (same_varnames) {
      # We can call the parallel C++ wrapper once with the list of prepared dfs
      # and a single covariates vector (the common varnames)
      prepared_list <- prepared
      raw_list <- logisregRcpp(
        data = prepared_list,
        event = event,
        covariates = metas[[1]]$varnames,
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
      # raw_list is a list of raw C results
      out <- vector("list", length(raw_list))
      for (i in seq_along(raw_list)) {
        out[[i]] <- postprocess_one(raw_list[[i]], metas[[i]])
      }
      return(out)
    } else {
      # Fallback: varnames differ; process sequentially to preserve correctness.
      # (We could implement a more advanced parallel version that accepts
      # per-dataset covariates, but for now we keep correctness.)
      out <- vector("list", length(data))
      for (i in seq_along(prepared)) {
        raw <- logisregRcpp(
          data = prepared[[i]],
          event = event,
          covariates = metas[[i]]$varnames,
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
        out[[i]] <- postprocess_one(raw, metas[[i]])
      }
      return(out)
    }
  }
  
  stop("Input 'data' must be either a data.frame or a list of data.frames.")
}

