#' @title Residuals for Parametric Regression Models for Failure Time Data
#' @description Obtains the response, deviance, dfbeta, and likelihood
#' displacement residuals for a parametric regression model for failure 
#' time data.
#'
#' @param object The output from the \code{phregr} call.
#' @param type The type of residuals desired, with options including
#'   \code{"response"}, \code{"deviance"}, \code{"dfbeta"},
#'   \code{"dfbetas"}, \code{"working"}, \code{"ldcase"}, 
#'   \code{"ldresp"}, \code{"ldshape"}, and \code{"matrix"}.
#' @param collapse Whether to collapse the residuals by \code{id}.
#' @param weighted Whether to compute weighted residuals.
#'
#' @details
#' The algorithms follow the \code{residuals.survreg} function in the 
#' \code{survival} package.
#'
#' @return 
#' Either a vector or a matrix of residuals, depending on the specified type:
#' 
#' * \code{response} residuals are on the scale of the original data. 
#' 
#' * \code{working} residuals are on the scale of the linear predictor. 
#' 
#' * \code{deviance} residuals are on the log-likelihood scale. 
#' 
#' * \code{dfbeta} residuals are returned as a matrix, where the 
#'   \eqn{i}-th row represents the approximate change in the model 
#'   coefficients resulting from the inclusion of subject \eqn{i}. 
#'   
#' * \code{dfbetas} residuals are similar to \code{dfbeta} residuals, but 
#'   each column is scaled by the standard deviation of the 
#'   corresponding coefficient.
#'   
#' * \code{matrix} residuals are a matrix of derivatives of the 
#'   log-likelihood function. Let \eqn{L} be the log-likelihood, \eqn{p} be 
#'   the linear predictor (\eqn{X\beta}), and \eqn{s} be \eqn{log(\sigma)}.
#'   Then the resulting matrix contains six columns: \eqn{L}, 
#'   \eqn{\partial L/\partial p}, \eqn{\partial^2 L/\partial p^2},
#'   \eqn{\partial L/\partial s}, \eqn{\partial^2 L/\partial s^2}, and
#'   \eqn{\partial L^2/\partial p\partial s}.
#'   
#' * \code{ldcase} residulas are likelihood displacement for case weight 
#'   perturbation.
#' 
#' * \code{ldresp} residuals are likelihood displacement for response value 
#'   perturbation. 
#' 
#' * \code{ldshape} residuals are likelihood displacement related to the 
#'   shape parameter.
#' 
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @references
#' Escobar, L. A. and Meeker, W. Q. 
#' Assessing influence in regression analysis with censored data. 
#' Biometrics 1992; 48:507-528.
#'
#' @examples
#'
#' library(dplyr)
#' 
#' fit1 <- liferegr(
#'   data = tobin %>% mutate(time = ifelse(durable>0, durable, NA)),
#'   time = "time", time2 = "durable",
#'   covariates = c("age", "quant"), dist = "normal")
#'
#'  resid <- residuals_liferegr(fit1, type = "response")
#'
#' @export
residuals_liferegr <- function(
    object, type=c("response", "deviance", "dfbeta", "dfbetas",
                   "working", "ldcase", "ldresp", "ldshape", "matrix"),
    collapse=FALSE, weighted=(type %in% c("dfbeta", "dfbetas"))) {
  
  p = object$p
  beta = object$beta
  
  data = object$data
  stratum = object$stratum
  time = object$time
  time2 = object$time2
  event = object$event
  covariates = object$covariates
  weight = object$weight
  offset = object$offset
  id = object$id
  dist = object$dist
  param = object$param

  rownames(data) = NULL
  
  elements = c(stratum, covariates, weight, offset, id)
  elements = unique(elements[elements != "" & elements != "none"])
  if (!(length(elements) == 0)) {
    mf = model.frame(formula(paste("~", paste(elements, collapse = "+"))),
                     data = data)
  } else {
    mf = model.frame(formula("~1"), data = data)
  }
  
  rownum = as.integer(rownames(mf))
  df = data[rownum,]
  
  nvar = length(covariates)
  if (missing(covariates) || is.null(covariates) || (nvar == 1 && (
    covariates[1] == "" || tolower(covariates[1]) == "none"))) {
    t1 = terms(formula("~1"))
    p3 = 0
  } else {
    t1 = terms(formula(paste("~", paste(covariates, collapse = "+"))))
    t2 = attr(t1, "factors")
    t3 = rownames(t2)
    p3 = length(t3)
  }
  
  if (p >= 1 && p3 >= 1) {
    mf = model.frame(t1, df)
    mm = model.matrix(t1, mf)
    colnames(mm) = make.names(colnames(mm))
    varnames = colnames(mm)[-1]
    for (i in 1:length(varnames)) {
      if (!(varnames[i] %in% names(df))) {
        df[,varnames[i]] = mm[,varnames[i]]
      }
    }
  } else {
    varnames = ""
  }

  type <- match.arg(type)
  
  if (type=='dfbeta' || type=='dfbetas') {
    if (missing(weighted))
      weighted <- TRUE  # different default for this case
  }
  
  n = nrow(data)

  vv <- drop(object$vbeta_naive)
  if (is.null(vv)) vv <- drop(object$vbeta)

  if (weight != "") {
    weights <- df[[weight]]
  } else {
    weights <- rep(1,n)
  }

  if (id != "") {
    idn <- df[[id]]
  } else {
    idn <- seq(1,n)
  }
  
  rr = residuals_liferegcpp(beta = beta,
                            vbeta = vv,
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
                            type = type,
                            collapse = collapse,
                            weighted = weighted)
  
  if (type=="response" || type=="deviance" || type=="working" ||
      type=="ldcase" || type=="ldresp" || type=="ldshape") {
    rr <- as.numeric(rr)
  } else if (type=="dfbeta" || type=="dfbetas") {
    colnames(rr) <- param
  } else {
    colnames(rr) <- c("g", "dg", "ddg", "ds", "dds", "dsg");
  }

  rr
}
