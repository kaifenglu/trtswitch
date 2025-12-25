#' @title Residuals for Parametric Regression Models for Failure Time Data
#' @description Obtains the response, martingale, deviance, dfbeta, and 
#' likelihood displacement residuals for a parametric regression model 
#' for failure time data.
#'
#' @param object The output from the \code{phregr} call.
#' @param type The type of residuals desired, with options including
#'   \code{"response"}, \code{"martingale"}, \code{"deviance"}, 
#'   \code{"dfbeta"}, \code{"dfbetas"}, \code{"working"}, \code{"ldcase"}, 
#'   \code{"ldresp"}, \code{"ldshape"}, and \code{"matrix"}.
#' @param collapse Whether to collapse the residuals by \code{id}.
#' @param weighted Whether to compute weighted residuals.
#'
#' @details
#' The algorithms follow the \code{residuals.survreg} function in the 
#' \code{survival} package, except for martingale residuals, which 
#' are defined only for event or right-censored data for exponential, 
#' weibull, lognormal, and loglogistic distributions.
#'
#' @return 
#' Either a vector or a matrix of residuals, depending on the specified type:
#' 
#' * \code{response} residuals are on the scale of the original data. 
#' 
#' * \code{martingale} residuals are event indicators minus the cumulative
#'   hazards for event or right-censored data.
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
#' resid <- residuals_liferegr(fit1, type = "response")
#' head(resid)
#'
#' @export
residuals_liferegr <- function(
    object, type=c("response", "martingale", "deviance", "dfbeta", "dfbetas",
                   "working", "ldcase", "ldresp", "ldshape", "matrix"),
    collapse=FALSE, weighted=(type %in% c("dfbeta", "dfbetas"))) {
  
  if (!inherits(object, "liferegr"))
    stop("object must be of class 'liferegr'");
  
  type <- match.arg(type)
  
  if (type=='dfbeta' || type=='dfbetas') {
    if (missing(weighted))
      weighted <- TRUE  # different default for this case
  }
  
  vv <- drop(object$vbeta_naive)
  if (is.null(vv)) vv <- drop(object$vbeta)
  
  rr = residuals_liferegRcpp(beta = object$beta,
                             vbeta = vv,
                             data = object$settings$data,
                             stratum = object$settings$stratum,
                             time = object$settings$time,
                             time2 = object$settings$time2,
                             event = object$settings$event,
                             covariates = object$settings$covariates,
                             weight = object$settings$weight,
                             offset = object$settings$offset,
                             id = object$settings$id,
                             dist = object$settings$dist,
                             type = type,
                             collapse = collapse,
                             weighted = weighted)
  
  if (type=="response" || type=="martingale" || type=="deviance" || 
      type=="working" || type=="ldcase" || type=="ldresp" || 
      type=="ldshape") {
    rr <- as.numeric(rr)
  } else if (type=="dfbeta" || type=="dfbetas") {
    colnames(rr) <- object$param
  } else {
    colnames(rr) <- c("g", "dg", "ddg", "ds", "dds", "dsg");
  }
  
  rr
}
