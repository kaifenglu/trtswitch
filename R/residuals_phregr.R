#' @title Residuals for Proportional Hazards Regression Models
#' @description Obtains the martingale, deviance, score, or Schoenfeld
#' residuals for a proportional hazards regression model.
#'
#' @param fit_phregr The output from the \code{phregr} call.
#' @param type The type of residuals desired, with options including
#'   \code{"martingale"}, \code{"deviance"}, \code{"score"},
#'   \code{"schoenfeld"}, \code{"dfbeta"}, \code{"dfbetas"}, and
#'   \code{"scaledsch"}.
#' @param collapse Whether to collapse the residuals by \code{id}.
#'   This is not applicable for Schoenfeld type residuals.
#' @param weighted Whether to compute weighted residuals.
#'
#' @details
#' For score and Schoenfeld type residuals, the proportional hazards model
#' must include at least one covariate. The algorithms for \code{deviance},
#' \code{dfbeta}, \code{dfbetas}, and \code{scaledsch} residuals follow
#' the \code{residuals.coxph} function in the \code{survival} package.
#'
#' @return For martingale and deviance residuals, the result is a vector
#' with one element corresponding to each subject (without \code{collapse}).
#' For score residuals, the result is a matrix where each row represents
#' a subject and each column corresponds to a variable. The row order
#' aligns with the input data used in the original fit. For Schoenfeld
#' residuals, the result is a matrix with one row for each event and
#' one column per variable. These rows are sorted by time within strata,
#' with the attributes \code{stratum} and \code{time} included.
#'
#' Score residuals represent each individual's contribution to the score
#' vector. Two commonly used transformations of this are \code{dfbeta},
#' which represents the approximate change in the coefficient vector
#' if the observation is excluded, and \code{dfbetas}, which gives the
#' approximate change in the coefficients scaled by the standard error
#' of the coefficients.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @references
#' Terry M. Therneau, Patricia M. Grambsch, and Thomas M. Fleming.
#' Martingale based residuals for survival models.
#' Biometrika 1990; 77:147-160.
#'
#' Patricia M. Grambsch and Terry M. Therneau.
#' Proportional hazards tests and diagnostics based on weighted residuals.
#' Biometrika 1994; 81:515-26.
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
#' ressco <- residuals_phregr(fit1, type = "score")
#'
#' # Example 2 with counting process data
#' fit2 <- phregr(data = heart %>% mutate(rx = as.numeric(transplant) - 1),
#'                time = "start", time2 = "stop", event = "event",
#'                covariates = c("rx", "age"), id = "id", robust = TRUE)
#'
#' resssch <- residuals_phregr(fit2, type = "scaledsch")
#'
#' @export
residuals_phregr <- function(
    fit_phregr, type=c("martingale", "deviance", "score", "schoenfeld",
                       "dfbeta", "dfbetas", "scaledsch"),
    collapse=FALSE, weighted=(type %in% c("dfbeta", "dfbetas"))) {

  p = fit_phregr$p
  beta = fit_phregr$beta
  residuals = fit_phregr$residuals

  data = fit_phregr$data
  stratum = fit_phregr$stratum
  time = fit_phregr$time
  time2 = fit_phregr$time2
  event = fit_phregr$event
  covariates = fit_phregr$covariates
  weight = fit_phregr$weight
  offset = fit_phregr$offset
  id = fit_phregr$id
  ties = fit_phregr$ties
  param = fit_phregr$param

  rownames(data) = NULL
  
  elements = c(stratum, time, event, covariates, weight, offset, id)
  elements = unique(elements[elements != "" & elements != "none"])
  mf = model.frame(formula(paste("~", paste(elements, collapse = "+"))),
                   data = data)
  rownum = as.integer(rownames(mf))
  df = data[rownum,]
  
  nvar = length(covariates)
  if (missing(covariates) || is.null(covariates) || (nvar == 1 && (
    covariates[1] == "" || tolower(covariates[1]) == "none"))) {
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
  otype <- type
  if (type=='dfbeta' || type=='dfbetas') {
    type <- 'score'
    if (missing(weighted))
      weighted <- TRUE  # different default for this case
  }
  if (type=='scaledsch') type<-'schoenfeld'

  n <- length(residuals)
  rr <- residuals

  vv <- drop(fit_phregr$vbeta_naive)
  if (is.null(vv)) vv <- drop(fit_phregr$vbeta)

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


  if (type == 'martingale') rr <- fit_phregr$residuals

  if (type=='schoenfeld') {
    if (p == 0) stop("covariates must be present for schoenfeld residuals")

    temp = residuals_phregcpp(p = p,
                              beta = beta,
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
                              type = "schoenfeld")

    if (p==1) {
      rr <- c(temp$resid)
    } else {
      rr <- temp$resid
    }
    if (weighted) rr <- rr * weights[temp$obs]
    if (length(unique(temp$stratumn)) > 1) {
      attr(rr, "stratum") <- temp$stratumn
    }
    attr(rr, "time") <- temp$time

    if (otype=='scaledsch') {
      ndead <- length(temp$obs)
      if (nvar==1) {
        rr <- rr * vv * ndead + beta
      } else {
        rr <- drop(rr %*% vv *ndead + rep(beta, each=nrow(rr)))
      }
    }

    if (is.matrix(rr)) colnames(rr) <- param

    return(rr)
  }

  if (type=='score') {
    if (p == 0) stop("covariates must be present for score residuals")

    temp = residuals_phregcpp(p = p,
                              beta = beta,
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
                              type = "score")


    if (p==1) {
      rr <- c(temp$resid)
    } else {
      rr <- temp$resid
    }

    if (otype=='dfbeta') {
      if (is.matrix(rr)) {
        rr <- rr %*% vv
      } else {
        rr <- rr * vv
      }
    }
    else if (otype=='dfbetas') {
      if (is.matrix(rr)) {
        rr <- (rr %*% vv) %*% diag(sqrt(1/diag(vv)))
      } else {
        rr <- rr * sqrt(vv)
      }
    }

    if (is.matrix(rr)) colnames(rr) <- param
  }

  #
  # Multiply up by case weights (which will be 1 for unweighted)
  #
  if (weighted) rr <- rr * weights

  status <- df[[event]]
  
  # Collapse if desired
  if (collapse) {
    rr <- drop(rowsum(rr, idn))
    if (type=='deviance') status <- drop(rowsum(status, idn))
  }

  # Deviance residuals are computed after collapsing occurs
  if (type=='deviance') {
    sign(rr) *sqrt(-2* (rr+ ifelse(status==0, 0, status*log(status-rr))))
  } else {
    rr
  }
}
