#' @title Residuals for Proportional Hazards Regression Models
#' @description Obtains the martingale, deviance, score, or Schoenfeld
#' residuals for a proportional hazards regression model.
#'
#' @param object The output from the \code{phregr} call.
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
#' head(ressco)
#'
#' # Example 2 with counting process data
#' fit2 <- phregr(data = heart %>% mutate(rx = as.numeric(transplant) - 1),
#'                time = "start", time2 = "stop", event = "event",
#'                covariates = c("rx", "age"), id = "id", robust = TRUE)
#'
#' resssch <- residuals_phregr(fit2, type = "scaledsch")
#' head(resssch)
#'
#' @export
residuals_phregr <- function(
    object, type=c("martingale", "deviance", "score", "schoenfeld",
                   "dfbeta", "dfbetas", "scaledsch"),
    collapse=FALSE, weighted=(type %in% c("dfbeta", "dfbetas"))) {
  
  p = object$p
  beta = object$beta
  residuals = object$residuals

  data = object$data
  stratum = object$stratum
  time = object$time
  time2 = object$time2
  event = object$event
  covariates = object$covariates
  weight = object$weight
  offset = object$offset
  id = object$id
  ties = object$ties
  param = object$param

  rownames(data) = NULL
  
  elements = c(stratum, time, event, covariates, weight, offset, id)
  elements = unique(elements[elements != "" & elements != "none"])
  fml = formula(paste("~", paste(elements, collapse = "+")))
  mf = model.frame(fml, data = data, na.action = na.omit)
  
  rownum = as.integer(rownames(mf))
  df = data[rownum,]
  
  nvar = length(covariates)
  if (missing(covariates) || is.null(covariates) || (nvar == 1 && (
    covariates[1] == "" || tolower(covariates[1]) == "none"))) {
    p3 = 0
  } else {
    fml1 = formula(paste("~", paste(covariates, collapse = "+")))
    p3 = length(rownames(attr(terms(fml1), "factors")))
  }

  if (p >= 1 && p3 >= 1) {
    mf1 <- model.frame(fml1, data = df, na.action = na.pass)
    mm <- model.matrix(fml1, mf1)
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

  n <- length(residuals)

  vv <- object$vbeta_naive
  if (is.null(vv)) vv <- object$vbeta

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
  
  
  if (p == 0) { # null Cox model case
    beta = 0
    vv = matrix(0,1,1)
  }

  temp = residuals_phregcpp(p = p,
                            beta = beta,
                            vbeta = vv,
                            resmart = residuals,
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
                            type = type, 
                            collapse = collapse,
                            weighted = weighted)
  

  if (type == "martingale" || type == "deviance") {
    rr <- temp$resid
  } else {
    if (p == 1) {
      rr = c(temp$resid)
    } else {
      rr = temp$resid
    }
    
    if (type == "schoenfeld" || type == "scaledsch") {
      attr(rr, "time") = temp$time
      if (length(temp) == 3) {
        attr(rr, "strata") = temp$strata
      }
    }
  }
  
  if (is.matrix(rr)) colnames(rr) <- param
  
  rr
}
