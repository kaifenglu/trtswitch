#' @title Assess Proportional Hazards Assumption Based on Scaled 
#' Schoenfeld Residuals
#' @description Obtains the scaled Schoenfeld residuals and tests the 
#' proportional hazards assumption using a score test for the interaction
#' between each covariate and a transformed time variable.
#'
#' @param object The output from the \code{phregr} call.
#' @param transform A character string indicating how survival times
#'   should be transformed before the test is performed. Supported values
#'   include "identity", "log", "rank", and "km" (default).
#'   
#' @details
#' This corresponds to the \code{cox.zph} function from the \code{survival} 
#' package with \code{terms = FALSE} and \code{global = TRUE}.
#' 
#' @return A list with the following components: 
#' 
#' * \code{table} A matrix with one row for each parameter and a final 
#'   row for the global test. The columns contain the score test 
#'   for adding the time-dependent term, the degrees of freedom, 
#'   and the two-sided p-value.
#' 
#' * \code{x} The transformed time values.
#' 
#' * \code{time} The original (untransformed) event times, with tied event
#'   times repeated.
#' 
#' * \code{strata} The stratum index for each event.
#' 
#' * \code{y} The matrix of scaled Schoenfeld residuals, with one column
#'   for each parameter and one row for each event. Column names correspond 
#'   to the parameter names.
#' 
#' * \code{var} An approximate covariance matrix of the scaled Schoenfeld 
#'   residuals, used to construct an approximate standard error band for 
#'   plots.
#'   
#' * \code{transform} the transformation applied to the time values.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @references
#'
#' Patricia M. Grambsch and Terry M. Therneau.
#' Proportional hazards tests and diagnostics based on weighted residuals.
#' Biometrika 1994; 81:515-26.
#'
#' @examples
#'
#' library(dplyr)
#' 
#' fit <- phregr(data = liver, time = "Time", event = "Status", 
#'               covariates = c("log(Bilirubin)", "log(Protime)", 
#'                              "log(Albumin)", "Age", "Edema"),
#'               ties = "breslow")
#'               
#' zph <- zph_phregr(fit, transform = "log")
#'   
#' zph$table
#'
#' @export
zph_phregr <- function(object, transform = "km") {
  
  if (!inherits(object, "phregr")) stop("object must be of class 'phregr'");
  
  p = object$p
  df = object$settings$data
  stratum = object$settings$stratum
  time = object$settings$time
  event = object$settings$event
  covariates = object$settings$covariates
  weight = object$settings$weight
  offset = object$settings$offset
  
  elements <- unique(c(stratum, time, event, covariates, weight, offset))
  elements <- elements[elements != ""]
  fml_all <- formula(paste("~", paste(elements, collapse = "+")))
  var_all <- all.vars(fml_all)
  rows_ok <- which(complete.cases(df[, var_all, drop = FALSE]))
  if (length(rows_ok) == 0) stop("No complete cases found for the specified variables.")
  df <- df[rows_ok, , drop = FALSE]
  
  misscovariates <- length(covariates) == 0 || 
    (length(covariates) == 1 && (covariates[1] == ""))
  
  if (!misscovariates) {
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
      varnames <- covariates
    } else {
      # FALLBACK (existing robust behavior): use model.frame + model.matrix on df
      mf <- model.frame(fml_cov, data = df, na.action = na.pass)
      mm <- model.matrix(fml_cov, mf)
      colnames(mm) = make.names(colnames(mm))
      varnames = colnames(mm)[-1]
      missing_cols <- setdiff(varnames, names(df))
      if (length(missing_cols) > 0) {
        for (vn in missing_cols) df[[vn]] <- mm[, vn, drop = TRUE]
      }
    }
  } else {
    varnames <- ""
  }
  
  zph <- zph_phregRcpp(p = p,
                       beta = object$beta,
                       vbeta = object$vbeta,
                       resmart = object$residuals,
                       data = df,
                       stratum = stratum,
                       time = time,
                       time2 = object$settings$time2,
                       event = event,
                       covariates = varnames,
                       weight = weight,
                       offset = offset,
                       ties = object$settings$ties,
                       transform = transform)
  
  rownames(zph$table) <- c(object$param, "GLOBAL")
  colnames(zph$table) <- c("chisq", "df", "p")
  colnames(zph$y) <- object$param
  
  class(zph) <- "cox.zph"
  zph
}
