#' @title Assess Proportional Hazards Assumption Based on Supremum Test
#' @description Obtains the standardized score processes and the simulated
#' distribution under the null hypothesis as well as the p-values for 
#' the supremum tests.
#'
#' @param object The output from the \code{phregr} call.
#' @param resample The number of simulation samples for the supremem test.
#' @param seed The random seed for the simulations.
#'   
#' @details
#' The supremum test corresponds to the ASSESS statement with \code{ph} 
#' option of SAS PROC PHREG.
#' 
#' @return A list with the following components: 
#' 
#' * \code{time} the unique event times.
#' 
#' * \code{score_t} the observed standardized score process.
#' 
#' * \code{score_t_list} a list of simulated standardized score processes
#'   under the null hypothesis.
#' 
#' * \code{max_abs_value} the supremum of the absolute value of the observed
#'   standardized score process for each covariate and the supremum of 
#'   the sum of absolute values of the observed standardized score processes 
#'   across all covariates.
#' 
#' * \code{p_value} the p-values for the supremum tests for each covariate
#'   and the global test.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @references
#' D. Y. Lin, L. J. Wei, and Z. Ying.
#' Checking the Cox model with cumulative sums of martingale-based
#' residuals. 
#' Biometrika 1993; 80:557-572.
#'
#' @examples
#'
#' fit <- phregr(data = liver, time = "Time", event = "Status", 
#'               covariates = c("log(Bilirubin)", "log(Protime)", 
#'                              "log(Albumin)", "Age", "Edema"),
#'               ties = "breslow")
#'               
#' aph <- assess_phregr(fit, resample = 1000, seed = 314159)
#'   
#' aph
#' 
#' plot(aph, nsim = 20)
#' 
#' @export
assess_phregr <- function(object, resample = 1000, seed = 12345) {
  
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
  
  aph <- assess_phregRcpp(p = p,
                          beta = object$beta,
                          vbeta = object$vbeta,
                          data = df,
                          stratum = stratum,
                          time = time,
                          time2 = object$settings$time2,
                          event = event,
                          covariates = varnames,
                          weight = weight,
                          offset = offset,
                          ties = object$settings$ties,
                          resample = resample,
                          seed = seed)
  
  class(aph) <- "assess_phregr"
  aph
}
