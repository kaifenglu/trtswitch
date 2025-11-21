#' @title Assess Proportional Hazards Assumption
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
#'   standardized score process for each covariate.
#' 
#' * \code{max_sum_abs_value} the supremum of the sum of absolute values
#'   of the observed standardized score processes across all covariates.
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
#' library(dplyr)
#' 
#' fit <- phregr(data = pbc, time = "Time", event = "Status", 
#'               covariates = c("log(Bilirubin)", "log(Protime)", 
#'                              "log(Albumin)", "Age", "Edema"),
#'               ties = "breslow")
#'               
#' aph <- assess_phregr(fit, resample = 1000, seed = 314159)
#'   
#' aph
#'
#' @export
assess_phregr <- function(object, resample = 1000, seed = 12345) {
  
  if (!inherits(object, "phregr"))
    stop("object must be of class 'phregr'");
  
  p = object$p
  beta = object$beta
  vbeta = object$vbeta
  data = object$settings$data
  stratum = object$settings$stratum
  time = object$settings$time
  time2 = object$settings$time2
  event = object$settings$event
  covariates = object$settings$covariates
  weight = object$settings$weight
  offset = object$settings$offset
  ties = object$settings$ties
  
  rownames(data) = NULL
  
  elements = c(stratum, time, event, covariates, weight, offset)
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
  
  aph <- assess_phregcpp(p = p,
                         beta = beta,
                         vbeta = vbeta,
                         data = df,
                         stratum = stratum,
                         time = time,
                         time2 = time2,
                         event = event,
                         covariates = varnames,
                         weight = weight,
                         offset = offset,
                         resample = resample,
                         seed = seed)
  
  aph$covariates <- varnames
  aph$resample <- resample
  aph$seed <- seed
  
  class(aph) <- "assess_phregr"
  aph
}
