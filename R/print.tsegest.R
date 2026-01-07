#' @title Print method for tsegest objects
#' @description Prints the concise information of a tsegest fit.
#'
#' @param x An object of class \code{tsegest}.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A printout from the fit of two-stage estimation with g-estimation.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.tsegest <- function(x, ...) {
  pvalue1 <- x$pvalue
  
  if (is.na(pvalue1)) {
    pvalue <- NA
  } else if (pvalue1 > 0.9999) {
    pvalue <- ">.9999" 
  } else if (pvalue1 < 0.0001) {
    pvalue <- "<.0001"
  } else {
    pvalue <- formatC(pvalue1, format = "f", digits = 4)
  }
  
  df0 <- x$event_summary[, -1]
  rownames(df0) <- c("Control", "Treatment")
  j0 <- grep("pct$", names(df0))
  df0[j0] <- lapply(df0[j0], formatC, format = "f", digits = 1)
  print(df0, ..., na.print = "", quote = FALSE)
  cat("\n")
  
  level <- paste0(100*(1 - x$settings$alpha), "%")
  
  if (x$settings$swtrt_control_only) {
    df1 <- data.frame(
      psi = c(x$psi, x$psi_CI[1], x$psi_CI[2]),
      surv_time_ratio = c(exp(-x$psi), exp(-x$psi_CI[2]), exp(-x$psi_CI[1])),
      hr = c(x$hr, x$hr_CI[1], x$hr_CI[2]),
      pvalue = c(pvalue, "", "")
    )
    
    j1 <- c(1,2,3)
    df1[j1] <- lapply(df1[j1], formatC, format = "f", digits = 3)
    
    df2 <- t(df1)
    
    colnames(df2) <- c("Estimate", paste("Lower", level), 
                       paste("Upper", level))
    
    rownames(df2) <- c("Causal parameter psi", "Causal survival time ratio", 
                       "Hazard ratio (HR)", 
                       paste0("P-value", " (", x$pvalue_type, ")"))
  } else {
    df1 <- data.frame(
      psi = c(x$psi, x$psi_CI[1], x$psi_CI[2]),
      surv_time_ratio = c(exp(-x$psi), exp(-x$psi_CI[2]), exp(-x$psi_CI[1])),
      psi_trt = c(x$psi_trt, x$psi_trt_CI[1], x$psi_trt_CI[2]),
      surv_time_ratio_trt = c(exp(-x$psi_trt), exp(-x$psi_trt_CI[2]), 
                              exp(-x$psi_trt_CI[1])),
      hr = c(x$hr, x$hr_CI[1], x$hr_CI[2]),
      pvalue = c(pvalue, "", "")
    )
    
    j1 <- c(1,2,3,4,5)
    df1[j1] <- lapply(df1[j1], formatC, format = "f", digits = 3)
    
    df2 <- t(df1)
    
    colnames(df2) <- c("Estimate", paste("Lower", level), 
                       paste("Upper", level))
    
    rownames(df2) <- c("Causal parameter psi for control arm", 
                       "Causal survival time ratio for control arm",
                       "Causal parameter psi for treatment arm",
                       "Causal survival time ratio for treatment arm",
                       "Hazard ratio (HR)", 
                       paste0("P-value", " (", x$pvalue_type, ")"))
  }
  
  print(df2, ..., na.print = "" , quote = FALSE )
  
  invisible(x)
}
