#' @title Print ipe Object
#' @description Prints the concise information of ipe fit.
#'
#' @param x The ipe object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A printout from the fit of a rank-preserving structural failure 
#' time model.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.ipe <- function(x, ...) {
  pvalue1 = if (x$settings$boot) x$cox_pvalue else x$logrank_pvalue
  if (pvalue1 > 0.9999) {
    pvalue = ">.9999" 
  } else if (pvalue1 < 0.0001) {
    pvalue = "<.0001"
  } else {
    pvalue = formatC(pvalue1, format = "f", digits = 4)
  }
  
  df1 <- data.frame(
    psi = c(x$psi, x$psi_CI[1], x$psi_CI[2]),
    surv_time_ratio = c(exp(-x$psi), exp(-x$psi_CI[2]), exp(-x$psi_CI[1])),
    hr = c(x$hr, x$hr_CI[1], x$hr_CI[2]),
    pvalue = c(pvalue, "", "")
  )
  
  j1 = c(1,2,3)
  df1[j1] <- lapply(df1[j1], formatC, format = "f", digits = 3)
  
  df2 <- t(df1)
  
  level = paste0(100*(1 - x$settings$alpha), "%")
  
  colnames(df2) <- c("Estimate", paste("Lower", level), paste("Upper", level))
  
  rownames(df2) = c("Causal parameter psi", "Causal survival time ratio", 
            "Hazard ratio (HR)", "P-value")
  
  print(df2, ..., na.print = "" , quote = FALSE )
  
  invisible(x)
}

