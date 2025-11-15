#' @title Print method for ipcw objects
#' @description Prints the concise information of an ipcw fit.
#'
#' @param x An object of class \code{ipcw}.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A printout from the fit of an inverse-probability of censoring
#' weights model.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.ipcw <- function(x, ...) {
  if (is.na(x$cox_pvalue)) {
    pvalue = NA
  } else if (x$cox_pvalue > 0.9999) {
    pvalue = ">.9999" 
  } else if (x$cox_pvalue < 0.0001) {
    pvalue = "<.0001"
  } else {
    pvalue = formatC(x$cox_pvalue, format = "f", digits = 4)
  }
  
  
  df0 <- x$event_summary[,-1]
  rownames(df0) <- c("Control", "Treatment")
  j0 <- grep("pct$", names(df0))
  df0[j0] <- lapply(df0[j0], formatC, format = "f", digits = 1)
  print(df0, ..., na.print = "", quote = FALSE)
  cat("\n")
  
  
  level = paste0(100*(1 - x$settings$alpha), "%")
  
  df1 <- x$weight_summary[,-1]
  j1 = c(2,3,4,5,6,7)
  df1[j1] <- lapply(df1[j1], formatC, format = "f", digits = 4)
  if (x$settings$swtrt_control_only) {
    rownames(df1) <- c("Control")
    cat("                                Weight summary\n")
  } else {
    rownames(df1) <- c("Control", "Treatment")
    cat("                                  Weight summary\n")
  }
  print(df1, ..., na.print = "" , quote = FALSE )
  cat("\n")
  
  df2 <- data.frame(
    hr = c(x$hr, x$hr_CI[1], x$hr_CI[2]),
    pvalue = c(pvalue, "", "")
  )
  
  df2[1] <- lapply(df2[1], formatC, format = "f", digits = 3)
  df3 <- t(df2)
  
  colnames(df3) <- c("Estimate", paste("Lower", level), 
                     paste("Upper", level))
  
  rownames(df3) = c("Hazard ratio (HR)", "P-value")
  
  print(df3, ..., na.print = "" , quote = FALSE )
  
  invisible(x)
}
