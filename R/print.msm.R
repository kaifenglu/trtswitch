#' @title Print msm Object
#' @description Prints the concise information of msm fit.
#'
#' @param x The msm object to print.
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
print.msm <- function(x, ...) {
  level = paste0(100*(1 - x$settings$alpha), "%")
  if (x$cox_pvalue > 0.9999) {
    pvalue = ">.9999" 
  } else if (x$cox_pvalue < 0.0001) {
    pvalue = "<.0001"
  } else {
    pvalue = formatC(x$cox_pvalue, format = "f", digits = 4)
  }

  if (x$settings$swtrt_control_only) {
    if (x$settings$stabilized_weights) {
      s0 <- summary(x$data_outcome$stabilized_weight[
        x$data_outcome$treated==0])
    } else {
      s0 <- summary(x$data_outcome$unstabilized_weight[
        x$data_outcome$treated==0])
    }
    
    df1 <- data.frame(s0 = as.vector(s0))
    df1[1] <- lapply(df1[1], formatC, format = "f", digits = 4)
    
    df2 <- t(df1)
    
    colnames(df2) = c("Min.", "1st Qu.", "Median", "Mean", 
                      "3rd Qu.", "Max.")
    rownames(df2) <- c("Weight summary for control arm")
    
    print(df2, ..., na.print = "" , quote = FALSE )
    cat("\n")
  } else {
    if (x$settings$stabilized_weights) {
      s0 <- summary(x$data_outcome$stabilized_weight[
        x$data_outcome$treated==0])
      s1 <- summary(x$data_outcome$stabilized_weight[
        x$data_outcome$treated==1])
    } else {
      s0 <- summary(x$data_outcome$unstabilized_weight[
        x$data_outcome$treated==0])
      s1 <- summary(x$data_outcome$unstabilized_weight[
        x$data_outcome$treated==1])
    }
    
    
    df1 <- data.frame(s0 = as.vector(s0), s1 = as.vector(s1))
    j1 = c(1,2)
    df1[j1] <- lapply(df1[j1], formatC, format = "f", digits = 4)
    
    df2 <- t(df1)
    
    colnames(df2) = c("Min.", "1st Qu.", "Median", "Mean", 
                      "3rd Qu.", "Max.")
    rownames(df2) <- c("Weight summary for control arm", 
                       "Weight summary for treatment arm")
    
    print(df2, ..., na.print = "" , quote = FALSE )
    cat("\n")
  }
  

  df3 <- data.frame(
    hr = c(x$hr, x$hr_CI[1], x$hr_CI[2]),
    pvalue = c(pvalue, "", "")
  )
  
  df3[1] <- lapply(df3[1], formatC, format = "f", digits = 3)
  df4 <- t(df3)
  
  colnames(df4) <- c("Estimate", paste("Lower", level), 
                     paste("Upper", level))
  
  rownames(df4) = c("Hazard ratio (HR)", "P-value")
  
  print(df4, ..., na.print = "" , quote = FALSE )
  
  invisible(x)
}

