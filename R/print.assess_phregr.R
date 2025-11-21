#' @title Print method for assess_phregr objects
#' @description Prints the concise information of an assess_phregr fit.
#'
#' @param x An object of class \code{assess_phregr}.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A printout from the fit of an assessment of proportional hazards
#' assumption of a Cox model.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.assess_phregr <- function(x, ...) {
  
  format_pvalue <- function(p) {
    # Handle the case of p < 0.0001
    ifelse(p < 0.0001,
           "<.0001",
           # Handle the case of p > 0.9999
           ifelse(p > 0.9999,
                  ">.9999",
                  # For all other cases, format to 4 decimal places
                  sprintf("%.4f", p)))
  }
  
  df <- data.frame(covariate = c(x$covariates, "Global test"),
                   max_abs_value = c(x$max_abs_value, 
                                     x$max_sum_abs_value),
                   resample = x$resample,
                   seed = x$seed,
                   p_value = format_pvalue(x$p_value))
  
  j0 <- 2
  df[j0] <- lapply(df[j0], formatC, format = "f", digits = 4)
  print(df, ..., na.print = "", quote = FALSE)
 
  invisible(x)
}
