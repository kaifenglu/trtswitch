#' @title Print phregr Object
#' @description Prints the concise information of phregr fit.
#'
#' @param x The phregr object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A printout from the fit of a Cox proportional hazards model.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.phregr <- function(x, ...) {
  lrchisq <- -2*(x$sumstat$loglik0 - x$sumstat$loglik1)
  degrees <- x$sumstat$p
  pvalue <- sapply(1:nrow(x$sumstat), function(i) {
    ifelse(degrees[i] > 0, 
           pchisq(lrchisq[i], degrees[i], 0, lower.tail = FALSE), 
           NA)
  })
  df1 <- cbind(x$sumstat[, c("n", "nevents", "loglik0", "loglik1")],
               lrchisq = lrchisq, df = degrees, pvalue = pvalue,
               x$sumstat[, c("scoretest", "niter", "ties")])
  print(df1, ..., na.print = "" , quote = FALSE )
  cat("\n")
  
  p <- x$p
  if (p > 0) {
    if (!x$settings$robust) {
      if (x$settings$plci) {
        df <- data.frame(param = x$param,
                         coef = x$parest$beta,
                         expcoef = x$parest$expbeta,
                         se = x$parest$sebeta,
                         z = x$parest$z,
                         lower = x$parest$lower,
                         upper = x$parest$upper,
                         p = x$parest$p,
                         method = x$parest$method)
        
        colnames(df) <- c("param", "coef", "exp(coef)", "se(coef)", "z",
                          paste("lower", 1-x$settings$alpha),
                          paste("upper", 1-x$settings$alpha), "p", "method")
        
      } else {
        df <- data.frame(param = x$param,
                         coef = x$parest$beta,
                         expcoef = x$parest$expbeta,
                         se = x$parest$sebeta,
                         z = x$parest$z,
                         p = x$parest$p)
        
        colnames(df) <- c("param", "coef", "exp(coef)", "se(coef)", "z", "p")
      }
    } else {
      if (x$settings$plci) {
        df <- data.frame(param = x$param,
                         coef = x$parest$beta,
                         expcoef = x$parest$expbeta,
                         nse = x$parest$sebeta_naive,
                         se = x$parest$sebeta,
                         z = x$parest$z,
                         lower = x$parest$lower,
                         upper = x$parest$upper,
                         p = x$parest$p,
                         method = x$parest$method)
        
        colnames(df) <- c("param", "coef", "exp(coef)", "se(coef)",
                          "robust se", "z", 
                          paste("lower", 1-x$settings$alpha),
                          paste("upper", 1-x$settings$alpha), 
                          "p", "method")
      } else {
        df <- data.frame(param = x$param,
                         coef = x$parest$beta,
                         expcoef = x$parest$expbeta,
                         nse = x$parest$sebeta_naive,
                         se = x$parest$sebeta,
                         z = x$parest$z,
                         p = x$parest$p)
        
        colnames(df) <- c("param", "coef", "exp(coef)", "se(coef)",
                          "robust se", "z", "p")
      }
    }
    
    print(df, ..., na.print = "" , quote = FALSE )
  }
  
  invisible(x)
}
