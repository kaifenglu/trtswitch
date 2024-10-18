#' @title Print logisregr Object
#' @description Prints the concise information of logisregr fit.
#'
#' @param x The logisregr object to print.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A printout from the fit of a logistic regression model.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
print.logisregr <- function(x, ...) {
  lrchisq = -2*(x$sumstat$loglik0 - x$sumstat$loglik1)
  degrees = x$sumstat$p - 1
  pvalue = sapply(1:nrow(x$sumstat), function(i) {
    ifelse(degrees[i] > 0, 
           pchisq(lrchisq[i], degrees[i], 0, lower.tail = FALSE), 
           NA)
  })
  df1 <- cbind(x$sumstat[, c("n", "nevents", "loglik0", "loglik1")],
               lrchisq = lrchisq, df = degrees, pvalue = pvalue,
               x$sumstat[, c("niter", "firth", "flic")])

  p = x$p
  if (p > 0) {
    nreps = nrow(x$parest)/p
    
    if (!x$robust) {
      if (x$plci) {
        df = data.frame(param = rep(x$param, nreps),
                        coef = x$parest$beta,
                        expcoef = x$parest$expbeta,
                        se = x$parest$sebeta,
                        z = x$parest$z,
                        lower = x$parest$lower,
                        upper = x$parest$upper,
                        p = x$parest$p,
                        method = x$parest$method)
        
        if (nreps > 1) {
          df = cbind(df, x$parest[, (p+10):ncol(x$parest)])
          colnames(df) <- c("param", "coef", "exp(coef)", "se(coef)", "z",
                            paste("lower", 1-x$alpha),
                            paste("upper", 1-x$alpha), "p", "method",
                            colnames(x$parest)[(p+10):ncol(x$parest)])
        } else {
          colnames(df) <- c("param", "coef", "exp(coef)", "se(coef)", "z",
                            paste("lower", 1-x$alpha),
                            paste("upper", 1-x$alpha), "p", "method")
        }
      } else {
        df = data.frame(param = rep(x$param, nreps),
                        coef = x$parest$beta,
                        expcoef = x$parest$expbeta,
                        se = x$parest$sebeta,
                        z = x$parest$z,
                        p = x$parest$p)
        
        if (nreps > 1) {
          df = cbind(df, x$parest[, (p+10):ncol(x$parest)])
          colnames(df) <- c("param", "coef", "exp(coef)", "se(coef)", "z",
                            "p", colnames(x$parest)[(p+10):ncol(x$parest)])
        } else {
          colnames(df) <- c("param", "coef", "exp(coef)", "se(coef)", "z",
                            "p")
        }
      }
    } else {
      if (x$plci) {
        df = data.frame(param = rep(x$param, nreps),
                        coef = x$parest$beta,
                        expcoef = x$parest$expbeta,
                        nse = x$parest$sebeta_naive,
                        se = x$parest$sebeta,
                        z = x$parest$z,
                        lower = x$parest$lower,
                        upper = x$parest$upper,
                        p = x$parest$p,
                        method = x$parest$method)
        
        if (nreps > 1) {
          df = cbind(df, x$parest[, (2*p+11):ncol(x$parest)])
          colnames(df) <- c("param", "coef", "exp(coef)", "se(coef)",
                            "robust se", "z", paste("lower", 1-x$alpha),
                            paste("upper", 1-x$alpha), "p", "method",
                            colnames(x$parest)[(2*p+11):ncol(x$parest)])
        } else {
          colnames(df) <- c("param", "coef", "exp(coef)", "se(coef)",
                            "robust se", "z", paste("lower", 1-x$alpha),
                            paste("upper", 1-x$alpha), "p", "method")
        }
      } else {
        df = data.frame(param = rep(x$param, nreps),
                        coef = x$parest$beta,
                        expcoef = x$parest$expbeta,
                        nse = x$parest$sebeta_naive,
                        se = x$parest$sebeta,
                        z = x$parest$z,
                        p = x$parest$p)
        
        if (nreps > 1) {
          df = cbind(df, x$parest[, (2*p+11):ncol(x$parest)])
          colnames(df) <- c("param", "coef", "exp(coef)", "se(coef)",
                            "robust se", "z", "p",
                            colnames(x$parest)[(2*p+11):ncol(x$parest)])
        } else {
          colnames(df) <- c("param", "coef", "exp(coef)", "se(coef)",
                            "robust se", "z", "p")
        }
      }
    }
  }
  
  print(df1, ..., na.print = "" , quote = FALSE )
  cat("\n")
  print(df, ..., na.print = "" , quote = FALSE )
  invisible(x)
}
