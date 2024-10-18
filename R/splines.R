#' @export
predict.bscpp <- function(object, newx, ...) {
  if(missing(newx)) return(object)
  a <- c(list(x = newx), attributes(object)[
    c("degree", "knots", "boundary_knots", "intercept", "warn_outside")])
  do.call("bscpp", a)
}

#' @export
predict.nscpp <- function(object, newx, ...) {
  if(missing(newx)) return(object)
  a <- c(list(x = newx), attributes(object)[
    c("knots", "boundary_knots", "intercept")])
  do.call("nscpp", a)
}

#' @export
makepredictcall.bscpp <- function(var, call) {
  if(as.character(call)[1L] == "bscpp" || 
     (is.call(call) && identical(eval(call[[1L]]), bscpp))) {
    at <- attributes(var)[c("degree", "knots", "boundary_knots", 
                            "intercept", "warn_outside")]
    call <- call[1L:2L]
    call[names(at)] <- at
  }
  call
}

#' @export
makepredictcall.nscpp <- function(var, call) {
  if(as.character(call)[1L] == "nscpp" || 
     (is.call(call) && identical(eval(call[[1L]]), nscpp))) {
    at <- attributes(var)[c("knots", "boundary_knots", "intercept")]
    call <- call[1L:2L]
    call[names(at)] <- at
  }
  call
}

