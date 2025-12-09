#include <Rcpp.h>
#include <R_ext/Applic.h>

using namespace Rcpp;

#ifndef __SPLINES__
#define __SPLINES__

NumericMatrix splineDesigncpp(
    NumericVector knots, 
    NumericVector x, 
    int ord, 
    IntegerVector derivs);

NumericMatrix bscpp(NumericVector x, int df, 
                    NumericVector knots, int degree, 
                    bool intercept, 
                    NumericVector boundary_knots, 
                    bool warn_outside);

NumericMatrix nscpp(NumericVector x, int df, 
                    NumericVector knots, bool intercept,
                    NumericVector boundary_knots);

#endif // __SPLINES__
