#include <Rcpp.h>
#include "splines.h"
#include "dataframe_list.h"

using namespace Rcpp;

//' @title B-Spline Basis for Polynomial Splines 
//' @description Computes the B-spline basis matrix for a given polynomial 
//' spline. 
//' 
//' @param x A numeric vector representing the predictor variable. 
//' @param df Degrees of freedom, specifying the number of columns in the 
//'   basis matrix. If \code{df} is provided, the function automatically 
//'   selects \code{df - degree - intercept} internal knots based on 
//'   appropriate quantiles of \code{x}, ignoring any missing values. 
//' @param knots A numeric vector specifying the internal breakpoints 
//'   that define the spline. If not provided, \code{df} must be specified. 
//' @param degree An integer specifying the degree of the piecewise 
//'   polynomial. The default value is \code{3}, which corresponds to 
//'   cubic splines. 
//' @param intercept A logical value indicating whether to include an 
//'   intercept in the basis. The default is \code{FALSE}. 
//' @param boundary_knots A numeric vector of length 2 specifying the 
//'   boundary points where the B-spline basis should be anchored. 
//'   If not supplied, the default is the range of non-missing values 
//'   in \code{x}. 
//' @param warn_outside A logical value indicating whether a warning 
//'   should be issued if any values of \code{x} fall outside the 
//'   specified boundary knots. 
//'
//' @return A list with a basis matrix with dimensions \code{c(length(x), df)}. 
//' If \code{df} is provided, the matrix will have \code{df} columns. 
//' Alternatively, if \code{knots} are supplied, the number of columns 
//' will be \code{length(knots) + degree + intercept}. The list also 
//' contains components corresponding to the input to the function, i.e., 
//' `degree`, `knots`, `boundary_knots`, and `intercept`. 
//' 
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' bsRcpp(women$height, df = 5)
//'
//' @export
// [[Rcpp::export]]
Rcpp::List bsRcpp(Rcpp::NumericVector x = NA_REAL, 
                  int df = NA_INTEGER, 
                  Rcpp::NumericVector knots = NA_REAL, 
                  int degree = 3, 
                  bool intercept = false, 
                  Rcpp::NumericVector boundary_knots = NA_REAL, 
                  bool warn_outside = true) {
  std::vector<double> xv = Rcpp::as<std::vector<double>>(x);
  std::vector<double> knotsv = Rcpp::as<std::vector<double>>(knots);
  std::vector<double> bkn = Rcpp::as<std::vector<double>>(boundary_knots);
  
  // call C++ implementation
  ListCpp out = bscpp(xv, df, knotsv, degree, intercept, bkn, warn_outside);
  
  // convert to R list
  return Rcpp::wrap(out);
}

//' @title Natural Cubic Spline Basis 
//' @description Computes the B-spline basis matrix for a natural cubic 
//' spline. 
//' 
//' @param x A numeric vector representing the predictor variable. 
//'   Missing values are allowed. 
//' @param df Degrees of freedom, specifying the number of columns in 
//'   the basis matrix. If \code{df} is provided, the function selects 
//'   \code{df - 1 - intercept} internal knots based on appropriate 
//'   quantiles of \code{x}, ignoring any missing values. 
//' @param knots A numeric vector specifying the internal breakpoints 
//'   that define the spline. If provided, the number of degrees of 
//'   freedom will be determined by the length of \code{knots}. 
//' @param intercept A logical value indicating whether to include an 
//'   intercept in the basis. The default is \code{FALSE}. 
//' @param boundary_knots A numeric vector of length 2 specifying the 
//'   boundary points where the natural boundary conditions are applied 
//'   and the B-spline basis is anchored. If not supplied, the default 
//'   is the range of non-missing values in \code{x}. 
//'
//' @return A list with a basis matrix with dimensions \code{c(length(x), df)}, 
//' where \code{df} is either provided directly or computed as 
//' \code{length(knots) + 1 + intercept} when \code{knots} are supplied. 
//' The list also contains components corresponding to the input to 
//' the function, i.e., `degree`, `knots`, `boundary_knots`, and `intercept`. 
//' 
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' nsRcpp(women$height, df = 5)
//'
//' @export
// [[Rcpp::export]]
Rcpp::List nsRcpp(Rcpp::NumericVector x = NA_REAL, 
                  int df = 0, 
                  Rcpp::NumericVector knots = NA_REAL, 
                  bool intercept = false,
                  Rcpp::NumericVector boundary_knots = NA_REAL) {
  std::vector<double> xv = Rcpp::as<std::vector<double>>(x);
  std::vector<double> knotsv = Rcpp::as<std::vector<double>>(knots);
  std::vector<double> bkn = Rcpp::as<std::vector<double>>(boundary_knots);
  
  ListCpp out = nscpp(xv, df, knotsv, intercept, bkn);
  return Rcpp::wrap(out);
}

