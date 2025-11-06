#include "utilities.h"
using namespace Rcpp;


//' @title B-Spline Design Matrix 
//' @description Computes the design matrix for B-splines based on the 
//' specified \code{knots} and evaluated at the values in \code{x}. 
//' 
//' @param knots A numeric vector specifying the positions of the knots, 
//'   including both boundary and internal knots. 
//' @param x A numeric vector of values where the B-spline functions 
//'   or their derivatives will be evaluated. The values of \code{x} 
//'   must lie within the range of the "inner" knots, i.e., between 
//'   \code{knots[ord]} and \code{knots[length(knots) - (ord - 1)]}. 
//' @param ord A positive integer indicating the order of the B-spline. 
//'   This corresponds to the number of coefficients in each piecewise 
//'   polynomial segment, where \code{ord = degree + 1}. 
//' @param derivs An integer vector specifying the order of derivatives 
//'   to be evaluated at the corresponding \code{x} values. Each value 
//'   must be between \code{0} and \code{ord - 1}, and the vector is 
//'   conceptually recycled to match the length of \code{x}. 
//'   The default is \code{0}, meaning the B-spline functions themselves 
//'   are evaluated.
//'   
//' @return A matrix with dimensions 
//' \code{c(length(x), length(knots) - ord)}. Each row corresponds 
//' to a value in \code{x} and contains the coefficients of the 
//' B-splines, or the specified derivatives, as defined by the 
//' \code{knots} and evaluated at that particular value of \code{x}. 
//' The total number of B-splines is \code{length(knots) - ord}, 
//' with each B-spline defined by a set of \code{ord} consecutive knots.
//' 
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' splineDesigncpp(knots = 1:10, x = 4:7)
//' splineDesigncpp(knots = 1:10, x = 4:7, derivs = 1)
//'
//' @export
// [[Rcpp::export]]
NumericMatrix splineDesigncpp(
    NumericVector knots = NA_REAL, 
    NumericVector x = NA_REAL, 
    int ord = 4, 
    IntegerVector derivs = IntegerVector::create(0)) {
  
  if (ord > 4) {
    stop("splines with ord > 4 are not implemented");
  }
  
  int nk = static_cast<int>(knots.size());
  if (nk <= 0) {
    stop("must have at least 'ord' knots");
  }
  if (!is_sorted(knots)) {
    knots.sort();
  }
  
  int nx = static_cast<int>(x.size());
  int nd = static_cast<int>(derivs.size());
  if (nd > nx) {
    stop("length of 'derivs' is larger than length of 'x'");
  }
  if (nd < 1) {
    stop("empty 'derivs'");
  }
  
  if (ord > nk || ord < 1) {
    stop("'ord' must be a positive integer, at most the number of knots");
  }
  if (nk < 2*ord - 1) {
    stop("need at least 2*ord - 1 knots");
  }
  
  int degree = ord - 1;
  
  bool need_outer = is_true(any((x < knots[ord-1]) | 
                            (knots[nk-degree-1] < x)));
  if (need_outer) {
    std::string errmsg = "the 'x' values must be in the range of ";
    errmsg = errmsg + "knots[ord], knots[nk - degree]";
    stop(errmsg);
  }
  
  NumericVector u = knots;
  auto f = [u](double y, int k, int ord, int deriv)->NumericVector {
    NumericVector s(ord);
    
    if (ord == 2) {
      if (deriv == 0) {
        s[0] = (u[k] - y)/(u[k] - u[k-1]);
        s[1] = (y - u[k-1])/(u[k] - u[k-1]);
      } else if (deriv == 1) {
        s[0] = (-1)/(u[k] - u[k-1]);
        s[1] = 1/(u[k] - u[k-1]);
      }
    } else if (ord == 3) {
      if (deriv == 0) {
        s[0] = pow(u[k] - y, 2)/
          ((u[k] - u[k-1])*(u[k] - u[k-2]));
        
        s[1] = (u[k] - y)*(y - u[k-2])/
          ((u[k] - u[k-1])*(u[k] - u[k-2])) + 
            (u[k+1] - y)*(y - u[k-1])/
              ((u[k+1] - u[k-1])*(u[k] - u[k-1]));
        
        s[2] = pow(y - u[k-1], 2)/
          ((u[k+1] - u[k-1])*(u[k] - u[k-1]));
      } else if (deriv == 1) {
        s[0] = 2*(y - u[k])/
          ((u[k] - u[k-1])*(u[k] - u[k-2]));
        
        s[1] = (-2*y + u[k] + u[k-2])/
          ((u[k] - u[k-1])*(u[k] - u[k-2])) + 
            (-2*y + u[k+1] + u[k-1])/
              ((u[k+1] - u[k-1])*(u[k] - u[k-1]));
        
        s[2] = 2*(y - u[k-1])/
          ((u[k+1] - u[k-1])*(u[k] - u[k-1]));
      } else if (deriv == 2) {
        s[0] = 2/((u[k] - u[k-1])*(u[k] - u[k-2]));
        
        s[1] = (-2)/((u[k] - u[k-1])*(u[k] - u[k-2])) + 
          (-2)/((u[k+1] - u[k-1])*(u[k] - u[k-1]));
        
        s[2] = 2/((u[k+1] - u[k-1])*(u[k] - u[k-1]));
      }
    } else if (ord == 4) {
      if (deriv == 0) {
        s[0] = pow(u[k] - y, 3)/
          ((u[k] - u[k-1])*(u[k] - u[k-2])*(u[k] - u[k-3]));
        
        s[1] = pow(u[k] - y, 2)*(y - u[k-3])/
          ((u[k] - u[k-1])*(u[k] - u[k-2])*(u[k] - u[k-3])) + 
            (u[k+1] - y)*(u[k] - y)*(y - u[k-2])/
              ((u[k+1] - u[k-2])*(u[k] - u[k-1])*(u[k] - u[k-2])) + 
                pow(u[k+1] - y, 2)*(y - u[k-1])/
                  ((u[k+1] - u[k-1])*(u[k+1] - u[k-2])*(u[k] - u[k-1]));
        
        s[2] = (u[k] - y)*pow(y - u[k-2], 2)/
          ((u[k+1] - u[k-2])*(u[k] - u[k-1])*(u[k] - u[k-2])) + 
            (u[k+1] - y)*(y - u[k-1])*(y - u[k-2])/
              ((u[k+1] - u[k-1])*(u[k+1] - u[k-2])*(u[k] - u[k-1])) + 
                (u[k+2] - y)*pow(y - u[k-1], 2)/
                  ((u[k+2] - u[k-1])*(u[k+1] - u[k-1])*(u[k] - u[k-1]));
        
        s[3] = pow(y - u[k-1], 3)/
          ((u[k+2] - u[k-1])*(u[k+1] - u[k-1])*(u[k] - u[k-1]));
      } else if (deriv == 1) {
        s[0] = -3*pow(y - u[k], 2)/
          ((u[k] - u[k-1])*(u[k] - u[k-2])*(u[k] - u[k-3]));
        
        s[1] = (y - u[k])*(3*y - u[k] - 2*u[k-3])/
          ((u[k] - u[k-1])*(u[k] - u[k-2])*(u[k] - u[k-3])) + 
            (y*(3*y - 2*u[k+1]) + u[k]*(-2*y + u[k+1]) + 
            u[k-2]*(-2*y + u[k] + u[k+1]))/
              ((u[k+1] - u[k-2])*(u[k] - u[k-1])*(u[k] - u[k-2])) + 
                (y - u[k+1])*(3*y - u[k+1] - 2*u[k-1])/
                  ((u[k+1] - u[k-1])*(u[k+1] - u[k-2])*(u[k] - u[k-1]));
        
        s[2] = -(y - u[k-2])*(3*y - 2*u[k] - u[k-2])/
          ((u[k+1] - u[k-2])*(u[k] - u[k-1])*(u[k] - u[k-2])) + 
            (u[k-1]*(2*y - u[k+1]) + u[k-2]*(2*y - u[k+1] - u[k-1]) + 
            y*(-3*y + 2*u[k+1]))/
              ((u[k+1] - u[k-1])*(u[k+1] - u[k-2])*(u[k] - u[k-1])) + 
                -(y - u[k-1])*(3*y - 2*u[k+2] - u[k-1])/
                ((u[k+2] - u[k-1])*(u[k+1] - u[k-1])*(u[k] - u[k-1]));
        
        s[3] = 3*pow(y - u[k-1], 2)/
          ((u[k+2] - u[k-1])*(u[k+1] - u[k-1])*(u[k] - u[k-1]));       
      } else if (deriv == 2) {
        s[0] = 6*(-y + u[k])/
          ((u[k] - u[k-1])*(u[k] - u[k-2])*(u[k] - u[k-3]));
        
        s[1] = 2*(3*y - 2*u[k] - u[k-3])/
          ((u[k] - u[k-1])*(u[k] - u[k-2])*(u[k] - u[k-3])) + 
            2*(3*y - u[k+1] - u[k] - u[k-2])/
              ((u[k+1] - u[k-2])*(u[k] - u[k-1])*(u[k] - u[k-2])) + 
                2*(3*y - 2*u[k+1] - u[k-1])/
                  ((u[k+1] - u[k-1])*(u[k+1] - u[k-2])*(u[k] - u[k-1]));
        
        s[2] = 2*(-3*y + u[k] + 2*u[k-2])/
          ((u[k+1] - u[k-2])*(u[k] - u[k-1])*(u[k] - u[k-2])) + 
            2*(-3*y + u[k+1] + u[k-1] + u[k-2])/
              ((u[k+1] - u[k-1])*(u[k+1] - u[k-2])*(u[k] - u[k-1])) + 
                2*(-3*y + u[k+2] + 2*u[k-1])/
                  ((u[k+2] - u[k-1])*(u[k+1] - u[k-1])*(u[k] - u[k-1]));
        
        s[3] = 6*(y - u[k-1])/
          ((u[k+2] - u[k-1])*(u[k+1] - u[k-1])*(u[k] - u[k-1]));
      } else if (deriv == 3) {
        s[0] = (-6)/((u[k] - u[k-1])*(u[k] - u[k-2])*(u[k] - u[k-3]));
        
        s[1] = 6/((u[k] - u[k-1])*(u[k] - u[k-2])*(u[k] - u[k-3])) + 
          6/((u[k+1] - u[k-2])*(u[k] - u[k-1])*(u[k] - u[k-2])) + 
          6/((u[k+1] - u[k-1])*(u[k+1] - u[k-2])*(u[k] - u[k-1]));
        
        s[2] = (-6)/((u[k+1] - u[k-2])*(u[k] - u[k-1])*(u[k] - u[k-2])) + 
          (-6)/((u[k+1] - u[k-1])*(u[k+1] - u[k-2])*(u[k] - u[k-1])) + 
          (-6)/((u[k+2] - u[k-1])*(u[k+1] - u[k-1])*(u[k] - u[k-1]));
        
        s[3] = 6/((u[k+2] - u[k-1])*(u[k+1] - u[k-1])*(u[k] - u[k-1]));
      }
    }
    
    return s;
  };
  

  int K = nk - 2*ord;
  
  IntegerVector derivs2(nx);
  for (int i=0; i<nx; ++i) {
    derivs2[i] = derivs[i % nd];
  }
  
  IntegerVector idx = findInterval3(x, u[Range(degree,K+ord)], 0, 0, 0);
  // include the maximum value in the last sub-interval
  idx = pmin(idx, K+1) + degree;
  
  NumericMatrix design(nx,K+ord);
  for (int i=0; i<nx; ++i) {
    int k = idx[i];
    NumericVector s = f(x[i], k, ord, derivs2[i]);
    for (int j=0; j<ord; ++j) {
      design(i, k-ord+j) = s[j];
    }
  }
  
  return design;
}


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
//' @return A matrix with dimensions \code{c(length(x), df)}. If 
//' \code{df} is provided, the matrix will have \code{df} columns. 
//' Alternatively, if \code{knots} are supplied, the number of columns 
//' will be \code{length(knots) + degree + intercept}. The matrix 
//' contains attributes that correspond to the arguments 
//' passed to the \code{bscpp} function. 
//' 
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' bscpp(women$height, df = 5)
//'
//' @export
// [[Rcpp::export]]
NumericMatrix bscpp(NumericVector x = NA_REAL, int df = NA_INTEGER, 
                    NumericVector knots = NA_REAL, int degree = 3, 
                    bool intercept = 0, 
                    NumericVector boundary_knots = NA_REAL, 
                    bool warn_outside = 1) {
  
  int m = static_cast<int>(x.size());
  
  int ord = 1 + degree;
  if (ord <= 1) {
    stop("'degree' must be a positive integer");
  }
  
  LogicalVector nax = is_na(x);
  bool nas = is_true(any(nax));
  NumericVector z = x[!nax];
  int n = static_cast<int>(z.size());
  
  LogicalVector outleft(n), outright(n), outside(n);
  if (is_false(any(is_na(boundary_knots)))) {
    boundary_knots.sort();
    outleft = z < boundary_knots[0];
    outright = z > boundary_knots[1];
    outside = outleft | outright;
  } else if (n == 1) {
    NumericVector boundary(2);
    boundary[0] = z[0]*7.0/8.0;
    boundary[1] = z[0]*9.0/8.0;
    boundary_knots = boundary;
  } else {
    NumericVector boundary(2);
    boundary[0] = min(z);
    boundary[1] = max(z);
    boundary_knots = boundary;
  }
  
  int K;
  bool mk_knots = (df != NA_INTEGER) && is_true(any(is_na(knots)));
  if (mk_knots) {
    K = df - ord + (1 - intercept);
    if (K < 0) {
      K = 0;
      warning("'df' is too small");
    }
    
    NumericVector knots1(K);
    NumericVector z1 = z[!outside];
    if (K > 0) {
      for (int k=0; k<K; ++k) {
        knots1[k] = quantilecpp(z1, (k+1.0)/(K+1.0));
      }
    }
    knots = knots1;
  } else {
    if (is_false(all(is_finite(knots)))) {
      stop("non-finite knots");
    }
    K = static_cast<int>(knots.size());
  }
  
  if (mk_knots && K > 0) {
    LogicalVector lrEq(2);
    lrEq[0] = is_true(any(boundary_knots == min(knots)));
    lrEq[1] = is_true(any(boundary_knots == max(knots)));
    if (is_true(any(lrEq))) {
      bool aE0 = 0, aE1 = 0;
      if (lrEq[0]) {
        double piv = boundary_knots[0];
        LogicalVector sub = knots == piv;
        aE0 = is_true(all(sub));
        if (aE0) {
          warning("all interior knots match left boundary knot");
        } else {
          NumericVector knots2 = knots[knots > piv];
          double shift = (min(knots2) - piv)/8;
          for (int i=0; i<K; ++i) {
            if (sub[i]) {
              knots[i] = knots[i] + shift;
            }
          }
        }
      }
      
      if (lrEq[1]) {
        double piv = boundary_knots[1];
        LogicalVector sub = knots == piv;
        aE1 = is_true(all(sub));
        if (aE1) {
          warning("all interior knots match right boundary knot");
        } else {
          NumericVector knots2 = knots[knots < piv];
          double shift = (piv - max(knots2))/8;
          for (int i=0; i<K; ++i) {
            if (sub[i]) {
              knots[i] = knots[i] - shift;
            }
          }
        }
      }
      
      if (!((lrEq[0] && aE0) || (lrEq[1] && aE1))) {
        warning("shoving 'interior' knots matching boundary knots to inside");
      }
    }
  }
  
  NumericVector u(K+2*ord);
  for (int k=0; k<ord; ++k) {
    u[k] = boundary_knots[0];
  }
  
  for (int k=ord; k<K+ord; ++k) {
    u[k] = knots[k-ord];
  }
  
  for (int k=K+ord; k<K+2*ord; ++k) {
    u[k] = boundary_knots[1];
  }
  
  IntegerVector deriv0(1); 
  deriv0[0] = 0;
  
  NumericMatrix design(n, K+ord);
  if (is_true(any(outside))) {
    if (warn_outside) {
      std::string warn_message = "some 'x' values beyond boundary knots ";
      warn_message = warn_message + "may cause ill-conditioned bases";
      warning(warn_message);
    }
    
    IntegerVector derivs = Range(0,degree);
    NumericVector scalef(ord); 
    scalef[0] = 1;
    for (int i=1; i<ord; ++i) {
      scalef[i] = scalef[i-1]*i;
    }
    
    double e = 0.25;
    if (is_true(any(outleft))) {
      double k_pivot = (1 - e)*boundary_knots[0] + e*u[ord];
      
      NumericVector zol = z[outleft];
      int r = static_cast<int>(zol.size());
      
      NumericMatrix zl(r,ord);
      for (int i=0; i<r; ++i) {
        zl(i,0) = 1.0;
        for (int j=1; j<=degree; ++j) {
          zl(i,j) = pow(zol[i] - k_pivot, j);
        }
      }
      
      NumericVector kp(ord,k_pivot);
      IntegerVector derivs = Range(0,ord-1);
      NumericMatrix tt = splineDesigncpp(u, kp, ord, derivs);
      int k = 0;
      for (int i=0; i<n; ++i) {
        if (outleft[i]) {
          for (int j=0; j<K+ord; ++j) {
            for (int l=0; l<ord; ++l) {
              design(i,j) += zl(k,l)*tt(l,j)/scalef[l];
            }
          }
          ++k;
        }
      }
    }
    
    if (is_true(any(outright))) {
      double k_pivot = (1 - e)*boundary_knots[1] + e*u[K+ord-1];
      
      NumericVector zor = z[outright];
      int r = static_cast<int>(zor.size());
      
      NumericMatrix zr(r,ord);
      for (int i=0; i<r; ++i) {
        zr(i,0) = 1.0;
        for (int j=1; j<=degree; ++j) {
          zr(i,j) = pow(zor[i] - k_pivot, j);
        }
      }
      
      NumericVector kp(ord, k_pivot);
      IntegerVector derivs = Range(0,ord-1);
      NumericMatrix tt = splineDesigncpp(u, kp, ord, derivs);
      int k = 0;
      for (int i=0; i<n; ++i) {
        if (outright[i]) {
          for (int j=0; j<K+ord; ++j) {
            for (int l=0; l<ord; ++l) {
              design(i,j) += zr(k,l)*tt(l,j)/scalef[l];
            }
          }
          ++k;
        }
      }
    }
    
    LogicalVector inside = !outside;
    if (is_true(any(inside))) {
      NumericMatrix v = splineDesigncpp(u, z[inside], ord, deriv0);
      int k = 0;
      for (int i=0; i<n; ++i) {
        if (inside[i]) {
          design(i,_) = v(k,_);
          ++k;
        }
      }
    }
  } else {
    design = splineDesigncpp(u, z, ord, deriv0);
  }
  
  if (!intercept) {
    design = design(_, Range(1,K+ord-1));
  }
  
  int ncol = design.ncol();
  
  NumericMatrix basis(m,ncol);
  if (nas) {
    basis.fill(NA_REAL);
    int j=0;
    for (int i=0; i<m; ++i) {
      if (!nax[i]) {
        basis(i,_) = design(j,_);
        ++j;
      }
    }
  } else {
    basis = design;
  }
  
  basis.attr("dimnames") = List::create(x.names(), seq(1,ncol));
  basis.attr("degree") = degree;
  basis.attr("knots") = knots;
  basis.attr("boundary_knots") = boundary_knots;
  basis.attr("intercept") = intercept;
  basis.attr("warn_outside") = warn_outside;
  basis.attr("class") = StringVector::create("bscpp", "basis", "matrix");
  return basis;
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
//' @return A matrix with dimensions \code{c(length(x), df)}, where 
//' \code{df} is either provided directly or computed as 
//' \code{length(knots) + 1 + intercept} when \code{knots} are supplied. 
//' The matrix contains attributes that correspond to the arguments 
//' passed to the \code{nscpp} function. 
//' 
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' nscpp(women$height, df = 5)
//'
//' @export
// [[Rcpp::export]]
NumericMatrix nscpp(NumericVector x = NA_REAL, int df = NA_INTEGER, 
                    NumericVector knots = NA_REAL, bool intercept = 0,
                    NumericVector boundary_knots = NA_REAL) {
  
  int m = static_cast<int>(x.size());
  
  LogicalVector nax = is_na(x);
  bool nas = is_true(any(nax));
  NumericVector z = x[!nax];
  int n = static_cast<int>(z.size());
  
  LogicalVector outleft(n), outright(n), outside(n);
  if (is_false(any(is_na(boundary_knots)))) {
    boundary_knots.sort();
    outleft = z < boundary_knots[0];
    outright = z > boundary_knots[1];
    outside = outleft | outright;
  } else if (n == 1) {
    NumericVector boundary(2);
    boundary[0] = z[0]*7.0/8.0;
    boundary[1] = z[0]*9.0/8.0;
    boundary_knots = boundary;
  } else {
    NumericVector boundary(2);
    boundary[0] = min(z);
    boundary[1] = max(z);
    boundary_knots = boundary;
  }
  
  int K;
  bool mk_knots = (df != NA_INTEGER) && is_true(any(is_na(knots)));
  if (mk_knots) {
    K = df - 1 - intercept;
    if (K < 0) {
      K = 0;
      warning("'df' is too small");
    }
    
    NumericVector knots1(K);
    NumericVector z1 = z[!outside];
    if (K > 0) {
      for (int k=0; k<K; ++k) {
        knots1[k] = quantilecpp(z1, (k+1.0)/(K+1.0));
      }
    }
    knots = knots1;
  } else {
    if (is_false(all(is_finite(knots)))) {
      stop("non-finite knots");
    }
    K = static_cast<int>(knots.size());
  }
  
  if (mk_knots && K > 0) {
    LogicalVector lrEq(2);
    lrEq[0] = is_true(any(boundary_knots == min(knots)));
    lrEq[1] = is_true(any(boundary_knots == max(knots)));
    if (is_true(any(lrEq))) {
      if (lrEq[0]) {
        double piv = boundary_knots[0];
        LogicalVector sub = knots == piv;
        if (is_true(all(sub))) {
          stop("all interior knots match left boundary knot");
        } 
        
        NumericVector knots2 = knots[knots > piv];
        double shift = (min(knots2) - piv)/8;
        for (int i=0; i<K; ++i) {
          if (sub[i]) {
            knots[i] = knots[i] + shift;
          }
        }
      }
      
      if (lrEq[1]) {
        double piv = boundary_knots[1];
        LogicalVector sub = knots == piv;
        if (is_true(all(sub))) {
          stop("all interior knots match right boundary knot");
        }
        
        NumericVector knots2 = knots[knots < piv];
        double shift = (piv - max(knots2))/8;
        for (int i=0; i<K; ++i) {
          if (sub[i]) {
            knots[i] = knots[i] - shift;
          }
        }
      }
      
      warning("shoving 'interior' knots matching boundary knots to inside");
    }
  }
  
  NumericVector u(K+8);
  for (int k=0; k<4; ++k) {
    u[k] = boundary_knots[0];
  }
  
  for (int k=4; k<K+4; ++k) {
    u[k] = knots[k-4];
  }
  
  for (int k=K+4; k<K+8; ++k) {
    u[k] = boundary_knots[1];
  }
  
  IntegerVector deriv0(1);
  deriv0[0] = 0;
  
  NumericMatrix design(n, K+4);
  if (is_true(any(outside))) {
    
    if (is_true(any(outleft))) {
      double k_pivot = boundary_knots[0];
      
      NumericVector zol = z[outleft];
      int r = static_cast<int>(zol.size());
      
      NumericMatrix zl(r,2);
      for (int i=0; i<r; ++i) {
        zl(i,0) = 1.0;
        zl(i,1) = zol[i] - k_pivot;
      }
      
      NumericVector kp(2,k_pivot);
      IntegerVector derivs = Range(0,1);
      NumericMatrix tt = splineDesigncpp(u, kp, 4, derivs);
      int k = 0;
      for (int i=0; i<n; ++i) {
        if (outleft[i]) {
          for (int j=0; j<K+4; ++j) {
            for (int l=0; l<2; ++l) {
              design(i,j) += zl(k,l)*tt(l,j);
            }
          }
          ++k;
        }
      }
    }
    
    if (is_true(any(outright))) {
      double k_pivot = boundary_knots[1];
      
      NumericVector zor = z[outright];
      int r = static_cast<int>(zor.size());
      
      NumericMatrix zr(r,4);
      for (int i=0; i<r; ++i) {
        zr(i,0) = 1.0;
        zr(i,1) = zor[i] - k_pivot;
      }
      
      NumericVector kp(2, k_pivot);
      IntegerVector derivs = Range(0,1);
      NumericMatrix tt = splineDesigncpp(u, kp, 4, derivs);
      int k = 0;
      for (int i=0; i<n; ++i) {
        if (outright[i]) {
          for (int j=0; j<K+4; ++j) {
            for (int l=0; l<2; ++l) {
              design(i,j) += zr(k,l)*tt(l,j);
            }
          }
          ++k;
        }
      }
    }
    
    LogicalVector inside = !outside;
    if (is_true(any(inside))) {
      NumericMatrix v = splineDesigncpp(u, z[inside], 4, deriv0);
      int k = 0;
      for (int i=0; i<n; ++i) {
        if (inside[i]) {
          design(i,_) = v(k,_);
          ++k;
        }
      }
    }
  } else {
    design = splineDesigncpp(u, z, 4, deriv0);
  }
  
  IntegerVector deriv2(1);
  deriv2[0] = 2;
  NumericMatrix con = splineDesigncpp(u, boundary_knots, 4, deriv2);
  
  if (!intercept) {
    con = con(_, Range(1,K+3));
    design = design(_, Range(1,K+3));
  }


  NumericMatrix tcon = transpose(con);
  List w = qrcpp(tcon, 1e-12);
  NumericMatrix Q = as<NumericMatrix>(w["Q"]);

  int L = static_cast<int>(design.ncol());
  NumericMatrix basis2(n,L-2);
  for (int i=0; i<n; ++i) {
    for (int j=0; j<L-2; ++j) {
      for (int k=0; k<L; ++k) {
        basis2(i,j) += design(i,k)*Q(k,j+2);
      }
    }
  }
  
  NumericMatrix basis(m,L-2);
  if (nas) {
    basis.fill(NA_REAL);
    int j=0;
    for (int i=0; i<m; ++i) {
      if (!nax[i]) {
        basis(i,_) = basis2(j,_);
        ++j;
      }
    }
  } else {
    basis = basis2;
  }
  
  basis.attr("dimnames") = List::create(x.names(), seq(1,L-2));
  basis.attr("degree") = 3;
  basis.attr("knots") = knots;
  basis.attr("boundary_knots") = boundary_knots;
  basis.attr("intercept") = intercept;
  basis.attr("class") = StringVector::create("nscpp", "basis", "matrix");
  return basis;
}
