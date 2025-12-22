#ifndef __SPLINES__
#define __SPLINES__

#include "dataframe_list.h" // FlatMatrix, ListCpp, DataFrameCpp


// Evaluate B-spline design (or derivatives) at points x.
// - knots: full knot vector (including boundary/internals).
// - x: evaluation points.
// - ord: spline order (degree + 1).
// - derivs: per-x derivative order (recycled if shorter than x).
// Returns an nx x ncol matrix stored in column-major FlatMatrix:
// element (r,c) is at data[c*nrow + r].
FlatMatrix splineDesigncpp(
    const std::vector<double>& knots,
    const std::vector<double>& x,
    int ord = 4,
    const std::vector<int>& derivs = std::vector<int>{0});

// Compute B-spline basis.
// - x: input vector (may contain NaNs represented by std::nan("")). 
// - df: degrees of freedom (optional if knots provided; pass <=0 if not used).
// - knots: internal knots (if empty and df > 0, knots are chosen).
// - degree: polynomial degree (default 3).
// - intercept: whether to include intercept column.
// - boundary_knots: pair {left, right}; if empty, computed from x.
// Returns a ListCpp with elements:
//   "basis" -> FlatMatrix (m x ncol, column-major), "dimnames" -> ListCpp, etc.
ListCpp bscpp(
    const std::vector<double>& x,
    int df,
    const std::vector<double>& knots = {},
    int degree = 3,
    bool intercept = false,
    const std::vector<double>& boundary_knots = {},
    bool warn_outside = true);

// Compute natural cubic spline basis.
ListCpp nscpp(
    const std::vector<double>& x,
    int df,
    const std::vector<double>& knots = {},
    bool intercept = false,
    const std::vector<double>& boundary_knots = {});



#endif // __SPLINES__