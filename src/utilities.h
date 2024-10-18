#include <Rcpp.h>
#include <R_ext/Applic.h>

using namespace Rcpp;

#ifndef __UTILITIES__
#define __UTILITIES__

void set_seed(int seed);

IntegerVector which(const LogicalVector& vector);

IntegerVector findInterval3(NumericVector x,
                            NumericVector breaks);

double brent(const std::function<double(double)>& f,
             double x1, double x2, double tol);

bool hasVariable(DataFrame df, std::string varName);

double quantilecpp(const NumericVector& x, const double p);

double squantilecpp(const std::function<double(double)>& S, double p);

IntegerVector c_vectors_i(IntegerVector vec1, IntegerVector vec2);
NumericVector c_vectors(NumericVector vec1, NumericVector vec2);
NumericMatrix subset_matrix_by_row(NumericMatrix a, IntegerVector q);
NumericMatrix c_matrices(NumericMatrix a1, NumericMatrix a2);

List bygroup(DataFrame data, const StringVector& variables);

int cholesky2(NumericMatrix matrix, int n, double toler);
void chsolve2(NumericMatrix matrix, int n, NumericVector y);
void chinv2(NumericMatrix matrix, int n);
NumericMatrix invsympd(NumericMatrix matrix, int n, double toler);

DataFrame survsplit(NumericVector tstart,
                    NumericVector tstop,
                    NumericVector cut);

bool is_sorted(NumericVector x);

NumericVector house(NumericVector x);
void row_house(NumericMatrix A, NumericVector v);
List qrcpp(NumericMatrix x, double tol);

#endif // __UTILITIES__
