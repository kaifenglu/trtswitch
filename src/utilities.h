#include <Rcpp.h>
#include <R_ext/Applic.h>

using namespace Rcpp;

#ifndef __UTILITIES__
#define __UTILITIES__

void set_seed(int seed);

IntegerVector which(const LogicalVector& vector);

IntegerVector findInterval3(NumericVector x, NumericVector v, 
                            bool rightmost_closed, bool all_inside, 
                            bool left_open);

double brent(const std::function<double(double)>& f,
             double x1, double x2, double tol);

double bisect(const std::function<double(double)>& f,
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

NumericVector house(const NumericVector& x);
void row_house(NumericMatrix& A, const int i1, const int i2,
               const int j1, const int j2, const NumericVector& v);
List qrcpp(const NumericMatrix& X, double tol);

IntegerVector match3(
    const IntegerVector id1, const NumericVector v1,
    const IntegerVector id2, const NumericVector v2);

DataFrame untreated(
    const double psi,
    const IntegerVector& id,
    const NumericVector& time,
    const IntegerVector& event,
    const IntegerVector& treat,
    const NumericVector& rx,
    const NumericVector& censor_time,
    const bool recensor,
    const bool autoswitch);

DataFrame unswitched(
    const double psi,
    const int n,
    const IntegerVector& id,
    const NumericVector& time,
    const IntegerVector& event,
    const IntegerVector& treat,
    const NumericVector& rx,
    const NumericVector& censor_time,
    const bool recensor,
    const bool autoswitch);

std::string sanitize(const std::string& s);

double qtpwexpcpp1(const double p,
                   const NumericVector& piecewiseSurvivalTime,
                   const NumericVector& lambda,
                   const double lowerBound,
                   const bool lowertail,
                   const bool logp);

List getpsiest(const double target, const NumericVector& psi, 
               const NumericVector& Z, const int direction);

double getpsiend(const std::function<double(double)>& f,
                 const bool lowerend, const double initialend);

#endif // __UTILITIES__
