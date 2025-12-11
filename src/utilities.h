// utilities.h - helper declarations for trtswitch2
#ifndef __UTILITIES_H__
#define __UTILITIES_H__

#include <vector>
#include <string>
#include <functional>
#include <cstddef> // for std::size_t
#include <algorithm>
#include <Rcpp.h>

#include "dataframe_list.h" // FlatMatrix, IntMatrix, DataFrameCpp, ListCpp

// --------------------------- Distribution helpers --------------------------
double boost_pnorm(double q, double mean = 0.0, double sd = 1.0, bool lower_tail = true);
double boost_qnorm(double p, double mean = 0.0, double sd = 1.0, bool lower_tail = true);
double boost_dnorm(double x, double mean = 0.0, double sd = 1.0);

double boost_plogis(double q, double location = 0.0, double scale = 1.0, bool lower_tail = true);
double boost_qlogis(double p, double location = 0.0, double scale = 1.0, bool lower_tail = true);
double boost_dlogis(double x, double location = 0.0, double scale = 1.0);

double boost_pextreme(double q, double location = 0.0, double scale = 1.0, bool lower_tail = true);
double boost_qextreme(double p, double location = 0.0, double scale = 1.0, bool lower_tail = true);
double boost_dextreme(double x, double location = 0.0, double scale = 1.0);

double boost_pchisq(double q, double df, bool lower_tail = true);
double boost_qchisq(double p, double df, bool lower_tail = true);

// --------------------------- Small utilities --------------------------------

// seqcpp: integer inclusive sequence; inputs are std::size_t
std::vector<std::size_t> seqcpp(std::size_t start, std::size_t end);

// which: return indices of true values in vector<bool>
std::vector<std::size_t> which(const std::vector<bool>& vec);

// findInterval3: adapted helper (returns int indices following R-like convention)
std::vector<int> findInterval3(const std::vector<double>& x,
                               const std::vector<double>& v,
                               bool rightmost_closed = false,
                               bool all_inside = false,
                               bool left_open = false);

// subset_in_place: reorder/keep elements of v according to 'order' (indices)
template <typename T>
void subset_in_place(std::vector<T>& v, const std::vector<std::size_t>& order) {
  std::vector<T> temp_subset(order.size());
  for (std::size_t i = 0; i < order.size(); ++i) {
    std::size_t index = order[i];
    if (index < 0 || index >= v.size()) {
      throw std::out_of_range("Index in 'order' is out of bounds for the source vector.");
    }
    temp_subset[i] = v[index];
  }
  v.assign(temp_subset.begin(), temp_subset.end());
}

// subset: return a subset of v according to 'order' (indices)
template <typename T>
std::vector<T> subset(const std::vector<T>& v, const std::vector<std::size_t>& order) {
  std::vector<T> result(order.size());
  for (std::size_t i = 0; i < order.size(); ++i) {
    std::size_t index = order[i];
    if (index < 0 || index >= v.size()) {
      throw std::out_of_range("Index in 'order' is out of bounds for the source vector.");
    }
    result[i] = v[index];
  }
  return result;
}

// concat: concatenate two vectors
template <typename T>
std::vector<T> concat(const std::vector<T>& v1, const std::vector<T>& v2) {
  std::size_t total_size = v1.size() + v2.size();
  std::vector<T> result;
  result.reserve(total_size);
  std::copy(v1.begin(), v1.end(), std::back_inserter(result));
  std::copy(v2.begin(), v2.end(), std::back_inserter(result));
  return result;
}

// unique_sorted: return sorted unique values
template <typename T>
std::vector<T> unique_sorted(const std::vector<T>& v) {
  std::vector<T> w = v;
  std::sort(w.begin(), w.end());
  w.erase(std::unique(w.begin(), w.end()), w.end());
  return w;
}

// matchcpp: for each element of x find its index in table or -1 if not found
template <typename T>
std::vector<int> matchcpp(const std::vector<T>& x, const std::vector<T>& table) {
  std::vector<int> result(x.size());
  for (std::size_t i = 0; i < x.size(); ++i) {
    auto it = std::find(table.begin(), table.end(), x[i]);
    if (it != table.end()) {
      result[i] = static_cast<int>(std::distance(table.begin(), it));
    } else {
      result[i] = -1;
    }
  }
  return result;
}


// --------------------------- Root finders -----------------------------------
double brent(const std::function<double(double)>& f,
             double x1, double x2, double tol = 1e-8, std::size_t maxiter = 100);
double bisect(const std::function<double(double)>& f,
              double x1, double x2, double tol = 1e-8, std::size_t maxiter = 100);

// --------------------------- Quantiles -------------------------------------
double quantilecpp(const std::vector<double>& x, double p);
double squantilecpp(const std::function<double(double)>& S, double p, double tol);

// --------------------------- Matrix utilities (FlatMatrix) ------------------
std::vector<double> mat_vec_mult(const FlatMatrix& A, const std::vector<double>& x);
FlatMatrix mat_mat_mult(const FlatMatrix& A, const FlatMatrix& B);
FlatMatrix transpose(const FlatMatrix& A);

// --------------------------- Linear algebra helpers (FlatMatrix-backed) ----
// NOTE: size_t is used for matrix dimension parameters to avoid signed/unsigned conversions.
// Algorithmic behavior preserved.
int cholesky2(FlatMatrix& matrix, std::size_t n, double toler = 1e-12);
void chsolve2(FlatMatrix& matrix, std::size_t n, std::vector<double>& y);
void chinv2(FlatMatrix& matrix, std::size_t n);
FlatMatrix invsympd(const FlatMatrix& matrix, std::size_t n, double toler = 1e-12);

// Survival helpers
DataFrameCpp survsplit(const std::vector<double>& tstart,
                       const std::vector<double>& tstop,
                       const std::vector<double>& cut);

// QR and other helpers
ListCpp qrcpp(const FlatMatrix& X, double tol = 1e-07);

// --------------------------- Matching and other helpers ----------------------
std::vector<int> match3(const std::vector<int>& id1,
                        const std::vector<double>& v1,
                        const std::vector<int>& id2,
                        const std::vector<double>& v2);

// --------------------------- Counterfactual helpers --------------------------
DataFrameCpp untreated(double psi,
                       const std::vector<int>& id,
                       const std::vector<double>& time,
                       const std::vector<int>& event,
                       const std::vector<int>& treat,
                       const std::vector<double>& rx,
                       const std::vector<double>& censor_time,
                       bool recensor,
                       bool autoswitch);

DataFrameCpp unswitched(double psi,
                        const std::vector<int>& id,
                        const std::vector<double>& time,
                        const std::vector<int>& event,
                        const std::vector<int>& treat,
                        const std::vector<double>& rx,
                        const std::vector<double>& censor_time,
                        bool recensor,
                        bool autoswitch);


// Misc math helpers
std::string sanitize(const std::string& s);

double qtpwexpcpp1(const double p,
                   const std::vector<double>& piecewiseSurvivalTime,
                   const std::vector<double>& lambda,
                   const double lowerBound = 0.0,
                   const bool lowertail = true,
                   const bool logp = false);

// Root selection helpers
ListCpp getpsiest(double target,
                  const std::vector<double>& psi,
                  const std::vector<double>& Z,
                  int direction);

double getpsiend(const std::function<double(double)>& f,
                 bool lowerend,
                 double initialend);

ListCpp bygroup(const DataFrameCpp& data, const std::vectorstd::string& variables);

#endif // __UTILITIES_H__