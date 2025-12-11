#ifndef __UTILITIES__
#define __UTILITIES__

#include <vector>
#include <string>
#include <algorithm>
#include <stdexcept>
#include <limits>
#include <functional>

#include <Rcpp.h>

// Use the concrete types defined in dataframe_list.h
#include "dataframe_list.h"

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

std::vector<int> seqcpp(int start, int end);

std::vector<int> which(const std::vector<bool>& vec);

std::vector<int> findInterval3(const std::vector<double>& x,
                               const std::vector<double>& v,
                               bool rightmost_closed = false,
                               bool all_inside = false,
                               bool left_open = false);

double brent(const std::function<double(double)>& f,
             double x1, double x2, double tol = std::pow(std::numeric_limits<double>::epsilon(), 0.25), int maxiter = 1000);

double bisect(const std::function<double(double)>& f,
              double x1, double x2, double tol = std::pow(std::numeric_limits<double>::epsilon(), 0.25), int maxiter = 1000);

double quantilecpp(const std::vector<double>& x, double p);

double squantilecpp(const std::function<double(double)>& S, double p, double tol = 1e-6);

template <typename T>
void subset_in_place(std::vector<T>& v, const std::vector<int>& order) {
  std::vector<T> temp_subset(order.size());
  for (size_t i = 0; i < order.size(); ++i) {
    int index = order[i];
    if (index < 0 || index >= static_cast<int>(v.size())) {
      throw std::out_of_range("Index in 'order' is out of bounds for the source vector.");
    }
    temp_subset[i] = v[index];
  }
  v.assign(temp_subset.begin(), temp_subset.end());
}

template <typename T>
std::vector<T> subset(const std::vector<T>& v, const std::vector<int>& order) {
  std::vector<T> result(order.size());
  for (size_t i = 0; i < order.size(); ++i) {
    int index = order[i];
    if (index < 0 || index >= static_cast<int>(v.size())) {
      throw std::out_of_range("Index in 'order' is out of bounds for the source vector.");
    }
    result[i] = v[index];
  }
  return result;
}

template <typename T>
std::vector<T> concat(const std::vector<T>& v1, const std::vector<T>& v2) {
  size_t total_size = v1.size() + v2.size();
  std::vector<T> result;
  result.reserve(total_size);
  std::copy(v1.begin(), v1.end(), std::back_inserter(result));
  std::copy(v2.begin(), v2.end(), std::back_inserter(result));
  return result;
}

template <typename T>
std::vector<T> unique_sorted(const std::vector<T>& v) {
  std::vector<T> w = v;
  std::sort(w.begin(), w.end());
  w.erase(std::unique(w.begin(), w.end()), w.end());
  return w;
}

template <typename T>
std::vector<int> matchcpp(const std::vector<T>& x, const std::vector<T>& table) {
  std::vector<int> result(x.size());
  for (size_t i = 0; i < x.size(); ++i) {
    auto it = std::find(table.begin(), table.end(), x[i]);
    if (it != table.end()) {
      result[i] = static_cast<int>(std::distance(table.begin(), it));
    } else {
      result[i] = -1;
    }
  }
  return result;
}

// Use FlatMatrix-based signatures (efficient, column-major)
ListCpp bygroup(const DataFrameCpp& data, const std::vector<std::string>& variables);

// Matrix operations using FlatMatrix (efficient implementations in utilities.cpp)
std::vector<double> mat_vec_mult(const FlatMatrix& A, const std::vector<double>& x);
FlatMatrix mat_mat_mult(const FlatMatrix& A, const FlatMatrix& B);
FlatMatrix transpose(const FlatMatrix& A);

// Linear algebra helpers operating on FlatMatrix (declarations match implementations)
int cholesky2(FlatMatrix& matrix, int n, double toler = std::pow(std::numeric_limits<double>::epsilon(), 0.75));
void chsolve2(FlatMatrix& matrix, int n, std::vector<double>& y);
void chinv2(FlatMatrix& matrix, int n);
FlatMatrix invsympd(const FlatMatrix& matrix, int n, double toler = std::pow(std::numeric_limits<double>::epsilon(), 0.75));

// Survival, QR and other helpers (use FlatMatrix where appropriate)
DataFrameCpp survsplit(const std::vector<double>& tstart,
                       const std::vector<double>& tstop,
                       const std::vector<double>& cut);

// QR uses FlatMatrix as input
ListCpp qrcpp(const FlatMatrix& X, double tol = 1e-07);

std::vector<int> match3(const std::vector<int>& id1,
                        const std::vector<double>& v1,
                        const std::vector<int>& id2,
                        const std::vector<double>& v2);

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

std::string sanitize(const std::string& s);

double qtpwexpcpp1(const double p,
                   const std::vector<double>& piecewiseSurvivalTime,
                   const std::vector<double>& lambda,
                   const double lowerBound = 0.0,
                   const bool lowertail = true,
                   const bool logp = false);

ListCpp getpsiest(double target,
                  const std::vector<double>& psi,
                  const std::vector<double>& Z,
                  int direction = 0);

double getpsiend(const std::function<double(double)>& f,
                 bool lowerend,
                 double initialend);

#endif // __UTILITIES__