#ifndef __UTILITIES_H__
#define __UTILITIES_H__

struct FlatMatrix;
struct DataFrameCpp;
struct ListCpp;

#include <algorithm>   // std::sort, std::unique, std::find, std::copy
#include <cstddef>     // std::size_t
#include <functional>  // std::function
#include <iomanip>     // std::fixed, std::setprecision
#include <iostream>    // std::ostream, std::cout
#include <iterator>    // std::back_inserter, std::distance
#include <limits>      // std::numeric_limits
#include <sstream>     // std::ostringstream
#include <stdexcept>   // std::out_of_range
#include <string>      // std::string
#include <type_traits> // std::is_convertible
#include <utility>     // std::declval
#include <vector>      // std::vector

inline constexpr double NaN = std::numeric_limits<double>::quiet_NaN();

// --------------------------- Distribution helpers --------------------------
double boost_pnorm(double q, double mean = 0.0, double sd = 1.0, 
                   bool lower_tail = true);
double boost_qnorm(double p, double mean = 0.0, double sd = 1.0, 
                   bool lower_tail = true);
double boost_dnorm(double x, double mean = 0.0, double sd = 1.0);

double boost_plogis(double q, double location = 0.0, double scale = 1.0, 
                    bool lower_tail = true);
double boost_qlogis(double p, double location = 0.0, double scale = 1.0, 
                    bool lower_tail = true);
double boost_dlogis(double x, double location = 0.0, double scale = 1.0);

double boost_pextreme(double q, double location = 0.0, double scale = 1.0, 
                      bool lower_tail = true);
double boost_qextreme(double p, double location = 0.0, double scale = 1.0, 
                      bool lower_tail = true);
double boost_dextreme(double x, double location = 0.0, double scale = 1.0);

double boost_pchisq(double q, double df, bool lower_tail = true);
double boost_qchisq(double p, double df, bool lower_tail = true);

// --------------------------- Small utilities --------------------------------
inline double sq(double x) noexcept { return x * x; }

// seqcpp: inclusive sequence; inputs are int
std::vector<int> seqcpp(int start, int end);

// which: return indices of true values
std::vector<int> which(const std::vector<unsigned char>& vec);

// findInterval3: adapted helper (return indices following R-like convention)
std::vector<int> findInterval3(const std::vector<double>& x,
                               const std::vector<double>& v,
                               bool rightmost_closed = false,
                               bool all_inside = false,
                               bool left_open = false);

// --------------------------- Root finders -----------------------------------
double brent(const std::function<double(double)>& f,
             double x1, double x2, double tol = 1e-8, int maxiter = 100);
double bisect(const std::function<double(double)>& f,
              double x1, double x2, double tol = 1e-8, int maxiter = 100);

// --------------------------- Quantiles -------------------------------------
double quantilecpp(const std::vector<double>& x, double p);
double squantilecpp(const std::function<double(double)>& S, double p, double tol);

// subset_in_place: reorder/keep elements of v according to 'order' (indices)
template <typename T>
void subset_in_place(std::vector<T>& v, const std::vector<int>& order) {
  std::vector<T> temp_subset(order.size());
  int n = order.size();
  int nv = v.size();
  for (int i = 0; i < n; ++i) {
    int index = order[i];
    if (index < 0 || index >= nv) {
      throw std::out_of_range(
          "Index in 'order' is out of bounds for the source vector.");
    }
    temp_subset[i] = v[index];
  }
  v.assign(temp_subset.begin(), temp_subset.end());
}

// subset: return a subset of v according to 'order' (indices)
template <typename T>
std::vector<T> subset(const std::vector<T>& v, const std::vector<int>& order) {
  std::vector<T> result(order.size());
  int n = order.size();
  int nv = v.size();
  for (int i = 0; i < n; ++i) {
    int index = order[i];
    if (index < 0 || index >= nv) {
      throw std::out_of_range(
          "Index in 'order' is out of bounds for the source vector.");
    }
    result[i] = v[index];
  }
  return result;
}

// concat: concatenate two vectors
template <typename T>
std::vector<T> concat(const std::vector<T>& v1, const std::vector<T>& v2) {
  std::vector<T> result;
  result.reserve(v1.size() + v2.size());
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
std::vector<int> matchcpp(const std::vector<T>& x, const std::vector<T>& table, 
                          const int start_index = 0) {
  std::vector<int> result(x.size());
  int n = x.size();
  for (int i = 0; i < n; ++i) {
    auto it = std::find(table.begin(), table.end(), x[i]);
    if (it != table.end()) {
      result[i] = static_cast<int>(std::distance(table.begin(), it)) + start_index;
    } else {
      result[i] = -1;
    }
  }
  return result;
}

// bygroup: process grouping variables and return indices and lookup tables
ListCpp bygroup(const DataFrameCpp& data, const std::vector<std::string>& variables);

// --------------------------- Matrix utilities (FlatMatrix) ------------------
std::vector<double> mat_vec_mult(const FlatMatrix& A, const std::vector<double>& x);
FlatMatrix mat_mat_mult(const FlatMatrix& A, const FlatMatrix& B);
FlatMatrix transpose(const FlatMatrix& A);
double quadsym(const std::vector<double>& u, const FlatMatrix& v);

// --------------------------- Linear algebra helpers (FlatMatrix-backed) ----
int cholesky2(FlatMatrix& matrix, int n, double toler = 1e-12);
void chsolve2(FlatMatrix& matrix, int n, std::vector<double>& y);
void chinv2(FlatMatrix& matrix, int n);
FlatMatrix invsympd(const FlatMatrix& matrix, int n, double toler = 1e-12);

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
                 bool lowerend, double initialend);




#include <iostream>
#include <iomanip>
#include <sstream>
#include <type_traits>

// Print a std::vector<T> to std::cout.
// Requirements: T must be streamable via operator<< to std::ostream.
//
// Parameters:
//  - v: vector to print
//  - label: optional prefix printed before the vector
//  - precision: if >= 0, sets std::fixed and std::setprecision(precision) 
//    for floating values
//  - head: number of leading elements to show when truncated
//  - tail: number of trailing elements to show when truncated
//  - sep: separator between elements (default ", ")
//  - show_indices: if true prints each element as "idx: value"
//  - endline: whether to append a newline at the end (true by default)
template <typename T>
void print_vector(const std::vector<T>& v,
                  const std::string& label = "",
                  int precision = -1,
                  std::size_t head = 5,
                  std::size_t tail = 5,
                  const std::string& sep = ", ",
                  bool show_indices = false,
                  bool endline = true) {
  static_assert(
    std::is_convertible<decltype(
      std::declval<std::ostream&>() << std::declval<T>()), std::ostream&>::value,
                   "Type T must be streamable to std::ostream (operator<<)");
  std::ostringstream ss;
  if (!label.empty()) ss << label << ": ";
  
  std::size_t n = v.size();
  if (n == 0) {
    ss << "[]";
    if (endline) ss << '\n';
    std::cout << ss.str();
    return;
  }
  
  // Configure precision only if requested
  bool use_precision = (precision >= 0);
  if (use_precision) ss << std::fixed << std::setprecision(precision);
  
  ss << "[";
  auto print_elem = [&](std::size_t i) {
    if (show_indices) ss << i << ": ";
    ss << v[i];
  };
  
  if (n <= head + tail || head + tail == 0) {
    for (std::size_t i = 0; i < n; ++i) {
      if (i) ss << sep;
      print_elem(i);
    }
  } else {
    // print head
    for (std::size_t i = 0; i < head; ++i) {
      if (i) ss << sep;
      print_elem(i);
    }
    ss << sep << "..." << sep;
    // print tail
    for (std::size_t j = n - tail; j < n; ++j) {
      if (j != n - tail) ss << sep;
      print_elem(j);
    }
  }
  ss << "]";
  if (endline) ss << '\n';
  
  std::cout << ss.str();
}


#endif // __UTILITIES_H__