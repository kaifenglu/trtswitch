// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>    // RcppParallel::Worker, parallelFor
#include <RcppThread.h>      // RcppThread::Rcerr

#include "survival_analysis.h"
#include "utilities.h"
#include "dataframe_list.h"  // FlatMatrix, IntMatrix, DataFrameCpp, ListCpp
#include "thread_utils.h"    // push_thread_warning / drain_thread_warnings_to_R

#include <vector>
#include <string>
#include <numeric>   // iota, inner_product
#include <cmath>     // isnan, isinf, fabs, NaN, exp, log
#include <stdexcept> // exceptions
#include <algorithm> // sort, none_of, any_of
#include <random>
#include <cctype>


// Helper function to compute confidence interval for survival probability
std::vector<double> fsurvci(double surv, double sesurv, std::string& ct, double z) {
  double grad, hw, lower = NaN, upper = NaN;
  if (surv == 1.0 && sesurv == 0.0) {
    lower = upper = 1.0;
  } else if (ct == "plain" || ct == "linear") {
    lower = std::max(surv - z * sesurv, 0.0);
    upper = std::min(surv + z * sesurv, 1.0);
  } else if (ct == "log") {
    grad = 1.0 / surv;
    hw = z * grad * sesurv;
    lower = std::exp(std::log(surv) - hw);
    upper = std::min(std::exp(std::log(surv) + hw), 1.0);
  } else if (ct == "log-log" || ct == "loglog" || ct == "cloglog") {
    grad = 1.0 / (surv * std::log(surv));
    hw = z * grad * sesurv;
    lower = std::exp(-std::exp(std::log(-std::log(surv)) - hw));
    upper = std::exp(-std::exp(std::log(-std::log(surv)) + hw));
  } else if (ct == "logit") {
    grad = 1.0 / (surv * (1.0 - surv));
    hw = z * grad * sesurv;
    lower = boost_plogis(boost_qlogis(surv) - hw);
    upper = boost_plogis(boost_qlogis(surv) + hw);
  } else if (ct == "arcsin" || ct == "asin" || ct == "asinsqrt") {
    grad = 1.0 / (2.0 * std::sqrt(surv * (1.0 - surv)));
    hw = z * grad * sesurv;
    lower = sq(std::sin(std::asin(std::sqrt(surv)) - hw));
    upper = sq(std::sin(std::asin(std::sqrt(surv)) + hw));
  } else {
    throw std::invalid_argument("Unknown confidence type: " + ct);
  }

  return {lower, upper};
}


// Compute survival quantiles and confidence intervals
DataFrameCpp survQuantilecpp(const std::vector<double>& time,
                             const std::vector<int>& event,
                             const double cilevel,
                             const std::string& transform,
                             const std::vector<double>& probs) {
  
  // Basic input checks
  if (time.size() != event.size()) {
    throw std::invalid_argument("time and event must have the same length");
  }
  int n = static_cast<int>(time.size());
  for (int i = 0; i < n; ++i) {
    if (std::isnan(time[i])) throw std::invalid_argument("time must be provided");
  }
  for (int i = 0; i < n; ++i) {
    // event is integer; check not NaN is irrelevant; but check allowed values:
    if (!(event[i] == 0 || event[i] == 1)) {
      throw std::invalid_argument("event must be 1 or 0");
    }
  }
  
  int total_events = 0;
  for (int v : event) if (v == 1) ++total_events;
  if (total_events == 0) throw std::invalid_argument("at least 1 event is needed");
  
  if (!(cilevel > 0.0 && cilevel < 1.0)) {
    throw std::invalid_argument("cilevel must lie between 0 and 1");
  }
  
  std::string ct = transform;
  std::for_each(ct.begin(), ct.end(), [](char & c) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });
  
  if (!(ct == "linear" || ct == "plain" || ct == "log" ||
      ct == "loglog" || ct == "log-log" || ct == "cloglog" ||
      ct == "asinsqrt" || ct == "arcsin"|| ct == "asin" ||
      ct == "logit")) {
    throw std::invalid_argument("Invalid value for transform");
  }
  
  int code;
  if (ct == "plain" || ct == "linear") code = 1;
  else if (ct == "log") code = 2;
  else if (ct == "log-log" || ct == "loglog" || ct == "cloglog") code = 3;
  else if (ct == "logit") code = 4;
  else if (ct == "arcsin" || ct == "asin" || ct == "asinsqrt") code = 5;
  else throw std::invalid_argument("Unknown confidence type: " + ct);
  
  // probs checks
  for (double p : probs) {
    if (std::isnan(p)) throw std::invalid_argument("probs must be provided");
    if (!(p > 0.0 && p < 1.0)) 
      throw std::invalid_argument("Elements of probs must lie between 0 and 1");
  }

  // sort by time, and event with event in descending order
  std::vector<int> order(n);
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&](int i, int j) {
    if (time[i] != time[j]) return time[i] < time[j];
    return event[i] > event[j];
  });
  
  std::vector<double> time2(n);
  std::vector<int> event2(n);
  for (int i = 0; i < n; ++i) {
    time2[i] = time[order[i]];
    event2[i] = event[order[i]];
  }
  
  std::vector<double> time0, nrisk0, nevent0, surv0, sesurv0;
  time0.reserve(total_events);
  nrisk0.reserve(total_events);
  nevent0.reserve(total_events);
  surv0.reserve(total_events);
  sesurv0.reserve(total_events);
  
  double t = 0, nrisk = n, nevent = 0, surv = 1, vcumhaz = 0, sesurv;
  
  bool cache = false; // buffer for the current event time
  for (int j=0; j<n; ++j) {
    if ((j == 0 && event2[j] == 1) ||
        (j >= 1 && event2[j] == 1 && time2[j] > time2[j-1])) {
      // new event, add the info for the previous event
      if (cache) {
        surv *= (1.0 - nevent / nrisk);
        if (nrisk > nevent) {
          vcumhaz += nevent / (nrisk * (nrisk - nevent));
        } else {
          vcumhaz = NaN;
        }
        sesurv = surv * std::sqrt(vcumhaz);
        
        time0.push_back(t);
        nrisk0.push_back(nrisk);
        nevent0.push_back(nevent);
        surv0.push_back(surv);
        sesurv0.push_back(sesurv);
      }
      
      // update the buffer for the current event time
      t = time2[j];
      nrisk = n-j;
      nevent = 1;
      cache = true;
    } else if (j >= 1 && event2[j] == 1 && event2[j-1] == 1 && time2[j] == time2[j-1]) {
      // tied event
      ++nevent;
    } else if (j >= 1 && event2[j] == 0 && event2[j-1] == 1) {
      // new censoring, add the info for the previous event
      surv *= (1.0 - nevent / nrisk);
      if (nrisk > nevent) {
        vcumhaz += nevent / (nrisk * (nrisk - nevent));
      } else {
        vcumhaz = NaN;
      }
      sesurv = surv * std::sqrt(vcumhaz);
      
      time0.push_back(t);
      nrisk0.push_back(nrisk);
      nevent0.push_back(nevent);
      surv0.push_back(surv);
      sesurv0.push_back(sesurv);
      
      // empty the buffer for the current event time
      cache = false;
    }
  }
  
  // add the info for the last event
  if (cache) {
    surv *= (1.0 - nevent / nrisk);
    if (nrisk > nevent) {
      vcumhaz += nevent / (nrisk * (nrisk - nevent));
    } else {
      vcumhaz = NaN;
    }
    sesurv = surv * std::sqrt(vcumhaz);
    
    time0.push_back(t);
    nrisk0.push_back(nrisk);
    nevent0.push_back(nevent);
    surv0.push_back(surv);
    sesurv0.push_back(sesurv);
  }
  
  int n0 = time0.size();

  std::vector<double> z(n0, NaN), grad(n0, NaN);
  double zcrit = boost_qnorm((1.0 + cilevel) /2.0);
  
  int m = probs.size();
  std::vector<double> quantile(m), lower(m), upper(m);
  for (int j=0; j<m; ++j) {
    double q = 1.0 - probs[j];
    for (int i=0; i<n0; ++i) {
      if (nrisk0[i] > nevent0[i]) {
        switch (code) {
        case 1: // linear or plain
          z[i] = (surv0[i] - q) / sesurv0[i];
          break;
          
        case 2: // log
          grad[i] = 1.0 / surv0[i];
          z[i] = (std::log(surv0[i]) - std::log(q)) / (grad[i] * sesurv0[i]);
          break;
          
        case 3: // loglog / log-log / cloglog
          grad[i] = 1.0 / (surv0[i] * std::log(surv0[i]));
          z[i] = (std::log(-std::log(surv0[i])) - std::log(-std::log(q))) / (grad[i] * sesurv0[i]);
          break;
          
        case 4: // logit
          grad[i] = 1.0 / (surv0[i] * (1.0 - surv0[i]));
          z[i] = (boost_qlogis(surv0[i]) - boost_qlogis(q)) / (grad[i] * sesurv0[i]);
          break;
          
        case 5: // arcsin / asin / asinsqrt
          grad[i] = 1.0 / (2.0 * std::sqrt(surv0[i] * (1.0 - surv0[i])));
          z[i] = (std::asin(std::sqrt(surv0[i])) - std::asin(std::sqrt(q))) / (grad[i] * sesurv0[i]);
          break;
        }
      }
    }
    
   // find indices with |z| <= zcrit (we only need min and max)
    int imin = -1, imax = -1;
    for (int i = 0; i < n0; ++i) {
      if (!std::isnan(z[i]) && std::fabs(z[i]) <= zcrit) {
        if (imin == -1) imin = i;
        imax = i;
      }
    }
    if (imin == -1) {
      lower[j] = upper[j] = NaN;
    } else {
      lower[j] = time0[imin];
      if (imax < n0 - 1) upper[j] = time0[imax + 1];
      else upper[j] = NaN;
    }
    
    // quantile: first time where surv0 < q. surv0 is non-increasing, do binary search
    int first_lt = -1;
    int lo = 0, hi = n0 - 1;
    while (lo <= hi) {
      int mid = (lo + hi) >> 1;
      if (surv0[mid] < q) {
        first_lt = mid;
        hi = mid - 1;
      } else {
        lo = mid + 1;
      }
    }
    if (first_lt >= 0) quantile[j] = time0[first_lt];
    else quantile[j] = NaN;
  }
  
  DataFrameCpp result;
  result.push_back(probs, "prob");
  result.push_back(quantile, "quantile");
  result.push_back(lower, "lower");
  result.push_back(upper, "upper");
  result.push_back(cilevel, "cilevel");
  result.push_back(ct, "transform");
  return result;
}

//' @title Brookmeyer-Crowley Confidence Interval for Quantiles of
//' Right-Censored Time-to-Event Data
//' @description Obtains the Brookmeyer-Crowley confidence
//' interval for quantiles of right-censored time-to-event data.
//'
//' @param time The vector of possibly right-censored survival times.
//' @param event The vector of event indicators.
//' @param cilevel The confidence interval level. Defaults to 0.95.
//' @param transform The transformation of the survival function to use
//'   to construct the confidence interval. Options include
//'   "linear" (alternatively "plain"), "log",
//'   "loglog" (alternatively "log-log" or "cloglog"),
//'   "asinsqrt" (alternatively "asin" or "arcsin"), and "logit".
//'   Defaults to "loglog".
//'
//' @param probs The vector of probabilities to calculate the quantiles.
//'   Defaults to c(0.25, 0.5, 0.75).
//'
//' @return A data frame containing the estimated quantile and
//' confidence interval corresponding to each specified probability.
//' It includes the following variables:
//'
//' * \code{prob}: The probability to calculate the quantile.
//'
//' * \code{quantile}: The estimated quantile.
//'
//' * \code{lower}: The lower limit of the confidence interval.
//'
//' * \code{upper}: The upper limit of the confidence interval.
//'
//' * \code{cilevel}: The confidence interval level.
//'
//' * \code{transform}: The transformation of the survival function to use
//'   to construct the confidence interval.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' survQuantile(
//'   time = c(33.7, 3.9, 10.5, 5.4, 19.5, 23.8, 7.9, 16.9, 16.6,
//'            33.7, 17.1, 7.9, 10.5, 38),
//'   event = c(0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1),
//'   probs = c(0.25, 0.5, 0.75))
//'
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame survQuantile(const Rcpp::NumericVector& time,
                             const Rcpp::IntegerVector& event,
                             const double cilevel = 0.95,
                             const std::string& transform = "loglog",
                             const Rcpp::NumericVector& probs = 
                               Rcpp::NumericVector::create(0.25, 0.5, 0.75)) {
  
  std::vector<double> time_vec = Rcpp::as<std::vector<double>>(time);
  std::vector<int> event_vec = Rcpp::as<std::vector<int>>(event);
  std::vector<double> probs_vec = Rcpp::as<std::vector<double>>(probs);
  
  DataFrameCpp result = survQuantilecpp(time_vec, event_vec, cilevel, transform, probs_vec);
  
  return Rcpp::wrap(result);
}


// Compute Kaplan-Meier estimator
DataFrameCpp kmestcpp(const DataFrameCpp& data,
                      const std::vector<std::string>& stratum,
                      const std::string& time,
                      const std::string& time2,
                      const std::string& event,
                      const std::string& weight,
                      const std::string& conftype,
                      const double conflev,
                      const bool keep_censor) {
  
  int n = data.nrows();
  
  bool has_stratum;
  std::vector<int> stratumn(n);
  DataFrameCpp u_stratum;
  int p_stratum = stratum.size();
  if (p_stratum == 0 || (p_stratum == 1 && stratum[0] == "")) {
    has_stratum = false;
  } else {
    ListCpp out = bygroup(data, stratum);
    has_stratum = true;
    stratumn = out.get<std::vector<int>>("index");
    u_stratum = out.get<DataFrameCpp>("lookup");
  }
  
  if (!data.containElementNamed(time))
    throw std::invalid_argument("data must contain the time variable");
  
  std::vector<double> timen = data.get<double>(time);
  for (double v : timen) if (v < 0.0)
    throw std::invalid_argument("time must be nonnegative for each subject");
  
  bool has_time2 = !time2.empty() && data.containElementNamed(time2);
  std::vector<double> time2n(n);
  if (has_time2) {
    time2n = data.get<double>(time2);
    for (int i = 0; i < n; ++i) {
      if (time2n[i] <= timen[i]) {
        throw std::invalid_argument("time2 must be greater than time for each observation");
      }
    }
  }
  
  // unify right censored data with counting process data
  std::vector<double> tstartn(n), tstopn(n);
  if (!has_time2) {
    tstopn = std::move(timen);
    // tstartn remains zero-initialized
  } else {
    tstartn = std::move(timen);
    tstopn = std::move(time2n);
  }
  
  if (!data.containElementNamed(event))
    throw std::invalid_argument("data must contain the event variable");
  
  // event -> numeric 0/1
  std::vector<double> eventn(n);
  if (data.bool_cols.count(event)) {
    const std::vector<unsigned char>& vb = data.get<unsigned char>(event);
    for (int i = 0; i < n; ++i) eventn[i] = vb[i] ? 1.0 : 0.0;
  } else if (data.int_cols.count(event)) {
    const std::vector<int>& vi = data.get<int>(event);
    for (int i = 0; i < n; ++i) eventn[i] = static_cast<double>(vi[i]);
  } else if (data.numeric_cols.count(event)) {
    eventn = data.get<double>(event);
  } else {
    throw std::invalid_argument("event variable must be bool, integer or numeric");
  }
  for (double val : eventn) if (val != 0 && val != 1)
    throw std::invalid_argument("event must be 1 or 0 for each observation");
  
  std::vector<double> weightn(n, 1.0);
  if (!weight.empty() && data.containElementNamed(weight)) {
    weightn = data.get<double>(weight);
    for (double v : weightn) if (v <= 0.0)
      throw std::invalid_argument("weight must be greater than 0");
  }
  
  std::string ct = conftype;
  std::for_each(ct.begin(), ct.end(), [](char & c) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });
  
  if (!(ct == "none" || ct == "plain" || ct == "log" || ct == "log-log" || 
      ct == "logit" || ct == "arcsin")) {
    throw std::invalid_argument("conftype must be none, plain, log, log-log, logit, or arcsin");
  }
  
  if (conflev <= 0.0 || conflev >= 1.0) {
    throw std::invalid_argument("conflev must lie between 0 and 1");
  }
  
  // confidence interval for survival probability
  double z = boost_qnorm((1.0 + conflev) / 2.0);
  
  std::vector<int> stratum1(n), size1(n);
  std::vector<double> time1(n), nrisk1(n), nevent1(n), ncensor1(n);
  std::vector<double> surv1(n), sesurv1(n);
  std::vector<double> lower1(n), upper1(n);
  
  // sort by stopping time in descending order within each stratum
  std::vector<int> order(n);
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&](int i, int j) {
    if (stratumn[i] != stratumn[j]) return stratumn[i] > stratumn[j];
    if (tstopn[i] != tstopn[j]) return tstopn[i] > tstopn[j];
    return eventn[i] < eventn[j];
  });
  
  subset_in_place(stratumn, order);
  subset_in_place(tstartn, order);
  subset_in_place(tstopn, order);
  subset_in_place(eventn, order);
  subset_in_place(weightn, order);
  
  // sort by starting time in descending order within each stratum
  std::vector<int> order1(n);
  std::iota(order1.begin(), order1.end(), 0);
  std::sort(order1.begin(), order1.end(), [&](int i, int j) {
    if (stratumn[i] != stratumn[j]) return stratumn[i] > stratumn[j];
    return tstartn[i] > tstartn[j];
  });
  
  double nrisk = 0;   // number at risk
  double nriskw = 0;  // weighted # at risk
  double nriskw2 = 0; // weight^2 # at risk
  double nevent = 0;  // number of events
  double neventw = 0; // weighted # events
  double ncensor = 0; // number of censored
  double surv = 1;    // survival probability
  double sesurv = 0;  // standard error of survival probability
  double vcumhaz = 0; // cumulative hazard variance
  
  int istratum = stratumn.empty() ? 0 : stratumn[0]; // current stratum
  int i1 = 0;                 // index for removing out-of-risk subjects
  
  int stratum_size = 0;
  std::vector<int> stratum0(n), size00(n);
  std::vector<double> time0(n), nrisk0(n), nriskw0(n), nriskw20(n);
  std::vector<double> nevent0(n), neventw0(n), ncensor0(n);
  
  int index = 0;
  for (int i = 0; i < n; ) {
    // Reset when entering a new stratum
    if (stratumn[i] != istratum) {
      istratum = stratumn[i];
      i1 = i;
      stratum_size = 0;
      nrisk = nriskw = nriskw2 = 0;
    }
    
    const double dtime = tstopn[i];
    
    // Process all persons tied at this dtime
    while (i < n && tstopn[i] == dtime) {
      ++nrisk;
      nriskw += weightn[i];
      nriskw2 += weightn[i] * weightn[i];
      
      if (eventn[i] == 1) {
        ++nevent;
        neventw += weightn[i];
      } else {
        ++ncensor;
      }
      
      if (tstartn[i] == 0) ++stratum_size; // unique subjects in the stratum
      
      ++i;
      
      if (i < n && stratumn[i] != istratum) {
        size00[istratum] = stratum_size; // update size of the stratum
        break;
      }
    }
    
    // remove subjects no longer at risk
    for (; i1 < n; ++i1) {
      const int p1 = order1[i1];
      if (tstartn[p1] < dtime || stratumn[p1] != istratum) break;
      
      --nrisk;
      nriskw -= weightn[p1];
      nriskw2 -= weightn[p1] * weightn[p1];
    }
    
    if (nevent > 0 || keep_censor) {
      stratum0[index] = istratum;
      time0[index] = dtime;
      nrisk0[index] = nrisk;
      nriskw0[index] = nriskw;
      nriskw20[index] = nriskw2;
      nevent0[index] = nevent;
      neventw0[index] = neventw;
      ncensor0[index] = ncensor;
      
      nevent = neventw = 0;
      ++index;
    }
    
    ncensor = 0;
  }
  
  
  size00[istratum] = stratum_size; // update size of the last stratum
  
  std::vector<int> size0(index);
  std::vector<double> surv0(index), sesurv0(index);
  std::vector<double> lower0(index), upper0(index);
  
  if (index > 0) {
    istratum = stratum0[index - 1];
    for (int i = index - 1; i >= 0; --i) {
      double nevent_l = nevent0[i];
      double neventw_l = neventw0[i];
      double nriskw_l = nriskw0[i];
      double nriskw2_l = nriskw20[i];
      
      if (stratum0[i] != istratum) { // hit a new stratum
        // reset temporary variables
        istratum = stratum0[i];
        surv = 1.0; vcumhaz = 0.0; sesurv = 0.0;
      }
      
      if (nevent_l > 0) {
        double p = neventw_l / nriskw_l;
        surv *= (1.0 - p);
        double m = nriskw_l * nriskw_l / nriskw2_l;
        vcumhaz += p / (m * (1.0 - p));
        sesurv = surv * std::sqrt(vcumhaz);
      }
      
      size0[i] = size00[istratum];
      surv0[i] = surv;
      sesurv0[i] = sesurv;
      if (ct != "none") {
        std::vector<double> ci = fsurvci(surv, sesurv, ct, z);
        lower0[i] = ci[0];
        upper0[i] = ci[1];
      }
    }
  }
  
  // update the global index (reverse)
  for (int i = 0; i < index; ++i) {
    int k = index - 1 - i;
    stratum1[i] = stratum0[k];
    size1[i] = size0[k];
    time1[i] = time0[k];
    nrisk1[i] = nrisk0[k];
    nevent1[i] = nevent0[k];
    ncensor1[i] = ncensor0[k];
    surv1[i] = surv0[k];
    sesurv1[i] = sesurv0[k];
    if (ct != "none") {
      lower1[i] = lower0[k];
      upper1[i] = upper0[k];
    }
  }
  
  // subset to the used portion
  std::vector<int> sub(index);
  std::iota(sub.begin(), sub.end(), 0);
  subset_in_place(stratum1, sub);
  subset_in_place(size1, sub);
  subset_in_place(time1, sub);
  subset_in_place(nrisk1, sub);
  subset_in_place(nevent1, sub);
  subset_in_place(ncensor1, sub);
  subset_in_place(surv1, sub);
  subset_in_place(sesurv1, sub);
  if (ct != "none") {
    subset_in_place(lower1, sub);
    subset_in_place(upper1, sub);
  }
  
  DataFrameCpp result;
  result.push_back(size1, "size");
  result.push_back(time1, "time");
  result.push_back(nrisk1, "nrisk");
  result.push_back(nevent1, "nevent");
  result.push_back(ncensor1, "ncensor");
  result.push_back(surv1, "surv");
  result.push_back(sesurv1, "sesurv");
  
  if (ct != "none") {
    result.push_back(lower1, "lower");
    result.push_back(upper1, "upper");
    result.push_back(conflev, "conflev");
    result.push_back(ct, "conftype");
  }
  
  // add stratum lookup columns (if any)
  if (has_stratum) {
    for (int i = 0; i < p_stratum; ++i) {
      std::string s = stratum[i];
      if (u_stratum.int_cols.count(s)) {
        auto v = u_stratum.get<int>(s);
        subset_in_place(v, stratum1);
        result.push_back(v, s);
      } else if (u_stratum.numeric_cols.count(s)) {
        auto v = u_stratum.get<double>(s);
        subset_in_place(v, stratum1);
        result.push_back(v, s);
      } else if (u_stratum.string_cols.count(s)) {
        auto v = u_stratum.get<std::string>(s);
        subset_in_place(v, stratum1);
        result.push_back(v, s);
      } else {
        throw std::invalid_argument("unsupported type for stratum variable " + s);
      }
    }
  }
  
  return result;
}

// ----------------------------------------------------------------------------
// Parallel worker
// ----------------------------------------------------------------------------

struct kmestWorker : public RcppParallel::Worker {
  const std::vector<DataFrameCpp>* data_ptr;
  const std::vector<std::string>& stratum;
  const std::string& time;
  const std::string& time2;
  const std::string& event;
  const std::string& weight;
  const std::string& conftype;
  const double conflev;
  const bool keep_censor;
  
  std::vector<DataFrameCpp>* results;
  
  kmestWorker(const std::vector<DataFrameCpp>* data_ptr_,
              const std::vector<std::string>& stratum_,
              const std::string& time_,
              const std::string& time2_,
              const std::string& event_,
              const std::string& weight_,
              const std::string& conftype_,
              const double conflev_,
              const bool keep_censor_,
              std::vector<DataFrameCpp>* results_)
    : data_ptr(data_ptr_), stratum(stratum_), time(time_), 
      time2(time2_), event(event_), weight(weight_), conftype(conftype_), 
      conflev(conflev_), keep_censor(keep_censor_), results(results_) {}
  
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; ++i) {
      // Call the pure C++ function kmestcpp on data_ptr->at(i)
      DataFrameCpp out = kmestcpp(
        (*data_ptr)[i], stratum, time, time2, event, weight, 
        conftype, conflev, keep_censor);
      (*results)[i] = std::move(out);
    }
  }
};

// [[Rcpp::export]]
Rcpp::List kmestRcpp(SEXP data,
                     const std::vector<std::string>& stratum,
                     const std::string& time = "time",
                     const std::string& time2 = "",
                     const std::string& event = "event",
                     const std::string& weight = "",
                     const std::string& conftype = "log-log",
                     const double conflev = 0.95,
                     const bool keep_censor = false) {
  
  // Case A: single data.frame -> call kmestcpp on main thread
  if (Rf_inherits(data, "data.frame")) {
    Rcpp::DataFrame rdf(data);
    DataFrameCpp dfcpp = convertRDataFrameToCpp(rdf);
    
    // Call core C++ function directly on the DataFrameCpp
    DataFrameCpp cpp_result = kmestcpp(
      dfcpp, stratum, time, time2, event, weight, 
      conftype, conflev, keep_censor
    );
    
    thread_utils::drain_thread_warnings_to_R();
    return Rcpp::wrap(cpp_result);
  }
  
  // Case B: list of data.frames -> process in parallel
  if (TYPEOF(data) == VECSXP) {
    Rcpp::List lst(data);
    int m = lst.size();
    if (m == 0) return Rcpp::List(); // nothing to do
    
    // Convert each element to DataFrameCpp.
    std::vector<DataFrameCpp> data_vec;
    data_vec.reserve(m);
    for (int i = 0; i < m; ++i) {
      SEXP el = lst[i];
      if (!Rf_inherits(el, "data.frame")) {
        Rcpp::stop("When 'data' is a list, every element must be a data.frame.");
      }
      Rcpp::DataFrame rdf(el);
      data_vec.emplace_back(convertRDataFrameToCpp(rdf));
    }
    
    // Pre-allocate result vector of C++ objects (no R API used inside worker threads)
    std::vector<DataFrameCpp> results(m);
    
    // Build worker and run parallelFor across all indices [0, m)
    kmestWorker worker(
        &data_vec, stratum, time, time2, event, weight, 
        conftype, conflev, keep_censor, &results
    );
    
    // Execute parallelFor (this will schedule work across threads)
    RcppParallel::parallelFor(0, m, worker);
    
    // Drain thread-collected warnings (on main thread) into R's warning system
    thread_utils::drain_thread_warnings_to_R();
    
    // Convert C++ DataFrameCpp results back to R on the main thread
    Rcpp::List out(m);
    for (int i = 0; i < m; ++i) {
      out[i] = Rcpp::wrap(results[i]);
    }
    return out;
  }
  
  // Neither a data.frame nor a list: error
  Rcpp::stop("Input 'data' must be either a data.frame or a list of data.frames.");
  return R_NilValue; // unreachable
}


// Compute difference in KM estimates at a milestone time
DataFrameCpp kmdiffcpp(const DataFrameCpp& data,
                       const std::vector<std::string>& stratum,
                       const std::string& treat,
                       const std::string& time,
                       const std::string& time2,
                       const std::string& event,
                       const std::string& weight,
                       const double milestone,
                       const double survDiffH0,
                       const double conflev) {
  int n = data.nrows();

  bool has_stratum;
  std::vector<int> stratumn(n);
  DataFrameCpp u_stratum;
  int p_stratum = stratum.size();
  if (p_stratum == 0 || (p_stratum == 1 && stratum[0] == "")) {
    has_stratum = false;
  } else {
    ListCpp out = bygroup(data, stratum);
    has_stratum = true;
    stratumn = out.get<std::vector<int>>("index");
    u_stratum = out.get<DataFrameCpp>("lookup");
  }
  
  if (!data.containElementNamed(treat)) {
    throw std::invalid_argument("data must contain the treat variable");
  }
  
  // create the numeric treat variable
  std::vector<int> treatn(n);
  std::vector<int> treatwi;
  std::vector<double> treatwn;
  std::vector<std::string> treatwc;
  if (data.bool_cols.count(treat) || data.int_cols.count(treat)) {
    std::vector<int> treatv(n);
    if (data.bool_cols.count(treat)) {
      const std::vector<unsigned char>& treatvb = data.get<unsigned char>(treat);
      for (int i = 0; i < n; ++i) treatv[i] = treatvb[i] ? 1 : 0;
    } else treatv = data.get<int>(treat);
    treatwi = unique_sorted(treatv); // obtain unique treatment values
    if (treatwi.size() != 2)
      throw std::invalid_argument("treat must have two and only two distinct values");
    if (std::all_of(treatwi.begin(), treatwi.end(), [](int v) { return v == 0 || v == 1; })) {
      treatwi = {1, 0}; // special handling for 1/0 treatment coding
      for (int i = 0; i < n; ++i) treatn[i] = 2 - treatv[i];
    } else {
      treatn = matchcpp(treatv, treatwi, 1);
    }
  } else if (data.numeric_cols.count(treat)) {
    const std::vector<double>& treatv = data.get<double>(treat);
    treatwn = unique_sorted(treatv);
    if (treatwn.size() != 2)
      throw std::invalid_argument("treat must have two and only two distinct values");
    if (std::all_of(treatwn.begin(), treatwn.end(), [](double v) { return v == 0.0 || v == 1.0; })) {
      treatwn = {1.0, 0.0};
      for (int i = 0; i < n; ++i) treatn[i] = 2 - static_cast<int>(treatv[i]);
    } else {
      treatn = matchcpp(treatv, treatwn, 1);
    }
  } else if (data.string_cols.count(treat)) {
    const std::vector<std::string>& treatv = data.get<std::string>(treat);
    treatwc = unique_sorted(treatv);
    if (treatwc.size() != 2) 
      throw std::invalid_argument("treat must have two and only two distinct values");
    treatn = matchcpp(treatv, treatwc, 1);
  } else {
    throw std::invalid_argument("incorrect type for the treat variable in the input data");
  }
  
  if (!data.containElementNamed(time)) {
    throw std::invalid_argument("data must contain the time variable");
  }
  std::vector<double> timen = data.get<double>(time);
  for (double v : timen) if (v < 0.0)
    throw std::invalid_argument("time must be nonnegative for each subject");
  
  bool has_time2 = !time2.empty() && data.containElementNamed(time2);
  std::vector<double> time2n(n);
  if (has_time2) {
    time2n = data.get<double>(time2);
    for (int i = 0; i < n; ++i) {
      if (time2n[i] <= timen[i]) {
        throw std::invalid_argument("time2 must be greater than time for each observation");
      }
    }
  }
  
  // unify right censored data with counting process data
  std::vector<double> tstartn(n), tstopn(n);
  if (!has_time2) {
    tstopn = std::move(timen);
    // tstartn remains zero-initialized
  } else {
    tstartn = std::move(timen);
    tstopn = std::move(time2n);
  }
  
  if (!data.containElementNamed(event))
    throw std::invalid_argument("data must contain the event variable");
  
  // event -> numeric 0/1
  std::vector<double> eventn(n);
  if (data.bool_cols.count(event)) {
    const std::vector<unsigned char>& vb = data.get<unsigned char>(event);
    for (int i = 0; i < n; ++i) eventn[i] = vb[i] ? 1.0 : 0.0;
  } else if (data.int_cols.count(event)) {
    const std::vector<int>& vi = data.get<int>(event);
    for (int i = 0; i < n; ++i) eventn[i] = static_cast<double>(vi[i]);
  } else if (data.numeric_cols.count(event)) {
    eventn = data.get<double>(event);
  } else {
    throw std::invalid_argument("event variable must be bool, integer or numeric");
  }
  for (double val : eventn) if (val != 0 && val != 1)
    throw std::invalid_argument("event must be 1 or 0 for each observation");
  
  std::vector<double> weightn(n, 1.0);
  if (!weight.empty() && data.containElementNamed(weight)) {
    weightn = data.get<double>(weight);
    for (double v : weightn) if (v <= 0.0)
      throw std::invalid_argument("weight must be greater than 0");
  }
  
  if (std::isnan(milestone)) {
    throw std::invalid_argument("milestone must be provided");
  }
  
  if (milestone <= 0) {
    throw std::invalid_argument("milestone must be positive");
  }
  
  if (survDiffH0 <= -1 || survDiffH0 >= 1) {
    throw std::invalid_argument("survDiffH0 must lie between -1 and 1");
  }
  
  if (conflev <= 0 || conflev >= 1) {
    throw std::invalid_argument("conflev must lie between 0 and 1");
  }
  
  bool noerr = true;
  
  // sort by stratum
  std::vector<int> order(n);
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&](int i, int j) {
    return stratumn[i] < stratumn[j];
  });
  
  subset_in_place(stratumn, order);
  subset_in_place(treatn, order);
  subset_in_place(tstartn, order);
  subset_in_place(tstopn, order);
  subset_in_place(eventn, order);
  subset_in_place(weightn, order);
  
  // identify the locations of the unique values of stratum
  std::vector<int> idx(1,0);
  for (int i=1; i<n; ++i) {
    if (stratumn[i] != stratumn[i-1]) {
      idx.push_back(i);
    }
  }
  
  int nstrata = idx.size();
  idx.push_back(n);
  
  // whether the milestone exceeds the largest observed time
  for (int i=0; i<nstrata; ++i) {
    int n1 = idx[i+1] - idx[i];
    std::vector<int> q1 = seqcpp(idx[i], idx[i+1]-1);
    std::vector<int> treat1 = subset(treatn, q1);
    std::vector<double> tstop1 = subset(tstopn, q1);
    std::vector<double> tstop11, tstop12;
    for (int i=0; i<n1; ++i) {
      if (treat1[i] == 1) {
        tstop11.push_back(tstop1[i]);
      } else {
        tstop12.push_back(tstop1[i]);
      }
    }
    
    double max_time11 = *std::max_element(tstop11.begin(), tstop11.end());
    double max_time12 = *std::max_element(tstop12.begin(), tstop12.end());
    if (milestone > std::min(max_time11, max_time12)) {
      std::string stratumerr;
      if (!has_stratum) {
        stratumerr = "";
      } else {
        for (int j = 0; j < p_stratum; ++j) {
          std::string s = stratum[j];
          if (u_stratum.int_cols.count(s)) {
            auto v = u_stratum.get<int>(s);
            stratumerr += " " + s + " = " + std::to_string(v[i]);
          } else if (u_stratum.numeric_cols.count(s)) {
            auto v = u_stratum.get<double>(s);
            stratumerr += " " + s + " = " + std::to_string(v[i]);
          } else if (u_stratum.string_cols.count(s)) {
            auto v = u_stratum.get<std::string>(s);
            stratumerr += " " + s + " = " + v[i];
          } else {
            throw std::invalid_argument("unsupported type for stratum variable " + s);
          }
        }
      }
      
      int k = milestone > max_time11 ? 0 : 1;
      std::string treaterr;
      if (data.bool_cols.count(treat) || data.int_cols.count(treat)) {
        treaterr = " " + treat + " = " + std::to_string(treatwi[k]);
      } else if (data.numeric_cols.count(treat)) {
        treaterr = " " + treat + " = " + std::to_string(treatwn[k]);
      } else if (data.string_cols.count(treat)) {
        treaterr = " " + treat + " = " + treatwc[k];
      } else {
        throw std::invalid_argument("unsupported type for treat variable " + treat);
      }
      
      std::string str1 = "The milestone is larger than the largest observed time for";
      std::string errmsg = str1 + treaterr;
      if (!stratumerr.empty()) errmsg = errmsg + " " + stratumerr;
      
      if (noerr) {
        thread_utils::push_thread_warning(
          errmsg + "\nAdditional warning messages are suppressed.");
        noerr = false;
      }
      
      break;
    }
  }

  DataFrameCpp dfin;
  dfin.push_back(stratumn, "stratum");
  dfin.push_back(treatn, "treat");
  dfin.push_back(tstartn, "tstart");
  dfin.push_back(tstopn, "tstop");
  dfin.push_back(eventn, "event");
  dfin.push_back(weightn, "weight");
  
  std::vector<std::string> stratum_arg = {"stratum", "treat"};
  
  DataFrameCpp dfout = kmestcpp(dfin, stratum_arg, "tstart", "tstop", 
                                "event", "weight", "none", 0.95, false);
  
  std::vector<int> stratum2 = dfout.get<int>("stratum");
  std::vector<int> treat2 = dfout.get<int>("treat");
  std::vector<int> treatsize = dfout.get<int>("size");
  std::vector<double> time20 = dfout.get<double>("time");
  std::vector<double> survival2 = dfout.get<double>("surv");
  std::vector<double> stderr2 = dfout.get<double>("sesurv");
  
  int n2 = stratum2.size();
  
  // identify the locations of the unique values of stratum
  std::vector<int> idx2(1,0);
  for (int i=1; i<n2; ++i) {
    if (stratum2[i] != stratum2[i-1]) {
      idx2.push_back(i);
    }
  }
  
  nstrata = idx2.size();
  idx2.push_back(n2);
  
  std::vector<double> m(nstrata); // number of subjects in each stratum
  for (int i=0; i<nstrata; ++i) {
    int j1 = idx2[i], j2 = idx2[i+1] - 1;
    if (treat2[j1] != 1 || treat2[j2] != 2) {
      std::string stratumerr;
      if (!has_stratum) {
        stratumerr = "";
      } else {
        for (int j = 0; j < p_stratum; ++j) {
          std::string s = stratum[j];
          if (u_stratum.int_cols.count(s)) {
            auto v = u_stratum.get<int>(s);
            stratumerr += " " + s + " = " + std::to_string(v[i]);
          } else if (u_stratum.numeric_cols.count(s)) {
            auto v = u_stratum.get<double>(s);
            stratumerr += " " + s + " = " + std::to_string(v[i]);
          } else if (u_stratum.string_cols.count(s)) {
            auto v = u_stratum.get<std::string>(s);
            stratumerr += " " + s + " = " + v[i];
          } else {
            throw std::invalid_argument("unsupported type for stratum variable " + s);
          }
        }
      }
      
      int k = treat2[j1] != 1 ? 0 : 1;
      std::string treaterr;
      if (data.bool_cols.count(treat) || data.int_cols.count(treat)) {
        treaterr = " " + treat + " = " + std::to_string(treatwi[k]);
      } else if (data.numeric_cols.count(treat)) {
        treaterr = " " + treat + " = " + std::to_string(treatwn[k]);
      } else if (data.string_cols.count(treat)) {
        treaterr = " " + treat + " = " + treatwc[k];
      } else {
        throw std::invalid_argument("unsupported type for treat variable " + treat);
      }
      
      std::string str1 = "The data set does not contain";
      std::string errmsg = str1 + treaterr;
      if (!stratumerr.empty()) errmsg = errmsg + " " + stratumerr;
      
      if (noerr) {
        thread_utils::push_thread_warning(
          errmsg + "\nAdditional warning messages are suppressed.");
        noerr = false;
      }
      
      break;
    }
    
    m[i] = treatsize[j1] + treatsize[j2];
  }
  
  double M = std::accumulate(m.begin(), m.end(), 0.0);
  std::vector<double> p(nstrata);
  
  double surv1 = 0.0, surv2 = 0.0, vsurv1 = 0.0, vsurv2 = 0.0;
  for (int i=0; i<nstrata; ++i) {
    p[i] = m[i]/M; // fraction of subjects in the stratum
    std::vector<int> q = seqcpp(idx2[i], idx2[i+1]-1);
    std::vector<int> treatx = subset(treat2, q);
    std::vector<double> timex = subset(time20, q);
    std::vector<double> survivalx = subset(survival2, q);
    std::vector<double> stderrx = subset(stderr2, q);
    int nx = q.size();
    
    std::vector<double> surv(2), vsurv(2);
    for (int j=0; j<2; ++j) {
      
      std::vector<double> time0, survival0, stderr0;
      for (int k=0; k<nx; ++k) {
        if (treatx[k] == j+1) {
          time0.push_back(timex[k]);
          survival0.push_back(survivalx[k]);
          stderr0.push_back(stderrx[k]);
        }
      }
      int K = time0.size();
      
      // find the latest event time before milestone for each treat
      int k = 0;
      for (; k < K; ++k) {
        if (time0[k] > milestone) break;
      }
      
      if (k == 0) {
        surv[j] = 1;
        vsurv[j] = 0;
      } else {
        --k;
        surv[j] = survival0[k];
        vsurv[j] = stderr0[k] * stderr0[k];
      }
    }
    
    surv1 += p[i] * surv[0];
    surv2 += p[i] * surv[1];
    vsurv1 += p[i] * p[i] * vsurv[0];
    vsurv2 += p[i] * p[i] * vsurv[1];
  }
  
  double z = boost_qnorm((1.0 + conflev) / 2.0);
  
  double survDiff = surv1 - surv2;
  double sesurvDiff = std::sqrt(vsurv1 + vsurv2);
  double survDiffZ = (survDiff - survDiffH0)/sesurvDiff;
  double survDiffPValue = 2.0 * boost_pnorm(-std::fabs(survDiffZ));
  double lower = survDiff - z * sesurvDiff;
  double upper = survDiff + z * sesurvDiff;
  
  DataFrameCpp result;
  result.push_back(milestone, "milestone");
  result.push_back(survDiffH0, "survDiffH0");
  result.push_back(surv1, "surv1");
  result.push_back(surv2, "surv2");
  result.push_back(survDiff, "survDiff");
  result.push_back(vsurv1, "vsurv1");
  result.push_back(vsurv2, "vsurv2");
  result.push_back(sesurvDiff, "sesurvDiff");
  result.push_back(survDiffZ, "survDiffZ");
  result.push_back(survDiffPValue, "survDiffPValue");
  result.push_back(lower, "lower");
  result.push_back(upper, "upper");
  result.push_back(conflev, "conflev");
  
  return result;
}


// ----------------------------------------------------------------------------
// Parallel worker
// ----------------------------------------------------------------------------

struct kmdiffWorker : public RcppParallel::Worker {
  const std::vector<DataFrameCpp>* data_ptr;
  const std::vector<std::string>& stratum;
  const std::string& treat;
  const std::string& time;
  const std::string& time2;
  const std::string& event;
  const std::string& weight;
  const double milestone;
  const double survDiffH0;
  const double conflev;
  
  std::vector<DataFrameCpp>* results;
  
  kmdiffWorker(const std::vector<DataFrameCpp>* data_ptr_,
               const std::vector<std::string>& stratum_,
               const std::string& treat_,
               const std::string& time_,
               const std::string& time2_,
               const std::string& event_,
               const std::string& weight_,
               const double milestone_,
               const double survDiffH0_,
               const double conflev_,
               std::vector<DataFrameCpp>* results_)
    : data_ptr(data_ptr_), stratum(stratum_), treat(treat_), time(time_), 
      time2(time2_), event(event_), weight(weight_), milestone(milestone_), 
      survDiffH0(survDiffH0_), conflev(conflev_), results(results_) {}
  
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; ++i) {
      // Call the pure C++ function kmdiffcpp on data_ptr->at(i)
      DataFrameCpp out = kmdiffcpp(
        (*data_ptr)[i], stratum, treat, time, time2, event, weight, 
        milestone, survDiffH0, conflev);
      (*results)[i] = std::move(out);
    }
  }
};

// [[Rcpp::export]]
Rcpp::List kmdiffRcpp(SEXP data,
                      const std::vector<std::string>& stratum,
                      const std::string& treat = "treat",
                      const std::string& time = "time",
                      const std::string& time2 = "",
                      const std::string& event = "event",
                      const std::string& weight = "",
                      const double milestone = 0,
                      const double survDiffH0 = 0,
                      const double conflev = 0.95) {
  
  // Case A: single data.frame -> call kmdiffcpp on main thread
  if (Rf_inherits(data, "data.frame")) {
    Rcpp::DataFrame rdf(data);
    DataFrameCpp dfcpp = convertRDataFrameToCpp(rdf);
    
    // Call core C++ function directly on the DataFrameCpp
    DataFrameCpp cpp_result = kmdiffcpp(
      dfcpp, stratum, treat, time, time2, event, weight, 
      milestone, survDiffH0, conflev);
    
    thread_utils::drain_thread_warnings_to_R();
    return Rcpp::wrap(cpp_result);
  }
  
  // Case B: list of data.frames -> process in parallel
  if (TYPEOF(data) == VECSXP) {
    Rcpp::List lst(data);
    int m = lst.size();
    if (m == 0) return Rcpp::List(); // nothing to do
    
    // Convert each element to DataFrameCpp.
    std::vector<DataFrameCpp> data_vec;
    data_vec.reserve(m);
    for (int i = 0; i < m; ++i) {
      SEXP el = lst[i];
      if (!Rf_inherits(el, "data.frame")) {
        Rcpp::stop("When 'data' is a list, every element must be a data.frame.");
      }
      Rcpp::DataFrame rdf(el);
      
      DataFrameCpp dfcpp = convertRDataFrameToCpp(rdf);
      data_vec.push_back(std::move(dfcpp));
    }
    
    // Pre-allocate result vector of C++ objects (no R API used inside worker threads)
    std::vector<DataFrameCpp> results(m);
    
    // Build worker and run parallelFor across all indices [0, m)
    kmdiffWorker worker(
        &data_vec, stratum, treat, time, time2, event, weight, 
        milestone, survDiffH0, conflev, &results
    );
    
    // Execute parallelFor (this will schedule work across threads)
    RcppParallel::parallelFor(0, m, worker);
    
    // Drain thread-collected warnings (on main thread) into R's warning system
    thread_utils::drain_thread_warnings_to_R();
    
    // Convert C++ DataFrameCpp results back to R on the main thread
    Rcpp::List out(m);
    for (int i = 0; i < m; ++i) {
      out[i] = Rcpp::wrap(results[i]);
    }
    return out;
  }
  
  // Neither a data.frame nor a list: error
  Rcpp::stop("Input 'data' must be either a data.frame or a list of data.frames.");
  return R_NilValue; // unreachable
}


// Compute log-rank test statistic
DataFrameCpp lrtestcpp(const DataFrameCpp& data,
                       const std::vector<std::string>& stratum,
                       const std::string& treat,
                       const std::string& time,
                       const std::string& time2,
                       const std::string& event,
                       const std::string& weight,
                       const bool weight_readj,
                       const double rho1,
                       const double rho2) {
  int n = data.nrows();
  
  std::vector<int> stratumn(n);
  int p_stratum = stratum.size();
  if (!(p_stratum == 0 || (p_stratum == 1 && stratum[0] == ""))) {
    ListCpp out = bygroup(data, stratum);
    stratumn = out.get<std::vector<int>>("index");
  }
  
  // ---- treatment -> integer coding 1/2 ----
  if (!data.containElementNamed(treat)) {
    throw std::invalid_argument("data must contain the treat variable");
  }
  std::vector<int> treatn(n);
  if (data.bool_cols.count(treat) || data.int_cols.count(treat)) {
    std::vector<int> treatv(n);
    if (data.bool_cols.count(treat)) {
      const auto& vb = data.get<unsigned char>(treat);
      for (int i = 0; i < n; ++i) treatv[i] = vb[i] ? 1 : 0;
    } else treatv = data.get<int>(treat);
    auto treatwi = unique_sorted(treatv); // obtain unique treatment values
    if (treatwi.size() != 2)
      throw std::invalid_argument("treat must have two and only two distinct values");
    if (std::all_of(treatwi.begin(), treatwi.end(), [](int v) { return v == 0 || v == 1; })) {
      // special handling for 0/1 => map to 2,1 (so treated = 1, control = 2)
      for (int i = 0; i < n; ++i) treatn[i] = 2 - treatv[i];
    } else {
      treatn = matchcpp(treatv, treatwi, 1);
    }
  } else if (data.numeric_cols.count(treat)) {
    const auto& tv = data.get<double>(treat);
    auto treatwn = unique_sorted(tv);
    if (treatwn.size() != 2)
      throw std::invalid_argument("treat must have two and only two distinct values");
    if (std::all_of(treatwn.begin(), treatwn.end(), [](double v) { return v == 0.0 || v == 1.0; })) {
      for (int i = 0; i < n; ++i) treatn[i] = 2 - static_cast<int>(tv[i]);
    } else {
      treatn = matchcpp(tv, treatwn, 1);
    }
  } else if (data.string_cols.count(treat)) {
    const auto& tv = data.get<std::string>(treat);
    auto treatwc = unique_sorted(tv);
    if (treatwc.size() != 2) 
      throw std::invalid_argument("treat must have two and only two distinct values");
    treatn = matchcpp(tv, treatwc, 1);
  } else {
    throw std::invalid_argument("incorrect type for the treat variable in the input data");
  }
  
  // ---- time / time2 ----
  if (!data.containElementNamed(time)) {
    throw std::invalid_argument("data must contain the time variable");
  }
  std::vector<double> timen = data.get<double>(time);
  for (double v : timen) if (v < 0.0)
    throw std::invalid_argument("time must be nonnegative for each subject");
  
  bool has_time2 = !time2.empty() && data.containElementNamed(time2);
  std::vector<double> time2n(n);
  if (has_time2) {
    time2n = data.get<double>(time2);
    for (int i = 0; i < n; ++i) {
      if (time2n[i] <= timen[i]) {
        throw std::invalid_argument("time2 must be greater than time for each observation");
      }
    }
  }
  
  // unify right censored data with counting process data
  std::vector<double> tstartn(n), tstopn(n);
  if (!has_time2) {
    tstopn = std::move(timen);
    // tstartn remains zero-initialized
  } else {
    tstartn = std::move(timen);
    tstopn = std::move(time2n);
  }
  
  // ---- event -> numeric 0/1 ----
  if (!data.containElementNamed(event))
    throw std::invalid_argument("data must contain the event variable");
  std::vector<double> eventn(n);
  if (data.bool_cols.count(event)) {
    const auto& vb = data.get<unsigned char>(event);
    for (int i = 0; i < n; ++i) eventn[i] = vb[i] ? 1.0 : 0.0;
  } else if (data.int_cols.count(event)) {
    const auto& vi = data.get<int>(event);
    for (int i = 0; i < n; ++i) eventn[i] = static_cast<double>(vi[i]);
  } else if (data.numeric_cols.count(event)) {
    eventn = data.get<double>(event);
  } else {
    throw std::invalid_argument("event variable must be bool, integer or numeric");
  }
  for (double val : eventn) if (val != 0 && val != 1)
    throw std::invalid_argument("event must be 1 or 0 for each observation");
  
  std::vector<double> weightn(n, 1.0);
  if (!weight.empty() && data.containElementNamed(weight)) {
    weightn = data.get<double>(weight);
    for (double v : weightn) if (v <= 0.0)
      throw std::invalid_argument("weight must be greater than 0");
  }
  
  if (rho1 < 0) throw std::invalid_argument("rho1 must be non-negative");
  if (rho2 < 0) throw std::invalid_argument("rho2 must be non-negative");
  
  // sort by stopping time in descending order within each stratum
  std::vector<int> order(n);
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&](int i, int j) {
    if (stratumn[i] != stratumn[j]) return stratumn[i] < stratumn[j];
    if (tstopn[i] != tstopn[j]) return tstopn[i] > tstopn[j];
    return eventn[i] < eventn[j];
  });
  
  subset_in_place(stratumn, order);
  subset_in_place(treatn, order);
  subset_in_place(tstartn, order);
  subset_in_place(tstopn, order);
  subset_in_place(eventn, order);
  subset_in_place(weightn, order);
  
  // sort by starting time in descending order within each stratum
  std::vector<int> order1(n);
  std::iota(order1.begin(), order1.end(), 0);
  std::sort(order1.begin(), order1.end(), [&](int i, int j) {
    if (stratumn[i] != stratumn[j]) return stratumn[i] < stratumn[j];
    return tstartn[i] > tstartn[j];
  });
  
  double u = 0; // score
  double v = 0; // information (variance for score)
  
  double nrisk = 0, nrisk_1 = 0, nrisk_2 = 0;    // number at risk
  double nriskw = 0, nriskw_1 = 0, nriskw_2 = 0; // weighted # at risk
  double nriskw2_1 = 0, nriskw2_2 = 0;           // weight^2 # at risk
  double nevent = 0;                             // number of events
  double neventw = 0, neventw_1 = 0;             // weighted # events

  int istratum = stratumn[0]; // current stratum
  int i1 = 0;                 // index for removing out-of-risk subjects
  
  if (rho1 == 0 && rho2 == 0) {
    for (int i = 0; i < n; ) {
      // Reset when entering a new stratum
      if (stratumn[i] != istratum) {
        istratum = stratumn[i];
        i1 = i;
        
        nrisk = nrisk_1 = nrisk_2 = 0;
        nriskw = nriskw_1 = nriskw_2 = 0;
        nriskw2_1 = nriskw2_2 = 0;
      }
      
      const double dtime = tstopn[i];
      
      // Process all persons tied at this dtime
      for (; i < n && tstopn[i] == dtime && stratumn[i] == istratum; ++i) {
        ++nrisk;
        double w = weightn[i];
        nriskw += w;
        if (treatn[i] == 1) {
          ++nrisk_1; nriskw_1 += w; nriskw2_1 += w*w;
        } else {
          ++nrisk_2; nriskw_2 += w; nriskw2_2 += w*w;
        }
        
        if (eventn[i] == 1) {
          ++nevent; neventw += w;
          if (treatn[i] == 1) neventw_1 += w;
        }
      }
      
      // Remove subjects leaving risk set
      for (; i1<n; ++i1) {
        const int p1 = order1[i1];
        if (tstartn[p1] < dtime || stratumn[p1] != istratum) break;
        
        --nrisk;
        double w = weightn[p1]; 
        nriskw -= w;
        if (treatn[p1] == 1) {
          --nrisk_1; nriskw_1 -= w; nriskw2_1 -= w*w;
        } else {
          --nrisk_2; nriskw_2 -= w; nriskw2_2 -= w*w;
        }
      }
      
      // Add contributions for deaths at this time
      if (nevent > 0) {
        if (nrisk_1 > 0 && nrisk_2 > 0) {
          if (!weight_readj) {
            u += neventw_1 - neventw * (nriskw_1 / nriskw);
            double v1 = sq(nriskw_1 / nriskw) * nriskw2_2;
            double v2 = sq(nriskw_2 / nriskw) * nriskw2_1;
            v += nevent * (nrisk - nevent)/(nrisk * (nrisk - 1)) * (v1 + v2);
          } else {
            double neventw_2 = neventw - neventw_1;
            double neventwa_1 = neventw_1 * nrisk_1 / nriskw_1;
            double neventwa_2 = neventw_2 * nrisk_2 / nriskw_2;
            double neventwa = neventwa_1 + neventwa_2;
            u += neventwa_1 - neventwa * (nrisk_1 / nrisk);
            double v1 = sq(nrisk_1 / nrisk) * nriskw2_2 * sq(nrisk_2 / nriskw_2);
            double v2 = sq(nrisk_2 / nrisk) * nriskw2_1 * sq(nrisk_1 / nriskw_1);
            v += nevent * (nrisk - nevent)/(nrisk * (nrisk - 1)) * (v1 + v2);
          }
        }
        
        // Reset after processing deaths
        nevent = neventw = neventw_1 = 0;
      }
    }
  } else { // need to collect summaries for each event time and then iterate backward
    std::vector<int> stratum0;
    std::vector<double> time0;
    std::vector<double> nrisk0, nrisk_10, nrisk_20;
    std::vector<double> nriskw0, nriskw_10, nriskw_20;
    std::vector<double> nriskw2_10, nriskw2_20;
    std::vector<double> nevent0, neventw0, neventw_10;
    int d = n/4; // upper bound on number of unique event times
    stratum0.reserve(d); time0.reserve(d);
    nrisk0.reserve(d); nrisk_10.reserve(d); nrisk_20.reserve(d);
    nriskw0.reserve(d); nriskw_10.reserve(d); nriskw_20.reserve(d);
    nriskw2_10.reserve(d); nriskw2_20.reserve(d);
    nevent0.reserve(d); neventw0.reserve(d); neventw_10.reserve(d);
    
    for (int i = 0; i < n; ) {
      // Reset when entering a new stratum
      if (stratumn[i] != istratum) {
        istratum = stratumn[i];
        i1 = i;
        
        nrisk = nrisk_1 = nrisk_2 = 0;
        nriskw = nriskw_1 = nriskw_2 = 0;
        nriskw2_1 = nriskw2_2 = 0;
      }
      
      const double dtime = tstopn[i];
      
      // Process all persons tied at this dtime
      for (; i < n && tstopn[i] == dtime && stratumn[i] == istratum; ++i) {
        ++nrisk;
        double w = weightn[i]; 
        nriskw += w;
        if (treatn[i] == 1) {
          ++nrisk_1; nriskw_1 += w; nriskw2_1 += w*w;
        } else {
          ++nrisk_2; nriskw_2 += w; nriskw2_2 += w*w;
        }
        
        if (eventn[i] == 1) {
          ++nevent; neventw += w;
          if (treatn[i] == 1) neventw_1 += w;
        }
      }
      
      // Remove subjects leaving risk set
      for (; i1<n; ++i1) {
        const int p1 = order1[i1];
        if (tstartn[p1] < dtime || stratumn[p1] != istratum) break;
        
        --nrisk;
        double w = weightn[p1]; 
        nriskw -= w;
        if (treatn[p1] == 1) {
          --nrisk_1; nriskw_1 -= w; nriskw2_1 -= w*w;
        } else {
          --nrisk_2; nriskw_2 -= w; nriskw2_2 -= w*w;
        }
      }
      
      // add to the temporary storage
      if (nevent > 0) {
        stratum0.push_back(istratum);
        time0.push_back(dtime);
        nrisk0.push_back(nrisk);
        nrisk_10.push_back(nrisk_1);
        nrisk_20.push_back(nrisk_2);
        nriskw0.push_back(nriskw);
        nriskw_10.push_back(nriskw_1);
        nriskw_20.push_back(nriskw_2);
        nriskw2_10.push_back(nriskw2_1);
        nriskw2_20.push_back(nriskw2_2);
        nevent0.push_back(nevent);
        neventw0.push_back(neventw);
        neventw_10.push_back(neventw_1);

        // Reset after processing deaths
        nevent = neventw = neventw_1 = 0;
      }
    }
    
    int m = time0.size();
    if (m > 0) {
      double surv = 1.0;
      istratum = stratum0[m - 1];
      for (int i = m-1; i >= 0; --i) {
        nrisk = nrisk0[i];
        nrisk_1 = nrisk_10[i];
        nrisk_2 = nrisk_20[i];
        nriskw = nriskw0[i];
        nriskw_1 = nriskw_10[i];
        nriskw_2 = nriskw_20[i];
        nriskw2_1 = nriskw2_10[i];
        nriskw2_2 = nriskw2_20[i];
        nevent = nevent0[i];
        neventw = neventw0[i];
        neventw_1 = neventw_10[i];
        
        // Reset when entering a new stratum
        if (stratum0[i] != istratum) {
          istratum = stratum0[i];
          surv = 1;
        }
        
        // add to the main terms with S(t-)^rho1*(1 - S(t-))^rho2 weights
        double w = std::pow(surv, rho1) * std::pow(1 - surv, rho2);
        
        double neventwa = 0, neventwa_1 = 0, neventwa_2 = 0;
        if (weight_readj) {
          double neventw_2 = neventw - neventw_1;
          neventwa_1 = neventw_1 * nrisk_1 / nriskw_1;
          neventwa_2 = neventw_2 * nrisk_2 / nriskw_2;
          neventwa = neventwa_1 + neventwa_2;
        }
        
        if (nrisk_1 > 0 && nrisk_2 > 0) {
          if (!weight_readj) {
            u += w * (neventw_1 - neventw * (nriskw_1 / nriskw));
            double v1 = sq(nriskw_1 / nriskw) * nriskw2_2;
            double v2 = sq(nriskw_2 / nriskw) * nriskw2_1;
            v += w * w * nevent * (nrisk - nevent) / (nrisk * (nrisk - 1)) * (v1 + v2);
          } else {
            u += w * (neventwa_1 - neventwa * (nrisk_1 / nrisk));
            double v1 = sq(nrisk_1 / nrisk) * nriskw2_2 * sq(nrisk_2 / nriskw_2);
            double v2 = sq(nrisk_2 / nrisk) * nriskw2_1 * sq(nrisk_1 / nriskw_1);
            v += w * w * nevent * (nrisk - nevent) / (nrisk * (nrisk - 1)) * (v1 + v2);
          }
        }
        
        // update survival probability
        if (!weight_readj) surv *= (1.0 - neventw / nriskw);
        else surv *= (1.0 - neventwa / nrisk);
      }
    }
  }
  
  double z, p;
  if (v <= 0.0) {
    z = NaN; p = NaN;
  } else {
    z = u / std::sqrt(v);
    p = 2.0 * boost_pnorm(-std::fabs(z));
  }
  
  DataFrameCpp result;
  result.push_back(u, "uscore");
  result.push_back(v, "vscore");
  result.push_back(z, "logRankZ");
  result.push_back(p, "logRankPValue");
  result.push_back(weight_readj, "weight_readj");
  result.push_back(rho1, "rho1");
  result.push_back(rho2, "rho2");
  
  return result;
}


// ----------------------------------------------------------------------------
// Parallel worker
// ----------------------------------------------------------------------------

struct lrtestWorker : public RcppParallel::Worker {
  const std::vector<DataFrameCpp>* data_ptr;
  const std::vector<std::string>& stratum;
  const std::string& treat;
  const std::string& time;
  const std::string& time2;
  const std::string& event;
  const std::string& weight;
  const bool weight_readj;
  const double rho1;
  const double rho2;
  
  std::vector<DataFrameCpp>* results;
  
  lrtestWorker(const std::vector<DataFrameCpp>* data_ptr_,
               const std::vector<std::string>& stratum_,
               const std::string& treat_,
               const std::string& time_,
               const std::string& time2_,
               const std::string& event_,
               const std::string& weight_,
               const bool weight_readj_,
               const double rho1_,
               const double rho2_,
               std::vector<DataFrameCpp>* results_)
    : data_ptr(data_ptr_), stratum(stratum_), treat(treat_), time(time_), 
      time2(time2_), event(event_), weight(weight_), weight_readj(weight_readj_), 
      rho1(rho1_), rho2(rho2_), results(results_) {}
  
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; ++i) {
      // Call the pure C++ function lrtestcpp on data_ptr->at(i)
      DataFrameCpp out = lrtestcpp(
        (*data_ptr)[i], stratum, treat, time, time2, event, weight, 
        weight_readj, rho1, rho2);
      (*results)[i] = std::move(out);
    }
  }
};

// [[Rcpp::export]]
Rcpp::List lrtestRcpp(SEXP data,
                      const std::vector<std::string>& stratum,
                      const std::string& treat = "treat",
                      const std::string& time = "time",
                      const std::string& time2 = "",
                      const std::string& event = "event",
                      const std::string& weight = "",
                      const bool weight_readj = false,
                      const double rho1 = 0,
                      const double rho2 = 0) {
  
  // Case A: single data.frame -> call lrtestcpp on main thread
  if (Rf_inherits(data, "data.frame")) {
    Rcpp::DataFrame rdf(data);
    DataFrameCpp dfcpp = convertRDataFrameToCpp(rdf);
    
    // Call core C++ function directly on the DataFrameCpp
    DataFrameCpp cpp_result = lrtestcpp(
      dfcpp, stratum, treat, time, time2, event, weight, 
      weight_readj, rho1, rho2);
    
    thread_utils::drain_thread_warnings_to_R();
    return Rcpp::wrap(cpp_result);
  }
  
  // Case B: list of data.frames -> process in parallel
  if (TYPEOF(data) == VECSXP) {
    Rcpp::List lst(data);
    int m = lst.size();
    if (m == 0) return Rcpp::List(); // nothing to do
    
    // Convert each element to DataFrameCpp.
    std::vector<DataFrameCpp> data_vec;
    data_vec.reserve(m);
    for (int i = 0; i < m; ++i) {
      SEXP el = lst[i];
      if (!Rf_inherits(el, "data.frame")) {
        Rcpp::stop("When 'data' is a list, every element must be a data.frame.");
      }
      Rcpp::DataFrame rdf(el);
      
      DataFrameCpp dfcpp = convertRDataFrameToCpp(rdf);
      data_vec.push_back(std::move(dfcpp));
    }
    
    // Pre-allocate result vector of C++ objects (no R API used inside worker threads)
    std::vector<DataFrameCpp> results(m);
    
    // Build worker and run parallelFor across all indices [0, m)
    lrtestWorker worker(
        &data_vec, stratum, treat, time, time2, event, weight, 
        weight_readj, rho1, rho2, &results
    );
    
    // Execute parallelFor (this will schedule work across threads)
    RcppParallel::parallelFor(0, m, worker);
    
    // Drain thread-collected warnings (on main thread) into R's warning system
    thread_utils::drain_thread_warnings_to_R();
    
    // Convert C++ DataFrameCpp results back to R on the main thread
    Rcpp::List out(m);
    for (int i = 0; i < m; ++i) {
      out[i] = Rcpp::wrap(results[i]);
    }
    return out;
  }
  
  // Neither a data.frame nor a list: error
  Rcpp::stop("Input 'data' must be either a data.frame or a list of data.frames.");
  return R_NilValue; // unreachable
}
