#include <Rcpp.h>

#include "utilities.h"
#include "dataframe_list.h"

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/exponential_distribution.hpp>
#include <boost/random/bernoulli_distribution.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <random>
#include <stdexcept>
#include <vector>

// [[Rcpp::export]]
Rcpp::List tssimcpp(const bool tdxo, 
                    const bool coxo, 
                    const int allocation1,
                    const int allocation2,
                    const double p_X_1, 
                    const double p_X_0, 
                    const double rate_T, 
                    const double beta1, 
                    const double beta2, 
                    const double gamma0, 
                    const double gamma1, 
                    const double gamma2, 
                    const double gamma3, 
                    const double gamma4,
                    const double zeta0, 
                    const double zeta1, 
                    const double zeta2, 
                    const double zeta3, 
                    const double alpha0, 
                    const double alpha1, 
                    const double alpha2, 
                    const double theta1_1, 
                    const double theta1_0, 
                    const double theta2,
                    const double rate_C,
                    const std::vector<double>& accrualTime,
                    const std::vector<double>& accrualIntensity,
                    const double followupTime,
                    const bool fixedFollowup,
                    const double plannedTime,
                    const double days,
                    const int n, 
                    const int NSim, 
                    const int seed) {
  
  if (allocation1 < 1) 
    throw std::invalid_argument("allocation1 must be a positive integer");
  if (allocation2 < 1) 
    throw std::invalid_argument("allocation2 must be a positive integer");
  
  if (p_X_1 <= 0 || p_X_1 >= 1) 
    throw std::invalid_argument("p_X_1 must lie between 0 and 1");
  if (p_X_0 <= 0 || p_X_0 >= 1) 
    throw std::invalid_argument("p_X_0 must lie between 0 and 1");
  
  if (rate_T <= 0) throw std::invalid_argument("rate_T must be positive");
  if (rate_C < 0) throw std::invalid_argument("rate_C must be nonnegative");
  
  if (accrualTime[0] != 0) 
    throw std::invalid_argument("accrualTime must start with 0");
  if (accrualTime.size() > 1) {
    for (std::size_t i = 1; i < accrualTime.size(); ++i) {
      double prev = accrualTime[i - 1];
      double cur  = accrualTime[i];
      // If either element is NaN, skip the comparison
      if (std::isnan(prev) || std::isnan(cur)) continue;
      if (cur <= prev) {
        throw std::invalid_argument("accrualTime should be increasing");
      }
    }
  }
  
  // Check accrualIntensity has no missing values (NaN)
  if (std::any_of(accrualIntensity.begin(), accrualIntensity.end(),
                  [](double x){ return std::isnan(x); })) {
    throw std::invalid_argument("accrualIntensity must be provided");
  }
  
  // length check
  if (accrualTime.size() != accrualIntensity.size()) {
    throw std::invalid_argument(
        "accrualTime must have the same length as accrualIntensity");
  }
  
  // accrualIntensity: must be provided (no NaN) and non-negative
  for (std::size_t i = 0; i < accrualIntensity.size(); ++i) {
    double v = accrualIntensity[i];
    if (std::isnan(v)) {
      throw std::invalid_argument("accrualIntensity must be provided");
    }
    if (v < 0.0) {
      throw std::invalid_argument("accrualIntensity must be non-negative");
    }
  }
  
  // fixed follow-up checks
  if (fixedFollowup) {
    if (std::isnan(followupTime)) {
      throw std::invalid_argument("followupTime must be provided for fixed follow-up");
    }
    if (followupTime <= 0.0) {
      throw std::invalid_argument("followupTime must be positive for fixed follow-up");
    }
  }
  
  // plannedTime checks
  if (std::isnan(plannedTime)) 
    throw std::invalid_argument("plannedTime must be provided");
  if (plannedTime <= 0.0) throw std::invalid_argument("plannedTime must be positive");

  // days, n, NSim checks
  if (days <= 0.0) throw std::invalid_argument("days must be positive");
  if (n <= 0) throw std::invalid_argument("n must be positive");
  if (NSim <= 0) throw std::invalid_argument("NSim must be positive");
  
  // random number generator
  std::mt19937_64 rng;
  if (seed >= 0) {
    rng.seed(static_cast<std::uint64_t>(seed));
  } else {
    std::random_device rd;
    rng.seed(rd());
  }
  
  // distributions reused
  boost::random::uniform_real_distribution<double> unif01(0.0, 1.0);
  
  int maxFollowup = static_cast<int>(std::ceil(plannedTime / days));
  int K = n * maxFollowup;
  std::vector<DataFrameCpp> sims(NSim);
  
  for (int iter = 0; iter < NSim; ++iter) {
    std::vector<int> idx(K), trtrandx(K), tpointx(K);
    std::vector<int> bprogx(K), Lx(K), Llagx(K), Zx(K), Zlagx(K);
    std::vector<int> Ax(K), Alagx(K), Alag2x(K), eventx(K);
    std::vector<int> diedx(K), progressedx(K), xox(K);
    std::vector<double> tstartx(K), tstopx(K), timeOSx(K), timePDx(K);
    std::vector<double> xotimex(K, NA_REAL), censor_timex(K);
    std::vector<double> arrivalTimex(K);
    
    double b1 = allocation1, b2 = allocation2;
    double enrollt = 0;
    int k = 0;
    for (int i = 1; i <= n; ++i) { // subject index (1..n)
      int id = i;
      
      // generate accrual time
      double u = unif01(rng);
      enrollt = qtpwexpcpp1(u, accrualTime, accrualIntensity, enrollt, 1, 0);
      double arrivalTime = std::ceil(enrollt);
      
      // stratified block randomization
      // stratified block randomization
      u = unif01(rng);
      int trtrand;
      if (u <= b1 / (b1 + b2)) { trtrand = 1; b1 -= 1.0; }
      else { trtrand = 0; b2 -= 1.0; }
      if (b1 == 0.0 && b2 == 0.0) { b1 = allocation1; b2 = allocation2; }
      
      // baseline prognostic bprog (Bernoulli)
      double p_bprog = (trtrand == 1) ? p_X_1 : p_X_0;
      std::bernoulli_distribution bern_bprog(p_bprog);
      int bprog = bern_bprog(rng) ? 1 : 0;
      
      // event time T from exponential with rate = rate_T * exp(...)
      double rate_this = rate_T * std::exp(beta1 * trtrand + beta2 * bprog);
      std::exponential_distribution<double> expT(rate_this);
      double T = std::round(expT(rng));
      if (T == 0.0) T = 1.0;
      
      // follow-up and cycles as before...
      double fu = std::max(plannedTime - arrivalTime, 0.0);
      if (fixedFollowup) fu = std::min(followupTime, fu);
      int followup = static_cast<int>(std::floor(fu / days));
      
      int tpoint = followup; 
      int L = 0, Llag = 0, Z = 0, Zlag = 0, A = 0, Alag = 0, Alag2 = 0;
      for (int j = 1; j <= followup; ++j) { // j = cycle index
        tpoint = j; 
        double tstart = days * (j - 1);
        if (j == 1) {
          Llag = 0; Zlag = 0; Alag2 = 0; Alag = 0;
        } else {
          Llag = L; Zlag = Z; Alag2 = Alag; Alag = A;
        }
        
        // generate time-dependent covariate
        double probL = boost_plogis(gamma0 + gamma1 * Alag + gamma2 * Llag +
                                    gamma3 * bprog + gamma4 * trtrand);
        std::bernoulli_distribution bernL(probL);
        L = bernL(rng) ? 1 : 0;
        
        double tstop;
        int event;
        if (T <= days * j) { // died in cycle j, complete data for the subject
          tstop = T; event = 1;
          Z = INT_MIN; A = INT_MIN;
        } else { // alive at the end of cycle j, continue to the next cycle
          tstop = days * j; event = 0;
          
          // generate disease progression status
          if (Zlag == 0) {
            double probZ = boost_plogis(zeta0 + zeta1 * L + zeta2 * bprog + 
                                        zeta3 * trtrand);
            std::bernoulli_distribution bernZ(probZ);
            Z = bernZ(rng) ? 1 : 0;
          } else {
            Z = 1;
          }
          
          // generate treatment switching status
          if (Alag == 0) {
            bool condition_for_switch = ((tdxo == 0 && Z == 1 && Zlag == 0) ||
                                         (tdxo == 1 && Z == 1)) &&
                                         ((coxo == 1 && trtrand == 0) || (coxo == 0));
            if (condition_for_switch) {
              double probA = boost_plogis(alpha0 + alpha1 * L + alpha2 * bprog);
              std::bernoulli_distribution bernA(probA);
              A = bernA(rng) ? 1 : 0;
            } else {
              A = 0;
            }
          } else {
            A = 1;
          }
          
          // update survival time
          double theta1 = theta1_1 * trtrand + theta1_0 * (1 - trtrand);
          T = days * j + (T - days * j) * std::exp(-theta1 * (A - Alag) 
                                                     - theta2 * (L - Llag));
          T = std::round(T);
        }
        
        // add the data from the current cycle
        idx[k] = id;
        arrivalTimex[k] = arrivalTime;
        trtrandx[k] = trtrand;
        bprogx[k] = bprog;
        tpointx[k] = tpoint;
        tstartx[k] = tstart;
        tstopx[k] = tstop;
        Lx[k] = L;
        Llagx[k] = Llag;
        Zx[k] = Z;
        Zlagx[k] = Zlag;
        Ax[k] = A;
        Alagx[k] = Alag;
        Alag2x[k] = Alag2;
        eventx[k] = event;
        
        ++k;
        
        if (event == 1) break;
      }
      
      // generate random censoring due to dropout
      std::exponential_distribution<double> expC(rate_C);
      double C = std::round(expC(rng));
      if (C == 0.0) C = 1.0;
      double time = std::min({T, C, fu});
      
      int J;  // J is the number of treatment cycles
      if (time <= days*followup) {
        J = static_cast<int>(std::ceil(time / days));
        k = k - tpoint + J; // discard treatment cycles after censoring
        tstopx[k-1] = time; // update the ending time and event indicator
        eventx[k-1] = (T == time) ? 1 : 0;
        Zx[k-1] = INT_MIN;
        Ax[k-1] = INT_MIN;
      } else { // add one more record 
        J = followup + 1;
        idx[k] = id;
        arrivalTimex[k] = arrivalTime;
        trtrandx[k] = trtrand;
        bprogx[k] = bprog;
        tpointx[k] = followup + 1;
        tstartx[k] = days*followup;
        tstopx[k] = time;
        double probL_final = boost_plogis(gamma0 + gamma1 * Alag + gamma2 * Llag 
                                            + gamma3 * bprog + gamma4 * trtrand);
        std::bernoulli_distribution bernL_final(probL_final);
        Lx[k] = bernL_final(rng) ? 1 : 0;
        Llagx[k] = L;
        Zx[k] = INT_MIN;
        Zlagx[k] = Z;
        Ax[k] = INT_MIN;
        Alagx[k] = A;
        Alag2x[k] = Alag;
        eventx[k] = (T == time) ? 1 : 0;
        
        ++k;
      }
      
      // create subject-level survival time and death indicator
      for (int j = k - J; j < k; j++) {
        timeOSx[j] = tstopx[k-1];
        diedx[j] = eventx[k-1];
        censor_timex[j] = fu;
      }
      
      // progression and time to progression (if applicable)
      int pd = 0;
      double pd_time = NA_REAL;
      for (int j = k - J; j < k; j++) {
        if (Zx[j] == 1) { pd = 1; pd_time = tstopx[j]; break; }
      }

      // switching and time to switching (if applicable)
      int xo = 0;
      double xo_time = NaN;
      for (int j = k - J; j < k; j++) {
        if (Ax[j] == 1) { xo = 1; xo_time = tstopx[j]; break; }
      }
      
      for (int j = k - J; j < k; j++) {
        progressedx[j] = pd;
        timePDx[j] = pd_time;
        xox[j] = xo;
        xotimex[j] = xo_time;
      }
      
      // shift disease progression and alternative therapy status downward
      for (int j = k - J; j < k; j++) {
        Zx[j] = Zlagx[j];
        Ax[j] = Alagx[j];
        Alagx[j] = Alag2x[j];
      }
    }
    
    std::vector<int> sub = seqcpp(0, k - 1);
    subset_in_place(idx, sub);
    subset_in_place(arrivalTimex, sub);
    subset_in_place( trtrandx, sub);
    subset_in_place(bprogx, sub);
    subset_in_place(tpointx, sub);
    subset_in_place(tstartx, sub);
    subset_in_place(tstopx, sub);
    subset_in_place(Lx, sub);
    subset_in_place(Llagx, sub);
    subset_in_place(Zx, sub);
    subset_in_place(Ax, sub);
    subset_in_place(Alagx, sub);
    subset_in_place(eventx, sub);
    subset_in_place(diedx, sub);
    subset_in_place(progressedx, sub);
    subset_in_place(timeOSx, sub);
    subset_in_place(timePDx, sub);
    subset_in_place(xox, sub);
    subset_in_place(xotimex, sub);
    subset_in_place(censor_timex, sub);
    
    DataFrameCpp data;
    data.push_back(idx, "id");
    data.push_back(arrivalTimex, "arrival_time");
    data.push_back(trtrandx, "trtrand");
    data.push_back(bprogx, "bprog");
    data.push_back(tpointx, "tpoint");
    data.push_back(tstartx, "tstart");
    data.push_back(tstopx, "tstop");
    data.push_back(Lx, "L");
    data.push_back(Llagx, "Llag");
    data.push_back(Zx, "Z");
    data.push_back(Ax, "A");
    data.push_back(Alagx, "Alag");
    data.push_back(eventx, "event");
    data.push_back(timeOSx, "timeOS");
    data.push_back(diedx, "died");
    data.push_back(progressedx, "progressed");
    data.push_back(timePDx, "timePD");
    data.push_back(xox, "xo");
    data.push_back(xotimex, "xotime");
    data.push_back(censor_timex, "censor_time");
    
    sims[iter] = std::move(data);
  }
  
  Rcpp::List out(sims.size());
  for (std::size_t i = 0; i < sims.size(); ++i) {
    out[i] = convertDataFrameCppToR(sims[i]);
  }
  return out;
}