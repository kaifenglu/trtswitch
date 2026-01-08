#include <Rcpp.h>

#include <boost/random.hpp>

#include "survival_analysis.h"
#include "utilities.h"
#include "dataframe_list.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <random>
#include <stdexcept>
#include <vector>


//' @title Simulate Survival Data for Two-Stage Estimation with  
//' g-estimation
//' @description Obtains the simulated data for baseline prognosis, 
//' disease progression, treatment switching, death, and 
//' time-dependent covariates.
//'
//' @param n The total sample size for two treatment arms combined.
//' @param allocation1 The number of subjects in the active treatment group 
//'   in a randomization block.
//' @param allocation2 The number of subjects in the control group in
//'   a randomization block.
//' @param pbprog The probability of having poor prognosis at baseline.
//' @param trtlghr The treatment effect in terms of log hazard ratio.
//' @param bprogsl The poor prognosis effect in terms of log hazard ratio.
//' @param shape1 The shape parameter for the Weibull event distribution 
//'   for the first component.
//' @param scale1 The scale parameter for the Weibull event distribution 
//'   for the first component.
//' @param shape2 The shape parameter for the Weibull event distribution 
//'   for the second component.
//' @param scale2 The scale parameter for the Weibull event distribution 
//'   for the second component.
//' @param pmix The mixing probability of the first component Weibull 
//'   distribution.
//' @param admin The administrative censoring time.
//' @param pcatnotrtbprog The probability of developing metastatic disease
//'   on control treatment with poor baseline prognosis.
//' @param pcattrtbprog The probability of developing metastatic disease
//'   on active treatment with poor baseline prognosis.
//' @param pcatnotrt The probability of developing metastatic disease
//'   on control treatment with good baseline prognosis.
//' @param pcattrt The probability of developing metastatic disease
//'   on active treatment with good baseline prognosis.
//' @param catmult The impact of metastatic disease on shortening remaining 
//'   survival time.
//' @param tdxo Whether treatment crossover depends on time-dependent 
//'   covariates between disease progression and treatment switching.
//' @param ppoor The probability of switching for poor baseline prognosis
//'   with no metastatic disease.
//' @param pgood The probability of switching for good baseline prognosis
//'   with no metastatic disease.
//' @param ppoormet The probability of switching for poor baseline prognosis
//'   after developing metastatic disease.
//' @param pgoodmet The probability of switching for good baseline prognosis
//'   after developing metastatic disease.
//' @param xomult The direct effect of crossover on extending remaining 
//'   survival time.
//' @param milestone The milestone to calculate restricted mean survival 
//'   time.
//' @param seed The seed to reproduce the simulation results.
//'
//' @return A list with two data frames.
//' 
//' * \code{sumdata}: A summary data frame with the following variables:
//'
//'     - \code{simtrueconstmean}: The true control group restricted mean 
//'       survival time (RMST).
//'     
//'     - \code{simtrueconstlb}: The lower bound for control group RMST.
//'     
//'     - \code{simtrueconstub}: The upper bound for control group RMST.
//'     
//'     - \code{simtrueconstse}: The standard error for control group RMST.
//'     
//'     - \code{simtrueexpstmean}: The true experimental group restricted 
//'       mean survival time (RMST).
//'     
//'     - \code{simtrueexpstlb}: The lower bound for experimental group RMST.
//'     
//'     - \code{simtrueexpstub}: The upper bound for experimental group RMST.
//'     
//'     - \code{simtrueexpstse}: The standard error for experimental group 
//'       RMST.
//'     
//'     - \code{simtrue_coxwbprog_hr}: The treatment hazard ratio from the 
//'       Cox model adjusting for baseline prognosis.
//'     
//'     - \code{simtrue_cox_hr}: The treatment hazard ratio from the Cox 
//'       model without adjusting for baseline prognosis.
//'
//'     - \code{simtrue_aftwbprog_af}: The average acceleration factor from 
//'       the Weibull AFT model adjusting for baseline prognosis.
//'
//'     - \code{simtrue_aft_af}: The average acceleration factor from 
//'       the Weibull AFT model without adjusting for baseline prognosis.
//'
//' * \code{paneldata}: A counting process style subject-level data frame 
//'   with the following variables:
//'
//'     - \code{id}: The subject ID.
//'
//'     - \code{trtrand}: The randomized treatment arm.
//'
//'     - \code{bprog}: Whether the patient had poor baseline prognosis. 
//'     
//'     - \code{tstart}: The left end of time interval.
//'     
//'     - \code{tstop}: The right end of time interval.
//'     
//'     - \code{event}: Whether the patient died at the end of the interval. 
//'     
//'     - \code{timeOS}: The observed survival time.
//'     
//'     - \code{died}: Whether the patient died during the study. 
//'     
//'     - \code{progressed}: Whether the patient had disease progression. 
//'     
//'     - \code{timePFSobs}: The observed time of disease progression at 
//'       regular scheduled visits.
//'       
//'     - \code{progtdc}: The time-dependent covariate for progression.
//'       
//'     - \code{catevent}: Whether the patient developed metastatic disease.
//'     
//'     - \code{cattime}: When the patient developed metastatic disease.
//'     
//'     - \code{cattdc}: The time-dependent covariate for cat event.
//'     
//'     - \code{xo}: Whether the patient switched treatment. 
//'     
//'     - \code{xotime}: When the patient switched treatment.
//'     
//'     - \code{xotdc}: The time-dependent covariate for treatment 
//'       switching.
//'       
//'     - \code{censor_time}: The administrative censoring time.
//'     
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @references
//' NR Latimer, IR White, K Tilling, and U Siebert.
//' Improved two-stage estimation to adjust for treatment switching in 
//' randomised trials: g-estimation to address time-dependent confounding.
//' Statistical Methods in Medical Research. 2020;29(10):2900-2918.
//'
//' @examples
//'
//' sim1 <- tsegestsim(
//'   n = 500, allocation1 = 2, allocation2 = 1, pbprog = 0.5, 
//'   trtlghr = -0.5, bprogsl = 0.3, shape1 = 1.8, 
//'   scale1 = 360, shape2 = 1.7, scale2 = 688, 
//'   pmix = 0.5, admin = 5000, pcatnotrtbprog = 0.5, 
//'   pcattrtbprog = 0.25, pcatnotrt = 0.2, pcattrt = 0.1, 
//'   catmult = 0.5, tdxo = 1, ppoor = 0.1, pgood = 0.04, 
//'   ppoormet = 0.4, pgoodmet = 0.2, xomult = 1.4188308, 
//'   milestone = 546, seed = 2000)
//'
//' @export
// [[Rcpp::export]]
Rcpp::List tsegestsim(const int n = 500,
                      const int allocation1 = 2,
                      const int allocation2 = 1,
                      const double pbprog = 0.5,
                      const double trtlghr = -0.5,
                      const double bprogsl = 0.3,
                      const double shape1 = 1.8,
                      const double scale1 = 360,
                      const double shape2 = 1.7,
                      const double scale2 = 688,
                      const double pmix = 0.5,
                      const double admin = 5000,
                      const double pcatnotrtbprog = 0.5,
                      const double pcattrtbprog = 0.25,
                      const double pcatnotrt = 0.2,
                      const double pcattrt = 0.1,
                      const double catmult = 0.5,
                      const double tdxo = 1,
                      const double ppoor = 0.1,
                      const double pgood = 0.04,
                      const double ppoormet = 0.4,
                      const double pgoodmet = 0.2,
                      const double xomult = 1.4188308,
                      const double milestone = 546,
                      const int seed = 0) {
  
  // random number generator
  boost::random::mt19937_64 rng(seed);
  
  // distributions reused
  boost::random::uniform_real_distribution<double> unif(0.0, 1.0);
  boost::random::gamma_distribution<double> ga(5, 1.0);  // shape = 5, scale = 1
  boost::random::gamma_distribution<double> gb(10, 1.0); // shape = 10, scale = 1
  
  // survival function of the Weibull mixture
  auto S = [shape1, scale1, shape2, scale2, pmix](double t) -> double {
    double a1 = pmix * std::exp(-std::pow(t / scale1, shape1));
    double a2 = (1 - pmix) * std::exp(-std::pow(t / scale2, shape2));
    return a1 + a2;
  };
  
  double simtrueconstmean, simtrueconstlb, simtrueconstub, simtrueconstse;
  double simtrueexpstmean, simtrueexpstlb, simtrueexpstub, simtrueexpstse;
  double simtrue_coxwbprog_hr, simtrue_cox_hr;
  double simtrue_aftwbprog_af, simtrue_aft_af;
  
  std::vector<int> id(n), trtrand(n), bprog(n), dead(n), progressed(n);
  std::vector<double> timeOS(n), timeOS5(n);
  std::vector<double> timePFS(n), timePFSobs(n, NaN);
  std::vector<double> probcat(n);
  std::vector<int> cat2(n, INT_MIN), cat3(n, INT_MIN);
  std::vector<int> cat4(n, INT_MIN), cat5(n, INT_MIN);
  std::vector<int> catevent(n, INT_MIN);
  std::vector<double> catOSloss(n, NaN), cattime(n, NaN);
  
  double b1 = allocation1, b2 = allocation2;
  for (int i = 0; i < n; ++i) {
    id[i] = i + 1; // 1-based ID
    
    // generate treatment indicators using stratified block randomization
    if (unif(rng) <= b1 / (b1 + b2)) { trtrand[i] = 1; b1--; } 
    else { trtrand[i] = 0; b2--; }
    
    // start a new block after depleting the current block
    if (b1 == 0 && b2 == 0) { b1 = allocation1; b2 = allocation2; }
    
    // generate poor baseline prognosis indicators
    bprog[i] = (unif(rng) <= pbprog) ? 1 : 0;
    
    // generate survival times from Weibull mixture
    double eta = trtlghr * trtrand[i] + bprogsl * bprog[i];
    double v = std::pow(unif(rng), std::exp(-eta));
    timeOS[i] = squantilecpp(S, v, 1.0e-6);
    dead[i] = (timeOS[i] <= admin);
    timeOS[i] = std::min(timeOS[i], admin);
    timeOS[i] = std::round(timeOS[i]);
    if (timeOS[i] == 0) timeOS[i] = 1;
    
    // generate observed time to disease progression
    double g1 = ga(rng), g2 = gb(rng);
    double u = g1 / (g1 + g2);  // Beta(5, 10) distribution
    timePFS[i] = std::round(timeOS[i] * u);
    
    // scheduled visits are every 21 days
    int k = static_cast<int>(std::floor(timeOS[i] / 21));
    for (int j = 1; j <= k; ++j) {
      if (timePFS[i] < j*21 && timeOS[i] > j*21) {
        timePFSobs[i] = j*21;
        break;
      }
    }
    
    if (!std::isnan(timePFSobs[i])) {
      progressed[i] = 1;
    } else {
      timePFSobs[i] = timeOS[i];
      progressed[i] = 0;
    }
    
    // apply censoring to observed time to disease progression
    if (timePFSobs[i] > admin) {
      progressed[i] = 0;
      timePFSobs[i] = admin;
    }
    
    // generate catastrophic event that can happen after progression
    if (trtrand[i] == 0 && bprog[i] == 1) {
      probcat[i] = pcatnotrtbprog;
    } else if (trtrand[i] == 1 && bprog[i] == 1) {
      probcat[i] = pcattrtbprog;
    } else if (trtrand[i] == 0 && bprog[i] == 0) {
      probcat[i] = pcatnotrt;
    } else if (trtrand[i] == 1 && bprog[i] == 0) {
      probcat[i] = pcattrt;
    }
    
    if (progressed[i] == 1 && timeOS[i] > timePFSobs[i] + 21) {
      cat2[i] = (unif(rng) <= probcat[i]) ? 1 : 0;
    }
    
    if (cat2[i] == 0 &&
        progressed[i] == 1 && timeOS[i] > timePFSobs[i] + 42) {
      cat3[i] = (unif(rng) <= probcat[i]) ? 1 : 0;
    }
    
    if (cat2[i] == 0 && cat3[i] == 0 &&
        progressed[i] == 1 && timeOS[i] > timePFSobs[i] + 63) {
      cat4[i] = (unif(rng) <= probcat[i]) ? 1 : 0;
    }
    
    if (cat2[i] == 0 && cat3[i] == 0 && cat4[i] == 0 &&
        progressed[i] == 1 && timeOS[i] > timePFSobs[i] + 84) {
      cat5[i] = (unif(rng) <= probcat[i]) ? 1 : 0;
    }
    
    if (cat2[i] == 1) {
      cattime[i] = timePFSobs[i] + 21;
    } else if (cat3[i] == 1) {
      cattime[i] = timePFSobs[i] + 42;
    } else if (cat4[i] == 1) {
      cattime[i] = timePFSobs[i] + 63;
    } else if (cat5[i] == 1) {
      cattime[i] = timePFSobs[i] + 84;
    }
    
    // amend survival times
    timeOS5[i] = timeOS[i]; // original survival time
    if (cat2[i] == 1 || cat3[i] == 1 || cat4[i] == 1 || cat5[i] == 1) {
      catevent[i] = 1;
      catOSloss[i] = timeOS[i] - cattime[i];
      catOSloss[i] = std::round(catOSloss[i] * catmult);
      timeOS[i] = catOSloss[i] + cattime[i];
    }
  }
  
  // calculate HR and RMST with no switching
  std::vector<int> event(n);
  std::vector<double> time(n);
  for (int i = 0; i < n; ++i) {
    if (timeOS[i] > admin) {
      event[i] = 0;
      time[i] = admin;
    } else {
      event[i] = 1;
      time[i] = timeOS[i];
    }
  }
  
  DataFrameCpp a1;
  a1.push_back(time, "time");
  a1.push_back(event, "event");
  a1.push_back(trtrand, "trtrand");
  a1.push_back(bprog, "bprog");
  
  // RMST
  DataFrameCpp a2 = rmestcpp(a1, {"trtrand"}, "time", "event", milestone);
  std::vector<double> simtruestmean = a2.get<double>("rmst");
  std::vector<double> simtruestlb = a2.get<double>("lower");
  std::vector<double> simtruestub = a2.get<double>("upper");
  std::vector<double> simtruestse = a2.get<double>("stderr");
  simtrueconstmean = simtruestmean[0];
  simtrueexpstmean = simtruestmean[1];
  simtrueconstlb = simtruestlb[0];
  simtrueexpstlb = simtruestlb[1];
  simtrueconstub = simtruestub[0];
  simtrueexpstub = simtruestub[1];
  simtrueconstse = simtruestse[0];
  simtrueexpstse = simtruestse[1];
  
  // HR with bprog
  std::vector<std::string> covariates = {"trtrand", "bprog"};
  
  std::vector<double> init(1, NaN);
  ListCpp a3 = phregcpp(a1, {""}, "time", "", "event", covariates, 
                        "", "", "", "efron", init, 0, 0, 0, 0, 0);
  DataFrameCpp parest = a3.get<DataFrameCpp>("parest");
  std::vector<double> beta = parest.get<double>("beta");
  simtrue_coxwbprog_hr = std::exp(beta[0]);
  
  // HR without bprog
  a3 = phregcpp(a1, {""}, "time", "", "event", {"trtrand"}, 
                "", "", "", "efron", init, 0, 0, 0, 0, 0);
  parest = a3.get<DataFrameCpp>("parest");
  beta = parest.get<double>("beta");
  simtrue_cox_hr = std::exp(beta[0]);
  
  // acceleration factor from AFT with bprog
  ListCpp a4 = liferegcpp(a1, {""}, "time", "", "event", covariates,
                          "", "", "", "weibull", init, 0, 0);
  parest = a4.get<DataFrameCpp>("parest");
  beta = parest.get<double>("beta");
  simtrue_aftwbprog_af = std::exp(beta[1]);
  
  // acceleration factor from AFT without bprog
  ListCpp a5 = liferegcpp(a1, {""}, "time", "", "event", {"trtrand"},
                          "", "", "", "weibull", init, 0, 0);
  parest = a5.get<DataFrameCpp>("parest");
  beta = parest.get<double>("beta");
  simtrue_aft_af = std::exp(beta[1]);
  
  // apply switch and effect
  std::vector<int> xo1(n, INT_MIN), xo2(n, INT_MIN), xo3(n, INT_MIN);
  std::vector<int> xo4(n, INT_MIN), xo5(n, INT_MIN), xo6(n, INT_MIN);
  std::vector<int> xo(n, INT_MIN), xoprecat(n, INT_MIN);
  std::vector<double> xotime(n, NaN);
  std::vector<int> extra2v1(n), extra3v1(n), extra4v1(n), extra5v1(n);
  std::vector<int> extraobsv1(n);
  std::vector<double> xoOSgainobs(n, NaN), timeOS2(n);
  std::vector<int> extra2v2(n), extra3v2(n), extra4v2(n), extra5v2(n);
  std::vector<int> extraobsv2(n);
  for (int i = 0; i < n; ++i) {
    double p1, p2, p3, p4, p5, p6;
    
    // prob of switching depends on bprog for the first 2 visits after PD
    p1 = bprog[i] == 1 ? ppoor : pgood;
    p2 = p1;
    
    // for subsequent visits, prob of switching depends on bprog and 
    // metastatic disease status at the previous visit
    if (bprog[i] == 1 && cat2[i] == 1) {
      p3 = ppoormet;
    } else if (bprog[i] == 0 && cat2[i] == 1) {
      p3 = pgoodmet;
    } else if (bprog[i] == 1 && cat2[i] != 1) {
      p3 = ppoor;
    } else {
      p3 = pgood;
    }
    
    if (bprog[i] == 1 && (cat2[i] == 1 || cat3[i] == 1)) {
      p4 = ppoormet;
    } else if (bprog[i] == 0 && (cat2[i] == 1 || cat3[i] == 1)) {
      p4 = pgoodmet;
    } else if (bprog[i] == 1 && (cat2[i] != 1 && cat3[i] != 1)) {
      p4 = ppoor;
    } else {
      p4 = pgood;
    }
    
    if (bprog[i] == 1 && (cat2[i] == 1 || cat3[i] == 1 || cat4[i] == 1)) {
      p5 = ppoormet;
    } else if (bprog[i] == 0 && (cat2[i] == 1 || cat3[i] == 1 ||
      cat4[i] == 1)) {
      p5 = pgoodmet;
    } else if (bprog[i] == 1 && (cat2[i] != 1 && cat3[i] != 1 &&
      cat4[i] != 1)) {
      p5 = ppoor;
    } else {
      p5 = pgood;
    }
    
    if (bprog[i] == 1 && (cat2[i] == 1 || cat3[i] == 1 || cat4[i] == 1 ||
        cat5[i] == 1)) {
      p6 = ppoormet;
    } else if (bprog[i] == 0 && (cat2[i] == 1 || cat3[i] == 1 ||
      cat4[i] == 1 || cat5[i] == 1)) {
      p6 = pgoodmet;
    } else if (bprog[i] == 1 && (cat2[i] != 1 && cat3[i] != 1 &&
      cat4[i] != 1 && cat5[i] != 1)) {
      p6 = ppoor;
    } else {
      p6 = pgood;
    }
    
    // generate the indicator and time of treatment switching
    if (trtrand[i] == 0 &&
        timeOS[i] > timePFSobs[i] && progressed[i] == 1) {
      xo1[i] = (unif(rng) <= p1) ? 1 : 0;;
    }
    
    if (trtrand[i] == 0 && xo1[i] == 0 &&
        timeOS[i] > timePFSobs[i] + 21 && progressed[i] == 1 && 
        tdxo == 1) {
      xo2[i] = (unif(rng) <= p2) ? 1 : 0;
    }
    
    if (trtrand[i] == 0 && xo1[i] == 0 && xo2[i] == 0 &&
        timeOS[i] > timePFSobs[i] + 42 && progressed[i] == 1 && 
        tdxo == 1) {
      xo3[i] = (unif(rng) <= p3) ? 1 : 0;
    }
    
    if (trtrand[i] == 0 && xo1[i] == 0 && xo2[i] == 0 && xo3[i] == 0 &&
        timeOS[i] > timePFSobs[i] + 63 && progressed[i] == 1 && 
        tdxo == 1) {
      xo4[i] = (unif(rng) <= p4) ? 1 : 0;
    }
    
    if (trtrand[i] == 0 && xo1[i] == 0 && xo2[i] == 0 && xo3[i] == 0 &&
        xo4[i] == 0 &&
        timeOS[i] > timePFSobs[i] + 84 && progressed[i] == 1 && 
        tdxo == 1) {
      xo5[i] = (unif(rng) <= p5) ? 1 : 0;
    }
    
    if (trtrand[i] == 0 && xo1[i] == 0 && xo2[i] == 0 && xo3[i] == 0 &&
        xo4[i] == 0 && xo5[i] == 0 &&
        timeOS[i] > timePFSobs[i] + 105 && progressed[i] == 1 && 
        tdxo == 1) {
      xo6[i] = (unif(rng) <= p6) ? 1 : 0;
    }
    
    if (xo1[i] == 1) xotime[i] = timePFSobs[i];
    if (xo2[i] == 1) xotime[i] = timePFSobs[i] + 21;
    if (xo3[i] == 1) xotime[i] = timePFSobs[i] + 42;
    if (xo4[i] == 1) xotime[i] = timePFSobs[i] + 63;
    if (xo5[i] == 1) xotime[i] = timePFSobs[i] + 84;
    if (xo6[i] == 1) xotime[i] = timePFSobs[i] + 105;
    
    if (xo1[i] == 1 || xo2[i] == 1 || xo3[i] == 1 || xo4[i] == 1 ||
        xo5[i] == 1 || xo6[i] == 1) {
      xo[i] = 1;
    }
    
    // apply switch effect first through reducing prob. of catastrophic event
    if (xo[i] == 1 && (xotime[i] < cattime[i] || catevent[i] == INT_MIN)) {
      xoprecat[i] = 1;
    }
    
    // after treatment switching for control patients, prob of developing 
    // metastatic disease only depends on baseline prognosis
    if (trtrand[i] == 0 && xo[i] == 1 && bprog[i] == 1 && xoprecat[i] == 1) {
      probcat[i] = pcattrtbprog;
    }
    
    if (trtrand[i] == 0 && xo[i] == 1 && bprog[i] == 0 && xoprecat[i] == 1) {
      probcat[i] = pcattrt;
    }
    
    // regenerate metastatic disease status after treatment switching
    if (xo[i] == 1 && xoprecat[i] == 1 && xotime[i] == timePFSobs[i]) {
      cat2[i] = INT_MIN;
    }
    
    if (xo[i] == 1 && xoprecat[i] == 1 && (xotime[i] == timePFSobs[i] ||
        xotime[i] == timePFSobs[i] + 21)) {
      cat3[i] = INT_MIN;
    }
    
    if (xo[i] == 1 && xoprecat[i] == 1 && (xotime[i] == timePFSobs[i] ||
        xotime[i] == timePFSobs[i] + 21 || xotime[i] == timePFSobs[i] + 42)) {
      cat4[i] = INT_MIN;
    }
    
    if (xo[i] == 1 && xoprecat[i] == 1 && (xotime[i] == timePFSobs[i] ||
        xotime[i] == timePFSobs[i] + 21 || xotime[i] == timePFSobs[i] + 42 ||
        xotime[i] == timePFSobs[i] + 63)) {
      cat5[i] = INT_MIN;
    }
    
    if (progressed[i] == 1 && timeOS[i] > timePFSobs[i] + 21 &&
        xo[i] == 1 && xotime[i] == timePFSobs[i]) {
      cat2[i] = (unif(rng) <= probcat[i]) ? 1 : 0;
    }
    
    if (cat2[i] == 0 && progressed[i] == 1 &&
        timeOS[i] > timePFSobs[i] + 42 && xo[i] == 1 &&
        (xotime[i] == timePFSobs[i] || xotime[i] == timePFSobs[i] + 21)) {
      cat3[i] = (unif(rng) <= probcat[i]) ? 1 : 0;
    }
    
    if (cat2[i] == 0 && cat3[i] == 0 && progressed[i] == 1 &&
        timeOS[i] > timePFSobs[i] + 63 && xo[i] == 1 &&
        (xotime[i] == timePFSobs[i] || xotime[i] == timePFSobs[i] + 21 ||
        xotime[i] == timePFSobs[i] + 42)) {
      cat4[i] = (unif(rng) <= probcat[i]) ? 1 : 0;
    }
    
    if (cat2[i] == 0 && cat3[i] == 0 && cat4[i] == 0 && progressed[i] == 1 &&
        timeOS[i] > timePFSobs[i] + 84 && xo[i] == 1 &&
        (xotime[i] == timePFSobs[i] || xotime[i] == timePFSobs[i] + 21 ||
        xotime[i] == timePFSobs[i] + 42 || xotime[i] == timePFSobs[i] + 63)) {
      cat5[i] = (unif(rng) <= probcat[i]) ? 1 : 0;
    }
    
    if (xo[i] == 1 && xoprecat[i] == 1) {
      catevent[i] = INT_MIN;
      cattime[i] = NaN;
      catOSloss[i] = NaN;
    }
    
    if (cat2[i] == 1 || cat3[i] == 1 || cat4[i] == 1 || cat5[i] == 1) {
      catevent[i] = 1;
    }
    
    if (cat2[i] == 1 && xo[i] == 1 && xotime[i] == timePFSobs[i]) {
      cattime[i] = timePFSobs[i] + 21;
    }
    
    if (cat3[i] == 1 && xo[i] == 1 && (xotime[i] == timePFSobs[i] ||
        xotime[i] == timePFSobs[i] + 21)) {
      cattime[i] = timePFSobs[i] + 42;
    }
    
    if (cat4[i] == 1 && xo[i] == 1 && (xotime[i] == timePFSobs[i] ||
        xotime[i] == timePFSobs[i] + 21 ||
        xotime[i] == timePFSobs[i] + 42)) {
      cattime[i] = timePFSobs[i] + 63;
    }
    
    if (cat5[i] == 1 && xo[i] == 1 && (xotime[i] == timePFSobs[i] ||
        xotime[i] == timePFSobs[i] + 21 || xotime[i] == timePFSobs[i] + 42 ||
        xotime[i] == timePFSobs[i] + 63)) {
      cattime[i] = timePFSobs[i] + 84;
    }
    
    // modify survival time after developing metastatic disease
    if (catevent[i] == 1 && xoprecat[i] == 1) {
      catOSloss[i] = timeOS5[i] - cattime[i];
      catOSloss[i] = std::round(catOSloss[i] * catmult);
      timeOS[i] = catOSloss[i] + cattime[i];
    }
    
    // if the patient did not develop metastatic disease after all, use
    // original survival time
    if (catevent[i] != 1 && xoprecat[i] == 1) {
      timeOS[i] = timeOS5[i];
    }
    
    // apply chance of catastrophic event for people who don't have cat event
    // or have later cat event and now live past an additional observation.
    // Only replacing for extra observations - others don't change
    if (progressed[i] == 1 && timeOS[i] > timePFSobs[i] + 21 &&
        xo[i] == 1 && xotime[i] == timePFSobs[i] && cat2[i] == INT_MIN) {
      extra2v1[i] = 1;
    }
    
    if (cat2[i] == 0 && progressed[i] == 1 &&
        timeOS[i] > timePFSobs[i] + 42 &&
        xo[i] == 1 && (xotime[i] == timePFSobs[i] ||
        xotime[i] == timePFSobs[i] + 21) && cat3[i] == INT_MIN) {
      extra3v1[i] = 1;
    }
    
    if (cat2[i] == 0 && cat3[i] == 0 && progressed[i] == 1 &&
        timeOS[i] > timePFSobs[i] + 63 &&
        xo[i] == 1 && (xotime[i] == timePFSobs[i] ||
        xotime[i] == timePFSobs[i] + 21 ||
        xotime[i] == timePFSobs[i] + 42) && cat4[i] == INT_MIN) {
      extra4v1[i] = 1;
    }
    
    if (cat2[i] == 0 && cat3[i] == 0 && cat4[i] == 0 && progressed[i] == 1 &&
        timeOS[i] > timePFSobs[i] + 84 &&
        xo[i] == 1 && (xotime[i] == timePFSobs[i] ||
        xotime[i] == timePFSobs[i] + 21 || xotime[i] == timePFSobs[i] + 42 ||
        xotime[i] == timePFSobs[i] + 63) && cat5[i] == INT_MIN) {
      extra5v1[i] = 1;
    }
    
    if (extra2v1[i] == 1 || extra3v1[i] == 1 || extra4v1[i] == 1 ||
        extra5v1[i] == 1) {
      extraobsv1[i] = 1;
    }
    
    if (progressed[i] == 1 && timeOS[i] > timePFSobs[i] + 21 && 
        xo[i] == 1 && xotime[i] == timePFSobs[i] && cat2[i] == INT_MIN) {
      cat2[i] = (unif(rng) <= probcat[i]) ? 1 : 0;
    }
    
    if (cat2[i] == 0 && progressed[i] == 1 &&
        timeOS[i] > timePFSobs[i] + 42 && 
        xo[i] == 1 && (xotime[i] == timePFSobs[i] || 
        xotime[i] == timePFSobs[i] +21) && cat3[i] == INT_MIN) {
      cat3[i] = (unif(rng) <= probcat[i]) ? 1 : 0;
    }
    
    if (cat2[i] == 0 && cat3[i] == 0 && progressed[i] == 1 &&
        timeOS[i] > timePFSobs[i] + 63 && 
        xo[i] == 1 && (xotime[i] == timePFSobs[i] || 
        xotime[i] == timePFSobs[i] + 21 || 
        xotime[i] == timePFSobs[i] + 42) && cat4[i] == INT_MIN) {
      cat4[i] = (unif(rng) <= probcat[i]) ? 1 : 0;
    }
    
    if (cat2[i] == 0 && cat3[i] == 0 && cat4[i] == 0 && progressed[i] == 1 &&
        timeOS[i] > timePFSobs[i] + 84 && 
        xo[i] == 1 && (xotime[i] == timePFSobs[i] || 
        xotime[i] == timePFSobs[i] + 21 || xotime[i] == timePFSobs[i] + 42 ||
        xotime[i] == timePFSobs[i] + 63) && cat5[i] == INT_MIN) {
      cat5[i] = (unif(rng) <= probcat[i]) ? 1 : 0;
    }
    
    if (xo[i] == 1 && xoprecat[i] == 1 && extraobsv1[i] == 1) {
      cattime[i] = NaN;
    }
    
    if (cat2[i] == 1 && xo[i] == 1 && xotime[i] == timePFSobs[i] &&
        extraobsv1[i] == 1) {
      cattime[i] = timePFSobs[i] + 21;
    }
    
    if (cat3[i] == 1 && xo[i] == 1 && (xotime[i] == timePFSobs[i] ||
        xotime[i] == timePFSobs[i] + 21) && extraobsv1[i] == 1) {
      cattime[i] = timePFSobs[i] + 42;
    }
    
    if (cat4[i] == 1 && xo[i] == 1 && (xotime[i] == timePFSobs[i] ||
        xotime[i] == timePFSobs[i] + 21 ||
        xotime[i] == timePFSobs[i] + 42) && extraobsv1[i] == 1) {
      cattime[i] = timePFSobs[i] + 63;
    }
    
    if (cat5[i] == 1 && xo[i] == 1 && (xotime[i] == timePFSobs[i] ||
        xotime[i] == timePFSobs[i] + 21 || xotime[i] == timePFSobs[i] + 42 ||
        xotime[i] == timePFSobs[i] + 63) && extraobsv1[i] == 1) {
      cattime[i] = timePFSobs[i] + 84;
    }
    
    if (xo[i] == 1 && xoprecat[i] == 1 && extraobsv1[i] == 1) {
      catevent[i] = INT_MIN;
    }
    
    if ((cat2[i] == 1 || cat3[i] == 1 || cat4[i] == 1 || cat5[i] == 1) &&
        extraobsv1[i] == 1) {
      catevent[i] = 1;
    }
    
    // remember, these are only patients who hadn't had a cat event before,
    // and now may have done in their additional observations. So survival
    // times cannot go up, this only allows for possible additional events
    if (xo[i] == 1 && xoprecat[i] == 1 && extraobsv1[i] == 1) {
      catOSloss[i] = NaN;
    }
    
    if (catevent[i] == 1 && xoprecat[i] == 1 && extraobsv1[i] == 1) {
      catOSloss[i] = timeOS5[i] - cattime[i];
      catOSloss[i] = std::round(catOSloss[i]*catmult);
      timeOS[i] = catOSloss[i] + cattime[i];
    }
    
    if (catevent[i] != 1 && xoprecat[i] == 1 && extraobsv1[i] == 1) {
      timeOS[i] = timeOS5[i];
    }
    
    // Apply switch effect second through xomult
    if (trtrand[i] == 0 && xo[i] == 1) {
      xoOSgainobs[i] = timeOS[i] - xotime[i];
      xoOSgainobs[i] = std::round(xoOSgainobs[i] * xomult);
    }
    
    timeOS2[i] = xo[i] == 1 ? xoOSgainobs[i] + xotime[i] : timeOS[i];
    
    // apply chance of catastrophic event for xo patients who now live past
    // an additional observation. Only replacing for extra observations in
    // switchers, so no additional chance to switch and OS can't increase
    // as currently in these observations there will zero chance of cat event
    if (progressed[i] == 1 && timeOS2[i] > timePFSobs[i] + 21 &&
        xo[i] == 1 && xotime[i] == timePFSobs[i] && cat2[i] == INT_MIN) {
      extra2v2[i] = 1;
    }
    
    if (cat2[i] == 0 && progressed[i] == 1 &&
        timeOS2[i] > timePFSobs[i] + 42 &&
        xo[i] == 1 && (xotime[i] == timePFSobs[i] ||
        xotime[i] == timePFSobs[i] + 21) && cat3[i] == INT_MIN) {
      extra3v2[i] = 1;
    }
    
    if (cat2[i] == 0 && cat3[i] == 0 && progressed[i] == 1 &&
        timeOS2[i] > timePFSobs[i] + 63 &&
        xo[i] == 1 && (xotime[i] == timePFSobs[i] ||
        xotime[i] == timePFSobs[i] + 21 ||
        xotime[i] == timePFSobs[i] + 42) && cat4[i] == INT_MIN) {
      extra4v2[i] = 1;
    }
    
    if (cat2[i] == 0 && cat3[i] == 0 && cat4[i] == 0 && progressed[i] == 1 &&
        timeOS2[i] > timePFSobs[i] + 84 &&
        xo[i] == 1 && (xotime[i] == timePFSobs[i] ||
        xotime[i] == timePFSobs[i] + 21 || xotime[i] == timePFSobs[i] + 42 ||
        xotime[i] == timePFSobs[i] + 63) && cat5[i] == INT_MIN) {
      extra5v2[i] = 1;
    }
    
    if (extra2v2[i] == 1 || extra3v2[i] == 1 || extra4v2[i] == 1 ||
        extra5v2[i] == 1) {
      extraobsv2[i] = 1;
    }
    
    if (progressed[i] == 1 && timeOS2[i] > timePFSobs[i] + 21 &&
        xo[i] == 1 && xotime[i] == timePFSobs[i] && cat2[i] == INT_MIN) {
      cat2[i] = (unif(rng) <= probcat[i]) ? 1 : 0;
    }
    
    if (cat2[i] == 0 && progressed[i] == 1 &&
        timeOS2[i] > timePFSobs[i] + 42 && 
        xo[i] == 1 && (xotime[i] == timePFSobs[i] || 
        xotime[i] == timePFSobs[i] +21) && cat3[i] == INT_MIN) {
      cat3[i] = (unif(rng) <= probcat[i]) ? 1 : 0;
    }
    
    if (cat2[i] == 0 && cat3[i] == 0 && progressed[i] == 1 &&
        timeOS2[i] > timePFSobs[i] + 63 && 
        xo[i] == 1 && (xotime[i] == timePFSobs[i] || 
        xotime[i] == timePFSobs[i] + 21 ||
        xotime[i] == timePFSobs[i] + 42) && cat4[i] == INT_MIN) {
      cat4[i] = (unif(rng) <= probcat[i]) ? 1 : 0;
    }
    
    if (cat2[i] == 0 && cat3[i] == 0 && cat4[i] == 0 && progressed[i] == 1 &&
        timeOS2[i] > timePFSobs[i] + 84 && 
        xo[i] == 1 && (xotime[i] == timePFSobs[i] || 
        xotime[i] == timePFSobs[i] + 21 || xotime[i] == timePFSobs[i] + 42 ||
        xotime[i] == timePFSobs[i] + 63) && cat5[i] == INT_MIN) {
      cat5[i] = (unif(rng) <= probcat[i]) ? 1 : 0;
    }
    
    if (xo[i] == 1 && xoprecat[i] == 1 && extraobsv2[i] == 1) {
      cattime[i] = NaN;
    }
    
    if (cat2[i] == 1 && xo[i] == 1 && xotime[i] == timePFSobs[i] &&
        extraobsv2[i] == 1) {
      cattime[i] = timePFSobs[i] + 21;
    }
    
    if (cat3[i] == 1 && xo[i] == 1 && (xotime[i] == timePFSobs[i] ||
        xotime[i] == timePFSobs[i] + 21) && extraobsv2[i] == 1) {
      cattime[i] = timePFSobs[i] + 42;
    }
    
    if (cat4[i] == 1 && xo[i] == 1 && (xotime[i] == timePFSobs[i] ||
        xotime[i] == timePFSobs[i] + 21 ||
        xotime[i] == timePFSobs[i] + 42) && extraobsv2[i] == 1) {
      cattime[i] = timePFSobs[i] + 63;
    }
    
    if (cat5[i] == 1 && xo[i] == 1 && (xotime[i] == timePFSobs[i] ||
        xotime[i] == timePFSobs[i] + 21 || xotime[i] == timePFSobs[i] + 42 ||
        xotime[i] == timePFSobs[i] + 63) && extraobsv2[i] == 1) {
      cattime[i] = timePFSobs[i] + 84;
    }
    
    if (xo[i] == 1 && xoprecat[i] == 1 && extraobsv2[i] == 1) {
      catevent[i] = INT_MIN;
    }
    
    if ((cat2[i] == 1 || cat3[i] == 1 || cat4[i] == 1 || cat5[i] == 1) &&
        extraobsv2[i] == 1) {
      catevent[i] = 1;
    }
    
    // remember, these are only patients who hadn't had a cat event before,
    // and now may have done in their additional observations, which
    // previously they did not live past. So if they have an event,
    // cat time will be >OS5. Hence only makes sense to apply cat event
    // survival reduction to extended survival beyond this point, i.e., OS2
    if (xo[i] == 1 && xoprecat[i] == 1 && extraobsv2[i] == 1) {
      catOSloss[i] = NaN;
    }
    
    if (catevent[i] == 1 && xoprecat[i] == 1 && extraobsv2[i] == 1) {
      catOSloss[i] = timeOS2[i] - cattime[i];
      catOSloss[i] = std::round(catOSloss[i]*catmult);
      timeOS2[i] = catOSloss[i] + cattime[i];
    }
    
    // don't replace OS2 for those with an extraobsv2 who did not have a cat
    // event because OS has not changed for them
  }
  
  // apply censoring
  std::vector<int> died(n);
  for (int i = 0; i < n; ++i) {
    if (timeOS2[i] <= admin && dead[i] == 1) {
      died[i] = 1;
    } else {
      died[i] = 0;
    }
    
    if (timeOS2[i] > admin) timeOS2[i] = admin;
    
    if (xotime[i] >= admin) {
      xo[i] = INT_MIN;
      xotime[i] = NaN;
    }
    
    if (cattime[i] >= admin) {
      catevent[i] = INT_MIN;
      cattime[i] = NaN;
    }
  }
  
  // create panel
  std::vector<double> zero(n);
  double maxtime = *std::max_element(timeOS2.begin(), timeOS2.end());
  int kmax = static_cast<int>(std::ceil(maxtime / 21));
  std::vector<double> cut(kmax);
  for (int k = 0; k < kmax; ++k) {
    cut[k] = k * 21;
  }
  
  for (int i = 0; i < n; i++) {
    if (progressed[i] == INT_MIN) progressed[i] = 0;
    if (catevent[i] == INT_MIN) catevent[i] = 0;
    if (xo[i] == INT_MIN) xo[i] = 0;
  }
  
  DataFrameCpp a = survsplitcpp(zero, timeOS2, cut);
  int n2 = static_cast<int>(a.nrows());
  std::vector<int> q2 = a.get<int>("row");
  std::vector<double> tstart = a.get<double>("start");
  std::vector<double> tstop = a.get<double>("end");
  std::vector<int> censor = a.get<int>("censor");
  
  std::vector<int> id2 = subset(id, q2);
  std::vector<int> trtrand2 = subset(trtrand, q2);
  std::vector<int> bprog2 = subset(bprog, q2);
  std::vector<int> died2 = subset(died, q2);
  for (int i = 0; i < n2; ++i) if (censor[i] == 1) died2[i] = 0;
  
  std::vector<double> timeOS8 = subset(timeOS2, q2);
  std::vector<int> died8 = subset(died, q2);
  std::vector<int> progressed2 = subset(progressed, q2);
  std::vector<double> timePFSobs2 = subset(timePFSobs, q2);
  std::vector<int> catevent2 = subset(catevent, q2);
  std::vector<double> cattime2 = subset(cattime, q2);
  std::vector<int> xoo2 = subset(xo, q2);
  std::vector<double> xotime2 = subset(xotime, q2);
  
  // make time-dependent covariates for progression, cat event, and switch
  std::vector<int> progtdc(n2), cattdc(n2), xotdc(n2);
  for (int i = 0; i < n2; ++i) {
    if (progressed2[i] == 1 && tstart[i] >= timePFSobs2[i]) progtdc[i] = 1;
    if (catevent2[i] == 1 && tstart[i] >= cattime2[i]) cattdc[i] = 1;
    if (xoo2[i] == 1 && tstart[i] >= xotime2[i]) xotdc[i] = 1;
  }
  
  // create the lagged value of cattdc
  std::vector<int> idx(1, 0); // first observation within an id
  for (int i = 1; i < n2; ++i) {
    if (id2[i] != id2[i-1]) {
      idx.push_back(i);
    }
  }
  
  ListCpp result;
  
  DataFrameCpp sumstat;
  sumstat.push_back(simtrueconstmean, "simtrueconstmean");
  sumstat.push_back(simtrueconstlb, "simtrueconstlb");
  sumstat.push_back(simtrueconstub, "simtrueconstub");
  sumstat.push_back(simtrueconstse, "simtrueconstse");
  sumstat.push_back(simtrueexpstmean, "simtrueexpstmean");
  sumstat.push_back(simtrueexpstlb, "simtrueexpstlb");
  sumstat.push_back(simtrueexpstub, "simtrueexpstub");
  sumstat.push_back(simtrueexpstse, "simtrueexpstse");
  sumstat.push_back(simtrue_coxwbprog_hr, "simtrue_coxwbprog_hr");
  sumstat.push_back(simtrue_cox_hr, "simtrue_cox_hr");
  sumstat.push_back(simtrue_aftwbprog_af, "simtrue_aftwbprog_af");
  sumstat.push_back(simtrue_aft_af, "simtrue_aft_af");
  
  DataFrameCpp paneldata;
  paneldata.push_back(id2, "id");
  paneldata.push_back(trtrand2, "trtrand");
  paneldata.push_back(bprog2, "bprog");
  paneldata.push_back(tstart, "tstart");
  paneldata.push_back(tstop, "tstop");
  paneldata.push_back(died2, "event");
  paneldata.push_back(timeOS8, "timeOS");
  paneldata.push_back(died8, "died");
  paneldata.push_back(progressed2, "progressed");
  paneldata.push_back(timePFSobs2, "timePFSobs");
  paneldata.push_back(progtdc, "progtdc");
  paneldata.push_back(catevent2, "catevent");
  paneldata.push_back(cattime2, "cattime");
  paneldata.push_back(cattdc, "cattdc");
  paneldata.push_back(xoo2, "xo");
  paneldata.push_back(xotime2, "xotime");
  paneldata.push_back(xotdc, "xotdc");
  paneldata.push_back(admin, "censor_time");
  
  result.push_back(sumstat, "sumstat");
  result.push_back(paneldata, "paneldata");
  return Rcpp::wrap(result);
}
