#include <RcppThread.h>
#include <Rcpp.h>

#include <boost/random.hpp>

#include "survival_analysis.h"
#include "utilities.h"
#include "dataframe_list.h"
#include "thread_utils.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <random>
#include <stdexcept>
#include <vector>


// counterfactual untreated survival times and event indicators
DataFrameCpp f_untreated(
    const double psi,
    const int n,
    const std::vector<int>& id,
    const std::vector<double>& time,
    const std::vector<int>& event,
    const std::vector<int>& treat,
    const std::vector<double>& rx,
    const std::vector<double>& censor_time,
    const int recensor_type,
    const bool autoswitch) {
  
  std::vector<int> nonswitchers;
  nonswitchers.reserve(n);
  for (int i = 0; i < n; ++i) {
    if (treat[i] == 1) {
      if (rx[i] == 1.0) nonswitchers.push_back(i);
    } else {
      if (rx[i] == 0.0) nonswitchers.push_back(i);
    }
  }
  
  double a = std::exp(psi);
  std::vector<double> u_star(n);
  for (int i = 0; i < n; ++i) 
    u_star[i] = time[i] * ((1 - rx[i]) + rx[i] * a);
  std::vector<double> t_star = u_star;
  std::vector<int> d_star = event;
  
  if (recensor_type == 1) { // recensor all patients
    std::vector<double> c_star(n);
    double multiplier = std::min(1.0, a);
    for (int i=0; i<n; i++) 
      c_star[i] = censor_time[i] * multiplier;
    
    if (autoswitch) {
      std::vector<double> rx1, rx0;
      rx1.reserve(n); rx0.reserve(n);
      for (int i = 0; i < n; ++i) {
        if (treat[i] == 1) rx1.push_back(rx[i]);
        else rx0.push_back(rx[i]);
      }
      
      if (all_equal(rx1, 1.0, 0.0)) {
        for (int i = 0; i < n; ++i) {
          if (treat[i] == 1) c_star[i] = POS_INF;
        }
      }
      
      if (all_equal(rx0, 0.0, 0.0)) {
        for (int i = 0; i < n; ++i) {
          if (treat[i] == 0) c_star[i] = POS_INF;
        }
      }
    }
    
    for (int i = 0; i < n; ++i) {
      t_star[i] = std::min(u_star[i], c_star[i]);
      if (c_star[i] < u_star[i]) d_star[i] = 0;
    }
  } else if (recensor_type == 2) { // recensor only switchers
    std::vector<double> c_star(n);
    double multiplier = std::min(1.0, a);
    for (int i = 0; i < n; ++i) 
      c_star[i] = censor_time[i] * multiplier;
    
    if (autoswitch) {
      std::vector<double> rx1, rx0;
      rx1.reserve(n); rx0.reserve(n);
      for (int i = 0; i < n; ++i) {
        if (treat[i] == 1) rx1.push_back(rx[i]);
        else rx0.push_back(rx[i]);
      }
      
      if (all_equal(rx1, 1.0, 0.0)) {
        for (int i = 0; i < n; ++i) {
          if (treat[i] == 1) c_star[i] = POS_INF;
        }
      }
      
      if (all_equal(rx0, 0.0, 0.0)) {
        for (int i = 0; i < n; ++i) {
          if (treat[i] == 0) c_star[i] = POS_INF;
        }
      }
    }
    
    for (int i : nonswitchers) c_star[i] = POS_INF;
    
    for (int i = 0; i < n; ++i) {
      t_star[i] = std::min(u_star[i], c_star[i]);
      if (c_star[i] < u_star[i]) d_star[i] = 0;
    }
  } else if (recensor_type == 3) { // recensor only switchers with 
    // projected latent event time beyond the end of study
    std::vector<double> c_star = censor_time;
    
    if (autoswitch) {
      std::vector<double> rx1, rx0;
      rx1.reserve(n); rx0.reserve(n);
      for (int i = 0; i < n; ++i) {
        if (treat[i] == 1) rx1.push_back(rx[i]);
        else rx0.push_back(rx[i]);
      }
      
      if (all_equal(rx1, 1.0, 0.0)) {
        for (int i = 0; i < n; ++i) {
          if (treat[i] == 1) c_star[i] = POS_INF;
        }
      }
      
      if (all_equal(rx0, 0.0, 0.0)) {
        for (int i = 0; i < n; ++i) {
          if (treat[i] == 0) c_star[i] = POS_INF;
        }
      }
    }
    
    for (int i : nonswitchers) c_star[i] = POS_INF;
    
    for (int i = 0; i < n; ++i) {
      t_star[i] = std::min(u_star[i], c_star[i]);
      if (c_star[i] < u_star[i]) d_star[i] = 0;
    }
  }
  
  DataFrameCpp result;
  result.push_back(id, "uid");
  result.push_back(t_star, "t_star");
  result.push_back(d_star, "d_star");
  result.push_back(treat, "treated");
  return result;
}


// counterfactual unswitched survival times and event indicators
DataFrameCpp f_unswitched(
    const double psi,
    const int n,
    const std::vector<int>& id,
    const std::vector<double>& time,
    const std::vector<int>& event,
    const std::vector<int>& treat,
    const std::vector<double>& rx,
    const std::vector<double>& censor_time,
    const int recensor_type,
    const bool autoswitch) {
  
  std::vector<int> nonswitchers;
  nonswitchers.reserve(n);
  for (int i = 0; i < n; ++i) {
    if (treat[i] == 1) {
      if (rx[i] == 1.0) nonswitchers.push_back(i);
    } else {
      if (rx[i] == 0.0) nonswitchers.push_back(i);
    }
  }
  
  double a = std::exp(psi), a1 = std::exp(-psi);
  std::vector<double> u_star(n), t_star(n);
  std::vector<int> d_star(n);
  for (int i = 0; i < n; ++i) {
    if (treat[i] == 0) {
      u_star[i] = time[i] * ((1 - rx[i]) + rx[i] * a);
    } else {
      u_star[i] = time[i] * (rx[i] + (1 - rx[i]) * a1);
    }
    t_star[i] = u_star[i];
    d_star[i] = event[i];
  }
  
  if (recensor_type == 1) { // recensor all patients
    std::vector<double> c_star(n);
    double c0 = std::min(1.0, a), c1 = std::min(1.0, a1);
    for (int i = 0; i < n; ++i) {
      if (treat[i] == 0) c_star[i] = censor_time[i] * c0;
      else c_star[i] = censor_time[i] * c1;
    }
    
    if (autoswitch) {
      std::vector<double> rx1, rx0;
      rx1.reserve(n); rx0.reserve(n);
      for (int i = 0; i < n; ++i) {
        if (treat[i] == 1) rx1.push_back(rx[i]);
        else rx0.push_back(rx[i]);
      }
      
      if (all_equal(rx1, 1.0, 0.0)) {
        for (int i = 0; i < n; ++i) {
          if (treat[i] == 1) c_star[i] = POS_INF;
        }
      }
      
      if (all_equal(rx0, 0.0, 0.0)) {
        for (int i = 0; i < n; ++i) {
          if (treat[i] == 0) c_star[i] = POS_INF;
        }
      }
    }
    
    for (int i = 0; i < n; ++i) {
      t_star[i] = std::min(u_star[i], c_star[i]);
      if (c_star[i] < u_star[i]) d_star[i] = 0;
    }
  } else if (recensor_type == 2) { // recensor only switchers
    std::vector<double> c_star(n);
    double c0 = std::min(1.0, a), c1 = std::min(1.0, a1);
    for (int i = 0; i < n; ++i) {
      if (treat[i] == 0) c_star[i] = censor_time[i] * c0;
      else c_star[i] = censor_time[i] * c1;
    }
    
    if (autoswitch) {
      std::vector<double> rx1, rx0;
      rx1.reserve(n); rx0.reserve(n);
      for (int i = 0; i < n; ++i) {
        if (treat[i] == 1) rx1.push_back(rx[i]);
        else rx0.push_back(rx[i]);
      }
      
      if (all_equal(rx1, 1.0, 0.0)) {
        for (int i = 0; i < n; ++i) {
          if (treat[i] == 1) c_star[i] = POS_INF;
        }
      }
      
      if (all_equal(rx0, 0.0, 0.0)) {
        for (int i = 0; i < n; ++i) {
          if (treat[i] == 0) c_star[i] = POS_INF;
        }
      }
    }
    
    for (int i : nonswitchers) c_star[i] = POS_INF;
    
    for (int i = 0; i < n; ++i) {
      t_star[i] = std::min(u_star[i], c_star[i]);
      if (c_star[i] < u_star[i]) d_star[i] = 0;
    }
  } else if (recensor_type == 3) { // recensor only switchers with 
    // projected latent event time beyond the end of study
    std::vector<double> c_star = censor_time;
    
    if (autoswitch) {
      std::vector<double> rx1, rx0;
      rx1.reserve(n); rx0.reserve(n);
      for (int i = 0; i < n; ++i) {
        if (treat[i] == 1) rx1.push_back(rx[i]);
        else rx0.push_back(rx[i]);
      }
      
      if (all_equal(rx1, 1.0, 0.0)) {
        for (int i = 0; i < n; ++i) {
          if (treat[i] == 1) c_star[i] = POS_INF;
        }
      }
      
      if (all_equal(rx0, 0.0, 0.0)) {
        for (int i = 0; i < n; ++i) {
          if (treat[i] == 0) c_star[i] = POS_INF;
        }
      }
    }
    
    for (int i : nonswitchers) c_star[i] = POS_INF;
    
    for (int i = 0; i < n; ++i) {
      t_star[i] = std::min(u_star[i], c_star[i]);
      if (c_star[i] < u_star[i]) d_star[i] = 0;
    }
  }
  
  DataFrameCpp result;
  result.push_back(id, "uid");
  result.push_back(t_star, "t_star");
  result.push_back(d_star, "d_star");
  result.push_back(treat, "treated");
  return result;
}


double f_est_psi_rpsftm(
    const double psi,
    const int n,
    const std::vector<int>& id,
    const std::vector<int>& stratum,
    const std::vector<double>& time,
    const std::vector<int>& event,
    const std::vector<int>& treat,
    const std::vector<double>& rx,
    const std::vector<double>& censor_time,
    const double treat_modifier,
    const int recensor_type,
    const bool autoswitch,
    const double target) {
  
  DataFrameCpp Sstar = f_untreated(
    psi * treat_modifier, n, id, time, event, treat,
    rx, censor_time, recensor_type, autoswitch);
  
  Sstar.push_back(stratum, "ustratum");
  DataFrameCpp df = lrtestcpp(Sstar, {"ustratum"}, "treated", "t_star", "", "d_star");
  
  double z = df.get<double>("logRankZ")[0];
  return z - target;
}


//' @title Simulation Study to Evaluate Recensoring Rules in RPSFTM
//' 
//' @description 
//' Simulates datasets to evaluate the performance of various recensoring 
//' strategies under the Rank Preserving Structural Failure Time Model 
//' (RPSFTM) for handling treatment switching in survival analysis.
//'
//' @param nsim Number of simulated datasets.
//' @param n Number of subjects per simulation.
//' @param shape Shape parameter of the Weibull distribution for time to 
//'   death.
//' @param scale Scale parameter of the Weibull distribution for time to 
//'   death in the control group.
//' @param gamma Rate parameter of the exponential distribution for random 
//'   dropouts in the control group.
//' @param tfmin Minimum planned follow-up time (in days).
//' @param tfmax Maximum planned follow-up time (in days).
//' @param psi Log time ratio of death time for control vs experimental 
//'   treatment.
//' @param omega Log time ratio of dropout time for control vs experimental 
//'   treatment.
//' @param pswitch Probability of treatment switching at disease progression.
//' @param a Shape parameter 1 of the Beta distribution for time to disease 
//'   progression as a fraction of time to death.
//' @param b Shape parameter 2 of the Beta distribution for time to disease 
//'   progression.
//' @param low_psi Lower bound for the search interval of the causal 
//'   parameter \eqn{\psi}.
//' @param hi_psi Upper bound for the search interval of the causal 
//'   parameter \eqn{\psi}.
//' @param treat_modifier Sensitivity parameter modifying the constant 
//'   treatment effect assumption.
//' @param recensor_type Type of recensoring to apply:
//'   \itemize{
//'     \item 0: No recensoring
//'     \item 1: Recensor all control-arm subjects
//'     \item 2: Recensor only switchers in the control arm
//'     \item 3: Recensor only control-arm switchers whose counterfactual 
//'              survival exceeds the planned follow-up time
//'   }
//' @param admin_recensor_only Logical. If \code{TRUE}, recensoring is 
//'   applied only to administrative censoring times. 
//'   If \code{FALSE}, it is also applied to dropout times.
//' @param autoswitch Logical. If \code{TRUE}, disables recensoring in arms 
//'   without any treatment switching.
//' @param alpha Significance level for confidence interval calculation 
//'   (default is 0.05).
//' @param ties Method for handling tied event times in the Cox model. 
//'   Options are \code{"efron"} (default) or \code{"breslow"}.
//' @param tol Convergence tolerance for root-finding in estimation of 
//'   \eqn{\psi}.
//' @param boot Logical. If \code{TRUE}, bootstrap is used to estimate 
//'   the confidence interval for the hazard ratio. If \code{FALSE}, 
//'   the confidence interval is matched to the log-rank p-value.
//' @param n_boot Number of bootstrap samples, used only if 
//'   \code{boot = TRUE}.
//' @param seed Optional. Random seed for reproducibility. 
//'
//' @return A data frame summarizing the simulation results, including:
//' \itemize{
//'   \item \code{recensor_type}, \code{admin_recensor_only}: Settings 
//'         used in the simulation.
//'   \item Event rates: \code{p_event_1}, \code{p_dropout_1}, 
//'         \code{p_admin_censor_1}, \code{p_event_0}, 
//'         \code{p_dropout_0}, \code{p_admin_censor_0}.
//'   \item Progression and switching: \code{p_pd_0}, \code{p_swtrt_0},
//'         \code{p_recensored_0}.
//'   \item Causal parameter (\eqn{\psi}) estimates: \code{psi}, 
//'         \code{psi_est}, \code{psi_bias}, 
//'         \code{psi_se}, \code{psi_mse}.
//'   \item Log hazard ratio estimates: \code{loghr}, \code{loghr_est}, 
//'         \code{loghr_se}, \code{loghr_mse}.
//'   \item Hazard ratio metrics: \code{hr}, \code{hr_est} (geometric mean), 
//'         \code{hr_pctbias} (percent bias).
//'   \item Standard errors of log hazard ratio: \code{loghr_se_cox}, 
//'         \code{loghr_se_lr}, \code{loghr_se_boot}.
//'   \item Coverage probabilities: \code{hr_ci_cover_cox}, 
//'         \code{hr_ci_cover_lr}, \code{hr_ci_cover_boot}.
//' }
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' \donttest{
//' result <- recensor_sim_rpsftm(
//'   nsim = 10, n = 400, shape = 1.5, scale = exp(6.3169),
//'   gamma = 0.001, tfmin = 407.5, tfmax = 407.5,
//'   psi = log(0.5) / 1.5, omega = log(1), pswitch = 0.7,
//'   a = 2, b = 4, low_psi = -5, hi_psi = 5,
//'   treat_modifier = 1, recensor_type = 1,
//'   admin_recensor_only = TRUE, autoswitch = TRUE,
//'   alpha = 0.05, tol = 1e-6, boot = TRUE,
//'   n_boot = 10, seed = 314159)
//' }
//'
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame recensor_sim_rpsftm(
    const int nsim = 100,
    const int n = 400,
    const double shape = 1.5,
    const double scale = 553.9,
    const double gamma = 0.001,
    const double tfmin = 407.5,
    const double tfmax = 407.5,
    const double psi = -0.4621,
    const double omega = 0.0,
    const double pswitch = 0.7,
    const double a = 2,
    const double b = 4,
    const double low_psi = -5,
    const double hi_psi = 5,
    const double treat_modifier = 1,
    const int recensor_type = 1,
    const bool admin_recensor_only = true,
    const bool autoswitch = true,
    const double alpha = 0.05,
    const std::string ties = "efron",
    const double tol = 1.0e-6,
    const bool boot = true,
    const int n_boot = 100,
    const int seed = 0) {
  
  std::string str1 = "recensor_type = ";
  str1 += std::to_string(recensor_type) + ", admin_recensor_only = " + 
    std::to_string(admin_recensor_only) + "\n";
  
  boost::random::mt19937_64 rng(seed);
  boost::random::uniform_real_distribution<double> unif(0.0, 1.0);
  boost::random::exponential_distribution<double> expo(gamma);
  boost::random::weibull_distribution<double> weib(shape, scale);
  boost::random::gamma_distribution<double> g1(a, 1.0);
  boost::random::gamma_distribution<double> g2(b, 1.0);
  
  double twidth = tfmax - tfmin;
  double npergrp = static_cast<double>(n) / 2.0;
  
  std::vector<int> id = seqcpp(1, n), stratum(n, 1);
  std::vector<int> treat(n), event(n), dropout(n), admin_censor(n);
  std::vector<int> pd(n), swtrt_latent(n), swtrt(n);
  std::vector<double> survivalTime_latent(n), survivalTime(n);
  std::vector<double> dropoutTime_latent(n), dropoutTime(n);
  std::vector<double> followupTime(n), censorTime(n), time(n);
  std::vector<double> pd_time_latent(n), pd_time(n);
  std::vector<double> swtrt_time_latent(n), swtrt_time(n);
  std::vector<double> rx(n), censor_time(n);
  
  std::vector<double> p_event1(nsim), p_event0(nsim); 
  std::vector<double> p_dropout1(nsim), p_dropout0(nsim); 
  std::vector<double> p_admin_censor1(nsim), p_admin_censor0(nsim);
  std::vector<double> p_pd0(nsim), p_swtrt0(nsim), p_recensored0(nsim);
  std::vector<double> psihat(nsim), loghrhat(nsim), seloghrcox(nsim);
  std::vector<double> hrhat(nsim), hrlowercox(nsim), hruppercox(nsim);
  std::vector<double> seloghrlr(nsim), hrlowerlr(nsim), hrupperlr(nsim);
  std::vector<double> seloghrboot(nsim), hrlowerboot(nsim), hrupperboot(nsim);
  
  double zcrit = boost_qnorm(1.0 - alpha / 2.0);
  double tcrit = boost_qt(1.0 - alpha / 2.0, n_boot - 1);
  
  // treatment indicator  
  for (int i = 0; i < n; ++i) {
    if (id[i] <= npergrp) treat[i] = 1; 
    else treat[i] = 0;
  }
  
  auto f = [n, low_psi, hi_psi, treat_modifier, recensor_type, 
            autoswitch, alpha, ties, tol](
                std::vector<int>& idb, 
                std::vector<int>& stratumb, 
                std::vector<double>& timeb, 
                std::vector<int>& eventb, 
                std::vector<int>& treatb, 
                std::vector<double>& rxb, 
                std::vector<double>& censor_timeb) -> ListCpp {
                  std::vector<double> init(1, NaN);
                  
                  // obtain the estimate of psi
                  double target = 0.0;
                  auto g = [&target, n, idb, stratumb, timeb, eventb, 
                            treatb, rxb, censor_timeb, treat_modifier, 
                            recensor_type, autoswitch](double x) -> double {
                              return f_est_psi_rpsftm(
                                x, n, idb, stratumb, timeb, eventb, 
                                treatb, rxb, censor_timeb, 
                                treat_modifier, recensor_type, 
                                autoswitch, target);
                            };
                  
                  double psihat = NaN, loghrhat = NaN, seloghrcox = NaN;
                  if (g(low_psi) > 0 && g(hi_psi) < 0) {
                    psihat = brent(g, low_psi, hi_psi, tol);
                    
                    // run Cox model to obtain the hazard ratio estimate
                    DataFrameCpp data_outcome = f_unswitched(
                      psihat*treat_modifier, n, idb, timeb, eventb, treatb,
                      rxb, censor_timeb, recensor_type, autoswitch);
                    
                    data_outcome.push_back(stratumb, "ustratum");
                    
                    ListCpp fit_outcome = phregcpp(
                      data_outcome, {"ustratum"}, "t_star", "", "d_star", 
                      {"treated"}, "", "", "", ties, init, 
                      0, 0, 0, 0, 0, alpha);
                    
                    DataFrameCpp parest = fit_outcome.get<DataFrameCpp>("parest");
                    loghrhat = parest.get<double>("beta")[0];
                    seloghrcox = parest.get<double>("sebeta")[0];
                  }
                  
                  ListCpp out;
                  out.push_back(psihat, "psihat");
                  out.push_back(loghrhat, "loghrhat");
                  out.push_back(seloghrcox, "seloghrcox");
                  return out;
                };
  
  
  for (int iter = 0; iter < nsim; ++iter) {
    // data generation
    for (int i = 0; i < n; ++i) {
      // latent survival time
      survivalTime_latent[i] = weib(rng) * exp(-psi*treat[i]);

      // latent dropout time
      dropoutTime_latent[i] = expo(rng) * exp(-omega*treat[i]);
      
      // administrative censoring time
      followupTime[i] = unif(rng) * twidth + tfmin;
      
      // latent disease progression time
      double x1 = g1(rng), x2 = g2(rng);
      double u = x1 / (x1 + x2);
      pd_time_latent[i] = u * survivalTime_latent[i];
      
      if (treat[i] == 0) {
        // latent switch indicator
        swtrt_latent[i] = (unif(rng) <= pswitch) ? 1 : 0;
        
        // latent switch time
        if (swtrt_latent[i] == 1) swtrt_time_latent[i] = pd_time_latent[i];
        else swtrt_time_latent[i] = NaN;
      } else {
        swtrt_latent[i] = 0;
        swtrt_time_latent[i] = NaN;
      }
      
      // survival time
      if (swtrt_latent[i] == 1) {
        survivalTime[i] = swtrt_time_latent[i] + 
          exp(-psi)*(survivalTime_latent[i] - swtrt_time_latent[i]);
      } else {
        survivalTime[i] = survivalTime_latent[i];
      }
      
      // dropout time
      if (swtrt_latent[i] == 1 && swtrt_time_latent[i] <= dropoutTime_latent[i]) {
        dropoutTime[i] =  swtrt_time_latent[i] + 
          exp(-omega)*(dropoutTime_latent[i] - swtrt_time_latent[i]);
      } else {
        dropoutTime[i] = dropoutTime_latent[i];
      }
      
      // censoring time
      censorTime[i] = std::min(dropoutTime[i], followupTime[i]);
      
      // observed survival time
      time[i] = std::min(survivalTime[i], censorTime[i]);
      
      // event indicator
      event[i] = (time[i] == survivalTime[i]) ? 1 : 0;
      
      // dropout indicator
      dropout[i] = (time[i] == dropoutTime[i]) ? 1 : 0;
      
      // administrative censoring indicator
      admin_censor[i] = (time[i] == followupTime[i]) ? 1 : 0;
      
      // observed disease progression time
      pd_time[i] = std::min(pd_time_latent[i], censorTime[i]);
      
      // disease progression indicator
      pd[i] = (pd_time[i] == pd_time_latent[i]) ? 1 : 0;
      
      if (treat[i] == 0) {
        // observed switch indicator
        if (swtrt_latent[i] == 1 && swtrt_time_latent[i] <= censorTime[i]) {
          swtrt[i] = 1;
        } else {
          swtrt[i] = 0;
        }
        
        // observed switch time
        if (swtrt[i] == 1) {
          swtrt_time[i] = swtrt_time_latent[i];  
        } else if (swtrt_latent[i] == 1) {
          swtrt_time[i] = censorTime[i];
        } else {
          swtrt_time[i] = NaN;
        }
      } else {
        swtrt[i] = 0;
        swtrt_time[i] = NaN;
      }
      
      // proportion of time on treatment
      if (treat[i] == 0 && swtrt[i] == 1) {
        rx[i] = (time[i] - swtrt_time[i]) / time[i];
      } else if (treat[i] == 0 && swtrt[i] == 0) {
        rx[i] = 0.0;
      } else {
        rx[i] = 1.0;
      }
      
      // whether to incorporate random censoring
      if (admin_recensor_only) {
        censor_time[i] = followupTime[i];
      } else {
        censor_time[i] = (event[i] == 1) ? followupTime[i] : time[i];
      }
    }
    
    double sum_event_treat = 0.0;   // sum(event * treat)
    double sum_event_notreat = 0.0; // sum(event * (1 - treat))
    double sum_dropout_treat = 0.0;   // sum(event * treat)
    double sum_dropout_notreat = 0.0; // sum(event * (1 - treat))
    double sum_admin_treat = 0.0;   // sum(event * treat)
    double sum_admin_notreat = 0.0; // sum(event * (1 - treat))
    double sum_pd_notreat = 0.0; // sum(pd * (1 - treat))
    double sum_swtrt_notreat = 0.0; // sum(swtrt * (1 - treat))
    for (int i = 0; i < n; ++i) {
      if (treat[i] == 1) {
        sum_dropout_treat   += dropout[i];
        sum_admin_treat     += admin_censor[i];
      } else {
        sum_dropout_notreat += dropout[i];
        sum_admin_notreat   += admin_censor[i];
        sum_pd_notreat      += pd[i];
        sum_swtrt_notreat   += swtrt[i];
      }
    }
    
    p_event1[iter] = sum_event_treat / npergrp;
    p_event0[iter] = sum_event_notreat / npergrp;
    p_dropout1[iter] = sum_dropout_treat / npergrp;
    p_dropout0[iter] = sum_dropout_notreat / npergrp;
    p_admin_censor1[iter] = sum_admin_treat / npergrp;
    p_admin_censor0[iter] = sum_admin_notreat / npergrp;
    p_pd0[iter] = sum_pd_notreat / npergrp;
    p_swtrt0[iter] = sum_swtrt_notreat / npergrp;
    
    DataFrameCpp dt;
    dt.push_back(stratum, "stratum");
    dt.push_back(treat, "treat");
    dt.push_back(time, "time");
    dt.push_back(event, "event");
    
    DataFrameCpp lr = lrtestcpp(dt, {"stratum"}, "treat", "time", "", "event");
    double logRankZ = lr.get<double>("logRankZ")[0];
    
    ListCpp out = f(id, stratum, time, event, treat, rx, censor_time);
    
    psihat[iter] = out.get<double>("psihat");
    loghrhat[iter] = out.get<double>("loghrhat");
    seloghrcox[iter] = out.get<double>("seloghrcox");
    hrhat[iter] = std::exp(loghrhat[iter]);
    hrlowercox[iter] = std::exp(loghrhat[iter] - zcrit * seloghrcox[iter]);
    hruppercox[iter] = std::exp(loghrhat[iter] + zcrit * seloghrcox[iter]);
    
    seloghrlr[iter] = loghrhat[iter] / logRankZ;
    if (seloghrlr[iter] <= 0) {
      std::string str2 = "invalid standard error estimate for iter ";
      str2 += std::to_string(iter) + "\n";
      std::string errmsg = str1 + str2;
      thread_utils::push_thread_warning(errmsg);
    }
    
    hrlowerlr[iter] = std::exp(loghrhat[iter] - zcrit * seloghrlr[iter]);
    hrupperlr[iter] = std::exp(loghrhat[iter] + zcrit * seloghrlr[iter]);

    // bootstrap
    std::vector<int> idb(n), stratumb(n), treatb(n), eventb(n);
    std::vector<double> timeb(n), rxb(n), censor_timeb(n);

    // sort data by treatment group
    std::vector<int> idx0, idx1;
    idx0.reserve(n);
    idx1.reserve(n);
    for (int i = 0; i < n; ++i) {
      if (treat[i] == 0) {
        idx0.push_back(i);
      } else {
        idx1.push_back(i);
      }
    }  
    int n0 = static_cast<int>(idx0.size());
    int n1 = static_cast<int>(idx1.size());
    
    std::uniform_int_distribution<int> index_dist0(0, n0 - 1);
    std::uniform_int_distribution<int> index_dist1(0, n1 - 1);
    
    std::vector<int> order(n);
    for (int i = 0; i < n0; i++) {
      order[i] = idx0[i];
    }
    for (int i = 0; i < n1; i++){
      order[n0 + i] = idx1[i];
    }
    
    subset_in_place(id, order);
    subset_in_place(stratum, order);
    subset_in_place(time, order);
    subset_in_place(event, order);
    subset_in_place(treat, order);
    subset_in_place(rx, order);
    subset_in_place(censor_time, order);
    
    std::vector<double> loghrhats(n_boot);
    for (int k = 0; k < n_boot; ++k) {
      // sample the data with replacement by treatment group
      for (int i = 0; i < n; ++i) {
        int j = (treat[i] == 0) ? index_dist0(rng) : n0 + index_dist1(rng);
        idb[i] = id[j];
        stratumb[i] = stratum[j];
        timeb[i] = time[j];
        eventb[i] = event[j];
        treatb[i] = treat[j];
        rxb[i] = rx[j];
        censor_timeb[i] = censor_time[j];
      }
      
      ListCpp out = f(idb, stratumb, timeb, eventb, treatb, rxb, censor_timeb);
      loghrhats[k] = out.get<double>("loghrhat");
    }
    
    if (std::any_of(loghrhats.begin(), loghrhats.end(), 
                    [](bool x){ return std::isnan(x); }))  {
      std::string str2 = "invalid log HR estimate using boot for iter ";
      str2 += std::to_string(iter) + "\n" + 
        "bootstrap SD based on valid log HR estimates only" + "\n";
      std::string errmsg = str1 + str2;
      thread_utils::push_thread_warning(errmsg);
      
      std::vector<double> loghrhats1;
      loghrhats1.reserve(n_boot);
      std::size_t count = 0;
      for (int k = 0; k < n_boot; ++k) {
        if (!std::isnan(loghrhats[k])) {
          loghrhats1.push_back(loghrhats[k]);
          count++;
        }
      }
      double meanloghr1, sdloghr1;
      mean_sd(loghrhats1.data(), count, meanloghr1, sdloghr1);
      seloghrboot[iter] = sdloghr1;
    } else {
      double meanloghr, sdloghr;
      mean_sd(loghrhats.data(), n_boot, meanloghr, sdloghr);
      seloghrboot[iter] = sdloghr;
    }

    hrlowerboot[iter] = std::exp(loghrhat[iter] - tcrit * seloghrboot[iter]);
    hrupperboot[iter] = std::exp(loghrhat[iter] + tcrit * seloghrboot[iter]);
    
    // proportion of recensored events among control patients
    DataFrameCpp Sstar = f_untreated(
      psihat[iter] * treat_modifier, n, id, time, event, 
      treat, rx, censor_time, recensor_type, autoswitch);
    
    std::vector<int> d_star = Sstar.get<int>("d_star");
    std::vector<int> recensored(n);
    for (int i = 0; i < n; ++i) {
      recensored[i] = event[i] * (1 - d_star[i]);
    }
    double sum_recens_notreat = 0.0;
    for (int i = 0; i < n; ++i) {
      if (treat[i] == 0) sum_recens_notreat += recensored[i];
    }
    p_recensored0[iter] = sum_recens_notreat / npergrp;
  }
  
  double p_event_1 = mean_kahan(p_event1);
  double p_event_0 = mean_kahan(p_event0);
  double p_dropout_1 = mean_kahan(p_dropout1);
  double p_dropout_0 = mean_kahan(p_dropout0);
  double p_admin_censor_1 = mean_kahan(p_admin_censor1);
  double p_admin_censor_0 = mean_kahan(p_admin_censor0);
  double p_pd_0 = mean_kahan(p_pd0);
  double p_swtrt_0 = mean_kahan(p_swtrt0);
  double p_recensored_0 = mean_kahan(p_recensored0);
  
  double psi_est, psi_se;
  mean_sd(psihat.data(), nsim, psi_est, psi_se);
  double psi_bias = psi_est - psi;
  std::vector<double> psi_diff(nsim), psi_diff_squared(nsim);
  for (int i = 0; i < nsim; ++i) {
    psi_diff[i] = psihat[i] - psi;
    psi_diff_squared[i] = psi_diff[i] * psi_diff[i];
  }
  double psi_mse = mean_kahan(psi_diff_squared);
  
  double loghr = shape * psi;
  double loghr_est, loghr_se;
  mean_sd(loghrhat.data(), nsim, loghr_est, loghr_se);
  double loghr_bias = loghr_est - loghr;
  std::vector<double> loghr_diff(nsim), loghr_diff_squared(nsim);
  for (int i = 0; i < nsim; ++i) {
    loghr_diff[i] = loghrhat[i] - loghr;
    loghr_diff_squared[i] = loghr_diff[i] * loghr_diff[i];
  }
  double loghr_mse = mean_kahan(loghr_diff_squared);
  
  double loghr_se_cox = mean_kahan(seloghrcox);
  double loghr_se_lr = mean_kahan(seloghrlr);
  double loghr_se_boot = mean_kahan(seloghrboot);
  double hr = std::exp(loghr);
  double hr_est = std::exp(loghr_est);
  double hr_pctbias = 100.0 * (hr_est - hr) / hr;
  
  double hr_ci_cover_cox = 0.0, hr_ci_cover_lr = 0.0, hr_ci_cover_boot = 0.0;
  for (int i = 0; i < nsim; ++i) {
    hr_ci_cover_cox += (hrlowercox[i] < hr && hruppercox[i] > hr) ? 1.0 : 0.0;
    hr_ci_cover_lr += (hrlowerlr[i] < hr && hrupperlr[i] > hr) ? 1.0 : 0.0;
    hr_ci_cover_boot += (hrlowerboot[i] < hr && hrupperboot[i] > hr) ? 1.0 : 0.0;
  }
  hr_ci_cover_cox /= static_cast<double>(nsim);
  hr_ci_cover_lr /= static_cast<double>(nsim);
  hr_ci_cover_boot /= static_cast<double>(nsim);
  
  DataFrameCpp result;
  result.push_back(recensor_type, "recensor_type");
  result.push_back(admin_recensor_only, "admin_recensor_only");
  result.push_back(p_event_1, "p_event_1");
  result.push_back(p_dropout_1, "p_dropout_1");
  result.push_back(p_admin_censor_1, "p_admin_censor_1");
  result.push_back(p_event_0, "p_event_0");
  result.push_back(p_dropout_0, "p_dropout_0");
  result.push_back(p_admin_censor_0, "p_admin_censor_0");
  result.push_back(p_pd_0, "p_pd_0");
  result.push_back(p_swtrt_0, "p_swtrt_0");
  result.push_back(p_recensored_0, "p_recensored_0");
  result.push_back(psi, "psi");
  result.push_back(psi_est, "psi_est");
  result.push_back(psi_bias, "psi_bias");
  result.push_back(psi_se, "psi_se");
  result.push_back(psi_mse, "psi_mse");
  result.push_back(loghr, "loghr");
  result.push_back(loghr_est, "loghr_est");
  result.push_back(loghr_bias, "loghr_bias");
  result.push_back(loghr_se, "loghr_se");
  result.push_back(loghr_mse, "loghr_mse");
  result.push_back(hr, "hr");
  result.push_back(hr_est, "hr_est");
  result.push_back(hr_pctbias, "hr_pctbias");
  result.push_back(loghr_se_cox, "loghr_se_cox");
  result.push_back(loghr_se_lr, "loghr_se_lr");
  result.push_back(loghr_se_boot, "loghr_se_boot");
  result.push_back(hr_ci_cover_cox, "hr_ci_cover_cox");
  result.push_back(hr_ci_cover_lr, "hr_ci_cover_lr");
  result.push_back(hr_ci_cover_boot, "hr_ci_cover_boot");
  
  thread_utils::drain_thread_warnings_to_R();
  return Rcpp::wrap(result);
}
