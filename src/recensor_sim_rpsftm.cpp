#include "utilities.h"
#include "survival_analysis.h"

using namespace Rcpp;


// counterfactual untreated survival times and event indicators
DataFrame f_untreated(
    const double psi,
    const IntegerVector& id,
    const NumericVector& time,
    const IntegerVector& event,
    const IntegerVector& treat,
    const NumericVector& rx,
    const NumericVector& censor_time,
    const int recensor_type,
    const bool autoswitch) {
  
  IntegerVector nonswitchers = which(((treat == 1) & (rx == 1.0)) | 
    ((treat == 0) & (rx == 0.0)));
  
  double a = exp(psi);
  NumericVector u_star = time*((1 - rx) + rx*a);
  NumericVector t_star = clone(u_star);
  IntegerVector d_star = clone(event);
  
  if (recensor_type == 1) { // recensor all patients
    NumericVector c_star = censor_time*std::min(1.0, a);
    
    if (autoswitch) {
      NumericVector rx1 = rx[treat == 1];
      NumericVector rx0 = rx[treat == 0];
      if (is_true(all(rx1 == 1.0))) c_star[treat == 1] = R_PosInf;
      if (is_true(all(rx0 == 0.0))) c_star[treat == 0] = R_PosInf;
    }
    
    t_star = pmin(u_star, c_star);
    d_star[c_star < u_star] = 0;
  } else if (recensor_type == 2) { // recensor only switchers
    NumericVector c_star = censor_time*std::min(1.0, a);
    
    if (autoswitch) {
      NumericVector rx1 = rx[treat == 1];
      NumericVector rx0 = rx[treat == 0];
      if (is_true(all(rx1 == 1.0))) c_star[treat == 1] = R_PosInf;
      if (is_true(all(rx0 == 0.0))) c_star[treat == 0] = R_PosInf;
    }
    
    c_star[nonswitchers] = R_PosInf;
    
    t_star = pmin(u_star, c_star);
    d_star[c_star < u_star] = 0;
  } else if (recensor_type == 3) { // recensor only switchers with 
    // projected latent event time beyond the end of study
    NumericVector c_star = censor_time;
    
    if (autoswitch) {
      NumericVector rx1 = rx[treat == 1];
      NumericVector rx0 = rx[treat == 0];
      if (is_true(all(rx1 == 1.0))) c_star[treat == 1] = R_PosInf;
      if (is_true(all(rx0 == 0.0))) c_star[treat == 0] = R_PosInf;
    }
    
    c_star[nonswitchers] = R_PosInf;
    
    t_star = pmin(u_star, c_star);
    d_star[c_star < u_star] = 0;
  }
  
  DataFrame result = DataFrame::create(
    Named("uid") = id,
    Named("t_star") = t_star,
    Named("d_star") = d_star,
    Named("treated") = treat
  );
  
  return result;
}


// counterfactual unswitched survival times and event indicators
DataFrame f_unswitched(
    const double psi,
    const int n,
    const IntegerVector& id,
    const NumericVector& time,
    const IntegerVector& event,
    const IntegerVector& treat,
    const NumericVector& rx,
    const NumericVector& censor_time,
    const int recensor_type,
    const bool autoswitch) {
  
  IntegerVector nonswitchers = which(((treat == 1) & (rx == 1.0)) | 
    ((treat == 0) & (rx == 0.0)));
  
  double a = exp(psi), a1 = exp(-psi);
  NumericVector u_star(n), t_star(n);
  IntegerVector d_star(n);
  for (int i=0; i<n; i++) {
    if (treat[i] == 0) {
      u_star[i] = time[i]*((1 - rx[i]) + rx[i]*a);
    } else {
      u_star[i] = time[i]*(rx[i] + (1 - rx[i])*a1);
    }
    t_star[i] = u_star[i];
    d_star[i] = event[i];
  }
  
  if (recensor_type == 1) { // recensor all patients
    NumericVector c_star(n);
    double c0 = std::min(1.0, a), c1 = std::min(1.0, a1);
    for (int i=0; i<n; i++) {
      if (treat[i] == 0) {
        c_star[i] = censor_time[i]*c0;
      } else {
        c_star[i] = censor_time[i]*c1;
      }
    }
    
    if (autoswitch) {
      NumericVector rx1 = rx[treat == 1];
      NumericVector rx0 = rx[treat == 0];
      if (is_true(all(rx1 == 1.0))) c_star[treat == 1] = R_PosInf;
      if (is_true(all(rx0 == 0.0))) c_star[treat == 0] = R_PosInf;
    }
    
    t_star = pmin(u_star, c_star);
    d_star[c_star < u_star] = 0;
  } else if (recensor_type == 2) { // recensor only switchers
    NumericVector c_star(n);
    double c0 = std::min(1.0, a), c1 = std::min(1.0, a1);
    for (int i=0; i<n; i++) {
      if (treat[i] == 0) {
        c_star[i] = censor_time[i]*c0;
      } else {
        c_star[i] = censor_time[i]*c1;
      }
    }
    
    if (autoswitch) {
      NumericVector rx1 = rx[treat == 1];
      NumericVector rx0 = rx[treat == 0];
      if (is_true(all(rx1 == 1.0))) c_star[treat == 1] = R_PosInf;
      if (is_true(all(rx0 == 0.0))) c_star[treat == 0] = R_PosInf;
    }
    
    c_star[nonswitchers] = R_PosInf;
    
    t_star = pmin(u_star, c_star);
    d_star[c_star < u_star] = 0;
  } else if (recensor_type == 3) { // recensor only switchers with 
    // projected latent event time beyond the end of study
    NumericVector c_star = censor_time;
    
    if (autoswitch) {
      NumericVector rx1 = rx[treat == 1];
      NumericVector rx0 = rx[treat == 0];
      if (is_true(all(rx1 == 1.0))) c_star[treat == 1] = R_PosInf;
      if (is_true(all(rx0 == 0.0))) c_star[treat == 0] = R_PosInf;
    }
    
    c_star[nonswitchers] = R_PosInf;
    
    t_star = pmin(u_star, c_star);
    d_star[c_star < u_star] = 0;
  }
  
  DataFrame result = DataFrame::create(
    Named("uid") = id,
    Named("t_star") = t_star,
    Named("d_star") = d_star,
    Named("treated") = treat
  );
  
  return result;
}


double f_est_psi_rpsftm(
    const double psi,
    const IntegerVector& id,
    const IntegerVector& stratum,
    const NumericVector& time,
    const IntegerVector& event,
    const IntegerVector& treat,
    const NumericVector& rx,
    const NumericVector& censor_time,
    const double treat_modifier,
    const int recensor_type,
    const bool autoswitch,
    const double target) {
  
  DataFrame Sstar = f_untreated(psi*treat_modifier, id, time, event, treat, 
                                rx, censor_time, recensor_type, autoswitch);
  
  Sstar.push_back(stratum, "ustratum");
  DataFrame df = lrtest(Sstar, "", "ustratum", "treated", "t_star", "", 
                        "d_star", "", 0, 0);
  
  double z = df["logRankZ"];
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
//' @param seed Optional. Random seed for reproducibility. If not provided, 
//'   the global seed is used.
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
DataFrame recensor_sim_rpsftm(const int nsim = NA_INTEGER,
                              const int n = NA_INTEGER,
                              const double shape = NA_REAL,
                              const double scale = NA_REAL,
                              const double gamma = NA_REAL,
                              const double tfmin = NA_REAL,
                              const double tfmax = NA_REAL,
                              const double psi = NA_REAL,
                              const double omega = NA_REAL,
                              const double pswitch = NA_REAL,
                              const double a = NA_REAL,
                              const double b = NA_REAL,
                              const double low_psi = -1,
                              const double hi_psi = 1,
                              const double treat_modifier = 1,
                              const int recensor_type = 1,
                              const bool admin_recensor_only = true,
                              const bool autoswitch = true,
                              const double alpha = 0.05,
                              const std::string ties = "efron",
                              const double tol = 1.0e-6,
                              const bool boot = true,
                              const int n_boot = 1000,
                              const int seed = NA_INTEGER) {
  
  std::string str1 = "recensor_type = ";
  str1 += std::to_string(recensor_type) + ", admin_recensor_only = " + 
    std::to_string(admin_recensor_only) + "\n";
  
  if (seed != NA_INTEGER) set_seed(seed);
  
  IntegerVector id = seq_len(n), stratum(n, 1);
  IntegerVector treat(n), event(n), dropout(n), admin_censor(n);
  IntegerVector pd(n), swtrt_latent(n), swtrt(n);
  NumericVector survivalTime_latent(n), survivalTime(n);
  NumericVector dropoutTime_latent(n), dropoutTime(n);
  NumericVector followupTime(n), censorTime(n), time(n);
  NumericVector pd_time_latent(n), pd_time(n);
  NumericVector swtrt_time_latent(n), swtrt_time(n);
  NumericVector rx(n), censor_time(n);
  
  NumericVector p_event1(nsim), p_event0(nsim); 
  NumericVector p_dropout1(nsim), p_dropout0(nsim); 
  NumericVector p_admin_censor1(nsim), p_admin_censor0(nsim);
  NumericVector p_pd0(nsim), p_swtrt0(nsim), p_recensored0(nsim);
  NumericVector psihat(nsim), loghrhat(nsim), seloghrcox(nsim);
  NumericVector hrhat(nsim), hrlowercox(nsim), hruppercox(nsim);
  NumericVector seloghrlr(nsim), hrlowerlr(nsim), hrupperlr(nsim);
  NumericVector seloghrboot(nsim), hrlowerboot(nsim), hrupperboot(nsim);
  
  double zcrit = R::qnorm(1-alpha/2, 0, 1, 1, 0);
  double tcrit = R::qt(1-alpha/2, n_boot-1, 1, 0);
  
  // treatment indicator  
  treat[id <= n/2] = 1;
  for (int iter=0; iter<nsim; iter++) {
    // data generation
    for (int i=0; i<n; i++) {
      // latent survival time
      survivalTime_latent[i] = R::rweibull(shape, scale*exp(-psi*treat[i]));

      // latent dropout time
      dropoutTime_latent[i] = R::rweibull(1, 1/gamma*exp(-omega*treat[i]));
      
      // administrative censoring time
      followupTime[i] = R::runif(tfmin, tfmax);
      
      // latent disease progression time
      double u = R::rbeta(a, b);
      pd_time_latent[i] = u*survivalTime_latent[i];
      
      if (treat[i] == 0) {
        // latent switch indicator
        u = R::runif(0,1);
        swtrt_latent[i] = 1*(u <= pswitch);
        
        // latent switch time
        if (swtrt_latent[i] == 1) {
          swtrt_time_latent[i] = pd_time_latent[i];
        } else {
          swtrt_time_latent[i] = NA_REAL;
        }
      } else {
        swtrt_latent[i] = 0;
        swtrt_time_latent[i] = NA_REAL;
      }
      
      // survival time
      if (swtrt_latent[i] == 1) {
        survivalTime[i] = swtrt_time_latent[i] + 
          exp(-psi)*(survivalTime_latent[i] - swtrt_time_latent[i]);
      } else {
        survivalTime[i] = survivalTime_latent[i];
      }
      
      // dropout time
      if (swtrt_latent[i] == 1 && 
          swtrt_time_latent[i] <= dropoutTime_latent[i]) {
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
      event[i] = 1*(time[i] == survivalTime[i]);
      
      // dropout indicator
      dropout[i] = 1*(time[i] == dropoutTime[i]);
      
      // administrative censoring indicator
      admin_censor[i] = 1*(time[i] == followupTime[i]);
      
      // observed disease progression time
      pd_time[i] = std::min(pd_time_latent[i], censorTime[i]);
      
      // disease progression indicator
      pd[i] = 1*(pd_time[i] == pd_time_latent[i]);
      
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
          swtrt_time[i] = NA_REAL;
        }
      } else {
        swtrt[i] = 0;
        swtrt_time[i] = NA_REAL;
      }
      
      // proportion of time on treatment
      if (treat[i] == 0 && swtrt[i] == 1) {
        rx[i] = (time[i] - swtrt_time[i])/time[i];
      } else if (treat[i] == 0 && swtrt[i] == 0) {
        rx[i] = 0.0;
      } else {
        rx[i] = 1.0;
      }
      
      // whether to incorporate random censoring
      if (admin_recensor_only) {
        censor_time[i] = followupTime[i];
      } else {
        censor_time[i] = event[i] == 1 ? followupTime[i] : time[i];
      }
    }
    
    p_event1[iter] = sum(event*treat)/(n/2.0);
    p_event0[iter] = sum(event*(1-treat))/(n/2.0);
    p_dropout1[iter] = sum(dropout*treat)/(n/2.0);
    p_dropout0[iter] = sum(dropout*(1-treat))/(n/2.0);
    p_admin_censor1[iter] = sum(admin_censor*treat)/(n/2.0);
    p_admin_censor0[iter] = sum(admin_censor*(1-treat))/(n/2.0);
    p_pd0[iter] = sum(pd*(1-treat))/(n/2.0);
    p_swtrt0[iter] = sum(swtrt*(1-treat))/(n/2.0);
    
    DataFrame dt = DataFrame::create(
      Named("stratum") = stratum,
      Named("treat") = treat,
      Named("time") = time,
      Named("event") = event
    );
    
    DataFrame lr = lrtest(dt,"","stratum","treat","time","","event","",0,0);
    double logRankZ = lr["logRankZ"];
    
    auto f = [n, low_psi, hi_psi, treat_modifier, recensor_type, 
              autoswitch, alpha, ties, tol](
                  IntegerVector& idb, IntegerVector& stratumb, 
                  NumericVector& timeb, IntegerVector& eventb, 
                  IntegerVector& treatb, NumericVector& rxb, 
                  NumericVector& censor_timeb)->List {
                    NumericVector init(1, NA_REAL);
                    
                    // obtain the estimate of psi
                    double target = 0;
                    auto g = [&target, idb, stratumb, timeb, eventb, 
                              treatb, rxb, censor_timeb, treat_modifier, 
                              recensor_type, autoswitch](double x)->double {
                                return f_est_psi_rpsftm(
                                  x, idb, stratumb, timeb, eventb, 
                                  treatb, rxb, censor_timeb, 
                                  treat_modifier, recensor_type, 
                                  autoswitch, target);
                              };
                    
                    double psihat = NA_REAL;
                    double loghrhat = NA_REAL, seloghrcox = NA_REAL;
                    if (g(low_psi) > 0 && g(hi_psi) < 0) {
                      psihat = brent(g, low_psi, hi_psi, tol);
                      
                      // run Cox model to obtain the hazard ratio estimate
                      DataFrame data_outcome = f_unswitched(
                        psihat*treat_modifier, n, idb, timeb, eventb, treatb,
                        rxb, censor_timeb, recensor_type, autoswitch);
                      
                      data_outcome.push_back(stratumb, "ustratum");
                      
                      List fit_outcome = phregcpp(
                        data_outcome, "", "ustratum", "t_star", "", "d_star", 
                        "treated", "", "", "", ties, init, 
                        0, 0, 0, 0, 0, alpha, 50, 1.0e-9);
                      
                      DataFrame parest = DataFrame(fit_outcome["parest"]);
                      NumericVector beta = parest["beta"];
                      NumericVector sebeta = parest["sebeta"];
                      loghrhat = beta[0];
                      seloghrcox = sebeta[0];
                    }
                    
                    List out = List::create(
                      Named("psihat") = psihat,
                      Named("loghrhat") = loghrhat,
                      Named("seloghrcox") = seloghrcox);
                    
                    return out;
                  };
    
    List out = f(id, stratum, time, event, treat, rx, censor_time);
    
    psihat[iter] = out["psihat"];
    loghrhat[iter] = out["loghrhat"];
    seloghrcox[iter] = out["seloghrcox"];
    hrhat[iter] = exp(loghrhat[iter]);
    hrlowercox[iter] = exp(loghrhat[iter] - zcrit*seloghrcox[iter]);
    hruppercox[iter] = exp(loghrhat[iter] + zcrit*seloghrcox[iter]);
    
    seloghrlr[iter] = loghrhat[iter]/logRankZ;
    if (seloghrlr[iter] <= 0) {
      std::string str2 = "invalid standard error estimate for iter ";
      str2 += std::to_string(iter) + "\n";
      std::string errmsg = str1 + str2;
      Rcout << errmsg << "\n";
    }
    
    hrlowerlr[iter] = exp(loghrhat[iter] - zcrit*seloghrlr[iter]);
    hrupperlr[iter] = exp(loghrhat[iter] + zcrit*seloghrlr[iter]);

    // bootstrap
    IntegerVector idb(n), stratumb(n), treatb(n), eventb(n);
    NumericVector timeb(n), rxb(n), censor_timeb(n);

    // sort data by treatment group
    IntegerVector idx0 = which(treat == 0);
    IntegerVector idx1 = which(treat == 1);
    int n0 = static_cast<int>(idx0.size());
    int n1 = static_cast<int>(idx1.size());
    IntegerVector order(n);
    for (int i=0; i<n0; i++) {
      order[i] = idx0[i];
    }
    for (int i=0; i<n1; i++){
      order[n0+i] = idx1[i];
    }
    
    id = id[order];
    stratum = stratum[order];
    time = time[order];
    event = event[order];
    treat = treat[order];
    rx = rx[order];
    censor_time = censor_time[order];
    
    NumericVector loghrhats(n_boot);
    for (int k=0; k<n_boot; k++) {
      // sample the data with replacement by treatment group
      for (int i=0; i<n; i++) {
        double u = R::runif(0,1);
        int j = treat[i] == 0 ? static_cast<int>(std::floor(u*n0)) 
          : n0 + static_cast<int>(std::floor(u*n1));
        
        idb[i] = id[j];
        stratumb[i] = stratum[j];
        timeb[i] = time[j];
        eventb[i] = event[j];
        treatb[i] = treat[j];
        rxb[i] = rx[j];
        censor_timeb[i] = censor_time[j];
      }
      
      List out = f(idb, stratumb, timeb, eventb, treatb, rxb, censor_timeb);
      
      loghrhats[k] = out["loghrhat"];
    }
    
    if (is_true(any(is_na(loghrhats)))) {
      std::string str2 = "invalid log HR estimate using boot for iter ";
      str2 += std::to_string(iter) + "\n" + 
        "bootstrap SD based on valid log HR estimates only" + "\n";
      std::string errmsg = str1 + str2;
      Rcout << errmsg << "\n";
      LogicalVector index1 = !is_na(loghrhats);
      NumericVector loghrhats1 = loghrhats[index1];
      seloghrboot[iter] = sd(loghrhats1);
    } else {
      seloghrboot[iter] = sd(loghrhats);
    }

    hrlowerboot[iter] = exp(loghrhat[iter] - tcrit*seloghrboot[iter]);
    hrupperboot[iter] = exp(loghrhat[iter] + tcrit*seloghrboot[iter]);
    
    // proportion of recensored events among control patients
    DataFrame Sstar = f_untreated(psihat[iter]*treat_modifier, id, time, 
                                  event, treat, rx, censor_time, 
                                  recensor_type, autoswitch);
    
    IntegerVector d_star = Sstar["d_star"];
    IntegerVector recensored(n);
    for (int i=0; i<n; i++) {
      recensored[i] = event[i] * (1 - d_star[i]);
    }
    p_recensored0[iter] = sum(recensored*(1-treat))/(n/2.0);
  }
  
  double p_event_1 = mean(p_event1);
  double p_event_0 = mean(p_event0);
  double p_dropout_1 = mean(p_dropout1);
  double p_dropout_0 = mean(p_dropout0);
  double p_admin_censor_1 = mean(p_admin_censor1);
  double p_admin_censor_0 = mean(p_admin_censor0);
  double p_pd_0 = mean(p_pd0);
  double p_swtrt_0 = mean(p_swtrt0);
  double p_recensored_0 = mean(p_recensored0);
  double psi_est = mean(psihat);
  double psi_bias = psi_est - psi;
  double psi_se = sd(psihat);
  NumericVector psi_diff = psihat - psi;
  double psi_mse = mean(psi_diff*psi_diff);
  
  double loghr = shape*psi;
  double loghr_est = mean(loghrhat);
  double loghr_bias = loghr_est - loghr;
  double loghr_se = sd(loghrhat);
  NumericVector loghr_diff = loghrhat - loghr;
  double loghr_mse = mean(loghr_diff*loghr_diff);
  
  double loghr_se_cox = mean(seloghrcox);
  double loghr_se_lr = mean(seloghrlr);
  double loghr_se_boot = mean(seloghrboot);
  double hr = exp(loghr);
  double hr_est = exp(loghr_est);
  double hr_pctbias = 100*(hr_est - hr)/hr;
  double hr_ci_cover_cox = mean((hrlowercox < hr)*(hruppercox > hr));
  double hr_ci_cover_lr = mean((hrlowerlr < hr)*(hrupperlr > hr));
  double hr_ci_cover_boot = mean((hrlowerboot < hr)*(hrupperboot > hr));
  
  DataFrame result = DataFrame::create(
    Named("recensor_type") = recensor_type,
    Named("admin_recensor_only") = admin_recensor_only,
    Named("p_event_1") = p_event_1,
    Named("p_dropout_1") = p_dropout_1,
    Named("p_admin_censor_1") = p_admin_censor_1,
    Named("p_event_0") = p_event_0,
    Named("p_dropout_0") = p_dropout_0,
    Named("p_admin_censor_0") = p_admin_censor_0,
    Named("p_pd_0") = p_pd_0,
    Named("p_swtrt_0") = p_swtrt_0,
    Named("p_recensored_0") = p_recensored_0,
    Named("psi") = psi,
    Named("psi_est") = psi_est,
    Named("psi_bias") = psi_bias,
    Named("psi_se") = psi_se,
    Named("psi_mse") = psi_mse,
    Named("loghr") = loghr,
    Named("loghr_est") = loghr_est,
    Named("loghr_bias") = loghr_bias,
    Named("loghr_se") = loghr_se,
    Named("loghr_mse") = loghr_mse,
    Named("hr") = hr,
    Named("hr_est") = hr_est,
    Named("hr_pctbias") = hr_pctbias,
    Named("loghr_se_cox") = loghr_se_cox,
    Named("loghr_se_lr") = loghr_se_lr,
    Named("loghr_se_boot") = loghr_se_boot,
    Named("hr_ci_cover_cox") = hr_ci_cover_cox,
    Named("hr_ci_cover_lr") = hr_ci_cover_lr,
    Named("hr_ci_cover_boot") = hr_ci_cover_boot
  );
  
  return result;
}
