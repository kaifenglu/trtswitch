#include "utilities.h"
#include "survival_analysis.h"

using namespace Rcpp;


// [[Rcpp::export]]
List tsesimpcpp1(const DataFrame data,
                 const StringVector& stratum = "",
                 const std::string time = "time",
                 const std::string event = "event",
                 const std::string treat = "treat",
                 const std::string censor_time = "censor_time",
                 const std::string pd = "pd",
                 const std::string pd_time = "pd_time",
                 const std::string swtrt = "swtrt",
                 const std::string swtrt_time = "swtrt_time",
                 const StringVector& base_cov = "",
                 const StringVector& base2_cov = "",
                 const std::string aft_dist = "weibull",
                 const bool strata_main_effect_only = 1,
                 const bool recensor = 1,
                 const bool admin_recensor_only = 1,
                 const bool swtrt_control_only = 1,
                 const double alpha = 0.05,
                 const std::string ties = "efron",
                 const double offset = 1,
                 const bool boot = 1,
                 const int n_boot = 1000,
                 const int seed = NA_INTEGER) {
  
  int i, j, k, n = data.nrow();
  
  int p = static_cast<int>(base_cov.size());
  if (p == 1 && (base_cov[0] == "" || base_cov[0] == "none")) p = 0;
  
  int p2 = static_cast<int>(base2_cov.size());
  if (p2 == 1 && (base2_cov[0] == "" || base2_cov[0] == "none")) p2 = 0;
  
  int p_stratum = static_cast<int>(stratum.size());
  
  IntegerVector stratumn(n);
  IntegerVector d(p_stratum);
  IntegerMatrix stratan(n,p_stratum);
  if (p_stratum == 1 && (stratum[0] == "" || stratum[0] == "none")) {
    stratumn.fill(1);
    d[0] = 1;
    stratan(_,0) = stratumn;
  } else {
    List out = bygroup(data, stratum);
    stratumn = out["index"];
    d = out["nlevels"];
    stratan = as<IntegerMatrix>(out["indices"]);
  }
  
  IntegerVector stratumn_unique = unique(stratumn);
  int nstrata = static_cast<int>(stratumn_unique.size());
  
  NumericVector timenz = data[time];
  NumericVector timen = clone(timenz);
  
  IntegerVector eventnz = data[event];
  IntegerVector eventn = clone(eventnz);

  IntegerVector treatnz = data[treat];
  IntegerVector treatn = clone(treatnz);
  
  NumericVector censor_timenz = data[censor_time];
  NumericVector censor_timen = clone(censor_timenz);
  
  if (!admin_recensor_only) {
    for (i=0; i<n; i++) {
      if (eventn[i] == 0) { // use the actual censoring time for dropouts
        censor_timen[i] = timen[i];
      }
    }
  }
  
  IntegerVector pdnz = data[pd];
  IntegerVector pdn = clone(pdnz);
  
  NumericVector pd_timenz = data[pd_time];
  NumericVector pd_timen = clone(pd_timenz);
  
  IntegerVector swtrtnz = data[swtrt];
  IntegerVector swtrtn = clone(swtrtnz);
  
  NumericVector swtrt_timenz = data[swtrt_time];
  NumericVector swtrt_timen = clone(swtrt_timenz);
  
  // covariates for the Cox model containing treat and base_cov
  StringVector covariates(p+1);
  NumericMatrix zn(n,p);
  covariates[0] = "treated";
  for (j=0; j<p; j++) {
    String zj = base_cov[j];
    NumericVector u = data[zj];
    covariates[j+1] = zj;
    zn(_,j) = u;
  }
  
  // covariates for the accelerated failure time model for control with pd
  // including stratum and base2_cov
  int q; // number of columns corresponding to the strata effects
  if (strata_main_effect_only) {
    q = sum(d - 1);
  } else {
    q = nstrata - 1;
  }
  
  StringVector covariates_aft(q+p2+1);
  NumericMatrix zn_aft(n,q+p2);
  covariates_aft[0] = "swtrt";
  if (strata_main_effect_only) {
    k = 0;
    for (i=0; i<p_stratum; i++) {
      for (j=0; j<d[i]-1; j++) {
        covariates_aft[k+j+1] = "stratum_" + std::to_string(i+1) +
          "_level_" + std::to_string(j+1);
        zn_aft(_,k+j) = 1.0*(stratan(_,i) == j+1);
      }
      k += d[i]-1;
    }
  } else {
    for (j=0; j<nstrata-1; j++) {
      covariates_aft[j+1] = "stratum_" + std::to_string(j+1);
      zn_aft(_,j) = 1.0*(stratumn == j+1);
    }
  }
  
  for (j=0; j<p2; j++) {
    String zj = base2_cov[j];
    NumericVector u = data[zj];
    covariates_aft[q+j+1] = zj;
    zn_aft(_,q+j) = u;
  }
  
  double zcrit = R::qnorm(1-alpha/2, 0, 1, 1, 0);
  
  auto f = [n, q, p, p2, covariates, covariates_aft, aft_dist, 
            recensor, swtrt_control_only, alpha, zcrit, ties, offset](
                IntegerVector stratumb, NumericVector timeb,
                IntegerVector eventb, IntegerVector treatb,
                NumericVector censor_timeb, 
                IntegerVector pdb, NumericVector pd_timeb, 
                IntegerVector swtrtb, NumericVector swtrt_timeb, 
                NumericMatrix zb, NumericMatrix zb_aft)->List {
                  int h, i, j;
                  
                  // order data by treat
                  IntegerVector order = seq(0, n-1);
                  std::sort(order.begin(), order.end(), [&](int i, int j) {
                    return (treatb[i] < treatb[j]);
                  });
                  
                  stratumb = stratumb[order];
                  timeb = timeb[order];
                  eventb = eventb[order];
                  treatb = treatb[order];
                  censor_timeb = censor_timeb[order];
                  pdb = pdb[order];
                  pd_timeb = pd_timeb[order];
                  swtrtb = swtrtb[order];
                  swtrt_timeb = swtrt_timeb[order];
                  zb = subset_matrix_by_row(zb, order);
                  zb_aft = subset_matrix_by_row(zb_aft, order);
                  
                  // time and event adjusted for treatment switching
                  NumericVector t_star = clone(timeb);
                  IntegerVector d_star = clone(eventb);
                  
                  double psi0hat = 0, psi0lower = 0, psi0upper = 0;
                  double psi1hat = 0, psi1lower = 0, psi1upper = 0;
                  
                  DataFrame data_outcome;
                  
                  // treat arms that include patients who switched treatment
                  IntegerVector treats(1);
                  treats[0] = 0;
                  if (!swtrt_control_only) {
                    treats.push_back(1);
                  }
                  
                  int K = static_cast<int>(treats.size());
                  for (h=0; h<K; h++) {
                    // post progression data
                    IntegerVector l = which((treatb == h) & (pdb == 1));
                    NumericVector time2 = timeb[l] - pd_timeb[l] + offset;
                    IntegerVector event2 = eventb[l];
                    IntegerVector swtrt2 = swtrtb[l];
                    
                    DataFrame data1 = DataFrame::create(
                      Named("time") = time2,
                      Named("event") = event2,
                      Named("swtrt") = swtrt2);
                    
                    for (j=0; j<q+p2; j++) {
                      String zj = covariates_aft[j+1];
                      NumericVector u = zb_aft(_,j);
                      data1.push_back(u[l], zj);
                    }
                    
                    List fit1 = liferegcpp(
                      data1, "", "", "time", "", "event", 
                      covariates_aft, "", "", "", aft_dist, 0, 0, alpha);
                    
                    DataFrame parest1 = DataFrame(fit1["parest"]);
                    NumericVector beta1 = parest1["beta"];
                    NumericVector sebeta1 = parest1["sebeta"];
                    double psihat = -beta1[1];
                    double psilower = -(beta1[1] + zcrit*sebeta1[1]);
                    double psiupper = -(beta1[1] - zcrit*sebeta1[1]);
                    
                    double a = exp(psihat);
                    for (i=0; i<n; i++) {
                      if (treatb[i] == h) {
                        double b2, u_star, c_star;
                        if (swtrtb[i] == 1) {
                          b2 = pdb[i] == 1 ? pd_timeb[i] : swtrt_timeb[i];
                          b2 = b2 - offset;
                          u_star = b2 + (timeb[i] - b2)*a;
                        } else {
                          u_star = timeb[i];
                        }
                        
                        if (recensor) {
                          c_star = censor_timeb[i]*std::min(1.0, a);
                          t_star[i] = std::min(u_star, c_star);
                          d_star[i] = c_star < u_star ? 0 : eventb[i];
                        } else {
                          t_star[i] = u_star;
                          d_star[i] = eventb[i];
                        }
                      }
                    }
                    
                    // update treatment-specific causal parameter estimates
                    if (h == 0) {
                      psi0hat = psihat;
                      psi0lower = psilower;
                      psi0upper = psiupper;
                    } else {
                      psi1hat = psihat;
                      psi1lower = psilower;
                      psi1upper = psiupper;
                    }
                  }
                  
                  // Cox model for hypothetical treatment effect estimate
                  data_outcome = DataFrame::create(
                    Named("stratum") = stratumb,
                    Named("t_star") = t_star,
                    Named("d_star") = d_star,
                    Named("treated") = treatb);
                  
                  for (j=0; j<p; j++) {
                    String zj = covariates[j+1];
                    NumericVector u = zb(_,j);
                    data_outcome.push_back(u, zj);
                  }
                  
                  List fit_outcome = phregcpp(
                    data_outcome, "", "stratum", "t_star", "", "d_star",
                    covariates, "", "", "", ties, 0, 0, 0, 0, 0, alpha);
                  
                  DataFrame parest = DataFrame(fit_outcome["parest"]);
                  NumericVector beta = parest["beta"];
                  NumericVector sebeta = parest["sebeta"];
                  NumericVector z = parest["z"];
                  double hrhat = exp(beta[0]);
                  double hrlower = exp(beta[0] - zcrit*sebeta[0]);
                  double hrupper = exp(beta[0] + zcrit*sebeta[0]);
                  double pvalue = 2*(1 - R::pnorm(fabs(z[0]), 0, 1, 1, 0));
                  
                  List out;
                    out = List::create(
                      Named("psihat") = psi0hat,
                      Named("psilower") = psi0lower,
                      Named("psiupper") = psi0upper,
                      Named("psi1hat") = psi1hat,
                      Named("psi1lower") = psi1lower,
                      Named("psi1upper") = psi1upper,
                      Named("hrhat") = hrhat,
                      Named("hrlower") = hrlower,
                      Named("hrupper") = hrupper,
                      Named("pvalue") = pvalue);

                  return out;
                };
  
  List out = f(stratumn, timen, eventn, treatn, censor_timen,
               pdn, pd_timen, swtrtn, swtrt_timen, zn, zn_aft);
  
  double psihat = out["psihat"];
  double psilower = out["psilower"];
  double psiupper = out["psiupper"];
  double psi1hat = out["psi1hat"];
  double psi1lower = out["psi1lower"];
  double psi1upper = out["psi1upper"];
  double hrhat = out["hrhat"];
  double hrlower = out["hrlower"];
  double hrupper = out["hrupper"];
  double pvalue = out["pvalue"];
  String psi_CI_type = "AFT model";
  
  // construct the confidence interval for HR
  String hr_CI_type;
  NumericVector hrhats(n_boot), psihats(n_boot), psi1hats(n_boot);
  if (!boot) { // use Cox model to construct CI for HR if no boot
    hr_CI_type = "Cox model";
  } else { // bootstrap the entire process to construct CI for HR
    if (seed != NA_INTEGER) set_seed(seed);
    
    IntegerVector stratumb(n), eventb(n), treatb(n), pdb(n), swtrtb(n);
    NumericVector timeb(n), censor_timeb(n), pd_timeb(n), swtrt_timeb(n);
    NumericMatrix zb(n,p), zb_aft(n,q+p2);
    
    // sort data by treatment group
    IntegerVector idx0 = which(treatn == 0);
    IntegerVector idx1 = which(treatn == 1);
    int n0 = static_cast<int>(idx0.size());
    int n1 = static_cast<int>(idx1.size());
    IntegerVector order(n);
    for (i=0; i<n0; i++) {
      order[i] = idx0[i];
    }
    for (i=0; i<n1; i++){
      order[n0+i] = idx1[i];
    }
    
    stratumn = stratumn[order];
    timen = timen[order];
    eventn = eventn[order];
    treatn = treatn[order];
    censor_timen = censor_timen[order];
    pdn = pdn[order];
    pd_timen = pd_timen[order];
    swtrtn = swtrtn[order];
    swtrt_timen = swtrt_timen[order];
    zn = subset_matrix_by_row(zn, order);
    zn_aft = subset_matrix_by_row(zn_aft, order);
    
    for (k=0; k<n_boot; k++) {
      // sample the data with replacement by treatment group
      for (i=0; i<n; i++) {
        double u = R::runif(0,1);
        if (i < n0) {
          j = static_cast<int>(std::floor(u*n0));
        } else {
          j = n0 + static_cast<int>(std::floor(u*n1));
        }
        
        stratumb[i] = stratumn[j];
        timeb[i] = timen[j];
        eventb[i] = eventn[j];
        treatb[i] = treatn[j];
        censor_timeb[i] = censor_timen[j];
        pdb[i] = pdn[j];
        pd_timeb[i] = pd_timen[j];
        swtrtb[i] = swtrtn[j];
        swtrt_timeb[i] = swtrt_timen[j];
        zb(i,_) = zn(j,_);
        zb_aft(i,_) = zn_aft(j,_);
      }
      
      out = f(stratumb, timeb, eventb, treatb, censor_timeb,
              pdb, pd_timeb, swtrtb, swtrt_timeb, zb, zb_aft);
      
      hrhats[k] = out["hrhat"];
      psihats[k] = out["psihat"];
      psi1hats[k] = out["psi1hat"];
    }
    
    // obtain bootstrap confidence interval for HR
    double loghr = log(hrhat);
    NumericVector loghrs = log(hrhats);
    double sdloghr = sd(loghrs);
    double tcrit = R::qt(1-alpha/2, n_boot-1, 1, 0);
    hrlower = exp(loghr - tcrit*sdloghr);
    hrupper = exp(loghr + tcrit*sdloghr);
    hr_CI_type = "bootstrap";
    pvalue = 2*(1 - R::pt(fabs(loghr/sdloghr), n_boot-1, 1, 0));
    
    // obtain bootstrap confidence interval for psi
    double sdpsi = sd(psihats);
    psilower = psihat - tcrit*sdpsi;
    psiupper = psihat + tcrit*sdpsi;
    psi_CI_type = "bootstrap";
    
    double sdpsi1 = sd(psi1hats);
    psi1lower = psi1hat - tcrit*sdpsi1;
    psi1upper = psi1hat + tcrit*sdpsi1;
  }
  
  List result = List::create(
    Named("psi") = psihat,
    Named("psi_lower") = psilower, 
    Named("psi_upper") = psiupper,
    Named("psi_CI_type") = psi_CI_type,
    Named("cox_pvalue") = pvalue,
    Named("hr") = hrhat,
    Named("hr_lower") = hrlower, 
    Named("hr_upper") = hrupper,
    Named("hr_CI_type") = hr_CI_type);
  
  if (!swtrt_control_only) {
    result.push_back(psi1hat, "psi_trt");
    result.push_back(psi1lower, "psi_trt_lower");
    result.push_back(psi1upper, "psi_trt_upper");
  }
  
  if (boot) {
    result.push_back(hrhats, "hr_boots");
    result.push_back(psihats, "psi_boots");
    if (!swtrt_control_only) {
      result.push_back(psi1hats, "psi_trt_boots");
    }
  }
  
  return result;
}


//' @title Simulate Survival Data for Simple Two-Stage Estimation Method
//' @description Obtains the simulated data for baseline prognosis, 
//' disease progression, treatment switching, and death.
//'
//' @param n The total sample size for two treatment arms combined.
//' @param allocation1 The number of subjects in the active treatment group 
//'   in a randomization block.
//' @param allocation2 The number of subjects in the control group in
//'   a randomization block.
//' @param rate1 The rate parameter for survival for the treatment group.
//' @param rate2 The rate parameter for survival for the control group.
//' @param p1 The parameter for disease progression for the treatment group.
//' @param p2 The parameter for disease progression for the control group.
//' @param followupTime The follow-up time for fixed follow-up design.
//' @param pswtrt1 The probability of treatment switching for the treatment 
//'   group.
//' @param pswtrt2 The probability of treatment switching for the control
//'   group.
//' @param mult1 The multiplier for post-switching survival time for the 
//'   treatment group.
//' @param mult2 The multiplier for post-switching survival time for the 
//'   control group.
//' @param gamma1 The rate parameter for dropout for the treatment group. 
//' @param gamma2 The rate parameter for dropout for the control group.
//' @param maxNumberOfIterations The number of simulation iterations.
//'   Defaults to 1000.
//' @param seed The seed to reproduce the simulation results.
//'   The seed from the environment will be used if left unspecified.
//'
//' @return A data frame of estimated hazard ratios and confidence intervals
//'   for simple two-stage estimation methods with and without recensoring
//'   and with or without recensoring for administrative censoring times only.
//' 
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' sim1 <- tsesimpsim(
//'   n = 1000, allocation1 = 1, allocation2 = 1,  
//'   rate1 = 0.0012, rate2 = 0.002, p1 = 0.6971, p2 = 0.6225, 
//'   followupTime = 2000, pswtrt1 = 0.3, pswtrt2 = 0.5,
//'   mult1 = 0.9, mult2 = 0.5, gamma1 = 0.002, gamma2 = 0.002,
//'   maxNumberOfIterations = 1000, seed = 2000)
//'
//' @export
// [[Rcpp::export]]
DataFrame tsesimpsim(const int n = NA_INTEGER,
                     const int allocation1 = 1,
                     const int allocation2 = 1,
                     const double rate1 = NA_REAL,
                     const double rate2 = NA_REAL,
                     const double p1 = NA_REAL,
                     const double p2 = NA_REAL,
                     const double followupTime = NA_REAL,
                     const double pswtrt1 = NA_REAL,
                     const double pswtrt2 = NA_REAL,
                     const double mult1 = NA_REAL,
                     const double mult2 = NA_REAL,
                     const double gamma1 = NA_REAL,
                     const double gamma2 = NA_REAL,
                     const int maxNumberOfIterations = NA_INTEGER,
                     const int seed = NA_INTEGER) {
  
  int i, j, k, iter;
  double u;
  
  if (seed != NA_INTEGER) {
    set_seed(seed);
  }
  
  IntegerVector treat(n), event(n), dropoutEvent(n);
  IntegerVector pd(n), swtrt(n);
  NumericVector survivalTime(n), dropoutTime(n);
  NumericVector pdTime(n), swtrtTime(n);
  NumericVector time(n), pd_time(n), swtrt_time(n);
  
  int N = 3*maxNumberOfIterations;
  LogicalVector recensor(N), admin_recensor_only(N);
  NumericVector hr(N), hr_lower(N), hr_upper(N);

  for (iter=0; iter<maxNumberOfIterations; iter++) {
    int b1 = allocation1, b2 = allocation2;
    for (i=0; i<n; i++) {

      // generate treatment indicators using stratified block randomization
      u = R::runif(0,1);
      if (u <= b1/(b1+b2+0.0)) {
        treat[i] = 1;
        b1--;
      } else {
        treat[i] = 0;
        b2--;
      }
      
      // start a new block after depleting the current block
      if (b1+b2==0) {
        b1 = allocation1;
        b2 = allocation2;
      }
      
      // generate survival times from exponential distribution
      u = R::runif(0,1);
      double rate = treat[i] == 1 ? rate1 : rate2;
      survivalTime[i] = -log(u)/rate;
      survivalTime[i] = std::round(survivalTime[i]);
      if (survivalTime[i] == 0) survivalTime[i] = 1;
      
      // generate disease progression time - scheduled visits are every 21 days
      k = static_cast<int>(std::floor(survivalTime[i] / 21));
      
      NumericVector breaks(k);
      double p = treat[i] == 1 ? p1 : p2;
      for (j=0; j<k; j++) {
        breaks[j] = 1 - pow(p, j+1);
      }
      
      u = R::runif(0,1);
      NumericVector u1(1); 
      u1[0] = u;
      
      IntegerVector j1 = findInterval3(u1, breaks);
      j = j1[0];
      
      if (j < k) {
        pdTime[i] = (j+1)*21;
        pd[i] = 1;
      } else {
        pdTime[i] = survivalTime[i];
        pd[i] = 0;
      }
      
      // generate treatment switch indicators
      if (pd[i] == 1) {
        double pswtrt = treat[i] == 1 ? pswtrt1 : pswtrt2;
        swtrt[i] = static_cast<int>(R::rbinom(1, pswtrt));
        if (swtrt[i] == 1) swtrtTime[i] = pdTime[i];
      }
      
      if (swtrt[i] == 1) {
        // generate survival time after treatment switch
        double mult = treat[i] == 1 ? mult1 : mult2;
        u = R::runif(0,1);
        double survivalTime2 = -log(u)/(rate*mult);
        survivalTime2 = std::round(survivalTime2);
        if (survivalTime2 == 0) survivalTime2 = 1;
        survivalTime[i] = swtrtTime[i] + survivalTime2;
      }
      
      // generate dropout time
      double gamma = treat[i] == 1 ? gamma1 : gamma2;
      u = R::runif(0,1);
      dropoutTime[i] = -log(u)/(gamma);
      dropoutTime[i] = std::round(dropoutTime[i]);
      if (dropoutTime[i] == 0) dropoutTime[i] = 1;
      
      // censoring
      if (survivalTime[i] <= dropoutTime[i] && 
          survivalTime[i] <= followupTime) {
        time[i] = survivalTime[i];
        event[i] = 1;
        dropoutEvent[i] = 0;
      } else if (dropoutTime[i] < survivalTime[i] &&
        dropoutTime[i] <= followupTime) {
        time[i] = dropoutTime[i];
        event[i] = 0;
        dropoutEvent[i] = 1;
      } else {
        time[i] = followupTime;
        event[i] = 0;
        dropoutEvent[i] = 0;
      }
      
      if (pd[i] == 1 && pdTime[i] <= time[i]) {
        pd_time[i] = pdTime[i];
      } else {
        pd[i] = 0;
        pd_time[i] = NA_REAL;
      }

      if (swtrt[i] == 1 && swtrtTime[i] <= time[i]) {
        swtrt_time[i] = swtrtTime[i];
      } else {
        swtrt[i] = 0;
        swtrt_time[i] = NA_REAL;
      }
    }
    
    DataFrame data1 = DataFrame::create(
      Named("time") = time,
      Named("event") = event,
      Named("treat") = treat,
      Named("censor_time") = followupTime,
      Named("pd") = pd,
      Named("pd_time") = pd_time,
      Named("swtrt") = swtrt,
      Named("swtrt_time") = swtrt_time
    );
    
    List out1 = tsesimpcpp1(
      data1, "", "time", "event", "treat", "censor_time", "pd", 
      "pd_time", "swtrt", "swtrt_time", "", "", "weibull", 1, 
      0, 1, 0, 0.05, "efron", 1, 0, 1000, NA_INTEGER);

    List out2 = tsesimpcpp1(
      data1, "", "time", "event", "treat", "censor_time", "pd", 
      "pd_time", "swtrt", "swtrt_time", "", "",  "weibull", 1, 
      1, 1, 0, 0.05, "efron", 1, 0, 1000, NA_INTEGER);
    
    List out3 = tsesimpcpp1(
      data1, "", "time", "event", "treat", "censor_time", "pd", 
      "pd_time", "swtrt", "swtrt_time", "", "",  "weibull", 1, 
      1, 0, 0, 0.05, "efron", 1, 0, 1000, NA_INTEGER);
    
    recensor[iter] = 0;
    admin_recensor_only[iter] = 1;
    hr[iter] = out1["hr"];
    hr_lower[iter] = out1["hr_lower"];
    hr_upper[iter] = out1["hr_upper"];
    
    recensor[iter+maxNumberOfIterations] = 1;
    admin_recensor_only[iter+maxNumberOfIterations] = 1;
    hr[iter+maxNumberOfIterations] = out2["hr"];
    hr_lower[iter+maxNumberOfIterations] = out2["hr_lower"];
    hr_upper[iter+maxNumberOfIterations] = out2["hr_upper"];
    
    recensor[iter+2*maxNumberOfIterations] = 1;
    admin_recensor_only[iter+2*maxNumberOfIterations] = 0;
    hr[iter+2*maxNumberOfIterations] = out3["hr"];
    hr_lower[iter+2*maxNumberOfIterations] = out3["hr_lower"];
    hr_upper[iter+2*maxNumberOfIterations] = out3["hr_upper"];
  }
  
  DataFrame result = DataFrame::create(
    Named("recensor") = recensor,
    Named("admin_recensor_only") = admin_recensor_only,
    Named("hr") = hr,
    Named("hr_lower") = hr_lower,
    Named("hr_upper") = hr_upper
  );
  
  return result;
}
