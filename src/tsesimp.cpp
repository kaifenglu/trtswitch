#include "utilities.h"
#include "survival_analysis.h"

using namespace Rcpp;


// [[Rcpp::export]]
List tsesimpcpp(const DataFrame data,
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

  int i, j, k, l, n = data.nrow();

  int p = static_cast<int>(base_cov.size());
  if (p == 1 && (base_cov[0] == "" || base_cov[0] == "none")) p = 0;

  int p2 = static_cast<int>(base2_cov.size());
  if (p2 == 1 && (base2_cov[0] == "" || base2_cov[0] == "none")) p2 = 0;

  int p_stratum = static_cast<int>(stratum.size());
  
  bool has_stratum;
  IntegerVector stratumn(n);
  DataFrame u_stratum;
  IntegerVector d(p_stratum);
  IntegerMatrix stratan(n,p_stratum);
  if (p_stratum == 1 && (stratum[0] == "" || stratum[0] == "none")) {
    has_stratum = 0;
    stratumn.fill(1);
    d[0] = 1;
    stratan(_,0) = stratumn;
  } else {
    List out = bygroup(data, stratum);
    has_stratum = 1;
    stratumn = out["index"];
    u_stratum = DataFrame(out["lookup"]);
    d = out["nlevels"];
    stratan = as<IntegerMatrix>(out["indices"]);
  }
  
  IntegerVector stratumn_unique = unique(stratumn);
  int nstrata = static_cast<int>(stratumn_unique.size());
  
  bool has_time = hasVariable(data, time);
  bool has_event = hasVariable(data, event);
  bool has_treat = hasVariable(data, treat);
  bool has_censor_time = hasVariable(data, censor_time);
  bool has_pd = hasVariable(data, pd);
  bool has_pd_time = hasVariable(data, pd_time);
  bool has_swtrt = hasVariable(data, swtrt);
  bool has_swtrt_time = hasVariable(data, swtrt_time);

  if (!has_time) {
    stop("data must contain the time variable");
  }

  if (TYPEOF(data[time]) != INTSXP && TYPEOF(data[time]) != REALSXP) {
    stop("time must take numeric values");
  }

  NumericVector timenz = data[time];
  NumericVector timen = clone(timenz);
  if (is_true(any(timen <= 0.0))) {
    stop("time must be positive");
  }

  if (!has_event) {
    stop("data must contain the event variable");
  }

  if (TYPEOF(data[event]) != INTSXP && TYPEOF(data[event]) != LGLSXP) {
    stop("event must take integer or logical values");
  }

  IntegerVector eventnz = data[event];
  IntegerVector eventn = clone(eventnz);
  if (is_true(any((eventn != 1) & (eventn != 0)))) {
    stop("event must be 1 or 0");
  }

  if (is_true(all(eventn == 0))) {
    stop("at least 1 event is needed");
  }

  if (!has_treat) {
    stop("data must contain the treat variable");
  }

  // create the numeric treat variable
  if (!has_treat) stop("data must contain the treat variable");
  IntegerVector treatn(n);
  IntegerVector treatwi;
  NumericVector treatwn;
  StringVector treatwc;
  if (TYPEOF(data[treat]) == LGLSXP || TYPEOF(data[treat]) == INTSXP) {
    IntegerVector treatv = data[treat];
    treatwi = unique(treatv);
    if (treatwi.size() != 2) {
      stop("treat must have two and only two distinct values");
    }
    // special handling for 1/0 treatment coding
    if (is_true(all((treatwi == 0) | (treatwi == 1)))) {
      treatwi = IntegerVector::create(1,0);
      treatn = 2 - treatv;
    } else {
      treatwi.sort();
      treatn = match(treatv, treatwi);
    }
  } else if (TYPEOF(data[treat]) == REALSXP) {
    NumericVector treatv = data[treat];
    treatwn = unique(treatv);
    if (treatwn.size() != 2) {
      stop("treat must have two and only two distinct values");
    }
    // special handling for 1/0 treatment coding
    if (is_true(all((treatwn == 0) | (treatwn == 1)))) {
      treatwn = NumericVector::create(1,0);
      treatn = 2 - as<IntegerVector>(treatv);
    } else {
      treatwn.sort();
      treatn = match(treatv, treatwn);
    }
  } else if (TYPEOF(data[treat]) == STRSXP) {
    StringVector treatv = data[treat];
    treatwc = unique(treatv);
    if (treatwc.size() != 2) {
      stop("treat must have two and only two distinct values");
    }
    treatwc.sort();
    treatn = match(treatv, treatwc);
  } else {
    stop("incorrect type for the treat variable in the input data");
  }
  
  treatn = 2 - treatn; // use the 1/0 treatment coding
  
  if (!has_censor_time) {
    stop("data must contain the censor_time variable");
  }

  if (TYPEOF(data[censor_time]) != INTSXP &&
      TYPEOF(data[censor_time]) != REALSXP) {
    stop("censor_time must take numeric values");
  }

  NumericVector censor_timenz = data[censor_time];
  NumericVector censor_timen = clone(censor_timenz);
  if (is_true(any(censor_timen < timen))) {
    stop("censor_time must be greater than or equal to time");
  }

  if (!admin_recensor_only) {
    for (i=0; i<n; i++) {
      if (eventn[i] == 0) { // use the actual censoring time for dropouts
        censor_timen[i] = timen[i];
      }
    }
  }

  if (!has_pd) {
    stop("data must contain the pd variable");
  }

  if (TYPEOF(data[pd]) != INTSXP && TYPEOF(data[pd]) != LGLSXP) {
    stop("pd must take integer or logical values");
  }

  IntegerVector pdnz = data[pd];
  IntegerVector pdn = clone(pdnz);
  if (is_true(any((pdn != 1) & (pdn != 0)))) {
    stop("pd must be 1 or 0");
  }

  if (!has_pd_time) {
    stop("data must contain the pd_time variable");
  }

  if (TYPEOF(data[pd_time]) != INTSXP && TYPEOF(data[pd_time]) != REALSXP) {
    stop("pd_time must take numeric values");
  }

  NumericVector pd_timenz = data[pd_time];
  NumericVector pd_timen = clone(pd_timenz);
  for (i=0; i<n; i++) {
    if (pdn[i] == 1 && pd_timen[i] < 0.0) {
      stop("pd_time must be nonnegative");
    }
    
    if (pdn[i] == 1 && pd_timen[i] > timen[i]) {
      stop("pd_time must be less than or equal to time");
    }
  }

  if (!has_swtrt) {
    stop("data must contain the swtrt variable");
  }

  if (TYPEOF(data[swtrt]) != INTSXP && TYPEOF(data[swtrt]) != LGLSXP) {
    stop("swtrt must take integer or logical values");
  }

  IntegerVector swtrtnz = data[swtrt];
  IntegerVector swtrtn = clone(swtrtnz);
  if (is_true(any((swtrtn != 1) & (swtrtn != 0)))) {
    stop("swtrt must be 1 or 0");
  }

  if (is_false(any((pdn == 1) & (swtrtn == 1) & (treatn == 0)))) {
    stop("at least 1 pd and swtrt is needed in the control group");
  }

  if (!swtrt_control_only) {
    if (is_false(any((pdn == 1) & (swtrtn == 1) & (treatn == 1)))) {
      stop("at least 1 pd and swtrt is needed in the treatment group");
    }
  }

  if (!has_swtrt_time) {
    stop("data must contain the swtrt_time variable");
  }

  if (TYPEOF(data[swtrt_time]) != INTSXP &&
      TYPEOF(data[swtrt_time]) != REALSXP) {
    stop("swtrt_time must take numeric values");
  }

  NumericVector swtrt_timenz = data[swtrt_time];
  NumericVector swtrt_timen = clone(swtrt_timenz);
  for (i=0; i<n; i++) {
    if (swtrtn[i] == 1 && swtrt_timen[i] < 0.0) {
      stop("swtrt_time must be nonnegative");
    }
    
    if (swtrtn[i] == 1 && swtrt_timen[i] > timen[i]) {
      stop("swtrt_time must be less than or equal to time");
    }
  }

  // covariates for the Cox model containing treat and base_cov
  StringVector covariates(p+1);
  NumericMatrix zn(n,p);
  covariates[0] = "treated";
  for (j=0; j<p; j++) {
    String zj = base_cov[j];
    if (!hasVariable(data, zj)) {
      stop("data must contain the variables in base_cov");
    }
    if (zj == treat) {
      stop("treat should be excluded from base_cov");
    }
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
    if (!hasVariable(data, zj)) {
      stop("data must contain the variables in base2_cov");
    }
    if (zj == treat) {
      stop("treat should be excluded from base2_cov");
    }
    NumericVector u = data[zj];
    covariates_aft[q+j+1] = zj;
    zn_aft(_,q+j) = u;
  }

  std::string dist = aft_dist;
  std::for_each(dist.begin(), dist.end(), [](char & c) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });

  if ((dist == "log-logistic") || (dist == "llogistic")) {
    dist = "loglogistic";
  } else if  ((dist == "log-normal") || (dist == "lnormal")) {
    dist = "lognormal";
  }

  if (!((dist == "exponential") || (dist == "weibull") ||
      (dist == "lognormal") || (dist == "loglogistic"))) {
    stop("aft_dist must be exponential, weibull, lognormal, or loglogistic");
  }

  if (alpha <= 0.0 || alpha >= 0.5) {
    stop("alpha must lie between 0 and 0.5");
  }

  if (ties != "efron" && ties != "breslow") {
    stop("ties must be efron or breslow");
  }
  
  if (offset < 0.0) {
    stop("offset must be nonnegative");
  }
  
  if (n_boot < 100) {
    stop("n_boot must be greater than or equal to 100");
  }
  

  DataFrame lr = lrtest(data, "", stratum, treat, time, event, 0, 0);
  double logRankPValue = as<double>(lr["logRankPValue"]);
  
  double zcrit = R::qnorm(1-alpha/2, 0, 1, 1, 0);

  k = -1; // indicate the observed data
  auto f = [&k, data, has_stratum, stratum, p_stratum, u_stratum, 
            treat, treatwi, treatwn, treatwc, 
            n, q, p, p2, covariates, covariates_aft, dist, 
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
                  NumericVector time_ts = clone(timeb);
                  IntegerVector event_ts = clone(eventb);
                  
                  double psi0hat = 0, psi0lower = 0, psi0upper = 0;
                  double psi1hat = 0, psi1lower = 0, psi1upper = 0;
                  
                  // initialize data_aft
                  List data_aft(2);
                  if (k == -1) {
                    for (h=0; h<2; h++) {
                      List data_aftx = List::create(
                        Named("data") = R_NilValue,
                        Named(treat) = R_NilValue
                      );
                      
                      if (TYPEOF(data[treat]) == LGLSXP ||
                          TYPEOF(data[treat]) == INTSXP) {
                        data_aftx[treat] = treatwi[1-h];
                      } else if (TYPEOF(data[treat]) == REALSXP) {
                        data_aftx[treat] = treatwn[1-h];
                      } else if (TYPEOF(data[treat]) == STRSXP) {
                        data_aftx[treat] = treatwc[1-h];
                      }
                      
                      data_aft[h] = data_aftx;
                    }
                  }
                  
                  // initialize fit_aft
                  List fit_aft(2);
                  if (k == -1) {
                    for (h=0; h<2; h++) {
                      List fit_aftx = List::create(
                        Named("fit") = R_NilValue,
                        Named(treat) = R_NilValue
                      );
                      
                      if (TYPEOF(data[treat]) == LGLSXP ||
                          TYPEOF(data[treat]) == INTSXP) {
                        fit_aftx[treat] = treatwi[1-h];
                      } else if (TYPEOF(data[treat]) == REALSXP) {
                        fit_aftx[treat] = treatwn[1-h];
                      } else if (TYPEOF(data[treat]) == STRSXP) {
                        fit_aftx[treat] = treatwc[1-h];
                      }
                      
                      fit_aft[h] = fit_aftx;
                    }
                  }
                  
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
                    IntegerVector sub = which((treatb == h) & (pdb == 1));
                    int nsub = static_cast<int>(sub.size());
                    NumericVector time2(nsub); // post-pd survival time
                    IntegerVector event2(nsub);
                    IntegerVector swtrt2(nsub);
                    for (i=0; i<nsub; i++) {
                      j = sub[i];
                      time2[i] = timeb[j] - pd_timeb[j] + offset;
                      event2[i] = eventb[j];
                      swtrt2[i] = swtrtb[j];
                    }
                    
                    DataFrame data1 = DataFrame::create(
                      Named("time") = time2,
                      Named("event") = event2,
                      Named("swtrt") = swtrt2);
                    
                    for (j=0; j<q+p2; j++) {
                      String zj = covariates_aft[j+1];
                      NumericVector u = zb_aft(_,j);
                      data1.push_back(u[sub], zj);
                    }
                    
                    List fit1 = liferegcpp(
                      data1, "", "", "time", "", "event", 
                      covariates_aft, "", "", "", dist, 0, 0, alpha);
                    
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
                          time_ts[i] = std::min(u_star, c_star);
                          event_ts[i] = c_star < u_star ? 0 : eventb[i];
                        } else {
                          time_ts[i] = u_star;
                          event_ts[i] = eventb[i];
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
                    
                    // update data_aft and fit_aft
                    if (k == -1) {
                      IntegerVector stratum2 = stratumb[sub];
                      
                      if (has_stratum) {
                        for (i=0; i<p_stratum; i++) {
                          String s = stratum[i];
                          if (TYPEOF(data[s]) == INTSXP) {
                            IntegerVector stratumwi = u_stratum[s];
                            data1.push_back(stratumwi[stratum2-1], s);
                          } else if (TYPEOF(data[s]) == REALSXP) {
                            NumericVector stratumwn = u_stratum[s];
                            data1.push_back(stratumwn[stratum2-1], s);
                          } else if (TYPEOF(data[s]) == STRSXP) {
                            StringVector stratumwc = u_stratum[s];
                            data1.push_back(stratumwc[stratum2-1], s);
                          }
                        }
                      }
                      
                      List data_aftx = List::create(
                        Named("data") = data1,
                        Named(treat) = R_NilValue
                      );
                      
                      if (TYPEOF(data[treat]) == LGLSXP ||
                          TYPEOF(data[treat]) == INTSXP) {
                        data_aftx[treat] = treatwi[1-h];
                      } else if (TYPEOF(data[treat]) == REALSXP) {
                        data_aftx[treat] = treatwn[1-h];
                      } else if (TYPEOF(data[treat]) == STRSXP) {
                        data_aftx[treat] = treatwc[1-h];
                      }
                      
                      data_aft[h] = data_aftx;
                      
                      List fit_aftx = List::create(
                        Named("fit") = fit1,
                        Named(treat) = R_NilValue
                      );
                      
                      if (TYPEOF(data[treat]) == LGLSXP ||
                          TYPEOF(data[treat]) == INTSXP) {
                        fit_aftx[treat] = treatwi[1-h];
                      } else if (TYPEOF(data[treat]) == REALSXP) {
                        fit_aftx[treat] = treatwn[1-h];
                      } else if (TYPEOF(data[treat]) == STRSXP) {
                        fit_aftx[treat] = treatwc[1-h];
                      }
                      
                      fit_aft[h] = fit_aftx;
                    }
                  }
                  
                  // Cox model for hypothetical treatment effect estimate
                  data_outcome = DataFrame::create(
                    Named("time") = time_ts,
                    Named("event") = event_ts,
                    Named("treated") = treatb);

                  if (has_stratum) {
                    for (i=0; i<p_stratum; i++) {
                      String s = stratum[i];
                      if (TYPEOF(data[s]) == INTSXP) {
                        IntegerVector stratumwi = u_stratum[s];
                        data_outcome.push_back(stratumwi[stratumb-1], s);
                      } else if (TYPEOF(data[s]) == REALSXP) {
                        NumericVector stratumwn = u_stratum[s];
                        data_outcome.push_back(stratumwn[stratumb-1], s);
                      } else if (TYPEOF(data[s]) == STRSXP) {
                        StringVector stratumwc = u_stratum[s];
                        data_outcome.push_back(stratumwc[stratumb-1], s);
                      }
                    }
                  }
                  
                  for (j=0; j<p; j++) {
                    String zj = covariates[j+1];
                    NumericVector u = zb(_,j);
                    data_outcome.push_back(u, zj);
                  }
                  
                  List fit_outcome = phregcpp(
                    data_outcome, "", stratum, "time", "", "event",
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
                  if (k == -1) {
                    out = List::create(
                      Named("data_aft") = data_aft,
                      Named("fit_aft") = fit_aft,
                      Named("data_outcome") = data_outcome,
                      Named("fit_outcome") = fit_outcome,
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
                  } else {
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
                  }

                  return out;
                };

  List out = f(stratumn, timen, eventn, treatn, censor_timen,
               pdn, pd_timen, swtrtn, swtrt_timen, zn, zn_aft);

  List data_aft = out["data_aft"];
  List fit_aft = out["fit_aft"];
  DataFrame data_outcome = DataFrame(out["data_outcome"]);
  List fit_outcome = out["fit_outcome"];
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
  
  IntegerVector treated = data_outcome["treated"];
  if (TYPEOF(data[treat]) == LGLSXP || TYPEOF(data[treat]) == INTSXP) {
    data_outcome.push_back(treatwi[1-treated], treat);
  } else if (TYPEOF(data[treat]) == REALSXP) {
    data_outcome.push_back(treatwn[1-treated], treat);
  } else if (TYPEOF(data[treat]) == STRSXP) {
    data_outcome.push_back(treatwc[1-treated], treat);
  }
  
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
        for (l=0; l<p; l++) {
          zb(i,l) = zn(j,l);
        }
        for (l=0; l<q+p2; l++) {
          zb_aft(i,l) = zn_aft(j,l);
        }
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

  List settings = List::create(
    Named("aft_dist") = aft_dist,
    Named("strata_main_effect_only") = strata_main_effect_only,
    Named("recensor") = recensor,
    Named("admin_recensor_only") = admin_recensor_only,
    Named("swtrt_control_only") = swtrt_control_only,
    Named("alpha") = alpha,
    Named("ties") = ties,
    Named("offset") = offset,
    Named("boot") = boot,
    Named("n_boot") = n_boot,
    Named("seed") = seed);

  List result = List::create(
    Named("psi") = psihat,
    Named("psi_CI") = NumericVector::create(psilower, psiupper),
    Named("psi_CI_type") = psi_CI_type,
    Named("logrank_pvalue") = 2*std::min(logRankPValue, 1-logRankPValue),
    Named("cox_pvalue") = pvalue,
    Named("hr") = hrhat,
    Named("hr_CI") = NumericVector::create(hrlower, hrupper),
    Named("hr_CI_type") = hr_CI_type,
    Named("data_aft") = data_aft,
    Named("fit_aft") = fit_aft,
    Named("data_outcome") = data_outcome,
    Named("fit_outcome") = fit_outcome,
    Named("settings") = settings);

  if (!swtrt_control_only) {
    result.push_back(psi1hat, "psi_trt");
    NumericVector psi1_CI = NumericVector::create(psi1lower, psi1upper);
    result.push_back(psi1_CI, "psi_trt_CI");
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
