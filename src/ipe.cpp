#include "utilities.h"
#include "survival_analysis.h"

using namespace Rcpp;


double est_psi_ipe(
    const double psi,
    const int n,
    const int q,
    const int p,
    const NumericVector& time,
    const IntegerVector& event,
    const IntegerVector& treat,
    const NumericVector& rx,
    const NumericVector& censor_time,
    const StringVector& covariates_aft,
    const NumericMatrix& zb_aft,
    const std::string dist,
    const double treat_modifier,
    const bool recensor,
    const bool autoswitch,
    const double alpha) {

  DataFrame Tstar = unswitched(psi*treat_modifier, n, time, event, treat,
                               rx, censor_time, recensor, autoswitch);

  NumericVector t_star = Tstar["t_star"];
  IntegerVector d_star = Tstar["d_star"];

  DataFrame data = DataFrame::create(
    Named("time") = t_star,
    Named("event") = d_star,
    Named("treated") = treat);

  for (int j=0; j<q+p; j++) {
    String zj = covariates_aft[j+1];
    NumericVector u = zb_aft(_,j);
    data.push_back(u, zj);
  }

  List fit = liferegcpp(
    data, "", "", "time", "", "event",
    covariates_aft, "", "", "", dist, 0, 0, alpha);
  
  DataFrame parest = DataFrame(fit["parest"]);
  NumericVector beta = parest["beta"];
  double psihat = -beta[1]/treat_modifier;

  return psihat;
}


// [[Rcpp::export]]
List ipecpp(const DataFrame data,
            const StringVector& stratum = "",
            const std::string time = "time",
            const std::string event = "event",
            const std::string treat = "treat",
            const std::string rx = "rx",
            const std::string censor_time = "censor_time",
            const StringVector& base_cov = "",
            const std::string aft_dist = "weibull",
            const bool strata_main_effect_only = 1,
            const double treat_modifier = 1,
            const bool recensor = 1,
            const bool admin_recensor_only = 1,
            const bool autoswitch = 1,
            const double alpha = 0.05,
            const std::string ties = "efron",
            const double tol = 1.0e-6,
            const bool boot = 0,
            const int n_boot = 1000,
            const int seed = NA_INTEGER) {

  int i, j, k, l, n = data.nrow();
  int p = static_cast<int>(base_cov.size());
  if (p == 1 && (base_cov[0] == "" || base_cov[0] == "none")) p = 0;

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
  bool has_rx = hasVariable(data, rx);
  bool has_censor_time = hasVariable(data, censor_time);

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
  
  if (!has_rx) {
    stop("data must contain the rx variable");
  }

  if (TYPEOF(data[rx]) != INTSXP && TYPEOF(data[rx]) != REALSXP) {
    stop("rx must take numeric values");
  }

  NumericVector rxnz = data[rx];
  NumericVector rxn = clone(rxnz);
  if (is_true(any((rxn < 0.0) | (rxn > 1.0)))) {
    stop("rx must take values between 0 and 1");
  }

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

  // covariates for the accelerated failure time model
  // including treat, stratum, and base_cov
  int q; // number of columns corresponding to the strata effects
  if (strata_main_effect_only) {
    q = sum(d - 1);
  } else {
    q = nstrata - 1;
  }

  StringVector covariates_aft(q+p+1);
  NumericMatrix zn_aft(n,q+p);
  covariates_aft[0] = "treated";
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

  for (j=0; j<p; j++) {
    String zj = base_cov[j];
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
  
  if (treat_modifier <= 0.0) {
    stop("treat_modifier must be positive");
  }
  
  if (alpha <= 0.0 || alpha >= 0.5) {
    stop("alpha must lie between 0 and 0.5");
  }

  if (ties != "efron" && ties != "breslow") {
    stop("ties must be efron or breslow");
  }
  
  if (tol <= 0.0) {
    stop("tol must be positive");
  }
  
  if (n_boot < 100) {
    stop("n_boot must be greater than or equal to 100");
  }
  

  DataFrame lr = lrtest(data, "", stratum, treat, time, event, 0, 0);
  double logRankPValue = as<double>(lr["logRankPValue"]);
  
  double zcrit = R::qnorm(1-alpha/2, 0, 1, 1, 0);

  k = -1;
  auto f = [&k, data, has_stratum, stratum, p_stratum, u_stratum, 
            n, q, p, covariates, covariates_aft, dist,
            treat_modifier, recensor, autoswitch, alpha, ties, tol](
                IntegerVector stratumb, NumericVector timeb,
                IntegerVector eventb, IntegerVector treatb,
                NumericVector rxb, NumericVector censor_timeb,
                NumericMatrix zb, NumericMatrix zb_aft)->List {
                  int i, j;
                  
                  // estimate psi
                  auto g = [n, q, p, timeb, eventb, treatb, rxb, 
                            censor_timeb, covariates_aft, zb_aft, 
                            dist, treat_modifier, recensor, 
                            autoswitch, alpha](double psi)->double{
                              double psinew = est_psi_ipe(
                                psi, n, q, p, timeb, eventb, treatb, rxb,
                                censor_timeb, covariates_aft, zb_aft, 
                                dist, treat_modifier, recensor, 
                                autoswitch, alpha);
                              return psinew - psi;
                            };

                  double psihat = brent(g, -3, 3, tol);

                  // construct the counterfactual survival times
                  DataFrame Sstar = untreated(
                    psihat*treat_modifier, timeb, eventb, treatb,
                    rxb, censor_timeb, recensor, autoswitch);

                  NumericVector t_star = Sstar["t_star"];
                  IntegerVector d_star = Sstar["d_star"];
                  
                  // obtain the Kaplan-Meier estimates
                  DataFrame kmstar;
                  if (k == -1) {
                    DataFrame kmdata = DataFrame::create(
                      Named("treat") = treatb,
                      Named("time") = t_star,
                      Named("event") = d_star);

                    kmstar = kmest(kmdata, "", "treat", "time",
                                   "event", "log-log", 1-alpha);
                  }
                  
                  // run Cox model to obtain the hazard ratio estimate
                  DataFrame Tstar = unswitched(
                    psihat*treat_modifier, n, timeb, eventb, treatb,
                    rxb, censor_timeb, recensor, autoswitch);
                  
                  t_star = Tstar["t_star"];
                  d_star = Tstar["d_star"];
                  
                  DataFrame data_outcome = DataFrame::create(
                    Named("time") = t_star,
                    Named("event") = d_star,
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
                  NumericVector z = parest["z"];
                  double hrhat = exp(beta[0]/treat_modifier);
                  double pvalue = 2*(1 - R::pnorm(fabs(z[0]), 0, 1, 1, 0));
                  
                  List out;
                  if (k == -1) {
                    out = List::create(
                      Named("Sstar") = Sstar,
                      Named("kmstar") = kmstar,
                      Named("data_outcome") = data_outcome,
                      Named("fit_outcome") = fit_outcome,
                      Named("psihat") = psihat,
                      Named("hrhat") = hrhat,
                      Named("pvalue") = pvalue);
                  } else {
                    out = List::create(
                      Named("psihat") = psihat,
                      Named("hrhat") = hrhat,
                      Named("pvalue") = pvalue);
                  }

                  return out;
                };

  List out = f(stratumn, timen, eventn, treatn, rxn, censor_timen, zn,
               zn_aft);

  DataFrame Sstar = DataFrame(out["Sstar"]);
  DataFrame kmstar = DataFrame(out["kmstar"]);
  DataFrame data_outcome = DataFrame(out["data_outcome"]);
  List fit_outcome = out["fit_outcome"];
  double psihat = out["psihat"];
  double zipe = R::qnorm(logRankPValue, 0, 1, 1, 0);
  double sepsi = psihat/zipe;
  double psilower = psihat - zcrit*sepsi;
  double psiupper = psihat + zcrit*sepsi;
  String psi_CI_type = "log-rank p-value";
  double hrhat = out["hrhat"];
  double pvalue = out["pvalue"];
  
  IntegerVector treated = data_outcome["treated"];
  if (TYPEOF(data[treat]) == LGLSXP || TYPEOF(data[treat]) == INTSXP) {
    data_outcome.push_back(treatwi[1-treated], treat);
  } else if (TYPEOF(data[treat]) == REALSXP) {
    data_outcome.push_back(treatwn[1-treated], treat);
  } else if (TYPEOF(data[treat]) == STRSXP) {
    data_outcome.push_back(treatwc[1-treated], treat);
  }
  
  // construct the confidence interval for HR
  double hrlower, hrupper;
  NumericVector hrhats(n_boot), psihats(n_boot);
  String hr_CI_type;
  if (!boot) { // use log-rank p-value to construct CI for HR if no boot
    double loghr = log(hrhat);
    double zcox = R::qnorm(logRankPValue, 0, 1, 1, 0);
    double seloghr = loghr/zcox;
    hrlower = exp(loghr - zcrit*seloghr);
    hrupper = exp(loghr + zcrit*seloghr);
    hr_CI_type = "log-rank p-value";
  } else { // bootstrap the entire process to construct CI for HR
    if (seed != NA_INTEGER) set_seed(seed);

    IntegerVector stratumb(n), treatb(n), eventb(n);
    NumericVector timeb(n), rxb(n), censor_timeb(n);
    NumericMatrix zb(n,p), zb_aft(n,q+p);

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
    rxn = rxn[order];
    censor_timen = censor_timen[order];
    zn = subset_matrix_by_row(zn, order);
    zn_aft = subset_matrix_by_row(zn_aft, order);

    for (k=0; k<n_boot; k++) {
      // sample the data with replacement by treatment group
      for (i=0; i<n; i++) {
        double u = R::runif(0,1);
        if (treatn[i] == 0) {
          j = static_cast<int>(std::floor(u*n0));
        } else {
          j = n0 + static_cast<int>(std::floor(u*n1));
        }

        stratumb[i] = stratumn[j];
        timeb[i] = timen[j];
        eventb[i] = eventn[j];
        treatb[i] = treatn[j];
        rxb[i] = rxn[j];
        censor_timeb[i] = censor_timen[j];
        for (l=0; l<p; l++) {
          zb(i,l) = zn(j,l);
        }
        for (l=0; l<q+p; l++) {
          zb_aft(i,l) = zn_aft(j,l);
        }
      }

      List out = f(stratumb, timeb, eventb, treatb, rxb, censor_timeb, zb,
                   zb_aft);
      hrhats[k] = out["hrhat"];
      psihats[k] = out["psihat"];
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
  }

  List settings = List::create(
    Named("aft_dist") = aft_dist,
    Named("strata_main_effect_only") = strata_main_effect_only,
    Named("treat_modifer") = treat_modifier,
    Named("recensor") = recensor,
    Named("admin_recensor_only") = admin_recensor_only,
    Named("autoswitch") = autoswitch,
    Named("alpha") = alpha,
    Named("ties") = ties,
    Named("tol") = tol,
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
    Named("Sstar") = Sstar,
    Named("kmstar") = kmstar,
    Named("data_outcome") = data_outcome,
    Named("fit_outcome") = fit_outcome,
    Named("settings") = settings);

  if (boot) {
    result.push_back(hrhats, "hr_boots");
    result.push_back(psihats, "psi_boots");
  }
  
  return result;
}
