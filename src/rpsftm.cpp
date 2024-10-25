#include "utilities.h"
#include "survival_analysis.h"

using namespace Rcpp;


DataFrame untreated(
    const double psi,
    const NumericVector& time,
    const IntegerVector& event,
    const IntegerVector& treat,
    const NumericVector& rx,
    const NumericVector& censor_time,
    const bool recensor,
    const bool autoswitch) {

  NumericVector u = time*((1 - rx) + rx*exp(psi));
  NumericVector t_star = clone(u);
  IntegerVector d_star = clone(event);

  if (recensor) {
    NumericVector c_star = pmin(censor_time, censor_time*exp(psi));

    if (autoswitch) {
      NumericVector rx1 = rx[treat == 1];
      NumericVector rx0 = rx[treat == 0];
      if (is_true(all(rx1 == 1.0))) c_star[treat == 1] = R_PosInf;
      if (is_true(all(rx0 == 0.0))) c_star[treat == 0] = R_PosInf;
    }

    t_star = pmin(u, c_star);
    d_star[c_star < u] = 0;
  }

  DataFrame result = DataFrame::create(
    Named("t_star") = t_star,
    Named("d_star") = d_star,
    Named("treat") = treat
  );

  return result;
}


double est_eqn(
    const double psi,
    const IntegerVector& stratum,
    const NumericVector& time,
    const IntegerVector& event,
    const IntegerVector& treat,
    const NumericVector& rx,
    const NumericVector& censor_time,
    const double treat_modifier,
    const bool recensor,
    const bool autoswitch,
    double target = 0) {

  DataFrame Sstar = untreated(psi*treat_modifier, time, event, treat, rx,
                              censor_time, recensor, autoswitch);

  NumericVector t_star = Sstar["t_star"];
  IntegerVector d_star = Sstar["d_star"];

  DataFrame data = DataFrame::create(
    Named("stratum") = stratum,
    Named("treat") = treat,
    Named("time") = t_star,
    Named("event") = d_star);

  DataFrame df = lrtest(data, "", "stratum", "treat", "time", "event", 0, 0);

  double result = as<double>(df["logRankZ"]) - target;
  return result;
}


// [[Rcpp::export]]
List rpsftmcpp(const DataFrame data,
               const StringVector& stratum = "",
               const std::string time = "time",
               const std::string event = "event",
               const std::string treat = "treat",
               const std::string rx = "rx",
               const std::string censor_time = "censor_time",
               const StringVector& base_cov = "",
               const double low_psi = -1,
               const double hi_psi = 1,
               const int n_eval_z = 100,
               const double treat_modifier = 1,
               const bool recensor = 1,
               const bool admin_recensor_only = 1,
               const bool autoswitch = 1,
               const bool gridsearch = 0,
               const double alpha = 0.05,
               const std::string ties = "efron",
               const double tol = 1.0e-6,
               const bool boot = 0,
               const int n_boot = 1000,
               const int seed = NA_INTEGER) {

  int i, j, k, l, n = data.nrow();
  int p = static_cast<int>(base_cov.size());
  if (p == 1 && (base_cov[0] == "" || base_cov[0] == "none")) p = 0;

  IntegerVector stratumn(n);
  int p_stratum = static_cast<int>(stratum.size());
  if (p_stratum == 1 && (stratum[0] == "" || stratum[0] == "none")) {
    stratumn.fill(1);
  } else {
    List out = bygroup(data, stratum);
    stratumn = out["index"];
  }

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
  if (TYPEOF(data[treat]) == LGLSXP || TYPEOF(data[treat]) == INTSXP) {
    IntegerVector treatv = data[treat];
    IntegerVector treatwi = unique(treatv);
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
    NumericVector treatwn = unique(treatv);
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
    StringVector treatwc = unique(treatv);
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

  // covariates for the Cox proportional hazards model
  // containing treat and base_cov
  StringVector covariates(p+1);
  NumericMatrix zn(n,p+1);
  covariates[0] = "treat";
  zn(_,0) = treatn;
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
    zn(_,j+1) = u;
  }

  if (low_psi >= hi_psi) {
    stop("low_psi must be less than hi_psi");
  }

  if (n_eval_z < 2) {
    stop("n_eval_z must be greater than or equal to 2");
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
  
  if (n_boot < 100) {
    stop("n_boot must be greater than or equal to 100");
  }

  DataFrame lr = lrtest(data, "", stratum, treat, time, event, 0, 0);
  double logRankPValue = as<double>(lr["logRankPValue"]);

  // evaluate the log-rank test statistic at each psi
  double step_psi = (hi_psi - low_psi)/(n_eval_z - 1);
  NumericVector psi(n_eval_z), Z(n_eval_z);
  for (i=0; i<n_eval_z; i++) {
    psi[i] = low_psi + i*step_psi;
    Z[i] = est_eqn(psi[i], stratumn, timen, eventn, treatn,
                   rxn, censor_timen, treat_modifier,
                   recensor, autoswitch, 0);
  }

  DataFrame eval_z = DataFrame::create(
    Named("psi") = psi,
    Named("Z") = Z);

  double zcrit = R::qnorm(1-alpha/2, 0, 1, 1, 0);

  k = -1;
  auto f = [&k, n, p, covariates, low_psi, hi_psi, n_eval_z, psi, 
            treat_modifier, recensor, autoswitch, gridsearch, 
            alpha, zcrit, ties, tol](
                IntegerVector stratumb, NumericVector timeb,
                IntegerVector eventb, IntegerVector treatb,
                NumericVector rxb, NumericVector censor_timeb,
                NumericMatrix zb)->List {
                  int i, j;

                  // obtain the estimate and confidence interval of psi
                  double psihat, psilower = 0, psiupper = 0;
                  String psi_CI_type;

                  if (gridsearch) {
                    NumericVector Z(n_eval_z);
                    for (i=0; i<n_eval_z; i++) {
                      Z[i] = est_eqn(psi[i], stratumb, timeb, eventb, 
                                     treatb, rxb, censor_timeb, 
                                     treat_modifier, recensor, 
                                     autoswitch, 0);
                    }

                    auto g = [psi, Z](double target)->double{
                      NumericVector Z1 = Z - target;
                      NumericVector Zsq = Z1*Z1;
                      return psi[which_min(Zsq)];
                    };

                    psihat = g(0);
                    psi_CI_type = "grid search";

                    if (k == -1) {
                      psilower = g(zcrit);
                      psiupper = g(-zcrit);
                    }
                  } else {
                    double target = 0;
                    auto g = [&target, stratumb, timeb, eventb, treatb, 
                              rxb, censor_timeb, treat_modifier, recensor, 
                              autoswitch](double x)->double {
                                return est_eqn(x, stratumb, timeb, eventb,
                                               treatb, rxb, censor_timeb,
                                               treat_modifier, recensor,
                                               autoswitch, target);
                              };

                    psihat = brent(g, low_psi, hi_psi, tol);
                    psi_CI_type = "root finding";

                    if (k == -1) {
                      target = zcrit;
                      psilower = brent(g, low_psi, psihat, tol);
                      target = -zcrit;
                      psiupper = brent(g, psihat, hi_psi, tol);
                    }
                  }

                  // construct the counter-factual untreated survival times
                  double psi1 = psihat*treat_modifier;
                  DataFrame Sstar = untreated(psi1, timeb, eventb, treatb,
                                              rxb, censor_timeb, recensor,
                                              autoswitch);

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
                  NumericVector rx1 = rxb[treatb == 1];
                  bool switch1 = is_false(all(rx1 == 1));
                  NumericVector timewb(n);
                  IntegerVector eventwb(n);
                  for (i=0; i<n; i++) {
                    if (treatb[i] == 1) {
                      if (switch1) {
                        // counter-factual survival time on treatment
                        double t_tilde, c_tilde;
                        t_tilde = timeb[i]*((1-rxb[i])*exp(-psi1) + rxb[i]);
                        if (recensor) {
                          c_tilde = std::min(censor_timeb[i],
                                             censor_timeb[i]*exp(-psi1));
                          timewb[i] = std::min(t_tilde, c_tilde);
                          eventwb[i] = c_tilde < t_tilde ? 0 : eventb[i];
                        } else {
                          timewb[i] = t_tilde;
                          eventwb[i] = eventb[i];
                        }
                      } else {
                        timewb[i] = timeb[i];
                        eventwb[i] = eventb[i];
                      }
                    } else {
                      timewb[i] = t_star[i];
                      eventwb[i] = d_star[i];
                    }
                  }

                  DataFrame data_outcome = DataFrame::create(
                    Named("stratum") = stratumb,
                    Named("time") = timewb,
                    Named("event") = eventwb,
                    Named("treat") = treatb);

                  for (j=0; j<p; j++) {
                    String zj = covariates[j+1];
                    NumericVector u = zb(_,j+1);
                    data_outcome.push_back(u, zj);
                  }

                  List fit_outcome = phregcpp(
                    data_outcome, "", "stratum", "time", "", "event", 
                    covariates, "", "", "", ties, 0, 0, 0, 0, 0, alpha);

                  DataFrame parest = DataFrame(fit_outcome["parest"]);
                  NumericVector beta = parest["beta"];
                  NumericVector z = parest["z"];
                  double hrhat = exp(beta[0]);
                  double pvalue = 2*(1 - R::pnorm(fabs(z[0]), 0, 1, 1, 0));
                  
                  List out;
                  if (k == -1) {
                    out = List::create(
                      Named("psihat") = psihat,
                      Named("psilower") = psilower,
                      Named("psiupper") = psiupper,
                      Named("psi_CI_type") = psi_CI_type,
                      Named("Sstar") = Sstar,
                      Named("kmstar") = kmstar,
                      Named("data_outcome") = data_outcome,
                      Named("fit_outcome") = fit_outcome,
                      Named("hrhat") = hrhat,
                      Named("pvalue") = pvalue);
                  } else {
                    out = List::create(
                      Named("psihat") = psihat,
                      Named("psi_CI_type") = psi_CI_type,
                      Named("hrhat") = hrhat,
                      Named("pvalue") = pvalue);
                  }

                  return out;
                };

  List out = f(stratumn, timen, eventn, treatn, rxn, censor_timen, zn);

  DataFrame Sstar = DataFrame(out["Sstar"]);
  DataFrame kmstar = DataFrame(out["kmstar"]);
  DataFrame data_outcome = DataFrame(out["data_outcome"]);
  List fit_outcome = out["fit_outcome"];
  double psihat = out["psihat"];
  double psilower = out["psilower"];
  double psiupper = out["psiupper"];
  String psi_CI_type = out["psi_CI_type"];
  double hrhat = out["hrhat"];
  double pvalue = out["pvalue"];

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
    NumericMatrix zb(n,p+1);

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
        for (l=0; l<p+1; l++) {
          zb(i,l) = zn(j,l);
        }
      }

      List out = f(stratumb, timeb, eventb, treatb, rxb, censor_timeb, zb);
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
    Named("low_psi") = low_psi,
    Named("hi_psi") = hi_psi,
    Named("n_eval_z") = n_eval_z,
    Named("treat_modifer") = treat_modifier,
    Named("recensor") = recensor,
    Named("admin_recensor_only") = admin_recensor_only,
    Named("autoswitch") = autoswitch,
    Named("gridsearch") =gridsearch,
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
    Named("eval_z") = eval_z,
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
