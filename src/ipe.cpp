#include "utilities.h"
#include "survival_analysis.h"

using namespace Rcpp;


List est_psi_ipe(
    const double psi,
    const int n,
    const int q,
    const int p,
    const IntegerVector& id,
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
  NumericVector init(1, NA_REAL);
  DataFrame df = unswitched(psi*treat_modifier, n, id, time, event, treat,
                            rx, censor_time, recensor, autoswitch);

  for (int j=0; j<q+p; j++) {
    String zj = covariates_aft[j+1];
    NumericVector u = zb_aft(_,j);
    df.push_back(u, zj);
  }

  List fit = liferegcpp(df, "", "", "t_star", "", "d_star",
                        covariates_aft, "", "", "", dist, init, 
                        0, 0, alpha, 50, 1.0e-9);
  
  DataFrame sumstat = DataFrame(fit["sumstat"]);
  bool fail = sumstat["fail"];

  DataFrame parest = DataFrame(fit["parest"]);
  NumericVector beta = parest["beta"];
  double psinew = -beta[1]/treat_modifier;
  
  List out = List::create(
    Named("data_aft") = df,
    Named("fit_aft") = fit,
    Named("psinew") = psinew,
    Named("fail") = fail
  );
  
  return out;
}


// [[Rcpp::export]]
List ipecpp(const DataFrame data,
            const std::string id = "id",
            const StringVector& stratum = "",
            const std::string time = "time",
            const std::string event = "event",
            const std::string treat = "treat",
            const std::string rx = "rx",
            const std::string censor_time = "censor_time",
            const StringVector& base_cov = "",
            const std::string aft_dist = "weibull",
            const bool strata_main_effect_only = 1,
            const double low_psi = -2,
            const double hi_psi = 2,
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

  int i, j, k, n = data.nrow();
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

  bool has_id = hasVariable(data, id);
  bool has_time = hasVariable(data, time);
  bool has_event = hasVariable(data, event);
  bool has_treat = hasVariable(data, treat);
  bool has_rx = hasVariable(data, rx);
  bool has_censor_time = hasVariable(data, censor_time);

  if (!has_id) {
    stop("data must contain the id variable");
  }
  
  IntegerVector idn(n);
  IntegerVector idwi;
  NumericVector idwn;
  StringVector idwc;
  if (TYPEOF(data[id]) == INTSXP) {
    IntegerVector idv = data[id];
    idwi = unique(idv);
    idwi.sort();
    idn = match(idv, idwi);
  } else if (TYPEOF(data[id]) == REALSXP) {
    NumericVector idv = data[id];
    idwn = unique(idv);
    idwn.sort();
    idn = match(idv, idwn);
  } else if (TYPEOF(data[id]) == STRSXP) {
    StringVector idv = data[id];
    idwc = unique(idv);
    idwc.sort();
    idn = match(idv, idwc);
  } else {
    stop("incorrect type for the id variable in the input data");
  }
  
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

  // covariates for the Cox proportional hazards model 
  // containing treat and base_cov
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
  
  if (low_psi >= hi_psi) {
    stop("low_psi must be less than hi_psi");
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
  double logRankPValue = lr["logRankPValue"];
  
  double zcrit = R::qnorm(1-alpha/2, 0, 1, 1, 0);

  k = -1;
  auto f = [&k, n, q, p, covariates, covariates_aft, dist, low_psi, hi_psi,
            treat_modifier, recensor, autoswitch, alpha, ties, tol](
                IntegerVector& idb, 
                IntegerVector& stratumb, NumericVector& timeb,
                IntegerVector& eventb, IntegerVector& treatb,
                NumericVector& rxb, NumericVector& censor_timeb,
                NumericMatrix& zb, NumericMatrix& zb_aft)->List {
                  int j;
                  bool fail = 0; // whether any model fails to converge
                  NumericVector init(1, NA_REAL);
                  
                  // estimate psi
                  auto g = [n, q, p, idb, timeb, eventb, treatb, rxb, 
                            censor_timeb, covariates_aft, zb_aft, 
                            dist, treat_modifier, recensor, 
                            autoswitch, alpha](double psi)->double{
                              List out_aft = est_psi_ipe(
                                psi, n, q, p, idb, timeb, eventb, treatb, 
                                rxb, censor_timeb, covariates_aft, zb_aft, 
                                dist, treat_modifier, recensor, autoswitch, 
                                alpha);
                              
                              double psinew = out_aft["psinew"];
                              return psinew - psi;
                            };

                  double psihat = brent(g, low_psi, hi_psi, tol);

                  // obtain the Kaplan-Meier estimates
                  DataFrame Sstar, kmstar, data_aft;
                  List fit_aft;
                  if (k == -1) {
                    // construct the counterfactual survival times
                    Sstar = untreated(
                      psihat*treat_modifier, idb, timeb, eventb, treatb,
                      rxb, censor_timeb, recensor, autoswitch);
                    
                    kmstar = kmest(Sstar, "", "treated", "t_star",
                                   "d_star", "log-log", 1-alpha, 1);
                    
                    Sstar.push_back(stratumb, "ustratum");
                    
                    for (j=0; j<p; j++) {
                      String zj = covariates[j+1];
                      NumericVector u = zb(_,j);
                      Sstar.push_back(u, zj);
                    }
                    
                    List out_aft = est_psi_ipe(
                      psihat, n, q, p, idb, timeb, eventb, treatb, rxb,
                      censor_timeb, covariates_aft, zb_aft, dist, 
                      treat_modifier, recensor, autoswitch, alpha);
                    
                    bool fail_aft = out_aft["fail"];
                    if (fail_aft) fail = 1;
                    
                    data_aft = DataFrame(out_aft["data_aft"]);
                    data_aft.push_back(stratumb, "ustratum");
                    
                    fit_aft = out_aft["fit_aft"];
                  }
                  
                  // run Cox model to obtain the hazard ratio estimate
                  DataFrame data_outcome = unswitched(
                    psihat*treat_modifier, n, idb, timeb, eventb, treatb,
                    rxb, censor_timeb, recensor, autoswitch);

                  data_outcome.push_back(stratumb, "ustratum");
                  
                  for (j=0; j<p; j++) {
                    String zj = covariates[j+1];
                    NumericVector u = zb(_,j);
                    data_outcome.push_back(u, zj);
                  }

                  List fit_outcome = phregcpp(
                    data_outcome, "", "ustratum", "t_star", "", "d_star", 
                    covariates, "", "", "", ties, init, 
                    0, 0, 0, 0, 0, alpha, 50, 1.0e-9);

                  DataFrame sumstat_cox = DataFrame(fit_outcome["sumstat"]);
                  bool fail_cox = sumstat_cox["fail"];
                  if (fail_cox == 1) fail = 1;
                  
                  DataFrame parest = DataFrame(fit_outcome["parest"]);
                  NumericVector beta = parest["beta"];
                  NumericVector pval = parest["p"];
                  double hrhat = exp(beta[0]/treat_modifier);
                  double pvalue = pval[0];
                  
                  List out;
                  if (k == -1) {
                    out = List::create(
                      Named("Sstar") = Sstar,
                      Named("kmstar") = kmstar,
                      Named("data_aft") = data_aft,
                      Named("fit_aft") = fit_aft,
                      Named("data_outcome") = data_outcome,
                      Named("fit_outcome") = fit_outcome,
                      Named("psihat") = psihat,
                      Named("hrhat") = hrhat,
                      Named("pvalue") = pvalue,
                      Named("fail") = fail);
                  } else {
                    out = List::create(
                      Named("psihat") = psihat,
                      Named("hrhat") = hrhat,
                      Named("pvalue") = pvalue,
                      Named("fail") = fail);
                  }

                  return out;
                };

  List out = f(idn, stratumn, timen, eventn, treatn, rxn, censor_timen, 
               zn, zn_aft);

  DataFrame Sstar = DataFrame(out["Sstar"]);
  DataFrame kmstar = DataFrame(out["kmstar"]);
  DataFrame data_aft = DataFrame(out["data_aft"]);
  List fit_aft = out["fit_aft"];
  DataFrame data_outcome = DataFrame(out["data_outcome"]);
  List fit_outcome = out["fit_outcome"];
  
  IntegerVector uid = Sstar["uid"];
  if (TYPEOF(data[id]) == INTSXP) {
    Sstar.push_front(idwi[uid-1], id);
  } else if (TYPEOF(data[id]) == REALSXP) {
    Sstar.push_front(idwn[uid-1], id);
  } else if (TYPEOF(data[id]) == STRSXP) {
    Sstar.push_front(idwc[uid-1], id);
  }
  
  uid = data_aft["uid"];
  if (TYPEOF(data[id]) == INTSXP) {
    data_aft.push_front(idwi[uid-1], id);
  } else if (TYPEOF(data[id]) == REALSXP) {
    data_aft.push_front(idwn[uid-1], id);
  } else if (TYPEOF(data[id]) == STRSXP) {
    data_aft.push_front(idwc[uid-1], id);
  }
  
  uid = data_outcome["uid"];
  if (TYPEOF(data[id]) == INTSXP) {
    data_outcome.push_front(idwi[uid-1], id);
  } else if (TYPEOF(data[id]) == REALSXP) {
    data_outcome.push_front(idwn[uid-1], id);
  } else if (TYPEOF(data[id]) == STRSXP) {
    data_outcome.push_front(idwc[uid-1], id);
  }
  
  IntegerVector treated = Sstar["treated"];
  if (TYPEOF(data[treat]) == LGLSXP || TYPEOF(data[treat]) == INTSXP) {
    Sstar.push_back(treatwi[1-treated], treat);
  } else if (TYPEOF(data[treat]) == REALSXP) {
    Sstar.push_back(treatwn[1-treated], treat);
  } else if (TYPEOF(data[treat]) == STRSXP) {
    Sstar.push_back(treatwc[1-treated], treat);
  }
  
  treated = kmstar["treated"];
  if (TYPEOF(data[treat]) == LGLSXP || TYPEOF(data[treat]) == INTSXP) {
    kmstar.push_back(treatwi[1-treated], treat);
  } else if (TYPEOF(data[treat]) == REALSXP) {
    kmstar.push_back(treatwn[1-treated], treat);
  } else if (TYPEOF(data[treat]) == STRSXP) {
    kmstar.push_back(treatwc[1-treated], treat);
  }
  
  treated = data_aft["treated"];
  if (TYPEOF(data[treat]) == LGLSXP || TYPEOF(data[treat]) == INTSXP) {
    data_aft.push_back(treatwi[1-treated], treat);
  } else if (TYPEOF(data[treat]) == REALSXP) {
    data_aft.push_back(treatwn[1-treated], treat);
  } else if (TYPEOF(data[treat]) == STRSXP) {
    data_aft.push_back(treatwc[1-treated], treat);
  }
  
  treated = data_outcome["treated"];
  if (TYPEOF(data[treat]) == LGLSXP || TYPEOF(data[treat]) == INTSXP) {
    data_outcome.push_back(treatwi[1-treated], treat);
  } else if (TYPEOF(data[treat]) == REALSXP) {
    data_outcome.push_back(treatwn[1-treated], treat);
  } else if (TYPEOF(data[treat]) == STRSXP) {
    data_outcome.push_back(treatwc[1-treated], treat);
  }
  
  
  if (has_stratum) {
    IntegerVector ustratum = Sstar["ustratum"];
    for (i=0; i<p_stratum; i++) {
      String s = stratum[i];
      if (TYPEOF(data[s]) == INTSXP) {
        IntegerVector stratumwi = u_stratum[s];
        Sstar.push_back(stratumwi[ustratum-1], s);
      } else if (TYPEOF(data[s]) == REALSXP) {
        NumericVector stratumwn = u_stratum[s];
        Sstar.push_back(stratumwn[ustratum-1], s);
      } else if (TYPEOF(data[s]) == STRSXP) {
        StringVector stratumwc = u_stratum[s];
        Sstar.push_back(stratumwc[ustratum-1], s);
      }
    }
    
    ustratum = data_aft["ustratum"];
    for (i=0; i<p_stratum; i++) {
      String s = stratum[i];
      if (TYPEOF(data[s]) == INTSXP) {
        IntegerVector stratumwi = u_stratum[s];
        data_aft.push_back(stratumwi[ustratum-1], s);
      } else if (TYPEOF(data[s]) == REALSXP) {
        NumericVector stratumwn = u_stratum[s];
        data_aft.push_back(stratumwn[ustratum-1], s);
      } else if (TYPEOF(data[s]) == STRSXP) {
        StringVector stratumwc = u_stratum[s];
        data_aft.push_back(stratumwc[ustratum-1], s);
      }
    }
    
    ustratum = data_outcome["ustratum"];
    for (i=0; i<p_stratum; i++) {
      String s = stratum[i];
      if (TYPEOF(data[s]) == INTSXP) {
        IntegerVector stratumwi = u_stratum[s];
        data_outcome.push_back(stratumwi[ustratum-1], s);
      } else if (TYPEOF(data[s]) == REALSXP) {
        NumericVector stratumwn = u_stratum[s];
        data_outcome.push_back(stratumwn[ustratum-1], s);
      } else if (TYPEOF(data[s]) == STRSXP) {
        StringVector stratumwc = u_stratum[s];
        data_outcome.push_back(stratumwc[ustratum-1], s);
      }
    }
  }
  
  
  double psihat = out["psihat"];
  double zipe = R::qnorm(logRankPValue, 0, 1, 1, 0);
  double sepsi = psihat/zipe;
  double psilower = psihat - zcrit*sepsi;
  double psiupper = psihat + zcrit*sepsi;
  String psi_CI_type = "log-rank p-value";
  double hrhat = out["hrhat"];
  double pvalue = out["pvalue"];
  bool fail = out["fail"];
  
  // construct the confidence interval for HR
  double hrlower, hrupper;
  NumericVector hrhats(n_boot), psihats(n_boot);
  LogicalVector fails(n_boot);
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

    IntegerVector idb(n), stratumb(n), treatb(n), eventb(n);
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

    idn = idn[order];
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

        idb[i] = idn[j];
        stratumb[i] = stratumn[j];
        timeb[i] = timen[j];
        eventb[i] = eventn[j];
        treatb[i] = treatn[j];
        rxb[i] = rxn[j];
        censor_timeb[i] = censor_timen[j];
        zb(i,_) = zn(j,_);
        zb_aft(i,_) = zn_aft(j,_);
      }

      List out = f(idb, stratumb, timeb, eventb, treatb, rxb, censor_timeb, 
                   zb, zb_aft);
      
      fails[k] = out["fail"];
      hrhats[k] = out["hrhat"];
      psihats[k] = out["psihat"];
    }

    // obtain bootstrap confidence interval for HR
    double loghr = log(hrhat);
    LogicalVector ok = 1 - fails;
    int n_ok = sum(ok);
    NumericVector loghrs = log(hrhats[ok]);
    double sdloghr = sd(loghrs);
    double tcrit = R::qt(1-alpha/2, n_ok-1, 1, 0);
    hrlower = exp(loghr - tcrit*sdloghr);
    hrupper = exp(loghr + tcrit*sdloghr);
    hr_CI_type = "bootstrap";
    pvalue = 2*(1 - R::pt(fabs(loghr/sdloghr), n_ok-1, 1, 0));
    
    // obtain bootstrap confidence interval for psi
    NumericVector psihats1 = psihats[ok];
    double sdpsi = sd(psihats1);
    psilower = psihat - tcrit*sdpsi;
    psiupper = psihat + tcrit*sdpsi;
    psi_CI_type = "bootstrap";
  }

  List settings = List::create(
    Named("aft_dist") = aft_dist,
    Named("strata_main_effect_only") = strata_main_effect_only,
    Named("low_psi") = low_psi,
    Named("hi_psi") = hi_psi,
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
    Named("data_aft") = data_aft,
    Named("fit_aft") = fit_aft,
    Named("data_outcome") = data_outcome,
    Named("fit_outcome") = fit_outcome,
    Named("fail") = fail,
    Named("settings") = settings);

  if (boot) {
    result.push_back(fails, "fail_boots");
    result.push_back(hrhats, "hr_boots");
    result.push_back(psihats, "psi_boots");
  }
  
  return result;
}
