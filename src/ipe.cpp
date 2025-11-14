#include "utilities.h"
#include "survival_analysis.h"

using namespace Rcpp;


List est_psi_ipe(
    const double psi,
    const int n,
    const int q,
    const int p,
    const IntegerVector& idb,
    const NumericVector& timeb,
    const IntegerVector& eventb,
    const IntegerVector& treatb,
    const NumericVector& rxb,
    const NumericVector& censor_timeb,
    const StringVector& covariates_aft,
    const NumericMatrix& z_aftb,
    const std::string dist,
    const double treat_modifier,
    const bool recensor,
    const bool autoswitch,
    const double alpha) {
  
  NumericVector init(1, NA_REAL);
  DataFrame df = unswitched(psi*treat_modifier, n, idb, timeb, eventb, treatb,
                            rxb, censor_timeb, recensor, autoswitch);
  
  for (int j=0; j<q+p; ++j) {
    String zj = covariates_aft[j+1];
    NumericVector u = z_aftb(_,j);
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
            const bool strata_main_effect_only = true,
            const double low_psi = -2,
            const double hi_psi = 2,
            const double treat_modifier = 1,
            const bool recensor = true,
            const bool admin_recensor_only = true,
            const bool autoswitch = true,
            const std::string root_finding = "brent",
            const double alpha = 0.05,
            const std::string ties = "efron",
            const double tol = 1.0e-6,
            const bool boot = false,
            const int n_boot = 1000,
            const int seed = NA_INTEGER) {
  
  int k, n = data.nrow();
  int p = static_cast<int>(base_cov.size());
  if (p == 1 && (base_cov[0] == "" || base_cov[0] == "none")) p = 0;
  
  int p_stratum = static_cast<int>(stratum.size());
  
  bool has_stratum;
  IntegerVector stratumn(n);
  DataFrame u_stratum;
  IntegerVector d(p_stratum);
  IntegerMatrix stratan(n,p_stratum);
  List levels(p_stratum);
  if (p_stratum == 1 && (stratum[0] == "" || stratum[0] == "none")) {
    has_stratum = 0;
    stratumn.fill(0);
    d[0] = 1;
    stratan(_,0) = stratumn;
    levels[0] = 1;
  } else {
    List out = bygroup(data, stratum);
    has_stratum = 1;
    stratumn = out["index"];
    u_stratum = DataFrame(out["lookup"]);
    d = out["nlevels"];
    stratan = as<IntegerMatrix>(out["indices"]);
    levels = out["lookups"];
  }
  
  IntegerVector stratumn_unique = unique(stratumn);
  int nstrata = static_cast<int>(stratumn_unique.size());
  
  bool has_id = hasVariable(data, id);
  bool has_time = hasVariable(data, time);
  bool has_event = hasVariable(data, event);
  bool has_treat = hasVariable(data, treat);
  bool has_rx = hasVariable(data, rx);
  bool has_censor_time = hasVariable(data, censor_time);
  
  
  if (!has_id) stop("data must contain the id variable");
  
  IntegerVector idn(n);
  IntegerVector idwi;
  NumericVector idwn;
  StringVector idwc;
  
  SEXP col_id = data[id];
  SEXPTYPE type_id = TYPEOF(col_id);
  
  if (type_id == INTSXP) {
    IntegerVector idv = col_id;
    idwi = unique(idv);
    idwi.sort();
    idn = match(idv, idwi) - 1;
  } else if (type_id == REALSXP) {
    NumericVector idv = col_id;
    idwn = unique(idv);
    idwn.sort();
    idn = match(idv, idwn) - 1;
  } else if (type_id == STRSXP) {
    StringVector idv = col_id;
    idwc = unique(idv);
    idwc.sort();
    idn = match(idv, idwc) - 1;
  } else {
    stop("incorrect type for the id variable in the input data");
  }
  
  
  if (!has_time) stop("data must contain the time variable");
  
  SEXP col_time = data[time];
  SEXPTYPE type_time = TYPEOF(col_time);
  
  if (type_time != INTSXP && type_time != REALSXP) {
    stop("time must take numeric values");
  }
  
  NumericVector timenz = col_time;
  NumericVector timen = clone(timenz);
  if (is_true(any(timen <= 0.0))) {
    stop("time must be positive");
  }
  
  
  if (!has_event) stop("data must contain the event variable");
  
  SEXP col_event = data[event];
  SEXPTYPE type_event = TYPEOF(col_event);
  
  IntegerVector eventn(n);
  if (type_event == LGLSXP || type_event == INTSXP) {
    IntegerVector eventnz = col_event;
    if (is_true(any((eventnz != 1) & (eventnz != 0)))) {
      stop("event must be 1 or 0 for each subject");
    } else {
      eventn = clone(eventnz);
    }
  } else if (type_event == REALSXP) {
    NumericVector eventnz = col_event;
    if (is_true(any((eventnz != 1) & (eventnz != 0)))) {
      stop("event must be 1 or 0 for each subject");
    } else {
      NumericVector eventnz2 = clone(eventnz);
      eventn = as<IntegerVector>(eventnz2);
    }
  } else {
    stop("event must take logical, integer, or real values");
  }
  
  if (is_true(all(eventn == 0))) {
    stop("at least 1 event is needed");
  }
  
  
  if (!has_treat) stop("data must contain the treat variable");
  
  IntegerVector treatn(n);
  IntegerVector treatwi;
  NumericVector treatwn;
  StringVector treatwc;
  
  SEXP col_treat = data[treat];
  SEXPTYPE type_treat = TYPEOF(col_treat);
  
  if (type_treat == LGLSXP || type_treat == INTSXP) {
    IntegerVector treatv = col_treat;
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
  } else if (type_treat == REALSXP) {
    NumericVector treatv = col_treat;
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
  } else if (type_treat == STRSXP) {
    StringVector treatv = col_treat;
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
  
  
  if (!has_rx) stop("data must contain the rx variable");
  
  SEXP col_rx = data[rx];
  SEXPTYPE type_rx = TYPEOF(col_rx);
  
  if (type_rx != INTSXP && type_rx != REALSXP) {
    stop("rx must take numeric values");
  }
  
  NumericVector rxnz = col_rx;
  NumericVector rxn = clone(rxnz);
  if (is_true(any((rxn < 0.0) | (rxn > 1.0)))) {
    stop("rx must take values between 0 and 1");
  }
  
  
  if (!has_censor_time) stop("data must contain the censor_time variable");

  SEXP col_censor_time = data[censor_time];
  SEXPTYPE type_censor_time = TYPEOF(col_censor_time);
  
  if (type_censor_time != INTSXP && type_censor_time != REALSXP) {
    stop("censor_time must take numeric values");
  }
  
  NumericVector censor_timenz = col_censor_time;
  NumericVector censor_timen = clone(censor_timenz);
  if (is_true(any(censor_timen < timen))) {
    stop("censor_time must be greater than or equal to time");
  }
  
  if (!admin_recensor_only) {
    for (int i=0; i<n; ++i) {
      if (eventn[i] == 0) { // use the actual censoring time for dropouts
        censor_timen[i] = timen[i];
      }
    }
  }
  

  int q; // number of columns corresponding to the strata effects
  if (strata_main_effect_only) {
    q = sum(d - 1);
  } else {
    q = nstrata - 1;
  }
  
  // covariates for the accelerated failure time model
  // including treat, stratum, and base_cov
  StringVector covariates_aft(q+p+1);
  NumericMatrix z_aftn(n,q+p);
  covariates_aft[0] = "treated";
  if (strata_main_effect_only) {
    k = 0;
    for (int i=0; i<p_stratum; ++i) {
      SEXP col_level = levels[i];
      SEXPTYPE type_level = TYPEOF(col_level);
      
      int di = d[i]-1;
      for (int j=0; j<di; ++j) {
        covariates_aft[k+j+1] = as<std::string>(stratum[i]);
        
        if (type_level == STRSXP) {
          StringVector u = col_level;
          std::string label = sanitize(as<std::string>(u[j]));
          covariates_aft[k+j+1] += label;
        } else if (type_level == REALSXP) {
          NumericVector u = col_level;
          covariates_aft[k+j+1] += std::to_string(u[j]);
        } else if (type_level == INTSXP || type_level == LGLSXP) {
          IntegerVector u = col_level;
          covariates_aft[k+j+1] += std::to_string(u[j]);
        }
        
        z_aftn(_,k+j) = (stratan(_,i) == j);
      }
      k += di;
    }
  } else {
    for (int j=0; j<nstrata-1; ++j) {
      // locate the first observation in the stratum
      int first_k = 0;
      for (; first_k<n; ++first_k) {
        if (stratumn[first_k] == j) break;
      }
      covariates_aft[j+1] = "";
      for (int i=0; i<p_stratum; ++i) {
        SEXP col_level = levels[i];
        SEXPTYPE type_level = TYPEOF(col_level);
        
        IntegerVector q_col = stratan(_,i);
        int l = q_col[first_k];
        
        covariates_aft[j+1] += as<std::string>(stratum[i]);
        
        if (type_level == STRSXP) {
          StringVector u = col_level;
          std::string label = sanitize(as<std::string>(u[l]));
          covariates_aft[j+1] += label;
        } else if (type_level == REALSXP) {
          NumericVector u = col_level;
          covariates_aft[j+1] += std::to_string(u[l]);
        } else if (type_level == INTSXP || type_level == LGLSXP) {
          IntegerVector u = col_level;
          covariates_aft[j+1] += std::to_string(u[l]);
        }
        
        if (i < p_stratum-1) {
          covariates_aft[j+1] += ".";
        }
      }
      z_aftn(_,j) = (stratumn == j);
    }
  }
  
  // covariates for the Cox model containing treat and base_cov
  StringVector covariates(p+1);
  NumericMatrix zn(n,p);
  covariates[0] = "treated";
  for (int j=0; j<p; ++j) {
    String zj = base_cov[j];
    
    if (!hasVariable(data, zj)) {
      stop("data must contain the variables in base_cov");
    }
    
    if (zj == treat) {
      stop("treat should be excluded from base_cov");
    }
    
    NumericVector u = data[zj];
    
    covariates_aft[q+j+1] = zj;
    z_aftn(_,q+j) = u;
    
    covariates[j+1] = zj;
    zn(_,j) = u;
  }
  

  std::string dist = aft_dist;
  std::for_each(dist.begin(), dist.end(), [](char & c) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });
  
  if (dist == "log-logistic" || dist == "llogistic") {
    dist = "loglogistic";
  } else if  (dist == "log-normal" || dist == "lnormal") {
    dist = "lognormal";
  }
  
  if (!(dist == "exponential" || dist == "weibull" ||
      dist == "lognormal" || dist == "loglogistic")) {
    stop("aft_dist must be exponential, weibull, lognormal, or loglogistic");
  }
  
  if (low_psi >= hi_psi) {
    stop("low_psi must be less than hi_psi");
  }
  
  if (treat_modifier <= 0.0) {
    stop("treat_modifier must be positive");
  }
  
  std::string rooting = root_finding;
  std::for_each(rooting.begin(), rooting.end(), [](char & c) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });
  
  if (rooting == "uniroot" || rooting.find("br", 0) == 0) {
    rooting = "brent";
  } else if (rooting.find("bi", 0) == 0) {
    rooting = "bisection";
  }
  
  if (!(rooting == "brent" || rooting == "bisection")) {
    stop("root_finding must be brent or bisection");
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
  
  
  // exclude observations with missing values
  LogicalVector sub(n,1);
  for (int i=0; i<n; ++i) {
    if (idn[i] == NA_INTEGER || stratumn[i] == NA_INTEGER || 
        std::isnan(timen[i]) || eventn[i] == NA_INTEGER || 
        treatn[i] == NA_INTEGER || std::isnan(rxn[i]) ||
        std::isnan(censor_timen[i])) {
      sub[i] = 0;
    }
    for (int j=0; j<p; ++j) {
      if (std::isnan(zn(i,j))) sub[i] = 0;
    }
  }
  
  IntegerVector order = which(sub);
  idn = idn[order];
  stratumn = stratumn[order];
  timen = timen[order];
  eventn = eventn[order];
  treatn = treatn[order];
  rxn = rxn[order];
  censor_timen = censor_timen[order];
  if (p > 0) zn = subset_matrix_by_row(zn, order);
  z_aftn = subset_matrix_by_row(z_aftn, order);
  n = sum(sub);
  if (n == 0) stop("no observations left after removing missing values");
  
  
  // summarize number of deaths and switches by treatment arm
  IntegerVector treat_out = IntegerVector::create(0, 1);
  NumericVector n_total(2);
  NumericVector n_event(2);
  NumericVector n_switch(2);
  for (int i = 0; i < n; ++i) {
    int g = treatn[i];
    n_total[g]++;
    if (eventn[i] == 1) n_event[g]++;
    if ((treatn[i] == 0 && rxn[i] > 0) || (treatn[i] == 1 && rxn[i] < 1)) {
      n_switch[g]++;
    }
  }
  
  // Compute percentages
  NumericVector pct_event(2);
  NumericVector pct_switch(2);
  for (int g = 0; g < 2; g++) {
    pct_event[g] = 100.0 * n_event[g] / n_total[g];
    pct_switch[g] = 100.0 * n_switch[g] / n_total[g];
  }
  
  // Combine count and percentage
  List event_summary = List::create(
    _["treated"] = treat_out,
    _["n"] = n_total,
    _["event_n"] = n_event,
    _["event_pct"] = pct_event,
    _["switch_n"] = n_switch,
    _["switch_pct"] = pct_switch
  );
  
  
  // ITT analysis log-rank test
  DataFrame lr = lrtest(data, "", stratum, treat, time, "", event, "", 0,0,0);
  double logRankZ = lr["logRankZ"];
  double logRankPValue = lr["logRankPValue"];
  double zcrit = R::qnorm(1-alpha/2, 0, 1, 1, 0);
  
  k = -1; // indicate the observed data
  auto f = [&k, n, q, p, covariates, covariates_aft, dist, low_psi, hi_psi,
            treat_modifier, recensor, autoswitch, rooting, alpha, ties, tol](
                IntegerVector& idb, IntegerVector& stratumb, 
                NumericVector& timeb, IntegerVector& eventb, 
                IntegerVector& treatb, NumericVector& rxb, 
                NumericVector& censor_timeb,
                NumericMatrix& zb, NumericMatrix& z_aftb)->List {
                  bool fail = false; // whether any model fails to converge
                  NumericVector init(1, NA_REAL);
                  
                  // estimate psi
                  auto g = [n, q, p, idb, timeb, eventb, treatb, rxb, 
                            censor_timeb, covariates_aft, z_aftb, 
                            dist, treat_modifier, recensor, 
                            autoswitch, alpha](double psi)->double{
                              List out_aft = est_psi_ipe(
                                psi, n, q, p, idb, timeb, eventb, treatb, 
                                rxb, censor_timeb, covariates_aft, z_aftb, 
                                dist, treat_modifier, recensor, autoswitch, 
                                alpha);
                              
                              bool fail = out_aft["fail"];
                              if (!fail) {
                                double psinew = out_aft["psinew"];
                                return psinew - psi;
                              } else {
                                return NA_REAL;
                              }
                            };
                  
                  double psihat = NA_REAL;
                  double psilo = getpsiend(g, 1, low_psi);
                  double psihi = getpsiend(g, 0, hi_psi);
                  if (!std::isnan(psilo) && !std::isnan(psihi)) {
                    if (rooting == "brent") {
                      psihat = brent(g, psilo, psihi, tol);  
                    } else {
                      psihat = bisect(g, psilo, psihi, tol);
                    }
                  }
                  
                  // obtain the Kaplan-Meier estimates
                  List Sstar, kmstar, data_aft, data_outcome;
                  List fit_aft, fit_outcome;
                  NumericVector res_aft;
                  double hrhat = NA_REAL, pvalue = NA_REAL;
                  List km_outcome, lr_outcome;
                  
                  bool psimissing = std::isnan(psihat);
                  
                  if (!psimissing) {
                    if (k == -1) {
                      // construct the counterfactual survival times
                      Sstar = untreated(
                        psihat*treat_modifier, idb, timeb, eventb, treatb,
                        rxb, censor_timeb, recensor, autoswitch);
                      
                      kmstar = kmest(Sstar, "", "treated", "t_star", "",
                                     "d_star", "", "log-log", 1-alpha, 1);
                      
                      Sstar.push_back(stratumb, "ustratum");
                      
                      for (int j=0; j<p; ++j) {
                        String zj = covariates[j+1];
                        NumericVector u = zb(_,j);
                        Sstar.push_back(u, zj);
                      }
                      
                      List out_aft = est_psi_ipe(
                        psihat, n, q, p, idb, timeb, eventb, treatb, rxb,
                        censor_timeb, covariates_aft, z_aftb, dist, 
                        treat_modifier, recensor, autoswitch, alpha);
                      
                      bool fail_aft = out_aft["fail"];
                      if (fail_aft) fail = true;
                      
                      data_aft = out_aft["data_aft"];
                      data_aft.push_back(stratumb, "ustratum");
                      
                      fit_aft = out_aft["fit_aft"];
                      
                      DataFrame parest = DataFrame(fit_aft["parest"]);
                      NumericVector beta = parest["beta"];
                      NumericMatrix vbeta(q+p+1, q+p+1);
                      NumericMatrix rr = residuals_liferegcpp(
                        beta, vbeta, data_aft, "", "t_star", "", "d_star",
                        covariates_aft, "", "", "", dist, "deviance", 0, 0);
                      res_aft = rr(_,0);
                    }
                    
                    // run Cox model to obtain the hazard ratio estimate
                    data_outcome = unswitched(
                      psihat*treat_modifier, n, idb, timeb, eventb, treatb,
                      rxb, censor_timeb, recensor, autoswitch);
                    
                    data_outcome.push_back(stratumb, "ustratum");
                    
                    for (int j=0; j<p; ++j) {
                      String zj = covariates[j+1];
                      NumericVector u = zb(_,j);
                      data_outcome.push_back(u, zj);
                    }
                    
                    // generate KM estimate and log-rank test
                    if (k == -1) {
                      km_outcome = kmest(
                        data_outcome, "", "treated", "t_star", "", 
                        "d_star", "", "log-log", 1-alpha, 1);
                      
                      lr_outcome = lrtest(
                        data_outcome, "", "ustratum", "treated", "t_star", 
                        "", "d_star", "", 0, 0, 0);
                    }
                    
                    // fit the outcome model
                    fit_outcome = phregcpp(
                      data_outcome, "", "ustratum", "t_star", "", "d_star", 
                      covariates, "", "", "", ties, init, 
                      0, 0, 0, 0, 0, alpha, 50, 1.0e-9);
                    
                    DataFrame sumstat = DataFrame(fit_outcome["sumstat"]);
                    bool fail_cox = sumstat["fail"];
                    if (fail_cox) fail = true;
                    
                    DataFrame parest = DataFrame(fit_outcome["parest"]);
                    NumericVector beta = parest["beta"];
                    NumericVector pval = parest["p"];
                    hrhat = exp(beta[0]);
                    pvalue = pval[0];
                  }
                  
                  List out;
                  if (k == -1) {
                    out = List::create(
                      Named("Sstar") = Sstar,
                      Named("kmstar") = kmstar,
                      Named("data_aft") = data_aft,
                      Named("fit_aft") = fit_aft,
                      Named("res_aft") = res_aft,
                      Named("data_outcome") = data_outcome,
                      Named("km_outcome") = km_outcome,
                      Named("lr_outcome") = lr_outcome,
                      Named("fit_outcome") = fit_outcome,
                      Named("psihat") = psihat,
                      Named("hrhat") = hrhat,
                      Named("pvalue") = pvalue,
                      Named("fail") = fail,
                      Named("psimissing") = psimissing);
                  } else {
                    out = List::create(
                      Named("psihat") = psihat,
                      Named("hrhat") = hrhat,
                      Named("pvalue") = pvalue,
                      Named("fail") = fail,
                      Named("psimissing") = psimissing);
                  }
                  
                  return out;
                };
  
  List out = f(idn, stratumn, timen, eventn, treatn, rxn, censor_timen, 
               zn, z_aftn);
  
  List Sstar = out["Sstar"];
  List kmstar = out["kmstar"];
  List data_aft = out["data_aft"];
  List fit_aft = out["fit_aft"];
  NumericVector res_aft = out["res_aft"];
  List data_outcome = out["data_outcome"];
  List km_outcome = out["km_outcome"];
  List lr_outcome = out["lr_outcome"];
  List fit_outcome = out["fit_outcome"];
  
  double psihat = out["psihat"];
  double sepsi = logRankZ != 0 ? psihat/logRankZ : NA_REAL;
  double psilower = psihat - zcrit*sepsi;
  double psiupper = psihat + zcrit*sepsi;
  String psi_CI_type = "log-rank p-value";
  double hrhat = out["hrhat"];
  double pvalue = out["pvalue"];
  bool fail = out["fail"];
  bool psimissing = out["psimissing"];
  
  double hrlower = NA_REAL, hrupper = NA_REAL;
  NumericVector hrhats(n_boot), psihats(n_boot);
  LogicalVector fails(n_boot);
  List fail_boots_data;
  String hr_CI_type;
  
  if (!psimissing) {
    IntegerVector uid = Sstar["uid"];
    if (type_id == INTSXP) {
      Sstar.push_front(idwi[uid], id);
    } else if (type_id == REALSXP) {
      Sstar.push_front(idwn[uid], id);
    } else if (type_id == STRSXP) {
      Sstar.push_front(idwc[uid], id);
    }
    
    uid = data_aft["uid"];
    if (type_id == INTSXP) {
      data_aft.push_front(idwi[uid], id);
    } else if (type_id == REALSXP) {
      data_aft.push_front(idwn[uid], id);
    } else if (type_id == STRSXP) {
      data_aft.push_front(idwc[uid], id);
    }
    
    uid = data_outcome["uid"];
    if (type_id == INTSXP) {
      data_outcome.push_front(idwi[uid], id);
    } else if (type_id == REALSXP) {
      data_outcome.push_front(idwn[uid], id);
    } else if (type_id == STRSXP) {
      data_outcome.push_front(idwc[uid], id);
    }
    
    
    IntegerVector treated = event_summary["treated"];
    if (type_treat == LGLSXP || type_treat == INTSXP) {
      event_summary.push_back(treatwi[1-treated], treat);
    } else if (type_treat == REALSXP) {
      event_summary.push_back(treatwn[1-treated], treat);
    } else if (type_treat == STRSXP) {
      event_summary.push_back(treatwc[1-treated], treat);
    }
    
    treated = Sstar["treated"];
    if (type_treat == LGLSXP || type_treat == INTSXP) {
      Sstar.push_back(treatwi[1-treated], treat);
    } else if (type_treat == REALSXP) {
      Sstar.push_back(treatwn[1-treated], treat);
    } else if (type_treat == STRSXP) {
      Sstar.push_back(treatwc[1-treated], treat);
    }
    
    treated = kmstar["treated"];
    if (type_treat == LGLSXP || type_treat == INTSXP) {
      kmstar.push_back(treatwi[1-treated], treat);
    } else if (type_treat == REALSXP) {
      kmstar.push_back(treatwn[1-treated], treat);
    } else if (type_treat == STRSXP) {
      kmstar.push_back(treatwc[1-treated], treat);
    }
    
    treated = data_aft["treated"];
    if (type_treat == LGLSXP || type_treat == INTSXP) {
      data_aft.push_back(treatwi[1-treated], treat);
    } else if (type_treat == REALSXP) {
      data_aft.push_back(treatwn[1-treated], treat);
    } else if (type_treat == STRSXP) {
      data_aft.push_back(treatwc[1-treated], treat);
    }
    
    treated = data_outcome["treated"];
    if (type_treat == LGLSXP || type_treat == INTSXP) {
      data_outcome.push_back(treatwi[1-treated], treat);
    } else if (type_treat == REALSXP) {
      data_outcome.push_back(treatwn[1-treated], treat);
    } else if (type_treat == STRSXP) {
      data_outcome.push_back(treatwc[1-treated], treat);
    }
    
    treated = km_outcome["treated"];
    if (type_treat == LGLSXP || type_treat == INTSXP) {
      km_outcome.push_back(treatwi[1-treated], treat);
    } else if (type_treat == REALSXP) {
      km_outcome.push_back(treatwn[1-treated], treat);
    } else if (type_treat == STRSXP) {
      km_outcome.push_back(treatwc[1-treated], treat);
    }
    
    
    if (has_stratum) {
      IntegerVector ustratum = Sstar["ustratum"];
      for (int i=0; i<p_stratum; ++i) {
        std::string s = as<std::string>(stratum[i]);
        SEXP col_stratum = u_stratum[s];
        SEXPTYPE type_stratum = TYPEOF(col_stratum);
        if (type_stratum == INTSXP) {
          IntegerVector v = col_stratum;
          Sstar.push_back(v[ustratum], s);
        } else if (type_stratum == REALSXP) {
          NumericVector v = col_stratum;
          Sstar.push_back(v[ustratum], s);
        } else if (type_stratum == STRSXP) {
          StringVector v = col_stratum;
          Sstar.push_back(v[ustratum], s);
        }
      }
      
      ustratum = data_aft["ustratum"];
      for (int i=0; i<p_stratum; ++i) {
        std::string s = as<std::string>(stratum[i]);
        SEXP col_stratum = u_stratum[s];
        SEXPTYPE type_stratum = TYPEOF(col_stratum);
        if (type_stratum == INTSXP) {
          IntegerVector v = col_stratum;
          data_aft.push_back(v[ustratum], s);
        } else if (type_stratum == REALSXP) {
          NumericVector v = col_stratum;
          data_aft.push_back(v[ustratum], s);
        } else if (type_stratum == STRSXP) {
          StringVector v = col_stratum;
          data_aft.push_back(v[ustratum], s);
        }
      }
      
      ustratum = data_outcome["ustratum"];
      for (int i=0; i<p_stratum; ++i) {
        std::string s = as<std::string>(stratum[i]);
        SEXP col_stratum = u_stratum[s];
        SEXPTYPE type_stratum = TYPEOF(col_stratum);
        if (type_stratum == INTSXP) {
          IntegerVector v = col_stratum;
          data_outcome.push_back(v[ustratum], s);
        } else if (type_stratum == REALSXP) {
          NumericVector v = col_stratum;
          data_outcome.push_back(v[ustratum], s);
        } else if (type_stratum == STRSXP) {
          StringVector v = col_stratum;
          data_outcome.push_back(v[ustratum], s);
        }
      }
    }
    
    // construct the confidence interval for HR
    if (!boot) { // use log-rank p-value to construct CI for HR if no boot
      double loghr = log(hrhat);
      double seloghr = logRankZ != 0 ? loghr/logRankZ : NA_REAL;
      hrlower = exp(loghr - zcrit*seloghr);
      hrupper = exp(loghr + zcrit*seloghr);
      hr_CI_type = "log-rank p-value";
    } else { // bootstrap the entire process to construct CI for HR
      if (seed != NA_INTEGER) set_seed(seed);
      
      IntegerVector idb(n), stratumb(n), treatb(n), eventb(n);
      NumericVector timeb(n), rxb(n), censor_timeb(n);
      NumericMatrix zb(n,p), z_aftb(n,q+p);
      
      int B = n*n_boot;
      IntegerVector boot_indexc(B);
      IntegerVector idc(B), stratumc(B), treatc(B), eventc(B);
      NumericVector timec(B), rxc(B), censor_timec(B);
      NumericMatrix z_aftc(B,q+p);
      int index1 = 0; // current index for fail boots data
      
      // sort data by treatment group, stratum and id
      IntegerVector order = seq(0, n-1);
      if (has_stratum) {
        std::sort(order.begin(), order.end(), [&](int i, int j) {
          return std::tie(treatn[i], stratumn[i], idn[i]) <
            std::tie(treatn[j], stratumn[j], idn[j]);
        });
      } else {
        std::sort(order.begin(), order.end(), [&](int i, int j) {
          return std::tie(treatn[i], idn[i]) < std::tie(treatn[j], idn[j]);
        });
      }
      
      idn = idn[order];
      stratumn = stratumn[order];
      timen = timen[order];
      eventn = eventn[order];
      treatn = treatn[order];
      rxn = rxn[order];
      censor_timen = censor_timen[order];
      zn = subset_matrix_by_row(zn, order);
      z_aftn = subset_matrix_by_row(z_aftn, order);
      
      IntegerVector tsx(1,0); // first observation within each treat/stratum
      for (int i=1; i<n; ++i) {
        if (treatn[i] != treatn[i-1] || stratumn[i] != stratumn[i-1]) {
          tsx.push_back(i);
        }
      }
      
      int ntss = static_cast<int>(tsx.size());
      tsx.push_back(n);
      
      for (k=0; k<n_boot; ++k) {
        // sample the data with replacement by treatment group and stratum
        for (int h=0; h<ntss; ++h) {
          for (int i=tsx[h]; i<tsx[h+1]; ++i) {
            double u = R::runif(0,1);
            int j = tsx[h] + static_cast<int>(std::floor(u*(tsx[h+1]-tsx[h])));
            
            idb[i] = idn[j];
            stratumb[i] = stratumn[j];
            timeb[i] = timen[j];
            eventb[i] = eventn[j];
            treatb[i] = treatn[j];
            rxb[i] = rxn[j];
            censor_timeb[i] = censor_timen[j];
            zb(i,_) = zn(j,_);
            z_aftb(i,_) = z_aftn(j,_);
          }
        }
        
        List out = f(idb, stratumb, timeb, eventb, treatb, rxb, censor_timeb, 
                     zb, z_aftb);
        
        fails[k] = out["fail"];
        hrhats[k] = out["hrhat"];
        psihats[k] = out["psihat"];
        
        if (fails[k]) {
          for (int i=0; i<n; ++i) {
            int j = index1 + i;
            boot_indexc[j] = k+1;
            idc[j] = idb[i];
            stratumc[j] = stratumb[i];
            timec[j] = timeb[i];
            eventc[j] = eventb[i];
            treatc[j] = treatb[i];
            rxc[j] = rxb[i];
            censor_timec[j] = censor_timeb[i];
            z_aftc(j,_) = z_aftb(i,_);
          }
          index1 += n;
        }
      }
      
      if (is_true(any(fails))) {
        IntegerVector sub = seq(0,index1-1);
        boot_indexc = boot_indexc[sub];
        idc = idc[sub];
        stratumc = stratumc[sub];
        timec = timec[sub];
        eventc = eventc[sub];
        treatc = treatc[sub];
        rxc = rxc[sub];
        censor_timec = censor_timec[sub];
        z_aftc = subset_matrix_by_row(z_aftc,sub);
        
        fail_boots_data = List::create(
          Named("boot_index") = boot_indexc,
          Named("time") = timec,
          Named("event") = eventc,
          Named("treated") = treatc,
          Named("rx") = rxc,
          Named("censor_time") = censor_timec
        );
        
        for (int j=0; j<q+p; ++j) {
          String zj = covariates_aft[j+1];
          NumericVector u = z_aftc(_,j);
          fail_boots_data.push_back(u, zj);
        }
        
        if (type_id == INTSXP) {
          fail_boots_data.push_back(idwi[idc], id);
        } else if (type_id == REALSXP) {
          fail_boots_data.push_back(idwn[idc], id);
        } else if (type_id == STRSXP) {
          fail_boots_data.push_back(idwc[idc], id);
        }
        
        if (type_treat == LGLSXP || type_treat == INTSXP) {
          fail_boots_data.push_back(treatwi[1-treatc], treat);
        } else if (type_treat == REALSXP) {
          fail_boots_data.push_back(treatwn[1-treatc], treat);
        } else if (type_treat == STRSXP) {
          fail_boots_data.push_back(treatwc[1-treatc], treat);
        }
        
        if (has_stratum) {
          for (int i=0; i<p_stratum; ++i) {
            std::string s = as<std::string>(stratum[i]);
            SEXP col_stratum = u_stratum[s];
            SEXPTYPE type_stratum = TYPEOF(col_stratum);
            if (type_stratum == INTSXP) {
              IntegerVector v = col_stratum;
              fail_boots_data.push_back(v[stratumc], s);
            } else if (type_stratum == REALSXP) {
              NumericVector v = col_stratum;
              fail_boots_data.push_back(v[stratumc], s);
            } else if (type_stratum == STRSXP) {
              StringVector v = col_stratum;
              fail_boots_data.push_back(v[stratumc], s);
            }
          }
        }
      }
      
      // obtain bootstrap confidence interval for HR
      double loghr = log(hrhat);
      LogicalVector ok = (!fails) & !is_na(hrhats);
      int n_ok = sum(ok);
      NumericVector subset_hrhats = hrhats[ok];
      NumericVector loghrs = log(subset_hrhats);
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
  }

  List result = List::create(
    Named("psi") = psihat,
    Named("psi_CI") = NumericVector::create(psilower, psiupper),
    Named("psi_CI_type") = psi_CI_type,
    Named("logrank_pvalue") = logRankPValue,
    Named("cox_pvalue") = pvalue,
    Named("hr") = hrhat,
    Named("hr_CI") = NumericVector::create(hrlower, hrupper),
    Named("hr_CI_type") = hr_CI_type,
    Named("event_summary") = as<DataFrame>(event_summary),
    Named("Sstar") = as<DataFrame>(Sstar),
    Named("kmstar") = as<DataFrame>(kmstar),
    Named("data_aft") = as<DataFrame>(data_aft),
    Named("fit_aft") = fit_aft,
    Named("res_aft") = res_aft,
    Named("data_outcome") = as<DataFrame>(data_outcome),
    Named("km_outcome") = as<DataFrame>(km_outcome),
    Named("lr_outcome") = as<DataFrame>(lr_outcome),
    Named("fit_outcome") = fit_outcome,
    Named("fail") = fail,
    Named("psimissing") = psimissing);
  
  if (boot) {
    result.push_back(fails, "fail_boots");
    result.push_back(hrhats, "hr_boots");
    result.push_back(psihats, "psi_boots");
    if (is_true(any(fails))) {
      result.push_back(as<DataFrame>(fail_boots_data), "fail_boots_data");
    }
  }
  
  return result;
}
