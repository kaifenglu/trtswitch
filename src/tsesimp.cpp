#include "utilities.h"
#include "survival_analysis.h"

using namespace Rcpp;


// [[Rcpp::export]]
List tsesimpcpp(const DataFrame data,
                const std::string id = "id",
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
                const bool strata_main_effect_only = true,
                const bool recensor = true,
                const bool admin_recensor_only = true,
                const bool swtrt_control_only = true,
                const double alpha = 0.05,
                const std::string ties = "efron",
                const double offset = 1,
                const bool boot = true,
                const int n_boot = 1000,
                const int seed = NA_INTEGER) {
  
  int k, n = data.nrow();
  
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
  bool has_censor_time = hasVariable(data, censor_time);
  bool has_pd = hasVariable(data, pd);
  bool has_pd_time = hasVariable(data, pd_time);
  bool has_swtrt = hasVariable(data, swtrt);
  bool has_swtrt_time = hasVariable(data, swtrt_time);
  
  
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
  
  
  if (!has_pd) stop("data must contain the pd variable"); 
  
  SEXP col_pd = data[pd];
  SEXPTYPE type_pd = TYPEOF(col_pd);
  
  IntegerVector pdn(n);
  if (type_pd == LGLSXP || type_pd == INTSXP) {
    IntegerVector pdnz = col_pd;
    if (is_true(any((pdnz != 1) & (pdnz != 0)))) {
      stop("pd must be 1 or 0 for each subject");
    } else {
      pdn = clone(pdnz);
    }
  } else if (type_pd == REALSXP) {
    NumericVector pdnz = col_pd;
    if (is_true(any((pdnz != 1) & (pdnz != 0)))) {
      stop("pd must be 1 or 0 for each subject");
    } else {
      NumericVector pdnz2 = clone(pdnz);
      pdn = as<IntegerVector>(pdnz2);
    }
  } else {
    stop("pd must take logical, integer, or real values");
  }
  
  if (!has_pd_time) {
    stop("data must contain the pd_time variable"); 
  }
  
  
  SEXP col_pd_time = data[pd_time];
  SEXPTYPE type_pd_time = TYPEOF(col_pd_time);
  if (type_pd_time != INTSXP && type_pd_time != REALSXP) {
    stop("pd_time must take numeric values");
  }
  
  NumericVector pd_timenz = col_pd_time;
  NumericVector pd_timen = clone(pd_timenz);
  for (int i=0; i<n; ++i) {
    if (pdn[i] == 1 && std::isnan(pd_timen[i])) {
      stop("pd_time must not be missing when pd=1");
    }
    if (pdn[i] == 1 && pd_timen[i] < 0.0) {
      stop("pd_time must be nonnegative when pd=1");
    }
  }
  
  
  if (!has_swtrt) stop("data must contain the swtrt variable");
  
  SEXP col_swtrt = data[swtrt];
  SEXPTYPE type_swtrt = TYPEOF(col_swtrt);
  
  IntegerVector swtrtn(n);
  if (type_swtrt == LGLSXP || type_swtrt == INTSXP) {
    IntegerVector swtrtnz = col_swtrt;
    if (is_true(any((swtrtnz != 1) & (swtrtnz != 0)))) {
      stop("swtrt must be 1 or 0 for each subject");
    } else {
      swtrtn = clone(swtrtnz);
    }
  } else if (type_swtrt == REALSXP) {
    NumericVector swtrtnz = col_swtrt;
    if (is_true(any((swtrtnz != 1) & (swtrtnz != 0)))) {
      stop("swtrt must be 1 or 0 for each subject");
    } else {
      NumericVector swtrtnz2 = clone(swtrtnz);
      swtrtn = as<IntegerVector>(swtrtnz2);
    }
  } else {
    stop("swtrt must take logical, integer, or real values");
  }
  
  
  if (is_false(any((pdn == 1) & (swtrtn == 1) & (treatn == 0)))) {
    stop("at least 1 pd and swtrt is needed in the control group");
  }
  
  if (!swtrt_control_only) {
    if (is_false(any((pdn == 1) & (swtrtn == 1) & (treatn == 1)))) {
      stop("at least 1 pd and swtrt is needed in the treatment group");
    }
  }
  
  if (!has_swtrt_time) stop("data must contain the swtrt_time variable"); 
  
  SEXP col_swtrt_time = data[swtrt_time];
  SEXPTYPE type_swtrt_time = TYPEOF(col_swtrt_time);
  if (type_swtrt_time != INTSXP && type_swtrt_time != REALSXP) {
    stop("swtrt_time must take numeric values");
  }
  
  NumericVector swtrt_timenz = col_swtrt_time;
  NumericVector swtrt_timen = clone(swtrt_timenz);
  for (int i=0; i<n; ++i) {
    if (swtrtn[i] == 1 && std::isnan(swtrt_timen[i])) {
      stop("swtrt_time must not be missing when swtrt=1");
    }
    if (swtrtn[i] == 1 && swtrt_timen[i] < 0.0) {
      stop("swtrt_time must be nonnegative when swtrt=1");
    }
    if (swtrtn[i] == 1 && swtrt_timen[i] > timen[i]) {
      stop("swtrt_time must be less than or equal to time");
    }
  }
  
  // if the patient switched before pd, set pd time equal to switch time
  for (int i=0; i<n; ++i) {
    if (pdn[i] == 1 && swtrtn[i] == 1 && swtrt_timen[i] < pd_timen[i]) {
      pd_timen[i] = swtrt_timen[i];
    }
    
    if (pdn[i] == 0 && swtrtn[i] == 1) {
      pdn[i] = 1; pd_timen[i] = swtrt_timen[i];
    }
  }
  
  // make sure offset is less than or equal to observed time variables
  for (int i=0; i<n; ++i) {
    if (pdn[i] == 1 && pd_timen[i] < offset) {
      stop("pd_time must be great than or equal to offset");
    }
    if (swtrtn[i] == 1 && swtrt_timen[i] < offset) {
      stop("swtrt_time must be great than or equal to offset");
    }
  }
  
  // ensure pd time < os time
  for (int i=0; i<n; ++i) {
    if (pdn[i] == 1 && pd_timen[i] == timen[i]) {
      timen[i] = timen[i] + 1.0e-8;
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
    covariates[j+1] = zj;
    zn(_,j) = u;
  }
  
  
  int q; // number of columns corresponding to the strata effects
  if (strata_main_effect_only) {
    q = sum(d - 1);
  } else {
    q = nstrata - 1;
  }
  
  // covariates for the accelerated failure time model for control with pd
  // including swtrt, stratum and base2_cov
  StringVector covariates_aft(q+p2+1);
  NumericMatrix z_aftn(n,q+p2);
  covariates_aft[0] = "swtrt";
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
  
  for (int j=0; j<p2; ++j) {
    String zj = base2_cov[j];
    
    if (!hasVariable(data, zj)) {
      stop("data must contain the variables in base2_cov");
    }
    
    if (zj == treat) {
      stop("treat should be excluded from base2_cov");
    }
    
    NumericVector u = data[zj];
    covariates_aft[q+j+1] = zj;
    z_aftn(_,q+j) = u;
  }
  
  std::string dist = aft_dist;
  std::for_each(dist.begin(), dist.end(), [](char & c) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });
  
  if (dist == "log-logistic" || dist == "llogistic") {
    dist = "loglogistic";
  } else if (dist == "log-normal" || dist == "lnormal") {
    dist = "lognormal";
  }
  
  if (!(dist == "exponential" || dist == "weibull" ||
      dist == "lognormal" || dist == "loglogistic")) {
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
  
  
  // exclude observations with missing values
  LogicalVector sub(n,1);
  for (int i=0; i<n; ++i) {
    if (idn[i] == NA_INTEGER || stratumn[i] == NA_INTEGER || 
        std::isnan(timen[i]) || eventn[i] == NA_INTEGER || 
        treatn[i] == NA_INTEGER || std::isnan(censor_timen[i]) || 
        pdn[i] == NA_INTEGER || swtrtn[i] == NA_INTEGER) {
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
  censor_timen = censor_timen[order];
  pdn = pdn[order];
  pd_timen = pd_timen[order];
  swtrtn = swtrtn[order];
  swtrt_timen = swtrt_timen[order];
  if (p > 0) zn = subset_matrix_by_row(zn, order);
  z_aftn = subset_matrix_by_row(z_aftn, order);
  n = sum(sub);
  if (n == 0) stop("no observations left after removing missing values");
  
  // summarize number of deaths and switches by treatment arm
  IntegerVector treat_out = IntegerVector::create(0, 1);
  NumericVector n_total(2);
  NumericVector n_event(2);
  NumericVector n_pd(2);
  NumericVector n_switch(2);
  for (int i = 0; i < n; ++i) {
    int g = treatn[i];
    n_total[g]++;
    if (eventn[i] == 1) n_event[g]++;
    if (pdn[i] == 1) n_pd[g]++;
    if (swtrtn[i] == 1) n_switch[g]++;
  }
  
  // Compute percentages
  NumericVector pct_event(2);
  NumericVector pct_pd(2);
  NumericVector pct_switch(2);
  for (int g = 0; g < 2; g++) {
    pct_event[g] = 100.0 * n_event[g] / n_total[g];
    pct_pd[g] = 100.0 * n_pd[g] / n_total[g];
    pct_switch[g] = 100.0 * n_switch[g] / n_total[g];
  }
  
  // Combine count and percentage
  List event_summary = List::create(
    _["treated"] = treat_out,
    _["n"] = n_total,
    _["event_n"] = n_event,
    _["event_pct"] = pct_event,
    _["pd_n"] = n_pd,
    _["pd_pct"] = pct_pd,
    _["switch_n"] = n_switch,
    _["switch_pct"] = pct_switch
  );
  
  DataFrame lr = lrtest(data, "", stratum, treat, time, "", event, "", 0,0,0);
  double logRankPValue = lr["logRankPValue"];
  double zcrit = R::qnorm(1-alpha/2, 0, 1, 1, 0);
  
  k = -1; // indicate the observed data
  auto f = [&k, has_stratum, stratum, p_stratum, u_stratum, 
            type_treat, treat, treatwi, treatwn, treatwc, 
            type_id, id, idwi, idwn, idwc,
            n, q, p, p2, covariates, covariates_aft, dist, 
            recensor, swtrt_control_only, alpha, zcrit, ties, offset](
                IntegerVector& idb, IntegerVector& stratumb, 
                NumericVector& timeb, IntegerVector& eventb, 
                IntegerVector& treatb, NumericVector& censor_timeb, 
                IntegerVector& pdb, NumericVector& pd_timeb, 
                IntegerVector& swtrtb, NumericVector& swtrt_timeb, 
                NumericMatrix& zb, NumericMatrix& z_aftb)->List {
                  bool fail = false; // whether any model fails to converge
                  NumericVector init(1, NA_REAL);
                  
                  // time and event adjusted for treatment switching
                  NumericVector t_star = clone(timeb);
                  IntegerVector d_star = clone(eventb);
                  
                  double psi0hat = NA_REAL;
                  double psi0lower = NA_REAL, psi0upper = NA_REAL;
                  double psi1hat = NA_REAL;
                  double psi1lower = NA_REAL, psi1upper = NA_REAL;
                  
                  // initialize data_aft and fit_aft
                  List data_aft(2), fit_aft(2), res_aft(2);
                  if (k == -1) {
                    for (int h=0; h<2; ++h) {
                      List data_x = List::create(
                        Named("data") = R_NilValue,
                        Named(treat) = R_NilValue
                      );
                      
                      if (type_treat == LGLSXP ||
                          type_treat == INTSXP) {
                        data_x[treat] = treatwi[1-h];
                      } else if (type_treat == REALSXP) {
                        data_x[treat] = treatwn[1-h];
                      } else if (type_treat == STRSXP) {
                        data_x[treat] = treatwc[1-h];
                      }
                      
                      data_aft[h] = data_x;
                      
                      List fit_x = List::create(
                        Named("fit") = R_NilValue,
                        Named(treat) = R_NilValue
                      );
                      
                      if (type_treat == LGLSXP ||
                          type_treat == INTSXP) {
                        fit_x[treat] = treatwi[1-h];
                      } else if (type_treat == REALSXP) {
                        fit_x[treat] = treatwn[1-h];
                      } else if (type_treat == STRSXP) {
                        fit_x[treat] = treatwc[1-h];
                      }
                      
                      fit_aft[h] = fit_x;
                      
                      List res_x = List::create(
                        Named("res") = R_NilValue,
                        Named(treat) = R_NilValue
                      );
                      
                      if (type_treat == LGLSXP ||
                          type_treat == INTSXP) {
                        res_x[treat] = treatwi[1-h];
                      } else if (type_treat == REALSXP) {
                        res_x[treat] = treatwn[1-h];
                      } else if (type_treat == STRSXP) {
                        res_x[treat] = treatwc[1-h];
                      }
                      
                      res_aft[h] = res_x;
                    }
                  }
                  
                  List data_outcome;
                  List fit_outcome;
                  double hrhat = NA_REAL, hrlower = NA_REAL, 
                    hrupper = NA_REAL, pvalue = NA_REAL;
                  List km_outcome, lr_outcome;
                  
                  bool psimissing = false;
                  
                  // treat arms that include patients who switched treatment
                  IntegerVector treats(1);
                  treats[0] = 0;
                  if (!swtrt_control_only) {
                    treats.push_back(1);
                  }
                  
                  int K = static_cast<int>(treats.size());
                  for (int h=0; h<K; ++h) {
                    // post progression data
                    IntegerVector l = which((treatb == h) & (pdb == 1));
                    IntegerVector id2 = idb[l];
                    NumericVector time2 = timeb[l] - pd_timeb[l] + offset;
                    IntegerVector event2 = eventb[l];
                    IntegerVector swtrt2 = swtrtb[l];
                    
                    List data1 = List::create(
                      Named("pps") = time2,
                      Named("event") = event2,
                      Named("swtrt") = swtrt2);
                    
                    for (int j=0; j<q+p2; ++j) {
                      String zj = covariates_aft[j+1];
                      NumericVector u = z_aftb(_,j);
                      data1.push_back(u[l], zj);
                    }
                    
                    List fit1 = liferegcpp(
                      data1, "", "", "pps", "", "event", 
                      covariates_aft, "", "", "", dist, init, 
                      0, 0, alpha, 50, 1.0e-9);
                    
                    DataFrame sumstat1 = DataFrame(fit1["sumstat"]);
                    bool fail1 = sumstat1["fail"];
                    if (fail1) fail = true;
                    
                    DataFrame parest1 = DataFrame(fit1["parest"]);
                    NumericVector beta1 = parest1["beta"];
                    double psihat = -beta1[1];
                    double psilower = NA_REAL, psiupper = NA_REAL;
                    if (k == -1) {
                      NumericVector sebeta1 = parest1["sebeta"];
                      psilower = -beta1[1] - zcrit*sebeta1[1];
                      psiupper = -beta1[1] + zcrit*sebeta1[1];
                    }
                    
                    NumericVector res;
                    if (k == -1) {
                      NumericMatrix vbeta1(q+p2+1, q+p2+1);
                      NumericMatrix rr = residuals_liferegcpp(
                        beta1, vbeta1, data1, "", "pps", "", "event",
                        covariates_aft, "", "", "", dist, "deviance", 0, 0);
                      res = rr(_,0);
                    }
                    
                    // update treatment-specific causal parameter estimates
                    if (h == 0) {
                      psi0hat = psihat;
                      if (k == -1) {
                        psi0lower = psilower;
                        psi0upper = psiupper;
                      }
                    } else {
                      psi1hat = psihat;
                      if (k == -1) {
                        psi1lower = psilower;
                        psi1upper = psiupper;
                      }
                    }
                    
                    if (!std::isnan(psihat)) {
                      // calculate counter-factual survival times
                      double a = exp(psihat);
                      double c0 = std::min(1.0, a);
                      for (int i=0; i<n; ++i) {
                        if (treatb[i] == h) {
                          double b2, u_star, c_star;
                          if (swtrtb[i] == 1) {
                            b2 = pd_timeb[i] - offset;
                            u_star = b2 + (timeb[i] - b2)*a;
                          } else {
                            u_star = timeb[i];
                          }
                          
                          if (recensor) {
                            c_star = censor_timeb[i]*c0;
                            t_star[i] = std::min(u_star, c_star);
                            d_star[i] = c_star < u_star ? 0 : eventb[i];
                          } else {
                            t_star[i] = u_star;
                            d_star[i] = eventb[i];
                          }
                        }
                      }
                      
                      // update data_aft, fit_aft, and res_aft
                      if (k == -1) {
                        IntegerVector stratum2 = stratumb[l];
                        
                        if (has_stratum) {
                          for (int i=0; i<p_stratum; ++i) {
                            std::string s = as<std::string>(stratum[i]);
                            SEXP col_stratum = u_stratum[s];
                            SEXPTYPE type_stratum = TYPEOF(col_stratum);
                            if (type_stratum == INTSXP) {
                              IntegerVector v = col_stratum;
                              data1.push_back(v[stratum2], s);
                            } else if (type_stratum == REALSXP) {
                              NumericVector v = col_stratum;
                              data1.push_back(v[stratum2], s);
                            } else if (type_stratum == STRSXP) {
                              StringVector v = col_stratum;
                              data1.push_back(v[stratum2], s);
                            }
                          }
                        }
                        
                        if (type_id == INTSXP) {
                          data1.push_front(idwi[id2], id);
                        } else if (type_id == REALSXP) {
                          data1.push_front(idwn[id2], id);
                        } else if (type_id == STRSXP) {
                          data1.push_front(idwc[id2], id);
                        }
                        
                        List data_x = data_aft[h];
                        data_x["data"] = as<DataFrame>(data1);
                        data_aft[h] = data_x;
                        
                        List fit_x = fit_aft[h];
                        fit_x["fit"] = fit1;
                        fit_aft[h] = fit_x;
                        
                        List res_x = res_aft[h];
                        res_x["res"] = res;
                        res_aft[h] = res_x;
                      }
                    } else {
                      psimissing = 1;
                    }
                  }
                  
                  
                  if (!psimissing) {
                    // Cox model for hypothetical treatment effect estimate
                    data_outcome = List::create(
                      Named("uid") = idb,
                      Named("t_star") = t_star,
                      Named("d_star") = d_star,
                      Named("treated") = treatb);
                    
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
                    NumericVector sebeta = parest["sebeta"];
                    NumericVector pval = parest["p"];
                    hrhat = exp(beta[0]);
                    hrlower = exp(beta[0] - zcrit*sebeta[0]);
                    hrupper = exp(beta[0] + zcrit*sebeta[0]);
                    pvalue = pval[0];
                  }
                  
                  List out;
                  if (k == -1) {
                    out = List::create(
                      Named("data_aft") = data_aft,
                      Named("fit_aft") = fit_aft,
                      Named("res_aft") = res_aft,
                      Named("data_outcome") = data_outcome,
                      Named("km_outcome") = km_outcome,
                      Named("lr_outcome") = lr_outcome,
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
                      Named("pvalue") = pvalue,
                      Named("fail") = fail,
                      Named("psimissing") = psimissing);
                  } else {
                    out = List::create(
                      Named("psihat") = psi0hat,
                      Named("psi1hat") = psi1hat,
                      Named("hrhat") = hrhat,
                      Named("hrlower") = hrlower,
                      Named("hrupper") = hrupper,
                      Named("pvalue") = pvalue,
                      Named("fail") = fail,
                      Named("psimissing") = psimissing);
                  }
                  
                  return out;
                };
  
  List out = f(idn, stratumn, timen, eventn, treatn, censor_timen,
               pdn, pd_timen, swtrtn, swtrt_timen, zn, z_aftn);
  
  List data_aft = out["data_aft"];
  List fit_aft = out["fit_aft"];
  List res_aft = out["res_aft"];
  List data_outcome = out["data_outcome"];
  List km_outcome = out["km_outcome"];
  List lr_outcome = out["lr_outcome"];
  List fit_outcome = out["fit_outcome"];
  
  double psihat = out["psihat"];
  double psilower = out["psilower"];
  double psiupper = out["psiupper"];
  double psi1hat = out["psi1hat"];
  double psi1lower = out["psi1lower"];
  double psi1upper = out["psi1upper"];
  String psi_CI_type = "AFT model";
  double hrhat = out["hrhat"];
  double hrlower = out["hrlower"];
  double hrupper = out["hrupper"];
  double pvalue = out["pvalue"];
  bool fail = out["fail"];
  bool psimissing = out["psimissing"];
  
  NumericVector hrhats(n_boot), psihats(n_boot), psi1hats(n_boot);
  LogicalVector fails(n_boot);
  DataFrame fail_boots_data;
  String hr_CI_type;
  
  if (!psimissing) {
    IntegerVector uid = data_outcome["uid"];
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
      IntegerVector ustratum = data_outcome["ustratum"];
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
    if (!boot) { // use Cox model to construct CI for HR if no boot
      hr_CI_type = "Cox model";
    } else { // bootstrap the entire process to construct CI for HR
      if (seed != NA_INTEGER) set_seed(seed);
      
      IntegerVector idb(n), stratumb(n), treatb(n), eventb(n);
      IntegerVector pdb(n), swtrtb(n);
      NumericVector timeb(n), censor_timeb(n), pd_timeb(n), swtrt_timeb(n);
      NumericMatrix zb(n,p), z_aftb(n,q+p2);
      
      int B = n*n_boot;
      IntegerVector boot_indexc(B);
      IntegerVector idc(B), stratumc(B), treatc(B), eventc(B);
      IntegerVector pdc(B), swtrtc(B);
      NumericVector timec(B), censor_timec(B), pd_timec(B), swtrt_timec(B);
      NumericMatrix z_aftc(B,q+p2);
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
      censor_timen = censor_timen[order];
      pdn = pdn[order];
      pd_timen = pd_timen[order];
      swtrtn = swtrtn[order];
      swtrt_timen = swtrt_timen[order];
      zn = subset_matrix_by_row(zn, order);
      z_aftn = subset_matrix_by_row(z_aftn, order);
      
      IntegerVector tsx(1,0); // first observation within each treat/stratum
      for (int i=1; i<n; ++i) {
        if (treatn[i] != treatn[i-1] || stratumn[i] != stratumn[i-1]) {
          tsx.push_back(i);
        }
      }
      
      int ntss = static_cast<int>(tsx.size());
      tsx.push_back(n); // add the end index
      
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
            censor_timeb[i] = censor_timen[j];
            pdb[i] = pdn[j];
            pd_timeb[i] = pd_timen[j];
            swtrtb[i] = swtrtn[j];
            swtrt_timeb[i] = swtrt_timen[j];
            zb(i,_) = zn(j,_);
            z_aftb(i,_) = z_aftn(j,_);
          }
        }
        
        List out = f(idb, stratumb, timeb, eventb, treatb, censor_timeb,
                     pdb, pd_timeb, swtrtb, swtrt_timeb, zb, z_aftb);
        
        fails[k] = out["fail"];
        hrhats[k] = out["hrhat"];
        psihats[k] = out["psihat"];
        psi1hats[k] = out["psi1hat"];
        
        if (fails[k]) {
          for (int i=0; i<n; ++i) {
            int j = index1 + i;
            boot_indexc[j] = k+1;
            idc[j] = idb[i];
            stratumc[j] = stratumb[i];
            timec[j] = timeb[i];
            eventc[j] = eventb[i];
            treatc[j] = treatb[i];
            censor_timec[j] = censor_timeb[i];
            pdc[j] = pdb[i];
            pd_timec[j] = pd_timeb[i];
            swtrtc[j] = swtrtb[i];
            swtrt_timec[j] = swtrt_timeb[i];
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
        censor_timec = censor_timec[sub];
        pdc = pdc[sub];
        pd_timec = pd_timec[sub];
        swtrtc = swtrtc[sub];
        swtrt_timec = swtrt_timec[sub];
        z_aftc = subset_matrix_by_row(z_aftc,sub);
        
        fail_boots_data = List::create(
          Named("boot_index") = boot_indexc,
          Named("time") = timec,
          Named("event") = eventc,
          Named("treated") = treatc,
          Named("censor_time") = censor_timec,
          Named("pd") = pdc,
          Named("pd_time") = pd_timec,
          Named("swtrt") = swtrtc,
          Named("swtrt_time") = swtrt_timec
        );
        
        for (int j=0; j<q+p2; ++j) {
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
      
      NumericVector psi1hats1 = psi1hats[ok];
      double sdpsi1 = sd(psi1hats1);
      psi1lower = psi1hat - tcrit*sdpsi1;
      psi1upper = psi1hat + tcrit*sdpsi1;
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
    Named("data_aft") = data_aft,
    Named("fit_aft") = fit_aft,
    Named("res_aft") = res_aft,
    Named("data_outcome") = as<DataFrame>(data_outcome),
    Named("km_outcome") = as<DataFrame>(km_outcome),
    Named("lr_outcome") = as<DataFrame>(lr_outcome),
    Named("fit_outcome") = fit_outcome,
    Named("fail") = fail,
    Named("psimissing") = psimissing);
  
  if (!swtrt_control_only) {
    result.push_back(psi1hat, "psi_trt");
    NumericVector psi1_CI = NumericVector::create(psi1lower, psi1upper);
    result.push_back(psi1_CI, "psi_trt_CI");
  }
  
  if (boot) {
    result.push_back(fails, "fail_boots");
    result.push_back(hrhats, "hr_boots");
    result.push_back(psihats, "psi_boots");
    if (!swtrt_control_only) {
      result.push_back(psi1hats, "psi_trt_boots");
    }
    if (is_true(any(fails))) {
      result.push_back(as<DataFrame>(fail_boots_data), "fail_boots_data");
    }
  }
  
  return result;
}
