#include "utilities.h"
#include "survival_analysis.h"
#include "logistic_regression.h"
#include "splines.h"

using namespace Rcpp;


List est_psi_tsegest(int n2, int q, int p2, int nids2, 
                     IntegerVector idx2, IntegerVector stratum3, 
                     IntegerVector os3, NumericVector os_time3, 
                     NumericVector censor_time3,  
                     IntegerVector swtrt3, NumericVector swtrt_time3, 
                     IntegerVector id2, IntegerVector cross2, 
                     NumericVector tstart2, NumericVector tstop2,
                     StringVector covariates_lgs, NumericMatrix z_lgs2, 
                     int ns_df, NumericMatrix s, bool firth, bool flic, 
                     bool recensor, double alpha, 
                     std::string ties, double offset, double x) {
  
  NumericVector init(1, NA_REAL);
  
  double a = exp(x);
  double c0 = std::min(1.0, a);
  
  // counterfactual survival times and event indicators
  NumericVector t_star(nids2);
  IntegerVector d_star(nids2);
  for (int i=0; i<nids2; ++i) {
    double b2, u_star, c_star;
    if (swtrt3[i] == 1) {
      b2 = swtrt_time3[i] - offset;
      u_star = b2 + (os_time3[i] - b2)*a;
    } else {
      u_star = os_time3[i];
    }
    
    if (recensor) {
      c_star = censor_time3[i]*c0;
      t_star[i] = std::min(u_star, c_star);
      d_star[i] = c_star < u_star ? 0 : os3[i];
    } else {
      t_star[i] = u_star;
      d_star[i] = os3[i];
    }
  }
  
  // martingale residuals from the null model
  DataFrame data1 = DataFrame::create(
    Named("t_star") = t_star,
    Named("d_star") = d_star,
    Named("ustratum") = stratum3
  );
  
  List fit1 = phregcpp(
    data1, "", "ustratum", "t_star", "", "d_star",
    "", "", "", "", ties, init, 
    0, 0, 1, 0, 0, alpha, 50, 1.0e-9);

  // replicate counterfactual residuals within subjects
  NumericVector resid3 = fit1["residuals"];
  NumericVector resid(n2);
  for (int i=0; i<nids2; ++i) {
    for (int j=idx2[i]; j<idx2[i+1]; ++j) {
      resid[j] = resid3[i];
    }
  }
  // logistic regression switching model
  DataFrame data2 = DataFrame::create(
    Named("uid") = id2,
    Named("tstart") = tstart2,
    Named("tstop") = tstop2,
    Named("cross") = cross2,
    Named("counterfactual") = resid);
  
  for (int j=0; j<q+p2; ++j) {
    String zj = covariates_lgs[j+1];
    NumericVector u = z_lgs2(_,j);
    data2.push_back(u,zj);
  }
  for (int j=0; j<ns_df; ++j) {
    String zj = covariates_lgs[q+p2+j+1];
    NumericVector u = s(_,j);
    data2.push_back(u,zj);
  }
  
  List fit2 = logisregcpp(
    data2, "", "cross", covariates_lgs, "", "", 
    "", "uid", "logit", init, 
    1, firth, flic, 0, alpha, 50, 1.0e-9);
    
  DataFrame sumstat = DataFrame(fit2["sumstat"]);
  bool fail = sumstat["fail"];
  
  DataFrame parest = DataFrame(fit2["parest"]);
  NumericVector z = parest["z"];
  double z_counterfactual = z[1];
  
  List out = List::create(
    Named("data_nullcox") = data1,
    Named("fit_nullcox") = fit1,
    Named("data_logis") = data2,
    Named("fit_logis") = fit2,
    Named("z_counterfactual") = z_counterfactual,
    Named("fail") = fail);
  
  return out;
};


// [[Rcpp::export]]
List tsegestcpp(
    const DataFrame data,
    const std::string id = "id",
    const StringVector& stratum = "",
    const std::string tstart = "tstart",
    const std::string tstop = "tstop",
    const std::string event = "event",
    const std::string treat = "treat",
    const std::string censor_time = "censor_time",
    const std::string pd = "pd",
    const std::string pd_time = "pd_time",
    const std::string swtrt = "swtrt",
    const std::string swtrt_time = "swtrt_time",
    const StringVector& base_cov = "",
    const StringVector& conf_cov = "",
    const bool strata_main_effect_only = true,
    const int ns_df = 3,
    const bool firth = false,
    const bool flic = false,
    const double low_psi = -2,
    const double hi_psi = 2,
    const int n_eval_z = 101,
    const bool recensor = true,
    const bool admin_recensor_only = true,
    const bool swtrt_control_only = true,
    const bool gridsearch = true,
    const std::string root_finding = "brent",
    const double alpha = 0.05,
    const std::string ties = "efron",
    const double tol = 1.0e-6,
    const double offset = 1,
    const bool boot = true,
    const int n_boot = 1000,
    const int seed = NA_INTEGER) {
  
  int k, n = data.nrow();
  int p = static_cast<int>(base_cov.size());
  if (p == 1 && (base_cov[0] == "" || base_cov[0] == "none")) p = 0;
  
  int p2 = static_cast<int>(conf_cov.size());
  if (p2 == 1 && (conf_cov[0] == "" || conf_cov[0] == "none")) p2 = 0;
  
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
  bool has_tstart = hasVariable(data, tstart);
  bool has_tstop = hasVariable(data, tstop);
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
  
  
  
  if (!has_tstart) stop("data must contain the tstart variable");
  
  SEXP col_tstart = data[tstart];
  SEXPTYPE type_tstart = TYPEOF(col_tstart);
  if (type_tstart != INTSXP && type_tstart != REALSXP) {
    stop("tstart must take numeric values");
  }
  
  NumericVector tstartnz = col_tstart;
  NumericVector tstartn = clone(tstartnz);
  if (is_true(any(tstartn < 0.0))) {
    stop("tstart must be nonnegative for each observation");
  }
  
  
  if (!has_tstop) stop("data must contain the tstop variable");
  
  SEXP col_tstop = data[tstop];
  SEXPTYPE type_tstop = TYPEOF(col_tstop);
  if (type_tstop != INTSXP && type_tstop != REALSXP) {
    stop("tstop must take numeric values");
  }
  
  NumericVector tstopnz = col_tstop;
  NumericVector tstopn = clone(tstopnz);
  if (is_true(any(tstopn <= tstartn))) {
    stop("tstop must be greater than tstart for each observation");
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
  if (is_true(any(censor_timen < tstopn))) {
    stop("censor_time must be greater than or equal to tstop");
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
  
  StringVector covariates_lgs(q+p2+ns_df+1);
  NumericMatrix z_lgsn(n,q+p2);
  covariates_lgs[0] = "counterfactual";
  if (strata_main_effect_only) {
    int k = 0;
    for (int i=0; i<p_stratum; ++i) {
      SEXP col_level = levels[i];
      SEXPTYPE type_level = TYPEOF(col_level);
      
      int di = d[i]-1;
      for (int j=0; j<di; ++j) {
        covariates_lgs[k+j+1] = as<std::string>(stratum[i]);
        if (type_level == STRSXP) {
          StringVector u = col_level;
          std::string label = sanitize(as<std::string>(u[j]));
          covariates_lgs[k+j+1] += label;
        } else if (type_level == REALSXP) {
          NumericVector u = col_level;
          covariates_lgs[k+j+1] += std::to_string(u[j]);
        } else if (type_level == INTSXP || type_level == LGLSXP) {
          IntegerVector u = col_level;
          covariates_lgs[k+j+1] += std::to_string(u[j]);
        }
        z_lgsn(_,k+j) = (stratan(_,i) == j);
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
      covariates_lgs[j+1] = "";
      for (int i=0; i<p_stratum; ++i) {
        SEXP col_level = levels[i];
        SEXPTYPE type_level = TYPEOF(col_level);
        
        IntegerVector q_col = stratan(_,i);
        int l = q_col[first_k];
        
        covariates_lgs[j+1] += as<std::string>(stratum[i]);
        
        if (type_level == STRSXP) {
          StringVector u = col_level;
          std::string label = sanitize(as<std::string>(u[l]));
          covariates_lgs[j+1] += label;
        } else if (type_level == REALSXP) {
          NumericVector u = col_level;
          covariates_lgs[j+1] += std::to_string(u[l]);
        } else if (type_level == INTSXP || type_level == LGLSXP) {
          IntegerVector u = col_level;
          covariates_lgs[j+1] += std::to_string(u[l]);
        }
        
        if (i < p_stratum-1) {
          covariates_lgs[j+1] += ".";
        }
      }
      z_lgsn(_,j) = (stratumn == j);
    }
  }
  
  for (int j=0; j<p2; ++j) {
    String zj = conf_cov[j];
    
    if (!hasVariable(data, zj)) {
      stop("data must contain the variables in conf_cov");
    }
    
    if (zj == treat) {
      stop("treat should be excluded from conf_cov");
    }
    
    NumericVector u = data[zj];
    covariates_lgs[q+j+1] = zj;
    z_lgsn(_,q+j) = u;
  }
  
  if (ns_df < 0) {
    stop("ns_df must be a nonnegative integer");
  }
  
  for (int j=0; j<ns_df; ++j) {
    covariates_lgs[q+p2+j+1] = "ns" + std::to_string(j+1);
  }
  
  
  if (low_psi >= hi_psi) {
    stop("low_psi must be less than hi_psi");
  }
  
  if (n_eval_z < 2) {
    stop("n_eval_z must be greater than or equal to 2");
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
        std::isnan(tstartn[i]) || std::isnan(tstopn[i]) || 
        eventn[i] == NA_INTEGER || treatn[i] == NA_INTEGER || 
        std::isnan(censor_timen[i]) || pdn[i] == NA_INTEGER ||
        swtrtn[i] == NA_INTEGER) {
      sub[i] = 0;
    }
    for (int j=0; j<p; ++j) {
      if (std::isnan(zn(i,j))) sub[i] = 0;
    }
  }
  
  IntegerVector order = which(sub);
  idn = idn[order];
  stratumn = stratumn[order];
  tstartn = tstartn[order];
  tstopn = tstopn[order];
  eventn = eventn[order];
  treatn = treatn[order];
  censor_timen = censor_timen[order];
  pdn = pdn[order];
  pd_timen = pd_timen[order];
  swtrtn = swtrtn[order];
  swtrt_timen = swtrt_timen[order];
  if (p > 0) zn = subset_matrix_by_row(zn, order);
  z_lgsn = subset_matrix_by_row(z_lgsn, order);
  n = sum(sub);
  if (n == 0) stop("no observations left after removing missing values");
  
  // split at treatment switching into two observations if treatment 
  // switching occurs strictly between tstart and tstop for a subject
  LogicalVector tosplit(n);
  for (int i=0; i<n; ++i) {
    tosplit[i] = swtrtn[i] == 1 && swtrt_timen[i] > tstartn[i] && 
      swtrt_timen[i] < tstopn[i] ? 1 : 0;
  }
  
  k = sum(tosplit);
  if (k > 0) {
    // copy old matrices to new matrices
    NumericMatrix z1(n+k, zn.ncol());
    NumericMatrix z_lgsn1(n+k, z_lgsn.ncol());
    for (int i=0; i<n; ++i) {
      z1(i,_) = zn(i,_);
      z_lgsn1(i,_) = z_lgsn(i,_);
    }
    
    IntegerVector sub = which(tosplit);
    for (int i=0; i<k; ++i) {
      // append a new observation by changing tstart
      int l = sub[i];
      idn.push_back(idn[l]);
      stratumn.push_back(stratumn[l]);
      tstartn.push_back(swtrt_timen[l]);
      tstopn.push_back(tstopn[l]);
      eventn.push_back(eventn[l]);
      treatn.push_back(treatn[l]);
      censor_timen.push_back(censor_timen[l]);
      pdn.push_back(pdn[l]);
      pd_timen.push_back(pd_timen[l]);
      swtrtn.push_back(swtrtn[l]);
      swtrt_timen.push_back(swtrt_timen[l]);
      z1(n+i,_) = zn(l,_);
      z_lgsn1(n+i,_) = z_lgsn(l,_);
      
      // change tstop and event for the old observation
      tstopn[l] = swtrt_timen[l];
      eventn[l] = 0;
    }
    
    // update number of rows and old matrices
    n = n + k;
    zn = z1;
    z_lgsn = z_lgsn1;
  }
  
  // sort data by treatment group, id, and time
  order = seq(0, n-1);
  std::sort(order.begin(), order.end(), [&](int i, int j) {
    return std::tie(treatn[i], idn[i], tstopn[i]) <
      std::tie(treatn[j], idn[j], tstopn[j]);
  });
  
  idn = idn[order];
  stratumn = stratumn[order];
  tstartn = tstartn[order];
  tstopn = tstopn[order];
  eventn = eventn[order];
  treatn = treatn[order];
  censor_timen = censor_timen[order];
  pdn = pdn[order];
  pd_timen = pd_timen[order];
  swtrtn = swtrtn[order];
  swtrt_timen = swtrt_timen[order];
  zn = subset_matrix_by_row(zn, order);
  z_lgsn = subset_matrix_by_row(z_lgsn, order);
  
  IntegerVector idx(1,0); // first observation within an id
  for (int i=1; i<n; ++i) {
    if (idn[i] != idn[i-1]) {
      idx.push_back(i);
    }
  }
  
  int nids = static_cast<int>(idx.size());
  idx.push_back(n);
  
  IntegerVector idx1(nids); // last observation within an id
  for (int i=0; i<nids; ++i) {
    idx1[i] = idx[i+1]-1;
  }
  
  IntegerVector osn(n);
  NumericVector os_timen(n);
  for (int i=0; i<nids; ++i) {
    k = idx1[i];
    for (int j=idx[i]; j<idx[i+1]; ++j) {
      osn[j] = eventn[k];
      os_timen[j] = tstopn[k];
    }
  }
  
  if (is_true(any(ifelse(pdn == 1, pd_timen > os_timen, 0)))) {
    stop("pd_time must be less than or equal to os_time");
  }
  if (is_true(any(ifelse(swtrtn == 1, swtrt_timen > os_timen, 0)))) {
    stop("swtrt_time must be less than or equal to os_time");
  }
  
  if (!admin_recensor_only) { // use the actual censoring time for dropouts
    for (int i=0; i<nids; ++i) {
      if (eventn[idx1[i]] == 0) {
        for (int j=idx[i]; j<idx[i+1]; ++j) {
          censor_timen[j] = tstopn[idx1[i]];
        }
      }
    }
  }
  
  IntegerVector stratum1 = stratumn[idx1];
  IntegerVector treat1 = treatn[idx1];
  NumericVector tstopn1 = tstopn[idx1];
  IntegerVector event1 = eventn[idx1];
  
  DataFrame dt = DataFrame::create(
    Named("stratum") = stratum1,
    Named("treat") = treat1,
    Named("time") = tstopn1,
    Named("event") = event1);
  
  DataFrame lr = lrtest(dt,"","stratum","treat","time","","event","",0,0,0);
  double logRankPValue = lr["logRankPValue"];
  double zcrit = R::qnorm(1-alpha/2, 0, 1, 1, 0);
  
  double step_psi = (hi_psi - low_psi)/(n_eval_z - 1);
  NumericVector psi(n_eval_z);
  for (int i=0; i<n_eval_z; ++i) {
    psi[i] = low_psi + i*step_psi;
  }
  
  k = -1; // indicate the observed data
  auto f = [&k, has_stratum, stratum, p_stratum, u_stratum, 
            type_treat, treat, treatwi, treatwn, treatwc, 
            type_id, id, idwi, idwn, idwc,
            q, p, p2, covariates, covariates_lgs, ns_df, firth, flic, 
            low_psi, hi_psi, n_eval_z, psi, recensor, swtrt_control_only, 
            gridsearch, rooting, alpha, zcrit, ties, tol, offset] (
                IntegerVector& idb, IntegerVector& stratumb, 
                NumericVector& tstartb, NumericVector& tstopb, 
                IntegerVector& eventb, IntegerVector& treatb, 
                IntegerVector& osb, NumericVector& os_timeb, 
                NumericVector& censor_timeb, 
                IntegerVector& pdb, NumericVector& pd_timeb, 
                IntegerVector& swtrtb, NumericVector& swtrt_timeb, 
                NumericMatrix& zb, NumericMatrix& z_lgsb)->List {
                  int n = static_cast<int>(idb.size());
                  bool fail = false; // whether any model fails to converge
                  NumericVector init(1, NA_REAL);
                  
                  IntegerVector idx(1,0); // first observation within an id
                  for (int i=1; i<n; ++i) {
                    if (idb[i] != idb[i-1]) {
                      idx.push_back(i);
                    }
                  }
                  
                  int nids = static_cast<int>(idx.size());
                  idx.push_back(n);
                  
                  IntegerVector idx1(nids); // last observation within an id
                  for (int i=0; i<nids; ++i) {
                    idx1[i] = idx[i+1]-1;
                  }
                  
                  // one observation per subject
                  IntegerVector id1 = idb[idx1];
                  IntegerVector stratum1 = stratumb[idx1];
                  NumericVector time1 = tstopb[idx1];
                  IntegerVector event1 = eventb[idx1];
                  IntegerVector treat1 = treatb[idx1];
                  NumericVector censor_time1 = censor_timeb[idx1];
                  IntegerVector swtrt1 = swtrtb[idx1];
                  NumericVector swtrt_time1 = swtrt_timeb[idx1];
                  NumericMatrix z1 = subset_matrix_by_row(zb, idx1);
                  
                  // time and event adjusted for treatment switching
                  NumericVector t_star = clone(time1);
                  IntegerVector d_star = clone(event1);
                  
                  double psi0hat = NA_REAL;
                  double psi0lower = NA_REAL, psi0upper = NA_REAL;
                  double psi1hat = NA_REAL;
                  double psi1lower = NA_REAL, psi1upper = NA_REAL;
                  NumericVector psi0hat_vec, psi1hat_vec;
                  String psi_CI_type;
                  
                  // initialize data_switch, km_switch, eval_z, 
                  // data_nullcox, fit_nullcox, data_logis, fit_logis
                  List data_switch(2), km_switch(2), eval_z(2);
                  List data_nullcox(2), fit_nullcox(2);
                  List data_logis(2), fit_logis(2);
                  if (k == -1) {
                    for (int h=0; h<2; ++h) {
                      List data_x = List::create(
                        Named("data") = R_NilValue,
                        Named(treat) = R_NilValue
                      );
                      
                      if (type_treat == LGLSXP || type_treat == INTSXP) {
                        data_x[treat] = treatwi[1-h];
                      } else if (type_treat == REALSXP) {
                        data_x[treat] = treatwn[1-h];
                      } else if (type_treat == STRSXP) {
                        data_x[treat] = treatwc[1-h];
                      }
                      
                      data_switch[h] = data_x;
                      km_switch[h] = data_x;
                      eval_z[h] = data_x;
                      data_nullcox[h] = data_x;
                      data_logis[h] = data_x;
                      
                      List fit_x = List::create(
                        Named("fit") = R_NilValue,
                        Named(treat) = R_NilValue
                      );
                      
                      if (type_treat == LGLSXP || type_treat == INTSXP) {
                        fit_x[treat] = treatwi[1-h];
                      } else if (type_treat == REALSXP) {
                        fit_x[treat] = treatwn[1-h];
                      } else if (type_treat == STRSXP) {
                        fit_x[treat] = treatwc[1-h];
                      }
                      
                      fit_nullcox[h] = fit_x;
                      fit_logis[h] = fit_x;
                    }
                  }
                  
                  List data_outcome;
                  List fit_outcome;
                  double hrhat = NA_REAL, hrlower = NA_REAL, 
                    hrupper = NA_REAL, pvalue = NA_REAL;
                  
                  bool psimissing = false;
                  
                  // treat arms that include patients who switched treatment
                  IntegerVector treats(1);
                  treats[0] = 0;
                  if (!swtrt_control_only) {
                    treats.push_back(1);
                  }
                  
                  int K = static_cast<int>(treats.size());
                  for (int h=0; h<K; ++h) {
                    // post progression data up to switching for the treat
                    LogicalVector c1 = (treatb == h);
                    LogicalVector c2 = ifelse(pdb==1, tstopb >= pd_timeb, 0);
                    LogicalVector c3 = ifelse(
                      swtrtb == 1, tstartb < swtrt_timeb, tstopb < os_timeb);
                    IntegerVector l = which(c1 & c2 & c3);
                    
                    IntegerVector id2 = idb[l];
                    IntegerVector stratum2 = stratumb[l];
                    NumericVector tstart2 = tstartb[l];
                    NumericVector tstop2 = tstopb[l];
                    NumericVector pd_time2 = pd_timeb[l];
                    IntegerVector os2 = osb[l];
                    NumericVector os_time2 = os_timeb[l];
                    NumericVector censor_time2 = censor_timeb[l];
                    IntegerVector swtrt2 = swtrtb[l];
                    NumericVector swtrt_time2 = swtrt_timeb[l];
                    NumericMatrix z_lgs2 = subset_matrix_by_row(z_lgsb, l);
                    
                    int n2 = static_cast<int>(l.size());
                    
                    // treatment switching indicators
                    IntegerVector cross2(n2);
                    for (int i=0; i<n2; ++i) {
                      cross2[i] = swtrt2[i] && tstop2[i] >= swtrt_time2[i];
                    }
                    
                    // obtain natural cubic spline knots
                    NumericMatrix s(n2, ns_df);
                    if (ns_df > 0) {
                      NumericVector x = tstop2[cross2 == 1];
                      NumericVector knots(1, NA_REAL);
                      NumericVector boundary_knots(1, NA_REAL);
                      s = nscpp(x, ns_df, knots, 0, boundary_knots);
                      knots = s.attr("knots");
                      boundary_knots = s.attr("boundary_knots");
                      s = nscpp(tstop2, NA_INTEGER, knots, 0, 
                                boundary_knots);
                    }
                    
                    // re-baseline based on disease progression date
                    IntegerVector idx2(1,0);
                    for (int i=1; i<n2; ++i) {
                      if (id2[i] != id2[i-1]) {
                        idx2.push_back(i);
                      }
                    }
                    
                    int nids2 = static_cast<int>(idx2.size());
                    idx2.push_back(n2);
                    
                    IntegerVector id3(nids2);
                    IntegerVector stratum3(nids2);
                    IntegerVector os3(nids2);
                    NumericVector os_time3(nids2);
                    NumericVector censor_time3(nids2);
                    IntegerVector swtrt3(nids2);
                    NumericVector swtrt_time3(nids2, NA_REAL);
                    for (int i=0; i<nids2; ++i) {
                      int j = idx2[i];
                      double b2 = pd_time2[j] - offset;
                      id3[i] = id2[j];
                      stratum3[i] = stratum2[j];
                      os3[i] = os2[j];
                      os_time3[i] = os_time2[j] - b2;
                      censor_time3[i] = censor_time2[j] - b2;
                      swtrt3[i] = swtrt2[j];
                      if (swtrt3[i] == 1) {
                        swtrt_time3[i] = swtrt_time2[j] - b2;
                      }
                    }
                    
                    // obtain the estimate and confidence interval of psi
                    double psihat = NA_REAL;
                    double psilower = NA_REAL, psiupper = NA_REAL;
                    NumericVector psihat_vec;
                    NumericVector Z(n_eval_z, NA_REAL);
                    if (gridsearch) {
                      for (int i=0; i<n_eval_z; ++i) {
                        List out = est_psi_tsegest(
                          n2, q, p2, nids2, idx2, stratum3, 
                          os3, os_time3, censor_time3, 
                          swtrt3, swtrt_time3, id2, cross2, 
                          tstart2, tstop2, covariates_lgs, 
                          z_lgs2, ns_df, s, firth, flic, 
                          recensor, alpha, ties, offset, psi[i]);
                        
                        bool fail = out["fail"];
                        if (!fail) Z[i] = out["z_counterfactual"];
                      }
                      
                      List psihat_list = getpsiest(0, psi, Z, 0);
                      psihat = psihat_list["root"];
                      psihat_vec = psihat_list["roots"];
                      psi_CI_type = "grid search";
                      
                      if (k == -1) {
                        List psilower_list = getpsiest(zcrit, psi, Z, -1);
                        psilower = psilower_list["root"];
                        List psiupper_list = getpsiest(-zcrit, psi, Z, 1);
                        psiupper = psiupper_list["root"];
                      }
                    } else {
                      // z-stat for the slope of counterfactual survival time
                      // martingale residuals in logistic regression model
                      double target = 0;
                      auto g = [&target, n2, q, p2, nids2, idx2, stratum3, 
                                os3, os_time3, censor_time3, swtrt3, 
                                swtrt_time3, id2, cross2, tstart2, tstop2, 
                                covariates_lgs, z_lgs2, ns_df, s, 
                                firth, flic, recensor, alpha, ties, 
                                offset](double x)->double{
                                  List out = est_psi_tsegest(
                                    n2, q, p2, nids2, idx2, stratum3, 
                                    os3, os_time3, censor_time3, 
                                    swtrt3, swtrt_time3, id2, cross2, 
                                    tstart2, tstop2, covariates_lgs, 
                                    z_lgs2, ns_df, s, firth, flic, 
                                    recensor, alpha, ties, offset, x);
                                  
                                  bool fail = out["fail"];
                                  if (!fail) {
                                    double z = out["z_counterfactual"];
                                    return z - target;
                                  } else {
                                    return NA_REAL;
                                  }
                                };
                      
                      // causal parameter estimates
                      double psilo = getpsiend(g, 1, low_psi);
                      double psihi = getpsiend(g, 0, hi_psi);
                      if (!std::isnan(psilo) && !std::isnan(psihi)) {
                        if (rooting == "brent") {
                          psihat = brent(g, psilo, psihi, tol);
                        } else {
                          psihat = bisect(g, psilo, psihi, tol);
                        }
                      }
                      psi_CI_type = "root finding";

                      if (k == -1) {
                        target = zcrit;
                        psilo = getpsiend(g, 1, low_psi);
                        psihi = getpsiend(g, 0, hi_psi);
                        if (!std::isnan(psilo) && !std::isnan(psihi)) {
                          if (!std::isnan(psihat) && g(psihat) < 0) {
                            if (rooting == "brent") {
                              psilower = brent(g, psilo, psihat, tol);
                            } else {
                              psilower = bisect(g, psilo, psihat, tol);
                            }
                          } else {
                            if (rooting == "brent") {
                              psilower = brent(g, psilo, psihi, tol);
                            } else {
                              psilower = bisect(g, psilo, psihi, tol);
                            }
                          }
                        }
                        
                        target = -zcrit;
                        psilo = getpsiend(g, 1, low_psi);
                        psihi = getpsiend(g, 0, hi_psi);
                        if (!std::isnan(psilo) && !std::isnan(psihi)) {
                          if (!std::isnan(psihat) && g(psihat) > 0) {
                            if (rooting == "brent") {
                              psiupper = brent(g, psihat, psihi, tol);  
                            } else {
                              psiupper = bisect(g, psihat, psihi, tol);  
                            }
                          } else {
                            if (rooting == "brent") {
                              psiupper = brent(g, psilo, psihi, tol);  
                            } else {
                              psiupper = bisect(g, psilo, psihi, tol);
                            }
                          }
                        }
                      }
                    }
                    
 
                    if (!std::isnan(psihat)) {
                      // calculate counter-factual survival times
                      double a = exp(psihat);
                      double c0 = std::min(1.0, a);
                      for (int i=0; i<nids; ++i) {
                        if (treat1[i] == h) {
                          double b2, u_star, c_star;
                          if (swtrt1[i] == 1) {
                            b2 = swtrt_time1[i] - offset;
                            u_star = b2 + (time1[i] - b2)*a;
                          } else {
                            u_star = time1[i];
                          }
                          
                          if (recensor) {
                            c_star = censor_time1[i]*c0;
                            t_star[i] = std::min(u_star, c_star);
                            d_star[i] = c_star < u_star ? 0 : event1[i];
                          } else {
                            t_star[i] = u_star;
                            d_star[i] = event1[i];
                          }
                        }
                      }
                      
                      if (k == -1) {
                        // obtain data and KM plot for time to switch
                        NumericVector swtrt_timen4(nids2);
                        for (int i=0; i<nids2; ++i) {
                          swtrt_timen4[i] = swtrt3[i] == 1 ? 
                          swtrt_time3[i] : os_time3[i];
                        }
                        
                        List data1 = List::create(
                          Named("swtrt") = swtrt3,
                          Named("swtrt_time") = swtrt_timen4);
                        
                        if (type_id == INTSXP) {
                          data1.push_front(idwi[id3], id);
                        } else if (type_id == REALSXP) {
                          data1.push_front(idwn[id3], id);
                        } else if (type_id == STRSXP) {
                          data1.push_front(idwc[id3], id);
                        }
                        
                        if (has_stratum) {
                          for (int i=0; i<p_stratum; ++i) {
                            std::string s = as<std::string>(stratum[i]);
                            SEXP col_stratum = u_stratum[s];
                            SEXPTYPE type_stratum = TYPEOF(col_stratum);
                            if (type_stratum == INTSXP) {
                              IntegerVector v = col_stratum;
                              data1.push_back(v[stratum3], s);
                            } else if (type_stratum == REALSXP) {
                              NumericVector v = col_stratum;
                              data1.push_back(v[stratum3], s);
                            } else if (type_stratum == STRSXP) {
                              StringVector v = col_stratum;
                              data1.push_back(v[stratum3], s);
                            }
                          }
                        }
                        
                        List km1 = kmest(data1,"","","swtrt_time","",
                                         "swtrt","","log-log",1-alpha,1);
                        
                        // obtain the Wald statistics for the coefficient of 
                        // the counterfactual in the logistic regression 
                        // switching model at a sequence of psi values
                        if (!gridsearch) {
                          for (int i=0; i<n_eval_z; ++i) {
                            List out = est_psi_tsegest(
                              n2, q, p2, nids2, idx2, stratum3, 
                              os3, os_time3, censor_time3, 
                              swtrt3, swtrt_time3, id2, cross2, 
                              tstart2, tstop2, covariates_lgs, 
                              z_lgs2, ns_df, s, firth, flic, 
                              recensor, alpha, ties, offset, psi[i]);
                            
                            bool fail = out["fail"];
                            if (!fail) Z[i] = out["z_counterfactual"];
                          }
                          
                          List psihat_list = getpsiest(0, psi, Z, 0);
                          psihat_vec = psihat_list["roots"];
                        }
                        
                        List data2 = List::create(
                          Named("psi") = psi,
                          Named("Z") = Z);
                        
                        // obtain data and fit for null Cox & logistic models
                        List out = est_psi_tsegest(
                          n2, q, p2, nids2, idx2, stratum3, 
                          os3, os_time3, censor_time3, 
                          swtrt3, swtrt_time3, id2, cross2, 
                          tstart2, tstop2, covariates_lgs, 
                          z_lgs2, ns_df, s, firth, flic, 
                          recensor, alpha, ties, offset, psihat);
                        
                        bool fail_lgs = out["fail"];
                        if (fail_lgs) fail = true;
                        
                        List data3 = out["data_nullcox"];
                        List data4 = out["data_logis"];
                        List fit3 = out["fit_nullcox"];
                        List fit4 = out["fit_logis"];
                        
                        if (type_id == INTSXP) {
                          data3.push_front(idwi[id3], id);
                          data4.push_front(idwi[id2], id);
                        } else if (type_id == REALSXP) {
                          data3.push_front(idwn[id3], id);
                          data4.push_front(idwn[id2], id);
                        } else if (type_id == STRSXP) {
                          data3.push_front(idwc[id3], id);
                          data4.push_front(idwc[id2], id);
                        }
                        
                        if (has_stratum) {
                          for (int i=0; i<p_stratum; ++i) {
                            std::string s = as<std::string>(stratum[i]);
                            SEXP col_stratum = u_stratum[s];
                            SEXPTYPE type_stratum = TYPEOF(col_stratum);
                            if (type_stratum == INTSXP) {
                              IntegerVector v = col_stratum;
                              data3.push_back(v[stratum3], s);
                              data4.push_back(v[stratum2], s);
                            } else if (type_stratum == REALSXP) {
                              NumericVector v = col_stratum;
                              data3.push_back(v[stratum3], s);
                              data4.push_back(v[stratum2], s);
                            } else if (type_stratum == STRSXP) {
                              StringVector v = col_stratum;
                              data3.push_back(v[stratum3], s);
                              data4.push_back(v[stratum2], s);
                            }
                          }
                        }
                        
                        // update the data and model fits
                        List data_1 = data_switch[h];
                        data_1["data"] = as<DataFrame>(data1);
                        data_switch[h] = data_1;
                        
                        List data_2 = clone(data_1);
                        data_2["data"] = as<DataFrame>(data2);
                        eval_z[h] = data_2;
                        
                        List data_3 = clone(data_1);
                        data_3["data"] = as<DataFrame>(data3);
                        data_nullcox[h] = data_3;
                        
                        List data_4 = clone(data_1);
                        data_4["data"] = as<DataFrame>(data4);
                        data_logis[h] = data_4;
                        
                        List data_5 = clone(data_1);
                        data_5["data"] = as<DataFrame>(km1);
                        km_switch[h] = data_5;
                        
                        List fit_3 = fit_nullcox[h];
                        fit_3["fit"] = fit3;
                        fit_nullcox[h] = fit_3;
                        
                        List fit_4 = clone(fit_3);
                        fit_4["fit"] = fit4;
                        fit_logis[h] = fit_4;
                      }
                    } else {
                      psimissing = 1;
                    }
                    
                    
                    // update treatment-specific causal parameter estimates
                    if (h == 0) {
                      psi0hat = psihat;
                      psi0hat_vec = psihat_vec;
                      if (k == -1) {
                        psi0lower = psilower;
                        psi0upper = psiupper;
                      }
                    } else {
                      psi1hat = psihat;
                      psi1hat_vec = psihat_vec;
                      if (k == -1) {
                        psi1lower = psilower;
                        psi1upper = psiupper;
                      }
                    }
                  }
                  
                  if (!psimissing) {
                    // Cox model for hypothetical treatment effect estimate
                    data_outcome = List::create(
                      Named("uid") = id1,
                      Named("t_star") = t_star,
                      Named("d_star") = d_star,
                      Named("treated") = treat1);
                    
                    data_outcome.push_back(stratum1, "ustratum");
                    
                    for (int j=0; j<p; ++j) {
                      String zj = covariates[j+1];
                      NumericVector u = z1(_,j);
                      data_outcome.push_back(u, zj);
                    }
                    
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
                      Named("data_switch") = data_switch,
                      Named("km_switch") = km_switch,
                      Named("eval_z") = eval_z,
                      Named("data_nullcox") = data_nullcox,
                      Named("fit_nullcox") = fit_nullcox,
                      Named("data_logis") = data_logis,
                      Named("fit_logis") = fit_logis,
                      Named("data_outcome") = data_outcome,
                      Named("fit_outcome") = fit_outcome,
                      Named("psihat") = psi0hat,
                      Named("psihat_vec") = psi0hat_vec,
                      Named("psilower") = psi0lower,
                      Named("psiupper") = psi0upper,
                      Named("psi1hat") = psi1hat,
                      Named("psi1hat_vec") = psi1hat_vec,
                      Named("psi1lower") = psi1lower,
                      Named("psi1upper") = psi1upper,
                      Named("psi_CI_type") = psi_CI_type,
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
  
  List out = f(idn, stratumn, tstartn, tstopn, eventn, treatn,
               osn, os_timen, censor_timen, pdn, pd_timen, 
               swtrtn, swtrt_timen, zn, z_lgsn);
  
  List data_switch = out["data_switch"];
  List km_switch = out["km_switch"];
  List eval_z = out["eval_z"];
  List data_nullcox = out["data_nullcox"];
  List fit_nullcox = out["fit_nullcox"];
  List data_logis = out["data_logis"];
  List fit_logis = out["fit_logis"];
  List data_outcome = out["data_outcome"];
  List fit_outcome = out["fit_outcome"];
  
  double psihat = out["psihat"];
  NumericVector psihat_vec = out["psihat_vec"];
  double psilower = out["psilower"];
  double psiupper = out["psiupper"];
  double psi1hat = out["psi1hat"];
  NumericVector psi1hat_vec = out["psi1hat_vec"];
  double psi1lower = out["psi1lower"];
  double psi1upper = out["psi1upper"];
  String psi_CI_type = out["psi_CI_type"];
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
    
    IntegerVector treated = data_outcome["treated"];
    if (type_treat == LGLSXP || type_treat == INTSXP) {
      data_outcome.push_back(treatwi[1-treated], treat);
    } else if (type_treat == REALSXP) {
      data_outcome.push_back(treatwn[1-treated], treat);
    } else if (type_treat == STRSXP) {
      data_outcome.push_back(treatwc[1-treated], treat);
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
      
      IntegerVector nobs = diff(idx);
      int N = max(nobs)*nids;
      
      int B = N*n_boot;
      IntegerVector boot_indexc(B);
      IntegerVector oidc(B);
      IntegerVector idc(B), stratumc(B), treatc(B), eventc(B);
      IntegerVector osc(B), pdc(B), swtrtc(B);
      NumericVector tstartc(B), tstopc(B), os_timec(B), censor_timec(B);
      NumericVector pd_timec(B), swtrt_timec(B);
      NumericMatrix z_lgsc(B,q+p2);
      int index1 = 0; // current index for fail boots data
      
      if (has_stratum) {
        // sort data by treatment group, stratum, id, and time
        IntegerVector order = seq(0, n-1);
        std::sort(order.begin(), order.end(), [&](int i, int j) {
          return std::tie(treatn[i], stratumn[i], idn[i], tstopn[i]) <
            std::tie(treatn[j], stratumn[j], idn[j], tstopn[j]);
        });
        
        idn = idn[order];
        stratumn = stratumn[order];
        tstartn = tstartn[order];
        tstopn = tstopn[order];
        eventn = eventn[order];
        treatn = treatn[order];
        censor_timen = censor_timen[order];
        pdn = pdn[order];
        pd_timen = pd_timen[order];
        swtrtn = swtrtn[order];
        swtrt_timen = swtrt_timen[order];
        zn = subset_matrix_by_row(zn, order);
        z_lgsn = subset_matrix_by_row(z_lgsn, order);
      }
      
      IntegerVector idx(1,0); // first observation within an id
      for (int i=1; i<n; ++i) {
        if (idn[i] != idn[i-1]) {
          idx.push_back(i);
        }
      }
      
      int nids = static_cast<int>(idx.size());
      idx.push_back(n);
      
      IntegerVector idx1(nids); // last observation within an id
      for (int i=0; i<nids; ++i) {
        idx1[i] = idx[i+1]-1;
      }
      
      IntegerVector treat1 = treatn[idx1];
      IntegerVector stratum1 = stratumn[idx1];
      
      IntegerVector tsx(1,0); // first id within each treat/stratum
      for (int i=1; i<nids; ++i) {
        if (treat1[i] != treat1[i-1] || stratum1[i] != stratum1[i-1]) {
          tsx.push_back(i);
        }
      }
      
      int ntss = static_cast<int>(tsx.size());
      tsx.push_back(nids); // add the end index
      
      for (k=0; k<n_boot; ++k) {
        IntegerVector oidb(N);
        IntegerVector idb(N), stratumb(N), treatb(N), eventb(N);
        IntegerVector osb(N), pdb(N), swtrtb(N);
        NumericVector tstartb(N), tstopb(N), os_timeb(N), censor_timeb(N);
        NumericVector pd_timeb(N), swtrt_timeb(N);
        NumericMatrix zb(N,p), z_lgsb(N,q+p2);
        
        // sample the subject-level data with replacement by treat/stratum
        int l = 0; // current index for bootstrap data
        for (int h=0; h<ntss; ++h) {
          for (int r=tsx[h]; r<tsx[h+1]; ++r) {
            double u = R::runif(0,1);
            int i = tsx[h] + static_cast<int>(std::floor(u*(tsx[h+1]-tsx[h])));
            
            // create unique ids for bootstrap data sets
            int oidb1 = idn[idx[i]];
            int idb1 = oidb1 + r*nids;
            
            for (int j=idx[i]; j<idx[i+1]; ++j) {
              oidb[l] = oidb1;
              idb[l] = idb1;
              stratumb[l] = stratumn[j];
              tstartb[l] = tstartn[j];
              tstopb[l] = tstopn[j];
              eventb[l] = eventn[j];
              treatb[l] = treatn[j];
              osb[l] = osn[j];
              os_timeb[l] = os_timen[j];
              censor_timeb[l] = censor_timen[j];
              pdb[l] = pdn[j];
              pd_timeb[l] = pd_timen[j];
              swtrtb[l] = swtrtn[j];
              swtrt_timeb[l] = swtrt_timen[j];
              zb(l,_) = zn(j,_);
              z_lgsb(l,_) = z_lgsn(j,_);
              ++l;
            }
          }
        }
        
        IntegerVector sub = Range(0,l-1);
        oidb = oidb[sub];
        idb = idb[sub];
        stratumb = stratumb[sub];
        tstartb = tstartb[sub];
        tstopb = tstopb[sub];
        eventb = eventb[sub];
        treatb = treatb[sub];
        osb = osb[sub];
        os_timeb = os_timeb[sub];
        censor_timeb = censor_timeb[sub];
        pdb = pdb[sub];
        pd_timeb = pd_timeb[sub];
        swtrtb = swtrtb[sub];
        swtrt_timeb = swtrt_timeb[sub];
        zb = subset_matrix_by_row(zb, sub);
        z_lgsb = subset_matrix_by_row(z_lgsb, sub);
        
        out = f(idb, stratumb, tstartb, tstopb, eventb, treatb,
                osb, os_timeb, censor_timeb, pdb, pd_timeb, 
                swtrtb, swtrt_timeb, zb, z_lgsb);
        
        fails[k] = out["fail"];
        hrhats[k] = out["hrhat"];
        psihats[k] = out["psihat"];
        psi1hats[k] = out["psi1hat"];
        
        if (fails[k]) {
          for (int i=0; i<l; ++i) {
            int j = index1 + i;
            boot_indexc[j] = k+1;
            oidc[j] = oidb[i];
            idc[j] = idb[i];
            stratumc[j] = stratumb[i];
            tstartc[j] = tstartb[i];
            tstopc[j] = tstopb[i];
            eventc[j] = eventb[i];
            treatc[j] = treatb[i];
            osc[j] = osb[i];
            os_timec[j] = os_timeb[i];
            censor_timec[j] = censor_timeb[i];
            pdc[j] = pdb[i];
            pd_timec[j] = pd_timeb[i];
            swtrtc[j] = swtrtb[i];
            swtrt_timec[j] = swtrt_timeb[i];
            z_lgsc(j,_) = z_lgsb(i,_);
          }
          index1 += l;
        }
      }
      
      if (is_true(any(fails))) {
        IntegerVector sub = seq(0,index1-1);
        boot_indexc = boot_indexc[sub];
        oidc = oidc[sub];
        idc = idc[sub];
        stratumc = stratumc[sub];
        tstartc = tstartc[sub];
        tstopc = tstopc[sub];
        eventc = eventc[sub];
        treatc = treatc[sub];
        osc = osc[sub];
        os_timec = os_timec[sub];
        censor_timec = censor_timec[sub];
        pdc = pdc[sub];
        pd_timec = pd_timec[sub];
        swtrtc = swtrtc[sub];
        swtrt_timec = swtrt_timec[sub];
        z_lgsc = subset_matrix_by_row(z_lgsc,sub);
        
        fail_boots_data = List::create(
          Named("boot_index") = boot_indexc,
          Named("uid") = idc,
          Named("tstart") = tstartc,
          Named("tstop") = tstopc,
          Named("event") = eventc,
          Named("treated") = treatc,
          Named("os") = osc,
          Named("os_time") = os_timec,
          Named("censor_time") = censor_timec,
          Named("pd") = pdc,
          Named("pd_time") = pd_timec,
          Named("swtrt") = swtrtc,
          Named("swtrt_time") = swtrt_timec
        );
        
        for (int j=0; j<q+p2; ++j) {
          String zj = covariates_lgs[j+1];
          NumericVector u = z_lgsc(_,j);
          fail_boots_data.push_back(u, zj);
        }
        
        if (type_id == INTSXP) {
          fail_boots_data.push_back(idwi[oidc], id);
        } else if (type_id == REALSXP) {
          fail_boots_data.push_back(idwn[oidc], id);
        } else if (type_id == STRSXP) {
          fail_boots_data.push_back(idwc[oidc], id);
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
  
  List settings = List::create(
    Named("strata_main_effect_only") = strata_main_effect_only,
    Named("ns_df") = ns_df,
    Named("firth") = firth,
    Named("flic") = flic,
    Named("low_psi") = low_psi,
    Named("hi_psi") = hi_psi,
    Named("n_eval_z") = n_eval_z,
    Named("recensor") = recensor,
    Named("admin_recensor_only") = admin_recensor_only,
    Named("swtrt_control_only") = swtrt_control_only,
    Named("gridsearch") = gridsearch,
    Named("root_finding") = root_finding,
    Named("alpha") = alpha,
    Named("ties") = ties,
    Named("tol") = tol,
    Named("offset") = offset,
    Named("boot") = boot,
    Named("n_boot") = n_boot,
    Named("seed") = seed);
  
  List result = List::create(
    Named("psi") = psihat,
    Named("psi_roots") = psihat_vec,
    Named("psi_CI") = NumericVector::create(psilower, psiupper),
    Named("psi_CI_type") = psi_CI_type,
    Named("logrank_pvalue") = logRankPValue,
    Named("cox_pvalue") = pvalue,
    Named("hr") = hrhat,
    Named("hr_CI") = NumericVector::create(hrlower, hrupper),
    Named("hr_CI_type") = hr_CI_type,
    Named("data_switch") = data_switch,
    Named("km_switch") = km_switch,
    Named("eval_z") = eval_z,
    Named("data_nullcox") = data_nullcox,
    Named("fit_nullcox") = fit_nullcox,
    Named("data_logis") = data_logis,
    Named("fit_logis") = fit_logis,
    Named("data_outcome") = as<DataFrame>(data_outcome),
    Named("fit_outcome") = fit_outcome,
    Named("fail") = fail,
    Named("psimissing") = psimissing,
    Named("settings") = settings);
  
  if (!swtrt_control_only) {
    result.push_back(psi1hat, "psi_trt");
     result.push_back(psi1hat_vec, "psi_trt_roots");
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
