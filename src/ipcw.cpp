#include "utilities.h"
#include "survival_analysis.h"
#include "logistic_regression.h"
#include "splines.h"

using namespace Rcpp;

// [[Rcpp::export]]
List ipcwcpp(
    const DataFrame data,
    const std::string id = "id",
    const StringVector& stratum = "",
    const std::string tstart = "tstart",
    const std::string tstop = "tstop",
    const std::string event = "event",
    const std::string treat = "treat",
    const std::string swtrt = "swtrt",
    const std::string swtrt_time = "swtrt_time",
    const StringVector& base_cov = "",
    const StringVector& numerator = "",
    const StringVector& denominator = "",
    const bool logistic_switching_model = false,
    const bool strata_main_effect_only = true,
    const int ns_df = 3,
    const bool firth = false,
    const bool flic = false,
    const bool stabilized_weights = true,
    const double trunc = 0,
    const bool trunc_upper_only = true,
    const bool swtrt_control_only = true,
    const double alpha = 0.05,
    const std::string ties = "efron",
    const bool boot = false,
    const int n_boot = 1000,
    const int seed = NA_INTEGER) {
  
  int k, n = data.nrow();
  
  int p = static_cast<int>(base_cov.size());
  if (p == 1 && (base_cov[0] == "" || base_cov[0] == "none")) p = 0;
  
  int p1 = static_cast<int>(numerator.size());
  if (p1 == 1 && (numerator[0] == "" || numerator[0] == "none")) p1 = 0;
  
  int p2 = static_cast<int>(denominator.size());
  if (p2 == 1 && (denominator[0] == "" || denominator[0] == "none")) {
    stop("covariates for the switch model must be provided");
  }
  
  if (p1 > 0) {
    if (p == 0 || is_true(any(is_na(match(numerator, base_cov))))) {
      stop("numerator must be a subset of base_cov");
    }
  }
  
  if (p > 0) {
    if (p2 == 0 || is_true(any(is_na(match(base_cov, denominator))))) {
      stop("base_cov must be a subset of denominator");
    }
  }
  
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

  
  // covariates for the Cox model containing treat and base_cov
  StringVector covariates(p+1);
  NumericMatrix zn(n,p);
  covariates[0] = "treated";
  for (int j=0; j<p; ++j) {
    String zj = base_cov[j];
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
  
  // covariates for the logistic regression switch model for denominator
  // including stratum, denominator, and ns_df spline terms
  StringVector covariates_lgs_den(q+p2+ns_df);
  NumericMatrix z_lgs_denn(n,q+p2);
  if (strata_main_effect_only) {
    k = 0;
    for (int i=0; i<p_stratum; ++i) {
      SEXP col_level = levels[i];
      SEXPTYPE type_level = TYPEOF(col_level);
      
      int di = d[i]-1;
      for (int j=0; j<di; ++j) {
        covariates_lgs_den[k+j] = as<std::string>(stratum[i]);
        
        if (type_level == STRSXP) {
          StringVector u = col_level;
          std::string label = sanitize(as<std::string>(u[j]));
          covariates_lgs_den[k+j] += label;
        } else if (type_level == REALSXP) {
          NumericVector u = col_level;;
          covariates_lgs_den[k+j] += std::to_string(u[j]);
        } else if (type_level == INTSXP || type_level == LGLSXP) {
          IntegerVector u = col_level;
          covariates_lgs_den[k+j] += std::to_string(u[j]);
        }
        
        z_lgs_denn(_,k+j) = (stratan(_,i) == j);
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
      covariates_lgs_den[j] = "";
      for (int i=0; i<p_stratum; ++i) {
        SEXP col_level = levels[i];
        SEXPTYPE type_level = TYPEOF(col_level);
        
        IntegerVector q_col = stratan(_,i);
        int l = q_col[first_k];
        
        covariates_lgs_den[j] += as<std::string>(stratum[i]);
        
        if (type_level == STRSXP) {
          StringVector u = col_level;
          std::string label = sanitize(as<std::string>(u[l]));
          covariates_lgs_den[j] += label;
        } else if (type_level == REALSXP) {
          NumericVector u = col_level;
          covariates_lgs_den[j] += std::to_string(u[l]);
        } else if (type_level == INTSXP || type_level == LGLSXP) {
          IntegerVector u = col_level;
          covariates_lgs_den[j] += std::to_string(u[l]);
        }
        
        if (i < p_stratum-1) {
          covariates_lgs_den[j] += ".";
        }
      }
      z_lgs_denn(_,j) = (stratumn == j);
    }
  }
  
  
  NumericMatrix z_cox_denn(n,p2);
  for (int j=0; j<p2; ++j) {
    String zj = denominator[j];
    
    if (!hasVariable(data, zj)) {
      stop("data must contain the variables in denominator");
    }
    
    if (zj == treat) {
      stop("treat should be excluded from denominator");
    }
    
    NumericVector u = data[zj];
    covariates_lgs_den[q+j] = zj;
    z_lgs_denn(_,q+j) = u;
    z_cox_denn(_,j) = u;
  }
  
  if (ns_df < 0) {
    stop("ns_df must be a nonnegative integer");
  }
  
  for (int j=0; j<ns_df; ++j) {
    covariates_lgs_den[q+p2+j] = "ns" + std::to_string(j+1);
  }
  
  // covariates for the logistic regression switch model for numerator
  // including stratum, numerator, and ns_df spline terms
  StringVector covariates_lgs_num(q+p1+ns_df);
  for (int j=0; j<q; ++j) {
    covariates_lgs_num[j] = covariates_lgs_den[j];
  }
  for (int j=0; j<p1; ++j) {
    String zj = numerator[j];
    covariates_lgs_num[q+j] = zj;
  }
  for (int j=0; j<ns_df; ++j) {
    covariates_lgs_num[q+p1+j] = "ns" + std::to_string(j+1);
  }
  
  if (trunc < 0.0 || trunc >= 0.5) {
    stop("trunc must lie in [0, 0.5)");
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
  
  // exclude observations with missing values
  LogicalVector sub(n,1);
  for (int i=0; i<n; ++i) {
    if (idn[i] == NA_INTEGER || stratumn[i] == NA_INTEGER || 
        std::isnan(tstartn[i]) || std::isnan(tstopn[i]) || 
        eventn[i] == NA_INTEGER || treatn[i] == NA_INTEGER || 
        swtrtn[i] == NA_INTEGER) {
      sub[i] = 0;
    }
    for (int j=0; j<q+p2; ++j) {
      if (std::isnan(z_lgs_denn(i,j))) sub[i] = 0;
    }
  }
  
  IntegerVector order = which(sub);
  idn = idn[order];
  stratumn = stratumn[order];
  tstartn = tstartn[order];
  tstopn = tstopn[order];
  eventn = eventn[order];
  treatn = treatn[order];
  swtrtn = swtrtn[order];
  swtrt_timen = swtrt_timen[order];
  if (p > 0) zn = subset_matrix_by_row(zn, order);
  z_cox_denn = subset_matrix_by_row(z_cox_denn, order);
  z_lgs_denn = subset_matrix_by_row(z_lgs_denn, order);
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
    NumericMatrix zn1(n+k, zn.ncol());
    NumericMatrix z_cox_denn1(n+k, z_cox_denn.ncol());
    NumericMatrix z_lgs_denn1(n+k, z_lgs_denn.ncol());
    for (int i=0; i<n; ++i) {
      zn1(i,_) = zn(i,_);
      z_cox_denn1(i,_) = z_cox_denn(i,_);
      z_lgs_denn1(i,_) = z_lgs_denn(i,_);
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
      swtrtn.push_back(swtrtn[l]);
      swtrt_timen.push_back(swtrt_timen[l]);
      zn1(n+i,_) = zn(l,_);
      z_cox_denn1(n+i,_) = z_cox_denn(l,_);
      z_lgs_denn1(n+i,_) = z_lgs_denn(l,_);
      
      // change tstop and event for the old observation
      tstopn[l] = swtrt_timen[l];
      eventn[l] = 0;
    }
    
    // update number of rows and old matrices
    n = n + k;
    zn = zn1;
    z_cox_denn = z_cox_denn1;
    z_lgs_denn = z_lgs_denn1;
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
  swtrtn = swtrtn[order];
  swtrt_timen = swtrt_timen[order];
  zn = subset_matrix_by_row(zn, order);
  z_cox_denn = subset_matrix_by_row(z_cox_denn, order);
  z_lgs_denn = subset_matrix_by_row(z_lgs_denn, order);
  
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
    int k = idx1[i];
    for (int j=idx[i]; j<idx[i+1]; ++j) {
      osn[j] = eventn[k];
      os_timen[j] = tstopn[k];
    }
  }
  
  if (is_true(any(ifelse(swtrtn == 1, swtrt_timen > os_timen, 0)))) {
    stop("swtrt_time must be less than or equal to os_time");
  }
  
  IntegerVector stratumn1 = stratumn[idx1];
  IntegerVector treatn1 = treatn[idx1];
  NumericVector tstopn1 = tstopn[idx1];
  IntegerVector eventn1 = eventn[idx1];
  
  DataFrame dt = DataFrame::create(
    Named("stratum") = stratumn1,
    Named("treat") = treatn1,
    Named("time") = tstopn1,
    Named("event") = eventn1);
  
  DataFrame lr = lrtest(dt,"","stratum","treat","time","","event","",0,0,0);
  double logRankPValue = lr["logRankPValue"];
  
  double zcrit = R::qnorm(1-alpha/2, 0, 1, 1, 0);
  
  k = -1; // indicate the observed data
  auto f = [&k, has_stratum, stratum, p_stratum, u_stratum, 
            type_treat, treat, treatwi, treatwn, treatwc,
            type_id, id, idwi, idwn, idwc,
            q, p, p1, p2, covariates, numerator, denominator, 
            covariates_lgs_num, covariates_lgs_den, 
            logistic_switching_model, ns_df, firth, flic, 
            stabilized_weights, trunc, trunc_upper_only, 
            swtrt_control_only, alpha, zcrit, ties](
                IntegerVector& idb, IntegerVector& stratumb,
                NumericVector& tstartb, NumericVector& tstopb,
                IntegerVector& eventb, IntegerVector& treatb,
                NumericVector& os_timeb, 
                IntegerVector& swtrtb, NumericVector& swtrt_timeb,
                NumericMatrix& zb, NumericMatrix& z_cox_denb,
                NumericMatrix& z_lgs_denb)->List {
                  bool fail = false; // whether any model fails to converge
                  NumericVector init(1, NA_REAL);
                  
                  // exclude observations after treatment switch
                  LogicalVector c1 = ifelse(
                    swtrtb == 1, tstartb < swtrt_timeb, 1);
                  IntegerVector l = which(c1);
                  IntegerVector id1 = idb[l];
                  IntegerVector stratum1 = stratumb[l];
                  NumericVector tstart1 = tstartb[l];
                  NumericVector tstop1 = tstopb[l];
                  IntegerVector event1 = eventb[l];
                  IntegerVector treat1 = treatb[l];
                  NumericVector os_time1 = os_timeb[l];
                  IntegerVector swtrt1 = swtrtb[l];
                  NumericVector swtrt_time1 = swtrt_timeb[l];
                  NumericMatrix z1 = subset_matrix_by_row(zb, l);
                  NumericMatrix z_cox_den1 = 
                    subset_matrix_by_row(z_cox_denb, l);
                  NumericMatrix z_lgs_den1 = 
                    subset_matrix_by_row(z_lgs_denb, l);
                  int n1 = static_cast<int>(l.size());
                  
                  // set up crossover indicators
                  IntegerVector cross1(n1);
                  for (int i=0; i<n1; ++i) {
                    if (i == n1-1 || id1[i] != id1[i+1]) {
                      if (swtrt1[i] == 1 && tstop1[i] >= swtrt_time1[i]) {
                        cross1[i] = 1;
                      }
                    }
                  }
                  
                  // initialize data_switch and fit_switch
                  List data_switch(2), fit_switch(2);
                  if (k == -1) {
                    for (int h=0; h<2; ++h) {
                      List data_x = List::create(
                        Named("data") = R_NilValue,
                        Named(treat) = R_NilValue
                      );
                      
                      List fit_x = List::create(
                        Named("fit_den") = R_NilValue,
                        Named("fit_num") = R_NilValue,        
                        Named(treat) = R_NilValue
                      );
                      
                      if (type_treat == LGLSXP || type_treat == INTSXP) {
                        data_x[treat] = treatwi[1-h];
                        fit_x[treat] = treatwi[1-h];
                      } else if (type_treat == REALSXP) {
                        data_x[treat] = treatwn[1-h];
                        fit_x[treat] = treatwn[1-h];
                      } else if (type_treat == STRSXP) {
                        data_x[treat] = treatwc[1-h];
                        fit_x[treat] = treatwc[1-h];
                      }
                      
                      data_switch[h] = data_x;
                      fit_switch[h] = fit_x;
                    }
                  }
                  
                  List data_outcome;
                  
                  // treat arms with patients who switched treatment
                  IntegerVector treats(1);
                  treats[0] = 0;
                  if (!swtrt_control_only) {
                    treats.push_back(1);
                  }
                  
                  int K = static_cast<int>(treats.size());
                  if (!logistic_switching_model) { // time-dependent Cox
                    // obtain the unique event times across treatment groups
                    NumericVector cut = tstop1[event1 == 1];
                    cut = unique(cut);
                    cut.sort();
                    
                    // replicate event times within each subject
                    IntegerVector id2, stratum2, event2, treat2, cross2;
                    NumericVector tstart2, tstop2;
                    NumericMatrix z2, z_cox_den2;
                    int n2;
                    if (!swtrt_control_only) {
                      DataFrame a = survsplit(tstart1, tstop1, cut);
                      IntegerVector censor = a["censor"];
                      IntegerVector l = a["row"];
                      id2 = id1[l];
                      stratum2 = stratum1[l];
                      tstart2 = a["start"];
                      tstop2 = a["end"];
                      event2 = event1[l];
                      treat2 = treat1[l];
                      cross2 = cross1[l];
                      z2 = subset_matrix_by_row(z1, l);
                      z_cox_den2 = subset_matrix_by_row(z_cox_den1, l);
                      n2 = static_cast<int>(l.size());
                      for (int i=0; i<n2; ++i) {
                        if (censor[i] == 1) {
                          event2[i] = 0;
                          cross2[i] = 0;
                        }
                      }
                    } else {
                      // extract data for the control group
                      IntegerVector l = which(treat1 == 0);
                      IntegerVector id10 = id1[l];
                      IntegerVector stratum10 = stratum1[l];
                      NumericVector tstart10 = tstart1[l];
                      NumericVector tstop10 = tstop1[l];
                      IntegerVector event10 = event1[l];
                      IntegerVector treat10 = treat1[l];
                      IntegerVector cross10 = cross1[l];
                      NumericMatrix z10 = subset_matrix_by_row(z1, l);
                      NumericMatrix z_cox_den10 = 
                        subset_matrix_by_row(z_cox_den1, l);
                      
                      // replicate event times within each subject
                      DataFrame a = survsplit(tstart10, tstop10, cut);
                      IntegerVector censor = a["censor"];
                      l = a["row"];
                      IntegerVector id20 = id10[l];
                      IntegerVector stratum20 = stratum10[l];
                      NumericVector tstart20 = a["start"];
                      NumericVector tstop20 = a["end"];
                      IntegerVector event20 = event10[l];
                      IntegerVector treat20 = treat10[l];
                      IntegerVector cross20 = cross10[l];
                      NumericMatrix z20 = subset_matrix_by_row(z10, l);
                      NumericMatrix z_cox_den20 = 
                        subset_matrix_by_row(z_cox_den10, l);
                      int n20 = static_cast<int>(l.size());
                      for (int i=0; i<n20; ++i) {
                        if (censor[i] == 1) {
                          event20[i] = 0;
                          cross20[i] = 0;
                        }
                      }
                      
                      // extract data for the active group
                      l = which(treat1 == 1);
                      IntegerVector id21 = id1[l];
                      IntegerVector stratum21 = stratum1[l];
                      NumericVector tstart21 = tstart1[l];
                      NumericVector tstop21 = tstop1[l];
                      IntegerVector event21 = event1[l];
                      IntegerVector treat21 = treat1[l];
                      IntegerVector cross21 = cross1[l];
                      NumericMatrix z21 = subset_matrix_by_row(z1, l);
                      NumericMatrix z_cox_den21 = 
                        subset_matrix_by_row(z_cox_den1, l);
                      int n21 = static_cast<int>(l.size());
                      
                      // combine weighted control with unweighted active data
                      id2 = c_vectors_i(id20, id21);
                      stratum2 = c_vectors_i(stratum20, stratum21);
                      tstart2 = c_vectors(tstart20, tstart21);
                      tstop2 = c_vectors(tstop20, tstop21);
                      event2 = c_vectors_i(event20, event21);
                      treat2 = c_vectors_i(treat20, treat21);
                      cross2 = c_vectors_i(cross20, cross21);
                      z2 = c_matrices(z20, z21);
                      z_cox_den2 = c_matrices(z_cox_den20, z_cox_den21);
                      n2 = n20 + n21;
                    }
                    
                    // initialize weights
                    NumericVector w2(n2, 1.0), sw2(n2, 1.0);
                    
                    // fit the switching models by treatment group
                    for (int h=0; h<K; ++h) {
                      IntegerVector l = which(treat2 == h);
                      IntegerVector id3 = id2[l];
                      IntegerVector stratum3 = stratum2[l];
                      NumericVector tstart3 = tstart2[l];
                      NumericVector tstop3 = tstop2[l];
                      IntegerVector cross3 = cross2[l];
                      NumericMatrix z_cox_den3 = 
                        subset_matrix_by_row(z_cox_den2, l);
                      int n3 = static_cast<int>(l.size());
                      
                      // prepare the data for fitting the switching model
                      List data1 = List::create(
                        Named("uid") = id3,
                        Named("ustratum") = stratum3,
                        Named("tstart") = tstart3,
                        Named("tstop") = tstop3,
                        Named("cross") = cross3);
                      
                      for (int j=0; j<p2; ++j) {
                        String zj = denominator[j];
                        NumericVector u = z_cox_den3(_,j);
                        data1.push_back(u,zj);
                      }
                      
                      // fit the denominator model for crossover
                      List fit_den = phregcpp(
                        data1, "", "ustratum", "tstart", "tstop", "cross",
                        denominator, "", "", "uid", ties, init, 
                        0, 1, 0, 0, 0, alpha, 50, 1.0e-9);
                      
                      // obtain the survival probabilities for crossover
                      DataFrame sumstat_den = DataFrame(fit_den["sumstat"]);
                      bool fail_den = sumstat_den["fail"];
                      if (fail_den) fail = true;
                      
                      DataFrame parest_den = DataFrame(fit_den["parest"]);
                      NumericVector beta_den = parest_den["beta"];
                      NumericMatrix vbeta_den(p2, p2);
                      DataFrame basehaz_den = DataFrame(fit_den["basehaz"]);
                      
                      DataFrame km_den = survfit_phregcpp(
                        p2, beta_den, vbeta_den, basehaz_den, data1,
                        denominator, "ustratum", "", "uid", "tstart", 
                        "tstop", 0, "log-log", 1-alpha);
                      
                      NumericVector surv_den = km_den["surv"];
                      
                      List fit_num = phregcpp(
                        data1, "", "ustratum", "tstart", "tstop", "cross", 
                        numerator, "", "", "uid", ties, init, 
                        0, 1, 0, 0, 0, alpha, 50, 1.0e-9);
                      
                      int p10 = p1 > 0 ? p1 : 1;
                      NumericVector beta_num(p10);
                      NumericMatrix vbeta_num(p10,p10);
                      if (p1 > 0) {
                        DataFrame sumstat_num= DataFrame(fit_num["sumstat"]);
                        bool fail_num = sumstat_num["fail"];
                        if (fail_num) fail = true;
                        
                        DataFrame parest_num = DataFrame(fit_num["parest"]);
                        beta_num = parest_num["beta"];
                      }
                      DataFrame basehaz_num = DataFrame(fit_num["basehaz"]);
                      
                      DataFrame km_num = survfit_phregcpp(
                        p1, beta_num, vbeta_num, basehaz_num, data1,
                        numerator, "ustratum", "", "uid", "tstart", 
                        "tstop", 0, "log-log", 1-alpha);
                      
                      NumericVector surv_num = km_num["surv"];
                      
                      // update data_switch and fit_switch
                      if (k == -1) {
                        IntegerVector uid = data1["uid"];
                        if (type_id == INTSXP) {
                          data1.push_front(idwi[uid], id);
                        } else if (type_id == REALSXP) {
                          data1.push_front(idwn[uid], id);
                        } else if (type_id == STRSXP) {
                          data1.push_front(idwc[uid], id);
                        }
                        
                        
                        if (has_stratum) {
                          IntegerVector ustratum = data1["ustratum"];
                          for (int i=0; i<p_stratum; ++i) {
                            std::string s = as<std::string>(stratum[i]);
                            SEXP col_stratum = u_stratum[s];
                            SEXPTYPE type_stratum = TYPEOF(col_stratum);
                            if (type_stratum == INTSXP) {
                              IntegerVector v = col_stratum;
                              data1.push_back(v[ustratum], s);
                            } else if (type_stratum == REALSXP) {
                              NumericVector v = col_stratum;
                              data1.push_back(v[ustratum], s);
                            } else if (type_stratum == STRSXP) {
                              StringVector v = col_stratum;
                              data1.push_back(v[ustratum], s);
                            }
                          }
                        }
                        
                        List data_x = data_switch[h];
                        data_x["data"] = as<DataFrame>(data1);
                        data_switch[h] = data_x;
                        
                        List fit_x = fit_switch[h];
                        fit_x["fit_den"] = fit_den;
                        fit_x["fit_num"]= fit_num;
                        fit_switch[h] = fit_x;
                      }
                      
                      // unstabilized and stabilized weights
                      NumericVector w = 1.0/surv_den;
                      NumericVector sw = surv_num/surv_den;
                      
                      // weights in the outcome model by matching id and time
                      IntegerVector id4 = km_den["uid"];
                      NumericVector time4 = km_den["time"];
                      IntegerVector sub = match3(id3, tstop3, id4, time4);
                      NumericVector w3 = w[sub];
                      NumericVector sw3 = sw[sub];
                      
                      // fill in missing weights with LOCF starting at 1
                      IntegerVector idx3(1,0);
                      for (int i=1; i<n3; ++i) {
                        if (id3[i] != id3[i-1]) {
                          idx3.push_back(i);
                        }
                      }
                      
                      int nids3 = static_cast<int>(idx3.size());
                      idx3.push_back(n3);
                      
                      for (int i=0; i<nids3; ++i) {
                        if (std::isnan(w3[idx3[i]])) {
                          w3[idx3[i]] = 1.0;
                          sw3[idx3[i]] = 1.0;
                        }
                        for (int j=idx3[i]+1; j<idx3[i+1]; ++j) {
                          if (std::isnan(w3[j])) {
                            w3[j] = w3[j-1];
                            sw3[j] = sw3[j-1];
                          }
                        }
                      }
                      
                      // truncate the weights if requested
                      if (trunc > 0.0) {
                        int m = w3.size();
                        
                        // truncated unstabilized weights
                        if (trunc_upper_only) {
                          double upper = quantilecpp(w3, 1-trunc);
                          for (int i=0; i<m; ++i) {
                            if (w3[i] > upper) w3[i] = upper;
                          }
                        } else {
                          double lower = quantilecpp(w3, trunc);
                          double upper = quantilecpp(w3, 1-trunc);
                          for (int i=0; i<m; ++i) {
                            if (w3[i] < lower) {
                              w3[i] = lower;
                            } else if (w3[i] > upper) {
                              w3[i] = upper;
                            }
                          }
                        }
                        
                        // truncated stabilized weights
                        if (trunc_upper_only) {
                          double upper = quantilecpp(sw3, 1-trunc);
                          for (int i=0; i<m; ++i) {
                            if (sw3[i] > upper) sw3[i] = upper;
                          }
                        } else {
                          double lower = quantilecpp(sw3, trunc);
                          double upper = quantilecpp(sw3, 1-trunc);
                          for (int i=0; i<m; ++i) {
                            if (sw3[i] < lower) {
                              sw3[i] = lower;
                            } else if (sw3[i] > upper) {
                              sw3[i] = upper;
                            }
                          }
                        }
                      }
                      
                      // fill in the weights
                      w2[l] = w3;
                      sw2[l] = sw3;
                    }
                    
                    // prepare data for the outcome model
                    data_outcome = List::create(
                      Named("uid") = id2,
                      Named("tstart") = tstart2,
                      Named("tstop") = tstop2,
                      Named("event") = event2,
                      Named("treated") = treat2,
                      Named("unstabilized_weight") = w2,
                      Named("stabilized_weight") = sw2);
                    
                    data_outcome.push_back(stratum2, "ustratum");
                    
                    for (int j=0; j<p; ++j) {
                      NumericVector u = z2(_,j);
                      String zj = covariates[j+1];
                      data_outcome.push_back(u,zj);
                    }
                  } else { // logistic regression switching model
                    NumericVector w1(n1, 1.0), sw1(n1, 1.0);
                    
                    // fit the switching models by treatment group
                    for (int h=0; h<K; ++h) {
                      LogicalVector c2 = ifelse(
                        swtrt1 == 1, tstart1 < swtrt_time1, tstop1<os_time1);
                      IntegerVector l = which(c2 & (treat1 == h));
                      IntegerVector id2 = id1[l];
                      IntegerVector stratum2 = stratum1[l];
                      NumericVector tstart2 = tstart1[l];
                      NumericVector tstop2 = tstop1[l];
                      IntegerVector cross2 = cross1[l];
                      NumericMatrix z_lgs_den2 = 
                        subset_matrix_by_row(z_lgs_den1, l);
                      int n2 = static_cast<int>(l.size());
                      
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
                      
                      // prepare the data for fitting the switching model
                      List data1 = List::create(
                        Named("uid") = id2,
                        Named("ustratum") = stratum2,
                        Named("tstart") = tstart2,
                        Named("tstop") = tstop2,
                        Named("cross") = cross2);
                      
                      for (int j=0; j<q+p2; ++j) {
                        String zj = covariates_lgs_den[j];
                        NumericVector u = z_lgs_den2(_,j);
                        data1.push_back(u,zj);
                      }
                      for (int j=0; j<ns_df; ++j) {
                        String zj = covariates_lgs_den[q+p2+j];
                        NumericVector u = s(_,j);
                        data1.push_back(u,zj);
                      }
                      
                      List fit_den = logisregcpp(
                        data1, "", "cross", covariates_lgs_den, "", "", 
                        "", "uid", "logit", init, 
                        0, firth, flic, 0, alpha, 50, 1.0e-9);
                      
                      DataFrame sumstat_den = DataFrame(fit_den["sumstat"]);
                      bool fail_den = sumstat_den["fail"];
                      if (fail_den) fail = true;
                      
                      DataFrame f_den = DataFrame(fit_den["fitted"]);
                      NumericVector h_den = f_den["fitted_values"];
                      
                      List fit_num = logisregcpp(
                        data1, "", "cross", covariates_lgs_num, "", "", 
                        "", "uid", "logit", init, 
                        0, firth, flic, 0, alpha, 50, 1.0e-9);
                      
                      DataFrame sumstat_num = DataFrame(fit_num["sumstat"]);
                      bool fail_num = sumstat_num["fail"];
                      if (fail_num) fail = true;
                      
                      DataFrame f_num = DataFrame(fit_num["fitted"]);
                      NumericVector h_num = f_num["fitted_values"];
                      
                      // update data_switch and fit_switch
                      if (k == -1) {
                        IntegerVector uid = data1["uid"];
                        if (type_id == INTSXP) {
                          data1.push_front(idwi[uid], id);
                        } else if (type_id == REALSXP) {
                          data1.push_front(idwn[uid], id);
                        } else if (type_id == STRSXP) {
                          data1.push_front(idwc[uid], id);
                        }
                        
                        if (has_stratum) {
                          IntegerVector ustratum = data1["ustratum"];
                          for (int i=0; i<p_stratum; ++i) {
                            std::string s = as<std::string>(stratum[i]);
                            SEXP col_stratum = u_stratum[s];
                            SEXPTYPE type_stratum = TYPEOF(col_stratum);
                            if (type_stratum == INTSXP) {
                              IntegerVector v = col_stratum;
                              data1.push_back(v[ustratum], s);
                            } else if (type_stratum == REALSXP) {
                              NumericVector v = col_stratum;
                              data1.push_back(v[ustratum], s);
                            } else if (type_stratum == STRSXP) {
                              StringVector v = col_stratum;
                              data1.push_back(v[ustratum], s);
                            }
                          }
                        }
                        
                        List data_x = data_switch[h];
                        data_x["data"] = as<DataFrame>(data1);
                        data_switch[h] = data_x;
                        
                        List fit_x = fit_switch[h];
                        fit_x["fit_den"] = fit_den;
                        fit_x["fit_num"] = fit_num;
                        fit_switch[h] = fit_x;
                      }
                      
                      // convert to probability of observed response 
                      NumericVector o_den(n2), o_num(n2);
                      for (int i=0; i<n2; ++i) {
                        o_den[i] = cross2[i] == 1 ? h_den[i] : 1 - h_den[i];
                        o_num[i] = cross2[i] == 1 ? h_num[i] : 1 - h_num[i];
                      }
                      
                      // obtain cumulative products within a subject
                      l = which(treat1 == h);
                      IntegerVector id3 = id1[l];
                      IntegerVector swtrt3 = swtrt1[l];
                      int n3 = static_cast<int>(l.size());
                      
                      IntegerVector idx3(1,0);
                      for (int i=1; i<n3; ++i) {
                        if (id3[i] != id3[i-1]) {
                          idx3.push_back(i);
                        }
                      }
                      
                      IntegerVector swtrt3u = swtrt3[idx3];
                      
                      int nids3 = static_cast<int>(idx3.size());
                      idx3.push_back(n3);
                      
                      NumericVector p_den(n3, 1.0), p_num(n3, 1.0);
                      
                      int m = 0; // index for id2
                      for (int i=0; i<nids3; ++i) {
                        int r = m - idx3[i] - 1;
                        for (int j=idx3[i]+1; j<idx3[i+1]; ++j) {
                          p_den[j] = p_den[j-1]*o_den[r+j];
                          p_num[j] = p_num[j-1]*o_num[r+j];
                        }
                        
                        int ni = idx3[i+1] - idx3[i];
                        int mi = swtrt3u[i] == 1 ? ni : ni - 1;
                        m += mi;
                      }
                      
                      // unstabilized and stabilized weights
                      NumericVector w3 = 1.0/p_den;
                      NumericVector sw3 = p_num/p_den;
                      
                      // truncate the weights if requested
                      if (trunc > 0.0) {
                        // truncated unstabilized weights
                        if (trunc_upper_only) {
                          double upper = quantilecpp(w3, 1-trunc);
                          for (int i=0; i<n3; ++i) {
                            if (w3[i] > upper) w3[i] = upper;
                          }
                        } else {
                          double lower = quantilecpp(w3, trunc);
                          double upper = quantilecpp(w3, 1-trunc);
                          for (int i=0; i<n3; ++i) {
                            if (w3[i] < lower) {
                              w3[i] = lower;
                            } else if (w3[i] > upper) {
                              w3[i] = upper;
                            }
                          }
                        }
                        
                        // truncated stabilized weights
                        if (trunc_upper_only) {
                          double upper = quantilecpp(sw3, 1-trunc);
                          for (int i=0; i<n3; ++i) {
                            if (sw3[i] > upper) sw3[i] = upper;
                          }
                        } else {
                          double lower = quantilecpp(sw3, trunc);
                          double upper = quantilecpp(sw3, 1-trunc);
                          for (int i=0; i<n3; ++i) {
                            if (sw3[i] < lower) {
                              sw3[i] = lower;
                            } else if (sw3[i] > upper) {
                              sw3[i] = upper;
                            }
                          }
                        }
                      }
                      
                      // fill in the weights
                      w1[l] = w3;
                      sw1[l] = sw3;
                    }
                    
                    // prepare data for the outcome model
                    data_outcome = List::create(
                      Named("uid") = id1,
                      Named("tstart") = tstart1,
                      Named("tstop") = tstop1,
                      Named("event") = event1,
                      Named("treated") = treat1,
                      Named("unstabilized_weight") = w1,
                      Named("stabilized_weight") = sw1);
                    
                    data_outcome.push_back(stratum1, "ustratum");
                    
                    for (int j=0; j<p; ++j) {
                      String zj = covariates[j+1];
                      NumericVector u = z1(_,j);
                      data_outcome.push_back(u,zj);
                    }
                  }
                  

                  std::string weight_variable = stabilized_weights ? 
                  "stabilized_weight" : "unstabilized_weight";
                  
                  // generate weighted KM estimate and log-rank test
                  List km_outcome, lr_outcome;
                  if (k == -1) {
                    km_outcome = kmest(
                      data_outcome, "", "treated", "tstart", "tstop", 
                      "event", weight_variable, "log-log", 1-alpha, 1);
                    
                    lr_outcome = lrtest(
                      data_outcome, "", "ustratum", "treated", "tstart", 
                      "tstop", "event", weight_variable, 0, 0, 0);
                  }
                  
                  // fit the outcome model with weights
                  List fit_outcome = phregcpp(
                    data_outcome, "", "ustratum", "tstart", "tstop",
                    "event", covariates, weight_variable, "",
                    "uid", ties, init, 1, 0, 0, 0, 0, alpha, 50, 1.0e-9);
                  
                  DataFrame sumstat = DataFrame(fit_outcome["sumstat"]);
                  bool fail_cox = sumstat["fail"];
                  if (fail_cox) fail = true;
                  
                  DataFrame parest = DataFrame(fit_outcome["parest"]);
                  NumericVector beta = parest["beta"];
                  NumericVector sebeta = parest["sebeta"];
                  NumericVector pval = parest["p"];
                  double hrhat = exp(beta[0]);
                  double hrlower = exp(beta[0] - zcrit*sebeta[0]);
                  double hrupper = exp(beta[0] + zcrit*sebeta[0]);
                  double pvalue = pval[0];
                  
                  List out;
                  if (k == -1) {
                    out = List::create(
                      Named("data_switch") = data_switch,
                      Named("fit_switch") = fit_switch,
                      Named("data_outcome") = data_outcome,
                      Named("km_outcome") = km_outcome,
                      Named("lr_outcome") = lr_outcome,
                      Named("fit_outcome") = fit_outcome,
                      Named("hrhat") = hrhat,
                      Named("hrlower") = hrlower,
                      Named("hrupper") = hrupper,
                      Named("pvalue") = pvalue,
                      Named("fail") = fail);
                  } else {
                    out = List::create(
                      Named("hrhat") = hrhat,
                      Named("hrlower") = hrlower,
                      Named("hrupper") = hrupper,
                      Named("pvalue") = pvalue,
                      Named("fail") = fail);
                  }
                  
                  return out;
                };
  
  List out = f(idn, stratumn, tstartn, tstopn, eventn, treatn, os_timen,
               swtrtn, swtrt_timen, zn, z_cox_denn, z_lgs_denn);
  
  List data_switch = out["data_switch"];
  List fit_switch = out["fit_switch"];
  List data_outcome = out["data_outcome"];
  List km_outcome = out["km_outcome"];
  List lr_outcome = out["lr_outcome"];
  List fit_outcome = out["fit_outcome"];
  
  double hrhat = out["hrhat"];
  double hrlower = out["hrlower"];
  double hrupper = out["hrupper"];
  double pvalue = out["pvalue"];
  bool fail = out["fail"];
  
  
  NumericVector hrhats(n_boot);
  LogicalVector fails(n_boot);
  List fail_boots_data;
  String hr_CI_type;
  
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
        IntegerVector stratumwi = col_stratum;
        data_outcome.push_back(stratumwi[ustratum], s);
      } else if (type_stratum == REALSXP) {
        NumericVector stratumwn = col_stratum;
        data_outcome.push_back(stratumwn[ustratum], s);
      } else if (type_stratum == STRSXP) {
        StringVector stratumwc = col_stratum;
        data_outcome.push_back(stratumwc[ustratum], s);
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
    IntegerVector osc(B), swtrtc(B);
    NumericVector tstartc(B), tstopc(B), os_timec(B), swtrt_timec(B);
    NumericMatrix z_lgs_denc(B,q+p2);
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
      swtrtn = swtrtn[order];
      swtrt_timen = swtrt_timen[order];
      zn = subset_matrix_by_row(zn, order);
      z_cox_denn = subset_matrix_by_row(z_cox_denn, order);
      z_lgs_denn = subset_matrix_by_row(z_lgs_denn, order);
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
      IntegerVector osb(N), swtrtb(N);
      NumericVector tstartb(N), tstopb(N), os_timeb(N), swtrt_timeb(N);
      NumericMatrix zb(N,p), z_cox_denb(N,p2), z_lgs_denb(N,q+p2);
      
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
            swtrtb[l] = swtrtn[j];
            swtrt_timeb[l] = swtrt_timen[j];
            zb(l,_) = zn(j,_);
            z_cox_denb(l,_) = z_cox_denn(j,_);
            z_lgs_denb(l,_) = z_lgs_denn(j,_);
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
      swtrtb = swtrtb[sub];
      swtrt_timeb = swtrt_timeb[sub];
      zb = subset_matrix_by_row(zb, sub);
      z_cox_denb = subset_matrix_by_row(z_cox_denb, sub);
      z_lgs_denb = subset_matrix_by_row(z_lgs_denb, sub);
      
      out = f(idb, stratumb, tstartb, tstopb, eventb, treatb, os_timeb,
              swtrtb, swtrt_timeb, zb, z_cox_denb, z_lgs_denb);
      
      fails[k] = out["fail"];
      hrhats[k] = out["hrhat"];
      
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
          swtrtc[j] = swtrtb[i];
          swtrt_timec[j] = swtrt_timeb[i];
          z_lgs_denc(j,_) = z_lgs_denb(i,_);
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
      swtrtc = swtrtc[sub];
      swtrt_timec = swtrt_timec[sub];
      z_lgs_denc = subset_matrix_by_row(z_lgs_denc,sub);
      
      fail_boots_data = List::create(
        Named("boot_index") = boot_indexc,
        Named("uid") = idc,
        Named("tstart") = tstartc,
        Named("tstop") = tstopc,
        Named("event") = eventc,
        Named("treated") = treatc,
        Named("os") = osc,
        Named("os_time") = os_timec,
        Named("swtrt") = swtrtc,
        Named("swtrt_time") = swtrt_timec
      );
      
      for (int j=0; j<q+p2; ++j) {
        String zj = covariates_lgs_den[j];
        NumericVector u = z_lgs_denc(_,j);
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
            IntegerVector stratumwi = col_stratum;
            fail_boots_data.push_back(stratumwi[stratumc], s);
          } else if (type_stratum == REALSXP) {
            NumericVector stratumwn = col_stratum;
            fail_boots_data.push_back(stratumwn[stratumc], s);
          } else if (type_stratum == STRSXP) {
            StringVector stratumwc = col_stratum;
            fail_boots_data.push_back(stratumwc[stratumc], s);
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
  }
  
  List settings = List::create(
    Named("logistic_switching_model") = logistic_switching_model,
    Named("strata_main_effect_only") = strata_main_effect_only,
    Named("ns_df") = ns_df,
    Named("firth") = firth,
    Named("flic") = flic,
    Named("stabilized_weights") = stabilized_weights,
    Named("trunc") = trunc,
    Named("trunc_upper_only") = trunc_upper_only,
    Named("swtrt_control_only") = swtrt_control_only,
    Named("alpha") = alpha,
    Named("ties") = ties,
    Named("boot") = boot,
    Named("n_boot") = n_boot,
    Named("seed") = seed);
  
  List result = List::create(
    Named("logrank_pvalue") = logRankPValue,
    Named("cox_pvalue") = pvalue,
    Named("hr") = hrhat,
    Named("hr_CI") = NumericVector::create(hrlower, hrupper),
    Named("hr_CI_type") = hr_CI_type,
    Named("data_switch") = data_switch,
    Named("fit_switch") = fit_switch,
    Named("data_outcome") = as<DataFrame>(data_outcome),
    Named("km_outcome") = as<DataFrame>(km_outcome),
    Named("lr_outcome") = as<DataFrame>(lr_outcome),
    Named("fit_outcome") = fit_outcome,
    Named("fail") = fail,
    Named("settings") = settings);
  
  if (boot) {
    result.push_back(fails, "fail_boots");
    result.push_back(hrhats, "hr_boots"); 
    if (is_true(any(fails))) {
      result.push_back(as<DataFrame>(fail_boots_data), "fail_boots_data");
    }
  }
  
  return result;
}
