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
    const std::string swtrt_time_lower = "",
    const std::string swtrt_time_upper = "",
    const StringVector& base_cov = "",
    const StringVector& numerator = "",
    const StringVector& denominator = "",
    const bool logistic_switching_model = 0,
    const bool strata_main_effect_only = 1,
    const bool firth = 0,
    const bool flic = 0,
    const int ns_df = 3,
    const bool relative_time = 0,
    const bool stabilized_weights = 1,
    const double trunc = 0,
    const bool trunc_upper_only = 1,
    const bool swtrt_control_only = 1,
    const double alpha = 0.05,
    const std::string ties = "efron",
    const bool boot = 1,
    const int n_boot = 1000,
    const int seed = NA_INTEGER) {
  
  int i, j, k, n = data.nrow();
  
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
    
    if (is_true(any(is_na(match(numerator, denominator))))) {
      stop("numerator must be a subset of denominator");
    }
  }
  
  bool has_id = hasVariable(data, id);
  bool has_tstart = hasVariable(data, tstart);
  bool has_tstop = hasVariable(data, tstop);
  bool has_event = hasVariable(data, event);
  bool has_treat = hasVariable(data, treat);
  bool has_swtrt = hasVariable(data, swtrt);
  bool has_swtrt_time = hasVariable(data, swtrt_time);
  bool has_swtrt_time_lower = hasVariable(data, swtrt_time_lower);
  bool has_swtrt_time_upper = hasVariable(data, swtrt_time_upper);
  
  // create the numeric id variable
  if (!has_id) stop("data must contain the id variable");
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
  
  if (!has_tstart) {
    stop("data must contain the tstart variable"); 
  }
  
  if (TYPEOF(data[tstart]) != INTSXP && TYPEOF(data[tstart]) != REALSXP) {
    stop("tstart must take numeric values");
  }
  
  NumericVector tstartnz = data[tstart];
  NumericVector tstartn = clone(tstartnz);
  if (is_true(any(tstartn < 0.0))) {
    stop("tstart must be nonnegative for each observation");
  }
  
  if (!has_tstop) {
    stop("data must contain the tstop variable"); 
  }
  
  if (TYPEOF(data[tstop]) != INTSXP && TYPEOF(data[tstop]) != REALSXP) {
    stop("tstop must take numeric values");
  }
  
  NumericVector tstopnz = data[tstop];
  NumericVector tstopn = clone(tstopnz);
  if (is_true(any(tstopn <= tstartn))) {
    stop("tstop must be greater than tstart for each observation");
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
    stop("event must be 1 or 0 for each subject");
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
  
  if (!has_swtrt) {
    stop("data must contain the swtrt variable"); 
  }
  
  if (TYPEOF(data[swtrt]) != INTSXP && TYPEOF(data[swtrt]) != LGLSXP) {
    stop("swtrt must take integer or logical values");
  }
  
  IntegerVector swtrtnz = data[swtrt];
  IntegerVector swtrtn = clone(swtrtnz);
  if (is_true(any((swtrtn != 1) & (swtrtn != 0)))) {
    stop("swtrt must be 1 or 0 for each subject");
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
  }
  
  if (!has_swtrt_time_lower) {
    stop("data must contain the swtrt_time_lower variable"); 
  }
  
  if (TYPEOF(data[swtrt_time_lower]) != INTSXP &&
      TYPEOF(data[swtrt_time_lower]) != REALSXP) {
    stop("swtrt_time_lower must take numeric values");
  }
  
  NumericVector swtrt_time_lowernz = data[swtrt_time_lower];
  NumericVector swtrt_time_lowern = clone(swtrt_time_lowernz);
  if (is_true(any(swtrt_time_lowern < 0.0))) {
    stop("swtrt_time_lower must be nonnegative");
  }
  
  if (!has_swtrt_time_upper) {
    stop("data must contain the swtrt_time_upper variable"); 
  }
  
  if (TYPEOF(data[swtrt_time_upper]) != INTSXP &&
      TYPEOF(data[swtrt_time_upper]) != REALSXP) {
    stop("swtrt_time_upper must take numeric values");
  }
  
  NumericVector swtrt_time_uppernz = data[swtrt_time_upper];
  NumericVector swtrt_time_uppern = clone(swtrt_time_uppernz);
  if (is_true(any(swtrt_time_uppern < 0.0))) {
    stop("swtrt_time_upper must be nonnegative");
  }
  
  for (i=0; i<n; i++) {
    if (swtrtn[i] == 1 && swtrt_timen[i] < swtrt_time_lowern[i]) {
      stop("swtrt_time must be greater than or equal to swtrt_time_lower");
    }
    
    if (swtrtn[i] == 1 && swtrt_timen[i] > swtrt_time_uppern[i]) {
      stop("swtrt_time must be less than or equal to swtrt_time_upper");
    }
  }

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
  
  NumericMatrix zn_cox_den(n,p2);
  for (j=0; j<p2; j++) {
    String zj = denominator[j];
    if (!hasVariable(data, zj)) {
      stop("data must contain the variables in denominator");
    }
    if (zj == treat) {
      stop("treat should be excluded from denominator");
    }
    NumericVector u = data[zj];
    zn_cox_den(_,j) = u;
  }
  
  int q; // number of columns corresponding to the strata effects
  if (strata_main_effect_only) {
    q = sum(d - 1);
  } else {
    q = nstrata - 1;
  }
  
  if (ns_df < 0) {
    stop("ns_df must be a nonnegative integer");
  }
  
  StringVector covariates_lgs_den(q+p2+ns_df);
  NumericMatrix zn_lgs_den(n,q+p2);
  if (strata_main_effect_only) {
    k = 0;
    for (i=0; i<p_stratum; i++) {
      for (j=0; j<d[i]-1; j++) {
        covariates_lgs_den[k+j] = "stratum_" + std::to_string(i+1) +
          "_level_" + std::to_string(j+1);
        zn_lgs_den(_,k+j) = 1.0*(stratan(_,i) == j+1);
      }
      k += d[i]-1;
    }
  } else {
    for (j=0; j<nstrata-1; j++) {
      covariates_lgs_den[j] = "stratum_" + std::to_string(j+1);
      zn_lgs_den(_,j) = 1.0*(stratumn == j+1);
    }
  }
  
  for (j=0; j<p2; j++) {
    String zj = denominator[j];
    covariates_lgs_den[q+j] = zj;
    NumericVector u = data[zj];
    zn_lgs_den(_,q+j) = u;
  }
  
  for (j=0; j<ns_df; j++) {
    covariates_lgs_den[q+p2+j] = "s" + std::to_string(j+1);
  }
  
  StringVector covariates_lgs_num(q+p1+ns_df);
  for (j=0; j<q; j++) {
    covariates_lgs_num[j] = covariates_lgs_den[j];
  }
  for (j=0; j<p1; j++) {
    String zj = numerator[j];
    covariates_lgs_num[q+j] = zj;
  }
  for (j=0; j<ns_df; j++) {
    covariates_lgs_num[q+p1+j] = "s" + std::to_string(j+1);
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
  
  
  // order data by treat, id, and time
  IntegerVector order = seq(0, n-1);
  std::sort(order.begin(), order.end(), [&](int i, int j) {
    return ((treatn[i] < treatn[j]) ||
            ((treatn[i] == treatn[j]) && (idn[i] < idn[j])) ||
            ((treatn[i] == treatn[j]) && (idn[i] == idn[j]) &&
            (tstopn[i] < tstopn[j])));
  });
  
  idn = idn[order];
  stratumn = stratumn[order];
  tstartn = tstartn[order];
  tstopn = tstopn[order];
  eventn = eventn[order];
  treatn = treatn[order];
  swtrtn = swtrtn[order];
  swtrt_timen = swtrt_timen[order];
  swtrt_time_lowern = swtrt_time_lowern[order];
  swtrt_time_uppern = swtrt_time_uppern[order];
  zn = subset_matrix_by_row(zn, order);
  zn_cox_den = subset_matrix_by_row(zn_cox_den, order);
  zn_lgs_den = subset_matrix_by_row(zn_lgs_den, order);
  
  IntegerVector idx(1,0); // first observation within an id
  for (i=1; i<n; i++) {
    if (idn[i] != idn[i-1]) {
      idx.push_back(i);
    }
  }
  
  int nids = static_cast<int>(idx.size());
  idx.push_back(n);
  
  IntegerVector idx1(nids); // last observation within an id
  for (i=0; i<nids; i++) {
    idx1[i] = idx[i+1]-1;
  }
  
  IntegerVector stratumn1 = stratumn[idx1];
  IntegerVector treatn1 = treatn[idx1];
  NumericVector tstopn1 = tstopn[idx1];
  IntegerVector eventn1 = eventn[idx1];
  
  DataFrame lrdata = DataFrame::create(
    Named("stratum") = stratumn1,
    Named("treat") = treatn1,
    Named("time") = tstopn1,
    Named("event") = eventn1);
  
  DataFrame lr = lrtest(lrdata, "", "stratum", "treat", "time", "event", 
                        0, 0);
  double logRankPValue = lr["logRankPValue"];
  
  double zcrit = R::qnorm(1-alpha/2, 0, 1, 1, 0);
  
  k = -1; // indicate the observed data
  auto f = [&k, data, has_stratum, stratum, p_stratum, u_stratum, 
            treat, treatwi, treatwn, treatwc, 
            q, p, p1, p2, base_cov, numerator, denominator, 
            covariates, covariates_lgs_num, covariates_lgs_den, 
            logistic_switching_model, firth, flic, ns_df, 
            relative_time, stabilized_weights, trunc, trunc_upper_only, 
            swtrt_control_only, alpha, zcrit, ties](
                IntegerVector idb, IntegerVector stratumb,
                NumericVector tstartb, NumericVector tstopb,
                IntegerVector eventb, IntegerVector treatb,
                IntegerVector swtrtb, NumericVector swtrt_timeb,
                NumericVector swtrt_time_lowerb, 
                NumericVector swtrt_time_upperb,
                NumericMatrix zb, NumericMatrix zb_cox_den,
                NumericMatrix zb_lgs_den)->List {
                  int h, i, j, n = static_cast<int>(idb.size());
                  
                  // order data by treat, id, and time
                  IntegerVector order = seq(0, n-1);
                  std::sort(order.begin(), order.end(), [&](int i, int j) {
                    return ((treatb[i] < treatb[j]) ||
                            ((treatb[i]==treatb[j]) && (idb[i] < idb[j])) ||
                            ((treatb[i]==treatb[j]) && (idb[i]== idb[j]) &&
                            (tstopb[i] < tstopb[j])));
                  });
                  
                  idb = idb[order];
                  stratumb = stratumb[order];
                  tstartb = tstartb[order];
                  tstopb = tstopb[order];
                  eventb = eventb[order];
                  treatb = treatb[order];
                  swtrtb = swtrtb[order];
                  swtrt_timeb = swtrt_timeb[order];
                  swtrt_time_lowerb = swtrt_time_lowerb[order];
                  swtrt_time_upperb = swtrt_time_upperb[order];
                  zb = subset_matrix_by_row(zb, order);
                  zb_cox_den = subset_matrix_by_row(zb_cox_den, order);
                  zb_lgs_den = subset_matrix_by_row(zb_lgs_den, order);
                  
                  // exclude observations after treatment switch
                  IntegerVector l = which((swtrtb == 0) |
                    (tstartb < swtrt_timeb));
                  IntegerVector id1 = idb[l];
                  IntegerVector stratum1 = stratumb[l];
                  NumericVector tstart1 = tstartb[l];
                  NumericVector tstop1 = tstopb[l];
                  IntegerVector event1 = eventb[l];
                  IntegerVector treat1 = treatb[l];
                  IntegerVector swtrt1 = swtrtb[l];
                  NumericVector swtrt_time1 = swtrt_timeb[l];
                  NumericVector swtrt_time_lower1 = swtrt_time_lowerb[l];
                  NumericVector swtrt_time_upper1 = swtrt_time_upperb[l];
                  NumericMatrix z1 = subset_matrix_by_row(zb, l);
                  NumericMatrix z1_cox_den = 
                    subset_matrix_by_row(zb_cox_den, l);
                  NumericMatrix z1_lgs_den = 
                    subset_matrix_by_row(zb_lgs_den, l);
                  
                  // set up crossover and event indicators
                  int n1 = static_cast<int>(l.size());
                  IntegerVector cross1(n1);
                  for (i=0; i<n1; i++) {
                    if (i == n1-1 || id1[i] != id1[i+1]) {
                      if (swtrt1[i] == 1 && tstop1[i] >= swtrt_time1[i]) {
                        cross1[i] = 1;
                        event1[i] = 0;
                        tstop1[i] = swtrt_time1[i];
                      }
                    }
                  }
                  
                  // initialize data_switch and fit_switch
                  List data_switch(2), fit_switch(2);
                  if (k == -1) {
                    for (h=0; h<2; h++) {
                      List data_x = List::create(
                        Named("data") = R_NilValue,
                        Named(treat) = R_NilValue
                      );
                      
                      if (TYPEOF(data[treat]) == LGLSXP ||
                          TYPEOF(data[treat]) == INTSXP) {
                        data_x[treat] = treatwi[1-h];
                      } else if (TYPEOF(data[treat]) == REALSXP) {
                        data_x[treat] = treatwn[1-h];
                      } else if (TYPEOF(data[treat]) == STRSXP) {
                        data_x[treat] = treatwc[1-h];
                      }
                      
                      data_switch[h] = data_x;
                      
                      List fit_x = List::create(
                        Named("fit_den") = R_NilValue,
                        Named("fit_num") = R_NilValue,        
                        Named(treat) = R_NilValue
                      );
                      
                      if (TYPEOF(data[treat]) == LGLSXP ||
                          TYPEOF(data[treat]) == INTSXP) {
                        fit_x[treat] = treatwi[1-h];
                      } else if (TYPEOF(data[treat]) == REALSXP) {
                        fit_x[treat] = treatwn[1-h];
                      } else if (TYPEOF(data[treat]) == STRSXP) {
                        fit_x[treat] = treatwc[1-h];
                      }
                      
                      fit_switch[h] = fit_x;
                    }
                  }
                  
                  DataFrame data_outcome;
                  
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
                    NumericVector swtrt_time_lower2, swtrt_time_upper2;
                    NumericMatrix z2, z2_cox_den;
                    int n2;
                    if (!swtrt_control_only) {
                      DataFrame a = survsplit(tstart1, tstop1, cut);
                      IntegerVector censor = a["censor"];
                      l = a["row"];
                      id2 = id1[l];
                      stratum2 = stratum1[l];
                      tstart2 = a["start"];
                      tstop2 = a["end"];
                      event2 = event1[l];
                      treat2 = treat1[l];
                      cross2 = cross1[l];
                      swtrt_time_lower2 = swtrt_time_lower1[l];
                      swtrt_time_upper2 = swtrt_time_upper1[l];
                      z2 = subset_matrix_by_row(z1, l);
                      z2_cox_den = subset_matrix_by_row(z1_cox_den, l);
                      n2 = static_cast<int>(l.size());
                      for (i=0; i<n2; i++) {
                        if (censor[i] == 1) {
                          event2[i] = 0;
                          cross2[i] = 0;
                        }
                      }
                    } else {
                      // extract data for the control group
                      l = which(treat1 == 0);
                      IntegerVector id0 = id1[l];
                      IntegerVector stratum0 = stratum1[l];
                      NumericVector tstart0 = tstart1[l];
                      NumericVector tstop0 = tstop1[l];
                      IntegerVector event0 = event1[l];
                      IntegerVector treat0 = treat1[l];
                      IntegerVector cross0 = cross1[l];
                      NumericVector swtrt_time_lower0 = swtrt_time_lower1[l];
                      NumericVector swtrt_time_upper0 = swtrt_time_upper1[l];
                      NumericMatrix z0 = subset_matrix_by_row(z1, l);
                      NumericMatrix z0_cox_den = 
                        subset_matrix_by_row(z1_cox_den, l);
                      
                      // replicate event times within each subject
                      DataFrame a = survsplit(tstart0, tstop0, cut);
                      IntegerVector censor = a["censor"];
                      l = a["row"];
                      IntegerVector id20 = id0[l];
                      IntegerVector stratum20 = stratum0[l];
                      NumericVector tstart20 = a["start"];
                      NumericVector tstop20 = a["end"];
                      IntegerVector event20 = event0[l];
                      IntegerVector treat20 = treat0[l];
                      IntegerVector cross20 = cross0[l];
                      NumericVector swtrt_time_lower20 = swtrt_time_lower0[l];
                      NumericVector swtrt_time_upper20 = swtrt_time_upper0[l];
                      NumericMatrix z20 = subset_matrix_by_row(z0, l);
                      NumericMatrix z20_cox_den = 
                        subset_matrix_by_row(z0_cox_den, l);
                      int n20 = static_cast<int>(l.size());
                      for (i=0; i<n20; i++) {
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
                      NumericVector swtrt_time_lower21 = swtrt_time_lower1[l];
                      NumericVector swtrt_time_upper21 = swtrt_time_upper1[l];
                      NumericMatrix z21 = subset_matrix_by_row(z1, l);
                      NumericMatrix z21_cox_den = 
                        subset_matrix_by_row(z1_cox_den, l);
                      int n21 = static_cast<int>(l.size());
                      
                      // combine weighted control with unweighted active data
                      id2 = c_vectors_i(id20, id21);
                      stratum2 = c_vectors_i(stratum20, stratum21);
                      tstart2 = c_vectors(tstart20, tstart21);
                      tstop2 = c_vectors(tstop20, tstop21);
                      event2 = c_vectors_i(event20, event21);
                      treat2 = c_vectors_i(treat20, treat21);
                      cross2 = c_vectors_i(cross20, cross21);
                      swtrt_time_lower2 = c_vectors(swtrt_time_lower20, 
                                                    swtrt_time_lower21);
                      swtrt_time_upper2 = c_vectors(swtrt_time_upper20, 
                                                    swtrt_time_upper21);
                      z2 = c_matrices(z20, z21);
                      z2_cox_den = c_matrices(z20_cox_den, z21_cox_den);
                      n2 = n20 + n21;
                    }

                    // initialize weights
                    NumericVector w2(n2, NA_REAL), sw2(n2, NA_REAL);
                    
                    // fit the switching models by treatment group
                    for (h=0; h<K; h++) {
                      l = which((treat2 == h) & 
                        (tstop2 >= swtrt_time_lower2) & 
                        (tstop2 <= swtrt_time_upper2));
                      IntegerVector id3 = id2[l];
                      IntegerVector stratum3 = stratum2[l];
                      NumericVector tstart3 = tstart2[l];
                      NumericVector tstop3 = tstop2[l];
                      IntegerVector cross3 = cross2[l];
                      NumericMatrix z3_cox_den = 
                        subset_matrix_by_row(z2_cox_den, l);
                      int n3 = static_cast<int>(l.size());
                      
                      // prepare the data for fitting the switching model
                      DataFrame data1 = DataFrame::create(
                        Named("id") = id3,
                        Named("stratum") = stratum3,
                        Named("tstart") = tstart3,
                        Named("tstop") = tstop3,
                        Named("cross") = cross3);
                      
                      for (j=0; j<p2; j++) {
                        String zj = denominator[j];
                        NumericVector u = z3_cox_den(_,j);
                        data1.push_back(u,zj);
                      }
                      
                      // fit the denominator model for crossover
                      List fit_den = phregcpp(
                        data1, "", "stratum", "tstart", "tstop", "cross",
                        denominator, "", "", "id", ties, 
                        1, 1, 0, 0, 0, alpha);
                      
                      // obtain the survival probabilities for crossover
                      DataFrame parest_den = DataFrame(fit_den["parest"]);
                      NumericVector beta_den = parest_den["beta"];
                      NumericMatrix vbeta_den(p2, p2);
                      DataFrame basehaz_den = DataFrame(fit_den["basehaz"]);
                      
                      DataFrame km_den = survfit_phregcpp(
                        p2, beta_den, vbeta_den, basehaz_den, data1,
                        denominator, "stratum", "", "id", "tstart", "tstop", 
                        0, "log-log", 1-alpha);
                      
                      NumericVector surv_den = km_den["surv"];
                      int m = km_den.nrows();
                      
                      List fit_num = phregcpp(
                        data1, "", "stratum", "tstart", "tstop", "cross", 
                        numerator, "", "", "id", ties, 
                        1, 1, 0, 0, 0, alpha);
                      
                      NumericVector surv_num(m);
                      if (p1 > 0) {
                        DataFrame parest_num = DataFrame(fit_num["parest"]);
                        NumericVector beta_num = parest_num["beta"];
                        NumericMatrix vbeta_num(p1, p1);
                        DataFrame basehaz_num = DataFrame(fit_num["basehaz"]);
                        
                        DataFrame km_num = survfit_phregcpp(
                          p1, beta_num, vbeta_num, basehaz_num, data1,
                          numerator, "stratum", "", "id", "tstart", "tstop", 
                          0, "log-log", 1-alpha);
                        
                        surv_num = km_num["surv"];
                      } else {
                        NumericVector beta_num(1);
                        NumericMatrix vbeta_num(1,1);
                        DataFrame basehaz_num = DataFrame(fit_num["basehaz"]);
                        
                        DataFrame km_num = survfit_phregcpp(
                          p1, beta_num, vbeta_num, basehaz_num, data1,
                          numerator, "stratum", "", "id", "tstart", "tstop", 
                          0, "log-log", 1-alpha);
                        
                        surv_num = km_num["surv"];
                      }
                      
                      // unstabilized and stabilized weights
                      NumericVector w = 1.0/surv_den;
                      NumericVector sw = surv_num/surv_den;
                      
                      // truncate the weights if requested
                      if (trunc > 0.0) {
                        // truncated unstabilized weights
                        if (trunc_upper_only) {
                          double upper = quantilecpp(w, 1-trunc);
                          for (i=0; i<m; i++) {
                            if (w[i] > upper) w[i] = upper;
                          }
                        } else {
                          double lower = quantilecpp(w, trunc);
                          double upper = quantilecpp(w, 1-trunc);
                          for (i=0; i<m; i++) {
                            if (w[i] < lower) {
                              w[i] = lower;
                            } else if (w[i] > upper) {
                              w[i] = upper;
                            }
                          }
                        }
                        
                        // truncated stabilized weights
                        if (trunc_upper_only) {
                          double upper = quantilecpp(sw, 1-trunc);
                          for (i=0; i<m; i++) {
                            if (sw[i] > upper) sw[i] = upper;
                          }
                        } else {
                          double lower = quantilecpp(sw, trunc);
                          double upper = quantilecpp(sw, 1-trunc);
                          for (i=0; i<m; i++) {
                            if (sw[i] < lower) {
                              sw[i] = lower;
                            } else if (sw[i] > upper) {
                              sw[i] = upper;
                            }
                          }
                        }
                      }
                      
                      // update data_switch and fit_switch
                      if (k == -1) {
                        List data_x = List::create(
                          Named("data") = data1,
                          Named(treat) = R_NilValue
                        );
                        
                        if (TYPEOF(data[treat]) == LGLSXP ||
                            TYPEOF(data[treat]) == INTSXP) {
                          data_x[treat] = treatwi[1-h];
                        } else if (TYPEOF(data[treat]) == REALSXP) {
                          data_x[treat] = treatwn[1-h];
                        } else if (TYPEOF(data[treat]) == STRSXP) {
                          data_x[treat] = treatwc[1-h];
                        }
                        
                        data_switch[h] = data_x;
                        
                        List fit_x = List::create(
                          Named("fit_den") = fit_den,
                          Named("fit_num") = fit_num,
                          Named(treat) = R_NilValue
                        );
                        
                        if (TYPEOF(data[treat]) == LGLSXP ||
                            TYPEOF(data[treat]) == INTSXP) {
                          fit_x[treat] = treatwi[1-h];
                        } else if (TYPEOF(data[treat]) == REALSXP) {
                          fit_x[treat] = treatwn[1-h];
                        } else if (TYPEOF(data[treat]) == STRSXP) {
                          fit_x[treat] = treatwc[1-h];
                        }
                        
                        fit_switch[h] = fit_x;
                      }
                      
                      // weights in the outcome model by matching id and time
                      IntegerVector idx = km_den["id"];
                      NumericVector timex = km_den["time"];
                      double limit = 2*max(timex) + 1;
                      NumericVector idtimex(m);
                      NumericVector idtime3(n3);
                      for (i=0; i<m; i++) {
                        idtimex[i] = idx[i]*limit + timex[i];
                      }
                      for (i=0; i<n3; i++) {
                        idtime3[i] = id3[i]*limit + tstop3[i];
                      }
                      IntegerVector sub = match(idtime3, idtimex) - 1;
                      NumericVector w3 = w[sub];
                      NumericVector sw3 = sw[sub];
                      w2[l] = w3;
                      sw2[l] = sw3;
                    }
                    
                    // fill in missing weights with LOCF starting at 1
                    IntegerVector idx(1,0);
                    for (i=1; i<n2; i++) {
                      if (id2[i] != id2[i-1]) {
                        idx.push_back(i);
                      }
                    }
                    
                    int nids2 = static_cast<int>(idx.size());
                    idx.push_back(n2);
                    
                    for (i=0; i<nids2; i++) {
                      if (std::isnan(w2[idx[i]])) {
                        w2[idx[i]] = 1.0;
                        sw2[idx[i]] = 1.0;
                      }
                      for (j=idx[i]+1; j<idx[i+1]; j++) {
                        if (std::isnan(w2[j])) {
                          w2[j] = w2[j-1];
                          sw2[j] = sw2[j-1];
                        }
                      }
                    }
                    
                    // prepare data for the outcome model
                    data_outcome = DataFrame::create(
                      Named("uid") = id2,
                      Named("tstart") = tstart2,
                      Named("tstop") = tstop2,
                      Named("event") = event2,
                      Named("treated") = treat2,
                      Named("unstabilized_weight") = w2,
                      Named("stabilized_weight") = sw2);
                    
                    if (has_stratum) {
                      for (i=0; i<p_stratum; i++) {
                        String s = stratum[i];
                        if (TYPEOF(data[s]) == INTSXP) {
                          IntegerVector stratumwi = u_stratum[s];
                          data_outcome.push_back(stratumwi[stratum2-1], s);
                        } else if (TYPEOF(data[s]) == REALSXP) {
                          NumericVector stratumwn = u_stratum[s];
                          data_outcome.push_back(stratumwn[stratum2-1], s);
                        } else if (TYPEOF(data[s]) == STRSXP) {
                          StringVector stratumwc = u_stratum[s];
                          data_outcome.push_back(stratumwc[stratum2-1], s);
                        }
                      }
                    }
                    
                    for (j=0; j<p; j++) {
                      NumericVector u = z2(_,j);
                      String zj = base_cov[j];
                      data_outcome.push_back(u,zj);
                    }
                  } else { // logistic regression switching model
                    // initialize weights
                    NumericVector w1(n1, 1.0), sw1(n1, 1.0);
                    
                    // fit the switching models by treatment group
                    for (h=0; h<K; h++) {
                      l = which((treat1 == h) & 
                        (tstop1 >= swtrt_time_lower1) &
                        (tstop1 <= swtrt_time_upper1));
                      IntegerVector id2 = id1[l];
                      IntegerVector stratum2 = stratum1[l];
                      NumericVector tstop2 = tstop1[l];
                      NumericVector swtrt_time_lower2 = swtrt_time_lower1[l];
                      IntegerVector cross2 = cross1[l];
                      NumericMatrix z2_lgs_den = 
                        subset_matrix_by_row(z1_lgs_den, l);
                      int n2 = static_cast<int>(l.size());
                      
                      // obtain natural cubic spline knots
                      NumericMatrix s(n2, ns_df);
                      if (ns_df > 0) {
                        NumericVector x0(n2);
                        if (relative_time) {
                          x0 = tstop2 - swtrt_time_lower2;
                        } else {
                          x0 = tstop2;
                        }
                        NumericVector x = x0[cross2 == 1];
                        NumericVector knots(1, NA_REAL);
                        NumericVector boundary_knots(1, NA_REAL);
                        s = nscpp(x, ns_df, knots, 0, boundary_knots);
                        knots = s.attr("knots");
                        boundary_knots = s.attr("boundary_knots");
                        s = nscpp(x0, NA_INTEGER, knots, 0, boundary_knots);
                      }
                      
                      // prepare the data for fitting the switching model
                      DataFrame data1 = DataFrame::create(
                        Named("id") = id2,
                        Named("cross") = cross2);
                      
                      for (j=0; j<q+p2; j++) {
                        String zj = covariates_lgs_den[j];
                        NumericVector u = z2_lgs_den(_,j);
                        data1.push_back(u,zj);
                      }
                      for (j=0; j<ns_df; j++) {
                        String zj = covariates_lgs_den[q+p2+j];
                        NumericVector u = s(_,j);
                        data1.push_back(u,zj);
                      }
                      
                      List fit_den = logisregcpp(
                        data1, "", "cross", covariates_lgs_den, "", "", 
                        "", "id", "logit", 1, firth, 0, flic, 0, alpha);
                      
                      DataFrame f_den = DataFrame(fit_den["fitted"]);
                      NumericVector h_den = f_den["fitted_values"];
                      
                      // convert to probability of remaining uncensored
                      NumericVector s_den = 1.0 - h_den;
                      
                      // replace missing probabilities with 1 within subjects
                      l = which(treat1 == h);
                      IntegerVector id3 = id1[l];
                      NumericVector tstop3 = tstop1[l];
                      int n3 = static_cast<int>(l.size());
                      
                      // match on id and time
                      double limit = 2*max(tstop3) + 1;
                      NumericVector idtime2(n2);
                      NumericVector idtime3(n3);
                      for (i=0; i<n2; i++) {
                        idtime2[i] = id2[i]*limit + tstop2[i];
                      }
                      for (i=0; i<n3; i++) {
                        idtime3[i] = id3[i]*limit + tstop3[i];
                      }
                      IntegerVector sub = match(idtime2, idtime3) - 1;
                      NumericVector pstay_den(n3, 1.0);
                      pstay_den[sub] = s_den;
                      
                      // obtain cumulative products within a subject
                      IntegerVector idx(1,0);
                      for (i=1; i<n3; i++) {
                        if (id3[i] != id3[i-1]) {
                          idx.push_back(i);
                        }
                      }
                      
                      int nids3 = static_cast<int>(idx.size());
                      idx.push_back(n3);
                      
                      NumericVector surv_den(n3);
                      for (i=0; i<nids3; i++) {
                        surv_den[idx[i]] = pstay_den[idx[i]];
                        for (j=idx[i]+1; j<idx[i+1]; j++) {
                          surv_den[j] = surv_den[j-1]*pstay_den[j];
                        }
                      }
                      
                      List fit_num = logisregcpp(
                        data1, "", "cross", covariates_lgs_num, "", "", 
                        "", "id", "logit", 1, firth, 0, flic, 0, alpha);
                      
                      DataFrame f_num = DataFrame(fit_num["fitted"]);
                      NumericVector h_num = f_num["fitted_values"];
                      
                      // convert to probability of remaining uncensored
                      NumericVector s_num = 1.0 - h_num;
                      NumericVector pstay_num(n3, 1.0);
                      pstay_num[sub] = s_num;
                      
                      NumericVector surv_num(n3);
                      for (i=0; i<nids3; i++) {
                        surv_num[idx[i]] = pstay_num[idx[i]];
                        for (j=idx[i]+1; j<idx[i+1]; j++) {
                          surv_num[j] = surv_num[j-1]*pstay_num[j];
                        }
                      }
                      
                      // unstabilized and stabilized weights
                      NumericVector w = 1.0/surv_den;
                      NumericVector sw = surv_num/surv_den;

                      // truncate the weights if requested
                      if (trunc > 0.0) {
                        // truncated unstabilized weights
                        if (trunc_upper_only) {
                          double upper = quantilecpp(w, 1-trunc);
                          for (i=0; i<n3; i++) {
                            if (w[i] > upper) w[i] = upper;
                          }
                        } else {
                          double lower = quantilecpp(w, trunc);
                          double upper = quantilecpp(w, 1-trunc);
                          for (i=0; i<n3; i++) {
                            if (w[i] < lower) {
                              w[i] = lower;
                            } else if (w[i] > upper) {
                              w[i] = upper;
                            }
                          }
                        }
                        
                        // truncated stabilized weights
                        if (trunc_upper_only) {
                          double upper = quantilecpp(sw, 1-trunc);
                          for (i=0; i<n3; i++) {
                            if (sw[i] > upper) sw[i] = upper;
                          }
                        } else {
                          double lower = quantilecpp(sw, trunc);
                          double upper = quantilecpp(sw, 1-trunc);
                          for (i=0; i<n3; i++) {
                            if (sw[i] < lower) {
                              sw[i] = lower;
                            } else if (sw[i] > upper) {
                              sw[i] = upper;
                            }
                          }
                        }
                      }
                      
                      // update data_switch and fit_switch
                      if (k == -1) {
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
                        
                        List data_x = List::create(
                          Named("data") = data1,
                          Named(treat) = R_NilValue
                        );
                        
                        if (TYPEOF(data[treat]) == LGLSXP ||
                            TYPEOF(data[treat]) == INTSXP) {
                          data_x[treat] = treatwi[1-h];
                        } else if (TYPEOF(data[treat]) == REALSXP) {
                          data_x[treat] = treatwn[1-h];
                        } else if (TYPEOF(data[treat]) == STRSXP) {
                          data_x[treat] = treatwc[1-h];
                        }
                        
                        data_switch[h] = data_x;
                        
                        List fit_x = List::create(
                          Named("fit_den") = fit_den,
                          Named("fit_num") = fit_num,
                          Named(treat) = R_NilValue
                        );
                        
                        if (TYPEOF(data[treat]) == LGLSXP ||
                            TYPEOF(data[treat]) == INTSXP) {
                          fit_x[treat] = treatwi[1-h];
                        } else if (TYPEOF(data[treat]) == REALSXP) {
                          fit_x[treat] = treatwn[1-h];
                        } else if (TYPEOF(data[treat]) == STRSXP) {
                          fit_x[treat] = treatwc[1-h];
                        }
                        
                        fit_switch[h] = fit_x;
                      }
                      
                      w1[l] = w;
                      sw1[l] = sw;
                    }
                    
                    // prepare data for the outcome model
                    data_outcome = DataFrame::create(
                      Named("uid") = id1,
                      Named("tstart") = tstart1,
                      Named("tstop") = tstop1,
                      Named("event") = event1,
                      Named("treated") = treat1,
                      Named("unstabilized_weight") = w1,
                      Named("stabilized_weight") = sw1);
                    
                    if (has_stratum) {
                      for (i=0; i<p_stratum; i++) {
                        String s = stratum[i];
                        if (TYPEOF(data[s]) == INTSXP) {
                          IntegerVector stratumwi = u_stratum[s];
                          data_outcome.push_back(stratumwi[stratum1-1], s);
                        } else if (TYPEOF(data[s]) == REALSXP) {
                          NumericVector stratumwn = u_stratum[s];
                          data_outcome.push_back(stratumwn[stratum1-1], s);
                        } else if (TYPEOF(data[s]) == STRSXP) {
                          StringVector stratumwc = u_stratum[s];
                          data_outcome.push_back(stratumwc[stratum1-1], s);
                        }
                      }
                    }
                    
                    for (j=0; j<p; j++) {
                      NumericVector u = z1(_,j);
                      String zj = base_cov[j];
                      data_outcome.push_back(u,zj);
                    }
                  }
                  
                  // fit the outcome model with weights
                  List fit_outcome;
                  if (stabilized_weights) {
                    fit_outcome = phregcpp(
                      data_outcome, "", stratum, "tstart", "tstop",
                      "event", covariates, "stabilized_weight", "",
                      "uid", ties, 1, 1, 0, 0, 0, alpha);
                  } else {
                    fit_outcome = phregcpp(
                      data_outcome, "", stratum, "tstart", "tstop",
                      "event", covariates, "unstabilized_weight", "",
                      "uid", ties, 1, 1, 0, 0, 0, alpha);
                  }
                  
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
                      Named("data_switch") = data_switch,
                      Named("fit_switch") = fit_switch,
                      Named("data_outcome") = data_outcome,
                      Named("fit_outcome") = fit_outcome,
                      Named("hrhat") = hrhat,
                      Named("hrlower") = hrlower,
                      Named("hrupper") = hrupper,
                      Named("pvalue") = pvalue);
                  } else {
                    out = List::create(
                      Named("hrhat") = hrhat,
                      Named("hrlower") = hrlower,
                      Named("hrupper") = hrupper,
                      Named("pvalue") = pvalue);
                  }
                  
                  return out;
                };
  
  List out = f(idn, stratumn, tstartn, tstopn, eventn, treatn,
               swtrtn, swtrt_timen, swtrt_time_lowern, 
               swtrt_time_uppern, zn, zn_cox_den, zn_lgs_den);
  
  List data_switch = out["data_switch"];
  List fit_switch = out["fit_switch"];
  DataFrame data_outcome = DataFrame(out["data_outcome"]);
  List fit_outcome = out["fit_outcome"];
  double hrhat = out["hrhat"];
  double hrlower = out["hrlower"];
  double hrupper = out["hrupper"];
  double pvalue = out["pvalue"];
  
  IntegerVector uid = data_outcome["uid"];
  if (TYPEOF(data[id]) == INTSXP) {
    data_outcome.push_front(idwi[uid-1], id);
  } else if (TYPEOF(data[id]) == REALSXP) {
    data_outcome.push_front(idwn[uid-1], id);
  } else if (TYPEOF(data[id]) == STRSXP) {
    data_outcome.push_front(idwc[uid-1], id);
  }
  
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
  NumericVector hrhats(n_boot);
  if (!boot) { // use Cox model to construct CI for HR if no boot
    hr_CI_type = "Cox model";
  } else { // bootstrap the entire process to construct CI for HR
    if (seed != NA_INTEGER) set_seed(seed);
    
    int n0 = sum(treatn1==0);
    int n1 = sum(treatn1==1);
    IntegerVector nobs = diff(idx);
    int N = max(nobs)*nids;
    
    for (k=0; k<n_boot; k++) {
      IntegerVector idb(N), stratumb(N), treatb(N), eventb(N), swtrtb(N);
      NumericVector tstartb(N), tstopb(N), swtrt_timeb(N);
      NumericVector swtrt_time_lowerb(N), swtrt_time_upperb(N);
      NumericMatrix zb(N, p), zb_cox_den(N, p2), zb_lgs_den(N, q+p2);
      
      // sample the subject-level data with replacement by treatment group
      int l = 0;
      for (int h=0; h<nids; h++) {
        double u = R::runif(0,1);
        if (h < n0) {
          i = static_cast<int>(std::floor(u*n0));
        } else {
          i = n0 + static_cast<int>(std::floor(u*n1));
        }
        
        // create unique ids for bootstrap data sets
        int idb1 = idn[idx[i]] + h*nids;
        
        for (j=idx[i]; j<idx[i+1]; j++) {
          int r = l + j - idx[i];
          
          idb[r] = idb1;
          stratumb[r] = stratumn[j];
          tstartb[r] = tstartn[j];
          tstopb[r] = tstopn[j];
          eventb[r] = eventn[j];
          treatb[r] = treatn[j];
          swtrtb[r] = swtrtn[j];
          swtrt_timeb[r] = swtrt_timen[j];
          swtrt_time_lowerb[r] = swtrt_time_lowern[j];
          swtrt_time_upperb[r] = swtrt_time_uppern[j];
          zb(r,_) = zn(j,_);
          zb_cox_den(r,_) = zn_cox_den(j,_);
          zb_lgs_den(r,_) = zn_lgs_den(j,_);
        }
        
        l += idx[i+1] - idx[i];
      }
      
      IntegerVector sub = Range(0,l-1);
      idb = idb[sub];
      stratumb = stratumb[sub];
      tstartb = tstartb[sub];
      tstopb = tstopb[sub];
      eventb = eventb[sub];
      treatb = treatb[sub];
      swtrtb = swtrtb[sub];
      swtrt_timeb = swtrt_timeb[sub];
      swtrt_time_lowerb = swtrt_time_lowerb[sub];
      swtrt_time_upperb = swtrt_time_upperb[sub];
      zb = subset_matrix_by_row(zb, sub);
      zb_cox_den = subset_matrix_by_row(zb_cox_den, sub);
      zb_lgs_den = subset_matrix_by_row(zb_lgs_den, sub);
      
      out = f(idb, stratumb, tstartb, tstopb, eventb, treatb,
              swtrtb, swtrt_timeb, swtrt_time_lowerb, 
              swtrt_time_upperb, zb, zb_cox_den, zb_lgs_den);
      
      hrhats[k] = out["hrhat"];
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
  }
  
  List settings = List::create(
    Named("logistic_switching_model") = logistic_switching_model,
    Named("strata_main_effect_only") = strata_main_effect_only,
    Named("firth") = firth,
    Named("flic") = flic,
    Named("ns_df") = ns_df,
    Named("relative_time") = relative_time,
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
    Named("logrank_pvalue") = 2*std::min(logRankPValue, 1-logRankPValue),
    Named("cox_pvalue") = pvalue,
    Named("hr") = hrhat,
    Named("hr_CI") = NumericVector::create(hrlower, hrupper),
    Named("hr_CI_type") = hr_CI_type,
    Named("data_switch") = data_switch,
    Named("fit_switch") = fit_switch,
    Named("data_outcome") = data_outcome,
    Named("fit_outcome") = fit_outcome,
    Named("settings") = settings);
  
  if (boot) result.push_back(hrhats, "hr_boots");
  
  return result;
}
