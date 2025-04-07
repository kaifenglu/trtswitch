#include "utilities.h"
#include "survival_analysis.h"
#include "logistic_regression.h"

using namespace Rcpp;


List est_psi_tsegest(int n2, int q, int p2, int nids2, 
          IntegerVector idx2, IntegerVector stratumn3, 
          IntegerVector osn3, NumericVector os_timen3, 
          NumericVector censor_timen3,  
          IntegerVector swtrtn3, NumericVector swtrt_timen3, 
          IntegerVector idn2, IntegerVector y, 
          NumericVector tstartn2, NumericVector tstopn2,
          StringVector covariates_lgs, NumericMatrix zn_lgs2, 
          bool firth, bool flic, bool recensor, double alpha, 
          std::string ties, double offset, double psi) {
  
  int i, j; 
  double a = exp(psi);
  
  // counterfactual survival times and event indicators
  NumericVector t_star(nids2);
  IntegerVector d_star(nids2);
  for (i=0; i<nids2; i++) {
    double b2, u_star, c_star;
    if (swtrtn3[i] == 1) {
      b2 = swtrt_timen3[i] - offset;
      u_star = b2 + (os_timen3[i] - b2)*a;
    } else {
      u_star = os_timen3[i];
    }
    
    if (recensor) {
      c_star = censor_timen3[i]*std::min(1.0, a);
      t_star[i] = std::min(u_star, c_star);
      d_star[i] = c_star < u_star ? 0 : osn3[i];
    } else {
      t_star[i] = u_star;
      d_star[i] = osn3[i];
    }
  }
  
  // martingale residuals from the null model
  DataFrame data1 = DataFrame::create(
    Named("t_star") = t_star,
    Named("d_star") = d_star,
    Named("ustratum") = stratumn3
  );
  
  List fit1 = phregcpp(
    data1, "", "ustratum", "t_star", "", "d_star",
    "", "", "", "", ties, 0, 0, 1, 0, 0, alpha);
  
  // replicate counterfactual residuals within subjects
  NumericVector resid3 = fit1["residuals"];
  NumericVector resid(n2);
  for (i=0; i<nids2; i++) {
    for (j=idx2[i]; j<idx2[i+1]; j++) {
      resid[j] = resid3[i];
    }
  }
  
  // logistic regression switching model
  DataFrame data2 = DataFrame::create(
    Named("uid") = idn2,
    Named("tstart") = tstartn2,
    Named("tstop") = tstopn2,
    Named("cross") = y,
    Named("counterfactual") = resid);
  
  for (int j=0; j<q+p2; j++) {
    String zj = covariates_lgs[j+1];
    NumericVector u = zn_lgs2(_,j);
    data2.push_back(u, zj);
  }
  
  List fit2 = logisregcpp(
    data2, "", "cross", covariates_lgs, "", "", 
    "", "uid", "logit", 1, firth, 0, flic, 0, alpha);
  
  DataFrame parest = DataFrame(fit2["parest"]);
  NumericVector z = parest["z"];
  double z_counterfactual = z[1];
  
  List out = List::create(
    Named("data_nullcox") = data1,
    Named("fit_nullcox") = fit1,
    Named("data_logis") = data2,
    Named("fit_logis") = fit2,
    Named("z_counterfactual") = z_counterfactual);
  
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
    const std::string swtrt_time_upper = "",
    const StringVector& base_cov = "",
    const StringVector& conf_cov = "",
    const double low_psi = -2,
    const double hi_psi = 2,
    const int n_eval_z = 101,
    const bool strata_main_effect_only = 1,
    const bool firth = 0,
    const bool flic = 0,
    const bool recensor = 1,
    const bool admin_recensor_only = 1,
    const bool swtrt_control_only = 1,
    const double alpha = 0.05,
    const std::string ties = "efron",
    const double tol = 1.0e-6,
    const double offset = 1,
    const bool boot = 1,
    const int n_boot = 1000,
    const int seed = NA_INTEGER) {
  
  int i, j, k, n = data.nrow();
  
  int p = static_cast<int>(base_cov.size());
  if (p == 1 && (base_cov[0] == "" || base_cov[0] == "none")) p = 0;
  
  int p2 = static_cast<int>(conf_cov.size());
  if (p2 == 1 && (conf_cov[0] == "" || conf_cov[0] == "none")) p2 = 0;
  
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
  
  if (!has_censor_time) {
    stop("data must contain the censor_time variable"); 
  }
  
  if (TYPEOF(data[censor_time]) != INTSXP &&
      TYPEOF(data[censor_time]) != REALSXP) {
    stop("censor_time must take numeric values");
  }
  
  NumericVector censor_timenz = data[censor_time];
  NumericVector censor_timen = clone(censor_timenz);
  if (is_true(any(censor_timen < tstopn))) {
    stop("censor_time must be greater than or equal to tstop");
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
    stop("swtrt must be 1 or 0 for each subject");
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
  
  int q; // number of columns corresponding to the strata effects
  if (strata_main_effect_only) {
    q = sum(d - 1);
  } else {
    q = nstrata - 1;
  }
  
  StringVector covariates_lgs(q+p2+1);
  NumericMatrix zn_lgs(n,q+p2);
  covariates_lgs[0] = "counterfactual";
  if (strata_main_effect_only) {
    k = 0;
    for (i=0; i<p_stratum; i++) {
      for (j=0; j<d[i]-1; j++) {
        covariates_lgs[k+j+1] = "stratum_" + std::to_string(i+1) +
          "_level_" + std::to_string(j+1);
        zn_lgs(_,k+j) = 1.0*(stratan(_,i) == j+1);
      }
      k += d[i]-1;
    }
  } else {
    for (j=0; j<nstrata-1; j++) {
      covariates_lgs[j+1] = "stratum_" + std::to_string(j+1);
      zn_lgs(_,j) = 1.0*(stratumn == j+1);
    }
  }
  
  for (j=0; j<p2; j++) {
    String zj = conf_cov[j];
    if (!hasVariable(data, zj)) {
      stop("data must contain the variables in conf_cov");
    }
    if (zj == treat) {
      stop("treat should be excluded from conf_cov");
    }
    
    NumericVector u = data[zj];
    covariates_lgs[q+j+1] = zj;
    zn_lgs(_,q+j) = u;
  }
  
  if (low_psi >= hi_psi) {
    stop("low_psi must be less than hi_psi");
  }
  
  if (n_eval_z < 2) {
    stop("n_eval_z must be greater than or equal to 2");
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
  censor_timen = censor_timen[order];
  pdn = pdn[order];
  pd_timen = pd_timen[order];
  swtrtn = swtrtn[order];
  swtrt_timen = swtrt_timen[order];
  swtrt_time_uppern = swtrt_time_uppern[order];
  zn = subset_matrix_by_row(zn, order);
  zn_lgs = subset_matrix_by_row(zn_lgs, order);
  
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
  
  IntegerVector osn(n);
  NumericVector os_timen(n);
  for (i=0; i<nids; i++) {
    for (j=idx[i]; j<idx[i+1]; j++) {
      osn[j] = eventn[idx1[i]];
      os_timen[j] = tstopn[idx1[i]];
    }
  }
  
  if (is_true(any((pdn == 1) & (pd_timen > os_timen)))) {
    stop("pd_time must be less than or equal to os_time");
  }
  
  if (!admin_recensor_only) { // use the actual censoring time for dropouts
    for (i=0; i<nids; i++) {
      if (eventn[idx1[i]] == 0) {
        for (j=idx[i]; j<idx[i+1]; j++) {
          censor_timen[j] = tstopn[idx1[i]];
        }
      }
    }
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
  
  k = -1;
  auto f = [&k, data, has_stratum, stratum, p_stratum, u_stratum, 
            treat, treatwi, treatwn, treatwc, id, idwi, idwn, idwc,
            q, p, p2, covariates, covariates_lgs,
            low_psi, hi_psi, n_eval_z, firth, flic, recensor, 
            swtrt_control_only, alpha, zcrit, ties, tol, offset] (
                IntegerVector& idb, IntegerVector& stratumb, 
                NumericVector& tstartb, NumericVector& tstopb, 
                IntegerVector& eventb, IntegerVector& treatb, 
                IntegerVector& osb, NumericVector& os_timeb, 
                NumericVector& censor_timeb, 
                IntegerVector& pdb, NumericVector& pd_timeb, 
                IntegerVector& swtrtb, NumericVector& swtrt_timeb, 
                NumericVector& swtrt_time_upperb,
                NumericMatrix& zb, NumericMatrix& zb_lgs)->List {
                  int h, i, j, n = static_cast<int>(idb.size());
                  
                  // order data by treat, id, and time
                  IntegerVector order = seq(0, n-1);
                  std::sort(order.begin(), order.end(), [&](int i, int j) {
                    return ((treatb[i] < treatb[j]) ||
                            ((treatb[i] == treatb[j]) && 
                            (idb[i] < idb[j])) ||
                            ((treatb[i] == treatb[j]) && 
                            (idb[i] == idb[j]) &&
                            (tstopb[i] < tstopb[j])));
                  });
                  
                  idb = idb[order];
                  stratumb = stratumb[order];
                  tstartb = tstartb[order];
                  tstopb = tstopb[order];
                  eventb = eventb[order];
                  treatb = treatb[order];
                  osb = osb[order];
                  os_timeb = os_timeb[order];
                  censor_timeb = censor_timeb[order];
                  pdb = pdb[order];
                  pd_timeb = pd_timeb[order];
                  swtrtb = swtrtb[order];
                  swtrt_timeb = swtrt_timeb[order];
                  swtrt_time_upperb = swtrt_time_upperb[order];
                  zb = subset_matrix_by_row(zb, order);
                  zb_lgs = subset_matrix_by_row(zb_lgs, order);
                  
                  IntegerVector idx(1,0); // first observation within an id
                  for (i=1; i<n; i++) {
                    if (idb[i] != idb[i-1]) {
                      idx.push_back(i);
                    }
                  }
                  
                  int nids = static_cast<int>(idx.size());
                  idx.push_back(n);
                  
                  IntegerVector idx1(nids); // last observation within an id
                  for (i=0; i<nids; i++) {
                    idx1[i] = idx[i+1]-1;
                  }
                  
                  // one observation per subject
                  IntegerVector idn1 = idb[idx1];
                  IntegerVector stratumn1 = stratumb[idx1];
                  NumericVector timen1 = tstopb[idx1];
                  IntegerVector eventn1 = eventb[idx1];
                  IntegerVector treatn1 = treatb[idx1];
                  NumericVector censor_timen1 = censor_timeb[idx1];
                  IntegerVector swtrtn1 = swtrtb[idx1];
                  NumericVector swtrt_timen1 = swtrt_timeb[idx1];
                  NumericMatrix zn1 = subset_matrix_by_row(zb, idx1);
                  
                  // time and event adjusted for treatment switching
                  NumericVector t_star = clone(timen1);
                  IntegerVector d_star = clone(eventn1);
                  
                  double psi0hat = 0, psi0lower = 0, psi0upper = 0;
                  double psi1hat = 0, psi1lower = 0, psi1upper = 0;
                  
                  // initialize data_switch, km_switch, eval_z, 
                  // data_nullcox, fit_nullcox, data_logis, fit_logis
                  List data_switch(2), km_switch(2), eval_z(2);
                  List data_nullcox(2), fit_nullcox(2);
                  List data_logis(2), fit_logis(2);
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
                      km_switch[h] = data_x;
                      eval_z[h] = data_x;
                      data_nullcox[h] = data_x;
                      data_logis[h] = data_x;
                      
                      List fit_x = List::create(
                        Named("fit") = R_NilValue,
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
                      
                      fit_nullcox[h] = fit_x;
                      fit_logis[h] = fit_x;
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
                    // post progression data up to switching for the treat
                    IntegerVector l = which(
                      (treatb == h) & (pdb == 1) & (tstopb >= pd_timeb) & 
                        ((swtrtb != 1) | (tstopb <= swtrt_timeb)) & 
                        (tstopb <= swtrt_time_upperb));
                    
                    IntegerVector idn2 = idb[l];
                    IntegerVector stratumn2 = stratumb[l];
                    NumericVector tstartn2 = tstartb[l];
                    NumericVector tstopn2 = tstopb[l];
                    NumericVector pd_timen2 = pd_timeb[l];
                    IntegerVector osn2 = osb[l];
                    NumericVector os_timen2 = os_timeb[l];
                    NumericVector censor_timen2 = censor_timeb[l];
                    IntegerVector swtrtn2 = swtrtb[l];
                    NumericVector swtrt_timen2 = swtrt_timeb[l];
                    NumericMatrix zn_lgs2 = subset_matrix_by_row(zb_lgs, l);
                    
                    int n2 = static_cast<int>(l.size());
                    
                    // treatment switching indicators
                    IntegerVector y(n2);
                    for (i=0; i<n2; i++) {
                      y[i] = swtrtn2[i] * (tstopn2[i] == swtrt_timen2[i]);
                    }
                    
                    // re-baseline based on disease progression date
                    IntegerVector idx2(1,0);
                    for (i=1; i<n2; i++) {
                      if (idn2[i] != idn2[i-1]) {
                        idx2.push_back(i);
                      }
                    }
                    
                    int nids2 = static_cast<int>(idx2.size());
                    idx2.push_back(n2);
                    
                    IntegerVector idn3(nids2);
                    IntegerVector stratumn3(nids2);
                    IntegerVector osn3(nids2);
                    NumericVector os_timen3(nids2);
                    NumericVector censor_timen3(nids2);
                    IntegerVector swtrtn3(nids2);
                    NumericVector swtrt_timen3(nids2);
                    for (i=0; i<nids2; i++) {
                      j = idx2[i];
                      double b2 = pd_timen2[j] - offset;
                      idn3[i] = idn2[j];
                      stratumn3[i] = stratumn2[j];
                      osn3[i] = osn2[j];
                      os_timen3[i] = os_timen2[j] - b2;
                      censor_timen3[i] = censor_timen2[j] - b2;
                      swtrtn3[i] = swtrtn2[j];
                      if (swtrtn3[i] == 1) {
                        swtrt_timen3[i] = swtrt_timen2[j] - b2;
                      }
                    }
                    

                    // z-stat for the slope of counterfactual survival time
                    // martingale residuals in the logistic regression model
                    double target = 0;
                    auto g = [&target, n2, q, p2, nids2, idx2, stratumn3, 
                              osn3, os_timen3, censor_timen3, swtrtn3, 
                              swtrt_timen3, idn2, y, tstartn2, tstopn2, 
                              covariates_lgs, zn_lgs2, firth, flic, 
                              recensor, alpha, ties, 
                              offset](double psi)->double{
                                List out = est_psi_tsegest(
                                  n2, q, p2, nids2, idx2, stratumn3, 
                                  osn3, os_timen3, censor_timen3, 
                                  swtrtn3, swtrt_timen3, idn2, y, 
                                  tstartn2, tstopn2, 
                                  covariates_lgs, zn_lgs2, firth, flic, 
                                  recensor, alpha, ties, offset, psi);
                                
                                double z = out["z_counterfactual"];
                                return z - target;
                              };
                    
                    // causal parameter estimates
                    double psihat = brent(g, low_psi, hi_psi, tol);
                    double psilower = 0, psiupper = 0;
                    if (k == -1) {
                      target = zcrit;
                      psilower = brent(g, low_psi, psihat, tol);
                      target = -zcrit;
                      psiupper = brent(g, psihat, hi_psi, tol);
                    }
                    
                    // counter-factual survival times and event indicators
                    double a = exp(psihat);
                    for (i=0; i<nids; i++) {
                      if (treatn1[i] == h) {
                        double b2, u_star, c_star;
                        if (swtrtn1[i] == 1) {
                          b2 = swtrt_timen1[i] - offset;
                          u_star = b2 + (timen1[i] - b2)*a;
                        } else {
                          u_star = timen1[i];
                        }
                        
                        if (recensor) {
                          c_star = censor_timen1[i]*std::min(1.0, a);
                          t_star[i] = std::min(u_star, c_star);
                          d_star[i] = c_star < u_star ? 0 : eventn1[i];
                        } else {
                          t_star[i] = u_star;
                          d_star[i] = eventn1[i];
                        }
                      }
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
                    
                    if (k == -1) {
                      // obtain data and Kaplan-Meier plot for time to switch
                      NumericVector swtrt_timen4(nids2);
                      for (i=0; i<nids2; i++) {
                        if (swtrtn3[i] == 1) {
                          swtrt_timen4[i] = swtrt_timen3[i];
                        } else {
                          swtrt_timen4[i] = os_timen3[i];
                        }
                      }
                      
                      DataFrame data1 = DataFrame::create(
                        Named("swtrt") = swtrtn3,
                        Named("swtrt_time") = swtrt_timen4);
                      
                      if (TYPEOF(data[id]) == INTSXP) {
                        data1.push_front(idwi[idn3-1], id);
                      } else if (TYPEOF(data[id]) == REALSXP) {
                        data1.push_front(idwn[idn3-1], id);
                      } else if (TYPEOF(data[id]) == STRSXP) {
                        data1.push_front(idwc[idn3-1], id);
                      }
                      
                      if (has_stratum) {
                        for (i=0; i<p_stratum; i++) {
                          String s = stratum[i];
                          if (TYPEOF(data[s]) == INTSXP) {
                            IntegerVector stratumwi = u_stratum[s];
                            data1.push_back(stratumwi[stratumn3-1], s);
                          } else if (TYPEOF(data[s]) == REALSXP) {
                            NumericVector stratumwn = u_stratum[s];
                            data1.push_back(stratumwn[stratumn3-1], s);
                          } else if (TYPEOF(data[s]) == STRSXP) {
                            StringVector stratumwc = u_stratum[s];
                            data1.push_back(stratumwc[stratumn3-1], s);
                          }
                        }
                      }
                      
                      DataFrame km1 = kmest(data1, "", "", "swtrt_time", 
                                            "swtrt", "log-log", 1-alpha, 1);
                      
                      // obtain the Wald statistics for the coefficient of 
                      // the counterfactual in the logistic regression 
                      // switching model at a sequence of psi values
                      double step_psi = (hi_psi - low_psi)/(n_eval_z - 1);
                      NumericVector psi(n_eval_z), Z(n_eval_z);
                      for (i=0; i<n_eval_z; i++) {
                        psi[i] = low_psi + i*step_psi;
                        List out = est_psi_tsegest(
                          n2, q, p2, nids2, idx2, stratumn3, 
                          osn3, os_timen3, censor_timen3, 
                          swtrtn3, swtrt_timen3, idn2, y, 
                          tstartn2, tstopn2, 
                          covariates_lgs, zn_lgs2, firth, flic, 
                          recensor, alpha, ties, offset, psi[i]);
                        
                        Z[i] = out["z_counterfactual"];
                      }
                      
                      DataFrame data2 = DataFrame::create(
                        Named("psi") = psi,
                        Named("Z") = Z);
                      
                      // obtain data and fit for null Cox and logistic models
                      List out = est_psi_tsegest(
                        n2, q, p2, nids2, idx2, stratumn3, 
                        osn3, os_timen3, censor_timen3, 
                        swtrtn3, swtrt_timen3, idn2, y, 
                        tstartn2, tstopn2,
                        covariates_lgs, zn_lgs2, firth, flic, 
                        recensor, alpha, ties, offset, psihat);
                      
                      DataFrame data3 = DataFrame(out["data_nullcox"]);
                      DataFrame data4 = DataFrame(out["data_logis"]);
                      List fit3 = out["fit_nullcox"];
                      List fit4 = out["fit_logis"];
                      
                      if (TYPEOF(data[id]) == INTSXP) {
                        data3.push_front(idwi[idn3-1], id);
                        data4.push_front(idwi[idn2-1], id);
                      } else if (TYPEOF(data[id]) == REALSXP) {
                        data3.push_front(idwn[idn3-1], id);
                        data4.push_front(idwn[idn2-1], id);
                      } else if (TYPEOF(data[id]) == STRSXP) {
                        data3.push_front(idwc[idn3-1], id);
                        data4.push_front(idwc[idn2-1], id);
                      }
                      
                      if (has_stratum) {
                        for (i=0; i<p_stratum; i++) {
                          String s = stratum[i];
                          if (TYPEOF(data[s]) == INTSXP) {
                            IntegerVector stratumwi = u_stratum[s];
                            data3.push_back(stratumwi[stratumn3-1], s);
                            data4.push_back(stratumwi[stratumn2-1], s);
                          } else if (TYPEOF(data[s]) == REALSXP) {
                            NumericVector stratumwn = u_stratum[s];
                            data3.push_back(stratumwn[stratumn3-1], s);
                            data4.push_back(stratumwn[stratumn2-1], s);
                          } else if (TYPEOF(data[s]) == STRSXP) {
                            StringVector stratumwc = u_stratum[s];
                            data3.push_back(stratumwc[stratumn3-1], s);
                            data4.push_back(stratumwc[stratumn2-1], s);
                          }
                        }
                      }
                      
                      // update the data and model fits
                      List data_1 = data_switch[h];
                      data_1["data"] = data1;
                      data_switch[h] = data_1;
                      
                      List data_2 = clone(data_1);
                      data_2["data"] = data2;
                      eval_z[h] = data_2;
                      
                      List data_3 = clone(data_1);
                      data_3["data"] = data3;
                      data_nullcox[h] = data_3;
                      
                      List data_4 = clone(data_1);
                      data_4["data"] = data4;
                      data_logis[h] = data_4;
                      
                      List data_5 = clone(data_1);
                      data_5["data"] = km1;
                      km_switch[h] = data_5;
                      
                      List fit_3 = fit_nullcox[h];
                      fit_3["fit"] = fit3;
                      fit_nullcox[h] = fit_3;
                      
                      List fit_4 = clone(fit_3);
                      fit_4["fit"] = fit4;
                      fit_logis[h] = fit_4;
                    }
                  }
                  
                  // Cox model for hypothetical treatment effect estimate
                  data_outcome = DataFrame::create(
                    Named("uid") = idn1,
                    Named("t_star") = t_star,
                    Named("d_star") = d_star,
                    Named("treated") = treatn1);
                  
                  data_outcome.push_back(stratumn1, "ustratum");
                  
                  for (j=0; j<p; j++) {
                    String zj = covariates[j+1];
                    NumericVector u = zn1(_,j);
                    data_outcome.push_back(u, zj);
                  }
                  
                  List fit_outcome = phregcpp(
                    data_outcome, "", "ustratum", "t_star", "", "d_star",
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
                      Named("psi1hat") = psi1hat,
                      Named("hrhat") = hrhat,
                      Named("hrlower") = hrlower,
                      Named("hrupper") = hrupper,
                      Named("pvalue") = pvalue);
                  }
                  
                  return out;
                };
  
  List out = f(idn, stratumn, tstartn, tstopn, eventn, treatn,
               osn, os_timen, censor_timen, pdn, pd_timen, 
               swtrtn, swtrt_timen, swtrt_time_uppern, 
               zn, zn_lgs);
  
  List data_switch = out["data_switch"];
  List km_switch = out["km_switch"];
  List eval_z = out["eval_z"];
  List data_nullcox = out["data_nullcox"];
  List fit_nullcox = out["fit_nullcox"];
  List data_logis = out["data_logis"];
  List fit_logis = out["fit_logis"];
  DataFrame data_outcome = DataFrame(out["data_outcome"]);
  List fit_outcome = out["fit_outcome"];
  
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
  
  if (has_stratum) {
    IntegerVector ustratum = data_outcome["ustratum"];
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
  double psilower = out["psilower"];
  double psiupper = out["psiupper"];
  double psi1hat = out["psi1hat"];
  double psi1lower = out["psi1lower"];
  double psi1upper = out["psi1upper"];
  double hrhat = out["hrhat"];
  double hrlower = out["hrlower"];
  double hrupper = out["hrupper"];
  double pvalue = out["pvalue"];
  String psi_CI_type = "logistic model";
  
  // construct the confidence interval for HR
  String hr_CI_type;
  NumericVector hrhats(n_boot), psihats(n_boot), psi1hats(n_boot);
  if (!boot) { // use Cox model to construct CI for HR if no boot
    hr_CI_type = "Cox model";
  } else { // bootstrap the entire process to construct CI for HR
    if (seed != NA_INTEGER) set_seed(seed);
    
    int n0 = sum(treatn1 == 0);
    int n1 = sum(treatn1 == 1);
    IntegerVector nobs = diff(idx);
    int N = max(nobs)*nids;
    
    for (k=0; k<n_boot; k++) {
      IntegerVector idb(N), stratumb(N), eventb(N), treatb(N);
      IntegerVector osb(N), pdb(N), swtrtb(N);
      NumericVector tstartb(N), tstopb(N), os_timeb(N), censor_timeb(N);
      NumericVector pd_timeb(N), swtrt_timeb(N), swtrt_time_upperb(N);
      NumericMatrix zb(N, p), zb_lgs(N, q+p2);
      
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
          osb[r] = osn[j];
          os_timeb[r] = os_timen[j];
          censor_timeb[r] = censor_timen[j];
          pdb[r] = pdn[j];
          pd_timeb[r] = pd_timen[j];
          swtrtb[r] = swtrtn[j];
          swtrt_timeb[r] = swtrt_timen[j];
          swtrt_time_upperb[r] = swtrt_time_uppern[j];
          zb(r,_) = zn(j,_);
          zb_lgs(r,_) = zn_lgs(j,_);
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
      osb = osb[sub];
      os_timeb = os_timeb[sub];
      censor_timeb = censor_timeb[sub];
      pdb = pdb[sub];
      pd_timeb = pd_timeb[sub];
      swtrtb = swtrtb[sub];
      swtrt_timeb = swtrt_timeb[sub];
      swtrt_time_upperb = swtrt_time_upperb[sub];
      zb = subset_matrix_by_row(zb, sub);
      zb_lgs = subset_matrix_by_row(zb_lgs, sub);
      
      out = f(idb, stratumb, tstartb, tstopb, eventb, treatb,
              osb, os_timeb, censor_timeb, pdb, pd_timeb, 
              swtrtb, swtrt_timeb, swtrt_time_upperb, 
              zb, zb_lgs);
      
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
  
  List analysis_switch = List::create(
    Named("data_switch") = data_switch,
    Named("km_switch") = km_switch,
    Named("eval_z") = eval_z,
    Named("data_nullcox") = data_nullcox,
    Named("fit_nullcox") = fit_nullcox,
    Named("data_logis") = data_logis,
    Named("fit_logis") = fit_logis
  );
  
  List settings = List::create(
    Named("low_psi") = low_psi,
    Named("hi_psi") = hi_psi,
    Named("n_eval_z") = n_eval_z,
    Named("strata_main_effect_only") = strata_main_effect_only,
    Named("firth") = firth,
    Named("flic") = flic,
    Named("recensor") = recensor,
    Named("admin_recensor_only") = admin_recensor_only,
    Named("swtrt_control_only") = swtrt_control_only,
    Named("alpha") = alpha,
    Named("ties") = ties,
    Named("tol") = tol,
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
    Named("analysis_switch") = analysis_switch,
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
