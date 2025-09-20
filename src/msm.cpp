#include "utilities.h"
#include "survival_analysis.h"
#include "logistic_regression.h"
#include "splines.h"

using namespace Rcpp;

// [[Rcpp::export]]
List msmcpp(
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
    const bool strata_main_effect_only = 1,
    const bool firth = 0,
    const bool flic = 0,
    const int ns_df = 3,
    const bool stabilized_weights = 1,
    const double trunc = 0,
    const bool trunc_upper_only = 1,
    const bool swtrt_control_only = 1,
    const bool treat_alt_interaction = 0,
    const double alpha = 0.05,
    const std::string ties = "efron",
    const bool boot = 0,
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
  }
  
  if (p > 0) {
    if (p2 == 0 || is_true(any(is_na(match(base_cov, denominator))))) {
      stop("base_cov must be a subset of denominator");
    }
  }
  
  bool has_id = hasVariable(data, id);
  bool has_tstart = hasVariable(data, tstart);
  bool has_tstop = hasVariable(data, tstop);
  bool has_event = hasVariable(data, event);
  bool has_treat = hasVariable(data, treat);
  bool has_swtrt = hasVariable(data, swtrt);
  bool has_swtrt_time = hasVariable(data, swtrt_time);
  
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
  List levels(p_stratum);
  if (p_stratum == 1 && (stratum[0] == "" || stratum[0] == "none")) {
    has_stratum = 0;
    stratumn.fill(1);
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
  
  IntegerVector eventn(n);
  if (TYPEOF(data[event]) == LGLSXP || TYPEOF(data[event]) == INTSXP) {
    IntegerVector eventnz = data[event];
    if (is_true(any((eventnz != 1) & (eventnz != 0)))) {
      stop("event must be 1 or 0 for each subject");
    } else {
      eventn = clone(eventnz);
    }
  } else if (TYPEOF(data[event]) == REALSXP) {
    NumericVector eventnz = data[event];
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
  
  IntegerVector swtrtn(n);
  if (TYPEOF(data[swtrt]) == LGLSXP || TYPEOF(data[swtrt]) == INTSXP) {
    IntegerVector swtrtnz = data[swtrt];
    if (is_true(any((swtrtnz != 1) & (swtrtnz != 0)))) {
      stop("swtrt must be 1 or 0 for each subject");
    } else {
      swtrtn = clone(swtrtnz);
    }
  } else if (TYPEOF(data[swtrt]) == REALSXP) {
    NumericVector swtrtnz = data[swtrt];
    if (is_true(any((swtrtnz != 1) & (swtrtnz != 0)))) {
      stop("swtrt must be 1 or 0 for each subject");
    } else {
      NumericVector swtrtnz2 = clone(swtrtnz);
      swtrtn = as<IntegerVector>(swtrtnz2);
    }
  } else {
    stop("swtrt must take logical, integer, or real values");
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
  
  StringVector covariates(p+2);
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
    covariates[j+1] = zj;
    NumericVector u = data[zj];
    zn(_,j) = u;
  }
  covariates[p+1] = "crossed";
  if (!swtrt_control_only && treat_alt_interaction) {
    covariates.push_back("treated_crossed");
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
    for (i=0; i<p_stratum; ++i) {
      int di = d[i]-1;
      for (j=0; j<di; ++j) {
        covariates_lgs_den[k+j] = as<std::string>(stratum[i]);
        if (TYPEOF(levels[i]) == STRSXP) {
          StringVector u = levels[i];
          std::string label = sanitize(as<std::string>(u[j]));
          covariates_lgs_den[k+j] += label;
        } else if (TYPEOF(levels[i]) == REALSXP) {
          NumericVector u = levels[i];
          covariates_lgs_den[k+j] += std::to_string(u[j]);
        } else if (TYPEOF(levels[i]) == INTSXP 
                     || TYPEOF(levels[i]) == LGLSXP) {
          IntegerVector u = levels[i];
          covariates_lgs_den[k+j] += std::to_string(u[j]);
        }
        zn_lgs_den(_,k+j) = (stratan(_,i) == j+1);
      }
      k += di;
    }
  } else {
    for (j=0; j<nstrata-1; ++j) {
      // locate the first observation in the stratum
      int first_k = 0;
      for (; first_k<n; ++first_k) {
        if (stratumn[first_k] == j+1) break;
      }
      covariates_lgs_den[j] = "";
      for (i=0; i<p_stratum; ++i) {
        IntegerVector q_col = stratan(_,i);
        int l = q_col[first_k] - 1;
        covariates_lgs_den[j] += as<std::string>(stratum[i]);
        if (TYPEOF(levels[i]) == STRSXP) {
          StringVector u = levels[i];
          std::string label = sanitize(as<std::string>(u[l]));
          covariates_lgs_den[j] += label;
        } else if (TYPEOF(levels[i]) == REALSXP) {
          NumericVector u = levels[i];
          covariates_lgs_den[j] += std::to_string(u[l]);
        } else if (TYPEOF(levels[i]) == INTSXP 
                     || TYPEOF(levels[i]) == LGLSXP) {
          IntegerVector u = levels[i];
          covariates_lgs_den[j] += std::to_string(u[l]);
        }
        if (i < p_stratum-1) {
          covariates_lgs_den[j] += ".";
        }
      }
      zn_lgs_den(_,j) = (stratumn == j+1);
    }
  }
  
  for (j=0; j<p2; j++) {
    String zj = denominator[j];
    covariates_lgs_den[q+j] = zj;
    NumericVector u = data[zj];
    zn_lgs_den(_,q+j) = u;
  }
  
  for (j=0; j<ns_df; j++) {
    covariates_lgs_den[q+p2+j] = "ns" + std::to_string(j+1);
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
  
  
  // split at treatment switching into two observations if treatment 
  // switching occurs strictly between tstart and tstop for a subject
  LogicalVector tosplit(n);
  for (i=0; i<n; i++) {
    tosplit[i] = swtrtn[i] == 1 && swtrt_timen[i] > tstartn[i] && 
      swtrt_timen[i] < tstopn[i] ? 1 : 0;
  }
  
  k = sum(tosplit);
  if (k > 0) {
    // copy old matrices to new matrices
    NumericMatrix zn1(n+k, zn.ncol());
    NumericMatrix zn_lgs_den1(n+k, zn_lgs_den.ncol());
    for (i=0; i<n; i++) {
      zn1(i,_) = zn(i,_);
      zn_lgs_den1(i,_) = zn_lgs_den(i,_);
    }
    
    IntegerVector sub = which(tosplit);
    for (i=0; i<k; i++) {
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
      zn_lgs_den1(n+i,_) = zn_lgs_den(l,_);
      
      // change tstop and event for the old observation
      tstopn[l] = swtrt_timen[l];
      eventn[l] = 0;
    }
    
    // update number of rows and old matrices
    n = n + k;
    zn = zn1;
    zn_lgs_den = zn_lgs_den1;
  }
  
  
  // sort data by treatment group, stratum, id, and time
  IntegerVector order = seq(0, n-1);
  if (has_stratum) {
    std::sort(order.begin(), order.end(), [&](int i, int j) {
      return ((treatn[i] < treatn[j]) ||
              ((treatn[i] == treatn[j]) && (stratumn[i] < stratumn[j])) ||
              ((treatn[i] == treatn[j]) && (stratumn[i] == stratumn[j]) &&
              (idn[i] < idn[j])) ||
              ((treatn[i] == treatn[j]) && (stratumn[i] == stratumn[j]) &&
              (idn[i] == idn[j]) && (tstopn[i] < tstopn[j])));
    });
  } else {
    std::sort(order.begin(), order.end(), [&](int i, int j) {
      return ((treatn[i] < treatn[j]) ||
              ((treatn[i] == treatn[j]) && (idn[i] < idn[j])) ||
              ((treatn[i] == treatn[j]) && (idn[i] == idn[j]) && 
              (tstopn[i] < tstopn[j])));
    });
  }
  
  idn = idn[order];
  stratumn = stratumn[order];
  tstartn = tstartn[order];
  tstopn = tstopn[order];
  eventn = eventn[order];
  treatn = treatn[order];
  swtrtn = swtrtn[order];
  swtrt_timen = swtrt_timen[order];
  zn = subset_matrix_by_row(zn, order);
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
  
  IntegerVector osn(n);
  NumericVector os_timen(n);
  for (i=0; i<nids; i++) {
    k = idx1[i];
    for (j=idx[i]; j<idx[i+1]; j++) {
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
  
  DataFrame lrdata = DataFrame::create(
    Named("stratum") = stratumn1,
    Named("treat") = treatn1,
    Named("time") = tstopn1,
    Named("event") = eventn1);
  
  DataFrame lr = lrtest(lrdata, "", "stratum", "treat", "time", "event",0,0);
  double logRankPValue = lr["logRankPValue"];
  
  double zcrit = R::qnorm(1-alpha/2, 0, 1, 1, 0);
  
  k = -1; // indicate the observed data
  auto f = [&k, data, has_stratum, stratum, p_stratum, u_stratum, 
            treat, treatwi, treatwn, treatwc, id, idwi, idwn, idwc,
            q, p, p2, base_cov, numerator, denominator, 
            covariates, covariates_lgs_num, covariates_lgs_den, 
            firth, flic, ns_df, 
            stabilized_weights, trunc, trunc_upper_only, 
            swtrt_control_only, treat_alt_interaction, alpha, zcrit, ties](
                IntegerVector& idb, IntegerVector& stratumb,
                NumericVector& tstartb, NumericVector& tstopb,
                IntegerVector& eventb, IntegerVector& treatb,
                NumericVector& os_timeb, 
                IntegerVector& swtrtb, NumericVector& swtrt_timeb,
                NumericMatrix& zb, 
                NumericMatrix& zb_lgs_den)->List {
                  int h, i, j, n = static_cast<int>(idb.size());
                  bool fail = 0; // whether any model fails to converge
                  NumericVector init(1, NA_REAL);
                  
                  // exclude observations after treatment switch for logistic
                  LogicalVector c1 = ifelse(
                    swtrtb == 1, tstartb < swtrt_timeb, tstopb < os_timeb);
                  IntegerVector l = which(c1);
                  IntegerVector id1 = idb[l];
                  IntegerVector stratum1 = stratumb[l];
                  NumericVector tstart1 = tstartb[l];
                  NumericVector tstop1 = tstopb[l];
                  IntegerVector treat1 = treatb[l];
                  IntegerVector swtrt1 = swtrtb[l];
                  NumericVector swtrt_time1 = swtrt_timeb[l];
                  NumericMatrix z1_lgs_den = 
                    subset_matrix_by_row(zb_lgs_den, l);
                  int n1 = static_cast<int>(l.size());
                  
                  // set up crossover indicators
                  IntegerVector cross1(n1);
                  for (i=0; i<n1; i++) {
                    if (i == n1-1 || id1[i] != id1[i+1]) {
                      if (swtrt1[i] == 1 && tstop1[i] >= swtrt_time1[i]) {
                        cross1[i] = 1;
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
                  
                  // initialize weights
                  NumericVector wb(n, 1.0), swb(n, 1.0);
                  
                  // fit the switching models by treatment group
                  for (h=0; h<K; h++) {
                    IntegerVector l = which(treat1 == h);
                    IntegerVector id2 = id1[l];
                    IntegerVector stratum2 = stratum1[l];
                    NumericVector tstart2 = tstart1[l];
                    NumericVector tstop2 = tstop1[l];
                    IntegerVector cross2 = cross1[l];
                    NumericMatrix z2_lgs_den = 
                      subset_matrix_by_row(z1_lgs_den, l);
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
                    DataFrame data1 = DataFrame::create(
                      Named("uid") = id2,
                      Named("ustratum") = stratum2,
                      Named("tstart") = tstart2,
                      Named("tstop") = tstop2,
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
                      "", "uid", "logit", init, 
                      0, firth, flic, 0, alpha, 50, 1.0e-9);
                    
                    DataFrame sumstat_den = DataFrame(fit_den["sumstat"]);
                    bool fail_den = sumstat_den["fail"];
                    if (fail_den == 1) fail = 1;
                    
                    DataFrame f_den = DataFrame(fit_den["fitted"]);
                    NumericVector h_den = f_den["fitted_values"];
                    
                    List fit_num = logisregcpp(
                      data1, "", "cross", covariates_lgs_num, "", "", 
                      "", "uid", "logit", init, 
                      0, firth, flic, 0, alpha, 50, 1.0e-9);
                    
                    DataFrame sumstat_num = DataFrame(fit_num["sumstat"]);
                    bool fail_num = sumstat_num["fail"];
                    if (fail_num == 1) fail = 1;
                    
                    DataFrame f_num = DataFrame(fit_num["fitted"]);
                    NumericVector h_num = f_num["fitted_values"];
                    
                    // convert to probability of observed response 
                    NumericVector o_den(n2), o_num(n2);
                    for (i=0; i<n2; i++) {
                      o_den[i] = cross2[i] == 1 ? h_den[i] : 1 - h_den[i];
                      o_num[i] = cross2[i] == 1 ? h_num[i] : 1 - h_num[i];
                    }
                    
                    // obtain cumulative products within a subject
                    IntegerVector idx2(1,0);
                    for (i=1; i<n2; i++) {
                      if (id2[i] != id2[i-1]) {
                        idx2.push_back(i);
                      }
                    }
                    
                    idx2.push_back(n2);
                    
                    l = which(treatb == h);
                    IntegerVector id3 = idb[l];
                    IntegerVector swtrt3 = swtrtb[l];
                    int n3 = static_cast<int>(l.size());
                    
                    IntegerVector idx3(1,0);
                    for (i=1; i<n3; i++) {
                      if (id3[i] != id3[i-1]) {
                        idx3.push_back(i);
                      }
                    }
                    
                    IntegerVector swtrt3u = swtrt3[idx3];
                    
                    int nids3 = static_cast<int>(idx3.size());
                    idx3.push_back(n3);
                    
                    NumericVector p_den(n3, 1.0), p_num(n3, 1.0);
                    
                    int m = 0; // cum obs prior to current id2
                    int v = 0; // index for current unique id2
                    for (i=0; i<nids3; i++) {
                      int mi = swtrt3u[i] == 1 ? idx2[v+1] - idx2[v] : 
                      idx3[i+1] - idx3[i] - 1;
                      int r = m - idx3[i] - 1;
                      if (swtrt3u[i] == 1) {
                        // cum prod before switch
                        int jj = std::min(idx3[i]+mi, idx3[i+1]-1);
                        for (j=idx3[i]+1; j<=jj; j++) {
                          p_den[j] = p_den[j-1]*o_den[r+j];
                          p_num[j] = p_num[j-1]*o_num[r+j];
                        }
                        // LOCF after switch
                        for (j=jj+1; j<idx3[i+1]; j++) {
                          p_den[j] = p_den[j-1];
                          p_num[j] = p_num[j-1];
                        }
                      } else {
                        for (j=idx3[i]+1; j<idx3[i+1]; j++) {
                          p_den[j] = p_den[j-1]*o_den[r+j];
                          p_num[j] = p_num[j-1]*o_num[r+j];
                        }
                      }
                      
                      m += mi;
                      if (mi > 0) v++;
                    }
                    
                    // unstabilized and stabilized weights
                    NumericVector w = 1.0/p_den;
                    NumericVector sw = p_num/p_den;
                    
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
                      IntegerVector uid = data1["uid"];
                      if (TYPEOF(data[id]) == INTSXP) {
                        data1.push_front(idwi[uid-1], id);
                      } else if (TYPEOF(data[id]) == REALSXP) {
                        data1.push_front(idwn[uid-1], id);
                      } else if (TYPEOF(data[id]) == STRSXP) {
                        data1.push_front(idwc[uid-1], id);
                      }
                      
                      if (has_stratum) {
                        IntegerVector ustratum = data1["ustratum"];
                        for (i=0; i<p_stratum; i++) {
                          String s = stratum[i];
                          if (TYPEOF(data[s]) == INTSXP) {
                            IntegerVector stratumwi = u_stratum[s];
                            data1.push_back(stratumwi[ustratum-1], s);
                          } else if (TYPEOF(data[s]) == REALSXP) {
                            NumericVector stratumwn = u_stratum[s];
                            data1.push_back(stratumwn[ustratum-1], s);
                          } else if (TYPEOF(data[s]) == STRSXP) {
                            StringVector stratumwc = u_stratum[s];
                            data1.push_back(stratumwc[ustratum-1], s);
                          }
                        }
                      }
                      
                      List data_x = data_switch[h];
                      data_x["data"] = data1;
                      data_switch[h] = data_x;
                      
                      List fit_x = fit_switch[h];
                      fit_x["fit_den"] = fit_den;
                      fit_x["fit_num"] = fit_num;
                      fit_switch[h] = fit_x;
                    }
                    
                    wb[l] = w;
                    swb[l] = sw;
                  }
                  
                  // set up time-dependent treatment switching indicators
                  IntegerVector crossb(n);
                  for (i=0; i<n; i++) {
                    if (swtrtb[i] == 1 && tstartb[i] >= swtrt_timeb[i]) {
                      crossb[i] = 1;
                    } else {
                      crossb[i] = 0;
                    }
                  }
                  
                  // prepare data for the outcome model
                  data_outcome = DataFrame::create(
                    Named("uid") = idb,
                    Named("tstart") = tstartb,
                    Named("tstop") = tstopb,
                    Named("event") = eventb,
                    Named("treated") = treatb,
                    Named("crossed") = crossb,
                    Named("unstabilized_weight") = wb,
                    Named("stabilized_weight") = swb);
                  
                  if (!swtrt_control_only && treat_alt_interaction) {
                    IntegerVector treat_cross = treatb*crossb;
                    data_outcome.push_back(treat_cross, "treated_crossed");
                  }
                  
                  data_outcome.push_back(stratumb, "ustratum");
                  
                  for (j=0; j<p; j++) {
                    NumericVector u = zb(_,j);
                    String zj = base_cov[j];
                    data_outcome.push_back(u,zj);
                  }
                  
                  // fit the outcome model with weights
                  List fit_outcome;
                  if (stabilized_weights) {
                    fit_outcome = phregcpp(
                      data_outcome, "", "ustratum", "tstart", "tstop",
                      "event", covariates, "stabilized_weight", "",
                      "uid", ties, init, 1, 0, 0, 0, 0, alpha, 50, 1.0e-9);
                  } else {
                    fit_outcome = phregcpp(
                      data_outcome, "", "ustratum", "tstart", "tstop",
                      "event", covariates, "unstabilized_weight", "",
                      "uid", ties, init, 1, 0, 0, 0, 0, alpha, 50, 1.0e-9);
                  }
                  
                  DataFrame sumstat_cox = DataFrame(fit_outcome["sumstat"]);
                  bool fail_cox = sumstat_cox["fail"];
                  if (fail_cox == 1) fail = 1;
                  
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
               swtrtn, swtrt_timen, zn, zn_lgs_den);
  
  List data_switch = out["data_switch"];
  List fit_switch = out["fit_switch"];
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
  
  
  double hrhat = out["hrhat"];
  double hrlower = out["hrlower"];
  double hrupper = out["hrupper"];
  double pvalue = out["pvalue"];
  bool fail = out["fail"];
  
  // construct the confidence interval for HR
  String hr_CI_type;
  NumericVector hrhats(n_boot);
  LogicalVector fails(n_boot);
  DataFrame fail_boots_data;
  if (!boot) { // use Cox model to construct CI for HR if no boot
    hr_CI_type = "Cox model";
  } else { // bootstrap the entire process to construct CI for HR
    if (seed != NA_INTEGER) set_seed(seed);
    
    IntegerVector nobs = diff(idx);
    int N = max(nobs)*nids;
    
    int B = N*n_boot;
    IntegerVector boot_indexc(B);
    IntegerVector oidc(B);
    IntegerVector idc(B), stratumc(B), eventc(B), treatc(B);
    IntegerVector osc(B), swtrtc(B);
    NumericVector tstartc(B), tstopc(B), os_timec(B), swtrt_timec(B);
    NumericMatrix zc_lgs_den(B, q+p2);
    int index1 = 0;
    
    IntegerVector tsx(1,0); // first id within each treat/stratum
    for (i=1; i<nids; i++) {
      if ((treatn1[i] != treatn1[i-1]) || 
          ((treatn1[i] == treatn1[i-1]) && 
          (stratumn1[i] != stratumn1[i-1]))) {
        tsx.push_back(i);
      }
    }
    
    int ntss = static_cast<int>(tsx.size());
    tsx.push_back(nids);
    
    for (k=0; k<n_boot; k++) {
      IntegerVector oidb(N);
      IntegerVector idb(N), stratumb(N), eventb(N), treatb(N);
      IntegerVector osb(N), swtrtb(N);
      NumericVector tstartb(N), tstopb(N), os_timeb(N), swtrt_timeb(N);
      NumericMatrix zb(N, p), zb_lgs_den(N, q+p2);
      
      // sample the subject-level data with replacement by treat/stratum
      int l = 0;
      for (int h=0; h<ntss; h++) {
        for (int r=tsx[h]; r<tsx[h+1]; r++) {
          double u = R::runif(0,1);
          i = tsx[h] + static_cast<int>(std::floor(u*(tsx[h+1]-tsx[h])));
          
          // create unique ids for bootstrap data sets
          int oidb1 = idn[idx[i]];
          int idb1 = idn[idx[i]] + r*nids;
          
          for (j=idx[i]; j<idx[i+1]; j++) {
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
            zb_lgs_den(l,_) = zn_lgs_den(j,_);
            l++;
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
      zb_lgs_den = subset_matrix_by_row(zb_lgs_den, sub);
      
      out = f(idb, stratumb, tstartb, tstopb, eventb, treatb, os_timeb,
              swtrtb, swtrt_timeb, zb, zb_lgs_den);
      
      fails[k] = out["fail"];
      hrhats[k] = out["hrhat"];
      
      if (fails[k]) {
        for (i=0; i<l; i++) {
          j = index1 + i;
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
          zc_lgs_den(j,_) = zb_lgs_den(i,_);
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
      zc_lgs_den = subset_matrix_by_row(zc_lgs_den,sub);
      
      fail_boots_data = DataFrame::create(
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
      
      for (j=0; j<q+p2; j++) {
        String zj = covariates_lgs_den[j];
        NumericVector u = zc_lgs_den(_,j);
        fail_boots_data.push_back(u, zj);
      }
      
      if (TYPEOF(data[id]) == INTSXP) {
        fail_boots_data.push_back(idwi[oidc-1], id);
      } else if (TYPEOF(data[id]) == REALSXP) {
        fail_boots_data.push_back(idwn[oidc-1], id);
      } else if (TYPEOF(data[id]) == STRSXP) {
        fail_boots_data.push_back(idwc[oidc-1], id);
      }
      
      if (TYPEOF(data[treat]) == LGLSXP || TYPEOF(data[treat]) == INTSXP) {
        fail_boots_data.push_back(treatwi[1-treatc], treat);
      } else if (TYPEOF(data[treat]) == REALSXP) {
        fail_boots_data.push_back(treatwn[1-treatc], treat);
      } else if (TYPEOF(data[treat]) == STRSXP) {
        fail_boots_data.push_back(treatwc[1-treatc], treat);
      }
      
      if (has_stratum) {
        for (i=0; i<p_stratum; i++) {
          String s = stratum[i];
          if (TYPEOF(data[s]) == INTSXP) {
            IntegerVector stratumwi = u_stratum[s];
            fail_boots_data.push_back(stratumwi[stratumc-1], s);
          } else if (TYPEOF(data[s]) == REALSXP) {
            NumericVector stratumwn = u_stratum[s];
            fail_boots_data.push_back(stratumwn[stratumc-1], s);
          } else if (TYPEOF(data[s]) == STRSXP) {
            StringVector stratumwc = u_stratum[s];
            fail_boots_data.push_back(stratumwc[stratumc-1], s);
          }
        }
      }
    }
    
    // obtain bootstrap confidence interval for HR
    double loghr = log(hrhat);
    LogicalVector ok = (1 - fails) & !is_na(hrhats);
    int n_ok = sum(ok);
    NumericVector loghrs = log(hrhats[ok]);
    double sdloghr = sd(loghrs);
    double tcrit = R::qt(1-alpha/2, n_ok-1, 1, 0);
    hrlower = exp(loghr - tcrit*sdloghr);
    hrupper = exp(loghr + tcrit*sdloghr);
    hr_CI_type = "bootstrap";
    pvalue = 2*(1 - R::pt(fabs(loghr/sdloghr), n_ok-1, 1, 0));
  }
  
  List settings = List::create(
    Named("strata_main_effect_only") = strata_main_effect_only,
    Named("firth") = firth,
    Named("flic") = flic,
    Named("ns_df") = ns_df,
    Named("stabilized_weights") = stabilized_weights,
    Named("trunc") = trunc,
    Named("trunc_upper_only") = trunc_upper_only,
    Named("swtrt_control_only") = swtrt_control_only,
    Named("treat_alt_interaction") = treat_alt_interaction,
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
    Named("fail") = fail,
    Named("settings") = settings);
  
  if (boot) {
    result.push_back(fails, "fail_boots");
    result.push_back(hrhats, "hr_boots"); 
    if (is_true(any(fails))) {
      result.push_back(fail_boots_data, "fail_boots_data");
    }
  }
  
  return result;
}
