#include <RcppParallel.h>
#include <RcppThread.h>
#include <Rcpp.h>

#include <boost/random.hpp>

#include "survival_analysis.h"
#include "logistic_regression.h"
#include "splines.h"
#include "utilities.h"
#include "dataframe_list.h"
#include "thread_utils.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <functional>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <vector>

// Helper: estimate psi using counterfactual residuals and logistic regression
ListCpp est_psi_tsegest(
    int n2, int q, int p2, int nids2, 
    const std::vector<int>& idx2, const std::vector<int>& stratum3,
    const std::vector<int>& os3, const std::vector<double>& os_time3,
    const std::vector<double>& censor_time3,
    const std::vector<int>& swtrt3, const std::vector<double>& swtrt_time3,
    const std::vector<int>& id2, const std::vector<int>& cross2,
    const std::vector<double>& tstart2, const std::vector<double>& tstop2,
    const std::vector<std::string>& covariates_lgs,
    const FlatMatrix& z_lgs2,
    int ns_df, const FlatMatrix& s, bool firth, bool flic,
    bool recensor, double alpha,
    const std::string& ties, double offset, double x) {
  
  double a = std::exp(x);
  double c0 = std::min(1.0, a);
  
  // Counterfactual times and events at subject-level (nids2)
  std::vector<double> t_star(nids2);
  std::vector<int> d_star(nids2);
  
  for (int i = 0; i < nids2; ++i) {
    double u_star;
    if (swtrt3[i] == 1) {
      double b2 = swtrt_time3[i] - offset;
      u_star = b2 + (os_time3[i] - b2) * a;
    } else {
      u_star = os_time3[i];
    }
    
    if (recensor) {
      double c_star = censor_time3[i] * c0;
      t_star[i] = std::min(u_star, c_star);
      d_star[i] = (c_star < u_star) ? 0 : os3[i];
    } else {
      t_star[i] = u_star;
      d_star[i] = os3[i];
    }
  }
  
  // Build DataFrameCpp for null Cox (subject-level)
  DataFrameCpp dn;
  dn.push_back(t_star, "t_star");
  dn.push_back(d_star, "d_star");
  dn.push_back(stratum3, "ustratum");
  
  std::vector<double> init(1, NaN);
  ListCpp fn = phregcpp(
    dn, {"ustratum"}, "t_star", "", "d_star", 
    {""}, "", "", "", ties, init, 0, 0, 1, 0, 0, alpha);
  
  // Extract residuals (subject-level) and expand to observation-level
  std::vector<double> resid3 = fn.get<std::vector<double>>("residuals");
  std::vector<double> resid(n2);
  for (int i = 0; i < nids2; ++i) {
    int start = idx2[i], end = idx2[i+1];
    std::fill(resid.begin() + start, resid.begin() + end, resid3[i]);
  }
  
  // Build logistic regression dataset (observation-level)
  DataFrameCpp dl;
  dl.push_back(id2, "uid");
  dl.push_back(tstart2, "tstart");
  dl.push_back(tstop2, "tstop");
  dl.push_back(cross2, "cross");
  dl.push_back(resid, "counterfactual");
  
  // Append covariates (q + p2) from z_lgs2 (FlatMatrix gives column-major)
  for (int j = 0; j < q + p2; ++j) {
    std::vector<double> col = flatmatrix_get_column(z_lgs2, j);
    dl.push_back(std::move(col), covariates_lgs[j + 1]);
  }
  // Append spline columns (ns_df)
  for (int j = 0; j < ns_df; ++j) {
    std::vector<double> col = flatmatrix_get_column(s, j);
    dl.push_back(std::move(col), covariates_lgs[q + p2 + j + 1]);
  }
  
  // Fit logistic regression
  ListCpp fl = logisregcpp(
    dl, "cross", covariates_lgs, "", "", "", 
    "uid", "logit", init, 1, firth, flic, 0, alpha);
  
  DataFrameCpp sumstat = fl.get<DataFrameCpp>("sumstat");
  bool fail = sumstat.get<unsigned char>("fail")[0];
  
  DataFrameCpp parest = fl.get<DataFrameCpp>("parest");
  std::vector<double> z = parest.get<double>("z");
  double z_counterfactual = z[1];
  
  ListCpp out;
  out.push_back(std::move(dn), "data_nullcox");
  out.push_back(std::move(fn), "fit_nullcox");
  out.push_back(std::move(dl), "data_logis");
  out.push_back(std::move(fl), "fit_logis");
  out.push_back(z_counterfactual, "z_counterfactual");
  out.push_back(fail, "fail");
  return out;
};


// [[Rcpp::export]]
Rcpp::List tsegestcpp(const Rcpp::DataFrame& df,
                      const std::string& id,
                      const std::vector<std::string>& stratum,
                      const std::string& tstart,
                      const std::string& tstop,
                      const std::string& event,
                      const std::string& treat,
                      const std::string& censor_time,
                      const std::string& pd,
                      const std::string& pd_time,
                      const std::string& swtrt,
                      const std::string& swtrt_time,
                      const std::vector<std::string>& base_cov,
                      const std::vector<std::string>& conf_cov,
                      bool strata_main_effect_only,
                      int ns_df,
                      bool firth,
                      bool flic,
                      double low_psi,
                      double hi_psi,
                      int n_eval_z,
                      bool recensor,
                      bool admin_recensor_only,
                      bool swtrt_control_only,
                      bool gridsearch,
                      const std::string& root_finding,
                      double alpha,
                      const std::string& ties,
                      double tol,
                      double offset,
                      bool boot,
                      int n_boot,
                      int seed) {
  
  DataFrameCpp data = convertRDataFrameToCpp(df);  

  int n = static_cast<int>(data.nrows());
  int p = static_cast<int>(base_cov.size());
  if (p == 1 && base_cov[0] == "") p = 0;
  
  
  int p2 = static_cast<int>(conf_cov.size());
  if (p2 == 1 && conf_cov[0] == "") p2 = 0;
  
  // process stratification variables
  int p_stratum = static_cast<int>(stratum.size());
  bool has_stratum = false;
  std::vector<int> stratumn(n);
  DataFrameCpp u_stratum;
  std::vector<int> d(p_stratum);
  IntMatrix stratan(n, p_stratum);
  ListCpp levels;
  if (!(p_stratum == 0 || (p_stratum == 1 && stratum[0] == ""))) {
    ListCpp out = bygroup(data, stratum);
    has_stratum = true;
    stratumn = out.get<std::vector<int>>("index");
    u_stratum = out.get<DataFrameCpp>("lookup");
    d = out.get<std::vector<int>>("nlevels");
    stratan = out.get<IntMatrix>("indices");
    levels = out.get_list("lookups_per_variable");
  }
  std::vector<int> stratumn_unique = unique_sorted(stratumn);
  int nstrata = static_cast<int>(stratumn_unique.size());
  
  // create the numeric id variable
  if (id.empty() || !data.containElementNamed(id))
    throw std::invalid_argument("data must contain the id variable");
  std::vector<int> idn(n);
  std::vector<int> idwi;
  std::vector<double> idwn;
  std::vector<std::string> idwc;
  if (data.int_cols.count(id)) {
    auto v = data.get<int>(id);
    idwi = unique_sorted(v);
    idn = matchcpp(v, idwi);
  } else if (data.numeric_cols.count(id)) {
    auto v = data.get<double>(id);
    idwn = unique_sorted(v);
    idn = matchcpp(v, idwn);
  } else if (data.string_cols.count(id)) {
    auto v = data.get<std::string>(id);
    idwc = unique_sorted(v);
    idn = matchcpp(v, idwc);
  } else throw std::invalid_argument(
      "incorrect type for the id variable in data");
  
  // --- tstart / tstop existence and checks ---
  if (tstart.empty() || !data.containElementNamed(tstart))
    throw std::invalid_argument("data must contain the tstart variable");
  std::vector<double> tstartn(n);
  if (data.int_cols.count(tstart)) {
    const std::vector<int>& vi = data.get<int>(tstart);
    for (int i = 0; i < n; ++i) tstartn[i] = static_cast<double>(vi[i]);
  } else if (data.numeric_cols.count(tstart)) {
    tstartn = data.get<double>(tstart);
  } else {
    throw std::invalid_argument("tstart variable must be integer or numeric");
  }
  for (int i = 0; i < n; ++i) {
    if (!std::isnan(tstartn[i]) && tstartn[i] < 0.0)
      throw std::invalid_argument("tstart must be nonnegative");
  }
  
  if (tstop.empty() || !data.containElementNamed(tstop))
    throw std::invalid_argument("data must contain the tstop variable");
  std::vector<double> tstopn(n);
  if (data.int_cols.count(tstop)) {
    const std::vector<int>& vi = data.get<int>(tstop);
    for (int i = 0; i < n; ++i) tstopn[i] = static_cast<double>(vi[i]);
  } else if (data.numeric_cols.count(tstop)) {
    tstopn = data.get<double>(tstop);
  } else {
    throw std::invalid_argument("tstop variable must be integer or numeric");
  }
  for (int i = 0; i < n; ++i) {
    if (!std::isnan(tstartn[i]) && !std::isnan(tstopn[i]) && tstopn[i] <= tstartn[i])
      throw std::invalid_argument("tstop must be greater than tstart");
  }
  
  // --- event variable ---
  if (event.empty() || !data.containElementNamed(event)) {
    throw std::invalid_argument("data must contain the event variable");
  }
  std::vector<int> eventn(n);
  if (data.bool_cols.count(event)) {
    const std::vector<unsigned char>& vb = data.get<unsigned char>(event);
    for (int i = 0; i < n; ++i) eventn[i] = vb[i] ? 1 : 0;
  } else if (data.int_cols.count(event)) {
    eventn = data.get<int>(event);
  } else if (data.numeric_cols.count(event)) {
    const std::vector<double>& vd = data.get<double>(event);
    for (int i = 0; i < n; ++i) eventn[i] = static_cast<int>(vd[i]);
  } else {
    throw std::invalid_argument("event variable must be bool, integer or numeric");
  }
  for (double val : eventn) if (val != 0 && val != 1)
    throw std::invalid_argument("event must be 1 or 0 for each observation");
  if (std::all_of(eventn.begin(), eventn.end(), [](int x){ return x == 0; })) {
    throw std::invalid_argument("at least 1 event is needed");
  }
  
  // create the numeric treat variable
  if (treat.empty() || !data.containElementNamed(treat))
    throw std::invalid_argument("data must contain the treat variable");
  std::vector<int> treatn(n);
  std::vector<int> treatwi;
  std::vector<double> treatwn;
  std::vector<std::string> treatwc;
  if (data.bool_cols.count(treat) || data.int_cols.count(treat)) {
    std::vector<int> treatv(n);
    if (data.bool_cols.count(treat)) {
      const std::vector<unsigned char>& treatvb = data.get<unsigned char>(treat);
      for (int i = 0; i < n; ++i) treatv[i] = treatvb[i] ? 1 : 0;
    } else treatv = data.get<int>(treat);
    treatwi = unique_sorted(treatv); // obtain unique treatment values
    if (treatwi.size() != 2)
      throw std::invalid_argument("treat must have two and only two distinct values");
    if (std::all_of(treatwi.begin(), treatwi.end(), [](int v) {
      return v == 0 || v == 1; })) {
      treatwi = {1, 0}; // special handling for 1/0 treatment coding
      for (int i = 0; i < n; ++i) treatn[i] = 2 - treatv[i];
    } else {
      treatn = matchcpp(treatv, treatwi, 1);
    }
  } else if (data.numeric_cols.count(treat)) {
    const std::vector<double>& treatv = data.get<double>(treat);
    treatwn = unique_sorted(treatv);
    if (treatwn.size() != 2)
      throw std::invalid_argument("treat must have two and only two distinct values");
    if (std::all_of(treatwn.begin(), treatwn.end(), [](double v) {
      return v == 0.0 || v == 1.0; })) {
      treatwn = {1.0, 0.0};
      for (int i = 0; i < n; ++i) treatn[i] = 2 - static_cast<int>(treatv[i]);
    } else {
      treatn = matchcpp(treatv, treatwn, 1);
    }
  } else if (data.string_cols.count(treat)) {
    const std::vector<std::string>& treatv = data.get<std::string>(treat);
    treatwc = unique_sorted(treatv);
    if (treatwc.size() != 2)
      throw std::invalid_argument("treat must have two and only two distinct values");
    treatn = matchcpp(treatv, treatwc, 1);
  } else {
    throw std::invalid_argument(
        "incorrect type for the treat variable in the input data");
  }
  for (int i = 0; i < n; ++i) {
    treatn[i] = 2 - treatn[i]; // convert to 1/0 coding
  }
  
  // --- censor_time variable ---
  if (censor_time.empty() || !data.containElementNamed(censor_time))
    throw std::invalid_argument("data must contain the censor_time variable");
  std::vector<double> censor_timen(n);
  if (data.int_cols.count(censor_time)) {
    const std::vector<int>& vi = data.get<int>(censor_time);
    for (int i = 0; i < n; ++i) censor_timen[i] = static_cast<double>(vi[i]);
  } else if (data.numeric_cols.count(censor_time)) {
    censor_timen = data.get<double>(censor_time);
  } else {
    throw std::invalid_argument("censor_time variable must be integer or numeric");
  }
  for (double v : censor_timen) if (v < 0.0 || std::isnan(v))
    throw std::invalid_argument("censor_time cannot be missing");
  for (int i = 0; i < n; ++i) {
    if (censor_timen[i] < tstopn[i]) throw std::invalid_argument(
        "censor_time must be greater than or equal to tstop");
  }
  
  // --- pd variable ---
  if (pd.empty() || !data.containElementNamed(pd)) {
    throw std::invalid_argument("data must contain the pd variable");
  }
  std::vector<int> pdn(n);
  if (data.bool_cols.count(pd)) {
    const std::vector<unsigned char>& vb = data.get<unsigned char>(pd);
    for (int i = 0; i < n; ++i) pdn[i] = vb[i] ? 1 : 0;
  } else if (data.int_cols.count(pd)) {
    pdn = data.get<int>(pd);
  } else if (data.numeric_cols.count(pd)) {
    const std::vector<double>& vd = data.get<double>(pd);
    for (int i = 0; i < n; ++i) pdn[i] = static_cast<int>(vd[i]);
  } else {
    throw std::invalid_argument("pd variable must be bool, integer or numeric");
  }
  for (double val : eventn) if (val != 0 && val != 1)
    throw std::invalid_argument("event must be 1 or 0 for each observation");
  
  // --- pd_time variable ---
  if (pd_time.empty() || !data.containElementNamed(pd_time))
    throw std::invalid_argument("data must contain the pd_time variable");
  std::vector<double> pd_timen(n);
  if (data.int_cols.count(pd_time)) {
    const std::vector<int>& vi = data.get<int>(pd_time);
    for (int i = 0; i < n; ++i) pd_timen[i] = static_cast<double>(vi[i]);
  } else if (data.numeric_cols.count(pd_time)) {
    pd_timen = data.get<double>(pd_time);
  } else {
    throw std::invalid_argument("pd_time variable must be integer or numeric");
  }
  
  // check consistency between pd and pd_time
  for (int i = 0; i < n; ++i) {
    if (pdn[i] == 1 && std::isnan(pd_timen[i])) {
      throw std::runtime_error("pd_time must not be missing when pd = 1");
    }
    if (pdn[i] == 1 && pd_timen[i] < 0.0) {
      throw std::runtime_error("pd_time must be nonnegative when pd = 1");
    }
  }
  
  // --- swtrt variable ---
  if (swtrt.empty() || !data.containElementNamed(swtrt)) {
    throw std::invalid_argument("data must contain the swtrt variable");
  }
  std::vector<int> swtrtn(n);
  if (data.bool_cols.count(swtrt)) {
    const std::vector<unsigned char>& vb = data.get<unsigned char>(swtrt);
    for (int i = 0; i < n; ++i) swtrtn[i] = vb[i] ? 1 : 0;
  } else if (data.int_cols.count(swtrt)) {
    swtrtn = data.get<int>(swtrt);
  } else if (data.numeric_cols.count(swtrt)) {
    const std::vector<double>& vd = data.get<double>(swtrt);
    for (int i = 0; i < n; ++i) swtrtn[i] = static_cast<int>(vd[i]);
  } else {
    throw std::invalid_argument("swtrt variable must be bool, integer or numeric");
  }
  for (double val : swtrtn) if (val != 0 && val != 1)
    throw std::invalid_argument("swtrt must be 1 or 0 for each observation");
  
  // check presence of at least 1 pd and swtrt in each group
  bool found_control = false;
  for (int i = 0; i < n; ++i) {
    if (pdn[i] == 1 && swtrtn[i] == 1 && treatn[i] == 0) {
      found_control = true;
      break;
    }
  }
  if (!found_control) {
    throw std::runtime_error(
        "at least 1 pd and swtrt is needed in the control group");
  }
  
  if (!swtrt_control_only)  {
    bool found_treated = false;
    for (int i = 0; i < n; ++i) {
      if (pdn[i] == 1 && swtrtn[i] == 1 && treatn[i] == 1) {
        found_treated = true;
        break;
      }
    }
    if (!found_treated) {
      throw std::runtime_error(
          "at least 1 pd and swtrt is needed in the treatment group");
    }
  }

  // --- swtrt_time variable ---
  if (swtrt_time.empty() || !data.containElementNamed(swtrt_time))
    throw std::invalid_argument("data must contain the swtrt_time variable");
  std::vector<double> swtrt_timen(n);
  if (data.int_cols.count(swtrt_time)) {
    const std::vector<int>& vi = data.get<int>(swtrt_time);
    for (int i = 0; i < n; ++i) swtrt_timen[i] = static_cast<double>(vi[i]);
  } else if (data.numeric_cols.count(swtrt_time)) {
    swtrt_timen = data.get<double>(swtrt_time);
  } else {
    throw std::invalid_argument("swtrt_time variable must be integer or numeric");
  }
  
  // check consistency between swtrt and swtrt_time
  for (int i = 0; i < n; ++i) {
    if (swtrtn[i] == 1 && std::isnan(swtrt_timen[i])) {
      throw std::runtime_error("swtrt_time must not be missing when swtrt=1");
    }
    if (swtrtn[i] == 1 && swtrt_timen[i] < 0.0) {
      throw std::runtime_error("swtrt_time must be nonnegative when swtrt=1");
    }
  }
      
  // adjust pd and pd_time based on swtrt and swtrt_time:
  // if the patient switched before pd, set pd time equal to switch time
  for (int i = 0; i < n; ++i) {
    if (pdn[i] == 1 && swtrtn[i] == 1 && swtrt_timen[i] < pd_timen[i]) {
      pd_timen[i] = swtrt_timen[i];
    }
    
    if (pdn[i] == 0 && swtrtn[i] == 1) {
      pdn[i] = 1; 
      pd_timen[i] = swtrt_timen[i];
    }
  }
  
  // make sure offset is less than or equal to observed time variables
  for (int i = 0; i < n; ++i) {
    if (pdn[i] == 1 && pd_timen[i] < offset) {
      throw std::runtime_error("pd_time must be great than or equal to offset");
    }
    if (swtrtn[i] == 1 && swtrt_timen[i] < offset) {
      throw std::runtime_error("swtrt_time must be great than or equal to offset");
    }
  }
  
  // covariates for the Cox model containing treat and base_cov
  std::vector<std::string> covariates(p + 1);
  FlatMatrix zn(n, p);
  covariates[0] = "treated";
  for (int j = 0; j < p; ++j) {
    const std::string& zj = base_cov[j];
    if (!data.containElementNamed(zj))
      throw std::invalid_argument("data must contain the variables in base_cov");
    if (zj == treat)
      throw std::invalid_argument("treat should be excluded from base_cov");
    covariates[j + 1] = zj;
    double* zn_col = zn.data_ptr() + j * n;
    if (data.bool_cols.count(zj)) {
      const std::vector<unsigned char>& vb = data.get<unsigned char>(zj);
      for (int i = 0; i < n; ++i) zn_col[i] = vb[i] ? 1.0 : 0.0;
    } else if (data.int_cols.count(zj)) {
      const std::vector<int>& vi = data.get<int>(zj);
      for (int i = 0; i < n; ++i) zn_col[i] = static_cast<double>(vi[i]);
    } else if (data.numeric_cols.count(zj)) {
      const std::vector<double>& vd = data.get<double>(zj);
      std::memcpy(zn_col, vd.data(), n * sizeof(double));
    } else {
      throw std::invalid_argument("covariates must be bool, integer or numeric");
    }
  }
  
  // number of columns corresponding to the strata effects
  int q = 0;
  if (has_stratum) {
    if (strata_main_effect_only) {
      q = 0;
      for (int i = 0; i < p_stratum; ++i) q += d[i] - 1;
    } else {
      q = nstrata - 1;
    }
  }
  
  // covariates for the pooled logistic regression model for control with pd
  // including counterfactual, stratum, conf_cov, and ns_df spline terms
  std::vector<std::string> covariates_lgs(q + p2 + ns_df + 1);
  FlatMatrix z_lgsn(n, q + p2);
  covariates_lgs[0] = "counterfactual";
  if (has_stratum) {
    if (strata_main_effect_only) {
      int k = 0;
      for (int i = 0; i < p_stratum; ++i) {
        const std::string& s = stratum[i];
        int di = d[i] - 1;
        
        if (u_stratum.string_cols.count(s)) {
          auto u = levels.get<std::vector<std::string>>(s);
          for (int j = 0; j < di; ++j) {
            covariates_lgs[k + j + 1] = s + sanitize(u[j]);
          }
        } else if (u_stratum.numeric_cols.count(s)) {
          auto u = levels.get<std::vector<double>>(s);
          for (int j = 0; j < di; ++j) {
            covariates_lgs[k + j + 1] = s + std::to_string(u[j]);
          }
        } else if (u_stratum.int_cols.count(s)) {
          auto u = levels.get<std::vector<int>>(s);
          for (int j = 0; j < di; ++j) {
            covariates_lgs[k + j + 1] = s + std::to_string(u[j]);
          }
        } else if (u_stratum.bool_cols.count(s)) {
          auto u = levels.get<std::vector<unsigned char>>(s);
          for (int j = 0; j < di; ++j) {
            covariates_lgs[k + j + 1] = s + std::to_string(u[j]);
          }
        }
        
        for (int j = 0; j < di; ++j) {
          const int* stratan_col = stratan.data_ptr() + i * n;
          double* z_lgsn_col = z_lgsn.data_ptr() + (k + j) * n;
          for (int r = 0; r < n; ++r) {
            z_lgsn_col[r] = stratan_col[r] == j ? 1.0 : 0.0;
          }
        }
        
        k += di;
      }
    } else {
      for (int j = 0; j < nstrata - 1; ++j) {
        // locate the first observation in the stratum
        int first_k = 0;
        for (; first_k<n; ++first_k) {
          if (stratumn[first_k] == j) break;
        }
        
        covariates_lgs[j + 1] = "";
        
        for (int i = 0; i < p_stratum; ++i) {
          const std::string& s = stratum[i];
          
          std::vector<int> q_col = intmatrix_get_column(stratan, i);
          int l = q_col[first_k];
          
          if (u_stratum.string_cols.count(s)) {
            auto u = levels.get<std::vector<std::string>>(s);
            covariates_lgs[j + 1] += s + sanitize(u[l]);
          } else if (u_stratum.numeric_cols.count(s)) {
            auto u = levels.get<std::vector<double>>(s);
            covariates_lgs[j + 1] += s + std::to_string(u[l]);
          } else if (u_stratum.int_cols.count(s)) {
            auto u = levels.get<std::vector<int>>(s);
            covariates_lgs[j + 1] += s + std::to_string(u[l]);
          } else if (u_stratum.bool_cols.count(s)) {
            auto u = levels.get<std::vector<unsigned char>>(s);
            covariates_lgs[j + 1] += s + std::to_string(u[l]);
          }
          
          if (i < p_stratum - 1) {
            covariates_lgs[j + 1] += ".";
          }
        }
        
        double* z_lgsn_col = z_lgsn.data_ptr() + j * n;
        for (int r = 0; r < n; ++r) {
          z_lgsn_col[r] = stratumn[r] == j ? 1.0 : 0.0;
        }
      }
    }
  }
  
  // append confounding covariates
  for (int j = 0; j < p2; ++j) {
    const std::string& zj = conf_cov[j];
    if (!data.containElementNamed(zj))
      throw std::invalid_argument("data must contain the variables in conf_cov");
    if (zj == treat)
      throw std::invalid_argument("treat should be excluded from conf_cov");
    covariates_lgs[q + j + 1] = zj;
    double* z_lgsn_col = z_lgsn.data_ptr() + (q + j) * n;
    if (data.bool_cols.count(zj)) {
      const std::vector<unsigned char>& vb = data.get<unsigned char>(zj);
      for (int i = 0; i < n; ++i)
        z_lgsn_col[i] = vb[i] ? 1.0 : 0.0;
    } else if (data.int_cols.count(zj)) {
      const std::vector<int>& vi = data.get<int>(zj);
      for (int i = 0; i < n; ++i)
        z_lgsn_col[i] = static_cast<double>(vi[i]);
    } else if (data.numeric_cols.count(zj)) {
      const std::vector<double>& vd = data.get<double>(zj);
      std::memcpy(z_lgsn_col, vd.data(), n * sizeof(double));
    } else {
      throw std::invalid_argument("covariates must be bool, integer or numeric");
    }
  }
  
  if (ns_df < 0) {
    throw std::invalid_argument("ns_df must be a nonnegative integer");
  }
  
  for (int j = 0; j < ns_df; ++j) {
    covariates_lgs[q + p2 + j + 1] = "ns" + std::to_string(j+1);
  }
  
  if (low_psi >= hi_psi) {
    throw std::invalid_argument("low_psi must be less than hi_psi");
  }
  
  if (n_eval_z < 2) {
    throw std::invalid_argument("n_eval_z must be greater than or equal to 2");
  }
  
  std::string rooting = root_finding;
  std::for_each(rooting.begin(), rooting.end(), [](char & c) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });
  if (rooting == "uniroot" || rooting.find("br", 0) == 0) rooting = "brent";
  else if (rooting.find("bi", 0) == 0) rooting = "bisection";
  if (!(rooting == "brent" || rooting == "bisection"))
    throw std::invalid_argument("root_finding must be brent or bisection");
  
  if (alpha <= 0.0 || alpha >= 0.5)
    throw std::invalid_argument("alpha must lie between 0 and 0.5");
  if (ties != "efron" && ties != "breslow")
    throw std::invalid_argument("ties must be efron or breslow");
  if (tol <= 0.0) throw std::invalid_argument("tol must be positive");
  if (offset < 0.0) throw std::invalid_argument("offset must be nonnegative");
  if (n_boot < 100)
    throw std::invalid_argument("n_boot must be greater than or equal to 100");

  // exclude observations with missing values
  std::vector<unsigned char> sub(n,1);
  for (int i = 0; i < n; ++i) {
    if (idn[i] == INT_MIN || stratumn[i] == INT_MIN ||
        std::isnan(tstartn[i]) || std::isnan(tstopn[i]) ||
        eventn[i] == INT_MIN || treatn[i] == INT_MIN ||
        std::isnan(censor_timen[i]) || pdn[i] == INT_MIN ||
        swtrtn[i] == INT_MIN) {
      sub[i] = 0; continue;
    }
    for (int j = 0; j < p; ++j) {
      if (std::isnan(zn(i,j))) { sub[i] = 0; break; }
    }
  }
  
  std::vector<int> keep = which(sub);
  if (keep.empty())
    throw std::invalid_argument("no observations without missing values");
  subset_in_place(idn, keep);
  subset_in_place(stratumn, keep);
  subset_in_place(tstartn, keep);
  subset_in_place(tstopn, keep);
  subset_in_place(eventn, keep);
  subset_in_place(treatn, keep);
  subset_in_place(censor_timen, keep);
  subset_in_place(pdn, keep);
  subset_in_place(pd_timen, keep);
  subset_in_place(swtrtn, keep);
  subset_in_place(swtrt_timen, keep);
  subset_in_place_flatmatrix(zn, keep);
  subset_in_place_flatmatrix(z_lgsn, keep);
  n = static_cast<int>(keep.size());
  
  // split at treatment switching into two observations if treatment
  // switching occurs strictly between tstart and tstop for a subject
  std::vector<unsigned char> tosplit(n);
  for (int i = 0; i < n; ++i) {
    tosplit[i] = swtrtn[i] == 1 && swtrt_timen[i] > tstartn[i] &&
      swtrt_timen[i] < tstopn[i] ? 1 : 0;
  }
  
  int k = std::accumulate(tosplit.begin(), tosplit.end(), 0);
  if (k > 0) {
    std::vector<int> sub = which(tosplit);
    for (int i = 0; i < k; ++i) {
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
      
      // change tstop and event for the old observation
      tstopn[l] = swtrt_timen[l];
      eventn[l] = 0;
    }
    
    // append new rows to the covariate matrices
    FlatMatrix zn_new = subset_flatmatrix(zn, sub);
    FlatMatrix z_lgsn_new = subset_flatmatrix(z_lgsn, sub);
    append_flatmatrix(zn, zn_new);
    append_flatmatrix(z_lgsn, z_lgsn_new);
    
    // update number of rows
    n = n + k;
  }
  
  // sort data by treatment group, id, and time
  std::vector<int> order = seqcpp(0, n - 1);
  std::sort(order.begin(), order.end(), [&](int i, int j) {
    return std::tie(treatn[i], idn[i], tstopn[i]) <
      std::tie(treatn[j], idn[j], tstopn[j]);
  });
  
  subset_in_place(idn, order);
  subset_in_place(stratumn, order);
  subset_in_place(tstartn, order);
  subset_in_place(tstopn, order);
  subset_in_place(eventn, order);
  subset_in_place(treatn, order);
  subset_in_place(censor_timen, order);
  subset_in_place(pdn, order);
  subset_in_place(pd_timen, order);
  subset_in_place(swtrtn, order);
  subset_in_place(swtrt_timen, order);
  subset_in_place_flatmatrix(zn, order);
  subset_in_place_flatmatrix(z_lgsn, order);
  
  // identify first and last observation within an id
  std::vector<int> idx(1,0); // first observation within an id
  for (int i = 1; i < n; ++i) {
    if (idn[i] != idn[i-1]) {
      idx.push_back(i);
    }
  }
  int nids = static_cast<int>(idx.size());
  idx.push_back(n);
  
  std::vector<int> idx1(nids); // last observation within an id
  for (int i = 0; i < nids; ++i) {
    idx1[i] = idx[i+1] - 1;
  }
  
  // create os_time variable
  std::vector<int> osn(n);
  std::vector<double> os_timen(n);
  for (int i = 0; i < nids; ++i) {
    int k = idx1[i];
    int ev = eventn[k];
    double ts = tstopn[k];
    int start = idx[i], end = idx[i+1];
    std::fill(osn.begin() + start, osn.begin() + end, ev);
    std::fill(os_timen.begin() + start, os_timen.begin() + end, ts);
  }
  
  for (int i = 0; i < n; ++i) {
    if (pdn[i] == 1 && pd_timen[i] > os_timen[i]) {
      throw std::invalid_argument("pd_time must be less than or equal to os_time");
    }
    if (swtrtn[i] == 1 && swtrt_timen[i] > os_timen[i]) {
      throw std::invalid_argument("swtrt_time must be less than or equal to os_time");
    }
  }
  
  // adjust censor_time for dropouts if recensoring is requested
  if (!admin_recensor_only) { // use the actual censoring time for dropouts
    for (int i = 0; i < nids; ++i) {
      int k = idx1[i];
      int ev = eventn[k];
      double ts = tstopn[k];
      if (ev == 0) {
        int start = idx[i], end = idx[i+1];
        std::fill(censor_timen.begin() + start, censor_timen.begin() + end, ts);
      }
    }
  }
  
  // subset to one observation per id for event summary
  std::vector<int> treatn1 = subset(treatn, idx1);
  std::vector<int> eventn1 = subset(eventn, idx1);
  std::vector<int> pdn1 = subset(pdn, idx1);
  std::vector<int> swtrtn1 = subset(swtrtn, idx1);
  
  // summarize number of deaths and switches by treatment arm
  std::vector<int> treat_out = {0, 1};
  std::vector<double> n_total(2);
  std::vector<double> n_event(2);
  std::vector<double> n_pd(2);
  std::vector<double> n_switch(2);
  for (int i = 0; i < nids; ++i) {
    int g = treatn1[i];
    ++n_total[g];
    if (eventn1[i] == 1) ++n_event[g];
    if (pdn1[i] == 1) ++n_pd[g];
    if (swtrtn1[i] == 1) ++n_switch[g];
  }
  
  // Compute percentages
  std::vector<double> pct_event(2);
  std::vector<double> pct_pd(2);
  std::vector<double> pct_switch(2);
  for (int g = 0; g < 2; g++) {
    pct_event[g] = 100.0 * n_event[g] / n_total[g];
    pct_pd[g] = 100.0 * n_pd[g] / n_total[g];
    pct_switch[g] = 100.0 * n_switch[g] / n_total[g];
  }
  
  // Combine count and percentage
  DataFrameCpp event_summary;
  event_summary.push_back(std::move(treat_out), "treated");
  event_summary.push_back(n_total, "n");
  event_summary.push_back(std::move(n_event), "event_n");
  event_summary.push_back(std::move(pct_event), "event_pct");
  event_summary.push_back(std::move(n_pd), "pd_n");
  event_summary.push_back(std::move(pct_pd), "pd_pct");
  event_summary.push_back(std::move(n_switch), "switch_n");
  event_summary.push_back(std::move(pct_switch), "switch_pct");
  
  double zcrit = boost_qnorm(1.0 - alpha / 2.0);
  
  // evaluation points for psi
  double step_psi = (hi_psi - low_psi) / (n_eval_z - 1);
  std::vector<double> psi(n_eval_z);
  for (int i = 0; i < n_eval_z; ++i) {
    psi[i] = low_psi + i * step_psi;
  }
  
  auto f = [data, has_stratum, stratum, p_stratum, u_stratum, 
            treat, treatwi, treatwn, treatwc, id, idwi, idwn, idwc,
            p, p2, q, covariates, covariates_lgs, ns_df, firth, flic, 
            low_psi, hi_psi, n_eval_z, psi, recensor, swtrt_control_only, 
            gridsearch, rooting, alpha, zcrit, ties, tol, offset] (
                const std::vector<int>& idb,
                const std::vector<int>& stratumb,
                const std::vector<double>& tstartb,
                const std::vector<double>& tstopb,
                const std::vector<int>& eventb,
                const std::vector<int>& treatb,
                const std::vector<int>& osb,
                const std::vector<double>& os_timeb,
                const std::vector<double>& censor_timeb,
                const std::vector<int>& pdb,
                const std::vector<double>& pd_timeb,
                const std::vector<int>& swtrtb,
                const std::vector<double>& swtrt_timeb,
                const FlatMatrix& zb,
                const FlatMatrix& z_lgsb, int k) -> ListCpp {
                  int n = static_cast<int>(idb.size());
                  bool fail = false; // whether any model fails to converge
                  std::vector<double> init(1, NaN);
                  
                  std::vector<int> idx(1, 0); // first observation within an id
                  for (int i = 1; i < n; ++i) {
                    if (idb[i] != idb[i-1]) {
                      idx.push_back(i);
                    }
                  }
                  int nids = static_cast<int>(idx.size());
                  idx.push_back(n);
                  
                  std::vector<int> idx1(nids); // last observation within an id
                  for (int i = 0; i < nids; ++i) {
                    idx1[i] = idx[i+1] - 1;
                  }
                  
                  // one observation per subject
                  std::vector<int> id1 = subset(idb, idx1);
                  std::vector<int> stratum1 = subset(stratumb, idx1);
                  std::vector<double> time1 = subset(tstopb, idx1);
                  std::vector<int> event1 = subset(eventb, idx1);
                  std::vector<int> treat1 = subset(treatb, idx1);
                  std::vector<double> censor_time1 = subset(censor_timeb, idx1);
                  std::vector<int> swtrt1 = subset(swtrtb, idx1);
                  std::vector<double> swtrt_time1 = subset(swtrt_timeb, idx1);
                  FlatMatrix z1 = subset_flatmatrix(zb, idx1);

                  // time and event adjusted for treatment switching
                  std::vector<double> t_star = time1;
                  std::vector<int> d_star = event1;
                  
                  double psi0hat = NaN, psi0lower = NaN, psi0upper = NaN;
                  double psi1hat = NaN, psi1lower = NaN, psi1upper = NaN;
                  std::vector<double> psi0hat_vec, psi1hat_vec;
                  std::string psi_CI_type;
                  
                  // initialize data_switch, km_switch, eval_z, 
                  // data_nullcox, fit_nullcox, data_logis, fit_logis
                  std::vector<ListPtr> data_switch(2), km_switch(2), eval_z(2);
                  std::vector<ListPtr> data_nullcox(2), fit_nullcox(2);
                  std::vector<ListPtr> data_logis(2), fit_logis(2);
                  if (k == -1) {
                    DataFrameCpp nulldata;
                    ListCpp nullfit;
                    auto make_group = [&](int h) -> std::array<ListPtr, 7> {
                      ListPtr ds = std::make_shared<ListCpp>();
                      ListPtr km = std::make_shared<ListCpp>();
                      ListPtr ez = std::make_shared<ListCpp>();
                      ListPtr dn = std::make_shared<ListCpp>();
                      ListPtr dl = std::make_shared<ListCpp>();
                      ListPtr fn = std::make_shared<ListCpp>();
                      ListPtr fl = std::make_shared<ListCpp>();
                      
                      // push the common placeholders
                      ds->push_back(nulldata, "data");
                      km->push_back(nulldata, "data");
                      ez->push_back(nulldata, "data");
                      dn->push_back(nulldata, "data");
                      dl->push_back(nulldata, "data");
                      fn->push_back(nullfit, "fit");
                      fl->push_back(nullfit, "fit");
                      
                      // push the treatment value to each group element
                      if (data.bool_cols.count(treat) || 
                          data.int_cols.count(treat)) {
                        int v = treatwi[1 - h];
                        for (auto ptr : {ds, km, ez, dn, dl, fn, fl}) 
                          ptr->push_back(v, treat);
                      } else if (data.numeric_cols.count(treat)) {
                        double v = treatwn[1 - h];
                        for (auto ptr : {ds, km, ez, dn, dl, fn, fl}) 
                          ptr->push_back(v, treat);
                      } else if (data.string_cols.count(treat)) {
                        std::string v = treatwc[1 - h];
                        for (auto ptr : {ds, km, ez, dn, dl, fn, fl}) 
                          ptr->push_back(v, treat);
                      }
                      
                      return {ds, km, ez, dn, dl, fn, fl};
                    };
                    
                    for (int h = 0; h < 2; ++h) {
                      auto grp = make_group(h);
                      data_switch[h]  = grp[0];
                      km_switch[h]    = grp[1];
                      eval_z[h]       = grp[2];
                      data_nullcox[h] = grp[3];
                      data_logis[h]   = grp[4];
                      fit_nullcox[h]  = grp[5];
                      fit_logis[h]    = grp[6];
                    }
                  }
                  
                  DataFrameCpp data_outcome, km_outcome, lr_outcome;
                  ListCpp fit_outcome;
                  double hrhat = NaN, hrlower = NaN, hrupper = NaN, pvalue = NaN;
                  
                  bool psimissing = false;
                  
                  // # arms that include patients who switched treatment
                  int K = swtrt_control_only ? 1 : 2;
                  for (int h = 0; h < K; ++h) {
                    // post progression data up to switching for the treat
                    std::vector<int> l;
                    l.reserve(n);
                    for (int i = 0; i < n; ++i) {
                      if (treatb[i] != h) continue;
                      if (pdb[i] != 1 || tstopb[i] < pd_timeb[i]) continue;
                      if (swtrtb[i] == 1) {
                        if (tstartb[i] < swtrt_timeb[i]) l.push_back(i);
                      } else {
                        if (tstopb[i] < os_timeb[i]) l.push_back(i);
                      }
                    }
                    
                    std::vector<int> id2 = subset(idb, l);
                    std::vector<int> stratum2 = subset(stratumb, l);
                    std::vector<double> tstart2 = subset(tstartb, l);
                    std::vector<double> tstop2 = subset(tstopb, l);
                    std::vector<double> pd_time2 = subset(pd_timeb, l);
                    std::vector<int> os2 = subset(osb, l);
                    std::vector<double> os_time2 = subset(os_timeb, l);
                    std::vector<double> censor_time2 = subset(censor_timeb, l);
                    std::vector<int> swtrt2 = subset(swtrtb, l);
                    std::vector<double> swtrt_time2 = subset(swtrt_timeb, l);
                    FlatMatrix z_lgs2 = subset_flatmatrix(z_lgsb, l);
                    int n2 = static_cast<int>(l.size());
                    
                    // treatment switching indicators
                    std::vector<int> cross2(n2);
                    for (int i = 0; i < n2; ++i) {
                      if (swtrt2[i] == 1 && tstop2[i] >= swtrt_time2[i]) {
                        cross2[i] = 1;                        
                      }
                    }
                    
                    // obtain natural cubic spline knots
                    FlatMatrix s(n2, ns_df);
                    if (ns_df > 0) {
                      std::vector<double> x;
                      x.reserve(n2);
                      for (int i = 0; i < n2; ++i) {
                        if (cross2[i] == 1) x.push_back(tstop2[i]);
                      }
                      ListCpp out = nscpp(x, ns_df);
                      auto knots = out.get<std::vector<double>>("knots");
                      auto b_knots = out.get<std::vector<double>>("boundary_knots");
                      ListCpp out2 = nscpp(tstop2, ns_df, knots, 0, b_knots);
                      s = out2.get<FlatMatrix>("basis");
                    }
                    
                    // re-baseline based on disease progression date
                    std::vector<int> idx2(1, 0);
                    for (int i = 1; i < n2; ++i) {
                      if (id2[i] != id2[i-1]) {
                        idx2.push_back(i);
                      }
                    }
                    int nids2 = static_cast<int>(idx2.size());
                    idx2.push_back(n2);
                    
                    std::vector<int> id3(nids2), stratum3(nids2);
                    std::vector<int> os3(nids2), swtrt3(nids2);
                    std::vector<double> os_time3(nids2), censor_time3(nids2);
                    std::vector<double> swtrt_time3(nids2, NaN);
                    for (int i = 0; i < nids2; ++i) {
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
                    double psihat = NaN, psilower = NaN, psiupper = NaN;
                    std::vector<double> psihat_vec;
                    std::vector<double> Z(n_eval_z, NaN);
                    if (gridsearch) {
                      for (int i = 0; i < n_eval_z; ++i) {
                        ListCpp out = est_psi_tsegest(
                          n2, q, p2, nids2, idx2, stratum3, 
                          os3, os_time3, censor_time3, 
                          swtrt3, swtrt_time3, id2, cross2, 
                          tstart2, tstop2, covariates_lgs, 
                          z_lgs2, ns_df, s, firth, flic, 
                          recensor, alpha, ties, offset, psi[i]);
                        
                        bool fail = out.get<bool>("fail");
                        if (!fail) Z[i] = out.get<double>("z_counterfactual");
                      }
                      
                      ListCpp psihat_list = getpsiest(0, psi, Z, 0);
                      psihat = psihat_list.get<double>("selected_root");
                      psihat_vec = psihat_list.get<std::vector<double>>("all_roots");
                      psi_CI_type = "grid search";
                      
                      if (k == -1) {
                        ListCpp psilower_list = getpsiest(zcrit, psi, Z, -1);
                        psilower = psilower_list.get<double>("selected_root");
                        ListCpp psiupper_list = getpsiest(-zcrit, psi, Z, 1);
                        psiupper = psiupper_list.get<double>("selected_root");
                      }
                    } else {
                      // z-stat for the slope of counterfactual survival time
                      // martingale residuals in logistic regression model
                      double target = 0.0;
                      auto g = [&target, n2, q, p2, nids2, idx2, stratum3, 
                                os3, os_time3, censor_time3, swtrt3, 
                                swtrt_time3, id2, cross2, tstart2, tstop2, 
                                covariates_lgs, z_lgs2, ns_df, s, 
                                firth, flic, recensor, alpha, ties, 
                                offset](double x) -> double{
                                  ListCpp out = est_psi_tsegest(
                                    n2, q, p2, nids2, idx2, stratum3, 
                                    os3, os_time3, censor_time3, 
                                    swtrt3, swtrt_time3, id2, cross2, 
                                    tstart2, tstop2, covariates_lgs, 
                                    z_lgs2, ns_df, s, firth, flic, 
                                    recensor, alpha, ties, offset, x);
                                  
                                  bool fail = out.get<bool>("fail");
                                  if (!fail) {
                                    double z = out.get<double>("z_counterfactual");
                                    return z - target;
                                  } else {
                                    return NaN;
                                  }
                                };
                      
                      // causal parameter estimates
                      double psilo = getpsiend(g, true, low_psi);
                      double psihi = getpsiend(g, false, hi_psi);
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
                        psilo = getpsiend(g, true, low_psi);
                        psihi = getpsiend(g, false, hi_psi);
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
                        psilo = getpsiend(g, true, low_psi);
                        psihi = getpsiend(g, false, hi_psi);
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
                    
                    if (k == -1) {
                      // obtain data and KM plot for time to switch
                      std::vector<double> swtrt_timen4(nids2);
                      for (int i = 0; i < nids2; ++i) {
                        swtrt_timen4[i] = (swtrt3[i] == 1) ? 
                        swtrt_time3[i] : os_time3[i];
                      }
                      
                      DataFrameCpp ds;
                      ds.push_back(swtrt3, "swtrt");
                      ds.push_back(swtrt_timen4, "swtrt_time");
                      if (data.int_cols.count(id)) {
                        ds.push_front(subset(idwi, id3), id);
                      } else if (data.numeric_cols.count(id)) {
                        ds.push_front(subset(idwn, id3), id);
                      } else if (data.string_cols.count(id)) {
                        ds.push_front(subset(idwc, id3), id);
                      }
                      
                      if (has_stratum) {
                        for (int i = 0; i < p_stratum; ++i) {
                          const std::string& s = stratum[i];
                          if (data.bool_cols.count(s)) {
                            auto v = u_stratum.get<unsigned char>(s);
                            ds.push_back(subset(v, stratum3), s);
                          } else if (data.int_cols.count(s)) {
                            auto v = u_stratum.get<int>(s);
                            ds.push_back(subset(v, stratum3), s);
                          } else if (data.numeric_cols.count(s)) {
                            auto v = u_stratum.get<double>(s);
                            ds.push_back(subset(v, stratum3), s);
                          } else if (data.string_cols.count(s)) {
                            auto v = u_stratum.get<std::string>(s);
                            ds.push_back(subset(v, stratum3), s);
                          }
                        }
                      }
                      
                      DataFrameCpp km = kmestcpp(
                        ds, {""}, "swtrt_time", "", "swtrt", "", 
                        "log-log", 1.0 - alpha, 1);
                      
                      // obtain the Wald statistics for the coefficient of 
                      // the counterfactual in the logistic regression 
                      // switching model at a sequence of psi values
                      if (!gridsearch) {
                        for (int i = 0; i < n_eval_z; ++i) {
                          ListCpp out = est_psi_tsegest(
                            n2, q, p2, nids2, idx2, stratum3, 
                            os3, os_time3, censor_time3, 
                            swtrt3, swtrt_time3, id2, cross2, 
                            tstart2, tstop2, covariates_lgs, 
                            z_lgs2, ns_df, s, firth, flic, 
                            recensor, alpha, ties, offset, psi[i]);
                          
                          bool fail = out.get<bool>("fail");
                          if (!fail) Z[i] = out.get<double>("z_counterfactual");
                        }
                        
                        ListCpp psihat_list = getpsiest(0, psi, Z, 0);
                        psihat_vec = psihat_list.get<std::vector<double>>("all_roots");
                      }
                      
                      DataFrameCpp ez;
                      ez.push_back(psi, "psi");
                      ez.push_back(Z, "Z");
                      
                      ListPtr& dsp = data_switch[h];
                      dsp->get<DataFrameCpp>("data") = ds;
                      
                      ListPtr& kmp = km_switch[h];
                      kmp->get<DataFrameCpp>("data") = km;
                      
                      ListPtr& ezp = eval_z[h];
                      ezp->get<DataFrameCpp>("data") = ez;
                    }
                    
                    if (!std::isnan(psihat)) {
                      // calculate counter-factual survival times
                      double a = std::exp(psihat);
                      double c0 = std::min(1.0, a);
                      for (int i = 0; i < nids; ++i) {
                        if (treat1[i] == h) {
                          double u_star;
                          if (swtrt1[i] == 1) {
                            double b2 = swtrt_time1[i] - offset;
                            u_star = b2 + (time1[i] - b2) * a;
                          } else {
                            u_star = time1[i];
                          }
                          
                          if (recensor) {
                            double c_star = censor_time1[i] * c0;
                            t_star[i] = std::min(u_star, c_star);
                            d_star[i] = c_star < u_star ? 0 : event1[i];
                          } else {
                            t_star[i] = u_star;
                            d_star[i] = event1[i];
                          }
                        }
                      }
                      
                      if (k == -1) {
                        // obtain data and fit for null Cox & logistic models
                        ListCpp out = est_psi_tsegest(
                          n2, q, p2, nids2, idx2, stratum3, 
                          os3, os_time3, censor_time3, 
                          swtrt3, swtrt_time3, id2, cross2, 
                          tstart2, tstop2, covariates_lgs, 
                          z_lgs2, ns_df, s, firth, flic, 
                          recensor, alpha, ties, offset, psihat);
                        if (out.get<bool>("fail")) fail = true;
                        
                        DataFrameCpp dn = out.get<DataFrameCpp>("data_nullcox");
                        DataFrameCpp dl = out.get<DataFrameCpp>("data_logis");
                        ListCpp fn = out.get_list("fit_nullcox");
                        ListCpp fl = out.get_list("fit_logis");
                        
                        if (data.int_cols.count(id)) {
                          dn.push_front(subset(idwi, id3), id);
                          dl.push_front(subset(idwi, id2), id);
                        } else if (data.numeric_cols.count(id)) {
                          dn.push_front(subset(idwn, id3), id);
                          dl.push_front(subset(idwn, id2), id);
                        } else if (data.string_cols.count(id)) {
                          dn.push_front(subset(idwc, id3), id);
                          dl.push_front(subset(idwc, id2), id);
                        }
                        
                        if (has_stratum) {
                          for (int i=0; i<p_stratum; ++i) {
                            const std::string& s = stratum[i];
                            if (data.bool_cols.count(s)) {
                              auto v = u_stratum.get<unsigned char>(s);
                              dn.push_back(subset(v, stratum3), s);
                              dl.push_back(subset(v, stratum2), s);
                            } else if (data.int_cols.count(s)) {
                              auto v = u_stratum.get<int>(s);
                              dn.push_back(subset(v, stratum3), s);
                              dl.push_back(subset(v, stratum2), s);
                            } else if (data.numeric_cols.count(s)) {
                              auto v = u_stratum.get<double>(s);
                              dn.push_back(subset(v, stratum3), s);
                              dl.push_back(subset(v, stratum2), s);
                            } else if (data.string_cols.count(s)) {
                              auto v = u_stratum.get<std::string>(s);
                              dn.push_back(subset(v, stratum3), s);
                              dl.push_back(subset(v, stratum2), s);
                            }
                          }
                        }
                        
                        // update the data and model fits
                        ListPtr& dnp = data_nullcox[h];
                        dnp->get<DataFrameCpp>("data") = dn;
                        
                        ListPtr& dlp = data_logis[h];
                        dlp->get<DataFrameCpp>("data") = dl;
                        
                        ListPtr& fnp = fit_nullcox[h];
                        fnp->get_list("fit") = fn;

                        ListPtr& flp = fit_logis[h];
                        flp->get_list("fit") = fl;
                      }
                    } else {
                      psimissing = true;
                      fail = true;
                    }
                  }
                  
                  if (!psimissing) {
                    // Cox model for hypothetical treatment effect estimate
                    data_outcome.push_back(id1, "uid");
                    data_outcome.push_back(t_star, "t_star");
                    data_outcome.push_back(d_star, "d_star");
                    data_outcome.push_back(treat1, "treated");
                    data_outcome.push_back(stratum1, "ustratum");
                    
                    for (int j = 0; j < p; ++j) {
                      const std::string& zj = covariates[j+1];
                      std::vector<double> u = flatmatrix_get_column(z1, j);
                      data_outcome.push_back(std::move(u), zj);
                    }
                    
                    // generate KM estimate and log-rank test
                    if (k == -1) {
                      km_outcome = kmestcpp(
                        data_outcome, {"treated"}, "t_star", "", 
                        "d_star", "", "log-log", 1.0 - alpha, 1);
                      
                      lr_outcome = lrtestcpp(
                        data_outcome, {"ustratum"}, "treated", "t_star", 
                        "", "d_star");
                    }
                    
                    // fit the outcome model
                    fit_outcome = phregcpp(
                      data_outcome, {"ustratum"}, "t_star", "", "d_star", 
                      covariates, "", "", "", ties, init, 0, 0, 0, 0, 0, alpha);
                    
                    DataFrameCpp sumstat = fit_outcome.get<DataFrameCpp>("sumstat");
                    if (sumstat.get<unsigned char>("fail")[0]) fail = true;
                    
                    DataFrameCpp parest = fit_outcome.get<DataFrameCpp>("parest");
                    double beta0 = parest.get<double>("beta")[0];
                    double sebeta0 = parest.get<double>("sebeta")[0];
                    hrhat = std::exp(beta0);
                    if (k == -1) {
                      hrlower = std::exp(beta0 - zcrit * sebeta0);
                      hrupper = std::exp(beta0 + zcrit * sebeta0);
                      pvalue = parest.get<double>("p")[0];
                    }
                  }
                  
                  ListCpp out;
                  if (k == -1) {
                    out.push_back(std::move(data_switch), "data_switch");
                    out.push_back(std::move(km_switch), "km_switch");
                    out.push_back(std::move(eval_z), "eval_z");
                    out.push_back(std::move(data_nullcox), "data_nullcox");
                    out.push_back(std::move(fit_nullcox), "fit_nullcox");
                    out.push_back(std::move(data_logis), "data_logis");
                    out.push_back(std::move(fit_logis), "fit_logis");
                    out.push_back(std::move(data_outcome), "data_outcome");
                    out.push_back(std::move(km_outcome), "km_outcome");
                    out.push_back(std::move(lr_outcome), "lr_outcome");
                    out.push_back(std::move(fit_outcome), "fit_outcome");
                    out.push_back(psi0hat, "psihat");
                    out.push_back(std::move(psi0hat_vec), "psihat_vec");
                    out.push_back(psi0lower, "psilower");
                    out.push_back(psi0upper, "psiupper");
                    out.push_back(psi1hat, "psi1hat");
                    out.push_back(std::move(psi1hat_vec), "psi1hat_vec");
                    out.push_back(psi1lower, "psi1lower");
                    out.push_back(psi1upper, "psi1upper");
                    out.push_back(psi_CI_type, "psi_CI_type");
                    out.push_back(hrhat, "hrhat");
                    out.push_back(hrlower, "hrlower");
                    out.push_back(hrupper, "hrupper");
                    out.push_back(pvalue, "pvalue");
                    out.push_back(fail, "fail");
                    out.push_back(psimissing, "psimissing");
                  } else {
                    out.push_back(psi0hat, "psihat");
                    out.push_back(psi1hat, "psi1hat");
                    out.push_back(hrhat, "hrhat");
                    out.push_back(fail, "fail");
                  }
                  
                  return out;
                };

  ListCpp out = f(idn, stratumn, tstartn, tstopn, eventn, treatn,
                  osn, os_timen, censor_timen, pdn, pd_timen, 
                  swtrtn, swtrt_timen, zn, z_lgsn, -1);
  
  auto data_switch = out.get<std::vector<ListPtr>>("data_switch");
  auto km_switch = out.get<std::vector<ListPtr>>("km_switch");
  auto eval_z = out.get<std::vector<ListPtr>>("eval_z");
  auto data_nullcox = out.get<std::vector<ListPtr>>("data_nullcox");
  auto fit_nullcox = out.get<std::vector<ListPtr>>("fit_nullcox");
  auto data_logis = out.get<std::vector<ListPtr>>("data_logis");
  auto fit_logis = out.get<std::vector<ListPtr>>("fit_logis");
  DataFrameCpp data_outcome = out.get<DataFrameCpp>("data_outcome");
  DataFrameCpp km_outcome = out.get<DataFrameCpp>("km_outcome");
  DataFrameCpp lr_outcome = out.get<DataFrameCpp>("lr_outcome");
  ListCpp fit_outcome = out.get_list("fit_outcome");
  double psihat = out.get<double>("psihat");
  std::vector<double> psihat_vec = out.get<std::vector<double>>("psihat_vec");
  double psilower = out.get<double>("psilower");
  double psiupper = out.get<double>("psiupper");
  double psi1hat = out.get<double>("psi1hat");
  std::vector<double> psi1hat_vec = out.get<std::vector<double>>("psi1hat_vec");
  double psi1lower = out.get<double>("psi1lower");
  double psi1upper = out.get<double>("psi1upper");
  std::string psi_CI_type = out.get<std::string>("psi_CI_type");
  double hrhat = out.get<double>("hrhat");
  double hrlower = out.get<double>("hrlower");
  double hrupper = out.get<double>("hrupper");
  double pvalue = out.get<double>("pvalue");
  bool fail = out.get<bool>("fail");
  bool psimissing = out.get<bool>("psimissing");
  
  std::vector<double> hrhats(n_boot), psihats(n_boot), psi1hats(n_boot);
  std::vector<unsigned char> fails(n_boot);
  DataFrameCpp fail_boots_data;
  std::string hr_CI_type;
  
  if (!psimissing) {
    // summarize number of deaths by treatment arm in the outcome data
    std::vector<int> treated = data_outcome.get<int>("treated");
    std::vector<int> event_out = data_outcome.get<int>("d_star");
    std::vector<double> n_event_out(2);
    // note: outcome data excludes data after switch
    for (int i = 0; i < static_cast<int>(treated.size()); ++i) {
      int g = treated[i];
      if (event_out[i] == 1) ++n_event_out[g];
    }
    std::vector<double> pct_event_out(2);
    for (int g = 0; g < 2; g++) {
      pct_event_out[g] = 100.0 * n_event_out[g] / n_total[g];
    }
    event_summary.push_back(std::move(n_event_out), "event_out_n");
    event_summary.push_back(std::move(pct_event_out), "event_out_pct");
    
    std::vector<int> uid = data_outcome.get<int>("uid");
    if (data.int_cols.count(id)) {
      data_outcome.push_front(subset(idwi, uid), id);
    } else if (data.numeric_cols.count(id)) {
      data_outcome.push_front(subset(idwn, uid), id);
    } else if (data.string_cols.count(id)) {
      data_outcome.push_front(subset(idwc, uid), id);
    }
    
    treated = event_summary.get<int>("treated");
    std::vector<int> nottreated(treated.size());
    std::transform(treated.begin(), treated.end(), nottreated.begin(),
                   [](int value) { return 1 - value; });
    if (data.bool_cols.count(treat) || data.int_cols.count(treat)) {
      event_summary.push_back(subset(treatwi, nottreated), treat);
    } else if (data.numeric_cols.count(treat)) {
      event_summary.push_back(subset(treatwn, nottreated), treat);
    } else if (data.string_cols.count(treat)) {
      event_summary.push_back(subset(treatwc, nottreated), treat);
    }
    
    treated = data_outcome.get<int>("treated");
    nottreated.resize(treated.size());
    std::transform(treated.begin(), treated.end(), nottreated.begin(),
                   [](int value) { return 1 - value; });
    if (data.bool_cols.count(treat) || data.int_cols.count(treat)) {
      data_outcome.push_back(subset(treatwi, nottreated), treat);
    } else if (data.numeric_cols.count(treat)) {
      data_outcome.push_back(subset(treatwn, nottreated), treat);
    } else if (data.string_cols.count(treat)) {
      data_outcome.push_back(subset(treatwc, nottreated), treat);
    }
    
    treated = km_outcome.get<int>("treated");
    nottreated.resize(treated.size());
    std::transform(treated.begin(), treated.end(), nottreated.begin(),
                   [](int value) { return 1 - value; });
    if (data.bool_cols.count(treat) || data.int_cols.count(treat)) {
      km_outcome.push_back(subset(treatwi, nottreated), treat);
    } else if (data.numeric_cols.count(treat)) {
      km_outcome.push_back(subset(treatwn, nottreated), treat);
    } else if (data.string_cols.count(treat)) {
      km_outcome.push_back(subset(treatwc, nottreated), treat);
    }
    
    if (has_stratum) {
      std::vector<int> ustratum = data_outcome.get<int>("ustratum");
      for (int i = 0; i < p_stratum; ++i) {
        const std::string& s = stratum[i];
        if (data.bool_cols.count(s)) {
          auto v = u_stratum.get<unsigned char>(s);
          data_outcome.push_back(subset(v, ustratum), s);
        } else if (data.int_cols.count(s)) {
          auto v = u_stratum.get<int>(s);
          data_outcome.push_back(subset(v, ustratum), s);
        } else if (data.numeric_cols.count(s)) {
          auto v = u_stratum.get<double>(s);
          data_outcome.push_back(subset(v, ustratum), s);
        } else if (data.string_cols.count(s)) {
          auto v = u_stratum.get<std::string>(s);
          data_outcome.push_back(subset(v, ustratum), s);
        }
      }
    }
    
    // construct the confidence interval for HR
    if (!boot) { // use Cox model to construct CI for HR if no boot
      hr_CI_type = "Cox model";
    } else { // bootstrap the entire process to construct CI for HR
      if (has_stratum) {
        // sort data by treatment group, stratum, id, and time
        std::vector<int> order = seqcpp(0, n-1);
        std::sort(order.begin(), order.end(), [&](int i, int j) {
          return std::tie(treatn[i], stratumn[i], idn[i], tstopn[i]) <
            std::tie(treatn[j], stratumn[j], idn[j], tstopn[j]);
        });
        
        subset_in_place(idn, order);
        subset_in_place(stratumn, order);
        subset_in_place(tstartn, order);
        subset_in_place(tstopn, order);
        subset_in_place(eventn, order);
        subset_in_place(treatn, order);
        subset_in_place(osn, order);
        subset_in_place(os_timen, order);
        subset_in_place(censor_timen, order);
        subset_in_place(pdn, order);
        subset_in_place(pd_timen, order);
        subset_in_place(swtrtn, order);
        subset_in_place(swtrt_timen, order);
        subset_in_place_flatmatrix(zn, order);
        subset_in_place_flatmatrix(z_lgsn, order);
      }
      
      std::vector<int> idx(1, 0); // first observation within an id
      for (int i = 1; i < n; ++i) {
        if (idn[i] != idn[i-1]) {
          idx.push_back(i);
        }
      }
      
      int nids = static_cast<int>(idx.size());
      idx.push_back(n);
      
      std::vector<int> idx1(nids); // last observation within an id
      for (int i = 0; i < nids; ++i) {
        idx1[i] = idx[i+1] - 1;
      }
      
      std::vector<int> treat1 = subset(treatn, idx1);
      std::vector<int> stratum1 = subset(stratumn, idx1);
      
      std::vector<int> tsx(1, 0); // first id within each treat/stratum
      for (int i = 1; i < nids; ++i) {
        if (treat1[i] != treat1[i-1] || stratum1[i] != stratum1[i-1]) {
          tsx.push_back(i);
        }
      }
      
      int ntss = static_cast<int>(tsx.size());
      tsx.push_back(nids); // add the end index
      
      // Before running the parallel loop: pre-generate deterministic seeds
      std::vector<uint64_t> seeds(n_boot);
      boost::random::mt19937_64 master_rng(static_cast<uint64_t>(seed));
      for (int k = 0; k < n_boot; ++k) seeds[k] = master_rng();
      
      // We'll collect failure bootstrap data per-worker and merge via Worker::join.
      struct BootstrapWorker : public RcppParallel::Worker {
        // references to read-only inputs (no mutation)
        const int n;
        const int nids;
        const int ntss;
        const std::vector<int>& idx;
        const std::vector<int>& tsx;
        const std::vector<int>& idn;
        const std::vector<int>& stratumn;
        const std::vector<double>& tstartn;
        const std::vector<double>& tstopn;
        const std::vector<int>& eventn;
        const std::vector<int>& treatn;
        const std::vector<int>& osn;
        const std::vector<double>& os_timen;
        const std::vector<double>& censor_timen;
        const std::vector<int>& pdn;
        const std::vector<double>& pd_timen;
        const std::vector<int>& swtrtn;
        const std::vector<double>& swtrt_timen;
        const FlatMatrix& zn;
        const FlatMatrix& z_lgsn;
        const std::vector<uint64_t>& seeds;
        
        // function f and other params that f needs are captured from outer scope
        // capture them by reference here so worker can call f(...)
        std::function<ListCpp(const std::vector<int>&, 
                              const std::vector<int>&,
                              const std::vector<double>&, 
                              const std::vector<double>&,
                              const std::vector<int>&, 
                              const std::vector<int>&,
                              const std::vector<int>&,
                              const std::vector<double>&,
                              const std::vector<double>&,
                              const std::vector<int>&,
                              const std::vector<double>&,
                              const std::vector<int>&,
                              const std::vector<double>&,
                              const FlatMatrix&,
                              const FlatMatrix&, int)> f;
        
        // result references (each iteration writes unique index into these)
        std::vector<unsigned char>& fails_out;
        std::vector<double>& hrhats_out;
        
        // Per-worker storage for failed-boot data (to be merged in join)
        std::vector<int> boot_indexc_local;
        std::vector<int> oidc_local;
        std::vector<int> idc_local;
        std::vector<int> stratumc_local;
        std::vector<double> tstartc_local;
        std::vector<double> tstopc_local;
        std::vector<int> eventc_local;
        std::vector<int> treatc_local;
        std::vector<int> osc_local;
        std::vector<double> os_timec_local;
        std::vector<double> censor_timec_local;
        std::vector<int> pdc_local;
        std::vector<double> pd_timec_local;
        std::vector<int> swtrtc_local;
        std::vector<double> swtrt_timec_local;
        
        // store column-wise z_lgsc_local: outer vector length == z_lgsn.ncol
        // each inner vector stores column data across failed boots
        std::vector<std::vector<double>> z_lgsc_local;
        
        int index1_local = 0; // number of rows stored so far for z_lgsc_local
        
        // constructor
        BootstrapWorker(int n_, int nids_, int ntss_,
                        const std::vector<int>& idx_,
                        const std::vector<int>& tsx_,
                        const std::vector<int>& idn_,
                        const std::vector<int>& stratumn_,
                        const std::vector<double>& tstartn_,
                        const std::vector<double>& tstopn_,
                        const std::vector<int>& eventn_,
                        const std::vector<int>& treatn_,
                        const std::vector<int>& osn_,
                        const std::vector<double>& os_timen_,
                        const std::vector<double>& censor_timen_,
                        const std::vector<int>& pdn_,
                        const std::vector<double>& pd_timen_,
                        const std::vector<int>& swtrtn_,
                        const std::vector<double>& swtrt_timen_,
                        const FlatMatrix& zn_,
                        const FlatMatrix& z_lgsn_,
                        const std::vector<uint64_t>& seeds_,
                        decltype(f) f_,
                        std::vector<unsigned char>& fails_out_,
                        std::vector<double>& hrhats_out_) :
          n(n_), nids(nids_), ntss(ntss_), idx(idx_), tsx(tsx_), idn(idn_),
          stratumn(stratumn_), tstartn(tstartn_), tstopn(tstopn_),
          eventn(eventn_), treatn(treatn_), osn(osn_), os_timen(os_timen_),
          censor_timen(censor_timen_), pdn(pdn_), pd_timen(pd_timen_),
          swtrtn(swtrtn_), swtrt_timen(swtrt_timen_), zn(zn_),
          z_lgsn(z_lgsn_), seeds(seeds_), f(std::move(f_)),
          fails_out(fails_out_), hrhats_out(hrhats_out_) {
          
          // heuristic reservation to reduce reallocations:
          int ncols_lgs = z_lgsn.ncol;
          z_lgsc_local.resize(ncols_lgs);
          
          // reserve some capacity per column. This is a heuristic; adjust if needed.
          std::size_t per_col_reserve = 10 * n;
          for (int col = 0; col < ncols_lgs; ++col)
            z_lgsc_local[col].reserve(per_col_reserve);
          
          // Reserve scalar buffers heuristically (reduce reallocs)
          boot_indexc_local.reserve(per_col_reserve);
          oidc_local.reserve(per_col_reserve);
          idc_local.reserve(per_col_reserve);
          stratumc_local.reserve(per_col_reserve);
          tstartc_local.reserve(per_col_reserve);
          tstopc_local.reserve(per_col_reserve);
          eventc_local.reserve(per_col_reserve);
          treatc_local.reserve(per_col_reserve);
          osc_local.reserve(per_col_reserve);
          os_timec_local.reserve(per_col_reserve);
          censor_timec_local.reserve(per_col_reserve);
          pdc_local.reserve(per_col_reserve);
          pd_timec_local.reserve(per_col_reserve);
          swtrtc_local.reserve(per_col_reserve);
          swtrt_timec_local.reserve(per_col_reserve);
        }
        
        // operator() processes a range of bootstrap iterations [begin, end)
        void operator()(std::size_t begin, std::size_t end) {
          // per-worker reusable buffers (avoid reallocation per iteration)
          std::vector<int> oidb, idb, stratumb, eventb, treatb, osb, pdb, swtrtb;
          std::vector<double> tstartb, tstopb, os_timeb, censor_timeb;
          std::vector<double> pd_timeb, swtrt_timeb;
          std::vector<std::vector<double>> zb_cols, z_lgsb_cols;
          int ncols_z = zn.ncol, ncols_lgs = z_lgsn.ncol;
          zb_cols.resize(ncols_z); z_lgsb_cols.resize(ncols_lgs);
          
          oidb.reserve(n); idb.reserve(n); stratumb.reserve(n); 
          tstartb.reserve(n); tstopb.reserve(n); eventb.reserve(n); 
          treatb.reserve(n); osb.reserve(n); os_timeb.reserve(n); 
          censor_timeb.reserve(n); pdb.reserve(n); pd_timeb.reserve(n);
          swtrtb.reserve(n); swtrt_timeb.reserve(n);
          for (int col = 0; col < ncols_z; ++col) zb_cols[col].reserve(n);
          for (int col = 0; col < ncols_lgs; ++col) z_lgsb_cols[col].reserve(n);
          
          for (std::size_t k = begin; k < end; ++k) {
            // Reset per-iteration buffers so each k starts fresh.
            // (They were reserved earlier for performance but must be emptied.)
            oidb.clear(); idb.clear(); stratumb.clear(); tstartb.clear();
            tstopb.clear(); eventb.clear(); treatb.clear(); osb.clear();
            os_timeb.clear(); censor_timeb.clear(); pdb.clear(); 
            pd_timeb.clear(); swtrtb.clear(); swtrt_timeb.clear();
            for (int col = 0; col < ncols_z; ++col) zb_cols[col].clear();
            for (int col = 0; col < ncols_lgs; ++col) z_lgsb_cols[col].clear();
            
            // deterministic RNG per-iteration
            std::mt19937_64 rng(seeds[k]);
            
            // sample by treatment/stratum blocks
            for (int h = 0; h < ntss; ++h) {
              int start = tsx[h], end = tsx[h + 1];
              int len = end - start;
              boost::random::uniform_int_distribution<int> index_dist(0, len - 1);
              for (int r = start; r < end; ++r) {
                int i = start + index_dist(rng);
                int oidb1 = idn[idx[i]];     // original id
                int idb1 = oidb1 + r * nids; // unique id
                
                int start1 = idx[i], end1 = idx[i+1]; // within-id block
                int len1 = end1 - start1;
                
                std::vector<int> oidn1(len1, oidb1);
                std::vector<int> idn1(len1, idb1);
                std::vector<int> stratumn1 = subset(stratumn, start1, end1);
                std::vector<double> tstartn1 = subset(tstartn, start1, end1);
                std::vector<double> tstopn1 = subset(tstopn, start1, end1);
                std::vector<int> eventn1 = subset(eventn, start1, end1);
                std::vector<int> treatn1 = subset(treatn, start1, end1);
                std::vector<int> osn1 = subset(osn, start1, end1);
                std::vector<double> os_timen1 = subset(os_timen, start1, end1);
                std::vector<double> censor_timen1 = subset(censor_timen, start1, end1);
                std::vector<int> pdn1 = subset(pdn, start1, end1);
                std::vector<double> pd_timen1 = subset(pd_timen, start1, end1);
                std::vector<int> swtrtn1 = subset(swtrtn, start1, end1);
                std::vector<double> swtrt_timen1 = subset(swtrt_timen, start1, end1);
                FlatMatrix zn1 = subset_flatmatrix(zn, start1, end1);
                FlatMatrix z_lgsn1 = subset_flatmatrix(z_lgsn, start1, end1);
                
                append(oidb, oidn1);
                append(idb, idn1);
                append(stratumb, stratumn1);
                append(tstartb, tstartn1);
                append(tstopb, tstopn1);
                append(eventb, eventn1);
                append(treatb, treatn1);
                append(osb, osn1);
                append(os_timeb, os_timen1);
                append(censor_timeb, censor_timen1);
                append(pdb, pdn1);
                append(pd_timeb, pd_timen1);
                append(swtrtb, swtrtn1);
                append(swtrt_timeb, swtrt_timen1);
                append_flatmatrix(zb_cols, zn1);
                append_flatmatrix(z_lgsb_cols, z_lgsn1);
              }
            } // end block sampling
            
            // call the (thread-safe) per-iteration function f
            FlatMatrix zb = cols_to_flatmatrix(zb_cols);
            FlatMatrix z_lgsb = cols_to_flatmatrix(z_lgsb_cols);
            
            ListCpp out = f(idb, stratumb, tstartb, tstopb, eventb, treatb, 
                            osb, os_timeb, censor_timeb, pdb, pd_timeb,
                            swtrtb, swtrt_timeb, zb, z_lgsb, 
                            static_cast<int>(k));
            
            // write results
            fails_out[k] = out.get<bool>("fail");
            hrhats_out[k] = out.get<double>("hrhat");
            
            // existing code that collects failure data when fails_out[k] is true...
            if (fails_out[k]) {
              int l = static_cast<int>(idb.size());
              append(boot_indexc_local, std::vector<int>(l, k+1));
              append(oidc_local, oidb);
              append(idc_local, idb);
              append(stratumc_local, stratumb);
              append(tstartc_local, tstartb);
              append(tstopc_local, tstopb);
              append(eventc_local, eventb);
              append(treatc_local, treatb);
              append(osc_local, osb);
              append(os_timec_local, os_timeb);
              append(censor_timec_local, censor_timeb);
              append(pdc_local, pdb);
              append(pd_timec_local, pd_timeb);
              append(swtrtc_local, swtrtb);
              append(swtrt_timec_local, swtrt_timeb);
              append_flatmatrix(z_lgsc_local, z_lgsb);
              index1_local += l;
            }
          } // end for k
        } // end operator()
        
        // join merges other (worker copy) into this instance (called by parallelFor)
        void join(const BootstrapWorker& other) {
          // append other's local vectors into this worker's storage
          append(boot_indexc_local, other.boot_indexc_local);
          append(oidc_local, other.oidc_local);
          append(idc_local, other.idc_local);
          append(stratumc_local, other.stratumc_local);
          append(tstartc_local, other.tstartc_local);
          append(tstopc_local, other.tstopc_local);
          append(eventc_local, other.eventc_local);
          append(treatc_local, other.treatc_local);
          append(osc_local, other.osc_local);
          append(os_timec_local, other.os_timec_local);
          append(censor_timec_local, other.censor_timec_local);
          append(pdc_local, other.pdc_local);
          append(pd_timec_local, other.pd_timec_local);
          append(swtrtc_local, other.swtrtc_local);
          append(swtrt_timec_local, other.swtrt_timec_local);
          append_flatmatrix(z_lgsc_local, other.z_lgsc_local);
          index1_local += other.index1_local;
        }
      }; // end BootstrapWorker
      
      // Instantiate the Worker with references to inputs and outputs
      BootstrapWorker worker(
          n, nids, ntss, idx, tsx, idn, stratumn, tstartn, tstopn, 
          eventn, treatn, osn, os_timen, censor_timen, pdn, pd_timen, 
          swtrtn, swtrt_timen, zn, z_lgsn, seeds,
          // bind f into std::function (capture the f we already have)
          std::function<ListCpp(const std::vector<int>&, 
                                const std::vector<int>&,
                                const std::vector<double>&, 
                                const std::vector<double>&,
                                const std::vector<int>&, 
                                const std::vector<int>&,
                                const std::vector<int>&,
                                const std::vector<double>&,
                                const std::vector<double>&,
                                const std::vector<int>&,
                                const std::vector<double>&,
                                const std::vector<int>&,
                                const std::vector<double>&,
                                const FlatMatrix&,
                                const FlatMatrix&, int)>(f),
                                fails, hrhats
      );
      
      // Run the parallel loop over bootstrap iterations
      RcppParallel::parallelFor(0, n_boot, worker, 1 /*grain size*/);
      
      // After parallelFor returns, worker contains merged per-worker failure data.
      // Move them into the outer-scope vectors used later in the function:
      std::vector<int> oidc = std::move(worker.oidc_local);
      std::vector<int> stratumc = std::move(worker.stratumc_local);
      std::vector<int> treatc = std::move(worker.treatc_local);
      
      // assemble the failed bootstrap data into a DataFrame
      if (worker.index1_local > 0) {
        fail_boots_data.push_back(std::move(worker.boot_indexc_local), "boot_index");
        fail_boots_data.push_back(std::move(worker.idc_local), "uid");
        fail_boots_data.push_back(std::move(worker.tstartc_local), "tstart");
        fail_boots_data.push_back(std::move(worker.tstopc_local), "tstop");
        fail_boots_data.push_back(std::move(worker.eventc_local), "event");
        fail_boots_data.push_back(treatc, "treated");
        fail_boots_data.push_back(std::move(worker.osc_local), "os");
        fail_boots_data.push_back(std::move(worker.os_timec_local), "os_time");
        fail_boots_data.push_back(std::move(worker.censor_timec_local), "censor_time");
        fail_boots_data.push_back(std::move(worker.pdc_local), "pd");
        fail_boots_data.push_back(std::move(worker.pd_timec_local), "pd_time");
        fail_boots_data.push_back(std::move(worker.swtrtc_local), "swtrt");
        fail_boots_data.push_back(std::move(worker.swtrt_timec_local), "swtrt_time");
        
        int ncols_lgs = worker.z_lgsc_local.size();
        for (int j = 0; j < ncols_lgs; ++j) {
          const std::string& zj = covariates_lgs[j];
          fail_boots_data.push_back(std::move(worker.z_lgsc_local[j]), zj);
        }
        
        if (data.int_cols.count(id)) {
          fail_boots_data.push_back(subset(idwi, oidc), id);
        } else if (data.numeric_cols.count(id)) {
          fail_boots_data.push_back(subset(idwn, oidc), id);
        } else if (data.string_cols.count(id)) {
          fail_boots_data.push_back(subset(idwc, oidc), id);
        }
        
        std::vector<int> nottreatc(treatc.size());
        std::transform(treatc.begin(), treatc.end(), nottreatc.begin(),
                       [](int value) { return 1 - value; });
        if (data.bool_cols.count(treat) || data.int_cols.count(treat)) {
          fail_boots_data.push_back(subset(treatwi, nottreatc), treat);
        } else if (data.numeric_cols.count(treat)) {
          fail_boots_data.push_back(subset(treatwn, nottreatc), treat);
        } else if (data.string_cols.count(treat)) {
          fail_boots_data.push_back(subset(treatwc, nottreatc), treat);
        }
        
        if (has_stratum) {
          for (int i = 0; i < p_stratum; ++i) {
            const std::string& s = stratum[i];
            if (data.bool_cols.count(s)) {
              auto v = u_stratum.get<unsigned char>(s);
              fail_boots_data.push_back(subset(v, stratumc), s);
            } else if (data.int_cols.count(s)) {
              auto v = u_stratum.get<int>(s);
              fail_boots_data.push_back(subset(v, stratumc), s);
            } else if (data.numeric_cols.count(s)) {
              auto v = u_stratum.get<double>(s);
              fail_boots_data.push_back(subset(v, stratumc), s);
            } else if (data.string_cols.count(s)) {
              auto v = u_stratum.get<std::string>(s);
              fail_boots_data.push_back(subset(v, stratumc), s);
            }
          }
        }
      }
      
      // retrieve the bootstrap results
      fails = worker.fails_out;
      hrhats = worker.hrhats_out;
      
      // obtain bootstrap confidence interval for HR
      double loghr = std::log(hrhat);
      std::vector<int> ok;
      ok.reserve(n_boot);
      int n_ok = 0;
      for (int k = 0; k < n_boot; ++k) {
        if (!fails[k] && !std::isnan(hrhats[k])) {
          ok.push_back(k);
          ++n_ok;
        }
      }
      std::vector<double> hrhats1 = subset(hrhats, ok);
      std::vector<double> loghrs(n_ok);
      std::transform(hrhats1.begin(), hrhats1.end(), loghrs.begin(),
                     [](double value) { return std::log(value); });
      double meanloghr, sdloghr;
      mean_sd(loghrs.data(), n_ok, meanloghr, sdloghr);
      double tcrit = boost_qt(1.0 - alpha / 2.0, n_ok - 1);
      hrlower = std::exp(loghr - tcrit * sdloghr);
      hrupper = std::exp(loghr + tcrit * sdloghr);
      hr_CI_type = "bootstrap";
      pvalue = 2.0 * (1.0 - boost_pt(std::fabs(loghr / sdloghr), n_ok - 1));
      
      // obtain bootstrap confidence interval for psi
      std::vector<double> psihats1 = subset(psihats, ok);
      double meanpsi, sdpsi;
      mean_sd(psihats1.data(), n_ok, meanpsi, sdpsi);
      psilower = psihat - tcrit * sdpsi;
      psiupper = psihat + tcrit * sdpsi;
      psi_CI_type = "bootstrap";
      
      std::vector<double> psi1hats1 = subset(psi1hats, ok);
      double meanpsi1, sdpsi1;
      mean_sd(psi1hats1.data(), n_ok, meanpsi1, sdpsi1);
      psi1lower = psi1hat - tcrit * sdpsi1;
      psi1upper = psi1hat + tcrit * sdpsi1;
    }
  }
  
  ListCpp result;
  std::string pvalue_type = boot ? "bootstrap" : "Cox model";
  std::vector<double> psi_CI = {psilower, psiupper};
  std::vector<double> hr_CI = {hrlower, hrupper};
  result.push_back(psihat, "psi");
  result.push_back(std::move(psihat_vec), "psi_roots");
  result.push_back(std::move(psi_CI), "psi_CI");
  result.push_back(psi_CI_type, "psi_CI_type");
  result.push_back(pvalue, "pvalue");
  result.push_back(pvalue_type, "pvalue_type");
  result.push_back(hrhat, "hr");
  result.push_back(std::move(hr_CI), "hr_CI");
  result.push_back(hr_CI_type, "hr_CI_type");
  result.push_back(std::move(event_summary), "event_summary");
  result.push_back(std::move(data_switch), "data_switch");
  result.push_back(std::move(km_switch), "km_switch");
  result.push_back(std::move(eval_z), "eval_z");
  result.push_back(std::move(data_nullcox), "data_nullcox");
  result.push_back(std::move(fit_nullcox), "fit_nullcox");
  result.push_back(std::move(data_logis), "data_logis");
  result.push_back(std::move(fit_logis), "fit_logis");
  result.push_back(std::move(data_outcome), "data_outcome");
  result.push_back(std::move(km_outcome), "km_outcome");
  result.push_back(std::move(lr_outcome), "lr_outcome");
  result.push_back(std::move(fit_outcome), "fit_outcome");
  result.push_back(fail, "fail");
  result.push_back(psimissing, "psimissing");
  
  if (!swtrt_control_only) {
    result.push_back(psi1hat, "psi_trt");
     result.push_back(std::move(psi1hat_vec), "psi_trt_roots");
    std::vector<double> psi1_CI = {psi1lower, psi1upper};
    result.push_back(std::move(psi1_CI), "psi_trt_CI");
  }
  
  if (boot) {
    result.push_back(fails, "fail_boots");
    result.push_back(std::move(hrhats), "hr_boots");
    result.push_back(std::move(psihats), "psi_boots");
    if (!swtrt_control_only) {
      result.push_back(std::move(psi1hats), "psi_trt_boots");
    }
    if (std::any_of(fails.begin(), fails.end(), [](bool x){ return x; })) {
      result.push_back(std::move(fail_boots_data), "fail_boots_data");
    }
  }
  
  thread_utils::drain_thread_warnings_to_R();
  return Rcpp::wrap(result);
}
