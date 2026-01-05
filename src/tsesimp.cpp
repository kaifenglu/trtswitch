#include <RcppParallel.h>
#include <RcppThread.h>
#include <Rcpp.h>

#include "survival_analysis.h"
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

// [[Rcpp::export]]
Rcpp::List tsesimpcpp(const Rcpp::DataFrame& df,
                      const std::string& id,
                      const std::vector<std::string>& stratum,
                      const std::string& time,
                      const std::string& event,
                      const std::string& treat,
                      const std::string& censor_time,
                      const std::string& pd,
                      const std::string& pd_time,
                      const std::string& swtrt,
                      const std::string& swtrt_time,
                      const std::vector<std::string>& base_cov,
                      const std::vector<std::string>& base2_cov,
                      const std::string& aft_dist,
                      const bool strata_main_effect_only,
                      const bool recensor,
                      const bool admin_recensor_only,
                      const bool swtrt_control_only,
                      const double alpha,
                      const std::string& ties,
                      const double offset,
                      const bool boot,
                      const int n_boot,
                      const int seed) {
  
  DataFrameCpp data = convertRDataFrameToCpp(df);
  
  int n = static_cast<int>(data.nrows());
  int p = static_cast<int>(base_cov.size());
  if (p == 1 && base_cov[0] == "") p = 0;
  
  int p2 = static_cast<int>(base2_cov.size());
  if (p2 == 1 && base2_cov[0] == "") p2 = 0;
  
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
  
  // --- time existence and checks ---
  if (time.empty() || !data.containElementNamed(time))
    throw std::invalid_argument("data must contain the time variable");
  std::vector<double> timen(n);
  if (data.int_cols.count(time)) {
    const std::vector<int>& vi = data.get<int>(time);
    for (int i = 0; i < n; ++i) timen[i] = static_cast<double>(vi[i]);
  } else if (data.numeric_cols.count(time)) {
    timen = data.get<double>(time);
  } else {
    throw std::invalid_argument("time variable must be integer or numeric");
  }
  for (int i = 0; i < n; ++i) {
    if (!std::isnan(timen[i]) && timen[i] < 0.0)
      throw std::invalid_argument("time must be nonnegative");
  }
  
  // --- event variable ---
  if (event.empty() || !data.containElementNamed(event)) {
    throw std::invalid_argument("data must contain the event variable");
  }
  std::vector<double> eventn(n);
  if (data.bool_cols.count(event)) {
    const std::vector<unsigned char>& vb = data.get<unsigned char>(event);
    for (int i = 0; i < n; ++i) eventn[i] = vb[i] ? 1.0 : 0.0;
  } else if (data.int_cols.count(event)) {
    const std::vector<int>& vi = data.get<int>(event);
    for (int i = 0; i < n; ++i) eventn[i] = static_cast<double>(vi[i]);
  } else if (data.numeric_cols.count(event)) {
    eventn = data.get<double>(event);
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
    if (censor_timen[i] < timen[i]) throw std::invalid_argument(
        "censor_time must be greater than or equal to time");
  }
  if (!admin_recensor_only) { // use the actual censoring time for dropouts
    for (int i = 0; i < n; ++i) if (eventn[i] == 0) censor_timen[i] = timen[i];
  }
  
  // --- pd variable ---
  if (pd.empty() || !data.containElementNamed(pd)) {
    throw std::invalid_argument("data must contain the pd variable");
  }
  std::vector<double> pdn(n);
  if (data.bool_cols.count(pd)) {
    const std::vector<unsigned char>& vb = data.get<unsigned char>(pd);
    for (int i = 0; i < n; ++i) pdn[i] = vb[i] ? 1.0 : 0.0;
  } else if (data.int_cols.count(pd)) {
    const std::vector<int>& vi = data.get<int>(pd);
    for (int i = 0; i < n; ++i) pdn[i] = static_cast<double>(vi[i]);
  } else if (data.numeric_cols.count(pd)) {
    pdn = data.get<double>(pd);
  } else {
    throw std::invalid_argument("pd variable must be bool, integer or numeric");
  }
  for (double val : pdn) if (val != 0 && val != 1)
    throw std::invalid_argument("pd must be 1 or 0 for each observation");

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
  std::vector<double> swtrtn(n);
  if (data.bool_cols.count(swtrt)) {
    const std::vector<unsigned char>& vb = data.get<unsigned char>(swtrt);
    for (int i = 0; i < n; ++i) swtrtn[i] = vb[i] ? 1.0 : 0.0;
  } else if (data.int_cols.count(swtrt)) {
    const std::vector<int>& vi = data.get<int>(swtrt);
    for (int i = 0; i < n; ++i) swtrtn[i] = static_cast<double>(vi[i]);
  } else if (data.numeric_cols.count(swtrt)) {
    swtrtn = data.get<double>(swtrt);
  } else {
    throw std::invalid_argument("swtrt variable must be bool, integer or numeric");
  }
  for (double val : swtrtn) if (val != 0 && val != 1)
    throw std::invalid_argument("swtrt must be 1 or 0 for each observation");
  
  // check that at least 1 pd and swtrt in each group
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
    if (swtrtn[i] == 1 && swtrt_timen[i] > timen[i]) {
      throw std::runtime_error("swtrt_time must be less than or equal to time");
    }
  }
  
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
  
  // ensure pd time < os time so that post-progression survival (pps) is positive
  for (int i = 0; i < n; ++i) {
    if (pdn[i] == 1 && pd_timen[i] == timen[i]) {
      timen[i] = timen[i] + 1.0e-8;
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
  
  int qp2 = q + p2;
  
  // covariates for the accelerated failure time model for control with pd
  // including swtrt, stratum and base2_cov
  std::vector<std::string> covariates_aft(qp2 + 1);
  FlatMatrix z_aftn(n, qp2);
  covariates_aft[0] = "swtrt";
  if (has_stratum) {
    if (strata_main_effect_only) {
      int k = 0;
      for (int i = 0; i < p_stratum; ++i) {
        const std::string& s = stratum[i];
        int di = d[i] - 1;
        
        if (u_stratum.string_cols.count(s)) {
          auto u = levels.get<std::vector<std::string>>(s);
          for (int j = 0; j < di; ++j) {
            covariates_aft[k + j + 1] = s + sanitize(u[j]);
          }
        } else if (u_stratum.numeric_cols.count(s)) {
          auto u = levels.get<std::vector<double>>(s);
          for (int j = 0; j < di; ++j) {
            covariates_aft[k + j + 1] = s + std::to_string(u[j]);
          }
        } else if (u_stratum.int_cols.count(s)) {
          auto u = levels.get<std::vector<int>>(s);
          for (int j = 0; j < di; ++j) {
            covariates_aft[k + j + 1] = s + std::to_string(u[j]);
          }
        } else if (u_stratum.bool_cols.count(s)) {
          auto u = levels.get<std::vector<unsigned char>>(s);
          for (int j = 0; j < di; ++j) {
            covariates_aft[k + j + 1] = s + std::to_string(u[j]);
          }
        }
        
        for (int j = 0; j < di; ++j) {
          const int* stratan_col = stratan.data_ptr() + i * n;
          double* z_aftn_col = z_aftn.data_ptr() + (k + j) * n;
          for (int r = 0; r < n; ++r) {
            z_aftn_col[r] = stratan_col[r] == j ? 1.0 : 0.0;
          }
        }
        
        k += di;
      }
    } else {
      for (int j = 0; j < nstrata - 1; ++j) {
        // locate the first observation in the stratum
        int first_k = 0;
        for (; first_k < n; ++first_k) {
          if (stratumn[first_k] == j) break;
        }
        
        covariates_aft[j + 1] = "";
        
        for (int i = 0; i < p_stratum; ++i) {
          const std::string& s = stratum[i];
          
          std::vector<int> q_col = intmatrix_get_column(stratan, i);
          int l = q_col[first_k];
          
          if (u_stratum.string_cols.count(s)) {
            auto u = levels.get<std::vector<std::string>>(s);
            covariates_aft[j + 1] += s + sanitize(u[l]);
          } else if (u_stratum.numeric_cols.count(s)) {
            auto u = levels.get<std::vector<double>>(s);
            covariates_aft[j + 1] += s + std::to_string(u[l]);
          } else if (u_stratum.int_cols.count(s)) {
            auto u = levels.get<std::vector<int>>(s);
            covariates_aft[j + 1] += s + std::to_string(u[l]);
          } else if (u_stratum.bool_cols.count(s)) {
            auto u = levels.get<std::vector<unsigned char>>(s);
            covariates_aft[j + 1] += s + std::to_string(u[l]);
          }
          
          if (i < p_stratum - 1) {
            covariates_aft[j + 1] += ".";
          }
        }
        
        double* z_aftn_col = z_aftn.data_ptr() + j * n;
        for (int r = 0; r < n; ++r) {
          z_aftn_col[r] = stratumn[r] == j ? 1.0 : 0.0;
        }
      }
    }
  }
  
  for (int j = 0; j < p2; ++j) {
    const std::string& zj = base2_cov[j];
    if (!data.containElementNamed(zj))
      throw std::invalid_argument("data must contain the variables in base2_cov");
    if (zj == treat)
      throw std::invalid_argument("treat should be excluded from base2_cov");
    covariates_aft[q + j + 1] = zj;
    double* z_aftn_col = z_aftn.data_ptr() + (q + j) * n;
    if (data.bool_cols.count(zj)) {
      const std::vector<unsigned char>& vb = data.get<unsigned char>(zj);
      for (int i = 0; i < n; ++i) z_aftn_col[i] = vb[i] ? 1.0 : 0.0;
    } else if (data.int_cols.count(zj)) {
      const std::vector<int>& vi = data.get<int>(zj);
      for (int i = 0; i < n; ++i) z_aftn_col[i] = static_cast<double>(vi[i]);
    } else if (data.numeric_cols.count(zj)) {
      const std::vector<double>& vd = data.get<double>(zj);
      std::memcpy(z_aftn_col, vd.data(), n * sizeof(double));
    } else {
      throw std::invalid_argument("base2_cov must be bool, integer or numeric");
    }
  }
  
  std::string dist = aft_dist;
  std::for_each(dist.begin(), dist.end(), [](char & c) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });
  if (dist == "log-logistic" || dist == "llogistic") dist = "loglogistic";
  else if  (dist == "log-normal" || dist == "lnormal") dist = "lognormal";
  if (!(dist == "exponential" || dist == "weibull" || dist == "lognormal" ||
      dist == "loglogistic")) throw std::invalid_argument(
        "aft_dist must be exponential, weibull, lognormal, or loglogistic");
  
  if (alpha <= 0.0 || alpha >= 0.5)
    throw std::invalid_argument("alpha must lie between 0 and 0.5");
  if (ties != "efron" && ties != "breslow")
    throw std::invalid_argument("ties must be efron or breslow");  
  if (offset < 0.0) 
    throw std::invalid_argument("offset must be nonnegative");
  if (n_boot < 100)
    throw std::invalid_argument("n_boot must be greater than or equal to 100");
  
  // exclude observations with missing values
  std::vector<unsigned char> sub(n,1);
  for (int i = 0; i < n; ++i) {
    if (idn[i] == INT_MIN || stratumn[i] == INT_MIN || 
        std::isnan(timen[i]) || std::isnan(eventn[i]) || 
        treatn[i] == INT_MIN || std::isnan(censor_timen[i]) || 
        std::isnan(pdn[i]) || std::isnan(swtrtn[i])) {
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
  subset_in_place(timen, keep);
  subset_in_place(eventn, keep);
  subset_in_place(treatn, keep);
  subset_in_place(censor_timen, keep);
  subset_in_place(pdn, keep);
  subset_in_place(pd_timen, keep);
  subset_in_place(swtrtn, keep);
  subset_in_place(swtrt_timen, keep);
  subset_in_place_flatmatrix(zn, keep);
  subset_in_place_flatmatrix(z_aftn, keep);
  n = static_cast<int>(keep.size());
  
  // summarize number of deaths and switches by treatment arm
  std::vector<int> treat_out = {0, 1};
  std::vector<double> n_total(2);
  std::vector<double> n_event(2);
  std::vector<double> n_pd(2);
  std::vector<double> n_switch(2);
  for (int i = 0; i < n; ++i) {
    int g = treatn[i];
    ++n_total[g];
    if (eventn[i] == 1) ++n_event[g];
    if (pdn[i] == 1) ++n_pd[g];
    if (swtrtn[i] == 1) ++n_switch[g];
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
  
  auto f = [data, has_stratum, stratum, p_stratum, u_stratum, 
            treat, treatwi, treatwn, treatwc, id, idwi, idwn, idwc,
            n, p, qp2, covariates, covariates_aft, dist, 
            recensor, swtrt_control_only, alpha, zcrit, ties, offset](
                const std::vector<int>& idb, 
                const std::vector<int>& stratumb, 
                const std::vector<double>& timeb, 
                const std::vector<double>& eventb, 
                const std::vector<int>& treatb, 
                const std::vector<double>& censor_timeb, 
                const std::vector<double>& pdb, 
                const std::vector<double>& pd_timeb, 
                const std::vector<double>& swtrtb, 
                const std::vector<double>& swtrt_timeb, 
                const FlatMatrix& zb, 
                const FlatMatrix& z_aftb, int k) -> ListCpp {
                  bool fail = false; // whether any model fails to converge
                  std::vector<double> init(1, NaN);
                  
                  // time and event adjusted for treatment switching
                  std::vector<double> t_star = timeb;
                  std::vector<double> d_star = eventb;
                  
                  double psi0hat = NaN, psi0lower = NaN, psi0upper = NaN;
                  double psi1hat = NaN, psi1lower = NaN, psi1upper = NaN;
                  
                  // initialize data_aft and fit_aft
                  std::vector<ListPtr> data_aft(2), fit_aft(2), res_aft(2);
                  if (k == -1) {
                    DataFrameCpp nulldata;
                    ListCpp nullfit;
                    std::vector<double> nullres;
                    for (int h = 0; h < 2; ++h) {
                      ListPtr data_x = std::make_shared<ListCpp>();
                      ListPtr fit_x  = std::make_shared<ListCpp>();
                      ListPtr res_x  = std::make_shared<ListCpp>();
                      data_x->push_back(nulldata, "data");
                      fit_x->push_back(nullfit, "fit");
                      res_x->push_back(nullres, "res");
                      if (data.bool_cols.count(treat) || 
                          data.int_cols.count(treat)) {
                        data_x->push_back(treatwi[1 - h], treat);
                        fit_x->push_back(treatwi[1 - h], treat);
                        res_x->push_back(treatwi[1 - h], treat);
                      } else if (data.numeric_cols.count(treat)) {
                        data_x->push_back(treatwn[1 - h], treat);
                        fit_x->push_back(treatwn[1 - h], treat);
                        res_x->push_back(treatwn[1 - h], treat);
                      } else if (data.string_cols.count(treat)) {
                        data_x->push_back(treatwc[1 - h], treat);
                        fit_x->push_back(treatwc[1 - h], treat);
                        res_x->push_back(treatwc[1 - h], treat);
                      }
                      data_aft[h] = std::move(data_x);
                      fit_aft[h]  = std::move(fit_x);
                      res_aft[h]  = std::move(res_x);
                    }
                  }
                  
                  
                  bool psimissing = false;
                  
                  // # arms that include patients who switched treatment
                  int K = swtrt_control_only ? 1 : 2;
                  for (int h = 0; h < K; ++h) {
                    // post progression data
                    std::vector<int> l; 
                    l.reserve(n);
                    for (int i = 0; i < n; ++i) 
                      if (treatb[i] == h && pdb[i] == 1) l.push_back(i);
                    
                    int m = static_cast<int>(l.size());
                    std::vector<int> id2(m);
                    std::vector<double> time2(m), event2(m), swtrt2(m);
                    for (int i = 0; i < m; ++i) {
                      int j = l[i];
                      id2[i] = idb[j];
                      time2[i] = timeb[j] - pd_timeb[j] + offset;
                      event2[i] = eventb[j];
                      swtrt2[i] = swtrtb[j];
                    }
                    
                    DataFrameCpp data1;
                    data1.push_back(id2, "uid");
                    data1.push_back(std::move(time2), "pps");
                    data1.push_back(std::move(event2), "event");
                    data1.push_back(std::move(swtrt2), "swtrt");
                    
                    for (int j = 0; j < qp2; ++j) {
                      const std::string& zj = covariates_aft[j+1];
                      std::vector<double> u = flatmatrix_get_column(z_aftb, j);
                      data1.push_back(subset(u, l), zj);
                    }
                    
                    ListCpp fit1 = liferegcpp(
                      data1, {""}, "pps", "", "event", covariates_aft, 
                      "", "", "", dist, init, 0, 0, alpha);
                    
                    DataFrameCpp sumstat1 = fit1.get<DataFrameCpp>("sumstat");
                    if (sumstat1.get<unsigned char>("fail")[0]) fail = true;
                    
                    DataFrameCpp parest1 = fit1.get<DataFrameCpp>("parest");
                    std::vector<double> beta1 = parest1.get<double>("beta");
                    double psihat = -beta1[1];
                    double psilower = NaN, psiupper = NaN;
                    if (k == -1) {
                      std::vector<double> sebeta1 = parest1.get<double>("sebeta");
                      psilower = -beta1[1] - zcrit * sebeta1[1];
                      psiupper = -beta1[1] + zcrit * sebeta1[1];
                    }
                    
                    std::vector<double> res1;
                    if (k == -1) {
                      FlatMatrix vbeta1(qp2 + 1, qp2 + 1);
                      FlatMatrix rr = residuals_liferegcpp(
                        beta1, vbeta1, data1, {""}, "pps", "", "event",
                        covariates_aft, "", "", "", dist, "deviance");
                      res1 = flatmatrix_get_column(rr, 0);
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
                    
                    // update data_aft, fit_aft, and res_aft
                    if (k == -1) {
                      std::vector<int> stratum2 = subset(stratumb, l);
                      if (has_stratum) {
                        for (int i = 0; i < p_stratum; ++i) {
                          const std::string& s = stratum[i];
                          if (data.bool_cols.count(s)) {
                            auto v = u_stratum.get<unsigned char>(s);
                            data1.push_back(subset(v, stratum2), s);
                          } else if (data.int_cols.count(s)) {
                            auto v = u_stratum.get<int>(s);
                            data1.push_back(subset(v, stratum2), s);
                          } else if (data.numeric_cols.count(s)) {
                            auto v = u_stratum.get<double>(s);
                            data1.push_back(subset(v, stratum2), s);
                          } else if (data.string_cols.count(s)) {
                            auto v = u_stratum.get<std::string>(s);
                            data1.push_back(subset(v, stratum2), s);
                          }
                        }
                      }
                      
                      if (data.int_cols.count(id)) {
                        data1.push_front(subset(idwi, id2), id);
                      } else if (data.numeric_cols.count(id)) {
                        data1.push_front(subset(idwn, id2), id);
                      } else if (data.string_cols.count(id)) {
                        data1.push_front(subset(idwc, id2), id);
                      }
                      
                      ListPtr& data_x = data_aft[h];
                      data_x->get<DataFrameCpp>("data") = data1;
                      
                      ListPtr& fit_x = fit_aft[h];
                      fit_x->get_list("fit") = fit1;
                      
                      ListPtr& res_x = res_aft[h];
                      res_x->get<std::vector<double>>("res") = res1;
                    }
                    
                    if (!std::isnan(psihat)) {
                      // calculate counter-factual survival times
                      double a = std::exp(psihat);
                      double c0 = std::min(1.0, a);
                      for (int i = 0; i < n; ++i) {
                        if (treatb[i] == h) {
                          double b2, u_star, c_star;
                          if (swtrtb[i] == 1) {
                            b2 = pd_timeb[i] - offset;
                            u_star = b2 + (timeb[i] - b2) * a;
                          } else {
                            u_star = timeb[i];
                          }
                          
                          if (recensor) {
                            c_star = censor_timeb[i] * c0;
                            t_star[i] = std::min(u_star, c_star);
                            d_star[i] = c_star < u_star ? 0 : eventb[i];
                          } else {
                            t_star[i] = u_star;
                            d_star[i] = eventb[i];
                          }
                        }
                      }
                    } else {
                      psimissing = true;
                      fail = true;
                    }
                  }
                  
                  
                  DataFrameCpp data_outcome, km_outcome, lr_outcome;
                  ListCpp fit_outcome;
                  double hrhat = NaN, hrlower = NaN, hrupper = NaN, pvalue = NaN;
                  
                  if (!psimissing) {
                    // Cox model for hypothetical treatment effect estimate
                    data_outcome.push_back(idb, "uid");
                    data_outcome.push_back(std::move(t_star), "t_star");
                    data_outcome.push_back(std::move(d_star), "d_star");
                    data_outcome.push_back(treatb, "treated");
                    data_outcome.push_back(stratumb, "ustratum");
                    
                    for (int j = 0; j < p; ++j) {
                      const std::string& zj = covariates[j+1];
                      std::vector<double> u = flatmatrix_get_column(zb, j);
                      data_outcome.push_back(std::move(u), zj);
                    }
                    
                    // generate KM estimate and log-rank test
                    if (k == -1) {
                      km_outcome = kmestcpp(
                        data_outcome, {"treated"}, "t_star", "", "d_star", 
                        "", "log-log", 1.0 - alpha, 1);
                      
                      lr_outcome = lrtestcpp(
                        data_outcome, {"ustratum"}, "treated", "t_star", "", "d_star");
                    }
                    
                    // fit the outcome model
                    fit_outcome = phregcpp(
                      data_outcome, {"ustratum"}, "t_star", "", "d_star", covariates, 
                      "", "", "", ties, init, 0, 0, 0, 0, 0, alpha);
                    
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
                    out.push_back(std::move(data_aft), "data_aft");
                    out.push_back(std::move(fit_aft), "fit_aft");
                    out.push_back(std::move(res_aft), "res_aft");
                    out.push_back(std::move(data_outcome), "data_outcome");
                    out.push_back(std::move(km_outcome), "km_outcome");
                    out.push_back(std::move(lr_outcome), "lr_outcome");
                    out.push_back(std::move(fit_outcome), "fit_outcome");
                    out.push_back(psi0hat, "psihat");
                    out.push_back(psi0lower, "psilower");
                    out.push_back(psi0upper, "psiupper");
                    out.push_back(psi1hat, "psi1hat");
                    out.push_back(psi1lower, "psi1lower");
                    out.push_back(psi1upper, "psi1upper");
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
  
  ListCpp out = f(idn, stratumn, timen, eventn, treatn, censor_timen,
                  pdn, pd_timen, swtrtn, swtrt_timen, zn, z_aftn, -1);
  
  auto data_aft = out.get<std::vector<ListPtr>>("data_aft");
  auto fit_aft = out.get<std::vector<ListPtr>>("fit_aft");
  auto res_aft = out.get<std::vector<ListPtr>>("res_aft");
  DataFrameCpp data_outcome = out.get<DataFrameCpp>("data_outcome");
  DataFrameCpp km_outcome = out.get<DataFrameCpp>("km_outcome");
  DataFrameCpp lr_outcome = out.get<DataFrameCpp>("lr_outcome");
  ListCpp fit_outcome = out.get_list("fit_outcome");
  double psihat = out.get<double>("psihat");
  double psilower = out.get<double>("psilower");
  double psiupper = out.get<double>("psiupper");
  double psi1hat = out.get<double>("psi1hat");
  double psi1lower = out.get<double>("psi1lower");
  double psi1upper = out.get<double>("psi1upper");
  double hrhat = out.get<double>("hrhat");
  double hrlower = out.get<double>("hrlower");
  double hrupper = out.get<double>("hrupper");
  double pvalue = out.get<double>("pvalue");
  bool fail = out.get<bool>("fail");
  bool psimissing = out.get<bool>("psimissing");
  std::string psi_CI_type = "AFT model";
  
  std::vector<double> hrhats(n_boot), psihats(n_boot), psi1hats(n_boot);
  std::vector<unsigned char> fails(n_boot);
  DataFrameCpp fail_boots_data;
  std::string hr_CI_type;
  
  if (!psimissing) {
    // summarize number of deaths by treatment arm in the outcome data
    std::vector<int> treated = data_outcome.get<int>("treated");
    std::vector<double> event_out = data_outcome.get<double>("d_star");
    std::vector<double> n_event_out(2);
    for (int i = 0; i < n; ++i) {
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
      // sort data by treatment group, stratum and id
      std::vector<int> order = seqcpp(0, n-1);
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
      
      subset_in_place(idn, order);
      subset_in_place(stratumn, order);
      subset_in_place(timen, order);
      subset_in_place(eventn, order);
      subset_in_place(treatn, order);
      subset_in_place(censor_timen, order);
      subset_in_place(pdn, order);
      subset_in_place(pd_timen, order);
      subset_in_place(swtrtn, order);
      subset_in_place(swtrt_timen, order);
      subset_in_place_flatmatrix(zn, order);
      subset_in_place_flatmatrix(z_aftn, order);
      
      std::vector<int> tsx(1,0); // first observation within each treat/stratum
      for (int i = 1; i < n; ++i) {
        if (treatn[i] != treatn[i-1] || stratumn[i] != stratumn[i-1]) {
          tsx.push_back(i);
        }
      }
      
      int ntss = static_cast<int>(tsx.size());
      tsx.push_back(n); // add the end index
      
      // Before running the parallel loop: pre-generate deterministic seeds
      std::vector<uint64_t> seeds(n_boot);
      std::mt19937_64 master_rng(static_cast<uint64_t>(seed)); // user-provided seed
      for (int k = 0; k < n_boot; ++k) seeds[k] = master_rng();
      
      // We'll collect failure bootstrap data per-worker and merge via Worker::join.
      struct BootstrapWorker : public RcppParallel::Worker {
        // references to read-only inputs (no mutation)
        const int n;
        const int ntss;
        const std::vector<int>& tsx;
        const std::vector<int>& idn;
        const std::vector<int>& stratumn;
        const std::vector<double>& timen;
        const std::vector<double>& eventn;
        const std::vector<int>& treatn;
        const std::vector<double>& censor_timen;
        const std::vector<double>& pdn;
        const std::vector<double>& pd_timen;
        const std::vector<double>& swtrtn;
        const std::vector<double>& swtrt_timen;
        const FlatMatrix& zn;
        const FlatMatrix& z_aftn;
        const std::vector<uint64_t>& seeds;
        // function f and other params that f needs are captured from outer scope
        // capture them by reference here so worker can call f(...)
        std::function<ListCpp(const std::vector<int>&, const std::vector<int>&,
                              const std::vector<double>&, const std::vector<double>&,
                              const std::vector<int>&, const std::vector<double>&,
                              const std::vector<double>&, const std::vector<double>&, 
                              const std::vector<double>&, const std::vector<double>&, 
                              const FlatMatrix&, const FlatMatrix&, int)> f;
        
        // result references (each iteration writes unique index into these)
        std::vector<unsigned char>& fails_out;
        std::vector<double>& hrhats_out;
        std::vector<double>& psihats_out;
        std::vector<double>& psi1hats_out;
        
        // Per-worker storage for failed-boot data (to be merged in join)
        std::vector<int> boot_indexc_local;
        std::vector<int> oidc_local;
        std::vector<int> idc_local;
        std::vector<int> stratumc_local;
        std::vector<int> treatc_local;
        std::vector<double> timec_local;
        std::vector<double> eventc_local;
        std::vector<double> censor_timec_local;
        std::vector<double> pdc_local;
        std::vector<double> pd_timec_local;
        std::vector<double> swtrtc_local;
        std::vector<double> swtrt_timec_local;
        // store column-wise z_aftc_local: outer vector length == z_aftn.ncol
        // each inner vector stores column data across failed boots
        std::vector<std::vector<double>> z_aftc_local;
        int index1_local = 0; // number of rows stored so far for z_aftc_local
        
        // constructor
        BootstrapWorker(int n_, int ntss_,
                        const std::vector<int>& tsx_,
                        const std::vector<int>& idn_,
                        const std::vector<int>& stratumn_,
                        const std::vector<double>& timen_,
                        const std::vector<double>& eventn_,
                        const std::vector<int>& treatn_,
                        const std::vector<double>& censor_timen_,
                        const std::vector<double>& pdn_,
                        const std::vector<double>& pd_timen_,
                        const std::vector<double>& swtrtn_,
                        const std::vector<double>& swtrt_timen_,
                        const FlatMatrix& zn_,
                        const FlatMatrix& z_aftn_,
                        const std::vector<uint64_t>& seeds_,
                        decltype(f) f_,
                        std::vector<unsigned char>& fails_out_,
                        std::vector<double>& hrhats_out_,
                        std::vector<double>& psihats_out_,
                        std::vector<double>& psi1hats_out_) :
          n(n_), ntss(ntss_), tsx(tsx_), idn(idn_), 
          stratumn(stratumn_), timen(timen_), eventn(eventn_), treatn(treatn_), 
          censor_timen(censor_timen_), pdn(pdn_), pd_timen(pd_timen_),
          swtrtn(swtrtn_), swtrt_timen(swtrt_timen_), zn(zn_), z_aftn(z_aftn_),
          seeds(seeds_), f(std::move(f_)),
          fails_out(fails_out_), hrhats_out(hrhats_out_), 
          psihats_out(psihats_out_), psi1hats_out(psi1hats_out_) {
          // heuristic reservation to reduce reallocations:
          int ncols_aft = z_aftn.ncol;
          z_aftc_local.resize(static_cast<std::size_t>(ncols_aft));
          // reserve some capacity per column. This is a heuristic; adjust if needed.
          std::size_t per_col_reserve = static_cast<std::size_t>(10 * n);
          for (int col = 0; col < ncols_aft; ++col) 
            z_aftc_local[col].reserve(per_col_reserve);
          // Reserve scalar buffers heuristically (reduce reallocs)
          boot_indexc_local.reserve(per_col_reserve);
          oidc_local.reserve(per_col_reserve);
          idc_local.reserve(per_col_reserve);
          stratumc_local.reserve(per_col_reserve);
          treatc_local.reserve(per_col_reserve);
          timec_local.reserve(per_col_reserve);
          eventc_local.reserve(per_col_reserve);
          censor_timec_local.reserve(per_col_reserve);
          pdc_local.reserve(per_col_reserve);
          pd_timec_local.reserve(per_col_reserve);
          swtrtc_local.reserve(per_col_reserve);
          swtrt_timec_local.reserve(per_col_reserve);
        }
        
        // operator() processes a range of bootstrap iterations [begin, end)
        void operator()(std::size_t begin, std::size_t end) {
          // per-worker reusable buffers (avoid reallocation per iteration)
          std::vector<int> oidb(n), idb(n), stratumb(n), treatb(n);
          std::vector<double> timeb(n), eventb(n), censor_timeb(n);
          std::vector<double> pdb(n), pd_timeb(n), swtrtb(n), swtrt_timeb(n);
          FlatMatrix zb(n, zn.ncol), z_aftb(n, z_aftn.ncol);
          
          for (std::size_t k = begin; k < end; ++k) {
            // deterministic RNG per-iteration
            std::mt19937_64 rng(seeds[k]);
            
            // sample by treatment/stratum blocks
            for (int h = 0; h < ntss; ++h) {
              int start = tsx[h], end = tsx[h + 1];
              int len = end - start;
              std::uniform_int_distribution<int> index_dist(0, len - 1);
              
              std::vector<int> indices(len);
              for (int i = start; i < end; ++i) {
                int j = start + index_dist(rng);
                indices[i - start] = j;
                oidb[i] = idn[j];
                idb[i] = idn[j] + i * n; // make unique ids within bootstrap
                stratumb[i] = stratumn[j];
                timeb[i] = timen[j];
                eventb[i] = eventn[j];
                treatb[i] = treatn[j];
                censor_timeb[i] = censor_timen[j];
                pdb[i] = pdn[j];
                pd_timeb[i] = pd_timen[j];
                swtrtb[i] = swtrtn[j];
                swtrt_timeb[i] = swtrt_timen[j];
              }
              // copy covariates using sampled indices
              for (int l = 0; l < zn.ncol; ++l) {
                const double* zn_col = zn.data_ptr() + l * n;
                double* zb_col = zb.data_ptr() + l * n;
                for (int i = start; i < end; ++i) {
                  int j = indices[i - start];
                  zb_col[i] = zn_col[j];
                }
              }
              for (int l = 0; l < z_aftn.ncol; ++l) {
                const double* z_aftn_col = z_aftn.data_ptr() + l * n;
                double* z_aftb_col = z_aftb.data_ptr() + l * n;
                for (int i = start; i < end; ++i) {
                  int j = indices[i - start];
                  z_aftb_col[i] = z_aftn_col[j];
                }
              }
            } // end block sampling
            
            // call the (thread-safe) per-iteration function f
            ListCpp out = f(idb, stratumb, timeb, eventb, treatb, censor_timeb,
                            pdb, pd_timeb, swtrtb, swtrt_timeb, zb, z_aftb,
                            static_cast<int>(k));
            
            // write results
            fails_out[k] = out.get<bool>("fail");
            hrhats_out[k] = out.get<double>("hrhat");
            psihats_out[k] = out.get<double>("psihat");
            psi1hats_out[k] = out.get<double>("psi1hat");
            
            // existing code that collects failure data when fails_out[k] is true...
            if (fails_out[k]) {
              // collect into boot_indexc_local, idc_local, ... and z_aftc_local
              append(boot_indexc_local, std::vector<int>(n, static_cast<int>(k + 1)));
              append(oidc_local, oidb);
              append(idc_local, idb);
              append(stratumc_local, stratumb);
              append(treatc_local, treatb);
              append(timec_local, timeb);
              append(eventc_local, eventb);
              append(censor_timec_local, censor_timeb);
              append(pdc_local, pdb);
              append(pd_timec_local, pd_timeb);
              append(swtrtc_local, swtrtb);
              append(swtrt_timec_local, swtrt_timeb);
              append_flatmatrix(z_aftc_local, z_aftb);
              index1_local += n;
            }
          } // end for k
        } // end operator()
        
        // join merges other (worker copy) into this instance (called by parallelFor)
        void join(const BootstrapWorker& other) {
          // append other's local vectors into this worker's storage
          append(boot_indexc_local, other.boot_indexc_local);
          append(idc_local, other.idc_local);
          append(stratumc_local, other.stratumc_local);
          append(treatc_local, other.treatc_local);
          append(timec_local, other.timec_local);
          append(eventc_local, other.eventc_local);
          append(censor_timec_local, other.censor_timec_local);
          append(pdc_local, other.pdc_local);
          append(pd_timec_local, other.pd_timec_local);
          append(swtrtc_local, other.swtrtc_local);
          append(swtrt_timec_local, other.swtrt_timec_local);
          append_flatmatrix(z_aftc_local, other.z_aftc_local);
          index1_local += other.index1_local;
        }
      }; // end BootstrapWorker
      
      // Instantiate the Worker with references to inputs and outputs
      BootstrapWorker worker(
          n, ntss, tsx, idn, stratumn, timen, eventn, treatn, censor_timen,
          pdn, pd_timen, swtrtn, swtrt_timen, zn, z_aftn, seeds,
          // bind f into std::function (capture the f we already have)
          std::function<ListCpp(const std::vector<int>&, const std::vector<int>&,
                                const std::vector<double>&, const std::vector<double>&,
                                const std::vector<int>&, const std::vector<double>&,
                                const std::vector<double>&, const std::vector<double>&, 
                                const std::vector<double>&, const std::vector<double>&, 
                                const FlatMatrix&, const FlatMatrix&, int)>(f),
                                fails, hrhats, psihats, psi1hats
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
        fail_boots_data.push_back(std::move(worker.timec_local), "time");
        fail_boots_data.push_back(std::move(worker.eventc_local), "event");
        fail_boots_data.push_back(treatc, "treated");
        fail_boots_data.push_back(std::move(worker.censor_timec_local), "censor_time");
        fail_boots_data.push_back(std::move(worker.pdc_local), "pd");
        fail_boots_data.push_back(std::move(worker.pd_timec_local), "pd_time");
        fail_boots_data.push_back(std::move(worker.swtrtc_local), "swtrt");
        fail_boots_data.push_back(std::move(worker.swtrt_timec_local), "swtrt_time");
        
        int ncols_aft = worker.z_aftc_local.size();
        for (int j = 0; j < ncols_aft; ++j) {
          const std::string& zj = covariates_aft[j+1];
          std::vector<double> u = std::move(worker.z_aftc_local[j]);
          fail_boots_data.push_back(std::move(u), zj);
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
      psihats = worker.psihats_out;
      psi1hats = worker.psi1hats_out;
      
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
      
      if (!swtrt_control_only) {
        std::vector<double> psi1hats1 = subset(psi1hats, ok);
        double meanpsi1, sdpsi1;
        mean_sd(psihats1.data(), n_ok, meanpsi1, sdpsi1);
        psi1lower = psi1hat - tcrit * sdpsi1;
        psi1upper = psi1hat + tcrit * sdpsi1;
      }
    }
  }
  
  ListCpp result;
  std::string pvalue_type = boot ? "bootstrap" : "Cox model";
  std::vector <double> psi_CI = {psilower, psiupper};
  std::vector<double> hr_CI = {hrlower, hrupper};
  result.push_back(psihat, "psi");
  result.push_back(std::move(psi_CI), "psi_CI");
  result.push_back(psi_CI_type, "psi_CI_type");
  result.push_back(pvalue, "pvalue");
  result.push_back(pvalue_type, "pvalue_type");
  result.push_back(hrhat, "hr");
  result.push_back(std::move(hr_CI), "hr_CI");
  result.push_back(hr_CI_type, "hr_CI_type");
  result.push_back(std::move(event_summary), "event_summary");
  result.push_back(std::move(data_aft), "data_aft");
  result.push_back(std::move(fit_aft), "fit_aft");
  result.push_back(std::move(res_aft), "res_aft");
  result.push_back(std::move(data_outcome), "data_outcome");
  result.push_back(std::move(km_outcome), "km_outcome");
  result.push_back(std::move(lr_outcome), "lr_outcome");
  result.push_back(std::move(fit_outcome), "fit_outcome");
  result.push_back(fail, "fail");
  result.push_back(psimissing, "psimissing");
  
  if (!swtrt_control_only) {
    result.push_back(psi1hat, "psi_trt");
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
    if (std::any_of(fails.begin(), fails.end(), [](double x){ return x; })) {
      result.push_back(std::move(fail_boots_data), "fail_boots_data");
    }
  }
  
  thread_utils::drain_thread_warnings_to_R();
  return Rcpp::wrap(result);
}
