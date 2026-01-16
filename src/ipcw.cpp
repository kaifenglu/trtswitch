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
#include <cstdio>
#include <cstdlib>
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
Rcpp::List ipcwcpp(const Rcpp::DataFrame df,
                   const std::string& id,
                   const std::vector<std::string>& stratum,
                   const std::string& tstart,
                   const std::string& tstop,
                   const std::string& event,
                   const std::string& treat,
                   const std::string& swtrt,
                   const std::string& swtrt_time,
                   const std::vector<std::string>& base_cov,
                   const std::vector<std::string>& numerator,
                   const std::vector<std::string>& denominator,
                   const bool logistic_switching_model,
                   const bool strata_main_effect_only,
                   const int ns_df,
                   const bool firth,
                   const bool flic,
                   const bool stabilized_weights,
                   const double trunc,
                   const bool trunc_upper_only,
                   const bool swtrt_control_only,
                   const double alpha,
                   const std::string& ties,
                   const bool boot,
                   const int n_boot,
                   const int seed) {
  
  DataFrameCpp data = convertRDataFrameToCpp(df);
  
  int n = static_cast<int>(data.nrows());
  int p = static_cast<int>(base_cov.size());
  if (p == 1 && base_cov[0] == "") p = 0;
  
  int p1 = static_cast<int>(numerator.size());
  if (p1 == 1 && numerator[0] == "") p1 = 0;
  
  int p2 = static_cast<int>(denominator.size());
  if (p2 == 1 && denominator[0] == "") {
    throw std::invalid_argument("covariates for the switch model must be provided");
  }
  
  if (p1 > 0) {
    if (p == 0) throw std::invalid_argument("numerator must be a subset of base_cov");
    
    std::unordered_set<std::string> base_set;
    base_set.reserve(base_cov.size());
    for (const auto& s : base_cov) base_set.insert(s);
    
    for (const auto& name : numerator) {
      if (base_set.find(name) == base_set.end()) {
        throw std::invalid_argument("numerator must be a subset of base_cov");
      }
    }
  }
  
  if (p > 0) {
    if (p2 == 0) 
      throw std::invalid_argument("base_cov must be a subset of denominator");
    
    std::unordered_set<std::string> denom_set;
    denom_set.reserve(denominator.size());
    for (const auto& s : denominator) denom_set.insert(s);
    
    for (const auto& name : base_cov) {
      if (denom_set.find(name) == denom_set.end()) {
        throw std::invalid_argument("base_cov must be a subset of denominator");
      }
    }
  }
  
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
  
  // covariates for the logistic regression switch model for denominator
  // including stratum, denominator, and ns_df spline terms
  std::vector<std::string> covariates_lgs_den(q + p2 + ns_df);
  FlatMatrix z_lgs_denn(n, q + p2);
  if (has_stratum) {
    if (strata_main_effect_only) {
      int k = 0;
      for (int i = 0; i < p_stratum; ++i) {
        const std::string& s = stratum[i];
        int di = d[i] - 1;
        
        if (u_stratum.string_cols.count(s)) {
          auto u = levels.get<std::vector<std::string>>(s);
          for (int j = 0; j < di; ++j) {
            covariates_lgs_den[k + j] = s + sanitize(u[j]);
          }
        } else if (u_stratum.numeric_cols.count(s)) {
          auto u = levels.get<std::vector<double>>(s);
          for (int j = 0; j < di; ++j) {
            covariates_lgs_den[k + j] = s + std::to_string(u[j]);
          }
        } else if (u_stratum.int_cols.count(s)) {
          auto u = levels.get<std::vector<int>>(s);
          for (int j = 0; j < di; ++j) {
            covariates_lgs_den[k + j] = s + std::to_string(u[j]);
          }
        } else if (u_stratum.bool_cols.count(s)) {
          auto u = levels.get<std::vector<unsigned char>>(s);
          for (int j = 0; j < di; ++j) {
            covariates_lgs_den[k + j] = s + std::to_string(u[j]);
          }
        }
        
        for (int j = 0; j < di; ++j) {
          const int* stratan_col = stratan.data_ptr() + i * n;
          double* z_lgs_denn_col = z_lgs_denn.data_ptr() + (k + j) * n;
          for (int r = 0; r < n; ++r) {
            z_lgs_denn_col[r] = stratan_col[r] == j ? 1.0 : 0.0;
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
        
        covariates_lgs_den[j] = "";
        
        for (int i = 0; i < p_stratum; ++i) {
          const std::string& s = stratum[i];
          
          std::vector<int> q_col = intmatrix_get_column(stratan, i);
          int l = q_col[first_k];
          
          if (u_stratum.string_cols.count(s)) {
            auto u = levels.get<std::vector<std::string>>(s);
            covariates_lgs_den[j] += s + sanitize(u[l]);
          } else if (u_stratum.numeric_cols.count(s)) {
            auto u = levels.get<std::vector<double>>(s);
            covariates_lgs_den[j] += s + std::to_string(u[l]);
          } else if (u_stratum.int_cols.count(s)) {
            auto u = levels.get<std::vector<int>>(s);
            covariates_lgs_den[j] += s + std::to_string(u[l]);
          } else if (u_stratum.bool_cols.count(s)) {
            auto u = levels.get<std::vector<unsigned char>>(s);
            covariates_lgs_den[j] += s + std::to_string(u[l]);
          }
          
          if (i < p_stratum - 1) {
            covariates_lgs_den[j] += ".";
          }
        }
        
        double* z_lgs_denn_col = z_lgs_denn.data_ptr() + j * n;
        for (int r = 0; r < n; ++r) {
          z_lgs_denn_col[r] = stratumn[r] == j ? 1.0 : 0.0;
        }
      }
    }
  }
  
  FlatMatrix z_cox_denn(n, p2);
  for (int j = 0; j < p2; ++j) {
    const std::string& zj = denominator[j];
    if (!data.containElementNamed(zj))
      throw std::invalid_argument("data must contain the variables in denominator");
    if (zj == treat)
      throw std::invalid_argument("treat should be excluded from denominator");
    covariates_lgs_den[q + j] = zj;
    double* z_lgs_denn_col = z_lgs_denn.data_ptr() + (q + j) * n;
    double* z_cox_denn_col = z_cox_denn.data_ptr() + j * n;
    if (data.bool_cols.count(zj)) {
      const std::vector<unsigned char>& vb = data.get<unsigned char>(zj);
      for (int i = 0; i < n; ++i) 
        z_lgs_denn_col[i] = z_cox_denn_col[i] = vb[i] ? 1.0 : 0.0;
    } else if (data.int_cols.count(zj)) {
      const std::vector<int>& vi = data.get<int>(zj);
      for (int i = 0; i < n; ++i) 
        z_lgs_denn_col[i] = z_cox_denn_col[i] = static_cast<double>(vi[i]);
    } else if (data.numeric_cols.count(zj)) {
      const std::vector<double>& vd = data.get<double>(zj);
      std::memcpy(z_lgs_denn_col, vd.data(), n * sizeof(double));
      std::memcpy(z_cox_denn_col, vd.data(), n * sizeof(double));
    } else {
      throw std::invalid_argument("covariates must be bool, integer or numeric");
    }
  }
  
  if (ns_df < 0) {
    throw std::invalid_argument("ns_df must be a nonnegative integer");
  }
  
  for (int j = 0; j < ns_df; ++j) {
    covariates_lgs_den[q + p2 + j] = "ns" + std::to_string(j+1);
  }
  
  // covariates for the logistic regression switch model for numerator
  // including stratum, numerator, and ns_df spline terms
  std::vector<std::string> covariates_lgs_num(q + p1 + ns_df);
  for (int j = 0; j < q; ++j) {
    covariates_lgs_num[j] = covariates_lgs_den[j];
  }
  for (int j = 0; j < p1; ++j) {
    covariates_lgs_num[q + j] = numerator[j];
  }
  for (int j=0; j<ns_df; ++j) {
    covariates_lgs_num[q + p1 + j] = "ns" + std::to_string(j+1);
  }
  
  if (trunc < 0.0 || trunc >= 0.5) {
    throw std::invalid_argument("trunc must lie in [0, 0.5)");
  }
  
  if (alpha <= 0.0 || alpha >= 0.5) {
    throw std::invalid_argument("alpha must lie between 0 and 0.5");
  }
  
  if (ties != "efron" && ties != "breslow") {
    throw std::invalid_argument("ties must be efron or breslow");
  }
  
  if (n_boot < 100) {
    throw std::invalid_argument("n_boot must be greater than or equal to 100");
  }
  
  // exclude observations with missing values
  std::vector<unsigned char> sub(n,1);
  for (int i = 0; i < n; ++i) {
    if (idn[i] == INT_MIN || stratumn[i] == INT_MIN || 
        std::isnan(tstartn[i]) || std::isnan(tstopn[i]) || 
        eventn[i] == INT_MIN || treatn[i] == INT_MIN || 
        swtrtn[i] == INT_MIN) {
      sub[i] = 0; continue;
    }
    for (int j = 0; j < q + p2; ++j) {
      if (std::isnan(z_lgs_denn(i,j))) { sub[i] = 0; break; }
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
  subset_in_place(swtrtn, keep);
  subset_in_place(swtrt_timen, keep);
  subset_in_place_flatmatrix(zn, keep);
  subset_in_place_flatmatrix(z_cox_denn, keep);
  subset_in_place_flatmatrix(z_lgs_denn, keep);
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
      swtrtn.push_back(swtrtn[l]);
      swtrt_timen.push_back(swtrt_timen[l]);
      
      // change tstop and event for the old observation
      tstopn[l] = swtrt_timen[l];
      eventn[l] = 0;
    }
    
    // append new rows to the covariate matrices
    FlatMatrix zn_new = subset_flatmatrix(zn, sub);
    FlatMatrix z_cox_denn_new = subset_flatmatrix(z_cox_denn, sub);
    FlatMatrix z_lgs_denn_new = subset_flatmatrix(z_lgs_denn, sub);
    zn = concat_flatmatrix(zn, zn_new);
    z_cox_denn = concat_flatmatrix(z_cox_denn, z_cox_denn_new);
    z_lgs_denn = concat_flatmatrix(z_lgs_denn, z_lgs_denn_new);
    
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
  subset_in_place(swtrtn, order);
  subset_in_place(swtrt_timen, order);
  subset_in_place_flatmatrix(zn, order);
  subset_in_place_flatmatrix(z_cox_denn, order);
  subset_in_place_flatmatrix(z_lgs_denn, order);
  
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
    if (swtrtn[i] == 1 && swtrt_timen[i] > os_timen[i]) {
      throw std::invalid_argument("swtrt_time must be less than or equal to os_time");
    }
  }
  
  // subset to one observation per id for event summary
  std::vector<int> treatn1 = subset(treatn, idx1);
  std::vector<int> eventn1 = subset(eventn, idx1);
  std::vector<int> swtrtn1 = subset(swtrtn, idx1);
  
  // summarize number of deaths and switches by treatment arm
  std::vector<int> treat_out = {0, 1};
  std::vector<double> n_total(2);
  std::vector<double> n_event(2);
  std::vector<double> n_switch(2);
  for (int i = 0; i < nids; ++i) {
    int g = treatn1[i];
    ++n_total[g];
    if (eventn1[i] == 1) ++n_event[g];
    if (swtrtn1[i] == 1) ++n_switch[g];
  }
  
  // Compute percentages
  std::vector<double> pct_event(2);
  std::vector<double> pct_switch(2);
  for (int g = 0; g < 2; g++) {
    pct_event[g] = 100.0 * n_event[g] / n_total[g];
    pct_switch[g] = 100.0 * n_switch[g] / n_total[g];
  }
  
  // Combine count and percentage
  DataFrameCpp event_summary;
  event_summary.push_back(std::move(treat_out), "treated");
  event_summary.push_back(n_total, "n");
  event_summary.push_back(std::move(n_event), "event_n");
  event_summary.push_back(std::move(pct_event), "event_pct");
  event_summary.push_back(std::move(n_switch), "switch_n");
  event_summary.push_back(std::move(pct_switch), "switch_pct");
  
  double zcrit = boost_qnorm(1.0 - alpha / 2.0);
  
  auto f = [data, has_stratum, stratum, p_stratum, u_stratum, 
            treat, treatwi, treatwn, treatwc, id, idwi, idwn, idwc,
            p, p1, p2, q, covariates, numerator, denominator, 
            covariates_lgs_num, covariates_lgs_den, 
            logistic_switching_model, ns_df, firth, flic, 
            stabilized_weights, trunc, trunc_upper_only, 
            swtrt_control_only, alpha, zcrit, ties](
                const std::vector<int>& idb, 
                const std::vector<int>& stratumb, 
                const std::vector<double>& tstartb,
                const std::vector<double>& tstopb,
                const std::vector<int>& eventb, 
                const std::vector<int>& treatb, 
                const std::vector<double>& os_timeb, 
                const std::vector<int>& swtrtb, 
                const std::vector<double>& swtrt_timeb, 
                const FlatMatrix& zb,
                const FlatMatrix& z_cox_denb,
                const FlatMatrix& z_lgs_denb, int k) -> ListCpp {
                  // the total number of rows change across bootstrap samples 
                  // because bootstrap is done at the subject level and 
                  // different subjects generally have different number of rows
                  int n = static_cast<int>(idb.size());
                  bool fail = false; // whether any model fails to converge
                  std::vector<double> init(1, NaN);
                  
                  int n1;
                  std::vector<int> id1, stratum1, treat1;
                  std::vector<int> event1, swtrt1, cross1;
                  std::vector<double> tstart1, tstop1, os_time1, swtrt_time1;
                  FlatMatrix z1, z_cox_den1, z_lgs_den1;
                  if (!swtrt_control_only) {
                    // exclude observations after switch
                    std::vector<int> l; 
                    l.reserve(n);
                    for (int i = 0; i < n; ++i) {
                      if (swtrtb[i] == 1 && tstartb[i] >= swtrt_timeb[i]) {
                        continue;
                      }
                      l.push_back(i);
                    }
                    
                    id1 = subset(idb, l);
                    stratum1 = subset(stratumb, l);
                    tstart1 = subset(tstartb, l);
                    tstop1 = subset(tstopb, l);
                    event1 = subset(eventb, l);
                    treat1 = subset(treatb, l);
                    os_time1 = subset(os_timeb, l);
                    swtrt1 = subset(swtrtb, l);
                    swtrt_time1 = subset(swtrt_timeb, l);
                    z1 = subset_flatmatrix(zb, l);
                    z_cox_den1 = subset_flatmatrix(z_cox_denb, l);
                    z_lgs_den1 = subset_flatmatrix(z_lgs_denb, l);
                    n1 = static_cast<int>(l.size());
                    
                    // set up crossover indicators
                    cross1 = std::vector<int>(n1);
                    for (int i = 0; i < n1; ++i) {
                      if (i == n1 - 1 || id1[i] != id1[i+1]) {
                        if (swtrt1[i] == 1 && tstop1[i] >= swtrt_time1[i]) {
                          cross1[i] = 1;
                        }
                      }
                    }
                  } else {
                    // for control group, exclude observations after switch
                    std::vector<int> l;
                    l.reserve(n);
                    for (int i = 0; i < n; ++i) {
                      if (treatb[i] == 1) continue;
                      if (swtrtb[i] == 1 && tstartb[i] >= swtrt_timeb[i]) {
                        continue;
                      }
                      l.push_back(i);
                    }
                    
                    std::vector<int> id10 = subset(idb, l);
                    std::vector<int> stratum10 = subset(stratumb, l);
                    std::vector<double> tstart10 = subset(tstartb, l);
                    std::vector<double> tstop10 = subset(tstopb, l);
                    std::vector<int> event10 = subset(eventb, l);
                    std::vector<int> treat10 = subset(treatb, l);
                    std::vector<double> os_time10 = subset(os_timeb, l);
                    std::vector<int> swtrt10 = subset(swtrtb, l);
                    std::vector<double> swtrt_time10 = subset(swtrt_timeb, l);
                    FlatMatrix z10 = subset_flatmatrix(zb, l);
                    FlatMatrix z_cox_den10 = subset_flatmatrix(z_cox_denb, l);
                    FlatMatrix z_lgs_den10 = subset_flatmatrix(z_lgs_denb, l);
                    int n10 = static_cast<int>(l.size());
                    
                    // set up crossover indicators for control
                    std::vector<int> cross10(n10);
                    for (int i = 0; i < n10; ++i) {
                      if (i == n10 - 1 || id10[i] != id10[i+1]) {
                        if (swtrt10[i]== 1 && tstop10[i] >= swtrt_time10[i]) {
                          cross10[i] = 1;
                        }
                      }
                    }
                    
                    // extract data for the active group
                    int start = 0;
                    for (; start < n; ++start) {
                      if (treatb[start] == 1) break;
                    }
                    
                    std::vector<int> id11 = subset(idb, start, n);
                    std::vector<int> stratum11 = subset(stratumb, start, n);
                    std::vector<double> tstart11 = subset(tstartb, start, n);
                    std::vector<double> tstop11 = subset(tstopb, start, n);
                    std::vector<int> event11 = subset(eventb, start, n);
                    std::vector<int> treat11 = subset(treatb, start, n);
                    std::vector<double> os_time11 = subset(os_timeb, start, n);
                    std::vector<int> swtrt11 = subset(swtrtb, start, n);
                    std::vector<double> swtrt_time11 = subset(swtrt_timeb, start, n);
                    FlatMatrix z11 = subset_flatmatrix(zb, start, n);
                    FlatMatrix z_cox_den11 = subset_flatmatrix(z_cox_denb, start, n);
                    FlatMatrix z_lgs_den11 = subset_flatmatrix(z_lgs_denb, start, n);
                    int n11 = n - start;
                    
                    // no crossover in active group
                    std::vector<int> cross11(n11);
                    
                    // combine control and active group data
                    id1 = concat(id10, id11);
                    stratum1 = concat(stratum10, stratum11);
                    tstart1 = concat(tstart10, tstart11);
                    tstop1 = concat(tstop10, tstop11);
                    event1 = concat(event10, event11);
                    treat1 = concat(treat10, treat11);
                    os_time1 = concat(os_time10, os_time11);
                    swtrt1 = concat(swtrt10, swtrt11);
                    swtrt_time1 = concat(swtrt_time10, swtrt_time11);
                    cross1 = concat(cross10, cross11);
                    z1 = concat_flatmatrix(z10, z11);
                    z_cox_den1 = concat_flatmatrix(z_cox_den10, z_cox_den11);
                    z_lgs_den1 = concat_flatmatrix(z_lgs_den10, z_lgs_den11);
                    n1 = n10 + n11;
                  }
                  
                  // initialize data_switch and fit_switch
                  std::vector<ListPtr> data_switch(2), fit_switch(2);
                  if (k == -1) {
                    DataFrameCpp nulldata;
                    ListCpp nullfit;
                    for (int h = 0; h < 2; ++h) {
                      ListPtr data_x = std::make_shared<ListCpp>();
                      ListPtr fit_x  = std::make_shared<ListCpp>();
                      data_x->push_back(nulldata, "data");
                      fit_x->push_back(nullfit, "fit_den");
                      fit_x->push_back(nullfit, "fit_num");
                      if (data.bool_cols.count(treat) || 
                          data.int_cols.count(treat)) {
                        data_x->push_back(treatwi[1 - h], treat);
                        fit_x->push_back(treatwi[1 - h], treat);
                      } else if (data.numeric_cols.count(treat)) {
                        data_x->push_back(treatwn[1 - h], treat);
                        fit_x->push_back(treatwn[1 - h], treat);
                      } else if (data.string_cols.count(treat)) {
                        data_x->push_back(treatwc[1 - h], treat);
                        fit_x->push_back(treatwc[1 - h], treat);
                      }
                      data_switch[h] = std::move(data_x);
                      fit_switch[h]  = std::move(fit_x);
                    }
                  }
                  
                  DataFrameCpp data_outcome;
                  
                  // # arms that include patients who switched treatment
                  int K = swtrt_control_only ? 1 : 2;
                  std::vector<int> w_treated(K), w_n(K);
                  std::vector<double> w_min(K), w_Q1(K), w_median(K);
                  std::vector<double> w_mean(K), w_Q3(K), w_max(K);
                  
                  if (!logistic_switching_model) { // time-dependent Cox
                    // obtain the unique event times across treatment groups
                    
                    // unique event times
                    std::vector<double> cut;
                    cut.reserve(n1);
                    for (int i = 0; i < n1; ++i) {
                      if (event1[i] == 1) cut.push_back(tstop1[i]);
                    }
                    cut = unique_sorted(cut);
                    
                    // replicate event times within each subject
                    std::vector<int> id2, stratum2, event2, treat2, cross2;
                    std::vector<double> tstart2, tstop2;
                    FlatMatrix z2, z_cox_den2;
                    int n2;
                    if (!swtrt_control_only) {
                      DataFrameCpp a = survsplitcpp(tstart1, tstop1, cut);
                      std::vector<int> censor = a.get<int>("censor");
                      std::vector<int> l = a.get<int>("row");
                      id2 = subset(id1, l);
                      stratum2 = subset(stratum1, l);
                      tstart2 = a.get<double>("start");
                      tstop2 = a.get<double>("end");
                      event2 = subset(event1, l);
                      treat2 = subset(treat1, l);
                      cross2 = subset(cross1, l);
                      z2 = subset_flatmatrix(z1, l);
                      z_cox_den2 = subset_flatmatrix(z_cox_den1, l);
                      n2 = static_cast<int>(l.size());
                      for (int i = 0; i < n2; ++i) {
                        if (censor[i] == 1) {
                          event2[i] = 0;
                          cross2[i] = 0;
                        }
                      }
                    } else {
                      // extract data for the control group
                      int end = 0;
                      for (; end < n1; ++end) {
                        if (treat1[end] == 1) break;
                      }
                      std::vector<int> id10 = subset(id1, 0, end);
                      std::vector<int> stratum10 = subset(stratum1, 0, end);
                      std::vector<double> tstart10 = subset(tstart1, 0, end);
                      std::vector<double> tstop10 = subset(tstop1, 0, end);
                      std::vector<int> event10 = subset(event1, 0, end);
                      std::vector<int> treat10 = subset(treat1, 0, end);
                      std::vector<int> cross10 = subset(cross1, 0, end);
                      FlatMatrix z10 = subset_flatmatrix(z1, 0, end);
                      FlatMatrix z_cox_den10 = subset_flatmatrix(z_cox_den1, 0, end);
                      
                      // replicate event times within each subject
                      DataFrameCpp a = survsplitcpp(tstart10, tstop10, cut);
                      std::vector<int> censor = a.get<int>("censor");
                      std::vector<int> l = a.get<int>("row");
                      std::vector<int> id20 = subset(id10, l);
                      std::vector<int> stratum20 = subset(stratum10, l);
                      std::vector<double> tstart20 = a.get<double>("start");
                      std::vector<double> tstop20 = a.get<double>("end");
                      std::vector<int> event20 = subset(event10, l);
                      std::vector<int> treat20 = subset(treat10, l);
                      std::vector<int> cross20 = subset(cross10, l);
                      FlatMatrix z20 = subset_flatmatrix(z10, l);
                      FlatMatrix z_cox_den20 = subset_flatmatrix(z_cox_den10, l);
                      int n20 = static_cast<int>(l.size());
                      for (int i = 0; i < n20; ++i) {
                        if (censor[i] == 1) {
                          event20[i] = 0;
                          cross20[i] = 0;
                        }
                      }
                      
                      // extract data for the active group
                      int start = end;
                      std::vector<int> id21 = subset(id1, start, n1);
                      std::vector<int> stratum21 = subset(stratum1, start, n1);
                      std::vector<double> tstart21 = subset(tstart1, start, n1);
                      std::vector<double> tstop21 = subset(tstop1, start, n1);
                      std::vector<int> event21 = subset(event1, start, n1);
                      std::vector<int> treat21 = subset(treat1, start, n1);
                      std::vector<int> cross21 = subset(cross1, start, n1);
                      FlatMatrix z21 = subset_flatmatrix(z1, start, n1);
                      FlatMatrix z_cox_den21 = 
                        subset_flatmatrix(z_cox_den1, start, n1);
                      int n21 = n1 - start;
                      
                      // combine weighted control with unweighted active data
                      id2 = concat(id20, id21);
                      stratum2 = concat(stratum20, stratum21);
                      tstart2 = concat(tstart20, tstart21);
                      tstop2 = concat(tstop20, tstop21);
                      event2 = concat(event20, event21);
                      treat2 = concat(treat20, treat21);
                      cross2 = concat(cross20, cross21);
                      z2 = concat_flatmatrix(z20, z21);
                      z_cox_den2 = concat_flatmatrix(z_cox_den20, z_cox_den21);
                      n2 = n20 + n21;
                    }
                    
                    // initialize weights
                    std::vector<double> w2(n2, 1.0), sw2(n2, 1.0);
                    
                    // fit the switching models by treatment group
                    for (int h = 0; h < K; ++h) {
                      int mid = 0;
                      for (; mid < n2; ++mid) {
                        if (treat2[mid] == 1) break;
                      }
                      
                      int start, end;
                      if (h == 0) {
                        start = 0; end = mid;
                      } else {
                        start = mid; end = n2;
                      }
                      
                      std::vector<int> id3 = subset(id2, start, end);
                      std::vector<int> stratum3 = subset(stratum2, start, end);
                      std::vector<double> tstart3 = subset(tstart2, start, end);
                      std::vector<double> tstop3 = subset(tstop2, start, end);
                      std::vector<int> cross3 = subset(cross2, start, end);
                      FlatMatrix z_cox_den3 = 
                        subset_flatmatrix(z_cox_den2, start, end);
                      int n3 = end - start;
                      
                      // prepare the data for fitting the switching model
                      DataFrameCpp data1;
                      data1.push_back(id3, "uid");
                      data1.push_back(std::move(stratum3), "ustratum");
                      data1.push_back(std::move(tstart3), "tstart");
                      data1.push_back(tstop3, "tstop");
                      data1.push_back(std::move(cross3), "cross");
                      
                      for (int j = 0; j < p2; ++j) {
                        const std::string& zj = denominator[j];
                        std::vector<double> u = flatmatrix_get_column(z_cox_den3, j);
                        data1.push_back(std::move(u), zj);
                      }
                      
                      // fit the denominator model for crossover
                      ListCpp fit_den = phregcpp(
                        data1, {"ustratum"}, "tstart", "tstop", "cross", 
                        denominator, "", "", "uid", ties, init, 0, 1, 0, 0, 0, alpha);
                      
                      // obtain the survival probabilities for crossover
                      DataFrameCpp sumstat_den = fit_den.get<DataFrameCpp>("sumstat");
                      if (sumstat_den.get<unsigned char>("fail")[0]) fail = true;
                      
                      DataFrameCpp parest_den = fit_den.get<DataFrameCpp>("parest");
                      std::vector<double> beta_den = parest_den.get<double>("beta");
                      FlatMatrix vbeta_den = fit_den.get<FlatMatrix>("vbeta");
                      DataFrameCpp basehaz_den = fit_den.get<DataFrameCpp>("basehaz");
                      
                      DataFrameCpp km_den = survfit_phregcpp(
                        p2, beta_den, vbeta_den, basehaz_den, data1, 
                        denominator, {"ustratum"}, "", "uid", "tstart", 
                        "tstop", 0, "log-log", 1.0 - alpha);
                      std::vector<double> surv_den = km_den.get<double>("surv");
                      
                      ListCpp fit_num = phregcpp(
                        data1, {"ustratum"}, "tstart", "tstop", "cross", 
                        numerator, "", "", "uid", ties, init, 0, 1, 0, 0, 0, alpha);
                      
                      int p10 = (p1 > 0) ? p1 : 1;
                      std::vector<double> beta_num(p10);
                      FlatMatrix vbeta_num(p10,p10);
                      if (p1 > 0) {
                        DataFrameCpp sumstat_num = 
                          fit_num.get<DataFrameCpp>("sumstat");
                        if (sumstat_num.get<unsigned char>("fail")[0]) fail = true;
                        
                        DataFrameCpp parest_num = fit_num.get<DataFrameCpp>("parest");
                        beta_num = parest_num.get<double>("beta");
                        vbeta_num = fit_num.get<FlatMatrix>("vbeta");
                      }
                      DataFrameCpp basehaz_num = fit_num.get<DataFrameCpp>("basehaz");
                      
                      DataFrameCpp km_num = survfit_phregcpp(
                        p1, beta_num, vbeta_num, basehaz_num, data1, 
                        numerator, {"ustratum"}, "", "uid", "tstart", 
                        "tstop", 0, "log-log", 1.0 - alpha);
                      std::vector<double> surv_num = km_num.get<double>("surv");
                      
                      // update data_switch and fit_switch
                      if (k == -1) {
                        std::vector<int> uid = data1.get<int>("uid");
                        if (data.int_cols.count(id)) {
                          data1.push_front(subset(idwi, uid), id);
                        } else if (data.numeric_cols.count(id)) {
                          data1.push_front(subset(idwn, uid), id);
                        } else if (data.string_cols.count(id)) {
                          data1.push_front(subset(idwc, uid), id);
                        }
                        
                        if (has_stratum) {
                          std::vector<int> ustratum = data1.get<int>("ustratum");
                          for (int i = 0; i < p_stratum; ++i) {
                            const std::string& s = stratum[i];
                            if (data.bool_cols.count(s)) {
                              auto v = u_stratum.get<unsigned char>(s);
                              data1.push_back(subset(v, ustratum), s);
                            } else if (data.int_cols.count(s)) {
                              auto v = u_stratum.get<int>(s);
                              data1.push_back(subset(v, ustratum), s);
                            } else if (data.numeric_cols.count(s)) {
                              auto v = u_stratum.get<double>(s);
                              data1.push_back(subset(v, ustratum), s);
                            } else if (data.string_cols.count(s)) {
                              auto v = u_stratum.get<std::string>(s);
                              data1.push_back(subset(v, ustratum), s);
                            }
                          }
                        }
                        
                        ListPtr& data_x = data_switch[h];
                        data_x->get<DataFrameCpp>("data") = data1;
                        
                        ListPtr& fit_x = fit_switch[h];
                        fit_x->get_list("fit_den") = fit_den;
                        fit_x->get_list("fit_num") = fit_num;
                      }
                      
                      // unstabilized and stabilized weights
                      int m3 = surv_den.size();
                      std::vector<double> w(m3), sw(m3);
                      for (int i = 0; i < m3; ++i) {
                        if (surv_den[i] == 0.0) {
                          w[i] = NaN;
                          sw[i] = NaN;
                        } else {
                          w[i] = 1.0 / surv_den[i];
                          sw[i] = surv_num[i] / surv_den[i];
                        }
                      }
                      
                      // weights in the outcome model by matching id and time
                      std::vector<int> id4 = km_den.get<int>("uid");
                      std::vector<double> time4 = km_den.get<double>("time");
                      std::vector<int> sub = match3(id3, tstop3, id4, time4);
                      std::vector<double> w3 = subset(w, sub);
                      std::vector<double> sw3 = subset(sw, sub);
                      
                      // fill in missing weights with LOCF starting at 1
                      std::vector<int> idx3(1, 0);
                      for (int i = 1; i < n3; ++i) {
                        if (id3[i] != id3[i-1]) {
                          idx3.push_back(i);
                        }
                      }
                      
                      int nids3 = static_cast<int>(idx3.size());
                      idx3.push_back(n3);
                      
                      for (int i = 0; i < nids3; ++i) {
                        int j1 = idx3[i], j2 = idx3[i+1];
                        if (std::isnan(w3[j1])) {
                          w3[j1] = 1.0;
                          sw3[j1] = 1.0;
                        }
                        for (int j = j1 + 1; j < j2; ++j) {
                          if (std::isnan(w3[j])) {
                            w3[j] = w3[j-1];
                            sw3[j] = sw3[j-1];
                          }
                        }
                      }
                      
                      // truncate the weights if requested
                      if (trunc > 0.0) {
                        truncate_in_place(w3, trunc_upper_only, trunc);
                        truncate_in_place(sw3, trunc_upper_only, trunc);
                      }
                      
                      // summarize weights for the treatment arm
                      if (k == -1) {
                        w_treated[h] = h;
                        w_n[h] = n3;
                        if (stabilized_weights) {
                          w_min[h] = *std::min_element(sw3.begin(), sw3.end());
                          w_Q1[h] = quantilecpp(sw3, 0.25);
                          w_median[h] = quantilecpp(sw3, 0.5);
                          w_mean[h] = mean_kahan(sw3);
                          w_Q3[h] = quantilecpp(sw3, 0.75);
                          w_max[h] = *std::max_element(sw3.begin(), sw3.end());
                        } else {
                          w_min[h] = *std::min_element(w3.begin(), w3.end());
                          w_Q1[h] = quantilecpp(w3, 0.25);
                          w_median[h] = quantilecpp(w3, 0.5);
                          w_mean[h] = mean_kahan(w3);
                          w_Q3[h] = quantilecpp(w3, 0.75);
                          w_max[h] = *std::max_element(w3.begin(), w3.end());
                        }
                      }
                      
                      // fill in the weights
                      std::memcpy(w2.data() + start, w3.data(), n3 * sizeof(double));
                      std::memcpy(sw2.data() + start, sw3.data(), n3 * sizeof(double));
                    }
                    
                    // prepare data for the outcome model
                    data_outcome.push_back(std::move(id2), "uid");
                    data_outcome.push_back(std::move(tstart2), "tstart");
                    data_outcome.push_back(std::move(tstop2), "tstop");
                    data_outcome.push_back(std::move(event2), "event");
                    data_outcome.push_back(std::move(treat2), "treated");
                    data_outcome.push_back(std::move(w2), "unstabilized_weight");
                    data_outcome.push_back(std::move(sw2), "stabilized_weight");
                    data_outcome.push_back(std::move(stratum2), "ustratum");
                    
                    for (int j = 0; j < p; ++j) {
                      const std::string& zj = covariates[j+1];
                      std::vector<double> u = flatmatrix_get_column(z2, j);
                      data_outcome.push_back(std::move(u), zj);
                    }
                  } else { // logistic regression switching model
                    std::vector<double> w1(n1, 1.0), sw1(n1, 1.0);
                    
                    // fit the switching models by treatment group
                    for (int h = 0; h < K; ++h) {
                      std::vector<int> l;
                      l.reserve(n1);
                      for (int i = 0; i < n1; ++i) {
                        if (treat1[i] != h) {
                          continue; 
                        }
                        if (swtrt1[i] == 1 && tstart1[i] >= swtrt_time1[i]) {
                          continue;
                        }
                        if (swtrt1[i] == 0 && tstop1[i] >= os_time1[i]) {
                          continue;
                        }
                        l.push_back(i);
                      }
                      
                      std::vector<int> id2 = subset(id1, l);
                      std::vector<int> stratum2 = subset(stratum1, l);
                      std::vector<double> tstart2 = subset(tstart1, l);
                      std::vector<double> tstop2 = subset(tstop1, l);
                      std::vector<int> cross2 = subset(cross1, l);
                      FlatMatrix z_lgs_den2 = subset_flatmatrix(z_lgs_den1, l);
                      int n2 = static_cast<int>(l.size());
                      
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
                      
                      // prepare the data for fitting the switching model
                      DataFrameCpp data1;
                      data1.push_back(std::move(id2), "uid");
                      data1.push_back(std::move(stratum2), "ustratum");
                      data1.push_back(std::move(tstart2), "tstart");
                      data1.push_back(std::move(tstop2), "tstop");
                      data1.push_back(cross2, "cross");
                      
                      for (int j = 0; j < q + p2; ++j) {
                        const std::string& zj = covariates_lgs_den[j];
                        std::vector<double> u = flatmatrix_get_column(z_lgs_den2, j);
                        data1.push_back(std::move(u), zj);
                      }
                      for (int j = 0; j < ns_df; ++j) {
                        const std::string& zj = covariates_lgs_den[q + p2 + j];
                        std::vector<double> u = flatmatrix_get_column(s, j);
                        data1.push_back(std::move(u), zj);
                      }
                      
                      ListCpp fit_den = logisregcpp(
                        data1, "cross", covariates_lgs_den, "", "", "", 
                        "uid", "logit", init, 0, firth, flic, 0, alpha);
                      DataFrameCpp sumstat_den = fit_den.get<DataFrameCpp>("sumstat");
                      if (sumstat_den.get<unsigned char>("fail")[0]) fail = true;
                      DataFrameCpp f_den = fit_den.get<DataFrameCpp>("fitted");
                      std::vector<double> h_den = f_den.get<double>("fitted_values");
                      
                      ListCpp fit_num = logisregcpp(
                        data1, "cross", covariates_lgs_num, "", "", "", 
                        "uid", "logit", init, 0, firth, flic, 0, alpha);
                      DataFrameCpp sumstat_num = fit_num.get<DataFrameCpp>("sumstat");
                      if (sumstat_num.get<unsigned char>("fail")[0]) fail = true;
                      DataFrameCpp f_num = fit_num.get<DataFrameCpp>("fitted");
                      std::vector<double> h_num = f_num.get<double>("fitted_values");
                      
                      // update data_switch and fit_switch
                      if (k == -1) {
                        std::vector<int> uid = data1.get<int>("uid");
                        if (data.int_cols.count(id)) {
                          data1.push_front(subset(idwi, uid), id);
                        } else if (data.numeric_cols.count(id)) {
                          data1.push_front(subset(idwn, uid), id);
                        } else if (data.string_cols.count(id)) {
                          data1.push_front(subset(idwc, uid), id);
                        }
                        
                        if (has_stratum) {
                          std::vector<int> ustratum = data1.get<int>("ustratum");
                          for (int i = 0; i < p_stratum; ++i) {
                            const std::string& s = stratum[i];
                            if (data.bool_cols.count(s)) {
                              auto v = u_stratum.get<unsigned char>(s);
                              data1.push_back(subset(v, ustratum), s);
                            } else if (data.int_cols.count(s)) {
                              auto v = u_stratum.get<int>(s);
                              data1.push_back(subset(v, ustratum), s);
                            } else if (data.numeric_cols.count(s)) {
                              auto v = u_stratum.get<double>(s);
                              data1.push_back(subset(v, ustratum), s);
                            } else if (data.string_cols.count(s)) {
                              auto v = u_stratum.get<std::string>(s);
                              data1.push_back(subset(v, ustratum), s);
                            }
                          }
                        }
                        
                        ListPtr& data_x = data_switch[h];
                        data_x->get<DataFrameCpp>("data") = data1;
                        
                        ListPtr& fit_x = fit_switch[h];
                        fit_x->get_list("fit_den") = fit_den;
                        fit_x->get_list("fit_num") = fit_num;
                      }
                      
                      // convert to probability of observed response 
                      std::vector<double> o_den(n2), o_num(n2);
                      for (int i = 0; i < n2; ++i) {
                        o_den[i] = cross2[i] == 1 ? h_den[i] : 1 - h_den[i];
                        o_num[i] = cross2[i] == 1 ? h_num[i] : 1 - h_num[i];
                      }
                      
                      // obtain cumulative products within a subject
                      int mid = 0;
                      for (; mid < n1; ++mid) {
                        if (treat1[mid] == 1) break;
                      }
                      
                      int start, end;
                      if (h == 0) {
                        start = 0; end = mid;
                      } else {
                        start = mid; end = n1;
                      }
                      
                      // extract data for the treatment group
                      std::vector<int> id3 = subset(id1, start, end);
                      std::vector<int> swtrt3 = subset(swtrt1, start, end);
                      int n3 = end - start;
                      
                      // indices for each subject
                      std::vector<int> idx3(1, 0);
                      for (int i = 1; i < n3; ++i) {
                        if (id3[i] != id3[i-1]) {
                          idx3.push_back(i);
                        }
                      }
                      
                      // extract switch indicators for each subject
                      std::vector<int> swtrt3u = subset(swtrt3, idx3);
                      
                      int nids3 = static_cast<int>(idx3.size());
                      idx3.push_back(n3);
                      
                      // cumulative products of the observed probabilities
                      std::vector<double> p_den(n3, 1.0), p_num(n3, 1.0);
                      
                      int m = 0; // index for id2
                      for (int i = 0; i < nids3; ++i) {
                        int r = m - idx3[i] - 1;
                        int j1 = idx3[i], j2 = idx3[i+1];
                        for (int j = j1 + 1; j < j2; ++j) {
                          p_den[j] = p_den[j-1] * o_den[r + j];
                          p_num[j] = p_num[j-1] * o_num[r + j];
                        }
                        
                        int ni = j2 - j1;
                        int mi = (swtrt3u[i] == 1) ? ni : ni - 1;
                        m += mi;
                      }
                      
                      // unstabilized and stabilized weights
                      int m3 = p_den.size();
                      std::vector<double> w3(m3), sw3(m3);
                      for (int i = 0; i < m3; ++i) {
                        if (p_den[i] == 0.0) {
                          w3[i] = NaN;
                          sw3[i] = NaN;
                        } else {
                          w3[i] = 1.0 / p_den[i];
                          sw3[i] = p_num[i] / p_den[i];
                        }
                      }
                      
                      // truncate the weights if requested
                      if (trunc > 0.0) {
                        truncate_in_place(w3, trunc_upper_only, trunc);
                        truncate_in_place(sw3, trunc_upper_only, trunc);
                      }
                      
                      // summarize weights for the treatment arm
                      if (k == -1) {
                        w_treated[h] = h;
                        w_n[h] = n3;
                        if (stabilized_weights) {
                          w_min[h] = *std::min_element(sw3.begin(), sw3.end());
                          w_Q1[h] = quantilecpp(sw3, 0.25);
                          w_median[h] = quantilecpp(sw3, 0.5);
                          w_mean[h] = mean_kahan(sw3);
                          w_Q3[h] = quantilecpp(sw3, 0.75);
                          w_max[h] = *std::max_element(sw3.begin(), sw3.end());
                        } else {
                          w_min[h] = *std::min_element(w3.begin(), w3.end());
                          w_Q1[h] = quantilecpp(w3, 0.25);
                          w_median[h] = quantilecpp(w3, 0.5);
                          w_mean[h] = mean_kahan(w3);
                          w_Q3[h] = quantilecpp(w3, 0.75);
                          w_max[h] = *std::max_element(w3.begin(), w3.end());
                        }
                      }
                      
                      // fill in the weights
                      std::memcpy(w1.data() + start, w3.data(), n3 * sizeof(double));
                      std::memcpy(sw1.data() + start, sw3.data(), n3 * sizeof(double));
                    }
                    
                    // prepare data for the outcome model
                    data_outcome.push_back(std::move(id1), "uid");
                    data_outcome.push_back(std::move(tstart1), "tstart");
                    data_outcome.push_back(std::move(tstop1), "tstop");
                    data_outcome.push_back(std::move(event1), "event");
                    data_outcome.push_back(std::move(treat1), "treated");
                    data_outcome.push_back(std::move(w1), "unstabilized_weight");
                    data_outcome.push_back(std::move(sw1), "stabilized_weight");
                    data_outcome.push_back(std::move(stratum1), "ustratum");
                    
                    for (int j = 0; j < p; ++j) {
                      const std::string& zj = covariates[j+1];
                      std::vector<double> u = flatmatrix_get_column(z1, j);
                      data_outcome.push_back(std::move(u), zj);
                    }
                  }
                  
                  std::string weight_variable = stabilized_weights ? 
                  "stabilized_weight" : "unstabilized_weight";
                  
                  DataFrameCpp weight_summary, km_outcome, lr_outcome;
                  if (k == -1) {
                    weight_summary.push_back(std::move(w_treated), "treated");
                    weight_summary.push_back(std::move(w_n), "N");
                    weight_summary.push_back(std::move(w_min), "Min");
                    weight_summary.push_back(std::move(w_Q1), "Q1");
                    weight_summary.push_back(std::move(w_median), "Median");
                    weight_summary.push_back(std::move(w_mean), "Mean");
                    weight_summary.push_back(std::move(w_Q3), "Q3");
                    weight_summary.push_back(std::move(w_max), "Max");
                    
                    // generate weighted KM estimate and log-rank test
                    km_outcome = kmestcpp(
                      data_outcome, {"treated"}, "tstart", "tstop", "event", 
                      weight_variable, "log-log", 1.0 - alpha, 1);
                    
                    lr_outcome = lrtestcpp(
                      data_outcome, {"ustratum"}, "treated", "tstart", "tstop", 
                      "event", weight_variable);
                  }
                  
                  // fit the outcome model with weights
                  ListCpp fit_outcome = phregcpp(
                    data_outcome, {"ustratum"}, "tstart", "tstop", "event", 
                    covariates, weight_variable, "", "uid", ties, init, 
                    1, 0, 0, 0, 0, alpha);
                  
                  DataFrameCpp sumstat = fit_outcome.get<DataFrameCpp>("sumstat");
                  if (sumstat.get<unsigned char>("fail")[0]) fail = true;
                  
                  DataFrameCpp parest = fit_outcome.get<DataFrameCpp>("parest");
                  double beta0 = parest.get<double>("beta")[0];
                  double sebeta0 = parest.get<double>("sebeta")[0];
                  double hrhat = std::exp(beta0);
                  double hrlower = NaN, hrupper = NaN, pvalue = NaN;
                  if (k == -1) {
                    hrlower = std::exp(beta0 - zcrit * sebeta0);
                    hrupper = std::exp(beta0 + zcrit * sebeta0);
                    pvalue = parest.get<double>("p")[0];
                  }
                  
                  ListCpp out;
                  if (k == -1) {
                    out.push_back(std::move(data_switch), "data_switch");
                    out.push_back(std::move(fit_switch), "fit_switch");
                    out.push_back(std::move(data_outcome), "data_outcome");
                    out.push_back(std::move(weight_summary), "weight_summary");
                    out.push_back(std::move(km_outcome), "km_outcome");
                    out.push_back(std::move(lr_outcome), "lr_outcome");
                    out.push_back(std::move(fit_outcome), "fit_outcome");
                    out.push_back(hrhat, "hrhat");
                    out.push_back(hrlower, "hrlower");
                    out.push_back(hrupper, "hrupper");
                    out.push_back(pvalue, "pvalue");
                    out.push_back(fail, "fail");
                  } else {
                    out.push_back(hrhat, "hrhat");
                    out.push_back(fail, "fail");
                  }
                  
                  return out;
                };
  
  ListCpp out = f(idn, stratumn, tstartn, tstopn, eventn, treatn, os_timen,
                  swtrtn, swtrt_timen, zn, z_cox_denn, z_lgs_denn, -1);
  
  auto data_switch = out.get<std::vector<ListPtr>>("data_switch");
  auto fit_switch = out.get<std::vector<ListPtr>>("fit_switch");
  DataFrameCpp data_outcome = out.get<DataFrameCpp>("data_outcome");
  DataFrameCpp weight_summary = out.get<DataFrameCpp>("weight_summary");
  DataFrameCpp km_outcome = out.get<DataFrameCpp>("km_outcome");
  DataFrameCpp lr_outcome = out.get<DataFrameCpp>("lr_outcome");
  ListCpp fit_outcome = out.get_list("fit_outcome");
  double hrhat = out.get<double>("hrhat");
  double hrlower = out.get<double>("hrlower");
  double hrupper = out.get<double>("hrupper");
  double pvalue = out.get<double>("pvalue");
  bool fail = out.get<bool>("fail");
  
  std::vector<double> hrhats(n_boot);
  std::vector<unsigned char> fails(n_boot);
  DataFrameCpp fail_boots_data;
  std::string hr_CI_type;
  
  // summarize number of deaths by treatment arm in the outcome data
  std::vector<int> treated = data_outcome.get<int>("treated");
  std::vector<int> event_out = data_outcome.get<int>("event");
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
  
  treated = weight_summary.get<int>("treated");
  nottreated.resize(treated.size());
  std::transform(treated.begin(), treated.end(), nottreated.begin(),
                 [](int value) { return 1 - value; });
  if (data.bool_cols.count(treat) || data.int_cols.count(treat)) {
    weight_summary.push_back(subset(treatwi, nottreated), treat);
  } else if (data.numeric_cols.count(treat)) {
    weight_summary.push_back(subset(treatwn, nottreated), treat);
  } else if (data.string_cols.count(treat)) {
    weight_summary.push_back(subset(treatwc, nottreated), treat);
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
      subset_in_place(swtrtn, order);
      subset_in_place(swtrt_timen, order);
      subset_in_place_flatmatrix(zn, order);
      subset_in_place_flatmatrix(z_cox_denn, order);
      subset_in_place_flatmatrix(z_lgs_denn, order);
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
    boost::random::mt19937_64 master_rng(static_cast<uint64_t>(seed)); // user-provided seed
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
      const std::vector<int>& swtrtn;
      const std::vector<double>& swtrt_timen;
      const FlatMatrix& zn;
      const FlatMatrix& z_cox_denn;
      const FlatMatrix& z_lgs_denn;
      const std::vector<uint64_t>& seeds;
      
      // function f and other params that f needs are captured from outer scope
      // capture them by reference here so worker can call f(...)
      std::function<ListCpp(const std::vector<int>&, 
                            const std::vector<int>&,
                            const std::vector<double>&, 
                            const std::vector<double>&,
                            const std::vector<int>&, 
                            const std::vector<int>&,
                            const std::vector<double>&, 
                            const std::vector<int>&,
                            const std::vector<double>&, 
                            const FlatMatrix&,
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
      std::vector<int> swtrtc_local;
      std::vector<double> swtrt_timec_local;
      
      // store column-wise z_lgs_denc_local: outer vector length == z_lgs_denn.ncol
      // each inner vector stores column data across failed boots
      std::vector<std::vector<double>> z_lgs_denc_local;
      
      int index1_local = 0; // number of rows stored so far for z_lgs_denc_local
      
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
                      const std::vector<int>& swtrtn_,
                      const std::vector<double>& swtrt_timen_,
                      const FlatMatrix& zn_,
                      const FlatMatrix& z_cox_denn_,
                      const FlatMatrix& z_lgs_denn_,
                      const std::vector<uint64_t>& seeds_,
                      decltype(f) f_,
                      std::vector<unsigned char>& fails_out_,
                      std::vector<double>& hrhats_out_) :
        n(n_), nids(nids_), ntss(ntss_), idx(idx_), tsx(tsx_), idn(idn_), 
        stratumn(stratumn_), tstartn(tstartn_), tstopn(tstopn_),
        eventn(eventn_), treatn(treatn_), osn(osn_), os_timen(os_timen_),
        swtrtn(swtrtn_), swtrt_timen(swtrt_timen_), zn(zn_),
        z_cox_denn(z_cox_denn_), z_lgs_denn(z_lgs_denn_),
        seeds(seeds_), f(std::move(f_)),
        fails_out(fails_out_), hrhats_out(hrhats_out_) {
        
        // heuristic reservation to reduce reallocations:
        int ncols_lgs = z_lgs_denn.ncol;
        z_lgs_denc_local.resize(ncols_lgs);
        
        // reserve some capacity per column. This is a heuristic; adjust if needed.
        std::size_t per_col_reserve = 10 * n;
        for (int col = 0; col < ncols_lgs; ++col) 
          z_lgs_denc_local[col].reserve(per_col_reserve);
        
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
        swtrtc_local.reserve(per_col_reserve);
        swtrt_timec_local.reserve(per_col_reserve);
      }
      
      // operator() processes a range of bootstrap iterations [begin, end)
      void operator()(std::size_t begin, std::size_t end) {
        // per-worker reusable buffers (avoid reallocation per iteration)
        std::vector<int> oidb, idb, stratumb, eventb, treatb, osb, swtrtb;
        std::vector<double> tstartb, tstopb, os_timeb, swtrt_timeb;
        std::vector<std::vector<double>> zb_cols, z_cox_denb_cols, z_lgs_denb_cols;
        int ncols_z = zn.ncol, ncols_cox = z_cox_denn.ncol;
        int ncols_lgs = z_lgs_denn.ncol;
        zb_cols.resize(ncols_z); z_cox_denb_cols.resize(ncols_cox);
        z_lgs_denb_cols.resize(ncols_lgs);
        
        oidb.reserve(n); idb.reserve(n); stratumb.reserve(n); 
        tstartb.reserve(n); tstopb.reserve(n); eventb.reserve(n); 
        treatb.reserve(n); osb.reserve(n); os_timeb.reserve(n); 
        swtrtb.reserve(n); swtrt_timeb.reserve(n);
        for (int col = 0; col < ncols_z; ++col) zb_cols[col].reserve(n);
        for (int col = 0; col < ncols_cox; ++col) z_cox_denb_cols[col].reserve(n);
        for (int col = 0; col < ncols_lgs; ++col) z_lgs_denb_cols[col].reserve(n);
        
        for (std::size_t k = begin; k < end; ++k) {
          // Reset per-iteration buffers so each k starts fresh.
          // (They were reserved earlier for performance but must be emptied.)
          oidb.clear(); idb.clear(); stratumb.clear(); tstartb.clear();
          tstopb.clear(); eventb.clear(); treatb.clear(); osb.clear();
          os_timeb.clear(); swtrtb.clear(); swtrt_timeb.clear();
          for (int col = 0; col < ncols_z; ++col) zb_cols[col].clear();
          for (int col = 0; col < ncols_cox; ++col) z_cox_denb_cols[col].clear();
          for (int col = 0; col < ncols_lgs; ++col) z_lgs_denb_cols[col].clear();
          
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
              std::vector<int> swtrtn1 = subset(swtrtn, start1, end1);
              std::vector<double> swtrt_timen1 = subset(swtrt_timen, start1, end1);
              FlatMatrix zn1 = subset_flatmatrix(zn, start1, end1);
              FlatMatrix z_cox_denn1 = subset_flatmatrix(z_cox_denn, start1, end1);
              FlatMatrix z_lgs_denn1 = subset_flatmatrix(z_lgs_denn, start1, end1);
              
              append(oidb, oidn1);
              append(idb, idn1);
              append(stratumb, stratumn1);
              append(tstartb, tstartn1);
              append(tstopb, tstopn1);
              append(eventb, eventn1);
              append(treatb, treatn1);
              append(osb, osn1);
              append(os_timeb, os_timen1);
              append(swtrtb, swtrtn1);
              append(swtrt_timeb, swtrt_timen1);
              append_flatmatrix(zb_cols, zn1);
              append_flatmatrix(z_cox_denb_cols, z_cox_denn1);
              append_flatmatrix(z_lgs_denb_cols, z_lgs_denn1);
            }
          } // end block sampling
          
          // call the (thread-safe) per-iteration function f
          FlatMatrix zb = cols_to_flatmatrix(zb_cols);
          FlatMatrix z_cox_denb = cols_to_flatmatrix(z_cox_denb_cols);
          FlatMatrix z_lgs_denb = cols_to_flatmatrix(z_lgs_denb_cols);
          
          ListCpp out = f(idb, stratumb, tstartb, tstopb, eventb, treatb, os_timeb, 
                          swtrtb, swtrt_timeb, zb, z_cox_denb, z_lgs_denb,
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
            append(swtrtc_local, swtrtb);
            append(swtrt_timec_local, swtrt_timeb);
            append_flatmatrix(z_lgs_denc_local, z_lgs_denb);
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
        append(swtrtc_local, other.swtrtc_local);
        append(swtrt_timec_local, other.swtrt_timec_local);
        append_flatmatrix(z_lgs_denc_local, other.z_lgs_denc_local);
        index1_local += other.index1_local;
      }
    }; // end BootstrapWorker
    
    // Instantiate the Worker with references to inputs and outputs
    BootstrapWorker worker(
        n, nids, ntss, idx, tsx, idn, stratumn, tstartn, tstopn, eventn, treatn, 
        osn, os_timen, swtrtn, swtrt_timen, zn, z_cox_denn, z_lgs_denn, seeds,
        // bind f into std::function (capture the f we already have)
        std::function<ListCpp(const std::vector<int>&, 
                              const std::vector<int>&,
                              const std::vector<double>&, 
                              const std::vector<double>&,
                              const std::vector<int>&, 
                              const std::vector<int>&,
                              const std::vector<double>&, 
                              const std::vector<int>&,
                              const std::vector<double>&, 
                              const FlatMatrix&,
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
      fail_boots_data.push_back(std::move(worker.swtrtc_local), "swtrt");
      fail_boots_data.push_back(std::move(worker.swtrt_timec_local), "swtrt_time");
      
      int ncols_lgs = worker.z_lgs_denc_local.size();
      for (int j = 0; j < ncols_lgs; ++j) {
        const std::string& zj = covariates_lgs_den[j];
        fail_boots_data.push_back(std::move(worker.z_lgs_denc_local[j]), zj);
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
  }
  
  ListCpp result;
  std::string pvalue_type = boot ? "bootstrap" : "Cox model";
  std::vector<double> hr_CI = {hrlower, hrupper};
  result.push_back(pvalue, "pvalue");
  result.push_back(pvalue_type, "pvalue_type");
  result.push_back(hrhat, "hr");
  result.push_back(std::move(hr_CI), "hr_CI");
  result.push_back(hr_CI_type, "hr_CI_type");
  result.push_back(std::move(event_summary), "event_summary");
  result.push_back(std::move(data_switch), "data_switch");
  result.push_back(std::move(fit_switch), "fit_switch");
  result.push_back(std::move(data_outcome), "data_outcome");
  result.push_back(std::move(weight_summary), "weight_summary");
  result.push_back(std::move(km_outcome), "km_outcome");
  result.push_back(std::move(lr_outcome), "lr_outcome");
  result.push_back(std::move(fit_outcome), "fit_outcome");
  result.push_back(fail, "fail");
  
  if (boot) {
    result.push_back(fails, "fail_boots");
    result.push_back(std::move(hrhats), "hr_boots"); 
    if (std::any_of(fails.begin(), fails.end(), [](bool x){ return x; })) {
      result.push_back(std::move(fail_boots_data), "fail_boots_data");
    }
  }
  
  thread_utils::drain_thread_warnings_to_R();
  return Rcpp::wrap(result);
}
