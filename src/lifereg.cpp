// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>    // RcppParallel::Worker, parallelFor
#include <RcppThread.h>      // RcppThread::Rcerr
#include <Rcpp.h>

#include "logistic_regression.h"
#include "utilities.h"      // boost_pnorm, boost_plogis, etc.
#include "dataframe_list.h"  // FlatMatrix, IntMatrix, DataFrameCpp, ListCpp
#include "thread_utils.h"    // push_thread_warning / drain_thread_warnings_to_R

#include <vector>
#include <string>
#include <numeric>   // iota, inner_product
#include <cmath>     // isnan, isinf, fabs, NaN, exp, log
#include <stdexcept> // exceptions
#include <algorithm> // sort, none_of, any_of
#include <cstring>   // memcpy

struct aftparams {
  int dist_code; // 1: exponential, 2: weibull, 3: lognormal, 4: normal, 
  // 5: loglogistic, 6: logistic
  std::vector<int> strata;
  std::vector<double> tstart;
  std::vector<double> tstop;
  std::vector<int> status;
  std::vector<double> weight;
  std::vector<double> offset;
  FlatMatrix z;
  int nstrata;
};


// all-in-one function for log-likelihood, score, and information matrix for the AFT model
ListCpp f_der_1(int p, const std::vector<double>& par, void* ex) {
  aftparams *param = (aftparams *) ex;
  const int n = param->z.nrow;
  const int nvar = param->z.ncol;
  const int dist_code = param->dist_code;
  
  const std::vector<int>& strata = param->strata; 
  const std::vector<double>& tstart = param->tstart; 
  const std::vector<double>& tstop = param->tstop; 
  const std::vector<int>& status = param->status; 
  const std::vector<double>& weight = param->weight; 
  const std::vector<double>& offset = param->offset; 
  const double* zdata = param->z.data_ptr(); // column-major: zdata[col * n + row]
  
  // compute linear predictor eta efficiently using column-major storage:
  std::vector<double> eta = offset; // initialize with offset
  // add contributions of each coefficient times column
  for (int i = 0; i < nvar; ++i) {
    double beta = par[i];
    if (beta == 0.0) continue;
    const double* col = zdata + i * n;
    for (int r = 0; r < n; ++r) {
      eta[r] += beta * col[r];
    }
  }
  
  // Precompute sigma and logsigma per person (exponential has sigma = 1)
  std::vector<double> sigma(n, 1.0); 
  std::vector<double> logsigma(n, 0.0); 
  if (dist_code != 1) { // all except exponential 
    for (int person = 0; person < n; ++person) { 
      int k = strata[person] + nvar; // index in par for log(sigma) 
      sigma[person] = std::exp(par[k]); 
      logsigma[person] = par[k]; // par[k] == log(sigma) 
    } 
  }
  
  // Initialize accumulators
  double loglik = 0.0;
  std::vector<double> score(p);
  FlatMatrix imat(p,p);
  
  // Main loop over persons 
  for (int person = 0; person < n; ++person) {
    const double wt = weight[person]; 
    const double s = sigma[person]; 
    const double inv_s = 1.0 / s; 
    const double logsig = logsigma[person]; 
    const double eta_p = eta[person]; 
    const double tstart_p = tstart[person]; 
    const double tstop_p = tstop[person]; 
    const int st = status[person]; 
    const int k = strata[person] + nvar;
    
    // helper to get z_{j}(person) / sigma without allocating
    auto z = [&](int j)->double {
      return zdata[j * n + person] * inv_s;
    };
    
    switch (st) {
    case 1: // event
      switch (dist_code) {
      case 1: case 2: { // exponential / weibull
      double u = (std::log(tstop_p) - eta_p) * inv_s;
      double eu = std::exp(u);
      loglik += wt * (u - eu - logsig);
      
      double c1 = -wt * (1.0 - eu);
      for (int i = 0; i < nvar; ++i) score[i] += c1 * z(i);
      if (dist_code == 2)
        score[k] += wt * ((1.0 - eu) * (-u) - 1.0);
      
      c1 = wt * eu;
      for (int j = 0; j < nvar; ++j) {
        double zj = z(j);
        for (int i = j; i < nvar; ++i) {
          imat(i, j) += c1 * z(i) * zj;
        }
      }
      if (dist_code == 2) { // weibull
        double c2 = wt * (eu * u - (1.0 - eu));
        for (int j = 0; j < nvar; ++j) imat(k, j) += c2 * z(j);
        imat(k, k) += c2 * u;
      }
      break;
    }
      case 3: case 4: { // lognormal / normal
        double u;
        if (dist_code == 3) u = (std::log(tstop_p) - eta_p) * inv_s;
        else u = (tstop_p - eta_p) * inv_s;
        loglik += wt * (std::log(boost_dnorm(u)) - logsig);
        
        double c1 = wt * u;
        for (int i = 0; i < nvar; ++i) score[i] += c1 * z(i);
        score[k] += wt * (u * u - 1.0);
        
        // information: beta-beta
        for (int j = 0; j < nvar; ++j) {
          double zj = z(j);
          for (int i = j; i < nvar; ++i) {
            imat(i, j) += wt * z(i) * zj;
          }
        }
        double c2 = wt * 2.0 * u;
        for (int j = 0; j < nvar; ++j) imat(k, j) += c2 * z(j);
        imat(k, k) += c2 * u;
        break;
      }
      case 5: case 6: { // loglogistic / logistic
        double u;
        if (dist_code == 5) u = (std::log(tstop_p) - eta_p) * inv_s;
        else u = (tstop_p - eta_p) * inv_s;
        loglik += wt * (std::log(boost_dlogis(u)) - logsig);
        
        double c = 1.0 - 2.0 * boost_plogis(u, 0.0, 1.0, 0);
        double c1 = wt * c;
        for (int i = 0; i < nvar; ++i) score[i] += c1 * z(i);
        score[k] += wt * (c * u - 1.0);
        
        c1 = wt * 2.0 * boost_dlogis(u);
        double c2 = wt * (2.0 * boost_dlogis(u) * u + 
                          1.0 - 2.0 * boost_plogis(u, 0.0, 1.0, 0));
        for (int j = 0; j < nvar; ++j) {
          double zj = z(j);
          for (int i = j; i < nvar; ++i) {
            imat(i, j) += c1 * z(i) * zj;
          }
        }
        for (int j = 0; j < nvar; ++j) imat(k, j) += c2 * z(j);
        imat(k, k) += c2 * u;
        break;
      }
      }
      break;
      
    case 3: // interval censoring
      switch (dist_code) {
      case 1: case 2: {
        double u1 = (std::log(tstart_p) - eta_p) * inv_s;
        double u2 = (std::log(tstop_p)  - eta_p) * inv_s;
        double e_u1 = std::exp(u1);
        double e_u2 = std::exp(u2);
        double q1 = std::exp(-e_u1);
        double q2 = std::exp(-e_u2);
        double d1 = e_u1 * q1;
        double d2 = e_u2 * q2;
        double num = d1 - d2;
        double den = q1 - q2;
        double tmp = num / den;
        double ddu = d1 * u1 - d2 * u2;
        double term = ddu / den;
        
        loglik += wt * std::log(den);
        
        double c1 = wt * tmp;
        for (int i = 0; i < nvar; ++i) score[i] += c1 * z(i);
        if (dist_code == 2) {
          score[k] += wt * term;
        }
        
        c1 = wt * (tmp * tmp + (d1 * (1.0 - e_u1) - d2 * (1.0 - e_u2)) / den);
        for (int j = 0; j < nvar; ++j) {
          double zj = z(j);
          for (int i = j; i < nvar; ++i) {
            imat(i, j) += c1 * z(i) * zj;
          }
        }
        
        if (dist_code == 2) {
          double du1 = d1 * (1.0 + (1.0 - e_u1) * u1);
          double du2 = d2 * (1.0 + (1.0 - e_u2) * u2);
          double c2 = wt * (tmp * term +  (du1 - du2) / den);
          for (int j = 0; j < nvar; ++j) imat(k, j) += c2 * z(j);
          imat(k, k) += wt * (term * term + (du1 * u1 - du2 * u2) / den);
        }
        break;
      }
      case 3: case 4: {
        double u1, u2;
        if (dist_code == 3) {
          u1 = (std::log(tstart_p) - eta_p) * inv_s;
          u2 = (std::log(tstop_p)  - eta_p) * inv_s;
        } else {
          u1 = (tstart_p - eta_p) * inv_s;
          u2 = (tstop_p  - eta_p) * inv_s;
        }
        double q1 = boost_pnorm(u1, 0.0, 1.0, 0);
        double q2 = boost_pnorm(u2, 0.0, 1.0, 0);
        double d1 = boost_dnorm(u1);
        double d2 = boost_dnorm(u2);
        double num = d1 - d2;
        double den = q1 - q2;
        double tmp = num / den;
        double ddu = d1 * u1 - d2 * u2;
        double term = ddu / den;
        
        loglik += wt * std::log(den);
        
        double c1 = wt * tmp;
        for (int i = 0; i < nvar; ++i) score[i] += c1 * z(i);
        score[k] += wt * term;
        
        c1 = wt * ( tmp * tmp - term );
        for (int j = 0; j < nvar; ++j) {
          double zj = z(j);
          for (int i = j; i < nvar; ++i) {
            imat(i, j) += c1 * z(i) * zj;
          }
        }
        
        double du1 = d1 * (1.0 - u1 * u1);
        double du2 = d2 * (1.0 - u2 * u2);
        double c2 = wt * ( tmp * term + (du1 - du2) / den );
        for (int j = 0; j < nvar; ++j) imat(k, j) += c2 * z(j);
        imat(k, k) += wt * ( term * term + (du1 * u1 - du2 * u2) / den );
        break;
      }
      case 5: case 6: {
        double u1, u2;
        if (dist_code == 5) {
          u1 = (std::log(tstart_p) - eta_p) * inv_s;
          u2 = (std::log(tstop_p)  - eta_p) * inv_s;
        } else {
          u1 = (tstart_p - eta_p) * inv_s;
          u2 = (tstop_p  - eta_p) * inv_s;
        }
        double q1 = boost_plogis(u1, 0.0, 1.0, 0);
        double q2 = boost_plogis(u2, 0.0, 1.0, 0);
        double d1 = boost_dlogis(u1);
        double d2 = boost_dlogis(u2);
        double num = d1 - d2;
        double den = q1 - q2;
        double tmp = num / den;
        double ddu = d1 * u1 - d2 * u2;
        double term = ddu / den;
        
        loglik += wt * std::log(den);
        
        double c1 = wt * tmp;
        for (int i = 0; i < nvar; ++i) score[i] += c1 * z(i);
        score[k] += wt * term;
        
        c1 = wt * ( tmp * tmp + (d1 * (2.0 * q1 - 1.0) - d2 * (2.0 * q2 - 1.0)) / den );
        for (int j = 0; j < nvar; ++j) {
          double zj = z(j);
          for (int i = j; i < nvar; ++i) {
            imat(i, j) += c1 * z(i) * zj;
          }
        }
        
        double du1 = d1 * (1.0 + (2.0 * q1 - 1.0) * u1);
        double du2 = d2 * (1.0 + (2.0 * q2 - 1.0) * u2);
        double c2 = wt * ( tmp * term + (du1 - du2) / den );
        for (int j = 0; j < nvar; ++j) imat(k, j) += c2 * z(j);
        imat(k, k) += wt * ( term * term + (du1 * u1 - du2 * u2) / den );
        break;
      }
      } // dist_code
      break;
      
    case 2: // left censoring
      switch (dist_code) {
      case 1: case 2: {
        double u2 = (std::log(tstop_p) - eta_p) * inv_s;
        double e_u2 = std::exp(u2);
        double q2 = std::exp(-e_u2);
        double d2 = e_u2 * q2;
        double num = -d2;
        double den = 1.0 - q2;
        double tmp = num / den;
        double term = tmp * u2;
        
        loglik += wt * std::log(den);
        
        double c1 = wt * tmp;
        for (int i = 0; i < nvar; ++i) score[i] += c1 * z(i);
        if (dist_code == 2) score[k] += c1 * u2;
        
        c1 = wt * ( tmp * tmp - d2 * (1.0 - e_u2) / den );
        for (int j = 0; j < nvar; ++j) {
          double zj = z(j);
          for (int i = j; i < nvar; ++i) {
            imat(i, j) += c1 * z(i) * zj;
          }
        }
        
        if (dist_code == 2) {
          double du2 = d2 * (1.0 + (1.0 - e_u2) * u2);
          double c2 = wt * ( tmp * term - du2 / den );
          for (int j = 0; j < nvar; ++j) imat(k, j) += c2 * z(j);
          imat(k, k) += c2 * u2;
        }
        break;
      }
      case 3: case 4: {
        double u2;
        if (dist_code == 3) u2 = (std::log(tstop_p) - eta_p) * inv_s;
        else u2 = (tstop_p - eta_p) * inv_s;
        double q2 = boost_pnorm(u2, 0.0, 1.0, 0); 
        double d2 = boost_dnorm(u2);
        double num = -d2;
        double den = 1.0 - q2;
        double tmp = num / den;
        double term = tmp * u2;
        
        loglik += wt * std::log(den);
        
        double c1 = wt * tmp;
        for (int i = 0; i < nvar; ++i) score[i] += c1 * z(i);
        score[k] += c1 * u2;
        
        c1 = wt * ( tmp * tmp - term );
        for (int j = 0; j < nvar; ++j) {
          double zj = z(j);
          for (int i = j; i < nvar; ++i) {
            imat(i, j) += c1 * z(i) * zj;
          }
        }
        
        double du2 = d2 * (1.0 - u2 * u2);
        double c2 = wt * ( tmp * term - du2 / den );
        for (int j = 0; j < nvar; ++j) imat(k, j) += c2 * z(j);
        imat(k, k) += c2 * u2;
        break;
      }
      case 5: case 6: {
        double u2;
        if (dist_code == 5) u2 = (std::log(tstop_p) - eta_p) * inv_s;
        else u2 = (tstop_p - eta_p) * inv_s;
        double q2 = boost_plogis(u2, 0.0, 1.0, 0); 
        double d2 = boost_dlogis(u2);
        double num = -d2;
        double den = 1.0 - q2;
        double tmp = num / den;
        
        loglik += wt * std::log(den);
        
        double c1 = wt * tmp;
        for (int i = 0; i < nvar; ++i) score[i] += c1 * z(i);
        score[k] += c1 * u2;
        
        c1 = wt * d2;
        for (int j = 0; j < nvar; ++j) {
          double zj = z(j);
          for (int i = j; i < nvar; ++i) {
            imat(i, j) += c1 * z(i) * zj;
          }
        }
        
        double c2 = wt * (d2 * u2 - q2);
        for (int j = 0; j < nvar; ++j) imat(k, j) += c2 * z(j);
        imat(k, k) += c2 * u2;
        break;
      }
      }
      break;
      
    case 0: // right censoring
      switch (dist_code) {
      case 1: case 2: {
        double u1 = (std::log(tstart_p) - eta_p) * inv_s;
        double e_u1 = std::exp(u1);
        loglik += wt * (-e_u1);
        
        double c1 = wt * e_u1;
        for (int i = 0; i < nvar; ++i) score[i] += c1 * z(i);
        if (dist_code == 2) score[k] += c1 * u1;
        
        for (int j = 0; j < nvar; ++j) {
          double zj = z(j);
          for (int i = j; i < nvar; ++i) {
            imat(i, j) += c1 * z(i) * zj;
          }
        }
        
        if (dist_code == 2) {
          double c2 = wt * e_u1 * (1.0 + u1);
          for (int j = 0; j < nvar; ++j) imat(k, j) += c2 * z(j);
          imat(k, k) += c2 * u1;
        }
        break;
      }
      case 3: case 4: {
        double u1;
        if (dist_code == 3) u1 = (std::log(tstart_p) - eta_p) * inv_s;
        else u1 = (tstart_p - eta_p) * inv_s;
        double d1 = boost_dnorm(u1);
        double q1 = boost_pnorm(u1, 0.0, 1.0, 0);
        double tmp = d1 / q1;
        double term = tmp * u1;
        
        loglik += wt * std::log(q1);
        
        double c1 = wt * tmp;
        for (int i = 0; i < nvar; ++i) score[i] += c1 * z(i);
        score[k] += c1 * u1;
        
        c1 = wt * ( tmp * tmp - term );
        for (int j = 0; j < nvar; ++j) {
          double zj = z(j);
          for (int i = j; i < nvar; ++i) {
            imat(i, j) += c1 * z(i) * zj;
          }
        }
        
        double c2 = wt * ( tmp * term + tmp * (1.0 - u1 * u1) );
        for (int j = 0; j < nvar; ++j) imat(k, j) += c2 * z(j);
        imat(k, k) += c2 * u1;
        break;
      }
      case 5: case 6: {
        double u1;
        if (dist_code == 5) u1 = (std::log(tstart_p) - eta_p) * inv_s;
        else u1 = (tstart_p - eta_p) * inv_s;
        double q1 = boost_plogis(u1, 0.0, 1.0, 0);
        double d1 = boost_dlogis(u1); 
        double tmp  = d1 / q1;
        
        loglik += wt * std::log(q1);
        
        double c1 = wt * tmp;
        for (int i = 0; i < nvar; ++i) score[i] += c1 * z(i);
        score[k] += c1 * u1;
        
        c1 = wt * d1;
        for (int j = 0; j < nvar; ++j) {
          double zj = z(j);
          for (int i = j; i < nvar; ++i) {
            imat(i, j) += c1 * z(i) * zj;
          }
        }
        
        double c2 = wt * (1.0 - q1 + d1 * u1);
        for (int j = 0; j < nvar; ++j) imat(k, j) += c2 * z(j);
        imat(k, k) += c2 * u1;
        break;
      }
      }
      break;
      
    default:
      throw std::runtime_error("Unknown status: " + std::to_string(st));
    } // switch(status)
  } // person loop
  
  // mirror lower triangle to upper triangle (imat is symmetric) 
  for (int i = 0; i < p - 1; ++i) 
    for (int j = i + 1; j < p; ++j) 
      imat(i, j) = imat(j, i);
  
  // Build result 
  ListCpp result; 
  result.push_back(loglik, "loglik"); 
  result.push_back(std::move(score), "score"); 
  result.push_back(std::move(imat), "imat"); 
  return result; 
}


// score residual matrix
FlatMatrix f_ressco_1(int p, const std::vector<double>& par, void *ex) {
  aftparams *param = (aftparams *) ex;
  const int n = param->z.nrow;
  const int nvar = param->z.ncol;
  const int dist_code = param->dist_code;
  
  const std::vector<int>& strata = param->strata; 
  const std::vector<double>& tstart = param->tstart; 
  const std::vector<double>& tstop = param->tstop; 
  const std::vector<int>& status = param->status; 
  const std::vector<double>& offset = param->offset; 
  double* zdata = param->z.data_ptr(); // column-major: zdata[col * n + row]
  
  // compute linear predictor eta efficiently using column-major storage:
  std::vector<double> eta = offset; // initialize with offset
  // add contributions of each coefficient times column
  for (int i = 0; i < nvar; ++i) {
    double beta = par[i];
    if (beta == 0.0) continue;
    const double* col = zdata + i * n;
    for (int r = 0; r < n; ++r) {
      eta[r] += beta * col[r];
    }
  }
  
  // Precompute sigma per person (exponential has sigma = 1)
  std::vector<double> sigma(n, 1.0); 
  if (dist_code != 1) { // all except exponential 
    for (int person = 0; person < n; ++person) { 
      int k = strata[person] + nvar; // index in par for log(sigma) 
      sigma[person] = std::exp(par[k]); 
    } 
  }
  
  // Main loop to compute residuals
  FlatMatrix resid(n, p);
  double* rdata = resid.data_ptr(); // column-major: rdata[col * n + row]
  for (int person = 0; person < n; ++person) {
    const double s = sigma[person];
    const double inv_s = 1.0 / s;
    const double eta_p = eta[person];
    const double tstart_p = tstart[person];
    const double tstop_p = tstop[person];
    const int st = status[person];
    const int k = strata[person] + nvar;

    // helper to get z_{j}(person) / sigma without allocating
    auto z = [&](int j)->double {
      return zdata[j * n + person] * inv_s;
    };
    
    // helper to get resid_{j}(person) without allocating
    auto set_resid = [&](int person, int j, double val) {
      rdata[j * n + person] = val;
    };
    
    switch (st) {

    case 1: // event
      switch (dist_code) {
      case 1: case 2: { // exponential / weibull
        double u = (std::log(tstop_p) - eta_p) * inv_s;
        double c1 = -(1.0 - std::exp(u)); // -f'/f
        for (int j = 0; j < nvar; ++j) {
          set_resid(person, j, c1 * z(j));
        }
        if (dist_code == 2) {
          set_resid(person, k, c1 * u - 1.0);
        }
        break;
    }
      case 3: case 4: { // lognormal / normal
        double u;
        if (dist_code == 3)
          u = (std::log(tstop_p) - eta_p) * inv_s;
        else
          u = (tstop_p - eta_p) * inv_s;
        double c1 = u;
        for (int j = 0; j < nvar; ++j) {
          set_resid(person, j, c1 * z(j));
        }
        set_resid(person, k, c1 * u - 1.0);
        break;
      }
      case 5: case 6: { // loglogistic / logistic
        double u;
        if (dist_code == 5)
          u = (std::log(tstop_p) - eta_p) * inv_s;
        else
          u = (tstop_p - eta_p) * inv_s;
        double c1 = 1.0 - 2.0 * boost_plogis(u, 0.0, 1.0, 0);
        for (int j = 0; j < nvar; ++j) {
          set_resid(person, j, c1 * z(j));
        }
        set_resid(person, k, c1 * u - 1.0);
        break;
      }
      }
      break;
      
    case 3: // interval censoring
      switch (dist_code) {
      case 1: case 2: { // exponential / weibull
        double u1 = (std::log(tstart_p) - eta_p) * inv_s;
        double u2 = (std::log(tstop_p) - eta_p) * inv_s;
        double e_u1 = std::exp(u1);
        double e_u2 = std::exp(u2);
        double q1 = std::exp(-e_u1);
        double q2 = std::exp(-e_u2);
        double d1 = e_u1 * q1;
        double d2 = e_u2 * q2;
        double c1 = (d1 - d2) / (q1 - q2);
        for (int j = 0; j < nvar; ++j) {
          set_resid(person, j, c1 * z(j));
        }
        if (dist_code == 2) {
          set_resid(person, k, (d1 * u1 - d2 * u2) / (q1 - q2));
        }
        break;
      }
      case 3: case 4: { // lognormal / normal
        double u1, u2;
        if (dist_code == 3) {
          u1 = (std::log(tstart_p) - eta_p) * inv_s;
          u2 = (std::log(tstop_p) - eta_p) * inv_s;
        } else {
          u1 = (tstart_p - eta_p) * inv_s;
          u2 = (tstop_p - eta_p) * inv_s;
        }
        double q1 = boost_pnorm(u1, 0.0, 1.0, 0); 
        double q2 = boost_pnorm(u2, 0.0, 1.0, 0);
        double d1 = boost_dnorm(u1); 
        double d2 = boost_dnorm(u2);
        double c1 = (d1 - d2) / (q1 - q2);
        for (int j = 0; j < nvar; ++j) {
          set_resid(person, j, c1 * z(j));
        }
        set_resid(person, k, (d1 * u1 - d2 * u2) / (q1 - q2));
        break;
      }
      case 5: case 6: { // loglogistic / logistic
        double u1, u2;
        if (dist_code == 5) {
          u1 = (std::log(tstart_p) - eta_p) * inv_s;
          u2 = (std::log(tstop_p) - eta_p) * inv_s;
        } else {
          u1 = (tstart_p - eta_p) * inv_s;
          u2 = (tstop_p - eta_p) * inv_s;
        }
        double q1 = boost_plogis(u1, 0.0, 1.0, 0); 
        double q2 = boost_plogis(u2, 0.0, 1.0, 0);
        double d1 = boost_dlogis(u1); 
        double d2 = boost_dlogis(u2);
        double c1 = (d1 - d2) / (q1 - q2);
        for (int j = 0; j < nvar; ++j) {
          set_resid(person, j, c1 * z(j));
        }
        set_resid(person, k, (d1 * u1 - d2 * u2) / (q1 - q2));
        break;
      }
      }
      break;
      
    case 2: // left censoring
      switch (dist_code) {
      case 1: case 2: { // exponential / weibull
        double u2 = (std::log(tstop_p) - eta_p) * inv_s;
        double e_u2 = std::exp(u2);
        double q2 = std::exp(-e_u2);
        double d2 = e_u2 * q2;
        double c1 = -d2 / (1.0 - q2);
        for (int j = 0; j < nvar; ++j) {
          set_resid(person, j, c1 * z(j));
        }
        if (dist_code == 2) set_resid(person, k, c1 * u2);
        break;
      }
      case 3: case 4: { // lognormal / normal
        double u2;
        if (dist_code == 3)
          u2 = (std::log(tstop_p) - eta_p) * inv_s;
        else
          u2 = (tstop_p - eta_p) * inv_s;
        double q2 = boost_pnorm(u2, 0.0, 1.0, 0);
        double d2 = boost_dnorm(u2);
        double c1 = -d2 / (1.0 - q2);
        for (int j = 0; j < nvar; ++j) {
          set_resid(person, j, c1 * z(j));
        }
        set_resid(person, k, c1 * u2);
        break;
      }
      case 5: case 6: { // loglogistic / logistic
        double u2;
        if (dist_code == 5)
          u2 = (std::log(tstop_p) - eta_p) * inv_s;
        else
          u2 = (tstop_p - eta_p) * inv_s;
        double q2 = boost_plogis(u2, 0.0, 1.0, 0);
        double d2 = boost_dlogis(u2);
        double c1 = -d2 / (1.0 - q2);
        for (int j = 0; j < nvar; ++j) {
          set_resid(person, j, c1 * z(j));
        }
        set_resid(person, k, c1 * u2);
        break;
      }
      }
      break;
      
    case 0: // right censoring
      switch (dist_code) {
      case 1: case 2: { // exponential / weibull
        double u1 = (std::log(tstart_p) - eta_p) * inv_s;
        double c1 = std::exp(u1);
        for (int j = 0; j < nvar; ++j) {
          set_resid(person, j, c1 * z(j));
        }
        if (dist_code == 2) set_resid(person, k, c1 * u1);
        break;
      }
      case 3: case 4: { // lognormal / normal
        double u1;
        if (dist_code == 3)
          u1 = (std::log(tstart_p) - eta_p) * inv_s;
        else
          u1 = (tstart_p - eta_p) * inv_s;
        double q1 = boost_pnorm(u1, 0.0, 1.0, 0);
        double d1 = boost_dnorm(u1);
        double c1 = d1 / q1;
        for (int j = 0; j < nvar; ++j) {
          set_resid(person, j, c1 * z(j));
        }
        set_resid(person, k, c1 * u1);
        break;
      }
      case 5: case 6: { // loglogistic / logistic
        double u1;
        if (dist_code == 5)
          u1 = (std::log(tstart_p) - eta_p) * inv_s;
        else
          u1 = (tstart_p - eta_p) * inv_s;
        double c1 = boost_plogis(u1);
        for (int j = 0; j < nvar; ++j) {
          set_resid(person, j, c1 * z(j));
        }
        set_resid(person, k, c1 * u1);
        break;
      }
      }
      break;
      
    default:
      throw std::runtime_error("Unknown status: " + std::to_string(st));
    }
  }

  return resid;
}


// substitute information matrix guaranteed to be positive definite
FlatMatrix f_jj_1(int p, const std::vector<double>& par, void *ex) {
  aftparams *param = (aftparams *) ex;
  const int n = param->z.nrow;

  FlatMatrix resid = f_ressco_1(p, par, param);
  FlatMatrix jj(p,p);
  
  // Fast access pointers
  const double* rdata = resid.data_ptr();     // length n * p, column-major
  const double* wptr  = param->weight.data(); // length n
  double* jjdata      = jj.data_ptr();        // length p * p, column-major
  
  // Compute jj(a,b) = sum_{person=0..n-1} w[person] * resid(person,a) * resid(person,b)
  // We compute only lower triangle (b <= a) and mirror to upper triangle for efficiency.
  for (int a = 0; a < p; ++a) {
    const double* colA = rdata + a * n; // resid[:, a]
    for (int b = 0; b <= a; ++b) {
      const double* colB = rdata + b * n; // resid[:, b]
      double sum = 0.0;
      // accumulate dot product of colA and colB, weighted by wptr
      for (int person = 0; person < n; ++person) {
        sum += wptr[person] * colA[person] * colB[person];
      }
      // store at (row=a, col=b) and (row=b, col=a)
      // column-major index: data[col * nrow + row], here nrow == p for jj
      jjdata[ b * p + a ] = sum; // (a, b)
      if (a != b) jjdata[ a * p + b ] = sum; // (b, a) mirror
    }
  }
  
  return jj;
}


// underlying optimization algorithm for lifereg
ListCpp liferegloop(int p, const std::vector<double>& par, void *ex,
                    int maxiter, double eps,
                    const std::vector<int>& colfit, int ncolfit) {
  aftparams *param = (aftparams *) ex;

  int iter = 0, halving = 0;
  bool fail = false;

  int nstrata = param->nstrata;
  int nvar = param->z.ncol;
  int nsub = param->z.nrow;

  FlatMatrix z1 = param->z;
  std::vector<double> mu(nvar), sigma(nvar);
  FlatMatrix z2(nsub, nvar);

  // helpers: mean and sample sd
  auto mean_vec = [&](const double* colptr, int n) -> double {
    double s = 0.0;
    for (int i = 0; i < n; ++i) s += colptr[i];
    return s / static_cast<double>(n);
  };
  auto sd_vec = [&](const double* colptr, int n, double m) -> double {
    double s = 0.0;
    for (int i = 0; i < n; ++i) {
      double d = colptr[i] - m;
      s += d * d;
    }
    return std::sqrt(s / static_cast<double>(n - 1));
  };
  
  // --- standardize z once --- (work column-by-column using column-major layout)
  for (int c = 0; c < nvar; ++c) {
    const double* colptr = z1.data_ptr() + c * nsub;
    // compute mean and sd
    double m = mean_vec(colptr, nsub);
    double s = sd_vec(colptr, nsub, m);
    
    // check if column is indicator (only 0 or 1) -> keep m=0, s=1
    bool all_zero_or_one = true;
    for (int r = 0; r < nsub; ++r) {
      double v = colptr[r];
      if (!(v == 0.0 || v == 1.0)) { all_zero_or_one = false; break; }
    }
    if (all_zero_or_one) { m = 0.0; s = 1.0; }
    
    mu[c] = m;
    sigma[c] = s;
    
    // fill standardized column into z2
    for (int r = 0; r < nsub; ++r) {
      z2(r, c) = (colptr[r] - m) / s;
    }
  }

  // --- initial beta ---
  std::vector<double> beta(p), newbeta(p);
  beta[0] = par[0];
  for (int i = 1; i < nvar; ++i) {
    beta[i] = par[i] * sigma[i];
    beta[0] += par[i] * mu[i];
  }
  if (param->dist_code != 1)
    std::copy(par.begin() + nvar, par.end(), beta.begin() + nvar);

  // local aftparams using standardized covariates z2  
  aftparams para = {param->dist_code, param->strata, param->tstart,
                    param->tstop, param->status, param->weight,
                    param->offset, z2, nstrata};

  ListCpp der = f_der_1(p, beta, &para);
  double loglik = der.get<double>("loglik");
  double newlk = 0;
  std::vector<double> u = der.get<std::vector<double>>("score");
  FlatMatrix imat = der.get<FlatMatrix>("imat");
  FlatMatrix jj; // will be used if needed
  std::vector<double> u1(ncolfit);
  FlatMatrix imat1(ncolfit, ncolfit);
  FlatMatrix jj1(ncolfit, ncolfit);
  
  // fill u1 with selected components
  for (int i = 0; i < ncolfit; ++i) u1[i] = u[colfit[i]];
  
  // fill imat1 from imat using colfit indices
  for (int i = 0; i < ncolfit; ++i) {
    for (int j = 0; j < ncolfit; ++j) {
      imat1(i, j) = imat(colfit[i], colfit[j]);
    }
  }
  
  // --- first step: solve system using imat1 (cholesky) or fallback to jj1 ---
  if (cholesky2(imat1, ncolfit) < 0) {
    jj = f_jj_1(p, beta, &para); // substitute information matrix
    for (int i = 0; i < ncolfit; ++i)
      for (int j = 0; j < ncolfit; ++j)
        jj1(i, j) = jj(colfit[i], colfit[j]);
    cholesky2(jj1, ncolfit);
    chsolve2(jj1, ncolfit, u1);
  } else {
    chsolve2(imat1, ncolfit, u1);
  }

  // construct update vector u (length p) from solved u1
  std::fill(u.begin(), u.end(), 0.0);
  for (int i = 0; i < ncolfit; ++i) u[colfit[i]] = u1[i];
  
  // newbeta = beta + u
  for (int i = 0; i < p; ++i) newbeta[i] = beta[i] + u[i];
  
  // --- main iteration ---
  for (iter = 0; iter < maxiter; ++iter) {
    der = f_der_1(p, newbeta, &para);
    newlk = der.get<double>("loglik");
    
    fail = std::isnan(newlk) || std::isinf(newlk);
    if (!fail && halving == 0 && std::fabs(1.0 - loglik / newlk) < eps) break;
    
    if (fail || newlk < loglik) {
      ++halving;
      for (int i = 0; i < p; ++i) newbeta[i] = 0.5 * (beta[i] + newbeta[i]);
      
      // special handling of sigmas
      if (halving == 1 && param->dist_code != 1) {
        for (int i = 0; i < nstrata; ++i) {
          int idx = nvar + i;
          if (beta[idx] - newbeta[idx] > 1.1) newbeta[idx] = beta[idx] - 1.1;
        }
      }
      continue;
    }
    
    // --- update step: accept newbeta and compute next increment ---
    halving = 0;
    beta = newbeta;         // copy accepted parameters
    loglik = newlk;
    u = der.get<std::vector<double>>("score");
    imat = der.get<FlatMatrix>("imat");
    
    // extract relevant components for solving
    for (int i = 0; i < ncolfit; ++i) u1[i] = u[colfit[i]];
    
    for (int i = 0; i < ncolfit; ++i)
      for (int j = 0; j < ncolfit; ++j)
        imat1(i, j) = imat(colfit[i], colfit[j]);
    
    if (cholesky2(imat1, ncolfit) < 0) {
      jj = f_jj_1(p, beta, &para);
      for (int i = 0; i < ncolfit; ++i)
        for (int j = 0; j < ncolfit; ++j)
          jj1(i, j) = jj(colfit[i], colfit[j]);
      cholesky2(jj1, ncolfit);
      chsolve2(jj1, ncolfit, u1);
    } else {
      chsolve2(imat1, ncolfit, u1);
    }
    
    std::fill(u.begin(), u.end(), 0.0);
    for (int i = 0; i < ncolfit; ++i) u[colfit[i]] = u1[i];
    for (int i = 0; i < p; ++i) newbeta[i] = beta[i] + u[i];
  }
  
  if (iter == maxiter) fail = true;
  
  // --- rescale back (undo standardization) ---
  for (int i = 1; i < nvar; ++i) {
    newbeta[i] /= sigma[i];
    newbeta[0] -= newbeta[i] * mu[i];
  }
  
  // --- rescale the information matrix accordingly ---
  imat = der.get<FlatMatrix>("imat");
  FlatMatrix jmat = imat; // copy
  
  // adjust the top-left nvar x nvar block
  for (int i = 0; i < nvar; ++i) {
    for (int j = 0; j <= i; ++j) {
      imat(i, j) = jmat(0,0) * mu[i] * mu[j]
      + jmat(0,j) * mu[i] * sigma[j]
      + jmat(i,0) * mu[j] * sigma[i]
      + jmat(i,j) * sigma[i] * sigma[j];
      if (i != j) imat(j, i) = imat(i, j); // symmetric
    }
  }
  
  // adjust remaining rows/cols that involve shape/scale parameters
  for (int i = nvar; i < p; ++i) {
    for (int j = 0; j < nvar; ++j) {
      imat(i, j) = jmat(i,0) * mu[j] + jmat(i,j) * sigma[j];
      imat(j, i) = imat(i, j); // symmetric
    }
  }
  
  // compute variance matrix for the fitted parameters (subset colfit)
  for (int i = 0; i < ncolfit; ++i)
    for (int j = 0; j < ncolfit; ++j)
      imat1(i, j) = imat(colfit[i], colfit[j]);
  
  FlatMatrix var1 = invsympd(imat1, ncolfit); // inverse of submatrix
  FlatMatrix var(p, p); // zero-initialized
  // place var1 into the appropriate locations in the full variance matrix
  for (int i = 0; i < ncolfit; ++i)
    for (int j = 0; j < ncolfit; ++j)
      var(colfit[i], colfit[j]) = var1(i, j);
  
  // Build and return result as ListCpp
  ListCpp result;
  result.push_back(std::move(newbeta), "coef");
  result.push_back(iter, "iter");
  result.push_back(std::move(var), "var");
  result.push_back(newlk, "loglik");
  result.push_back(fail, "fail");
  return result;
}


// confidence limit of profile likelihood method
double liferegplloop(int p, const std::vector<double>& par, void *ex,
                     int maxiter, double eps,
                     int k, int direction, double l0) {
  aftparams *param = (aftparams *) ex;
  int iter = 0;
  bool fail = false;
  
  // use std::vector for parameter vectors
  std::vector<double> beta = par;
  std::vector<double> newbeta(p, 0.0);
  double loglik = 0.0;
  double newlk = 0.0;
  
  // containers for score, delta and matrices (FlatMatrix)
  std::vector<double> u(p), delta(p);
  FlatMatrix imat(p, p), jj(p, p), v(p, p);
  
  // --- first step ---
  // Evaluate derivatives / information at initial beta
  ListCpp der = f_der_1(p, beta, param);
  loglik = der.get<double>("loglik");
  u = der.get<std::vector<double>>("score");
  imat = der.get<FlatMatrix>("imat");
  
  // compute inverse of imat (with cholesky fallback to jj)
  // test cholesky on a copy
  jj = imat; // copy
  if (cholesky2(jj, p) < 0) {
    // fallback: compute substitute information jj and invert
    jj = f_jj_1(p, beta, param);
    v = invsympd(jj, p); // inv of jj
  } else {
    v = invsympd(imat, p); // inv of imat
  }
  
  // Lagrange multiplier method used in SAS PROC LOGISTIC
  double w = -quadsym(u, v);
  double underroot = -2.0 * (l0 - loglik + 0.5 * w) / v(k, k);
  double lambda = (underroot < 0.0 ? 0.0 : direction * std::sqrt(underroot));
  u[k] += lambda;
  delta = mat_vec_mult(v, u);
  for (int i = 0; i < p; ++i) newbeta[i] = beta[i] + delta[i];
  
  // --- main iteration ---
  for (iter = 0; iter < maxiter; ++iter) {
    der = f_der_1(p, newbeta, param);
    newlk = der.get<double>("loglik");
    
    // check convergence
    fail = std::isnan(newlk) || std::isinf(newlk);
    if (!fail && std::fabs(newlk - l0) < eps && w < eps) break;
    
    // update step: accept newbeta and recompute
    beta = newbeta;
    loglik = newlk;
    u = der.get<std::vector<double>>("score");
    imat = der.get<FlatMatrix>("imat");
    jj = imat; // copy for cholesky test
    
    if (cholesky2(jj, p) < 0) {
      jj = f_jj_1(p, beta, param);
      v = invsympd(jj, p);
    } else {
      v = invsympd(imat, p);
    }
    
    // Lagrange multiplier step again
    w = -quadsym(u, v);
    underroot = -2.0 * (l0 - loglik + 0.5 * w) / v(k, k);
    lambda = underroot < 0.0 ? 0.0 : direction * std::sqrt(underroot);
    u[k] += lambda;
    delta = mat_vec_mult(v, u);
    for (int i = 0; i < p; ++i) newbeta[i] = beta[i] + delta[i];
  }
  
  if (iter == maxiter) fail = true;
  if (fail) thread_utils::push_thread_warning(
      "liferegplloop did not converge within the maximum iterations.");
  
  return newbeta[k];
}


// main liferegcpp function
ListCpp liferegcpp(const DataFrameCpp& data,
                   const std::vector<std::string>& stratum,
                   const std::string time,
                   const std::string time2,
                   const std::string event,
                   const std::vector<std::string>& covariates,
                   const std::string weight,
                   const std::string offset,
                   const std::string id,
                   const std::string dist,
                   const std::vector<double>& init,
                   const bool robust,
                   const bool plci,
                   const double alpha,
                   const int maxiter,
                   const double eps) {
  
  // --- sizes and distribution normalization ---
  int n = data.nrows();
  int nvar = static_cast<int>(covariates.size()) + 1;
  if (nvar == 2 && covariates[0] == "") nvar = 1;
  
  std::string dist1 = dist;
  std::for_each(dist1.begin(), dist1.end(), [](char &c){ 
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });
  if (dist1 == "log-logistic" || dist1 == "llogistic") dist1 = "loglogistic";
  else if (dist1 == "log-normal" || dist1 == "lnormal") dist1 = "lognormal";
  else if (dist1 == "gaussian") dist1 = "normal";
  
  int dist_code = 0;
  if (dist1 == "exponential") dist_code = 1;
  else if (dist1 == "weibull") dist_code = 2;
  else if (dist1 == "lognormal") dist_code = 3;
  else if (dist1 == "normal") dist_code = 4;
  else if (dist1 == "loglogistic") dist_code = 5;
  else if (dist1 == "logistic") dist_code = 6;
  else throw std::invalid_argument("invalid distribution: " + dist1);
  
  // --- handle strata (bygroup) ---
  std::vector<int> stratumn(n);
  int p_stratum = stratum.size();
  if (!(p_stratum == 0 || (p_stratum == 1 && stratum[0] == ""))) {
    ListCpp out = bygroup(data, stratum);
    stratumn = out.get<std::vector<int>>("index");
  }
  
  std::vector<int> stratumn1 = unique_sorted(stratumn);
  int nstrata = static_cast<int>(stratumn1.size());
  int p = (dist_code == 1) ? nvar : (nvar + nstrata);
  if (dist_code == 1 && nstrata > 1) {
    throw std::invalid_argument("Stratification is not valid with the exponential distribution");
  }
  
  // --- time / time2 existence and checks ---
  if (!data.containElementNamed(time)) 
    throw std::invalid_argument("data must contain the time variable");
  std::vector<double> timen = data.get<double>(time);
  if ((dist_code == 1 || dist_code == 2 || dist_code == 3 || dist_code == 5)) {
    for (int i = 0; i < n; ++i) if (!std::isnan(timen[i]) && timen[i] <= 0.0)
      throw std::invalid_argument("time must be positive for the " + dist1 + " distribution");
  }
  
  bool has_time2 = !time2.empty() && data.containElementNamed(time2);
  std::vector<double> time2n(n);
  if (has_time2) {
    time2n = data.get<double>(time2);
    if ((dist_code == 1 || dist_code == 2 || dist_code == 3 || dist_code == 5)) {
      for (int i = 0; i < n; ++i) if (!std::isnan(time2n[i]) && time2n[i] <= 0.0)
        throw std::invalid_argument("time2 must be positive for the " + dist1 + " distribution");
    }
  }
  
  // --- event variable ---
  bool has_event = !event.empty() && data.containElementNamed(event);
  if (!has_time2 && !has_event) {
    throw std::invalid_argument("data must contain the event variable for right censored data"); 
  }
  std::vector<double> eventn(n);
  if (has_event) {
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
  }
  
  // --- build design matrix zn (n x nvar) column-major FlatMatrix ---
  FlatMatrix zn(n, nvar);
  // intercept
  for (int i = 0; i < n; ++i) zn.data[i] = 1.0;
  // covariates
  for (int j = 0; j < nvar - 1; ++j) {
    const std::string& zj = covariates[j];
    if (!data.containElementNamed(zj)) 
      throw std::invalid_argument("data must contain the variables in covariates");
    if (data.bool_cols.count(zj)) {
      const std::vector<unsigned char>& vb = data.get<unsigned char>(zj);
      int off = (j + 1) * n;
      for (int i = 0; i < n; ++i) zn.data[off + i] = vb[i] ? 1.0 : 0.0;
    } else if (data.int_cols.count(zj)) {
      const std::vector<int>& vi = data.get<int>(zj);
      int off = (j + 1) * n;
      for (int i = 0; i < n; ++i) zn.data[off + i] = static_cast<double>(vi[i]);
    } else if (data.numeric_cols.count(zj)) {
      const std::vector<double>& vd = data.get<double>(zj);
      int off = (j + 1) * n;
      for (int i = 0; i < n; ++i) zn.data[off + i] = vd[i];
    } else {
      throw std::invalid_argument("covariates must be bool, integer or numeric");
    }
  }
  
  // --- weight and offset ---
  std::vector<double> weightn(n, 1.0);
  if (!weight.empty() && data.containElementNamed(weight)) {
    weightn = data.get<double>(weight);
    for (double w : weightn) if (std::isnan(w) || w <= 0.0) 
      throw std::invalid_argument("weight must be greater than 0");
  }
  std::vector<double> offsetn(n, 0.0);
  if (!offset.empty() && data.containElementNamed(offset)) 
    offsetn = data.get<double>(offset);
  
  // --- id mapping ---
  bool has_id = !id.empty() && data.containElementNamed(id);
  std::vector<int> idn(n);
  if (!has_id) {
    std::iota(idn.begin(), idn.end(), 0);
  } else {
    if (data.int_cols.count(id)) {
      auto v = data.get<int>(id);
      auto w = unique_sorted(v);
      idn = matchcpp(v, w);
    } else if (data.numeric_cols.count(id)) {
      auto v = data.get<double>(id);
      auto w = unique_sorted(v);
      idn = matchcpp(v, w);
    } else if (data.string_cols.count(id)) {
      auto v = data.get<std::string>(id);
      auto w = unique_sorted(v);
      idn = matchcpp(v, w);
    } else throw std::invalid_argument("incorrect type for the id variable in the input data");
  }
  
  // --- exclude observations with missing values ---
  std::vector<unsigned char> sub(n, 1);
  for (int i = 0; i < n; ++i) {
    if (stratumn[i] == INT_MIN || idn[i] == INT_MIN ||
        (std::isnan(timen[i]) && std::isnan(time2n[i])) ||
        std::isnan(weightn[i]) || std::isnan(offsetn[i])) {
      sub[i] = 0;
      continue;
    }
    // check covariates columns
    for (int j = 0; j < nvar - 1; ++j) {
      double v = zn.data[(j + 1) * n + i];
      if (std::isnan(v)) { sub[i] = 0; break; }
    }
  }
  std::vector<int> keep = which(sub);
  if (keep.empty()) throw std::invalid_argument("no observations without missing values");
  
  subset_in_place(stratumn, keep);
  subset_in_place(timen, keep);
  subset_in_place(time2n, keep);
  subset_in_place(eventn, keep);
  subset_in_place(weightn, keep);
  subset_in_place(offsetn, keep);
  subset_in_place(idn, keep);
  subset_in_place_flatmatrix(zn, keep);
  n = keep.size();
  
  // sumstat data set
  double nobs, nevents;
  double loglik0, loglik1;
  int niter;
  bool fail;
  
  // parest data set
  std::vector<std::string> par(p);
  std::vector<double> b(p), seb(p), rseb(p);
  std::vector<double> z(p), expbeta(p);
  FlatMatrix vb(p, p), rvb(p, p);
  std::vector<double> lb(p), ub(p), prob(p);
  std::vector<std::string> clparm(p);
  
  // linear predictor and fitted values for all observations
  std::vector<double> linear_predictors(n);
  
  double zcrit = boost_qnorm(1.0 - alpha / 2.0);
  double xcrit = zcrit * zcrit;
  
  // unify right censored data with interval censored data
  std::vector<double> tstart(n), tstop(n);
  if (!has_time2) {
    for (int i = 0; i < n; ++i) {
      tstart[i] = timen[i];
      tstop[i] = eventn[i] == 1 ? tstart[i] : NaN;
    }
  } else {
    tstart = timen;
    tstop = time2n;
  }
  
  std::vector<int> status(n);
  for (int i = 0; i < n; ++i) {
    if (!std::isnan(tstart[i]) && !std::isnan(tstop[i]) && tstart[i] == tstop[i]) {
      status[i] = 1; // event
    } else if (!std::isnan(tstart[i]) && !std::isnan(tstop[i]) && tstart[i] < tstop[i]) {
      status[i] = 3; // interval censoring
    } else if (std::isnan(tstart[i]) && !std::isnan(tstop[i])) {
      status[i] = 2; // left censoring
    } else if (!std::isnan(tstart[i]) && std::isnan(tstop[i])) {
      status[i] = 0; // right censoring
    } else {
      status[i] = -1; // exclude the observation
    }
  }
  
  nobs = n;
  nevents = 0;
  for (int i = 0; i < n; ++i) if (status[i] == 1) ++nevents;
  if (nevents == 0) {
    for (int i = 0; i < p; ++i) {
      if (i == 0) par[i] = "(Intercept)";
      else if (i < nvar) par[i] = covariates[i-1];
      else {
        if (nstrata == 1) par[i] = "Log(scale)";
        else par[i] = std::string("Log(scale ") + std::to_string(i - nvar + 1) + ")";
      }
      
      b[i] = NaN;
      seb[i] = NaN;
      rseb[i] = NaN;
      z[i] = NaN;
      expbeta[i] = NaN;
      lb[i] = NaN;
      ub[i] = NaN;
      prob[i] = NaN;
      clparm[i] = "Wald";
    }
    
    for (int j = 0; j < p; ++j) {
      for (int i = 0; i < p; ++i) {
        vb(i,j) = NaN;
        rvb(i,j) = NaN;
      }
    }
    
    for (int i = 0; i < n; ++i) {
      linear_predictors[i] = offsetn[i];
    }
    
    loglik0 = NaN;
    loglik1 = NaN;
    niter = 0;
    fail = true;
  } else { // nevents > 0
    // exclude invalid status
    std::vector<unsigned char> good(n, 1);
    for (int i = 0; i < n; ++i) if (status[i] == -1) good[i] = 0;
    std::vector<int> q = which(good);
    int n1 = q.size();
    
    if (n1 < n) {
      subset_in_place(stratumn, q);
      subset_in_place(tstart, q);
      subset_in_place(tstop, q);
      subset_in_place(status, q);
      subset_in_place(weightn, q);
      subset_in_place(offsetn, q);
      subset_in_place(idn, q);
      subset_in_place_flatmatrix(zn, q);
    }
    
    // intercept only model
    std::vector<double> time0(n1);
    for (int i = 0; i < n1; ++i) {
      if (status[i] == 0 || status[i] == 1) { // right censoring or event
        time0[i] = tstart[i];
      } else if (status[i] == 2) { // left censoring
        time0[i] = tstop[i];
      } else if (status[i] == 3) { // interval censoring
        time0[i] = 0.5 * (tstart[i] + tstop[i]);
      }
    }
    
    std::vector<double> y0 = time0;
    if (dist_code == 1 || dist_code == 2 || dist_code == 3 || dist_code == 5) {
      for (int i = 0; i < n1; ++i) y0[i] = std::log(y0[i]);
    }
    
    // helpers: mean and sample sd
    auto mean_vec = [&](const double* colptr, int n) -> double {
      double s = 0.0;
      for (int i = 0; i < n; ++i) s += colptr[i];
      return s / static_cast<double>(n);
    };
    auto sd_vec = [&](const double* colptr, int n, double m) -> double {
      double s = 0.0;
      for (int i = 0; i < n; ++i) {
        double d = colptr[i] - m;
        s += d * d;
      }
      return std::sqrt(s / static_cast<double>(n - 1));
    };
    
    const double* yptr = y0.data();
    double int0 = mean_vec(yptr, n1);
    double logsig0 = std::log(sd_vec(yptr, n1, int0));
    
    std::vector<double> bint0(p);
    int ncolfit0 = (dist_code == 1) ? 1 : (nstrata + 1);
    std::vector<int> colfit0(ncolfit0); // indices of parameters to be fitted
    if (dist_code == 1) {
      bint0[0] = int0;
      colfit0[0] = 0;
    } else {
      bint0[0] = int0;
      for (int i = 0; i < nstrata; ++i) bint0[nvar + i] = logsig0;
      colfit0[0] = 0;
      for (int i = 0; i < nstrata; ++i) colfit0[i + 1] = nvar + i;
    }
    
    // parameter estimates and standard errors for the null model
    aftparams param = {dist_code, stratumn, tstart, tstop, status, weightn, offsetn, zn, nstrata};
    ListCpp outint = liferegloop(p, bint0, &param, maxiter, eps, colfit0, ncolfit0);
    
    std::vector<double> bint = outint.get<std::vector<double>>("coef");
    FlatMatrix vbint = outint.get<FlatMatrix>("var");
    
    ListCpp out;
    
    if (nvar > 1) {
      std::vector<int> colfit = seqcpp(0, p-1);
      if (static_cast<int>(init.size()) == p && 
          std::none_of(init.begin(), init.end(), [](double val){ 
            return std::isnan(val); })) {
        out = liferegloop(p, init, &param, maxiter, eps, colfit, p);
      } else {
        out = liferegloop(p, bint, &param, maxiter, eps, colfit, p);
      }
      
      fail = out.get<bool>("fail");
      
      if (fail) {
        // obtain initial values for model parameters using OLS
        std::vector<double> y1(n1);
        for (int i = 0; i < n1; ++i) y1[i] = y0[i] - offsetn[i];
        
        FlatMatrix v1(nvar, nvar); // X'WX
        const double *zdata = zn.data_ptr();
        double *vptr = v1.data_ptr();
        for (int j = 0; j < nvar; ++j) {
          const double* zj = zdata + j * n1;          // pointer to Z(:,j)
          for (int k = j; k < nvar; ++k) {
            const double* zk = zdata + k * n1;        // pointer to Z(:,k)
            double sum = 0.0;
            // inner loop reads zj[i] and zk[i] contiguously
            for (int i = 0; i < n1; ++i) {
              sum += weightn[i] * zj[i] * zk[i];
            }
            // write into v1(j,k) and mirror
            vptr[k * nvar + j] = sum;                
            if (k != j) {
              vptr[j * nvar + k] = sum;            
            }
          }
        }
        
        std::vector<double> u1(nvar, 0.0); // X'Wy
        for (int j = 0; j < nvar; ++j) {
          const double* zj = zdata + j * n1;          // pointer to Z(:,j)
          double sum = 0.0;
          // inner loop reads zj[i] contiguously
          for (int i = 0; i < n1; ++i) {
            sum += weightn[i] * zj[i] * y1[i];
          }
          u1[j] = sum;
        }
        
        cholesky2(v1, nvar);
        chsolve2(v1, nvar, u1);
        
        std::vector<double> binit(p);
        for (int j = 0; j < nvar; ++j) binit[j] = u1[j];
        
        if (dist_code != 1) {
          double ssum = 0.0;
          double wsum = 0.0;
          for (int i = 0; i < n1; ++i) {
            double pred = 0.0;
            for (int j = 0; j < nvar; ++j) pred += zn(i,j) * u1[j];
            double r = y1[i] - pred;
            ssum += weightn[i] * r * r;
            wsum += weightn[i];
          }
          double s = 0.5 * std::log(ssum / wsum * n1 / std::max(1, n1 - nvar));  // log(sigma)
          for (int j = nvar; j < p; ++j) binit[j] = s;
        }
        
        // fit the model using the initial values
        out = liferegloop(p, binit, &param, maxiter, eps, colfit, p);
        fail = out.get<bool>("fail");
      }
      
      if (fail) thread_utils::push_thread_warning("The algorithm in liferegr did not converge");
      
      b = out.get<std::vector<double>>("coef");
      vb = out.get<FlatMatrix>("var");
    } else {
      // intercept-only
      b = bint;
      vb = vbint;
      out = outint;
    }
    
    // compute standard errors
    for (int j = 0; j < p; ++j) {
      seb[j] = std::sqrt(vb(j,j));
    }
    
    // fill parest outputs
    for (int i = 0; i < p; ++i) {
      if (i == 0) par[i] = "(Intercept)";
      else if (i < nvar) par[i] = covariates[i-1];
      else {
        if (nstrata == 1) par[i] = "Log(scale)";
        else par[i] = std::string("Log(scale ") + std::to_string(i - nvar + 1) + ")";
      }
    }
    
    // fill linear predictors
    for (int i = 0; i < nvar; ++i) {
      double beta = b[i];
      if (beta == 0.0) continue;
      int off = i * n1;
      for (int r = 0; r < n1; ++r) {
        linear_predictors[q[r]] += beta * zn.data[off + r];
      }
    }
    
    niter = out.get<int>("iter");
    fail = out.get<bool>("fail");
    
    // robust variance estimates
    if (robust) {
      FlatMatrix ressco = f_ressco_1(p, b, &param);
      
      int nr; // number of rows in the score residual matrix
      if (!has_id) { // no clustering, just weight the score residuals
        const double *wptr = weightn.data();
        double *resscoptr = ressco.data_ptr();
        
        for (int j = 0; j < p; ++j) {
          double* colptr = resscoptr + j * n1;
          for (int i = 0; i < n1; ++i) {
            colptr[i] *= wptr[i];
          }
        }
        nr = n1;
      } else { // sum up the score residuals by id
        std::vector<int> order = seqcpp(0, n1-1);
        std::sort(order.begin(), order.end(), [&](int i, int j) {
          return idn[i] < idn[j];
        });
        
        std::vector<int> id1 = subset(idn, order);
        std::vector<int> idx(1,0);
        for (int i = 1; i < n1; ++i) {
          if (id1[i] != id1[i-1]) {
            idx.push_back(i);
          }
        }
        int nids = idx.size();
        idx.push_back(n1);
        
        FlatMatrix ressco1(nids, p); // score residuals summed by id
        for (int j = 0; j < p; ++j) {
          const int coloff = j * nids;
          for (int i = 0; i < nids; ++i) {
            double sum = 0.0;
            for (int k = idx[i]; k < idx[i+1]; ++k) {
              int row = order[k];
              sum  += weightn[row] * ressco(row,j);
            }
            ressco1.data[coloff + i] = sum;
          }
        }
        
        ressco = std::move(ressco1);  // update the score residuals
        nr = nids;
      }
      
      FlatMatrix D = mat_mat_mult(ressco, vb); // DFBETA
      
      const double* Dptr = D.data_ptr();
      double* rvbptr = rvb.data_ptr();
      for (int j = 0; j < p; ++j) {
        const double* Dj = Dptr + j * nr; // pointer to D(:,j)
        for (int k = 0; k <= j; ++k) {
          const double* Dk = Dptr + k * nr; // pointer to D(:,k)
          double sum = 0.0;
          for (int i = 0; i < nr; ++i) {
            sum += Dj[i] * Dk[i];
          }
          rvbptr[k * p + j] = sum;
          if (j != k) rvbptr[j * p + k] = sum;
        }
      }
      
      for (int i = 0; i < p; ++i) {
        rseb[i] = std::sqrt(rvb(i,i));
      }
    }
    
    // profile likelihood confidence interval for regression coefficients
    if (plci) {
      double lmax = out.get<double>("loglik");
      double l0 = lmax - 0.5 * xcrit;
      
      for (int k = 0; k < p; ++k) {
        lb[k] = liferegplloop(p, b, &param, maxiter, eps, k, -1, l0);
        ub[k] = liferegplloop(p, b, &param, maxiter, eps, k, 1, l0);
        
        std::vector<int> colfit1(p-1);
        for (int i = 0, j = 0; i < p; ++i) {
          if (i == k) continue;
          colfit1[j++] = i;
        }
        
        std::vector<double> b0(p);
        ListCpp out0 = liferegloop(p, b0, &param, maxiter, eps, colfit1, p-1);
        double lmax0 = out0.get<double>("loglik");
        prob[k] = boost_pchisq(-2.0 * (lmax0 - lmax), 1, 0);
        clparm[k] = "PL";
      }
    } else {
      for (int k = 0; k < p; ++k) {
        if (!robust) {
          lb[k] = b[k] - zcrit * seb[k];
          ub[k] = b[k] + zcrit * seb[k];
          prob[k] = boost_pchisq(sq(b[k] / seb[k]), 1, 0);
        } else {
          lb[k] = b[k] - zcrit * rseb[k];
          ub[k] = b[k] + zcrit * rseb[k];
          prob[k] = boost_pchisq(sq(b[k] / rseb[k]), 1, 0);
        }
        clparm[k] = "Wald";
      }
    }
    
    // log-likelihoods
    loglik0 = outint.get<double>("loglik");
    loglik1 = out.get<double>("loglik");
    
    // compute exp(beta)
    for (int i = 0; i < p; ++i) {
      expbeta[i] = std::exp(b[i]);
    }
    
    // compute z statistics
    if (robust) {
      for (int i = 0; i < p; ++i) {
        if (rseb[i] == 0) {
          z[i] = NaN;
        } else {
          z[i] = b[i]/rseb[i];
        }
      }
    } else {
      for (int i = 0; i < p; ++i) {
        if (seb[i] == 0) {
          z[i] = NaN;
        } else {
          z[i] = b[i]/seb[i];
        }
      }
    }
  }
  
  DataFrameCpp sumstat;
  sumstat.push_back(nobs, "n");
  sumstat.push_back(nevents, "nevents");
  sumstat.push_back(loglik0, "loglik0");
  sumstat.push_back(loglik1, "loglik1");
  sumstat.push_back(niter, "niter");
  sumstat.push_back(dist1, "dist");
  sumstat.push_back(p, "p");
  sumstat.push_back(nvar - 1, "nvar");
  sumstat.push_back(robust, "robust");
  sumstat.push_back(fail, "fail");
  
  std::vector sebeta = robust ? rseb : seb;
  FlatMatrix vbeta = robust ? rvb : vb;
  
  DataFrameCpp parest;
  parest.push_back(std::move(par), "param");
  parest.push_back(std::move(b), "beta");
  parest.push_back(std::move(sebeta), "sebeta");
  parest.push_back(std::move(z), "z");
  parest.push_back(std::move(expbeta), "expbeta");
  parest.push_back(std::move(vbeta), "vbeta");
  parest.push_back(std::move(lb), "lower");
  parest.push_back(std::move(ub), "upper");
  parest.push_back(std::move(prob), "p");
  parest.push_back(std::move(clparm), "method");
  if (robust) {
    parest.push_back(std::move(seb), "sebeta_naive");
    parest.push_back(std::move(vb), "vbeta_naive");
  }
  
  ListCpp result;
  result.push_back(std::move(sumstat), "sumstat");
  result.push_back(std::move(parest), "parest");
  result.push_back(std::move(linear_predictors), "linear_predictors");
  
  return result;
}


// [[Rcpp::export]]
Rcpp::List liferegRcpp(const Rcpp::DataFrame& data,
                       const std::vector<std::string>& stratum,
                       const std::string time,
                       const std::string time2,
                       const std::string event,
                       const std::vector<std::string>& covariates,
                       const std::string weight,
                       const std::string offset,
                       const std::string id,
                       const std::string dist,
                       const std::vector<double>& init,
                       const bool robust,
                       const bool plci,
                       const double alpha,
                       const int maxiter,
                       const double eps) {
  
  DataFrameCpp dfcpp = convertRDataFrameToCpp(data);
  
  ListCpp cpp_result = liferegcpp(
    dfcpp, stratum, time, time2, event, covariates, weight, offset, id, 
    dist, init, robust, plci, alpha, maxiter, eps
  );
  
  thread_utils::drain_thread_warnings_to_R();
  return Rcpp::wrap(cpp_result);
}


// first and second derivatives of log likelihood with respect to eta and tau
ListCpp f_ld_1(std::vector<double>& eta, std::vector<double>& sigma, void *ex) {
  aftparams *param = (aftparams *) ex;
  const int n = param->tstart.size();
  const int dist_code = param->dist_code;
  
  const std::vector<double>& tstart = param->tstart; 
  const std::vector<double>& tstop = param->tstop; 
  const std::vector<int>& status = param->status; 
  
  std::vector<double> g(n), dg(n), ddg(n), ds(n), dds(n), dsg(n);
  
  // Main loop to compute  derivatives
  for (int person = 0; person < n; ++person) {
    const double s = sigma[person];
    const double inv_s = 1.0 / s;
    const double logsig = std::log(s);
    const double eta_p = eta[person];
    const double tstart_p = tstart[person]; 
    const double tstop_p = tstop[person];
    const int st = status[person];
    
    double vg = NaN, vdg = NaN, vddg = NaN, vds = NaN, vdds = NaN, vdsg = NaN;
    
    switch (st) {
    
    case 1: // event
      switch (dist_code) {
      case 1: case 2: { // exponential / weibull
      double u = (std::log(tstop_p) - eta_p) * inv_s;
      double eu = std::exp(u);
      double su = s * u;
      vg = u - eu - logsig;        // log f(u) - log(s)
      vdg = -(1.0 - eu) * inv_s;   // -f'(u)/f(u) * 1/s
      vddg = -eu * inv_s * inv_s;  // (f''(u)/f(u) - (f'(u)/f(u))^2) * 1/s^2
      vds = -1.0 + vdg * su;
      vdds = vddg * su * su - vdg * su;
      vdsg = vddg * su - vdg;
      break;
    }
      case 3: case 4: { // lognormal / normal
        double u;
        if (dist_code == 3) u = (std::log(tstop_p) - eta_p) * inv_s;
        else u = (tstop_p - eta_p) * inv_s;
        double su = s * u;
        double d = boost_dnorm(u);
        vg = std::log(d) - logsig;
        vdg = u * inv_s;
        vddg = -inv_s * inv_s;
        vds = -1.0 + vdg * su;
        vdds = vddg * su * su - vdg * su;
        vdsg = vddg * su - vdg;
        break;
      }
      case 5: case 6: { // loglogistic / logistic
        double u;
        if (dist_code == 5) u = (std::log(tstop_p) - eta_p) * inv_s;
        else u = (tstop_p - eta_p) * inv_s;
        double su = s * u;
        double q = boost_plogis(u, 0.0, 1.0, 0);
        double d = boost_dlogis(u);
        vg = std::log(d) - logsig;
        vdg = (1.0 - 2.0 * q) * inv_s;
        vddg = -2.0 * d * inv_s * inv_s;
        vds = -1.0 + vdg * su;
        vdds = vddg * su * su - vdg * su;
        vdsg = vddg * su - vdg;
        break;
      }
      }
      break;
      
    case 3: // interval censoring
      switch (dist_code) {
      case 1: case 2: { // exponential / weibull
        double u1 = (std::log(tstart_p) - eta_p) * inv_s;
        double u2 = (std::log(tstop_p) - eta_p) * inv_s;
        double e_u1 = std::exp(u1); 
        double e_u2 = std::exp(u2);
        double q1 = std::exp(-e_u1); 
        double q2 = std::exp(-e_u2);
        double d1 = e_u1 * q1; 
        double d2 = e_u2 * q2;
        double r1 = 1.0 - e_u1;
        double r2 = 1.0 - e_u2;
        double den = q1 - q2;
        vg = std::log(den);
        vdg = (d1 - d2) / den * inv_s;
        vddg = -(d1 * r1 - d2 * r2) / den * inv_s * inv_s - vdg * vdg;
        vds = (u1 * d1 - u2 * d2) / den;
        vdds = (u2 * u2 * d2 * r2 - u1 * u1 * d1 * r1) / den - vds * (1.0 + vds);
        vdsg = (u2 * d2 * r2 - u1 * d1 * r1) / den * inv_s - vdg * (1.0 + vds);
        break;
      }
      case 3: case 4: { // lognormal / normal
        double u1, u2;
        if (dist_code == 3) {
          u1 = (std::log(tstart_p) - eta_p) * inv_s;
          u2 = (std::log(tstop_p) - eta_p) * inv_s;
        } else {
          u1 = (tstart_p - eta_p) * inv_s;
          u2 = (tstop_p - eta_p) * inv_s;
        }
        double q1 = boost_pnorm(u1, 0.0, 1.0, 0); 
        double q2 = boost_pnorm(u2, 0.0, 1.0, 0);
        double d1 = boost_dnorm(u1); 
        double d2 = boost_dnorm(u2);
        double r1 = -u1;
        double r2 = -u2;
        double den = q1 - q2;
        vg = std::log(den);
        vdg = (d1 - d2) / den * inv_s;
        vddg = -(d1 * r1 - d2 * r2) / den * inv_s * inv_s - vdg * vdg;
        vds = (u1 * d1 - u2 * d2) / den;
        vdds = (u2 * u2 * d2 * r2 - u1 * u1 * d1 * r1) / den - vds * (1.0 + vds);
        vdsg = (u2 * d2 * r2 - u1 * d1 * r1) / den * inv_s - vdg * (1.0 + vds);
        break;
      }
      case 5: case 6: { // loglogistic / logistic
        double u1, u2;
        if (dist_code == 5) {
          u1 = (std::log(tstart_p) - eta_p) * inv_s;
          u2 = (std::log(tstop_p) - eta_p) * inv_s;
        } else {
          u1 = (tstart_p - eta_p) * inv_s;
          u2 = (tstop_p - eta_p) * inv_s;
        }
        double d1 = boost_dlogis(u1); 
        double d2 = boost_dlogis(u2);
        double q1 = boost_plogis(u1, 0.0, 1.0, 0); 
        double q2 = boost_plogis(u2, 0.0, 1.0, 0);
        double r1 = 2.0 * q1 - 1.0;
        double r2 = 2.0 * q2 - 1.0;
        double den = q1 - q2;
        vg = std::log(den);
        vdg = (d1 - d2) / den * inv_s;
        vddg = -(d1 * r1 - d2 * r2) / den * inv_s * inv_s - vdg * vdg;
        vds = (u1 * d1 - u2 * d2) / den;
        vdds = (u2 * u2 * d2 * r2 - u1 * u1 * d1 * r1) / den - vds * (1.0 + vds);
        vdsg = (u2 * d2 * r2 - u1 * d1 * r1) / den * inv_s - vdg * (1.0 + vds);
        break;
      }
      }
      break;
      
    case 2: // left censoring
      switch (dist_code) {
      case 1: case 2: { // exponential / weibull
        double u2 = (std::log(tstop_p) - eta_p) * inv_s;
        double e_u2 = std::exp(u2);
        double q2 = std::exp(-e_u2);
        double d2 = e_u2 * q2;
        double r2 = 1.0 - e_u2;
        double den = 1.0 - q2;
        double su = s * u2;
        vg = std::log(den);
        vdg = -d2 / den * inv_s;
        vddg = -vdg * r2 * inv_s - vdg * vdg;
        vds = vdg * su;
        vdds = vddg * su * su - vds;
        vdsg = vddg * su - vdg;
        break;
      }
      case 3: case 4: { // lognormal / normal
        double u2;
        if (dist_code == 3) u2 = (std::log(tstop_p) - eta_p) * inv_s;
        else u2 = (tstop_p - eta_p) * inv_s;
        double q2 = boost_pnorm(u2, 0.0, 1.0, 0);
        double d2 = boost_dnorm(u2);
        double r2 = -u2;
        double den = 1.0 - q2;
        double su = s * u2;
        vg = std::log(den);
        vdg = -d2 / den * inv_s;
        vddg = -vdg * r2 * inv_s - vdg * vdg;
        vds = vdg * su;
        vdds = vddg * su * su - vds;
        vdsg = vddg * su - vdg;
        break;
      }
      case 5: case 6: { // loglogistic / logistic
        double u2;
        if (param->dist_code==5) u2 = (std::log(tstop_p) - eta_p) * inv_s;
        else u2 = (tstop_p - eta_p) * inv_s;
        double q2 = boost_plogis(u2, 0.0, 1.0, 0);
        double d2 = boost_dlogis(u2);
        double den = 1.0 - q2;
        double su = s * u2;
        vg = std::log(den);
        vdg = -q2 * inv_s;
        vddg = -d2 * inv_s * inv_s;
        vds = vdg * su;
        vdds = vddg * su * su - vds;
        vdsg = vddg * su - vdg;
        break;
      }
      }
      break;
      
    case 0: // right censoring
      switch (dist_code) {
      case 1: case 2: { // exponential / weibull
        double u1 = (std::log(tstart_p) - eta_p) * inv_s;
        double e_u1 = std::exp(u1);
        double su = s * u1;
        vg = -e_u1;
        vdg = e_u1 * inv_s;
        vddg = -e_u1 * inv_s * inv_s;
        vds = vdg * su;
        vdds = vddg * su * su - vds;
        vdsg = vddg * su - vdg;
        break;
      }
      case 3: case 4: { // lognormal / normal
        double u1;
        if (dist_code==3) u1 = (std::log(tstart_p) - eta_p) * inv_s;
        else u1 = (tstart_p - eta_p) * inv_s;
        double q1 = boost_pnorm(u1, 0.0, 1.0, 0);
        double d1 = boost_dnorm(u1);
        double r1 = -u1;
        double su = s * u1;
        vg = std::log(q1);
        vdg = d1 / q1 * inv_s;
        vddg = -vdg * r1 * inv_s - vdg * vdg;
        vds = vdg * su;
        vdds = vddg * su * su - vds;
        vdsg = vddg * su - vdg;
        break;
      }
      case 5: case 6: { // loglogistic / logistic
        double u1;
        if (dist_code == 5) u1 = (std::log(tstart_p) - eta_p) * inv_s;
        else u1 = (tstart_p - eta_p) * inv_s;
        double q1 = boost_plogis(u1, 0.0, 1.0, 0);
        double d1 = boost_dlogis(u1);
        double su = s * u1;
        vg = std::log(q1);
        vdg = (1.0 - q1) * inv_s;
        vddg = -d1 * inv_s * inv_s;
        vds = vdg * su;
        vdds = vddg * su * su - vds;
        vdsg = vddg * su - vdg;
        break;
      }
      }
      break;
      
    default:
      throw std::runtime_error("Unknown status: " + std::to_string(st));
    }
    
    g[person] = vg;
    dg[person] = vdg;
    ddg[person] = vddg;
    ds[person] = vds;
    dds[person] = vdds;
    dsg[person] = vdsg;
  }
  
  ListCpp result;
  result.push_back(std::move(g), "g");
  result.push_back(std::move(dg), "dg");
  result.push_back(std::move(ddg), "ddg");
  result.push_back(std::move(ds), "ds");
  result.push_back(std::move(dds), "dds");
  result.push_back(std::move(dsg), "dsg");
  
  return result;
}


// residuals of the AFT model
FlatMatrix residuals_liferegcpp(const std::vector<double>& beta,
                                const FlatMatrix& vbeta,
                                const DataFrameCpp& data,
                                const std::vector<std::string>& stratum,
                                const std::string time,
                                const std::string time2,
                                const std::string event,
                                const std::vector<std::string>& covariates,
                                const std::string weight,
                                const std::string offset,
                                const std::string id,
                                const std::string dist,
                                const std::string type,
                                bool collapse,
                                bool weighted) {
  // --- sizes and distribution normalization ---
  int n = data.nrows();
  int nvar = static_cast<int>(covariates.size()) + 1;
  if (nvar == 2 && covariates[0] == "") nvar = 1;

  std::string dist1 = dist;
  std::for_each(dist1.begin(), dist1.end(), [](char &c){
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });
  if (dist1 == "log-logistic" || dist1 == "llogistic") dist1 = "loglogistic";
  else if (dist1 == "log-normal" || dist1 == "lnormal") dist1 = "lognormal";
  else if (dist1 == "gaussian") dist1 = "normal";

  int dist_code = 0;
  if (dist1 == "exponential") dist_code = 1;
  else if (dist1 == "weibull") dist_code = 2;
  else if (dist1 == "lognormal") dist_code = 3;
  else if (dist1 == "normal") dist_code = 4;
  else if (dist1 == "loglogistic") dist_code = 5;
  else if (dist1 == "logistic") dist_code = 6;
  else throw std::invalid_argument("invalid distribution: " + dist1);

  // --- handle strata (bygroup) ---
  std::vector<int> stratumn(n);
  int p_stratum = stratum.size();
  if (!(p_stratum == 0 || (p_stratum == 1 && stratum[0] == ""))) {
    ListCpp out = bygroup(data, stratum);
    stratumn = out.get<std::vector<int>>("index");
  }

  std::vector<int> stratumn1 = unique_sorted(stratumn);
  int nstrata = static_cast<int>(stratumn1.size());
  int p = (dist_code == 1) ? nvar : (nvar + nstrata);
  if (dist_code == 1 && nstrata > 1) {
    throw std::invalid_argument("Stratification is not valid with the exponential distribution");
  }

  // --- time / time2 existence and checks ---
  if (!data.containElementNamed(time))
    throw std::invalid_argument("data must contain the time variable");
  std::vector<double> timen = data.get<double>(time);
  if ((dist_code == 1 || dist_code == 2 || dist_code == 3 || dist_code == 5)) {
    for (int i = 0; i < n; ++i) if (!std::isnan(timen[i]) && timen[i] <= 0.0)
      throw std::invalid_argument("time must be positive for the " + dist1 + " distribution");
  }

  bool has_time2 = !time2.empty() && data.containElementNamed(time2);
  std::vector<double> time2n(n);
  if (has_time2) {
    time2n = data.get<double>(time2);
    if ((dist_code == 1 || dist_code == 2 || dist_code == 3 || dist_code == 5)) {
      for (int i = 0; i < n; ++i) if (!std::isnan(time2n[i]) && time2n[i] <= 0.0)
        throw std::invalid_argument("time2 must be positive for the " + dist1 + " distribution");
    }
  }

  // --- event variable ---
  bool has_event = !event.empty() && data.containElementNamed(event);
  if (!has_time2 && !has_event) {
    throw std::invalid_argument("data must contain the event variable for right censored data");
  }
  std::vector<double> eventn(n);
  if (has_event) {
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
  }

  // --- build design matrix zn (n x nvar) column-major FlatMatrix ---
  FlatMatrix zn(n, nvar);
  // intercept
  for (int i = 0; i < n; ++i) zn.data[i] = 1.0;
  // covariates
  for (int j = 0; j < nvar - 1; ++j) {
    const std::string& zj = covariates[j];
    if (!data.containElementNamed(zj))
      throw std::invalid_argument("data must contain the variables in covariates");
    if (data.bool_cols.count(zj)) {
      const std::vector<unsigned char>& vb = data.get<unsigned char>(zj);
      int off = (j + 1) * n;
      for (int i = 0; i < n; ++i) zn.data[off + i] = vb[i] ? 1.0 : 0.0;
    } else if (data.int_cols.count(zj)) {
      const std::vector<int>& vi = data.get<int>(zj);
      int off = (j + 1) * n;
      for (int i = 0; i < n; ++i) zn.data[off + i] = static_cast<double>(vi[i]);
    } else if (data.numeric_cols.count(zj)) {
      const std::vector<double>& vd = data.get<double>(zj);
      int off = (j + 1) * n;
      for (int i = 0; i < n; ++i) zn.data[off + i] = vd[i];
    } else {
      throw std::invalid_argument("covariates must be bool, integer or numeric");
    }
  }

  // --- weight and offset ---
  std::vector<double> weightn(n, 1.0);
  if (!weight.empty() && data.containElementNamed(weight)) {
    weightn = data.get<double>(weight);
    for (double w : weightn) if (std::isnan(w) || w <= 0.0)
      throw std::invalid_argument("weight must be greater than 0");
  }
  std::vector<double> offsetn(n, 0.0);
  if (!offset.empty() && data.containElementNamed(offset))
    offsetn = data.get<double>(offset);

  // --- id mapping ---
  bool has_id = !id.empty() && data.containElementNamed(id);
  std::vector<int> idn(n);
  if (!has_id) {
    std::iota(idn.begin(), idn.end(), 0);
  } else {
    if (data.int_cols.count(id)) {
      auto v = data.get<int>(id);
      auto w = unique_sorted(v);
      idn = matchcpp(v, w);
    } else if (data.numeric_cols.count(id)) {
      auto v = data.get<double>(id);
      auto w = unique_sorted(v);
      idn = matchcpp(v, w);
    } else if (data.string_cols.count(id)) {
      auto v = data.get<std::string>(id);
      auto w = unique_sorted(v);
      idn = matchcpp(v, w);
    } else throw std::invalid_argument("incorrect type for the id variable in the input data");
  }

  // --- exclude observations with missing values ---
  std::vector<unsigned char> sub(n, 1);
  for (int i = 0; i < n; ++i) {
    if (stratumn[i] == INT_MIN || idn[i] == INT_MIN ||
        (std::isnan(timen[i]) && std::isnan(time2n[i])) ||
        std::isnan(weightn[i]) || std::isnan(offsetn[i])) {
      sub[i] = 0;
      continue;
    }
    // check covariates columns
    for (int j = 0; j < nvar - 1; ++j) {
      double v = zn.data[(j + 1) * n + i];
      if (std::isnan(v)) { sub[i] = 0; break; }
    }
  }
  std::vector<int> keep = which(sub);
  if (keep.empty()) throw std::invalid_argument("no observations without missing values");

  subset_in_place(stratumn, keep);
  subset_in_place(timen, keep);
  subset_in_place(time2n, keep);
  subset_in_place(eventn, keep);
  subset_in_place(weightn, keep);
  subset_in_place(offsetn, keep);
  subset_in_place(idn, keep);
  subset_in_place_flatmatrix(zn, keep);
  n = static_cast<int>(keep.size());

  // --- tstart / tstop and status ---
  std::vector<double> tstart(n), tstop(n);
  if (!has_time2) {
    for (int i = 0; i < n; ++i) {
      tstart[i] = timen[i];
      tstop[i] = (eventn[i] == 1) ? tstart[i] : NaN;
    }
  } else {
    tstart = timen;
    tstop = time2n;
  }
  std::vector<int> status(n);
  for (int i = 0; i < n; ++i) {
    if (!std::isnan(tstart[i]) && !std::isnan(tstop[i]) && tstart[i] == tstop[i]) status[i] = 1;
    else if (!std::isnan(tstart[i]) && !std::isnan(tstop[i]) && tstart[i] < tstop[i]) status[i] = 3;
    else if (std::isnan(tstart[i]) && !std::isnan(tstop[i])) status[i] = 2;
    else if (!std::isnan(tstart[i]) && std::isnan(tstop[i])) status[i] = 0;
    else status[i] = -1;
  }

  // exclude invalid status
  std::vector<unsigned char> good(n, 1);
  for (int i = 0; i < n; ++i) if (status[i] == -1) good[i] = 0;
  std::vector<int> q = which(good);
  int n1 = static_cast<int>(q.size());
  if (n1 == 0) throw std::invalid_argument("no valid records after status filtering");

  if (n1 < n) {
    subset_in_place(stratumn, q);
    subset_in_place(tstart, q);
    subset_in_place(tstop, q);
    subset_in_place(status, q);
    subset_in_place(weightn, q);
    subset_in_place(offsetn, q);
    subset_in_place(idn, q);
    subset_in_place_flatmatrix(zn, q);
  }

  // --- compute eta (linear predictor) ---
  std::vector<double> eta = offsetn; // initialize with offset
  // add contributions of each coefficient times column
  const double* zdata = zn.data_ptr();
  for (int i = 0; i < nvar; ++i) {
    double b = beta[i];
    if (b == 0.0) continue;
    const double* col = zdata + i * n1;
    for (int r = 0; r < n1; ++r) eta[r] += b * col[i];
  }

  // --- compute sigma per observation ---
  std::vector<double> s(n1, 1.0);
  if (dist_code != 1) {
    for (int i = 0; i < n1; ++i) {
      int k = stratumn[i] + nvar;
      s[i] = std::exp(beta[k]);
    }
  }

  // --- map type string to code ---
  int K = 1;
  if (type == "dfbeta" || type == "dfbetas") K = p;
  else if (type == "matrix") K = 6;

  int type_code = 0;
  if (type == "response") type_code = 1;
  else if (type == "martingale") type_code = 2;
  else if (type == "deviance") type_code = 3;
  else if (type == "working") type_code = 4;
  else if (type == "dfbeta") type_code = 5;
  else if (type == "dfbetas") type_code = 6;
  else if (type == "ldcase") type_code = 7;
  else if (type == "ldresp") type_code = 8;
  else if (type == "ldshape") type_code = 9;
  else if (type == "matrix") type_code = 10;
  else throw std::invalid_argument("invalid type of residuals: " + type);

  // --- prepare result matrix rr (n1 x K) ---
  FlatMatrix rr(n1, K);
  double* rdata = rr.data_ptr();

  // Helper lambdas for rr indexing: rr(i,k) -> rdata[k * n1 + i]

  // --- simple types ---
  if (type_code == 1) { // response
    std::vector<double> yhat0(n1);
    for (int i = 0; i < n1; ++i) {
      if (status[i] == 0 || status[i] == 1) {
        if (dist_code == 1 || dist_code == 2 || dist_code == 3 || dist_code == 5) 
          yhat0[i] = std::log(tstart[i]);
        else yhat0[i] = tstart[i];
      } else if (status[i] == 2) {
        if (dist_code == 1 || dist_code == 2 || dist_code == 3 || dist_code == 5) 
          yhat0[i] = std::log(tstop[i]);
        else yhat0[i] = tstop[i];
      } else { // interval
        if (dist_code == 1 || dist_code == 2) {
          double width = (std::log(tstop[i]) - std::log(tstart[i])) / s[i];
          yhat0[i] = std::log(tstart[i]) - s[i] * std::log(width / (std::exp(width) - 1.0));
        } else if (dist_code == 3 || dist_code == 5) {
          yhat0[i] = 0.5 * (std::log(tstart[i]) + std::log(tstop[i]));
        } else {
          yhat0[i] = 0.5 * (tstart[i] + tstop[i]);
        }
      }
      double val;
      if (dist_code == 1 || dist_code == 2 || dist_code == 3 || dist_code == 5) 
        val = std::exp(yhat0[i]) - std::exp(eta[i]);
      else val = yhat0[i] - eta[i];
      rdata[i] = val;
    }
  } else if (type_code == 2) { // martingale
    if (dist_code == 4 || dist_code == 6) 
      throw std::invalid_argument("incorrect distribution for martingale residuals: " + dist1);
    for (int i = 0; i < n1; ++i) {
      if (status[i] == 0 || status[i] == 1) {
        double y = (std::log(tstart[i]) - eta[i]) / s[i];
        double val = 0.0;
        double evt = status[i] == 1 ? 1.0 : 0.0;
        if (dist_code == 1 || dist_code == 2) val = evt - std::exp(y);
        else if (dist_code == 3) val = evt + std::log(boost_pnorm(y,0,1,0));
        else if (dist_code == 5) val = evt + std::log(boost_plogis(y,0,1,0));
        rdata[i] = val;
      } else {
        rdata[i] = NaN;
      }
    }
  } else {
    // --- complex types using f_ld_1 ---
    aftparams param = { dist_code, stratumn, tstart, tstop, status, weightn, offsetn, zn, nstrata };
    ListCpp der = f_ld_1(eta, s, &param);
    std::vector<double> g   = der.get<std::vector<double>>("g");
    std::vector<double> dg  = der.get<std::vector<double>>("dg");
    std::vector<double> ddg = der.get<std::vector<double>>("ddg");
    std::vector<double> dsv = der.get<std::vector<double>>("ds");
    std::vector<double> ddsv= der.get<std::vector<double>>("dds");
    std::vector<double> dsg = der.get<std::vector<double>>("dsg");

    if (type_code == 3) { // deviance
      std::vector<double> loglik(n1);
      for (int i = 0; i < n1; ++i) {
        switch (status[i]) {
        case 0: case 2:
          loglik[i] = 0.0;
          break;
        case 1:
          if (dist_code == 1 || dist_code == 2) loglik[i] = -std::log(s[i]) - 1.0;
          else if (dist_code == 3 || dist_code == 4) 
            loglik[i] = -std::log(std::sqrt(2.0*M_PI) * s[i]);
          else loglik[i] = -std::log(4.0 * s[i]);
          break;
        default: { // interval censored
            double width;
            if (dist_code == 1 || dist_code == 2) {
              width = (std::log(tstop[i]) - std::log(tstart[i])) / s[i];
              loglik[i] = - width/(std::exp(width)-1.0) + std::log(1.0 - std::exp(-width));
            } else if (dist_code == 3) {
              width = (std::log(tstop[i]) - std::log(tstart[i])) / s[i];
              loglik[i] = std::log(2.0 * boost_pnorm(width/2.0) - 1.0);
            } else if (dist_code == 4) {
              width = (tstop[i] - tstart[i]) / s[i];
              loglik[i] = std::log(2.0 * boost_pnorm(width/2.0) - 1.0);
            } else if (dist_code == 5) {
              width = (std::log(tstop[i]) - std::log(tstart[i])) / s[i];
              loglik[i] = std::log((std::exp(width/2.0)-1.0)/(std::exp(width/2.0)+1.0));
            } else {
              width = (tstop[i] - tstart[i]) / s[i];
              loglik[i] = std::log((std::exp(width/2.0)-1.0)/(std::exp(width/2.0)+1.0));
            }
          }
        }
        double val = -dg[i] / ddg[i];
        double dev = 0.0;
        if (std::isfinite(loglik[i]) && std::isfinite(g[i])) {
          double tmp = 2.0 * (loglik[i] - g[i]);
          dev = (tmp > 0.0) ? std::sqrt(tmp) : 0.0;
          if (val < 0) dev = -dev;
        } else dev = NaN;
        rdata[i] = dev;
      }
    } else if (type_code == 4) { // working
      for (int i = 0; i < n1; ++i) rdata[i] = -dg[i] / ddg[i];
    } else if (type_code == 5 || type_code == 6 || type_code == 7) { // dfbeta, dfbetas, ldcase
      // vbeta is p x p FlatMatrix (column-major)
      const double* vptr = vbeta.data_ptr();
      for (int i = 0; i < n1; ++i) {
        // compute score vector
        std::vector<double> score(p, 0.0);
        for (int j = 0; j < nvar; ++j) score[j] = dg[i] * zn(j, i); // zn(j,i) helper via operator()
        for (int j = nvar; j < p; ++j) score[j] = (stratumn[i] == j - nvar) ? dsv[i] : 0.0;
        // resid = score * vbeta  (1 x p) * (p x p) -> vector length p
        std::vector<double> resid(p, 0.0);
        for (int k = 0; k < p; ++k) {
          double sum = 0.0;
          const double* colvk = vptr + k * p; // vbeta(:,k)
          for (int j = 0; j < p; ++j) sum += score[j] * colvk[j];
          resid[k] = sum;
        }
        if (type_code == 6) {
          for (int k = 0; k < p; ++k) {
            double denom = vbeta(k,k);
            if (denom <= 0.0) rdata[k * n1 + i] = NaN;
            else rdata[k * n1 + i] = resid[k] / std::sqrt(denom);
          }
        } else if (type_code == 7) {
          double acc = 0.0;
          for (int k = 0; k < p; ++k) acc += resid[k] * score[k];
          rdata[i] = acc;
        } else {
          for (int k = 0; k < p; ++k) rdata[k * n1 + i] = resid[k];
        }
      }
    } else if (type_code == 8) { // ldresp
      const double* vptr = vbeta.data_ptr();
      for (int i = 0; i < n1; ++i) {
        std::vector<double> rscore(p, 0.0);
        for (int j = 0; j < nvar; ++j) rscore[j] = -ddg[i] * zn(j, i) * s[i];
        for (int j = nvar; j < p; ++j) rscore[j] = (stratumn[i] == j - nvar) ? -dsg[i] * s[i] : 0.0;
        std::vector<double> temp(p, 0.0);
        for (int k = 0; k < p; ++k) {
          const double* colvk = vptr + k * p;
          double sum = 0.0;
          for (int j = 0; j < p; ++j) sum += rscore[j] * colvk[j];
          temp[k] = sum;
        }
        double acc = 0.0;
        for (int k = 0; k < p; ++k) acc += temp[k] * rscore[k];
        rdata[0 * n1 + i] = acc;
      }
    } else if (type_code == 9) { // ldshape
      const double* vptr = vbeta.data_ptr();
      for (int i = 0; i < n1; ++i) {
        std::vector<double> sscore(p, 0.0);
        for (int j = 0; j < nvar; ++j) sscore[j] = dsg[i] * zn(j, i);
        for (int j = nvar; j < p; ++j) sscore[j] = (stratumn[i] == j - nvar) ? ddsv[i] : 0.0;
        std::vector<double> temp(p, 0.0);
        for (int k = 0; k < p; ++k) {
          const double* colvk = vptr + k * p;
          double sum = 0.0;
          for (int j = 0; j < p; ++j) sum += sscore[j] * colvk[j];
          temp[k] = sum;
        }
        double acc = 0.0;
        for (int k = 0; k < p; ++k) acc += temp[k] * sscore[k];
        rdata[0 * n1 + i] = acc;
      }
    } else if (type_code == 10) { // matrix
      for (int i = 0; i < n1; ++i) {
        rdata[i] = g[i];
        rdata[1 * n1 + i] = dg[i];
        rdata[2 * n1 + i] = ddg[i];
        rdata[3 * n1 + i] = dsv[i];
        rdata[4 * n1 + i] = ddsv[i];
        rdata[5 * n1 + i] = dsg[i];
      }
    }
  } // end complex types

  // --- apply case weights if requested ---
  if (weighted) {
    for (int k = 0; k < K; ++k) {
      double* col = rdata + k * n1;
      for (int i = 0; i < n1; ++i) {
        col[i] *= weightn[i];
      }
    }
  }

  // --- collapse by id if requested ---
  if (collapse) {
    // order by id
    std::vector<int> order = seqcpp(0, n1 - 1);
    std::sort(order.begin(), order.end(), [&](int a, int b){ return idn[a] < idn[b]; });
    std::vector<int> id1 = subset(idn, order);
    std::vector<int> idx;
    idx.push_back(0);
    for (int i = 1; i < n1; ++i) if (id1[i] != id1[i-1]) idx.push_back(i);
    int nids = static_cast<int>(idx.size());
    idx.push_back(n1);

    FlatMatrix rr1(nids, K);
    double* rr1ptr = rr1.data_ptr();
    for (int grp = 0; grp < nids; ++grp) {
      for (int k = 0; k < K; ++k) {
        double acc = 0.0;
        const double* rrcol = rdata + k * n1;
        for (int j = idx[grp]; j < idx[grp+1]; ++j) {
          acc += rrcol[order[j]];
        }
        rr1ptr[ k * nids + grp ] = acc;
      }
    }
    return rr1;
  }

  return rr;
}
