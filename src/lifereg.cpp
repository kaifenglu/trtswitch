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
  const std::vector<double>& zdata = param->z.data; // column-major: zdata[col * n + row]
  
  // compute linear predictor eta efficiently using column-major storage:
  std::vector<double> eta = offset; // initialize with offset
  // add contributions of each coefficient times column
  for (int i = 0; i < nvar; ++i) {
    double beta = par[i];
    if (beta == 0.0) continue;
    int off = i * n;
    for (int r = 0; r < n; ++r) {
      eta[r] += beta * zdata[off + r];
    }
  }
  
  // Precompute sigma and logsigma per person (exponential has sigma = 1)
  std::vector<double> sigma(n, 1.0); 
  std::vector<double> logsigma(n, 0.0); 
  if (dist_code != 1) { // all except exponential 
    for (int person = 0; person < n; ++person) { 
      int k = strata[person] + nvar; // index in par for log(sigma) 
      double s = std::exp(par[k]); 
      sigma[person] = s; 
      logsigma[person] = par[k]; // par[k] == log(s) 
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
          double zi = z(i);
          imat(i, j) += c1 * zi * zj;
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
            double zi = z(i);
            imat(i, j) += wt * zi * zj;
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
            double zi = z(i);
            imat(i, j) += c1 * zi * zj;
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
        double e_u1 = std::exp(u1), e_u2 = std::exp(u2);
        double q1 = std::exp(-e_u1), q2 = std::exp(-e_u2);
        double d1 = e_u1 * q1, d2 = e_u2 * q2;
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
            double zi = z(i);
            imat(i, j) += c1 * zi * zj;
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
      double d1 = boost_dnorm(u1), d2 = boost_dnorm(u2);
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
          double zi = z(i);
          imat(i, j) += c1 * zi * zj;
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
      double d1 = boost_dlogis(u1), d2 = boost_dlogis(u2);
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
          double zi = z(i);
          imat(i, j) += c1 * zi * zj;
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
            double zi = z(i);
            imat(i, j) += c1 * zi * zj;
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
            double zi = z(i);
            imat(i, j) += c1 * zi * zj;
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
            double zi = z(i);
            imat(i, j) += c1 * zi * zj;
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
            double zi = z(i);
            imat(i, j) += c1 * zi * zj;
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
            double zi = z(i);
            imat(i, j) += c1 * zi * zj;
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
            double zi = z(i);
            imat(i, j) += c1 * zi * zj;
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
  result.push_back(score, "score"); 
  result.push_back(imat, "imat"); 
  return result; 
}


// score residual matrix
NumericMatrix f_ressco_1(int p, const NumericVector& par, void *ex) {
  aftparams *param = (aftparams *) ex;
  int n = param->z.nrow;
  int nvar = param->z.ncol;

  // compute eta
  NumericVector eta(n);
  for (int person = 0; person < n; ++person) {
    double val = param->offset[person];
    for (int i = 0; i < nvar; ++i) val += par[i] * param->z(person, i);
    eta[person] = val;
  }

  // compute sigma
  NumericVector sig(n, 1.0);
  if (param->dist_code != 1) { // not exponential
    for (int person = 0; person < n; ++person) {
      int k = param->strata[person] + nvar;
      sig[person] = std::exp(par[k]);
    }
  }

  // Main loop to compute residuals
  NumericMatrix resid(n, p);
  for (int person = 0; person < n; ++person) {
    double sigma = sig[person];
    NumericVector z = param->z(person, _) / sigma;
    int k = param->strata[person] + nvar;
    double eta_p = eta[person];

    double u, u1, u2, c1, w1, w2, q1, q2, d1, d2;

    switch (param->status[person]) {

    case 1: // event
      switch (param->dist_code) {
      case 1: case 2: // exponential / weibull
        u = (std::log(param->tstop[person]) - eta_p) / sigma;
        c1 = -(1 - std::exp(u));
        for (int i = 0; i < nvar; ++i) resid(person, i) = c1 * z[i];
        if (param->dist_code == 2) resid(person, k) = (1 - std::exp(u)) * (-u) - 1;
        break;
      case 3: case 4: // lognormal / normal
        u = (param->dist_code == 3) ?
        (std::log(param->tstop[person]) - eta_p) / sigma :
        (param->tstop[person] - eta_p) / sigma;
        for (int i = 0; i < nvar; ++i) resid(person, i) = u * z[i];
        resid(person, k) = u * u - 1;
        break;
      case 5: case 6: // loglogistic / logistic
        u = (param->dist_code == 5) ?
        (std::log(param->tstop[person]) - eta_p) / sigma :
        (param->tstop[person] - eta_p) / sigma;
        c1 = 1 - 2 * boost_plogis(u, 0, 1, 0, 0);
        for (int i = 0; i < nvar; ++i) resid(person, i) = c1 * z[i];
        resid(person, k) = c1 * u - 1;
        break;
      }
      break;

    case 3: // interval censoring
      switch (param->dist_code) {
      case 1: case 2: // exponential / weibull
        u1 = (std::log(param->tstart[person]) - eta_p) / sigma;
        u2 = (std::log(param->tstop[person]) - eta_p) / sigma;
        w1 = std::exp(u1); w2 = std::exp(u2);
        q1 = std::exp(-w1); q2 = std::exp(-w2);
        d1 = w1 * q1; d2 = w2 * q2;
        c1 = (d1 - d2) / (q1 - q2);
        for (int i = 0; i < nvar; ++i) resid(person, i) = c1 * z[i];
        if (param->dist_code == 2)
          resid(person, k) = (d1 * u1 - d2 * u2) / (q1 - q2);
        break;
      case 3: case 4: // lognormal / normal
        u1 = (param->dist_code == 3) ?
        (std::log(param->tstart[person]) - eta_p) / sigma :
        (param->tstart[person] - eta_p) / sigma;
        u2 = (param->dist_code == 3) ?
        (std::log(param->tstop[person]) - eta_p) / sigma :
          (param->tstop[person] - eta_p) / sigma;
        d1 = boost_dnorm(u1); d2 = boost_dnorm(u2);
        q1 = boost_pnorm(u1, 0, 1, 0, 0); q2 = boost_pnorm(u2, 0, 1, 0, 0);
        c1 = (d1 - d2) / (q1 - q2);
        for (int i = 0; i < nvar; ++i) resid(person, i) = c1 * z[i];
        resid(person, k) = (d1 * u1 - d2 * u2) / (q1 - q2);
        break;
      case 5: case 6: // loglogistic / logistic
        u1 = (param->dist_code == 5) ?
        (std::log(param->tstart[person]) - eta_p) / sigma :
        (param->tstart[person] - eta_p) / sigma;
        u2 = (param->dist_code == 5) ?
        (std::log(param->tstop[person]) - eta_p) / sigma :
          (param->tstop[person] - eta_p) / sigma;
        d1 = boost_dlogis(u1); d2 = boost_dlogis(u2);
        q1 = boost_plogis(u1, 0, 1, 0, 0); q2 = boost_plogis(u2, 0, 1, 0, 0);
        c1 = (d1 - d2) / (q1 - q2);
        for (int i = 0; i < nvar; ++i) resid(person, i) = c1 * z[i];
        resid(person, k) = (d1 * u1 - d2 * u2) / (q1 - q2);
        break;
      }
      break;

    case 2: // left censoring
      switch (param->dist_code) {
      case 1: case 2: // exponential / weibull
        u = (std::log(param->tstop[person]) - eta_p) / sigma;
        w2 = std::exp(u); q2 = std::exp(-w2); d2 = w2 * q2;
        c1 = -d2 / (1 - q2);
        for (int i = 0; i < nvar; ++i) resid(person, i) = c1 * z[i];
        if (param->dist_code == 2) resid(person, k) = c1 * u;
        break;
      case 3: case 4: // lognormal / normal
        u = (param->dist_code == 3) ?
        (std::log(param->tstop[person]) - eta_p) / sigma :
        (param->tstop[person] - eta_p) / sigma;
        c1 = -boost_dnorm(u) / boost_pnorm(u);
        for (int i = 0; i < nvar; ++i) resid(person, i) = c1 * z[i];
        resid(person, k) = c1 * u;
        break;
      case 5: case 6: // loglogistic / logistic
        u = (param->dist_code == 5) ?
        (std::log(param->tstop[person]) - eta_p) / sigma :
        (param->tstop[person] - eta_p) / sigma;
        c1 = -boost_plogis(u, 0, 1, 0, 0);
        for (int i = 0; i < nvar; ++i) resid(person, i) = c1 * z[i];
        resid(person, k) = c1 * u;
        break;
      }
      break;

    case 0: // right censoring
      switch (param->dist_code) {
      case 1: case 2: // exponential / weibull
        u = (std::log(param->tstart[person]) - eta_p) / sigma;
        c1 = std::exp(u);
        for (int i = 0; i < nvar; ++i) resid(person, i) = c1 * z[i];
        if (param->dist_code == 2) resid(person, k) = c1 * u;
        break;
      case 3: case 4: // lognormal / normal
        u = (param->dist_code == 3) ?
        (std::log(param->tstart[person]) - eta_p) / sigma :
        (param->tstart[person] - eta_p) / sigma;
        c1 = boost_dnorm(u) / boost_pnorm(u, 0, 1, 0, 0);
        for (int i = 0; i < nvar; ++i) resid(person, i) = c1 * z[i];
        resid(person, k) = c1 * u;
        break;
      case 5: case 6: // loglogistic / logistic
        u = (param->dist_code == 5) ?
        (std::log(param->tstart[person]) - eta_p) / sigma :
        (param->tstart[person] - eta_p) / sigma;
        c1 = boost_plogis(u);
        for (int i = 0; i < nvar; ++i) resid(person, i) = c1 * z[i];
        resid(person, k) = c1 * u;
        break;
      }
      break;

    default:
      throw std::runtime_error("Unknown status: " + std::to_string(param->status[person]));
    }
  }

  return resid;
}

// 
// // substitute information matrix guaranteed to be positive definite
// NumericMatrix f_jj_1(int p, const NumericVector& par, void *ex) {
//   aftparams *param = (aftparams *) ex;
//   int n = param->z.nrow();
//   
//   NumericMatrix resid = f_ressco_1(p, par, param);
//   NumericMatrix jj(p,p);
//   for (int person = 0; person < n; ++person) {
//     double w = param->weight[person];
//     for (int i = 0; i < p; ++i) {
//       for (int j = 0; j < p; ++j) {
//         jj(i,j) += w * resid(person,i) * resid(person,j);
//       }
//     }
//   }
//   
//   return jj;
// }
// 
// 
// // underlying optimization algorithm for lifereg
// List liferegloop(int p, const NumericVector& par, void *ex,
//                  int maxiter, double eps,
//                  const IntegerVector& colfit, int ncolfit) {
//   aftparams *param = (aftparams *) ex;
//   
//   int iter, halving = 0;
//   bool fail = false;
//   double toler = 1e-12;
//   
//   int nstrata = param->nstrata;
//   int nvar = param->z.ncol();
//   int nsub = param->z.nrow();
//   
//   NumericMatrix z1 = param->z;
//   NumericVector mu(nvar), sigma(nvar);
//   NumericMatrix z2(nsub, nvar);
//   
//   // --- standardize z once ---
//   for (int i = 0; i < nvar; ++i) {
//     NumericVector u = z1(_, i);
//     double s = sd(u), m = mean(u);
//     if (is_true(all((u == 0) | (u == 1)))) {
//       m = 0; s = 1;
//     }
//     mu[i] = m;
//     sigma[i] = s;
//     for (int k = 0; k < nsub; ++k) z2(k, i) = (u[k] - m) / s;
//   }
//   
//   // --- initial beta ---
//   NumericVector beta(p), newbeta(p);
//   beta[0] = par[0];
//   for (int i = 1; i < nvar; ++i) {
//     beta[i] = par[i] * sigma[i];
//     beta[0] += par[i] * mu[i];
//   }
//   if (param->dist_code != 1)
//     std::copy(par.begin() + nvar, par.end(), beta.begin() + nvar);
//   
//   aftparams para = {param->dist_code, param->strata, param->tstart,
//                     param->tstop, param->status, param->weight,
//                     param->offset, z2, nstrata};
//   
//   
//   List der = f_der_1(p, beta, &para);
//   double loglik = der["loglik"];
//   double newlk = 0;
//   NumericVector u = der["score"];
//   NumericMatrix imat = as<NumericMatrix>(der["imat"]);
//   NumericMatrix jj;
//   NumericVector u1(ncolfit);
//   NumericMatrix imat1(ncolfit, ncolfit);
//   NumericMatrix jj1(ncolfit, ncolfit);
//   
//   for (int i = 0; i < ncolfit; ++i) u1[i] = u[colfit[i]];
//   
//   for (int i = 0; i < ncolfit; ++i)
//     for (int j = 0; j < ncolfit; ++j)
//       imat1(i, j) = imat(colfit[i], colfit[j]);
//   
//   // --- first step ---
//   if (cholesky2(imat1, ncolfit, toler) < 0) {
//     jj = f_jj_1(p, beta, &para);
//     for (int i = 0; i < ncolfit; ++i)
//       for (int j = 0; j < ncolfit; ++j)
//         jj1(i, j) = jj(colfit[i], colfit[j]);
//     cholesky2(jj1, ncolfit, toler);
//     chsolve2(jj1, ncolfit, u1);
//   } else {
//     chsolve2(imat1, ncolfit, u1);
//   }
//   
//   u.fill(0.0);
//   for (int i = 0; i < ncolfit; ++i) u[colfit[i]] = u1[i];
//   newbeta = beta + u;
//   
//   // --- main iteration ---
//   for (iter = 0; iter < maxiter; ++iter) {
//     der = f_der_1(p, newbeta, &para);
//     newlk = der["loglik"];
//     
//     fail = std::isnan(newlk) || std::isinf(newlk);
//     if (!fail && halving == 0 && fabs(1 - loglik / newlk) < eps) break;
//     
//     if (fail || newlk < loglik) {
//       ++halving;
//       for (int i = 0; i < p; ++i)
//         newbeta[i] = 0.5 * (beta[i] + newbeta[i]);
//       
//       // special handling of sigmas
//       if (halving == 1 && param->dist_code != 1) {
//         for (int i = 0; i < nstrata; ++i)
//           if (beta[nvar + i] - newbeta[nvar + i] > 1.1)
//             newbeta[nvar + i] = beta[nvar + i] - 1.1;
//       }
//       continue;
//     }
//     
//     // --- update ---
//     halving = 0;
//     beta = clone(newbeta);
//     loglik = newlk;
//     u = der["score"];
//     imat = as<NumericMatrix>(der["imat"]);
//     
//     for (int i = 0; i < ncolfit; ++i) u1[i] = u[colfit[i]];
//     
//     for (int i = 0; i < ncolfit; ++i)
//       for (int j = 0; j < ncolfit; ++j)
//         imat1(i, j) = imat(colfit[i], colfit[j]);
//     
//     if (cholesky2(imat1, ncolfit, toler) < 0) {
//       jj = f_jj_1(p, beta, &para);
//       for (int i = 0; i < ncolfit; ++i)
//         for (int j = 0; j < ncolfit; ++j)
//           jj1(i, j) = jj(colfit[i], colfit[j]);
//       cholesky2(jj1, ncolfit, toler);
//       chsolve2(jj1, ncolfit, u1);
//     } else {
//       chsolve2(imat1, ncolfit, u1);
//     }
//     
//     u.fill(0.0);
//     for (int i = 0; i < ncolfit; ++i) u[colfit[i]] = u1[i];
//     newbeta = beta + u;
//   }
//   
//   if (iter == maxiter) fail = true;
//   
//   // --- rescale back ---
//   for (int i = 1; i < nvar; ++i) {
//     newbeta[i] /= sigma[i];
//     newbeta[0] -= newbeta[i] * mu[i];
//   }
//   
//   
//   // rescale the information matrix accordingly
//   imat = as<NumericMatrix>(der["imat"]);
//   NumericMatrix jmat = clone(imat);
//   for (int i = 0; i < nvar; ++i)
//     for (int j = 0; j < nvar; ++j)
//       imat(i,j) = jmat(0,0)*mu[i]*mu[j] + jmat(0,j)*mu[i]*sigma[j] +
//         jmat(i,0)*mu[j]*sigma[i] + jmat(i,j)*sigma[i]*sigma[j];
//   
//   for (int i = nvar; i < p; ++i) {
//     for (int j = 0; j < nvar; ++j) {
//       imat(i,j) = jmat(i,0)*mu[j] + jmat(i,j)*sigma[j];
//       imat(j,i) = imat(i,j);
//     }
//   }
//   
//   // compute variance matrix
//   for (int i = 0; i < ncolfit; ++i)
//     for (int j = 0; j < ncolfit; ++j)
//       imat1(i, j) = imat(colfit[i], colfit[j]);
//   
//   NumericMatrix var1 = invsympd(imat1, ncolfit, toler);
//   NumericMatrix var(p, p);
//   for (int i = 0; i < ncolfit; ++i)
//     for (int j = 0; j < ncolfit; ++j)
//       var(colfit[i], colfit[j]) = var1(i, j);
//   
//   return List::create(
//     Named("coef") = newbeta,
//     Named("iter") = iter,
//     Named("var") = var,
//     Named("loglik") = newlk,
//     Named("fail") = fail);
// }
// 
// 
// // confidence limit of profile likelihood method
// double liferegplloop(int p, const NumericVector& par, void *ex,
//                      int maxiter, double eps,
//                      int k, int which, double l0) {
//   aftparams *param = (aftparams *) ex;
//   
//   int iter;
//   bool fail = false;
//   double toler = 1e-12;
//   
//   NumericVector beta(p), newbeta(p);
//   double loglik, newlk;
//   NumericVector u(p), delta(p);
//   NumericMatrix imat(p, p), jj(p, p), v(p, p);
//   
//   // --- first step ---
//   beta = clone(par);
//   
//   List der = f_der_1(p, beta, param);
//   loglik = der["loglik"];
//   u = der["score"];
//   imat = as<NumericMatrix>(der["imat"]);
//   
//   // compute inverse of inv (with cholesky fallback to jj)
//   jj = clone(imat); // test cholesky on a copy (cholesky2 overwrites)
//   
//   if (cholesky2(jj, p, toler) < 0) {
//     jj = f_jj_1(p, beta, param);
//     v = invsympd(jj, p, toler); // inv of jj
//   } else {
//     v = invsympd(imat, p, toler); // inv of imat
//   }
//   
//   // Lagrange multiplier method as used in SAS PROC LOGISTIC
//   double w = 0;
//   for (int i=0; i<p; ++i)
//     for (int j=0; j<p; ++j)
//       w -= u[i]*v(i,j)*u[j];
//   
//   double underroot = -2*(l0 - loglik + 0.5*w)/v(k,k);
//   double lambda = underroot < 0.0 ? 0.0 : which*std::sqrt(underroot);
//   u[k] += lambda;
//   
//   delta.fill(0.0);
//   for (int i=0; i<p; ++i)
//     for (int j=0; j<p; ++j)
//       delta[i] += v(i,j)*u[j];
//   
//   // update beta
//   newbeta = beta + delta;
//   
//   // --- main iteration ---
//   for (iter=0; iter<maxiter; ++iter) {
//     der = f_der_1(p, newbeta, param);
//     newlk = der["loglik"];
//     
//     // check convergence
//     fail = std::isnan(newlk) || std::isinf(newlk);
//     if (!fail && fabs(newlk - l0) < eps && w < eps) break;
//     
//     // update step
//     beta = clone(newbeta);
//     loglik = newlk;
//     u = der["score"];
//     imat = as<NumericMatrix>(der["imat"]);
//     jj = clone(imat);
//     
//     if (cholesky2(jj, p, toler) < 0) {
//       jj = f_jj_1(p, beta, param);
//       v = invsympd(jj, p, toler);
//     } else {
//       v = invsympd(imat, p, toler);
//     }
//     
//     // Lagrange multiplier method as used in SAS PROC LOGISTIC
//     double w = 0;
//     for (int i=0; i<p; ++i)
//       for (int j=0; j<p; ++j)
//         w -= u[i]*v(i,j)*u[j];
//     
//     double underroot = -2*(l0 - loglik + 0.5*w)/v(k,k);
//     double lambda = underroot < 0.0 ? 0.0 : which*std::sqrt(underroot);
//     u[k] += lambda;
//     
//     delta.fill(0.0);
//     for (int i=0; i<p; ++i)
//       for (int j=0; j<p; ++j)
//         delta[i] += v(i,j)*u[j];
//     
//     // update beta
//     newbeta = beta + delta;
//   }
//   
//   if (iter == maxiter) fail = true;
//   if (fail) warning("The algorithm in liferegplloop did not converge");
//   
//   return newbeta[k];
// }
// 
// 
// 
// 
// 
// // [[Rcpp::export]]
// List liferegcpp(const DataFrame data,
//                 const StringVector& rep = "",
//                 const StringVector& stratum = "",
//                 const std::string time = "time",
//                 const std::string time2 = "",
//                 const std::string event = "event",
//                 const StringVector& covariates = "",
//                 const std::string weight = "",
//                 const std::string offset = "",
//                 const std::string id = "",
//                 const std::string dist = "weibull",
//                 const NumericVector& init = NA_REAL,
//                 const bool robust = false,
//                 const bool plci = false,
//                 const double alpha = 0.05,
//                 const int maxiter = 50,
//                 const double eps = 1.0e-9) {
//   
//   int n = data.nrows();
//   int nvar = static_cast<int>(covariates.size()) + 1;
//   if (nvar == 2 && (covariates[0] == "" || covariates[0] == "none")) {
//     nvar = 1;
//   }
//   
//   std::string dist1 = dist;
//   std::for_each(dist1.begin(), dist1.end(), [](char & c) {
//     c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
//   });
//   
//   if (dist1 == "log-logistic" || dist1 == "llogistic") {
//     dist1 = "loglogistic";
//   } else if (dist1 == "log-normal" || dist1 == "lnormal") {
//     dist1 = "lognormal";
//   } else if (dist1 == "gaussian") {
//     dist1 = "normal";
//   }
//   
//   int dist_code;
//   if (dist1 == "exponential") dist_code = 1;
//   else if (dist1 == "weibull") dist_code = 2;
//   else if (dist1 == "lognormal") dist_code = 3;
//   else if (dist1 == "normal") dist_code = 4;
//   else if (dist1 == "loglogistic") dist_code = 5;
//   else if (dist1 == "logistic") dist_code = 6;
//   else throw std::invalid_argument("invalid distribution: " + dist1);
//   
//   bool has_rep;
//   IntegerVector repn(n);
//   DataFrame u_rep;
//   int p_rep = static_cast<int>(rep.size());
//   if (p_rep == 1 && (rep[0] == "" || rep[0] == "none")) {
//     has_rep = 0;
//     repn.fill(0);
//   } else {
//     List out = bygroup(data, rep);
//     has_rep = 1;
//     repn = out["index"];
//     u_rep = DataFrame(out["lookup"]);
//   }
//   
//   IntegerVector stratumn(n);
//   int p_stratum = static_cast<int>(stratum.size());
//   if (p_stratum == 1 && (stratum[0] == "" || stratum[0] == "none")) {
//     stratumn.fill(0);
//   } else {
//     List out = bygroup(data, stratum);
//     stratumn = out["index"];
//   }
//   
//   IntegerVector stratumn1 = unique(stratumn);
//   int nstrata = static_cast<int>(stratumn1.size());
//   int p = dist_code == 1 ? nvar : (nvar+nstrata);
//   
//   if (dist_code == 1 && nstrata > 1) {
//     throw std::invalid_argument("Stratification is not valid with the exponential distribution");
//   }
//   
//   bool has_time = hasVariable(data, time);
//   if (!has_time) throw std::invalid_argument("data must contain the time variable");
//   NumericVector timenz = data[time];
//   NumericVector timen = clone(timenz);
//   for (int i=0; i<n; ++i) {
//     if (!std::isnan(timen[i]) && (dist_code == 1 || dist_code == 2 ||
//         dist_code == 3 || dist_code == 5) && timen[i] <= 0) {
//       std::string str1 = "time must be positive for each subject for the";
//       std::string str2 = "distribution";
//       std::string errmsg = str1 + " " + dist1 + " " + str2;
//       throw std::invalid_argument(errmsg);
//     }
//   }
//   
//   bool has_time2 = hasVariable(data, time2);
//   NumericVector time2n(n);
//   if (has_time2) {
//     NumericVector time2nz = data[time2];
//     time2n = clone(time2nz);
//     for (int i=0; i<n; ++i) {
//       if (!std::isnan(time2n[i]) && (dist_code == 1 || dist_code == 2 ||
//           dist_code == 3 || dist_code == 5) && time2n[i] <= 0) {
//         std::string str1 = "time2 must be positive for each subject for the";
//         std::string str2 = "distribution";
//         std::string errmsg = str1 + " " + dist1 + " " + str2;
//         throw std::invalid_argument(errmsg);
//       }
//     }
//   }
//   
//   bool has_event = hasVariable(data, event);
//   if (!has_time2 && !has_event) {
//     throw std::invalid_argument("data must contain the event variable for right censored data");
//   }
//   
//   IntegerVector eventn(n);
//   if (has_event) {
//     IntegerVector eventnz = data[event];
//     eventn = clone(eventnz);
//     if (is_true(any((eventn != 1) & (eventn != 0)))) {
//       throw std::invalid_argument("event must be 1 or 0 for each subject");
//     }
//   }
//   
//   NumericMatrix zn(n,nvar);
//   for (int i=0; i<n; ++i) {
//     zn(i,0) = 1; // intercept
//   }
//   
//   for (int j=0; j<nvar-1; ++j) {
//     String zj = covariates[j];
//     if (!hasVariable(data, zj)) {
//       throw std::invalid_argument("data must contain the variables in covariates");
//     }
//     NumericVector u = data[zj];
//     for (int i=0; i<n; ++i) {
//       zn(i,j+1) = u[i];
//     }
//   }
//   
//   bool has_weight = hasVariable(data, weight);
//   NumericVector weightn(n, 1.0);
//   if (has_weight) {
//     NumericVector weightnz = data[weight];
//     weightn = clone(weightnz);
//     if (is_true(any(weightn <= 0))) {
//       throw std::invalid_argument("weight must be greater than 0");
//     }
//   }
//   
//   bool has_offset = hasVariable(data, offset);
//   NumericVector offsetn(n);
//   if (has_offset) {
//     NumericVector offsetnz = data[offset];
//     offsetn = clone(offsetnz);
//   }
//   
//   // create the numeric id variable
//   bool has_id = hasVariable(data, id);
//   IntegerVector idn(n);
//   if (!has_id) {
//     idn = seq(0, n - 1);
//   } else {
//     SEXP col = data[id];
//     SEXPTYPE col_type = TYPEOF(col);
//     if (col_type == INTSXP) {
//       IntegerVector v = col;
//       IntegerVector w = unique(v);
//       w.sort();
//       idn = match(v, w) - 1;
//     } else if (col_type == REALSXP) {
//       NumericVector v = col;
//       NumericVector w = unique(v);
//       w.sort();
//       idn = match(v, w) - 1;
//     } else if (col_type == STRSXP) {
//       StringVector v = col;
//       StringVector w = unique(v);
//       w.sort();
//       idn = match(v, w) - 1;
//     } else {
//       throw std::invalid_argument("incorrect type for the id variable in the input data");
//     }
//   }
//   
//   // sort the data by rep
//   if (has_rep) {
//     IntegerVector order = seq(0, n-1);
//     std::stable_sort(order.begin(), order.end(), [&](int i, int j) {
//       return repn[i] < repn[j];
//     });
//     
//     repn = repn[order];
//     stratumn = stratumn[order];
//     timen = timen[order];
//     time2n = time2n[order];
//     eventn = eventn[order];
//     weightn = weightn[order];
//     offsetn = offsetn[order];
//     idn = idn[order];
//     zn = subset_matrix_by_row(zn, order);
//   }
//   
//   // exclude observations with missing values
//   LogicalVector sub(n,1);
//   for (int i=0; i<n; ++i) {
//     if (repn[i] == NA_INTEGER || stratumn[i] == NA_INTEGER ||
//         (std::isnan(timen[i]) && std::isnan(time2n[i])) ||
//         std::isnan(weightn[i]) || std::isnan(offsetn[i]) ||
//         idn[i] == NA_INTEGER) {
//       sub[i] = 0;
//     }
//     for (int j=0; j<nvar-1; ++j) {
//       if (std::isnan(zn(i,j+1))) sub[i] = 0;
//     }
//   }
//   
//   IntegerVector order = which(sub);
//   repn = repn[order];
//   stratumn = stratumn[order];
//   timen = timen[order];
//   time2n = time2n[order];
//   eventn = eventn[order];
//   weightn = weightn[order];
//   offsetn = offsetn[order];
//   idn = idn[order];
//   zn = subset_matrix_by_row(zn, order);
//   n = sum(sub);
//   if (n == 0) throw std::invalid_argument("no observation is left after removing missing values");
//   
//   // identify the locations of the unique values of rep
//   IntegerVector idx(1,0);
//   for (int i=1; i<n; ++i) {
//     if (repn[i] != repn[i-1]) {
//       idx.push_back(i);
//     }
//   }
//   
//   int nreps = static_cast<int>(idx.size());
//   idx.push_back(n);
//   
//   // variables in the output data sets
//   // sumstat data set
//   IntegerVector rep01 = seq(0,nreps-1);
//   IntegerVector nobs(nreps), nevents(nreps);
//   NumericVector loglik0(nreps), loglik1(nreps);
//   IntegerVector niter(nreps);
//   LogicalVector fails(nreps);
//   
//   // parest data set
//   IntegerVector rep0(nreps*p);
//   StringVector par0(nreps*p);
//   NumericVector beta0(nreps*p), sebeta0(nreps*p), rsebeta0(nreps*p);
//   NumericVector z0(nreps*p), expbeta0(nreps*p);
//   NumericMatrix vbeta0(nreps*p,p), rvbeta0(nreps*p,p);
//   NumericVector lb0(nreps*p), ub0(nreps*p), prob0(nreps*p);
//   StringVector clparm0(nreps*p);
//   
//   NumericVector linear_predictors(n);
//   
//   double zcrit = boost_qnorm(1-alpha/2);
//   double xcrit = zcrit * zcrit;
//   int bign0 = 0;
//   for (int h=0; h<nreps; ++h) {
//     IntegerVector q1 = Range(idx[h], idx[h+1]-1);
//     int n1 = static_cast<int>(q1.size());
//     
//     IntegerVector stratum1 = stratumn[q1];
//     NumericVector time1 = timen[q1];
//     NumericVector time21 = time2n[q1];
//     IntegerVector event1 = eventn[q1];
//     NumericVector weight1 = weightn[q1];
//     NumericVector offset1 = offsetn[q1];
//     IntegerVector id1 = idn[q1];
//     NumericMatrix z1 = subset_matrix_by_row(zn, q1);
//     
//     // unify right censored data with interval censored data
//     NumericVector tstart(n1), tstop(n1);
//     if (!has_time2) {
//       tstart = time1;
//       for (int i=0; i<n1; ++i) {
//         tstop[i] = event1[i] == 1 ? tstart[i] : NA_REAL;
//       }
//     } else {
//       tstart = time1;
//       tstop = time21;
//     }
//     
//     IntegerVector status(n1);
//     for (int i=0; i<n1; ++i) {
//       if (!std::isnan(tstart[i]) && !std::isnan(tstop[i]) &&
//           tstart[i] == tstop[i]) {
//         status[i] = 1; // event
//       } else if (!std::isnan(tstart[i]) && !std::isnan(tstop[i]) &&
//         tstart[i] < tstop[i]) {
//         status[i] = 3; // interval censoring
//       } else if (std::isnan(tstart[i]) && !std::isnan(tstop[i])) {
//         status[i] = 2; // left censoring
//       } else if (!std::isnan(tstart[i]) && std::isnan(tstop[i])) {
//         status[i] = 0; // right censoring
//       } else {
//         status[i] = -1; // exclude the observation
//       }
//     }
//     
//     nobs[h] = n1;
//     nevents[h] = sum(status == 1);
//     
//     if (nevents[h] == 0) {
//       for (int i=0; i<p; ++i) {
//         int k = h*p+i;
//         rep0[k] = h;
//         
//         if (i==0) {
//           par0[k] = "(Intercept)";
//         } else if (i < nvar) {
//           par0[k] = covariates[i-1];
//         } else {
//           if (nstrata == 1) {
//             par0[k] = "Log(scale)";
//           } else {
//             std::string str1 = "Log(scale ";
//             std::string str2 = ")";
//             par0[k] = str1 + std::to_string(i-nvar+1) + str2;
//           }
//         }
//         
//         beta0[k] = NA_REAL;
//         sebeta0[k] = 0;
//         rsebeta0[k] = 0;
//         z0[k] = NA_REAL;
//         expbeta0[k] = NA_REAL;
//         for (int j=0; j<p; ++j) {
//           vbeta0(k,j) = 0;
//           rvbeta0(k,j) = 0;
//         }
//         lb0[k] = NA_REAL;
//         ub0[k] = NA_REAL;
//         prob0[k] = NA_REAL;
//         clparm0[k] = "Wald";
//       }
//       
//       for (int i=0; i<n1; ++i) {
//         linear_predictors[bign0+i] = offset1[i];
//       }
//       
//       bign0 += n1;
//       
//       continue;
//     }
//     
//     // exclude records with invalid status
//     IntegerVector q2 = which(status != -1);
//     int n2 = static_cast<int>(q2.size());
//     
//     if (n2 < n1) {
//       stratum1 = stratum1[q2];
//       tstart = tstart[q2];
//       tstop = tstop[q2];
//       status = status[q2];
//       weight1 = weight1[q2];
//       offset1 = offset1[q2];
//       id1 = id1[q2];
//       z1 = subset_matrix_by_row(z1,q2);
//     }
//     
//     // intercept only model
//     NumericVector time0(n2);
//     for (int i=0; i<n2; ++i) {
//       if (status[i] == 0 || status[i] == 1) { // right censoring or event
//         time0[i] = tstart[i];
//       } else if (status[i] == 2) { // left censoring
//         time0[i] = tstop[i];
//       } else if (status[i] == 3) { // interval censoring
//         time0[i] = (tstart[i] + tstop[i])/2;
//       }
//     }
//     
//     NumericVector y0 = clone(time0);
//     if (dist_code == 1 || dist_code == 2 ||
//         dist_code == 3 || dist_code == 5) {
//       y0 = std::log(y0);
//     }
//     
//     double int0 = mean(y0);
//     double logsig0 = std::log(sd(y0));
//     
//     NumericVector bint0(p);
//     int ncolfit0 = dist_code == 1 ? 1 : nstrata + 1;
//     IntegerVector colfit0(ncolfit0);
//     if (dist_code == 1) {
//       bint0[0] = int0;
//       ncolfit0 = 1;
//       colfit0[0] = 0;
//     } else {
//       bint0[0] = int0;
//       for (int i=0; i<nstrata; ++i) {
//         bint0[nvar+i] = logsig0;
//       }
//       
//       colfit0[0] = 0;
//       for (int i=0; i<nstrata; ++i) {
//         colfit0[i+1] = nvar+i;
//       }
//     }
//     
//     // parameter estimates and standard errors for the null model
//     aftparams param = {dist_code, stratum1, tstart, tstop, status, weight1,
//                        offset1, z1, nstrata};
//     
//     List outint = liferegloop(p, bint0, &param, maxiter, eps,
//                               colfit0, ncolfit0);
//     NumericVector bint = outint["coef"];
//     NumericMatrix vbint = as<NumericMatrix>(outint["var"]);
//     
//     NumericVector b(p);
//     NumericMatrix vb(p,p);
//     List out;
//     
//     if (nvar > 1) {
//       IntegerVector colfit = seq(0,p-1);
//       if (is_false(any(is_na(init))) && init.size() == p) {
//         out = liferegloop(p, init, &param, maxiter, eps, colfit, p);
//       } else {
//         out = liferegloop(p, bint, &param, maxiter, eps, colfit, p);
//       }
//       
//       bool fail = out["fail"];
//       if (fail) {
//         // obtain initial values for model parameters using OLS
//         NumericVector y1 = y0 - offset1;
//         NumericMatrix v1(nvar,nvar);  // XWX matrix
//         NumericVector u1(nvar);       // XWY vector
//         for (int i=0; i<n2; ++i) {
//           for (int j=0; j<nvar; ++j) {
//             for (int k=0; k<nvar; ++k) {
//               v1(j,k) += weight1[i]*z1(i,j)*z1(i,k);
//             }
//             u1[j] += weight1[i]*z1(i,j)*y1[i];
//           }
//         }
//         
//         double toler = 1e-12;
//         cholesky2(v1, nvar, toler);
//         chsolve2(v1, nvar, u1);
//         
//         NumericVector binit(p);
//         for (int j=0; j<nvar; ++j) {
//           binit[j] = u1[j];
//         }
//         
//         if (dist_code != 1) {
//           double s = 0;
//           for (int i=0; i<n2; ++i) {
//             double r = y1[i] - std::inner_product(
//               z1(i, _).begin(), z1(i, _).end(), u1.begin(), 0.0);
//             s += weight1[i]*r*r;
//           }
//           s = 0.5*std::log(s/sum(weight1)*n2/(n2-nvar));  // log(sigma)
//           
//           for (int j=nvar; j<p; ++j) {
//             binit[j] = s;
//           }
//         }
//         
//         // fit the model using the initial values
//         out = liferegloop(p, binit, &param, maxiter, eps, colfit, p);
//         fail = out["fail"];
//       }
//       
//       if (fail) warning("The algorithm in liferegr did not converge");
//       
//       b = out["coef"];
//       vb = as<NumericMatrix>(out["var"]);
//     } else {
//       b = bint;
//       vb = vbint;
//       out = outint;
//     }
//     
//     NumericVector seb(p);
//     for (int j=0; j<p; ++j) {
//       seb[j] = std::sqrt(vb(j,j));
//     }
//     
//     for (int i=0; i<p; ++i) {
//       rep0[h*p+i] = h;
//       
//       if (i==0) {
//         par0[h*p+i] = "(Intercept)";
//       } else if (i < nvar) {
//         par0[h*p+i] = covariates[i-1];
//       } else {
//         if (nstrata == 1) {
//           par0[h*p+i] = "Log(scale)";
//         } else {
//           std::string str1 = "Log(scale ";
//           std::string str2 = ")";
//           par0[h*p+i] = str1 + std::to_string(i-nvar+1) + str2;
//         }
//       }
//       
//       beta0[h*p+i] = b[i];
//       sebeta0[h*p+i] = seb[i];
//       for (int j=0; j<p; ++j) {
//         vbeta0(h*p+i,j) = vb(i,j);
//       }
//     }
//     
//     for (int i=0; i<n2; ++i) {
//       linear_predictors[bign0+q2[i]] = offset1[i];
//       for (int j=0; j<nvar; ++j) {
//         linear_predictors[bign0+q2[i]] += b[j]*z1(i,j);
//       }
//     }
//     
//     bign0 += n1;
//     
//     niter[h] = out["iter"];
//     fails[h] = out["fail"];
//     
//     // robust variance estimates
//     NumericVector rseb(p);  // robust standard error for betahat
//     if (robust) {
//       NumericMatrix ressco = f_ressco_1(p, b, &param);
//       
//       int nr; // number of rows in the score residual matrix
//       if (!has_id) {
//         for (int i=0; i<n2; ++i) {
//           for (int j=0; j<p; ++j) {
//             ressco(i,j) = weight1[i]*ressco(i,j);
//           }
//         }
//         nr = n2;
//       } else { // need to sum up score residuals by id
//         IntegerVector order = seq(0, n2-1);
//         std::sort(order.begin(), order.end(), [&](int i, int j) {
//           return id1[i] < id1[j];
//         });
//         
//         IntegerVector id2 = id1[order];
//         IntegerVector idx(1,0);
//         for (int i=1; i<n2; ++i) {
//           if (id2[i] != id2[i-1]) {
//             idx.push_back(i);
//           }
//         }
//         
//         int nids = static_cast<int>(idx.size());
//         idx.push_back(n2);
//         
//         NumericVector weight2 = weight1[order];
//         
//         NumericMatrix ressco2(nids,p);
//         for (int i=0; i<nids; ++i) {
//           for (int j=0; j<p; ++j) {
//             for (int k=idx[i]; k<idx[i+1]; ++k) {
//               ressco2(i,j) += weight2[k]*ressco(order[k],j);
//             }
//           }
//         }
//         
//         ressco = ressco2;  // update the score residuals
//         nr = nids;
//       }
//       
//       NumericMatrix D(nr,p); // DFBETA
//       for (int i=0; i<nr; ++i) {
//         for (int j=0; j<p; ++j) {
//           for (int k=0; k<p; ++k) {
//             D(i,j) += ressco(i,k)*vb(k,j);
//           }
//         }
//       }
//       
//       NumericMatrix rvb(p,p); // robust variance matrix for betahat
//       for (int j=0; j<p; ++j) {
//         for (int k=0; k<p; ++k) {
//           for (int i=0; i<nr; ++i) {
//             rvb(j,k) += D(i,j)*D(i,k);
//           }
//         }
//       }
//       
//       for (int i=0; i<p; ++i) {
//         rseb[i] = std::sqrt(rvb(i,i));
//       }
//       
//       for (int i=0; i<p; ++i) {
//         rsebeta0[h*p+i] = rseb[i];
//         for (int j=0; j<p; ++j) {
//           rvbeta0(h*p+i,j) = rvb(i,j);
//         }
//       }
//     }
//     
//     // profile likelihood confidence interval for regression coefficients
//     NumericVector lb(p), ub(p), prob(p);
//     StringVector clparm(p);
//     
//     if (plci) {
//       double lmax = out["loglik"];
//       double l0 = lmax - 0.5*xcrit;
//       
//       for (int k=0; k<p; ++k) {
//         lb[k] = liferegplloop(p, b, &param, maxiter, eps, k, -1, l0);
//         ub[k] = liferegplloop(p, b, &param, maxiter, eps, k, 1, l0);
//         
//         IntegerVector colfit1(p-1);
//         for (int i=0; i<k; ++i) {
//           colfit1[i] = i;
//         }
//         for (int i=k+1; i<p; ++i) {
//           colfit1[i-1] = i;
//         }
//         
//         NumericVector b0(p);
//         List out0 = liferegloop(p, b0, &param, maxiter, eps, colfit1, p-1);
//         double lmax0 = out0["loglik"];
//         prob[k] = boost_pchisq(-2*(lmax0 - lmax), 1, 0, 0);
//         clparm[k] = "PL";
//       }
//     } else {
//       for (int k=0; k<p; ++k) {
//         if (!robust) {
//           lb[k] = b[k] - zcrit*seb[k];
//           ub[k] = b[k] + zcrit*seb[k];
//           prob[k] = boost_pchisq(pow(b[k]/seb[k], 2), 1, 0, 0);
//         } else {
//           lb[k] = b[k] - zcrit*rseb[k];
//           ub[k] = b[k] + zcrit*rseb[k];
//           prob[k] = boost_pchisq(pow(b[k]/rseb[k], 2), 1, 0, 0);
//         }
//         clparm[k] = "Wald";
//       }
//     }
//     
//     for (int i=0; i<p; ++i) {
//       lb0[h*p+i] = lb[i];
//       ub0[h*p+i] = ub[i];
//       prob0[h*p+i] = prob[i];
//       clparm0[h*p+i] = clparm[i];
//     }
//     
//     // log-likelihoods
//     loglik0[h] = outint["loglik"];
//     loglik1[h] = out["loglik"];
//   }
//   
//   expbeta0 = std::exp(beta0);
//   if (!robust) z0 = beta0/sebeta0;
//   else z0 = beta0/rsebeta0;
//   
//   List sumstat = List::create(
//     _["n"] = nobs,
//     _["nevents"] = nevents,
//     _["loglik0"] = loglik0,
//     _["loglik1"] = loglik1,
//     _["niter"] = niter,
//     _["dist"] = dist1,
//     _["p"] = p,
//     _["nvar"] = nvar-1,
//     _["robust"] = robust,
//     _["fail"] = fails);
//   
//   List parest = List::create(
//     _["param"] = par0,
//     _["beta"] = beta0,
//     _["sebeta"] = robust ? rsebeta0 : sebeta0,
//     _["z"] = z0,
//     _["expbeta"] = expbeta0,
//     _["vbeta"] = robust ? rvbeta0 : vbeta0,
//     _["lower"] = lb0,
//     _["upper"] = ub0,
//     _["p"] = prob0,
//     _["method"] = clparm0);
//   
//   if (robust) {
//     parest.push_back(sebeta0, "sebeta_naive");
//     parest.push_back(vbeta0, "vbeta_naive");
//   }
//   
//   if (has_rep) {
//     for (int i=0; i<p_rep; ++i) {
//       std::string s = as<std::string>(rep[i]);
//       SEXP col = u_rep[s];
//       SEXPTYPE col_type = TYPEOF(col);
//       if (col_type == INTSXP) {
//         IntegerVector v = col;
//         sumstat.push_back(v[rep01], s);
//         parest.push_back(v[rep0], s);
//       } else if (col_type == REALSXP) {
//         NumericVector v = col;
//         sumstat.push_back(v[rep01], s);
//         parest.push_back(v[rep0], s);
//       } else if (col_type == STRSXP) {
//         StringVector v = col;
//         sumstat.push_back(v[rep01], s);
//         parest.push_back(v[rep0], s);
//       } else {
//         throw std::invalid_argument("Unsupported type for rep variable" + s);
//       }
//     }
//   }
//   
//   List result = List::create(
//     _["sumstat"] = as<DataFrame>(sumstat),
//     _["parest"] = as<DataFrame>(parest),
//     _["linear_predictors"] = linear_predictors);
//   
//   return result;
// }
// 
// // first and second derivatives of log likelihood with respect to eta and tau
// List f_ld_1(NumericVector eta, NumericVector sig, void *ex) {
//   aftparams *param = (aftparams *) ex;
//   int n = param->tstart.size();
//   
//   NumericVector g(n), dg(n), ddg(n), ds(n), dds(n), dsg(n);
//   
//   // Main loop to compute  derivatives
//   for (int i = 0; i < n; ++i) {
//     double sigma = sig[i];
//     double u, u1, u2, d1, d2, q1, q2, w1, w2;
//     
//     switch (param->status[i]) {
//     
//     case 1: // event
//       switch (param->dist_code) {
//       case 1: case 2: // exponential / weibull
//         u = (std::log(param->tstop[i]) - eta[i])/sigma;
//         g[i] = u - std::exp(u) - std::log(sigma);
//         dg[i] = -(1-std::exp(u))/sigma;
//         ddg[i] = -std::exp(u)/(sigma*sigma);
//         ds[i] = -1 + dg[i]*(sigma*u);
//         dds[i] = ddg[i]*pow(sigma*u,2) - dg[i]*(sigma*u);
//         dsg[i] = ddg[i]*(sigma*u) - dg[i];
//         break;
//       case 3: case 4: // lognormal / normal
//         u = (param->dist_code==3) ?
//         (std::log(param->tstop[i]) - eta[i])/sigma
//         : (param->tstop[i] - eta[i])/sigma;
//         g[i] = boost_dnorm(u,0,1,1) - std::log(sigma);
//         dg[i] = u/sigma;
//         ddg[i] = -1/(sigma*sigma);
//         ds[i] = -1 + dg[i]*(sigma*u);
//         dds[i] = ddg[i]*pow(sigma*u,2) - dg[i]*(sigma*u);
//         dsg[i] = ddg[i]*(sigma*u) - dg[i];
//         break;
//       case 5: case 6: // loglogistic / logistic
//         u = (param->dist_code==5)
//         ? (std::log(param->tstop[i]) - eta[i])/sigma
//         : (param->tstop[i] - eta[i])/sigma;
//         g[i] = boost_dlogis(u,0,1,1) - std::log(sigma);
//         dg[i] = (1 - 2*boost_plogis(u,0,1,0,0))/sigma;
//         ddg[i] = -2*boost_dlogis(u)/(sigma*sigma);
//         ds[i] = -1 + dg[i]*(sigma*u);
//         dds[i] = ddg[i]*pow(sigma*u,2) - dg[i]*(sigma*u);
//         dsg[i] = ddg[i]*(sigma*u) - dg[i];
//         break;
//       }
//       break;
//       
//     case 3: // interval censoring
//       switch (param->dist_code) {
//       case 1: case 2: // exponential / weibull
//         u1 = (std::log(param->tstart[i]) - eta[i])/sigma;
//         u2 = (std::log(param->tstop[i]) - eta[i])/sigma;
//         w1 = std::exp(u1); w2 = std::exp(u2);
//         q1 = std::exp(-w1); q2 = std::exp(-w2);
//         d1 = w1*q1; d2 = w2*q2;
//         g[i] = std::log(q1 - q2);
//         dg[i] = (d1 - d2)/(q1 - q2)/sigma;
//         ddg[i] = -(d1*(1-w1) - d2*(1-w2))/(q1 - q2)/(sigma*sigma) -
//           pow(dg[i],2);
//         ds[i] = (u1*d1 - u2*d2)/(q1 - q2);
//         dds[i] = (u2*u2*(1-w2)*d2 - u1*u1*(1-w1)*d1)/(q1 - q2) -
//           ds[i]*(1+ds[i]);
//         dsg[i] = (u2*(1-w2)*d2 - u1*(1-w1)*d1)/((q1 - q2)*sigma) -
//           dg[i]*(1+ds[i]);
//         break;
//       case 3: case 4: // lognormal / normal
//         u1 = (param->dist_code==3)
//         ? (std::log(param->tstart[i]) - eta[i])/sigma
//         : (param->tstart[i] - eta[i])/sigma;
//         u2 = (param->dist_code==3)
//           ? (std::log(param->tstop[i]) - eta[i])/sigma
//         : (param->tstop[i] - eta[i])/sigma;
//         d1 = boost_dnorm(u1); d2 = boost_dnorm(u2);
//         q1 = boost_pnorm(u1,0,1,0,0); q2 = boost_pnorm(u2,0,1,0,0);
//         g[i] = std::log(q1 - q2);
//         dg[i] = (d1 - d2)/(q1 - q2)/sigma;
//         ddg[i] = (d1*u1 - d2*u2)/(q1 - q2)/(sigma*sigma) - pow(dg[i],2);
//         ds[i] = (u1*d1 - u2*d2)/(q1 - q2);
//         dds[i] = (-u2*u2*u2*d2 + u1*u1*u1*d1)/(q1 - q2) - ds[i]*(1+ds[i]);
//         dsg[i] = (-u2*u2*d2 + u1*u1*d1)/((q1 - q2)*sigma) - dg[i]*(1+ds[i]);
//         break;
//       case 5: case 6: // loglogistic / logistic
//         u1 = (param->dist_code==5)
//         ? (std::log(param->tstart[i]) - eta[i])/sigma
//         : (param->tstart[i] - eta[i])/sigma;
//         u2 = (param->dist_code==5) ?
//         (std::log(param->tstop[i]) - eta[i])/sigma
//         : (param->tstop[i] - eta[i])/sigma;
//         d1 = boost_dlogis(u1); d2 = boost_dlogis(u2);
//         q1 = boost_plogis(u1,0,1,0,0); q2 = boost_plogis(u2,0,1,0,0);
//         g[i] = std::log(q1 - q2);
//         dg[i] = (d1 - d2)/(q1 - q2)/sigma;
//         ddg[i] = -(d1*(2*q1-1) - d2*(2*q2-1))/(q1 - q2)/(sigma*sigma) -
//           pow(dg[i],2);
//         ds[i] = (u1*d1 - u2*d2)/(q1 - q2);
//         dds[i] = (u2*u2*d2*(2*q2-1) - u1*u1*d1*(2*q1-1))/(q1 - q2) -
//           ds[i]*(1+ds[i]);
//         dsg[i] = (u2*d2*(2*q2-1) - u1*d1*(2*q1-1))/((q1 - q2)*sigma) -
//           dg[i]*(1+ds[i]);
//         break;
//       }
//       break;
//       
//     case 2: // left censoring
//       switch (param->dist_code) {
//       case 1: case 2: // exponential / weibull
//         u = (std::log(param->tstop[i]) - eta[i])/sigma;
//         g[i] = std::log(1 - std::exp(-std::exp(u)));
//         dg[i] = -std::exp(u - std::exp(u))/(1 - std::exp(-std::exp(u)))/sigma;
//         ddg[i] = (1 - std::exp(u) - std::exp(-std::exp(u)))*std::exp(u - std::exp(u))/
//           pow((1 - std::exp(-std::exp(u)))*sigma, 2);
//         ds[i] = dg[i]*(sigma*u);
//         dds[i] = ddg[i]*pow(sigma*u,2) - ds[i];
//         dsg[i] = ddg[i]*(sigma*u) - dg[i];
//         break;
//       case 3: case 4: // lognormal / normal
//         u = (param->dist_code==3)
//         ? (std::log(param->tstop[i]) - eta[i])/sigma
//         : (param->tstop[i] - eta[i])/sigma;
//         g[i] = boost_pnorm(u,0,1,1,1);
//         dg[i] = -boost_dnorm(u)/boost_pnorm(u)/sigma;
//         ddg[i] = -u*boost_dnorm(u)/(boost_pnorm(u)*sigma*sigma) -
//           pow(dg[i],2);
//         ds[i] = dg[i]*(sigma*u);
//         dds[i] = ddg[i]*pow(sigma*u,2) - ds[i];
//         dsg[i] = ddg[i]*(sigma*u) - dg[i];
//         break;
//       case 5: case 6: // loglogistic / logistic
//         u = (param->dist_code==5)
//         ? (std::log(param->tstop[i]) - eta[i])/sigma
//         : (param->tstop[i] - eta[i])/sigma;
//         g[i] = boost_plogis(u,0,1,1,1);
//         dg[i] = -boost_plogis(u,0,1,0,0)/sigma;
//         ddg[i] = -boost_dlogis(u)/(sigma*sigma);
//         ds[i] = dg[i]*(sigma*u);
//         dds[i] = ddg[i]*pow(sigma*u,2) - ds[i];
//         dsg[i] = ddg[i]*(sigma*u) - dg[i];
//         break;
//       }
//       break;
//       
//     case 0: // right censoring
//       switch (param->dist_code) {
//       case 1: case 2: // exponential / weibull
//         u = (std::log(param->tstart[i]) - eta[i])/sigma;
//         g[i] = -std::exp(u);
//         dg[i] = std::exp(u)/sigma;
//         ddg[i] = -std::exp(u)/(sigma*sigma);
//         ds[i] = dg[i]*(sigma*u);
//         dds[i] = ddg[i]*pow(sigma*u,2) - ds[i];
//         dsg[i] = ddg[i]*(sigma*u) - dg[i];
//         break;
//       case 3: case 4: // lognormal / normal
//         u = (param->dist_code==3)
//         ? (std::log(param->tstart[i]) - eta[i])/sigma
//         : (param->tstart[i] - eta[i])/sigma;
//         g[i] = boost_pnorm(u,0,1,0,1);
//         dg[i] = boost_dnorm(u)/boost_pnorm(u,0,1,0,0)/sigma;
//         ddg[i] = u*boost_dnorm(u)/(boost_pnorm(u,0,1,0,0)*sigma*sigma) -
//           pow(dg[i],2);
//         ds[i] = dg[i]*(sigma*u);
//         dds[i] = ddg[i]*pow(sigma*u,2) - ds[i];
//         dsg[i] = ddg[i]*(sigma*u) - dg[i];
//         break;
//       case 5: case 6: // loglogistic / logistic
//         u = (param->dist_code==5)
//         ? (std::log(param->tstart[i]) - eta[i])/sigma
//         : (param->tstart[i] - eta[i])/sigma;
//         g[i] = boost_plogis(u,0,1,0,1);
//         dg[i] = boost_plogis(u)/sigma;
//         ddg[i] = -boost_dlogis(u)/(sigma*sigma);
//         ds[i] = dg[i]*(sigma*u);
//         dds[i] = ddg[i]*pow(sigma*u,2) - ds[i];
//         dsg[i] = ddg[i]*(sigma*u) - dg[i];
//         break;
//       }
//       break;
//       
//     default:
//       throw std::runtime_error("Unknown status: " + std::to_string(param->status[i]));
//     }
//   }
//   
//   return List::create(
//     Named("g") = g,
//     Named("dg") = dg,
//     Named("ddg") = ddg,
//     Named("ds") = ds,
//     Named("dds") = dds,
//     Named("dsg") = dsg
//   );
// }
// 
// // [[Rcpp::export]]
// NumericMatrix residuals_liferegcpp(const NumericVector& beta,
//                                    const NumericMatrix& vbeta,
//                                    DataFrame data,
//                                    const StringVector& stratum = "",
//                                    const std::string time = "time",
//                                    const std::string time2 = "",
//                                    const std::string event = "event",
//                                    const StringVector& covariates = "",
//                                    const std::string weight = "",
//                                    const std::string offset = "",
//                                    const std::string id = "",
//                                    const std::string dist = "weibull",
//                                    const std::string type = "response",
//                                    const bool collapse = false,
//                                    const bool weighted = false) {
//   
//   int n = data.nrows();
//   int nvar = static_cast<int>(covariates.size()) + 1;
//   if (nvar == 2 && (covariates[0] == "" || covariates[0] == "none")) {
//     nvar = 1;
//   }
//   
//   std::string dist1 = dist;
//   std::for_each(dist1.begin(), dist1.end(), [](char & c) {
//     c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
//   });
//   
//   if (dist1 == "log-logistic" || dist1 == "llogistic") {
//     dist1 = "loglogistic";
//   } else if  (dist1 == "log-normal" || dist1 == "lnormal") {
//     dist1 = "lognormal";
//   } else if (dist1 == "gaussian") {
//     dist1 = "normal";
//   }
//   
//   int dist_code;
//   if (dist1 == "exponential") dist_code = 1;
//   else if (dist1 == "weibull") dist_code = 2;
//   else if (dist1 == "lognormal") dist_code = 3;
//   else if (dist1 == "normal") dist_code = 4;
//   else if (dist1 == "loglogistic") dist_code = 5;
//   else if (dist1 == "logistic") dist_code = 6;
//   else throw std::invalid_argument("invalid distribution: " + dist1);
//   
//   IntegerVector stratumn(n);
//   int p_stratum = static_cast<int>(stratum.size());
//   if (p_stratum == 1 && (stratum[0] == "" || stratum[0] == "none")) {
//     stratumn.fill(0);
//   } else {
//     List out = bygroup(data, stratum);
//     stratumn = out["index"];
//   }
//   
//   IntegerVector stratumn1 = unique(stratumn);
//   int nstrata = static_cast<int>(stratumn1.size());
//   int p = dist_code == 1 ? nvar : (nvar+nstrata);
//   
//   if (dist_code == 1 && nstrata > 1) {
//     throw std::invalid_argument("Stratification is not valid with the exponential distribution");
//   }
//   
//   bool has_time = hasVariable(data, time);
//   if (!has_time) throw std::invalid_argument("data must contain the time variable");
//   NumericVector timenz = data[time];
//   NumericVector timen = clone(timenz);
//   for (int i=0; i<n; ++i) {
//     if (!std::isnan(timen[i]) && (dist_code == 1 || dist_code == 2 ||
//         dist_code == 3 || dist_code == 5) && timen[i] <= 0) {
//       std::string str1 = "time must be positive for each subject for the";
//       std::string str2 = "distribution";
//       std::string errmsg = str1 + " " + dist1 + " " + str2;
//       throw std::invalid_argument(errmsg);
//     }
//   }
//   
//   bool has_time2 = hasVariable(data, time2);
//   NumericVector time2n(n);
//   if (has_time2) {
//     NumericVector time2nz = data[time2];
//     time2n = clone(time2nz);
//     for (int i=0; i<n; ++i) {
//       if (!std::isnan(time2n[i]) && (dist_code == 1 || dist_code == 2 ||
//           dist_code == 3 || dist_code == 5) && time2n[i] <= 0) {
//         std::string str1 = "time2 must be positive for each subject for the";
//         std::string str2 = "distribution";
//         std::string errmsg = str1 + " " + dist1 + " " + str2;
//         throw std::invalid_argument(errmsg);
//       }
//     }
//   }
//   
//   bool has_event = hasVariable(data, event);
//   if (!has_time2 && !has_event) {
//     throw std::invalid_argument("data must contain the event variable for right censored data");
//   }
//   
//   IntegerVector eventn(n);
//   if (has_event) {
//     IntegerVector eventnz = data[event];
//     eventn = clone(eventnz);
//     if (is_true(any((eventn != 1) & (eventn != 0)))) {
//       throw std::invalid_argument("event must be 1 or 0 for each subject");
//     }
//   }
//   
//   NumericMatrix zn(n,nvar);
//   for (int i=0; i<n; ++i) {
//     zn(i,0) = 1; // intercept
//   }
//   
//   for (int j=0; j<nvar-1; ++j) {
//     String zj = covariates[j];
//     if (!hasVariable(data, zj)) {
//       throw std::invalid_argument("data must contain the variables in covariates");
//     }
//     NumericVector u = data[zj];
//     for (int i=0; i<n; ++i) {
//       zn(i,j+1) = u[i];
//     }
//   }
//   
//   bool has_weight = hasVariable(data, weight);
//   NumericVector weightn(n, 1.0);
//   if (has_weight) {
//     NumericVector weightnz = data[weight];
//     weightn = clone(weightnz);
//     if (is_true(any(weightn <= 0))) {
//       throw std::invalid_argument("weight must be greater than 0");
//     }
//   }
//   
//   bool has_offset = hasVariable(data, offset);
//   NumericVector offsetn(n);
//   if (has_offset) {
//     NumericVector offsetnz = data[offset];
//     offsetn = clone(offsetnz);
//   }
//   
//   // create the numeric id variable
//   bool has_id = hasVariable(data, id);
//   IntegerVector idn(n);
//   if (!has_id) {
//     idn = seq(0, n - 1);
//   } else {
//     SEXP col = data[id];
//     SEXPTYPE col_type = TYPEOF(col);
//     if (col_type == INTSXP) {
//       IntegerVector v = col;
//       IntegerVector w = unique(v);
//       w.sort();
//       idn = match(v, w) - 1;
//     } else if (col_type == REALSXP) {
//       NumericVector v = col;
//       NumericVector w = unique(v);
//       w.sort();
//       idn = match(v, w) - 1;
//     } else if (col_type == STRSXP) {
//       StringVector v = col;
//       StringVector w = unique(v);
//       w.sort();
//       idn = match(v, w) - 1;
//     } else {
//       throw std::invalid_argument("incorrect type for the id variable in the input data");
//     }
//   }
//   
//   // exclude observations with missing values
//   LogicalVector sub(n,1);
//   for (int i=0; i<n; ++i) {
//     if (stratumn[i] == NA_INTEGER || idn[i] == NA_INTEGER ||
//         (std::isnan(timen[i]) && std::isnan(time2n[i])) ||
//         std::isnan(weightn[i]) || std::isnan(offsetn[i])) {
//       sub[i] = 0;
//     }
//     for (int j=0; j<nvar-1; ++j) {
//       if (std::isnan(zn(i,j+1))) sub[i] = 0;
//     }
//   }
//   
//   IntegerVector q1 = which(sub);
//   stratumn = stratumn[q1];
//   timen = timen[q1];
//   time2n = time2n[q1];
//   eventn = eventn[q1];
//   weightn = weightn[q1];
//   offsetn = offsetn[q1];
//   idn = idn[q1];
//   zn = subset_matrix_by_row(zn, q1);
//   int n1 = sum(sub);
//   if (n1 == 0) throw std::invalid_argument("no observation is left after removing missing values");
//   
//   // unify right censored data with interval censored data
//   NumericVector tstart(n1), tstop(n1);
//   if (!has_time2) {
//     tstart = timen;
//     for (int i=0; i<n1; ++i) {
//       tstop[i] = eventn[i] == 1 ? tstart[i] : NA_REAL;
//     }
//   } else {
//     tstart = timen;
//     tstop = time2n;
//   }
//   
//   IntegerVector status(n1);
//   for (int i=0; i<n1; ++i) {
//     if (!std::isnan(tstart[i]) && !std::isnan(tstop[i]) &&
//         tstart[i] == tstop[i]) {
//       status[i] = 1; // event
//     } else if (!std::isnan(tstart[i]) && !std::isnan(tstop[i]) &&
//       tstart[i] < tstop[i]) {
//       status[i] = 3; // interval censoring
//     } else if (std::isnan(tstart[i]) && !std::isnan(tstop[i])) {
//       status[i] = 2; // left censoring
//     } else if (!std::isnan(tstart[i]) && std::isnan(tstop[i])) {
//       status[i] = 0; // right censoring
//     } else {
//       status[i] = -1; // exclude the observation
//     }
//   }
//   
//   // exclude records with invalid status
//   IntegerVector q2 = which(status != -1);
//   int n2 = static_cast<int>(q2.size());
//   
//   if (n2 < n1) {
//     stratumn = stratumn[q2];
//     tstart = tstart[q2];
//     tstop = tstop[q2];
//     status = status[q2];
//     weightn = weightn[q2];
//     offsetn = offsetn[q2];
//     idn = idn[q2];
//     zn = subset_matrix_by_row(zn,q2);
//   }
//   
//   NumericVector eta(n2);
//   for (int i = 0; i < n2; ++i) {
//     eta[i] = offsetn[i];
//     for (int j=0; j<nvar; ++j) {
//       eta[i] += beta[j]*zn(i,j);
//     }
//   }
//   
//   NumericVector sig(n2, 1.0);
//   if (dist_code != 1) {
//     for (int i = 0; i < n2; ++i) {
//       int j = stratumn[i] + nvar;
//       sig[i] = std::exp(beta[j]);
//     }
//   }
//   
//   int K = 1;
//   if (type == "dfbeta" || type == "dfbetas") {
//     K = p;
//   } else if (type == "matrix") {
//     K = 6;
//   }
//   
//   // Map type to integer code
//   int type_code;
//   if (type == "response") type_code = 1;
//   else if (type == "martingale") type_code = 2;
//   else if (type == "deviance") type_code = 3;
//   else if (type == "working") type_code = 4;
//   else if (type == "dfbeta") type_code = 5;
//   else if (type == "dfbetas") type_code = 6;
//   else if (type == "ldcase") type_code = 7;
//   else if (type == "ldresp") type_code = 8;
//   else if (type == "ldshape") type_code = 9;
//   else if (type == "matrix") type_code = 10;
//   else throw std::invalid_argument("invalid type of residuals" + type);
//   
//   // rr: residual matrix
//   NumericMatrix rr(n2, K);
//   
//   switch (type_code) {
//   
//   case 1: { // response
//     NumericVector yhat0(n2);
//     for (int i = 0; i < n2; ++i) {
//       switch (status[i]) {
//       case 0: case 1: // right-censored or event
//         switch (dist_code) {
//         case 1: case 2: case 3: case 5: yhat0[i] = std::log(tstart[i]); break;
//         default: yhat0[i] = tstart[i];
//         }
//         break;
//       case 2: // left-censored
//         switch (dist_code) {
//         case 1: case 2: case 3: case 5: yhat0[i] = std::log(tstop[i]); break;
//         default: yhat0[i] = tstop[i];
//         }
//         break;
//       default: // interval-censored
//         switch (dist_code) {
//         case 1: case 2: {
//           double width = (std::log(tstop[i]) - std::log(tstart[i])) / sig[i];
//           yhat0[i] = std::log(tstart[i]) - sig[i] * std::log(width / (std::exp(width) - 1));
//           break;
//         }
//         case 3: case 5:
//           yhat0[i] = 0.5 * (std::log(tstart[i]) + std::log(tstop[i])); break;
//         default: yhat0[i] = 0.5 * (tstart[i] + tstop[i]);
//         }
//       }
//       
//       switch (dist_code) {
//       case 1: case 2: case 3: case 5:
//         rr(i,0) = std::exp(yhat0[i]) - std::exp(eta[i]); break;
//       default: rr(i,0) = yhat0[i] - eta[i];
//       }
//     }
//     break;
//   }
//     
//   case 2: { // martingale
//     if (dist_code == 4 || dist_code == 6)
//       throw std::invalid_argument("incorrect type of distribution" + dist1 +
//         " for martingale residuals");
//     for (int i = 0; i < n2; ++i) {
//       if (status[i] == 0 || status[i] == 1) {
//         double y = (std::log(tstart[i]) - eta[i]) / sig[i];
//         switch (dist_code) {
//         case 1: case 2: rr(i,0) = (status[i] == 1) - std::exp(y); break;
//         case 3: rr(i,0) = (status[i] == 1) + boost_pnorm(y,0,1,0,1); break;
//         case 5: rr(i,0) = (status[i] == 1) + boost_plogis(y,0,1,0,1); break;
//         }
//       } else rr(i,0) = NA_REAL;
//     }
//     break;
//   }
//     
//   default: { // other types: deviance, working, dfbeta, ld*, matrix
//     aftparams param = {dist_code, stratumn, tstart, tstop, status, weightn,
//                        offsetn, zn, nstrata};
//     List der = f_ld_1(eta, sig, &param);
//     NumericVector g = der["g"], dg = der["dg"], ddg = der["ddg"];
//     NumericVector ds = der["ds"], dds = der["dds"], dsg = der["dsg"];
//     
//     switch (type_code) {
//     
//     case 3: { // deviance
//       NumericVector loglik(n2);
//       for (int i = 0; i < n2; ++i) {
//         switch (status[i]) {
//         case 0: case 2: loglik[i] = 0; break; // right or left censored
//         case 1:  // event
//           switch (dist_code) {
//           case 1: case 2: loglik[i] = -std::log(sig[i]) - 1; break;
//           case 3: case 4: loglik[i] = -std::log(std::sqrt(2*M_PI)*sig[i]); break;
//           default: loglik[i] = -std::log(4*sig[i]);
//           }
//           break;
//         default: { // interval censored
//             double width;
//             switch (dist_code) {
//             case 1: case 2: width = (std::log(tstop[i]) - std::log(tstart[i])) / sig[i];
//               loglik[i] = - width/(std::exp(width)-1) + std::log(1 - std::exp(-width)); break;
//             case 3: width = (std::log(tstop[i]) - std::log(tstart[i])) / sig[i];
//               loglik[i] = std::log(2*boost_pnorm(width/2) - 1); break;
//             case 4: width = (tstop[i] - tstart[i]) / sig[i];
//               loglik[i] = std::log(2*boost_pnorm(width/2) - 1); break;
//             case 5: width = (std::log(tstop[i]) - std::log(tstart[i])) / sig[i];
//               loglik[i] = std::log((std::exp(width/2)-1)/(std::exp(width/2)+1)); break;
//             default: width = (tstop[i] - tstart[i]) / sig[i];
//             loglik[i] = std::log((std::exp(width/2)-1)/(std::exp(width/2)+1));
//             }
//           }
//         }
//         double val = -dg[i]/ddg[i];
//         rr(i,0) = ((val>0) - (val<0))*std::sqrt(2*(loglik[i] - g[i]));
//       }
//       break;
//     }
//       
//     case 4: { // working
//       for (int i=0; i<n2; ++i) rr(i,0) = -dg[i]/ddg[i];
//       break;
//     }
//       
//     case 5: case 6: case 7: { // dfbeta, dfbetas, ldcase
//       for (int i=0; i<n2; ++i) {
//       NumericVector score(p), resid(p);
//       for (int j=0; j<nvar; ++j) score[j] = dg[i]*zn(i,j);
//       for (int j=nvar; j<p; ++j) score[j] = stratumn[i]==j-nvar ? ds[i] : 0;
//       for (int k=0; k<p; ++k) for (int j=0; j<p; ++j)
//         resid[k] += score[j]*vbeta(j,k);
//       if (type_code==6) for (int k=0; k<p; ++k) resid[k] /= std::sqrt(vbeta(k,k));
//       if (type_code==7) for (int k=0; k<p; ++k) rr(i,0) += resid[k]*score[k];
//       else for (int k=0; k<p; ++k) rr(i,k) = resid[k];
//     }
//       break;
//     }
//       
//     case 8: { // ldresp
//       for (int i=0; i<n2; ++i) {
//       NumericVector rscore(p), temp(p);
//       for (int j=0; j<nvar; ++j) rscore[j] = -ddg[i]*zn(i,j)*sig[i];
//       for (int j=nvar; j<p; ++j)
//         rscore[j] = stratumn[i]==j-nvar ? -dsg[i]*sig[i] : 0;
//       for (int k=0; k<p; ++k) for (int j=0; j<p; ++j)
//         temp[k] += rscore[j]*vbeta(j,k);
//       for (int k=0; k<p; ++k) rr(i,0) += temp[k]*rscore[k];
//     }
//       break;
//     }
//       
//     case 9: { // ldshape
//       for (int i=0; i<n2; ++i) {
//       NumericVector sscore(p), temp(p);
//       for (int j=0; j<nvar; ++j) sscore[j] = dsg[i]*zn(i,j);
//       for (int j=nvar; j<p; ++j) sscore[j] = stratumn[i]==j-nvar ? dds[i] : 0;
//       for (int k=0; k<p; ++k) for (int j=0; j<p; ++j)
//         temp[k] += sscore[j]*vbeta(j,k);
//       for (int k=0; k<p; ++k) rr(i,0) += temp[k]*sscore[k];
//     }
//       break;
//     }
//       
//     case 10: { // matrix
//       for (int i=0; i<n2; ++i) {
//       rr(i,0)=g[i]; rr(i,1)=dg[i]; rr(i,2)=ddg[i];
//       rr(i,3)=ds[i]; rr(i,4)=dds[i]; rr(i,5)=dsg[i];
//     }
//       break;
//     }
//       
//     }
//   } // end default
//   } // end switch (type_code)
//   
//   // case weights
//   if (weighted) {
//     for (int i = 0; i < n2; ++i) {
//       for (int k=0; k<K; ++k) {
//         rr(i,k) *= weightn[i];
//       }
//     }
//   }
//   
//   // collapse if needed
//   if (collapse) {
//     // order data by id
//     IntegerVector order = seq(0, n2-1);
//     std::sort(order.begin(), order.end(), [&](int i, int j) {
//       return idn[i] < idn[j];
//     });
//     
//     IntegerVector id2 = idn[order];
//     IntegerVector idx(1,0);
//     for (int i=1; i<n2; ++i) {
//       if (id2[i] != id2[i-1]) {
//         idx.push_back(i);
//       }
//     }
//     
//     int nids = static_cast<int>(idx.size());
//     idx.push_back(n2);
//     
//     // collapse over id
//     NumericMatrix rr2(nids,K);
//     for (int i=0; i<nids; ++i) {
//       for (int k=0; k<K; ++k) {
//         for (int j=idx[i]; j<idx[i+1]; ++j) {
//           rr2(i,k) += rr(order[j],k);
//         }
//       }
//     }
//     
//     rr = rr2;
//   }
//   
//   return rr;
// }
// 
