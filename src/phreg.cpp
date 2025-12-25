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

struct coxparams {
  int nused;
  std::vector<int> strata;
  std::vector<double> tstart;
  std::vector<double> tstop;
  std::vector<int> event;
  std::vector<double> weight;
  std::vector<double> offset;
  FlatMatrix z;
  std::vector<int> order1;
  int method; // 1: breslow, 2: efron
};

// all-in-one function for log-likelihood, score, and information matrix
// for the Cox model with or without Firth's correction
ListCpp f_der_2(int p, const std::vector<double>& par, void* ex, bool firth) {
  coxparams* param = (coxparams*) ex;
  
  const int nused = param->nused;
  const int method = param->method;
  
  // Precompute eta and exp(eta)
  std::vector<double> eta(nused);
  std::vector<double> exp_eta(nused);
  for (int person = 0; person < nused; ++person) {
    double val = param->offset[person];
    for (int i = 0; i < p; ++i) {
      val += par[i] * param->z(person,i);
    }
    eta[person] = val;
    exp_eta[person] = std::exp(val);
  }
  
  double loglik = 0.0;        // log-likelihood
  std::vector<double> u(p);         // score vector
  FlatMatrix imat(p,p);    // information matrix
  FlatMatrix dimat(p,p*p); // tensor for third order derivatives
  std::vector<double> a(p);         // s1(beta,k,t)
  std::vector<double> a2(p);        // sum of w*exp(zbeta)*z for the deaths
  FlatMatrix cmat(p,p);    // s2(beta,k,t)
  FlatMatrix cmat2(p,p);   // sum of w*exp(zbeta)*z*z' for the deaths
  FlatMatrix dmat(p,p*p);  // q2(beta,k,t)
  FlatMatrix dmat2(p,p*p); // sum of w*exp(zbeta)*z*z*z' for the deaths
  double denom = 0.0;         // s0(beta,k,t)
  double denom2 = 0.0;        // sum of weighted risks for deaths
  double deadwt = 0.0;        // sum of weights for the deaths
  int ndead = 0;              // number of deaths at this time point
  
  
  int istrata = param->strata[0];
  int i1 = 0; // index for removing out-of-risk subjects
  
  // Loop through subjects
  for (int person = 0; person < nused; ) {
    // Reset when entering a new stratum
    if (param->strata[person] != istrata) {
      istrata = param->strata[person];
      i1 = person;
      denom = 0.0;
      
      for (int i = 0; i < p; ++i) {
        a[i] = 0.0;
        for (int j = 0; j <= i; ++j) {
          cmat(i,j) = 0.0;
          if (firth) {
            for (int k = 0; k <= j; ++k) {
              dmat(i,k*p+j) = 0.0;
            }
          }
        }
      }
    }
    
    const double dtime = param->tstop[person];
    
    // Process all persons tied at this dtime
    for (; person < nused && param->tstop[person] == dtime &&
         param->strata[person] == istrata; ++person) {
      
      const double w = param->weight[person];
      const double r = w * exp_eta[person];
      
      if (param->event[person] == 0) {
        denom += r;
        
        for (int i = 0; i < p; ++i) {
          const double zi = param->z(person,i);
          a[i] += r * zi;
          for (int j = 0; j <= i; ++j) {
            const double zj = param->z(person,j);
            cmat(i,j) += r * zi * zj;
            if (firth) {
              for (int k = 0; k <= j; ++k) {
                dmat(i,k*p+j) += r * zi * zj * param->z(person,k);
              }
            }
          }
        }
      } else {
        ++ndead;
        deadwt += w;
        denom2 += r;
        loglik += w * eta[person];
        for (int i = 0; i < p; ++i) {
          const double zi = param->z(person,i);
          a2[i] += r * zi;
          u[i] += w * zi;
          for (int j = 0; j <= i; ++j) {
            const double zj = param->z(person,j);
            cmat2(i,j) += r * zi * zj;
            if (firth) {
              for (int k = 0; k <= j; ++k) {
                dmat2(i,k*p+j) += r * zi * zj * param->z(person,k);
              }
            }
          }
        }
      }
    }
    
    // Remove subjects leaving risk set
    for (; i1 < nused; ++i1) {
      const int p1 = param->order1[i1];
      if (param->tstart[p1] < dtime || param->strata[p1] != istrata) break;
      
      const double r = param->weight[p1] * exp_eta[p1];
      denom -= r;
      for (int i = 0; i < p; ++i) {
        const double zi = param->z(p1,i);
        a[i] -= r * zi;
        for (int j = 0; j <= i; ++j) {
          const double zj = param->z(p1,j);
          cmat(i,j) -= r * zi * zj;
          if (firth) {
            for (int k = 0; k <= j; ++k) {
              dmat(i,k*p+j) -= r * zi * zj * param->z(p1,k);
            }
          }
        }
      }
    }
    
    // Add contributions for deaths at this time
    if (ndead > 0) {
      if (method == 0 || ndead == 1) { // Breslow or single event
        denom += denom2;
        loglik -= deadwt * std::log(denom);
        for (int i = 0; i < p; ++i) {
          a[i] += a2[i];
          const double xbar = a[i] / denom;
          u[i] -= deadwt * xbar;
          for (int j = 0; j <= i; ++j) {
            cmat(i,j) += cmat2(i,j);
            imat(i,j) += deadwt * (cmat(i,j) - xbar * a[j]) / denom;
            if (firth) {
              for (int k = 0; k <= j; ++k) {
                dmat(i,k*p+j) += dmat2(i,k*p+j);
                dimat(i,k*p+j) += deadwt * (dmat(i,k*p+j) -
                  (cmat(i,j)*a[k] + cmat(i,k)*a[j] + cmat(j,k)*a[i]) / denom +
                  2.0 * a[i] * a[j] * a[k] / (denom * denom)) / denom;
              }
            }
          }
        }
      } else { // Efron method
        const double meanwt = deadwt / ndead;
        const double increment = denom2 / ndead;
        for (int l = 0; l < ndead; ++l) {
          denom += increment;
          loglik -= meanwt * std::log(denom);
          for (int i = 0; i < p; ++i) {
            a[i] += a2[i] / ndead;
            const double xbar = a[i] / denom;
            u[i] -= meanwt * xbar;
            for (int j = 0; j <= i; ++j) {
              cmat(i,j) += cmat2(i,j) / ndead;
              imat(i,j) += meanwt * (cmat(i,j) - xbar * a[j]) / denom;
              if (firth) {
                for (int k = 0; k <= j; ++k) {
                  dmat(i,k*p+j) += dmat2(i,k*p+j) / ndead;
                  dimat(i,k*p+j) += meanwt * (dmat(i,k*p+j) -
                    (cmat(i,j)*a[k] + cmat(i,k)*a[j] + cmat(j,k)*a[i])/denom +
                    2.0 * a[i] * a[j] * a[k] / (denom * denom)) / denom;
                }
              }
            }
          }
        }
      }
      
      // Reset after processing deaths
      ndead = 0;
      deadwt = denom2 = 0.0;
      for (int i = 0; i < p; ++i) {
        a2[i] = 0;
        for (int j = 0; j <= i; ++j) {
          cmat2(i,j) = 0;
          if (firth) {
            for (int k = 0; k <= j; ++k) {
              dmat2(i,k*p+j) = 0;
            }
          }
        }
      }
    }
  }
  
  
  // fill the symmetric elements of the information matrix
  for (int i = 0; i < p - 1; ++i)
    for (int j = i+1; j < p; ++j)
      imat(i,j) = imat(j,i);
  
  
  // fill the symmetric elements of the tensor array
  if (firth) {
    for (int i = 0; i < p-1; ++i)
      for (int j = i+1; j < p; ++j)
        for (int k = 0; k <= i; ++k)
          dimat(i,k*p+j) = dimat(j,k*p+i);
    
    for (int j = 0; j < p-1; ++j)
      for (int i = j; i < p; ++i)
        for (int k = j+1; k < p; ++k)
          dimat(i,k*p+j) = dimat(i,j*p+k);
    
    for (int i = 0; i < p-1; ++i)
      for (int j = i+1; j < p; ++j)
        for (int k = i+1; k < p; ++k)
          dimat(i,k*p+j) = dimat(k,i*p+j);
  }
  
  
  ListCpp result;
  
  // Firth adjustment
  if (p > 0 && firth) {
    // obtain the determinant of information matrix
    FlatMatrix imat0 = imat;
    double toler = 1e-12;
    cholesky2(imat0, p, toler);
    
    double v = 0;
    for (int i=0; i<p; ++i) {
      v += std::log(imat0(i,i));
    }
    
    // penalized log-likelihood adjustment
    double penloglik = loglik + 0.5*v;
    
    // compute the bias adjustment to the score function
    FlatMatrix y(p,p);
    std::vector<double> g(p);
    
    for (int k = 0; k < p; ++k) {
      // partial derivative of the information matrix w.r.t. beta[k]
      for (int i = 0; i < p; ++i) {
        for (int j = 0; j < p; ++j) {
          y(i,j) = dimat(i,k*p+j);
        }
      }
      
      // solve(imat, y)
      for (int h = 0; h < p; ++h) {
        for (int i = 0; i < p; ++i) {
          double temp = y(i,h);
          for (int j = 0; j < i; ++j)
            temp -= y(j,h) * imat0(j,i);
          y(i,h) = temp;
        }
        
        for (int i = p-1; i >= 0; --i) {
          if (imat0(i,i) == 0) y(i,h) = 0;
          else {
            double temp = y(i,h) / imat0(i,i);
            for (int j = i+1; j < p; ++j)
              temp -= y(j,h)*imat0(i,j);
            y(i,h) = temp;
          }
        }
      }
      
      // trace
      for (int i = 0; i < p; ++i) g[k] += y(i,i);
      
      g[k] = u[k] + 0.5 * g[k];
    }
    
    result.push_back(penloglik, "loglik");
    result.push_back(g, "score");
    result.push_back(imat, "imat");
    result.push_back(loglik, "regloglik");
    result.push_back(u, "regscore");
  } else {
    result.push_back(loglik, "loglik");
    if (p > 0) {
      result.push_back(u, "score");
      result.push_back(imat, "imat");
    }
  }
  
  return result;
}


// // underlying optimization algorithm for phreg
// List phregloop(int p, const NumericVector& par, void *ex,
//                int maxiter, double eps, bool firth,
//                const IntegerVector& colfit, int ncolfit) {
//   coxparams *param = (coxparams *) ex;
//   
//   int iter, halving = 0;
//   bool fail = false;
//   double toler = 1e-12;
//   
//   NumericVector beta(p), newbeta(p);
//   double loglik, newlk = 0;
//   NumericVector u(p);
//   NumericMatrix imat(p,p);
//   NumericVector u1(ncolfit);
//   NumericMatrix imat1(ncolfit, ncolfit);
//   NumericMatrix z1 = param->z;
//   
//   // --- first step ---
//   beta = clone(par);
//   
//   List der = f_der_2(p, beta, param, firth);
//   loglik = der["loglik"];
//   u = der["score"];
//   imat = as<NumericMatrix>(der["imat"]);
//   
//   for (int i = 0; i < ncolfit; ++i) u1[i] = u[colfit[i]];
//   
//   for (int i = 0; i < ncolfit; ++i)
//     for (int j = 0; j < ncolfit; ++j)
//       imat1(i,j) = imat(colfit[i], colfit[j]);
//   
//   cholesky2(imat1, ncolfit, toler);
//   chsolve2(imat1, ncolfit, u1);
//   
//   u.fill(0.0);
//   for (int i = 0; i < ncolfit; ++i) u[colfit[i]] = u1[i];
//   newbeta = beta + u;
//   
//   // --- main iteration ---
//   for (iter = 0; iter < maxiter; ++iter) {
//     der = f_der_2(p, newbeta, param, firth);
//     newlk = der["loglik"];
//     
//     fail = std::isnan(newlk) || std::isinf(newlk);
//     if (!fail && halving == 0 && fabs(1 - (loglik/newlk)) < eps) break;
//     
//     if (fail || newlk < loglik) {
//       ++halving; // adjust step size if likelihood decreases
//       for (int i = 0; i < p; ++i) {
//         newbeta[i] = 0.5 * (beta[i] + newbeta[i]);
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
//         imat1(i,j) = imat(colfit[i], colfit[j]);
//     
//     cholesky2(imat1, ncolfit, toler);
//     chsolve2(imat1, ncolfit, u1);
//     
//     u.fill(0.0);
//     for (int i = 0; i < ncolfit; ++i) u[colfit[i]] = u1[i];
//     newbeta = beta + u;
//   }
//   
//   if (iter == maxiter) fail = true;
//   
//   // --- final variance calculation ---
//   imat = as<NumericMatrix>(der["imat"]);
//   for (int i = 0; i < ncolfit; ++i)
//     for (int j = 0; j < ncolfit; ++j)
//       imat1(i,j) = imat(colfit[i], colfit[j]);
//   
//   NumericMatrix var1 = invsympd(imat1, ncolfit, toler);
//   NumericMatrix var(p,p);
//   for (int i = 0; i < ncolfit; ++i)
//     for (int j = 0; j < ncolfit; ++j)
//       var(colfit[i], colfit[j]) = var1(i,j);
//   
//   List result = List::create(
//     Named("coef") = newbeta,
//     Named("iter") = iter,
//     Named("var") = var,
//     Named("loglik") = newlk,
//     Named("fail") = fail);
//   
//   if (firth) {
//     double regloglik = as<double>(der["regloglik"]);
//     result.push_back(regloglik, "regloglik");
//   }
//   
//   return result;
// }
// 
// 
// // confidence limit of profile likelihood method
// double phregplloop(int p, const NumericVector& par, void *ex,
//                    int maxiter, double eps, bool firth,
//                    int k, int which, double l0) {
//   coxparams *param = (coxparams *) ex;
//   
//   int iter;
//   bool fail = false;
//   double toler = 1e-12;
//   
//   NumericVector beta(p), newbeta(p);
//   double loglik, newlk;
//   NumericVector u(p);
//   NumericVector delta(p);
//   NumericMatrix imat(p,p);
//   NumericMatrix v(p,p);
//   
//   // --- first step ---
//   beta = clone(par);
//   
//   List der = f_der_2(p, beta, param, firth);
//   loglik = der["loglik"];
//   u = der["score"];
//   imat = as<NumericMatrix>(der["imat"]);
//   
//   // Lagrange multiplier method as used in SAS PROC LOGISTIC
//   v = invsympd(imat, p, toler);
//   
//   double w = 0.0;
//   for (int i = 0; i < p; ++i)
//     for (int j = 0; j < p; ++j)
//       w -= u[i] * v(i,j) * u[j];
//   
//   double underroot = -2*(l0 - loglik + 0.5*w)/v(k,k);
//   double lambda = underroot < 0.0 ? 0.0 : which*std::sqrt(underroot);
//   u[k] += lambda;
//   
//   delta.fill(0.0);
//   for (int i = 0; i < p; ++i)
//     for (int j = 0; j < p; ++j)
//       delta[i] += v(i,j) * u[j];
//   
//   // update beta
//   newbeta = beta + delta;
//   
//   // --- main iteration ---
//   for (iter = 0; iter < maxiter; ++iter) {
//     der = f_der_2(p, newbeta, param, firth);
//     newlk = der["loglik"];
//     
//     fail = std::isnan(newlk) || std::isinf(newlk);
//     if (!fail && fabs(newlk - l0) < eps && w < eps) break;
//     
//     beta = clone(newbeta);
//     loglik = newlk;
//     u = as<NumericVector>(der["score"]);
//     imat = as<NumericMatrix>(der["imat"]);
//     
//     // Lagrange multiplier method as used in SAS PROC LOGISTIC
//     v = invsympd(imat, p, toler);
//     
//     double w = 0.0;
//     for (int i = 0; i < p; ++i)
//       for (int j = 0; j < p; ++j)
//         w -= u[i] * v(i,j) * u[j];
//     
//     double underroot = -2*(l0 - loglik + 0.5*w)/v(k,k);
//     double lambda = underroot < 0.0 ? 0.0 : which*std::sqrt(underroot);
//     u[k] += lambda;
//     
//     delta.fill(0.0);
//     for (int i = 0; i < p; ++i)
//       for (int j = 0; j < p; ++j)
//         delta[i] += v(i,j) * u[j];
//     
//     // update beta
//     newbeta = beta + delta;
//   }
//   
//   if (iter == maxiter) fail = true;
//   if (fail) warning("The algorithm in phregplloop did not converge");
//   
//   return newbeta[k];
// }
// 
// 
// // baseline hazard estimates
// List f_basehaz(int p, const NumericVector& par, void *ex) {
//   coxparams *param = (coxparams *) ex;
//   
//   const int nused = param->nused;
//   const int method = param->method;
//   
//   // Precompute exp(eta)
//   NumericVector exp_eta(nused);
//   for (int person = 0; person < nused; ++person) {
//     double val = param->offset[person];
//     for (int i = 0; i < p; ++i) {
//       val += par[i] * param->z(person,i);
//     }
//     exp_eta[person] = std::exp(val);
//   }
//   
//   NumericVector a(p);         // s1(beta,k,t)
//   NumericVector a2(p);        // sum of w*exp(zbeta)*z for the deaths
//   double deadwt = 0.0;        // sum of weights for the deaths
//   double denom = 0.0;         // s0(beta,k,t)
//   double denom2 = 0.0;        // sum of weighted risks for the deaths
//   int natrisk = 0;            // number at risk at this time point
//   int ndead = 0;              // number of deaths at this time point
//   int ncens = 0;              // number of censored at this time point
//   
//   // locate the first observation within each stratum
//   IntegerVector istratum(1,0);
//   for (int i = 1; i < nused; ++i) {
//     if (param->strata[i] != param->strata[i-1]) {
//       istratum.push_back(i);
//     }
//   }
//   
//   int nstrata = static_cast<int>(istratum.size());
//   istratum.push_back(nused);
//   
//   // add time 0 to each stratum
//   int J = nstrata;
//   for (int i = 0; i < nstrata; ++i) {
//     IntegerVector idx = seq(istratum[i], istratum[i+1] - 1);
//     NumericVector utime = param->tstop[idx];
//     utime = unique(utime);
//     J += static_cast<int>(utime.size());
//   }
//   
//   IntegerVector stratum(J);
//   NumericVector time(J), nrisk(J), nevent(J), ncensor(J), haz(J), varhaz(J);
//   NumericMatrix gradhaz(J,p);
//   
//   int istrata = param->strata[0];
//   int i1 = 0; // index for removing out-of-risk subjects
//   int j = J;  // index the unique time in ascending order
//   
//   // Loop through subjects
//   for (int person = 0; person < nused; ) {
//     if (param->strata[person] != istrata) { // hit a new stratum
//       // add time 0 at the start of a new stratum
//       j--;
//       stratum[j] = istrata;
//       time[j] = 0.0;
//       nrisk[j] = natrisk;
//       
//       istrata = param->strata[person]; // reset temporary variables
//       i1 = person;
//       natrisk = 0;
//       denom = 0.0;
//       a.fill(0.0);
//     }
//     
//     const double dtime = param->tstop[person];
//     
//     // Process all persons tied at this dtime
//     bool first = true;
//     for (; person < nused && param->tstop[person] == dtime &&
//          param->strata[person] == istrata; ++person) {
//       
//       if (first) { // first incidence at this time
//         j--;
//         stratum[j] = param->strata[person];
//         time[j] = dtime;
//         first = false;
//       }
//       
//       const double w = param->weight[person];
//       const double r = w * exp_eta[person];
//       
//       ++natrisk;
//       if (param->event[person] == 0) {
//         ++ncens;
//         denom += r;
//         for (int i = 0; i < p; ++i) {
//           a[i] += r * param->z(person,i);
//         }
//       } else {
//         ++ndead;
//         deadwt += w;
//         denom2 += r;
//         for (int i = 0; i < p; ++i) {
//           a2[i] += r * param->z(person,i);
//         }
//       }
//     }
//     
//     // Remove subjects leaving risk set
//     for (; i1 < nused; ++i1) {
//       const int p1 = param->order1[i1];
//       if (param->tstart[p1] < dtime || param->strata[p1] != istrata) break;
//       
//       const double r = param->weight[p1] * exp_eta[p1];
//       
//       natrisk--;
//       denom -= r;
//       for (int i = 0; i < p; ++i) {
//         a[i] -= r * param->z(p1,i);
//       }
//     }
//     
//     // Add contributions for deaths at this time
//     nrisk[j] = natrisk;
//     nevent[j] = ndead;
//     ncensor[j] = ncens;
//     ncens = 0; // reset for the next time point
//     if (ndead > 0) {
//       if (method == 0 || ndead == 1) { // Breslow or single event
//         denom += denom2;
//         const double temp = deadwt / denom;
//         haz[j] = temp;
//         varhaz[j] = temp / denom;
//         for (int i = 0; i < p; ++i) {
//           a[i] += a2[i];
//           gradhaz(j,i) = temp * a[i] / denom;
//         }
//       } else { // Efron method
//         const double meanwt = deadwt / ndead;
//         for (int k = 0; k < ndead; ++k) {
//           denom += denom2 / ndead;
//           const double temp = meanwt / denom;
//           haz[j] += temp;
//           varhaz[j] += temp / denom;
//           for (int i = 0; i < p; ++i) {
//             a[i] += a2[i] / ndead;
//             gradhaz(j,i) += temp * a[i] / denom;
//           }
//         }
//       }
//       
//       // reset for the next death time
//       ndead = 0;
//       deadwt = denom2 = 0.0;
//       a2.fill(0.0);
//     }
//   }
//   
//   // add time 0 for the first stratum
//   stratum[0] = istrata;
//   time[0] = 0.0;
//   nrisk[0] = natrisk;
//   
//   List result = List::create(
//     _["stratum"] = stratum,
//     _["time"] = time,
//     _["nrisk"] = nrisk,
//     _["nevent"] = nevent,
//     _["ncensor"] = ncensor,
//     _["haz"] = haz,
//     _["varhaz"] = varhaz
//   );
//   
//   if (p > 0) result.push_back(gradhaz, "gradhaz");
//   
//   return result;
// }
// 
// 
// 
// // martingale residuals
// NumericVector f_resmart(int p, const NumericVector& par, void *ex) {
//   coxparams *param = (coxparams *) ex;
//   
//   const int nused = param->nused;
//   const int method = param->method;
//   const int n = static_cast<int>(param->tstop.size());
//   
//   // Precompute sexp(eta)
//   NumericVector exp_eta(nused);
//   for (int person = 0; person < nused; ++person) {
//     double val = param->offset[person];
//     for (int i = 0; i < p; ++i) {
//       val += par[i] * param->z(person,i);
//     }
//     exp_eta[person] = std::exp(val);
//   }
//   
//   double denom = 0.0;         // s0(beta,k,t)
//   double denom2 = 0.0;        // sum of weighted risks for deaths
//   double deadwt = 0.0;        // sum of weights for the deaths
//   int ndead = 0;              // number of deaths at this time point
//   
//   // initialize the residuals to the event indicators
//   NumericVector resid(n);
//   for (int person = 0; person < nused; ++person) {
//     resid[person] = param->event[person];
//   }
//   
//   int istrata = param->strata[0];
//   int i1 = 0; // index for removing out-of-risk subjects
//   int j0 = 0; // first person in the stratum
//   
//   // Loop through subjects
//   for (int person = 0; person < nused; ) {
//     if (param->strata[person] != istrata) { // hit a new stratum
//       istrata = param->strata[person]; // reset temporary variables
//       i1 = person;
//       j0 = person;
//       denom = 0.0;
//     }
//     
//     const double dtime = param->tstop[person];
//     
//     // process all persons tied at this dtime
//     int j1 = person;   // first person in the stratum with the tied time
//     for (; person < nused && param->tstop[person] == dtime &&
//          param->strata[person] == istrata; ++person) {
//       
//       const double w = param->weight[person];
//       const double r = w * exp_eta[person];
//       
//       if (param->event[person] == 0) {
//         denom += r;
//       } else {
//         ++ndead;
//         deadwt += w;
//         denom2 += r;
//       }
//     }
//     
//     int j2 = person - 1; // last person in the stratum with the tied time
//     
//     // Remove subjects leaving risk set
//     for (; i1 < nused; ++i1) {
//       const int p1 = param->order1[i1];
//       if (param->tstart[p1] < dtime || param->strata[p1] != istrata) break;
//       denom -= param->weight[p1] * exp_eta[p1];
//     }
//     
//     // Add contributions for deaths at this time
//     if (ndead > 0) {
//       denom += denom2;
//       
//       for (int j = j0; j <= j2; ++j) {
//         if (param->tstart[j] < dtime) {
//           double hazard;
//           if (method == 0 || ndead == 1) {
//             hazard = deadwt / denom;
//           } else {
//             hazard = 0.0;
//             const double meanwt = deadwt / ndead;
//             if (j < j1 || param->event[j] == 0) {
//               for (int i = 0; i < ndead; ++i) {
//                 hazard += meanwt /(denom - (i + 0.0)/ndead * denom2);
//               }
//             } else {
//               for (int i = 0; i < ndead; ++i) {
//                 hazard += (1 - (i + 0.0)/ndead) * meanwt /
//                   (denom - (i + 0.0)/ndead * denom2);
//               }
//             }
//           }
//           resid[j] -= hazard * exp_eta[j];
//         }
//       }
//       
//       // reset for the next death time
//       ndead = 0;
//       deadwt = denom2 = 0.0;
//     }
//   }
//   
//   return resid;
// }
// 
// 
// // score residual matrix
// NumericMatrix f_ressco_2(int p, const NumericVector& par, void *ex) {
//   coxparams *param = (coxparams *) ex;
//   
//   const int nused = param->nused;
//   const int method = param->method;
//   const int n = static_cast<int>(param->tstart.size());
//   
//   // Precompute exp(eta)
//   NumericVector exp_eta(nused);
//   for (int person = 0; person < nused; ++person) {
//     double val = param->offset[person];
//     for (int i = 0; i < p; ++i) {
//       val += par[i] * param->z(person,i);
//     }
//     exp_eta[person] = std::exp(val);
//   }
//   
//   NumericMatrix resid(n,p);   // residual matrix
//   NumericVector a(p);         // s1(beta,k,t)
//   NumericVector a2(p);        // sum of w*exp(zbeta)*z for the deaths
//   double denom = 0.0;         // s0(beta,k,t)
//   double denom2 = 0.0;        // sum of weighted risks for deaths
//   double deadwt = 0.0;        // sum of weights for the deaths
//   int ndead = 0;              // number of deaths at this time point
//   double cumhaz = 0.0;        // cumulative hazard
//   
//   NumericVector xhaz(p), mh1(p), mh2(p), mh3(p); // temp vectors
//   
//   int istrata = param->strata[0];
//   int i1 = 0; // index for removing out-of-risk subjects
//   
//   // Loop through subjects
//   for (int person = 0; person < nused; ) {
//     // Reset when entering a new stratum
//     if (param->strata[person] != istrata) {
//       istrata = param->strata[person];
//       
//       // first obs of a new stratum, finish off the prior stratum
//       for (; i1 < nused && param->order1[i1] < person; ++i1) {
//         const int p1 = param->order1[i1];
//         for (int i = 0; i < p; ++i) {
//           resid(p1,i) -= exp_eta[p1] * (param->z(p1,i) * cumhaz - xhaz[i]);
//         }
//       }
//       
//       denom = 0.0; // reset temporary variables
//       cumhaz = 0.0;
//       a.fill(0.0); xhaz.fill(0.0);
//     }
//     
//     const double dtime = param->tstop[person];
//     
//     // process all persons tied at this dtime
//     for (; person < nused && param->tstop[person] == dtime &&
//          param->strata[person] == istrata; ++person) {
//       
//       // initialize residuals to score[i] * (x[i] * cumhaz - xhaz), before
//       // updating cumhaz and xhaz
//       for (int i = 0; i < p; ++i) {
//         resid(person,i) = exp_eta[person] *
//           (param->z(person,i) * cumhaz - xhaz[i]);
//       }
//       
//       const double w = param->weight[person];
//       const double r = w * exp_eta[person];
//       
//       if (param->event[person] == 0) {
//         denom += r;
//         for (int i = 0; i < p; ++i) {
//           a[i] += r * param->z(person,i);
//         }
//       } else {
//         ++ndead;
//         deadwt += w;
//         denom2 += r;
//         for (int i = 0; i < p; ++i) {
//           a2[i] += r * param->z(person,i);
//         }
//       }
//     }
//     
//     // Remove subjects leaving risk set
//     for (; i1 < nused; ++i1) {
//       const int p1 = param->order1[i1];
//       if (param->tstart[p1] < dtime || param->strata[p1] != istrata) break;
//       
//       const double r = param->weight[p1] * exp_eta[p1];
//       denom -= r;
//       for (int i = 0; i < p; ++i) {
//         // finish the residual by subtracting score[i] * (x[i] * cumhaz - xhaz)
//         resid(p1,i) -= exp_eta[p1] * (param->z(p1,i) * cumhaz - xhaz[i]);
//         a[i] -= r * param->z(p1,i);
//       }
//     }
//     
//     // Add contributions for deaths at this time
//     if (ndead > 0) {
//       if (method == 0 || ndead == 1) { // Breslow or single event
//         denom += denom2;
//         const double hazard = deadwt/denom;
//         cumhaz += hazard;
//         for (int i = 0; i < p; ++i) {
//           a[i] += a2[i];
//           const double xbar = a[i]/denom;
//           xhaz[i] += xbar*hazard;
//           for (int j = person-1; j >= person - ndead; j--) {
//             resid(j,i) += param->z(j,i) - xbar;
//           }
//         }
//       } else {  // Efron method
//         for (int i = 0; i < p; ++i) {
//           mh1[i] = 0.0;
//           mh2[i] = 0.0;
//           mh3[i] = 0.0;
//         }
//         
//         const double meanwt = deadwt / ndead;
//         const double increment = denom2 / ndead;
//         
//         for (int k = 0; k < ndead; ++k) {
//           denom += increment;
//           const double hazard = meanwt/denom;
//           cumhaz += hazard;
//           const double downwt = (ndead - k - 1.0)/ndead;
//           for (int i = 0; i < p; ++i) {
//             a[i] += a2[i] / ndead;
//             const double xbar = a[i] / denom;
//             xhaz[i] += xbar * hazard;
//             mh1[i]  += hazard * downwt;
//             mh2[i]  += xbar * hazard * downwt;
//             mh3[i]  += xbar / ndead;
//           }
//         }
//         
//         for (int j = person-1; j >= person - ndead; j--) {
//           for (int i = 0; i < p; ++i) {
//             resid(j,i) += (param->z(j,i) - mh3[i]) +
//               exp_eta[j] * (param->z(j,i) * mh1[i] - mh2[i]);
//           }
//         }
//       }
//       
//       // Reset after processing deaths
//       ndead = 0;
//       deadwt = denom2 = 0.0;
//       a2.fill(0.0);
//     }
//   }
//   
//   // finish those remaining in the final stratum
//   for (; i1 < nused; ++i1) {
//     const int p1 = param->order1[i1];
//     for (int i = 0; i < p; ++i)
//       resid(p1,i) -= exp_eta[p1] * (param->z(p1,i) * cumhaz - xhaz[i]);
//   }
//   
//   return resid;
// }
// 
// 
// // schoenfeld residuals
// List f_ressch(int p, const NumericVector& par, void *ex) {
//   coxparams *param = (coxparams *) ex;
//   
//   const int nused = param->nused;
//   const int method = param->method;
//   
//   // Precompute exp(eta)
//   NumericVector exp_eta(nused);
//   for (int person = 0; person < nused; ++person) {
//     double val = param->offset[person];
//     for (int i = 0; i < p; ++i) {
//       val += par[i] * param->z(person, i);
//     }
//     exp_eta[person] = std::exp(val);
//   }
//   
//   NumericVector xbar(p);      // weighted mean covariate at this time
//   NumericVector a(p);         // s1(beta,k,t)
//   NumericVector a2(p);        // sum of w*exp(zbeta)*z for the deaths
//   double denom = 0.0;         // s0(beta,k,t)
//   double denom2 = 0.0;        // sum of weighted risks for the deaths
//   int ndead = 0;              // number of deaths at this time point
//   
//   int nevent = sum(param->event);   // total number of events
//   NumericMatrix resid(nevent, p);   // residual matrix
//   IntegerVector index(nevent);      // index of residuals
//   
//   int istrata = param->strata[0];
//   int i1 = 0; // index for removing out-of-risk subjects
//   int j = nevent;  // index the events in descending order
//   
//   // Loop through subjects
//   for (int person = 0; person < nused; ) {
//     if (param->strata[person] != istrata) { // hit a new stratum
//       istrata = param->strata[person]; // reset temporary variables
//       i1 = person;
//       denom = 0.0;
//       a.fill(0.0);
//     }
//     
//     const double dtime = param->tstop[person];
//     
//     // process all persons tied at this dtime
//     for (; person < nused && param->tstop[person] == dtime &&
//          param->strata[person] == istrata; ++person) {
//       
//       const double r = param->weight[person] * exp_eta[person];
//       
//       if (param->event[person] == 0) {
//         denom += r;
//         for (int i=0; i<p; ++i) {
//           a[i] += r * param->z(person,i);
//         }
//       } else {
//         j--;
//         resid(j,_) = param->z(person,_);
//         index[j] = person;
//         
//         ++ndead;
//         denom2 += r;
//         for (int i=0; i<p; ++i) {
//           a2[i] += r * param->z(person,i);
//         }
//       }
//     }
//     
//     // remove subjects no longer at risk
//     for (; i1 < nused; ++i1) {
//       const int p1 = param->order1[i1];
//       if (param->tstart[p1] < dtime || param->strata[p1] != istrata) break;
//       const double r = param->weight[p1] * exp_eta[p1];
//       denom -= r;
//       for (int i = 0; i < p; ++i) {
//         a[i] -= r * param->z(p1,i);
//       }
//     }
//     
//     // add to the main terms
//     if (ndead > 0) {
//       if (method == 0 || ndead == 1) {
//         denom += denom2;
//         for (int i = 0; i < p; ++i) {
//           a[i] += a2[i];
//           xbar[i] = a[i] / denom;
//         }
//       } else {
//         xbar.fill(0.0);
//         for (int k = 0; k < ndead; ++k) {
//           denom += denom2 / ndead;
//           for (int i = 0; i < p; ++i) {
//             a[i] += a2[i] / ndead;
//             xbar[i] += a[i] / denom;
//           }
//         }
//         xbar = xbar / ndead;
//       }
//       
//       for (int k = 0; k < ndead; ++k) {
//         for (int i = 0; i < p; ++i) {
//           resid(j+k,i) -= xbar[i];
//         }
//       }
//       
//       // reset for the next death time
//       ndead = 0;
//       denom2 = 0.0;
//       a2.fill(0.0);
//     }
//   }
//   
//   return List::create(
//     Named("resid") = resid,
//     Named("index") = index);
// }
// 
// 
// // [[Rcpp::export]]
// List phregcpp(const DataFrame data,
//               const StringVector& rep = "",
//               const StringVector& stratum = "",
//               const std::string time = "time",
//               const std::string time2 = "",
//               const std::string event = "event",
//               const StringVector& covariates = "",
//               const std::string weight = "",
//               const std::string offset = "",
//               const std::string id = "",
//               const std::string ties = "efron",
//               const NumericVector& init = NA_REAL,
//               const bool robust = false,
//               const bool est_basehaz = true,
//               const bool est_resid = true,
//               const bool firth = false,
//               const bool plci = false,
//               const double alpha = 0.05,
//               const int maxiter = 50,
//               const double eps = 1.0e-9) {
//   
//   int n = data.nrows();
//   int p = static_cast<int>(covariates.size());
//   if (p == 1 && (covariates[0] == "" || covariates[0] == "none")) {
//     p = 0;
//   }
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
//   bool has_stratum;
//   IntegerVector stratumn(n);
//   DataFrame u_stratum;
//   int p_stratum = static_cast<int>(stratum.size());
//   if (p_stratum == 1 && (stratum[0] == "" || stratum[0] == "none")) {
//     has_stratum = 0;
//     stratumn.fill(0);
//   } else {
//     List out = bygroup(data, stratum);
//     has_stratum = 1;
//     stratumn = out["index"];
//     u_stratum = DataFrame(out["lookup"]);
//   }
//   
//   bool has_time = hasVariable(data, time);
//   if (!has_time) stop("data must contain the time variable");
//   NumericVector timenz = data[time];
//   NumericVector timen = clone(timenz);
//   if (is_true(any(timen < 0))) {
//     stop("time must be nonnegative for each observation");
//   }
//   
//   bool has_time2 = hasVariable(data, time2);
//   NumericVector time2n(n);
//   if (has_time2) {
//     NumericVector time2nz = data[time2];
//     time2n = clone(time2nz);
//     if (is_true(any(time2n <= timen))) {
//       stop("time2 must be greater than time for each observation");
//     }
//   }
//   
//   bool has_event = hasVariable(data, event);
//   if (!has_event) stop("data must contain the event variable");
//   IntegerVector eventnz = data[event];
//   IntegerVector eventn = clone(eventnz);
//   if (is_true(any((eventn != 1) & (eventn != 0)))) {
//     stop("event must be 1 or 0 for each observation");
//   }
//   
//   NumericMatrix zn(n,p);
//   if (p > 0) {
//     for (int j=0; j<p; ++j) {
//       String zj = covariates[j];
//       if (!hasVariable(data, zj)) {
//         stop("data must contain the variables in covariates");
//       }
//       NumericVector u = data[zj];
//       for (int i=0; i<n; ++i) {
//         zn(i,j) = u[i];
//       }
//     }
//   }
//   
//   bool has_weight = hasVariable(data, weight);
//   NumericVector weightn(n, 1.0);
//   if (has_weight) {
//     NumericVector weightnz = data[weight];
//     weightn = clone(weightnz);
//     if (is_true(any(weightn <= 0))) {
//       stop("weight must be greater than 0");
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
//       stop("incorrect type for the id variable in the input data");
//     }
//   }
//   
//   if (robust && has_time2 && !has_id) {
//     stop("id is needed for counting process data with robust variance");
//   }
//   
//   std::string meth = ties;
//   std::for_each(meth.begin(), meth.end(), [](char & c) {
//     c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
//   });
//   
//   int method = meth == "efron" ? 1 : 0;
//   
//   // unify right censored data with counting process data
//   NumericVector tstartn(n), tstopn(n);
//   if (!has_time2) {
//     tstopn = timen;
//   } else {
//     tstartn = timen;
//     tstopn = time2n;
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
//     tstartn = tstartn[order];
//     tstopn = tstopn[order];
//     eventn = eventn[order];
//     weightn = weightn[order];
//     offsetn = offsetn[order];
//     idn = idn[order];
//     if (p > 0) zn = subset_matrix_by_row(zn, order);
//   }
//   
//   // exclude observations with missing values
//   LogicalVector sub(n,1);
//   for (int i=0; i<n; ++i) {
//     if (repn[i] == NA_INTEGER || stratumn[i] == NA_INTEGER ||
//         std::isnan(tstartn[i]) || std::isnan(tstopn[i]) ||
//         eventn[i] == NA_INTEGER || std::isnan(weightn[i]) ||
//         std::isnan(offsetn[i]) || idn[i] == NA_INTEGER) {
//       sub[i] = 0;
//     }
//     for (int j=0; j<p; ++j) {
//       if (std::isnan(zn(i,j))) sub[i] = 0;
//     }
//   }
//   
//   IntegerVector order = which(sub);
//   repn = repn[order];
//   stratumn = stratumn[order];
//   tstartn = tstartn[order];
//   tstopn = tstopn[order];
//   eventn = eventn[order];
//   weightn = weightn[order];
//   offsetn = offsetn[order];
//   idn = idn[order];
//   if (p > 0) zn = subset_matrix_by_row(zn, order);
//   n = sum(sub);
//   if (n == 0) stop("no observations without missing values");
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
//   NumericMatrix loglik(nreps,2);
//   NumericMatrix regloglik(nreps,2);
//   NumericVector scoretest(nreps);
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
//   // baseline hazards data set
//   int N = 2*n; // account for additional time 0 rows
//   IntegerVector drep(N), dstratum(N);
//   NumericVector dtime(N), dnrisk(N), dnevent(N), dncensor(N);
//   NumericVector dhaz(N), dvarhaz(N);
//   NumericMatrix dgradhaz(N,p);
//   
//   // martingale residuals
//   NumericVector resmart(n);
//   
//   // linear predictors
//   NumericVector linear_predictors(n);
//   
//   int n0 = 0; // number of rows in the baseline hazard data set
//   int bign0 = 0; // number of elements in the martingale residuals vector
//   double toler = 1e-12;
//   double zcrit = boost_qnorm(1-alpha/2,0,1,1,0);
//   double xcrit = zcrit * zcrit;
//   
//   for (int h=0; h<nreps; ++h) {
//     IntegerVector q1 = Range(idx[h], idx[h+1]-1);
//     int n1 = static_cast<int>(q1.size());
//     
//     IntegerVector stratum1 = stratumn[q1];
//     NumericVector tstart1 = tstartn[q1];
//     NumericVector tstop1 = tstopn[q1];
//     IntegerVector event1 = eventn[q1];
//     NumericVector weight1 = weightn[q1];
//     NumericVector offset1 = offsetn[q1];
//     IntegerVector id1 = idn[q1];
//     
//     NumericMatrix z1(n1,p);
//     if (p > 0) z1 = subset_matrix_by_row(zn, q1);
//     
//     nobs[h] = n1;
//     nevents[h] = sum(event1);
//     
//     if (nevents[h] == 0) {
//       if (p > 0) {
//         for (int i=0; i<p; ++i) {
//           int k = h*p+i;
//           rep0[k] = h;
//           par0[k] = covariates[i];
//           beta0[k] = NaN;
//           sebeta0[k] = 0;
//           rsebeta0[k] = 0;
//           z0[k] = NaN;
//           expbeta0[k] = NaN;
//           for (int j=0; j<p; ++j) {
//             vbeta0(k,j) = 0;
//             rvbeta0(k,j) = 0;
//           }
//           lb0[k] = NaN;
//           ub0[k] = NaN;
//           prob0[k] = NaN;
//           clparm0[k] = "Wald";
//         }
//       }
//       
//       // baseline hazard
//       if (est_basehaz) {
//         drep[n0] = h;
//         dstratum[n0] = 0;
//         dtime[n0] = 0;
//         dnrisk[n0] = n1;
//         dnevent[n0] = 0;
//         dncensor[n0] = 0;
//         dhaz[n0] = 0;
//         dvarhaz[n0] = 0;
//         if (p > 0) {
//           for (int i=0; i<p; ++i) {
//             dgradhaz(n0,i) = 0;
//           }
//         }
//         ++n0;
//       }
//       
//       // martingale residuals
//       if (est_resid) {
//         for (int i=0; i<n1; ++i) {
//           resmart[bign0+i] = 0;
//         }
//       }
//       
//       // linear predictors
//       for (int i=0; i<n1; ++i) {
//         linear_predictors[bign0+i] = offset1[i];
//       }
//       
//       bign0 += n1;
//       
//       continue;
//     }
//     
//     // sort by stratum
//     IntegerVector order0 = seq(0, n1-1);
//     std::sort(order0.begin(), order0.end(), [&](int i, int j) {
//       return stratum1[i] < stratum1[j];
//     });
//     
//     IntegerVector stratum1z = stratum1[order0];
//     NumericVector tstart1z = tstart1[order0];
//     NumericVector tstop1z = tstop1[order0];
//     IntegerVector event1z = event1[order0];
//     
//     // locate the first observation within each stratum
//     IntegerVector istratum(1,0);
//     for (int i=1; i<n1; ++i) {
//       if (stratum1z[i] != stratum1z[i-1]) {
//         istratum.push_back(i);
//       }
//     }
//     
//     int nstrata = static_cast<int>(istratum.size());
//     istratum.push_back(n1);
//     
//     // ignore subjects not at risk for any event time
//     IntegerVector ignore1z(n1);
//     for (int i=0; i<nstrata; ++i) {
//       IntegerVector q0 = Range(istratum[i], istratum[i+1]-1);
//       NumericVector tstart0 = tstart1z[q0];
//       NumericVector tstop0 = tstop1z[q0];
//       IntegerVector event0 = event1z[q0];
//       NumericVector etime = tstop0[event0==1];
//       etime = unique(etime);
//       etime.sort();
//       IntegerVector index1 = findInterval3(tstart0, etime, 0, 0, 0);
//       IntegerVector index2 = findInterval3(tstop0, etime, 0, 0, 0);
//       for (int j=istratum[i]; j<istratum[i+1]; ++j) {
//         int j0 = j-istratum[i];
//         if (index1[j0] == index2[j0]) { // no event in (tstart, tstop]
//           ignore1z[j] = 1;
//         } else {
//           ignore1z[j] = 0;
//         }
//       }
//     }
//     
//     IntegerVector ignore1(n1); // back to the original order
//     for (int i=0; i<n1; ++i) {
//       ignore1[order0[i]] = ignore1z[i];
//     }
//     
//     int nused = n1 - sum(ignore1); // number of used observations
//     
//     // sort by stopping time in descending order within each stratum
//     IntegerVector order2 = seq(0, n1-1);
//     std::sort(order2.begin(), order2.end(), [&](int i, int j) {
//       if (ignore1[i] != ignore1[j]) return ignore1[i] < ignore1[j];
//       if (stratum1[i] != stratum1[j]) return stratum1[i] < stratum1[j];
//       if (tstop1[i] != tstop1[j]) return tstop1[i] > tstop1[j];
//       return event1[i] < event1[j];
//     });
//     
//     IntegerVector stratum1a = stratum1[order2];
//     NumericVector tstart1a = tstart1[order2];
//     NumericVector tstop1a = tstop1[order2];
//     IntegerVector event1a = event1[order2];
//     NumericVector weight1a = weight1[order2];
//     NumericVector offset1a = offset1[order2];
//     IntegerVector id1a = id1[order2];
//     IntegerVector ignore1a = ignore1[order2];
//     NumericMatrix z1a;
//     if (p > 0) z1a = subset_matrix_by_row(z1, order2);
//     
//     // sort by starting time in descending order within each stratum
//     IntegerVector order1a = seq(0, n1-1);
//     std::sort(order1a.begin(), order1a.end(), [&](int i, int j) {
//       if (ignore1a[i] != ignore1a[j]) return ignore1a[i] < ignore1a[j];
//       if (stratum1a[i] != stratum1a[j]) return stratum1a[i] < stratum1a[j];
//       return tstart1a[i] > tstart1a[j];
//     });
//     
//     coxparams param = {nused, stratum1a, tstart1a, tstop1a, event1a,
//                        weight1a, offset1a, z1a, order1a, method};
//     
//     NumericVector bint(p);
//     List derint = f_der_2(p, bint, &param, firth);
//     
//     NumericVector b(p);
//     NumericMatrix vb(p,p);
//     List out;
//     
//     if (p > 0) {
//       IntegerVector colfit = seq(0,p-1);
//       if (is_false(any(is_na(init))) && init.size() == p) {
//         out = phregloop(p, init, &param, maxiter, eps, firth, colfit, p);
//       } else {
//         out = phregloop(p, bint, &param, maxiter, eps, firth, colfit, p);
//       }
//       
//       bool fail = out["fail"];
//       if (fail) warning("The algorithm in phregr did not converge");
//       
//       b = out["coef"];
//       vb = as<NumericMatrix>(out["var"]);
//       
//       NumericVector seb(p);
//       for (int j=0; j<p; ++j) {
//         seb[j] = std::sqrt(vb(j,j));
//       }
//       
//       for (int i=0; i<p; ++i) {
//         int k = h*p+i;
//         rep0[k] = h;
//         par0[k] = covariates[i];
//         beta0[k] = b[i];
//         sebeta0[k] = seb[i];
//         for (int j=0; j<p; ++j) {
//           vbeta0(k,j) = vb(i,j);
//         }
//       }
//       
//       // score statistic
//       
//       NumericVector scorebint = firth ? derint["regscore"] : derint["score"];
//       NumericMatrix infobint = as<NumericMatrix>(derint["imat"]);
//       NumericMatrix vbint = invsympd(infobint, p, toler);
//       for (int i=0; i<p; ++i) {
//         for (int j=0; j<p; ++j) {
//           scoretest[h] += scorebint[i]*vbint(i,j)*scorebint[j];
//         }
//       }
//       
//       niter[h] = out["iter"];
//       fails[h] = out["fail"];
//       
//       // robust variance estimates
//       NumericVector rseb(p);  // robust standard error for betahat
//       if (robust) {
//         NumericMatrix ressco = f_ressco_2(p, b, &param);
//         
//         int nr; // number of rows in the score residual matrix
//         if (!has_id) {
//           for (int i=0; i<n1; ++i) {
//             for (int j=0; j<p; ++j) {
//               ressco(i,j) = weight1a[i]*ressco(i,j);
//             }
//           }
//           nr = n1;
//         } else { // need to sum up score residuals by id
//           IntegerVector order = seq(0, n1-1);
//           std::sort(order.begin(), order.end(), [&](int i, int j) {
//             return id1a[i] < id1a[j];
//           });
//           
//           IntegerVector id2 = id1a[order];
//           IntegerVector idx(1,0);
//           for (int i=1; i<n1; ++i) {
//             if (id2[i] != id2[i-1]) {
//               idx.push_back(i);
//             }
//           }
//           
//           int nids = static_cast<int>(idx.size());
//           idx.push_back(n1);
//           
//           NumericVector weight2 = weight1a[order];
//           
//           NumericMatrix ressco2(nids,p);
//           for (int i=0; i<nids; ++i) {
//             for (int j=0; j<p; ++j) {
//               for (int k=idx[i]; k<idx[i+1]; ++k) {
//                 ressco2(i,j) += weight2[k]*ressco(order[k],j);
//               }
//             }
//           }
//           
//           ressco = ressco2;  // update the score residuals
//           nr = nids;
//         }
//         
//         NumericMatrix D(nr,p); // DFBETA
//         for (int i=0; i<nr; ++i) {
//           for (int j=0; j<p; ++j) {
//             for (int k=0; k<p; ++k) {
//               D(i,j) += ressco(i,k)*vb(k,j);
//             }
//           }
//         }
//         
//         NumericMatrix rvb(p,p); // robust variance matrix for betahat
//         for (int j=0; j<p; ++j) {
//           for (int k=0; k<p; ++k) {
//             for (int i=0; i<nr; ++i) {
//               rvb(j,k) += D(i,j)*D(i,k);
//             }
//           }
//         }
//         
//         for (int i=0; i<p; ++i) {
//           rseb[i] = std::sqrt(rvb(i,i));
//         }
//         
//         for (int i=0; i<p; ++i) {
//           int k = h*p+i;
//           rsebeta0[k] = rseb[i];
//           for (int j=0; j<p; ++j) {
//             rvbeta0(k,j) = rvb(i,j);
//           }
//         }
//       }
//       
//       // profile likelihood confidence interval for regression coefficients
//       NumericVector lb(p), ub(p), prob(p);
//       StringVector clparm(p);
//       
//       if (plci) {
//         double lmax = out["loglik"];
//         double l0 = lmax - 0.5*xcrit;
//         
//         for (int k=0; k<p; ++k) {
//           lb[k] = phregplloop(p, b, &param, maxiter, eps, firth, k, -1, l0);
//           ub[k] = phregplloop(p, b, &param, maxiter, eps, firth, k, 1, l0);
//           
//           IntegerVector colfit1(p-1);
//           for (int i = 0, j = 0; i < p; ++i) {
//             if (i == k) continue;
//             colfit1[j++] = i;
//           }
//           
//           NumericVector b0(p);
//           List out0 = phregloop(p, b0, &param, maxiter,eps,firth,colfit1,p-1);
//           double lmax0 = out0["loglik"];
//           prob[k] = boost_pchisq(-2*(lmax0 - lmax), 1, 0, 0);
//           clparm[k] = "PL";
//         }
//       } else {
//         for (int k=0; k<p; ++k) {
//           if (!robust) {
//             lb[k] = b[k] - zcrit*seb[k];
//             ub[k] = b[k] + zcrit*seb[k];
//             prob[k] = boost_pchisq(std::pow(b[k]/seb[k], 2), 1, 0, 0);
//           } else {
//             lb[k] = b[k] - zcrit*rseb[k];
//             ub[k] = b[k] + zcrit*rseb[k];
//             prob[k] = boost_pchisq(std::pow(b[k]/rseb[k], 2), 1, 0, 0);
//           }
//           clparm[k] = "Wald";
//         }
//       }
//       
//       for (int i=0; i<p; ++i) {
//         int k = h*p+i;
//         lb0[k] = lb[i];
//         ub0[k] = ub[i];
//         prob0[k] = prob[i];
//         clparm0[k] = clparm[i];
//       }
//     }
//     
//     // log-likelihoods
//     if (p > 0 && firth) {
//       loglik(h,0) = derint["loglik"];
//       loglik(h,1) = out["loglik"];
//       regloglik(h,0) = derint["regloglik"];
//       regloglik(h,1) = out["regloglik"];
//     } else {
//       loglik(h,0) = derint["loglik"];
//       loglik(h,1) = p > 0 ? out["loglik"] : derint["loglik"];
//     }
//     
//     // estimate baseline hazard
//     if (est_basehaz) {
//       // prepare the data for estimating baseline hazards at all time points
//       
//       // sort by stopping time in descending order within each stratum
//       IntegerVector order3 = seq(0, n1-1);
//       std::sort(order3.begin(), order3.end(), [&](int i, int j) {
//         if (stratum1[i] != stratum1[j]) return stratum1[i] > stratum1[j];
//         return tstop1[i] > tstop1[j];
//       });
//       
//       IntegerVector stratum1b = stratum1[order3];
//       NumericVector tstart1b = tstart1[order3];
//       NumericVector tstop1b = tstop1[order3];
//       IntegerVector event1b = event1[order3];
//       NumericVector weight1b = weight1[order3];
//       NumericVector offset1b = offset1[order3];
//       NumericMatrix z1b(n1,p);
//       if (p > 0) z1b = subset_matrix_by_row(z1, order3);
//       
//       // sort by starting time in descending order within each stratum
//       IntegerVector order1b = seq(0, n1-1);
//       std::sort(order1b.begin(), order1b.end(), [&](int i, int j) {
//         if (stratum1b[i] != stratum1b[j]) return stratum1b[i] > stratum1b[j];
//         return tstart1b[i] > tstart1b[j];
//       });
//       
//       coxparams paramb = {n1, stratum1b, tstart1b, tstop1b, event1b,
//                           weight1b, offset1b, z1b, order1b, method};
//       
//       List basehaz1 = f_basehaz(p, b, &paramb);
//       
//       IntegerVector dstratum1 = basehaz1["stratum"];
//       NumericVector dtime1 = basehaz1["time"];
//       NumericVector dnrisk1 = basehaz1["nrisk"];
//       NumericVector dnevent1 = basehaz1["nevent"];
//       NumericVector dncensor1 = basehaz1["ncensor"];
//       NumericVector dhaz1 = basehaz1["haz"];
//       NumericVector dvarhaz1 = basehaz1["varhaz"];
//       int J = static_cast<int>(dstratum1.size());
//       
//       // add to output data frame
//       for (int j=0; j<J; ++j) {
//         int k = n0 + j;
//         drep[k] = h;
//         dstratum[k] = dstratum1[j];
//         dtime[k] = dtime1[j];
//         dnrisk[k] = dnrisk1[j];
//         dnevent[k] = dnevent1[j];
//         dncensor[k] = dncensor1[j];
//         dhaz[k] = dhaz1[j];
//         dvarhaz[k] = dvarhaz1[j];
//         
//         if (p > 0) {
//           NumericMatrix dgradhaz1 = basehaz1["gradhaz"];
//           for (int i=0; i<p; ++i) {
//             dgradhaz(k,i) = dgradhaz1(j,i);
//           }
//         }
//       }
//       
//       n0 += J;
//     }
//     
//     // martingale residuals
//     if (est_resid) {
//       NumericVector resid = f_resmart(p, b, &param);
//       
//       for (int i=0; i<n1; ++i) {
//         resmart[bign0 + order2[i]] = resid[i];
//       }
//     }
//     
//     // linear predictors
//     for (int i=0; i<n1; ++i) {
//       linear_predictors[bign0 + order2[i]] = offset1a[i];
//       if (p > 0) {
//         for (int j=0; j<p; ++j) {
//           linear_predictors[bign0 + order2[i]] += b[j]*z1a(i,j);
//         }
//       }
//     }
//     
//     bign0 += n1;
//   }
//   
//   if (est_basehaz) {
//     IntegerVector sub = Range(0, n0 - 1);
//     drep = drep[sub];
//     dstratum = dstratum[sub];
//     dtime = dtime[sub];
//     dnrisk = dnrisk[sub];
//     dnevent = dnevent[sub];
//     dncensor = dncensor[sub];
//     dhaz = dhaz[sub];
//     dvarhaz = dvarhaz[sub];
//     if (p > 0) dgradhaz = subset_matrix_by_row(dgradhaz, sub);
//   }
//   
//   // prepare the output data sets
//   List sumstat = List::create(
//     _["n"] = nobs,
//     _["nevents"] = nevents,
//     _["loglik0"] = loglik(_,0),
//     _["loglik1"] = loglik(_,1),
//     _["scoretest"] = scoretest,
//     _["niter"] = niter,
//     _["ties"] = meth,
//     _["p"] = p,
//     _["robust"] = robust,
//     _["firth"] = firth,
//     _["fail"] = fails);
//   
//   if (p > 0 && firth) {
//     sumstat.push_back(regloglik(_,0), "loglik0_unpenalized");
//     sumstat.push_back(regloglik(_,1), "loglik1_unpenalized");
//   }
//   
//   if (has_rep) {
//     for (int i = 0; i < p_rep; ++i) {
//       std::string s = as<std::string>(rep[i]);
//       SEXP col = u_rep[s];
//       SEXPTYPE col_type = TYPEOF(col);
//       if (col_type == INTSXP) {
//         IntegerVector v = col;
//         sumstat.push_back(v[rep01], s);
//       } else if (col_type == REALSXP) {
//         NumericVector v = col;
//         sumstat.push_back(v[rep01], s);
//       } else if (col_type == STRSXP) {
//         StringVector v = col;
//         sumstat.push_back(v[rep01], s);
//       } else {
//         stop("Unsupported type for rep variable" + s);
//       }
//     }
//   }
//   
//   List result = List::create(
//     _["sumstat"] = as<DataFrame>(sumstat)
//   );
//   
//   
//   if (p > 0) {
//     expbeta0 = std::exp(beta0);
//     if (!robust) z0 = beta0/sebeta0;
//     else z0 = beta0/rsebeta0;
//     
//     List parest = List::create(
//       _["param"] = par0,
//       _["beta"] = beta0,
//       _["sebeta"] = robust ? rsebeta0 : sebeta0,
//       _["z"] = z0,
//       _["expbeta"] = expbeta0,
//       _["vbeta"] = robust ? rvbeta0 : vbeta0,
//       _["lower"] = lb0,
//       _["upper"] = ub0,
//       _["p"] = prob0,
//       _["method"] = clparm0);
//     
//     if (robust) {
//       parest.push_back(sebeta0, "sebeta_naive");
//       parest.push_back(vbeta0, "vbeta_naive");
//     }
//     
//     if (has_rep) {
//       for (int i=0; i<p_rep; ++i) {
//         std::string s = as<std::string>(rep[i]);
//         SEXP col = u_rep[s];
//         SEXPTYPE col_type = TYPEOF(col);
//         if (col_type == INTSXP) {
//           IntegerVector v = col;
//           parest.push_back(v[rep0], s);
//         } else if (col_type == REALSXP) {
//           NumericVector v = col;
//           parest.push_back(v[rep0], s);
//         } else if (col_type == STRSXP) {
//           StringVector v = col;
//           parest.push_back(v[rep0], s);
//         } else {
//           stop("Unsupported type for rep variable" + s);
//         }
//       }
//     }
//     
//     result.push_back(as<DataFrame>(parest), "parest");
//   }
//   
//   
//   if (est_basehaz) {
//     List basehaz = List::create(
//       _["time"] = dtime,
//       _["nrisk"] = dnrisk,
//       _["nevent"] = dnevent,
//       _["ncensor"] = dncensor,
//       _["haz"] = dhaz,
//       _["varhaz"] = dvarhaz
//     );
//     
//     if (p > 0) {
//       basehaz.push_back(dgradhaz, "gradhaz");
//     }
//     
//     if (has_stratum) {
//       for (int i = 0; i < p_stratum; ++i) {
//         std::string s = as<std::string>(stratum[i]);
//         SEXP col = u_stratum[s];
//         SEXPTYPE col_type = TYPEOF(col);
//         if (col_type == INTSXP) {
//           IntegerVector v = col;
//           basehaz.push_back(v[dstratum], s);
//         } else if (col_type == REALSXP) {
//           NumericVector v = col;
//           basehaz.push_back(v[dstratum], s);
//         } else if (col_type == STRSXP) {
//           StringVector v = col;
//           basehaz.push_back(v[dstratum], s);
//         } else {
//           stop("Unsupported type for stratum variable" + s);
//         }
//       }
//     }
//     
//     if (has_rep) {
//       for (int i = 0; i < p_rep; ++i) {
//         std::string s = as<std::string>(rep[i]);
//         SEXP col = u_rep[s];
//         SEXPTYPE col_type = TYPEOF(col);
//         if (col_type == INTSXP) {
//           IntegerVector v = col;
//           basehaz.push_back(v[drep], s);
//         } else if (col_type == REALSXP) {
//           NumericVector v = col;
//           basehaz.push_back(v[drep], s);
//         } else if (col_type == STRSXP) {
//           StringVector v = col;
//           basehaz.push_back(v[drep], s);
//         } else {
//           stop("Unsupported type for rep variable" + s);
//         }
//       }
//     }
//     
//     result.push_back(as<DataFrame>(basehaz), "basehaz");
//   }
//   
//   
//   if (est_resid) {
//     result.push_back(resmart, "residuals");
//   }
//   
//   result.push_back(linear_predictors, "linear_predictors");
//   
//   return result;
// }
// 
// 
// // [[Rcpp::export]]
// DataFrame survfit_phregcpp(const int p,
//                            const NumericVector& beta,
//                            const NumericMatrix& vbeta,
//                            DataFrame basehaz,
//                            DataFrame newdata,
//                            const StringVector& covariates = "",
//                            const StringVector& stratum = "",
//                            const std::string offset = "",
//                            const std::string id = "",
//                            const std::string tstart = "",
//                            const std::string tstop = "",
//                            const bool sefit = true,
//                            const String conftype = "log-log",
//                            const double conflev = 0.95) {
//   
//   int n0 = basehaz.nrows();
//   int n = newdata.nrows();
//   int nvar = static_cast<int>(covariates.size());
//   
//   std::string ct = conftype;
//   std::for_each(ct.begin(), ct.end(), [](char & c) {
//     c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
//   });
//   
//   if (!(ct=="none" || ct=="plain" || ct=="log" || ct=="log-log" ||
//       ct=="logit" || ct=="arcsin")) {
//     stop("conftype must be none, plain, log, log-log, logit, or arcsin");
//   }
//   
//   if (conflev <= 0 || conflev >= 1) {
//     stop("conflev must lie between 0 and 1");
//   }
//   
//   double zcrit = boost_qnorm((1+conflev)/2,0,1,1,0);
//   
//   IntegerVector stratumn0(n0);
//   DataFrame u_stratum0;
//   IntegerVector nlevels;
//   List lookups;
//   int p_stratum = static_cast<int>(stratum.size());
//   if (p_stratum == 1 && (stratum[0] == "" || stratum[0] == "none")) {
//     stratumn0.fill(0);
//   } else {
//     List out = bygroup(basehaz, stratum);
//     stratumn0 = out["index"];
//     u_stratum0 = DataFrame(out["lookup"]);
//     nlevels = out["nlevels"];
//     lookups = out["lookups"];
//   }
//   
//   bool nullmodel = (p == 0 || (nvar == 1 &&
//                     (covariates[0] == "" || covariates[0] == "none")));
//   
//   if (!nullmodel && nvar != p) {
//     stop("incorrect number of covariates for the Cox model");
//   }
//   
//   NumericMatrix zn(n,p);
//   for (int j=0; j<p; ++j) {
//     String zj = covariates[j];
//     if (!hasVariable(newdata, zj)) {
//       stop("newdata must contain the variables in covariates");
//     }
//     NumericVector u = newdata[zj];
//     for (int i=0; i<n; ++i) {
//       zn(i,j) = u[i];
//     }
//   }
//   
//   bool has_stratum;
//   IntegerVector stratumn(n);
//   if (p_stratum == 1 && (stratum[0] == "" || stratum[0] == "none")) {
//     has_stratum = false;
//     stratumn.fill(0);
//   } else {
//     has_stratum = true;
//     // match stratum in newdata to stratum in basehaz
//     int orep = u_stratum0.nrow();
//     for (int i=0; i<p_stratum; ++i) {
//       orep /= nlevels[i];
//       std::string s = as<std::string>(stratum[i]);
//       IntegerVector idx;
//       SEXP col = newdata[s];
//       SEXPTYPE col_type = TYPEOF(col);
//       if (col_type == LGLSXP || col_type == INTSXP) {
//         IntegerVector v = col;
//         IntegerVector w = lookups[i];
//         idx = match(v, w) - 1;
//       } else if (col_type == REALSXP) {
//         NumericVector v = col;
//         NumericVector w = lookups[i];
//         idx = match(v, w) - 1;
//       } else if (col_type == STRSXP) {
//         StringVector v = col;
//         StringVector w = lookups[i];
//         idx = match(v, w) - 1;
//       } else {
//         stop("Unsupported type for stratum variable: " + s);
//       }
//       
//       stratumn = stratumn + idx * orep;
//     }
//   }
//   
//   bool has_offset = hasVariable(newdata, offset);
//   NumericVector offsetn(n);
//   if (has_offset) {
//     NumericVector offsetnz = newdata[offset];
//     offsetn = clone(offsetnz);
//   }
//   
//   NumericVector time0 = basehaz["time"];
//   NumericVector nrisk0 = basehaz["nrisk"];
//   NumericVector nevent0 = basehaz["nevent"];
//   NumericVector ncensor0 = basehaz["ncensor"];
//   NumericVector haz0 = basehaz["haz"];
//   NumericVector vhaz0 = basehaz["varhaz"];
//   NumericMatrix ghaz0(n0,p);
//   for (int j=0; j<p; ++j) {
//     std::string col_name = "gradhaz";
//     if (p>1) col_name += "." + std::to_string(j+1);
//     NumericVector u = basehaz[col_name];
//     ghaz0(_,j) = u;
//   }
//   
//   // create the numeric id variable
//   bool has_id = hasVariable(newdata, id);
//   IntegerVector idn(n);
//   IntegerVector idwi;
//   NumericVector idwn;
//   StringVector idwc;
//   if (!has_id) {
//     idn = seq(0,n-1);
//   } else { // input data has the counting process style of input
//     SEXP col = newdata[id];
//     SEXPTYPE col_type = TYPEOF(col);
//     if (col_type == INTSXP) {
//       IntegerVector idv = col;
//       idwi = unique(idv);
//       idwi.sort();
//       idn = match(idv, idwi) - 1;
//     } else if (col_type == REALSXP) {
//       NumericVector idv = col;
//       idwn = unique(idv);
//       idwn.sort();
//       idn = match(idv, idwn) - 1;
//     } else if (col_type == STRSXP) {
//       StringVector idv = col;
//       idwc = unique(idv);
//       idwc.sort();
//       idn = match(idv, idwc) - 1;
//     } else {
//       stop("incorrect type for the id variable in newdata");
//     }
//   }
//   
//   // unify right-censoring data with counting process data
//   NumericVector tstartn(n), tstopn(n);
//   if (!has_id) { // right-censored data
//     tstartn.fill(0.0);
//     double maxt0 = max(time0) + 1.0;
//     tstopn.fill(maxt0);
//   } else {
//     bool has_tstart = hasVariable(newdata, tstart);
//     if (!has_tstart) stop("newdata must contain the tstart variable");
//     NumericVector tstartnz = newdata[tstart];
//     tstartn = clone(tstartnz);
//     if (is_true(any(tstartn < 0))) {
//       stop("tstart must be nonnegative for each observation");
//     }
//     
//     bool has_tstop = hasVariable(newdata, tstop);
//     if (!has_tstop) stop("newdata must contain the tstop variable");
//     NumericVector tstopnz = newdata[tstop];
//     tstopn = clone(tstopnz);
//     if (is_true(any(tstopn <= tstartn))) {
//       stop("tstop must be greater than tstart for each observation");
//     }
//   }
//   
//   // order data by id and tstop, assuming consecutive intervals
//   IntegerVector order = seq(0, n-1);
//   std::sort(order.begin(), order.end(), [&](int i, int j) {
//     if (idn[i] != idn[j]) return idn[i] < idn[j];
//     return tstopn[i] < tstopn[j];
//   });
//   
//   idn = idn[order];
//   stratumn = stratumn[order];
//   tstartn = tstartn[order];
//   tstopn = tstopn[order];
//   offsetn = offsetn[order];
//   if (p > 0) zn = subset_matrix_by_row(zn, order);
//   
//   // exclude observations with missing values
//   LogicalVector sub(n,1);
//   for (int i=0; i<n; ++i) {
//     if (stratumn[i] == NA_INTEGER || idn[i] == NA_INTEGER ||
//         std::isnan(tstartn[i]) || std::isnan(tstopn[i]) ||
//         std::isnan(offsetn[i])) {
//       sub[i] = 0;
//     }
//     for (int j=0; j<p; ++j) {
//       if (std::isnan(zn(i,j))) sub[i] = 0;
//     }
//   }
//   
//   order = which(sub);
//   stratumn = stratumn[order];
//   offsetn = offsetn[order];
//   idn = idn[order];
//   tstartn = tstartn[order];
//   tstopn = tstopn[order];
//   if (p > 0) zn = subset_matrix_by_row(zn, order);
//   n = sum(sub);
//   if (n == 0) stop("no observations left after removing missing values");
//   
//   // risk score
//   NumericVector risk(n);
//   for (int i=0; i<n; ++i) {
//     double val = offsetn[i];
//     for (int j=0; j<p; ++j) {
//       val += beta[j]*zn(i,j);
//     }
//     risk[i] = std::exp(val);
//   }
//   
//   // count number of observations for each id
//   IntegerVector idx(1,0);
//   for (int i=1; i<n; ++i) {
//     if (idn[i] != idn[i-1]) {
//       idx.push_back(i);
//     }
//   }
//   
//   int nids = static_cast<int>(idx.size());
//   idx.push_back(n);
//   
//   int N = nids*n0; // upper bound on the number of rows in the output
//   NumericVector time(N);
//   NumericVector nrisk(N), nevent(N), ncensor(N);
//   NumericVector cumhaz(N), vcumhaz(N), secumhaz(N);
//   IntegerVector strata(N);
//   NumericMatrix z(N,p);
//   IntegerVector ids(N);
//   
//   // process by id
//   int l = 0;
//   for (int h=0; h<nids; ++h) {
//     IntegerVector q1 = Range(idx[h], idx[h+1]-1);
//     int n1 = static_cast<int>(q1.size());
//     
//     IntegerVector id2 = idn[q1];
//     IntegerVector stratum2 = stratumn[q1];
//     NumericVector tstart2 = tstartn[q1];
//     NumericVector tstop2 = tstopn[q1];
//     NumericVector risk2 = risk[q1];
//     NumericMatrix z2 = subset_matrix_by_row(zn, q1);
//     tstop2.push_front(tstart2[0]);
//     
//     // match the stratum in basehaz
//     IntegerVector idx1 = which(stratumn0 == stratum2[0]);
//     NumericVector time01 = time0[idx1];
//     
//     // left-open and right-closed intervals containing the event time
//     IntegerVector idx2 = findInterval3(time01, tstop2, 0, 0, 1);
//     IntegerVector sub = which((idx2 >= 1) & (idx2 <= n1));
//     int m1 = sub.size();
//     
//     if (m1 != 0) {
//       IntegerVector idx3 = idx1[sub];
//       NumericVector time1 = time0[idx3];
//       NumericVector nrisk1 = nrisk0[idx3];
//       NumericVector nevent1 = nevent0[idx3];
//       NumericVector ncensor1 = ncensor0[idx3];
//       NumericVector haz1 = haz0[idx3];
//       
//       IntegerVector idx4 = idx2[sub];
//       idx4 = idx4 - 1; // change to 0-1 indexing
//       
//       // cumulative hazards
//       for (int i=0; i<m1; ++i) {
//         int r = l + i;
//         time[r] = time1[i];
//         nrisk[r] = nrisk1[i];
//         nevent[r] = nevent1[i];
//         ncensor[r] = ncensor1[i];
//         
//         int k = idx4[i];
//         ids[r] = id2[k];
//         strata[r] = stratum2[k];
//         for (int j=0; j<p; ++j) {
//           z(r,j) = z2(k,j);
//         }
//         
//         if (i==0) {
//           cumhaz[r] = haz1[i]*risk2[k];
//         } else {
//           cumhaz[r] = cumhaz[r-1] + haz1[i]*risk2[k];
//         }
//       }
//       
//       if (sefit) {
//         NumericVector vhaz1 = vhaz0[idx3];
//         NumericMatrix ghaz1(m1,p);
//         for (int j=0; j<p; ++j) {
//           for (int i=0; i<m1; ++i) {
//             ghaz1(i,j) = ghaz0(idx3[i],j);
//           }
//         }
//         
//         NumericMatrix a(m1,p);
//         for (int j=0; j<p; ++j) {
//           for (int i=0; i<m1; ++i) {
//             int k = idx4[i];
//             if (i==0) {
//               a(i,j) = (haz1[i]*z2(k,j) - ghaz1(i,j))*risk2[k];
//             } else {
//               a(i,j) = a(i-1,j) + (haz1[i]*z2(k,j) - ghaz1(i,j))*risk2[k];
//             }
//           }
//         }
//         
//         // calculate the first component of variance
//         for (int i=0; i<m1; ++i) {
//           int r = l + i;
//           int k = idx4[i];
//           if (i==0) {
//             vcumhaz[r] = vhaz1[i]*risk2[k]*risk2[k];
//           } else {
//             vcumhaz[r] = vcumhaz[r-1] + vhaz1[i]*risk2[k]*risk2[k];
//           }
//         }
//         
//         // add the second component of variance
//         for (int i=0; i<m1; ++i) {
//           int r = l + i;
//           for (int j=0; j<p; ++j) {
//             for (int k=0; k<p; ++k) {
//               vcumhaz[r] += a(i,j)*vbeta(j,k)*a(i,k);
//             }
//           }
//           secumhaz[r] = std::sqrt(vcumhaz[r]);
//         }
//       }
//       
//       l += m1;
//     }
//   }
//   
//   IntegerVector sub2 = Range(0,l-1);
//   time = time[sub2];
//   nrisk = nrisk[sub2];
//   nevent = nevent[sub2];
//   ncensor = ncensor[sub2];
//   cumhaz = cumhaz[sub2];
//   secumhaz = secumhaz[sub2];
//   strata = strata[sub2];
//   z = subset_matrix_by_row(z, sub2);
//   ids = ids[sub2];
//   
//   NumericVector surv = std::exp(-cumhaz);
//   
//   DataFrame result = DataFrame::create(
//     _["time"] = time,
//     _["nrisk"] = nrisk,
//     _["nevent"] = nevent,
//     _["ncensor"] = ncensor,
//     _["cumhaz"] = cumhaz,
//     _["surv"] = surv);
//   
//   if (sefit) {
//     NumericVector sesurv = surv*secumhaz;
//     
//     NumericVector lower(l), upper(l);
//     for (int i=0; i<l; ++i) {
//       NumericVector ci = fsurvci(surv[i], sesurv[i], ct, zcrit);
//       lower[i] = ci[0];
//       upper[i] = ci[1];
//     }
//     
//     result.push_back(sesurv, "sesurv");
//     result.push_back(lower, "lower");
//     result.push_back(upper, "upper");
//     result.push_back(conflev, "conflev");
//     result.push_back(ct, "conftype");
//   }
//   
//   for (int j=0; j<p; ++j) {
//     NumericVector u = z(_,j);
//     String zj = covariates[j];
//     result.push_back(u, zj);
//   }
//   
//   if (has_stratum) {
//     for (int i=0; i<p_stratum; ++i) {
//       std::string s = as<std::string>(stratum[i]);
//       SEXP col = u_stratum0[s];
//       SEXPTYPE col_type = TYPEOF(col);
//       if (col_type == INTSXP) {
//         IntegerVector v = col;
//         result.push_back(v[strata], s);
//       } else if (col_type == REALSXP) {
//         NumericVector v = col;
//         result.push_back(v[strata], s);
//       } else if (col_type == STRSXP) {
//         StringVector v = col;
//         result.push_back(v[strata], s);
//       } else {
//         stop("Unsupported type for stratum variable: " + s);
//       }
//     }
//   }
//   
//   if (has_id) {
//     SEXP col = newdata[id];
//     SEXPTYPE col_type = TYPEOF(col);
//     if (col_type == INTSXP) {
//       result.push_back(idwi[ids], id);
//     } else if (col_type == REALSXP) {
//       result.push_back(idwn[ids], id);
//     } else if (col_type == STRSXP) {
//       result.push_back(idwc[ids], id);
//     } else {
//       stop("incorrect type for the id variable in newdata");
//     }
//   }
//   
//   return result;
// }
// 
// 
// // [[Rcpp::export]]
// List residuals_phregcpp(const int p,
//                         const NumericVector& beta,
//                         const NumericMatrix& vbeta,
//                         const NumericVector& resmart,
//                         DataFrame data,
//                         const StringVector& stratum = "",
//                         const std::string time = "time",
//                         const std::string time2 = "",
//                         const std::string event = "event",
//                         const StringVector& covariates = "",
//                         const std::string weight = "",
//                         const std::string offset = "",
//                         const std::string id = "",
//                         const std::string ties = "efron",
//                         const std::string type = "schoenfeld",
//                         const bool collapse = false,
//                         const bool weighted = false) {
//   
//   int n = data.nrows();
//   
//   IntegerVector stratumn(n);
//   DataFrame u_stratum;
//   int p_stratum = static_cast<int>(stratum.size());
//   if (p_stratum == 1 && (stratum[0] == "" || stratum[0] == "none")) {
//     stratumn.fill(0);
//   } else {
//     List out = bygroup(data, stratum);
//     stratumn = out["index"];
//     u_stratum = out["lookup"];
//   }
//   
//   bool has_time = hasVariable(data, time);
//   if (!has_time) stop("data must contain the time variable");
//   NumericVector timenz = data[time];
//   NumericVector timen = clone(timenz);
//   if (is_true(any(timen < 0))) {
//     stop("time must be nonnegative for each observation");
//   }
//   
//   bool has_time2 = hasVariable(data, time2);
//   NumericVector time2n(n);
//   if (has_time2) {
//     NumericVector time2nz = data[time2];
//     time2n = clone(time2nz);
//     if (is_true(any(time2n <= timen))) {
//       stop("time2 must be greater than time for each observation");
//     }
//   }
//   
//   bool has_event = hasVariable(data, event);
//   if (!has_event) stop("data must contain the event variable");
//   IntegerVector eventnz = data[event];
//   IntegerVector eventn = clone(eventnz);
//   if (is_true(any((eventn != 1) & (eventn != 0)))) {
//     stop("event must be 1 or 0 for each observation");
//   }
//   
//   NumericMatrix zn(n,p);
//   if (p > 0) {
//     for (int j=0; j<p; ++j) {
//       String zj = covariates[j];
//       if (!hasVariable(data, zj)) {
//         stop("data must contain the variables in covariates");
//       }
//       NumericVector u = data[zj];
//       for (int i=0; i<n; ++i) {
//         zn(i,j) = u[i];
//       }
//     }
//   }
//   
//   bool has_weight = hasVariable(data, weight);
//   NumericVector weightn(n, 1.0);
//   if (has_weight) {
//     NumericVector weightnz = data[weight];
//     weightn = clone(weightnz);
//     if (is_true(any(weightn <= 0))) {
//       stop("weight must be greater than 0");
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
//       stop("incorrect type for the id variable in the input data");
//     }
//   }
//   
//   std::string meth = ties;
//   std::for_each(meth.begin(), meth.end(), [](char & c) {
//     c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
//   });
//   
//   int method = meth == "efron" ? 1 : 0;
//   
//   // unify right censored data with counting process data
//   NumericVector tstartn(n), tstopn(n);
//   if (!has_time2) {
//     tstopn = timen;
//   } else {
//     tstartn = timen;
//     tstopn = time2n;
//   }
//   
//   // exclude observations with missing values
//   LogicalVector sub(n,1);
//   for (int i=0; i<n; ++i) {
//     if (stratumn[i] == NA_INTEGER || idn[i] == NA_INTEGER ||
//         std::isnan(tstartn[i]) || std::isnan(tstopn[i]) ||
//         eventn[i] == NA_INTEGER || std::isnan(weightn[i]) ||
//         std::isnan(offsetn[i])) {
//       sub[i] = 0;
//     }
//     for (int j=0; j<p; ++j) {
//       if (std::isnan(zn(i,j))) sub[i] = 0;
//     }
//   }
//   
//   IntegerVector order = which(sub);
//   stratumn = stratumn[order];
//   tstartn = tstartn[order];
//   tstopn = tstopn[order];
//   eventn = eventn[order];
//   weightn = weightn[order];
//   offsetn = offsetn[order];
//   idn = idn[order];
//   if (p > 0) zn = subset_matrix_by_row(zn, order);
//   n = sum(sub);
//   if (n == 0) stop("no observations left after removing missing values");
//   
//   // sort by stratum
//   IntegerVector order0 = seq(0, n-1);
//   std::sort(order0.begin(), order0.end(), [&](int i, int j) {
//     return stratumn[i] < stratumn[j];
//   });
//   
//   IntegerVector stratum1z = stratumn[order0];
//   NumericVector tstart1z = tstartn[order0];
//   NumericVector tstop1z = tstopn[order0];
//   IntegerVector event1z = eventn[order0];
//   
//   // locate the first observation within each stratum
//   IntegerVector istratum(1,0);
//   for (int i=1; i<n; ++i) {
//     if (stratum1z[i] != stratum1z[i-1]) {
//       istratum.push_back(i);
//     }
//   }
//   
//   int nstrata = static_cast<int>(istratum.size());
//   istratum.push_back(n);
//   
//   // ignore subjects not at risk for any event time
//   IntegerVector ignore1z(n);
//   for (int i=0; i<nstrata; ++i) {
//     IntegerVector q0 = Range(istratum[i], istratum[i+1]-1);
//     NumericVector tstart0 = tstart1z[q0];
//     NumericVector tstop0 = tstop1z[q0];
//     IntegerVector event0 = event1z[q0];
//     NumericVector etime = tstop0[event0==1];
//     etime = unique(etime);
//     etime.sort();
//     IntegerVector index1 = findInterval3(tstart0, etime, 0, 0, 0);
//     IntegerVector index2 = findInterval3(tstop0, etime, 0, 0, 0);
//     for (int j=istratum[i]; j<istratum[i+1]; ++j) {
//       int j0 = j - istratum[i];
//       if (index1[j0] == index2[j0]) { // no event in (tstart, tstop]
//         ignore1z[j] = 1;
//       } else {
//         ignore1z[j] = 0;
//       }
//     }
//   }
//   
//   IntegerVector ignore(n);
//   for (int i=0; i<n; ++i) {
//     ignore[order0[i]] = ignore1z[i];
//   }
//   
//   int nused = n - sum(ignore);
//   
//   order = seq(0, n-1);
//   IntegerVector idx(1,0);
//   int nids = n;
//   if (has_id) { // collapse over id
//     std::sort(order.begin(), order.end(), [&](int i, int j) {
//       return idn[i] < idn[j];
//     });
//     
//     IntegerVector id2 = idn[order];
//     for (int i=1; i<n; ++i) {
//       if (id2[i] != id2[i-1]) {
//         idx.push_back(i);
//       }
//     }
//     
//     nids = static_cast<int>(idx.size());
//     idx.push_back(n);
//   }
//   
//   
//   List result;
//   if (type == "martingale") {
//     NumericVector rr = clone(resmart);
//     if (weighted) rr = rr*weightn;
//     
//     if (collapse) { // collapse over id
//       NumericVector rr2(nids);
//       for (int i=0; i<nids; ++i) {
//         for (int j=idx[i]; j<idx[i+1]; ++j) {
//           rr2[i] += rr[order[j]];
//         }
//       }
//       
//       rr = rr2;
//     }
//     
//     result = List::create(Named("resid") = rr);
//   } else if (type == "deviance") {
//     NumericVector rr = clone(resmart);
//     IntegerVector status = clone(eventn);
//     int m = n;
//     
//     if (weighted) rr = rr*weightn;
//     
//     if (collapse) { // collapse over id
//       NumericVector rr2(nids);
//       IntegerVector status2(nids);
//       for (int i=0; i<nids; ++i) {
//         for (int j=idx[i]; j<idx[i+1]; ++j) {
//           int k = order[j];
//           rr2[i] += rr[k];
//           status2[i] += eventn[k];
//         }
//       }
//       
//       rr = rr2;
//       status = status2;
//       m = nids;
//     }
//     
//     for (int i=0; i<m; ++i) {
//       double temp = status[i] == 0 ? 0 : status[i]*std::log(status[i] - rr[i]);
//       rr[i] = ((rr[i]>0) - (rr[i]<0))*std::sqrt(-2*(rr[i] + temp));
//     }
//     
//     result = List::create(Named("resid") = rr);
//   } else if (p == 0) {
//     stop("covariates must be present for score and schoenfeld residuals");
//   } else {
//     NumericMatrix rr(n,p);
//     if (type == "score" || type == "dfbeta" || type == "dfbetas") {
//       // sort by stopping time in descending order within each stratum
//       IntegerVector order2 = seq(0, n-1);
//       std::sort(order2.begin(), order2.end(), [&](int i, int j) {
//         if (ignore[i] != ignore[j]) return ignore[i] < ignore[j];
//         if (stratumn[i] != stratumn[j]) return stratumn[i] < stratumn[j];
//         if (tstopn[i] != tstopn[j]) return tstopn[i] > tstopn[j];
//         return eventn[i] < eventn[j];
//       });
//       
//       IntegerVector stratum1 = stratumn[order2];
//       NumericVector tstart1 = tstartn[order2];
//       NumericVector tstop1 = tstopn[order2];
//       IntegerVector event1 = eventn[order2];
//       NumericVector weight1 = weightn[order2];
//       NumericVector offset1 = offsetn[order2];
//       IntegerVector ignore1 = ignore[order2];
//       NumericMatrix z1 = subset_matrix_by_row(zn, order2);
//       
//       // sort by starting time in descending order within each stratum
//       IntegerVector order1 = seq(0, n-1);
//       std::sort(order1.begin(), order1.end(), [&](int i, int j) {
//         if (ignore1[i] != ignore1[j]) return ignore1[i] < ignore1[j];
//         if (stratum1[i] != stratum1[j]) return stratum1[i] < stratum1[j];
//         return tstart1[i] > tstart1[j];
//       });
//       
//       coxparams param = {nused, stratum1, tstart1, tstop1, event1,
//                          weight1, offset1, z1, order1, method};
//       
//       NumericMatrix ressco = f_ressco_2(p, beta, &param);
//       
//       NumericMatrix score(n,p);
//       for (int i=0; i<n; ++i) {
//         score(order2[i],_) = ressco(i,_); // original order
//       }
//       
//       if (type == "dfbeta" || type == "dfbetas") {
//         for (int i=0; i<n; ++i) {
//           for (int k=0; k<p; ++k) {
//             for (int j=0; j<p; ++j) {
//               rr(i,k) += score(i,j)*vbeta(j,k);
//             }
//             if (type == "dfbetas") {
//               rr(i,k) /= std::sqrt(vbeta(k,k));
//             }
//           }
//         }
//       } else {
//         rr = score;
//       }
//       
//       if (weighted) {
//         for (int i=0; i<n; ++i) {
//           for (int k=0; k<p; ++k) {
//             rr(i,k) = rr(i,k)*weightn[i];
//           }
//         }
//       }
//       
//       if (collapse) { // collapse over id
//         NumericMatrix rr2(nids,p);
//         for (int i=0; i<nids; ++i) {
//           for (int k=0; k<p; ++k) {
//             for (int j=idx[i]; j<idx[i+1]; ++j) {
//               rr2(i,k) += rr(order[j],k);
//             }
//           }
//         }
//         
//         rr = rr2;
//       }
//       
//       result = List::create(Named("resid") = rr);
//     } else if (type == "schoenfeld" || type == "scaledsch") {
//       // sort by stopping time in descending order within each stratum
//       IntegerVector order2 = seq(0, n-1);
//       std::sort(order2.begin(), order2.end(), [&](int i, int j) {
//         if (ignore[i] != ignore[j]) return ignore[i] < ignore[j];
//         if (stratumn[i] != stratumn[j]) return stratumn[i] > stratumn[j];
//         if (tstopn[i] != tstopn[j]) return tstopn[i] > tstopn[j];
//         if (eventn[i] != eventn[j]) return eventn[i] < eventn[j];
//         return idn[i] > idn[j];
//       });
//       
//       IntegerVector stratum1 = stratumn[order2];
//       NumericVector tstart1 = tstartn[order2];
//       NumericVector tstop1 = tstopn[order2];
//       IntegerVector event1 = eventn[order2];
//       NumericVector weight1 = weightn[order2];
//       NumericVector offset1 = offsetn[order2];
//       IntegerVector id1 = idn[order2];
//       IntegerVector ignore1 = ignore[order2];
//       NumericMatrix z1 = subset_matrix_by_row(zn, order2);
//       
//       // sort by starting time in descending order within each stratum
//       IntegerVector order1 = seq(0, n-1);
//       std::sort(order1.begin(), order1.end(), [&](int i, int j) {
//         if (ignore1[i] != ignore1[j]) return ignore1[i] < ignore1[j];
//         if (stratum1[i] != stratum1[j]) return stratum1[i] > stratum1[j];
//         return tstart1[i] > tstart1[j];
//       });
//       
//       coxparams param = {nused, stratum1, tstart1, tstop1, event1,
//                          weight1, offset1, z1, order1, method};
//       
//       List out = f_ressch(p, beta, &param);
//       rr = as<NumericMatrix>(out["resid"]);
//       IntegerVector index = out["index"];
//       IntegerVector stratum2 = stratum1[index];
//       NumericVector time3 = tstop1[index];
//       int ndead = static_cast<int>(index.size());
//       
//       if (weighted) {
//         for (int i=0; i<ndead; ++i) {
//           for (int k=0; k<p; ++k) {
//             rr(i,k) = rr(i,k)*weightn[i];
//           }
//         }
//       }
//       
//       if (type == "scaledsch") {
//         NumericMatrix rr2(ndead,p);
//         for (int i=0; i<ndead; ++i) {
//           for (int k=0; k<p; ++k) {
//             for (int j=0; j<p; ++j) {
//               rr2(i,k) += rr(i,j)*vbeta(j,k);
//             }
//             rr2(i,k) = rr2(i,k)*ndead + beta[k];
//           }
//         }
//         rr = rr2;
//       }
//       
//       result = List::create(
//         Named("resid") = rr,
//         Named("time") = time3);
//       
//       IntegerVector stratum3 = unique(stratum2);
//       if (stratum3.size() > 1) {
//         List strata(p_stratum);
//         for (int i = 0; i < p_stratum; ++i) {
//           std::string s = as<std::string>(stratum[i]);
//           SEXP col = u_stratum[s];
//           SEXPTYPE col_type = TYPEOF(col);
//           if (col_type == INTSXP) {
//             IntegerVector v = col;
//             strata[i] = v[stratum2];
//           } else if (col_type == REALSXP) {
//             NumericVector v = col;
//             strata[i] = v[stratum2];
//           } else if (col_type == STRSXP) {
//             StringVector v = col;
//             strata[i] = v[stratum2];
//           } else {
//             stop("Unsupported type for stratum variable: " + s);
//           }
//         }
//         strata.attr("names") = stratum;
//         result.push_back(as<DataFrame>(strata), "strata");
//       }
//     } else {
//       stop("unknown type of residuals");
//     }
//   }
//   
//   return result;
// }
// 
// 
// // function for individual contributions to score and information matrix
// // for the Cox model
// List f_der_i_2(int p, const NumericVector& par, void* ex) {
//   coxparams* param = (coxparams*) ex;
//   
//   const int nused = param->nused;
//   const int method = param->method;
//   
//   // Precompute eta and exp(eta)
//   NumericVector eta(nused);
//   NumericVector exp_eta(nused);
//   for (int person = 0; person < nused; ++person) {
//     double val = param->offset[person];
//     for (int i = 0; i < p; ++i) {
//       val += par[i] * param->z(person,i);
//     }
//     eta[person] = val;
//     exp_eta[person] = std::exp(val);
//   }
//   
//   NumericMatrix u(nused,p);         // score vector for each individual
//   NumericMatrix imat(nused,p*p);    // information matrix for each individual
//   NumericVector a(p);               // s1(beta,k,t)
//   NumericVector a2(p);              // sum of w*exp(zbeta)*z for the deaths
//   NumericMatrix cmat(p,p);          // s2(beta,k,t)
//   NumericMatrix cmat2(p,p);         // sum of w*exp(zbeta)*z*z' for the deaths
//   double denom = 0.0;               // s0(beta,k,t)
//   double denom2 = 0.0;              // sum of weighted risks for deaths
//   int ndead = 0;                    // number of deaths at this time point
//   
//   
//   int istrata = param->strata[0];
//   int i1 = 0; // index for removing out-of-risk subjects
//   
//   // Loop through subjects
//   for (int person = 0; person < nused; ) {
//     // Reset when entering a new stratum
//     if (param->strata[person] != istrata) {
//       istrata = param->strata[person];
//       i1 = person;
//       denom = 0.0;
//       
//       for (int i = 0; i < p; ++i) {
//         a[i] = 0.0;
//         for (int j = 0; j <= i; ++j) {
//           cmat(i,j) = 0.0;
//         }
//       }
//     }
//     
//     const double dtime = param->tstop[person];
//     
//     // Process all persons tied at this dtime
//     int person1 = person;
//     
//     for (; person < nused && param->tstop[person] == dtime &&
//          param->strata[person] == istrata; ++person) {
//       
//       const double w = param->weight[person];
//       const double r = w * exp_eta[person];
//       
//       if (param->event[person] == 0) {
//         denom += r;
//         
//         for (int i = 0; i < p; ++i) {
//           const double zi = param->z(person,i);
//           a[i] += r * zi;
//           for (int j = 0; j <= i; ++j) {
//             const double zj = param->z(person,j);
//             cmat(i,j) += r * zi * zj;
//           }
//         }
//       } else {
//         ++ndead;
//         denom2 += r;
//         for (int i = 0; i < p; ++i) {
//           const double zi = param->z(person,i);
//           a2[i] += r * zi;
//           u(person,i) = w * zi;
//           for (int j = 0; j <= i; ++j) {
//             const double zj = param->z(person,j);
//             cmat2(i,j) += r * zi * zj;
//           }
//         }
//       }
//     }
//     
//     // Remove subjects leaving risk set
//     for (; i1 < nused; ++i1) {
//       const int p1 = param->order1[i1];
//       if (param->tstart[p1] < dtime || param->strata[p1] != istrata) break;
//       
//       const double r = param->weight[p1] * exp_eta[p1];
//       denom -= r;
//       for (int i = 0; i < p; ++i) {
//         const double zi = param->z(p1,i);
//         a[i] -= r * zi;
//         for (int j = 0; j <= i; ++j) {
//           const double zj = param->z(p1,j);
//           cmat(i,j) -= r * zi * zj;
//         }
//       }
//     }
//     
//     // Add contributions for deaths at this time
//     if (ndead > 0) {
//       if (method == 0 || ndead == 1) { // Breslow or single event
//         denom += denom2;
//         for (int i = 0; i < p; ++i) {
//           a[i] += a2[i];
//           const double xbar = a[i] / denom;
//           for (int human = person1; human < person; ++human) {
//             if (param->event[human] == 1)
//               u(human,i) -= param->weight[human] * xbar;
//           }
//           
//           for (int j = 0; j <= i; ++j) {
//             cmat(i,j) += cmat2(i,j);
//             
//             for (int human = person1; human < person; ++human) {
//               if (param->event[human] == 1)
//                 imat(human, i*p+j) = param->weight[human] *
//                   (cmat(i,j) - xbar * a[j]) / denom;
//             }
//           }
//         }
//       } else { // Efron method
//         const double increment = denom2 / ndead;
//         for (int l = 0; l < ndead; ++l) {
//           denom += increment;
//           for (int i = 0; i < p; ++i) {
//             a[i] += a2[i] / ndead;
//             const double xbar = a[i] / denom;
//             for (int human = person1; human < person; ++human) {
//               if (param->event[human] == 1)
//                 u(human,i) -= param->weight[human]/ndead * xbar;
//             }
//             
//             for (int j = 0; j <= i; ++j) {
//               cmat(i,j) += cmat2(i,j) / ndead;
//               for (int human = person1; human < person; ++human) {
//                 if (param->event[human] == 1)
//                   imat(human,i*p+j) +=  param->weight[human]/ndead *
//                     (cmat(i,j) - xbar * a[j]) / denom;
//               }
//             }
//           }
//         }
//       }
//       
//       // Reset after processing deaths
//       ndead = 0;
//       denom2 = 0.0;
//       for (int i = 0; i < p; ++i) {
//         a2[i] = 0;
//         for (int j = 0; j <= i; ++j) {
//           cmat2(i,j) = 0;
//         }
//       }
//     }
//   }
//   
//   // fill the symmetric elements of the information matrix
//   for (int person = 0; person < nused; ++person)
//     for (int i = 0; i < p - 1; ++i)
//       for (int j = i+1; j < p; ++j)
//         imat(person,i*p+j) = imat(person,j*p+i);
//   
//   List result = List::create(
//     _["score_i"] = u,
//     _["imat_i"] = imat
//   );
//   
//   return result;
// }
// 
// 
// // [[Rcpp::export]]
// List assess_phregcpp(const int p,
//                      const NumericVector& beta,
//                      const NumericMatrix& vbeta,
//                      DataFrame data,
//                      const StringVector& stratum = "",
//                      const std::string time = "time",
//                      const std::string time2 = "",
//                      const std::string event = "event",
//                      const StringVector& covariates = "",
//                      const std::string weight = "",
//                      const std::string offset = "",
//                      const std::string ties = "efron",
//                      const int resample = 1000,
//                      const std::uint32_t seed) {
//   
//   td::mt19937_64 rng(seed);                      // choose generator and seed
//   std::normal_distribution<double> dist(0.0, 1.0);
//   
//   if (seed != NA_INTEGER) set_seed(seed);
//   
//   if (p <= 0) {
//     stop("covariates must be present to test proportional hazards");
//   }
//   
//   int n = data.nrows();
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
//   bool has_time = hasVariable(data, time);
//   if (!has_time) stop("data must contain the time variable");
//   NumericVector timenz = data[time];
//   NumericVector timen = clone(timenz);
//   if (is_true(any(timen < 0))) {
//     stop("time must be nonnegative for each observation");
//   }
//   
//   bool has_time2 = hasVariable(data, time2);
//   NumericVector time2n(n);
//   if (has_time2) {
//     NumericVector time2nz = data[time2];
//     time2n = clone(time2nz);
//     if (is_true(any(time2n <= timen))) {
//       stop("time2 must be greater than time for each observation");
//     }
//   }
//   
//   bool has_event = hasVariable(data, event);
//   if (!has_event) stop("data must contain the event variable");
//   IntegerVector eventnz = data[event];
//   IntegerVector eventn = clone(eventnz);
//   if (is_true(any((eventn != 1) & (eventn != 0)))) {
//     stop("event must be 1 or 0 for each observation");
//   }
//   
//   NumericMatrix zn(n,p);
//   if (p > 0) {
//     for (int j=0; j<p; ++j) {
//       String zj = covariates[j];
//       if (!hasVariable(data, zj)) {
//         stop("data must contain the variables in covariates");
//       }
//       NumericVector u = data[zj];
//       for (int i=0; i<n; ++i) {
//         zn(i,j) = u[i];
//       }
//     }
//   }
//   
//   bool has_weight = hasVariable(data, weight);
//   NumericVector weightn(n, 1.0);
//   if (has_weight) {
//     NumericVector weightnz = data[weight];
//     weightn = clone(weightnz);
//     if (is_true(any(weightn <= 0))) {
//       stop("weight must be greater than 0");
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
//   std::string meth = ties;
//   std::for_each(meth.begin(), meth.end(), [](char & c) {
//     c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
//   });
//   
//   int method = meth == "efron" ? 1 : 0;
//   
//   // unify right censored data with counting process data
//   NumericVector tstartn(n), tstopn(n);
//   if (!has_time2) {
//     tstopn = timen;
//   } else {
//     tstartn = timen;
//     tstopn = time2n;
//   }
//   
//   // exclude observations with missing values
//   LogicalVector sub(n,1);
//   for (int i=0; i<n; ++i) {
//     if (stratumn[i] == NA_INTEGER || std::isnan(tstartn[i]) ||
//         std::isnan(tstopn[i]) || eventn[i] == NA_INTEGER ||
//         std::isnan(weightn[i]) || std::isnan(offsetn[i])) {
//       sub[i] = 0;
//     }
//     for (int j=0; j<p; ++j) {
//       if (std::isnan(zn(i,j))) sub[i] = 0;
//     }
//   }
//   
//   IntegerVector order = which(sub);
//   stratumn = stratumn[order];
//   tstartn = tstartn[order];
//   tstopn = tstopn[order];
//   eventn = eventn[order];
//   weightn = weightn[order];
//   offsetn = offsetn[order];
//   zn = subset_matrix_by_row(zn, order);
//   n = sum(sub);
//   if (n == 0) stop("no observations left after removing missing values");
//   
//   // sort by stratum
//   IntegerVector order0 = seq(0, n-1);
//   std::sort(order0.begin(), order0.end(), [&](int i, int j) {
//     return stratumn[i] < stratumn[j];
//   });
//   
//   IntegerVector stratum1z = stratumn[order0];
//   NumericVector tstart1z = tstartn[order0];
//   NumericVector tstop1z = tstopn[order0];
//   IntegerVector event1z = eventn[order0];
//   
//   // locate the first observation within each stratum
//   IntegerVector istratum(1,0);
//   for (int i=1; i<n; ++i) {
//     if (stratum1z[i] != stratum1z[i-1]) {
//       istratum.push_back(i);
//     }
//   }
//   
//   int nstrata = static_cast<int>(istratum.size());
//   istratum.push_back(n);
//   
//   // ignore subjects not at risk for any event time
//   IntegerVector ignore1z(n);
//   for (int i=0; i<nstrata; ++i) {
//     IntegerVector q0 = Range(istratum[i], istratum[i+1]-1);
//     NumericVector tstart0 = tstart1z[q0];
//     NumericVector tstop0 = tstop1z[q0];
//     IntegerVector event0 = event1z[q0];
//     NumericVector etime = tstop0[event0==1];
//     etime = unique(etime);
//     etime.sort();
//     IntegerVector index1 = findInterval3(tstart0, etime, 0, 0, 0);
//     IntegerVector index2 = findInterval3(tstop0, etime, 0, 0, 0);
//     for (int j=istratum[i]; j<istratum[i+1]; ++j) {
//       int j0 = j - istratum[i];
//       if (index1[j0] == index2[j0]) { // no event in (tstart, tstop]
//         ignore1z[j] = 1;
//       } else {
//         ignore1z[j] = 0;
//       }
//     }
//   }
//   
//   IntegerVector ignoren(n); // final ignore vector
//   for (int i=0; i<n; ++i) {
//     ignoren[order0[i]] = ignore1z[i];
//   }
//   
//   int nused = n - sum(ignoren);
//   
//   // sort by stopping time in descending order within each stratum
//   IntegerVector order1 = seq(0, n-1);
//   std::sort(order1.begin(), order1.end(), [&](int i, int j) {
//     if (ignoren[i] != ignoren[j]) return ignoren[i] < ignoren[j];
//     if (stratumn[i] != stratumn[j]) return stratumn[i] < stratumn[j];
//     if (tstopn[i] != tstopn[j]) return tstopn[i] > tstopn[j];
//     return eventn[i] < eventn[j];
//   });
//   
//   IntegerVector stratum1 = stratumn[order1];
//   NumericVector tstart1 = tstartn[order1];
//   NumericVector tstop1 = tstopn[order1];
//   IntegerVector event1 = eventn[order1];
//   NumericVector weight1 = weightn[order1];
//   NumericVector offset1 = offsetn[order1];
//   IntegerVector ignore1 = ignoren[order1];
//   NumericMatrix z1 = subset_matrix_by_row(zn, order1);
//   
//   // sort by starting time in descending order within each stratum
//   IntegerVector order10 = seq(0, n-1);
//   std::sort(order10.begin(), order10.end(), [&](int i, int j) {
//     if (ignore1[i] != ignore1[j]) return ignore1[i] < ignore1[j];
//     if (stratum1[i] != stratum1[j]) return stratum1[i] < stratum1[j];
//     return tstart1[i] > tstart1[j];
//   });
//   
//   coxparams param = {nused, stratum1, tstart1, tstop1, event1,
//                      weight1, offset1, z1, order10, method};
//   
//   List out = f_der_i_2(p, beta, &param);
//   NumericMatrix score_i = as<NumericMatrix>(out["score_i"]);
//   NumericMatrix imat_i = as<NumericMatrix>(out["imat_i"]);
//   
//   // order data by ascending tstop
//   IntegerVector order2 = seq(0, n-1);
//   std::sort(order2.begin(), order2.end(), [&](int i, int j) {
//     return tstop1[i] < tstop1[j];
//   });
//   
//   NumericVector tstop2 = tstop1[order2];
//   IntegerVector event2 = event1[order2];
//   NumericMatrix score_i2 = subset_matrix_by_row(score_i, order2);
//   NumericMatrix imat_i2 = subset_matrix_by_row(imat_i, order2);
//   
//   // only consider event times
//   IntegerVector q = which(event2 == 1);
//   NumericVector tstop3 = tstop2[q];
//   NumericMatrix score_i3 = subset_matrix_by_row(score_i2, q);
//   NumericMatrix imat_i3 = subset_matrix_by_row(imat_i2, q);
//   int nq = static_cast<int>(q.size());
//   
//   // identify unique event times
//   NumericVector t = unique(tstop3);
//   t.sort();
//   
//   int nt = static_cast<int>(t.size());
//   NumericMatrix score_i4(nt, p);
//   NumericMatrix imat_i4(nt, p*p);
//   
//   tstop3.push_front(-1.0); // to facilitate the loop
//   int k = -1;
//   for (int i=1; i<=nq; ++i) {
//     if (tstop3[i] != tstop3[i-1]) k++;
//     
//     for (int j=0; j<p; ++j) {
//       score_i4(k,j) += score_i3(i-1,j);
//       for (int l=0; l<p; ++l) {
//         imat_i4(k,j*p+l) += imat_i3(i-1,j*p+l);
//       }
//     }
//   }
//   
//   // cumulative sums of score and information processes
//   for (int i=1; i<nt; ++i) {
//     for (int j=0; j<p; ++j) {
//       score_i4(i,j) += score_i4(i-1,j);
//       for (int l=0; l<p; ++l) {
//         imat_i4(i,j*p+l) += imat_i4(i-1,j*p+l);
//       }
//     }
//   }
//   
//   // standardize the score processes
//   for (int i=0; i<nt; ++i) {
//     for (int j=0; j<p; ++j) {
//       score_i4(i,j) *= std::sqrt(vbeta(j,j));
//     }
//   }
//   
//   // resampling to approximate the null distribution
//   List G_list(resample);
//   for (int r=0; r<resample; ++r) {
//     NumericVector G(nq);
//     for (int i=0; i<nq; ++i) {
//       G[i] = dist(rng);
//     }
//     G_list[r] = G;
//   }
//   
//   // sum(u*G)
//   List U_list(resample);
//   for (int r=0; r<resample; ++r) {
//     NumericVector U(p);
//     NumericVector G = G_list[r];
//     for (int i=0; i<nq; ++i) {
//       for (int j=0; j<p; ++j) {
//         U[j] += score_i3(i,j)*G[i];
//       }
//     }
//     U_list[r] = U;
//   }
//   
//   // sum(I(X <= t)*U*G)
//   List score_i_list(resample);
//   for (int r=0; r<resample; ++r) {
//     NumericMatrix score_i4(nt, p);
//     NumericVector G = G_list[r];
//     
//     int k = -1;
//     for (int i=1; i<=nq; ++i) {
//       if (tstop3[i] != tstop3[i-1]) k++;
//       
//       for (int j=0; j<p; ++j) {
//         score_i4(k,j) += score_i3(i-1,j)*G[i-1];
//       }
//     }
//     
//     // cumulative sums of score processes
//     for (int i=1; i<nt; ++i) {
//       for (int j=0; j<p; ++j) {
//         score_i4(i,j) += score_i4(i-1,j);
//       }
//     }
//     
//     score_i_list[r] = score_i4;
//   }
//   
//   // piece together the result
//   List score_i_list_2(resample);
//   for (int r=0; r<resample; ++r) {
//     NumericMatrix score_i4 = score_i_list[r];
//     NumericMatrix score_i4_2(nt, p);
//     NumericVector U = U_list[r];
//     for (int i=0; i<nt; ++i) {
//       for (int j=0; j<p; ++j) {
//         score_i4_2(i,j) = score_i4(i,j);
//         for (int l=0; l<p; ++l) {
//           for (int s=0; s<p; ++s) {
//             score_i4_2(i,j) -= imat_i4(i,j*p+l)*vbeta(l,s)*U[s];
//           }
//         }
//       }
//     }
//     score_i_list_2[r] = score_i4_2;
//   }
//   
//   // standardize the resampled score processes
//   for (int r=0; r<resample; ++r) {
//     NumericMatrix score_i4 = score_i_list_2[r];
//     for (int i=0; i<nt; ++i) {
//       for (int j=0; j<p; ++j) {
//         score_i4(i,j) *= std::sqrt(vbeta(j,j));
//       }
//     }
//     score_i_list_2[r] = score_i4;
//   }
//   
//   // calculate individual p-values for proportional hazards assumption test
//   NumericVector p_values(p+1);
//   NumericVector maxabs(p+1);
//   for (int j=0; j<p; ++j) {
//     NumericVector obs(nt);
//     for (int i=0; i<nt; ++i) {
//       obs[i] = fabs(score_i4(i,j));
//     }
//     double tobs = max(obs);
//     maxabs[j] = tobs;
//     
//     int count = 0;
//     for (int r=0; r<resample; ++r) {
//       NumericMatrix score_i4 = score_i_list_2[r];
//       NumericVector sim(nt);
//       for (int i=0; i<nt; ++i) {
//         sim[i] = fabs(score_i4(i,j));
//       }
//       double tsim = max(sim);
//       if (tsim >= tobs) count++;
//     }
//     p_values[j] = static_cast<double>(count) / static_cast<double>(resample);
//   }
//   
//   // p-value for the global test
//   NumericVector obs(nt);
//   for (int i=0; i<nt; ++i) {
//     for (int j=0; j<p; ++j) {
//       obs[i] += fabs(score_i4(i,j));
//     }
//   }
//   double tobs = max(obs);
//   
//   int count = 0;
//   for (int r=0; r<resample; ++r) {
//     NumericMatrix score_i4 = score_i_list_2[r];
//     NumericVector sim(nt);
//     for (int i=0; i<nt; ++i) {
//       for (int j=0; j<p; ++j) {
//         sim[i] += fabs(score_i4(i,j));
//       }
//     }
//     double tsim = max(sim);
//     if (tsim >= tobs) count++;
//   }
//   p_values[p] = static_cast<double>(count) / static_cast<double>(resample);
//   
//   maxabs[p] = tobs;
//   List result = List::create(
//     Named("time") = t,
//     Named("score_t") = score_i4,
//     Named("score_t_list") = score_i_list_2,
//     Named("max_abs_value") = maxabs,
//     Named("p_value") = p_values
//   );
//   
//   return result;
// }
// 
// 
// // [[Rcpp::export]]
// List zph_phregcpp(const int p,
//                   const NumericVector& beta,
//                   const NumericMatrix& vbeta,
//                   const NumericVector& resmart,
//                   DataFrame data,
//                   const StringVector& stratum = "",
//                   const std::string time = "time",
//                   const std::string time2 = "",
//                   const std::string event = "event",
//                   const StringVector& covariates = "",
//                   const std::string weight = "",
//                   const std::string offset = "",
//                   const std::string ties = "efron",
//                   const std::string transform= "km") {
//   
//   if (p <= 0) {
//     stop("covariates must be present to test proportional hazards");
//   }
//   
//   int n = data.nrows();
//   
//   bool has_stratum;
//   IntegerVector stratumn(n);
//   int p_stratum = static_cast<int>(stratum.size());
//   if (p_stratum == 1 && (stratum[0] == "" || stratum[0] == "none")) {
//     has_stratum = false;
//     stratumn.fill(0);
//   } else {
//     has_stratum = true;
//     List out = bygroup(data, stratum);
//     stratumn = out["index"];
//   }
//   
//   bool has_time = hasVariable(data, time);
//   if (!has_time) stop("data must contain the time variable");
//   NumericVector timenz = data[time];
//   NumericVector timen = clone(timenz);
//   if (is_true(any(timen < 0))) {
//     stop("time must be nonnegative for each observation");
//   }
//   
//   bool has_time2 = hasVariable(data, time2);
//   NumericVector time2n(n);
//   if (has_time2) {
//     NumericVector time2nz = data[time2];
//     time2n = clone(time2nz);
//     if (is_true(any(time2n <= timen))) {
//       stop("time2 must be greater than time for each observation");
//     }
//   }
//   
//   bool has_event = hasVariable(data, event);
//   if (!has_event) stop("data must contain the event variable");
//   IntegerVector eventnz = data[event];
//   IntegerVector eventn = clone(eventnz);
//   if (is_true(any((eventn != 1) & (eventn != 0)))) {
//     stop("event must be 1 or 0 for each observation");
//   }
//   
//   NumericMatrix zn(n,p);
//   if (p > 0) {
//     for (int j=0; j<p; ++j) {
//       String zj = covariates[j];
//       if (!hasVariable(data, zj)) {
//         stop("data must contain the variables in covariates");
//       }
//       NumericVector u = data[zj];
//       for (int i=0; i<n; ++i) {
//         zn(i,j) = u[i];
//       }
//     }
//   }
//   
//   bool has_weight = hasVariable(data, weight);
//   NumericVector weightn(n, 1.0);
//   if (has_weight) {
//     NumericVector weightnz = data[weight];
//     weightn = clone(weightnz);
//     if (is_true(any(weightn <= 0))) {
//       stop("weight must be greater than 0");
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
//   std::string meth = ties;
//   std::for_each(meth.begin(), meth.end(), [](char & c) {
//     c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
//   });
//   
//   int method = meth == "efron" ? 1 : 0;
//   
//   // unify right censored data with counting process data
//   NumericVector tstartn(n), tstopn(n);
//   if (!has_time2) {
//     tstopn = timen;
//   } else {
//     tstartn = timen;
//     tstopn = time2n;
//   }
//   
//   // exclude observations with missing values
//   LogicalVector sub(n,1);
//   for (int i=0; i<n; ++i) {
//     if (stratumn[i] == NA_INTEGER || std::isnan(tstartn[i]) ||
//         std::isnan(tstopn[i]) || eventn[i] == NA_INTEGER ||
//         std::isnan(weightn[i]) || std::isnan(offsetn[i])) {
//       sub[i] = 0;
//     }
//     for (int j=0; j<p; ++j) {
//       if (std::isnan(zn(i,j))) sub[i] = 0;
//     }
//   }
//   
//   IntegerVector order = which(sub);
//   stratumn = stratumn[order];
//   tstartn = tstartn[order];
//   tstopn = tstopn[order];
//   eventn = eventn[order];
//   weightn = weightn[order];
//   offsetn = offsetn[order];
//   zn = subset_matrix_by_row(zn, order);
//   n = sum(sub);
//   if (n == 0) stop("no observations left after removing missing values");
//   
//   // sort by stratum
//   IntegerVector order0 = seq(0, n-1);
//   std::sort(order0.begin(), order0.end(), [&](int i, int j) {
//     return stratumn[i] < stratumn[j];
//   });
//   
//   IntegerVector stratum1z = stratumn[order0];
//   NumericVector tstart1z = tstartn[order0];
//   NumericVector tstop1z = tstopn[order0];
//   IntegerVector event1z = eventn[order0];
//   
//   // locate the first observation within each stratum
//   IntegerVector istratum(1,0);
//   for (int i=1; i<n; ++i) {
//     if (stratum1z[i] != stratum1z[i-1]) {
//       istratum.push_back(i);
//     }
//   }
//   
//   int nstrata = static_cast<int>(istratum.size());
//   istratum.push_back(n);
//   
//   // ignore subjects not at risk for any event time
//   IntegerVector ignore1z(n);
//   for (int i=0; i<nstrata; ++i) {
//     IntegerVector q0 = Range(istratum[i], istratum[i+1]-1);
//     NumericVector tstart0 = tstart1z[q0];
//     NumericVector tstop0 = tstop1z[q0];
//     IntegerVector event0 = event1z[q0];
//     NumericVector etime = tstop0[event0==1];
//     etime = unique(etime);
//     etime.sort();
//     IntegerVector index1 = findInterval3(tstart0, etime, 0, 0, 0);
//     IntegerVector index2 = findInterval3(tstop0, etime, 0, 0, 0);
//     for (int j=istratum[i]; j<istratum[i+1]; ++j) {
//       int j0 = j - istratum[i];
//       if (index1[j0] == index2[j0]) { // no event in (tstart, tstop]
//         ignore1z[j] = 1;
//       } else {
//         ignore1z[j] = 0;
//       }
//     }
//   }
//   
//   IntegerVector ignoren(n); // final ignore vector
//   for (int i=0; i<n; ++i) {
//     ignoren[order0[i]] = ignore1z[i];
//   }
//   
//   int nused = n - sum(ignoren);
//   
//   // sort by stopping time in descending order within each stratum
//   IntegerVector order1 = seq(0, n-1);
//   std::sort(order1.begin(), order1.end(), [&](int i, int j) {
//     if (ignoren[i] != ignoren[j]) return ignoren[i] < ignoren[j];
//     if (stratumn[i] != stratumn[j]) return stratumn[i] < stratumn[j];
//     if (tstopn[i] != tstopn[j]) return tstopn[i] > tstopn[j];
//     return eventn[i] < eventn[j];
//   });
//   
//   IntegerVector stratum1 = stratumn[order1];
//   NumericVector tstart1 = tstartn[order1];
//   NumericVector tstop1 = tstopn[order1];
//   IntegerVector event1 = eventn[order1];
//   NumericVector weight1 = weightn[order1];
//   NumericVector offset1 = offsetn[order1];
//   IntegerVector ignore1 = ignoren[order1];
//   NumericMatrix z1 = subset_matrix_by_row(zn, order1);
//   
//   // sort by starting time in descending order within each stratum
//   IntegerVector order10 = seq(0, n-1);
//   std::sort(order10.begin(), order10.end(), [&](int i, int j) {
//     if (ignore1[i] != ignore1[j]) return ignore1[i] < ignore1[j];
//     if (stratum1[i] != stratum1[j]) return stratum1[i] < stratum1[j];
//     return tstart1[i] > tstart1[j];
//   });
//   
//   coxparams param = {nused, stratum1, tstart1, tstop1, event1,
//                      weight1, offset1, z1, order10, method};
//   
//   List out = f_der_i_2(p, beta, &param);
//   NumericMatrix score_i = as<NumericMatrix>(out["score_i"]);
//   NumericMatrix imat_i = as<NumericMatrix>(out["imat_i"]);
//   
//   // order data by ascending tstop
//   IntegerVector order2 = seq(0, n-1);
//   std::sort(order2.begin(), order2.end(), [&](int i, int j) {
//     return tstop1[i] < tstop1[j];
//   });
//   
//   IntegerVector stratum2 = stratum1[order2];
//   NumericVector tstop2 = tstop1[order2];
//   IntegerVector event2 = event1[order2];
//   NumericMatrix score_i2 = subset_matrix_by_row(score_i, order2);
//   NumericMatrix imat_i2 = subset_matrix_by_row(imat_i, order2);
//   
//   // only consider event times
//   IntegerVector q = which(event2 == 1);
//   IntegerVector stratum3 = stratum2[q];
//   NumericVector tstop3 = tstop2[q];
//   NumericMatrix score_i3 = subset_matrix_by_row(score_i2, q);
//   NumericMatrix imat_i3 = subset_matrix_by_row(imat_i2, q);
//   int nq = static_cast<int>(q.size());
//   
//   // identify unique event times
//   NumericVector t = unique(tstop3);
//   t.sort();
//   
//   int nt = static_cast<int>(t.size());
//   NumericMatrix score_i4(nt, p);
//   NumericMatrix imat_i4(nt, p*p);
//   
//   NumericVector tstop3a = clone(tstop3);
//   tstop3a.push_front(-1.0); // to facilitate the loop
//   int k = -1;
//   for (int i=1; i<=nq; ++i) {
//     if (tstop3a[i] != tstop3a[i-1]) k++;
//     
//     for (int j=0; j<p; ++j) {
//       score_i4(k,j) += score_i3(i-1,j);
//       for (int l=0; l<p; ++l) {
//         imat_i4(k,j*p+l) += imat_i3(i-1,j*p+l);
//       }
//     }
//   }
//   
//   // transformed time points
//   NumericVector g(nt);
//   if (transform == "identity") {
//     for (int i=0; i<nt; ++i) {
//       g[i] = t[i];
//     }
//   } else if (transform == "log") {
//     for (int i=0; i<nt; ++i) {
//       g[i] = std::log(t[i]);
//     }
//   } else if (transform == "rank") {
//     for (int i=0; i<nt; ++i) {
//       g[i] = static_cast<double>(i+1);
//     }
//   } else if (transform == "km") {
//     DataFrame temp = DataFrame::create(
//       Named("tstart") = tstart1,
//       Named("tstop") = tstop1,
//       Named("event") = event1
//     );
//     
//     DataFrame df_km = kmest(temp, "", "", "tstart", "tstop", "event",
//                             "", "none", 0.95, false);
//     
//     NumericVector surv = df_km["surv"];
//     surv.push_front(1.0); // for left-continuous step function
//     for (int i=0; i<nt; ++i) {
//       g[i] = 1.0 - surv[i];
//     }
//   } else {
//     Function f(transform);
//     for (int i=0; i<nt; ++i) {
//       g[i] = as<double>(f(t[i]));
//     }
//   }
//   
//   // score for theta
//   NumericVector u_theta(p);
//   for (int i=0; i<nt; ++i) {
//     for (int j=0; j<p; ++j) {
//       u_theta[j] += score_i4(i,j)*g[i];
//     }
//   }
//   
//   // covariance between u_theta and u_beta
//   NumericMatrix imat_theta_beta(p, p);
//   for (int i=0; i<nt; ++i) {
//     for (int j=0; j<p; ++j) {
//       for (int l=0; l<p; ++l) {
//         imat_theta_beta(j,l) += imat_i4(i,j*p+l)*g[i];
//       }
//     }
//   }
//   
//   // covariance for u_theta
//   NumericMatrix imat_theta(p, p);
//   for (int i=0; i<nt; ++i) {
//     for (int j=0; j<p; ++j) {
//       for (int l=0; l<p; ++l) {
//         imat_theta(j,l) += imat_i4(i,j*p+l)*g[i]*g[i];
//       }
//     }
//   }
//   
//   // conditional variance for u_theta given u_beta
//   NumericMatrix vtheta(p, p);
//   for (int j=0; j<p; ++j) {
//     for (int l=0; l<p; ++l) {
//       vtheta(j,l) = imat_theta(j,l);
//       for (int s=0; s<p; ++s) {
//         for (int r=0; r<p; ++r) {
//           vtheta(j,l) -= imat_theta_beta(j,s)*vbeta(s,r)*imat_theta_beta(l,r);
//         }
//       }
//     }
//   }
//   
//   // individual score test for theta = 0
//   NumericVector score_test(p);
//   for (int j=0; j<p; ++j) {
//     score_test[j] = u_theta[j]*u_theta[j]/vtheta(j,j);
//   }
//   
//   // global score test for theta = 0
//   double score_test_global = 0.0;
//   NumericMatrix vtheta_inv = invsympd(vtheta, p, 1e-12);
//   for (int j=0; j<p; ++j) {
//     for (int l=0; l<p; ++l) {
//       score_test_global += u_theta[j]*vtheta_inv(j,l)*u_theta[l];
//     }
//   }
//   
//   score_test.push_back(score_test_global);
//   NumericVector df(p+1,1.0);
//   df[p] = p; // degrees of freedom for the global test
//   
//   // p-values
//   NumericVector p_values(p+1);
//   for (int j=0; j<=p; ++j) {
//     p_values[j] = 1.0 - boost_pchisq(score_test[j], df[j], true, false);
//   }
//   
//   // assemble the table
//   NumericMatrix table(p+1,3);
//   for (int j=0; j<=p; ++j) {
//     table(j,0) = score_test[j];
//     table(j,1) = df[j];
//     table(j,2) = p_values[j];
//   }
//   
//   StringVector covariates_vec(p+1);
//   for (int j=0; j<p; ++j) covariates_vec[j] = covariates[j];
//   covariates_vec[p] = "GLOBAL";
//   rownames(table) = covariates_vec;
//   colnames(table) = StringVector::create("chisq", "df", "p");
//   
//   // obtain scaled schoenfeld residuals
//   List resid_list = residuals_phregcpp(
//     p, beta, vbeta, resmart, data, stratum, time, time2, event,
//     covariates, weight, offset, "", ties, "scaledsch", false, true);
//   
//   NumericMatrix sresid = as<NumericMatrix>(resid_list["resid"]);
//   colnames(sresid) = covariates;
//   
//   // repeat the g values according to the number of events
//   NumericVector g_rep(nq);
//   int idx = 0;
//   for (int i=0; i<nt; ++i) {
//     int count = 0;
//     for (int j=0; j<nq; ++j) {
//       if (tstop3[j] == t[i]) {
//         count++;
//       }
//     }
//     for (int j=0; j<count; ++j) {
//       g_rep[idx] = g[i];
//       idx++;
//     }
//   }
//   
//   List result = List::create(
//     Named("table") = table,
//     Named("x") = g_rep,
//     Named("time") = tstop3,
//     Named("y") = sresid,
//     Named("var") = vbeta*nq,
//     Named("transform") = transform
//   );
//   
//   if (has_stratum) {
//     result.push_back(stratum3, "strata");
//   }
//   
//   return result;
// }