#include "utilities.h"
#include "survival_analysis.h"

using namespace Rcpp;


//' @title Simulate Survival Data for Two-Stage Estimation Method Using 
//' g-estimation
//' @description Obtains the simulated data for baseline prognosis, 
//' disease progression, treatment switching, death, and 
//' time-dependent covariates.
//'
//' @param n The total sample size for two treatment arms combined.
//' @param allocation1 The number of subjects in the active treatment group 
//'   in a randomization block.
//' @param allocation2 The number of subjects in the control group in
//'   a randomization block.
//' @param pbprog The probability of having poor prognosis at baseline.
//' @param trtlghr The treatment effect in terms of log hazard ratio.
//' @param bprogsl The poor prognosis effect in terms of log hazard ratio.
//' @param shape1 The shape parameter for the Weibull event distribution 
//'   for the first component.
//' @param scale1 The scale parameter for the Weibull event distribution 
//'   for the first component.
//' @param shape2 The shape parameter for the Weibull event distribution 
//'   for the second component.
//' @param scale2 The scale parameter for the Weibull event distribution 
//'   for the second component.
//' @param pmix The mixing probability of the first component Weibull 
//'   distribution.
//' @param admin The administrative censoring time.
//' @param pcatnotrtbprog The probability of developing metastatic disease
//'   on control treatment with poor baseline prognosis.
//' @param pcattrtbprog The probability of developing metastatic disease
//'   on active treatment with poor baseline prognosis.
//' @param pcatnotrt The probability of developing metastatic disease
//'   on control treatment with good baseline prognosis.
//' @param pcattrt The probability of developing metastatic disease
//'   on active treatment with good baseline prognosis.
//' @param catmult The impact of metastatic disease on shortening remaining 
//'   survival time.
//' @param tdxo Whether treatment crossover depends on time-dependent 
//'   covariates between disease progression and treatment switching.
//' @param ppoor The probability of switching for poor baseline prognosis
//'   with no metastatic disease.
//' @param pgood The probability of switching for good baseline prognosis
//'   with no metastatic disease.
//' @param ppoormet The probability of switching for poor baseline prognosis
//'   after developing metastatic disease.
//' @param pgoodmet The probability of switching for good baseline prognosis
//'   after developing metastatic disease.
//' @param xomult The direct effect of crossover on extending remaining 
//'   survival time.
//' @param milestone The milestone to calculate restricted mean survival 
//'   time.
//' @param swtrt_control_only Whether treatment switching occurred only in
//'   the control group.
//' @param outputRawDataset Whether to output the raw data set
//' @param seed The seed to reproduce the simulation results.
//'   The seed from the environment will be used if left unspecified.
//'
//' @return A list with two data frames
//' 
//' * \code{sumdata}: A data frame with the following variables:
//'
//'     - \code{simtrueconstmean}: The true control group restricted mean 
//'       survival time (RMST).
//'     
//'     - \code{simtrueconstlb}: The lower bound for control group RMST.
//'     
//'     - \code{simtrueconstub}: The upper bound for control group RMST.
//'     
//'     - \code{simtrueconstse}: The standard error for control group RMST.
//'     
//'     - \code{simtrueexpstmean}: The true experimental group restricted 
//'       mean survival time (RMST).
//'     
//'     - \code{simtrueexpstlb}: The lower bound for experimental group RMST.
//'     
//'     - \code{simtrueexpstub}: The upper bound for experimental group RMST.
//'     
//'     - \code{simtrueexpstse}: The standard error for experimental group 
//'       RMST.
//'     
//'     - \code{simtrue_coxwbprog_hr}: The treatment hazard ratio from the 
//'       Cox model adjusting for baseline prognosis.
//'     
//'     - \code{simtrue_cox_hr}: The treatment hazard ratio from the Cox 
//'       model without adjusting for baseline prognosis.
//'
//' * \code{paneldata}: A counting process style data frame with the 
//'   following variables:
//'
//'     - \code{id}: The subject ID.
//'
//'     - \code{trtrand}: The randomized treatment arm.
//'
//'     - \code{bprog}: Whether the patient had poor baseline prognosis. 
//'     
//'     - \code{tstart}: The left end of time interval.
//'     
//'     - \code{tstop}: The right end of time interval.
//'     
//'     - \code{died}: Whether the patient died. 
//'     
//'     - \code{progressed}: Whether the patient had disease progression. 
//'     
//'     - \code{timePFSobs}: The observed time of disease progression at 
//'       regular scheduled visits.
//'       
//'     - \code{progtdc}: The time-dependent covariate for progression.
//'       
//'     - \code{catevent}: Whether the patient developed metastatic disease.
//'     
//'     - \code{cattime}: When the patient developed metastatic disease.
//'     
//'     - \code{cattdc}: The time-dependent covariate for cat event.
//'     
//'     - \code{catlag}: The lagged value of \code{cattdc}.
//'     
//'     - \code{xo}: Whether the patient switched treatment. 
//'     
//'     - \code{xotime}: When the patient switched treatment.
//'     
//'     - \code{xotdc}: The time-dependent covariate for treatment 
//'       switching.
//'       
//'     - \code{xotime_upper}: The upper bound of treatment switching time.
//'     
//'     - \code{censor_time}: The administrative censoring time.
//'     
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' sim1 <- tsegestsim(
//'   n = 500, allocation1 = 2, allocation2 = 1, pbprog = 0.5, 
//'   trtlghr = -0.5, bprogsl = 0.3, shape1 = 1.8, 
//'   scale1 = 0.000025, shape2 = 1.7, scale2 = 0.000015, 
//'   pmix = 0.5, admin = 5000, pcatnotrtbprog = 0.5, 
//'   pcattrtbprog = 0.25, pcatnotrt = 0.2, pcattrt = 0.1, 
//'   catmult = 0.5, tdxo = 1, ppoor = 0.1, pgood = 0.04, 
//'   ppoormet = 0.4, pgoodmet = 0.2, xomult = 1.4188308, 
//'   milestone = 546, outputRawDataset = 1, seed = 2000)
//'
//' @export
// [[Rcpp::export]]
List tsegestsim(const int n = 500,
                const int allocation1 = 2,
                const int allocation2 = 1,
                const double pbprog = 0.5,
                const double trtlghr = -0.5,
                const double bprogsl = 0.3,
                const double shape1 = 1.8,
                const double scale1 = 0.000025,
                const double shape2 = 1.7,
                const double scale2 = 0.000015,
                const double pmix = 0.5,
                const double admin = 5000,
                const double pcatnotrtbprog = 0.5,
                const double pcattrtbprog = 0.25,
                const double pcatnotrt = 0.2,
                const double pcattrt = 0.1,
                const double catmult = 0.5,
                const double tdxo = 1,
                const double ppoor = 0.1,
                const double pgood = 0.04,
                const double ppoormet = 0.4,
                const double pgoodmet = 0.2,
                const double xomult = 1.4188308,
                const double milestone = 546,
                const bool swtrt_control_only = 1,
                const bool outputRawDataset = 1,
                const int seed = NA_INTEGER) {
  
  int i, j, k;

  if (seed != NA_INTEGER) {
    set_seed(seed);
  }

  // survival function of the Weibull mixture
  auto S = [shape1, scale1, shape2, scale2, pmix](double t)->double {
    double a1 = pmix*exp(-scale1*pow(t,shape1));
    double a2 = (1-pmix)*exp(-scale2*pow(t,shape2));
    return a1+a2;
  };

  double eta, u, v;
  double simtrueconstmean, simtrueconstlb, simtrueconstub, simtrueconstse;
  double simtrueexpstmean, simtrueexpstlb, simtrueexpstub, simtrueexpstse;
  double simtrue_coxwbprog_hr, simtrue_cox_hr;

  IntegerVector id(n), trtrand(n), bprog(n), dead(n), progressed(n);
  NumericVector timeOS(n), timeOS5(n);
  NumericVector timePFS(n), timePFSobs(n);
  NumericVector probcat(n);
  IntegerVector cat2(n, NA_INTEGER), cat3(n, NA_INTEGER);
  IntegerVector cat4(n, NA_INTEGER), cat5(n, NA_INTEGER);
  IntegerVector catevent(n, NA_INTEGER);
  NumericVector catOSloss(n, NA_REAL), cattime(n, NA_REAL);

  int b1 = allocation1, b2 = allocation2;
  for (i=0; i<n; i++) {
    id[i] = i+1;

    // generate treatment indicators using stratified block randomization
    u = R::runif(0,1);
    if (u <= b1/(b1+b2+0.0)) {
      trtrand[i] = 1;
      b1--;
    } else {
      trtrand[i] = 0;
      b2--;
    }

    // start a new block after depleting the current block
    if (b1+b2==0) {
      b1 = allocation1;
      b2 = allocation2;
    }

    // generate poor baseline prognosis indicators
    bprog[i] = static_cast<int>(R::rbinom(1, pbprog));

    // generate survival times from Weibull mixture
    eta = trtlghr*trtrand[i] + bprogsl*bprog[i];
    u = R::runif(0,1);
    v = pow(u, exp(-eta));
    timeOS[i] = squantilecpp(S, v);
    timeOS[i] = std::round(timeOS[i]);
    if (timeOS[i] == 0) timeOS[i] = 1;

    dead[i] = 1;

    // generate observed time to disease progression
    u = R::rbeta(5,10);
    timePFS[i] = std::round(timeOS[i]*u);

    // scheduled visits are every 21 days
    k = static_cast<int>(std::floor(timeOS[i] / 21));
    for (j=1; j<=k; j++) {
      if (timePFS[i] <= j*21) {
        timePFSobs[i] = j*21;
        break;
      }
    }

    if (j<=k) {
      progressed[i] = 1;
    } else {
      timePFSobs[i] = timeOS[i];
      progressed[i] = 0;
    }

    // apply censoring to observed time to disease progression
    if (timePFSobs[i] > admin) {
      progressed[i] = 0;
      timePFSobs[i] = admin;
    }

    // generate catastrophic event that can happen after progression
    if (trtrand[i] == 0 && bprog[i] == 1) {
      probcat[i] = pcatnotrtbprog;
    } else if (trtrand[i] == 1 && bprog[i] == 1) {
      probcat[i] = pcattrtbprog;
    } else if (trtrand[i] == 0 && bprog[i] == 0) {
      probcat[i] = pcatnotrt;
    } else if (trtrand[i] == 1 && bprog[i] == 0) {
      probcat[i] = pcattrt;
    }

    if (progressed[i] == 1 && timeOS[i] > timePFSobs[i] + 21) {
      cat2[i] = static_cast<int>(R::rbinom(1, probcat[i]));
    }

    if (cat2[i] == 0 &&
        progressed[i] == 1 && timeOS[i] > timePFSobs[i] + 42) {
      cat3[i] = static_cast<int>(R::rbinom(1, probcat[i]));
    }

    if (cat2[i] == 0 && cat3[i] == 0 &&
        progressed[i] == 1 && timeOS[i] > timePFSobs[i] + 63) {
      cat4[i] = static_cast<int>(R::rbinom(1, probcat[i]));
    }

    if (cat2[i] == 0 && cat3[i] == 0 && cat4[i] == 0 &&
        progressed[i] == 1 && timeOS[i] > timePFSobs[i] + 84) {
      cat5[i] = static_cast<int>(R::rbinom(1, probcat[i]));
    }

    if (cat2[i] == 1) {
      cattime[i] = timePFSobs[i] + 21;
    } else if (cat3[i] == 1) {
      cattime[i] = timePFSobs[i] + 42;
    } else if (cat4[i] == 1) {
      cattime[i] = timePFSobs[i] + 63;
    } else if (cat5[i] == 1) {
      cattime[i] = timePFSobs[i] + 84;
    }

    // amend survival times
    timeOS5[i] = timeOS[i]; // original survival time
    if (cat2[i] == 1 || cat3[i] == 1 || cat4[i] == 1 || cat5[i] == 1) {
      catevent[i] = 1;
      catOSloss[i] = timeOS[i] - cattime[i];
      catOSloss[i] = std::round(catOSloss[i]*catmult);
      timeOS[i] = catOSloss[i] + cattime[i];
    }
  }

  // calculate HR and RMST with no switching
  IntegerVector event(n);
  NumericVector time(n);
  for (i=0; i<n; i++) {
    if (timeOS[i] > admin) {
      event[i] = 0;
      time[i] = admin;
    } else {
      event[i] = 1;
      time[i] = timeOS[i];
    }
  }

  DataFrame a1 = DataFrame::create(
    Named("time") = time,
    Named("event") = event,
    Named("trtrand") = trtrand,
    Named("bprog") = bprog);

  // RMST
  DataFrame a2 = rmest(a1, "trtrand", "", "time", "event", milestone,
                       0.95, 0);
  NumericVector simtruestmean = a2["rmst"];
  NumericVector simtruestlb = a2["lower"];
  NumericVector simtruestub = a2["upper"];
  NumericVector simtruestse = a2["stderr"];
  simtrueconstmean = simtruestmean[0];
  simtrueexpstmean = simtruestmean[1];
  simtrueconstlb = simtruestlb[0];
  simtrueexpstlb = simtruestlb[1];
  simtrueconstub = simtruestub[0];
  simtrueexpstub = simtruestub[1];
  simtrueconstse = simtruestse[0];
  simtrueexpstse = simtruestse[1];

  // HR with bprog
  StringVector covariates(2);
  covariates[0] = "trtrand";
  covariates[1] = "bprog";

  List a3 = phregcpp(a1, "", "", "time", "", "event", covariates,
                     "", "", "", "efron", 0, 0, 0, 0, 0, 0.05);
  DataFrame parest = DataFrame(a3["parest"]);
  NumericVector beta = parest["beta"];
  simtrue_coxwbprog_hr = exp(beta[0]);

  // HR without bprog
  a3 = phregcpp(a1, "", "", "time", "", "event", "trtrand",
                "", "", "", "efron", 0, 0, 0, 0, 0, 0.05);
  parest = DataFrame(a3["parest"]);
  beta = parest["beta"];
  simtrue_cox_hr = exp(beta[0]);

  // apply switch and effect
  IntegerVector xo1(n, NA_INTEGER), xo2(n, NA_INTEGER), xo3(n, NA_INTEGER);
  IntegerVector xo4(n, NA_INTEGER), xo5(n, NA_INTEGER), xo6(n, NA_INTEGER);
  IntegerVector xo(n, NA_INTEGER), xoprecat(n, NA_INTEGER);
  NumericVector xotime(n, NA_REAL);
  IntegerVector extra2v1(n), extra3v1(n), extra4v1(n), extra5v1(n);
  IntegerVector extraobsv1(n);
  NumericVector xoOSgainobs(n, NA_REAL), timeOS2(n);
  IntegerVector extra2v2(n), extra3v2(n), extra4v2(n), extra5v2(n);
  IntegerVector extraobsv2(n);
  for (i=0; i<n; i++) {
    double p1, p2, p3, p4, p5, p6;

    // prob of switching depends on bprog for the first 2 visits after PD
    p1 = bprog[i] == 1 ? ppoor : pgood;
    p2 = p1;

    // for subsequent visits, prob of switching depends on bprog and 
    // metastatic disease status at the previous visit
    if (bprog[i] == 1 && cat2[i] == 1) {
      p3 = ppoormet;
    } else if (bprog[i] == 0 && cat2[i] == 1) {
      p3 = pgoodmet;
    } else if (bprog[i] == 1 && cat2[i] != 1) {
      p3 = ppoor;
    } else {
      p3 = pgood;
    }

    if (bprog[i] == 1 && (cat2[i] == 1 || cat3[i] == 1)) {
      p4 = ppoormet;
    } else if (bprog[i] == 0 && (cat2[i] == 1 || cat3[i] == 1)) {
      p4 = pgoodmet;
    } else if (bprog[i] == 1 && (cat2[i] != 1 && cat3[i] != 1)) {
      p4 = ppoor;
    } else {
      p4 = pgood;
    }

    if (bprog[i] == 1 && (cat2[i] == 1 || cat3[i] == 1 || cat4[i] == 1)) {
      p5 = ppoormet;
    } else if (bprog[i] == 0 && (cat2[i] == 1 || cat3[i] == 1 ||
      cat4[i] == 1)) {
      p5 = pgoodmet;
    } else if (bprog[i] == 1 && (cat2[i] != 1 && cat3[i] != 1 &&
      cat4[i] != 1)) {
      p5 = ppoor;
    } else {
      p5 = pgood;
    }

    if (bprog[i] == 1 && (cat2[i] == 1 || cat3[i] == 1 || cat4[i] == 1 ||
        cat5[i] == 1)) {
      p6 = ppoormet;
    } else if (bprog[i] == 0 && (cat2[i] == 1 || cat3[i] == 1 ||
      cat4[i] == 1 || cat5[i] == 1)) {
      p6 = pgoodmet;
    } else if (bprog[i] == 1 && (cat2[i] != 1 && cat3[i] != 1 &&
      cat4[i] != 1 && cat5[i] != 1)) {
      p6 = ppoor;
    } else {
      p6 = pgood;
    }

    // generate the indicator and time of treatment switching
    if (swtrt_control_only) {
      if (trtrand[i] == 0 &&
          timeOS[i] > timePFSobs[i] && progressed[i] == 1) {
        xo1[i] = static_cast<int>(R::rbinom(1, p1));
      }
      
      if (trtrand[i] == 0 && xo1[i] == 0 &&
          timeOS[i] > timePFSobs[i] + 21 && progressed[i] == 1 && 
          tdxo == 1) {
        xo2[i] = static_cast<int>(R::rbinom(1, p2));
      }
      
      if (trtrand[i] == 0 && xo1[i] == 0 && xo2[i] == 0 &&
          timeOS[i] > timePFSobs[i] + 42 && progressed[i] == 1 && 
          tdxo == 1) {
        xo3[i] = static_cast<int>(R::rbinom(1, p3));
      }
      
      if (trtrand[i] == 0 && xo1[i] == 0 && xo2[i] == 0 && xo3[i] == 0 &&
          timeOS[i] > timePFSobs[i] + 63 && progressed[i] == 1 && 
          tdxo == 1) {
        xo4[i] = static_cast<int>(R::rbinom(1, p4));
      }
      
      if (trtrand[i] == 0 && xo1[i] == 0 && xo2[i] == 0 && xo3[i] == 0 &&
          xo4[i] == 0 &&
          timeOS[i] > timePFSobs[i] + 84 && progressed[i] == 1 && 
          tdxo == 1) {
        xo5[i] = static_cast<int>(R::rbinom(1, p5));
      }
      
      if (trtrand[i] == 0 && xo1[i] == 0 && xo2[i] == 0 && xo3[i] == 0 &&
          xo4[i] == 0 && xo5[i] == 0 &&
          timeOS[i] > timePFSobs[i] + 105 && progressed[i] == 1 && 
          tdxo == 1) {
        xo6[i] = static_cast<int>(R::rbinom(1, p6));
      }
    } else {
      if (timeOS[i] > timePFSobs[i] && progressed[i] == 1) {
        xo1[i] = static_cast<int>(R::rbinom(1, p1));
      }
      
      if (xo1[i] == 0 &&
          timeOS[i] > timePFSobs[i] + 21 && progressed[i] == 1 && 
          tdxo == 1) {
        xo2[i] = static_cast<int>(R::rbinom(1, p2));
      }
      
      if (xo1[i] == 0 && xo2[i] == 0 &&
          timeOS[i] > timePFSobs[i] + 42 && progressed[i] == 1 && 
          tdxo == 1) {
        xo3[i] = static_cast<int>(R::rbinom(1, p3));
      }
      
      if (xo1[i] == 0 && xo2[i] == 0 && xo3[i] == 0 &&
          timeOS[i] > timePFSobs[i] + 63 && progressed[i] == 1 && 
          tdxo == 1) {
        xo4[i] = static_cast<int>(R::rbinom(1, p4));
      }
      
      if (xo1[i] == 0 && xo2[i] == 0 && xo3[i] == 0 && xo4[i] == 0 &&
          timeOS[i] > timePFSobs[i] + 84 && progressed[i] == 1 && 
          tdxo == 1) {
        xo5[i] = static_cast<int>(R::rbinom(1, p5));
      }
      
      if (xo1[i] == 0 && xo2[i] == 0 && xo3[i] == 0 && xo4[i] == 0 && 
          xo5[i] == 0 &&
          timeOS[i] > timePFSobs[i] + 105 && progressed[i] == 1 && 
          tdxo == 1) {
        xo6[i] = static_cast<int>(R::rbinom(1, p6));
      }
    }

    if (xo1[i] == 1) xotime[i] = timePFSobs[i];
    if (xo2[i] == 1) xotime[i] = timePFSobs[i] + 21;
    if (xo3[i] == 1) xotime[i] = timePFSobs[i] + 42;
    if (xo4[i] == 1) xotime[i] = timePFSobs[i] + 63;
    if (xo5[i] == 1) xotime[i] = timePFSobs[i] + 84;
    if (xo6[i] == 1) xotime[i] = timePFSobs[i] + 105;

    if (xo1[i] == 1 || xo2[i] == 1 || xo3[i] == 1 || xo4[i] == 1 ||
        xo5[i] == 1 || xo6[i] == 1) {
      xo[i] = 1;
    }

    // apply switch effect first through reducing prob. of catastrophic event
    if (xo[i] == 1 && (xotime[i] < cattime[i] || 
        catevent[i] == NA_INTEGER)) {
      xoprecat[i] = 1;
    }

    // after treatment switching for control patients, prob of developing 
    // metastatic disease only depends on baseline prognosis
    if (trtrand[i] == 0 && xo[i] == 1 && bprog[i] == 1 && xoprecat[i] == 1) {
      probcat[i] = pcattrtbprog;
    }
    
    if (trtrand[i] == 0 && xo[i] == 1 && bprog[i] == 0 && xoprecat[i] == 1) {
      probcat[i] = pcattrt;
    }
    
    if (!swtrt_control_only) {
      if (trtrand[i] == 1 && xo[i] == 1 && bprog[i] == 1 && 
          xoprecat[i] == 1) {
        probcat[i] = pcatnotrtbprog;
      }
      
      if (trtrand[i] == 1 && xo[i] == 1 && bprog[i] == 0 && 
          xoprecat[i] == 1) {
        probcat[i] = pcatnotrt;
      }
    }
    
    // regenerate metastatic disease status after treatment switching
    if (xo[i] == 1 && xoprecat[i] == 1 && xotime[i] == timePFSobs[i]) {
      cat2[i] = NA_INTEGER;
    }
    
    if (xo[i] == 1 && xoprecat[i] == 1 && (xotime[i] == timePFSobs[i] ||
        xotime[i] == timePFSobs[i] + 21)) {
      cat3[i] = NA_INTEGER;
    }
    
    if (xo[i] == 1 && xoprecat[i] == 1 && (xotime[i] == timePFSobs[i] ||
        xotime[i] == timePFSobs[i] + 21 ||
        xotime[i] == timePFSobs[i] + 42)) {
      cat4[i] = NA_INTEGER;
    }
    
    if (xo[i] == 1 && xoprecat[i] == 1 && (xotime[i] == timePFSobs[i] ||
        xotime[i] == timePFSobs[i] + 21 || xotime[i] == timePFSobs[i] + 42 ||
        xotime[i] == timePFSobs[i] + 63)) {
      cat5[i] = NA_INTEGER;
    }
    
    if (progressed[i] == 1 && timeOS[i] > timePFSobs[i] + 21 &&
        xo[i] == 1 && xotime[i] == timePFSobs[i]) {
      cat2[i] = static_cast<int>(R::rbinom(1, probcat[i]));
    }
    
    if (cat2[i] == 0 && progressed[i] == 1 &&
        timeOS[i] > timePFSobs[i] + 42 && xo[i] == 1 &&
        (xotime[i] == timePFSobs[i] || xotime[i] == timePFSobs[i] + 21)) {
      cat3[i] = static_cast<int>(R::rbinom(1, probcat[i]));
    }
    
    if (cat2[i] == 0 && cat3[i] == 0 && progressed[i] == 1 &&
        timeOS[i] > timePFSobs[i] + 63 && xo[i] == 1 &&
        (xotime[i] == timePFSobs[i] || xotime[i] == timePFSobs[i] + 21 ||
        xotime[i] == timePFSobs[i] + 42)) {
      cat4[i] = static_cast<int>(R::rbinom(1, probcat[i]));
    }
    
    if (cat2[i] == 0 && cat3[i] == 0 && cat4[i] == 0 && progressed[i] == 1 &&
        timeOS[i] > timePFSobs[i] + 84 && xo[i] == 1 &&
        (xotime[i] == timePFSobs[i] || xotime[i] == timePFSobs[i] + 21 ||
        xotime[i] == timePFSobs[i] + 42 ||
        xotime[i] == timePFSobs[i] + 63)) {
      cat5[i] = static_cast<int>(R::rbinom(1, probcat[i]));
    }
    
    if (xo[i] == 1 && xoprecat[i] == 1) {
      catevent[i] = NA_INTEGER;
      cattime[i] = NA_REAL;
      catOSloss[i] = NA_REAL;
    }
    
    if (cat2[i] == 1 || cat3[i] == 1 || cat4[i] == 1 || cat5[i] == 1) {
      catevent[i] = 1;
    }
    
    if (cat2[i] == 1 && xo[i] == 1 && xotime[i] == timePFSobs[i]) {
      cattime[i] = timePFSobs[i] + 21;
    }
    
    if (cat3[i] == 1 && xo[i] == 1 && (xotime[i] == timePFSobs[i] ||
        xotime[i] == timePFSobs[i] + 21)) {
      cattime[i] = timePFSobs[i] + 42;
    }
    
    if (cat4[i] == 1 && xo[i] == 1 && (xotime[i] == timePFSobs[i] ||
        xotime[i] == timePFSobs[i] + 21 ||
        xotime[i] == timePFSobs[i] + 42)) {
      cattime[i] = timePFSobs[i] + 63;
    }
    
    if (cat5[i] == 1 && xo[i] == 1 && (xotime[i] == timePFSobs[i] ||
        xotime[i] == timePFSobs[i] + 21 || xotime[i] == timePFSobs[i] + 42 ||
        xotime[i] == timePFSobs[i] + 63)) {
      cattime[i] = timePFSobs[i] + 84;
    }
    
    // modify survival time after developing metastatic disease
    if (catevent[i] == 1 && xo[i] == 1 && xoprecat[i] == 1) {
      catOSloss[i] = timeOS5[i] - cattime[i];
      catOSloss[i] = std::round(catOSloss[i]*catmult);
    }
    
    if (catevent[i] == 1 && xoprecat[i] == 1) {
      timeOS[i] = catOSloss[i] + cattime[i];
    }
    
    // if the patient did not develop metastatic disease after all, use
    // original survival time
    if (catevent[i] != 1 && xoprecat[i] == 1) {
      timeOS[i] = timeOS5[i];
    }
    
    // apply chance of catastrophic event for people who don't have cat event
    // or have later cat event and now live past an additional observation.
    // Only replacing for extra observations - others don't change
    if (progressed[i] == 1 && timeOS[i] > timePFSobs[i] + 21 &&
        xo[i] == 1 && xotime[i] == timePFSobs[i] && cat2[i] == NA_INTEGER) {
      extra2v1[i] = 1;
    }
    
    if (cat2[i] == 0 && progressed[i] == 1 &&
        timeOS[i] > timePFSobs[i] + 42 &&
        xo[i] == 1 && (xotime[i] == timePFSobs[i] ||
        xotime[i] == timePFSobs[i] + 21) && cat3[i] == NA_INTEGER) {
      extra3v1[i] = 1;
    }
    
    if (cat2[i] == 0 && cat3[i] == 0 && progressed[i] == 1 &&
        timeOS[i] > timePFSobs[i] + 63 &&
        xo[i] == 1 && (xotime[i] == timePFSobs[i] ||
        xotime[i] == timePFSobs[i] + 21 ||
        xotime[i] == timePFSobs[i] + 42) && cat4[i] == NA_INTEGER) {
      extra4v1[i] = 1;
    }
    
    if (cat2[i] == 0 && cat3[i] == 0 && cat4[i] == 0 && progressed[i] == 1 &&
        timeOS[i] > timePFSobs[i] + 84 &&
        xo[i] == 1 && (xotime[i] == timePFSobs[i] ||
        xotime[i] == timePFSobs[i] + 21 || xotime[i] == timePFSobs[i] + 42 ||
        xotime[i] == timePFSobs[i] + 63) && cat5[i] == NA_INTEGER) {
      extra5v1[i] = 1;
    }
    
    if (extra2v1[i] == 1 || extra3v1[i] == 1 || extra4v1[i] == 1 ||
        extra5v1[i] == 1) {
      extraobsv1[i] = 1;
    }
    
    if (progressed[i] == 1 && timeOS[i] > timePFSobs[i] + 21 && xo[i] == 1 &&
        xotime[i] == timePFSobs[i] && cat2[i] == NA_INTEGER) {
      cat2[i] = static_cast<int>(R::rbinom(1, probcat[i]));
    }
    
    if (cat2[i] == 0 && progressed[i] == 1 &&
        timeOS[i] > timePFSobs[i] + 42 && xo[i] == 1 &&
        (xotime[i] == timePFSobs[i] || xotime[i] == timePFSobs[i] +21) &&
        cat3[i] == NA_INTEGER) {
      cat3[i] = static_cast<int>(R::rbinom(1, probcat[i]));
    }
    
    if (cat2[i] == 0 && cat3[i] == 0 && progressed[i] == 1 &&
        timeOS[i] > timePFSobs[i] + 63 && xo[i] == 1 &&
        (xotime[i] == timePFSobs[i] || xotime[i] == timePFSobs[i] + 21 ||
        xotime[i] == timePFSobs[i] + 42) && cat4[i] == NA_INTEGER) {
      cat4[i] = static_cast<int>(R::rbinom(1, probcat[i]));
    }
    
    if (cat2[i] == 0 && cat3[i] == 0 && cat4[i] == 0 && progressed[i] == 1 &&
        timeOS[i] > timePFSobs[i] + 84 && xo[i] == 1 &&
        (xotime[i] == timePFSobs[i] || xotime[i] == timePFSobs[i] + 21 ||
        xotime[i] == timePFSobs[i] + 42 ||
        xotime[i] == timePFSobs[i] + 63) && cat5[i] == NA_INTEGER) {
      cat5[i] = static_cast<int>(R::rbinom(1, probcat[i]));
    }
    
    if (xo[i] == 1 && xoprecat[i] == 1 && extraobsv1[i] == 1) {
      cattime[i] = NA_REAL;
    }
    
    if (cat2[i] == 1 && xo[i] == 1 && xotime[i] == timePFSobs[i] &&
        extraobsv1[i] == 1) {
      cattime[i] = timePFSobs[i] + 21;
    }
    
    if (cat3[i] == 1 && xo[i] == 1 && (xotime[i] == timePFSobs[i] ||
        xotime[i] == timePFSobs[i] + 21) && extraobsv1[i] == 1) {
      cattime[i] = timePFSobs[i] + 42;
    }
    
    if (cat4[i] == 1 && xo[i] == 1 && (xotime[i] == timePFSobs[i] ||
        xotime[i] == timePFSobs[i] + 21 ||
        xotime[i] == timePFSobs[i] + 42) && extraobsv1[i] == 1) {
      cattime[i] = timePFSobs[i] + 63;
    }
    
    if (cat5[i] == 1 && xo[i] == 1 && (xotime[i] == timePFSobs[i] ||
        xotime[i] == timePFSobs[i] + 21 || xotime[i] == timePFSobs[i] + 42 ||
        xotime[i] == timePFSobs[i] + 63) && extraobsv1[i] == 1) {
      cattime[i] = timePFSobs[i] + 84;
    }
    
    if (xo[i] == 1 && xoprecat[i] == 1 && extraobsv1[i] == 1) {
      catevent[i] = NA_INTEGER;
    }
    
    if ((cat2[i] == 1 || cat3[i] == 1 || cat4[i] == 1 || cat5[i] == 1) &&
        extraobsv1[i] == 1) {
      catevent[i] = 1;
    }
    
    // remember, these are only patients who hadn't had a cat event before,
    // and now may have done in their additional observations. So survival
    // times cannot go up, this only allows for possible additional events
    if (xo[i] == 1 && xoprecat[i] == 1 && extraobsv1[i] == 1) {
      catOSloss[i] = NA_REAL;
    }
    
    if (catevent[i] == 1 && xo[i] == 1 && xoprecat[i] == 1 &&
        extraobsv1[i] == 1) {
      catOSloss[i] = timeOS5[i] - cattime[i];
      catOSloss[i] = std::round(catOSloss[i]*catmult);
    }
    
    if (catevent[i] == 1 && xoprecat[i] == 1 && extraobsv1[i] == 1) {
      timeOS[i] = catOSloss[i] + cattime[i];
    }
    
    if (catevent[i] != 1 && xoprecat[i] == 1 && extraobsv1[i] == 1) {
      timeOS[i] = timeOS5[i];
    }
    
    // Apply switch effect second through xomult
    if (xo[i] == 1) {
      xoOSgainobs[i] = timeOS[i] - xotime[i];
      xoOSgainobs[i] = std::round(xoOSgainobs[i]*xomult);
    }
    
    timeOS2[i] = xo[i] == 1 ? xoOSgainobs[i] + xotime[i] : timeOS[i];
    
    // apply chance of catastrophic event for xo patients who now live past
    // an additional observation. Only replacing for extra observations in
    // switchers, so no additional chance to switch and OS can't increase
    // as currently in these observations there will zero chance of cat event
    if (progressed[i] == 1 && timeOS2[i] > timePFSobs[i] + 21 &&
        xo[i] == 1 && xotime[i] == timePFSobs[i] && cat2[i] == NA_INTEGER) {
      extra2v2[i] = 1;
    }
    
    if (cat2[i] == 0 && progressed[i] == 1 &&
        timeOS2[i] > timePFSobs[i] + 42 &&
        xo[i] == 1 && (xotime[i] == timePFSobs[i] ||
        xotime[i] == timePFSobs[i] + 21) && cat3[i] == NA_INTEGER) {
      extra3v2[i] = 1;
    }
    
    if (cat2[i] == 0 && cat3[i] == 0 && progressed[i] == 1 &&
        timeOS2[i] > timePFSobs[i] + 63 &&
        xo[i] == 1 && (xotime[i] == timePFSobs[i] ||
        xotime[i] == timePFSobs[i] + 21 ||
        xotime[i] == timePFSobs[i] + 42) && cat4[i] == NA_INTEGER) {
      extra4v2[i] = 1;
    }
    
    if (cat2[i] == 0 && cat3[i] == 0 && cat4[i] == 0 && progressed[i] == 1 &&
        timeOS2[i] > timePFSobs[i] + 84 &&
        xo[i] == 1 && (xotime[i] == timePFSobs[i] ||
        xotime[i] == timePFSobs[i] + 21 || xotime[i] == timePFSobs[i] + 42 ||
        xotime[i] == timePFSobs[i] + 63) && cat5[i] == NA_INTEGER) {
      extra5v2[i] = 1;
    }
    
    if (extra2v2[i] == 1 || extra3v2[i] == 1 || extra4v2[i] == 1 ||
        extra5v2[i] == 1) {
      extraobsv2[i] = 1;
    }
    
    if (progressed[i] == 1 && timeOS2[i] > timePFSobs[i] + 21 &&
        xo[i] == 1 && xotime[i] == timePFSobs[i] && cat2[i] == NA_INTEGER) {
      cat2[i] = static_cast<int>(R::rbinom(1, probcat[i]));
    }
    
    if (cat2[i] == 0 && progressed[i] == 1 &&
        timeOS2[i] > timePFSobs[i] + 42 && xo[i] == 1 &&
        (xotime[i] == timePFSobs[i] || xotime[i] == timePFSobs[i] +21) &&
        cat3[i] == NA_INTEGER) {
      cat3[i] = static_cast<int>(R::rbinom(1, probcat[i]));
    }
    
    if (cat2[i] == 0 && cat3[i] == 0 && progressed[i] == 1 &&
        timeOS2[i] > timePFSobs[i] + 63 && xo[i] == 1 &&
        (xotime[i] == timePFSobs[i] || xotime[i] == timePFSobs[i] + 21 ||
        xotime[i] == timePFSobs[i] + 42) && cat4[i] == NA_INTEGER) {
      cat4[i] = static_cast<int>(R::rbinom(1, probcat[i]));
    }
    
    if (cat2[i] == 0 && cat3[i] == 0 && cat4[i] == 0 && progressed[i] == 1 &&
        timeOS2[i] > timePFSobs[i] + 84 && xo[i] == 1 &&
        (xotime[i] == timePFSobs[i] || xotime[i] == timePFSobs[i] + 21 ||
        xotime[i] == timePFSobs[i] + 42 ||
        xotime[i] == timePFSobs[i] + 63) && cat5[i] == NA_INTEGER) {
      cat5[i] = static_cast<int>(R::rbinom(1, probcat[i]));
    }
    
    if (xo[i] == 1 && xoprecat[i] == 1 && extraobsv2[i] == 1) {
      cattime[i] = NA_REAL;
    }
    
    if (cat2[i] == 1 && xo[i] == 1 && xotime[i] == timePFSobs[i] &&
        extraobsv2[i] == 1) {
      cattime[i] = timePFSobs[i] + 21;
    }
    
    if (cat3[i] == 1 && xo[i] == 1 && (xotime[i] == timePFSobs[i] ||
        xotime[i] == timePFSobs[i] + 21) && extraobsv2[i] == 1) {
      cattime[i] = timePFSobs[i] + 42;
    }
    
    if (cat4[i] == 1 && xo[i] == 1 && (xotime[i] == timePFSobs[i] ||
        xotime[i] == timePFSobs[i] + 21 ||
        xotime[i] == timePFSobs[i] + 42) && extraobsv2[i] == 1) {
      cattime[i] = timePFSobs[i] + 63;
    }
    
    if (cat5[i] == 1 && xo[i] == 1 && (xotime[i] == timePFSobs[i] ||
        xotime[i] == timePFSobs[i] + 21 || xotime[i] == timePFSobs[i] + 42 ||
        xotime[i] == timePFSobs[i] + 63) && extraobsv2[i] == 1) {
      cattime[i] = timePFSobs[i] + 84;
    }
    
    if (xo[i] == 1 && xoprecat[i] == 1 && extraobsv2[i] == 1) {
      catevent[i] = NA_INTEGER;
    }
    
    if ((cat2[i] == 1 || cat3[i] == 1 || cat4[i] == 1 || cat5[i] == 1) &&
        extraobsv2[i] == 1) {
      catevent[i] = 1;
    }
    
    // remember, these are only patients who hadn't had a cat event before,
    // and now may have done in their additional observations, which
    // previously they did not live past. So if they have an event,
    // cat time will be >OS5. Hence only makes sense to apply cat event
    // survival reduction to extended survival beyond this point, i.e., OS2
    if (xo[i] == 1 && xoprecat[i] == 1 && extraobsv2[i] == 1) {
      catOSloss[i] = NA_REAL;
    }
    
    if (catevent[i] == 1 && xo[i] == 1 && xoprecat[i] == 1 &&
        extraobsv2[i] == 1) {
      catOSloss[i] = timeOS2[i] - cattime[i];
      catOSloss[i] = std::round(catOSloss[i]*catmult);
    }
    
    if (catevent[i] == 1 && xoprecat[i] == 1 && extraobsv2[i] == 1) {
      timeOS2[i] = catOSloss[i] + cattime[i];
    }
    
    // don't replace OS2 for those with an extraobsv2 who did not have a cat
    // event because OS has not changed for them
  }
  
  // apply censoring
  IntegerVector died(n);
  for (i=0; i<n; i++) {
    if (timeOS2[i] <= admin && dead[i] == 1) {
      died[i] = 1;
    } else {
      died[i] = 0;
    }
    
    if (timeOS2[i] > admin) timeOS2[i] = admin;
    
    if (xotime[i] >= admin) {
      xoOSgainobs[i] = NA_REAL;
      xo[i] = NA_INTEGER;
      xotime[i] = NA_REAL;
    }
    
    if (cattime[i] >= admin) {
      catevent[i] = NA_INTEGER;
      cattime[i] = NA_REAL;
    }
  }
  
  // create panel
  NumericVector zero(n);
  int kmax = static_cast<int>(std::ceil(max(timeOS2)/21));
  NumericVector cut(kmax);
  for (k=0; k<kmax; k++) {
    cut[k] = k*21;
  }
  
  DataFrame a = survsplit(zero, timeOS2, cut);
  int n2 = a.nrow();
  IntegerVector q2 = a["row"];
  NumericVector tstart = a["start"];
  NumericVector tstop = a["end"];
  IntegerVector censor = a["censor"];
  
  IntegerVector id2 = id[q2];
  IntegerVector trtrand2 = trtrand[q2];
  IntegerVector bprog2 = bprog[q2];
  IntegerVector died2 = died[q2];
  died2[censor == 1] = 0;
  
  IntegerVector progressed2 = progressed[q2];
  NumericVector timePFSobs2 = timePFSobs[q2];
  IntegerVector catevent2 = catevent[q2];
  NumericVector cattime2 = cattime[q2];
  IntegerVector xoo2 = xo[q2];
  NumericVector xotime2 = xotime[q2];
  
  for (i=0; i<n2; i++) {
    if (progressed2[i] == NA_INTEGER) progressed2[i] = 0;
    if (catevent2[i] == NA_INTEGER) catevent2[i] = 0;
    if (xoo2[i] == NA_INTEGER) xoo2[i] = 0;
  }
  
  // make time-dependent covariates for progression, cat event, and switch
  IntegerVector progtdc(n2), cattdc(n2), xotdc(n2);
  for (i=0; i<n2; i++) {
    if (progressed2[i] == 1 && tstart[i] >= timePFSobs2[i]) progtdc[i] = 1;
    if (catevent2[i] == 1 && tstart[i] >= cattime2[i]) cattdc[i] = 1;
    if (xoo2[i] == 1 && tstart[i] >= xotime2[i]) xotdc[i] = 1;
  }
  
  IntegerVector idx(1,0); // first observation within an id
  for (i=1; i<n2; i++) {
    if (id2[i] != id2[i-1]) {
      idx.push_back(i);
    }
  }
  
  int nids = static_cast<int>(idx.size());
  idx.push_back(n2);
  
  IntegerVector catlag(n2);
  for (i=0; i<nids; i++) {
    catlag[idx[i]] = 0;
    for (j=idx[i]+1; j<idx[i+1]; j++) {
      catlag[j] = cattdc[j-1];
    }
  }
  
  NumericVector xotime_upper(n2);
  for (i=0; i<n2; i++) {
    if (progressed2[i] == 1) {
      xotime_upper[i] = timePFSobs2[i] + 105;
    } else {
      xotime_upper[i] = 1e8;
    }
  }
  
  NumericVector censor_time(n2, admin);
  
  List result;
  
  DataFrame sumstat = DataFrame::create(
    Named("simtrueconstmean") = simtrueconstmean,
    Named("simtrueconstlb") = simtrueconstlb,
    Named("simtrueconstub") = simtrueconstub,
    Named("simtrueconstse") = simtrueconstse,
    Named("simtrueexpstmean") = simtrueexpstmean,
    Named("simtrueexpstlb") = simtrueexpstlb,
    Named("simtrueexpstub") = simtrueexpstub,
    Named("simtrueexpstse") = simtrueexpstse,
    Named("simtrue_coxwbprog_hr") = simtrue_coxwbprog_hr,
    Named("simtrue_cox_hr") = simtrue_cox_hr);
  
  if (outputRawDataset) {
    DataFrame paneldata = DataFrame::create(
      Named("id") = id2,
      Named("trtrand") = trtrand2,
      Named("bprog") = bprog2,
      Named("tstart") = tstart,
      Named("tstop") = tstop,
      Named("died") = died2,
      Named("progressed") = progressed2,
      Named("timePFSobs") = timePFSobs2,
      Named("progtdc") = progtdc,
      Named("catevent") = catevent2,
      Named("cattime") = cattime2,
      Named("cattdc") = cattdc,
      Named("catlag") = catlag,
      Named("xo") = xoo2,
      Named("xotime") = xotime2,
      Named("xotdc") = xotdc,
      Named("xotime_upper") = xotime_upper,
      Named("censor_time") = censor_time);
    
    result = List::create(Named("sumstat") = sumstat,
                          Named("paneldata") = paneldata);
  } else {
    result = List::create(Named("sumstat") = sumstat);
  }
  
  return result;
}
