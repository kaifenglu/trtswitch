#include "utilities.h"
using namespace Rcpp;

//' @title Simulate Data for Treatment Switching
//' @description Simulates data for studies involving treatment switching, 
//' incorporating time-dependent confounding. The generated data can be used 
//' to evaluate methods for handling treatment switching in survival 
//' analysis.
//'
//' @param tdxo Logical indicator for timing of treatment switching:
//'   \itemize{
//'     \item 1: Treatment switching can occur at or after disease 
//'              progression.
//'     \item 0: Treatment switching is restricted to the time of disease 
//'              progression.
//'   }
//' @param coxo Logical indicator for arm-specific treatment switching:
//'   \itemize{
//'     \item 1: Treatment switching occurs only in the control arm.
//'     \item 0: Treatment switching is allowed in both arms.
//'   }
//' @param allocation1 Number of subjects in the active treatment group in
//'   a randomization block. Defaults to 1 for equal randomization.
//' @param allocation2 Number of subjects in the control group in
//'   a randomization block. Defaults to 1 for equal randomization.
//' @param p_X_1 Probability of poor baseline prognosis in the experimental 
//'   arm.
//' @param p_X_0 Probability of poor baseline prognosis in the control arm.
//' @param rate_T Baseline hazard rate for time to death.
//' @param beta1 Log hazard ratio for randomized treatment (\code{R}).
//' @param beta2 Log hazard ratio for baseline covariate (\code{X}).
//' @param gamma0 Intercept for the time-dependent covariate model 
//'   (\code{L}).
//' @param gamma1 Coefficient for lagged treatment switching (\code{Alag}) 
//'   in the \code{L} model.
//' @param gamma2 Coefficient for lagged \code{L} (\code{Llag}) in the 
//'   \code{L} model.
//' @param gamma3 Coefficient for baseline covariate (\code{X}) in the 
//'   \code{L} model.
//' @param gamma4 Coefficient for randomized treatment (\code{R}) in the 
//'   \code{L} model.
//' @param zeta0 Intercept for the disease progression model (\code{Z}).
//' @param zeta1 Coefficient for \code{L} in the \code{Z} model.
//' @param zeta2 Coefficient for baseline covariate (\code{X}) in the 
//'   \code{Z} model.
//' @param zeta3 Coefficient for randomized treatment (\code{R}) in the 
//'   \code{Z} model.
//' @param alpha0 Intercept for the treatment switching model (\code{A}).
//' @param alpha1 Coefficient for \code{L} in the \code{A} model.
//' @param alpha2 Coefficient for baseline covariate (\code{X}) in the 
//'   \code{A} model.
//' @param theta1_1 Negative log time ratio for \code{A} for the 
//'   experimental arm.
//' @param theta1_0 Negative log time ratio for \code{A} for the control arm.
//' @param theta2 Negative log time ratio for \code{L}.
//' @param rate_C Hazard rate for random (dropout) censoring.
//' @param accrualTime A vector that specifies the starting time of 
//'   piecewise Poisson enrollment time intervals. Must start with 0, 
//'   e.g., \code{c(0, 3)} breaks the time axis into 2 accrual intervals: 
//'   [0, 3) and [3, Inf).
//' @param accrualIntensity A vector of accrual intensities. One for 
//'   each accrual time interval.
//' @param followupTime	Follow-up time for a fixed follow-up design.
//' @param fixedFollowup Whether a fixed follow-up design is used. 
//'   Defaults to 0 for variable follow-up.
//' @param plannedTime The calendar time for the analysis.
//' @param days Number of days in each treatment cycle.
//' @param n Number of subjects per simulation.
//' @param NSim Number of simulated datasets.
//' @param seed Random seed for reproducibility.
//'
//' @details 
//' The simulation algorithm is adapted from Xu et al. (2022) and simulates
//' disease progression status while incorporating the multiplicative 
//' effects of both baseline and time-dependent covariates on survival time. 
//' The design options \code{tdxo} and \code{coxo} specify
//' the timing of treatment switching and the study arm eligibility 
//' for switching, respectively. Time is measured in days. 
//' 
//' In a fixed follow-up design, all subjects share the same follow-up
//' duration. In contrast, under a variable follow-up design, follow-up
//' time also depends on each subject's enrollment date.  
//' The number of treatment cycles for a subject is determined by
//' dividing the follow-up time by the number of days in each cycle. 
//' \enumerate{
//'   \item At randomization, subjects are assigned to treatment arms 
//'   using block randomization, with \code{allocation1} patients 
//'   assigned to active treatment and \code{allocation2} to control
//'   within each randomized block. A baseline covariate is 
//'   also generated for each subject:
//'   \deqn{X_i \sim \mbox{Bernoulli}(p_{X_1} R_i + p_{X_0} (1-R_i))}
//'         
//'   \item The initial survival time is drawn
//'   from an exponential distribution with hazard:
//'   \deqn{rate_T \exp(\beta_1 R_i + \beta_2 X_i)}
//'   We define the event indicator at cycle \eqn{j} as
//'   \deqn{Y_{i,j} = I(T_i \leq j\times days)}
//'         
//'   \item The initial states are set to
//'   \eqn{L_{i,0} = 0}, \eqn{Z_{i,0} = 0}, \eqn{A_{i,0} = 0},
//'   \eqn{Y_{i,0} = 0}. For each treatment cycle \eqn{j=1,\ldots,J},
//'   we set \eqn{tstart = (j-1) \times days}.
//'         
//'   \item Generate time-dependent covariates:
//'   \deqn{\mbox{logit} P(L_{i,j}=1|\mbox{history}) = 
//'   \gamma_0 + \gamma_1 A_{i,j-1} + \gamma_2 L_{i,j-1} + 
//'   \gamma_3 X_i + \gamma_4 R_i}
//'      
//'   \item If \eqn{T_i \leq j \times days}, set \eqn{tstop = T_i} and 
//'   \eqn{Y_{i,j} = 1}, which completes data generation
//'   for subject \eqn{i}.
//'           
//'   \item If \eqn{T_i > j \times days}, set \eqn{tstop = j\times days},  
//'   \eqn{Y_{i,j} = 0}, and perform the following before proceeding to 
//'   the next cycle for the subject.
//'   
//'   \item Generate disease progression status: 
//'   If \eqn{Z_{i,j-1} = 0},  
//'   \deqn{\mbox{logit} P(Z_{i,j}=1 | \mbox{history}) = \zeta_0 + 
//'   \zeta_1 L_{i,j} + \zeta_2 X_i + \zeta_3 R_i}
//'   Otherwise, set \eqn{Z_{i,j} = 1}.
//'     
//'   \item Generate alternative therapy status:     
//'   If \eqn{A_{i,j-1} = 0}, determine switching eligibility 
//'   based on design options.        
//'   If switching is allowed:
//'   \deqn{\mbox{logit} P(A_{i,j} = 1 | \mbox{history}) = \alpha_0 + 
//'   \alpha_1 L_{i,j} + \alpha_2 X_i}       
//'   If switching is now allowed, set \eqn{A_{i,j} = 0}.      
//'   If \eqn{A_{i,j-1} = 1}, set \eqn{A_{i,j} = 1} (once switched to 
//'   alternative therapy, remain on alternative therapy).
//'          
//'   \item Update survival time based on changes in alternative 
//'   therapy status and time-dependent covariates:
//'   \deqn{T_i = j\times days + (T_i - j\times days) \exp\{
//'   -(\theta_{1,1}R_i + \theta_{1,0}(1-R_i))(A_{i,j} - A_{i,j-1}) 
//'   -\theta_2 (L_{i,j} - L_{i,j-1})\}}   
//' }
//' 
//' Additional random censoring times are generated from an exponential 
//' distribution with hazard rate \eqn{rate_C}.
//' 
//' An extra record is generated when the minimum of the latent survival 
//' time, the random censoring time, and the administrative censoring time 
//' is greater than the number of regular treatment cycles times 
//' days per cycle.
//' 
//' Finally we apply the lag function so that \eqn{Z_{i,j}} and 
//' \eqn{A_{i,j}} represent the PD status and alternative therapy status 
//' at the start of cycle \eqn{j} (and thus remain appplicable for the 
//' entire cycle \eqn{j}) for subject \eqn{i}.
//' 
//' To estimate the true treatment effect in a Cox marginal 
//' structural model, one can set \eqn{\alpha_0 = -\infty}, which 
//' effectively forces \eqn{A_{i,j} = 0} (disabling treatment switching). 
//' The coefficient for the randomized treatment can then be estimated 
//' using a Cox proportional hazards model.
//'
//' @return
//' A list of data frames, each containing simulated longitudinal covariate, 
//' pd status, alternative therapy status, and event history data with the 
//' following variables:
//'
//'  * \code{id}: Subject identifier.
//'  
//'  * \code{arrivalTime}: The enrollment time for the subject.
//'  
//'  * \code{trtrand}: Randomized treatment assignment (0 = control, 
//'    1 = experimental)
//'
//'  * \code{bprog}: Baseline prognosis (0 = good, 1 = poor).
//'
//'  * \code{tpoint}: Treatment cycle index.
//'
//'  * \code{tstart}: Start day of the treatment cycle.
//'
//'  * \code{tstop}: End day of the treatment cycle.
//'
//'  * \code{L}: Time-dependent covariate at \code{tstart} predicting 
//'    survival and switching; affected by treatment switching.
//'
//'  * \code{Llag}: Lagged value of \code{L}.
//'
//'  * \code{Z}: Disease progression status at \code{tstart}.
//'
//'  * \code{A}: Treatment switching status at \code{tstart}.
//'
//'  * \code{Alag}: Lagged value of \code{A}.
//'
//'  * \code{event}: Death indicator at \code{tstop}.
//'
//'  * \code{timeOS}: Observed time to death or censoring.
//'     
//'  * \code{died}: Indicator of death by end of follow-up.
//'     
//'  * \code{progressed}: Indicator of disease progression by end of 
//'    follow-up.
//'       
//'  * \code{timePD}: Observed time to progression or censoring.
//'       
//'  * \code{xo}: Indicator for whether treatment switching occurred.
//'     
//'  * \code{xotime}: Time of treatment switching (if applicable).
//'     
//'  * \code{censor_time}: Administrative censoring time.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @references
//' 
//' Jessica G. Young, and Eric J. Tchetgen Tchetgen. 
//' Simulation from a known Cox MSM using standard parametric models 
//' for the g-formula.
//' Statistics in Medicine. 2014;33(6):1001-1014. 
//' 
//' NR Latimer, IR White, K Tilling, and U Siebert.
//' Improved two-stage estimation to adjust for treatment switching in 
//' randomised trials: g-estimation to address time-dependent confounding.
//' Statistical Methods in Medical Research. 2020;29(10):2900-2918.
//' 
//' Jing Xu, Guohui Liu, and Bingxia Wang.
//' Bias and type I error control in correcting treatment effect
//' for treatment switching using marginal structural models 
//' in Phse III oncology trials.
//' Journal of Biopharmaceutical Statistics. 2022;32(6):897-914.
//'
//' @examples
//'
//' library(dplyr)
//' 
//' simulated.data <- tssim(
//'   tdxo = 1, coxo = 1, allocation1 = 1, allocation2 = 1,
//'   p_X_1 = 0.3, p_X_0 = 0.3, 
//'   rate_T = 0.002, beta1 = -0.5, beta2 = 0.3, 
//'   gamma0 = 0.3, gamma1 = -0.9, gamma2 = 0.7, gamma3 = 1.1, gamma4 = -0.8,
//'   zeta0 = -3.5, zeta1 = 0.5, zeta2 = 0.2, zeta3 = -0.4, 
//'   alpha0 = 0.5, alpha1 = 0.5, alpha2 = 0.4, 
//'   theta1_1 = -0.4, theta1_0 = -0.4, theta2 = 0.2,
//'   rate_C = 0.0000855, accrualIntensity = 20/30, 
//'   fixedFollowup = FALSE, plannedTime = 1350, days = 30,
//'   n = 500, NSim = 100, seed = 314159)
//'   
//' simulated.data[[1]] %>% filter(id == 1)
//'
//' @export
// [[Rcpp::export]]
List tssim(const bool tdxo = 0, 
           const bool coxo = 1, 
           const int allocation1 = 1,
           const int allocation2 = 1,
           const double p_X_1 = NA_REAL, 
           const double p_X_0 = NA_REAL, 
           const double rate_T = NA_REAL, 
           const double beta1 = NA_REAL, 
           const double beta2 = NA_REAL, 
           const double gamma0 = NA_REAL, 
           const double gamma1 = NA_REAL, 
           const double gamma2 = NA_REAL, 
           const double gamma3 = NA_REAL, 
           const double gamma4 = NA_REAL,
           const double zeta0 = NA_REAL, 
           const double zeta1 = NA_REAL, 
           const double zeta2 = NA_REAL, 
           const double zeta3 = NA_REAL, 
           const double alpha0 = NA_REAL, 
           const double alpha1 = NA_REAL, 
           const double alpha2 = NA_REAL, 
           const double theta1_1 = NA_REAL, 
           const double theta1_0 = NA_REAL, 
           const double theta2 = NA_REAL,
           const double rate_C = NA_REAL,
           const NumericVector& accrualTime = 0,
           const NumericVector& accrualIntensity = NA_REAL,
           const double followupTime = NA_REAL,
           const bool fixedFollowup = 0,
           const double plannedTime = NA_REAL,
           const double days = NA_INTEGER,
           const int n = NA_INTEGER, 
           const int NSim = 1000, 
           const int seed = NA_INTEGER) {
  
  if (allocation1 < 1) {
    stop("allocation1 must be a positive integer");
  }
  
  if (allocation2 < 1) {
    stop("allocation2 must be a positive integer");
  }
  
  if (p_X_1 <= 0 || p_X_1 >= 1) {
    stop("p_X_1 must lie between 0 and 1");
  }
  
  if (p_X_0 <= 0 || p_X_0 >= 1) {
    stop("p_X_0 must lie between 0 and 1");
  }
  
  if (rate_T <= 0) {
    stop("rate_T must be positive");
  }
  
  if (rate_C < 0) {
    stop("rate_C must be nonnegative");
  }
  
  if (accrualTime[0] != 0) {
    stop("accrualTime must start with 0");
  }
  
  if (accrualTime.size() > 1 && is_true(any(diff(accrualTime) <= 0))) {
    stop("accrualTime should be increasing");
  }
  
  if (is_true(any(is_na(accrualIntensity)))) {
    stop("accrualIntensity must be provided");
  }
  
  if (accrualTime.size() != accrualIntensity.size()) {
    stop("accrualTime must have the same length as accrualIntensity");
  }
  
  if (is_true(any(accrualIntensity < 0))) {
    stop("accrualIntensity must be non-negative");
  }
  
  if (fixedFollowup) {
    if (R_isnancpp(followupTime)) {
      stop("followupTime must be provided for fixed follow-up");
    }
    
    if (followupTime <= 0) {
      stop("followupTime must be positive for fixed follow-up");
    }
  }
  
  if (R_isnancpp(plannedTime)) {
    stop("plannedTime must be provided");
  }
  
  if (plannedTime <= 0) {
    stop("plannedTime must be positive");
  }
  
  if (days <= 0) {
    stop("days must be positive");
  }
  
  if (n <= 0) {
    stop("n must be positive");
  }
  
  if (NSim <= 0) {
    stop("NSim must be positive");
  }
  
  if (seed != NA_INTEGER) set_seed(seed);
  int maxFollowup = static_cast<int>(std::ceil(plannedTime/days));
  int K = n*maxFollowup;
  List result(NSim);
  
  for (int iter=0; iter<NSim; iter++) {
    IntegerVector idx(K), trtrandx(K), bprogx(K), tpointx(K);
    IntegerVector Lx(K), Llagx(K), Zx(K), Zlagx(K);
    IntegerVector Ax(K), Alagx(K), Alag2x(K), eventx(K);
    IntegerVector diedx(K), progressedx(K), xox(K);
    NumericVector tstartx(K), tstopx(K), timeOSx(K), timePDx(K);
    NumericVector xotimex(K, NA_REAL), censor_timex(K);
    NumericVector arrivalTimex(K);
    
    int id, trtrand, bprog, tpoint, L, Llag, Z, Zlag, A, Alag, Alag2, event;
    double arrivalTime, tstart, tstop;
    double b1 = allocation1, b2 = allocation2;
    double enrollt = 0;
    int k = 0;
    for (int i=1; i<=n; i++) { // i = subject index
      id = i;
      
      // generate accrual time
      double u = R::runif(0,1);
      enrollt = qtpwexpcpp1(u, accrualTime, accrualIntensity, enrollt, 1, 0);
      arrivalTime = std::ceil(enrollt);
      
      // stratified block randomization
      u = R::runif(0,1);
      if (u <= b1/(b1+b2)) {
        trtrand = 1;
        b1--;
      } else {
        trtrand = 0;
        b2--;
      }
      if (b1 == 0 && b2 == 0) {
        b1 = allocation1; b2 = allocation2;
      }
      
      bprog = (R::runif(0,1) <= p_X_1*trtrand + p_X_0*(1-trtrand));
      double T = R::rexp(1.0/(rate_T*exp(beta1*trtrand + beta2*bprog)));
      T = std::round(T);
      if (T == 0) T = 1;
      
      double fu = std::max(plannedTime - arrivalTime, 0.0);
      if (fixedFollowup) fu = std::min(followupTime, fu);
      int followup = static_cast<int>(std::floor(fu/days));
      
      tpoint = followup; L = 0; Llag = 0; Z = 0; A = 0; Alag = 0;
      for (int j=1; j<=followup; j++) { // j = cycle index
        tpoint = j; tstart = days*(j-1);
        if (j == 1) {
          Llag = 0; Zlag = 0; Alag2 = 0; Alag = 0;
        } else {
          Llag = L; Zlag = Z; Alag2 = Alag; Alag = A;
        }
        
        // generate time-dependent covariate
        L = (R::runif(0,1) <= R::plogis(gamma0 + gamma1*Alag + gamma2*Llag + 
          gamma3*bprog + gamma4*trtrand, 0, 1, 1, 0));
        
        if (T <= days*j) { // died in cycle j, complete data for the subject
          tstop = T; event = 1;
          Z = NA_INTEGER; A = NA_INTEGER;
        } else { // alive at the end of cycle j, continue to the next cycle
          tstop = days*j; event = 0;
          
          // generate disease progression status
          if (Zlag == 0) {
            Z = (R::runif(0,1) <= R::plogis(zeta0 + zeta1*L + zeta2*bprog + 
              zeta3*trtrand, 0, 1, 1, 0));
          } else {
            Z = 1;
          }
          
          // generate treatment switching status
          if (Alag == 0) {
            if (((tdxo == 0 && Z == 1 && Zlag == 0) || 
                (tdxo == 1 && Z == 1)) && 
                ((coxo == 1 && trtrand == 0) || (coxo == 0))) {
              A = (R::runif(0,1) <= R::plogis(alpha0 + alpha1*L + 
                alpha2*bprog, 0, 1, 1, 0));
            } else {
              A = 0;
            }
          } else {
            A = 1;
          }
          
          // update survival time
          double theta1 = theta1_1*trtrand + theta1_0*(1-trtrand);
          T = days*j + (T - days*j)*exp(-theta1*(A-Alag) - theta2*(L-Llag));
          T = std::round(T);
        }
        
        // add the data from the current cycle
        idx[k] = id;
        arrivalTimex[k] = arrivalTime;
        trtrandx[k] = trtrand;
        bprogx[k] = bprog;
        tpointx[k] = tpoint;
        tstartx[k] = tstart;
        tstopx[k] = tstop;
        Lx[k] = L;
        Llagx[k] = Llag;
        Zx[k] = Z;
        Zlagx[k] = Zlag;
        Ax[k] = A;
        Alagx[k] = Alag;
        Alag2x[k] = Alag2;
        eventx[k] = event;
        
        k++;
        
        if (event == 1) break;
      }
      
      // generate random censoring due to dropout
      double C = R::rexp(1.0/rate_C);
      C = std::round(C);
      if (C == 0) C = 1;
      double time = std::min({T, C, fu});
      
      int J;  // J is the number of treatment cycles
      if (time <= days*followup) {
        J = static_cast<int>(std::ceil(time/days));
        k = k - tpoint + J; // discard treatment cycles after censoring
        tstopx[k-1] = time; // update the ending time and event indicator
        eventx[k-1] = T == time ? 1 : 0;
        Zx[k-1] = NA_INTEGER;
        Ax[k-1] = NA_INTEGER;
      } else { // add one more record 
        J = followup + 1;
        idx[k] = id;
        arrivalTimex[k] = arrivalTime;
        trtrandx[k] = trtrand;
        bprogx[k] = bprog;
        tpointx[k] = followup + 1;
        tstartx[k] = days*followup;
        tstopx[k] = time;
        Lx[k] = (R::runif(0,1) <= R::plogis(gamma0 + gamma1*Alag + 
          gamma2*Llag + gamma3*bprog + gamma4*trtrand, 0, 1, 1, 0));
        Llagx[k] = L;
        Zx[k] = NA_INTEGER;
        Zlagx[k] = Z;
        Ax[k] = NA_INTEGER;
        Alagx[k] = A;
        Alag2x[k] = Alag;
        eventx[k] = T == time ? 1 : 0;
        
        k++;
      }
      
      // create subject-level survival time and death indicator
      for (int j=k-J; j<k; j++) {
        timeOSx[j] = tstopx[k-1];
        diedx[j] = eventx[k-1];
        censor_timex[j] = fu;
      }
      
      // progression and time to progression (if applicable)
      int pd = 0;
      double pd_time = NA_REAL;
      for (int j=k-J; j<k; j++) {
        if (Zx[j] == 1) {
          pd = 1; pd_time = tstopx[j]; break;
        }
      }

      // switching and time to switching (if applicable)
      int xo = 0;
      double xo_time = NA_REAL;
      for (int j=k-J; j<k; j++) {
        if (Ax[j] == 1) {
          xo = 1; xo_time = tstopx[j]; break;
        }
      }
      
      for (int j=k-J; j<k; j++) {
        progressedx[j] = pd;
        timePDx[j] = pd_time;
        xox[j] = xo;
        xotimex[j] = xo_time;
      }
      
      // shift disease progression and alternative therapy status downward
      for (int j=k-J; j<k; j++) {
        Zx[j] = Zlagx[j];
        Ax[j] = Alagx[j];
        Alagx[j] = Alag2x[j];
      }
    }
    
    IntegerVector sub = seq(0, k-1);
    idx = idx[sub];
    arrivalTimex = arrivalTimex[sub];
    trtrandx = trtrandx[sub];
    bprogx = bprogx[sub];
    tpointx = tpointx[sub];
    Lx = Lx[sub];
    Llagx = Llagx[sub];
    Zx = Zx[sub];
    Ax = Ax[sub];
    Alagx = Alagx[sub];
    eventx = eventx[sub];
    diedx = diedx[sub];
    progressedx = progressedx[sub];
    xox = xox[sub];
    tstartx = tstartx[sub];
    tstopx = tstopx[sub];
    timeOSx = timeOSx[sub];
    timePDx = timePDx[sub];
    xotimex = xotimex[sub];
    censor_timex = censor_timex[sub];
    
    DataFrame data = DataFrame::create(
      _["id"] = idx,
      _["arrival_time"] = arrivalTimex,
      _["trtrand"] = trtrandx,
      _["bprog"] = bprogx,
      _["tpoint"] = tpointx,
      _["tstart"] = tstartx,
      _["tstop"] = tstopx,
      _["L"] = Lx,
      _["Llag"] = Llagx,
      _["Z"] = Zx,
      _["A"] = Ax,
      _["Alag"] = Alagx,
      _["event"] = eventx,
      _["timeOS"] = timeOSx,
      _["died"] = diedx,
      _["progressed"] = progressedx,
      _["timePD"] = timePDx,
      _["xo"] = xox,
      _["xotime"] = xotimex,
      _["censor_time"] = censor_timex
    );
    
    result[iter] = data;
  }
  
  return result;
}
