#' @title Simulate Data for Treatment Switching
#' @description Simulates data for studies involving treatment switching, 
#' incorporating time-dependent confounding. The generated data can be used 
#' to evaluate methods for handling treatment switching in survival 
#' analysis.
#'
#' @param tdxo Logical indicator for timing of treatment switching:
#'   \itemize{
#'     \item 1: Treatment switching can occur at or after disease 
#'              progression.
#'     \item 0: Treatment switching is restricted to the time of disease 
#'              progression.
#'   }
#' @param coxo Logical indicator for arm-specific treatment switching:
#'   \itemize{
#'     \item 1: Treatment switching occurs only in the control arm.
#'     \item 0: Treatment switching is allowed in both arms.
#'   }
#' @param allocation1 Number of subjects in the active treatment group in
#'   a randomization block. Defaults to 1 for equal randomization.
#' @param allocation2 Number of subjects in the control group in
#'   a randomization block. Defaults to 1 for equal randomization.
#' @param p_X_1 Probability of poor baseline prognosis in the experimental 
#'   arm.
#' @param p_X_0 Probability of poor baseline prognosis in the control arm.
#' @param rate_T Baseline hazard rate for time to death.
#' @param beta1 Log hazard ratio for randomized treatment (\code{R}).
#' @param beta2 Log hazard ratio for baseline covariate (\code{X}).
#' @param gamma0 Intercept for the time-dependent covariate model 
#'   (\code{L}).
#' @param gamma1 Coefficient for lagged treatment switching (\code{Alag}) 
#'   in the \code{L} model.
#' @param gamma2 Coefficient for lagged \code{L} (\code{Llag}) in the 
#'   \code{L} model.
#' @param gamma3 Coefficient for baseline covariate (\code{X}) in the 
#'   \code{L} model.
#' @param gamma4 Coefficient for randomized treatment (\code{R}) in the 
#'   \code{L} model.
#' @param zeta0 Intercept for the disease progression model (\code{Z}).
#' @param zeta1 Coefficient for \code{L} in the \code{Z} model.
#' @param zeta2 Coefficient for baseline covariate (\code{X}) in the 
#'   \code{Z} model.
#' @param zeta3 Coefficient for randomized treatment (\code{R}) in the 
#'   \code{Z} model.
#' @param alpha0 Intercept for the treatment switching model (\code{A}).
#' @param alpha1 Coefficient for \code{L} in the \code{A} model.
#' @param alpha2 Coefficient for baseline covariate (\code{X}) in the 
#'   \code{A} model.
#' @param theta1_1 Negative log time ratio for \code{A} for the 
#'   experimental arm.
#' @param theta1_0 Negative log time ratio for \code{A} for the control arm.
#' @param theta2 Negative log time ratio for \code{L}.
#' @param rate_C Hazard rate for random (dropout) censoring.
#' @param accrualTime A vector that specifies the starting time of 
#'   piecewise Poisson enrollment time intervals. Must start with 0, 
#'   e.g., \code{c(0, 3)} breaks the time axis into 2 accrual intervals: 
#'   [0, 3) and [3, Inf).
#' @param accrualIntensity A vector of accrual intensities. One for 
#'   each accrual time interval.
#' @param followupTime	Follow-up time for a fixed follow-up design.
#' @param fixedFollowup Whether a fixed follow-up design is used. 
#'   Defaults to 0 for variable follow-up.
#' @param plannedTime The calendar time for the analysis.
#' @param days Number of days in each treatment cycle.
#' @param n Number of subjects per simulation.
#' @param NSim Number of simulated datasets.
#' @param seed Random seed for reproducibility.
#'
#' @details 
#' The simulation algorithm is adapted from Xu et al. (2022) and simulates
#' disease progression status while incorporating the multiplicative 
#' effects of both baseline and time-dependent covariates on survival time. 
#' The design options \code{tdxo} and \code{coxo} specify
#' the timing of treatment switching and the study arm eligibility 
#' for switching, respectively. Time is measured in days. 
#' 
#' In a fixed follow-up design, all subjects share the same follow-up
#' duration. In contrast, under a variable follow-up design, follow-up
#' time also depends on each subject's enrollment date.  
#' The number of treatment cycles for a subject is determined by
#' dividing the follow-up time by the number of days in each cycle. 
#' \enumerate{
#'   \item At randomization, subjects are assigned to treatment arms 
#'   using block randomization, with \code{allocation1} patients 
#'   assigned to active treatment and \code{allocation2} to control
#'   within each randomized block. A baseline covariate is 
#'   also generated for each subject:
#'   \deqn{X_i \sim \mbox{Bernoulli}(p_1 R_i + p_0 (1-R_i))}
#'         
#'   \item The initial survival time is drawn
#'   from an exponential distribution with hazard:
#'   \deqn{\lambda_T \exp(\beta_1 R_i + \beta_2 X_i)}
#'   We define the event indicator at cycle \eqn{j} as
#'   \deqn{Y_{i,j} = I(T_i \leq j\times days)}
#'         
#'   \item The initial states are set to
#'   \eqn{L_{i,0} = 0}, \eqn{Z_{i,0} = 0}, \eqn{A_{i,0} = 0},
#'   \eqn{Y_{i,0} = 0}. For each treatment cycle \eqn{j=1,\ldots,J},
#'   we set \eqn{tstart = (j-1) \times days}.
#'         
#'   \item Generate time-dependent covariates:
#'   \deqn{\mbox{logit} P(L_{i,j}=1|\mbox{history}) = 
#'   \gamma_0 + \gamma_1 A_{i,j-1} + \gamma_2 L_{i,j-1} + 
#'   \gamma_3 X_i + \gamma_4 R_i}
#'      
#'   \item If \eqn{T_i \leq j \times days}, set \eqn{tstop = T_i} and 
#'   \eqn{Y_{i,j} = 1}, which completes data generation
#'   for subject \eqn{i}.
#'           
#'   \item If \eqn{T_i > j \times days}, set \eqn{tstop = j\times days},  
#'   \eqn{Y_{i,j} = 0}, and perform the following before proceeding to 
#'   the next cycle for the subject.
#'   
#'   \item Generate disease progression status: 
#'   If \eqn{Z_{i,j-1} = 0},  
#'   \deqn{\mbox{logit} P(Z_{i,j}=1 | \mbox{history}) = \zeta_0 + 
#'   \zeta_1 L_{i,j} + \zeta_2 X_i + \zeta_3 R_i}
#'   Otherwise, set \eqn{Z_{i,j} = 1}.
#'     
#'   \item Generate alternative therapy status:     
#'   If \eqn{A_{i,j-1} = 0}, determine switching eligibility 
#'   based on design options.        
#'   If switching is allowed:
#'   \deqn{\mbox{logit} P(A_{i,j} = 1 | \mbox{history}) = \alpha_0 + 
#'   \alpha_1 L_{i,j} + \alpha_2 X_i}       
#'   If switching is now allowed, set \eqn{A_{i,j} = 0}.      
#'   If \eqn{A_{i,j-1} = 1}, set \eqn{A_{i,j} = 1} (once switched to 
#'   alternative therapy, remain on alternative therapy).
#'          
#'   \item Update survival time based on changes in alternative 
#'   therapy status and time-dependent covariates:
#'   \deqn{T_i = j\times days + (T_i - j\times days) \exp\{
#'   -(\theta_{1,1}R_i + \theta_{1,0}(1-R_i))(A_{i,j} - A_{i,j-1}) 
#'   -\theta_2 (L_{i,j} - L_{i,j-1})\}}   
#' }
#' 
#' Additional random censoring times are generated from an exponential 
#' distribution with hazard rate \eqn{\lambda_C}.
#' 
#' An extra record is generated when the minimum of the latent survival 
#' time, the random censoring time, and the administrative censoring time 
#' is greater than the number of regular treatment cycles times 
#' days per cycle.
#' 
#' Finally we apply the lag function so that \eqn{Z_{i,j}} and 
#' \eqn{A_{i,j}} represent the PD status and alternative therapy status 
#' at the start of cycle \eqn{j} (and thus remain appplicable for the 
#' entire cycle \eqn{j}) for subject \eqn{i}.
#' 
#' To estimate the true treatment effect in a Cox marginal 
#' structural model, one can set \eqn{\alpha_0 = -\infty}, which 
#' effectively forces \eqn{A_{i,j} = 0} (disabling treatment switching). 
#' The coefficient for the randomized treatment can then be estimated 
#' using a Cox proportional hazards model.
#'
#' @return
#' A list of data frames, each containing simulated longitudinal covariate, 
#' pd status, alternative therapy status, and event history data with the 
#' following variables:
#'
#'  * \code{id}: Subject identifier.
#'  
#'  * \code{arrival_time}: The enrollment time for the subject.
#'  
#'  * \code{trtrand}: Randomized treatment assignment (0 = control, 
#'    1 = experimental)
#'
#'  * \code{bprog}: Baseline prognosis (0 = good, 1 = poor).
#'
#'  * \code{tpoint}: Treatment cycle index.
#'
#'  * \code{tstart}: Start day of the treatment cycle.
#'
#'  * \code{tstop}: End day of the treatment cycle.
#'
#'  * \code{L}: Time-dependent covariate at \code{tstart} predicting 
#'    survival and switching; affected by treatment switching.
#'
#'  * \code{Llag}: Lagged value of \code{L}.
#'
#'  * \code{Z}: Disease progression status at \code{tstart}.
#'
#'  * \code{A}: Treatment switching status at \code{tstart}.
#'
#'  * \code{Alag}: Lagged value of \code{A}.
#'
#'  * \code{event}: Death indicator at \code{tstop}.
#'
#'  * \code{timeOS}: Observed time to death or censoring.
#'     
#'  * \code{died}: Indicator of death by end of follow-up.
#'     
#'  * \code{progressed}: Indicator of disease progression by end of 
#'    follow-up.
#'       
#'  * \code{timePD}: Observed time to progression or censoring.
#'       
#'  * \code{xo}: Indicator for whether treatment switching occurred.
#'     
#'  * \code{xotime}: Time of treatment switching (if applicable).
#'     
#'  * \code{censor_time}: Administrative censoring time.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @references
#' 
#' Jessica G. Young, and Eric J. Tchetgen Tchetgen. 
#' Simulation from a known Cox MSM using standard parametric models 
#' for the g-formula.
#' Statistics in Medicine. 2014;33(6):1001-1014. 
#' 
#' NR Latimer, IR White, K Tilling, and U Siebert.
#' Improved two-stage estimation to adjust for treatment switching in 
#' randomised trials: g-estimation to address time-dependent confounding.
#' Statistical Methods in Medical Research. 2020;29(10):2900-2918.
#' 
#' Jing Xu, Guohui Liu, and Bingxia Wang.
#' Bias and type I error control in correcting treatment effect
#' for treatment switching using marginal structural models 
#' in Phse III oncology trials.
#' Journal of Biopharmaceutical Statistics. 2022;32(6):897-914.
#'
#' @examples
#'
#' library(dplyr)
#' 
#' simulated.data <- tssim(
#'   tdxo = 1, coxo = 1, allocation1 = 1, allocation2 = 1,
#'   p_X_1 = 0.3, p_X_0 = 0.3, 
#'   rate_T = 0.002, beta1 = -0.5, beta2 = 0.3, 
#'   gamma0 = 0.3, gamma1 = -0.9, gamma2 = 0.7, gamma3 = 1.1, gamma4 = -0.8,
#'   zeta0 = -3.5, zeta1 = 0.5, zeta2 = 0.2, zeta3 = -0.4, 
#'   alpha0 = 0.5, alpha1 = 0.5, alpha2 = 0.4, 
#'   theta1_1 = -0.4, theta1_0 = -0.4, theta2 = 0.2,
#'   rate_C = 0.0000855, accrualIntensity = 20/30, 
#'   fixedFollowup = FALSE, plannedTime = 1350, days = 30,
#'   n = 500, NSim = 100, seed = 314159)
#'   
#' simulated.data[[1]] %>% filter(id == 1)
#'
#' @export
tssim <- function(
    tdxo = FALSE, coxo = TRUE, allocation1 = 1, allocation2 = 1,
    p_X_1, p_X_0, rate_T, beta1, beta2, gamma0, gamma1, 
    gamma2, gamma3, gamma4, zeta0, zeta1, zeta2, zeta3, 
    alpha0, alpha1, alpha2, theta1_1, theta1_0, theta2, 
    rate_C, accrualTime = 0, accrualIntensity = NA, followupTime = NA, 
    fixedFollowup = FALSE, plannedTime = NA, days = NA, 
    n = NA, NSim = 1000, seed = 0) {
  
  tssimcpp(tdxo, coxo, allocation1, allocation2,
           p_X_1, p_X_0, rate_T, beta1, beta2, gamma0, gamma1,
           gamma2, gamma3, gamma4, zeta0, zeta1, zeta2, zeta3,
           alpha0, alpha1, alpha2, theta1_1, theta1_0, theta2,
           rate_C, accrualTime, accrualIntensity, followupTime,
           fixedFollowup, plannedTime, days, n, NSim, seed)
}
