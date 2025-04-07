#' @title The Two-Stage Estimation (TSE) Method Using g-estimation  
#' for Treatment Switching
#' @description Obtains the causal parameter estimate of the logistic
#' regression switching model and the hazard ratio estimate of 
#' the Cox model to adjust for treatment switching.
#'
#' @param data The input data frame that contains the following variables:
#'
#'   * \code{id}: The id to identify observations belonging to the same
#'     subject for counting process data with time-dependent covariates.
#'
#'   * \code{stratum}: The stratum.
#'
#'   * \code{tstart}: The starting time of the time interval for
#'     counting-process data with time-dependent covariates.
#'
#'   * \code{tstop}: The stopping time of the time interval for
#'     counting-process data with time-dependent covariates.
#'
#'   * \code{event}: The event indicator, 1=event, 0=no event.
#'
#'   * \code{treat}: The randomized treatment indicator, 1=treatment,
#'     0=control.
#'
#'   * \code{censor_time}: The administrative censoring time. It should
#'     be provided for all subjects including those who had events.
#'
#'   * \code{pd}: The disease progression indicator, 1=PD, 0=no PD.
#'
#'   * \code{pd_time}: The time from randomization to PD.
#'
#'   * \code{swtrt}: The treatment switch indicator, 1=switch, 0=no switch.
#'
#'   * \code{swtrt_time}: The time from randomization to treatment switch.
#'   
#'   * \code{swtrt_time_upper}: The upper bound of treatment switching time.
#'
#'   * \code{base_cov}: The baseline covariates (excluding treat).
#'
#'   * \code{conf_cov}: The confounding variables for predicting
#'     treatment switching (excluding treat).
#'
#' @param id The name of the id variable in the input data.
#' @param stratum The name(s) of the stratum variable(s) in the input data.
#' @param tstart The name of the tstart variable in the input data.
#' @param tstop The name of the tstop variable in the input data.
#' @param event The name of the event variable in the input data.
#' @param treat The name of the treatment variable in the input data.
#' @param censor_time The name of the censor_time variable in the input data.
#' @param pd The name of the pd variable in the input data.
#' @param pd_time The name of the pd_time variable in the input data.
#' @param swtrt The name of the swtrt variable in the input data.
#' @param swtrt_time The name of the swtrt_time variable in the input data.
#' @param swtrt_time_upper The name of the swtrt_time_upper variable in the 
#'   input data.
#' @param base_cov The names of baseline covariates (excluding
#'   treat) in the input data for the Cox model.
#' @param conf_cov The names of confounding variables (excluding 
#'   treat) in the input data for the logistic regression switching model.
#' @param low_psi The lower limit of the causal parameter.
#' @param hi_psi The upper limit of the causal parameter.
#' @param n_eval_z The number of points between \code{low_psi} and 
#'   \code{hi_psi} (inclusive) at which to evaluate the Wald 
#'   statistics for the coefficient of the counterfactual in the logistic
#'   regression switching model.
#' @param strata_main_effect_only Whether to only include the strata main
#'   effects in the logistic regression switching model. Defaults to 
#'   \code{TRUE}, otherwise all possible strata combinations will be 
#'   considered in the switching model.
#' @param firth Whether the Firth's bias reducing penalized likelihood
#'   should be used. The default is \code{FALSE}.
#' @param flic Whether to apply intercept correction to obtain more
#'   accurate predicted probabilities. The default is \code{FALSE}.
#' @param recensor Whether to apply recensoring to counterfactual
#'   survival times. Defaults to \code{TRUE}.
#' @param admin_recensor_only Whether to apply recensoring to administrative
#'   censoring times only. Defaults to \code{TRUE}. If \code{FALSE},
#'   recensoring will be applied to the actual censoring times for dropouts.
#' @param swtrt_control_only Whether treatment switching occurred only in
#'   the control group. The default is \code{TRUE}.
#' @param alpha The significance level to calculate confidence intervals. 
#'   The default value is 0.05.
#' @param ties The method for handling ties in the Cox model, either
#'   "breslow" or "efron" (default).
#' @param tol The desired accuracy (convergence tolerance) for \code{psi}. 
#' @param offset The offset to calculate the time to event, PD, and 
#'   treatment switch. We can set \code{offset} equal to 1 (default), 
#'   1/30.4375, or 1/365.25 if the time unit is day, month, or year.
#' @param boot Whether to use bootstrap to obtain the confidence
#'   interval for hazard ratio. Defaults to \code{TRUE}.
#' @param n_boot The number of bootstrap samples.
#' @param seed The seed to reproduce the bootstrap results. The default is 
#'   missing, in which case, the seed from the environment will be used.
#'
#' @details We use the following steps to obtain the hazard ratio estimate
#' and confidence interval had there been no treatment switching:
#'
#' * Use a pooled logistic regression switching model to estimate 
#'   the causal parameter \eqn{\psi} based on the patients in the 
#'   control group  who had disease progression:
#'   \deqn{\textrm{logit}(p(E_{ik})) = \alpha U_{i,\psi} + 
#'   \sum_{j} \beta_j x_{ijk}}
#'   where \eqn{E_{ik}} is the observed switch indicator for individual 
#'   \eqn{i} at observation \eqn{k}, 
#'   \deqn{U_{i,\psi} = T_{C_i} + e^{\psi}T_{E_i}}
#'   is the counterfactual survival time for individual \eqn{i} given a 
#'   specific value for \eqn{\psi}, and \eqn{x_{ijk}} are the confounders
#'   for individual \eqn{i} at observation \eqn{k}. 
#'   When applied from a secondary baseline, \eqn{U_{i,\psi}} 
#'   refers to post-secondary baseline counterfactual survival, where 
#'   \eqn{T_{C_i}} corresponds to the time spent after the secondary baseline 
#'   on control treatment, and \eqn{T_{E_i}} corresponds to the time spent 
#'   after the secondary baseline on the experimental treatment.
#'   
#' * Search for \eqn{\psi} such that the estimate of \eqn{\alpha} is close
#'   to zero. This will be the estimate of the caual parameter. The 
#'   confidence interval for \eqn{\psi} can be obtained as the value of 
#'   \eqn{\psi} such that the corresponding two-sided p-value for 
#'   testing \eqn{H_0:\alpha = 0} in the switching model is equal to the 
#'   nominal significance level.   
#'
#' * Derive the counterfactual survival times for control patients
#'   had there been no treatment switching.
#'
#' * Fit the Cox proportional hazards model to the observed survival times
#'   for the experimental group and the counterfactual survival times 
#'   for the control group to obtain the hazard ratio estimate.
#'   
#' * If bootstrapping is used, the confidence interval and corresponding 
#'   p-value for hazard ratio are calculated based on a t-distribution with 
#'   \code{n_boot - 1} degrees of freedom. 
#'
#' @return A list with the following components:
#'
#' * \code{psi}: The estimated causal parameter for the control group.
#'
#' * \code{psi_CI}: The confidence interval for \code{psi}.
#'
#' * \code{psi_CI_type}: The type of confidence interval for \code{psi},
#'   i.e., "logistic model" or "bootstrap".
#'   
#' * \code{logrank_pvalue}: The two-sided p-value of the log-rank test
#'   for the ITT analysis.
#'
#' * \code{cox_pvalue}: The two-sided p-value for treatment effect based on
#'   the Cox model applied to counterfactual unswitched survival times.
#'   If \code{boot} is \code{TRUE}, this value represents the 
#'   bootstrap p-value.
#'
#' * \code{hr}: The estimated hazard ratio from the Cox model.
#'
#' * \code{hr_CI}: The confidence interval for hazard ratio.
#'
#' * \code{hr_CI_type}: The type of confidence interval for hazard ratio,
#'   either "Cox model" or "bootstrap".
#'   
#' * \code{analysis_switch}: A list of data and analysis results related to  
#'   treatment switching. 
#'   
#'     - \code{data_switch}: The list of input data for the time from 
#'       secondary baseline to switch by treatment group. The variables 
#'       include \code{id}, \code{stratum} (if applicable), \code{swtrt}, 
#'       and \code{swtrt_time}. If \code{swtrt == 0}, then \code{swtrt_time} 
#'       is censored at the time from secondary baseline to either 
#'       death or censoring.
#'   
#'     - \code{km_switch}: The list of Kaplan-Meier plot data for the 
#'       time from secondary baseline to switch by treatment group.
#'   
#'     - \code{eval_z}: The list of data by treatment group containing 
#'       the Wald statistics for the coefficient of the counterfactual 
#'       in the logistic regression switching model, evaluated at 
#'       a sequence of \code{psi} values. Used to plot and check 
#'       if the range of \code{psi} values to search for the solution 
#'       and limits of confidence interval of \code{psi} need be modified.
#'   
#'     - \code{data_nullcox}: The list of input data for counterfactual 
#'       survival times for the null Cox model by treatment group.
#'   
#'     - \code{fit_nullcox}: The list of fitted null Cox models for 
#'       counterfactual survival times by treatment group, which contains
#'       the martingale residuals.
#' 
#'     - \code{data_logis}: The list of input data for pooled logistic 
#'       regression models for treatment switching using g-estimation.
#'   
#'     - \code{fit_logis}: The list of fitted pooled logistic regression 
#'       models for treatment switching using g-estimation.
#'   
#' * \code{data_outcome}: The input data for the outcome Cox model.
#'
#' * \code{fit_outcome}: The fitted outcome Cox model.
#'
#' * \code{settings}: A list with the following components:
#'
#'     - \code{low_psi}: The lower limit of the causal parameter.
#'     
#'     - \code{hi_psi}: The upper limit of the causal parameter.
#'     
#'     - \code{n_eval_z}: The number of points between \code{low_psi} and 
#'       \code{hi_psi} (inclusive) at which to evaluate the Wald statistics 
#'       for the coefficient for the counterfactual in the logistic 
#'       regression switching model.
#'
#'     - \code{strata_main_effect_only}: Whether to only include the strata
#'       main effects in the logistic regression switching model.
#'       
#'     - \code{firth}: Whether the Firth's penalized likelihood is used.
#'
#'     - \code{flic}: Whether to apply intercept correction.
#'
#'     - \code{recensor}: Whether to apply recensoring to counterfactual
#'       survival times.
#'
#'     - \code{admin_recensor_only}: Whether to apply recensoring to
#'       administrative censoring times only.
#'
#'     - \code{swtrt_control_only}: Whether treatment switching occurred
#'       only in the control group.
#'
#'     - \code{alpha}: The significance level to calculate confidence
#'       intervals.
#'
#'     - \code{ties}: The method for handling ties in the Cox model.
#'     
#'     - \code{tol}: The desired accuracy (convergence tolerance) 
#'       for \code{psi}.
#'
#'     - \code{offset}: The offset to calculate the time to event, PD, and
#'       treatment switch.
#'
#'     - \code{boot}: Whether to use bootstrap to obtain the confidence
#'       interval for hazard ratio.
#'
#'     - \code{n_boot}: The number of bootstrap samples.
#'
#'     - \code{seed}: The seed to reproduce the bootstrap results.
#'
#' * \code{psi_trt}: The estimated causal parameter for the experimental 
#'   group if \code{swtrt_control_only} is \code{FALSE}.
#'
#' * \code{psi_trt_CI}: The confidence interval for \code{psi_trt} if
#'   \code{swtrt_control_only} is \code{FALSE}.
#'
#' * \code{hr_boots}: The bootstrap hazard ratio estimates if \code{boot} is
#'   \code{TRUE}.
#'
#' * \code{psi_boots}: The bootstrap \code{psi} estimates if \code{boot} is 
#'   \code{TRUE}.
#' 
#' * \code{psi_trt_boots}: The bootstrap \code{psi_trt} estimates if 
#'   \code{boot} is \code{TRUE} and \code{swtrt_control_only} is 
#'   \code{FALSE}.
#' 
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @references
#' NR Latimer, IR White, K Tilling, and U Siebert.
#' Improved two-stage estimation to adjust for treatment switching in 
#' randomised trials: g-estimation to address time-dependent confounding.
#' Statistical Methods in Medical Research. 2020;29(10):2900-2918.
#'
#' @examples
#'
#' # Example 1: one-way treatment switching (control to active)
#' 
#' sim1 <- tsegestsim(
#'   n = 500, allocation1 = 2, allocation2 = 1, pbprog = 0.5, 
#'   trtlghr = -0.5, bprogsl = 0.3, shape1 = 1.8, 
#'   scale1 = 0.000025, shape2 = 1.7, scale2 = 0.000015, 
#'   pmix = 0.5, admin = 5000, pcatnotrtbprog = 0.5, 
#'   pcattrtbprog = 0.25, pcatnotrt = 0.2, pcattrt = 0.1, 
#'   catmult = 0.5, tdxo = 1, ppoor = 0.1, pgood = 0.04, 
#'   ppoormet = 0.4, pgoodmet = 0.2, xomult = 1.4188308, 
#'   milestone = 546, swtrt_control_only = TRUE,
#'   outputRawDataset = 1, seed = 2000)
#'   
#' fit1 <- tsegest(
#'   data = sim1$paneldata, id = "id", 
#'   tstart = "tstart", tstop = "tstop", event = "died", 
#'   treat = "trtrand", censor_time = "censor_time", 
#'   pd = "progressed", pd_time = "timePFSobs", swtrt = "xo", 
#'   swtrt_time = "xotime", swtrt_time_upper = "xotime_upper",
#'   base_cov = "bprog", conf_cov = "bprog*catlag", 
#'   low_psi = -3, hi_psi = 3, strata_main_effect_only = TRUE,
#'   recensor = TRUE, admin_recensor_only = TRUE, 
#'   swtrt_control_only = TRUE, alpha = 0.05, ties = "efron", 
#'   tol = 1.0e-6, boot = FALSE)
#'   
#' c(fit1$hr, fit1$hr_CI)
#'
#' # Example 2: two-way treatment switching
#'
#' sim2 <- tsegestsim(
#'   n = 500, allocation1 = 2, allocation2 = 1, pbprog = 0.5, 
#'   trtlghr = -0.5, bprogsl = 0.3, shape1 = 1.8, 
#'   scale1 = 0.000025, shape2 = 1.7, scale2 = 0.000015, 
#'   pmix = 0.5, admin = 5000, pcatnotrtbprog = 0.5, 
#'   pcattrtbprog = 0.25, pcatnotrt = 0.2, pcattrt = 0.1, 
#'   catmult = 0.5, tdxo = 1, ppoor = 0.1, pgood = 0.04, 
#'   ppoormet = 0.4, pgoodmet = 0.2, xomult = 1.4188308, 
#'   milestone = 546, swtrt_control_only = FALSE,
#'   outputRawDataset = 1, seed = 2000)
#'   
#' fit2 <- tsegest(
#'   data = sim2$paneldata, id = "id", 
#'   tstart = "tstart", tstop = "tstop", event = "died", 
#'   treat = "trtrand", censor_time = "censor_time", 
#'   pd = "progressed", pd_time = "timePFSobs", swtrt = "xo", 
#'   swtrt_time = "xotime", swtrt_time_upper = "xotime_upper",
#'   base_cov = "bprog", conf_cov = "bprog*catlag", 
#'   low_psi = -3, hi_psi = 3, strata_main_effect_only = TRUE,
#'   recensor = TRUE, admin_recensor_only = TRUE, 
#'   swtrt_control_only = FALSE, alpha = 0.05, ties = "efron", 
#'   tol = 1.0e-6, boot = FALSE)
#'   
#' c(fit2$hr, fit2$hr_CI)
#' 
#' @export
tsegest <- function(data, id = "id", stratum = "", 
                    tstart = "tstart", tstop = "tstop", event = "event",
                    treat = "treat", censor_time = "censor_time",
                    pd = "pd", pd_time = "pd_time",
                    swtrt = "swtrt", swtrt_time = "swtrt_time",
                    swtrt_time_upper = "", base_cov = "", conf_cov = "",
                    low_psi = -2, hi_psi = 2, n_eval_z = 101,
                    strata_main_effect_only = TRUE,
                    firth = FALSE, flic = FALSE,
                    recensor = TRUE, admin_recensor_only = TRUE,
                    swtrt_control_only = TRUE, alpha = 0.05, 
                    ties = "efron", tol = 1.0e-6, offset = 1, 
                    boot = TRUE, n_boot = 1000, seed = NA) {

  rownames(data) = NULL

  elements = c(id, stratum, tstart, tstop, event, treat, censor_time, 
               pd, swtrt)
  elements = unique(elements[elements != "" & elements != "none"])
  mf = model.frame(formula(paste("~", paste(elements, collapse = "+"))),
                   data = data)

  rownum = as.integer(rownames(mf))
  df = data[rownum,]

  nvar = length(base_cov)
  if (missing(base_cov) || is.null(base_cov) || (nvar == 1 && (
    base_cov[1] == "" || tolower(base_cov[1]) == "none"))) {
    p = 0
  } else {
    t1 = terms(formula(paste("~", paste(base_cov, collapse = "+"))))
    t2 = attr(t1, "factors")
    t3 = rownames(t2)
    p = length(t3)
  }

  if (p >= 1) {
    mm = model.matrix(t1, df)
    colnames(mm) = make.names(colnames(mm))
    varnames = colnames(mm)[-1]
    for (i in 1:length(varnames)) {
      if (!(varnames[i] %in% names(df))) {
        df[,varnames[i]] = mm[,varnames[i]]
      }
    }
  } else {
    varnames = ""
  }

  nvar2 = length(conf_cov)
  if (missing(conf_cov) || is.null(conf_cov) || (nvar2 == 1 && (
    conf_cov[1] == "" || tolower(conf_cov[1]) == "none"))) {
    p2 = 0
  } else {
    t1 = terms(formula(paste("~", paste(conf_cov, collapse = "+"))))
    t2 = attr(t1, "factors")
    t3 = rownames(t2)
    p2 = length(t3)
  }

  if (p2 >= 1) {
    mm2 = model.matrix(t1, df)
    colnames(mm2) = make.names(colnames(mm2))
    varnames2 = colnames(mm2)[-1]
    for (i in 1:length(varnames2)) {
      if (!(varnames2[i] %in% names(df))) {
        df[,varnames2[i]] = mm2[,varnames2[i]]
      }
    }
  } else {
    varnames2 = ""
  }
  
  if (missing(swtrt_time_upper) || is.null(swtrt_time_upper) || (
    swtrt_time_upper[1] == "" || tolower(swtrt_time_upper[1]) == "none")) {
    swtrt_time_upper = "swtrt_time_upper";
    df$swtrt_time_upper = 1.0e8;
  }

  out <- tsegestcpp(
    data = df, id = id, stratum = stratum, 
    tstart = tstart, tstop = tstop, event = event,
    treat = treat, censor_time = censor_time, pd = pd,
    pd_time = pd_time, swtrt = swtrt, swtrt_time = swtrt_time,
    swtrt_time_upper = swtrt_time_upper,
    base_cov = varnames, conf_cov = varnames2,
    low_psi = low_psi, hi_psi = hi_psi, n_eval_z = n_eval_z,
    strata_main_effect_only = strata_main_effect_only,
    firth = firth, flic = flic,
    recensor = recensor, admin_recensor_only = admin_recensor_only,
    swtrt_control_only = swtrt_control_only, alpha = alpha,
    ties = ties, tol = tol, offset = offset, 
    boot = boot, n_boot = n_boot, seed = seed)
  
  K = ifelse(swtrt_control_only, 1, 2)
  for (h in 1:K) {
    out$analysis_switch$data_logis[[h]]$data$uid <- NULL
    out$analysis_switch$data_nullcox[[h]]$data$ustratum <- NULL
  }
  
  out$data_outcome$uid <- NULL
  out$data_outcome$ustratum <- NULL
  
  df[, "tstart"] = df[, tstart]
  df[, "tstop"] = df[, tstop]
  
  if (p >= 1) {
    t1 = terms(formula(paste("~", paste(base_cov, collapse = "+"))))
    t2 = attr(t1, "factors")
    t3 = rownames(t2)
    
    add_vars <- setdiff(t3, varnames)
    if (length(add_vars) > 0) {
      out$data_outcome <- merge(out$data_outcome, df[, c(id, add_vars)], 
                                by = id, all.x = TRUE, sort = FALSE)
    }
    
    del_vars <- setdiff(varnames, t3)
    if (length(del_vars) > 0) {
      out$data_outcome[, del_vars] <- NULL
    }
  }
  
  if (p2 >= 1) {
    t1 = terms(formula(paste("~", paste(conf_cov, collapse = "+"))))
    t2 = attr(t1, "factors")
    t3 = rownames(t2)
    
    tem_vars <- c(pd_time, swtrt, swtrt_time, swtrt_time_upper)
    add_vars <- c(setdiff(t3, varnames2), tem_vars)
    if (length(add_vars) > 0) {
      for (h in 1:K) {
        out$analysis_switch$data_logis[[h]]$data <- 
          merge(out$analysis_switch$data_logis[[h]]$data, 
                df[, c(id, "tstart", "tstop", add_vars)], 
                by = c(id, "tstart", "tstop"), all.x = TRUE, sort = FALSE)
      }
    }
    
    del_vars <- setdiff(varnames2, t3)
    if (length(del_vars) > 0) {
      for (h in 1:K) {
        out$analysis_switch$data_logis[[h]]$data[, del_vars] <- NULL
      }
    }
    
    for (h in 1:K) {
      out$analysis_switch$data_logis[[h]]$data <- 
        out$analysis_switch$data_logis[[h]]$data[
          , !startsWith(names(out$analysis_switch$data_logis[[h]]$data), 
                        "stratum_")]
    }
  }
  
  out
}
