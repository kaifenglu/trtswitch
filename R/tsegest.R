#' @title Two-Stage Estimation with g-Estimation (TSEgest) for Treatment 
#' Switching
#' @description Estimates the causal parameter using g-estimation by 
#' fitting a pooled logistic regression switching model that includes 
#' counterfactual \emph{unswitched} survival times and time-dependent 
#' confounders as covariates. The adjusted hazard ratio is then obtained 
#' from the Cox model using counterfactual \emph{unswitched} survival 
#' times based on the estimated causal parameter.
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
#'   * \code{pd_time}: The time from randomization to disease progression.
#'
#'   * \code{swtrt}: The treatment switch indicator, 1=switch, 0=no switch.
#'
#'   * \code{swtrt_time}: The time from randomization to treatment switch.
#'   
#'   * \code{base_cov}: The baseline covariates (excluding treat).
#'
#'   * \code{conf_cov}: The confounding variables (excluding treat) for 
#'     predicting treatment switching.
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
#' @param base_cov The names of baseline covariates (excluding
#'   treat) in the input data for the Cox model.
#' @param conf_cov The names of confounding variables (excluding 
#'   treat) in the input data for the logistic regression switching model.
#' @param strata_main_effect_only Whether to only include the strata main
#'   effects in the logistic regression switching model. Defaults to 
#'   \code{TRUE}, otherwise all possible strata combinations will be 
#'   considered in the switching model.
#' @param ns_df Degrees of freedom for the natural cubic spline for 
#'   visit-specific intercepts of the pooled logistic regression model. 
#'   Defaults to 3 for two internal knots at the 33 and 67 percentiles
#'   of the treatment switching times.
#' @param firth Whether the Firth's bias reducing penalized likelihood
#'   should be used.
#' @param flic Whether to apply intercept correction to obtain more
#'   accurate predicted probabilities.
#' @param low_psi The lower limit of the causal parameter.
#' @param hi_psi The upper limit of the causal parameter.
#' @param n_eval_z The number of points between \code{low_psi} and 
#'   \code{hi_psi} (inclusive) at which to evaluate the Wald 
#'   statistics for the coefficient of the counterfactual in the logistic
#'   regression switching model.
#' @param recensor Whether to apply recensoring to counterfactual
#'   survival times. Defaults to \code{TRUE}.
#' @param admin_recensor_only Whether to apply recensoring to administrative
#'   censoring times only. Defaults to \code{TRUE}. If \code{FALSE},
#'   recensoring will be applied to the actual censoring times for dropouts.
#' @param swtrt_control_only Whether treatment switching occurred only in
#'   the control group. The default is \code{TRUE}.
#' @param gridsearch Whether to use grid search to estimate the causal
#'   parameter \code{psi}. Defaults to \code{TRUE}, otherwise, a root
#'   finding algorithm will be used.
#' @param root_finding Character string specifying the univariate 
#'   root-finding algorithm to use. Options are \code{"brent"} (default)
#'   for Brent's method, or \code{"bisection"} for the bisection method.
#' @param alpha The significance level to calculate confidence intervals. 
#' @param ties The method for handling ties in the Cox model, either
#'   "breslow" or "efron" (default).
#' @param tol The desired accuracy (convergence tolerance) for \code{psi}
#'    for the root finding algorithm. 
#' @param offset The offset to calculate the time from disease progression 
#'   to death or censoring, the time from disease progression to treatment 
#'   switch, and the time from treatment switch to death or censoring. 
#'   We can set \code{offset} equal to 0 (no offset), and 1 (default), 
#'   1/30.4375, or 1/365.25 if the time unit is day, month, or year, 
#'   respectively. 
#' @param boot Whether to use bootstrap to obtain the confidence
#'   interval for hazard ratio. Defaults to \code{TRUE}.
#' @param n_boot The number of bootstrap samples.
#' @param seed The seed to reproduce the bootstrap results. The default is 
#'   `NA`, in which case, the seed from the environment will be used.
#'
#' @details Assuming one-way switching from control to treatment, the 
#' hazard ratio and confidence interval under a no-switching scenario 
#' are obtained as follows:
#'
#' * Fit a pooled logistic regression switching model among control-arm 
#'   patients who experienced disease progression:
#'   \deqn{\text{logit}(p(E_{ik})) = \alpha U_{i,\psi} + \sum_{j} \beta_j 
#'   x_{ijk}} 
#'   where \eqn{E_{ik}} is the switch indicator for subject \eqn{i} at 
#'   observation \eqn{k}, 
#'   \deqn{U_{i,\psi} = T_{C_i} + e^{\psi}T_{E_i}} is the counterfactual 
#'   survival time given a specific \eqn{\psi}, and \eqn{x_{ijk}} 
#'   represents the time-dependent confounders. 
#'   Natural cubic splines of time can be included to model time-varying 
#'   baseline hazards. \eqn{U_{i,\psi}} is defined relative to the 
#'   secondary baseline at disease progression and represents 
#'   post-progression counterfactual survival, where \eqn{T_{C_i}} and 
#'   \eqn{T_{E_i}} correspond to  time spent after progression on control 
#'   and experimental treatments, respectively. 
#'   Martingale residuals may be used in place of counterfactual survival 
#'   times to account for censoring.
#'
#' * Identify the value of \eqn{\psi} for which the Z-statistic of 
#'   \eqn{\alpha} is approximately zero. This value is the causal 
#'   parameter estimate.
#'
#' * Compute counterfactual survival times for control patients using 
#'   the estimated  \eqn{\psi}.
#'
#' * Fit a Cox model to the observed survival times for the treatment group 
#'   and the counterfactual survival times for the control group to 
#'   estimate the hazard ratio.
#'
#' * When bootstrapping is used, derive the confidence interval and 
#'   p-value for the hazard ratio from a t-distribution with 
#'   \code{n_boot - 1} degrees of freedom.
#'
#' If treatment switching occurs before or in the absence of recorded disease 
#' progression, the patient is considered to have progressed at the time of 
#' treatment switching. 
#' 
#' If grid search is used to estimate \eqn{\psi}, the estimated \eqn{\psi} 
#' is the one with the smallest absolute value among those at which 
#' the Z-statistic is zero based on linear interpolation. 
#' If root finding is used, the estimated \eqn{\psi} is
#' the solution to the equation where the Z-statistic is zero.
#'
#' @return A list with the following components:
#'
#' * \code{psi}: The estimated causal parameter for the control group.
#' 
#' * \code{psi_roots}: Vector of \code{psi} values for the control group 
#'   at which the Z-statistic is zero, identified using grid search and 
#'   linear interpolation.
#'
#' * \code{psi_CI}: The confidence interval for \code{psi}.
#'
#' * \code{psi_CI_type}: The type of confidence interval for \code{psi},
#'   i.e., "grid search", "root finding", or "bootstrap".
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
#' * \code{event_summary}: A data frame containing the count and percentage
#'   of deaths, disease progressions, and switches by treatment arm.
#'   
#' * \code{data_switch}: The list of input data for the time from 
#'   disease progression to switch by treatment group. The variables 
#'   include \code{id}, \code{stratum}, \code{"swtrt"}, 
#'   and \code{"swtrt_time"}. If \code{swtrt == 0}, then \code{swtrt_time} 
#'   is censored at the time from disease progression to death or censoring.
#' 
#' * \code{km_switch}: The list of Kaplan-Meier plot data for the 
#'   time from disease progression to switch by treatment group.
#' 
#' * \code{eval_z}: The list of data by treatment group containing 
#'   the Wald statistics for the coefficient of the counterfactual 
#'   in the logistic regression switching model, evaluated at 
#'   a sequence of \code{psi} values. Used to plot and check 
#'   if the range of \code{psi} values to search for the solution 
#'   and limits of confidence interval of \code{psi} need be modified.
#' 
#' * \code{data_nullcox}: The list of input data for counterfactual 
#'   survival times for the null Cox model by treatment group.
#'   The variables include \code{id}, \code{stratum},
#'   \code{"t_star"} and \code{"d_star"}.
#' 
#' * \code{fit_nullcox}: The list of fitted null Cox models for 
#'   counterfactual survival times by treatment group, which contains
#'   the martingale residuals.
#' 
#' * \code{data_logis}: The list of input data for pooled logistic 
#'   regression models for treatment switching using g-estimation.
#'   The variables include \code{id}, \code{stratum}, 
#'   \code{"tstart"}, \code{"tstop"}, \code{"cross"}, 
#'   \code{"counterfactual"}, \code{conf_cov}, \code{ns}, 
#'   \code{pd_time}, \code{swtrt}, and \code{swtrt_time}.
#' 
#' * \code{fit_logis}: The list of fitted pooled logistic regression 
#'   models for treatment switching using g-estimation.
#' 
#' * \code{data_outcome}: The input data for the outcome Cox model
#'   of counterfactual unswitched survival times. 
#'   The variables include \code{id}, \code{stratum}, \code{"t_star"}, 
#'   \code{"d_star"}, \code{"treated"}, \code{base_cov} and \code{treat}.
#'
#' * \code{km_outcome}: The Kaplan-Meier estimates of the survival
#'   functions for the treatment and control groups based on the
#'   counterfactual unswitched survival times.
#'   
#' * \code{lr_outcome}: The log-rank test results for the treatment
#'   effect based on the counterfactual unswitched survival times.
#'   
#' * \code{fit_outcome}: The fitted outcome Cox model.
#'
#' * \code{fail}: Whether a model fails to converge.
#' 
#' * \code{psimissing}: Whether the `psi` parameter cannot be estimated.
#'
#' * \code{settings}: A list containing the input parameter values.
#'
#' * \code{psi_trt}: The estimated causal parameter for the experimental 
#'   group if \code{swtrt_control_only} is \code{FALSE}.
#'   
#' * \code{psi_trt_roots}: Vector of \code{psi_trt} values for the 
#'  experimental group at which the Z-statistic is zero, identified using
#'  grid search and linear interpolation, if \code{swtrt_control_only}
#'  is \code{FALSE}.
#'
#' * \code{psi_trt_CI}: The confidence interval for \code{psi_trt} if
#'   \code{swtrt_control_only} is \code{FALSE}.
#'
#' * \code{fail_boots}: The indicators for failed bootstrap samples
#'   if \code{boot} is \code{TRUE}.
#'
#' * \code{fail_boots_data}: The data for failed bootstrap samples
#'   if \code{boot} is \code{TRUE}.
#'
#' * \code{hr_boots}: The bootstrap hazard ratio estimates 
#'   if \code{boot} is \code{TRUE}.
#'
#' * \code{psi_boots}: The bootstrap \code{psi} estimates 
#'   if \code{boot} is \code{TRUE}.
#' 
#' * \code{psi_trt_boots}: The bootstrap \code{psi_trt} estimates 
#'   if \code{boot} is \code{TRUE} and \code{swtrt_control_only} is 
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
#' library(dplyr)
#' 
#' sim1 <- tsegestsim(
#'   n = 500, allocation1 = 2, allocation2 = 1, pbprog = 0.5, 
#'   trtlghr = -0.5, bprogsl = 0.3, shape1 = 1.8, 
#'   scale1 = 360, shape2 = 1.7, scale2 = 688, 
#'   pmix = 0.5, admin = 5000, pcatnotrtbprog = 0.5, 
#'   pcattrtbprog = 0.25, pcatnotrt = 0.2, pcattrt = 0.1, 
#'   catmult = 0.5, tdxo = 1, ppoor = 0.1, pgood = 0.04, 
#'   ppoormet = 0.4, pgoodmet = 0.2, xomult = 1.4188308, 
#'   milestone = 546, outputRawDataset = 1, seed = 2000)
#'   
#' data1 <- sim1$paneldata %>%
#'   mutate(visit7on = ifelse(progressed, tstop > timePFSobs + 105, 0))
#'
#' fit1 <- tsegest(
#'   data = data1, id = "id", 
#'   tstart = "tstart", tstop = "tstop", event = "event", 
#'   treat = "trtrand", censor_time = "censor_time", 
#'   pd = "progressed", pd_time = "timePFSobs", 
#'   swtrt = "xo", swtrt_time = "xotime", 
#'   base_cov = "bprog", 
#'   conf_cov = c("bprog*cattdc", "timePFSobs", "visit7on"), 
#'   ns_df = 3, low_psi = -1, hi_psi = 1, n_eval_z = 101,
#'   recensor = TRUE, admin_recensor_only = TRUE, 
#'   swtrt_control_only = TRUE, alpha = 0.05, 
#'   ties = "efron", tol = 1.0e-6, offset = 0, 
#'   boot = FALSE)
#'   
#' fit1
#' 
#' @export
tsegest <- function(data, id = "id", stratum = "", 
                    tstart = "tstart", tstop = "tstop", event = "event",
                    treat = "treat", censor_time = "censor_time",
                    pd = "pd", pd_time = "pd_time",
                    swtrt = "swtrt", swtrt_time = "swtrt_time",
                    base_cov = "", conf_cov = "",
                    strata_main_effect_only = TRUE, 
                    ns_df = 3, firth = FALSE, flic = FALSE,
                    low_psi = -2, hi_psi = 2, n_eval_z = 401,
                    recensor = TRUE, admin_recensor_only = TRUE,
                    swtrt_control_only = TRUE, 
                    gridsearch = TRUE, root_finding = "brent",
                    alpha = 0.05, ties = "efron", tol = 1.0e-6, offset = 1, 
                    boot = TRUE, n_boot = 1000, seed = NA) {

  rownames(data) = NULL

  elements = c(id, stratum, tstart, tstop, event, treat, censor_time, 
               pd, swtrt)
  elements = unique(elements[elements != "" & elements != "none"])
  fml = formula(paste("~", paste(elements, collapse = "+")))
  mf = model.frame(fml, data = data, na.action = na.omit)
  
  rownum = as.integer(rownames(mf))
  df = data[rownum,]

  nvar = length(base_cov)
  if (missing(base_cov) || is.null(base_cov) || (nvar == 1 && (
    base_cov[1] == "" || tolower(base_cov[1]) == "none"))) {
    p = 0
  } else {
    fml1 = formula(paste("~", paste(base_cov, collapse = "+")))
    vnames = rownames(attr(terms(fml1), "factors"))
    p = length(vnames)
  }

  if (p >= 1) {
    mf1 <- model.frame(fml1, data = df, na.action = na.pass)
    mm <- model.matrix(fml1, mf1)
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
    fml2 = formula(paste("~", paste(conf_cov, collapse = "+")))
    vnames2 = rownames(attr(terms(fml2), "factors"))
    p2 = length(vnames2)
  }

  if (p2 >= 1) {
    mf2 <- model.frame(fml2, data = df, na.action = na.pass)
    mm2 <- model.matrix(fml2, mf2)
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

  out <- tsegestcpp(
    data = df, id = id, stratum = stratum, 
    tstart = tstart, tstop = tstop, event = event,
    treat = treat, censor_time = censor_time, 
    pd = pd, pd_time = pd_time, 
    swtrt = swtrt, swtrt_time = swtrt_time,
    base_cov = varnames, conf_cov = varnames2,
    strata_main_effect_only = strata_main_effect_only,
    ns_df = ns_df, firth = firth, flic = flic,
    low_psi = low_psi, hi_psi = hi_psi, n_eval_z = n_eval_z, 
    recensor = recensor, admin_recensor_only = admin_recensor_only,
    swtrt_control_only = swtrt_control_only, 
    gridsearch = gridsearch, root_finding = root_finding, 
    alpha = alpha, ties = ties, tol = tol, offset = offset, 
    boot = boot, n_boot = n_boot, seed = seed)
  
  if (!out$psimissing) {
    K = ifelse(swtrt_control_only, 1, 2)
    for (h in 1:K) {
      out$data_logis[[h]]$data$uid <- NULL
      out$data_nullcox[[h]]$data$ustratum <- NULL
    }
    
    out$data_outcome$uid <- NULL
    out$data_outcome$ustratum <- NULL
    
    df[, "tstart"] = df[, tstart]
    df[, "tstop"] = df[, tstop]
    
    if (p >= 1) {
      add_vars <- setdiff(vnames, varnames)
      if (length(add_vars) > 0) {
        out$data_outcome <- merge(out$data_outcome, df[, c(id, add_vars)], 
                                  by = id, all.x = TRUE, sort = FALSE)
      }
      
      del_vars <- setdiff(varnames, vnames)
      if (length(del_vars) > 0) {
        out$data_outcome[, del_vars] <- NULL
      }
    }
    
    if (p2 >= 1) {
      tem_vars <- c(pd_time, swtrt, swtrt_time)
      add_vars <- c(setdiff(vnames2, varnames2), tem_vars)
      avars <- setdiff(add_vars, names(out$data_logis[[1]]$data))
      if (length(avars) > 0) {
        for (h in 1:K) {
          out$data_logis[[h]]$data <- 
            merge(out$data_logis[[h]]$data, 
                  df[, c(id, "tstart", "tstop", avars)], 
                  by = c(id, "tstart", "tstop"), all.x = TRUE, sort = FALSE)
        }
      }
      
      del_vars <- setdiff(varnames2, vnames2)
      if (length(del_vars) > 0) {
        for (h in 1:K) {
          out$data_logis[[h]]$data[, del_vars] <- NULL
        }
      }
    }
  }
  
  
  # convert treatment back to a factor variable if needed
  if (is.factor(data[[treat]])) {
    levs = levels(data[[treat]])
    
    out$event_summary[[treat]] <- factor(out$event_summary[[treat]], 
                                         levels = c(1,2), labels = levs)
    
    for (h in 1:2) {
      out$data_switch[[h]][[treat]] <- factor(
        out$data_switch[[h]][[treat]], levels = c(1,2), labels = levs)
      
      out$km_switch[[h]][[treat]] <- factor(
        out$km_switch[[h]][[treat]], levels = c(1,2), labels = levs)
      
      out$data_logis[[h]]$data[[treat]] <- factor(
        out$data_logis[[h]]$data[[treat]], levels = c(1,2), labels = levs)
      
      out$data_nullcox[[h]]$data[[treat]] <- factor(
        out$data_nullcox[[h]]$data[[treat]], levels = c(1,2), labels = levs)
    }
    
    out$data_outcome[[treat]] <- factor(out$data_outcome[[treat]], 
                                        levels = c(1,2), labels = levs)
    
    out$km_outcome[[treat]] <- factor(out$km_outcome[[treat]], 
                                      levels = c(1,2), labels = levs)
  }
  
  
  out$settings <- list(
    data = data, id = id, stratum = stratum, tstart = tstart, 
    tstop = tstop, event = event, treat = treat,
    censor_time = censor_time, pd = pd, pd_time = pd_time,
    swtrt = swtrt, swtrt_time = swtrt_time,
    base_cov = base_cov, conf_cov = conf_cov,
    strata_main_effect_only = strata_main_effect_only,
    ns_df = ns_df, firth = firth, flic = flic,
    low_psi = low_psi, hi_psi = hi_psi, n_eval = n_eval_z,
    recensor = recensor, admin_recensor_only = admin_recensor_only,
    swtrt_control_only = swtrt_control_only, 
    gridsearch = gridsearch, root_finding = root_finding,
    alpha = alpha, ties = ties, tol = tol, offset = offset,
    boot = boot, n_boot = n_boot, seed = seed
  )
  
  class(out) <- "tsegest"
  out
}
