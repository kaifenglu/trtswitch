#' @title The Simple Two-Stage Estimation (TSE) Method for Treatment 
#' Switching
#' @description Obtains the causal parameter estimate of the AFT model 
#' for switching after disease progression and the hazard ratio estimate 
#' of the outcome Cox model to adjust for treatment switching.
#'
#' @param data The input data frame that contains the following variables:
#'
#'   * \code{id}: The subject id.
#'
#'   * \code{stratum}: The stratum.
#'
#'   * \code{time}: The survival time for right censored data.
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
#'   * \code{base_cov}: The baseline covariates (excluding treat).
#'
#'   * \code{base2_cov}: The baseline and secondary baseline
#'     covariates (excluding swtrt).
#'
#' @param id The name of the id variable in the input data.
#' @param stratum The name(s) of the stratum variable(s) in the input data.
#' @param time The name of the time variable in the input data.
#' @param event The name of the event variable in the input data.
#' @param treat The name of the treatment variable in the input data.
#' @param censor_time The name of the censor_time variable in the input data.
#' @param pd The name of the pd variable in the input data.
#' @param pd_time The name of the pd_time variable in the input data.
#' @param swtrt The name of the swtrt variable in the input data.
#' @param swtrt_time The name of the swtrt_time variable in the input data.
#' @param base_cov The names of baseline covariates (excluding
#'   treat) in the input data for the outcome Cox model.
#' @param base2_cov The names of secondary baseline covariates
#'   (excluding swtrt) in the input data for the AFT model for 
#'   post-progression survival.
#' @param aft_dist The assumed distribution for time to event for the AFT
#'   model. Options include "exponential", "weibull" (default), 
#'   "loglogistic", and "lognormal".
#' @param strata_main_effect_only Whether to only include the strata main
#'   effects in the AFT model. Defaults to \code{TRUE}, otherwise all
#'   possible strata combinations will be considered in the AFT model.
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
#' * Fit an AFT model to post-progression survival data to estimate 
#'   the causal parameter \eqn{\psi} based on the patients 
#'   in the control group who had disease progression.
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
#'   i.e., "AFT model" or "bootstrap".
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
#' * \code{data_aft}: A list of input data for the AFT model by treatment 
#'   group.
#' 
#' * \code{fit_aft}: A list of fitted AFT models by treatment group.
#' 
#' * \code{data_outcome}: The input data for the outcome Cox model. 
#' 
#' * \code{fit_outcome}: The fitted outcome Cox model.
#'
#' * \code{settings}: A list with the following components:
#'
#'     - \code{aft_dist}: The distribution for time to event for the AFT
#'       model.
#'
#'     - \code{strata_main_effect_only}: Whether to only include the strata
#'       main effects in the AFT model.
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
#' Nicholas R Latimer, KR Abrams, PC Lambert, MK Crowther, AJ Wailoo,
#' JP Morden, RL Akehurst, and MJ Campbell.
#' Adjusting for treatment switching in randomised controlled trials - A
#' simulation study and a simplified two-stage method.
#' Statistical Methods in Medical Research. 2017;26(2):724-751.
#'
#' @examples
#'
#' library(dplyr)
#'
#' # the eventual survival time
#' shilong1 <- shilong %>%
#'   arrange(bras.f, id, tstop) %>%
#'   group_by(bras.f, id) %>%
#'   slice(n()) %>%
#'   select(-c("ps", "ttc", "tran"))
#'
#' # the last value of time-dependent covariates before pd
#' shilong2 <- shilong %>%
#'   filter(pd == 0 | tstart <= dpd) %>%
#'   arrange(bras.f, id, tstop) %>%
#'   group_by(bras.f, id) %>%
#'   slice(n()) %>%
#'   select(bras.f, id, ps, ttc, tran)
#'
#' # combine baseline and time-dependent covariates
#' shilong3 <- shilong1 %>%
#'   left_join(shilong2, by = c("bras.f", "id"))
#'
#' # apply the two-stage method
#' fit1 <- tsesimp(
#'   data = shilong3, id = "id", time = "tstop", event = "event",
#'   treat = "bras.f", censor_time = "dcut", pd = "pd",
#'   pd_time = "dpd", swtrt = "co", swtrt_time = "dco",
#'   base_cov = c("agerand", "sex.f", "tt_Lnum", "rmh_alea.c",
#'                 "pathway.f"),
#'   base2_cov = c("agerand", "sex.f", "tt_Lnum", "rmh_alea.c",
#'                 "pathway.f", "ps", "ttc", "tran"),
#'   aft_dist = "weibull", alpha = 0.05,
#'   recensor = TRUE, swtrt_control_only = FALSE, offset = 1,
#'   boot = FALSE)
#'
#' c(fit1$hr, fit1$hr_CI)
#'
#' @export
tsesimp <- function(data, id = "id", stratum = "", time = "time", 
                    event = "event", treat = "treat", 
                    censor_time = "censor_time",
                    pd = "pd", pd_time = "pd_time",
                    swtrt = "swtrt", swtrt_time = "swtrt_time",
                    base_cov = "", base2_cov = "",
                    aft_dist = "weibull", strata_main_effect_only = TRUE,
                    recensor = TRUE, admin_recensor_only = TRUE,
                    swtrt_control_only = TRUE, alpha = 0.05, 
                    ties = "efron", offset = 1, 
                    boot = TRUE, n_boot = 1000, seed = NA) {

  rownames(data) = NULL

  elements = c(stratum, time, event, treat, censor_time, pd, swtrt)
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

  nvar2 = length(base2_cov)
  if (missing(base2_cov) || is.null(base2_cov) || (nvar2 == 1 && (
    base2_cov[1] == "" || tolower(base2_cov[1]) == "none"))) {
    p2 = 0
  } else {
    t1 = terms(formula(paste("~", paste(base2_cov, collapse = "+"))))
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

  out <- tsesimpcpp(
    data = df, id = id, stratum = stratum, time = time, 
    event = event, treat = treat, censor_time = censor_time, 
    pd = pd, pd_time = pd_time, swtrt = swtrt, 
    swtrt_time = swtrt_time, base_cov = varnames, 
    base2_cov = varnames2, aft_dist = aft_dist,
    strata_main_effect_only = strata_main_effect_only,
    recensor = recensor, admin_recensor_only = admin_recensor_only,
    swtrt_control_only = swtrt_control_only, alpha = alpha,
    ties = ties, offset = offset, 
    boot = boot, n_boot = n_boot, seed = seed)
  
  out$data_outcome$uid <- NULL
  out$data_outcome$ustratum <- NULL
  
  
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
    t1 = terms(formula(paste("~", paste(base2_cov, collapse = "+"))))
    t2 = attr(t1, "factors")
    t3 = rownames(t2)
    
    K = ifelse(swtrt_control_only, 1, 2)
    add_vars <- setdiff(t3, varnames2)
    if (length(add_vars) > 0) {
      for (h in 1:K) {
        out$data_aft[[h]]$data <- merge(out$data_aft[[h]]$data, 
                                        df[, c(id, add_vars)], 
                                        by = id, all.x = TRUE, sort = FALSE)
      }
    }
    
    del_vars <- setdiff(varnames2, t3)
    if (length(del_vars) > 0) {
      for (h in 1:K) {
        out$data_aft[[h]]$data[, del_vars] <- NULL
      }
    }
    
    for (h in 1:K) {
      out$data_aft[[h]]$data <- out$data_aft[[h]]$data[
        , !startsWith(names(out$data_aft[[h]]$data), "stratum_")]
    }
  }

  out
}