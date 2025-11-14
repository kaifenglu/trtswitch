#' @title Simple Two-Stage Estimation (TSEsimp) for Treatment Switching
#' @description Estimates the causal parameter by fitting an accelerated 
#' failure time (AFT) model comparing post-progression survival between 
#' switchers and non-switchers, and derives the adjusted hazard ratio 
#' from the Cox model using counterfactual \emph{unswitched} survival 
#' times based on the estimated causal parameter.
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
#'   * \code{pd_time}: The time from randomization to disease progression.
#'
#'   * \code{swtrt}: The treatment switch indicator, 1=switch, 0=no switch.
#'
#'   * \code{swtrt_time}: The time from randomization to treatment switch.
#'
#'   * \code{base_cov}: The baseline covariates (excluding treat).
#'
#'   * \code{base2_cov}: The baseline and secondary baseline
#'     covariates (excluding treat).
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
#' @param base2_cov The names of baseline and secondary baseline covariates
#'   (excluding treat) in the input data for the AFT model for 
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
#' @param ties The method for handling ties in the Cox model, either
#'   "breslow" or "efron" (default).
#' @param offset The offset to calculate the time disease progression to 
#'   death or censoring. We can set \code{offset} equal to 0 (no offset), 
#'   and 1 (default), 1/30.4375, or 1/365.25 if the time unit is day, month, 
#'   or year, respectively.
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
#' * Estimate the causal parameter \eqn{\psi} by fitting an AFT model 
#'   comparing post-progression survival between switchers and non-switchers 
#'   in the control group who experienced disease progression.
#'
#' * Compute counterfactual survival times for control patients using 
#'   the estimated \eqn{\psi}.
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
#' * \code{event_summary}: A data frame containing the count and percentage
#'   of deaths, disease progressions, and switches by treatment arm.
#'   
#' * \code{data_aft}: A list of input data for the AFT model by treatment 
#'   group. The variables include \code{id}, \code{stratum}, \code{"pps"}, 
#'   \code{"event"}, \code{"swtrt"}, \code{base2_cov}, \code{pd_time}, 
#'   \code{swtrt_time}, and \code{time}.
#' 
#' * \code{fit_aft}: A list of fitted AFT models by treatment group.
#' 
#' * \code{res_aft}: A list of deviance residuals from the AFT models 
#'   by treatment group.
#'   
#' * \code{data_outcome}: The input data for the outcome Cox model 
#'   of counterfactual unswitched survival times.
#'   The variables include \code{id}, \code{stratum}, \code{"t_star"}, 
#'   \code{"d_star"}, \code{"treated"}, \code{base_cov}, and \code{treat}.
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
#' * \code{psi_trt_CI}: The confidence interval for \code{psi_trt} 
#'   if \code{swtrt_control_only} is \code{FALSE}.
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
#' # modify pd and dpd based on co and dco
#' shilong <- shilong %>%
#'   mutate(dpd = ifelse(co & !pd, dco, dpd),
#'          pd = ifelse(co & !pd, 1, pd)) %>%
#'   mutate(dpd = ifelse(pd & co & dco < dpd, dco, dpd))
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
#' fit1
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

  nvar2 = length(base2_cov)
  if (missing(base2_cov) || is.null(base2_cov) || (nvar2 == 1 && (
    base2_cov[1] == "" || tolower(base2_cov[1]) == "none"))) {
    p2 = 0
  } else {
    fml2 = formula(paste("~", paste(base2_cov, collapse = "+")))
    vnames2 = rownames(attr(terms(fml2), "factors"))
    p2 = length(vnames2)
  }

  if (p2 >= 1) {
    mf2 <- model.frame(fml2, data = df, na.action = na.pass)
    mm2 = model.matrix(fml2, mf2)
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
  
  
  if (!out$psimissing) {
    out$data_outcome$uid <- NULL
    out$data_outcome$ustratum <- NULL
    
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
      K = ifelse(swtrt_control_only, 1, 2)
      tem_vars <- c(pd_time, swtrt_time, time)
      add_vars <- c(setdiff(vnames2, varnames2), tem_vars)
      avars <- setdiff(add_vars, names(out$data_aft[[1]]$data))
      if (length(avars) > 0) {
        for (h in 1:K) {
          out$data_aft[[h]]$data <- merge(out$data_aft[[h]]$data, 
                                          df[, c(id, avars)], 
                                          by = id, all.x = TRUE, sort = FALSE)
        }
      }
      
      del_vars <- setdiff(varnames2, vnames2)
      if (length(del_vars) > 0) {
        for (h in 1:K) {
          out$data_aft[[h]]$data[, del_vars] <- NULL
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
      out$data_aft[[h]][[treat]] <- factor(
        out$data_aft[[h]][[treat]], levels = c(1,2), labels = levs)
    }

    out$data_outcome[[treat]] <- factor(out$data_outcome[[treat]], 
                                        levels = c(1,2), labels = levs)
    
    out$km_outcome[[treat]] <- factor(out$km_outcome[[treat]], 
                                      levels = c(1,2), labels = levs)
  }
  
  
  out$settings <- list(
    data = data, id = id, stratum = stratum, time = time, 
    event = event, treat = treat, censor_time = censor_time, 
    pd = pd, pd_time = pd_time, swtrt = swtrt, 
    swtrt_time = swtrt_time, base_cov = base_cov, 
    base2_cov = base2_cov, aft_dist = aft_dist,
    strata_main_effect_only = strata_main_effect_only,
    recensor = recensor, admin_recensor_only = admin_recensor_only,
    swtrt_control_only = swtrt_control_only,
    alpha = alpha, ties = ties, offset = offset,
    boot = boot, n_boot = n_boot, seed = seed
  )
  
  class(out) <- "tsesimp"
  out
}