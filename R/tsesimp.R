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
#' @param seed The seed to reproduce the bootstrap results.
#' @param nthreads The number of threads to use in bootstrapping (0 means 
#'   the default RcppParallel behavior)
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
#' * \code{pvalue}: The two-sided p-value.
#'
#' * \code{pvalue_type}: The type of two-sided p-value for treatment effect, 
#'   i.e., "Cox model" or "bootstrap".
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
                    swtrt_control_only = TRUE, 
                    alpha = 0.05, ties = "efron", offset = 1, 
                    boot = TRUE, n_boot = 1000, seed = 0, 
                    nthreads = 0) {

  # validate input
  if (!inherits(data, "data.frame")) {
    stop("Input 'data' must be a data frame");
  }
  
  if (inherits(data, "data.table") || inherits(data, "tbl") || 
      inherits(data, "tbl_df")) {
    df <- as.data.frame(data)
  } else {
    df <- data
  }
  
  for (nm in c(id, time, event, treat, censor_time, pd, pd_time, 
               swtrt, swtrt_time)) {
    if (!is.character(nm) || length(nm) != 1) {
      stop(paste(nm, "must be a single character string."));
    }
  }
  
  # Respect user-requested number of threads (best effort)
  if (nthreads > 0) {
    n_physical_cores <- parallel::detectCores(logical = FALSE)
    RcppParallel::setThreadOptions(min(nthreads, n_physical_cores))
  }
  
  # select complete cases for the relevant variables
  elements = unique(c(id, stratum, time, event, treat, censor_time, pd, swtrt))
  elements = elements[elements != ""]
  fml_all <- formula(paste("~", paste(elements, collapse = "+")))
  var_all <- all.vars(fml_all)
  rows_ok <- which(complete.cases(df[, var_all, drop = FALSE]))
  if (length(rows_ok) == 0) stop("No complete cases found for the specified variables.")
  df <- df[rows_ok, , drop = FALSE]
  
  misscovariates <- length(base_cov) == 0 || 
    (length(base_cov) == 1 && (base_cov[1] == ""))
  
  if (!misscovariates) {
    fml_cov <- as.formula(paste("~", paste(base_cov, collapse = "+")))
    vnames <- rownames(attr(terms(fml_cov), "factors"))
    
    # QUICK PATH: if all covariates present in df and are numeric, avoid model.matrix
    cov_present <- base_cov %in% names(df)
    all_numeric <- FALSE
    if (all(cov_present)) {
      all_numeric <- all(vapply(df[ base_cov ], is.numeric, logical(1)))
    }
    
    if (all_numeric) {
      # Build design columns directly from numeric covariates (intercept + columns)
      # This avoids model.matrix and is valid when covariates are simple numeric columns.
      varnames <- base_cov
    } else {
      # FALLBACK (existing robust behavior): use model.frame + model.matrix on df
      mf <- model.frame(fml_cov, data = df, na.action = na.pass)
      mm <- model.matrix(fml_cov, mf)
      colnames(mm) <- make.names(colnames(mm))
      varnames <- colnames(mm)[-1]
      missing_cols <- setdiff(varnames, names(df))
      if (length(missing_cols) > 0) {
        for (vn in missing_cols) df[[vn]] <- mm[, vn, drop = TRUE]
      }
    }
  } else {
    varnames <- ""
  }
  
  misscovariates2 <- length(base2_cov) == 0 || 
    (length(base2_cov) == 1 && (base2_cov[1] == ""))
  
  if (!misscovariates2) {
    fml_cov2 <- as.formula(paste("~", paste(base2_cov, collapse = "+")))
    vnames2 <- rownames(attr(terms(fml_cov2), "factors"))
    
    # QUICK PATH: if all covariates present in df and are numeric, avoid model.matrix
    cov_present2 <- base2_cov %in% names(df)
    all_numeric2 <- FALSE
    if (all(cov_present2)) {
      all_numeric2 <- all(vapply(df[ base2_cov ], is.numeric, logical(1)))
    }
    
    if (all_numeric2) {
      # Build design columns directly from numeric covariates (intercept + columns)
      # This avoids model.matrix and is valid when covariates are simple numeric columns.
      varnames2 <- base2_cov
    } else {
      # FALLBACK (existing robust behavior): use model.frame + model.matrix on df
      mf2 <- model.frame(fml_cov2, data = df, na.action = na.pass)
      mm2 <- model.matrix(fml_cov2, mf2)
      colnames(mm2) <- make.names(colnames(mm2))
      varnames2 <- colnames(mm2)[-1]
      missing_cols2 <- setdiff(varnames2, names(df))
      if (length(missing_cols2) > 0) {
        for (vn in missing_cols2) df[[vn]] <- mm2[, vn, drop = TRUE]
      }
    }
  } else {
    varnames2 <- ""
  }

  out <- tsesimpcpp(
    df = df, id = id, stratum = stratum, time = time, 
    event = event, treat = treat, 
    censor_time = censor_time, 
    pd = pd, pd_time = pd_time, 
    swtrt = swtrt, swtrt_time = swtrt_time, 
    base_cov = varnames, base2_cov = varnames2, 
    aft_dist = aft_dist, strata_main_effect_only = strata_main_effect_only,
    recensor = recensor, admin_recensor_only = admin_recensor_only,
    swtrt_control_only = swtrt_control_only, 
    alpha = alpha, ties = ties, offset = offset, 
    boot = boot, n_boot = n_boot, seed = seed)
  
  if (!out$psimissing) {
    out$data_outcome$uid <- NULL
    out$data_outcome$ustratum <- NULL
    
    if (!misscovariates) {
      add_vars <- setdiff(vnames, varnames)
      if (length(add_vars) > 0) {
        frame_df <- out$data_outcome
        idx <- match(frame_df[[id]], df[[id]])
        for (var in add_vars) frame_df[[var]] <- df[[var]][idx]
        out$data_outcome <- frame_df
      }
      
      del_vars <- setdiff(varnames, vnames)
      if (length(del_vars) > 0) {
        out$data_outcome[, del_vars] <- NULL
      }
    }
    
    if (!misscovariates2) {
      K = ifelse(swtrt_control_only, 1, 2)
      tem_vars <- c(pd_time, swtrt_time, time)
      add_vars <- c(setdiff(vnames2, varnames2), tem_vars)
      avars <- setdiff(add_vars, names(out$data_aft[[1]]$data))
      if (length(avars) > 0) {
        for (h in 1:K) {
          frame_df <- out$data_aft[[h]]$data
          idx <- match(frame_df[[id]], df[[id]])
          for (var in add_vars) frame_df[[var]] <- df[[var]][idx]
          out$data_aft[[h]]$data <- frame_df
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
    levs <- levels(data[[treat]])
    mf <- function(x) if (is.null(x)) x else factor(x, levels = c(1,2), labels = levs)
    
    # apply mf to a set of named containers that are data.frames with a column named `treat`
    for (nm in c("event_summary", "data_outcome", "km_outcome")) {
      out[[nm]][[treat]] <- mf(out[[nm]][[treat]])
    }
    
    # and for the list-of-lists
    out$data_aft <- lapply(out$data_aft, function(x) { x[[treat]] <- mf(x[[treat]]); x })
  }
  
  out$settings <- list(
    data = data, id = id, stratum = stratum, time = time, 
    event = event, treat = treat, 
    censor_time = censor_time, 
    pd = pd, pd_time = pd_time, 
    swtrt = swtrt, swtrt_time = swtrt_time, 
    base_cov = base_cov, base2_cov = base2_cov, 
    aft_dist = aft_dist, strata_main_effect_only = strata_main_effect_only,
    recensor = recensor, admin_recensor_only = admin_recensor_only,
    swtrt_control_only = swtrt_control_only,
    alpha = alpha, ties = ties, offset = offset,
    boot = boot, n_boot = n_boot, seed = seed
  )
  
  class(out) <- "tsesimp"
  out
}