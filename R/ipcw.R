#' @title Inverse Probability of Censoring Weights (IPCW) for Treatment 
#' Switching
#' @description Excludes data after treatment switching and fits a switching 
#' model to estimate the probability of not switching. The inverse of 
#' these probabilities (inverse probability of censoring weights) are then 
#' used as weights in a weighted Cox model to obtain the adjusted hazard 
#' ratio.
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
#'   * \code{swtrt}: The treatment switch indicator, 1=switch, 0=no switch.
#'
#'   * \code{swtrt_time}: The time from randomization to treatment switch.
#'
#'   * \code{base_cov}: The baseline covariates (excluding treat) used in
#'     the outcome model.
#'
#'   * \code{numerator}: The baseline covariates (excluding treat) used in 
#'     the numerator switching model for stabilized weights.
#'
#'   * \code{denominator}: The baseline (excluding treat) and time-dependent 
#'     covariates used in the denominator switching model.
#'
#' @param id The name of the id variable in the input data.
#' @param stratum The name(s) of the stratum variable(s) in the input data.
#' @param tstart The name of the tstart variable in the input data.
#' @param tstop The name of the tstop variable in the input data.
#' @param event The name of the event variable in the input data.
#' @param treat The name of the treatment variable in the input data.
#' @param swtrt The name of the swtrt variable in the input data.
#' @param swtrt_time The name of the swtrt_time variable in the input data.
#' @param base_cov The names of baseline covariates (excluding
#'   treat) in the input data for the Cox model.
#' @param numerator The names of baseline covariates 
#'   (excluding treat) in the input data for the numerator switching 
#'   model for stabilized weights.
#' @param denominator The names of baseline  (excluding treat) and 
#'   time-dependent covariates in the input data for the denominator 
#'   switching model.
#' @param logistic_switching_model Whether a pooled logistic regression 
#'   switching model is used.
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
#' @param stabilized_weights Whether to use the stabilized weights. 
#'   The default is \code{TRUE}.
#' @param trunc The truncation fraction of the weight distribution. 
#'   Defaults to 0 for no truncation in weights.
#' @param trunc_upper_only Whether to truncate the weights from the upper
#'   end of the weight distribution only. Defaults to \code{TRUE}, otherwise
#'   the weights will be truncated from both the lower and upper ends of
#'   the distribution.
#' @param swtrt_control_only Whether treatment switching occurred only in
#'   the control group. The default is \code{TRUE}.
#' @param alpha The significance level to calculate confidence intervals. 
#' @param ties The method for handling ties in the Cox model, either
#'   "breslow" or "efron" (default).
#' @param boot Whether to use bootstrap to obtain the confidence
#'   interval for hazard ratio. Defaults to \code{FALSE}.
#' @param n_boot The number of bootstrap samples.
#' @param seed The seed to reproduce the bootstrap results. The default is 
#'   `NA`, in which case, the seed from the environment will be used.
#'
#' @details The hazard ratio and confidence interval under a no-switching 
#' scenario are obtained as follows:
#'
#' * Exclude all observations after treatment switch.
#'
#' * Define the crossover and event indicators for the last time interval 
#'   of each subject.
#'
#' * For time-dependent Cox switching models, replicate unique event times 
#'   across treatment arms within each subject.
#'
#' * Fit the denominator switching model (and numerator model for 
#'   stabilized weights) to estimate inverse probability of censoring 
#'   weights. Either a Cox model with time-dependent covariates or 
#'   a pooled logistic regression model can be used.
#'   
#'     - For the pooled logistic regression model, the probability of 
#'       remaining uncensored (i.e., not switching) is calculated as 
#'       \eqn{1 - \hat{p}_{\text{switch}}} 
#'       and accumulated over time up to the start of each interval.
#'       
#'     - For the time-dependent Cox model, the probability of remaining 
#'       unswitched is derived from the estimated baseline hazard and 
#'       predicted risk score up to the end of each interval.
#'
#' * Fit a weighted Cox model to the outcome survival times (excluding 
#'   data after switching) to estimate the hazard ratio.
#'
#' * Construct the p-value and confidence interval for the hazard ratio 
#'   using either robust sandwich variance or bootstrapping. When 
#'   bootstrapping is used, the confidence interval and p-value are 
#'   based on a t-distribution with \code{n_boot - 1} degrees of freedom.
#'
#' @return A list with the following components:
#'
#' * \code{logrank_pvalue}:  The two-sided p-value of the log-rank test 
#'   for the ITT analysis.
#'
#' * \code{cox_pvalue}: The two-sided p-value for treatment effect based on
#'   the weighted Cox model excluding data after treatment switch. 
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
#' * \code{data_switch}: A list of input data for the switching models by 
#'   treatment group. The variables include \code{id}, \code{stratum}, 
#'   \code{"tstart"}, \code{"tstop"}, \code{"cross"}, \code{denominator}, 
#'   \code{swtrt}, and \code{swtrt_time}. For logistic switching models, 
#'   \code{stratum} variables are converted to dummy variables, and 
#'   natural cubic spline basis variables are created for the visit-specific 
#'   intercepts.
#'
#' * \code{fit_switch}: A list of fitted switching models for the
#'   denominator and numerator by treatment group.
#'
#' * \code{data_outcome}: The input data for the outcome Cox model 
#'   including the inverse probability of censoring weights.
#'   The variables include \code{id}, \code{stratum}, \code{"tstart"}, 
#'   \code{"tstop"}, \code{"event"}, \code{"treated"}, 
#'   \code{"unstablized_weight"}, \code{"stabilized_weight"}, 
#'   \code{base_cov}, and \code{treat}.
#'   
#' * \code{km_outcome}: The Kaplan-Meier estimates of the survival
#'   functions for the treatment and control groups based on the
#'   weighted outcome data.
#'   
#' * \code{lr_outcome}: The log-rank test results for the treatment
#'   effect based on the weighted outcome data.
#'
#' * \code{fit_outcome}: The fitted outcome Cox model.
#'
#' * \code{fail}: Whether a model fails to converge.
#'
#' * \code{settings}: A list with the following components:
#'
#'     - \code{logistic_switching_model}: Whether a pooled logistic 
#'       regression switching model is used.
#'       
#'     - \code{strata_main_effect_only}: Whether to only include the strata 
#'       main effects in the logistic regression switching model. 
#'       
#'     - \code{ns_df}: Degrees of freedom for the natural cubic spline.
#'   
#'     - \code{firth}: Whether the Firth's penalized likelihood is used.
#'       
#'     - \code{flic}: Whether to apply intercept correction.
#'       
#'     - \code{stabilized_weights}: Whether to use the stabilized weights.
#'
#'     - \code{trunc}: The truncation fraction of the weight distribution.
#'
#'     - \code{trunc_upper_only}: Whether to truncate the weights from the
#'       upper end of the distribution only.
#'
#'     - \code{swtrt_control_only} Whether treatment switching occurred 
#'       only in the control group.
#'
#'     - \code{alpha}: The significance level to calculate confidence
#'       intervals.
#'
#'     - \code{ties}: The method for handling ties in the Cox model.
#'
#'     - \code{boot}: Whether to use bootstrap to obtain the confidence
#'       interval for hazard ratio.
#'
#'     - \code{n_boot}: The number of bootstrap samples.
#'
#'     - \code{seed}: The seed to reproduce the bootstrap results.
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
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @references
#' James M. Robins and Dianne M. Finkelstein.
#' Correcting for noncompliance and dependent censoring in an AIDS clinical
#' trial with inverse probability of censoring weighted (IPCW) log-rank tests.
#' Biometrics. 2000;56(3):779-788.
#'
#' @examples
#'
#' # Example 1: pooled logistic regression switching model
#' library(dplyr)
#' 
#' sim1 <- tssim(
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
#' fit1 <- ipcw(
#'   sim1[[1]], id = "id", tstart = "tstart", 
#'   tstop = "tstop", event = "event", treat = "trtrand", 
#'   swtrt = "xo", swtrt_time = "xotime", 
#'   base_cov = "bprog", numerator = "bprog", 
#'   denominator = c("bprog", "L"),
#'   logistic_switching_model = TRUE, ns_df = 3,
#'   swtrt_control_only = TRUE, boot = FALSE)
#'   
#' c(fit1$hr, fit1$hr_CI) 
#' 
#' # Example 2: time-dependent covariates Cox switching model
#' 
#' fit2 <- ipcw(
#'   shilong, id = "id", tstart = "tstart", tstop = "tstop", 
#'   event = "event", treat = "bras.f", swtrt = "co", 
#'   swtrt_time = "dco", 
#'   base_cov = c("agerand", "sex.f", "tt_Lnum", "rmh_alea.c", 
#'                "pathway.f"),
#'   numerator = c("agerand", "sex.f", "tt_Lnum", "rmh_alea.c", 
#'                 "pathway.f"),
#'   denominator = c("agerand", "sex.f", "tt_Lnum", "rmh_alea.c",
#'                   "pathway.f", "ps", "ttc", "tran"),
#'   swtrt_control_only = FALSE, boot = FALSE)
#'
#' c(fit2$hr, fit2$hr_CI)
#'
#' @export
ipcw <- function(data, id = "id", stratum = "", tstart = "tstart",
                 tstop = "tstop", event = "event", treat = "treat",
                 swtrt = "swtrt", swtrt_time = "swtrt_time",
                 base_cov = "", numerator = "", denominator = "",
                 logistic_switching_model = FALSE, 
                 strata_main_effect_only = TRUE, 
                 ns_df = 3, firth = FALSE, flic = FALSE, 
                 stabilized_weights = TRUE, 
                 trunc = 0, trunc_upper_only = TRUE,
                 swtrt_control_only = TRUE, 
                 alpha = 0.05, ties = "efron", 
                 boot = FALSE, n_boot = 1000, seed = NA) {
  
  rownames(data) = NULL
  
  elements = c(id, stratum, tstart, tstop, event, treat, swtrt)
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
  
  nvar2 = length(numerator)
  if (missing(numerator) || is.null(numerator) || (nvar2 == 1 && (
    numerator[1] == "" || tolower(numerator[1]) == "none"))) {
    p2 = 0
  } else {
    fml2 = formula(paste("~", paste(numerator, collapse = "+")))
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
  
  nvar3 = length(denominator)
  if (missing(denominator) || is.null(denominator) || (nvar3 == 1 && (
    denominator[1] == "" || tolower(denominator[1]) == "none"))) {
    p3 = 0
  } else {
    fml3 = formula(paste("~", paste(denominator, collapse = "+")))
    vnames3 = rownames(attr(terms(fml3), "factors"))
    p3 = length(vnames3)
  }
  
  if (p3 >= 1) {
    mf3 <- model.frame(fml3, data = df, na.action = na.pass)
    mm3 <- model.matrix(fml3, mf3)
    colnames(mm3) = make.names(colnames(mm3))
    varnames3 = colnames(mm3)[-1]
    for (i in 1:length(varnames3)) {
      if (!(varnames3[i] %in% names(df))) {
        df[,varnames3[i]] = mm3[,varnames3[i]]
      }
    }
  } else {
    varnames3 = ""
  }
  
  out <- ipcwcpp(
    data = df, id = id, stratum = stratum, tstart = tstart,
    tstop = tstop, event = event, treat = treat, 
    swtrt = swtrt, swtrt_time = swtrt_time, 
    base_cov = varnames, numerator = varnames2, denominator = varnames3,
    logistic_switching_model = logistic_switching_model,
    strata_main_effect_only = strata_main_effect_only,
    ns_df = ns_df, firth = firth, flic = flic, 
    stabilized_weights = stabilized_weights, 
    trunc = trunc, trunc_upper_only = trunc_upper_only,
    swtrt_control_only = swtrt_control_only, 
    alpha = alpha, ties = ties, 
    boot = boot, n_boot = n_boot, seed = seed)
  
  # --- Update df ---
  df[, "tstart"] = df[, tstart]
  df[, "tstop"] = df[, tstop]
  
  df <- df[order(df[[id]]), ]          # Sort by id
  dfu <- df[!duplicated(df[[id]]), ]   # Keep the first row for each id
  
  # --- Update data_outcome ---
  out$data_outcome$uid <- NULL
  out$data_outcome$ustratum <- NULL
  
  if (p >= 1) {
    add_vars <- setdiff(vnames, varnames)
    if (length(add_vars) > 0) {
      out$data_outcome <- 
        merge(out$data_outcome, 
              dfu[, c(id, add_vars)], 
              by = id, all.x = TRUE, sort = FALSE)
    }
    
    del_vars <- setdiff(varnames, vnames)
    if (length(del_vars) > 0) {
      out$data_outcome[, del_vars] <- NULL
    }
  }
  
  # --- Update data_switch ---
  K = ifelse(swtrt_control_only, 1, 2)
  for (h in 1:K) {
    out$data_switch[[h]]$data$uid <- NULL
    out$data_switch[[h]]$data$ustratum <- NULL
  }
  
  if (p3 >= 1) {
    # exclude observations after treatment switch
    data1 <- df[!df[[swtrt]] | df[[tstart]] < df[[swtrt_time]], ]
    
    # sort by id and time
    data1 <- data1[order(data1[[id]], data1[[tstart]]), ]
    
    # identify the last obs within each id who switched 
    condition <- !duplicated(data1[[id]], fromLast = TRUE) & 
      data1[[swtrt]] & data1[[tstop]] >= data1[[swtrt_time]]
    
    # reset event and tstop at time of treatment switch
    data1[condition, event] <- 0
    data1[condition, tstop] <- data1[condition, swtrt_time]
    
    tem_vars <- c(swtrt, swtrt_time)
    add_vars <- c(setdiff(vnames3, varnames3), tem_vars)
    avars <- setdiff(add_vars, names(out$data_switch[[1]]$data))
    if (length(avars) > 0) {
      if (logistic_switching_model) {
        for (h in 1:K) {
          out$data_switch[[h]]$data <- 
            merge(out$data_switch[[h]]$data, 
                  data1[, c(id, "tstart", "tstop", avars)], 
                  by = c(id, "tstart", "tstop"), all.x = TRUE, sort = FALSE)
        }
      } else {
        # replicate event times within each subject
        cut <- sort(unique(data1$tstop[data1[[event]] == 1]))
        a1 <- survsplit(data1$tstart, data1$tstop, cut)
        data2 <- data1[a1$row+1,]
        data2$tstart = a1$start
        data2$tstop = a1$end
        
        for (h in 1:K) {
          out$data_switch[[h]]$data <- 
            merge(out$data_switch[[h]]$data, 
                  data2[, c(id, "tstart", "tstop", avars)], 
                  by = c(id, "tstart", "tstop"), all.x = TRUE, sort = FALSE)
        }
      }
    }
    
    del_vars <- setdiff(varnames3, vnames3)
    if (length(del_vars) > 0) {
      for (h in 1:K) {
        out$data_switch[[h]]$data[, del_vars] <- NULL
      }
    }
  }
  
  class(out) <- "ipcw"
  out
}
