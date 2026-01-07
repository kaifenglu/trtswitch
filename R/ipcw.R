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
#' @param seed The seed to reproduce the bootstrap results.
#' @param nthreads The number of threads to use in bootstrapping (0 means 
#'   the default RcppParallel behavior)
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
#'   of deaths and switches by treatment arm.
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
#' * \code{weight_summary}: A data frame summarizing the weights by
#'   treatment arm.
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
#' * \code{settings}: A list containing the input parameter values.
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
#'   tdxo = TRUE, coxo = TRUE, allocation1 = 1, allocation2 = 1,
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
#' fit1 
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
#' fit2
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
                 boot = FALSE, n_boot = 1000, seed = 0,
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
  
  for (nm in c(id, tstart, tstop, event, treat, swtrt, swtrt_time)) {
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
  elements = unique(c(id, stratum, tstart, tstop, event, treat, swtrt))
  elements = elements[elements != ""]
  fml_all <- formula(paste("~", paste(elements, collapse = "+")))
  var_all <- all.vars(fml_all)
  rows_ok <- which(complete.cases(df[, var_all, drop = FALSE]))
  if (length(rows_ok) == 0) 
    stop("No complete cases found for the specified variables.")
  df <- df[rows_ok, , drop = FALSE]

  # process covariate specifications
  res1 <- process_cov(base_cov, df)
  df <- res1$df
  vnames    <- res1$vnames
  varnames  <- res1$varnames
  
  res2 <- process_cov(numerator, df)
  df <- res2$df
  vnames2    <- res2$vnames
  varnames2  <- res2$varnames
  
  res3 <- process_cov(denominator, df)
  df <- res3$df
  vnames3    <- res3$vnames
  varnames3  <- res3$varnames
  
  # call the core cpp function
  out <- ipcwcpp(
    df = df, id = id, stratum = stratum, tstart = tstart,
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
  
  if (length(vnames) > 0) {
    add_vars <- setdiff(vnames, varnames)
    if (length(add_vars) > 0) {
      out$data_outcome <- merge_append(
        A = out$data_outcome, B = dfu, 
        by_vars = id, new_vars = add_vars,
        overwrite = FALSE, first_match = FALSE)
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
  
  if (length(vnames3) > 0) {
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
          out$data_switch[[h]]$data <- merge_append(
            A = out$data_switch[[h]]$data, B = data1, 
            by_vars = c(id, "tstart", "tstop"), new_vars = avars,
            overwrite = FALSE, first_match = FALSE)
        }
      } else {
        # replicate event times within each subject
        cut <- sort(unique(data1$tstop[data1[[event]] == 1]))
        a1 <- survsplit(data1$tstart, data1$tstop, cut)
        data2 <- data1[a1$row  + 1, ]
        data2$tstart = a1$start
        data2$tstop = a1$end
        for (h in 1:K) {
          out$data_switch[[h]]$data <- merge_append(
            A = out$data_switch[[h]]$data, B = data2, 
            by_vars = c(id, "tstart", "tstop"), new_vars = avars,
            overwrite = FALSE, first_match = FALSE)
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
  
  # convert treatment back to a factor variable if needed
  if (is.factor(data[[treat]])) {
    levs = levels(data[[treat]])
    mf <- function(x) factor(x, levels = c(1,2), labels = levs)
    
    # apply mf to a set of data.frames with a column named `treat`
    for (nm in c("event_summary", "weight_summary", "data_outcome", 
                 "km_outcome")) {
      out[[nm]][[treat]] <- mf(out[[nm]][[treat]])
    }
    
    # and for the list-of-lists
    out$data_switch <- lapply(out$data_switch, function(x) { 
      x[[treat]] <- mf(x[[treat]]); x 
    })
  }
  
  out$settings <- list(
    data = data, id = id, stratum = stratum, tstart = tstart, 
    tstop = tstop, event = event, treat = treat, swtrt = swtrt, 
    swtrt_time = swtrt_time, base_cov = base_cov, 
    numerator = numerator, denominator = denominator,
    logistic_switching_model = logistic_switching_model,
    strata_main_effect_only = strata_main_effect_only,
    ns_df = ns_df, firth = firth, flic = flic,
    stabilized_weights = stabilized_weights, 
    trunc = trunc, trunc_upper_only = trunc_upper_only,
    swtrt_control_only = swtrt_control_only, 
    alpha = alpha, ties = ties, boot = boot, 
    n_boot = n_boot, seed = seed
  )
  
  class(out) <- "ipcw"
  out
}
