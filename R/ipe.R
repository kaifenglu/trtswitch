#' @title Iterative Parameter Estimation (IPE) for Treatment Switching
#' @description Estimates the causal parameter by iteratively fitting an 
#' accelerated failure time (AFT) model to counterfactual 
#' \emph{unswitched} survival times, and derives the adjusted hazard 
#' ratio from the Cox model using counterfactual \emph{unswitched} 
#' survival times based on the estimated causal parameter.
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
#'   * \code{rx}: The proportion of time on active treatment.
#'
#'   * \code{censor_time}: The administrative censoring time. It should
#'     be provided for all subjects including those who had events.
#'
#'   * \code{base_cov}: The baseline covariates (excluding treat).
#'
#' @param id The name of the id variable in the input data.
#' @param stratum The name(s) of the stratum variable(s) in the input data.
#' @param time The name of the time variable in the input data.
#' @param event The name of the event variable in the input data.
#' @param treat The name of the treatment variable in the input data.
#' @param rx The name of the rx variable in the input data.
#' @param censor_time The name of the censor_time variable in the input data.
#' @param base_cov The names of baseline covariates (excluding
#'   treat) in the input data for the causal AFT model and the outcome 
#'   Cox model.
#' @param aft_dist The assumed distribution for time to event for the AFT
#'   model. Options include "exponential", "weibull" (default), 
#'   "loglogistic", and "lognormal".
#' @param strata_main_effect_only Whether to only include the strata main
#'   effects in the AFT model. Defaults to \code{TRUE}, otherwise all
#'   possible strata combinations will be considered in the AFT model.
#' @param low_psi The lower limit of the causal parameter.
#' @param hi_psi The upper limit of the causal parameter.
#' @param treat_modifier The optional sensitivity parameter for the
#'   constant treatment effect assumption.
#' @param recensor Whether to apply recensoring to counterfactual
#'   survival times. Defaults to \code{TRUE}.
#' @param admin_recensor_only Whether to apply recensoring to administrative
#'   censoring times only. Defaults to \code{TRUE}. If \code{FALSE},
#'   recensoring will be applied to the actual censoring times for dropouts.
#' @param autoswitch Whether to exclude recensoring for treatment arms
#'   with no switching. Defaults to \code{TRUE}.
#' @param root_finding Character string specifying the univariate 
#'   root-finding algorithm to use. Options are \code{"brent"} (default)
#'   for Brent's method, or \code{"bisection"} for the bisection method.
#' @param alpha The significance level to calculate confidence intervals.
#' @param ties The method for handling ties in the Cox model, either
#'   "breslow" or "efron" (default).
#' @param tol The desired accuracy (convergence tolerance) for \code{psi} 
#' for the root finding algorithm.
#' @param boot Whether to use bootstrap to obtain the confidence
#'   interval for hazard ratio. Defaults to \code{FALSE}, in which case,
#'   the confidence interval will be constructed to match the log-rank
#'   test p-value.
#' @param n_boot The number of bootstrap samples.
#' @param seed The seed to reproduce the bootstrap results.
#' @param nthreads The number of threads to use in bootstrapping (0 means 
#'   the default RcppParallel behavior)
#'
#' @details Assuming one-way switching from control to treatment, the 
#' hazard ratio and confidence interval under a no-switching scenario 
#' are obtained as follows:
#'
#' * Estimate the causal parameter \eqn{\psi} by iteratively fitting an 
#'   AFT model to the observed survival times for the treatment arm and 
#'   the counterfactual survival times for the control arm: 
#'   \deqn{U_{i,\psi} = T_{C_i} + e^{\psi}T_{E_i}}
#'
#' * Compute counterfactual survival times for control patients using 
#'   the estimated \eqn{\psi}.
#'
#' * Fit a Cox model to the observed survival times for the treatment group 
#'   and the counterfactual survival times for the control group to 
#'   estimate the hazard ratio.
#'
#' * Obtain the confidence interval for the hazard ratio using either 
#'   the ITT log-rank test p-value or bootstrap. When bootstrapping, 
#'   the interval and p-value are derived from a t-distribution 
#'   with \code{n_boot - 1} degrees of freedom.
#'
#' @return A list with the following components:
#'
#' * \code{psi}: The estimated causal parameter.
#'
#' * \code{psi_CI}: The confidence interval for \code{psi}.
#' 
#' * \code{psi_CI_type}: The type of confidence interval for \code{psi},
#'   i.e., "log-rank p-value" or "bootstrap".
#'
#' * \code{pvalue}: The two-sided p-value.
#'
#' * \code{pvalue_type}: The type of two-sided p-value for treatment effect, 
#'   i.e., "log-rank" or "bootstrap".
#'
#' * \code{hr}: The estimated hazard ratio from the Cox model.
#'
#' * \code{hr_CI}: The confidence interval for hazard ratio.
#'
#' * \code{hr_CI_type}: The type of confidence interval for hazard ratio,
#'   either "log-rank p-value" or "bootstrap".
#'
#' * \code{event_summary}: A data frame containing the count and percentage
#'   of deaths and switches by treatment arm.
#'
#' * \code{Sstar}: A data frame containing the counterfactual untreated
#'   survival times and event indicators for each treatment group.
#'   The variables include \code{id}, \code{stratum}, 
#'   \code{"t_star"}, \code{"d_star"}, \code{"treated"}, \code{base_cov}, 
#'   and \code{treat}.
#'   
#' * \code{kmstar}: A data frame containing the Kaplan-Meier estimates
#'   based on the counterfactual untreated survival times by treatment arm.
#'   
#' * \code{data_aft}: The input data for the AFT model for 
#'   estimating \code{psi}. The variables include \code{id}, \code{stratum}, 
#'   \code{"t_star"}, \code{"d_star"}, \code{"treated"}, \code{base_cov}, 
#'   and \code{treat}.
#' 
#' * \code{fit_aft}: The fitted AFT model for estimating \code{psi}.
#' 
#' * \code{res_aft}: The deviance residuals from the fitted AFT model.
#'
#' * \code{data_outcome}: The input data for the outcome Cox model of 
#'   counterfactual unswitched survival times.
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
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @references
#' Michael Branson and John Whitehead.
#' Estimating a treatment effect in survival studies in which patients
#' switch treatment.
#' Statistics in Medicine. 2002;21(17):2449-2463.
#'
#' Ian R White.
#' Letter to the Editor: Estimating treatment effects in randomized
#' trials with treatment switching.
#' Statistics in Medicine. 2006;25(9):1619-1622.
#'
#' @examples
#'
#' library(dplyr)
#'
#' # Example 1: one-way treatment switching (control to active)
#'
#' data <- immdef %>% mutate(rx = 1-xoyrs/progyrs)
#'
#' fit1 <- ipe(
#'   data, id = "id", time = "progyrs", event = "prog", treat = "imm", 
#'   rx = "rx", censor_time = "censyrs", aft_dist = "weibull",
#'   boot = FALSE)
#'
#' fit1
#'
#' # Example 2: two-way treatment switching (illustration only)
#'
#' # the eventual survival time
#' shilong1 <- shilong %>%
#'   arrange(bras.f, id, tstop) %>%
#'   group_by(bras.f, id) %>%
#'   slice(n()) %>%
#'   select(-c("ps", "ttc", "tran"))
#'
#' shilong2 <- shilong1 %>%
#'   mutate(rx = ifelse(co, ifelse(bras.f == "MTA", dco/ady, 
#'                                 1 - dco/ady),
#'                      ifelse(bras.f == "MTA", 1, 0)))
#'
#' fit2 <- ipe(
#'   shilong2, id = "id", time = "tstop", event = "event",
#'   treat = "bras.f", rx = "rx", censor_time = "dcut",
#'   base_cov = c("agerand", "sex.f", "tt_Lnum", "rmh_alea.c",
#'                "pathway.f"),
#'   aft_dist = "weibull", boot = FALSE)
#'
#' fit2
#'
#' @export
ipe <- function(data, id = "id", stratum = "", time = "time", 
                event = "event", treat = "treat", rx = "rx", 
                censor_time = "censor_time", base_cov = "", 
                aft_dist = "weibull", strata_main_effect_only = TRUE, 
                low_psi = -2, hi_psi = 2, treat_modifier = 1,
                recensor = TRUE, admin_recensor_only = TRUE,
                autoswitch = TRUE, root_finding = "brent",
                alpha = 0.05, ties = "efron", tol = 1.0e-6, 
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
  
  for (nm in c(id, time, event, treat, rx, censor_time)) {
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
  elements <- unique(c(id, stratum, time, event, treat, rx, censor_time, base_cov))
  elements <- elements[elements != ""]
  fml_all <- formula(paste("~", paste(elements, collapse = "+")))
  var_all <- all.vars(fml_all)
  rows_ok <- which(complete.cases(df[, var_all, drop = FALSE]))
  if (length(rows_ok) == 0) stop("No complete cases found for the specified variables.")
  df <- df[rows_ok, , drop = FALSE]
  
  # process covariate specifications
  res <- process_cov(base_cov, df)
  df <- res$df
  vnames <- res$vnames
  varnames <- res$varnames
  
  # call the core cpp function
  out <- ipecpp(
    df = df, id = id, stratum = stratum, time = time, 
    event = event, treat = treat, rx = rx, 
    censor_time = censor_time, base_cov = varnames, 
    aft_dist = aft_dist, strata_main_effect_only = strata_main_effect_only,
    low_psi = low_psi, hi_psi = hi_psi, treat_modifier = treat_modifier, 
    recensor = recensor, admin_recensor_only = admin_recensor_only, 
    autoswitch = autoswitch, root_finding = root_finding,
    alpha = alpha, ties = ties, tol = tol, 
    boot = boot, n_boot = n_boot, seed = seed)
  
  if (!out$psimissing) {
    out$Sstar$uid <- NULL
    out$Sstar$ustratum <- NULL
    out$data_aft$uid <- NULL
    out$data_aft$ustratum <- NULL
    out$data_outcome$uid <- NULL
    out$data_outcome$ustratum <- NULL
    
    if (length(vnames) > 0) {
      add_vars <- setdiff(vnames, varnames)
      if (length(add_vars) > 0) {
        for (frame_name in c("Sstar", "data_aft", "data_outcome")) {
          out[[frame_name]] <- merge_append(
            A = out[[frame_name]], B = df,
            by_vars = id, new_vars = add_vars,
            overwrite = FALSE, first_match = FALSE)
        }
      }
      
      del_vars <- setdiff(varnames, vnames)
      if (length(del_vars) > 0) {
        out$Sstar[, del_vars] <- NULL
        out$data_aft[, del_vars] <- NULL
        out$data_outcome[, del_vars] <- NULL
      }
    }
  }
  
  # convert treatment back to a factor variable if needed
  if (is.factor(data[[treat]])) {
    levs <- levels(data[[treat]])
    mf <- function(x) factor(x, levels = c(1,2), labels = levs)
    
    for (nm in c("event_summary", "Sstar", "kmstar", "data_aft", 
                 "data_outcome", "km_outcome")) {
      out[[nm]][[treat]] <- mf(out[[nm]][[treat]])
    }
  }
  
  out$settings <- list(
    data = data, id = id, stratum = stratum, time = time, 
    event = event, treat = treat, rx = rx, 
    censor_time = censor_time, base_cov = base_cov, 
    aft_dist = aft_dist, strata_main_effect_only = strata_main_effect_only,
    low_psi = low_psi, hi_psi = hi_psi, treat_modifier = treat_modifier, 
    recensor = recensor, admin_recensor_only = admin_recensor_only,
    autoswitch = autoswitch, root_finding = root_finding,
    alpha = alpha, ties = ties, tol = tol, 
    boot = boot, n_boot = n_boot, seed = seed
  )
  
  class(out) <- "ipe"
  out
}
