#' @title Rank Preserving Structural Failure Time Model (RPSFTM) for
#' Treatment Switching
#' @description Estimates the causal treatment effect parameter using 
#' g-estimation based on the log-rank test, Cox model, or parametric 
#' survival/accelerated failure time (AFT) model. The method uses 
#' counterfactual \emph{untreated} survival times to estimate the 
#' causal parameter and derives the adjusted hazard ratio from the 
#' Cox model using counterfactual \emph{unswitched} survival times.
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
#'   treat) in the input data for the outcome Cox model. 
#'   These covariates will also be used in the Cox model for estimating 
#'   \code{psi} when \code{psi_test = "phreg"} and in the parametric
#'   survival regression/AFT model for 
#'   estimating \code{psi} when \code{psi_test = "lifereg"}.
#' @param psi_test The survival function to calculate the Z-statistic, e.g., 
#'   "logrank" (default), "phreg", or "lifereg".
#' @param aft_dist The assumed distribution for time to event for the AFT 
#'   model when \code{psi_test = "lifereg"}. Options include "exponential", 
#'   "weibull" (default), "loglogistic", and "lognormal".
#' @param strata_main_effect_only Whether to only include the strata main
#'   effects in the AFT model. Defaults to \code{TRUE}, otherwise all
#'   possible strata combinations will be considered in the AFT model.
#' @param low_psi The lower limit of the causal parameter.
#' @param hi_psi The upper limit of the causal parameter.
#' @param n_eval_z The number of points between \code{low_psi} and 
#'   \code{hi_psi} (inclusive) at which to evaluate the Z-statistics.
#' @param treat_modifier The optional sensitivity parameter for the
#'   constant treatment effect assumption.
#' @param recensor Whether to apply recensoring to counterfactual
#'   survival times. Defaults to \code{TRUE}.
#' @param admin_recensor_only Whether to apply recensoring to administrative
#'   censoring times only. Defaults to \code{TRUE}. If \code{FALSE},
#'   recensoring will be applied to the actual censoring times for dropouts.
#' @param autoswitch Whether to exclude recensoring for treatment arms
#'   with no switching. Defaults to \code{TRUE}.
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
#'   for the root finding algorithm.
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
#' * Estimate the causal parameter \eqn{\psi} using g-estimation based on 
#'   the log-rank test (default), Cox model, or parametric survival/AFT 
#'   model, using counterfactual \emph{untreated} survival times for 
#'   both arms: 
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
#' If grid search is used to estimate \eqn{\psi}, the estimated \eqn{\psi} 
#' is the one with the smallest absolute value among those at which 
#' the Z-statistic is zero based on linear interpolation. 
#' If root finding is used, the estimated \eqn{\psi} is
#' the solution to the equation where the Z-statistic is zero.
#'
#' @return A list with the following components:
#'
#' * \code{psi}: The estimated causal parameter.
#' 
#' * \code{psi_roots}: Vector of \code{psi} values at which the Z-statistic 
#'   is zero, identified using grid search and linear interpolation.
#'
#' * \code{psi_CI}: The confidence interval for \code{psi}.
#' 
#' * \code{psi_CI_type}: The type of confidence interval for \code{psi},
#'   i.e., "grid search", "root finding", or "bootstrap".
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
#' * \code{eval_z}: A data frame containing the Z-statistics for treatment
#'   effect evaluated at a sequence of \code{psi} values. Used to plot and 
#'   check if the range of \code{psi} values to search for the solution and
#'   limits of confidence interval of \code{psi} need be modified.
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
#' James M. Robins and Anastasios A. Tsiatis.
#' Correcting for non-compliance in randomized trials using rank preserving
#' structural failure time models.
#' Communications in Statistics. 1991;20(8):2609-2631.
#'
#' Ian R. White, Adbel G. Babiker, Sarah Walker, and Janet H. Darbyshire.
#' Randomization-based methods for correcting for treatment changes:
#' Examples from the CONCORDE trial.
#' Statistics in Medicine. 1999;18(19):2617-2634.
#'
#' @examples
#'
#' library(dplyr)
#'
#' # Example 1: one-way treatment switching (control to active)
#'
#' data <- immdef %>% mutate(rx = 1-xoyrs/progyrs)
#'
#' fit1 <- rpsftm(
#'   data, id = "id", time = "progyrs", event = "prog", treat = "imm",
#'   rx = "rx", censor_time = "censyrs", boot = FALSE)
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
#' fit2 <- rpsftm(
#'   shilong2, id = "id", time = "tstop", event = "event",
#'   treat = "bras.f", rx = "rx", censor_time = "dcut",
#'   base_cov = c("agerand", "sex.f", "tt_Lnum", "rmh_alea.c",
#'                "pathway.f"),
#'   low_psi = -3, hi_psi = 3, boot = FALSE)
#'
#' fit2
#'
#' @export
rpsftm <- function(data, id = "id", stratum = "", time = "time", 
                   event = "event", treat = "treat", rx = "rx", 
                   censor_time = "censor_time", base_cov = "", 
                   psi_test = "logrank", aft_dist = "weibull", 
                   strata_main_effect_only = TRUE, 
                   low_psi = -2, hi_psi = 2, n_eval_z = 101, 
                   treat_modifier = 1, recensor = TRUE, 
                   admin_recensor_only = TRUE, autoswitch = TRUE, 
                   gridsearch = TRUE, root_finding = "brent",
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
  out <- rpsftmcpp(
    df = df, id = id, stratum = stratum, time = time, 
    event = event, treat = treat, rx = rx, 
    censor_time = censor_time, base_cov = varnames, 
    psi_test = psi_test, aft_dist = aft_dist,
    strata_main_effect_only = strata_main_effect_only, 
    low_psi = low_psi, hi_psi = hi_psi, n_eval_z = n_eval_z, 
    treat_modifier = treat_modifier, recensor = recensor, 
    admin_recensor_only = admin_recensor_only, autoswitch = autoswitch, 
    gridsearch = gridsearch, root_finding = root_finding,
    alpha = alpha, ties = ties, tol = tol, 
    boot = boot, n_boot = n_boot, seed = seed)
  
  if (!out$psimissing) {
    out$Sstar$uid <- NULL
    out$Sstar$ustratum <- NULL
    out$data_outcome$uid <- NULL
    out$data_outcome$ustratum <- NULL
    
    if (length(vnames) > 0) {
      add_vars <- setdiff(vnames, varnames)
      if (length(add_vars) > 0) {
        for (frame_name in c("Sstar", "data_outcome")) {
          out[[frame_name]] <- merge_append(
            A = out[[frame_name]],   B = df,
            by_vars = id, new_vars = add_vars,
            overwrite = FALSE, first_match = TRUE)
        }
      }
      
      del_vars <- setdiff(varnames, vnames)
      if (length(del_vars) > 0) {
        out$Sstar[, del_vars] <- NULL
        out$data_outcome[, del_vars] <- NULL
      }
    }
  }
  
  if (is.factor(data[[treat]])) {
    levs <- levels(data[[treat]])
    mf <- function(x) factor(x, levels = c(1,2), labels = levs)
    
    for (nm in c("event_summary", "Sstar", "kmstar", "data_outcome", 
                 "km_outcome")) {
      out[[nm]][[treat]] <- mf(out[[nm]][[treat]])
    }
  }
  
  out$settings <- list(
    data = data, id = id, stratum = stratum, time = time, 
    event = event, treat = treat, rx = rx, 
    censor_time = censor_time, base_cov = base_cov,
    psi_test = psi_test, aft_dist = aft_dist,
    strata_main_effect_only = strata_main_effect_only,
    low_psi = low_psi, hi_psi = hi_psi, n_eval_z = n_eval_z,
    treat_modifier = treat_modifier, recensor = recensor,
    admin_recensor_only = admin_recensor_only, autoswitch = autoswitch, 
    gridsearch = gridsearch, root_finding = root_finding,
    alpha = alpha, ties = ties, tol = tol, 
    boot = boot, n_boot = n_boot, seed = seed
  )
  
  class(out) <- "rpsftm"
  out
}
