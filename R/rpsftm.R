#' @title Rank Preserving Structural Failure Time Model (RPSFTM) for
#' Treatment Switching
#' @description Obtains the causal parameter estimate from the 
#' log-rank test and the hazard ratio estimate from the Cox model
#' to adjust for treatment switching.
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
#' @param low_psi The lower limit of the causal parameter.
#' @param hi_psi The upper limit of the causal parameter.
#' @param n_eval_z The number of points between \code{low_psi} and 
#'   \code{hi_psi} (inclusive) at which to evaluate the log-rank 
#'   Z-statistics.
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
#'   parameter \code{psi}. Defaults to \code{FALSE}, in which case, a root
#'   finding algorithm will be used.
#' @param alpha The significance level to calculate confidence intervals.
#' @param ties The method for handling ties in the Cox model, either
#'   "breslow" or "efron" (default).
#' @param tol The desired accuracy (convergence tolerance) for \code{psi}.
#' @param boot Whether to use bootstrap to obtain the confidence
#'   interval for hazard ratio. Defaults to \code{FALSE}, in which case,
#'   the confidence interval will be constructed to match the log-rank
#'   test p-value.
#' @param n_boot The number of bootstrap samples.
#' @param seed The seed to reproduce the bootstrap results. The default is 
#'   missing, in which case, the seed from the environment will be used.
#'
#' @details We use the following steps to obtain the hazard ratio estimate
#' and confidence interval had there been no treatment switching:
#'
#' * Use RPSFTM to estimate the causal parameter \eqn{\psi} based on the 
#'   log-rank test for counterfactual untreated survival times for 
#'   both arms: \deqn{U_{i,\psi} = T_{C_i} +  e^{\psi}T_{E_i}}
#'
#' * Fit the Cox proportional hazards model to the observed survival times
#'   for the experimental group and the counterfactual survival times
#'   for the control group to obtain the hazard ratio estimate.
#'
#' * Use either the log-rank test p-value for the intention-to-treat (ITT) 
#'   analysis or bootstrap to construct the confidence interval for 
#'   hazard ratio. If bootstrapping is used, the confidence interval 
#'   and corresponding p-value are calculated based on a t-distribution with 
#'   \code{n_boot - 1} degrees of freedom. 
#'
#' @return A list with the following components:
#'
#' * \code{psi}: The estimated causal parameter.
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
#'   either "log-rank p-value" or "bootstrap".
#'
#' * \code{eval_z}: A data frame containing the log-rank test Z-statistics
#'   evaluated at a sequence of \code{psi} values. Used to plot and check
#'   if the range of \code{psi} values to search for the solution and
#'   limits of confidence interval of \code{psi} need be modified.
#'
#' * \code{Sstar}: A data frame containing the counterfactual untreated
#'   survival times and event indicators for each treatment group.
#'
#' * \code{kmstar}: A data frame containing the Kaplan-Meier estimates
#'   based on the counterfactual untreated survival times by treatment arm.
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
#'       \code{hi_psi} (inclusive) at which to evaluate the log-rank 
#'       Z-statistics.
#'
#'     - \code{treat_modifier}: The sensitivity parameter for the
#'       constant treatment effect assumption.
#'
#'     - \code{recensor}: Whether to apply recensoring to counterfactual
#'       survival times.
#'
#'     - \code{admin_recensor_only}: Whether to apply recensoring to
#'       administrative censoring times only.
#'
#'     - \code{autoswitch}: Whether to exclude recensoring for treatment 
#'       arms with no switching.
#'
#'     - \code{gridsearch}: Whether to use grid search to estimate the 
#'       causal parameter \code{psi}.
#'
#'     - \code{alpha}: The significance level to calculate confidence
#'       intervals.
#'
#'     - \code{ties}: The method for handling ties in the Cox model.
#'
#'     - \code{tol}: The desired accuracy (convergence tolerance) 
#'       for \code{psi}.
#'
#'     - \code{boot}: Whether to use bootstrap to obtain the confidence
#'       interval for hazard ratio.
#'
#'     - \code{n_boot}: The number of bootstrap samples.
#'
#'     - \code{seed}: The seed to reproduce the bootstrap results.
#'
#' * \code{hr_boots}: The bootstrap hazard ratio estimates if \code{boot} is
#'   \code{TRUE}.
#'
#' * \code{psi_boots}: The bootstrap \code{psi} estimates if \code{boot} is 
#'   \code{TRUE}.
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
#' Statistics in Medicine. 1999;18:2617-2634.
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
#' c(fit1$hr, fit1$hr_CI)
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
#' c(fit2$hr, fit2$hr_CI)
#'
#' @export
rpsftm <- function(data, id = "id", stratum = "", 
                   time = "time", event = "event",
                   treat = "treat", rx = "rx", censor_time = "censor_time",
                   base_cov = "", low_psi = -2, hi_psi = 2,
                   n_eval_z = 101, treat_modifier = 1,
                   recensor = TRUE, admin_recensor_only = TRUE,
                   autoswitch = TRUE, gridsearch = FALSE,
                   alpha = 0.05, ties = "efron", tol = 1.0e-6,
                   boot = FALSE, n_boot = 1000, seed = NA) {

  rownames(data) = NULL

  elements = c(stratum, time, event, treat, rx, censor_time)
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

  out <- rpsftmcpp(
    data = df, id = id, stratum = stratum, 
    time = time, event = event,
    treat = treat, rx = rx, censor_time = censor_time,
    base_cov = varnames, low_psi = low_psi, hi_psi = hi_psi,
    n_eval_z = n_eval_z, treat_modifier = treat_modifier,
    recensor = recensor, admin_recensor_only = admin_recensor_only,
    autoswitch = autoswitch, gridsearch = gridsearch, 
    alpha = alpha, ties = ties, tol = tol, 
    boot = boot, n_boot = n_boot, seed = seed)
  
  out$Sstar$uid <- NULL
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
  
  out
}
