#' @name trtswitch-package
#' @aliases trtswitch-package
#' @title Treatment Switching
#'
#' @description Implements rank preserving structural failure time model 
#' (RPSFTM), iterative parameter estimation (IPE), inverse probability of 
#' censoring weights (IPCW), marginal structural model (MSM), simple 
#' two-stage estimation (TSEsimp), and improved two-stage estimation with 
#' g-estimation (TSEgest) methods for treatment switching in randomized 
#' clinical trials.
#'
#' @details To enable bootstrapping of the parameter estimates, we implements
#' many standard survival analysis methods in C++. These include but are not 
#' limited to Kaplan-Meier estimates of the survival curves, log-rank tests, 
#' accelerated failure time (AFT) models, and Cox proportional hazards 
#' models. 
#' 
#' All treatment switching adjustment methods allow treatment switching
#' in both treatment arms, stratification and covariates adjustment. 
#' For the AFT models, stratification factors are included as covariates 
#' (main effects only or all-way interactions) because SAS PROC LIFEREG 
#' does not have the strata statement. The RPSFTM, IPE and TSE methods 
#' can be used with or without recensoring. The IPCW and MSM methods can 
#' produce stabilized and truncated weights. 
#' 
#' The \code{treat} variable adopts a treatment-before-control order,  
#' except with 1/0 or TRUE/FALSE coding.
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
#' Michael Branson and John Whitehead.
#' Estimating a treatment effect in survival studies in which patients
#' switch treatment.
#' Statistics in Medicine. 2002;21:2449-2463.
#'
#' Ian R White.
#' Letter to the Editor: Estimating treatment effects in randomized
#' trials with treatment switching.
#' Statistics in Medicine. 2006;25:1619-1622.
#'
#' James M. Robins and Dianne M. Finkelstein.
#' Correcting for noncompliance and dependent censoring in an AIDS clinical
#' trial with inverse probability of censoring weighted (IPCW) log-rank tests.
#' Biometrics. 2000;56:779-788.
#'
#' Nicholas R. Latimer, Keith R. Abrams, Paul C. Lambert, Michael K. Crowther, 
#' Allan J. Wailoo, Jonathan P. Morden, Ron L. Akehurst, and 
#' Michael J. Campbell.
#' Adjusting for treatment switching in randomised controlled trials - A
#' simulation study and a simplified two-stage method.
#' Statistical Methods in Medical Research. 2017;26(2):724-751.
#' 
#' Nicholas R. Latimer, Ian R. White, Kate Tilling, and Ulrike Siebert. 
#' Improved two-stage estimation to adjust for treatment switching in 
#' randomised trials: g-estimation to address time-dependent confounding. 
#' Statistical Methods in Medical Research. 2020;29(10):2900-2918.
#' 
#' James M. Robins, Miguel Angel Hernan, and Babette Brumback. 
#' Marginal structural models and causal inference in epidemiology. 
#' Epidemiology. 2000;11(5):550-560.
#' 
#' Miguel Angel Hernan, Babette Brumback, and James M. Robins. 
#' Marginal structural modesl to estimate the causual effect of 
#' zidovudine on the survival of HIV-positive men. 
#' Epidemiology. 2000;11(5):561-570.
#'
#' @useDynLib trtswitch, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom stats formula makepredictcall model.frame model.matrix 
#' na.omit na.pass pchisq plogis qnorm rbinom rexp terms
#' @importFrom rlang .data
#' @importFrom data.table data.table dcast  fifelse nafill
#'   setDT setnames setorderv shift .N .SD :=
#' @importFrom ggplot2 aes annotate coord_cartesian element_blank 
#'  element_text facet_wrap geom_step geom_text geom_histogram geom_line
#'  geom_hline geom_vline ggplot ggplot_build labs margin 
#'  scale_x_continuous scale_y_continuous scale_y_discrete 
#'  theme theme_bw theme_minimal
#' @importFrom cowplot align_plots plot_grid
#'
NULL

