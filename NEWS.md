# trtswitch 0.1.6

- fixed an ASAN issue in survQuantile

# trtswitch 0.1.5

- updated the default search interval for psi to [-2,2] at 101 points in rpsftm
- added parameters and output low_psi = -2 and hi_psi = 2 to ipe
- renamed gest to est_psi_tsegest and updated the default search interval for psi to [-2,2] at 101 points in tsegest
- replaced dummy variables with original variables in data_outcome for rpsftm
- replaced dummy variables with original variables in data_outcome and data_aft for ipe
- removed dummy strata variables in data_aft for ipe
- replaced dummy variables with original variables in data_outcome and data_aft for tsesimp
- removed dummy strata variables in data_aft for tsesimp
- renamed time to pps in data_aft for tsesimp
- added pd_time, swtrt_time, and time to data_aft for tsesimp
- added tstart, tstop to data_logis for tsegest
- replaced dummy variables with original variables in data_outcome and data_logis for tsegest
- removed dummy strata variables in data_logis for tsegest
- added tstart, tstop to data_switch for ipcw
- replaced dummy variables with original variables in data_outcome and data_switch for ipcw
- removed dummy strata variables in data_switch for ipcw
- updated liferegr to use better starting values for model parameters
- updated documentations for rpsftm, ipe, tsesimp, tsegest, and ipcw to clarify the variables in output data frames
- added gridsearch parameter to tsegest
- updated unit tes
ts and vignettes to use functions from the survival package
- added maxiter and eps to logisregr, liferegr, and phregr
- added the special case for psilower and psiupper when no root exists when using uniroot finding in rpsftm and tsegest
- updated to shortening survival if patients switched from active treatment to control in tsegestsim
- added subject-level data to tsegestsim

# trtswitch 0.1.4

- Added the keep_censor parameter to the kmest function and added the ncensor variable to the output data frame.
- Added time 0 and the ncensor variable to the output data frame of the kmest function
- Added time 0 and the ncensor variable to the output baseline hazard data frame of the phregr function
- Added id to output data

# trtswitch 0.1.3

- Updated logistic regression and associated programs

# trtswitch 0.1.2

- Added probit and complementary log-log links to logistic regression
- Removed verbatim environment for ggplot in vignettes

# trtswitch 0.1.1

- Fixed clang-UBSAN issues.

- Updated description in trtswitch-package.

- Updates to survival.cpp:
    - Changed requirement of positive time to nonnegative time for kmest, kmtest, rmest, and rmdiff.
    
- Updates to tsegestsimp.cpp:
    - Changed the description of shape1, scale1, shape2, and scale2 parameters.
    
- Updates to utilities.cpp:
    - Moved untreated and unswitched functions here. These functions are shared by rpsftm and ipe.

- Updates to rpsftm:
    - Added treat variable to Sstar.
    - Changed the default value of admin_recensor_only from FALSE to TRUE.
    - Added data_outcome and fit_outcome to output.
    - Added raw stratum and treat information to data_outcome.
    
- Updates to ipe:
    - Added treat variable to Sstar.
    - Changed the default value of admin_recensor_only from FALSE to TRUE.
    - Modified Sstar and kmstar to be consistent with rpsftm.
    - Added data_aft and fit_aft to output.
    - Added data_outcome and fit_outcome to output.
    - Added raw stratum and treat information to data_outcome.

- Updates to tsesimp:
    - Changed the default value of admin_recensor_only from FALSE to TRUE.
    - Added data_aft, fit_aft, data_outcome, and fit_outcome to output.
    - Added psi1hats in bootstrap.
    - Added raw stratum information to data_outcome.
    - Added agerand, sex.f, tt_Lnum, rmh_alea.c, and pathway.f as baseline covariates in the example.

- Updates to tsegest:
    - Changed the default value of admin_recensor_only from FALSE to TRUE.
    - Added the parameter n_eval_z to evaluate the Wald statistics at a sequence of psi values.
    - Added data_switch, km_switch, eval_z, data_nullcox, fit_nullcox, data_logis, fit_logis, data_outcome, and fit_outcome to output.
    - Included data_switch, km_switch, eval_z, data_nullcox, fit_nullcox, data_logis, and fit_logis in the analysis_switch list.
    - Added raw stratum information to data_outcome.

- Updates to ipcw: 
    - Added relative_time parameter.
    - Added data_switch to output.
    - Added ns_df to output.
    - Renamed df_outcome to data_outcome.
    - Added raw stratum information to data_outcome.
    
- Added vignettes:
    - rpsftm.Rmd
    - ipe.Rmd
    - tsesimp.Rmd
    - tsegest.Rmd
    - ipcw.Rmd


# trtswitch 0.1.0

- Initial release.
