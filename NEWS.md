# trtswitch 0.1.1

- Fixed clang-UBSAN issues.

- Updated description in trtswitch-package.

- Updates to survival.cpp
    - Changed requirement of positive time to nonnegative time for kmest, kmtest, rmest, and rmdiff
    
- Updates to utilities.cpp
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

- Updates to test-rpsftm.R:
    - Removed redefinition of censyrs using progyrs for patients without an event

- Updates to test-ipe.R:
    - Removed redefinition of censyrs using progyrs for patients without an event
    
- Updates to test-tsesimp.R:
    - Replaced censor_time with dcut.
    - Added agerand, sex.f, tt_Lnum, rmh_alea.c, and pathway.f as baseline covariates

- Updates to test-tsegest.R:
    - Changed the default value of admin_recensor_only from FALSE to TRUE.

- Updates to test-ipcw.R:
    - Set relative_time = TRUE and replaced tstop with x0 = tstop - timePFSobs.
    - Changed df_outcome to data_outcome.
    
- Added vignettes:
    - rpsftm.Rmd
    - ipe.Rmd
    - tsesimp.Rmd
    - tsegest.Rmd
    - ipcw.Rmd


# trtswitch 0.1.0

- Initial release.
