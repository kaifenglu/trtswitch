# trtswitch 0.1.1

- Fixed clang-UBSAN issues.

- Updated description in trtswitch-package.

- Updates to rpsftm:
    - Added treat variable to Sstar.
    - Changed the default value of admin_recensor_only from FALSE to TRUE.
    - Added data_outcome and fit_outcome to output.

- Updates to ipe:
    - Added treat variable to Sstar.
    - Changed the default value of admin_recensor_only from FALSE to TRUE.
    - Modified Sstar and kmstar to be consistent with rpsftm.
    - Added data_outcome and fit_outcome to output.

- Updates to tsesimp:
    - Changed the default value of admin_recensor_only from FALSE to TRUE.
    - Added data_aft, fit_aft, data_outcome, and fit_outcome to output.
    - Added psi1hats in bootstrap.

- Updates to tsegest:
    - Changed the default value of admin_recensor_only from FALSE to TRUE.
    - Added data_nullcox, fit_nullcox, data_logis, fit_logis, data_outcome, and fit_outcome to output.

- Updates to ipcw: 
    - Added relative_time parameter.
    - Added data_switch to output.
    - Added ns_df to output.
    - Renamed df_outcome to data_outcome.

- Updates to test-tsesimp.R:
    - Replaced censor_time with dcut.

- Updates to test-tsegest.R:
    - Changed the default value of admin_recensor_only from FALSE to TRUE.

- Updates to test-ipcw.R:
    - Set relative_time = TRUE and replaced tstop with x0 = tstop - timePFSobs.
    - Changed df_outcome to data_outcome.
    
- Added vignettes:
    - rpsftm.Rmd
    - ipe.Rmd
    - tsesimp.Rmd


# trtswitch 0.1.0

- Initial release.
