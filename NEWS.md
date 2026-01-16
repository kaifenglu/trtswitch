# trtswitch 0.2.4

* updated the concat_flatmatrix function for empty flat matrices
* utilized column-major storage of flat matrices for cholesky decomposition
* simplified the bygroup function

# trtswitch 0.2.3

* added assess_phregr to assess the proportional hazards assumption for the Cox model using the supremum test of Lin et al. (1993)
* added an internal data pbc for testing assess_phregr
* updated the censoring of post switch survival times in ipcw and msm to be consistent with the swtrt_control_only parameter value
* added zph_phregr to compute the Schoenfeld residuals and test statistics for the proportional hazards assumption for the Cox model
* refactored the code in pure c++ and RcppParallel for core computations to improve computational speed

# trtswitch 0.2.2

* added a check for the class of object for residuals\_lieferegr.R, residuals\_phregr.R, and survfit\_phregr.R
* added psi\_roots for the vector of \code{psi} values at which the Z-statistic is zero, identified using grid search and linear interpolation for rpsftm.
* added psi\_roots and psi\_trt\_roots for the vector of psi values at which the Z-statistic is zero, identified using grid search and linear interpolation for tsegest
* separated NumericVector loghrs = log(hrhats\[ok]) into two statements: NumericVector subset\_hrhats = hrhats\[ok]; NumericVector loghrs = log(subset\_hrhats); to avoid ASAN errors for ipcw.cpp, ipe.cpp, msm.cpp, rpsftm.cpp, tsegest.cpp, and tsesimp.cpp
* updated the lrtest call to add the value for the new weight\_readj parameter for ipcw.cpp, ipe.cpp, msm.cpp, rpsftm.cpp, tsegest.cpp, and tsesimp.cpp
* added the weight\_readj parameter to lrtest to indicate whether the weight variable at each event time will be readjusted to be proportional to the number at risk by treatment group in survival\_analysis.
* updated the getpsiest function to output multiple roots depending on the new direction parameter value in utilities.cpp and utilities.h
* updated rpsftm and tsegest to use gridsearch by default
* added print method for rpsftm, ipe, tsesimp, tsegest, ipcw, and msm
* added class attribute for rpsftm, ipe, tsesimp, tsegest, ipcw, and msm objects
* added km\_outcome and lr\_outcome for rpsftm, ipe, tsesimp, tsegest, ipcw, and msm objects
* clarified the point estimate of psi when there are multiple roots and grid search is used
* initialized weights to 1.0 for the Cox switching model in ipcw
* added event summary and weight summary to output for ipcw and msm
* added event summary to output for rpsftm, ipe, tsesimp and tsegest
* added the plot method for rpsftm, ipe, tsesimp, tsegest, ipcw, and msm
* use settings to store input parameter values for logisregr, liferegr, phregr, rpsftm, ipe, tsesimp, tsegest, ipcw, and msm
* converted treatment back to a factor variable for output data sets in rpsftm, ipe, tsesimp, tsegest, ipcw, and msm
* add res\_aft to ipe and tsesimp



# trtswitch 0.2.1

* fixed the bootstrap procedures for IPCW, MSM, and TSEgest by reconstructing idx, idx1, treatn1, and stratumn1
* clarified the distinction between IPCW using pooled logistic regression and IPCW using a Cox model with time-dependent covariates
* excluded existing variables in out$data\_switch\[\[h]]$data from add\_vars for IPCW and MSM
* excluded existing variables in out$analysis\_switch$data\_logis\[\[h]]$data from add\_vars for TSEgest
* excluded existing variables in out$data\_aft\[\[h]]$data from add\_vars for TSEsimp
* replaced chained sorting using logical OR with multiple if–return statements to ensure stable sorting in survival\_analysis
* replaced chained sorting using logical OR with std::tie in IPCW, IPE, MSM, RPSFTM, TSEsimp, and TSEgest
* added id = "id" in test-ipe.R and test-rpsftm.R
* combined functions f and g into a single function in test-rpsftm.R
* added bootstrap testing for test-tsesimp.R
* removed the target argument from the definition of function f in test-tsegest.R
* added an option to use a Cox model with robust variance for constructing confidence intervals for the hazard ratio in IPCW and MSM
* removed robust variance from the switching model in IPCW and MSM to improve computational speed
* replaced geom\_density with geom\_histogram in the IPCW and MSM vignettes
* removed natural spline basis functions for time from the pooled logistic regression model equation in TSEgest and incorporated these effects into time-dependent confounders in the TSEgest vignette
* replaced event = "died" with event = "event" in bootstrap testing in the TSEgest vignette
* removed the analysis\_switch level from the output of TSEgest and moved its components to the top level of the output list
* updated getpsiest to return only the first root
* updated IPCW to apply weight truncation within each treatment group after computing weights for data\_outcome
* updated TSEsimp and TSEgest to use the switch time as the progression time when switching occurs before or in the absence of recorded progression
* updated kmest to add the weight parameter for adjusted Kaplan-Meier estimate by Xie \& Liu (2005)
* updated lrtest to add the weight parameter for adjusted log-rank test by Xie \& Liu (2005, 2011)
* updated ipcw and msm to add Kaplan-Meier estimates of the survival functions and the log-rank test for the treatment effect based on the weighted outcome data truncated at time of treatment switching

# trtswitch 0.2.0

* added MSM to package description
* added the plannedTime argument to the tssim function so that followupTime only refers to the follow-up time in a fixed follow-up design
* added an extra record for a subject to the tssim output dataset if the observed survival time exceeds the number of regular treatment cycles times days per cycle
* removed the getAccrualDurationFromN utility function
* changed the default value of the boot parameter to FALSE for ipcw and msm
* added the ns\_df parameter to the tsegest function
* added the preptdc function to prepare survival data with time-dependent covariates
* expanded the findInterval3 function definition to mimic R's findInterval() behavior with rightmost\_closed, all\_inside, and left\_open options
* updated the survfit\_phregr function to only include time points between the first tstart and the last tstop for each subject with counting process style of input
* added the exclusion of observations with missing values for survfit\_phregr and residuals\_phregr
* added the exclusion of observations with missing values for rpsftm, ipe, tsesimp, tsegest, ipcw, and msm
* changed to sort the input data by treatment, ID, and time variables for IPCW, MSM, and TSEgest, and by treatment, stratum, ID, and time variables for bootstrap samples
* changed the default value of treat\_alt\_interaction from FALSE to TRUE for MSM
* renamed Y to event for tssim output
* added the bisect utility function for root finding
* added the root\_finding parameter to rpsftm, ipe, tsegest to allow the user to choose between brent and bisect for root finding

# trtswitch 0.1.9

* added fail\_boots\_data to include the failed bootstrap sample data
* used strata variable names and values for dummy variables in logistic regression and life regression models
* updated the tssim function to allow both fixed and variable follow-up designs
* used stratified bootstrapping for the rpsftm, ipe, tsesimp, tsegest, ipcw, and msm functions
* allowed event, pd, and swtrt to take real values 1 or 0 in rpsftm, ipe, tsesimp, tsegest, ipcw, and msm functions
* added martingale residuals for lifereg for event or right-censored survival data
* added time-dependent covariate visit7on and used cattdc instead of catlag for tsegest sample call
* added interval splitting if treatment switching occurs between visits for ipcw, msm, and tsegest
* added a small number to os\_time if needed to ensure pd\_time is less than os\_time for tsesimp and tsegest
* added getpsiest and getpsiend utility functions for psi estimation of rpsftm, ipe, and tsegest
* updated the condition to include post progression data up to switching in tsegest
* updated the condition for time-dependent crossover indicator in tsegest
* updated the algorithm for generating timePFSobs in tsegestsim
* removed catlag and xotime\_upper from the output data set of tsegestsim
* applied weights to the next interval for ipcw and msm using logistic regression switching model
* updated tssim to use Z and A to denote the alternative therapy status at the start instead of end of interval
* updated logisregr, liferegr and phregr to allow zero events to run with trivial output

# trtswitch 0.1.8

* eliminated redundant code for the logistic regression switching model for ipcw
* added the marginal structural model (msm) method
* added na.action = na.pass for model frame construction involving covariates for all methods
* added the init parameter and the fail flag to the output to logisregr, liferegr, and phregr
* replaced the survreg initial value method with the OLS initial value method for liferegr
* added the fail flag to output of rpsftm, ipe, tsesimp, tsegest, ipcw, and msm

# trtswitch 0.1.7

* updated survival\_analysis to ignore intervals not at risk within each stratum without creating non overlapping times across strata
* updated documentation for the survsplit utility function
* updated tsegestsim to use the standard definition of weibull scale
* removed bc from logisregr
* removed swtrt\_time\_upper from tsegest
* removed subject-level adsl data from tsegestsim output
* removed swtrt\_time\_lower, swtrt\_time\_upper, and relative\_time from ipcw
* removed the robust option for logistic regression treatment switching model in ipcw
* added match3 utility function to match on id and value
* added residuals\_liferegr for residuals from parameteric regression models for failure time data
* added the psi\_test, aft\_dist, and strata\_main\_effect\_only parameters to rpsftm to allow parametric regression and cox regression models to estimate psi
* added baseline covariates to the output Sstar data set of the rpsftm and ipe functions
* added a check for switching time before progression time in tsesimp and tsegest
* added a check for switching time or progression time less than offset in tsesimp and tsegest
* added timeOS and died variables for subject-level data and used event for death status at the end of each time interval for tsegestsim
* added tssim for treatment switching data simulation

# trtswitch 0.1.6

* fixed an ASAN issue in survQuantile

# trtswitch 0.1.5

* updated the default search interval for psi to \[-2,2] at 101 points in rpsftm
* added parameters and output low\_psi = -2 and hi\_psi = 2 to ipe
* renamed gest to est\_psi\_tsegest and updated the default search interval for psi to \[-2,2] at 101 points in tsegest
* replaced dummy variables with original variables in data\_outcome for rpsftm
* replaced dummy variables with original variables in data\_outcome and data\_aft for ipe
* removed dummy strata variables in data\_aft for ipe
* replaced dummy variables with original variables in data\_outcome and data\_aft for tsesimp
* removed dummy strata variables in data\_aft for tsesimp
* renamed time to pps in data\_aft for tsesimp
* added pd\_time, swtrt\_time, and time to data\_aft for tsesimp
* added tstart, tstop to data\_logis for tsegest
* replaced dummy variables with original variables in data\_outcome and data\_logis for tsegest
* removed dummy strata variables in data\_logis for tsegest
* added tstart, tstop to data\_switch for ipcw
* replaced dummy variables with original variables in data\_outcome and data\_switch for ipcw
* removed dummy strata variables in data\_switch for ipcw
* updated liferegr to use better starting values for model parameters
* updated documentations for rpsftm, ipe, tsesimp, tsegest, and ipcw to clarify the variables in output data frames
* added gridsearch parameter to tsegest
* updated unit tests and vignettes to use functions from the survival package
* added maxiter and eps to logisregr, liferegr, and phregr
* added the special case for psilower and psiupper when no root exists when using uniroot finding in rpsftm and tsegest
* updated to shortening survival if patients switched from active treatment to control in tsegestsim
* added subject-level data to tsegestsim

# trtswitch 0.1.4

* Added the keep\_censor parameter to the kmest function and added the ncensor variable to the output data frame.
* Added time 0 and the ncensor variable to the output data frame of the kmest function
* Added time 0 and the ncensor variable to the output baseline hazard data frame of the phregr function
* Added id to output data

# trtswitch 0.1.3

* Updated logistic regression and associated programs

# trtswitch 0.1.2

* Added probit and complementary log-log links to logistic regression
* Removed verbatim environment for ggplot in vignettes

# trtswitch 0.1.1

* Fixed clang-UBSAN issues.
* Updated description in trtswitch-package.
* Updates to survival.cpp:

  * Changed requirement of positive time to nonnegative time for kmest, kmtest, rmest, and rmdiff.

* Updates to tsegestsimp.cpp:

  * Changed the description of shape1, scale1, shape2, and scale2 parameters.

* Updates to utilities.cpp:

  * Moved untreated and unswitched functions here. These functions are shared by rpsftm and ipe.

* Updates to rpsftm:

  * Added treat variable to Sstar.
  * Changed the default value of admin\_recensor\_only from FALSE to TRUE.
  * Added data\_outcome and fit\_outcome to output.
  * Added raw stratum and treat information to data\_outcome.

* Updates to ipe:

  * Added treat variable to Sstar.
  * Changed the default value of admin\_recensor\_only from FALSE to TRUE.
  * Modified Sstar and kmstar to be consistent with rpsftm.
  * Added data\_aft and fit\_aft to output.
  * Added data\_outcome and fit\_outcome to output.
  * Added raw stratum and treat information to data\_outcome.

* Updates to tsesimp:

  * Changed the default value of admin\_recensor\_only from FALSE to TRUE.
  * Added data\_aft, fit\_aft, data\_outcome, and fit\_outcome to output.
  * Added psi1hats in bootstrap.
  * Added raw stratum information to data\_outcome.
  * Added agerand, sex.f, tt\_Lnum, rmh\_alea.c, and pathway.f as baseline covariates in the example.

* Updates to tsegest:

  * Changed the default value of admin\_recensor\_only from FALSE to TRUE.
  * Added the parameter n\_eval\_z to evaluate the Wald statistics at a sequence of psi values.
  * Added data\_switch, km\_switch, eval\_z, data\_nullcox, fit\_nullcox, data\_logis, fit\_logis, data\_outcome, and fit\_outcome to output.
  * Included data\_switch, km\_switch, eval\_z, data\_nullcox, fit\_nullcox, data\_logis, and fit\_logis in the analysis\_switch list.
  * Added raw stratum information to data\_outcome.

* Updates to ipcw:

  * Added relative\_time parameter.
  * Added data\_switch to output.
  * Added ns\_df to output.
  * Renamed df\_outcome to data\_outcome.
  * Added raw stratum information to data\_outcome.

* Added vignettes:

  * rpsftm.Rmd
  * ipe.Rmd
  * tsesimp.Rmd
  * tsegest.Rmd
  * ipcw.Rmd



# trtswitch 0.1.0

* Initial release.
