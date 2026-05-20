# Inverse Probability of Censoring Weights

``` r
library(trtswitch)
library(dplyr, warn.conflicts = FALSE)
library(ggplot2)
```

## Introduction

The inverse probability of censoring weights (IPCW) method is a powerful
tool for adjusting survival analysis in the presence of treatment
switching. In scenarios where patients switch treatments, survival data
are artificially censored, leading to biased estimates if not properly
addressed. This bias arises because the decision to switch treatments
may depend on prognostic factors that also influence survival outcomes.

To mitigate this bias, researchers gather data on prognostic covariates
at baseline and at regular intervals until one of the following occurs:
a treatment switch, death, drop-out, or the end of the follow-up period.
A “switching model” is then developed to estimate weights for each
patient at each time point. These weights represent the inverse
probability of remaining unswitched over time. Each patient retains an
unstabilized weight greater than or equal to 1, which accounts for both
their own data and that of switchers with similar prognostic profiles.

Finally, an “outcome” model utilizing IPCW-weighted survival times is
employed to estimate the treatment effect, adjusted for switching.

## Estimation of weights

The switching model can be implemented using either pooled logistic
regression or proportional hazards regression with time-dependent
covariates.

### Pooled Logistic Regression Model

Let \\A\_{i,j}\\ denote the indicator of alternative therapy for subject
\\i\\ in treatment cycle \\j\\. We assume that once a patient switched
to an alternative therapy, he/she will stay on the alternative therapy.
Let \\X_i\\ denote the baseline covariates for subject \\i\\ and
\\L\_{i,j}\\ denote the time-dependent covariates for subject \\i\\ at
the start of cycle \\j\\. For the pooled logistic regression model for
treatment switching, the full data will be used for patients who didn’t
switch treatment, while the partial data up to and including the
treatment cycle when the patient switched treatment will be used for
patients who switched treatment. The logistic regression model is used
to estimate \\P(A_j = 0 \| \bar{A}\_{j-1} = \bar{0}, \bar{L}\_{j} =
\bar{l}\_{i,j}, X = x_i)\\ where \\x_i\\, \\l\_{i,j}\\, and \\a\_{i,j}\\
denote the observed value of baseline covariates, time-dependent
covariates at the start of cycle \\j\\, and alternative therapy status
for cycle \\j\\ of subject \\i\\, respectively. A subject is on
alternative therapy in cycle \\j\\ if and only if treatment switching
has occurred at the start of that cycle, which is equivalent to
switching having occurred by the end of cycle \\j−1\\. Consequently,
lagged values of time-dependent covariates are often used in the
logistic regression model. Natural cubic splines can be used to model
visit-specific intercepts for the pooled logistic regression model. The
unstabilized weights for subject \\i\\ for cycle \\t\\ are given by
\\w_i(t) = \prod\_{j=1}^{t} \frac{1}{P(A_j = 0 \| \bar{A}\_{j-1} =
\bar{0}, \bar{L}\_{j} = \bar{l}\_{i,j}, X = x_i)}\\ In the IPCW method,
all visits that occur after treatment switching are excluded from the
weighted outcome model. Accordingly, the weights are defined as the
inverse probability of not having switched at the start of each cycle.

When using pooled logistic regression for the switching model, the
following specifications can be made:

- **Visit-Specific Intercepts**: These can be modeled using a natural
  cubic spline with specified degrees of freedom, denoted as \\\nu\\
  (corresponding to the `ns_df` parameter of the `ipcw` function). The
  boundary and internal knots can be based on the range and percentiles
  of observed treatment switching times, respectively. Here, \\\nu\\
  equals the number of internal knots plus one; a value of \\\nu=0\\
  indicates a common intercept, while \\\nu=1\\ leads to a linear effect
  of visit. By default, \\\nu=3\\, which implies two internal knots for
  the cubic spline.

- **Stratification Factors**: If more than one stratification factor
  exists, they can be included in the logistic model either as main
  effects only (`strata_main_effect_only = TRUE`) or as all possible
  combinations of the levels of the stratification factors
  (`strata_main_effect_only = FALSE`).

- **Handling Nonconvergence**: If the pooled logistic regression model
  encounters issues such as complete or quasi-complete separation,
  Firth’s bias-reducing penalized likelihood method (`firth = TRUE`) can
  be applied alongside intercept correction (`flic = TRUE`) to yield
  more accurate predicted probabilities.

### Proportional Hazards Model

When using a Cox model with time-dependent covariates as the switching
model, adjacent observations within a subject may be combined as long as
the values of time-dependent covariates do not change over each combined
interval. This approach saves data storage but observations across
subjects will not be aligned at regular intervals. In this case, unique
event times across treatment arms should be replicated within each
subject, enabling the availability of weights at each event time.

In the time-dependent Cox switching model, survival times are
artificially censored at the point of treatment switching. The IPCW
weight is then defined as the probability of not having switched
treatment by the end of the time interval. This is to account for the
irregular visit structure for the time-dependent Cox switching model.

### Stabilized Weights

To address the high variability in “unstabilized” weights, researchers
often opt for “stabilized” weights. The switching model for stabilized
weights involves a model for the numerator of the weight in addition to
the model for the denominator. Typically, the numerator model mirrors
the denominator model but includes only prognostic baseline variables.
Any baseline variables included in the numerator must be part of the
weighted outcome model. Any baseline variables included in the weighted
outcome model must also be part of the denominator.

### Truncation of Weights

You can specify weight truncation using the `trunc` and
`trunc_upper_only` parameters of the `ipcw` function. The `trunc`
parameter is a number between 0 and 0.5 that represents the fraction of
the weight distribution to truncate. The `trunc_upper_only` parameter is
a flag that indicates whether to truncate weights only from the upper
end of the weight distribution.

## Estimation of Hazard Ratio

Once weights have been obtained for each event time in the data set, we
fit a (potentially stratified) Cox proportional hazards model to this
weighted data set to obtain an estimate of the hazard ratio. The
confidence interval for the hazard ratio can be obtained from the
weighted Cox model with robust sandwich variance estimates or derived by
bootstrapping the weight estimation and the subsequent model-fitting
process.

### Example for Pooled Logistic Regression Switching Model

First we prepare the data.

``` r
sim1 <- tssim(
  tdxo = 1, coxo = 1, allocation1 = 1, allocation2 = 1,
  p_X_1 = 0.3, p_X_0 = 0.3, 
  rate_T = 0.002, beta1 = -0.5, beta2 = 0.3, 
  gamma0 = 0.3, gamma1 = -0.9, gamma2 = 0.7, gamma3 = 1.1, gamma4 = -0.8,
  zeta0 = -3.5, zeta1 = 0.5, zeta2 = 0.2, zeta3 = -0.4, 
  alpha0 = 0.5, alpha1 = 0.5, alpha2 = 0.4, 
  theta1_1 = -0.4, theta1_0 = -0.4, theta2 = 0.2,
  rate_C = 0.0000855, accrualIntensity = 20/30, 
  fixedFollowup = 0, plannedTime = 1350, days = 30,
  n = 500, NSim = 100, seed = 314159)
```

Given that the observations are collected at regular intervals, we can
use a pooled logistic regression switching model. In this model, the
prognostic baseline variable is represented by `bprog`, while `L` serves
as a time-dependent covariate. Additionally, `ns1`, `ns2`, and `ns3` are
the three terms for the natural cubic spline with \\\nu=3\\ degrees of
freedom for visit-specific intercepts.

``` r
fit1 <- ipcw(
  sim1[[1]], id = "id", tstart = "tstart", 
  tstop = "tstop", event = "event", treat = "trtrand", 
  swtrt = "xo", swtrt_time = "xotime", 
  base_cov = "bprog", numerator = "bprog", 
  denominator = c("bprog", "L"),
  logistic_switching_model = TRUE, ns_df = 3,
  swtrt_control_only = TRUE, boot = FALSE)
```

The fits for the denominator and numerator switching models for the
control arm are as follows.

``` r
# denominator switching model fit
fit1$fit_switch[[1]]$fit_den$parest[, c("param", "beta", "sebeta", "z")]
#>         param       beta    sebeta           z
#> 1 (Intercept) -4.3938500 0.4364956 -10.0661945
#> 2       bprog  0.3602200 0.2453795   1.4680116
#> 3           L  0.3064046 0.3067097   0.9990053
#> 4         ns1 -0.8325131 0.6223178  -1.3377621
#> 5         ns2  2.1816608 0.8374703   2.6050606
#> 6         ns3  0.9189778 0.7471678   1.2299482

# numerator switching model fit
fit1$fit_switch[[1]]$fit_num$parest[, c("param", "beta", "sebeta", "z")]
#>         param       beta    sebeta          z
#> 1 (Intercept) -4.1873341 0.3784297 -11.065025
#> 2       bprog  0.4102063 0.2409982   1.702114
#> 3         ns1 -0.8233014 0.6204201  -1.327006
#> 4         ns2  2.2205113 0.8361129   2.655755
#> 5         ns3  0.9264175 0.7453640   1.242906
```

Since treatment switching is not allowed in the experimental arm, the
weights are identical to 1 for subjects in the experimental group. The
unstabilized and stablized weights for the control group are plotted
below.

``` r
# unstabilized weights
ggplot(fit1$data_outcome %>% filter(trtrand == 0), 
       aes(x = unstabilized_weight)) + 
  geom_histogram(fill="#77bd89", color="#1f6e34", alpha=0.8) +
  scale_x_continuous("unstabilized weights")
#> `stat_bin()` using `bins = 30`. Pick better value `binwidth`.
```

![](ipcw_files/figure-html/weights%20example%201-1.png)

``` r

# stabilized weights
ggplot(fit1$data_outcome %>% filter(trtrand == 0), 
       aes(x = stabilized_weight)) + 
  geom_histogram(fill="#77bd89", color="#1f6e34", alpha=0.8) +
  scale_x_continuous("stabilized weights")
#> `stat_bin()` using `bins = 30`. Pick better value `binwidth`.
```

![](ipcw_files/figure-html/weights%20example%201-2.png)

Now we fit a weighted outcome Cox model and compare the treatment hazard
ratio estimate with the reported.

``` r
fit1$fit_outcome$parest[, c("param", "beta", "sebeta", "z")]
#>     param       beta    sebeta         z
#> 1 treated -0.5919728 0.1118087 -5.294515
#> 2   bprog  0.1652406 0.1196620  1.380894
    
exp(fit1$fit_outcome$parest[1, c("beta", "lower", "upper")])
#>        beta     lower     upper
#> 1 0.5532348 0.4443629 0.6887811
```

### Example for Time-Dependent Cox Switching Model

Now we apply the IPCW method using a Cox proportional hazards model with
time-dependent covariates as the switching model.

``` r
fit2 <- ipcw(
  shilong, id = "id", tstart = "tstart", tstop = "tstop", 
  event = "event", treat = "bras.f", swtrt = "co", 
  swtrt_time = "dco", 
  base_cov = c("agerand", "sex.f", "tt_Lnum", "rmh_alea.c", 
               "pathway.f"),
  numerator = c("agerand", "sex.f", "tt_Lnum", "rmh_alea.c", 
                "pathway.f"),
  denominator = c("agerand", "sex.f", "tt_Lnum", "rmh_alea.c",
                  "pathway.f", "ps", "ttc", "tran"),
  swtrt_control_only = FALSE, boot = FALSE)
```

The fits for the denominator and numerator switching models for the
control arm are as follows.

``` r
# denominator switching model for the control group
fit2$fit_switch[[1]]$fit_den$parest[, c("param", "beta", "sebeta", "z")]
#>                    param         beta     sebeta           z
#> 1                agerand  0.007642073 0.01023088  0.74696131
#> 2            sex.fFemale -0.364211003 0.28631111 -1.27208128
#> 3                tt_Lnum  0.042406215 0.05314773  0.79789331
#> 4             rmh_alea.c -0.351616244 0.27416745 -1.28248719
#> 5            pathway.fHR -0.041768569 0.43444787 -0.09614173
#> 6 pathway.fPI3K.AKT.mTOR  0.267055320 0.43736340  0.61060281
#> 7                     ps  0.103478553 0.18544824  0.55799156
#> 8                    ttc -0.480378314 0.33557237 -1.43151926
#> 9                   tran  0.295828936 0.45806296  0.64582592

# numerator switching model for the control group
fit2$fit_switch[[1]]$fit_num$parest[, c("param", "beta", "sebeta", "z")]
#>                    param        beta      sebeta           z
#> 1                agerand  0.01059477 0.009501694  1.11503993
#> 2            sex.fFemale -0.32452112 0.283829120 -1.14336794
#> 3                tt_Lnum  0.05956492 0.051120376  1.16518947
#> 4             rmh_alea.c -0.29758402 0.266335628 -1.11732711
#> 5            pathway.fHR  0.02467202 0.429102680  0.05749678
#> 6 pathway.fPI3K.AKT.mTOR  0.28245603 0.434306994  0.65036030
```

The fits for the denominator and numerator switching models for the
experimental arm are as follows:

``` r
# denominator switching model for the experimental group
fit2$fit_switch[[2]]$fit_den$parest[, c("param", "beta", "sebeta", "z")]
#>                    param         beta     sebeta          z
#> 1                agerand -0.002639441 0.01796294 -0.1469381
#> 2            sex.fFemale  0.428028469 0.47154276  0.9077193
#> 3                tt_Lnum -0.152758042 0.10774494 -1.4177746
#> 4             rmh_alea.c  0.210365664 0.50502449  0.4165455
#> 5            pathway.fHR  1.774466817 1.04918789  1.6912765
#> 6 pathway.fPI3K.AKT.mTOR  0.859950234 1.08995673  0.7889765
#> 7                     ps  0.261142698 0.27003546  0.9670682
#> 8                    ttc -0.408175689 0.53635218 -0.7610218
#> 9                   tran  1.053347294 0.69246856  1.5211482

# numerator switching model for the experimental group
fit2$fit_switch[[2]]$fit_num$parest[, c("param", "beta", "sebeta", "z")]
#>                    param         beta     sebeta           z
#> 1                agerand  0.001625622 0.01812956  0.08966694
#> 2            sex.fFemale  0.476979559 0.45994387  1.03703863
#> 3                tt_Lnum -0.160020273 0.10615092 -1.50747886
#> 4             rmh_alea.c  0.154833260 0.48178554  0.32137382
#> 5            pathway.fHR  1.841535014 1.03960949  1.77137188
#> 6 pathway.fPI3K.AKT.mTOR  1.064637540 1.06675981  0.99801055
```

Below, the unstabilized and stabilized weights are plotted by treatment
group: the left panel displays data for the control group, while the
right panel shows data for the experimental group.

``` r
# unstabilized weights
ggplot(fit2$data_outcome, aes(x = unstabilized_weight)) + 
  geom_histogram(fill="#77bd89", color="#1f6e34", alpha=0.8) + 
  scale_x_continuous("unstabilized weights") + 
  facet_wrap(~treated) 
#> `stat_bin()` using `bins = 30`. Pick better value `binwidth`.
```

![](ipcw_files/figure-html/weights%20example%202-1.png)

``` r

# stabilized weights
ggplot(fit2$data_outcome, aes(x = stabilized_weight)) + 
  geom_histogram(fill="#77bd89", color="#1f6e34", alpha=0.8) + 
  scale_x_continuous("stabilized weights") + 
  facet_wrap(~treated)
#> `stat_bin()` using `bins = 30`. Pick better value `binwidth`.
```

![](ipcw_files/figure-html/weights%20example%202-2.png)

Finally, we fit the weighted outcome Cox model and compare the treatment
hazard ratio estimate with the reported.

``` r
fit2$fit_outcome$parest[, c("param", "beta", "sebeta", "z")]
#>                    param         beta     sebeta          z
#> 1                treated  0.356390611 0.25526832  1.3961412
#> 2                agerand -0.006047034 0.01018596 -0.5936636
#> 3            sex.fFemale -0.487409540 0.24458876 -1.9927716
#> 4                tt_Lnum  0.011244574 0.04147245  0.2711336
#> 5             rmh_alea.c  0.941651485 0.25315061  3.7197283
#> 6            pathway.fHR -0.127273307 0.35878557 -0.3547336
#> 7 pathway.fPI3K.AKT.mTOR -0.166035907 0.34997206 -0.4744262

c(fit2$hr, fit2$hr_CI)
#> [1] 1.4281653 0.8659517 2.3553924
```
