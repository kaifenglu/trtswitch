# trtswitch

`trtswitch` provides methods for treatment switching adjustment in
randomized clinical trials, including:

- rank preserving structural failure time model (RPSFTM)
- iterative parameter estimation (IPE)
- inverse probability of censoring weights (IPCW)
- marginal structural model (MSM)
- simple two-stage estimation (TSEsimp)
- improved two-stage estimation with g-estimation (TSEgest)

## Installation

Install from GitHub:

``` r
# install.packages("remotes")
remotes::install_github("kaifenglu/trtswitch")
```

## Example: RPSFTM

``` r
library(trtswitch)
library(dplyr)

# Build treatment exposure proportion used by rpsftm
# in the one-way switching example dataset.
data <- immdef %>%
  mutate(rx = 1 - xoyrs / progyrs)

fit <- rpsftm(
  data = data,
  id = "id",
  time = "progyrs",
  event = "prog",
  treat = "imm",
  rx = "rx",
  censor_time = "censyrs",
  boot = FALSE
)

fit

# Key estimates
fit$psi
fit$hr
fit$psi_CI
fit$hr_CI
```

## Documentation

- Website: <https://kaifenglu.github.io/trtswitch/>
- Function reference: <https://kaifenglu.github.io/trtswitch/reference/>
- Issues: <https://github.com/kaifenglu/trtswitch/issues>
