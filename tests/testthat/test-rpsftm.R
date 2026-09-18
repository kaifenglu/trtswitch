library(dplyr, warn.conflicts = FALSE)
library(survival)

testthat::test_that("rpsftm: control to active switch", {
  data1 <- immdef %>% mutate(rx = 1-xoyrs/progyrs)
  
  fit1 <- rpsftm(
    data1, id = "id", time = "progyrs", event = "prog", treat = "imm", 
    rx = "rx", censor_time = "censyrs", gridsearch = FALSE, boot = FALSE)
  
  # log-rank for ITT
  fit_lr <- survdiff(Surv(progyrs, prog) ~ imm, data = data1)
  z_lr <- sqrt(fit_lr$chisq)
  if (fit_lr$obs[2] < fit_lr$exp[2]) z_lr <- -z_lr
  
  f <- function(psi) {
    data2 <- data1 %>%
      mutate(u_star = xoyrs + (progyrs - xoyrs)*exp(psi),
             c_star = ifelse(imm == 0, pmin(censyrs, censyrs*exp(psi)), 1e10),
             t_star = pmin(u_star, c_star),
             d_star = ifelse(c_star < u_star, 0, prog))
    fit_lr <- survdiff(Surv(t_star, d_star) ~ imm, data = data2)
    z_lr <- sqrt(fit_lr$chisq)
    if (fit_lr$obs[2] < fit_lr$exp[2]) z_lr <- -z_lr
    z_lr
  }
  
  # psi based on log-rank test
  psi <- uniroot(f, c(-2,2), tol = 1e-6)$root
  
  data2 <- data1 %>%
    filter(imm == 0) %>%
    mutate(u_star = xoyrs + (progyrs - xoyrs)*exp(psi),
           c_star = pmin(censyrs, censyrs*exp(psi)),
           t_star = pmin(u_star, c_star),
           d_star = ifelse(c_star < u_star, 0, prog)) %>%
    select(-c("u_star", "c_star")) %>%
    bind_rows(data1 %>%
                filter(imm == 1) %>%
                mutate(t_star = progyrs, 
                       d_star = prog))
  
  fit <- coxph(Surv(t_star, d_star) ~ imm, data = data2)
  beta <- as.numeric(fit$coefficients[1])
  se <- beta/z_lr
  
  zcrit <- qnorm(0.975)
  hr1 <- exp(c(beta, beta - zcrit*se, beta + zcrit*se))
  testthat::expect_equal(hr1, c(fit1$hr, fit1$hr_CI))

  s <- summary(fit1)
  testthat::expect_s3_class(s, "summary.rpsftm")
  testthat::expect_true(is.call(fit1$call))
  testthat::expect_identical(s$call, fit1$call)
  testthat::expect_equal(s$n_input, nrow(data1))
  testthat::expect_equal(s$n_analyzed, sum(fit1$event_summary$n))
  testthat::expect_equal(
    s$estimates$estimate,
    c(fit1$psi, exp(-fit1$psi), fit1$hr)
  )
  testthat::expect_equal(nrow(s$population), 2L)
  testthat::expect_true(all(c(
    "model", "estimation", "diagnostics", "reporting_notes"
  ) %in% names(s)))
  testthat::expect_identical(s$pvalue_type, "ITT log-rank")
  testthat::expect_null(s$bootstrap)

  printed <- capture.output(print(s))
  testthat::expect_true(any(grepl("Model specification", printed)))
  testthat::expect_true(any(grepl("Diagnostics", printed)))

  summary_section <- function(start, end) {
    lines <- printed[(match(start, printed) + 1L):(match(end, printed) - 1L)]
    lines[nzchar(lines)]
  }

  model_lines <- summary_section("Model specification", "Estimation procedure")
  testthat::expect_true(all(nchar(model_lines) <= 80L))
  testthat::expect_false(any(grepl("^\\s", model_lines)))
  testthat::expect_true(any(grepl(
    "Structural model.*: RPSFTM; exposure measured by rx",
    model_lines
  )))
  testthat::expect_true(any(grepl(
    "Outcome model.*: Cox PH for counterfactual unswitched survival",
    model_lines
  )))

  estimation_lines <- summary_section(
    "Estimation procedure", "Treatment effect estimates (95% confidence intervals)"
  )
  diagnostic_lines <- summary_section("Diagnostics", "Reporting notes")
  for (lines in list(estimation_lines, diagnostic_lines)) {
    testthat::expect_true(all(nchar(lines) <= 80L))
    testthat::expect_false(any(grepl("^\\s", lines)))
    testthat::expect_true(all(grepl(": ", lines, fixed = TRUE)))
  }
  testthat::expect_true(any(grepl("Search method.*: ", estimation_lines)))
  testthat::expect_true(any(grepl("Model failure.*: ", diagnostic_lines)))

  plots <- plot(fit1, show_hr = FALSE, show_risk = FALSE)
  testthat::expect_named(plots, c("p_z", "p_kmstar", "p_km"))
  testthat::expect_s3_class(plots$p_kmstar, "ggplot")
  testthat::expect_identical(
    plots$p_kmstar$labels$title,
    "Kaplan-Meier Curves for Counterfactual Untreated Outcomes"
  )
  testthat::expect_equal(nrow(plots$p_kmstar$data), nrow(fit1$kmstar))
  testthat::expect_equal(
    plots$p_kmstar$data$month,
    plots$p_kmstar$data$time / 30.4375
  )
  testthat::expect_identical(plots$p_kmstar$labels$x, "Months")
  built_kmstar <- ggplot2::ggplot_build(plots$p_kmstar)
  testthat::expect_length(unique(built_kmstar$data[[1]]$group), 2L)
  testthat::expect_length(unique(built_kmstar$data[[1]]$colour), 2L)
  testthat::expect_identical(unique(built_kmstar$data[[1]]$linetype), 1L)

  min_surv_star <- data.table::data.table(plots$p_kmstar$data)[
    , min(get("surv")), by = "randomized_arm"][, get("V1")]
  expected_legend_position <- if (max(min_surv_star) < 0.5) {
    c(0.7, 0.85)
  } else {
    c(0.15, 0.25)
  }
  testthat::expect_equal(
    plots$p_kmstar$theme$legend.position.inside,
    expected_legend_position
  )
})
