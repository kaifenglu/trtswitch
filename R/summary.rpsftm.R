#' @title Summary method for rpsftm objects
#' @description Summarizes the specification, treatment effect estimates, and
#' diagnostics from a rank-preserving structural failure time model.
#'
#' @param object An object of class \code{rpsftm}.
#' @param ... Additional arguments passed to or from other methods.
#'
#' @return An object of class \code{summary.rpsftm} containing the analysis
#' population, model specification, estimation procedure, treatment effect
#' estimates, diagnostics, bootstrap summary, and reporting notes.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
summary.rpsftm <- function(object, ...) {
  if (!inherits(object, "rpsftm")) {
    stop("object must be of class 'rpsftm'")
  }

  settings <- object$settings
  event_summary <- object$event_summary
  treat_var <- settings$treat

  population <- event_summary
  arm <- if ("treated" %in% names(population)) population$treated else 0:1
  if (treat_var %in% names(population)) {
    arm <- as.character(population[[treat_var]])
  } else {
    arm <- ifelse(arm == 0, "Control", "Treatment")
  }
  population <- cbind(arm = arm, population[, setdiff(
    names(population), c("treated", treat_var)), drop = FALSE])

  n_analyzed <- if ("n" %in% names(event_summary)) {
    sum(event_summary$n, na.rm = TRUE)
  } else {
    NA_integer_
  }
  n_input <- if (!is.null(settings$data)) nrow(settings$data) else NA_integer_
  n_excluded <- if (is.na(n_input) || is.na(n_analyzed)) {
    NA_integer_
  } else {
    n_input - n_analyzed
  }

  g_test <- switch(
    settings$psi_test,
    logrank = "Log-rank test",
    phreg = "Cox PH model",
    lifereg = paste0("Parametric AFT (", settings$aft_dist, ")"),
    settings$psi_test
  )
  search_method <- if (isTRUE(settings$gridsearch)) {
    "Grid search with linear interpolation"
  } else {
    paste("Root finding (", settings$root_finding, ")", sep = "")
  }

  model <- data.frame(
    item = c(
      "Structural model", "G-estimation", "Strata", "Baseline covariates",
      "Treatment effect modifier", "Re-censoring",
      "Administrative censoring only", "Skip non-switching arms",
      "Outcome model", "Ties method"
    ),
    value = c(
      paste0("RPSFTM; exposure measured by ", settings$rx),
      g_test,
      if (length(settings$stratum) == 0 || all(settings$stratum == "")) {
        "None"
      } else {
        paste(settings$stratum, collapse = ", ")
      },
      if (length(settings$base_cov) == 0 || all(settings$base_cov == "")) {
        "None"
      } else {
        paste(settings$base_cov, collapse = ", ")
      },
      as.character(settings$treat_modifier),
      if (isTRUE(settings$recensor)) "Yes" else "No",
      if (isTRUE(settings$admin_recensor_only)) "Yes" else "No",
      if (isTRUE(settings$autoswitch)) "Yes" else "No",
      "Cox PH for counterfactual unswitched survival",
      settings$ties
    ),
    stringsAsFactors = FALSE
  )

  estimation <- data.frame(
    item = c(
      "Search method", "Search interval", "Number of Z evaluations",
      "Root-finding tolerance", "Number of zero-crossings", "Zero-crossings",
      "Selected root"
    ),
    value = c(
      search_method,
      paste0("[", settings$low_psi, ", ", settings$hi_psi, "]"),
      as.character(settings$n_eval_z),
      as.character(settings$tol),
      as.character(sum(is.finite(object$psi_roots))),
      if (any(is.finite(object$psi_roots))) {
        paste(format(object$psi_roots[is.finite(object$psi_roots)], digits = 6),
              collapse = ", ")
      } else {
        "None"
      },
      format(object$psi, digits = 6)
    ),
    stringsAsFactors = FALSE
  )

  conf_level <- 100 * (1 - settings$alpha)
  estimates <- data.frame(
    estimand = c(
      "Causal parameter psi", "Causal survival time ratio",
      "Adjusted hazard ratio"
    ),
    estimate = c(object$psi, exp(-object$psi), object$hr),
    lower = c(object$psi_CI[1], exp(-object$psi_CI[2]), object$hr_CI[1]),
    upper = c(object$psi_CI[2], exp(-object$psi_CI[1]), object$hr_CI[2]),
    ci_method = c(object$psi_CI_type, object$psi_CI_type, object$hr_CI_type),
    stringsAsFactors = FALSE
  )

  z_values <- object$eval_z$Z
  z_differences <- diff(z_values[is.finite(z_values)])
  z_direction <- if (length(z_differences) == 0) {
    "Not assessable"
  } else if (all(z_differences <= 0)) {
    "Monotonically decreasing"
  } else if (all(z_differences >= 0)) {
    "Monotonically increasing"
  } else {
    "Non-monotonic"
  }
  roots <- object$psi_roots[is.finite(object$psi_roots)]
  limits <- c(settings$low_psi, settings$hi_psi)
  interval_step <- diff(limits) / max(settings$n_eval_z - 1, 1)
  reported_limits <- c(object$psi, object$psi_CI)
  near_boundary <- any(
    is.finite(reported_limits) &
      (reported_limits <= limits[1] + interval_step |
         reported_limits >= limits[2] - interval_step)
  )

  diagnostics <- data.frame(
    diagnostic = c(
      "Model failure", "Psi not estimated", "Confidence limits incomplete",
      "Multiple zero-crossings", "Estimate or CI near search boundary",
      "Z-curve direction"
    ),
    result = c(
      as.character(isTRUE(object$fail)),
      as.character(isTRUE(object$psimissing) || !is.finite(object$psi)),
      as.character(any(!is.finite(object$psi_CI))),
      as.character(length(roots) > 1),
      as.character(near_boundary),
      z_direction
    ),
    stringsAsFactors = FALSE
  )

  bootstrap <- NULL
  if (isTRUE(settings$boot)) {
    failures <- as.logical(object$fail_boots)
    n_requested <- settings$n_boot
    n_failed <- sum(failures, na.rm = TRUE)
    bootstrap <- data.frame(
      requested = n_requested,
      successful = n_requested - n_failed,
      failed = n_failed,
      failure_pct = 100 * n_failed / n_requested,
      nonfinite_psi = sum(!is.finite(object$psi_boots)),
      nonfinite_hr = sum(!is.finite(object$hr_boots))
    )
  }

  reporting_notes <- c(
    paste0("The standard RPSFTM assumes a common treatment effect; ",
           "its plausibility must be justified outside the fitted object."),
    paste0("This analysis was fitted ",
           if (isTRUE(settings$recensor)) "with" else "without",
           " re-censoring."),
    "TSD 24 recommends companion analyses with and without re-censoring.",
    paste0("TSD 24 recommends sensitivity analyses over plausible ",
           "treatment-effect modifiers and alternative adjustment methods."),
    paste0("Review plot(object) to assess the g-estimation curve and ",
           "counterfactual survival estimates visually.")
  )

  out <- list(
    call = object$call,
    n_input = n_input,
    n_analyzed = n_analyzed,
    n_excluded = n_excluded,
    population = population,
    model = model,
    estimation = estimation,
    estimates = estimates,
    conf_level = conf_level,
    pvalue = object$pvalue,
    pvalue_type = if (identical(object$pvalue_type, "log-rank")) {
      "ITT log-rank"
    } else {
      object$pvalue_type
    },
    diagnostics = diagnostics,
    bootstrap = bootstrap,
    reporting_notes = reporting_notes
  )
  class(out) <- "summary.rpsftm"
  out
}

#' @title Print method for summary.rpsftm objects
#' @description Prints a detailed summary of a rpsftm fit.
#'
#' @param x An object of class \code{summary.rpsftm}.
#' @param digits The number of significant digits to print.
#' @param ... Additional arguments passed to \code{print.data.frame}.
#'
#' @return The input object, invisibly.
#'
#' @keywords internal
#'
#' @export
print.summary.rpsftm <- function(x, digits = max(3L, getOption("digits") - 3L),
                                 ...) {
  print_key_values <- function(labels, values) {
    label_width <- max(nchar(labels))
    labels <- format(labels, width = label_width, justify = "left")
    cat(paste0(labels, ": ", values, collapse = "\n"), "\n", sep = "")
  }

  cat("Rank Preserving Structural Failure Time Model\n\n")
  
  if(!is.null(cl <- x$call)) {
    cat("Call:\n")
    dput(cl)
    cat("\n")
  }
  
  cat("Analysis population\n")
  cat("Input:", x$n_input, " Analyzed:", x$n_analyzed,
      " Excluded for incomplete data:", x$n_excluded, "\n")
  print(x$population, row.names = FALSE, digits = digits, ...)

  cat("\nModel specification\n")
  print_key_values(x$model$item, x$model$value)

  cat("\nEstimation procedure\n")
  print_key_values(x$estimation$item, x$estimation$value)

  cat("\nTreatment effect estimates (", format(x$conf_level, trim = TRUE),
      "% confidence intervals)\n", sep = "")
  estimates <- x$estimates
  estimates[c("estimate", "lower", "upper")] <- lapply(
    estimates[c("estimate", "lower", "upper")], round, digits = digits)
  print(estimates, row.names = FALSE, ...)
  cat("P-value (", x$pvalue_type, "): ", format.pval(
    x$pvalue, digits = digits, eps = 10^-digits), "\n", sep = "")

  cat("\nDiagnostics\n")
  print_key_values(x$diagnostics$diagnostic, x$diagnostics$result)

  if (!is.null(x$bootstrap)) {
    cat("\nBootstrap performance\n")
    print(x$bootstrap, row.names = FALSE, digits = digits, ...)
  }

  cat("\nReporting notes\n")
  cat(paste0("* ", x$reporting_notes, collapse = "\n"), "\n")

  invisible(x)
}