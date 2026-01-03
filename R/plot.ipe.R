#' @title Plot method for ipe objects
#' @description Generate box plot for AFT model deviance residuals and 
#' Kaplan-Meier (KM) plot for potential outcomes of an ipe object.
#'
#' @param x An object of class \code{ipe}.
#' @param time_unit The time unit used in the input data.
#'   Options are "day" (default), "week", "month", or "year".
#' @param show_hr Logical; whether to show hazard ratio on the KM plot.
#'   Default is TRUE.
#' @param show_risk Logical; whether to show number at risk table
#'   below the KM plot. Default is TRUE.
#' @param ... Ensures that all arguments starting from "..." are named.
#' 
#' @return A list of two ggplot2 objects, one for box plot and the other 
#' for KM plot.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @method plot ipe
#' @export
plot.ipe <- function(x, time_unit = "day", 
                     show_hr = TRUE, show_risk = TRUE, ...) {
  if (!inherits(x, "ipe")) {
    stop("x must be of class 'ipe'")
  }
  
  alpha <- x$settings$alpha
  conflev <- 100*(1-alpha)
  treat_var <- x$settings$treat
  
  if (!is.null(x$data_outcome) && nrow(x$data_outcome) > 0) {
    # --- Deviance residuals plot for AFT models ---
    arm <- x$settings$data[[treat_var]]
    
    df1 <- data.frame(arm = x$data_aft[[treat_var]], 
                      res = x$res_aft)
    
    if (is.factor(arm)) {
      df1$arm <- factor(df1$arm, labels = levels(arm))
    } else if (is.numeric(arm) && all(arm %in% c(0, 1))) {
      df1$arm <- factor(df1$arm, levels = c(1, 0),
                        labels = c("Treatment", "Control"))
    } else {
      df1$arm <- factor(df1$arm)
    }
    
    
    p_res <- ggplot2::ggplot(df1, ggplot2::aes(x = .data$arm, 
                                               y = .data$res)) +
      ggplot2::geom_boxplot(fill="#77bd89", color="#1f6e34", alpha = 0.6) +
      ggplot2::scale_x_discrete(drop = FALSE) + 
      ggplot2::labs(x = NULL, y = "Deviance Residuals") +
      ggplot2::theme_bw()
    
    
    # --- Kaplan-Meier plot for counterfactual outcomes ---
    df <- x$km_outcome
    if (time_unit == "day") {
      df$month <- df$time / 30.4375
    } else if (time_unit == "week") {
      df$month <- df$time / 4.3482
    } else if (time_unit == "month") {
      df$month <- df$time
    } else if (time_unit == "year") {
      df$month <- df$time * 12
    } else {
      stop("time_unit must be one of 'day', 'week', 'month', or 'year'")
    }
    
    if (!is.factor(df[[treat_var]])) {
      if (is.numeric(df[[treat_var]]) && all(df[[treat_var]] %in% c(0, 1))) {
        df[[treat_var]] <- factor(df[[treat_var]], levels = c(1, 0),
                                  labels = c("Treatment", "Control"))
      } else {
        df[[treat_var]] <- factor(df[[treat_var]])
      }
    }
    
    min_surv <- data.table::data.table(df)[
      , min(get("surv")), by = treat_var][, get("V1")]
    
    p_km <- ggplot2::ggplot(df, ggplot2::aes(x = .data$month, 
                                             y = .data$surv, 
                                             group = .data[[treat_var]],
                                             colour = .data[[treat_var]])) +
      ggplot2::geom_step() +
      ggplot2::scale_x_continuous(n.breaks = 11) +
      ggplot2::scale_y_continuous(limits = c(0, 1)) +
      ggplot2::labs(
        x = "Months", 
        y = "Survival Probability",
        title = "Kaplan-Meier Curves for Counterfactual Outcomes") + 
      ggplot2::theme_bw() + 
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5),
        legend.title = ggplot2::element_blank(),
        panel.grid.minor.x = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(t = 2, r = 5, b = 0, l = 20))
    
    if (max(min_surv) < 0.5) {
      p_km <- p_km +
        ggplot2::theme(legend.position = c(0.7, 0.85))
    } else{
      p_km <- p_km +
        ggplot2::theme(legend.position = c(0.15, 0.25))
    }
    
    # add hazard ratio to plot
    if (show_hr) {
      if (max(min_surv) < 0.5) {
        p_km <- p_km  + 
          ggplot2::annotate(
            "text", x = 0.6*max(df$month), y = 0.7, hjust = 0,
            label = sprintf("HR = %.3f (%.0f%% CI: %.3f, %.3f)", 
                            x$hr, conflev, 
                            x$hr_CI[1], x$hr_CI[2]),
            size = 3.5, color = "black"
          )
      } else {
        p_km <- p_km  + 
          ggplot2::annotate(
            "text", x = 0, y = 0, hjust = 0,
            label = sprintf("HR = %.3f (%.0f%% CI: %.3f, %.3f)", 
                            x$hr, conflev, 
                            x$hr_CI[1], x$hr_CI[2]),
            size = 3.5, color = "black"
          )
      }
    }
    
    # add number at risk table below plot
    if (show_risk) {
      xbreaks <- ggplot2::ggplot_build(p_km)$layout$panel_params[[1]]$x$breaks
      xbreaks <- xbreaks[!is.na(xbreaks)]
      limits <- c(min(xbreaks), max(max(xbreaks), max(df$month)))
      
      # --- Compute number at risk at or before each tick ---
      tablist <- lapply(0:1, function(h) {
        t <- df$month[df$treated == h]
        n <- df$nrisk[df$treated == h]
        
        idx <- findInterval(xbreaks, t) + 1
        atrisk <- ifelse(idx <= length(n), n[idx], 0)
        
        df1 <- data.frame(time = xbreaks, atrisk = atrisk)
        df1[[treat_var]] <- df[[treat_var]][df$treated == h][1]
        df1
      })
      
      df_risk <- do.call(rbind, tablist)
      df_risk[[treat_var]] <- factor(df_risk[[treat_var]],
                                     levels = levels(df[[treat_var]]))
      
      # --- Create number at risk plot ---
      p_risk <- ggplot2::ggplot(df_risk, 
                                ggplot2::aes(x = .data$time, 
                                             y = .data[[treat_var]], 
                                             label = .data$atrisk, 
                                             colour = .data[[treat_var]])) +
        ggplot2::geom_text(size = 3.2, na.rm = TRUE) +
        ggplot2::scale_x_continuous(
          breaks = xbreaks, limits = range(xbreaks)) +
        ggplot2::scale_y_discrete(
          limits = rev(levels(df_risk[[treat_var]]))) +
        ggplot2::coord_cartesian(clip = "off") + 
        ggplot2::theme_minimal() +
        ggplot2::theme(
          axis.title = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_blank(),
          axis.ticks = ggplot2::element_blank(),
          legend.position = "none",
          panel.grid = ggplot2::element_blank(),
          plot.margin = ggplot2::margin(t = 6, r = 5, b = 0, l = 20)) + 
        ggplot2::annotate(
          "text", x = min(xbreaks), y = 3,
          label = "No. of Subjects at Risk",
          size = 4, hjust = 0.5)
      
      suppressMessages({ 
        p_km <- p_km +
          ggplot2::scale_x_continuous(
            breaks = xbreaks, limits = limits, expand = c(0.05, 0))
        
        p_risk <- p_risk +
          ggplot2::scale_x_continuous(
            breaks = xbreaks, limits = limits, expand = c(0.05, 0))
      })
      
      # 2. Align plots using align_plots()
      aligned <- cowplot::align_plots(p_km, p_risk, align = "v", axis = "lr")
      
      # 3. Combine with plot_grid()
      p_km <- cowplot::plot_grid(aligned[[1]], aligned[[2]], ncol = 1, 
                                 rel_heights = c(4, 0.6))    
    }
    
    list(p_res = p_res, p_km = p_km)
  } else {
    stop("No outcome data available to plot.")
  }
}
