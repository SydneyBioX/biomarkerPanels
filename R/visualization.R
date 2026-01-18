#' Plot Pareto Fronts
#'
#' Visualize trade-offs between two objectives for candidate biomarker panels.
#'
#' @param results Data frame with columns `objective_x`, `objective_y`, `label`.
#' @param xlab,ylab Axis labels.
#' @return A `ggplot` object.
#' @export
plot_pareto_front <- function(results,
                              xlab = "Objective X",
                              ylab = "Objective Y") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plot_pareto_front()", call. = FALSE)
  }

  ggplot2::ggplot(results, ggplot2::aes(x = objective_x, y = objective_y, label = label)) +
    ggplot2::geom_point(color = "#005B96") +
    ggplot2::geom_text(vjust = -0.4, size = 3) +
    ggplot2::labs(x = xlab, y = ylab) +
    ggplot2::theme_minimal(base_size = 12)
}

#' Plot Feature Inclusion Frequencies
#'
#' Creates a bar chart showing how often each feature appears across
#' Pareto-optimal solutions, highlighting highly stable features.
#'
#' @param stability_result Output from [`analyze_feature_stability()`].
#' @param top_n Integer; maximum number of features to display. Defaults to `20`.
#' @param highlight_threshold Numeric in \[0, 1\]; features appearing in more than
#'   this proportion of solutions are highlighted. Defaults to `0.5`.
#' @param bar_color Color for bars below the threshold. Defaults to `"#377EB8"`.
#' @param highlight_color Color for bars at or above the threshold. Defaults to `"#E41A1C"`.
#' @return A `ggplot` object.
#' @export
#' @seealso [analyze_feature_stability()]
#' @examples
#' # stability <- analyze_feature_stability(panel_result)
#' # plot_feature_stability(stability)
#' # plot_feature_stability(stability, top_n = 10, highlight_threshold = 0.7)
plot_feature_stability <- function(stability_result,
                                   top_n = 20L,
                                   highlight_threshold = 0.5,
                                   bar_color = "#377EB8",
                                   highlight_color = "#E41A1C") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plot_feature_stability()", call. = FALSE)
  }

  if (!inherits(stability_result, "FeatureStabilityResult")) {
    stop("`stability_result` must be output from analyze_feature_stability().",
         call. = FALSE)
  }

  if (!is.numeric(top_n) || length(top_n) != 1L || is.na(top_n) || top_n < 1L) {
    stop("`top_n` must be a positive integer.", call. = FALSE)
  }
  top_n <- as.integer(top_n)

  if (!is.numeric(highlight_threshold) || length(highlight_threshold) != 1L ||
      is.na(highlight_threshold) || highlight_threshold < 0 || highlight_threshold > 1) {
    stop("`highlight_threshold` must be a number between 0 and 1.", call. = FALSE)
  }

  freq_df <- stability_result$frequencies
  if (!nrow(freq_df)) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5,
                          label = "No features to display",
                          size = 5, color = "grey50") +
        ggplot2::theme_void()
    )
  }

  plot_data <- utils::head(freq_df, n = min(top_n, nrow(freq_df)))
  plot_data$stable <- plot_data$proportion >= highlight_threshold
  plot_data$feature <- factor(
    plot_data$feature,
    levels = rev(plot_data$feature)
  )

  mean_jac <- stability_result$mean_jaccard
  subtitle_text <- sprintf(
    "%d solutions analyzed (mean Jaccard: %s)",
    stability_result$n_solutions,
    if (is.na(mean_jac)) "N/A" else sprintf("%.2f", mean_jac)
  )

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .data$feature, y = .data$proportion, fill = .data$stable)
  ) +
    ggplot2::geom_col(show.legend = FALSE) +
    ggplot2::scale_fill_manual(
      values = c("FALSE" = bar_color, "TRUE" = highlight_color)
    ) +
    ggplot2::geom_hline(
      yintercept = highlight_threshold,
      linetype = "dashed",
      color = "grey40",
      linewidth = 0.5
    ) +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      labels = function(x) paste0(round(x * 100), "%"),
      expand = ggplot2::expansion(mult = c(0, 0.02))
    ) +
    ggplot2::labs(
      x = NULL,
      y = "Inclusion Frequency",
      title = "Feature Selection Stability",
      subtitle = subtitle_text
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )
}

#' Plot Per-Cohort Metric Comparison
#'
#' Creates a grouped bar chart comparing performance metrics across cohorts.
#' When a `method` column is present, bars are dodged by method for comparing
#' multiple approaches (e.g., MOO vs LASSO).
#'
#' @param cohort_metrics A `data.frame` from [evaluate_panel_by_cohort()], or
#'   any data frame with a `cohort` column and one or more metric columns.
#' @param metrics Character vector specifying which metrics to plot. Defaults
#'   to all numeric columns except count columns.
#' @param facet_by One of `"metric"` (facet by metric, x-axis = cohort) or
#'   `"cohort"` (facet by cohort, x-axis = metric). Defaults to `"metric"`.
#' @param show_values Logical; if `TRUE`, display metric values as text labels
#'   on bars. Defaults to `TRUE`.
#' @param value_size Numeric; text size for value labels. Defaults to `3`.
#' @param y_limits Numeric vector of length 2 specifying y-axis limits, or
#'   `NULL` for automatic limits. Defaults to `c(0, 1)`.
#' @param title Plot title. Defaults to `"Per-Cohort Performance Comparison"`.
#' @param subtitle Optional plot subtitle.
#' @return A `ggplot` object.
#' @export
#' @seealso [evaluate_panel_by_cohort()]
#' @examples
#' # cohort_metrics <- evaluate_panel_by_cohort(panel, x_list, y_list)
#' # plot_cohort_comparison(cohort_metrics)
#' #
#' # Compare methods:
#' # moo_metrics$method <- "MOO"
#' # lasso_metrics$method <- "LASSO"
#' # plot_cohort_comparison(rbind(moo_metrics, lasso_metrics))
plot_cohort_comparison <- function(cohort_metrics,
                                   metrics = NULL,
                                   facet_by = c("metric", "cohort"),
                                   show_values = TRUE,
                                   value_size = 3,
                                   y_limits = c(0, 1),
                                   title = "Per-Cohort Performance Comparison",
                                   subtitle = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plot_cohort_comparison()", call. = FALSE)
  }

  facet_by <- match.arg(facet_by)

  if (!is.data.frame(cohort_metrics)) {
    stop("`cohort_metrics` must be a data.frame.", call. = FALSE)
  }

  if (!"cohort" %in% names(cohort_metrics)) {
    stop("`cohort_metrics` must contain a 'cohort' column.", call. = FALSE)
  }

  count_cols <- c("n_samples", "n_positive", "n_negative")
  all_numeric <- names(cohort_metrics)[vapply(cohort_metrics, is.numeric, logical(1))]
  available_metrics <- setdiff(all_numeric, count_cols)

  if (is.null(metrics)) {
    metrics <- available_metrics
  } else {
    missing <- setdiff(metrics, available_metrics)
    if (length(missing)) {
      stop("Metric(s) not found in data: ", paste(missing, collapse = ", "),
           call. = FALSE)
    }
  }

  if (length(metrics) == 0L) {
    stop("No metrics to plot.", call. = FALSE)
  }

  has_method <- "method" %in% names(cohort_metrics)
  id_cols <- if (has_method) c("cohort", "method") else "cohort"

  plot_data <- stats::reshape(
    cohort_metrics[, c(id_cols, metrics), drop = FALSE],
    direction = "long",
    varying = metrics,
    v.names = "value",
    timevar = "metric",
    times = metrics,
    idvar = id_cols
  )
  plot_data$metric <- factor(plot_data$metric, levels = metrics)
  rownames(plot_data) <- NULL

  default_colors <- c(
    "sensitivity" = "#E41A1C",
    "specificity" = "#377EB8",
    "auc" = "#4DAF4A",
    "f1" = "#984EA3",
    "precision" = "#FF7F00",
    "npv" = "#A65628"
  )
  colors <- default_colors[metrics]
  missing_colors <- metrics[is.na(colors)]
  if (length(missing_colors)) {
    gray_vals <- grDevices::gray.colors(length(missing_colors), start = 0.3, end = 0.7)
    names(gray_vals) <- missing_colors
    colors[missing_colors] <- gray_vals
  }

  if (facet_by == "metric") {
    if (has_method) {
      p <- ggplot2::ggplot(plot_data,
                           ggplot2::aes(x = .data$cohort, y = .data$value, fill = .data$method)) +
        ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.7)
    } else {
      p <- ggplot2::ggplot(plot_data,
                           ggplot2::aes(x = .data$cohort, y = .data$value, fill = .data$metric)) +
        ggplot2::geom_col(width = 0.7, show.legend = FALSE) +
        ggplot2::scale_fill_manual(values = colors)
    }
    p <- p + ggplot2::facet_wrap(~ metric, scales = "fixed")
  } else {
    if (has_method) {
      p <- ggplot2::ggplot(plot_data,
                           ggplot2::aes(x = .data$metric, y = .data$value, fill = .data$method)) +
        ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.7)
    } else {
      p <- ggplot2::ggplot(plot_data,
                           ggplot2::aes(x = .data$metric, y = .data$value, fill = .data$metric)) +
        ggplot2::geom_col(width = 0.7) +
        ggplot2::scale_fill_manual(values = colors)
    }
    p <- p + ggplot2::facet_wrap(~ cohort)
  }

  if (show_values) {
    if (has_method) {
      p <- p + ggplot2::geom_text(
        ggplot2::aes(label = sprintf("%.2f", .data$value), group = .data$method),
        position = ggplot2::position_dodge(width = 0.8),
        vjust = -0.3,
        size = value_size
      )
    } else {
      p <- p + ggplot2::geom_text(
        ggplot2::aes(label = sprintf("%.2f", .data$value)),
        vjust = -0.3,
        size = value_size
      )
    }
  }

  p <- p +
    ggplot2::labs(
      x = if (facet_by == "metric") "Cohort" else "Metric",
      y = "Value",
      title = title,
      subtitle = subtitle
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid.major.x = ggplot2::element_blank(),
      legend.position = if (has_method) "bottom" else "none"
    )

  if (!is.null(y_limits)) {
    p <- p + ggplot2::coord_cartesian(ylim = y_limits)
  }

  p
}