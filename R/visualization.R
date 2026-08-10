#' Animate Pareto Front Evolution Alongside Feature Inclusion
#'
#' Two-panel animated GIF showing how the rank-1 (Pareto) population evolves
#' across NSGA generations: the left panel plots the Pareto front in
#' objective space, the right panel plots inclusion frequency of base
#' features in that generation's Pareto solutions. Both axes (Pareto-panel
#' x/y, frequency-panel features) are held fixed across frames so the
#' animation is comparable frame-to-frame.
#'
#' Requires the optimization to have been run with `record_history = TRUE`.
#' The set of features shown on the right panel is the top-N base features
#' by total appearances across the rank-1 population of all recorded
#' generations; features not present in a given generation drop to zero.
#'
#' @param optimization_result An `OptimizationResult` from [optimize_panel()].
#' @param x_metric,y_metric Objective names for the Pareto panel axes.
#'   Defaults `"sensitivity"` / `"specificity"`.
#' @param top_n_features Integer; number of features to display on the
#'   right panel. Default `15`.
#' @param generations One of: `NULL` (all recorded generations, default), a
#'   single integer `n` giving the number of evenly-spaced frames to render,
#'   or an integer vector of specific generation numbers.
#' @param fps Numeric; animation frames per second. Default `4`.
#' @param width,height Frame dimensions in pixels. Default `1000` x `500`.
#' @param point_size Point size in the Pareto panel. Default `3`.
#' @param viridis_option Viridis palette for coloring by panel size. Default
#'   `"plasma"`.
#' @param path Optional file path to write the GIF to. Default `NULL`
#'   (in-memory only).
#' @return A `magick-image` representing the animated GIF.
#' @export
#' @seealso [optimize_panel()], [nsga_history()], [hypervolume_history()]
#' @examples
#' \dontrun{
#' opt <- optimize_panel(x, y, record_history = TRUE)
#' anim <- plot_pareto_evolution(opt)
#' anim  # display in RStudio viewer
#' plot_pareto_evolution(opt, path = "evolution.gif")
#' }
plot_pareto_evolution <- function(optimization_result,
                                  x_metric = "sensitivity",
                                  y_metric = "specificity",
                                  top_n_features = 15L,
                                  generations = NULL,
                                  fps = 4,
                                  width = 1000,
                                  height = 500,
                                  point_size = 3,
                                  viridis_option = "plasma",
                                  path = NULL) {
  for (pkg in c("ggplot2", "patchwork", "magick")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required for plot_pareto_evolution(). ",
           "Install with install.packages('", pkg, "').", call. = FALSE)
    }
  }
  if (!inherits(optimization_result, "OptimizationResult")) {
    stop("`optimization_result` must be an OptimizationResult.", call. = FALSE)
  }

  hist <- nsga_history(optimization_result)
  if (!is.data.frame(hist) || !nrow(hist)) {
    stop("No NSGA history recorded. Re-run optimize_panel() with ",
         "record_history = TRUE.", call. = FALSE)
  }
  for (col in c(x_metric, y_metric, "is_pareto", "generation",
                "base_features", "n_features")) {
    if (!col %in% names(hist)) {
      stop("Column '", col, "' missing from nsga_history().", call. = FALSE)
    }
  }

  pareto_hist <- hist[hist$is_pareto, , drop = FALSE]
  if (!nrow(pareto_hist)) {
    stop("No rank-1 individuals in recorded history.", call. = FALSE)
  }

  all_gens <- sort(unique(pareto_hist$generation))
  if (is.null(generations)) {
    gens <- all_gens
  } else if (length(generations) == 1L && is.numeric(generations) &&
             !is.na(generations) && generations >= 1L &&
             generations == as.integer(generations)) {
    n_target <- as.integer(generations)
    if (n_target >= length(all_gens)) {
      gens <- all_gens
    } else {
      idx <- unique(round(seq.int(1, length(all_gens), length.out = n_target)))
      gens <- all_gens[idx]
    }
  } else {
    gens <- intersect(as.integer(generations), all_gens)
    if (!length(gens)) {
      stop("None of the requested generations were recorded.", call. = FALSE)
    }
  }

  all_feats <- unlist(pareto_hist$base_features, use.names = FALSE)
  if (!length(all_feats)) {
    stop("Recorded Pareto solutions contain no base features.", call. = FALSE)
  }
  feat_tab <- sort(table(all_feats), decreasing = TRUE)
  top_n <- min(as.integer(top_n_features), length(feat_tab))
  top_feats <- names(feat_tab)[seq_len(top_n)]

  bounded <- c("sensitivity", "specificity", "auc")
  axis_range <- function(metric) {
    if (metric %in% bounded) return(c(0, 1))
    rng <- range(pareto_hist[[metric]], finite = TRUE)
    pad <- max(diff(rng) * 0.05, .Machine$double.eps)
    c(rng[1L] - pad, rng[2L] + pad)
  }
  xlim <- axis_range(x_metric)
  ylim <- axis_range(y_metric)
  n_feat_range <- range(pareto_hist$n_features, finite = TRUE)

  frames <- magick::image_graph(width = width, height = height, res = 96)
  tryCatch({
    for (g in gens) {
      pareto_g <- pareto_hist[pareto_hist$generation == g, , drop = FALSE]

      p_pareto <- ggplot2::ggplot(
        pareto_g,
        ggplot2::aes(x = .data[[x_metric]],
                     y = .data[[y_metric]],
                     color = .data[["n_features"]])
      ) +
        ggplot2::geom_point(size = point_size, alpha = 0.85) +
        ggplot2::scale_color_viridis_c(
          name = "Panel size",
          option = viridis_option,
          limits = n_feat_range
        ) +
        ggplot2::coord_cartesian(xlim = xlim, ylim = ylim) +
        ggplot2::labs(
          x = x_metric, y = y_metric,
          title = sprintf("Pareto front - generation %d", g)
        ) +
        ggplot2::theme_bw(base_size = 12)

      feat_in_gen <- unlist(pareto_g$base_features, use.names = FALSE)
      counts <- as.integer(table(factor(feat_in_gen, levels = top_feats)))
      freq_df <- data.frame(
        feature = factor(top_feats, levels = rev(top_feats)),
        proportion = counts / max(nrow(pareto_g), 1L),
        stringsAsFactors = FALSE
      )

      p_freq <- ggplot2::ggplot(
        freq_df,
        ggplot2::aes(x = .data[["proportion"]], y = .data[["feature"]])
      ) +
        ggplot2::geom_col(fill = "#377EB8", width = 0.7) +
        ggplot2::scale_x_continuous(
          limits = c(0, 1),
          labels = function(z) paste0(round(z * 100), "%"),
          expand = ggplot2::expansion(mult = c(0, 0.02))
        ) +
        ggplot2::labs(
          x = "Inclusion frequency",
          y = NULL,
          title = sprintf("Top-%d Pareto features (%d solutions)",
                          length(top_feats), nrow(pareto_g))
        ) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
          panel.grid.major.y = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank()
        )

      print(patchwork::wrap_plots(p_pareto, p_freq, ncol = 2L))
    }
  }, finally = grDevices::dev.off())

  anim <- magick::image_animate(frames, fps = fps, optimize = FALSE)
  if (!is.null(path)) {
    magick::image_write(anim, path = path)
  }
  anim
}

#' Plot Pareto Front on Held-Out Data
#'
#' Fits and evaluates every Pareto-optimal solution from [optimize_panel()] on
#' held-out validation data and displays the trade-off between two metrics as a
#' scatter plot. The underlying data frame is accessible via `$data` on the
#' returned ggplot object.
#'
#' @param optimization_result An `OptimizationResult` from [optimize_panel()].
#' @param x Held-out validation data: a matrix, data.frame,
#'   `SummarizedExperiment`, or list of such objects. When `NULL` (default),
#'   automatically uses held-out data stored by
#'   [optimize_panel_transferable()].
#' @param y Held-out validation labels: a factor (or list of factors when `x`
#'   is a list). When `NULL` (default), automatically uses held-out labels
#'   stored by [optimize_panel_transferable()].
#' @param cohort Optional cohort factor for single-matrix `x`; passed to
#'   [evaluate_pareto_solutions()].
#' @param objectives Optional objective override; passed to
#'   [evaluate_pareto_solutions()].
#' @param x_metric Character; column name for the x-axis. Any column of
#'   [evaluate_pareto_solutions()] output, i.e. `"sensitivity"`,
#'   `"specificity"`, `"auc"`, `"n_features"`, `"n_base_features"`, or an
#'   `obj_`-prefixed optimisation objective such as `"obj_auc"`. Default
#'   `"sensitivity"`.
#' @param y_metric Character; column name for the y-axis. Same options as
#'   `x_metric`. Default `"specificity"`.
#' @param color_by Column name to map to point colour (e.g. `"n_features"`,
#'   `"auc"`, an `obj_` column), or `NULL` to disable. Default
#'   `"n_features"`.
#' @param point_size Numeric; point size. Default `4`.
#' @param point_alpha Numeric; point transparency. Default `0.7`.
#' @param viridis_option Character; viridis palette name. Default `"plasma"`.
#' @param title Character or `NULL`; plot title. Auto-generated when `NULL`.
#' @param xlab Character or `NULL`; x-axis label. Defaults to `x_metric`.
#' @param ylab Character or `NULL`; y-axis label. Defaults to `y_metric`.
#' @param color_label Character or `NULL`; color legend label. Defaults to
#'   `color_by`.
#' @param connect Logical; if `TRUE`, connect points with a line sorted by the
#'   x-axis metric. Default `FALSE`.
#' @param on_error One of `"warn"` or `"stop"`. Controls behaviour when a
#'   solution fails evaluation. `"warn"` (default) skips the solution with a
#'   warning; `"stop"` raises an error immediately.
#' @param regularized Logical; passed to [fit_panel()]. Default `TRUE`.
#' @param verbose Logical; print progress messages. Default `interactive()`.
#' @return A `ggplot` object. The underlying metrics data frame is accessible
#'   via `$data`.
#' @export
#' @seealso [evaluate_pareto_solutions()], [optimize_panel()], [fit_panel()],
#'   [evaluate_panel()]
#' @examples
#' \dontrun{
#' opt <- optimize_panel(x_train, y_train, objectives = define_objectives())
#' p <- plot_pareto_front(opt, x_test, y_test)
#' p          # display plot
#' p$data     # access metrics data.frame
#' }
plot_pareto_front <- function(optimization_result,
                              x = NULL,
                              y = NULL,
                              cohort = NULL,
                              objectives = NULL,
                              x_metric = "sensitivity",
                              y_metric = "specificity",
                              color_by = "n_features",
                              point_size = 4,
                              point_alpha = 0.7,
                              viridis_option = "plasma",
                              title = NULL,
                              xlab = NULL,
                              ylab = NULL,
                              color_label = NULL,
                              connect = FALSE,
                              on_error = c("warn", "stop"),
                              regularized = TRUE,
                              verbose = interactive()) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plot_pareto_front()", call. = FALSE)
  }
  plot_df <- evaluate_pareto_solutions(
    optimization_result,
    x = x,
    y = y,
    cohort = cohort,
    objectives = objectives,
    on_error = on_error,
    regularized = regularized,
    verbose = verbose
  )

  requested <- c(x_metric, y_metric, color_by)
  missing_cols <- setdiff(requested, names(plot_df))
  if (length(missing_cols)) {
    stop("Column(s) not available in the evaluated Pareto data: ",
         paste(missing_cols, collapse = ", "), ". Available: ",
         paste(names(plot_df), collapse = ", "), ".", call. = FALSE)
  }

  # Build plot
  if (is.null(title)) {
    title <- sprintf("Pareto Front (%d solutions)", nrow(plot_df))
  }

  aes_args <- list(x = as.name(x_metric), y = as.name(y_metric))
  if (!is.null(color_by)) aes_args$color <- as.name(color_by)
  p <- ggplot2::ggplot(plot_df, do.call(ggplot2::aes, aes_args)) +
    ggplot2::geom_point(size = point_size, alpha = point_alpha)

  if (connect) {
    sorted_df <- plot_df[order(plot_df[[x_metric]]), , drop = FALSE]
    p <- p + ggplot2::geom_line(
      data = sorted_df,
      ggplot2::aes(x = .data[[x_metric]], y = .data[[y_metric]]),
      inherit.aes = FALSE,
      alpha = 0.5
    )
  }

  if (!is.null(color_by) && is.numeric(plot_df[[color_by]])) {
    color_name <- if (!is.null(color_label)) color_label else color_by
    p <- p + ggplot2::scale_color_viridis_c(
      name = color_name, option = viridis_option
    )
  }

  bounded <- c("sensitivity", "specificity", "auc")
  if (x_metric %in% bounded && y_metric %in% bounded) {
    p <- p + ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
  }

  x_label <- if (!is.null(xlab)) xlab else x_metric
  y_label <- if (!is.null(ylab)) ylab else y_metric

  p + ggplot2::labs(x = x_label, y = y_label, title = title) +
    ggplot2::theme_bw(base_size = 12)
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

  top_n <- .validate_positive_integer(top_n, "top_n")

  .validate_probability(highlight_threshold, "highlight_threshold",
                        bounds = "closed")

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
