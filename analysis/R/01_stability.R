# Figure 1: network stability versus sample size.
#
# Inputs per cohort:
#   <cohort>/<suffix>/data/stability_df.csv
#   <cohort>/<suffix>/<cohort>_species.prevalence.<prev>.Control.tsv
#   <cohort>/<suffix>/<cohort>_species.prevalence.<prev>.<Disease>.tsv
# Outputs:
#   <results>/stab/stability_by_samplesize_<suffix>.png
#   <results>/stab/stability_overall_<suffix>.png
#   <results>/stab/stability_overall_den<suffix>.png
#   <results>/stab/stability_facet_<suffix>.png
#   <results>/stab/np_overall_<last suffix>.png
#   <results>/scripts/final_figures/fig_net_1.pdf

collect_stability <- function(cfg, suffixes) {
  all_results <- list()
  for (outdir_last in suffixes) {
    disease_cohorts <- get_disease_cohorts(network_dir(cfg))
    for (disease in names(disease_cohorts)) {
      if (disease %in% cfg$skip_diseases$stability) next

      disease_stability <- data.frame()
      for (cohort in disease_cohorts[[disease]]) {
        stability_file <- file.path(
          network_dir(cfg), cohort, outdir_last, "data", "stability_df.csv"
        )
        if (!file.exists(stability_file)) next

        df <- read.csv(stability_file)
        df$Cohort <- cohort
        df$Disease <- disease
        df$Method <- outdir_last
        df$Sample_Group <- ifelse(df$Sample_Number > 15, ">15", "≤15")
        disease_stability <- rbind(disease_stability, df)
      }

      if (nrow(disease_stability) > 0) {
        all_results[[paste(disease, outdir_last, sep = "|")]] <- disease_stability
      }
    }
  }
  all_results
}

stability_disease_order <- function(all_results) {
  unique(vapply(strsplit(names(all_results), "|", fixed = TRUE),
                function(x) x[1], character(1)))
}

make_stability_disease_plot <- function(data, disease) {
  ggplot(data, aes(x = Sample_Group, y = Stability, fill = Sample_Group)) +
    geom_boxplot(alpha = 0.7, width = 0.6, outlier.size = 0.8) +
    geom_jitter(width = 0.2, size = 1.2, alpha = 0.5, color = "black") +
    theme_bw() +
    labs(title = disease, x = "", y = "") +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 8),
      axis.text.y = element_text(size = 7),
      axis.title = element_text(size = 8),
      plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
      legend.position = "none",
      panel.grid.minor = element_blank()
    ) +
    scale_fill_manual(values = c(">15" = "#2E8B57", "≤15" = "#DC143C"))
}

make_stability_overall_box <- function(data, method) {
  ggplot(data, aes(x = Sample_Group, y = Stability, fill = Sample_Group)) +
    geom_boxplot(alpha = 0.7, width = 0.6) +
    geom_jitter(width = 0.2, size = 2, alpha = 0.4, color = "black") +
    theme_bw() +
    labs(
      title = paste("Overall Stability by Sample Size -", method),
      x = "Sample Size Group",
      y = "Stability"
    ) +
    theme(legend.position = "none") +
    scale_fill_manual(values = c(">15" = "#2E8B57", "≤15" = "#DC143C"))
}

make_stability_density <- function(data) {
  ggplot(data, aes(x = Stability)) +
    geom_density(fill = "lightblue", alpha = 0.5, color = "blue", linewidth = 1) +
    theme_bw() +
    labs(x = "Stability", y = "Density") +
    theme(
      panel.border = element_blank(),
      axis.line = element_line(color = "black"),
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
      axis.text.y = element_text(size = 12, face = "bold"),
      axis.title.x = element_text(size = 16, face = "bold"),
      axis.title.y = element_text(size = 16, face = "bold"),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 8),
      legend.position = "none"
    )
}

make_stability_facet_plot <- function(data, method) {
  ggplot(data, aes(x = Sample_Group, y = Stability, fill = Sample_Group)) +
    geom_boxplot(alpha = 0.7, width = 0.6, outlier.size = 0.5) +
    geom_jitter(width = 0.2, size = 0.8, alpha = 0.4, color = "black") +
    facet_wrap(~Disease, scales = "free_y", ncol = 6) +
    theme_bw() +
    labs(
      title = paste("Stability by Sample Size -", method),
      x = "Sample Size Group",
      y = "Stability"
    ) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
      axis.text.y = element_text(size = 6),
      axis.title = element_text(size = 10),
      plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
      strip.text = element_text(size = 7, face = "bold")
    ) +
    scale_fill_manual(
      values = c(">15" = "#2E8B57", "≤15" = "#DC143C"),
      name = "Sample Size"
    )
}

read_prevalence_dimensions <- function(path) {
  read_count <- data.table::fread(path, sep = "\t", header = TRUE)
  read_count <- as.data.frame(read_count)
  rownames(read_count) <- trimws(as.character(read_count[[1]]))
  read_count <- read_count[, -1]
  c(n = unname(nrow(read_count)), p = unname(ncol(read_count)))
}

collect_sample_size_data <- function(cfg, suffixes) {
  n_p_summary <- data.frame()
  for (outdir_last in suffixes) {
    disease_cohorts <- get_disease_cohorts(network_dir(cfg))
    for (disease in names(disease_cohorts)) {
      if (disease %in% cfg$skip_diseases$stability) next

      for (cohort in disease_cohorts[[disease]]) {
        control_data_path <- file.path(
          network_dir(cfg), cohort, outdir_last,
          paste0(cohort, "_species.prevalence.", cfg$prevalence, ".Control.tsv")
        )
        disease_name <- cohort_disease_token(cohort)
        disease_data_path <- file.path(
          network_dir(cfg), cohort, outdir_last,
          paste0(cohort, "_species.prevalence.", cfg$prevalence, ".", disease_name, ".tsv")
        )
        if (!file.exists(control_data_path) || !file.exists(disease_data_path)) next

        n_control <- read_prevalence_dimensions(control_data_path)
        n_disease <- read_prevalence_dimensions(disease_data_path)
        n_p_summary <- rbind(
          n_p_summary,
          data.frame(
            Cohort = cohort,
            Disease = disease,
            Group = c("Control", "Disease"),
            n = c(n_control["n"], n_disease["n"]),
            p = c(n_control["p"], n_disease["p"]),
            stringsAsFactors = FALSE
          )
        )
      }
    }
  }
  n_p_summary
}

run_stability <- function(cfg) {
  suffixes <- analysis_suffixes(cfg)
  all_results <- collect_stability(cfg, suffixes)
  output_dir <- ensure_dir(file.path(cfg$results_dir, cfg$paths$stability))
  diseases <- stability_disease_order(all_results)

  for (method in suffixes) {
    method_plots <- list()
    for (disease in diseases) {
      key <- paste(disease, method, sep = "|")
      if (!key %in% names(all_results)) next
      method_plots[[disease]] <- make_stability_disease_plot(all_results[[key]], disease)
    }

    if (length(method_plots) > 0) {
      n_plots <- length(method_plots)
      n_col <- 6
      n_row <- ceiling(n_plots / n_col)
      combined <- wrap_plots(method_plots, ncol = n_col, nrow = n_row) +
        plot_annotation(
          title = paste("Stability by Sample Size -", method),
          subtitle = paste("Total diseases:", n_plots,
                           "| Sample Size >15 (green) vs ≤15 (red)"),
          theme = theme(
            plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
            plot.subtitle = element_text(size = 12, hjust = 0.5)
          )
        )
      ggsave(
        file.path(output_dir, paste0("stability_by_samplesize_", method, ".png")),
        combined,
        width = n_col * 2.5,
        height = n_row * 2.5,
        dpi = 300,
        limitsize = FALSE
      )
      message("Saved: ", file.path(output_dir, paste0("stability_by_samplesize_", method, ".png")))
    }
  }

  # The overall plots are written per method, but the two cross-method
  # outputs below (np scatter and fig_net_1) historically use the last
  # method name and the last density plot that had data.
  p_overall <- NULL
  for (method in suffixes) {
    all_diseases_data <- data.frame()
    for (disease in diseases) {
      key <- paste(disease, method, sep = "|")
      if (!key %in% names(all_results)) next
      data <- all_results[[key]]
      data$Disease <- disease
      all_diseases_data <- rbind(all_diseases_data, data)
    }
    if (nrow(all_diseases_data) == 0) next

    p_box <- make_stability_overall_box(all_diseases_data, method)
    ggsave(file.path(output_dir, paste0("stability_overall_", method, ".png")),
           p_box, width = 6, height = 6, dpi = 300)

    p_overall <- make_stability_density(all_diseases_data)
    ggsave(file.path(output_dir, paste0("stability_overall_den", method, ".png")),
           p_overall, width = 6, height = 6, dpi = 300)

    group_gt15 <- all_diseases_data %>% filter(Sample_Group == ">15") %>% pull(Stability)
    group_le15 <- all_diseases_data %>% filter(Sample_Group == "≤15") %>% pull(Stability)
    if (length(group_gt15) > 0 && length(group_le15) > 0) {
      wilcox_result <- wilcox.test(group_gt15, group_le15)
      message(
        "Overall Wilcoxon (", method, "): >15 n=", length(group_gt15),
        " vs <=15 n=", length(group_le15),
        ", p=", signif(wilcox_result$p.value, 4),
        "; medians ", signif(median(group_gt15), 4), " / ", signif(median(group_le15), 4)
      )
    }
  }

  for (method in suffixes) {
    all_diseases_data <- data.frame()
    for (disease in diseases) {
      key <- paste(disease, method, sep = "|")
      if (!key %in% names(all_results)) next
      data <- all_results[[key]]
      data$Disease <- disease
      all_diseases_data <- rbind(all_diseases_data, data)
    }
    if (nrow(all_diseases_data) == 0) next

    p_facet <- make_stability_facet_plot(all_diseases_data, method)
    ggsave(file.path(output_dir, paste0("stability_facet_", method, ".png")),
           p_facet, width = 15, height = 18, dpi = 300, limitsize = FALSE)
  }

  n_p_summary <- collect_sample_size_data(cfg, suffixes)

  last_method <- tail(suffixes, 1)
  p_scatter <- ggplot(n_p_summary, aes(x = p, y = n)) +
    geom_point(size = 3, alpha = 0.7, color = "#db72fb") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 1.3,
                color = "black") +
    stat_function(fun = function(x) log(x), linetype = "dashed", color = "red",
                  linewidth = 1.3, n = 500) +
    labs(x = "Number of species (p)", y = "Sample size (n)") +
    theme_bw() +
    theme(
      panel.border = element_blank(),
      axis.line = element_line(color = "black"),
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
      axis.text.y = element_text(size = 12, face = "bold"),
      axis.title.x = element_text(size = 16, face = "bold"),
      axis.title.y = element_text(size = 16, face = "bold"),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 8),
      legend.position = "none"
    )
  ggsave(file.path(output_dir, paste0("np_overall_", last_method, ".png")),
         p_scatter, width = 6, height = 6, dpi = 300)

  if (is.null(p_overall)) {
    stop("no stability data available for fig_net_1")
  }
  p_merge <- p_scatter + p_overall +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(size = 16, face = "bold"))
  ggsave(file.path(cfg$results_dir, cfg$paths$final, "fig_net_1.pdf"),
         p_merge, width = 12, height = 6)

  invisible(list(scatter = p_scatter, overall = p_overall, merged = p_merge))
}
