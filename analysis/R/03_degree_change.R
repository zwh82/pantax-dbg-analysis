# Figure 3: mean degree change between Control and Disease networks.
#
# Degree is the normalized (transitive = FALSE) degree, matching the legacy
# script.  Delta = Disease mean degree - Control mean degree; the paired
# Wilcoxon test is computed across cohorts of the same disease.

degree_stats <- function(g) {
  if (is.null(g) || igraph::vcount(g) == 0) {
    return(c(mean = NA_real_, sd = NA_real_, n = 0))
  }
  d <- igraph::degree(g, mode = "all", normalized = TRUE)
  c(
    mean = mean(d, na.rm = TRUE),
    sd = stats::sd(d, na.rm = TRUE),
    n = length(d)
  )
}

load_network <- function(path) {
  if (file.exists(path)) readRDS(path) else NULL
}

collect_degree_change <- function(cfg, outdir_last) {
  disease_cohorts <- get_disease_cohorts(network_dir(cfg))
  records <- list()

  for (disease in names(disease_cohorts)) {
    if (disease %in% cfg$skip_diseases$degree) next
    for (cohort in disease_cohorts[[disease]]) {
      g_control_path <- file.path(
        network_dir(cfg), cohort, outdir_last, "g_Control_full_network.rds"
      )
      disease_net <- cohort_disease_token(cohort)
      g_disease_path <- file.path(
        network_dir(cfg), cohort, outdir_last,
        paste0("g_", disease_net, "_full_network.rds")
      )

      # ERAS has no standard disease network; fall back to the first
      # time-point network (g_D*_full_network.rds).
      if (!file.exists(g_disease_path)) {
        timepoints <- list.files(
          file.path(network_dir(cfg), cohort, outdir_last),
          pattern = "^g_D[0-9]+_full_network\\.rds$",
          full.names = TRUE
        )
        if (length(timepoints) > 0) {
          g_disease_path <- sort(timepoints)[1]
          message("No standard disease network for ", cohort,
                  "; using time point ", basename(g_disease_path))
        }
      }

      g_control <- load_network(g_control_path)
      g_disease <- load_network(g_disease_path)
      if (is.null(g_control) || is.null(g_disease)) {
        message("Skipping ", cohort, ": network files missing")
        next
      }

      sc <- degree_stats(g_control)
      sd_ <- degree_stats(g_disease)
      records[[cohort]] <- data.frame(
        Disease = disease,
        Cohort = cohort,
        Control_MeanDeg = sc["mean"], Control_SD = sc["sd"], Control_N = sc["n"],
        Disease_MeanDeg = sd_["mean"], Disease_SD = sd_["sd"], Disease_N = sd_["n"],
        stringsAsFactors = FALSE,
        row.names = NULL
      )
    }
  }
  bind_rows(records)
}

build_degree_box_plot <- function(deg_df, diseases_use) {
  if (length(diseases_use) == 0) return(NULL)

  box_w <- 0.60
  jit_w <- 0.15
  box_lw <- 0.9
  box_fill_alpha <- 0.60
  staple_w <- 0.5
  pt_size <- 2.6
  pt_alpha <- 0.9
  pt_stroke <- 0.9
  pt_fill <- "white"

  df <- deg_df %>%
    filter(Disease %in% diseases_use) %>%
    mutate(Direction = ifelse(Delta >= 0, "Increase", "Decrease"))

  ord <- df %>%
    group_by(Disease) %>%
    summarise(m = median(Delta), .groups = "drop") %>%
    arrange(desc(m)) %>%
    pull(Disease)
  df <- df %>% mutate(Disease = factor(Disease, levels = ord))

  n_dis <- nlevels(df$Disease)
  base_cols <- ggsci::pal_d3("category20")(20)
  dis_pal <- if (n_dis <= 20) base_cols[seq_len(n_dis)] else colorRampPalette(base_cols)(n_dis)
  names(dis_pal) <- levels(df$Disease)

  ggplot(df, aes(x = Disease, y = Delta)) +
    geom_hline(yintercept = 0, linewidth = 0.5, color = "grey55", linetype = "dashed") +
    geom_boxplot(
      aes(fill = Disease, color = Disease,
          group = interaction(Disease, Direction)),
      width = box_w, alpha = box_fill_alpha,
      linewidth = box_lw, staplewidth = staple_w,
      outlier.shape = NA, position = position_identity()
    ) +
    geom_jitter(
      aes(color = Disease), shape = 21, fill = pt_fill, stroke = pt_stroke,
      width = jit_w, height = 0, size = pt_size, alpha = pt_alpha
    ) +
    scale_fill_manual(values = dis_pal) +
    scale_color_manual(values = dis_pal) +
    scale_y_continuous(
      breaks = scales::pretty_breaks(6),
      expand = expansion(mult = c(0.03, 0.14))
    ) +
    labs(
      x = "Disease",
      y = expression(bold(Delta)~bold(" mean degree centrality (Disease - Control)"))
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", linewidth = 0.8),
      axis.text = element_text(color = "black", size = 13),
      axis.text.x = element_text(color = "black", size = 13, angle = 45, hjust = 1),
      axis.title.y = element_text(color = "black", size = 18, face = "bold"),
      axis.title.x = element_text(color = "black", size = 18, face = "bold"),
      axis.ticks = element_line(color = "black"),
      axis.ticks.length = unit(0.2, "cm")
    )
}

degree_figure_width <- function(diseases_use) {
  n_dis <- length(diseases_use)
  min(14, max(5, n_dis * 0.9 + 2))
}

run_degree_change <- function(cfg) {
  disease_map <- read_disease_map(cfg)
  suffixes <- analysis_suffixes(cfg)
  fig_height <- 7

  for (outdir_last in suffixes) {
    out_dir <- ensure_dir(file.path(
      cfg$results_dir, cfg$paths$degree, paste0("Degree_Change_", outdir_last)
    ))

    deg_df <- collect_degree_change(cfg, outdir_last) %>%
      mutate(
        Delta = Disease_MeanDeg - Control_MeanDeg,
        SE_Delta = sqrt(Disease_SD^2 / Disease_N + Control_SD^2 / Control_N),
        Direction = ifelse(Delta >= 0, "Increase", "Decrease")
      )

    disease_stats <- deg_df %>%
      replace_disease_name(disease_map) %>%
      group_by(Disease) %>%
      summarise(
        n_cohort = n(),
        n_up = sum(Delta > 0),
        n_down = sum(Delta < 0),
        mean_delta = mean(Delta),
        p_paired = if (n() >= 2) tryCatch(
          wilcox.test(Disease_MeanDeg, Control_MeanDeg, paired = TRUE)$p.value,
          error = function(e) NA_real_
        ) else NA_real_,
        .groups = "drop"
      ) %>%
      arrange(desc(n_cohort), Disease)

    overall_p <- tryCatch(
      wilcox.test(deg_df$Disease_MeanDeg, deg_df$Control_MeanDeg, paired = TRUE)$p.value,
      error = function(e) NA_real_
    )
    message(sprintf(
      "Overall paired Wilcoxon p = %.4g; up %d / down %d cohorts",
      overall_p, sum(deg_df$Delta > 0), sum(deg_df$Delta < 0)
    ))

    dir_stats <- deg_df %>%
      mutate(Direction = ifelse(Delta >= 0, "Increase", "Decrease")) %>%
      group_by(Disease, Direction) %>%
      summarise(
        n_cohort = n(),
        mean_absDel = mean(abs(Delta)),
        median_absDel = median(abs(Delta)),
        .groups = "drop"
      )

    write.csv(deg_df, file.path(out_dir, "Degree_Change_percohort.csv"), row.names = FALSE)
    write.csv(disease_stats, file.path(out_dir, "Degree_Change_disease_stats.csv"), row.names = FALSE)
    write.csv(dir_stats, file.path(out_dir, "Degree_Change_direction_stats.csv"), row.names = FALSE)

    multi_dis <- disease_stats %>% filter(n_cohort > 1) %>% pull(Disease)
    single_dis <- disease_stats %>% filter(n_cohort == 1) %>% pull(Disease)

    p_main <- build_degree_box_plot(deg_df, multi_dis)
    if (!is.null(p_main)) {
      ggsave(file.path(cfg$results_dir, cfg$paths$final, "fig_net_3.pdf"),
             p_main, width = 10, height = fig_height)
    }

    p_supp <- build_degree_box_plot(deg_df, single_dis)
    if (!is.null(p_supp)) {
      ggsave(file.path(out_dir, "Degree_Change_SUPP.png"),
             p_supp, width = degree_figure_width(single_dis), height = fig_height,
             limitsize = FALSE)
    }
  }

  invisible(TRUE)
}
