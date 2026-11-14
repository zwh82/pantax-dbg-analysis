# Consensus-hub abundance composition.
#
# Meta version recomputes consensus hubs from every cohort's network; the
# summary version plots an externally produced per-cohort hub table.  Both
# share the same layout and colours.

map_hub_abundance <- function(x) {
  if (is.null(x)) return(character(0))
  x <- tolower(as.character(x))
  dplyr::case_when(
    x == "high" ~ "abundant",
    x == "low" ~ "rare",
    x == "very_low" ~ "rare",
    TRUE ~ "Unknown_Other"
  )
}

consensus_hub_table <- function(rds_file, hub_cfg) {
  g <- as_undirected(readRDS(rds_file), mode = "collapse")
  aw <- if ("weight" %in% edge_attr_names(g)) abs(E(g)$weight) else rep(1, ecount(g))

  set.seed(hub_cfg$seed)
  cent <- data.frame(
    Node = V(g)$name,
    Degree = degree(g),
    Betweenness = betweenness(g, weights = aw),
    Closeness = closeness(g, weights = 1 / (aw + 1e-10)),
    Eigenvector = eigen_centrality(g, weights = aw)$vector,
    PageRank = page_rank(g, weights = aw)$vector,
    stringsAsFactors = FALSE
  )

  norm <- function(v) {
    r <- max(v, na.rm = TRUE) - min(v, na.rm = TRUE)
    if (r == 0) rep(0, length(v)) else (v - min(v, na.rm = TRUE)) / r
  }

  votes <- 0
  for (m in c("Degree", "Betweenness", "Closeness", "Eigenvector", "PageRank")) {
    votes <- votes +
      (norm(cent[[m]]) >= quantile(norm(cent[[m]]), 1 - hub_cfg$top_percent,
                                   na.rm = TRUE))
  }
  cent$Is_Hub <- (votes >= hub_cfg$vote_threshold) & (cent$Degree >= hub_cfg$min_degree)

  ab_raw <- if (!is.null(V(g)$abundance_category)) {
    V(g)$abundance_category
  } else {
    rep(NA, vcount(g))
  }
  cent$Abundance <- map_hub_abundance(ab_raw)
  cent[cent$Is_Hub, c("Node", "Abundance")]
}

collect_consensus_hub <- function(cfg) {
  network_suffix <- analysis_suffixes(cfg)[1]
  disease_cohorts <- get_disease_cohorts(network_dir(cfg))
  raw_list <- list()
  n_ok <- 0

  for (disease in names(disease_cohorts)) {
    if (disease %in% cfg$skip_diseases$hub) next

    for (cohort in disease_cohorts[[disease]]) {
      if (cohort %in% c("all_country") || !grepl("_", cohort)) next

      rc <- file.path(network_dir(cfg), cohort, network_suffix)
      ctrl_rds <- file.path(rc, "g_Control_full_network.rds")
      if (!file.exists(ctrl_rds)) next

      net_files <- list.files(rc, pattern = "^g_.*_full_network\\.rds$")
      dis_files <- setdiff(net_files, "g_Control_full_network.rds")
      if (length(dis_files) == 0) next
      case_rds <- file.path(rc, dis_files[1])

      message(sprintf("%s (Control vs %s)", cohort, disease))
      ht_c <- tryCatch(
        consensus_hub_table(ctrl_rds, cfg$hub),
        error = function(e) { message("  Control read failed"); NULL }
      )
      ht_d <- tryCatch(
        consensus_hub_table(case_rds, cfg$hub),
        error = function(e) { message("  Disease read failed"); NULL }
      )
      if (is.null(ht_c) || is.null(ht_d)) next
      n_ok <- n_ok + 1

      raw_list[[cohort]] <- bind_rows(
        ht_c %>% mutate(Group = "Control"),
        ht_d %>% mutate(Group = "Disease")
      ) %>%
        count(Group, Abundance) %>%
        mutate(Cohort = cohort, Disease = disease)
    }
  }

  message("Valid cohorts: ", n_ok)
  if (n_ok < 2) stop("fewer than 2 valid cohorts, cannot pair Control and Disease")

  master_df <- bind_rows(raw_list) %>%
    distinct() %>%
    filter(Abundance != "Unknown_Other") %>%
    mutate(
      Group = factor(Group, levels = c("Control", "Disease")),
      Abundance = factor(Abundance, levels = cfg$hub$levels)
    ) %>%
    complete(Cohort, Group, Abundance, fill = list(n = 0)) %>%
    group_by(Cohort, Group) %>%
    mutate(Total_Hub = sum(n)) %>%
    ungroup() %>%
    mutate(Proportion = ifelse(Total_Hub == 0, 0, n / Total_Hub)) %>%
    mutate(
      Group = factor(Group, levels = c("Control", "Disease")),
      Abundance = factor(Abundance, levels = cfg$hub$levels)
    )

  # complete() inserts rows with a missing Disease label; recover it by cohort.
  cohort_disease_map <- master_df %>%
    filter(!is.na(Disease)) %>%
    distinct(Cohort, Disease)
  master_df <- master_df %>%
    select(-Disease) %>%
    left_join(cohort_disease_map, by = "Cohort")

  master_df <- replace_disease_name2(master_df, read_disease_map(cfg), cols = c("Disease"))

  disease_counts <- master_df %>%
    distinct(Disease, Cohort) %>%
    count(Disease, name = "n_cohort") %>%
    arrange(desc(n_cohort))
  disease_order <- disease_counts$Disease

  cohort_order <- master_df %>%
    filter(Group == "Disease") %>%
    distinct(Cohort, Total_Hub) %>%
    arrange(desc(Total_Hub)) %>%
    pull(Cohort) %>%
    as.character()

  master_df <- master_df %>%
    mutate(Cohort = factor(Cohort, levels = cohort_order))

  list(
    master_df = master_df,
    disease_counts = disease_counts,
    disease_order = disease_order,
    cohort_order = cohort_order
  )
}

make_hub_positions <- function(master_df, diseases_use, disease_seq, cohort_order,
                               hub_cfg) {
  ds <- disease_seq[disease_seq %in% diseases_use]
  pos_list <- list()
  centers <- numeric(0)
  x_cursor <- 1

  for (d in ds) {
    ch <- master_df %>%
      filter(Disease == d) %>%
      distinct(Cohort) %>%
      pull(Cohort) %>%
      as.character()
    ch <- ch[order(match(ch, cohort_order))]
    this_x <- x_cursor + (seq_along(ch) - 1) * hub_cfg$step_within
    pos_list[[d]] <- data.frame(
      Cohort = ch, xpos = this_x, Disease = d, stringsAsFactors = FALSE
    )
    centers[d] <- mean(this_x)
    x_cursor <- max(this_x) + hub_cfg$step_between
  }

  list(pos = bind_rows(pos_list), centers = centers)
}

build_hub_grouped_plot <- function(master_df, diseases_use, disease_seq,
                                   cohort_order, hub_cfg, title = NULL) {
  if (length(diseases_use) == 0) return(NULL)

  pinfo <- make_hub_positions(
    master_df, diseases_use, disease_seq, cohort_order, hub_cfg
  )
  posmap <- pinfo$pos
  centers <- pinfo$centers

  df <- master_df %>%
    filter(Disease %in% diseases_use) %>%
    left_join(posmap %>% select(Cohort, xpos), by = "Cohort")
  tot <- df %>% distinct(Cohort, Group, Total_Hub, xpos)

  ggplot(df, aes(x = xpos, y = n, fill = Abundance)) +
    geom_bar(stat = "identity", width = hub_cfg$bar_width) +
    geom_text(
      aes(label = ifelse(n > 0, n, "")),
      position = position_stack(vjust = 0.5),
      size = 2.4, color = "white", fontface = "bold"
    ) +
    geom_text(
      data = tot, aes(x = xpos, y = Total_Hub, label = Total_Hub),
      inherit.aes = FALSE, vjust = -0.4, size = 2.6
    ) +
    facet_wrap(~ Group, ncol = 1) +
    scale_x_continuous(
      breaks = unname(centers), labels = names(centers),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    scale_fill_manual(values = unlist(hub_cfg$colors), drop = FALSE) +
    theme_bw(base_size = 13) +
    labs(title = title, x = "Disease", y = "Hub count", fill = "Abundance Type") +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 11, face = "bold"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      strip.text = element_text(face = "bold", size = 12),
      strip.background = element_rect(fill = "#f0f0f0"),
      legend.position = "bottom"
    )
}

hub_figure_width <- function(master_df, diseases_use) {
  n_bar <- master_df %>%
    filter(Disease %in% diseases_use) %>%
    distinct(Cohort) %>%
    nrow()
  min(16, max(6, n_bar * 0.5 + 2))
}

run_meta_consensus_hub <- function(cfg) {
  out_dir <- ensure_dir(file.path(cfg$results_dir, cfg$paths$hub))
  hub_data <- collect_consensus_hub(cfg)

  master_df <- hub_data$master_df
  disease_counts <- hub_data$disease_counts
  disease_order <- hub_data$disease_order
  cohort_order <- hub_data$cohort_order

  multi_dis <- disease_counts %>% filter(n_cohort > 1) %>% pull(Disease)
  single_dis <- disease_counts %>% filter(n_cohort == 1) %>% pull(Disease)
  main_seq <- disease_order[disease_order %in% multi_dis]
  supp_seq <- sort(single_dis)

  message("Main-figure diseases: ", paste(multi_dis, collapse = ", "))
  message("Supplementary diseases: ", paste(single_dis, collapse = ", "))

  p_main <- build_hub_grouped_plot(
    master_df, multi_dis, main_seq, cohort_order, cfg$hub
  )
  if (!is.null(p_main)) {
    ggsave(
      file.path(out_dir, "Meta_ConsensusHub_Abundance_MAIN.png"),
      p_main,
      width = hub_figure_width(master_df, multi_dis),
      height = cfg$hub$figure_height,
      limitsize = FALSE
    )
  }

  saveRDS(
    list(
      master_df = master_df,
      disease_counts = disease_counts,
      disease_order = disease_order,
      multi_dis = multi_dis,
      single_dis = single_dis,
      main_seq = main_seq,
      supp_seq = supp_seq,
      cohort_order = cohort_order,
      abund_levels = cfg$hub$levels,
      abund_colors = unlist(cfg$hub$colors),
      step_within = cfg$hub$step_within,
      step_between = cfg$hub$step_between,
      bar_w = cfg$hub$bar_width
    ),
    file = file.path(out_dir, "Meta_ConsensusHub_Abundance_plot_data.rds")
  )

  p_supp <- build_hub_grouped_plot(
    master_df, single_dis, supp_seq, cohort_order, cfg$hub
  )
  if (!is.null(p_supp)) {
    ggsave(
      file.path(out_dir, "Meta_ConsensusHub_Abundance_SUPP.png"),
      p_supp,
      width = hub_figure_width(master_df, single_dis),
      height = cfg$hub$figure_height,
      limitsize = FALSE
    )
  }

  invisible(TRUE)
}

run_hub_summary_plot <- function(cfg, in_file = NULL, out_dir = NULL) {
  if (is.null(in_file)) in_file <- cfg$hub_summary$input
  if (is.null(out_dir)) out_dir <- cfg$hub_summary$out_dir
  if (is.null(in_file) || !file.exists(in_file)) {
    stop("hub summary input not found: ", in_file,
         " (pass --input /path/to/ConsensusHub_master_summary.csv)")
  }
  out_dir <- ensure_dir(out_dir)

  raw <- read.csv(in_file, stringsAsFactors = FALSE, check.names = FALSE)
  master_df <- raw %>%
    mutate(Group = ifelse(Group == "Control", "Control", "Disease")) %>%
    transmute(
      Cohort, Disease, Group,
      rare = N_hub_rare,
      abundant = N_hub_abundant
    ) %>%
    pivot_longer(
      cols = c(rare, abundant), names_to = "Abundance", values_to = "n"
    ) %>%
    mutate(
      Group = factor(Group, levels = c("Control", "Disease")),
      Abundance = factor(Abundance, levels = cfg$hub$levels)
    ) %>%
    group_by(Cohort, Group) %>%
    mutate(Total_Hub = sum(n)) %>%
    ungroup()

  disease_counts <- master_df %>%
    distinct(Disease, Cohort) %>%
    count(Disease, name = "n_cohort") %>%
    arrange(desc(n_cohort))
  disease_order <- disease_counts$Disease

  cohort_order <- master_df %>%
    filter(Group == "Disease") %>%
    distinct(Cohort, Total_Hub) %>%
    arrange(desc(Total_Hub)) %>%
    pull(Cohort) %>%
    as.character()
  master_df <- master_df %>%
    mutate(Cohort = factor(Cohort, levels = cohort_order))

  multi_dis <- disease_counts %>% filter(n_cohort > 1) %>% pull(Disease)
  single_dis <- disease_counts %>% filter(n_cohort == 1) %>% pull(Disease)
  main_seq <- disease_order[disease_order %in% multi_dis]
  supp_seq <- sort(single_dis)

  message("Main-figure diseases: ", paste(multi_dis, collapse = ", "))
  message("Supplementary diseases: ", paste(single_dis, collapse = ", "))

  p_main <- build_hub_grouped_plot(
    master_df, multi_dis, main_seq, cohort_order, cfg$hub
  )
  if (!is.null(p_main)) {
    width <- hub_figure_width(master_df, multi_dis)
    ggsave(file.path(out_dir, "ConsensusHub_RareAbundant_MAIN.png"),
           p_main, width = width, height = cfg$hub$figure_height, limitsize = FALSE)
    ggsave(file.path(out_dir, "ConsensusHub_RareAbundant_MAIN.pdf"),
           p_main, width = width, height = cfg$hub$figure_height, limitsize = FALSE)
  }

  p_supp <- build_hub_grouped_plot(
    master_df, single_dis, supp_seq, cohort_order, cfg$hub
  )
  if (!is.null(p_supp)) {
    width <- hub_figure_width(master_df, single_dis)
    ggsave(file.path(out_dir, "ConsensusHub_RareAbundant_SUPP.png"),
           p_supp, width = width, height = cfg$hub$figure_height, limitsize = FALSE)
    ggsave(file.path(out_dir, "ConsensusHub_RareAbundant_SUPP.pdf"),
           p_supp, width = width, height = cfg$hub$figure_height, limitsize = FALSE)
  }

  invisible(TRUE)
}
