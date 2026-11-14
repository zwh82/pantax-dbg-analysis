# Per-cohort pipeline: input, prevalence filtering, networks, comparisons.

get_mode <- function(x) {
  x_valid <- x[x != "absent"]
  if (length(x_valid) == 0) return("absent")
  ux <- unique(x_valid)
  ux[which.max(tabulate(match(x_valid, ux)))]
}

read_matrix_file <- function(path) {
  df <- fread(path, sep = "\t", header = TRUE)
  df <- as.data.frame(df)
  rownames(df) <- trimws(as.character(df[[1]]))
  df[, -1, drop = FALSE]
}

read_metadata_file <- function(path) {
  metadata <- fread(path, header = TRUE)
  as.data.frame(metadata)
}

abundance_threshold_args <- function(cfg) {
  th <- cfg$abundance$thresholds
  list(
    high_threshold = th$high,
    low_threshold_min = th$low_min,
    low_threshold_max = th$low_max,
    very_low_min = th$very_low_min,
    very_low_max = th$very_low_max
  )
}

load_abundance <- function(row, cohort, cfg) {
  if (!cfg$input$merged) {
    relative_abund <- row$species_relative_abund
    if (cohort %in% c("PRJNA398089_598_CD", "PRJNA398089_375_UC")) {
      files <- unlist(str_split(relative_abund, ","))
      return(merge_and_filter_abundance(files, prevalence_threshold = 0))
    }
    abundance <- read_matrix_file(relative_abund)
    if (cfg$input$floor > 0) {
      abundance[abundance < cfg$input$floor] <- 0
    }
    return(abundance)
  }

  if (identical(cfg$input$data_type, "read_count")) {
    files <- unlist(strsplit(row$all_read_count, ","))
  } else {
    files <- unlist(strsplit(row$all_species_abund, ","))
  }
  merge_and_filter_abundance(files, prevalence_threshold = cfg$input$prevalence)
}

load_read_count <- function(row, cohort, cfg, common_samples) {
  if (cohort %in% c("PRJNA398089_598_CD", "PRJNA398089_375_UC")) {
    files <- unlist(str_split(row$species_abund, ","))
    read_count <- merge_and_filter_abundance(files, prevalence_threshold = 0)
  } else {
    read_count <- read_matrix_file(row$species_abund)
  }
  read_count[common_samples, , drop = FALSE]
}

subset_metadata <- function(metadata, cohort, cfg) {
  if (!cfg$input$merged) {
    if (cohort %in% c("PRJNA398089_598_CD", "PRJNA398089_375_UC")) {
      metadata <- subset(metadata, Cohort %in% c("PRJNA398089_363_Health", cohort))
      metadata$Cohort <- cohort
    } else if (cohort %in% c("PRJNA993675_474_UC", "PRJNA993675_631_CD")) {
      cohort <- "PRJNA993675_960_IBD"
      metadata <- subset(metadata, Cohort == cohort)
    } else if (!cfg$runs$network_only && !cfg$runs$core_network_only) {
      metadata <- subset(metadata, Cohort == cohort)
    }
  }
  metadata
}

determine_groups <- function(metadata, disease, cohort, cfg) {
  if (cfg$runs$eras) {
    metadata$Group <- ifelse(metadata$Group %in% "D0", "Control", metadata$Group)
    groups <- sort(unique(metadata$Group))
    return(list(metadata = metadata, groups = c("Control", setdiff(groups, "Control"))))
  }
  if (cohort == "PRJNA398089_IBD" && !cfg$runs$prjna398089_ibd_not) {
    groups <- sort(unique(metadata$Group))
    return(list(metadata = metadata, groups = c("Control", setdiff(groups, "Control"))))
  }
  if (cfg$runs$network_only || cfg$runs$core_network_only) {
    return(list(metadata = metadata, groups = "Control"))
  }
  if (disease == "IBD") {
    metadata <- metadata[!metadata$Group %in% "IndeterminateColitis", ]
  }
  groups <- sort(unique(metadata$Group))
  list(metadata = metadata, groups = c("Control", disease, setdiff(groups, c("Control", disease))))
}

compute_group_categories <- function(abundance, metadata, groups, cfg) {
  args <- abundance_threshold_args(cfg)
  categories <- list()
  for (grp in groups) {
    idx_grp <- which(metadata$Group == grp)
    if (length(idx_grp) > 0) {
      mat_grp <- abundance[idx_grp, , drop = FALSE]
      categories[[grp]] <- sapply(colnames(mat_grp), function(species) {
        do.call(assign_abundance_category, c(list(mat_grp[, species]), args))
      })
    } else {
      categories[[grp]] <- setNames(rep("absent", ncol(abundance)), colnames(abundance))
    }
  }
  categories
}

load_differential_markers <- function(maaslin_file) {
  if (!is.na(maaslin_file) && file.exists(maaslin_file)) {
    df_maaslin <- read.table(maaslin_file, sep = "\t", header = TRUE,
                             stringsAsFactors = FALSE)
  } else {
    warning(paste("Maaslin3 结果文件未找到:", maaslin_file, "已返回空数据框。"))
    df_maaslin <- data.frame(
      feature = character(), metadata = character(), value = character(),
      coef = numeric(), stderr = numeric(), pval = numeric(), qval = numeric(),
      stringsAsFactors = FALSE
    )
  }

  sig_df_list <- list()
  if ("metadata" %in% colnames(df_maaslin)) {
    comparisons <- unique(df_maaslin$metadata[grepl("Group|Disease", df_maaslin$metadata)])
    for (comp in comparisons) {
      sig_maaslin <- df_maaslin[
        df_maaslin$metadata == comp &
          df_maaslin$qval_joint < 0.25 &
          !is.na(df_maaslin$coef),
      ]
      if (nrow(sig_maaslin) > 0) {
        sig_maaslin$sig_marker <- ifelse(sig_maaslin$coef > 0, "risk", "protective")
        sig_df <- sig_maaslin[, c("feature", "sig_marker")]
        colnames(sig_df) <- c("species_taxid", "sig_marker")
        sig_df$node_label <- sig_df$species_taxid
        sig_df$color <- ifelse(sig_df$sig_marker == "risk", "red", "green")
        sig_df$comparison <- comp
        sig_df_list[[comp]] <- sig_df
      }
    }
  }

  if (length(sig_df_list) == 0) {
    sig_df_list[["Default"]] <- data.frame(
      species_taxid = character(), sig_marker = character(),
      node_label = character(), color = character(), comparison = character(),
      stringsAsFactors = FALSE
    )
  }
  sig_df_list
}

write_abundance_tables <- function(abundance_filt, metadata, groups, outdir,
                                   project_name, cfg) {
  write.table(
    abundance_filt,
    file = file.path(outdir, paste0(cfg$output$flag, project_name,
                                    "_species.prevalence.tsv")),
    sep = "\t", quote = FALSE, col.names = NA
  )

  for (grp in groups) {
    idx <- which(metadata$Group == grp)
    if (length(idx) > 0) {
      sub <- abundance_filt[idx, , drop = FALSE]
      out_file <- file.path(
        outdir,
        paste0(cfg$output$flag, project_name, "_species.prevalence.",
               cfg$input$prevalence, ".", grp, ".tsv")
      )
      write.table(sub, file = out_file, sep = "\t", quote = FALSE, col.names = NA)
      cat("[INFO]", grp, ":", nrow(sub), "samples,", ncol(sub), "species ->", out_file, "\n")
    }
  }
}

build_group_networks <- function(all_groups, abundance_filt, metadata,
                                 group_abundance_categories, sig_df_list,
                                 outdir, project_name, cfg) {
  args <- abundance_threshold_args(cfg)
  results <- list()
  all_metrics <- list()

  for (grp in all_groups) {
    idx <- which(metadata$Group == grp)
    if (length(idx) == 0) {
      cat("[WARNING] No samples for group:", grp, "\n")
      next
    }

    mat_grp <- abundance_filt[idx, , drop = FALSE]

    grp_categories <- if (!is.null(group_abundance_categories[[grp]])) {
      group_abundance_categories[[grp]]
    } else if (!is.null(group_abundance_categories[["All"]])) {
      group_abundance_categories[["All"]]
    } else {
      sapply(colnames(mat_grp), function(species) {
        do.call(assign_abundance_category, c(list(mat_grp[, species]), args))
      })
    }

    grp_sig_df <- NULL
    for (comp_name in names(sig_df_list)) {
      if (grepl(grp, comp_name, ignore.case = TRUE)) {
        grp_sig_df <- sig_df_list[[comp_name]]
        break
      }
    }
    if (is.null(grp_sig_df)) {
      grp_sig_df <- sig_df_list[[1]]
    }

    cat("\n", rep("-", 50), "\n", sep = "")
    cat("Building network for group:", grp, "\n")
    cat("Samples:", nrow(mat_grp), "\n")
    cat("Species:", ncol(mat_grp), "\n")
    cat(rep("-", 50), "\n", sep = "")

    set.seed(cfg$runtime$seed)
    out_file <- file.path(
      outdir,
      paste0(cfg$output$flag, project_name, "_species.prevalence.",
             cfg$input$prevalence, ".", grp, ".tsv")
    )

    result <- build_weighted_network(
      mat_grp, grp, colnames(mat_grp), grp_sig_df, outdir,
      grp_categories, out_file, cfg, cfg$output$flag
    )
    results[[grp]] <- result

    if (!cfg$runs$network_only || cfg$runs$core_network_only) {
      if (!is.null(result$graph) && vcount(result$graph) > 0) {
        all_metrics[[grp]] <- calculate_network_metrics(result$graph, grp)
      } else {
        cat("[WARNING] Network for group", grp, "is empty or NULL\n")
      }
    }
  }

  list(results = results, metrics = all_metrics)
}

stability_table <- function(graphs, group_names) {
  stability_df <- data.frame(
    Group = character(), Stability = numeric(), Sample_Number = numeric(),
    stringsAsFactors = FALSE
  )
  for (grp in group_names) {
    graph <- graphs[[grp]]
    stability_df <- rbind(stability_df, data.frame(
      Group = grp,
      Stability = ifelse(is.null(graph$Stability), NA, graph$Stability),
      Sample_Number = ifelse(is.null(graph$Samples_number), NA, graph$Samples_number),
      stringsAsFactors = FALSE
    ))
  }
  stability_df
}

write_network_metrics <- function(all_metrics, outdir, outflag) {
  if (length(all_metrics) == 0) return(invisible(NULL))
  metrics_df <- do.call(rbind, lapply(names(all_metrics), function(grp) {
    x <- all_metrics[[grp]]
    data.frame(
      Group = x$Group,
      Nodes = x$Nodes,
      Edges = x$Edges,
      Density = round(x$Density, 4),
      AvgDegree = round(x$AvgDegree, 4),
      GlobalClustering = round(x$GlobalClustering, 4),
      AvgPathLength = round(x$AvgPathLength, 4),
      PositiveEdges = x$PositiveEdges,
      NegativeEdges = x$NegativeEdges,
      PosNegRatio = round(x$PosNegRatio, 4),
      MaxDegree = x$MaxDegree,
      MinDegree = x$MinDegree,
      DegreeSD = round(x$DegreeSD, 4),
      stringsAsFactors = FALSE
    )
  }))
  write.csv(metrics_df, file = paste0(outdir, outflag, "network_metrics_summary.csv"),
            row.names = FALSE)
}

process_one_cohort <- function(i, config_data, cfg) {
  disease <- config_data$disease[i]
  row <- config_data[i, , drop = FALSE]
  project_name <- if (cfg$input$merged) row$project_name else row$cohort
  cohort <- if (cfg$input$merged) row$project_name else row$cohort

  outdir <- paste0(resolve_output_dir(cfg), project_name, "/",
                   resolve_output_suffix(cfg), "/")
  outdir <- gsub("%", "", outdir)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(paste0(outdir, "fig_rds"), recursive = TRUE, showWarnings = FALSE)
  dir.create(paste0(outdir, "data"), recursive = TRUE, showWarnings = FALSE)

  cat("\n", rep("=", 60), "\n", sep = "")
  cat("Processing:", i, "/", nrow(config_data), "-", disease, "\n")
  cat("Project:", project_name, "\n")
  cat("Output:", outdir, "\n")
  cat(rep("=", 60), "\n", sep = "")

  if (cfg$runs$network_only || cfg$runs$core_network_only) {
    output_net <- paste0(outdir, "g_Control_full_network.rds")
    if (file.exists(output_net)) {
      cat("Skipping - output already exists:", output_net, "\n")
      return(NULL)
    }
  }

  if (!cfg$input$merged) {
    abundance <- load_abundance(row, cohort, cfg)
    cat("[INFO] abundance samples:", nrow(abundance), "\n")
    if (cfg$input$floor > 0) {
      cat("[INFO] FilterLow ", cfg$input$floor, "\n")
    }

    metadata <- read_metadata_file(row$metadata_file)
    metadata <- subset_metadata(metadata, cohort, cfg)

    metadata$Sample_ID <- trimws(as.character(metadata$Sample_ID))
    metadata$Group <- trimws(as.character(metadata$Standard_Disease))

    grouped <- determine_groups(metadata, disease, cohort, cfg)
    metadata <- grouped$metadata
    all_groups <- grouped$groups
    cat("[INFO] Found groups:", paste(all_groups, collapse = ", "), "\n")
    cat("[INFO] Metadata samples:", nrow(metadata), "\n")

    common_samples <- intersect(rownames(abundance), metadata$Sample_ID)
    if (length(common_samples) == 0) {
      stop("No common samples found!\n")
    }
    abundance <- abundance[common_samples, , drop = FALSE]
    metadata <- metadata[metadata$Sample_ID %in% common_samples, ]
    metadata <- metadata[match(common_samples, metadata$Sample_ID), ]
    cat("[INFO] Row numbers match:", nrow(abundance), "\n")
    cat("Order check:", all(rownames(abundance) == metadata$Sample_ID), "\n")

    group_abundance_categories <- compute_group_categories(abundance, metadata, all_groups, cfg)

    prevalence_threshold <- if (cfg$input$skip_prevalence_filter) {
      0
    } else {
      cfg$input$prevalence
    }

    if (identical(cfg$input$data_type, "read_count")) {
      read_count <- load_read_count(row, cohort, cfg, common_samples)

      if (cfg$input$fill_zero_counts) {
        read_count[read_count == 0] <- 1e-8
      }

      if (prevalence_threshold > 0) {
        filtered <- filter_by_prevalence(read_count, prevalence_threshold,
                                         cfg$input$present)
        abundance_filt <- filtered$matrix
        keep_species <- filtered$keep
        for (grp in names(group_abundance_categories)) {
          group_abundance_categories[[grp]] <- group_abundance_categories[[grp]][keep_species]
        }
        cat("[INFO] Species after prevalence filter (threshold =", prevalence_threshold, "):",
            ncol(abundance_filt), "/", ncol(read_count), "\n")
      } else {
        abundance_filt <- read_count
        cat("[INFO] No prevalence filter applied (Notprefilter = TRUE)\n")
        cat("[INFO] Species:", ncol(abundance_filt), "\n")
      }
    } else {
      if (prevalence_threshold > 0) {
        filtered <- filter_by_prevalence(abundance, prevalence_threshold,
                                         cfg$input$present)
        abundance_filt <- filtered$matrix
        keep_species <- filtered$keep
        for (grp in names(group_abundance_categories)) {
          group_abundance_categories[[grp]] <- group_abundance_categories[[grp]][keep_species]
        }
        cat("[INFO] Species after prevalence filter (threshold =", prevalence_threshold, "):",
            ncol(abundance_filt), "/", ncol(abundance), "\n")
      } else {
        abundance_filt <- abundance
        cat("[INFO] No prevalence filter applied (Notprefilter = TRUE)\n")
        cat("[INFO] Species:", ncol(abundance_filt), "\n")
      }
    }

    cat("[INFO] Species after  filter:", ncol(abundance_filt), "\n")
  } else {
    abundance_filt <- load_abundance(row, NULL, cfg)

    abundance_files <- unlist(strsplit(row$all_species_abund, ","))
    cohort_cats_list <- map(abundance_files, function(path) {
      if (!file.exists(path)) return(NULL)
      cohort_name <- gsub(".*/meta_info/([^/]+)/.*", "\\1", path)
      df <- fread(path, header = TRUE, data.table = FALSE)
      rownames(df) <- df[[1]]
      df <- df[, -1, drop = FALSE]
      species_cats <- apply(df, 2, assign_abundance_category)
      data.frame(Species = names(species_cats),
                 Category = as.character(species_cats),
                 Cohort = cohort_name)
    }) %>% compact()

    all_cats_long <- bind_rows(cohort_cats_list)
    category_matrix <- all_cats_long %>%
      pivot_wider(names_from = Cohort, values_from = Category) %>%
      replace(is.na(.), "absent")

    cohort_cols <- setdiff(colnames(category_matrix), "Species")
    final_stratification <- category_matrix %>%
      rowwise() %>%
      mutate(
        Final_Category = get_mode(c_across(all_of(cohort_cols)))
      ) %>%
      ungroup()

    abundance_category <- setNames(final_stratification$Final_Category,
                                   final_stratification$Species)

    metadata <- read_metadata_file(row$metadata_file)
    colnames(metadata) <- trimws(colnames(metadata))
    metadata$Sample_ID <- trimws(as.character(metadata$Sample_ID))
    cat("[INFO] Metadata samples:", nrow(metadata), "\n")

    if (nrow(abundance_filt) != nrow(metadata)) {
      cat("Row numbers of abundance and metadata do NOT match!\n")
    }

    common_samples <- intersect(rownames(abundance_filt), metadata$Sample_ID)
    if (length(common_samples) == 0) {
      stop("No common samples found!\n")
    }
    abundance_filt <- abundance_filt[common_samples, , drop = FALSE]
    metadata <- metadata[metadata$Sample_ID %in% common_samples, ]
    metadata <- metadata[match(common_samples, metadata$Sample_ID), ]
    cat("[INFO] Row numbers match:", nrow(abundance_filt), "\n")
    cat("Order check:", all(rownames(abundance_filt) == metadata$Sample_ID), "\n")

    all_groups <- sort(unique(metadata$Group))
    group_abundance_categories <- list(All = abundance_category)
    abundance <- abundance_filt
  }

  write_abundance_tables(abundance_filt, metadata, all_groups, outdir,
                         project_name, cfg)
  abundance_filt <- as.matrix(abundance_filt)
  mode(abundance_filt) <- "numeric"

  sig_df_list <- load_differential_markers(row$maaslin3_res_file)

  built <- build_group_networks(
    all_groups, abundance_filt, metadata, group_abundance_categories,
    sig_df_list, outdir, project_name, cfg
  )
  results <- built$results
  all_metrics <- built$metrics

  if (cfg$runs$network_only || cfg$runs$core_network_only) {
    graphs <- list()
    for (grp in names(results)) {
      if (!is.null(results[[grp]]$graph) && vcount(results[[grp]]$graph) > 0) {
        graphs[[grp]] <- results[[grp]]$graph
      }
    }
    sorted_names <- c("Control")
    stability_df <- stability_table(graphs, sorted_names)
    for (grp in sorted_names) {
      graph_plot(graphs[[grp]], outdir, grp, unlist(cfg$abundance$colors),
                 unlist(cfg$abundance$labels), cfg$output$flag)
    }
    write.csv(stability_df, paste0(outdir, "data/stability_df.csv"), row.names = FALSE)
    return(results)
  }

  if (length(results) >= 2) {
    cat("\n", rep("=", 60), "\n", sep = "")
    cat("Multi-group Network Comparison\n")
    cat(rep("=", 60), "\n", sep = "")

    graphs <- list()
    for (grp in names(results)) {
      if (!is.null(results[[grp]]$graph) && vcount(results[[grp]]$graph) > 0) {
        graphs[[grp]] <- results[[grp]]$graph
      }
    }

    if (length(graphs) >= 2) {
      sorted_names <- sort(names(graphs))
      if (cfg$runs$eras) {
        sorted_names <- c("Control", "D4", "D8", "D42", "D180")
      } else if (cohort == "PRJNA398089_IBD" && !cfg$runs$prjna398089_ibd_not) {
        sorted_names <- c("Control", setdiff(sorted_names, c("Control")))
      } else {
        sorted_names <- c("Control", disease, setdiff(sorted_names, c("Control", disease)))
      }
      cat("[INFO] Group order:", paste(sorted_names, collapse = ", "), "\n")

      plot_multi_edges_connect(graphs, sorted_names, outdir, eras = cfg$runs$eras)
      plot_multi_edges(graphs, sorted_names, outdir, eras = cfg$runs$eras)

      stability_df <- stability_table(graphs, sorted_names)
      for (grp in sorted_names) {
        graph_plot(graphs[[grp]], outdir, grp, unlist(cfg$abundance$colors),
                   unlist(cfg$abundance$labels), cfg$output$flag)
      }
      write.csv(stability_df, paste0(outdir, "data/stability_df.csv"), row.names = FALSE)

      all_centrality <- list()
      for (grp in names(results)) {
        if (!is.null(results[[grp]]$hubs$all_metrics)) {
          all_centrality[[grp]] <- results[[grp]]$hubs$all_metrics
        }
      }
      if (length(all_centrality) >= 2) {
        plot_multi_abundance_centrality_boxplots(
          all_centrality, names(all_centrality), outdir, eras = cfg$runs$eras
        )
      }
    }
  }

  write_network_metrics(all_metrics, outdir, cfg$output$flag)
  results
}
