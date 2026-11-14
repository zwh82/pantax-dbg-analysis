# Network construction and topology metrics.

calculate_network_metrics <- function(g, group_name) {
  cat("\nCalculating network metrics for", group_name, "...\n")

  metrics <- list()
  metrics$Group <- group_name
  metrics$Nodes <- vcount(g)
  metrics$Edges <- ecount(g)
  metrics$Density <- edge_density(g)
  metrics$AvgDegree <- mean(degree(g))
  metrics$GlobalClustering <- transitivity(g, type = "global")

  g_dist <- g
  epsilon <- 1e-10
  E(g_dist)$distance <- 1 / (E(g_dist)$abs_weight + epsilon)

  if (is_connected(g_dist)) {
    metrics$AvgPathLength <- mean_distance(g_dist, weights = E(g_dist)$distance)
  } else {
    comps <- components(g_dist)
    path_lengths <- c()
    comp_weights <- c()

    for (i in 1:comps$no) {
      comp_nodes <- V(g_dist)[comps$membership == i]
      if (length(comp_nodes) > 1) {
        comp_graph <- induced_subgraph(g_dist, comp_nodes)
        comp_dist <- tryCatch({
          mean_distance(comp_graph, weights = E(comp_graph)$distance)
        }, error = function(e) NA)
        if (!is.na(comp_dist)) {
          path_lengths <- c(path_lengths, comp_dist)
          comp_weights <- c(comp_weights, length(comp_nodes))
        }
      }
    }

    if (length(path_lengths) > 0) {
      metrics$AvgPathLength <- weighted.mean(path_lengths, comp_weights, na.rm = TRUE)
      metrics$LargestComponentSize <- max(comps$csize)
    } else {
      metrics$AvgPathLength <- NA
      metrics$LargestComponentSize <- max(comps$csize)
    }
  }

  metrics$PositiveEdges <- sum(E(g)$weight > 0, na.rm = TRUE)
  metrics$NegativeEdges <- sum(E(g)$weight < 0, na.rm = TRUE)
  metrics$PosNegRatio <- ifelse(metrics$NegativeEdges > 0,
                                metrics$PositiveEdges / metrics$NegativeEdges,
                                Inf)

  deg_dist <- degree(g)
  metrics$MaxDegree <- max(deg_dist, na.rm = TRUE)
  metrics$MinDegree <- min(deg_dist, na.rm = TRUE)
  metrics$DegreeSD <- sd(deg_dist, na.rm = TRUE)

  metrics
}

identify_hubs <- function(g, group_name, outdir, top_percent = 0.1, min_degree = 5) {
  cat("\nIdentifying hub nodes for", group_name, "...\n")

  degree_cent <- degree(g)

  set.seed(123)
  betweenness_cent <- betweenness(g, weights = E(g)$abs_weight)

  g_dist <- g
  epsilon <- 1e-10
  E(g_dist)$distance <- 1 / (E(g_dist)$abs_weight + epsilon)
  set.seed(123)
  closeness_cent <- closeness(g_dist, weights = E(g_dist)$distance)

  set.seed(123)
  eigen_cent <- eigen_centrality(g, weights = E(g)$abs_weight)$vector

  set.seed(123)
  pagerank_cent <- page_rank(g, weights = E(g)$abs_weight)$vector

  centrality_df <- data.frame(
    Species = V(g)$name,
    AbundanceCategory = V(g)$abundance_category,
    Degree = degree_cent,
    Positive_Edges = V(g)$pos_edges,
    Negative_Edges = V(g)$neg_edges,
    Pos_Ratio = V(g)$pos_ratio,
    Betweenness = betweenness_cent,
    Closeness = closeness_cent,
    Eigenvector = eigen_cent,
    PageRank = pagerank_cent,
    stringsAsFactors = FALSE
  )

  normalize <- function(x) {
    if (length(x) == 0) return(x)
    if (max(x, na.rm = TRUE) - min(x, na.rm = TRUE) == 0) return(rep(0, length(x)))
    (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  }

  centrality_df$Degree_norm <- normalize(degree_cent)
  centrality_df$Betweenness_norm <- normalize(betweenness_cent)
  centrality_df$Closeness_norm <- normalize(closeness_cent)
  centrality_df$Eigenvector_norm <- normalize(eigen_cent)
  centrality_df$PageRank_norm <- normalize(pagerank_cent)

  centrality_df$CompositeScore <- rowMeans(
    centrality_df[, c("Degree_norm", "Betweenness_norm", "Closeness_norm",
                      "Eigenvector_norm", "PageRank_norm")],
    na.rm = TRUE
  )

  degree_threshold <- quantile(degree_cent, probs = 1 - top_percent, na.rm = TRUE)
  hubs_by_degree <- centrality_df[degree_cent >= degree_threshold & degree_cent >= min_degree, ]
  hubs_by_degree <- hubs_by_degree[order(hubs_by_degree$Degree, decreasing = TRUE), ]

  composite_threshold <- quantile(centrality_df$CompositeScore, probs = 1 - top_percent, na.rm = TRUE)
  hubs_by_composite <- centrality_df[centrality_df$CompositeScore >= composite_threshold, ]
  hubs_by_composite <- hubs_by_composite[order(hubs_by_composite$CompositeScore, decreasing = TRUE), ]

  for (col in c("Degree_norm", "Betweenness_norm", "Closeness_norm", "Eigenvector_norm", "PageRank_norm")) {
    threshold <- quantile(centrality_df[[col]], probs = 1 - top_percent, na.rm = TRUE)
    centrality_df[[paste0("IsTop_", col)]] <- ifelse(!is.na(centrality_df[[col]]),
                                                     centrality_df[[col]] >= threshold,
                                                     FALSE)
  }

  top_cols <- grep("IsTop_", names(centrality_df), value = TRUE)
  centrality_df$TopCount <- rowSums(centrality_df[, top_cols, drop = FALSE], na.rm = TRUE)
  hubs_by_consensus <- centrality_df[centrality_df$TopCount >= 3 & centrality_df$Degree >= min_degree, ]
  hubs_by_consensus <- hubs_by_consensus[order(hubs_by_consensus$TopCount,
                                               hubs_by_consensus$CompositeScore,
                                               decreasing = TRUE), ]

  write.csv(centrality_df, file = paste0(outdir, group_name, "_centrality_metrics.csv"),
            row.names = FALSE)

  cat("Hub nodes identified:\n")
  cat("  - By degree:", nrow(hubs_by_degree), "nodes\n")
  cat("  - By composite:", nrow(hubs_by_composite), "nodes\n")
  cat("  - By consensus:", nrow(hubs_by_consensus), "nodes\n")

  list(
    degree_hubs = if (nrow(hubs_by_degree) > 0) hubs_by_degree$Species else character(0),
    composite_hubs = if (nrow(hubs_by_composite) > 0) hubs_by_composite$Species else character(0),
    consensus_hubs = if (nrow(hubs_by_consensus) > 0) hubs_by_consensus$Species else character(0),
    all_metrics = centrality_df
  )
}

flashweave_build <- function(data_path, outdir, group_name, n_obs_min = NULL) {
  FlashWeave_other_para <- "sensitive=true"
  if (!is.null(n_obs_min)) {
    FlashWeave_other_para <- paste0(FlashWeave_other_para, ",n_obs_min=", n_obs_min)
  }
  L1 <- "using FlashWeave\n"
  L2 <- paste0('data_path = "', data_path, '"\n')
  L4 <- paste0(gsub(",$|,\\s+$", "",
                    paste0("netw_results = learn_network(data_path, ", FlashWeave_other_para)),
                    ")\n")
  L5 <- paste0('save_network("', outdir, '/network_FlashWeave_', group_name, '.gml", netw_results)')
  L <- paste0(L1, L2, L4, L5)
  openfile <- file(paste0(outdir, "/calculate_network.jl"), "wb")
  write(L, file = openfile)
  close(openfile)
  system(paste0("julia ", outdir, "/calculate_network.jl"))
  network <- read_graph(paste0(outdir, "network_FlashWeave_", group_name, ".gml"), format = "gml")
  set_vertex_attr(network, "name", value = V(network)$label)
}

build_weighted_network <- function(mat, group_name, taxa, sig_df, outdir,
                                   abundance_category, mat_file, cfg, outflag = "") {
  cat("\nBuilding full network:", group_name, "\n")
  cat("Number of taxa:", length(taxa), "\n")

  mat <- mat[, colSums(mat != 0) > 0, drop = FALSE]

  gml_files <- list.files(path = outdir, pattern = "\\.gml$", full.names = TRUE)
  gml_files <- gml_files[grepl(group_name, basename(gml_files))]

  net_files <- file.path(outdir, paste0("g_", group_name, "_full_network.rds"))
  if (length(gml_files) > 0) {
    g <- read_graph(gml_files[1], format = "gml")
    g <- set_vertex_attr(g, "name", value = V(g)$label)
    if (is.null(g$Stability)) {
      g$Stability <- 1
    }
    g$Samples_number <- nrow(mat)
    cat("GML file loaded:", gml_files[1], "\n")
  } else if (file.exists(net_files)) {
    g <- readRDS(net_files)
  } else {
    set.seed(cfg$runtime$seed)
    if (cfg$network$tool == "spiec-easi") {
      se_cfg <- cfg$network$spiec_easi
      if (se_cfg$relative_abundance) {
        se <- spiec.easi(
          mat,
          method = se_cfg$method,
          sel.criterion = se_cfg$sel.criterion,
          lambda.min.ratio = se_cfg$lambda.min.ratio,
          nlambda = se_cfg$nlambda,
          pulsar.select = se_cfg$pulsar_select,
          pulsar.params = se_cfg$pulsar,
          norm.params = se_cfg$norm.params
        )
      } else {
        se <- spiec.easi(
          mat,
          method = se_cfg$method,
          sel.criterion = se_cfg$sel.criterion,
          lambda.min.ratio = se_cfg$lambda.min.ratio,
          nlambda = se_cfg$nlambda,
          pulsar.select = se_cfg$pulsar_select,
          pulsar.params = se_cfg$pulsar
        )
      }

      stab <- getStability(se)
      cat(group_name, "network stability: ", stab, "\n")
      beta_sym <- as.matrix(symBeta(getOptBeta(se), mode = "maxabs"))
      refit <- getRefit(se)
      weighted_mat <- beta_sym * refit

      if (cfg$network$filter_by_weight) {
        weighted_mat[abs(weighted_mat) < 0.1] <- 0
      } else if (cfg$network$filter_by_confidence) {
        weighted_mat[getOptMerge(se) < 0.8] <- 0
      }

      g <- graph_from_adjacency_matrix(weighted_mat, mode = "undirected",
                                       weighted = TRUE, diag = FALSE)
      g$Stability <- stab
      g$Samples_number <- nrow(mat)
      if (length(V(g)) == length(taxa)) {
        V(g)$name <- taxa
      } else {
        V(g)$name <- colnames(mat)
      }

      g <- induced_subgraph(g, V(g)[degree(g) > 0])
      if (ecount(g) == 0) {
        cat("No edges detected in", group_name, "\n")
        return(list(graph = g, hubs = NULL))
      }
    } else if (cfg$network$tool == "flashweave") {
      cat("nrow", nrow(mat), "\n")
      if (nrow(mat) < 20) {
        n_obs_min <- nrow(mat) - 5
        g <- flashweave_build(mat_file, outdir, group_name, n_obs_min = n_obs_min)
      } else {
        g <- flashweave_build(mat_file, outdir, group_name)
      }
      g$Samples_number <- nrow(mat)
      g$Stability <- 1
    } else if (cfg$network$tool == "sparcc") {
      sparcc.amgut <- sparcc(mat)
      weight_mat <- sparcc.amgut$Cor
      weight_mat[abs(weight_mat) < 0.3] <- 0
      diag(weight_mat) <- 0
      g <- graph_from_adjacency_matrix(weight_mat, mode = "undirected",
                                       weighted = TRUE, diag = FALSE)
      g$Stability <- 1
      if (length(V(g)) == length(taxa)) {
        V(g)$name <- taxa
      } else {
        V(g)$name <- colnames(mat)
      }
      g <- induced_subgraph(g, V(g)[degree(g) > 0])
      if (ecount(g) == 0) {
        cat("No edges detected in", group_name, "\n")
        return(list(graph = g, hubs = NULL))
      }
    } else if (!is.null(cfg$network$graph_file)) {
      g <- read_graph(cfg$network$graph_file, format = "gml")
    } else {
      cat("No graph file")
      quit(1)
    }
  }

  abundance_category <- sapply(V(g)$name, function(species) {
    if (species %in% names(abundance_category)) {
      return(abundance_category[species])
    }
    "other"
  })

  V(g)$abundance_category <- abundance_category
  V(g)$abundance_color <- unlist(cfg$abundance$colors)[abundance_category]

  cat("Abundance category distribution:\n")
  cat("  High:", sum(abundance_category == "high"), "\n")
  cat("  Low:", sum(abundance_category == "low"), "\n")
  cat("  Very low:", sum(abundance_category == "very_low"), "\n")
  cat("  Other:", sum(abundance_category == "other"), "\n")

  V(g)$degree <- degree(g)

  V(g)$pos_edges <- 0
  V(g)$neg_edges <- 0
  V(g)$pos_ratio <- 0

  for (v in V(g)) {
    incident_edges <- incident(g, v)
    pos_count <- sum(E(g)$weight[incident_edges] > 0, na.rm = TRUE)
    neg_count <- sum(E(g)$weight[incident_edges] < 0, na.rm = TRUE)
    V(g)[v]$pos_edges <- pos_count
    V(g)[v]$neg_edges <- neg_count
    total_edges <- pos_count + neg_count
    V(g)[v]$pos_ratio <- ifelse(total_edges > 0, pos_count / total_edges, 0)
  }

  E(g)$abs_weight <- abs(E(g)$weight)

  hub_results <- identify_hubs(g, group_name, outdir,
                               top_percent = cfg$hub$top_percent,
                               min_degree = cfg$hub$min_degree)
  hub_diff_analysis <- analyze_hub_differential_overlap(
    hub_results, group_name, sig_df, outdir, outflag
  )

  V(g)$is_hub <- FALSE
  if (length(hub_results$consensus_hubs) > 0) {
    V(g)$is_hub <- V(g)$name %in% hub_results$consensus_hubs
  }

  V(g)$hub_color <- V(g)$abundance_color
  V(g)$hub_frame_color <- ifelse(V(g)$is_hub, "darkgoldenrod1", "gray50")
  V(g)$hub_frame_width <- ifelse(V(g)$is_hub, 2, 1)

  V(g)$size <- 3 + 2 * V(g)$degree / max(V(g)$degree)
  if (any(V(g)$is_hub)) {
    V(g)$size[V(g)$is_hub] <- V(g)$size[V(g)$is_hub] + 3
  }

  E(g)$color <- ifelse(E(g)$weight > 0, "#f783ac", "#3bc9db")
  E(g)$width <- 0.5 + 3 * E(g)$abs_weight / max(E(g)$abs_weight)

  sig_found <- sig_df[sig_df$species_taxid %in% V(g)$name, ]
  V(g)$label <- NA
  V(g)$label_color <- "black"
  V(g)$is_differential <- FALSE

  if (nrow(sig_found) > 0) {
    for (i in 1:nrow(sig_found)) {
      node_name <- as.character(sig_found$species_taxid[i])
      if (node_name %in% V(g)$name) {
        vertex_id <- which(V(g)$name == node_name)
        V(g)$label[vertex_id] <- as.character(sig_found$node_label[i])
        V(g)$label_color[vertex_id] <- as.character(sig_found$color[i])
        V(g)$is_differential[vertex_id] <- TRUE
      }
    }
  }

  list(graph = g, hubs = hub_results, hub_diff_analysis = hub_diff_analysis)
}
