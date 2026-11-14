# Node roles and targeted-removal robustness.

node_role_table <- function(g, nodes) {
  aw <- abs(igraph::E(g)$weight)
  if (!length(aw) || all(is.na(aw))) aw <- NULL
  distance_w <- if (is.null(aw)) NULL else 1 / (aw + 1e-10)
  dplyr::tibble(
    species = nodes$species,
    core_type = nodes$core_type,
    abundance_class = nodes$abundance_class,
    degree = igraph::degree(g),
    strength = igraph::strength(g, weights = aw),
    betweenness = igraph::betweenness(g, weights = distance_w, normalized = TRUE),
    eigenvector = igraph::eigen_centrality(g, weights = aw)$vector,
    pagerank = igraph::page_rank(g, weights = aw)$vector,
    coreness = igraph::coreness(g)
  ) %>%
    dplyr::mutate(
      is_degree_hub = degree >= quantile(degree, 0.9, na.rm = TRUE),
      is_bridge = betweenness >= quantile(betweenness, 0.9, na.rm = TRUE)
    )
}

summarise_node_roles <- function(x) {
  x %>%
    dplyr::mutate(node_group = dplyr::case_when(
      core_type == "core_rare" ~ "core_rare",
      core_type == "core_abundant" ~ "core_abundant",
      core_type == "noncore_rare" ~ "noncore_rare",
      core_type == "noncore_abundant" ~ "noncore_abundant",
      TRUE ~ "other"
    )) %>%
    dplyr::group_by(node_group) %>%
    dplyr::summarise(
      n_nodes = dplyr::n(),
      median_degree = median(degree, na.rm = TRUE),
      median_strength = median(strength, na.rm = TRUE),
      median_betweenness = median(betweenness, na.rm = TRUE),
      median_eigenvector = median(eigenvector, na.rm = TRUE),
      median_pagerank = median(pagerank, na.rm = TRUE),
      median_coreness = median(coreness, na.rm = TRUE),
      hub_fraction = mean(is_degree_hub),
      bridge_fraction = mean(is_bridge),
      .groups = "drop"
    )
}

largest_component_fraction <- function(g) {
  if (igraph::vcount(g) == 0) return(0)
  max(igraph::components(g)$csize) / igraph::vcount(g)
}

attack_summary <- function(g, nodes, type, n_perm) {
  targets <- nodes$species[nodes$core_type == type & nodes$is_abundance_eligible]
  if (!length(targets) || length(targets) >= igraph::vcount(g)) return(dplyr::tibble())
  attacked <- igraph::delete_vertices(g, targets)
  strata <- table(nodes$abundance_class[nodes$core_type == type & nodes$is_abundance_eligible])
  random_lcc <- replicate(n_perm, {
    selected <- unlist(purrr::map2(names(strata), as.integer(strata), function(s, n) {
      pool <- nodes$species[
        nodes$is_abundance_eligible &
          nodes$abundance_class == s &
          !nodes$is_core
      ]
      if (length(pool) < n) {
        pool <- nodes$species[
          nodes$is_abundance_eligible &
            nodes$abundance_class == s
        ]
      }
      sample(pool, n)
    }), use.names = FALSE)
    largest_component_fraction(igraph::delete_vertices(g, unique(selected)))
  })
  dplyr::tibble(
    target_type = type,
    n_removed = length(targets),
    observed_largest_component_fraction = largest_component_fraction(attacked),
    random_mean = mean(random_lcc),
    random_sd = sd(random_lcc),
    empirical_p = (sum(random_lcc <= largest_component_fraction(attacked)) + 1) /
      (n_perm + 1),
    n_perm = n_perm
  )
}
