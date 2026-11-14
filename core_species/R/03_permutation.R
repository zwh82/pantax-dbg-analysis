# Stratified label permutation for core-edge enrichment.
#
# Core labels are permuted within the broad abundance class (abundant / rare)
# used by the edge analysis, so the observed and random networks keep the same
# number of core rare and core abundant species in each class.

sample_core_labels <- function(nodes) {
  out <- rep("noncore", nrow(nodes))
  chosen <- integer()
  for (type in c("core_rare", "core_abundant")) {
    real <- which(nodes$core_type == type & nodes$is_abundance_eligible)
    # Permute core labels within the same broad abundance class used by the
    # edge analysis: high = abundant, low/very_low = rare.
    # Other abundance classes do not participate in the permutation.
    by_stratum <- table(nodes$abundance_class[real])
    picked <- unlist(purrr::map2(names(by_stratum), as.integer(by_stratum), function(s, n) {
      pool <- which(
        nodes$is_abundance_eligible &
          nodes$abundance_class == s &
          !seq_len(nrow(nodes)) %in% chosen
      )
      if (length(pool) < n) stop("Permutation pool is too small in stratum: ", s)
      sample(pool, n)
    }), use.names = FALSE)
    out[picked] <- type
    chosen <- c(chosen, picked)
  }
  out
}

apply_random_core <- function(nodes, labels) {
  nodes$core_type <- dplyr::case_when(
    labels == "core_rare" ~ "core_rare",
    labels == "core_abundant" ~ "core_abundant",
    nodes$is_abundance_eligible & nodes$abundance_class == "rare" ~ "noncore_rare",
    nodes$is_abundance_eligible & nodes$abundance_class == "abundant" ~ "noncore_abundant",
    TRUE ~ "noncore_other"
  )
  nodes$is_core <- grepl("^core_", nodes$core_type)
  nodes
}

permutation_enrichment <- function(g, nodes, observed, n_perm) {
  wanted <- c(
    "core_associated", "core_rare_associated", "core_abundant_associated",
    "core_rare_abundant", "core_rare_rare", "core_abundant_abundant"
  )
  random <- replicate(n_perm, {
    rn <- apply_random_core(nodes, sample_core_labels(nodes))
    rc <- edge_count_table(make_edge_df(g, rn))
    setNames(rc$observed, rc$metric)[wanted]
  })
  obs <- setNames(observed$observed, observed$metric)[wanted]
  dplyr::tibble(
    metric = wanted,
    observed = as.numeric(obs),
    random_mean = rowMeans(random),
    random_sd = apply(random, 1, sd),
    z_score = ifelse(random_sd == 0, NA_real_, (observed - random_mean) / random_sd),
    enrichment_ratio = ifelse(random_mean == 0, NA_real_, observed / random_mean),
    empirical_p = (rowSums(random >= observed) + 1) / (n_perm + 1),
    n_perm = n_perm
  )
}
