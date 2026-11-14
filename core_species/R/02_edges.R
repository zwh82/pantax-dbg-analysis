# Node annotation and edge tables.
#
# A species is core only when it belongs to the global core list AND its local
# abundance category is `high` (abundant) or `low`/`very_low` (rare).  The
# global list alone is not enough: a globally core species that is locally
# `other` is treated as non-core.

make_node_meta <- function(g, core_rare, core_abundant) {
  species <- as.character(igraph::V(g)$name)
  local_raw <- as.character(igraph::V(g)$abundance_category)
  if (length(species) != length(local_raw)) stop("Invalid abundance_category attribute.")
  dplyr::tibble(
    species = species,
    local_abundance_raw = local_raw,
    abundance_class = abundance_class(local_raw),
    is_abundance_eligible = is_abundance_eligible(abundance_class),
    core_type = dplyr::case_when(
      species %in% core_rare & is_abundance_eligible ~ "core_rare",
      species %in% core_abundant & is_abundance_eligible ~ "core_abundant",
      abundance_class == "rare" ~ "noncore_rare",
      abundance_class == "abundant" ~ "noncore_abundant",
      TRUE ~ "noncore_other"
    ),
    is_core = core_type %in% c("core_rare", "core_abundant")
  )
}

sort_pair <- function(a, b) {
  ifelse(a <= b, paste(a, b, sep = "__"), paste(b, a, sep = "__"))
}

make_edge_df <- function(g, nodes) {
  e <- igraph::as_edgelist(g, names = TRUE)
  if (!nrow(e)) return(dplyr::tibble())
  a <- as.character(e[, 1])
  b <- as.character(e[, 2])
  source <- pmin(a, b)
  target <- pmax(a, b)
  sm <- nodes[match(source, nodes$species), ]
  tm <- nodes[match(target, nodes$species), ]
  weight <- igraph::E(g)$weight
  if (is.null(weight)) weight <- rep(NA_real_, nrow(e))
  dplyr::tibble(
    source = source,
    target = target,
    edge_id = paste(source, target, sep = "__"),
    weight = as.numeric(weight),
    sign = dplyr::case_when(weight > 0 ~ "positive", weight < 0 ~ "negative", TRUE ~ "zero"),
    source_abundance = sm$abundance_class,
    target_abundance = tm$abundance_class,
    source_abundance_eligible = sm$is_abundance_eligible,
    target_abundance_eligible = tm$is_abundance_eligible,
    source_abundance_raw = sm$local_abundance_raw,
    target_abundance_raw = tm$local_abundance_raw,
    source_core_type = sm$core_type,
    target_core_type = tm$core_type,
    source_is_core = sm$is_core,
    target_is_core = tm$is_core
  ) %>%
    dplyr::mutate(
      edge_abundance_type = dplyr::case_when(
        source_abundance == "abundant" & target_abundance == "abundant" ~ "abundant-abundant",
        sort_pair(source_abundance, target_abundance) == "abundant__rare" ~ "abundant-rare",
        source_abundance == "rare" & target_abundance == "rare" ~ "rare-rare",
        TRUE ~ "other"
      ),
      edge_abundance_eligible = source_abundance_eligible & target_abundance_eligible,
      edge_core_type = sort_pair(source_core_type, target_core_type),
      is_core_associated = source_is_core | target_is_core,
      is_core_rare_associated =
        source_core_type == "core_rare" | target_core_type == "core_rare",
      is_core_abundant_associated =
        source_core_type == "core_abundant" | target_core_type == "core_abundant",
      is_core_rare_abundant =
        is_core_rare_associated &
          (source_abundance == "abundant" | target_abundance == "abundant"),
      is_core_rare_rare =
        is_core_rare_associated & edge_abundance_type == "rare-rare",
      is_core_abundant_abundant =
        is_core_abundant_associated & edge_abundance_type == "abundant-abundant"
    ) %>%
    dplyr::distinct(edge_id, .keep_all = TRUE)
}

# The 15 edge categories reported throughout the module.  Every category is
# restricted to edges whose two endpoints have an abundant/rare class.
edge_masks <- function(e) {
  eligible <- e$edge_abundance_eligible
  list(
    all_edges = eligible,
    core_associated = eligible & e$is_core_associated,
    noncore_noncore = eligible & !e$is_core_associated,
    core_rare_associated = eligible & e$is_core_rare_associated,
    noncore_rare_associated =
      eligible &
        !e$is_core_associated &
        (e$source_core_type == "noncore_rare" | e$target_core_type == "noncore_rare"),
    core_abundant_associated = eligible & e$is_core_abundant_associated,
    noncore_abundant_associated =
      eligible &
        !e$is_core_associated &
        (e$source_core_type == "noncore_abundant" |
          e$target_core_type == "noncore_abundant"),
    core_rare_abundant = eligible & e$is_core_rare_abundant,
    other_abundant_rare =
      eligible & e$edge_abundance_type == "abundant-rare" & !e$is_core_rare_associated,
    core_rare_rare = eligible & e$is_core_rare_rare,
    noncore_rare_rare =
      eligible & e$edge_abundance_type == "rare-rare" & !e$is_core_associated,
    core_abundant_abundant = eligible & e$is_core_abundant_abundant,
    abundant_rare = eligible & e$edge_abundance_type == "abundant-rare",
    rare_rare = eligible & e$edge_abundance_type == "rare-rare",
    abundant_abundant = eligible & e$edge_abundance_type == "abundant-abundant"
  )
}

edge_count_table <- function(e) {
  masks <- edge_masks(e)
  dplyr::tibble(
    metric = names(masks),
    observed = purrr::map_int(masks, ~ sum(.x, na.rm = TRUE))
  )
}
