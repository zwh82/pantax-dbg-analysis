# ERAS antibiotic recovery.

recovery_tables <- function(baseline, followups, orders) {
  categories <- edge_masks(baseline)
  recovery_matrix <- sapply(followups, function(f) {
    m <- match(baseline$edge_id, f$edge_id)
    present <- !is.na(m)
    strict <- rep(FALSE, nrow(baseline))
    strict[present] <-
      baseline$source_abundance[present] == f$source_abundance[m[present]] &
      baseline$target_abundance[present] == f$target_abundance[m[present]]
    strict
  })
  if (is.null(dim(recovery_matrix))) recovery_matrix <- matrix(recovery_matrix, ncol = 1)
  curve <- bind_rows_safe(purrr::imap(categories, function(mask, category) {
    bind_rows_safe(purrr::map(seq_along(followups), function(i) {
      dplyr::tibble(
        timepoint = names(followups)[i],
        time_order = orders[i],
        edge_category = category,
        n_baseline_edges = sum(mask),
        n_recovered_edges = sum(recovery_matrix[mask, i]),
        recovery_rate = safe_fraction(n_recovered_edges, n_baseline_edges)
      )
    }))
  }))
  first_index <- apply(recovery_matrix, 1, function(x) {
    z <- which(x)
    if (length(z)) z[1] else NA_integer_
  })
  first <- baseline %>%
    dplyr::transmute(
      edge_id, source, target, source_core_type, target_core_type,
      edge_abundance_type,
      first_recovery_timepoint = ifelse(
        is.na(first_index), NA_character_, names(followups)[first_index]
      ),
      first_recovery_order = ifelse(
        is.na(first_index), max(orders) + 1, orders[first_index]
      ),
      is_recovered = !is.na(first_index),
      censored = is.na(first_index)
    )
  list(curve = curve, first = first, matrix = recovery_matrix)
}
