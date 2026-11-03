# =============================================================================
# CAT Simulation - Custom Functions
# =============================================================================
# Description: Custom utility functions for CAT simulation
# =============================================================================

# --- A-Stratified B-Blocking (BAS) ---

#' Create A-Stratified B-Blocking structure
#' Following Chang, Qian, & Ying (2001), the BAS method uses a two-stage
#' stratification: (1) b-blocking first, (2) a-stratification within each
#' b-block, (3) recombine kth levels across all blocks into final strata.
#'
#' @param item_bank Data frame with item parameters (a, b columns)
#' @param n_a_strata Number of a-parameter strata (K levels)
#' @param n_b_blocks Number of b-parameter blocks (M blocks)
#' @return Item bank with stratum assignments (a_stratum column)
create_BAS_structure <- function(item_bank, n_a_strata = 3, n_b_blocks = NULL) {

  n_items <- nrow(item_bank)

  # Default: M = n_items / n_a_strata so each b-block has exactly n_a_strata items
  if (is.null(n_b_blocks)) {
    n_b_blocks <- floor(n_items / n_a_strata)
  }

  # Step 1: B-blocking — divide items into M blocks by b (ascending)
  item_bank <- item_bank[order(item_bank$b), ]
  items_per_block <- n_a_strata  # each block has K items (one per stratum)
  n_b_blocks_actual <- floor(n_items / items_per_block)
  n_used <- n_b_blocks_actual * items_per_block

  # Assign b-block membership
  item_bank$b_block <- NA
  item_bank$b_block[1:n_used] <- rep(1:n_b_blocks_actual, each = items_per_block)

  # Handle leftover items (if n_items not divisible by n_a_strata)
  if (n_used < n_items) {
    item_bank$b_block[(n_used + 1):n_items] <- n_b_blocks_actual
  }

  # Step 2: A-stratification within each b-block
  item_bank$a_stratum <- NA
  for (m in 1:n_b_blocks_actual) {
    block_idx <- which(item_bank$b_block == m)
    # Sort by a within this b-block, assign strata 1..K
    a_order <- order(item_bank$a[block_idx])
    item_bank$a_stratum[block_idx[a_order]] <- 1:length(block_idx)
  }

  # Handle leftover items — assign to closest stratum by a rank
  if (any(is.na(item_bank$a_stratum))) {
    leftover_idx <- which(is.na(item_bank$a_stratum))
    leftover_sorted <- leftover_idx[order(item_bank$a[leftover_idx])]
    item_bank$a_stratum[leftover_sorted] <- rep_len(1:n_a_strata, length(leftover_sorted))
  }

  # Step 3: Recombination is implicit — a_stratum now defines the final strata
  # Each stratum s contains items from all b-blocks at a-level s,
  # so each stratum covers the full range of b values.

  # Restore original order
  item_bank <- item_bank[order(item_bank$item_id), ]

  return(item_bank)
}

message("=== Custom CAT functions loaded ===")
