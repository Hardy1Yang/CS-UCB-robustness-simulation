# =============================================================================
# 01_generate_item_banks.R - Item Bank Generation
# =============================================================================
# Description: Generates the 15 robustness banks (3 sizes x 5 sparsity levels)
#              and the 45 calibration-error banks built from them.
#              Each bank has its own seed, so the banks do not depend on the
#              order of generation. (The original script also generated the
#              precision-study banks first; that block is removed here.)
# Usage: Rscript scripts/01_generate_item_banks.R   (from the repository root)
# =============================================================================

if (!exists("robustness_settings")) source("scripts/00_setup.R")

print_section("Item Bank Generation")

# --- Configuration ---

# Bank sizes (from paper)
bank_sizes <- item_bank_settings$bank_sizes

# Sparsity levels for Robustness study
# Paper: 10%, 30%, 50%
# Committee: add 3%, 5%
# sparsity_levels already loaded from config/simulation_settings.R

# Error levels for Robustness Scenario 2
error_levels <- error_settings$error_levels

# Parameter distributions (from paper)
# a ~ lognormal(0, 0.3) for Robustness
# b ~ N(0, 1) (committee feedback: use Normal, not Uniform)

# Base seed for reproducibility
base_seed <- common_settings$seed_base

# --- Generate Robustness Study Item Banks ---

print_section("Robustness Study Item Banks")

robustness_banks <- list()

for (K in bank_sizes) {
 for (rho in sparsity_levels) {

   bank_name <- paste0("K", K, "_rho", rho * 100)
   cat("Generating", bank_name, "...\n")

   # Seed format from spec (simulation_settings_summary.tex L266): seed = 42 + K + 100ρ
   seed <- 42 + K + round(100 * rho)

   set.seed(seed)

   # Discrimination: lognormal(0, 0.3) per paper Robustness section
   a <- rlnorm(K, meanlog = item_bank_settings$robustness$a_meanlog, sdlog = item_bank_settings$robustness$a_sdlog)
   a <- pmin(pmax(a, item_bank_settings$robustness$a_bounds[1]), item_bank_settings$robustness$a_bounds[2])

   # Difficulty with sparsity
   n_sparse <- round(K * rho)
   n_moderate <- K - n_sparse

   # Moderate items: N(0, 0.8) truncated to [-2, 2]
   b_moderate <- rnorm(n_moderate, mean = 0, sd = item_bank_settings$robustness$b_moderate_sd)
   b_moderate <- pmin(pmax(b_moderate, item_bank_settings$robustness$b_moderate_range[1]), item_bank_settings$robustness$b_moderate_range[2])

   # Sparse items: N(0, 1.5) truncated to |b| > 2
   if (n_sparse > 0) {
     b_sparse <- numeric(n_sparse)
     idx <- 1
     while (idx <= n_sparse) {
       candidate <- rnorm(1, mean = 0, sd = item_bank_settings$robustness$b_sparse_sd)
       if (abs(candidate) > item_bank_settings$robustness$b_sparse_threshold) {
         b_sparse[idx] <- candidate
         idx <- idx + 1
       }
     }
   } else {
     b_sparse <- numeric(0)
   }

   # Combine and shuffle with same permutation (sparse label must track b)
   b <- c(b_moderate, b_sparse)
   is_sparse_vec <- c(rep(FALSE, n_moderate), rep(TRUE, n_sparse))
   perm <- sample(K)
   b <- b[perm]
   is_sparse_vec <- is_sparse_vec[perm]

   # Verify sparse label consistency
   stopifnot(all(is_sparse_vec == (abs(b) > item_bank_settings$robustness$b_sparse_threshold)))

   # Create bank
   bank <- data.frame(
     item_id = 1:K,
     a = a,
     b = b,
     is_sparse = is_sparse_vec
   )

   robustness_banks[[bank_name]] <- bank

   cat("  Sparse items:", n_sparse, "(", rho * 100, "%)\n")
   cat("  Difficulty range: [",
       round(min(bank$b), 2), ",",
       round(max(bank$b), 2), "]\n")
 }
}

# Save Robustness banks
saveRDS(robustness_banks,
       "output/item_banks/robustness_banks.rds")
cat("\nSaved: output/item_banks/robustness_banks.rds\n")

# --- Generate Error Condition Banks (Robustness Scenario 2) ---

print_section("Error Condition Banks")

error_banks <- list()

for (K in bank_sizes) {
 for (rho in sparsity_levels) {

   base_name <- paste0("K", K, "_rho", rho * 100)
   base_bank <- robustness_banks[[base_name]]

   for (sigma_err in error_levels) {

     bank_name <- paste0(base_name, "_err", sigma_err * 100)
     cat("Generating", bank_name, "...\n")

     # Seed format: 99 + K + 100ρ + 1000σ_err
     seed <- 99 + K + round(100 * rho) + round(1000 * sigma_err)

     # Generate noisy parameters
     noisy <- generate_noisy_parameters(
       a_true = base_bank$a,
       b_true = base_bank$b,
       sigma_err = sigma_err,
       seed = seed
     )

     # Create error bank
     bank <- data.frame(
       item_id = 1:K,
       a_true = base_bank$a,
       b_true = base_bank$b,
       a_noisy = noisy$a_noisy,
       b_noisy = noisy$b_noisy,
       error_a = noisy$error_a,
       error_b = noisy$error_b,
       is_sparse = base_bank$is_sparse
     )

     error_banks[[bank_name]] <- list(
       bank = bank,
       sigma_err = sigma_err,
       rmse_a = noisy$rmse_a,
       rmse_b = noisy$rmse_b
     )

     cat("  RMSE(a):", round(noisy$rmse_a, 4),
         " RMSE(b):", round(noisy$rmse_b, 4), "\n")
   }
 }
}

# Save Error banks
saveRDS(error_banks,
       "output/item_banks/error_banks.rds")
cat("\nSaved: output/item_banks/error_banks.rds\n")

# --- Summary Statistics ---

print_section("Item Bank Summary")

cat("Robustness Study Banks:\n")
cat("  Total configurations:", length(robustness_banks), "\n")
cat("  Bank sizes:", paste(bank_sizes, collapse = ", "), "\n")
cat("  Sparsity levels:", paste(sparsity_levels * 100, "%", sep = "", collapse = ", "), "\n")

cat("\nError Condition Banks:\n")
cat("  Total configurations:", length(error_banks), "\n")
cat("  Error levels (σ_err):", paste(error_levels, collapse = ", "), "\n")

# --- Validation ---

print_section("Validation Checks")

# Check sparsity
cat("Sparsity check (Robustness K=300, ρ=30%):\n")
bank <- robustness_banks$K300_rho30
cat("  Sparse items:", sum(bank$is_sparse), "\n")
cat("  Moderate items:", sum(!bank$is_sparse), "\n")
cat("  Actual sparsity:", round(mean(bank$is_sparse) * 100, 1), "%\n")

# --- Complete ---

print_section("Item Bank Generation Complete")
cat("Generated banks saved to: output/item_banks/\n")
