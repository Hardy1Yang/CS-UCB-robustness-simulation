# =============================================================================
# CAT Simulation Settings
# =============================================================================
# Description: Centralized simulation parameters for Precision and Robustness
#              studies. Based on simulation_settings_summary.tex (committee-
#              approved) and the minutes of the proposal defense of 2026-01-13.
#
# Single Source of Truth: All scripts reference these settings via
# config variable names (e.g., precision_settings$T_max).
# Test scripts (*_test.R, run_full_simulation_test.R) may override values.
# =============================================================================

# --- Common Parameters ---

common_settings <- list(
  irt_model = "2PL",
  theta_range = c(-3, 3),
  # 2026-06-17 — Precision aligned to Robustness at the committee's request:
  # the two studies should share the test length, the ability grid and the
  # replications per ability level, so that they can be compared.
  # theta_grid_precision used to be by=0.5 (13 pts); now matches robustness.
  theta_grid_precision = seq(-3, 3, by = 0.2),   # 31 pts (aligned to robustness)
  theta_grid_robustness = seq(-3, 3, by = 0.2),  # 31 pts
  n_theta_precision = 31,                        # was 13; aligned to robustness
  n_theta_robustness = 31,
  seed_base = 20260113,                       # Proposal defense date
  alpha = 0.05                                # 95% CI
)

# --- Item Bank Settings ---

item_bank_settings <- list(
  bank_sizes = c(300, 600, 900),

  # Precision Study: a ~ LogN(0, 0.25), b ~ N(0, 1)
  precision = list(
    a_meanlog = 0,
    a_sdlog = 0.25,
    a_bounds = c(0.5, 2.5),
    b_mean = 0,
    b_sd = 1
  ),

  # Robustness Study: a ~ LogN(0, 0.3), b with sparsity model
  robustness = list(
    a_meanlog = 0,
    a_sdlog = 0.3,
    a_bounds = c(0.5, 2.5),
    # Moderate: N(0, 0.8) truncated [-2, 2]
    b_moderate_sd = 0.8,
    b_moderate_range = c(-2, 2),
    # Sparse: N(0, 1.5) truncated |b| > 2
    b_sparse_sd = 1.5,
    b_sparse_threshold = 2
  )
)

# --- Sparsity Levels ---
# Per committee: add extreme 3%, 5% to original 10%, 30%, 50%
sparsity_levels <- c(0.03, 0.05, 0.10, 0.30, 0.50)

# --- Error Model (Robustness Scenario 2) ---
# eps_a ~ N(0, sigma_err), eps_b ~ N(0, 2*sigma_err)
# a_noisy = max(a + eps_a, 0.1)
error_settings <- list(
  error_levels = c(0.05, 0.10, 0.20),
  a_lower_bound = 0.1,
  b_error_multiplier = 2,                    # b error SD = 2 * sigma_err
  seed_formula = "99 + K + 100*rho + 1000*sigma_err"
)

# --- Precision Study Settings ---

precision_settings <- list(
  # 2026-06-17 — aligned to Robustness at the committee's request (test length,
  # theta grid and replications should match). T_max 60→50, n_reps 500→100; theta grid
  # picked up from common_settings$n_theta_precision (now 31).
  T_max = 50,                                  # was 60; aligned to Robustness
  n_reps = 100,                                # was 500; aligned to Robustness

  # Crossed design: 3 x 4 x 3 x 31 = 1116 conditions (was 468 at 13 pts)
  initializations = c("staircase_3", "staircase_5", "censored_mle"),
  ci_methods = c("asymptotic", "azuma", "bernstein", "mixture"),

  # Selection: InfoMax ONLY
  selection_method = "infomax",

  # No operational controls
  use_stratification = FALSE,
  use_exposure_control = FALSE,

  # Trace-based CI evaluation at every trial t = 1, ..., T_max
  trace_based = TRUE
)

# --- Robustness Study Settings ---

robustness_settings <- list(
  T_max = 50,
  n_reps = 100,

  # 2x2x2 Factorial Design
  factorial = list(
    initialization = c("staircase", "censored_mle"),
    stratified = c(TRUE, FALSE),              # BAS: b-blocking + a-stratification
    exposure_control = c(TRUE, FALSE)         # Sympson-Hetter, r_max = 0.25
  ),

  # Staircase: 5 items (Robustness uses staircase-5 only)
  n_staircase = 5,

  # BAS stratification
  n_a_strata = 3,

  # Exposure control
  max_exposure_rate = 0.25,

  # Selection methods (all 8 methods enabled)
  selection_methods = c(
    "infomax",            # MFI (baseline)
    "csucb_theory",       # CS-UCB Azuma + theory H_lower (novel)
    "csucb_obsinfo",      # CS-UCB Azuma + observed info (novel)
    "csucb_bernstein",    # CS-UCB Bernstein + observed info (Issue #2)
    "csucb_mixture",      # CS-UCB Mixture martingale + observed info (Issue #2)
    "mpwi",               # Maximum Posterior Weighted Info (C++ or catR backend)
    "kl_information",     # Kullback-Leibler (C++ or catR backend)
    "progressive"         # Progressive explore-exploit (C++ or catR backend)
  )
)

# --- Traditional Method Backend ---
# Controls which implementation is used for MPWI, KL, Progressive
# "cpp"  = C++ accelerated (src/traditional_selection.cpp)
# "catR" = catR R package (gold standard, slower)
# Can be overridden before sourcing 00_setup.R:
#   traditional_backend <- "catR"; source("scripts/00_setup.R")
if (!exists("traditional_backend")) {
  traditional_backend <- "cpp"
}

# --- BAS Stratification ---
# Per Chang, Qian, & Ying (2001):
# 1. B-blocking: sort items by b, divide into M blocks
# 2. A-stratification: within blocks, sort by a, assign to S strata
# 3. Recombination: combine s-th items from all blocks into stratum s
# Implementation: create_BAS_structure() in cat_functions.R

stratification_settings <- list(
  method = "BAS",                             # B-blocking + A-stratification
  n_a_strata = 3,
  n_b_blocks = 3
)

# --- Exposure Control ---
# Sympson-Hetter method with r_max = 0.25

exposure_settings <- list(
  method = "sympson_hetter",
  max_rate = 0.25
)

# --- Content Balancing ---
# Note: content balancing not used in pipeline (assign_content_categories
# and related functions moved to utils/standalone/). BAS stratification
# serves as the item selection constraint mechanism instead.

# --- Evaluation Metrics ---

metrics_settings <- list(
  # Ability estimation
  ability = c("bias", "RMSE", "MAE", "coverage", "ci_width"),

  # Item bank usage (per committee)
  usage = c("exposure_rate", "overlap_rate", "chi_square_balance", "gini"),

  # Fairness (per committee)
  fairness = c("group_bias", "group_RMSE"),

  # Ability group boundaries for fairness analysis
  ability_groups = c(-1, 0, 1)                # Low/Mid-low/Mid-high/High
)

# --- Output Settings ---

output_settings <- list(
  save_results = TRUE,
  output_dir = "processed_data",
  formats = c("rds"),
  verbose = TRUE
)

# --- MLE Settings ---

mle_settings <- list(
  theta_bounds = c(-4, 4),                    # Wider than theta_range for numerical stability
  tol = 1e-6,
  max_iter = 50
)

# --- Staircase Settings ---

staircase_settings <- list(
  lengths = c(3, 5),                          # Precision: both; Robustness: 5 only
  step_size = 0.5,                            # Step up/down on correct/incorrect
  b_tolerance = 0.5                           # |b - theta_init| <= tolerance for initial item
)

# --- Calibration Settings (Sympson-Hetter) ---

calibration_settings <- list(
  n_examinees_per_iter = 1000,                # Pilot sample size per iteration
  n_calibration_iter = 10,                    # Default calibration iterations
  n_calibration_retry = 20                    # Retry iterations if first attempt fails
)

# --- Analysis Settings ---

analysis_settings <- list(
  # Theta region boundaries (used in 04, 06, 08 for ability group classification)
  # |theta| <= 1: Center, 1 < |theta| <= 2: Moderate, |theta| > 2: Extreme
  theta_region_breaks = c(1, 2),

  # Scoring weights for combined metrics
  combined_score_weights = c(coverage = 0.5, width = 0.5),
  fairness_composite_n = 4                    # (rmse + bias + coverage + width) / 4
)

# --- Visualization Defaults ---

viz_settings <- list(
  dpi = 300,
  default_width = 10,
  default_height = 6
)

# --- Parallel Computing Settings ---

parallel_settings <- list(
  # 2026-06-18 — reduced from detectCores()-1 (=9) to 5 after observing
  # cluster worker deaths and memory pressure during K=900 sparse banks.
  # Empirically 5 workers fit in available RAM without swap thrashing.
  n_cores = min(5L, max(1L, parallel::detectCores() - 1L)),
  use_parallel = TRUE
)

message("=== Simulation settings loaded ===")
