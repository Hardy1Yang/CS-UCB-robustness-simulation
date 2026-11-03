# =============================================================================
# 00_setup.R - Environment Setup and Validation
# =============================================================================
# Description: Sets up the simulation environment, loads all required packages
#              and custom functions, validates dependencies.
# Usage: source("scripts/00_setup.R") from the repository root.
#        Every script in scripts/ sources this file first.
# =============================================================================

cat("=== Robustness Simulation Setup ===\n")
cat("Starting at:", format(Sys.time()), "\n\n")

# --- Set Working Directory ---
# Assumes script is run from the repository root

# Try to find project root
possible_roots <- c(
  ".",
  "..",
  getwd()
)

# Add RStudio path if available
if (requireNamespace("rstudioapi", quietly = TRUE)) {
  tryCatch({
    if (rstudioapi::isAvailable()) {
      possible_roots <- c(possible_roots,
                          dirname(rstudioapi::getSourceEditorContext()$path))
    }
  }, error = function(e) NULL)
}

project_root <- NULL
for (root in possible_roots) {
 if (file.exists(file.path(root, "README.md")) &&
     dir.exists(file.path(root, "utils"))) {
   project_root <- normalizePath(root)
   break
 }
}

if (is.null(project_root)) {
 stop("Cannot find project root. Please set working directory to project root.")
}

setwd(project_root)
cat("Project root:", project_root, "\n\n")

# --- Load Configuration ---
cat("Loading configuration...\n")

# Load packages
if (file.exists("config/packages.R")) {
 source("config/packages.R")
} else {
 stop("config/packages.R not found")
}

# Load simulation settings
if (file.exists("config/simulation_settings.R")) {
 source("config/simulation_settings.R")
} else {
 stop("config/simulation_settings.R not found")
}

# --- Load Custom Functions ---
cat("\nLoading custom functions...\n")

utils_files <- c(
 "utils/helpers.R",
 "utils/irt_functions.R",
 "utils/mle_functions.R",
 "utils/cat_functions.R",
 "utils/staircase_functions.R",
 "utils/package_wrappers.R",
 "utils/exposure_control.R",
 "utils/ci_functions.R",
 "utils/error_model.R",
 "utils/robustness_functions.R",
 "utils/rcpp_loader.R",
 "utils/traditional_wrappers.R"
)
# Note: standalone/reference implementations in utils/standalone/ (not loaded here)
# Note: rcpp_loader.R is last so all pure R functions are available as fallbacks

for (f in utils_files) {
 if (file.exists(f)) {
   source(f)
   cat("  Loaded:", f, "\n")
 } else {
   warning("File not found: ", f)
 }
}

# --- §B-4: Validate C++ backend loaded — fail loudly otherwise ---
# Prevents the documented Table 7 MPWI RMSE 2.4 / KL RMSE 1.6 anomaly that
# happens when Rcpp compile silently fails and KL/MPWI/Progressive fall back
# to catR while still being labeled by their original method name.
required_cpp <- c("mpwi_select_cpp", "kl_select_cpp",
                  "progressive_select_cpp", "find_mle_censored_cpp",
                  "sh_calibrate_cpp")
missing_cpp <- required_cpp[!vapply(required_cpp, exists, logical(1))]
if (length(missing_cpp) > 0) {
  stop(sprintf(
    "[00_setup] FATAL: missing required C++ functions: %s.\n  Rebuild src/*.cpp before running simulation (do not let pipeline proceed with catR fallback).",
    paste(missing_cpp, collapse = ", ")))
}
cat("  All required C++ backends loaded:\n")
cat("    -", paste(required_cpp, collapse = "\n    - "), "\n")

# --- Validate Environment ---
cat("\nValidating environment...\n")

# Check required functions exist
required_functions <- c(
 # IRT
 "sigmoid", "p_fun", "info_fun",
 # MLE
 "find_mle_censored",
 # CAT
 "create_BAS_structure",
 # Staircase
 "run_staircase",
 # catR Wrappers
 "to_catR_bank", "catR_select",
 # Exposure
 "create_exposure_manager", "sympson_hetter_select",
 # CI
 "ci_asymptotic", "ci_azuma_hoeffding", "ci_bernstein", "ci_mixture",
 # Error
 "generate_noisy_parameters",
 # Helpers
 "print_section", "print_progress", "save_with_timestamp", "load_most_recent",
 # Precision
 # Robustness
 "run_robustness_test_factorial", "run_robustness_test_error_factorial",
 # Traditional method wrappers
 "traditional_select"
)

missing_functions <- required_functions[!sapply(required_functions, exists)]

if (length(missing_functions) > 0) {
 warning("Missing functions: ", paste(missing_functions, collapse = ", "))
} else {
 cat("  All required functions loaded successfully\n")
}

# --- Create Output Directories ---
cat("\nCreating output directories...\n")

output_dirs <- c(
  "output",
  "output/item_banks",
  "output/robustness",
  "output/robustness/analysis",
  "output/robustness_cpp",
  "output/visualization/robustness",
  "output/tables",
  "output/figures"
)

for (d in output_dirs) {
 if (!dir.exists(d)) {
   dir.create(d, recursive = TRUE)
   cat("  Created:", d, "\n")
 }
}

# --- Display Configuration Summary ---
cat("\n=== Configuration Summary ===\n")
cat("Item bank sizes:", paste(item_bank_settings$bank_sizes, collapse = ", "), "\n")
cat("Theta range:", paste(common_settings$theta_range, collapse = " to "), "\n")
cat("Robustness replications:", robustness_settings$n_reps, "\n")
cat("Max exposure rate:", exposure_settings$max_rate, "\n")
cat("Staircase lengths:", paste(staircase_settings$lengths, collapse = ", "), "\n")

# --- Helper Functions ---
# print_section, print_progress, save_with_timestamp, load_most_recent
# are now in utils/helpers.R (loaded above)

# --- Setup Complete ---
cat("\n=== Setup Complete ===\n")
cat("Finished at:", format(Sys.time()), "\n")
