# =============================================================================
# CAT Simulation - Package Configuration
# =============================================================================
# Description: Load and configure the packages used by the robustness study.
#              The original project also loaded mstR, mirtCAT, mirt, ltm, data.table,
#              kableExtra, here, yaml and jsonlite; the robustness code calls none of them.
# =============================================================================

# --- Package Installation Check ---
required_packages <- c(
  # CAT reference implementation (item-bank format and fallback selection)
  "catR",

  # Data manipulation and plotting
  "tidyverse",  # dplyr, ggplot2, tidyr, readr, forcats
  "scales",     # axis labels in scripts/05_figures.R

  # Parallel Computing
  "parallel",   # Parallel processing
  "foreach",    # Parallel loops
  "doParallel", # Parallel backend
  "doRNG",      # Reproducible parallel RNG

  # C++ Acceleration
  "Rcpp"        # R/C++ interface for hot-path functions
)

# Install missing packages
install_if_missing <- function(packages) {
  new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
  if (length(new_packages) > 0) {
    message("Installing missing packages: ", paste(new_packages, collapse = ", "))
    install.packages(new_packages, dependencies = TRUE)
  }
}

# Uncomment to auto-install missing packages
# install_if_missing(required_packages)

# --- Load Packages ---
load_packages <- function(packages) {
  for (pkg in packages) {
    if (requireNamespace(pkg, quietly = TRUE)) {
      library(pkg, character.only = TRUE)
      message(paste0("✓ Loaded: ", pkg))
    } else {
      stop(paste0("✗ Package not found: ", pkg, ". Install it with install.packages(\"", pkg, "\")."))
    }
  }
}

# Load all required packages
load_packages(required_packages)

# --- Global Options ---
options(
  scipen = 999,           # Disable scientific notation
  digits = 4,             # Default decimal places
  warn = 1,               # Show warnings immediately
  stringsAsFactors = FALSE
)

# --- Set Random Seed for Reproducibility ---
set.seed(20260113)  # Proposal defense date

message("\n=== Package configuration complete ===\n")
