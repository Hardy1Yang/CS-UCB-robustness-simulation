# =============================================================================
# Rcpp Loader - Compile and Load C++ Functions
# =============================================================================
# Description: Compiles C++ source files in src/ using Rcpp::sourceCpp().
#              Falls back gracefully to pure R if Rcpp is not available.
# =============================================================================

if (requireNamespace("Rcpp", quietly = TRUE)) {
  cpp_files <- list.files("src", pattern = "\\.cpp$", full.names = TRUE)
  if (length(cpp_files) > 0) {
    for (f in cpp_files) {
      # §B-3a (05_catR_package.md Action 1): fail-fast on compile error.
      # Silent warning() previously allowed Robust pipeline to proceed with
      # catR fallback while still labeling results as MPWI/KL/Progressive —
      # the documented root cause of Table 7 RMSE 1.6/2.4 anomalies.
      tryCatch({
        Rcpp::sourceCpp(f)
        message(paste0("  Compiled: ", f))
      }, error = function(e) {
        stop(sprintf(
          "[rcpp_loader] FATAL: failed to compile %s - %s\n  Aborting to prevent silent catR fallback which mislabels MFI/catR-MPWI results as C++ MPWI.",
          f, e$message))
      })
    }
  } else {
    stop("[rcpp_loader] FATAL: no .cpp files found in src/. Expected mle_fast.cpp, sh_calibration.cpp, traditional_selection.cpp, ucb_selection.cpp, irt_core.cpp.")
  }
} else {
  stop("[rcpp_loader] FATAL: Rcpp package not available. Robust pipeline depends on C++ backends — install Rcpp before running.")
}

# --- Override pure R functions with C++ versions ---

# Override find_mle_censored with C++ version (called directly, not via exists() check)
if (exists("find_mle_censored_cpp") && exists("find_mle_censored")) {
  find_mle_censored_r <- find_mle_censored  # Keep R version for testing
  find_mle_censored <- function(a_hist, b_hist, y_hist, theta_init = 0,
                                 theta_min = -4, theta_max = 4,
                                 tol = 1e-6, max_iter = 50) {
    find_mle_censored_cpp(a_hist, b_hist, as.integer(y_hist),
                           theta_init, theta_min, theta_max, tol, max_iter)
  }
  message("  Override: find_mle_censored -> C++ version")
}

# Report traditional selection C++ availability
if (exists("mpwi_select_cpp") && exists("kl_select_cpp") &&
    exists("progressive_select_cpp")) {
  message("  C++ available: mpwi_select_cpp, kl_select_cpp, progressive_select_cpp")
} else {
  message("  Traditional selection C++ not available - will use catR fallback")
}

message("=== Rcpp loader complete ===")
