# =============================================================================
# Traditional Method Wrappers - Switchable Backend (C++ / catR)
# =============================================================================
# Description: Unified interface for MPWI, KL, and Progressive item selection.
#              Automatically dispatches to C++ or catR based on
#              the `traditional_backend` setting in simulation_settings.R.
#
# Usage:
#   traditional_backend <- "cpp"   # Use C++ accelerated version
#   traditional_backend <- "catR"  # Use catR R package (gold standard)
#
# Dependencies:
#   - C++ backend: src/traditional_selection.cpp (compiled via rcpp_loader.R)
#   - catR backend: utils/package_wrappers.R (catR_select, to_catR_bank)
# =============================================================================

#' Select Item Using Traditional Method (Unified Backend)
#'
#' @description
#' Dispatches MPWI, KL, or Progressive item selection to either
#' C++ (fast) or catR (reference) implementation based on `traditional_backend`.
#'
#' @param method Selection method: "mpwi", "kl_information", "progressive"
#' @param a Full item bank discrimination parameters
#' @param b Full item bank difficulty parameters
#' @param available Vector of available item indices (1-indexed)
#' @param current_theta Current ability estimate
#' @param a_hist Discrimination parameters of administered items
#' @param b_hist Difficulty parameters of administered items
#' @param y_hist Response vector for administered items
#' @param n_admin Number of administered items
#' @param t Current trial number (for progressive, 1-indexed)
#' @param T_max Total test length (for progressive)
#' @param K Total items in bank (for catR fallback)
#' @param itemBank_catR Pre-computed catR bank matrix (for catR fallback)
#' @param backend "cpp" or "catR" (default: reads from traditional_backend)
#' @return Selected item index (1-indexed)
#' @export
traditional_select <- function(method, a, b, available, current_theta,
                                a_hist = NULL, b_hist = NULL,
                                y_hist = NULL, n_admin = 0L,
                                t = NULL, T_max = NULL,
                                K = NULL, itemBank_catR = NULL,
                                backend = NULL) {

  # Determine backend
  if (is.null(backend)) {
    backend <- if (exists("traditional_backend", envir = .GlobalEnv)) {
      get("traditional_backend", envir = .GlobalEnv)
    } else if (exists("traditional_backend", inherits = TRUE)) {
      traditional_backend
    } else {
      "cpp"  # Default to C++ if available
    }
  }

  # --- C++ Backend ---
  if (backend == "cpp") {
    if (method == "mpwi" && exists("mpwi_select_cpp")) {
      return(mpwi_select_cpp(a, b, as.integer(available),
                             a_hist, b_hist,
                             as.integer(y_hist),
                             as.integer(n_admin)))

    } else if (method == "kl_information" && exists("kl_select_cpp")) {
      return(kl_select_cpp(a, b, as.integer(available),
                           current_theta,
                           a_hist, b_hist,
                           as.integer(y_hist),
                           as.integer(n_admin)))

    } else if (method == "progressive" && exists("progressive_select_cpp")) {
      return(progressive_select_cpp(a, b, as.integer(available),
                                     current_theta,
                                     as.integer(t),
                                     as.integer(T_max)))
    }

    # §B-3b (05_catR_package.md Action 2): fail-fast instead of silent fallback.
    # Old behavior `message(...)` then proceed to catR caused Robust to label
    # catR-MPWI as C++-MPWI in results (root cause of Table 7 RMSE 2.4 anomaly).
    if (!exists("mpwi_select_cpp")) {
      stop(sprintf(
        "[traditional_select] FATAL: C++ backend missing for method '%s' (mpwi_select_cpp not found in global env). Rerun script/00_setup.R; if Rcpp compile failed, fix the .cpp file rather than silently falling back to catR (which would mislabel MFI/catR-MPWI as %s).",
        method, method))
    }
  }

  # --- catR Backend ---
  # Call catR::nextItem directly (not via catR_select) because
  # MPWI/KL need the response vector x and administered item params
  if (is.null(K)) K <- length(a)
  if (is.null(itemBank_catR)) itemBank_catR <- to_catR_bank(a, b)

  criterion <- switch(method,
    "mpwi" = "MPWI",
    "kl_information" = "KL",
    "progressive" = "progressive"
  )

  # Build 'out' from available (catR uses exclude list)
  out_items <- setdiff(1:K, available)
  if (length(out_items) == 0) out_items <- NULL

  result <- tryCatch({
    catR::nextItem(
      itemBank = itemBank_catR,
      model = NULL,
      theta = current_theta,
      out = out_items,
      x = if (!is.null(y_hist) && length(y_hist) > 0) y_hist else NULL,
      criterion = criterion,
      method = "BM",
      priorDist = "norm",
      priorPar = c(0, 1),
      D = 1,
      range = c(-4, 4),
      nAvailable = if (method == "progressive") length(available) else NULL,
      rule = "length",
      thr = if (!is.null(T_max)) T_max else 50
    )
  }, error = function(e) {
    # §B-3c (05_catR_package.md Action 2): fail-fast instead of silent MFI fallback.
    # Previous behavior masked MFI results as catR-MPWI/KL/Progressive, breaking
    # method attribution in Robust pipeline.
    stop(sprintf(
      "[traditional_select catR backend] FATAL: catR::nextItem failed for criterion '%s': %s\n  Pipeline aborted to prevent silently mislabeling MFI as %s.",
      criterion, e$message, criterion))
  })

  return(result$item)
}

#' Check Which Backend is Active
#'
#' @return List with backend name and C++ availability
#' @export
check_traditional_backend <- function() {
  backend <- if (exists("traditional_backend", envir = .GlobalEnv)) {
    get("traditional_backend", envir = .GlobalEnv)
  } else {
    "cpp"
  }

  cpp_available <- exists("mpwi_select_cpp") &&
                   exists("kl_select_cpp") &&
                   exists("progressive_select_cpp")

  list(
    backend = backend,
    cpp_available = cpp_available,
    effective = if (backend == "cpp" && cpp_available) "cpp" else "catR"
  )
}

message("=== Traditional method wrappers loaded ===")
