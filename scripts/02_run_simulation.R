# =============================================================================
# 02_run_simulation.R - Robustness Study Main Simulation
# =============================================================================
# Description: Robustness Study examining performance under practical CAT
#              conditions with operational controls
#
# Committee guidance:
#   - Manipulate exposure control and a-stratification (operational conditions)
#   - Use a 2×2×2 factorial design
#
# Factorial Design (2×2×2):
#   1. Initialization: Staircase vs Censored MLE
#   2. a-Stratified: Yes vs No
#   3. Exposure Control: Yes vs No
#
# Additional Factors:
#   - Bank sizes: K ∈ {300, 600, 900}
#   - Sparsity levels: {3%, 5%, 10%, 30%, 50%} (per committee feedback)
#   - Error levels: σ_err ∈ {0.05, 0.10, 0.20}
#   - Selection methods: MFI, CS-UCB-Theory, CS-UCB-ObsInfo,
#                        CS-UCB-Bernstein, CS-UCB-Mixture, MPWI, KL
#
# Prerequisites: Run scripts/01_generate_item_banks.R first
#
# Usage (from the repository root, in a fresh R process):
#   Rscript scripts/02_run_simulation.R [--cores=N] [--quick]
#   --cores=N  number of parallel workers (default: parallel_settings$n_cores).
#              Results do not depend on N: every test draws from its own
#              doRNG stream.
#   --quick    runs two short checks instead of the full design (see README):
#              bank K300_rho3 with factorial cell 1 (perfect parameters), and
#              bank K300_rho3_err5 with factorial cell 5 (calibration error).
#              Results go to output/quick/.
#
# Reproducibility: the doRNG streams are derived from the session's RNG state,
# which config/packages.R seeds with set.seed(20260113) and which the
# Sympson-Hetter calibrations and earlier loops then advance. A rerun therefore
# reproduces the archived results only when it runs in a fresh R process and
# from scratch, in the order of the loops below. Results cached in
# output/robustness_cpp/intermediate/ from an interrupted run are reused, but
# the banks computed after the restart then get different random streams.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
quick_mode <- "--quick" %in% args
cores_arg <- grep("^--cores=", args, value = TRUE)

source("scripts/00_setup.R")

if (length(cores_arg) > 0) {
  parallel_settings$n_cores <- as.integer(sub("^--cores=", "", cores_arg[1]))
  if (is.na(parallel_settings$n_cores) || parallel_settings$n_cores < 1) {
    stop("--cores must be a positive integer")
  }
}

# Results of this run: the full design writes to output/robustness/,
# the quick check to output/quick/.
data_dir <- if (quick_mode) "output/quick" else "output/robustness"
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

print_section("Robustness Study Simulation (2×2×2 Factorial)")

# --- Configuration ---

# Load item banks
if (!file.exists("output/item_banks/robustness_banks.rds")) {
  stop("Please run scripts/01_generate_item_banks.R first")
}
robustness_banks <- readRDS("output/item_banks/robustness_banks.rds")
cat("Loaded robustness item banks:", length(robustness_banks), "configurations\n")

if (!file.exists("output/item_banks/error_banks.rds")) {
  stop("Please run scripts/01_generate_item_banks.R first")
}
error_banks <- readRDS("output/item_banks/error_banks.rds")
cat("Loaded error condition banks:", length(error_banks), "configurations\n")

# Simulation parameters
robustness_config <- list(
  # Test configuration (from config/simulation_settings.R)
  T_max = robustness_settings$T_max,
  theta_range = common_settings$theta_range,
  mle_bounds = mle_settings$theta_bounds,
  n_theta = common_settings$n_theta_robustness,

  # Replications
  n_reps = robustness_settings$n_reps,

  # 2×2×2 Factorial Design (Committee Recommendation)
  factorial_design = expand.grid(
    initialization = robustness_settings$factorial$initialization,
    stratified = robustness_settings$factorial$stratified,
    exposure_control = robustness_settings$factorial$exposure_control,
    stringsAsFactors = FALSE
  ),

  # Staircase options (per committee: 3-5 items)
  n_staircase = robustness_settings$n_staircase,

  # Stratification settings
  n_strata = robustness_settings$n_a_strata,

  # Exposure control settings (Sympson-Hetter)
  max_exposure_rate = exposure_settings$max_rate,

  # Sparsity levels (per committee: add extreme 3%, 5%)
  sparsity_levels = sparsity_levels,

  # Error levels for Parameter Error condition
  error_levels = error_settings$error_levels,

  # Bank sizes
  bank_sizes = item_bank_settings$bank_sizes,

  # Selection methods (expanded per committee: add MPWI, KL, Progressive)
  selection_methods = robustness_settings$selection_methods,

  # Output paths (separated by backend for comparison)
  output_paths = {
    # Determine suffix based on backend
    has_traditional <- any(c("mpwi", "kl_information", "progressive") %in%
                          robustness_settings$selection_methods)
    suffix <- if (has_traditional && exists("traditional_backend") &&
                  traditional_backend == "catR") "_catR"
              else if (has_traditional) "_cpp"
              else ""
    list(
      data = if (quick_mode) "output/quick" else paste0("output/robustness", suffix),
      viz_factorial = paste0("output/visualization/robustness", suffix, "/factorial"),
      viz_sparsity = paste0("output/visualization/robustness", suffix, "/sparsity"),
      viz_error = paste0("output/visualization/robustness", suffix, "/error")
    )
  },

  # Confidence level
  alpha = common_settings$alpha,

  # Random seed base
  seed_base = common_settings$seed_base
)

# Auto-reduce config for catR backend
if (exists("traditional_backend") && traditional_backend == "catR") {
  cat("\n*** catR backend detected: using reduced configuration ***\n")
  robustness_config$bank_sizes <- 300
  robustness_config$n_reps <- 50
  robustness_config$factorial_design <- data.frame(
    initialization = "censored_mle",
    stratified = FALSE,
    exposure_control = FALSE,
    stringsAsFactors = FALSE
  )
  cat("  Bank sizes: K=300 only\n")
  cat("  Replications: 50\n")
  cat("  Factorial: 1 condition (censored_mle, no strat, no exposure)\n\n")
}

# Quick check: one bank per condition. The factorial design keeps all 8 rows,
# because the seeds use the row index f_idx; the loops below skip the other rows.
quick_cells <- c(perfect = 1L, error = 5L)
if (quick_mode) {
  robustness_config$bank_sizes <- 300
  robustness_config$sparsity_levels <- 0.03
  robustness_config$error_levels <- 0.05
  cat("\n*** --quick: bank K300_rho3, factorial cell", quick_cells[["perfect"]],
      "(perfect parameters); bank K300_rho3_err5, factorial cell",
      quick_cells[["error"]], "(calibration error) ***\n\n")
}

# Generate theta grid
theta_grid <- seq(robustness_config$theta_range[1],
                  robustness_config$theta_range[2],
                  length.out = robustness_config$n_theta)
cat("Theta grid:", round(theta_grid, 2), "\n")

# Display factorial design
cat("\n=== 2×2×2 Factorial Design ===\n")
print(robustness_config$factorial_design)
cat("Total factorial conditions:", nrow(robustness_config$factorial_design), "\n")

# --- Functions: see utils/robustness_functions.R ---
# run_robustness_test_factorial, run_robustness_test_error_factorial
# are loaded via 00_setup.R -> utils/robustness_functions.R
# Exposure control uses calibrated Sympson-Hetter from utils/exposure_control.R

# --- Condition 1: Perfect Parameter Knowledge (2×2×2 Factorial) ---

print_section("Condition 1: Perfect Parameters (2×2×2 Factorial)")

# Setup parallel cluster (shared by Conditions 1 and 2)
use_parallel <- exists("parallel_settings") &&
  parallel_settings$use_parallel &&
  requireNamespace("doRNG", quietly = TRUE)

if (use_parallel) {
  n_workers <- parallel_settings$n_cores
  cat("Parallel mode: ON (", n_workers, "cores)\n")
  cl <- makeCluster(n_workers)
  registerDoParallel(cl)

  # Set working directory on workers so they can find src/ files
  wd <- getwd()
  clusterExport(cl, "wd", envir = environment())
  clusterEvalQ(cl, setwd(wd))

  # Compile Rcpp on workers (shared libraries must be loaded per-process;
  # Rcpp functions cannot be serialized via clusterExport).
  # §B-3 hardening: workers must hard-fail on compile error to prevent silent
  # catR fallback from polluting MPWI/KL/Progressive results downstream.
  worker_compile_status <- clusterEvalQ(cl, {
    if (!requireNamespace("Rcpp", quietly = TRUE)) {
      stop("[worker] FATAL: Rcpp not available on parallel worker.")
    }
    cpp_files <- list.files("src", pattern = "\\.cpp$", full.names = TRUE)
    failed <- character(0)
    for (f in cpp_files) {
      tryCatch(Rcpp::sourceCpp(f),
               error = function(e) failed[length(failed) + 1L] <<- paste(f, e$message, sep = ": "))
    }
    required <- c("mpwi_select_cpp", "kl_select_cpp",
                  "progressive_select_cpp", "find_mle_censored_cpp",
                  "sh_calibrate_cpp")
    missing <- required[!vapply(required, exists, logical(1))]
    if (length(failed) > 0 || length(missing) > 0) {
      stop(sprintf("[worker %d] FATAL: cpp compile failures: %s; missing fns: %s",
                   Sys.getpid(),
                   paste(failed, collapse = "; "),
                   paste(missing, collapse = ", ")))
    }
    TRUE
  })
  if (!all(vapply(worker_compile_status, isTRUE, logical(1)))) {
    stopCluster(cl)
    stop("[03_robustness_study] FATAL: not all workers loaded C++ backend cleanly. Aborting.")
  }

  # Export R functions and data to workers
  # Note: do NOT export Rcpp functions (compiled locally on workers)
  export_fns <- c(
    "run_robustness_test_factorial",
    "run_robustness_test_error_factorial",
    "ci_asymptotic", "observed_info_2pl",
    "bernstein_score_bound", "bernstein_radius", "lr_mixture_radius",
    "sigmoid", "p_fun", "info_fun", "cW",
    "run_staircase", "to_catR_bank", "catR_select",
    "create_exposure_manager", "sympson_hetter_select",
    "create_BAS_structure",
    "generate_noisy_parameters",
    "traditional_select", "traditional_backend"
  )
  if (exists("find_mle_censored_r")) {
    export_fns <- c(export_fns, "find_mle_censored_r")
  } else {
    export_fns <- c(export_fns, "find_mle_censored")
  }
  clusterExport(cl, export_fns, envir = environment())

  # Apply find_mle_censored override on workers
  clusterEvalQ(cl, {
    if (exists("find_mle_censored_cpp")) {
      find_mle_censored <- function(a_hist, b_hist, y_hist,
                                     theta_init = 0,
                                     theta_min = -4,
                                     theta_max = 4,
                                     tol = 1e-6,
                                     max_iter = 50) {
        find_mle_censored_cpp(
          a_hist, b_hist, as.integer(y_hist),
          theta_init, theta_min, theta_max, tol, max_iter)
      }
    } else if (exists("find_mle_censored_r")) {
      find_mle_censored <- find_mle_censored_r
    }
  })
} else {
  cat("Parallel mode: OFF (sequential)\n")
}

factorial_results <- list()
item_traces_perfect <- list()
calibration_results <- list()
result_idx <- 1
trace_idx <- 1
# Note: trace data (robustness_traces_perfect) is NO LONGER accumulated globally.
# Each bank's trace data is saved to its intermediate file and loaded lazily
# by analysis scripts (05b_robustness_trace_analysis.R).
# This prevents memory exhaustion from ~6.2M rows per bank.

# Create intermediate directory for resume capability
# Uses backend-specific output path from robustness_config
intermediate_dir <- file.path(robustness_config$output_paths$data, "intermediate")
# A quick check always starts from scratch, as the comparison requires.
if (quick_mode) unlink(intermediate_dir, recursive = TRUE)
dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)
cat("Output directory:", robustness_config$output_paths$data, "\n")
cat("Backend:", if (exists("traditional_backend")) traditional_backend else "default", "\n")

for (K in robustness_config$bank_sizes) {
  for (rho in robustness_config$sparsity_levels) {
    bank_name <- paste0("K", K, "_rho", rho * 100)

    # Resume: check if this bank already has cached results
    intermediate_file <- file.path(intermediate_dir,
                                    paste0("cond1_", bank_name, ".rds"))
    if (file.exists(intermediate_file)) {
      cat("\n=== Bank:", bank_name, "=== [CACHED - skipping]\n")
      cached <- readRDS(intermediate_file)
      for (r in cached$factorial_results) {
        factorial_results[[result_idx]] <- r
        result_idx <- result_idx + 1
      }
      for (tr in cached$item_traces) {
        item_traces_perfect[[trace_idx]] <- tr
        trace_idx <- trace_idx + 1
      }
      # Trace data NOT loaded into memory — stays in intermediate file
      calibration_results[[bank_name]] <- cached$calibration
      rm(cached); gc()
      next
    }

    if (!bank_name %in% names(robustness_banks)) {
      cat("Skipping:", bank_name, "(not found)\n")
      next
    }

    item_bank <- robustness_banks[[bank_name]]
    cat("\n=== Bank:", bank_name, "===\n")

    # Pre-calibrate Sympson-Hetter for this bank (once per bank)
    cat("  Calibrating Sympson-Hetter (target rate =",
        robustness_config$max_exposure_rate, ")...")
    calibration_result <- calibrate_sympson_hetter(
      item_bank = data.frame(a = item_bank$a, b = item_bank$b),
      target_rate = robustness_config$max_exposure_rate,
      n_calibration = 50,
      n_examinees_per_iter = 1000,
      T_max = robustness_config$T_max,
      verbose = FALSE
    )
    # Guard: retry if calibration failed
    if (!calibration_result$converged ||
        calibration_result$max_exposure >
          robustness_config$max_exposure_rate) {
      cat("  WARNING: first calibration did not converge.",
          "Retrying with n_calibration=100...\n")
      calibration_result <- calibrate_sympson_hetter(
        item_bank = data.frame(a = item_bank$a,
                               b = item_bank$b),
        target_rate = robustness_config$max_exposure_rate,
        n_calibration = 100,
        n_examinees_per_iter = 1000,
        T_max = robustness_config$T_max,
        verbose = FALSE
      )
      cat("  Retry result: max_exp:",
          round(calibration_result$max_exposure, 3),
          ", converged:", calibration_result$converged, "\n")
      # Warn but continue if calibration is close (sparse banks may not converge exactly)
      if (calibration_result$max_exposure >
            robustness_config$max_exposure_rate) {
        cat("  WARNING: S-H above target (",
            round(calibration_result$max_exposure, 3),
            " > ", robustness_config$max_exposure_rate,
            "). Using best calibration found.\n")
      }
    }
    calibration_results[[bank_name]] <- calibration_result
    cat(" done (max_exp:",
        round(calibration_result$max_exposure, 3),
        ", converged:", calibration_result$converged, ")\n")

    # Record start indices for intermediate saving
    bank_result_start <- result_idx
    bank_trace_start <- trace_idx

    # Bank-local trace data list (written to disk per bank, then freed)
    bank_trace_data <- list()
    bank_tdata_idx <- 1

    # Loop through factorial design
    for (f_idx in 1:nrow(robustness_config$factorial_design)) {
      if (quick_mode && f_idx != quick_cells[["perfect"]]) next
      factorial_row <- robustness_config$factorial_design[f_idx, ]

      # Create factorial condition label
      f_label <- paste0(
        substr(factorial_row$initialization, 1, 5),
        "_", ifelse(factorial_row$stratified, "Str", "NoS"),
        "_", ifelse(factorial_row$exposure_control, "Exp", "NoE")
      )

      cat("  Factorial:", f_label, "\n")

      # Create exposure manager with calibrated params
      if (factorial_row$exposure_control) {
        exposure_manager <- create_exposure_manager(
          K, max_rate = robustness_config$max_exposure_rate
        )
        exposure_manager$control_params <- calibration_result$control_params
      } else {
        exposure_manager <- NULL
      }

      for (method in robustness_config$selection_methods) {
        cat("    Method:", method, "...")

        seed <- robustness_config$seed_base + K + round(rho * 1000) + f_idx * 100

        n_theta_local <- robustness_config$n_theta
        n_reps_local <- robustness_config$n_reps
        grid <- expand.grid(theta_idx = 1:n_theta_local, rep = 1:n_reps_local)

        if (use_parallel) {
          # Export loop-specific variables to workers
          clusterExport(cl, c("grid", "theta_grid", "item_bank",
                              "method", "robustness_config",
                              "factorial_row", "exposure_manager"),
                        envir = environment())
          par_results <- foreach(
            i = 1:nrow(grid),
            .packages = c("catR"),
            .noexport = ls()
          ) %dorng% {
            theta_true <- theta_grid[grid$theta_idx[i]]
            tryCatch(
              run_robustness_test_factorial(
                theta_true, item_bank, method, robustness_config,
                factorial_row, exposure_manager
              ),
              error = function(e) NULL
            )
          }
        } else {
          set.seed(seed)
          par_results <- vector("list", nrow(grid))
          for (i in 1:nrow(grid)) {
            theta_true <- theta_grid[grid$theta_idx[i]]
            par_results[[i]] <- tryCatch(
              run_robustness_test_factorial(
                theta_true, item_bank, method, robustness_config,
                factorial_row, exposure_manager
              ),
              error = function(e) NULL
            )
          }
        }

        # Post-process parallel results
        results_list <- vector("list", nrow(grid))
        for (i in 1:nrow(grid)) {
          test_result <- par_results[[i]]
          theta_true <- theta_grid[grid$theta_idx[i]]
          rep_num <- grid$rep[i]

          if (!is.null(test_result)) {
            # Aggregate exposure counts
            if (factorial_row$exposure_control && !is.null(exposure_manager)) {
              exposure_manager$exposure_counts[test_result$used_items] <-
                exposure_manager$exposure_counts[test_result$used_items] + 1
              exposure_manager$n_examinees <- exposure_manager$n_examinees + 1
            }

            # Save item-level trace
            item_traces_perfect[[trace_idx]] <- list(
              K = K, sparsity = rho, factorial_label = f_label,
              selection_method = method,
              theta_true = theta_true, rep = rep_num,
              used_items = test_result$used_items
            )
            trace_idx <- trace_idx + 1

            # Build per-trial trace (long format) — stored in bank-local list
            if (!is.null(test_result$theta_path)) {
              T_act <- test_result$T_actual
              bank_trace_data[[bank_tdata_idx]] <- data.frame(
                K = K, sparsity = rho,
                factorial_label = f_label,
                selection_method = method,
                theta_true = theta_true, rep = rep_num,
                t = 1:T_act,
                theta_hat_t = test_result$theta_path[2:(T_act + 1)],
                inst_regret_t = test_result$inst_regret_path,
                cum_regret_t = cumsum(test_result$inst_regret_path),
                selection_time_t = test_result$selection_time_path,
                stringsAsFactors = FALSE
              )
              bank_tdata_idx <- bank_tdata_idx + 1
            }

            results_list[[i]] <- data.frame(
              condition = "perfect_params",
              K = K,
              sparsity = rho,
              initialization = factorial_row$initialization,
              stratified = factorial_row$stratified,
              exposure_control = factorial_row$exposure_control,
              factorial_label = f_label,
              selection_method = method,
              sigma_err = NA,
              theta_true = theta_true,
              rep = rep_num,
              theta_final = test_result$theta_final,
              bias = test_result$bias,
              abs_error = test_result$abs_error,
              ci_coverage = test_result$ci_coverage,
              ci_width = test_result$ci_width,
              n_sparse_used = test_result$n_sparse_used,
              sparse_ratio = test_result$sparse_ratio,
              T_actual = test_result$T_actual,
              cum_regret = test_result$cum_regret,
              mean_selection_time = test_result$mean_selection_time
            )
          }
        }

        condition_results <- do.call(rbind, Filter(Negate(is.null), results_list))
        factorial_results[[result_idx]] <- condition_results
        result_idx <- result_idx + 1
        rm(par_results, results_list)

        cat(" RMSE =", round(sqrt(mean(condition_results$bias^2)), 3),
            ", Coverage =", round(mean(condition_results$ci_coverage) * 100, 1),
            "%\n")
      }
    }

    # Save intermediate results for this bank (resume capability)
    saveRDS(list(
      factorial_results = if (result_idx > bank_result_start)
        factorial_results[bank_result_start:(result_idx - 1)] else list(),
      item_traces = if (trace_idx > bank_trace_start)
        item_traces_perfect[bank_trace_start:(trace_idx - 1)] else list(),
      trace_data = bank_trace_data,
      calibration = calibration_results[[bank_name]]
    ), intermediate_file)
    cat("  [SAVED] Intermediate results:", bank_name, "\n")

    # Free bank-local trace data and par_results immediately
    rm(bank_trace_data); gc()
  }
}

# Combine factorial results
factorial_df <- do.call(rbind, factorial_results)
saveRDS(factorial_df, file.path(data_dir, "factorial_perfect_params.rds"))
saveRDS(item_traces_perfect, file.path(data_dir, "item_traces_perfect.rds"))
saveRDS(calibration_results, file.path(data_dir, "calibration_params.rds"))

# Save calibration diagnostics CSV
calib_diag <- do.call(rbind, lapply(names(calibration_results),
  function(bn) {
    cr <- calibration_results[[bn]]
    data.frame(
      bank = bn, max_exposure = cr$max_exposure,
      n_overexposed = cr$n_overexposed,
      iterations = cr$iterations,
      converged = cr$converged,
      stringsAsFactors = FALSE)
  }))
write.csv(calib_diag,
  file.path(data_dir, "calibration_diagnostics.csv"),
  row.names = FALSE)
cat("Saved calibration diagnostics CSV\n")

# Note: Per-trial trace data is now stored per-bank in intermediate files.
# Use 05b_robustness_trace_analysis.R to load and analyze traces lazily.
cat("\nSaved factorial results:", nrow(factorial_df), "rows\n")
cat("Saved item traces:", length(item_traces_perfect), "records\n")
cat("Trace data saved per-bank in:", intermediate_dir, "\n")
cat("Saved calibration params:", length(calibration_results), "banks\n")

# --- Condition 2: Parameter Estimation Error (2×2×2 Factorial) ---

print_section("Condition 2: Parameter Error (2×2×2 Factorial)")

error_factorial_results <- list()
item_traces_error <- list()
result_idx <- 1
trace_idx <- 1
# Trace data saved per-bank (same pattern as Condition 1)

for (K in robustness_config$bank_sizes) {
  for (rho in robustness_config$sparsity_levels) {
    for (sigma_err in robustness_config$error_levels) {

      bank_name <- paste0("K", K, "_rho", rho * 100, "_err", sigma_err * 100)

      # Resume: check if this bank already has cached results
      intermediate_file <- file.path(intermediate_dir,
                                      paste0("cond2_", bank_name, ".rds"))
      if (file.exists(intermediate_file)) {
        cat("\n=== Bank:", bank_name, "=== [CACHED - skipping]\n")
        cached <- readRDS(intermediate_file)
        for (r in cached$factorial_results) {
          error_factorial_results[[result_idx]] <- r
          result_idx <- result_idx + 1
        }
        for (tr in cached$item_traces) {
          item_traces_error[[trace_idx]] <- tr
          trace_idx <- trace_idx + 1
        }
        # Trace data NOT loaded into memory — stays in intermediate file
        rm(cached); gc()
        next
      }

      if (!bank_name %in% names(error_banks)) {
        cat("Skipping:", bank_name, "(not found)\n")
        next
      }

      error_bank <- error_banks[[bank_name]]
      K_bank <- nrow(error_bank$bank)
      cat("\n=== Bank:", bank_name, "===\n")

      # Calibrate S-H using noisy params (what the CAT algorithm sees)
      cat("  Calibrating S-H...")
      calibration_result_err <- calibrate_sympson_hetter(
        item_bank = data.frame(a = error_bank$bank$a_noisy,
                               b = error_bank$bank$b_noisy),
        target_rate = robustness_config$max_exposure_rate,
        n_calibration = 50,
        n_examinees_per_iter = 1000,
        T_max = robustness_config$T_max,
        verbose = FALSE
      )
      cat(" done (max_exp:",
          round(calibration_result_err$max_exposure, 3),
          ", converged:", calibration_result_err$converged,
          ")\n")

      # Guard: retry if calibration failed
      if (!calibration_result_err$converged ||
          calibration_result_err$max_exposure >
            robustness_config$max_exposure_rate) {
        cat("  WARNING: first calibration did not converge.",
            "Retrying with n_calibration=100...\n")
        calibration_result_err <- calibrate_sympson_hetter(
          item_bank = data.frame(
            a = error_bank$bank$a_noisy,
            b = error_bank$bank$b_noisy),
          target_rate = robustness_config$max_exposure_rate,
          n_calibration = 100,
          n_examinees_per_iter = 1000,
          T_max = robustness_config$T_max,
          verbose = FALSE
        )
        cat("  Retry result: max_exp:",
            round(calibration_result_err$max_exposure, 3),
            ", converged:",
            calibration_result_err$converged, "\n")
        # Warn but continue if calibration is close
        if (calibration_result_err$max_exposure >
              robustness_config$max_exposure_rate) {
          cat("  WARNING: S-H above target (",
              round(calibration_result_err$max_exposure, 3),
              " > ", robustness_config$max_exposure_rate,
              "). Using best calibration found.\n")
        }
      }

      # Record start indices for intermediate saving
      bank_result_start <- result_idx
      bank_trace_start <- trace_idx

      # Bank-local trace data list (written to disk per bank, then freed)
      bank_trace_data <- list()
      bank_tdata_idx <- 1

      # Loop through factorial design
      for (f_idx in 1:nrow(robustness_config$factorial_design)) {
        if (quick_mode && f_idx != quick_cells[["error"]]) next
        factorial_row <- robustness_config$factorial_design[f_idx, ]

        f_label <- paste0(
          substr(factorial_row$initialization, 1, 5),
          "_", ifelse(factorial_row$stratified, "Str", "NoS"),
          "_", ifelse(factorial_row$exposure_control, "Exp", "NoE")
        )

        cat("  Factorial:", f_label, "\n")

        if (factorial_row$exposure_control) {
          exposure_manager <- create_exposure_manager(
            K_bank, max_rate = robustness_config$max_exposure_rate
          )
          exposure_manager$control_params <- calibration_result_err$control_params
        } else {
          exposure_manager <- NULL
        }

        for (method in robustness_config$selection_methods) {
          cat("    Method:", method, "...")

          seed <- robustness_config$seed_base + K + round(rho * 1000) +
                  round(sigma_err * 10000) + f_idx * 100

          n_theta_local <- robustness_config$n_theta
          n_reps_local <- robustness_config$n_reps
          grid <- expand.grid(theta_idx = 1:n_theta_local, rep = 1:n_reps_local)

          if (use_parallel) {
            # Export loop-specific variables to workers (incl. sigma_err for
            # per-rep noise redraw, Patch 3 — Design A)
            clusterExport(cl, c("grid", "theta_grid", "error_bank",
                                "method", "robustness_config",
                                "factorial_row", "exposure_manager",
                                "sigma_err", "seed"),
                          envir = environment())
            par_results <- foreach(
              i = 1:nrow(grid),
              .packages = c("catR"),
              .noexport = ls()
            ) %dorng% {
              theta_true <- theta_grid[grid$theta_idx[i]]
              tryCatch(
                run_robustness_test_error_factorial(
                  theta_true, error_bank, method, robustness_config,
                  factorial_row, exposure_manager,
                  noise_seed = seed + grid$rep[i] * 1009L,
                  sigma_err = sigma_err
                ),
                error = function(e) NULL
              )
            }
          } else {
            set.seed(seed)
            par_results <- vector("list", nrow(grid))
            for (i in 1:nrow(grid)) {
              theta_true <- theta_grid[grid$theta_idx[i]]
              par_results[[i]] <- tryCatch(
                run_robustness_test_error_factorial(
                  theta_true, error_bank, method, robustness_config,
                  factorial_row, exposure_manager,
                  noise_seed = seed + grid$rep[i] * 1009L,
                  sigma_err = sigma_err
                ),
                error = function(e) NULL
              )
            }
          }

          # Post-process parallel results
          results_list <- vector("list", nrow(grid))
          for (i in 1:nrow(grid)) {
            test_result <- par_results[[i]]
            theta_true <- theta_grid[grid$theta_idx[i]]
            rep_num <- grid$rep[i]

            if (!is.null(test_result)) {
              if (factorial_row$exposure_control && !is.null(exposure_manager)) {
                exposure_manager$exposure_counts[test_result$used_items] <-
                  exposure_manager$exposure_counts[test_result$used_items] + 1
                exposure_manager$n_examinees <- exposure_manager$n_examinees + 1
              }

              item_traces_error[[trace_idx]] <- list(
                K = K, sparsity = rho, sigma_err = sigma_err,
                factorial_label = f_label, selection_method = method,
                theta_true = theta_true, rep = rep_num,
                used_items = test_result$used_items
              )
              trace_idx <- trace_idx + 1

              if (!is.null(test_result$theta_path)) {
                T_act <- test_result$T_actual
                bank_trace_data[[bank_tdata_idx]] <- data.frame(
                  K = K, sparsity = rho,
                  sigma_err = sigma_err,
                  factorial_label = f_label,
                  selection_method = method,
                  theta_true = theta_true, rep = rep_num,
                  t = 1:T_act,
                  theta_hat_t = test_result$theta_path[2:(T_act + 1)],
                  inst_regret_t = test_result$inst_regret_path,
                  cum_regret_t = cumsum(test_result$inst_regret_path),
                  selection_time_t = test_result$selection_time_path,
                  stringsAsFactors = FALSE
                )
                bank_tdata_idx <- bank_tdata_idx + 1
              }

              results_list[[i]] <- data.frame(
                condition = "param_error",
                K = K,
                sparsity = rho,
                initialization = factorial_row$initialization,
                stratified = factorial_row$stratified,
                exposure_control = factorial_row$exposure_control,
                factorial_label = f_label,
                selection_method = method,
                sigma_err = sigma_err,
                theta_true = theta_true,
                rep = rep_num,
                theta_final_noisy = test_result$theta_final_noisy,
                theta_final_oracle = test_result$theta_final_oracle,
                bias_noisy = test_result$bias_noisy,
                bias_oracle = test_result$bias_oracle,
                ci_coverage_noisy = test_result$ci_coverage_noisy,
                ci_coverage_oracle = test_result$ci_coverage_oracle,
                ci_width_noisy = test_result$ci_width_noisy,
                ci_width_oracle = test_result$ci_width_oracle,
                T_actual = test_result$T_actual,
                cum_regret = test_result$cum_regret,
                mean_selection_time = test_result$mean_selection_time
              )
            }
          }

          condition_results <- do.call(rbind, Filter(Negate(is.null), results_list))
          error_factorial_results[[result_idx]] <- condition_results
          result_idx <- result_idx + 1
          rm(par_results, results_list)

          cat(" RMSE(noisy) =", round(sqrt(mean(condition_results$bias_noisy^2)), 3),
              ", RMSE(oracle) =", round(sqrt(mean(condition_results$bias_oracle^2)), 3),
              "\n")
        }
      }

      # Save intermediate results for this bank (resume capability)
      saveRDS(list(
        factorial_results = if (result_idx > bank_result_start)
          error_factorial_results[bank_result_start:(result_idx - 1)] else list(),
        item_traces = if (trace_idx > bank_trace_start)
          item_traces_error[bank_trace_start:(trace_idx - 1)] else list(),
        trace_data = bank_trace_data
      ), intermediate_file)
      cat("  [SAVED] Intermediate results:", bank_name, "\n")

      # Free bank-local trace data immediately
      rm(bank_trace_data); gc()
    }
  }
}

# Combine error factorial results
error_factorial_df <- do.call(rbind, error_factorial_results)
saveRDS(error_factorial_df, file.path(data_dir, "factorial_param_error.rds"))
saveRDS(item_traces_error, file.path(data_dir, "item_traces_error.rds"))

# Note: Per-trial trace data is stored per-bank in intermediate files.
cat("\nSaved error factorial results:", nrow(error_factorial_df), "rows\n")
cat("Saved error item traces:", length(item_traces_error), "records\n")
cat("Error trace data saved per-bank in:", intermediate_dir, "\n")

# --- Save Combined Results ---

print_section("Saving Combined Results")

robustness_results <- list(
  factorial_perfect = factorial_df,
  factorial_error = error_factorial_df,
  config = robustness_config
)
saveRDS(robustness_results, file.path(data_dir, "robustness_results_full.rds"))

# --- Summary Statistics ---

print_section("Robustness Study Summary")

# Factorial Effect Summary (Perfect Parameters)
cat("=== 2×2×2 Factorial Effects (Perfect Parameters) ===\n\n")

factorial_summary <- factorial_df %>%
  group_by(initialization, stratified, exposure_control, selection_method) %>%
  summarise(
    n = n(),
    bias = mean(bias),
    RMSE = sqrt(mean(bias^2)),
    coverage = mean(ci_coverage) * 100,
    .groups = "drop"
  )
print(factorial_summary)

# Main Effects
cat("\n--- Main Effect: Initialization ---\n")
init_effect <- factorial_df %>%
  group_by(initialization) %>%
  summarise(
    RMSE = sqrt(mean(bias^2)),
    coverage = mean(ci_coverage) * 100,
    .groups = "drop"
  )
print(init_effect)

cat("\n--- Main Effect: Stratification ---\n")
strat_effect <- factorial_df %>%
  group_by(stratified) %>%
  summarise(
    RMSE = sqrt(mean(bias^2)),
    coverage = mean(ci_coverage) * 100,
    .groups = "drop"
  )
print(strat_effect)

cat("\n--- Main Effect: Exposure Control ---\n")
exp_effect <- factorial_df %>%
  group_by(exposure_control) %>%
  summarise(
    RMSE = sqrt(mean(bias^2)),
    coverage = mean(ci_coverage) * 100,
    .groups = "drop"
  )
print(exp_effect)

# Sparsity Effect
cat("\n=== Effect of Sparsity ===\n")
sparsity_summary <- factorial_df %>%
  group_by(sparsity, selection_method) %>%
  summarise(
    RMSE = sqrt(mean(bias^2)),
    coverage = mean(ci_coverage) * 100,
    .groups = "drop"
  )
print(sparsity_summary)

# Parameter Error Effect
cat("\n=== Effect of Parameter Error ===\n")
error_summary <- error_factorial_df %>%
  group_by(sigma_err, selection_method) %>%
  summarise(
    RMSE_noisy = sqrt(mean(bias_noisy^2)),
    RMSE_oracle = sqrt(mean(bias_oracle^2)),
    coverage_noisy = mean(ci_coverage_noisy) * 100,
    coverage_oracle = mean(ci_coverage_oracle) * 100,
    .groups = "drop"
  )
print(error_summary)

# Method Comparison across Factorial
cat("\n=== Method Comparison (Averaged over Factorial) ===\n")
method_summary <- factorial_df %>%
  group_by(selection_method) %>%
  summarise(
    RMSE = sqrt(mean(bias^2)),
    coverage = mean(ci_coverage) * 100,
    mean_sparse_ratio = mean(sparse_ratio) * 100,
    .groups = "drop"
  ) %>%
  arrange(RMSE)
print(method_summary)

# --- Cleanup parallel cluster ---
if (use_parallel && exists("cl")) {
  stopCluster(cl)
  cat("Parallel cluster stopped.\n")
}

# --- Complete ---

print_section("Robustness Study Complete")
cat("Results saved to:", data_dir, "\n")
cat("\nFiles created:\n")
cat("  - factorial_perfect_params.rds (2×2×2 factorial, perfect parameters)\n")
cat("  - factorial_param_error.rds (2×2×2 factorial, parameter error)\n")
cat("  - robustness_results_full.rds (combined)\n")
cat("\nFactorial Design:\n")
cat("  - Initialization: staircase vs censored_mle\n")
cat("  - Stratification: TRUE vs FALSE\n")
cat("  - Exposure Control: TRUE vs FALSE\n")
cat("\nNext steps:\n")
cat("  - For Analysis: source('script/05_analyze_robustness.R')\n")
cat("  - For Visualization: source('script/08_visualization.R')\n")
