# =============================================================================
# Exposure Control Functions
# =============================================================================
# Description: Comprehensive exposure control mechanisms for CAT simulation
#              Per committee feedback: limit how often each item is
#              administered across examinees (a maximum exposure rate).
# =============================================================================

# --- Exposure Control Manager Class ---

#' Create Exposure Control Manager
#'
#' @description
#' Creates an exposure control manager that tracks item usage across examinees
#' and controls exposure rates. This implements true cross-examinee exposure control.
#'
#' @param n_items Number of items in the bank
#' @param max_rate Maximum exposure rate (default = 0.20)
#' @param method Exposure control method ("sympson_hetter", "randomesque", "progressive")
#' @return Exposure control manager (environment)
#' @export
create_exposure_manager <- function(n_items, max_rate = 0.20, method = "sympson_hetter") {

  manager <- new.env()

  # Item tracking
  manager$n_items <- n_items
  manager$max_rate <- max_rate
  manager$method <- method

  # Exposure counts and control parameters
  manager$exposure_counts <- rep(0, n_items)      # Times each item was administered
  manager$selection_counts <- rep(0, n_items)     # Times each item was selected (before control)
  manager$n_examinees <- 0                        # Total examinees tested

  # Sympson-Hetter control parameters (calibrated iteratively)
  manager$control_params <- rep(1.0, n_items)     # Initial control parameters

  # Lock out items that exceed max rate
  manager$locked_items <- rep(FALSE, n_items)

  class(manager) <- c("exposure_manager", "environment")
  return(manager)
}

#' Update Exposure Manager After Test
#'
#' @description
#' Updates exposure counts after an examinee completes a test.
#'
#' @param manager Exposure manager object
#' @param administered_items Vector of administered item indices
#' @param selected_items Vector of originally selected items (before control)
#' @return Updated manager (invisibly)
#' @export
update_exposure_manager <- function(manager, administered_items, selected_items = NULL) {

  manager$n_examinees <- manager$n_examinees + 1

  # Update exposure counts
  manager$exposure_counts[administered_items] <-
    manager$exposure_counts[administered_items] + 1

  # Update selection counts (if provided)
  if (!is.null(selected_items)) {
    manager$selection_counts[selected_items] <-
      manager$selection_counts[selected_items] + 1
  }

  # Check for items exceeding max rate and lock them
  current_rates <- manager$exposure_counts / manager$n_examinees
  manager$locked_items <- current_rates >= manager$max_rate

  invisible(manager)
}

#' Get Current Exposure Rates
#'
#' @param manager Exposure manager object
#' @return Vector of current exposure rates
#' @export
get_exposure_rates <- function(manager) {
  if (manager$n_examinees == 0) {
    return(rep(0, manager$n_items))
  }
  manager$exposure_counts / manager$n_examinees
}

#' Get Available Items (Not Locked)
#'
#' @param manager Exposure manager object
#' @param administered Already administered items in current test
#' @return Vector of available item indices
#' @export
get_available_items <- function(manager, administered = integer(0)) {
  all_items <- 1:manager$n_items
  available <- setdiff(all_items, administered)
  available <- available[!manager$locked_items[available]]
  return(available)
}

# --- Theta Sampler for SH Calibration ---

#' Sample ability values for SH calibration / production
#'
#' Supports normal (clipped) and truncated normal (rejection-sampled).
#' Defaults align with Investigation H arm C: TN(0, 1.5, -3, 3).
#'
#' @param n Number of values to sample
#' @param dist list with type/mean/sd/lo/hi; type missing → "normal"
#' @return Numeric vector of length n
sample_theta_for_sh <- function(n, dist = list(type = "truncated_normal",
                                                  mean = 0, sd = 1.5,
                                                  lo = -3, hi = 3)) {
  type <- dist$type %||% "normal"
  mean <- dist$mean %||% 0
  sd   <- dist$sd   %||% 1
  lo   <- dist$lo   %||% -3
  hi   <- dist$hi   %||%  3

  if (type == "normal") {
    return(pmax(pmin(rnorm(n, mean, sd), hi), lo))
  }
  if (type == "uniform") {
    return(runif(n, lo, hi))
  }
  if (type == "truncated_normal") {
    out <- numeric(n)
    filled <- 0L
    while (filled < n) {
      need <- n - filled
      batch <- rnorm(ceiling(need * 1.3), mean, sd)
      keep <- batch[batch >= lo & batch <= hi]
      take <- min(length(keep), need)
      if (take > 0L) {
        out[(filled + 1L):(filled + take)] <- keep[seq_len(take)]
        filled <- filled + take
      }
    }
    return(out)
  }
  stop(sprintf("sample_theta_for_sh: unknown type '%s'", type))
}

# Local null-coalescing operator (avoid dependency on rlang)
`%||%` <- function(a, b) if (is.null(a)) b else a

# --- Sympson-Hetter Exposure Control ---

#' Sympson-Hetter Exposure Control (Enhanced)
#'
#' @description
#' Implements the Sympson-Hetter (1985) exposure control method.
#' Uses control parameters that are calibrated to achieve target exposure rates.
#'
#' @param manager Exposure manager object
#' @param selected_item Index of item selected by CAT algorithm
#' @param available_items Vector of currently available items
#' @param fallback_method How to handle rejected items ("random", "next_best", "info_weighted")
#' @param info_values Optional vector of information values for weighted fallback
#' @return Final administered item index
#' @export
sympson_hetter_select <- function(manager, selected_item, available_items,
                                  fallback_method = "random", info_values = NULL) {

  # Get control probability for selected item
  control_prob <- manager$control_params[selected_item]

  # Apply probabilistic control
  if (runif(1) <= control_prob) {
    # Accept the selected item
    return(selected_item)
  }

  # Item rejected - select alternative
  alternatives <- setdiff(available_items, selected_item)

  if (length(alternatives) == 0) {
    # No alternatives - must use selected item
    return(selected_item)
  }

  if (fallback_method == "random") {
    return(sample(alternatives, 1))

  } else if (fallback_method == "next_best" && !is.null(info_values)) {
    # Select item with next highest information
    alt_info <- info_values[alternatives]
    return(alternatives[which.max(alt_info)])

  } else if (fallback_method == "info_weighted" && !is.null(info_values)) {
    # Weighted random selection based on information
    alt_info <- info_values[alternatives]
    probs <- alt_info / sum(alt_info)
    return(sample(alternatives, 1, prob = probs))

  } else {
    return(sample(alternatives, 1))
  }
}

#' Calibrate Sympson-Hetter Parameters
#'
#' @description
#' Iteratively calibrates control parameters to achieve target exposure rates.
#' This should be run as a pre-simulation calibration step.
#'
#' @param item_bank Item bank data frame with a, b columns
#' @param target_rate Target maximum exposure rate
#' @param n_calibration Number of calibration iterations
#' @param n_examinees_per_iter Examinees per calibration iteration
#' @param selection_function Item selection function
#' @param theta_dist Distribution of abilities; list with fields type (one of
#'   "normal", "truncated_normal", "uniform"), mean, sd, lo, hi. Defaults match
#'   Investigation H arm C: TN(0, 1.5, -3, 3) matched with production.
#' @param T_max Test length
#' @param tolerance Convergence tolerance
#' @param verbose Print progress
#' @param backend "auto" routes to C++ sh_calibrate_cpp if available and
#'   selection_function is NULL (MFI default); "cpp" forces C++ (errors if
#'   unavailable); "r" forces pure-R implementation. C++ default for
#'   robustness re-runs per `script/test/results/benchmark_sh_calibration.json`
#'   (50-100x speedup; bit-identical reproducibility not guaranteed but final
#'   convergence outcome matches).
#' @return Calibrated control parameters
#' @export
calibrate_sympson_hetter <- function(item_bank, target_rate = 0.20,
                                     n_calibration = 50,
                                     n_examinees_per_iter = 1000,
                                     selection_function = NULL,
                                     theta_dist = list(
                                       type = "truncated_normal",
                                       mean = 0, sd = 1.5,
                                       lo = -3, hi = 3),
                                     T_max = 30,
                                     verbose = TRUE,
                                     backend = "auto") {
  if (!backend %in% c("auto", "cpp", "r")) {
    stop("calibrate_sympson_hetter: backend must be one of 'auto', 'cpp', 'r'")
  }

  use_cpp <- (backend %in% c("auto", "cpp")) &&
              is.null(selection_function) &&
              exists("sh_calibrate_cpp", mode = "function")
  if (backend == "cpp" && !use_cpp) {
    stop("calibrate_sympson_hetter: backend='cpp' requested but ",
         "sh_calibrate_cpp is not available or selection_function is not NULL ",
         "(C++ backend only supports default MFI selection).")
  }
  if (use_cpp) {
    type <- theta_dist$type %||% "normal"
    dist_code <- switch(type,
      truncated_normal = 0L,
      normal           = 1L,
      uniform          = 2L,
      stop(sprintf("calibrate_sympson_hetter: unknown theta_dist$type '%s'",
                   type)))
    res <- sh_calibrate_cpp(
      a = as.numeric(item_bank$a), b = as.numeric(item_bank$b),
      target_rate = target_rate,
      n_calibration = as.integer(n_calibration),
      n_examinees_per_iter = as.integer(n_examinees_per_iter),
      T_max = as.integer(T_max),
      dist_code = dist_code,
      dist_mean = theta_dist$mean %||% 0,
      dist_sd   = theta_dist$sd   %||% 1,
      dist_lo   = theta_dist$lo   %||% -3,
      dist_hi   = theta_dist$hi   %||%  3,
      verbose = isTRUE(verbose))
    return(res)
  }

  n_items <- nrow(item_bank)

  # Initialize control parameters
  control_params <- rep(1.0, n_items)

  # Default selection function (MFI)
  if (is.null(selection_function)) {
    selection_function <- function(theta, item_bank, available) {
      info <- sapply(available, function(i) {
        a <- item_bank$a[i]
        b <- item_bank$b[i]
        p <- 1 / (1 + exp(-a * (theta - b)))
        a^2 * p * (1 - p)
      })
      available[which.max(info)]
    }
  }

  max_rate_achieved <- NA
  converged <- FALSE

  for (iter in 1:n_calibration) {
    if (verbose) cat("Calibration iteration", iter, "of",
                     n_calibration, "\n")

    # Track selection and administration counts
    selection_counts <- rep(0, n_items)
    admin_counts <- rep(0, n_items)

    # Batch-sample ability values per iteration (TN by default)
    theta_pool <- sample_theta_for_sh(n_examinees_per_iter, theta_dist)

    # Simulate examinees
    for (e in 1:n_examinees_per_iter) {
      theta <- theta_pool[e]

      administered <- integer(0)

      for (t in 1:T_max) {
        available <- setdiff(1:n_items, administered)
        if (length(available) == 0) break

        # Select item
        selected <- selection_function(theta, item_bank, available)
        selection_counts[selected] <- selection_counts[selected] + 1

        # Apply exposure control
        if (runif(1) <= control_params[selected]) {
          administered <- c(administered, selected)
          admin_counts[selected] <- admin_counts[selected] + 1
        } else {
          # Select random alternative
          alternatives <- setdiff(available, selected)
          if (length(alternatives) > 0) {
            alt <- sample(alternatives, 1)
            administered <- c(administered, alt)
            admin_counts[alt] <- admin_counts[alt] + 1
          } else {
            administered <- c(administered, selected)
            admin_counts[selected] <- admin_counts[selected] + 1
          }
        }
      }
    }

    # Update control parameters (strict: no tolerance)
    for (i in 1:n_items) {
      if (selection_counts[i] > 0) {
        current_rate <- admin_counts[i] / n_examinees_per_iter

        if (current_rate > target_rate) {
          # Reduce control parameter
          control_params[i] <- control_params[i] *
            (target_rate / current_rate)
        } else if (current_rate < target_rate * 0.9 &&
                   control_params[i] < 1.0) {
          # Increase control parameter (but not above 1)
          adjustment <- min(1.0, control_params[i] * 1.1)
          control_params[i] <- adjustment
        }
      }
    }

    # Check convergence against the configured target itself.
    admin_rates <- admin_counts / n_examinees_per_iter
    max_rate_achieved <- max(admin_rates)

    if (verbose) {
      cat("  Max exposure rate:",
          round(max_rate_achieved, 3), "\n")
      cat("  Items exceeding target:",
          sum(admin_rates > target_rate), "\n")
    }

    if (max_rate_achieved <= target_rate) {
      if (verbose) cat("Converged at iteration", iter,
                       "(max_rate:", round(max_rate_achieved, 3), ")\n")
      converged <- TRUE
      break
    }
  }

  # Return diagnostics along with control parameters
  list(
    control_params = control_params,
    max_exposure = max_rate_achieved,
    n_overexposed = sum(admin_rates > target_rate),
    iterations = iter,
    converged = converged
  )
}

# --- Maximum Exposure Constraint (Hard Limit) ---

#' Apply Hard Exposure Limit
#'
#' @description
#' Enforces a hard maximum on item exposure (E_max = 1 per paper design).
#' Items that have reached the limit are excluded from selection.
#'
#' @param manager Exposure manager object
#' @param selected_item Index of selected item
#' @param available_items Vector of available items
#' @param max_uses Maximum times an item can be used (default = Inf for rate-based)
#' @return Administered item index (or NA if no items available)
#' @export
hard_exposure_limit <- function(manager, selected_item, available_items, max_uses = Inf) {

  # Check if selected item has reached limit
  if (manager$exposure_counts[selected_item] >= max_uses) {
    # Find alternatives that haven't reached limit
    alternatives <- setdiff(available_items, selected_item)
    valid_alts <- alternatives[manager$exposure_counts[alternatives] < max_uses]

    if (length(valid_alts) > 0) {
      return(sample(valid_alts, 1))
    } else {
      return(NA)  # No valid items available
    }
  }

  return(selected_item)
}

# --- Exposure Control Summary ---

#' Generate Exposure Control Report
#'
#' @param manager Exposure manager object
#' @return Data frame with exposure statistics
#' @export
exposure_report <- function(manager) {

  rates <- get_exposure_rates(manager)

  data.frame(
    metric = c("n_examinees", "n_items", "max_rate_target", "max_rate_achieved",
               "mean_rate", "median_rate", "sd_rate",
               "n_never_used", "n_overexposed", "n_locked"),
    value = c(
      manager$n_examinees,
      manager$n_items,
      manager$max_rate,
      max(rates),
      mean(rates),
      median(rates),
      sd(rates),
      sum(rates == 0),
      sum(rates > manager$max_rate),
      sum(manager$locked_items)
    )
  )
}

#' Print Exposure Manager Summary
#'
#' @param x Exposure manager object
#' @param ... Additional arguments
#' @export
print.exposure_manager <- function(x, ...) {
  rates <- get_exposure_rates(x)

  cat("=== Exposure Control Manager ===\n")
  cat("Method:", x$method, "\n")
  cat("Target max rate:", x$max_rate, "\n")
  cat("Examinees tested:", x$n_examinees, "\n")
  cat("Items in bank:", x$n_items, "\n")
  cat("\nExposure Statistics:\n")
  cat("  Achieved max rate:", round(max(rates), 4), "\n")
  cat("  Mean rate:", round(mean(rates), 4), "\n")
  cat("  Never used:", sum(rates == 0), "\n")
  cat("  Overexposed:", sum(rates > x$max_rate), "\n")
  cat("  Locked:", sum(x$locked_items), "\n")
}

message("=== Exposure control functions loaded ===")
