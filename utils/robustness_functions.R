# =============================================================================
# Robustness Study Functions
# =============================================================================
# Description: Functions for running robustness study simulations
# Source: Extracted from script/03_robustness_study.R
# Dependencies: find_mle_censored, ci_asymptotic, sigmoid, p_fun, info_fun,
#               cW, run_staircase, to_catR_bank, catR_select,
#               create_BAS_structure, sympson_hetter_select (all from utils/)
# =============================================================================

#' Run Single Robustness Test with Factorial Controls
#'
#' @description
#' Runs a single CAT test under the 2x2x2 factorial design:
#'   initialization (staircase/censored), stratified (yes/no),
#'   exposure_control (yes/no).
#' Implements all selection methods: infomax, 4 CS-UCB variants,
#' mpwi, kl_information, progressive.
#' Uses TRUE item parameters (perfect knowledge scenario).
#'
#' @param theta_true True ability
#' @param item_bank Item bank data frame with columns a, b
#' @param selection_method Selection method name
#' @param config Configuration list
#' @param factorial_row Row from factorial design data frame
#' @param exposure_manager Exposure manager (optional)
#' @return List with test results
run_robustness_test_factorial <- function(theta_true, item_bank,
                                           selection_method, config,
                                           factorial_row,
                                           exposure_manager = NULL) {

  K <- nrow(item_bank)
  a <- item_bank$a
  b <- item_bank$b
  T_max <- config$T_max
  theta_min <- config$theta_range[1]
  theta_max <- config$theta_range[2]
  mle_min <- config$mle_bounds[1]     # Censored MLE bounds [-4, 4]
  mle_max <- config$mle_bounds[2]
  alpha <- config$alpha

  # Extract factorial conditions
  use_staircase <- factorial_row$initialization == "staircase"
  use_stratified <- factorial_row$stratified
  use_exposure <- factorial_row$exposure_control

  # Setup BAS stratification if enabled
  if (use_stratified) {
    bank_for_bas <- data.frame(item_id = 1:K, a = a, b = b)
    bank_for_bas <- create_BAS_structure(bank_for_bas, n_a_strata = 3)
    strata <- bank_for_bas$a_stratum
    n_strata <- max(strata, na.rm = TRUE)
    strata_idx <- lapply(1:n_strata, function(s) which(strata == s))
    names(strata_idx) <- paste0("stratum_", 1:n_strata)
    stratum_names <- names(strata_idx)
    stage_size <- ceiling(T_max / n_strata)
  }

  # Pre-compute catR bank format for standard selection methods
  itemBank_catR <- to_catR_bank(a, b)

  # Track sparse item usage
  is_sparse <- if ("is_sparse" %in% names(item_bank)) {
    item_bank$is_sparse
  } else {
    rep(FALSE, K)
  }

  # Initialize storage
  chosen <- integer(T_max)
  y <- numeric(T_max)
  theta_hat <- numeric(T_max + 1)
  theta_hat[1] <- 0
  used_items <- integer(0)

  # Regret tracking
  p_oracle <- 1 / (1 + exp(-a * (theta_true - b)))
  I_oracle <- a^2 * p_oracle * (1 - p_oracle)
  I_star <- max(I_oracle)
  inst_regret <- numeric(T_max)
  selection_times <- numeric(T_max)

  # Variables for CS-UCB
  amax <- max(a)
  diam_Theta <- theta_max - theta_min
  # CS-UCB tracking (for curvature lower bound computation)
  r_history <- numeric(T_max)
  r_history[1] <- diam_Theta / 2  # Initial radius
  H_lower_prev <- 1e-3

  # Staircase initialization phase
  t_start <- 1
  if (use_staircase && config$n_staircase > 0) {
    n_staircase <- min(config$n_staircase, T_max)

    # Simple staircase: start with medium difficulty items
    medium_items <- which(abs(b) <= 0.5)
    if (length(medium_items) == 0) medium_items <- 1:K

    current_b <- 0
    for (t in 1:n_staircase) {
      # Find closest unused item to current target difficulty
      available <- setdiff(medium_items, used_items)
      if (length(available) == 0) available <- setdiff(1:K, used_items)

      distances <- abs(b[available] - current_b)
      it <- available[which.min(distances)]

      chosen[t] <- it
      used_items <- c(used_items, it)

      # Generate response
      prob_true <- 1 / (1 + exp(-a[it] * (theta_true - b[it])))
      y[t] <- rbinom(1, 1, prob_true)

      # Adjust target difficulty (up if correct, down if incorrect)
      current_b <- current_b + ifelse(y[t] == 1, 0.5, -0.5)
      current_b <- max(min(current_b, theta_max), theta_min)

      # Update MLE
      theta_hat[t + 1] <- find_mle_censored(
        a_hist = a[chosen[1:t]],
        b_hist = b[chosen[1:t]],
        y_hist = y[1:t],
        theta_init = theta_hat[t],
        theta_min = mle_min,
        theta_max = mle_max
      )
    }
    t_start <- n_staircase + 1
  }

  # Main adaptive phase
  for (t in t_start:T_max) {
    current_theta <- theta_hat[t]

    # Determine available items
    if (use_stratified) {
      # Stage-based administration (Chang, Qian, & Ying 2001)
      current_stage <- min(ceiling(t / stage_size), n_strata)
      stratum_name <- stratum_names[current_stage]
      available_pool <- strata_idx[[stratum_name]]
      available <- setdiff(available_pool, used_items)
      if (length(available) == 0) {
        available <- setdiff(1:K, used_items)
      }
    } else {
      available <- setdiff(1:K, used_items)
    }

    if (length(available) == 0) break

    # Select item based on method (timing tracked)
    t_start_time <- proc.time()[3]

    if (selection_method == "infomax") {
      if (exists("info_vec_cpp")) {
        infos <- info_vec_cpp(current_theta, a[available], b[available])
      } else {
        infos <- sapply(available, function(i) {
          p <- 1 / (1 + exp(-a[i] * (current_theta - b[i])))
          a[i]^2 * p * (1 - p)
        })
      }
      it <- available[which.max(infos)]

    } else if (selection_method %in% c("csucb_theory", "csucb_obsinfo",
                                       "csucb_bernstein", "csucb_mixture")) {
      # --- Curvature H_t ---
      if (t > 1 && length(used_items) > 0) {
        # Observed Fisher info (used by obsinfo, bernstein, mixture)
        if (exists("info_sum_cpp")) {
          H_obs <- info_sum_cpp(a[used_items], b[used_items], current_theta)
        } else {
          p_vals <- 1 / (1 + exp(-a[used_items] * (current_theta - b[used_items])))
          H_obs <- max(sum(a[used_items]^2 * p_vals * (1 - p_vals)), 1e-3)
        }

        if (selection_method == "csucb_theory") {
          # Theory curvature lower bound H_t^- (Lemma 2)
          if (exists("H_lower_cpp")) {
            H_t_minus <- H_lower_cpp(a[used_items], r_history, diam_Theta,
                                     length(used_items))
          } else {
            H_t_minus <- 0
            for (s in seq_along(used_items)) {
              W_s <- r_history[s] + diam_Theta
              H_t_minus <- H_t_minus + cW(a[used_items[s]], W_s)
            }
          }
          H_lower_prev <- max(H_t_minus, H_lower_prev, 1e-3)
          H_t <- H_lower_prev
        } else {
          H_t <- H_obs  # obsinfo, bernstein, mixture use observed info
        }
      } else {
        H_t <- 1e-3
        H_obs <- 1e-3
      }

      # --- Confidence Radius (method-specific) ---
      delta_t_val <- (6 * alpha) / (pi^2 * t^2)

      if (selection_method %in% c("csucb_theory", "csucb_obsinfo")) {
        # Azuma-Hoeffding: B_t = a_max * sqrt(2t * log(2/δ_t))
        B_t <- amax * sqrt(2 * t * log(2 / delta_t_val))
        r_t <- B_t / max(H_t, 1e-3)

      } else if (selection_method == "csucb_bernstein") {
        # Same Bernstein family as the precision study, with delta_t in place
        # of alpha to obtain a time-indexed exploration radius.
        r_t <- bernstein_radius(
          info = H_obs,
          alpha_level = delta_t_val,
          max_score = amax,
          denominator = H_t,
          min_info = 1e-3
        )

      } else if (selection_method == "csucb_mixture") {
        # Same Laplace-approximated LR / mixture family as the precision study,
        # with delta_t used as the exploration budget at time t.
        r_t <- lr_mixture_radius(
          info = H_obs,
          alpha_level = delta_t_val,
          tau_sq = 1,
          min_info = 1e-3
        )
      }

      # --- Radius capping (shared) ---
      r_t <- min(r_t, diam_Theta)
      max_radius <- min(current_theta - mle_min, mle_max - current_theta)
      r_t <- min(r_t, max(max_radius, 0.1))
      r_history[t] <- r_t

      # --- UCB selection: sup_{θ ∈ C_t} I_i(θ) ---
      if (exists("ucb_select_cpp")) {
        it <- ucb_select_cpp(a, b, as.integer(available), current_theta, r_t)
      } else {
        ucb <- sapply(available, function(i) {
          theta_lower <- current_theta - r_t
          theta_upper <- current_theta + r_t
          if (b[i] >= theta_lower && b[i] <= theta_upper) {
            return(a[i]^2 / 4)
          }
          p_lower <- 1 / (1 + exp(-a[i] * (theta_lower - b[i])))
          p_upper <- 1 / (1 + exp(-a[i] * (theta_upper - b[i])))
          max(a[i]^2 * p_lower * (1 - p_lower),
              a[i]^2 * p_upper * (1 - p_upper))
        })
        it <- available[which.max(ucb)]
      }

    } else if (selection_method %in% c("mpwi", "kl_information", "progressive")) {
      # Traditional methods via switchable backend (C++ or catR)
      it <- traditional_select(
        method = selection_method,
        a = a, b = b,
        available = available,
        current_theta = current_theta,
        a_hist = if (t > 1) a[chosen[1:(t-1)]] else numeric(0),
        b_hist = if (t > 1) b[chosen[1:(t-1)]] else numeric(0),
        y_hist = if (t > 1) y[1:(t-1)] else integer(0),
        n_admin = as.integer(t - 1),
        t = t, T_max = T_max,
        K = K, itemBank_catR = itemBank_catR
      )
    }

    # Apply Sympson-Hetter exposure control AFTER selection
    if (use_exposure && !is.null(exposure_manager)) {
      it <- sympson_hetter_select(exposure_manager, it, available)
    }

    selection_times[t] <- proc.time()[3] - t_start_time

    chosen[t] <- it
    used_items <- c(used_items, it)

    # Track regret
    inst_regret[t] <- I_star - info_fun(theta_true, a[it], b[it])

    # Generate response using TRUE parameters
    prob_true <- 1 / (1 + exp(-a[it] * (theta_true - b[it])))
    y[t] <- rbinom(1, 1, prob_true)

    # Update MLE
    theta_hat[t + 1] <- find_mle_censored(
      a_hist = a[chosen[1:t]],
      b_hist = b[chosen[1:t]],
      y_hist = y[1:t],
      theta_init = current_theta,
      theta_min = mle_min,
      theta_max = mle_max
    )
  }

  # Final results
  T_actual <- sum(chosen > 0)
  theta_final <- theta_hat[T_actual + 1]

  # Compute CI coverage (asymptotic)
  ci <- ci_asymptotic(theta_final, a[chosen[1:T_actual]],
                      b[chosen[1:T_actual]], alpha)
  covered <- theta_true >= ci$lower && theta_true <= ci$upper

  # Sparse item usage
  n_sparse_used <- sum(is_sparse[chosen[1:T_actual]])

  list(
    theta_true = theta_true,
    theta_final = theta_final,
    bias = theta_final - theta_true,
    abs_error = abs(theta_final - theta_true),
    T_actual = T_actual,
    ci_coverage = covered,
    ci_width = ci$width,
    n_sparse_used = n_sparse_used,
    sparse_ratio = n_sparse_used / T_actual,
    used_items = chosen[1:T_actual],
    cum_regret = sum(inst_regret[1:T_actual]),
    mean_selection_time = mean(selection_times[1:T_actual]),
    # Per-trial trace data
    theta_path = theta_hat[1:(T_actual + 1)],
    inst_regret_path = inst_regret[1:T_actual],
    selection_time_path = selection_times[1:T_actual]
  )
}

#' Run Robustness Test with Parameter Error
#'
#' @description
#' Variant of run_robustness_test_factorial for parameter error scenario.
#' Uses NOISY parameters for item selection and MLE,
#' but TRUE parameters for response generation.
#' Runs two parallel MLE streams (noisy and oracle).
#'
#' @param theta_true True ability
#' @param error_bank Error bank with true and noisy parameters
#' @param selection_method Selection method name
#' @param config Configuration list
#' @param factorial_row Factorial design row
#' @param exposure_manager Exposure manager (optional)
#' @param noise_seed Optional integer; when supplied with sigma_err, the noisy
#'   item parameters are redrawn per call (Design A — Cheng & Patton 2009 /
#'   Patton et al. 2013 standard). When NULL the function falls back to the
#'   pre-baked `bank$a_noisy` / `bank$b_noisy` (Design B; backward-compat for
#'   older callers and re-runnable cached outputs).
#' @param sigma_err Noise SD on `a` (b uses 2x). Required for per-rep redraw.
#' @param a_bounds Clip bounds for a after noise. Default c(0.1, 3.0).
#' @param b_bounds Clip bounds for b after noise. Default c(-4, 4).
#' @return List with test results including noisy and oracle estimates
run_robustness_test_error_factorial <- function(theta_true, error_bank,
                                                 selection_method, config,
                                                 factorial_row,
                                                 exposure_manager = NULL,
                                                 noise_seed = NULL,
                                                 sigma_err = NULL,
                                                 a_bounds = c(0.1, 3.0),
                                                 b_bounds = c(-4, 4)) {

  bank <- error_bank$bank
  K <- nrow(bank)

  # True parameters (for response generation)
  a_true <- bank$a_true
  b_true <- bank$b_true

  # === Patch 3 (08_calibration_error.md §6.3) — per-rep noise redraw ===
  # Design A: fresh noise per replication when noise_seed + sigma_err provided.
  # Design B fallback: use pre-baked noisy params when args are NULL (kept for
  # backward compatibility with older callers and cached runs).
  if (!is.null(noise_seed) && !is.null(sigma_err)) {
    set.seed(as.integer(noise_seed))
    err_a <- rnorm(K, 0, sigma_err)
    err_b <- rnorm(K, 0, 2 * sigma_err)
    a_noisy <- pmax(pmin(a_true + err_a, a_bounds[2]), a_bounds[1])
    b_noisy <- pmax(pmin(b_true + err_b, b_bounds[2]), b_bounds[1])
  } else {
    a_noisy <- bank$a_noisy
    b_noisy <- bank$b_noisy
  }

  T_max <- config$T_max
  theta_min <- config$theta_range[1]
  theta_max <- config$theta_range[2]
  mle_min <- config$mle_bounds[1]     # Censored MLE bounds [-4, 4]
  mle_max <- config$mle_bounds[2]
  alpha <- config$alpha

  # Extract factorial conditions
  use_staircase <- factorial_row$initialization == "staircase"
  use_stratified <- factorial_row$stratified
  use_exposure <- factorial_row$exposure_control

  # Pre-compute catR bank format (using noisy params for selection)
  itemBank_catR <- to_catR_bank(a_noisy, b_noisy)

  # Setup stratification based on NOISY parameters
  if (use_stratified) {
    bank_for_bas <- data.frame(item_id = 1:K,
                               a = a_noisy, b = b_noisy)
    bank_for_bas <- create_BAS_structure(bank_for_bas,
                                         n_a_strata = 3)
    strata <- bank_for_bas$a_stratum
    n_strata <- max(strata, na.rm = TRUE)
    strata_idx <- lapply(1:n_strata,
                         function(s) which(strata == s))
    names(strata_idx) <- paste0("stratum_", 1:n_strata)
    stratum_names <- names(strata_idx)
    stage_size <- ceiling(T_max / n_strata)
  }

  # Initialize storage
  chosen <- integer(T_max)
  y <- numeric(T_max)
  theta_hat <- numeric(T_max + 1)
  theta_oracle <- numeric(T_max + 1)
  theta_hat[1] <- 0
  theta_oracle[1] <- 0
  used_items <- integer(0)

  # Regret tracking (using TRUE params for oracle regret)
  p_oracle <- 1 / (1 + exp(-a_true * (theta_true - b_true)))
  I_oracle <- a_true^2 * p_oracle * (1 - p_oracle)
  I_star <- max(I_oracle)
  inst_regret <- numeric(T_max)
  selection_times <- numeric(T_max)

  # Variables for CS-UCB
  amax <- max(a_noisy)
  diam_Theta <- theta_max - theta_min
  # CS-UCB tracking (for curvature lower bound computation)
  r_history <- numeric(T_max)
  r_history[1] <- diam_Theta / 2  # Initial radius
  H_lower_prev <- 1e-3

  # Staircase initialization phase
  t_start <- 1
  if (use_staircase && config$n_staircase > 0) {
    n_staircase <- min(config$n_staircase, T_max)

    medium_items <- which(abs(b_noisy) <= 0.5)
    if (length(medium_items) == 0) medium_items <- 1:K

    current_b <- 0
    for (t in 1:n_staircase) {
      available <- setdiff(medium_items, used_items)
      if (length(available) == 0) available <- setdiff(1:K, used_items)

      distances <- abs(b_noisy[available] - current_b)
      it <- available[which.min(distances)]

      chosen[t] <- it
      used_items <- c(used_items, it)

      # Response uses TRUE parameters
      prob_true <- 1 / (1 + exp(-a_true[it] * (theta_true - b_true[it])))
      y[t] <- rbinom(1, 1, prob_true)

      current_b <- current_b + ifelse(y[t] == 1, 0.5, -0.5)
      current_b <- max(min(current_b, theta_max), theta_min)

      # Update MLE using NOISY parameters
      theta_hat[t + 1] <- find_mle_censored(
        a_hist = a_noisy[chosen[1:t]],
        b_hist = b_noisy[chosen[1:t]],
        y_hist = y[1:t],
        theta_init = theta_hat[t],
        theta_min = mle_min,
        theta_max = mle_max
      )

      # Oracle MLE using TRUE parameters
      theta_oracle[t + 1] <- find_mle_censored(
        a_hist = a_true[chosen[1:t]],
        b_hist = b_true[chosen[1:t]],
        y_hist = y[1:t],
        theta_init = theta_oracle[t],
        theta_min = mle_min,
        theta_max = mle_max
      )
    }
    t_start <- n_staircase + 1
  }

  # Main adaptive phase
  for (t in t_start:T_max) {
    current_theta <- theta_hat[t]

    # Determine available items
    if (use_stratified) {
      # Stage-based administration (Chang, Qian, & Ying 2001)
      current_stage <- min(ceiling(t / stage_size), n_strata)
      stratum_name <- stratum_names[current_stage]
      available_pool <- strata_idx[[stratum_name]]
      available <- setdiff(available_pool, used_items)
      if (length(available) == 0) {
        available <- setdiff(1:K, used_items)
      }
    } else {
      available <- setdiff(1:K, used_items)
    }

    if (length(available) == 0) break

    # Select item (timing tracked)
    t_start_time <- proc.time()[3]

    # Selection using NOISY parameters
    if (selection_method == "infomax") {
      if (exists("info_vec_cpp")) {
        infos <- info_vec_cpp(current_theta, a_noisy[available], b_noisy[available])
      } else {
        infos <- sapply(available, function(i) {
          p <- 1 / (1 + exp(-a_noisy[i] * (current_theta - b_noisy[i])))
          a_noisy[i]^2 * p * (1 - p)
        })
      }
      it <- available[which.max(infos)]

    } else if (selection_method %in% c("csucb_theory", "csucb_obsinfo",
                                       "csucb_bernstein", "csucb_mixture")) {
      # --- Curvature H_t ---
      if (t > 1 && length(used_items) > 0) {
        if (exists("info_sum_cpp")) {
          H_obs <- info_sum_cpp(a_noisy[used_items], b_noisy[used_items], current_theta)
        } else {
          p_vals <- 1 / (1 + exp(-a_noisy[used_items] *
                                  (current_theta - b_noisy[used_items])))
          H_obs <- max(sum(a_noisy[used_items]^2 * p_vals * (1 - p_vals)), 1e-3)
        }

        if (selection_method == "csucb_theory") {
          if (exists("H_lower_cpp")) {
            H_t_minus <- H_lower_cpp(a_noisy[used_items], r_history, diam_Theta,
                                     length(used_items))
          } else {
            H_t_minus <- 0
            for (s in seq_along(used_items)) {
              W_s <- r_history[s] + diam_Theta
              H_t_minus <- H_t_minus + cW(a_noisy[used_items[s]], W_s)
            }
          }
          H_lower_prev <- max(H_t_minus, H_lower_prev, 1e-3)
          H_t <- H_lower_prev
        } else {
          H_t <- H_obs
        }
      } else {
        H_t <- 1e-3
        H_obs <- 1e-3
      }

      # --- Confidence Radius (method-specific) ---
      delta_t_val <- (6 * alpha) / (pi^2 * t^2)

      if (selection_method %in% c("csucb_theory", "csucb_obsinfo")) {
        B_t <- amax * sqrt(2 * t * log(2 / delta_t_val))
        r_t <- B_t / max(H_t, 1e-3)

      } else if (selection_method == "csucb_bernstein") {
        r_t <- bernstein_radius(
          info = H_obs,
          alpha_level = delta_t_val,
          max_score = amax,
          denominator = H_t,
          min_info = 1e-3
        )

      } else if (selection_method == "csucb_mixture") {
        r_t <- lr_mixture_radius(
          info = H_obs,
          alpha_level = delta_t_val,
          tau_sq = 1,
          min_info = 1e-3
        )
      }

      # --- Radius capping (shared) ---
      r_t <- min(r_t, diam_Theta)
      max_radius <- min(current_theta - mle_min, mle_max - current_theta)
      r_t <- min(r_t, max(max_radius, 0.1))
      r_history[t] <- r_t

      # --- UCB selection: sup_{θ ∈ C_t} I_i(θ) ---
      if (exists("ucb_select_cpp")) {
        it <- ucb_select_cpp(a_noisy, b_noisy, as.integer(available),
                             current_theta, r_t)
      } else {
        ucb <- sapply(available, function(i) {
          theta_lower <- current_theta - r_t
          theta_upper <- current_theta + r_t
          if (b_noisy[i] >= theta_lower && b_noisy[i] <= theta_upper) {
            return(a_noisy[i]^2 / 4)
          }
          p_lower <- 1 / (1 + exp(-a_noisy[i] * (theta_lower - b_noisy[i])))
          p_upper <- 1 / (1 + exp(-a_noisy[i] * (theta_upper - b_noisy[i])))
          max(a_noisy[i]^2 * p_lower * (1 - p_lower),
              a_noisy[i]^2 * p_upper * (1 - p_upper))
        })
        it <- available[which.max(ucb)]
      }

    } else if (selection_method %in% c("mpwi", "kl_information", "progressive")) {
      # Traditional methods via switchable backend (using noisy parameters)
      it <- traditional_select(
        method = selection_method,
        a = a_noisy, b = b_noisy,
        available = available,
        current_theta = current_theta,
        a_hist = if (t > 1) a_noisy[chosen[1:(t-1)]] else numeric(0),
        b_hist = if (t > 1) b_noisy[chosen[1:(t-1)]] else numeric(0),
        y_hist = if (t > 1) y[1:(t-1)] else integer(0),
        n_admin = as.integer(t - 1),
        t = t, T_max = T_max,
        K = K, itemBank_catR = itemBank_catR
      )
    }

    # Apply Sympson-Hetter exposure control AFTER selection
    if (use_exposure && !is.null(exposure_manager)) {
      it <- sympson_hetter_select(exposure_manager, it, available)
    }

    selection_times[t] <- proc.time()[3] - t_start_time

    chosen[t] <- it
    used_items <- c(used_items, it)

    # Track regret (using TRUE params)
    inst_regret[t] <- I_star - info_fun(theta_true, a_true[it], b_true[it])

    # Response uses TRUE parameters
    prob_true <- 1 / (1 + exp(-a_true[it] * (theta_true - b_true[it])))
    y[t] <- rbinom(1, 1, prob_true)

    # MLE using NOISY parameters
    theta_hat[t + 1] <- find_mle_censored(
      a_hist = a_noisy[chosen[1:t]],
      b_hist = b_noisy[chosen[1:t]],
      y_hist = y[1:t],
      theta_init = current_theta,
      theta_min = mle_min,
      theta_max = mle_max
    )

    # Oracle MLE using TRUE parameters
    theta_oracle[t + 1] <- find_mle_censored(
      a_hist = a_true[chosen[1:t]],
      b_hist = b_true[chosen[1:t]],
      y_hist = y[1:t],
      theta_init = theta_oracle[t],
      theta_min = mle_min,
      theta_max = mle_max
    )
  }

  # Final results
  T_actual <- sum(chosen > 0)
  theta_final_noisy <- theta_hat[T_actual + 1]
  theta_final_oracle <- theta_oracle[T_actual + 1]

  # CI using noisy parameters
  ci_noisy <- ci_asymptotic(theta_final_noisy,
                            a_noisy[chosen[1:T_actual]],
                            b_noisy[chosen[1:T_actual]], alpha)
  covered_noisy <- theta_true >= ci_noisy$lower && theta_true <= ci_noisy$upper

  # CI using true parameters (oracle)
  ci_oracle <- ci_asymptotic(theta_final_oracle,
                              a_true[chosen[1:T_actual]],
                              b_true[chosen[1:T_actual]], alpha)
  covered_oracle <- theta_true >= ci_oracle$lower && theta_true <= ci_oracle$upper

  list(
    theta_true = theta_true,
    theta_final_noisy = theta_final_noisy,
    theta_final_oracle = theta_final_oracle,
    bias_noisy = theta_final_noisy - theta_true,
    bias_oracle = theta_final_oracle - theta_true,
    T_actual = T_actual,
    ci_coverage_noisy = covered_noisy,
    ci_coverage_oracle = covered_oracle,
    ci_width_noisy = ci_noisy$width,
    ci_width_oracle = ci_oracle$width,
    used_items = chosen[1:T_actual],
    cum_regret = sum(inst_regret[1:T_actual]),
    mean_selection_time = mean(selection_times[1:T_actual]),
    # Per-trial trace data
    theta_path = theta_hat[1:(T_actual + 1)],
    theta_oracle_path = theta_oracle[1:(T_actual + 1)],
    inst_regret_path = inst_regret[1:T_actual],
    selection_time_path = selection_times[1:T_actual]
  )
}

message("=== Robustness functions loaded ===")
