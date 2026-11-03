# =============================================================================
# Staircase Initialization Functions
# =============================================================================
# Description: Staircase procedure for CAT initialization
#              Shortened from 10 items to 3-5 items per committee feedback
# =============================================================================

#' Run Staircase Procedure
#'
#' @description
#' Implements the staircase initialization method for CAT.
#' Uses a simple up-down rule based on response correctness.
#' Per committee feedback: shortened to 3 or 5 items (from original 10).
#'
#' @param theta0 True ability (for simulation, generates responses)
#' @param a Vector of item discrimination parameters
#' @param b Vector of item difficulty parameters
#' @param n_staircase Number of staircase items (default = 5, can be 3)
#' @param step_size Step size for ability update (default = 0.5)
#' @param theta_init Initial ability estimate (default = 0)
#' @param theta_min Lower bound for ability
#' @param theta_max Upper bound for ability
#' @param b_tolerance Tolerance for finding items near target difficulty
#' @param verbose Print progress
#' @return List containing staircase results
#' @export
#'
#' @examples
#' bank <- make_item_bank(K = 300, seed = 123)
#' result <- run_staircase(theta0 = 1.0, a = bank$a, b = bank$b, n_staircase = 5)
run_staircase <- function(theta0, a, b, n_staircase = 5,
                          step_size = 0.5,
                          theta_init = 0,
                          theta_min = -4, theta_max = 4,
                          b_tolerance = 0.5,
                          verbose = FALSE) {

  K <- length(a)

  # Initialize storage
  chosen <- integer(n_staircase)
  y <- numeric(n_staircase)
  theta_current <- theta_init
  theta_path <- numeric(n_staircase + 1)
  theta_path[1] <- theta_current

  # Track used items
  used_items <- integer(0)

  for (t in 1:n_staircase) {
    if (verbose) cat("  Staircase step", t, "- target b:", round(theta_current, 2), "\n")

    # Find items close to current ability estimate (difficulty matching)
    available <- setdiff(1:K, used_items)
    b_diff <- abs(b[available] - theta_current)

    # Select item with difficulty closest to current theta
    # If multiple items within tolerance, pick randomly for diversity
    close_items <- available[b_diff <= b_tolerance]

    if (length(close_items) == 0) {
      # If no items within tolerance, pick closest
      it <- available[which.min(b_diff)]
    } else if (length(close_items) == 1) {
      it <- close_items
    } else {
      # Random selection among close items for item bank diversity
      it <- sample(close_items, 1)
    }

    chosen[t] <- it
    used_items <- c(used_items, it)

    # Generate response based on true ability
    prob_true <- 1 / (1 + exp(-a[it] * (theta0 - b[it])))
    y[t] <- rbinom(1, 1, prob_true)

    # Update theta estimate using staircase rule
    if (y[t] == 1) {
      # Correct response: increase difficulty
      theta_current <- min(theta_current + step_size, theta_max)
    } else {
      # Incorrect response: decrease difficulty
      theta_current <- max(theta_current - step_size, theta_min)
    }

    theta_path[t + 1] <- theta_current
  }

  # Final theta estimate from staircase
  theta_staircase <- theta_current

  return(list(
    method = paste0("Staircase-", n_staircase),
    theta_staircase = theta_staircase,
    theta_path = theta_path,
    chosen = chosen,
    y = y,
    n_items = n_staircase,
    a_hist = a[chosen],
    b_hist = b[chosen]
  ))
}

message("=== Staircase functions loaded ===")
