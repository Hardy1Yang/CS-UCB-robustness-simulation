# =============================================================================
# Maximum Likelihood Estimation Functions
# =============================================================================
# Description: Functions for computing Maximum Likelihood Estimates (MLE)
#              with censoring and bounds for numerical stability.
# Source: Migrated from OCO_CAT/finding/utils/mle_functions.R
# =============================================================================

#' Censored MLE Estimation
#'
#' @description
#' Computes MLE with bounds (censoring) for ability estimation.
#' Handles boundary cases (all correct/incorrect) appropriately.
#'
#' @param a_hist Vector of discrimination parameters for administered items
#' @param b_hist Vector of difficulty parameters for administered items
#' @param y_hist Vector of responses (0 or 1)
#' @param theta_init Initial value for theta (default = 0)
#' @param theta_min Lower bound for theta (default = -4)
#' @param theta_max Upper bound for theta (default = 4)
#' @param tol Tolerance for convergence (default = 1e-6)
#' @param max_iter Maximum iterations (default = 50)
#' @return Estimated theta value (censored to bounds)
#' @export
#' @examples
#' # Simple example with 3 items
#' a <- c(1, 1.2, 0.8)
#' b <- c(-1, 0, 1)
#' y <- c(1, 1, 0)
#' theta_hat <- find_mle_censored(a, b, y)
find_mle_censored <- function(a_hist, b_hist, y_hist,
                              theta_init = 0,
                              theta_min = -4, theta_max = 4,
                              tol = 1e-6, max_iter = 50) {

  # Handle empty case
  if (length(y_hist) == 0) return(theta_init)

  # Check if all responses are the same (boundary cases)
  if (all(y_hist == 1)) return(theta_max)
  if (all(y_hist == 0)) return(theta_min)

  # Score function for 2PL model
  score_fun_local <- function(theta) {
    p <- 1 / (1 + exp(-a_hist * (theta - b_hist)))
    sum(a_hist * (y_hist - p))
  }

  # Check boundary conditions
  score_at_min <- score_fun_local(theta_min)
  score_at_max <- score_fun_local(theta_max)

  if (score_at_min <= 0) return(theta_min)
  if (score_at_max >= 0) return(theta_max)

  # Try uniroot first
  result <- tryCatch({
    uniroot(score_fun_local,
            interval = c(theta_min, theta_max),
            tol = tol, maxiter = max_iter)$root
  }, error = function(e) {
    # Fallback to Newton-Raphson if uniroot fails
    theta <- theta_init

    for (iter in 1:max_iter) {
      p <- 1 / (1 + exp(-a_hist * (theta - b_hist)))
      score <- sum(a_hist * (y_hist - p))

      if (abs(score) < tol) break

      # Information (negative Hessian)
      info <- sum(a_hist^2 * p * (1 - p))

      if (info > 1e-8) {
        theta_new <- theta + score / info
        theta_new <- pmax(theta_min, pmin(theta_max, theta_new))
      } else {
        # Information too small, return current estimate
        break
      }

      if (abs(theta_new - theta) < tol) break
      theta <- theta_new
    }

    return(theta)
  })

  return(result)
}

message("=== MLE functions loaded ===")
