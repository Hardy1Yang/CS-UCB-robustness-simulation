# =============================================================================
# Parameter Error Model
# =============================================================================
# Description: Models for item parameter estimation error
#              Based on draft202601.tex Robustness Study Scenario 2:
#              "Item Parameter Estimation Error"
#
#              Error Model: Additive Gaussian noise
#              ã_i = max(a_i + ε_a, 0.1), ε_a ~ N(0, σ_err)
#              b̃_i = b_i + ε_b,           ε_b ~ N(0, 2σ_err)
# =============================================================================

# --- Error Generation Functions ---

#' Generate Noisy Item Parameters
#'
#' @description
#' Generates item parameters with calibration error.
#' Implements the error model from the paper:
#' - Error is additive Gaussian noise
#' - ε_a ~ N(0, σ_err), ε_b ~ N(0, 2σ_err)
#' - Errors are independent across items
#'
#' @param a_true True discrimination parameters
#' @param b_true True difficulty parameters
#' @param sigma_err Standard deviation of error (default = 0.1)
#' @param seed Random seed for reproducibility
#' @param a_bounds Bounds for a parameter (default c(0.1, 3.0))
#' @param b_bounds Bounds for b parameter (default c(-4, 4))
#' @return List with noisy parameters and error information
#' @export
generate_noisy_parameters <- function(a_true, b_true,
                                      sigma_err = 0.1,
                                      seed = NULL,
                                      a_bounds = c(0.1, 3.0),
                                      b_bounds = c(-4, 4)) {

  if (!is.null(seed)) set.seed(seed)

  n_items <- length(a_true)

  # Generate errors
  error_a <- rnorm(n_items, mean = 0, sd = sigma_err)
  error_b <- rnorm(n_items, mean = 0, sd = 2 * sigma_err)

  # Add errors to true parameters
  a_noisy <- a_true + error_a
  b_noisy <- b_true + error_b

  # Apply bounds (especially important for a > 0)
  a_noisy <- pmax(pmin(a_noisy, a_bounds[2]), a_bounds[1])
  b_noisy <- pmax(pmin(b_noisy, b_bounds[2]), b_bounds[1])

  list(
    a_noisy = a_noisy,
    b_noisy = b_noisy,
    a_true = a_true,
    b_true = b_true,
    error_a = error_a,
    error_b = error_b,
    sigma_err = sigma_err,
    rmse_a = sqrt(mean((a_noisy - a_true)^2)),
    rmse_b = sqrt(mean((b_noisy - b_true)^2)),
    correlation_a = cor(a_true, a_noisy),
    correlation_b = cor(b_true, b_noisy)
  )
}

message("=== Error model functions loaded ===")
