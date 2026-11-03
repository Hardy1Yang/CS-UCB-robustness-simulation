# =============================================================================
# Confidence Interval Methods
# =============================================================================
# Description: Multiple CI construction methods for CAT ability estimation
#              Based on draft202601.tex Section on CI methods:
#              - Asymptotic (Wald) CI
#              - Azuma-Hoeffding CI
#              - Bernstein CI
#              - Likelihood-Ratio / Mixture Confidence Sequences
# =============================================================================

# --- Asymptotic (Wald) CI ---

#' Asymptotic Confidence Interval
#'
#' @description
#' Standard asymptotic CI based on MLE and observed Fisher information.
#' CI = θ̂ ± z_{α/2} / √I(θ̂)
#'
#' @param theta_hat MLE estimate
#' @param a_hist Vector of administered item discriminations
#' @param b_hist Vector of administered item difficulties
#' @param alpha Confidence level (default = 0.05 for 95% CI)
#' @return List with lower, upper bounds and width
#' @export
ci_asymptotic <- function(theta_hat, a_hist, b_hist, alpha = 0.05) {

  # Calculate observed information
  info <- observed_info_2pl(theta_hat, a_hist, b_hist)

  if (info <= 0) {
    return(list(lower = -Inf, upper = Inf, width = Inf, se = NA))
  }

  se <- 1 / sqrt(info)
  z <- qnorm(1 - alpha/2)

  lower <- theta_hat - z * se
  upper <- theta_hat + z * se

  list(
    lower = lower,
    upper = upper,
    width = upper - lower,
    se = se,
    info = info,
    method = "Asymptotic"
  )
}

# --- Shared Helpers ---

#' Observed Fisher information under the 2PL model
#'
#' @keywords internal
observed_info_2pl <- function(theta_hat, a_hist, b_hist) {
  if (length(a_hist) == 0) {
    return(0)
  }

  p <- 1 / (1 + exp(-a_hist * (theta_hat - b_hist)))
  sum(a_hist^2 * p * (1 - p))
}

#' Bernstein score bound from the dissertation precision derivation
#'
#' @keywords internal
bernstein_score_bound <- function(info, alpha_level, max_score,
                                  min_info = 1e-12) {
  info_safe <- max(info, min_info)
  log_term <- log(2 / alpha_level)
  sqrt(2 * info_safe * log_term) + (max_score * log_term) / 3
}

#' Bernstein radius using a score bound divided by a curvature term
#'
#' @keywords internal
bernstein_radius <- function(info, alpha_level, max_score,
                             denominator = info,
                             min_info = 1e-12) {
  bound <- bernstein_score_bound(info, alpha_level, max_score, min_info)
  bound / max(denominator, min_info)
}

#' Laplace-approximated likelihood-ratio / mixture radius
#'
#' @keywords internal
lr_mixture_radius <- function(info, alpha_level, tau_sq = 1,
                              min_info = 1e-12) {
  info_safe <- max(info, min_info)
  sqrt((2 / info_safe) * log(sqrt(1 + info_safe * tau_sq) / alpha_level))
}

# --- Azuma-Hoeffding CI ---

#' Azuma-Hoeffding Confidence Interval
#'
#' @description
#' CI based on the dissertation's Azuma-Hoeffding plug-in radius:
#' a_max * sqrt(2 t log(2 / alpha)) / I_t(theta_hat).
#'
#' For 2PL model with bounded θ ∈ [θ_min, θ_max]:
#' |M_t - M_{t-1}| ≤ c_t where c_t depends on item parameters
#'
#' @param theta_hat MLE estimate
#' @param a_hist Vector of administered item discriminations
#' @param b_hist Vector of administered item difficulties
#' @param y_hist Vector of responses
#' @param alpha Confidence level
#' @param theta_min Lower bound for theta
#' @param theta_max Upper bound for theta
#' @return List with CI bounds
#' @export
ci_azuma_hoeffding <- function(theta_hat, a_hist, b_hist, y_hist,
                                alpha = 0.05,
                                theta_min = -4, theta_max = 4) {

  t <- length(y_hist)
  if (t == 0) {
    return(list(lower = theta_min, upper = theta_max, width = theta_max - theta_min,
                method = "Azuma-Hoeffding"))
  }

  # Dissertation precision-study formula:
  # r_t = a_max * sqrt(2 t log(2 / alpha)) / I_t(theta_hat)
  a_max <- max(a_hist)
  epsilon <- a_max * sqrt(2 * t * log(2 / alpha))

  # Calculate observed information for centering
  info <- observed_info_2pl(theta_hat, a_hist, b_hist)

  if (info <= 0) {
    return(list(lower = theta_min, upper = theta_max, width = theta_max - theta_min,
                method = "Azuma-Hoeffding"))
  }

  # CI radius
  radius <- epsilon / info

  lower <- max(theta_hat - radius, theta_min)
  upper <- min(theta_hat + radius, theta_max)

  list(
    lower = lower,
    upper = upper,
    width = upper - lower,
    epsilon = epsilon,
    info = info,
    method = "Azuma-Hoeffding"
  )
}

# --- Bernstein CI ---

#' Bernstein Confidence Interval
#'
#' @description
#' CI based on Bernstein's inequality, providing tighter bounds when
#' variance is small. Uses empirical variance estimate.
#'
#' @param theta_hat MLE estimate
#' @param a_hist Vector of administered item discriminations
#' @param b_hist Vector of administered item difficulties
#' @param y_hist Vector of responses
#' @param alpha Confidence level
#' @param theta_min Lower bound for theta
#' @param theta_max Upper bound for theta
#' @return List with CI bounds
#' @export
ci_bernstein <- function(theta_hat, a_hist, b_hist, y_hist,
                         alpha = 0.05,
                         theta_min = -4, theta_max = 4) {

  t <- length(y_hist)
  if (t == 0) {
    return(list(lower = theta_min, upper = theta_max, width = theta_max - theta_min,
                method = "Bernstein"))
  }

  # Observed information used as plug-in variance proxy
  var_score <- observed_info_2pl(theta_hat, a_hist, b_hist)

  # Maximum deviation (for bounded term in Bernstein)
  a_max <- max(a_hist)
  M <- a_max  # Maximum |X_i - E[X_i]|

  # Dissertation precision-study formula:
  # |S_t(θ0)| ≤ √(2 I_t log(2/α)) + a_max log(2/α) / 3
  epsilon <- bernstein_score_bound(var_score, alpha, M)
  info <- var_score

  if (info <= 0) {
    return(list(lower = theta_min, upper = theta_max, width = theta_max - theta_min,
                method = "Bernstein"))
  }

  # CI radius
  radius <- epsilon / info

  lower <- max(theta_hat - radius, theta_min)
  upper <- min(theta_hat + radius, theta_max)

  list(
    lower = lower,
    upper = upper,
    width = upper - lower,
    epsilon = epsilon,
    var_score = var_score,
    info = info,
    method = "Bernstein"
  )
}

# --- Confidence Sequence (Mixture Method) ---

#' Mixture Confidence Sequence
#'
#' @description
#' Anytime-valid confidence sequence using the dissertation's
#' Laplace-approximated likelihood-ratio / mixture construction.
#'
#' @param theta_hat MLE estimate
#' @param a_hist Vector of administered item discriminations
#' @param b_hist Vector of administered item difficulties
#' @param y_hist Vector of responses
#' @param alpha Confidence level
#' @param rho Legacy interface for the prior variance parameter τ²
#' @param theta_min Lower bound for theta
#' @param theta_max Upper bound for theta
#' @return List with CI bounds
#' @export
ci_mixture <- function(theta_hat, a_hist, b_hist, y_hist,
                       alpha = 0.05, rho = 1,
                       theta_min = -4, theta_max = 4) {

  t <- length(y_hist)
  if (t == 0) {
    return(list(lower = theta_min, upper = theta_max, width = theta_max - theta_min,
                method = "Mixture"))
  }

  # Observed information
  p <- 1 / (1 + exp(-a_hist * (theta_hat - b_hist)))
  info <- sum(a_hist^2 * p * (1 - p))

  if (info <= 0) {
    return(list(lower = theta_min, upper = theta_max, width = theta_max - theta_min,
                method = "Mixture"))
  }

  # Dissertation precision-study formula:
  # r_t = sqrt((2 / I_t) * log(sqrt(1 + I_t * τ²) / α))
  # Keep the legacy argument name `rho` for backward compatibility;
  # it corresponds to τ² in the current dissertation draft.
  tau_sq <- rho
  radius <- lr_mixture_radius(info, alpha, tau_sq)

  lower <- max(theta_hat - radius, theta_min)
  upper <- min(theta_hat + radius, theta_max)

  list(
    lower = lower,
    upper = upper,
    width = upper - lower,
    info = info,
    radius = radius,
    tau_sq = tau_sq,
    method = "Mixture"
  )
}

# --- Wrapper Function ---

#' Compute All CI Methods
#'
#' @description
#' Computes confidence intervals using all available methods.
#'
#' @param theta_hat MLE estimate
#' @param a_hist Vector of administered item discriminations
#' @param b_hist Vector of administered item difficulties
#' @param y_hist Vector of responses
#' @param alpha Confidence level
#' @param theta_min Lower bound for theta
#' @param theta_max Upper bound for theta
#' @return List with CIs from all methods
#' @export
compute_all_cis <- function(theta_hat, a_hist, b_hist, y_hist,
                            alpha = 0.05,
                            theta_min = -4, theta_max = 4) {

  list(
    asymptotic = ci_asymptotic(theta_hat, a_hist, b_hist, alpha),
    azuma = ci_azuma_hoeffding(theta_hat, a_hist, b_hist, y_hist, alpha, theta_min, theta_max),
    bernstein = ci_bernstein(theta_hat, a_hist, b_hist, y_hist, alpha, theta_min, theta_max),
    mixture = ci_mixture(theta_hat, a_hist, b_hist, y_hist, alpha, 1, theta_min, theta_max)
  )
}

#' Check CI Coverage
#'
#' @description
#' Checks if true theta is covered by each CI method.
#'
#' @param theta_true True ability value
#' @param cis List of CIs from compute_all_cis
#' @return Named logical vector indicating coverage
#' @export
check_coverage <- function(theta_true, cis) {
  sapply(cis, function(ci) {
    theta_true >= ci$lower && theta_true <= ci$upper
  })
}

message("=== CI methods loaded ===")
