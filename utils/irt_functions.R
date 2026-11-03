# =============================================================================
# IRT Model Functions
# =============================================================================
# Description: Core functions for Item Response Theory (IRT) models,
#              specifically the 2-Parameter Logistic (2PL) model.
# Source: Migrated from OCO_CAT/finding/utils/irt_functions.R
# =============================================================================

#' Sigmoid Function
#'
#' @param z Numeric vector of input values
#' @return Numeric vector of sigmoid values between 0 and 1
#' @export
#' @examples
#' sigmoid(0)  # Returns 0.5
#' sigmoid(c(-2, 0, 2))
sigmoid <- function(z) {
  1 / (1 + exp(-z))
}

#' 2PL Probability Function
#'
#' @param theta Ability parameter (can be vector)
#' @param a Discrimination parameter (can be vector)
#' @param b Difficulty parameter (can be vector)
#' @return Probability of correct response
#' @export
#' @examples
#' p_fun(theta = 0, a = 1, b = 0)  # Returns 0.5
p_fun <- function(theta, a, b) {
  sigmoid(a * (theta - b))
}

#' Fisher Information Function for 2PL Model
#'
#' @param theta Ability parameter
#' @param a Discrimination parameter
#' @param b Difficulty parameter
#' @return Fisher information value
#' @export
#' @examples
#' info_fun(theta = 0, a = 1, b = 0)
info_fun <- function(theta, a, b) {
  p <- p_fun(theta, a, b)
  (a^2) * p * (1 - p)
}

#' Clipped Weight Function cW (Worst-Case Logistic Variance)
#'
#' @description
#' Computes the worst-case Fisher information lower bound over
#' a bounded region of width W around the item difficulty.
#' Used in CS-UCB curvature lower bound H_t^-.
#' From Lemma 2: c_W(a) = a^2 * sigma(aW) * (1 - sigma(aW))
#'
#' @param a Discrimination parameter
#' @param W Width parameter (confidence radius + window)
#' @return Clipped weight value (worst-case information)
#' @export
cW <- function(a, W) {
  sigmoid_val <- sigmoid(a * W)
  a^2 * sigmoid_val * (1 - sigmoid_val)
}

message("=== IRT functions loaded ===")
