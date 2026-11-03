# =============================================================================
# package_wrappers.R - catR Package Wrapper Functions
# =============================================================================
# Description: Provides wrapper functions for catR package to standardize
#              item selection and ability estimation across the simulation.
#
# Strategy: Use catR for standard methods (MFI, MPWI, KL, Progressive)
#           while preserving custom implementations for novel contributions
#           (CS-UCB, confidence sequences, error model).
#
# Reference: Magis & Raiche (2012). catR: An R Package for Computerized
#            Adaptive Testing. Applied Psychological Measurement, 36(7).
# =============================================================================

#' Convert Item Parameters to catR Bank Format
#'
#' @description
#' Converts separate a and b vectors to catR's 4-column matrix format.
#' For 2PL model: columns are (a, b, c=0, d=1).
#'
#' @param a Vector of discrimination parameters
#' @param b Vector of difficulty parameters
#' @return Matrix with 4 columns (a, b, c, d) suitable for catR functions
#' @export
to_catR_bank <- function(a, b) {
  cbind(a = a, b = b, c = rep(0, length(a)), d = rep(1, length(a)))
}

#' Select Next Item Using catR
#'
#' @description
#' Wrapper for catR::nextItem() that handles format conversion and
#' available item filtering. Supports MFI, MPWI, KL, KLP, progressive,
#' and other criteria implemented in catR.
#'
#' The available item pool is specified directly; internally this is
#' converted to catR's 'out' parameter (items to exclude).
#'
#' @param itemBank_catR Item bank in catR format (4-column matrix from to_catR_bank)
#' @param theta Current ability estimate
#' @param available Vector of available item indices (1-based)
#' @param K Total number of items in the bank
#' @param criterion Selection criterion: "MFI", "MPWI", "KL", "KLP",
#'        "progressive", "MLWI", "MEI", "MEPV", "random", etc.
#' @param method Ability estimation method for posterior-based criteria (default "BM")
#' @param priorDist Prior distribution for Bayesian criteria (default "norm")
#' @param priorPar Prior parameters (default c(0, 1))
#' @param D Scaling constant (default 1.0)
#' @param range Theta range for integration (default c(-4, 4))
#' @param nAvailable Number of available items (for progressive criterion)
#' @param randomesque Number of top items to randomly choose from (default 1)
#' @return Index of selected item (1-based, referring to original bank)
#' @export
catR_select <- function(itemBank_catR, theta, available, K,
                        criterion = "MFI",
                        method = "BM",
                        priorDist = "norm",
                        priorPar = c(0, 1),
                        D = 1,
                        range = c(-4, 4),
                        nAvailable = NULL,
                        randomesque = 1) {

  if (length(available) == 0) return(NA)
  if (length(available) == 1) return(available)

  # Convert available items to excluded items (catR uses 'out' to exclude)
  out_items <- setdiff(1:K, available)
  if (length(out_items) == 0) out_items <- NULL

  # For progressive criterion, set nAvailable if not provided
  if (criterion == "progressive" && is.null(nAvailable)) {
    nAvailable <- length(available)
  }

  result <- tryCatch({
    catR::nextItem(
      itemBank = itemBank_catR,
      model = NULL,
      theta = theta,
      out = out_items,
      criterion = criterion,
      method = method,
      priorDist = priorDist,
      priorPar = priorPar,
      D = D,
      range = range,
      nAvailable = nAvailable,
      randomesque = randomesque
    )
  }, error = function(e) {
    warning("catR::nextItem failed (", criterion, "): ", e$message,
            ". Falling back to MFI.")
    # Fallback to simple MFI using available items
    a <- itemBank_catR[available, 1]
    b <- itemBank_catR[available, 2]
    p <- 1 / (1 + exp(-a * (theta - b)))
    infos <- a^2 * p * (1 - p)
    list(item = available[which.max(infos)])
  })

  return(result$item)
}

message("=== catR package wrappers loaded ===")
