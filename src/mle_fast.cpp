// =============================================================================
// MLE Fast - C++ Implementation
// =============================================================================
// Description: Fast Newton-Raphson MLE estimation for 2PL IRT model
//              with censoring (bounds enforcement)
// Replaces: utils/mle_functions.R::find_mle_censored
// =============================================================================

#include <Rcpp.h>
#include <cmath>
#include <algorithm>
using namespace Rcpp;

// [[Rcpp::export]]
double find_mle_censored_cpp(NumericVector a_hist, NumericVector b_hist,
                              IntegerVector y_hist,
                              double theta_init = 0.0,
                              double theta_min = -4.0,
                              double theta_max = 4.0,
                              double tol = 1e-6,
                              int max_iter = 50) {
  int n = y_hist.size();

  // Handle empty case
  if (n == 0) return theta_init;

  // Boundary check: all correct or all incorrect
  bool all_one = true, all_zero = true;
  for (int i = 0; i < n; i++) {
    if (y_hist[i] != 1) all_one = false;
    if (y_hist[i] != 0) all_zero = false;
  }
  if (all_one) return theta_max;
  if (all_zero) return theta_min;

  // Check score at boundaries
  double score_min = 0.0, score_max = 0.0;
  for (int i = 0; i < n; i++) {
    double p_min = 1.0 / (1.0 + std::exp(-a_hist[i] * (theta_min - b_hist[i])));
    double p_max = 1.0 / (1.0 + std::exp(-a_hist[i] * (theta_max - b_hist[i])));
    score_min += a_hist[i] * (y_hist[i] - p_min);
    score_max += a_hist[i] * (y_hist[i] - p_max);
  }
  if (score_min <= 0.0) return theta_min;
  if (score_max >= 0.0) return theta_max;

  // Newton-Raphson iteration
  double theta = theta_init;
  for (int iter = 0; iter < max_iter; iter++) {
    double score = 0.0, info = 0.0;
    for (int i = 0; i < n; i++) {
      double p = 1.0 / (1.0 + std::exp(-a_hist[i] * (theta - b_hist[i])));
      score += a_hist[i] * (y_hist[i] - p);
      info += a_hist[i] * a_hist[i] * p * (1.0 - p);
    }

    if (std::abs(score) < tol) break;

    if (info > 1e-8) {
      double theta_new = theta + score / info;
      // Enforce bounds
      theta_new = std::max(theta_min, std::min(theta_max, theta_new));
      if (std::abs(theta_new - theta) < tol) {
        theta = theta_new;
        break;
      }
      theta = theta_new;
    } else {
      // Information too small, return current estimate
      break;
    }
  }

  return std::max(theta_min, std::min(theta_max, theta));
}
