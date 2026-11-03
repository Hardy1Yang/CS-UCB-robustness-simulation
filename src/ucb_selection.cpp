// =============================================================================
// UCB Selection - C++ Implementation
// =============================================================================
// Description: Fast CS-UCB item selection via sup_{theta in C_t} I_i(theta)
// Used by: utils/robustness_functions.R (CS-UCB variants)
// =============================================================================

#include <Rcpp.h>
#include <cmath>
#include <algorithm>
using namespace Rcpp;

// [[Rcpp::export]]
int ucb_select_cpp(NumericVector a, NumericVector b,
                   IntegerVector available,
                   double current_theta, double r_t) {
  // For each available item, compute sup_{theta in [theta-r, theta+r]} I_i(theta)
  // where I_i(theta) = a_i^2 * p_i(theta) * (1 - p_i(theta))
  //
  // Maximum of I_i(theta) occurs at theta = b_i (where p = 0.5, info = a^2/4)
  // If b_i is within [theta-r, theta+r], max = a_i^2/4
  // Otherwise, evaluate at the closer boundary
  //
  // available: 1-indexed R item IDs
  // Returns: 1-indexed R item ID of selected item

  int n_avail = available.size();
  if (n_avail == 0) return NA_INTEGER;

  double best_val = -1.0;
  int best_idx = available[0];

  double lo = current_theta - r_t;
  double hi = current_theta + r_t;

  for (int k = 0; k < n_avail; k++) {
    int i = available[k] - 1;  // Convert to 0-indexed for C++ array access
    double val;

    if (b[i] >= lo && b[i] <= hi) {
      // b_i is within confidence set: info maximized at b_i
      val = a[i] * a[i] / 4.0;
    } else {
      // Evaluate at both boundaries, take maximum
      double p_lo = 1.0 / (1.0 + std::exp(-a[i] * (lo - b[i])));
      double p_hi = 1.0 / (1.0 + std::exp(-a[i] * (hi - b[i])));
      double info_lo = a[i] * a[i] * p_lo * (1.0 - p_lo);
      double info_hi = a[i] * a[i] * p_hi * (1.0 - p_hi);
      val = std::max(info_lo, info_hi);
    }

    if (val > best_val) {
      best_val = val;
      best_idx = available[k];  // Return 1-indexed
    }
  }

  return best_idx;
}
