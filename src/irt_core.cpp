// =============================================================================
// IRT Core Functions - C++ Implementation
// =============================================================================
// Description: Fast C++ implementations of IRT core computations
//              (sigmoid, Fisher information, curvature lower bound)
// Used by: utils/irt_functions.R, utils/robustness_functions.R,
//          utils/precision_functions.R
// =============================================================================

#include <Rcpp.h>
#include <cmath>
#include <algorithm>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector sigmoid_cpp(NumericVector z) {
  int n = z.size();
  NumericVector result(n);
  for (int i = 0; i < n; i++) {
    result[i] = 1.0 / (1.0 + std::exp(-z[i]));
  }
  return result;
}

// [[Rcpp::export]]
NumericVector info_vec_cpp(double theta, NumericVector a, NumericVector b) {
  // Vectorized Fisher information I_i(theta) = a_i^2 * p_i * (1-p_i)
  // for all items at a given theta
  int n = a.size();
  NumericVector result(n);
  for (int i = 0; i < n; i++) {
    double p = 1.0 / (1.0 + std::exp(-a[i] * (theta - b[i])));
    result[i] = a[i] * a[i] * p * (1.0 - p);
  }
  return result;
}

// [[Rcpp::export]]
double info_sum_cpp(NumericVector a, NumericVector b, double theta) {
  // Sum of Fisher information for administered items
  // Returns max(sum, 1e-3) for numerical stability
  int n = a.size();
  double total = 0.0;
  for (int i = 0; i < n; i++) {
    double p = 1.0 / (1.0 + std::exp(-a[i] * (theta - b[i])));
    total += a[i] * a[i] * p * (1.0 - p);
  }
  return std::max(total, 1e-3);
}

// [[Rcpp::export]]
double cW_cpp(double a, double W) {
  // Worst-case curvature lower bound (Lemma 2)
  // cW(a, W) = a^2 * sigma(aW) * (1 - sigma(aW))
  double sig = 1.0 / (1.0 + std::exp(-a * W));
  return a * a * sig * (1.0 - sig);
}

// [[Rcpp::export]]
double H_lower_cpp(NumericVector a_used, NumericVector r_history,
                   double diam_Theta, int n_used) {
  // Curvature lower bound sum for csucb_theory
  // H_t^- = sum_{s=1}^{t-1} cW(a_s, r_s + diam_Theta)
  double total = 0.0;
  for (int s = 0; s < n_used; s++) {
    double W_s = r_history[s] + diam_Theta;
    double sig = 1.0 / (1.0 + std::exp(-a_used[s] * W_s));
    total += a_used[s] * a_used[s] * sig * (1.0 - sig);
  }
  return total;
}
