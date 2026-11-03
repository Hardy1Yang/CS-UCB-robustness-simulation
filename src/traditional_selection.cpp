// =============================================================================
// Traditional Item Selection Methods - C++ Implementation
// =============================================================================
// Description: Fast C++ implementations of MPWI, KL, and Progressive item
//              selection methods, replacing catR::nextItem() for 2PL models.
//
// Mathematical Reference:
//   - MPWI: van der Linden (1998), Veerkamp & Berger (1997)
//   - KL:   Chang & Ying (1996)
//   - Progressive: Barrada et al. (2008, 2010)
//
// Integration: Uses same trapezoidal rule as catR:::integrate.catR():
//   sum((y[i] + y[i+1]) / 2 * dx) for equally spaced quadrature points
//
// Note: All item indices are 1-based (R convention)
// =============================================================================

#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <limits>
using namespace Rcpp;

// --- Constants ---
static const int N_QUAD = 33;           // Number of quadrature points (catR default)
static const double QUAD_LOWER = -4.0;
static const double QUAD_UPPER = 4.0;
static const double PROB_EPS = 1e-10;   // Clamp for log(0) prevention

// --- Helper: 2PL probability ---
static inline double prob_2pl(double theta, double a, double b) {
  double z = a * (theta - b);
  // Numerical stability for extreme z
  if (z > 35.0) return 1.0 - PROB_EPS;
  if (z < -35.0) return PROB_EPS;
  return 1.0 / (1.0 + std::exp(-z));
}

// --- Helper: Compute log-likelihood at quadrature points ---
// Returns log L(theta_q | x1..xk) for each quadrature point
// Uses log-sum-exp trick to avoid underflow
static std::vector<double> compute_likelihood_weights(
    NumericVector a_hist, NumericVector b_hist, IntegerVector y_hist,
    const std::vector<double>& quad_points, int n_admin) {

  int nq = quad_points.size();
  std::vector<double> log_lik(nq, 0.0);

  // Compute log-likelihood at each quadrature point
  for (int q = 0; q < nq; q++) {
    double theta_q = quad_points[q];
    double ll = 0.0;
    for (int j = 0; j < n_admin; j++) {
      double p = prob_2pl(theta_q, a_hist[j], b_hist[j]);
      p = std::max(PROB_EPS, std::min(1.0 - PROB_EPS, p));
      if (y_hist[j] == 1) {
        ll += std::log(p);
      } else {
        ll += std::log(1.0 - p);
      }
    }
    log_lik[q] = ll;
  }

  // Log-sum-exp: subtract max for numerical stability
  double max_ll = *std::max_element(log_lik.begin(), log_lik.end());
  std::vector<double> weights(nq);
  for (int q = 0; q < nq; q++) {
    weights[q] = std::exp(log_lik[q] - max_ll);
  }

  return weights;
}

// --- Helper: catR-style trapezoidal integration ---
// Replicates catR:::integrate.catR(x, y):
//   hauteur = x[2:n] - x[1:(n-1)]
//   base = mean of adjacent y pairs
//   sum(base * hauteur)
static double integrate_catR(const std::vector<double>& x,
                             const std::vector<double>& y) {
  double result = 0.0;
  for (size_t i = 0; i < x.size() - 1; i++) {
    double dx = x[i + 1] - x[i];
    double base = (y[i] + y[i + 1]) / 2.0;
    result += base * dx;
  }
  return result;
}

// =============================================================================
// MPWI: Maximum Posterior Weighted Information
// =============================================================================
// MPWI_i = ∫ I_i(θ) × L(θ|x) × π(θ) dθ
// where I_i(θ) = a_i² × p_i(θ) × (1 - p_i(θ))  [2PL Fisher info]
//       π(θ) = φ(θ; 0, 1)                         [N(0,1) prior]
// =============================================================================

// [[Rcpp::export]]
int mpwi_select_cpp(NumericVector a, NumericVector b,
                    IntegerVector available,
                    NumericVector a_hist, NumericVector b_hist,
                    IntegerVector y_hist, int n_admin) {

  int n_avail = available.size();
  if (n_avail == 0) return NA_INTEGER;
  if (n_avail == 1) return available[0];

  // Setup quadrature points
  std::vector<double> quad(N_QUAD);
  double dx = (QUAD_UPPER - QUAD_LOWER) / (N_QUAD - 1);
  for (int q = 0; q < N_QUAD; q++) {
    quad[q] = QUAD_LOWER + q * dx;
  }

  // Precompute likelihood weights (shared across all candidate items)
  std::vector<double> lik_weights = compute_likelihood_weights(
    a_hist, b_hist, y_hist, quad, n_admin);

  // Precompute prior density × likelihood (MPWI uses both)
  // π(θ) = 1/√(2π) × exp(-θ²/2)
  std::vector<double> weighted(N_QUAD);
  for (int q = 0; q < N_QUAD; q++) {
    double prior = std::exp(-0.5 * quad[q] * quad[q]) / std::sqrt(2.0 * M_PI);
    weighted[q] = lik_weights[q] * prior;
  }

  // Score each available item
  double best_score = -std::numeric_limits<double>::infinity();
  int best_item = available[0];

  for (int k = 0; k < n_avail; k++) {
    int idx = available[k] - 1;  // Convert to 0-indexed

    // Compute integrand: I_i(θ_q) × L(θ_q) × π(θ_q)
    std::vector<double> integrand(N_QUAD);
    for (int q = 0; q < N_QUAD; q++) {
      double p = prob_2pl(quad[q], a[idx], b[idx]);
      double info = a[idx] * a[idx] * p * (1.0 - p);
      integrand[q] = info * weighted[q];
    }

    double score = integrate_catR(quad, integrand);

    if (score > best_score) {
      best_score = score;
      best_item = available[k];  // Return 1-indexed
    }
  }

  return best_item;
}

// =============================================================================
// KL: Kullback-Leibler Information
// =============================================================================
// KL_i = ∫ KL_i(θ || θ̂) × L(θ|x) dθ
// where KL_i(θ||θ̂) = p_hat × log(p_hat/p_q) + (1-p_hat) × log((1-p_hat)/(1-p_q))
// Note: KL (not KLP) — no prior weighting
// =============================================================================

// [[Rcpp::export]]
int kl_select_cpp(NumericVector a, NumericVector b,
                  IntegerVector available,
                  double current_theta,
                  NumericVector a_hist, NumericVector b_hist,
                  IntegerVector y_hist, int n_admin) {

  int n_avail = available.size();
  if (n_avail == 0) return NA_INTEGER;
  if (n_avail == 1) return available[0];

  // Setup quadrature points
  std::vector<double> quad(N_QUAD);
  double dx_val = (QUAD_UPPER - QUAD_LOWER) / (N_QUAD - 1);
  for (int q = 0; q < N_QUAD; q++) {
    quad[q] = QUAD_LOWER + q * dx_val;
  }

  // Precompute likelihood weights (shared across all candidate items)
  std::vector<double> lik_weights = compute_likelihood_weights(
    a_hist, b_hist, y_hist, quad, n_admin);

  // Score each available item
  double best_score = -std::numeric_limits<double>::infinity();
  int best_item = available[0];

  for (int k = 0; k < n_avail; k++) {
    int idx = available[k] - 1;  // Convert to 0-indexed

    // p_hat: probability at current theta estimate (constant across quad points)
    double p_hat = prob_2pl(current_theta, a[idx], b[idx]);
    p_hat = std::max(PROB_EPS, std::min(1.0 - PROB_EPS, p_hat));

    // Compute integrand: KL_i(θ_q || θ̂) × L(θ_q)
    std::vector<double> integrand(N_QUAD);
    for (int q = 0; q < N_QUAD; q++) {
      double p_q = prob_2pl(quad[q], a[idx], b[idx]);
      p_q = std::max(PROB_EPS, std::min(1.0 - PROB_EPS, p_q));

      double kl = p_hat * std::log(p_hat / p_q) +
                  (1.0 - p_hat) * std::log((1.0 - p_hat) / (1.0 - p_q));

      integrand[q] = kl * lik_weights[q];
    }

    double score = integrate_catR(quad, integrand);

    if (score > best_score) {
      best_score = score;
      best_item = available[k];  // Return 1-indexed
    }
  }

  return best_item;
}

// =============================================================================
// Progressive: Information-weighted random selection
// =============================================================================
// score_i = I_i(θ̂) × w_t + U_i × (1 - w_t)
// where w_t = Σ_{s=1}^{t} s^AP / Σ_{s=1}^{T-1} s^AP, AP = 1
//       U_i ~ Uniform(0, max_info)
// =============================================================================

// [[Rcpp::export]]
int progressive_select_cpp(NumericVector a, NumericVector b,
                           IntegerVector available,
                           double current_theta,
                           int t, int T_max) {

  int n_avail = available.size();
  if (n_avail == 0) return NA_INTEGER;
  if (n_avail == 1) return available[0];

  // Compute Fisher information for available items
  std::vector<double> info(n_avail);
  double max_info = 0.0;

  for (int k = 0; k < n_avail; k++) {
    int idx = available[k] - 1;  // 0-indexed
    double p = prob_2pl(current_theta, a[idx], b[idx]);
    info[k] = a[idx] * a[idx] * p * (1.0 - p);
    if (info[k] > max_info) max_info = info[k];
  }

  // Compute progressive weight
  // AP = 1: w_t = sum(1:t) / sum(1:(T_max-1)) = t(t+1)/2 / (T_max-1)*T_max/2
  double wq = 0.0;
  if (t > 0 && T_max > 1) {
    double numerator = (double)t * (t + 1) / 2.0;
    double denominator = (double)(T_max - 1) * T_max / 2.0;
    wq = numerator / denominator;
    if (wq > 1.0) wq = 1.0;
  }

  // Generate random values using R's RNG (for reproducibility with set.seed)
  NumericVector rand_vals = Rcpp::runif(n_avail, 0.0, max_info);

  // Compute final scores
  double best_score = -std::numeric_limits<double>::infinity();
  int best_item = available[0];

  for (int k = 0; k < n_avail; k++) {
    double score = info[k] * wq + rand_vals[k] * (1.0 - wq);
    if (score > best_score) {
      best_score = score;
      best_item = available[k];  // Return 1-indexed
    }
  }

  return best_item;
}
