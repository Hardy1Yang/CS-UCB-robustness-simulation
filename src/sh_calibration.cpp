// =============================================================================
// Sympson-Hetter Calibration - C++ Implementation
// =============================================================================
// Description: Fast C++ implementation of Sympson-Hetter (1985) probabilistic
//              exposure-control parameter calibration. Inline MFI selection.
//              Default theta sampler: truncated normal TN(0, 1.5, -3, 3) — see
//              Dissertation-defense/analysis/01_sympson_hetter.md §6 Option 2.
//
// Replaces (when backend = "auto" / "cpp"):
//   utils/exposure_control.R::calibrate_sympson_hetter()  inner double loop
//
// Reference:
//   - Sympson & Hetter (1985), Proceedings of the 27th MTA
//   - Georgiadou, Triantafillou & Economides (2007), JTLA 5(8)
//   - Stocking & Lewis (1995a, 2000) for known design ceiling
//
// RNG: Uses R's RNG stream (Rcpp::runif / R::norm_rand), so a single set.seed()
// in R fully controls reproducibility. Bit-identical with the pure-R version is
// not guaranteed (the R version consumes the stream in a slightly different
// order), but statistical properties + final convergence outcome match.
// =============================================================================

#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <vector>
using namespace Rcpp;

namespace {

constexpr double PROB_EPS = 1e-10;

inline double prob_2pl(double theta, double a, double b) {
  const double z = a * (theta - b);
  if (z > 35.0) return 1.0 - PROB_EPS;
  if (z < -35.0) return PROB_EPS;
  return 1.0 / (1.0 + std::exp(-z));
}

// Sample one truncated-normal value via rejection. lo/hi must satisfy lo < hi.
inline double sample_tn(double mean, double sd, double lo, double hi) {
  for (int tries = 0; tries < 1000; ++tries) {
    const double x = mean + sd * R::norm_rand();
    if (x >= lo && x <= hi) return x;
  }
  // Pathological fall-through: clip a fresh normal draw
  return std::max(lo, std::min(hi, mean + sd * R::norm_rand()));
}

inline double sample_theta(int dist_code, double mean, double sd,
                           double lo, double hi) {
  // dist_code: 0 = truncated_normal (default), 1 = normal+clip, 2 = uniform
  if (dist_code == 1) {
    const double x = mean + sd * R::norm_rand();
    return std::max(lo, std::min(hi, x));
  }
  if (dist_code == 2) {
    return lo + (hi - lo) * R::unif_rand();
  }
  return sample_tn(mean, sd, lo, hi);
}

}  // namespace

// [[Rcpp::export]]
List sh_calibrate_cpp(NumericVector a, NumericVector b,
                      double target_rate,
                      int n_calibration,
                      int n_examinees_per_iter,
                      int T_max,
                      int dist_code,
                      double dist_mean,
                      double dist_sd,
                      double dist_lo,
                      double dist_hi,
                      bool verbose) {
  const int n_items = a.size();
  if (b.size() != n_items) stop("sh_calibrate_cpp: a and b length mismatch");
  if (T_max <= 0 || T_max > n_items) stop("sh_calibrate_cpp: invalid T_max");
  if (target_rate <= 0.0 || target_rate >= 1.0) {
    stop("sh_calibrate_cpp: target_rate must be in (0, 1)");
  }

  std::vector<double> control_params(n_items, 1.0);
  std::vector<int> admin_counts(n_items, 0);
  std::vector<int> selection_counts(n_items, 0);
  std::vector<char> is_administered(n_items, 0);  // per-session mask

  double max_rate_achieved = NA_REAL;
  bool converged = false;
  int iter_used = 0;

  for (int iter = 0; iter < n_calibration; ++iter) {
    iter_used = iter + 1;
    if (verbose) {
      Rcpp::Rcout << "Calibration iteration " << iter_used
                  << " of " << n_calibration << "\n";
    }
    std::fill(selection_counts.begin(), selection_counts.end(), 0);
    std::fill(admin_counts.begin(), admin_counts.end(), 0);

    for (int e = 0; e < n_examinees_per_iter; ++e) {
      const double theta = sample_theta(dist_code, dist_mean, dist_sd,
                                         dist_lo, dist_hi);
      // Reset per-session mask
      std::fill(is_administered.begin(), is_administered.end(), 0);
      int n_admin = 0;

      for (int t = 0; t < T_max; ++t) {
        // MFI selection over not-yet-administered items
        double best_info = -1.0;
        int selected = -1;
        for (int i = 0; i < n_items; ++i) {
          if (is_administered[i]) continue;
          const double p = prob_2pl(theta, a[i], b[i]);
          const double info_i = a[i] * a[i] * p * (1.0 - p);
          if (info_i > best_info) {
            best_info = info_i;
            selected = i;
          }
        }
        if (selected < 0) break;  // no items left

        selection_counts[selected] += 1;

        // SH gate
        const double u = R::unif_rand();
        int final_item = selected;
        if (u > control_params[selected]) {
          // Reject — pick uniform random alternative among remaining
          // (matches R's `sample(setdiff(available, selected), 1)`)
          int n_alt = 0;
          for (int i = 0; i < n_items; ++i) {
            if (!is_administered[i] && i != selected) ++n_alt;
          }
          if (n_alt > 0) {
            // Pick k-th alternative where k ~ Uniform{0, ..., n_alt-1}
            const int k = static_cast<int>(R::unif_rand() * n_alt);
            const int k_clamped = std::min(k, n_alt - 1);
            int seen = 0;
            for (int i = 0; i < n_items; ++i) {
              if (!is_administered[i] && i != selected) {
                if (seen == k_clamped) { final_item = i; break; }
                ++seen;
              }
            }
          }
          // else: no alternatives, use selected
        }

        admin_counts[final_item] += 1;
        is_administered[final_item] = 1;
        ++n_admin;
      }
    }

    // Update control parameters (matches R version L249-258)
    const double inv_n = 1.0 / static_cast<double>(n_examinees_per_iter);
    int n_overexposed = 0;
    double iter_max_rate = 0.0;
    for (int i = 0; i < n_items; ++i) {
      const double rate = admin_counts[i] * inv_n;
      if (rate > iter_max_rate) iter_max_rate = rate;
      if (rate > target_rate) ++n_overexposed;
      if (selection_counts[i] > 0) {
        if (rate > target_rate) {
          control_params[i] *= (target_rate / rate);
        } else if (rate < target_rate * 0.9 && control_params[i] < 1.0) {
          control_params[i] = std::min(1.0, control_params[i] * 1.1);
        }
      }
    }
    max_rate_achieved = iter_max_rate;

    if (verbose) {
      Rcpp::Rcout << "  Max exposure rate: "
                  << std::round(max_rate_achieved * 1000.0) / 1000.0
                  << "\n  Items exceeding target: " << n_overexposed << "\n";
    }

    if (max_rate_achieved <= target_rate) {
      if (verbose) {
        Rcpp::Rcout << "Converged at iteration " << iter_used
                    << " (max_rate: " << max_rate_achieved << ")\n";
      }
      converged = true;
      break;
    }
  }

  // Compute final overexposed count (using last iter's admin rates is the
  // R version's contract via `n_overexposed = sum(admin_rates > target)`)
  const double inv_n = 1.0 / static_cast<double>(n_examinees_per_iter);
  int final_overexposed = 0;
  for (int i = 0; i < n_items; ++i) {
    if (admin_counts[i] * inv_n > target_rate) ++final_overexposed;
  }

  return List::create(
    _["control_params"] = NumericVector(control_params.begin(),
                                          control_params.end()),
    _["max_exposure"] = max_rate_achieved,
    _["n_overexposed"] = final_overexposed,
    _["iterations"] = iter_used,
    _["converged"] = converged
  );
}
