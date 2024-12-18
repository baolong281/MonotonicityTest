#include <Rcpp.h>
#include <RcppEigen.h>
#include <tuple>

using namespace Rcpp;
using namespace Eigen;

// Compute b(s) and V(s) for a given s (x, h)
std::tuple<double, double> compute_b_V(const Eigen::VectorXd &Y,
                                       const Eigen::VectorXd &X,
                                       const Eigen::VectorXd &sigma_sq,
                                       double x_s, double h_s, int k) {
  int n = Y.size();
  double b_s = 0.0;
  double V_s = 0.0;

  // Precompute kernel weights for X_i and X_j
  Eigen::VectorXd K_i =
      (((X.array() - x_s) / h_s).abs().array() <= 1.0).cast<double>();

  for (int i = 0; i < n; ++i) {
    // Skip if K_i is zero
    if (K_i[i] == 0.0) {
      continue;
    }

    // Inner sum of variance
    double sum_Q = 0.0;

    for (int j = 0; j < n; ++j) {
      // Skip if K_j is zero
      if (K_i[j] == 0.0) {
        continue;
      }

      // Compute Q(X_i, X_j, s)
      double Q_ij = std::pow(std::abs(X[i] - X[j]), k);
      Q_ij *= K_i[i] * K_i[j];

      // Compute sign(X_j - X_i)
      double X_diff = X[j] - X[i];
      double sign_X = X_diff > 0 ? 1.0 : (X_diff < 0 ? -1.0 : 0.0);

      // Update b(s)
      b_s += (Y[i] - Y[j]) * sign_X * Q_ij;

      // Update sums for V(s)
      sum_Q += sign_X * Q_ij;
    }

    // Variance is the sum of sd^2 and inner sums squared
    V_s += sigma_sq[i] * sum_Q * sum_Q;
  }

  // b_s has a factor of 1/2
  b_s = 0.5 * b_s;

  return std::make_tuple(b_s, V_s);
}

// [[Rcpp::export]]
double adaptive_test_stat(const Eigen::VectorXd &X, const Eigen::VectorXd &Y,
                          const Eigen::VectorXd &sigma_sq,
                          const Eigen::VectorXd &x_vals,
                          const Eigen::VectorXd &h_vals, int k) {
  int num_x = x_vals.size();
  int num_h = h_vals.size();

  double T_max = -INFINITY;

  // Loop over all combinations of x and h in S_n
  for (int i = 0; i < num_x; ++i) {
    double x_s = x_vals[i];
    for (int j = 0; j < num_h; ++j) {
      double h_s = h_vals[j];

      // Compute b(s) and V(s)
      auto b_V = compute_b_V(Y, X, sigma_sq, x_s, h_s, k);
      double b_s = std::get<0>(b_V);
      double V_s = std::get<1>(b_V);

      // Avoid division by zero
      if (V_s <= 0.0) {
        continue;
      }

      double T_s = b_s / std::sqrt(V_s);

      // Update T_max
      if (T_s > T_max) {
        T_max = T_s;
      }
    }
  }

  return T_max;
}

// Estimate global variance using the Rice estimator
Eigen::VectorXd estimate_sigma_sq(const Eigen::VectorXd &Y) {
  int n = Y.size();
  Eigen::VectorXd sigma_sq(n);

  // Compute the global variance estimate
  double sigma2_hat = 0.0;
  for (int i = 0; i < n - 1; ++i) {
    double diff = Y[i + 1] - Y[i];
    sigma2_hat += diff * diff;
  }
  sigma2_hat = sigma2_hat / (2.0 * (n - 1));

  // Assign to all sigma_i^2
  sigma_sq.fill(sigma2_hat);

  return sigma_sq;
}

// [[Rcpp::export]]
double perform_adaptive_test(const Eigen::VectorXd &X, const Eigen::VectorXd &Y,
                             const Eigen::VectorXd &x_grid,
                             const Eigen::VectorXd &h_grid, int k) {
  Eigen::VectorXd sigma_sq = estimate_sigma_sq(Y);
  double T_stat = adaptive_test_stat(X, Y, sigma_sq, x_grid, h_grid, k);

  return T_stat;
}
