#include <Rcpp.h>
#include <RcppEigen.h>
#include <tuple>

using namespace Rcpp;
using namespace Eigen;

// Update P and w matrices using recursive least squares
// Update for a single x, y pair
// https://osquant.com/papers/recursive-least-squares-linear-regression/#:~:text=Rather%20than%20recalculating%20a%20least,greater%20emphasis%20on%20recent%20data.
std::tuple<MatrixXd, MatrixXd> rls_update_cpp(const MatrixXd &P,
                                              const MatrixXd &w,
                                              const VectorXd &x, double y) {
  // Calculate numerator and denominator
  MatrixXd Px = P * x; // Precompute since we use it multiple times
  MatrixXd xT_P = x.transpose() * P;

  double denom = 1 + (xT_P * x)(0, 0);
  double final_term = y - (x.transpose() * w)(0, 0);

  auto w_update = w + ((Px / denom) * final_term);
  auto P_update = P - ((Px * xT_P) / denom);

  return std::make_tuple(P_update, w_update);
}

// Calculate the hall statistic given a X and Y vector
// Uses recursive least squares instead of ordinary least squares at each step
// Also returns the "critical interval" where the t-statistic is the largest
// [[Rcpp::export]]
List get_hall_stat(const Eigen::VectorXd &x, const Eigen::VectorXd &y, int m) {
  int n = x.size();

  // If the length of the data is less than window size m, return NA
  if (n - m <= 0) {
    stop("m parameter must be less than the length of the data");
  }

  // initialize max_t_m and the critical interval
  auto max_t_m = -INFINITY;
  auto interval = std::make_tuple(-1, -1);

  // Start at each point
  for (int r = 0; r <= n - m; r++) {
    // Initialize P and w
    MatrixXd P = MatrixXd::Identity(2, 2) * 1e9;

    // W is a 2x1 matrix where row 1 is the intercept and row 2 is the slope
    MatrixXd w = MatrixXd::Zero(2, 1);

    VectorXd x_window(n - r);
    double curr_sum = 0.0;
    int curr_n = 0;
    int curr_ind = 0;

    // Calculate P and w for the first m points
    for (int ind = r; ind < r + m; ind++) {
      VectorXd x_vec(2); // Create a column vector of size 2
      x_vec(0) = 1.0;    // First index is 1
      x_vec(1) = x[ind]; // Second index is the value x[ind]

      auto new_state = rls_update_cpp(P, w, x_vec, y[ind]);
      P = std::get<0>(new_state);
      w = std::get<1>(new_state);

      x_window[curr_ind++] = x[ind];
      curr_sum += x[ind];
      curr_n++;
    }

    double curr_mean = curr_sum / curr_n;

    // Calculate P and w for the remaining points
    for (int s = r + m; s < n; s++) {
      VectorXd x_vec(2);
      x_vec(0) = 1.0;
      x_vec(1) = x[s];

      auto new_state = rls_update_cpp(P, w, x_vec, y[s]);

      P = std::get<0>(new_state);
      w = std::get<1>(new_state);

      x_window[curr_ind++] = x[s];
      curr_sum += x[s];
      curr_n++;
      curr_mean = curr_sum / curr_n;

      double q_sq =
          std::sqrt((x_window.head(curr_n).array() - curr_mean).square().sum());

      // Calculate the t-statistic
      double t_stat = -w(1, 0) * q_sq;

      // Update max_t_m and the critical interval if necessary
      if (t_stat > max_t_m) {
        max_t_m = t_stat;
        interval = std::make_tuple(r, s);
      }
    }
  }

  return List::create(Named("max_t_m") = max_t_m,
                      Named("interval") = NumericVector::create(
                          std::get<0>(interval), std::get<1>(interval)));
}
