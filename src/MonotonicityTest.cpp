#include "Matrix.h"
#include "Rcpp/vector/instantiation.h"
#include <Rcpp.h>
#include <iostream>
#include <tuple>

using namespace Rcpp;
using namespace matrix_utils;

// Update P and w matrices using recursive least squares
// Update for a single x, y pair
// https://osquant.com/papers/recursive-least-squares-linear-regression/#:~:text=Rather%20than%20recalculating%20a%20least,greater%20emphasis%20on%20recent%20data.
std::tuple<NumericMatrix, NumericMatrix>
rls_update_cpp(NumericMatrix &P, NumericMatrix &w, NumericMatrix &x, double y) {
  // Calculate numerator and denominator
  NumericMatrix numerator = P * x;

  NumericMatrix xT = transpose(x);

  double denom = 1 + as<double>(wrap(xT * P * x));
  double final_term = as<double>(y - xT * w);

  auto w_update = w + ((numerator / denom) * final_term);
  auto P_update = P - ((P * x * xT * P) / denom);

  return std::make_tuple(P_update, w_update);
}

// Calculate the hall statistic given a X and Y vector
// Uses recursive least squares instead of ordinary least squares at each step
// [[Rcpp::export]]
double get_t_from_data_cpp(NumericVector &x, const NumericVector &y, int m) {
  int n = x.length();

  // If the length of the data is less than window size m, return NA
  if (n - m <= 0) {
    stop("m parameter must be less than the length of the data");
  }

  auto max_t_m = -INFINITY;

  // Start at each point
  for (int r = 0; r <= n - m; r++) {
    // Initialize P and w
    NumericMatrix P = NumericMatrix(2);
    P.fill_diag(1e9);

    // W is a 2x1 matrix where row 1 is the intercept and row 2 is the slope
    NumericMatrix w = NumericMatrix(2, 1);

    NumericVector x_window(n - r);
    double curr_sum = 0.0;
    int curr_n = 0;
    int curr_ind = 0;

    // Calculate P and w for the first m points
    for (int ind = r; ind < r + m; ind++) {
      NumericVector x_vec = NumericVector::create(1, x[ind]);
      NumericMatrix x_colvec = create_colvec(x_vec);

      auto new_state = rls_update_cpp(P, w, x_colvec, y[ind]);
      P = std::get<0>(new_state);
      w = std::get<1>(new_state);

      x_window[curr_ind++] = x[ind];
      curr_sum += x[ind];
      curr_n++;
    }

    double curr_mean = curr_sum / curr_n;

    // Calculate P and w for the remaining points
    for (int s = r + m; s < n; s++) {
      NumericVector x_vec = NumericVector::create(1, x[s]);
      NumericMatrix x_colvec = create_colvec(x_vec);

      auto new_state = rls_update_cpp(P, w, x_colvec, y[s]);

      P = std::get<0>(new_state);
      w = std::get<1>(new_state);

      x_window[curr_ind++] = x[s];
      curr_sum += x[s];
      curr_n++;
      curr_mean = curr_sum / curr_n;

      NumericVector diff = x_window - curr_mean;
      double sum_sq = sum(diff * diff);
      double q_sq = sqrt(sum_sq);

      double t_stat = -w(1, 0) * q_sq; // t-statistic

      if (t_stat > max_t_m) {
        max_t_m = t_stat;
      }
    }
  }

  return max_t_m;
}
