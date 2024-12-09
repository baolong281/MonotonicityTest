#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

#include <vector>
using namespace Rcpp;

// Update P and w matrices using recursive least squares
// Update for a single x, y pair
// https://osquant.com/papers/recursive-least-squares-linear-regression/#:~:text=Rather%20than%20recalculating%20a%20least,greater%20emphasis%20on%20recent%20data.
List rls_update_cpp(arma::mat &P, arma::mat &w, arma::mat &x, float y) {
  arma::mat numerator = P * x;
  arma::mat denom = 1 + (x.t() * P * x);
  arma::mat final_term = (y - x.t() * w);

  arma::mat w_update = w + ((numerator / denom(0)) * final_term);

  arma::mat other_denom = 1 + (x.t() * P * x);
  arma::mat P_update = P - ((P * x * x.t() * P) / (other_denom(0)));
  return List::create(Rcpp::Named("P") = P_update, Rcpp::Named("w") = w_update);
}

// Calculate the hall statistic given a X and Y vector
// Uses recursive least squares instead of ordinary least squares at each step
// [[Rcpp::export]]
double get_t_from_data_cpp(const arma::vec &x, const arma::vec &y, int m) {
  int n = x.size();

  // If the length of the data is less than window size m, return NA
  if (n - m <= 0) {
    return NA_REAL;
  }

  auto max_t_m = -INFINITY;

  // Start at each point
  for (int r = 0; r <= n - m; r++) {
    // Initialize P and w
    arma::mat P = arma::eye<arma::mat>(2, 2) * 1e9;
    arma::colvec w = arma::zeros<arma::colvec>(2);

    arma::vec curr_x;
    double curr_sum = 0.0;
    int curr_n = 0;

    // Calculate P and w for the first m points
    for (int ind = r; ind < r + m; ind++) {
      arma::colvec x_vec = {1, x[ind]}; // Feature vector (1, x)
      List new_state = rls_update_cpp(P, w, x_vec, y[ind]);
      P = Rcpp::as<arma::mat>(new_state["P"]);
      w = Rcpp::as<arma::mat>(new_state["w"]);

      curr_x.insert_rows(curr_x.n_rows, arma::vec{x[ind]});
      curr_sum += x[ind];
      curr_n++;
    }

    double curr_mean = curr_sum / curr_n;

    // Calculate P and w for the remaining points
    for (int s = r + m; s < n; s++) {
      arma::colvec x_vec = {1, x[s]}; // f vector (1, x)
      List new_state = rls_update_cpp(P, w, x_vec, y[s]);
      P = Rcpp::as<arma::mat>(new_state["P"]);
      w = Rcpp::as<arma::mat>(new_state["w"]);

      curr_x.insert_rows(curr_x.n_rows, arma::vec{x[s]});
      curr_sum += x[s];
      curr_n++;
      curr_mean = curr_sum / curr_n;

      double q_sq = std::sqrt(
          arma::accu(arma::square(curr_x - curr_mean))); // sqrt of q stat

      double t_stat = -w(1) * q_sq; // t-statistic

      if (t_stat > max_t_m) {
        max_t_m = t_stat;
      }
    }
  }

  return max_t_m;
}
