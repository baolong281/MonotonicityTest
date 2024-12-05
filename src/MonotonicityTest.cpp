#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]


#include <vector>
using namespace Rcpp;

List rls_update_cpp(arma::mat &P, arma::mat &w, arma::mat &x, float y) {
  arma::mat numerator = P * x;
  arma::mat denom = 1 + (x.t() * P * x);
  arma::mat final_term = (y - x.t() * w);

  arma::mat w_update = w + ((numerator / denom(0)) * final_term);

  arma::mat other_denom = 1 + (x.t() * P * x);
  arma::mat P_update = P - ((P * x * x.t() * P) / (other_denom(0)));
  return List::create(Rcpp::Named("P") = P_update, Rcpp::Named("w") = w_update);
}

// [[Rcpp::export]]
double get_t_from_data_cpp(const arma::vec &x, const arma::vec &y, int m) {
  int n = x.size();

  if (n - m <= 0) {
    return NA_REAL;
  }

  auto max_t_m = -INFINITY;

  return(4.0);

  for (int r = 0; r <= n - m; r++) {
    // Initialize P and w
    arma::mat P = arma::eye<arma::mat>(2, 2) * 1e9;
    arma::colvec w = arma::zeros<arma::colvec>(2);

    arma::vec curr_x;
    double curr_sum = 0.0;
    int curr_n = 0;

    // calculate the initial window
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

