#include "RLSState.h"
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

void RLSState::update(double x, double y) {
  // Build regressor vector
  Eigen::VectorXd x_vec(2);
  x_vec << 1.0, x;

  auto new_state = rls_update_cpp(P, w, x_vec, y);
  P = std::get<0>(new_state);
  w = std::get<1>(new_state);

  // Update window stats
  x_window[curr_n++] = x;
  curr_sum += x;
}

double RLSState::mean() const { return curr_sum / curr_n; }

double RLSState::q_sq() const {
  return std::sqrt((x_window.head(curr_n).array() - mean()).square().sum());
}

double RLSState::t_stat() const { return -w(1, 0) * q_sq(); }
