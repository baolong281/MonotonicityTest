// RLSState.h
#ifndef RLSSTATE_H
#define RLSSTATE_H

#include <RcppEigen.h>
#include <tuple>

std::tuple<Eigen::MatrixXd, Eigen::MatrixXd>
rls_update_cpp(const Eigen::MatrixXd &P, const Eigen::MatrixXd &w,
               const Eigen::VectorXd &x_vec, double y);

// Struct to hold state of recursive least squares
struct RLSState {
  Eigen::MatrixXd P;        // covariance matrix
  Eigen::MatrixXd w;        // weights (intercept, slope)
  Eigen::VectorXd x_window; // window of x-values
  int curr_n = 0;
  double curr_sum = 0.0;

  RLSState(int max_window)
      : P(Eigen::MatrixXd::Identity(2, 2) * 1e9),
        w(Eigen::MatrixXd::Zero(2, 1)), x_window(max_window) {}

  // Update state with new observation (x, y)
  void update(double x, double y);

  // Mean of x values in the window
  double mean() const;

  // q statistic (scaled deviation from mean)
  double q_sq() const;

  // t statistic (based on slope coefficient and q_sq)
  double t_stat() const;
};

#endif // RLSSTATE_H
