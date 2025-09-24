#include "RLSState.h"
#include <Rcpp.h>
#include <RcppEigen.h>
#include <map>
#include <tuple>

using namespace Rcpp;
using namespace Eigen;

struct IntervalStat {
  double t_value;
  std::pair<int, int> interval;
};

// Calculate the hall statistic given a X and Y vector
// Uses recursive least squares instead of ordinary least squares at each step
// Also returns the "critical interval" where the t-statistic is the largest
// [[Rcpp::export]]
List get_hall_stat(const Eigen::VectorXd &x, const Eigen::VectorXd &y,
                   const Eigen::VectorXd &m) {
  int n = x.size();

  int min_window = m(0);

  // If the length of the data is less than window size m, return NA
  if (n - min_window <= 0) {
    stop("m parameter must be less than the length of the data");
  }

  std::map<int, IntervalStat> t_m_map;

  for (int i = 0; i < m.size(); i++) {
    t_m_map[m(i)] = IntervalStat{-INFINITY, std::make_pair(-1, -1)};
  }

  // Start at each point
  for (int r = 0; r <= n - min_window; r++) {

    RLSState state(n - r);

    // Calculate P and w for the first m points
    for (int ind = r; ind < r + min_window; ind++) {
      state.update(x[ind], y[ind]);
    }

    // Calculate P and w for the remaining points
    for (int s = r + min_window; s < n; s++) {
      state.update(x[s], y[s]);
      double t_stat = state.t_stat();

      for (auto &[window_size, stat] : t_m_map) {
        if (window_size <= s - r && t_stat > stat.t_value) {
          stat.t_value = t_stat;
          stat.interval = {r, s};
        }
      }
    }
  }

  NumericVector t_values(m.size());
  List intervals(m.size());

  for (int i = 0; i < m.size(); i++) {
    int window = m(i);
    auto &stat = t_m_map[window];
    t_values[i] = stat.t_value;
    intervals[i] =
        IntegerVector::create(stat.interval.first, stat.interval.second);
  }

  return List::create(Named("t_values") = t_values,
                      Named("intervals") = intervals);
}
