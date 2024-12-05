#' Do the monotonicity test
#'
#' Does the monotonicity test as described in (Hall and Heckman)
#'
#'
#' @param X Vector of x values
#' @param Y Vector of Y values
#' @param bandwidth Kernel bandwidth
#' @param boot_num Number of bootstrap samples
#' @param m M parameter described in (paper)
#' @param ncores Number of cores to run on
#' @return List with fields p, dist, stat
#' @export
monotonicity_test <- function(X,
                      Y,
                      bandwidth = bw.nrd(X) * (length(X) ^ -0.1),
                      boot_num = 200,
                      m = 10,
                      ncores = 1) {
  N <- length(X)

  # make sure data is sorted
  data_order <- order(X)
  X <- X[data_order]
  Y <- Y[data_order]

  t_stat <- get_t_from_data_cpp(X, Y, m)

  residuals <-
    calc_residuals_from_estimator(X, Y, bandwidth = bandwidth)

  boot_func <- function(i) {
    resampled_x_ind <- sample(1:N, size=N, replace=TRUE)
    resampled_y_ind <- sample(1:N, size=N, replace=TRUE)

    resampled_x <- X[resampled_x_ind]
    resampled_y <- X[resampled_x_ind]

    resample_order <- order(resampled_x)

    resampled_x <- resampled_x[resample_order]
    resampled_y <- resampled_y[resample_order]

    t_stat <- get_t_from_data_cpp(resampled_x, resampled_y, m)
    return(t_stat)
  }

  t_vals <-
    unlist(parallel::mclapply(1:boot_num, boot_func, mc.cores = ncores))

  p_val <- sum(t_vals >= t_stat) / length(t_vals)
  return(list(p = p_val, dist = t_vals, stat = t_stat))
}


# Get estimates with Nadaraya Watson kernel regression
watson_est <- function(X, Y, preds, bandwidth) {
  kernel_weights_x <-
    sapply(preds, function(x_i)
      dnorm((X - x_i) / bandwidth) / bandwidth)
  kernel_weights_y <- kernel_weights_x * Y

  return(sapply(1:length(preds), function(i)
    sum(
      kernel_weights_y[, i] / sum(kernel_weights_x[, i])
    )))
}

# Calculate the residuals after doing kernel regression
calc_residuals_from_estimator <- function(X, Y, bandwidth) {
  return(X - watson_est(X, Y, preds = X, bandwidth = bandwidth))
}

