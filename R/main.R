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

  if(length(X) != length(Y)) {
    stop("X and Y must be the same length.")
  }

  N <- length(X)

  # Data needs to be sorted for the test to work
  data_order <- order(X)
  X <- X[data_order]
  Y <- Y[data_order]

  # Get stat for actual dataset
  t_stat <- get_hall_stat(X, Y, m)

  residuals <-
    calc_residuals_from_estimator(X, Y, bandwidth = bandwidth)

  boot_func <- function(i) {
    # Resample each x and y value
    # We can mix them under the null
    resampled_x_ind <- sample(1:N, size=N, replace=TRUE)
    resampled_residuals_ind <- sample(1:N, size=N, replace=TRUE)

    # Order both
    resampled_x <- X[resampled_x_ind]
    resampled_residuals <- residuals[resampled_residuals_ind]

    resample_order <- order(resampled_x)

    resampled_x <- resampled_x[resample_order]
    resampled_residuals <- resampled_residuals[resample_order]

    t_stat <- get_hall_stat(resampled_x, resampled_residuals, m)
    return(t_stat)
  }

  # Do the bootstrap on multiple cores
  cl <- parallel::makeCluster(ncores)
  parallel::clusterExport(cl, c("X", "residuals", "N", "m", "get_hall_stat"), envir =
                  environment())

  t_vals <-
    unlist(parallel::parLapply(cl, 1:boot_num, boot_func))

  parallel::stopCluster(cl)

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
  return(Y - watson_est(X, Y, preds = X, bandwidth = bandwidth))
}

watson_visual <-
  function(X, Y, bandwidth = bw.nrd(X) * (length(X) ^ -0.1)) {
    # Set up the range for the x-axis
    x_range <-
      seq(min(X), max(X), length.out = 500)  # Fine grid for smooth curve

    # Compute watson estimates
    y_values <- watson_est(X, Y, x_range, bandwidth=bandwidth)

    plot(X, Y,
      main = paste("Nadaraya Watson", "bandwidth =", round(bandwidth, 4)),
      xlab = "X", ylab = "Y", pch = 16)

    # Add the smoothed curve
    lines(x_range, y_values, col = "green", lwd = 5.0)

    return(recordPlot())
  }

