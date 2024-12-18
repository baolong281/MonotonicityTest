#' Perform Monotonicity Test
#'
#' Performs a monotonicity test between the vectors \code{X} and \code{Y} as described in Hall and Heckman (2000).
#' This function uses a bootstrap approach to test for monotonicity in a nonparametric regression setting.
#'
#' The test evaluates the following hypotheses:
#'
#' \strong{\eqn{H_0}}: The regression function is monotonic
#' \itemize{
#'   \item \emph{Non-decreasing} if \code{negative = FALSE}
#'   \item \emph{Non-increasing} if \code{negative = TRUE}
#' }
#'
#' \strong{\eqn{H_A}}: The regression function is not monotonic
#'
#' @param X Numeric vector of predictor variable values. Must not contain missing or infinite values.
#' @param Y Numeric vector of response variable values. Must not contain missing or infinite values.
#' @param bandwidth Numeric value for the kernel bandwidth. Default is calculated as \code{bw.nrd(X) * (length(X) ^ -0.1)}.
#' @param boot_num Integer specifying the number of bootstrap samples. Default is \code{200}.
#' @param m Integer parameter used in the calculation of the test statistic. Corresponds to the minimum window size
#' to calculate the test statistic over. Default is \code{floor(0.05) * length(X)}.
#' @param ncores Integer specifying the number of cores to use for parallel processing. Default is \code{1}.
#' @param negative Logical value indicating whether to test for a monotonic decreasing (negative) relationship. Default is \code{FALSE}.
#' @return A list with the following components:
#' \describe{
#'   \item{\code{p}}{The p-value of the test.}
#'   \item{\code{dist}}{The distribution of test statistics from bootstrap samples.}
#'   \item{\code{stat}}{The test statistic calculated from the original data.}
#' }
#' @note For large datasets (e.g., \eqn{n \geq 6500}) this function may require significant computation time due to having to compute the statistic for every possible interval.
#' Consider reducing \code{boot_num}, using a subset of the data, or using parallel processing with \code{ncores} to improve performance.
#' @references
#'   Hall, P., & Heckman, N. E. (2000). Testing for monotonicity of a regression mean by calibrating for linear functions. \emph{The Annals of Statistics}, \strong{28}(1), 20–39.
#'
#' @examples
#' # Generate sample data
#' X <- runif(100)
#' Y <- X^2 + rnorm(100, sd = 0.1)
#'
#' # Perform monotonicity test
#' result <- monotonicity_test(X, Y, boot_num=1)
#' print(result$p)  # Display the p-value
#' @export
monotonicity_test <- function(X,
                      Y,
                      bandwidth = bw.nrd(X) * (length(X) ^ -0.1),
                      boot_num = 200,
                      m = floor(0.05 * length(X)),
                      ncores = 1,
                      negative=FALSE
                      ) {

  # Check if values are NA, NaN or infinite.
  if (any(!is.finite(X)) || any(!is.finite(Y))) {
    stop("X and Y must contain only finite values (no NA, NaN, or Inf).")
  }

  # Invert the data if testing for negativity
  if(negative) {
    Y <- -Y
  }

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

#' Generate Kernel Plot
#'
#' Creates a scatter plot of the input vectors `X` and `Y`, and overlays a Nadaraya-Watson kernel regression curve using the specified bandwidth.
#'
#' @param X Vector of x values.
#' @param Y Vector of y values.
#' @param bandwidth Kernel bandwidth used for the Nadaraya-Watson estimator. Default is calculated as \code{bw.nrd(X) * (length(X) ^ -0.1)}.
#' @return A recorded plot object containing the scatter plot with the kernel regression curve.
#' @references
#'   Nadaraya, E. A. (1964). On estimating regression. \emph{Theory of Probability and Its Applications}, \strong{9}(1), 141–142.
#'
#'   Watson, G. S. (1964). Smooth estimates of regression functions. \emph{Sankhyā: The Indian Journal of Statistics, Series A}, 359-372.
#' @export
create_kernel_plot <-
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
