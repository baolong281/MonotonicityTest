#' Perform Monotonicity Test
#'
#' Performs a monotonicity test between the vectors \eqn{X} and \eqn{Y}
#' as described in Hall and Heckman (2000).
#' This function uses a bootstrap approach to test for monotonicity
#' in a nonparametric regression setting.
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
#' @param X Numeric vector of predictor variable values.
#'          Must not contain missing or infinite values.
#' @param Y Numeric vector of response variable values.
#'          Must not contain missing or infinite values.
#' @param bandwidth Numeric value for the kernel bandwidth used in the
#'                 Nadaraya-Watson estimator. Default is calculated as
#'                 \code{bw.nrd(X) * (length(X) ^ -0.1)}.
#' @param boot_num Integer specifying the number of bootstrap samples.
#'                Default is \code{200}.
#' @param m Integer parameter used in the calculation of the test statistic.
#'          Corresponds to the minimum window size to calculate the test
#'          statistic over or a "smoothing" parameter. Lower values increase
#'          the sensitivity of the test to local deviations from monotonicity.
#'          Default is \code{floor(0.05 * length(X))}.
#' @param ncores Integer specifying the number of cores to use for parallel
#'              processing. Default is \code{1}.
#' @param negative Logical value indicating whether to test for a monotonic
#'                decreasing (negative) relationship. Default is \code{FALSE}.
#' @param seed Optional integer for setting the random seed. If NULL (default),
#'            the global random state is used.
#' @return A list with the following components:
#' \describe{
#'   \item{\code{p}}{The p-value of the test. A small p-value (e.g., < 0.05)
#'                   suggests evidence against the null hypothesis of
#'                   monotonicity.}
#'   \item{\code{dist}}{The distribution of test statistic under the null from
#'                      bootstrap samples. The length of \code{dist} is equal
#'                      to \code{boot_num}.}
#'   \item{\code{stat}}{The test statistic \eqn{T_m} calculated from the original data.}
#'   \item{\code{plot}}{A ggplot object with a scatter plot where the points of the
#'                      "critical interval" are highlighted. This critical interval
#'                       is the interval where \eqn{T_m} is greatest. }
#'   \item{\code{interval}}{Numeric vector containing the indices of the "critical interval".
#'                          The first index indicates where the interval starts, and
#'                          the second indicates where it ends in the sorted \code{X} vector.}
#' }
#' @note For large datasets (e.g., \eqn{n \geq 6500}) this function may require
#'       significant computation time due to having to compute the statistic
#'       for every possible interval. Consider reducing \code{boot_num}, using
#'       a subset of the data, or using parallel processing with \code{ncores}
#'       to improve performance.
#'
#'       In addition to this, a minimum of 300 observations is recommended for
#'       kernel estimates to be reliable.
#' @references
#'   Hall, P., & Heckman, N. E. (2000). Testing for monotonicity of a regression
#'   mean by calibrating for linear functions. \emph{The Annals of Statistics},
#'   \strong{28}(1), 20–39.
#'
#' @examples
#' # Example 1: Usage on monotonic increasing function
#' # Generate sample data
#' seed <- 42
#' set.seed(seed)
#'
#' X <- runif(500)
#' Y <- 4 * X + rnorm(500, sd = 1)
#' result <- monotonicity_test(X, Y, boot_num = 25, seed = seed)
#'
#' print(result$p)
#' print(result$stat)
#' print(result$dist)
#' print(result$interval)
#'
#' # Example 2: Usage on non-monotonic function
#' seed <- 42
#' set.seed(seed)
#'
#' X <- runif(500)
#' Y <- (X - 0.5) ^ 2 + rnorm(500, sd = 0.5)
#' result <- monotonicity_test(X, Y, boot_num = 25, seed = seed)
#'
#' print(result$p)
#' print(result$stat)
#' print(result$dist)
#' print(result$interval)
#'
#' @export
monotonicity_test <-
  function(X,
           Y,
           bandwidth = bw.nrd(X) * (length(X) ^ -0.1),
           boot_num = 200,
           m = floor(0.05 * length(X)),
           ncores = 1,
           negative = FALSE,
           seed = NULL) {
    validate_inputs(X, Y, m=m)
    N <- length(X)

    # Invert the data if testing for negativity
    if (negative) {
      Y <- -Y
    }

    # Data needs to be sorted for the test to work
    data_order <- order(X)
    X <- X[data_order]
    Y <- Y[data_order]

    # Get stat for actual dataset
    original_result <- get_hall_stat(X, Y, m)
    t_stat <- original_result$max_t_m
    crit_interval <- original_result$interval

    residuals <-
      calc_residuals_from_estimator(X, Y, bandwidth = bandwidth)

    boot_func <- function(i) {
      # Resample each x and y (residuals) independently
      resampled_x_ind <- sample(1:N, size = N, replace = TRUE)
      resampled_residuals_ind <- sample(1:N, size = N, replace = TRUE)

      # Order both
      resampled_x <- X[resampled_x_ind]
      resampled_residuals <- residuals[resampled_residuals_ind]
      resample_order <- order(resampled_x)
      resampled_x <- resampled_x[resample_order]
      resampled_residuals <- resampled_residuals[resample_order]

      t_stat <- get_hall_stat(resampled_x, resampled_residuals, m)$max_t_m
      return(t_stat)
    }

    # Do the bootstrap on multiple cores
    cl <- parallel::makeCluster(ncores)

    # Set seeds for reproducibility
    if (!is.null(seed)) {
      set.seed(seed)
      parallel::clusterSetRNGStream(cl, seed)
    }

    parallel::clusterExport(cl, c("X", "residuals", "N", "m", "get_hall_stat"),
                            envir = environment())

    t_vals <- unlist(parallel::parLapply(cl, 1:boot_num, boot_func))
    parallel::stopCluster(cl)

    plot <- plot_interval(X, Y, crit_interval, title = "Monotonicity Test: Critical Interval")
    p_val <- sum(t_vals >= t_stat) / length(t_vals)

    return(list(
      p = p_val,
      dist = t_vals,
      stat = t_stat,
      plot = plot,
      interval = crit_interval
    ))
  }

# Plot the "critical interval" which has the highest t-stat
plot_interval <- function(X, Y, interval, title) {
  plot_data <- data.frame(X = X, Y = Y)
  plot_data$InInterval <- ifelse(seq_along(X) >= interval[1] & seq_along(X) <= interval[2], "Yes", "No")

  interval_data <- plot_data[plot_data$InInterval == "Yes", ]

  lm_model <- lm(Y ~ X, data = interval_data)

  slope <- coef(lm_model)[2]
  intercept <- coef(lm_model)[1]

  plot <- ggplot(plot_data, aes(x = X, y = Y, color = .data$InInterval)) +
    geom_point(size = 2) +
    geom_abline(slope = slope, intercept = intercept, color = "red", linetype = "dashed", linewidth=1.5) +
    scale_color_manual(values = c("No" = "black", "Yes" = "red")) +
    labs(title = paste(title, " (Slope: ", round(slope, 3), ")", sep=""), x = "X", y = "Y", color = "Critical Interval")

  return(plot)
}

# Validate X, Y, and optionally, m inputs.
# Throws an error if invalid input otherwise returns nothing.
# Will only validate m if monotonicty=TRUE, otherwise ignore m validation
validate_inputs <- function(X, Y, m=NULL, monotonicity=TRUE) {
  # Check if values are NA, NaN or infinite.
  if (any(!is.finite(X)) || any(!is.finite(Y))) {
    stop("X and Y must contain only finite values (no NA, NaN, or Inf).")
  }

  if (!is.numeric(X) || !is.numeric(Y)) {
    stop("X and Y must be numeric vectors.")
  }

  if (length(X) != length(Y)) {
    stop("X and Y must be the same length.")
  }

  # Early return
  if(!monotonicity) return()

  # Check if m is valid
  if (!is.numeric(m) ||
      is.na(m) || m != as.integer(m) || m <= 0 || m >= length(X)) {
    stop("m must be a positive integer less than the length of the dataset.")
  }
}

#' Generate Kernel Plot
#'
#' Creates a scatter plot of the input vectors \eqn{X} and \eqn{Y}, and overlays
#' a Nadaraya-Watson kernel regression curve using the specified bandwidth.
#'
#' @param X Vector of x values.
#' @param Y Vector of y values.
#' @param bandwidth Kernel bandwidth used for the Nadaraya-Watson estimator. Can
#'                  be a single numeric value or a vector of bandwidths.
#'                  Default is calculated as
#'                  \code{bw.nrd(X) * (length(X) ^ -0.1)}.
#' @param nrows Number of rows in the facet grid if multiple bandwidths are provided.
#'              Does not do anything if only a single bandwidth value is provided.
#'              Default is \code{4}.
#' @return A ggplot object containing the scatter plot(s) with the kernel
#'         regression curve(s). If a vector of bandwidths is supplied, the plots
#'         are put into a grid using faceting.
#' @references
#'   Nadaraya, E. A. (1964). On estimating regression. \emph{Theory of
#'   Probability and Its Applications}, \strong{9}(1), 141–142.
#'
#'   Watson, G. S. (1964). Smooth estimates of regression functions.
#'   \emph{Sankhyā: The Indian Journal of Statistics, Series A}, 359-372.
#'
#' @examples
#' # Example 1: Basic plot on quadratic function
#' seed <- 42
#' set.seed(seed)
#' X <- runif(500)
#' Y <- X ^ 2 + rnorm(500, sd = 0.1)
#' plot <- create_kernel_plot(X, Y, bandwidth = bw.nrd(X) * (length(X) ^ -0.1))
#'
#' @export
create_kernel_plot <-
  function(X, Y, bandwidth = bw.nrd(X) * (length(X) ^ -0.1), nrows=4) {
    validate_inputs(X, Y, monotonicity=FALSE)

    if(!is.numeric(bandwidth)) {
      stop("'bandwidth' must be numeric value or vector")
    }

    # Check if nrows is an integer greater than zero
    if(nrows%%1!=0 | nrows<=0) {
      stop("'nrows' must be an integer greater than zero")
    }

    bandwidth_list <- as.list(bandwidth)

    # Set up the range for the x-axis
    x_range <- seq(min(X), max(X), length.out = 500)

    # Map the bandwidths into dataframes with their respective points
    # then concat with rbind
    plot_data <- do.call(rbind, Map(function(bw) {
      y_vals <- watson_est(X, Y, x_range, bandwidth = bw)
      data.frame(x = x_range,
                 y = y_vals,
                 bandwidth = sprintf("bandwidth = %.3f", bw))
    }, bandwidth_list))

    points_df <- data.frame(X = X, Y = Y)

    # Atleast 1 column
    n_cols <- max(ceiling(length(bandwidth_list) / nrows), 1)

    # Create facet wrap / grid plot
    plot <- ggplot() +
      geom_point(data = points_df, aes(x = .data$X, y = .data$Y)) +
      geom_line(
        data = plot_data,
        aes(x = .data$x, y = .data$y),
        color = "green",
        linewidth = 1.5
      ) +
      facet_wrap( ~ bandwidth, ncol = n_cols) +
      labs(title = "Nadaraya Watson Kernel Regression",
           x = "X",
           y = "Y")

    return(plot)
  }

# Get estimates with Nadaraya Watson kernel regression
watson_est <- function(X, Y, preds, bandwidth) {
  kernel_weights_x <- sapply(preds, function(x_i)
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
