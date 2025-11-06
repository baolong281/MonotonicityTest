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
#' @param check_m Boolean value indicating whether to run the test for many different
#'                         values of \code{m}. This produces extra plots when calling
#'                         \code{plot} and has a marginal impact on performance.
#'                         Default is \code{FALSE}.
#' @return A \code{monotonicity_result} object. Has associated `print`,
#' `summary`, and `plot` S3 functions.
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
#' print(result)
#'
#' # Example 2: Usage on non-monotonic function
#' seed <- 42
#' set.seed(seed)
#'
#' X <- runif(500)
#' Y <- (X - 0.5) ^ 2 + rnorm(500, sd = 0.5)
#' result <- monotonicity_test(X, Y, boot_num = 25, seed = seed)
#'
#' print(result)
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
           check_m = FALSE,
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

    m_grid <- c(m)

    if(check_m) {
      m_grid <- c(
        round(seq(from = 10, to = 0.80 * N, by =  10)),
        m_grid
      )
    }

    # Get stat for actual dataset
    original_result <- get_hall_stat(X, Y, m_grid)
    t_stats <- original_result$t_values
    crit_intervals <- original_result$intervals

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

      get_hall_stat(resampled_x, resampled_residuals, m_grid)$t_values
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

    t_vals <- parallel::parLapply(cl, 1:boot_num, boot_func)
    parallel::stopCluster(cl)

    crit_interval = crit_intervals[[length(m_grid)]]

    # Un-invert the data for our plot
    if (negative) {
      Y <- -Y
    }

    plot <- plot_interval(X, Y, crit_interval, title = "Monotonicity Test: Critical Interval")

    return(new_monotonicity_result(
      plot = plot,
      t_stats = t_stats,
      t_distributions = t_vals,
      m_grid = m_grid,
      intervals = crit_intervals,
      bandwidth=bandwidth,
      seed = seed,
      check_m = check_m
    ))
  }


