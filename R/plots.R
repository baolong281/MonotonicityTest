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
  function(X,
           Y,
           bandwidth = bw.nrd(X) * (length(X) ^ -0.1),
           nrows = 4) {
    validate_inputs(X, Y, monotonicity = FALSE)

    if (!is.numeric(bandwidth)) {
      stop("'bandwidth' must be numeric value or vector")
    }

    # Check if nrows is an integer greater than zero
    if (nrows %% 1 != 0 | nrows <= 0) {
      stop("'nrows' must be an integer greater than zero")
    }

    bandwidth_list <- as.list(bandwidth)
    sorted_bandwidths <- sort(bandwidth)

    # Set up the range for the x-axis
    x_range <- seq(min(X), max(X), length.out = 500)

    # Map the bandwidths into dataframes with their respective points
    # then concat with rbind
    plot_data <- do.call(rbind, Map(function(bw) {
      y_vals <- watson_est(X, Y, x_range, bandwidth = bw)
      data.frame(
        x = x_range,
        y = y_vals,
        bandwidth = factor(
          bw,
          levels = sorted_bandwidths,
          labels = sprintf("bandwidth = %.3f", sorted_bandwidths)
        )
      )
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
      facet_wrap(~ bandwidth, ncol = n_cols) +
      labs(title = "Nadaraya Watson Kernel Regression",
           x = "X",
           y = "Y")

    return(plot)
  }
