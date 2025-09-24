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
