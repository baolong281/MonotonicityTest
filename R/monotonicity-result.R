# File contains methods associated with the "monotonicity_result" class
new_monotonicity_result <- function(plot = ggplot(),
                                    t_stats = numeric(),
                                    t_distributions = list(),
                                    m_grid = numeric(),
                                    intervals = list(),
                                    bandwidth = double(),
                                    seed = NULL,
                                    check_m = FALSE
                                    ) {
  stopifnot(is.numeric(t_stats))
  stopifnot(is.numeric(m_grid))
  stopifnot(is.double(bandwidth))
  stopifnot(is.list(t_distributions))
  stopifnot(is.list(intervals))

  structure(
    list(
      plot = plot,
      t_stats = t_stats,
      t_distributions = t_distributions,
      m_grid = m_grid,
      intervals = intervals,
      bandwidth = bandwidth,
      seed = seed,
      check_m = check_m
    ),
    class = "monotonicity_result"
  )
}


calc_p <- function(obj, m_index) {
  t_stat <- obj$t_stats[m_index]
  t_vals <- sapply(obj$t_distributions, function(vec) vec[m_index])
  sum(t_vals >= t_stat) / length(t_vals)
}

#' @export
print.monotonicity_result <- function(x, ...) {
  index <- length(x$m_grid)
  p_val <- calc_p(x, index)
  t_stat <- x$t_stats[index]

  cat("\n")
  cat(sprintf("P-Value: %.3f\n", p_val))
  cat(sprintf("T-Statistic: %.3f\n\n", t_stat))
  cat("Call 'summary()' for more information.\n")
}

#' @export
summary.monotonicity_result <- function(object, ...) {
  x <- object

  index <- length(x$m_grid)
  p_val <- calc_p(x, index)
  t_stat <- x$t_stats[index]
  interval <- x$intervals[[index]]

  cat("\n")
  cat(sprintf("P-Value: %-10.3f     T-Statistic: %-10.3f\n", p_val, t_stat))
  cat(sprintf("Critical Interval: [%d, %d]   Random Seed: %s\n",
              interval[1], interval[2],
              if (is.null(x$seed)) "<none>" else as.character(x$seed)))
  cat(sprintf("Bandwidth Value: %.3f\n\n", object$bandwidth))

  dist <- x$t_distributions[[index]]

  dist_stats <- c(
    min(dist),
    quantile(dist),
    median(dist),
    quantile(dist, 0.75),
    max(dist)
  )

  cat(sprintf("Bootstrap Distribution (n=%d):\n", length(x$dist)))
  cat("      Min        1Q    Median        3Q       Max\n")
  cat(
    sprintf(
      "%10.3e %10.3e %10.3e %10.3e %10.3e",
      dist_stats[1],
      dist_stats[2],
      dist_stats[3],
      dist_stats[4],
      dist_stats[5]
    ),
    "\n\n"
  )

  if (isTRUE(x$check_m)) {
    cat("Statistics and P-Values for all m values:\n\n")

    m_vals <- head(x$m_grid, -1)
    t_vals <- head(x$t_stats, -1)
    p_vals <- vapply(seq_along(m_vals), function(i) calc_p(x, i), numeric(1))

    cat(sprintf("%-12s", ""), paste(sprintf("%-10d", m_vals), collapse = " "), "\n")
    cat(sprintf("%-12s", "T-Statistic"), paste(sprintf("%-10.3f", t_vals), collapse = " "), "\n")
    cat(sprintf("%-12s", "P-Value"), paste(sprintf("%-10.3f", p_vals), collapse = " "), "\n\n")
  }

}

#' @export
plot.monotonicity_result <- function(x, ...) {
  print(x$plot)
}
