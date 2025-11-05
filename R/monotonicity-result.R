# File contains methods associated with the "monotonicity_result" class
new_monotonicity_result <- function(plot = ggplot(),
                                    t_stats,
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

  m_index <- length(m_grid)
  t_stat <- t_stats[m_index]
  t_vals <- sapply(t_distributions, function(vec) vec[m_index])
  p <- sum(t_vals >= t_stat) / length(t_vals)

  structure(
    list(
      plot = plot,
      t_stats = t_stats,
      t_distributions = t_distributions,
      m_grid = m_grid,
      intervals = intervals,
      bandwidth = bandwidth,
      seed = seed,
      check_m = check_m,
      p = p
    ),
    class = "monotonicity_result"
  )
}


calc_p <- function(obj, m_index) {
  t_stat <- obj$t_stats[m_index]
  t_vals <- sapply(obj$t_distributions, function(vec) vec[m_index])
  sum(t_vals >= t_stat) / length(t_vals)
}

get_t <- function(obj, m_index) {
  t_vals <- sapply(obj$t_distributions, function(vec) vec[m_index])
  t_vals
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
  cat(sprintf("Bandwidth Value: %.3f\n", object$bandwidth))
  cat(sprintf("check_m: %s\n\n", if (isTRUE(x$check_m)) "TRUE" else "FALSE"))

  dist <- sapply(x$t_distributions, function(boot_s) boot_s[[index]])

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
}

#' @export
plot.monotonicity_result <- function(x, ...) {
  op <- par(ask = TRUE)
  on.exit(par(op))

  print(x$plot)

  if (isTRUE(x$check_m)) {
    m_vals <- head(x$m_grid, -1)
    t_vals <- head(x$t_stats, -1)
    p_vals <- vapply(seq_along(m_vals), function(i) calc_p(x, i), numeric(1))

    df <- data.frame(m = m_vals, t = t_vals, p = p_vals)

    p_m_p <- ggplot(df, aes(x = m, y = p)) +
      geom_line(color = "black") +
      geom_point(size = 1.5) +
      labs(x = "m", y = "p-value", title = "Monotonicity Test: Plot of m vs. p")

    print(p_m_p)

    p_m_t <- ggplot(df, aes(x = m, y = t)) +
      geom_line(color = "black") +
      geom_point(size = 1.5) +
      labs(x = "m", y = "t-statistic", title = "Monotonicity Test: Plot of m vs. t")

    print(p_m_t)
  }

}
