# File contains methods associated with the "monotonicity_result" class

new_monotonicity_result <- function(p = double(), dist=numeric(), stat=double(),
                                    plot = ggplot( ), interval = numeric(),
                                    bandwidth = double(),
                                    seed=NULL) {
  stopifnot(is.double(p))
  stopifnot(is.numeric(dist))
  stopifnot(is.double(stat))
  stopifnot(is.double(bandwidth))
  stopifnot(is.numeric(interval))
  stopifnot(length(interval) == 2)

  structure(
    list(
      p = p,
      dist = dist,
      stat = stat,
      plot = plot,
      interval = interval,
      bandwidth=bandwidth,
      seed = seed
    ),
    class = "monotonicity_result"
  )
}

#' @export
print.monotonicity_result <- function(x, ...) {
  cat("\n")
  cat(sprintf("P-Value: %.3f\n", x$p))
  cat(sprintf("T-Statistic: %.3f\n\n", x$stat))
  cat("Call 'summary()' for more information.\n")
}

#' @export
summary.monotonicity_result <- function(object, ...) {
  x <- object
  cat("\n")
  cat(sprintf("P-Value: %-10.3f     T-Statistic: %-10.3f\n", x$p, x$stat))
  cat(sprintf("Critical Interval: [%d, %d]   Random Seed: %s\n",
              x$interval[1], x$interval[2],
              if (is.null(x$seed)) "<none>" else as.character(x$seed)))
  cat(sprintf("Bandwidth Value: %.3f\n\n", object$bandwidth))

  dist <- x$dist
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
  print(x$plot)
}
