% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{create_kernel_plot}
\alias{create_kernel_plot}
\title{Generate Kernel Plot}
\usage{
create_kernel_plot(X, Y, bandwidth = bw.nrd(X) * (length(X)^-0.1), nrows = 4)
}
\arguments{
\item{X}{Vector of x values.}

\item{Y}{Vector of y values.}

\item{bandwidth}{Kernel bandwidth used for the Nadaraya-Watson estimator. Can
be a single numeric value or a vector of bandwidths.
Default is calculated as
\code{bw.nrd(X) * (length(X) ^ -0.1)}.}

\item{nrows}{Number of rows in the facet grid if multiple bandwidths are provided.
Does not do anything if only a single bandwidth value is provided.
Default is \code{4}.}
}
\value{
A ggplot object containing the scatter plot(s) with the kernel
        regression curve(s). If a vector of bandwidths is supplied, the plots
        are put into a grid using faceting.
}
\description{
Creates a scatter plot of the input vectors \eqn{X} and \eqn{Y}, and overlays
a Nadaraya-Watson kernel regression curve using the specified bandwidth.
}
\examples{
# Example 1: Basic plot on quadratic function
seed <- 42
set.seed(seed)
X <- runif(500)
Y <- X ^ 2 + rnorm(500, sd = 0.1)
plot <- create_kernel_plot(X, Y, bandwidth = bw.nrd(X) * (length(X) ^ -0.1))

}
\references{
Nadaraya, E. A. (1964). On estimating regression. \emph{Theory of
  Probability and Its Applications}, \strong{9}(1), 141–142.

  Watson, G. S. (1964). Smooth estimates of regression functions.
  \emph{Sankhyā: The Indian Journal of Statistics, Series A}, 359-372.
}
