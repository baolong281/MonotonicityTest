test_that("Check if main function works", {
  N <- 200
  X <- runif(N)
  Y <- rnorm(N)

  boot_num <- 1

  res <-
    MonotonicityTest::monotonicity_test(X,
                                        Y,
                                        boot_num = boot_num,
                                        ncores = 1)

  # Check all types match
  expect_equal(boot_num, length(res$dist))
  expect_type(res$p, "double")
  expect_type(res$stat, "double")
  expect_true(is.vector(res$dist))
})

test_that("create_kernel_plot generates a plot without errors", {
  X <- seq(0, 1, length.out = 50)
  Y <- sin(2 * pi * X) + rnorm(50, sd = 0.1)

  plot_obj <- create_kernel_plot(X, Y)

  # Check that the returned object is of class "recordedplot"
  expect_true(inherits(plot_obj, "recordedplot"))
})
