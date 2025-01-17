test_that("Check monotonicity_test works for valid inputs", {
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
  expect_true(is.vector(res$interval))
  expect_s3_class(res$plot, "ggplot")
})

test_that("Check if monotonicity_test input validation works", {
  # Valid inputs
  N <- 200
  X_valid <- runif(N)
  Y_valid <- rnorm(N)
  m_valid <- 50

  # Check if having non-finite values throws exception
  X_na <- X_valid
  X_na[1] <- NA
  expect_error(
    MonotonicityTest::monotonicity_test(X_na, Y_valid, m = m_valid),
    "X and Y must contain only finite values"
  )

  X_nan <- X_valid
  X_nan[1] <- NaN
  expect_error(
    MonotonicityTest::monotonicity_test(X_nan, Y_valid, m = m_valid),
    "X and Y must contain only finite values"
  )

  X_inf <- X_valid
  X_inf[1] <- Inf
  expect_error(
    MonotonicityTest::monotonicity_test(X_inf, Y_valid, m = m_valid),
    "X and Y must contain only finite values"
  )

  # Check if inequal lengths throws error
  Y_wrong_length <- rnorm(N + 1)
  expect_error(
    MonotonicityTest::monotonicity_test(X_valid, Y_wrong_length, m = m_valid),
    "X and Y must be the same length"
  )


  # Check if invalid M throws errors
  m_too_large <- N + 1
  expect_error(
    MonotonicityTest::monotonicity_test(X_valid, Y_valid, m = m_too_large),
    "m must be a positive integer less than the length of the dataset"
  )


  m_not_integer <- 50.5
  expect_error(
    MonotonicityTest::monotonicity_test(X_valid, Y_valid, m = m_not_integer),
    "m must be a positive integer less than the length of the dataset"
  )


  m_negative <- -10
  expect_error(
    MonotonicityTest::monotonicity_test(X_valid, Y_valid, m = m_negative),
    "m must be a positive integer less than the length of the dataset"
  )

})

test_that("Check if create_kernel_plot generates a plot without errors", {
  X <- 1:10
  Y <- X ^ 2

  # Normal plot with 1 bandwidth
  result <- create_kernel_plot(X, Y, bandwidth = 1.5)
  expect_s3_class(result, "ggplot")

  # Multiple bandwidths bandwidths
  result <- create_kernel_plot(X, Y, bandwidth = c(1, 2, 3, 4))
  expect_s3_class(result, "ggplot")
})

test_that("Check if create_kernel_plot input validation works", {
  # Valid inputs
  N <- 200
  X_valid <- runif(N)
  Y_valid <- rnorm(N)

  # Check if negative or non-integer nrows throws errors
  expect_error(
    MonotonicityTest::create_kernel_plot(X_valid, Y_valid, nrows=0),
    "'nrows' must be an integer greater than zero"
  )

  expect_error(
    MonotonicityTest::create_kernel_plot(X_valid, Y_valid, nrows=0.4),
    "'nrows' must be an integer greater than zero"
  )

  # Check if having non-finite values throws exception
  X_na <- X_valid
  X_na[1] <- NA
  expect_error(
    MonotonicityTest::create_kernel_plot(X_na, Y_valid),
    "X and Y must contain only finite values"
  )

  X_nan <- X_valid
  X_nan[1] <- NaN
  expect_error(
    MonotonicityTest::create_kernel_plot(X_nan, Y_valid),
    "X and Y must contain only finite values"
  )

  X_inf <- X_valid
  X_inf[1] <- Inf
  expect_error(
    MonotonicityTest::create_kernel_plot(X_inf, Y_valid),
    "X and Y must contain only finite values"
  )

  # Check if inequal lengths throws error
  Y_wrong_length <- rnorm(N + 1)
  expect_error(
    MonotonicityTest::create_kernel_plot(X_valid, Y_wrong_length),
    "X and Y must be the same length"
  )
})
