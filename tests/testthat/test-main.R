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

  expect_equal(boot_num, length(res$dist))
})
