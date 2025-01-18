# test script for swarmbHIVE.R - testcases are NOT comprehensive!

set.seed(42)

test_that("swarmbHIVE works for clustering", {
  data(iris)
  X <- as.matrix(iris[, 1:4])
  grid <- expand.grid(
    nAntibodies = c(10, 20),
    beta = c(3, 5),
    epsilon = c(0.01, 0.05)
  )
  
  # Tune hyperparameters
  res <- swarmbHIVE(X = X, task = "clustering", grid = grid, metric = "silhouette", verbose = FALSE)
  
  # Check results structure
  expect_named(res, c("best_params", "results"))
  expect_true(is.data.frame(res$results))
  expect_equal(ncol(res$results), ncol(grid) + 1) # grid columns + metric_value
  expect_true(is.data.frame(res$best_params))
  expect_equal(nrow(res$best_params), 1)
})

test_that("swarmbHIVE works for classification", {
  data(iris)
  X <- as.matrix(iris[, 1:4])
  y <- iris$Species
  grid <- expand.grid(
    nAntibodies = c(10, 20),
    beta = c(3, 5),
    epsilon = c(0.01, 0.05)
  )
  
  # Tune hyperparameters
  res <- swarmbHIVE(X = X, 
                    y = y, 
                    task = "classification", 
                    grid = grid, 
                    metric = "accuracy", 
                    verbose = FALSE)
  
  # Check results structure
  expect_named(res, c("best_params", "results"))
  expect_true(is.data.frame(res$results))
  expect_equal(ncol(res$results), ncol(grid) + 1)
  expect_true(is.data.frame(res$best_params))
  expect_equal(nrow(res$best_params), 1)
})

test_that("swarmbHIVE works for regression", {
  data(iris)
  X <- as.matrix(iris[, 2:4])
  y <- iris$Sepal.Length
  grid <- expand.grid(
    nAntibodies = c(10, 20),
    beta = c(3, 5),
    epsilon = c(0.01, 0.05)
  )
  
  # Tune hyperparameters
  res <- swarmbHIVE(X = X, 
                    y = y, 
                    task = "regression", 
                    grid = grid, 
                    metric = "rmse", 
                    verbose = FALSE)
  
  # Check results structure
  expect_named(res, c("best_params", "results"))
  expect_true(is.data.frame(res$results))
  expect_equal(ncol(res$results), ncol(grid) + 1)
  expect_true(is.data.frame(res$best_params))
  expect_equal(nrow(res$best_params), 1)
})

test_that("swarmbHIVE throws an error for invalid metric", {
  data(iris)
  X <- as.matrix(iris[, 1:4])
  grid <- expand.grid(
    nAntibodies = c(10, 20),
    beta = c(3, 5),
    epsilon = c(0.01, 0.05)
  )
  
  # Invalid metric
  expect_error(swarmbHIVE(X = X, task = "clustering", grid = grid, metric = "invalid_metric"))
})
