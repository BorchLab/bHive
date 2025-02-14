# test script for swarmbHIVE.R - testcases are NOT comprehensive!
set.seed(42)

test_that("swarmbHIVE works for clustering (silhouette)", {
  data(iris)
  X <- as.matrix(iris[, 1:4])
  grid <- expand.grid(
    nAntibodies = c(10, 20),
    beta        = c(3, 5),
    epsilon     = c(0.01, 0.05)
  )
  
  # Tune hyperparameters
  res <- swarmbHIVE(
    X     = X,
    task  = "clustering",
    grid  = grid,
    metric = "silhouette",
    maxIter = 10,
    verbose = FALSE
  )
  
  # Check results structure
  expect_named(res, c("best_params", "results"))
  expect_true(is.data.frame(res$results))
  expect_equal(ncol(res$results), ncol(grid) + 1) # grid cols + metric_value
  expect_true(is.data.frame(res$best_params))
  expect_equal(nrow(res$best_params), 1)
})

test_that("swarmbHIVE works for classification (accuracy)", {
  data(iris)
  X <- as.matrix(iris[, 1:4])
  y <- iris$Species
  
  grid <- expand.grid(
    nAntibodies = c(10, 20),
    beta        = c(3, 5),
    epsilon     = c(0.01, 0.05)
  )
  
  # Tune hyperparameters
  res <- swarmbHIVE(
    X     = X, 
    y     = y, 
    task  = "classification",
    grid  = grid, 
    metric = "accuracy",
    maxIter = 10,
    verbose = FALSE
  )
  
  # Check results structure
  expect_named(res, c("best_params", "results"))
  expect_true(is.data.frame(res$results))
  expect_equal(ncol(res$results), ncol(grid) + 1)
  expect_true(is.data.frame(res$best_params))
  expect_equal(nrow(res$best_params), 1)
})

test_that("swarmbHIVE works for regression (rmse)", {
  data(iris)
  X <- as.matrix(iris[, 2:4])
  y <- iris$Sepal.Length
  
  grid <- expand.grid(
    nAntibodies = c(10, 20),
    beta        = c(3, 5),
    epsilon     = c(0.01, 0.05)
  )
  
  # Tune hyperparameters
  res <- swarmbHIVE(
    X     = X, 
    y     = y, 
    task  = "regression",
    grid  = grid, 
    metric = "rmse",
    maxIter = 10,
    verbose = FALSE
  )
  
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
    beta        = c(3, 5),
    epsilon     = c(0.01, 0.05)
  )
  
  # Invalid metric
  expect_error(
    swarmbHIVE(
      X     = X,
      task  = "clustering",
      grid  = grid,
      metric = "invalid_metric"
    )
  )
})

test_that("swarmbHIVE works for classification (balanced_accuracy)", {
  data(iris)
  X <- as.matrix(iris[, 1:4])
  y <- iris$Species
  
  grid <- expand.grid(
    nAntibodies = 10,
    beta        = 3,
    epsilon     = 0.01
  )
  
  res <- swarmbHIVE(
    X         = X,
    y         = y,
    task      = "classification",
    grid      = grid,
    metric    = "balanced_accuracy",
    maxIter   = 5, 
    verbose   = FALSE
  )
  
  expect_named(res, c("best_params", "results"))
  expect_equal(nrow(res$best_params), 1)
  expect_true(all(res$results$metric_value <= 1))
})

test_that("swarmbHIVE works for classification (f1)", {
  data(iris)
  X <- as.matrix(iris[, 1:4])
  y <- iris$Species
  
  grid <- expand.grid(
    nAntibodies = 10,
    beta        = 3,
    epsilon     = 0.01
  )
  
  res <- swarmbHIVE(
    X      = X,
    y      = y,
    task   = "classification",
    grid   = grid,
    metric = "f1",
    maxIter= 5,
    verbose= FALSE
  )
  
  expect_named(res, c("best_params", "results"))
  expect_equal(nrow(res$best_params), 1)
  # F1 generally in [0,1]
  expect_true(all(res$results$metric_value >= 0 & res$results$metric_value <= 1))
})

test_that("swarmbHIVE works for classification (kappa)", {
  data(iris)
  X <- as.matrix(iris[, 1:4])
  y <- iris$Species
  
  grid <- expand.grid(
    nAntibodies = 10,
    beta        = 3,
    epsilon     = 0.01
  )
  
  res <- swarmbHIVE(
    X      = X,
    y      = y,
    task   = "classification",
    grid   = grid,
    metric = "kappa",
    maxIter= 5,
    verbose= FALSE
  )
  
  expect_named(res, c("best_params", "results"))
  expect_equal(nrow(res$best_params), 1)
  # Kappa can be negative in worst cases, but let's just check it doesn't break
  expect_true(all(!is.na(res$results$metric_value)))
})

test_that("swarmbHIVE works for regression (mae)", {
  data(iris)
  X <- as.matrix(iris[, 2:4])
  y <- iris$Sepal.Length
  
  grid <- expand.grid(
    nAntibodies = 10,
    beta        = 3,
    epsilon     = 0.01
  )
  
  res <- swarmbHIVE(
    X      = X,
    y      = y,
    task   = "regression",
    grid   = grid,
    metric = "mae",
    maxIter= 5,
    verbose= FALSE
  )
  
  expect_named(res, c("best_params", "results"))
  expect_equal(nrow(res$best_params), 1)
  # MAE >= 0
  expect_true(all(res$results$metric_value >= 0))
})

test_that("swarmbHIVE works for regression (r2)", {
  data(iris)
  X <- as.matrix(iris[, 2:4])
  y <- iris$Sepal.Length
  
  grid <- expand.grid(
    nAntibodies = 10,
    beta        = 3,
    epsilon     = 0.01
  )
  
  res <- swarmbHIVE(
    X      = X,
    y      = y,
    task   = "regression",
    grid   = grid,
    metric = "r2",
    maxIter= 5,
    verbose= FALSE
  )
  
  expect_named(res, c("best_params", "results"))
  expect_equal(nrow(res$best_params), 1)
  # R2 can be negative. Just ensure no NA
  expect_true(all(!is.na(res$results$metric_value)))
})

test_that("swarmbHIVE works for clustering (davies_bouldin)", {
  # skip if clusterCrit is not available
  skip_if_not_installed("clusterCrit")
  
  data(iris)
  X <- as.matrix(iris[, 1:4])
  grid <- expand.grid(
    nAntibodies = 5,
    beta        = 3,
    epsilon     = 0.01
  )
  
  res <- swarmbHIVE(
    X       = X,
    task    = "clustering",
    grid    = grid,
    metric  = "davies_bouldin",
    maxIter = 5,
    verbose = FALSE
  )
  
  # DB index is smaller=better
  # We won't check absolute range, just that we got a real number
  expect_named(res, c("best_params", "results"))
  expect_equal(nrow(res$best_params), 1)
  expect_true(all(!is.na(res$results$metric_value)))
})

test_that("swarmbHIVE works for clustering (calinski_harabasz)", {
  skip_if_not_installed("clusterCrit")
  
  data(iris)
  X <- as.matrix(iris[, 1:4])
  grid <- expand.grid(
    nAntibodies = 5,
    beta        = 3,
    epsilon     = 0.01
  )
  
  res <- swarmbHIVE(
    X       = X,
    task    = "clustering",
    grid    = grid,
    metric  = "calinski_harabasz",
    maxIter = 5,
    verbose = FALSE
  )
  
  expect_named(res, c("best_params", "results"))
  expect_equal(nrow(res$best_params), 1)
  expect_true(all(!is.na(res$results$metric_value)))
})
