# test script for honeycombHIVE.R - testcases are NOT comprehensive!

set.seed(42)
test_that("honeycombHIVE runs successfully for clustering task", {
  data(iris)
  X <- as.matrix(iris[, 1:4])
  
  # Run honeycombHIVE for clustering
  res <- honeycombHIVE(X = X, 
                       task = "clustering", 
                       layers = 3, 
                       nAntibodies = 15, 
                       beta = 5, 
                       maxIter = 5,
                       verbose = FALSE)
  
  # Test structure and results
  expect_type(res, "list")
  expect_length(res, 3) # 3 layers
  for (layer in res) {
    expect_named(layer, c("antibodies", "assignments", "task", "membership"))
    expect_equal(layer$task, "clustering")
    expect_true(is.data.frame(layer$antibodies))
    expect_true(is.numeric(layer$assignments))
    expect_equal(length(layer$membership), 150)
  }
})

test_that("honeycombHIVE runs successfully for classification task", {
  data(iris)
  X <- as.matrix(iris[, 1:4])
  y <- iris$Species
  
  # Run honeycombHIVE for classification
  res <- honeycombHIVE(X = X, 
                       y = as.factor(y), 
                       task = "classification", 
                       layers = 2, 
                       nAntibodies = 10, 
                       beta = 3, 
                       maxIter = 5, 
                       verbose = FALSE)
  
  # Test structure and results
  expect_type(res, "list")
  expect_length(res, 2) # 2 layers
  for (layer in res) {
    expect_named(layer, c("antibodies", "assignments", "task", "predictions", "membership" ))
    expect_equal(layer$task, "classification")
    expect_true(is.data.frame(layer$antibodies))
  }
})

test_that("honeycombHIVE runs successfully for regression task", {
  if (!requireNamespace("MASS", quietly = TRUE)) install.packages("MASS")
  library(MASS)
  data(Boston)
  X <- as.matrix(Boston[, -14])
  y <- Boston$medv
  
  # Run honeycombHIVE for regression
  res <- honeycombHIVE(X = X, 
                       y = y, 
                       task = "regression", 
                       layers = 3, 
                       nAntibodies = 15, 
                       beta = 5, 
                       maxIter = 5,
                       verbose = FALSE)
  
  # Test structure and results
  expect_type(res, "list")
  expect_length(res, 3) # 3 layers
  for (layer in res) {
    expect_named(layer, c("antibodies", "assignments", "task", "predictions", "membership"))
    expect_equal(layer$task, "regression")
    expect_true(is.data.frame(layer$antibodies))
    expect_true(is.numeric(layer$assignments))
    expect_equal(length(layer$membership), 506)
  }
})

test_that("honeycombHIVE handles different affinity functions", {
  data(iris)
  X <- as.matrix(iris[, 1:4])
  
  affinity_funcs <- c("gaussian", "laplace", "polynomial", "cosine")
  for (aff in affinity_funcs) {
    expect_silent(
      res <- honeycombHIVE(X = X, 
                           task = "clustering", 
                           layers = 2, 
                           nAntibodies = 10, 
                           affinityFunc = aff, 
                           maxIter = 5,
                           verbose = FALSE)
    )
    expect_length(res, 2)
    for (layer in res) {
      expect_equal(layer$task, "clustering")
    }
  }
})

test_that("honeycombHIVE handles different distance functions", {
  data(iris)
  X <- as.matrix(iris[, 1:4])
  
  dist_funcs <- c("euclidean", "manhattan", "minkowski")
  for (dist in dist_funcs) {
    expect_silent(
      res <- honeycombHIVE(X = X, 
                           task = "clustering", 
                           layers = 2, 
                           nAntibodies = 10, 
                           distFunc = dist, 
                           maxIter = 5,
                           verbose = FALSE)
    )
    expect_length(res, 2)
    for (layer in res) {
      expect_equal(layer$task, "clustering")
    }
  }
})

