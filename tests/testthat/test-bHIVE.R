# test script for bHIVE.R - testcases are NOT comprehensive!

set.seed(42)

test_that("bHIVE handles different affinity functions correctly", {
  data(iris)
  X <- as.matrix(iris[, 1:4])
  
  affinity_funcs <- c("gaussian", "laplace", "polynomial", "cosine")
  for (aff in affinity_funcs) {
    expect_silent(
      res <- bHIVE(X = X, 
                   task = "clustering", 
                   affinityFunc = aff, 
                   distFunc = "euclidean", 
                   nAntibodies = 10, 
                   maxIter = 5, 
                   verbose = FALSE)
    )
    expect_type(res, "list")
    expect_named(res, c("antibodies", "assignments", "task"))
  }
})

test_that("bHIVE handles different distance functions correctly", {
  data(iris)
  X <- as.matrix(iris[, 1:4])
  
  dist_funcs <- c("euclidean", "manhattan", "minkowski")
  for (dist in dist_funcs) {
    expect_silent(
      res <- bHIVE(X = X, 
                   task = "clustering", 
                   affinityFunc = "gaussian", 
                   distFunc = dist, 
                   nAntibodies = 10, 
                   maxIter = 5, 
                   verbose = FALSE)
    )
    expect_type(res, "list")
    expect_named(res, c("antibodies", "assignments", "task"))
  }
})

test_that("bHIVE works with different tasks", {
  data(iris)
  X <- as.matrix(iris[, 1:4])
  y_class <- iris$Species
  y_reg <- iris$Sepal.Length
  
  # Classification
  expect_silent(
    res_class <- bHIVE(X = X, 
                       y = y_class, 
                       task = "classification", 
                       affinityFunc = "gaussian", 
                       distFunc = "euclidean", 
                       nAntibodies = 10, 
                       maxIter = 5,
                       verbose = FALSE)
  )
  expect_type(res_class, "list")
  expect_named(res_class, c("antibodies", "assignments", "task"))
  expect_equal(length(res_class$assignments), nrow(X))
  
  # Regression
  expect_silent(
    res_reg <- bHIVE(X = X, 
                     y = y_reg, 
                     task = "regression", 
                     affinityFunc = "gaussian", 
                     distFunc = "euclidean", 
                     nAntibodies = 10, 
                     maxIter = 5, 
                     verbose = FALSE)
  )
  expect_type(res_reg, "list")
  expect_named(res_reg, c("antibodies", "assignments", "task"))
  expect_equal(length(res_reg$assignments), nrow(X))
  
  # Clustering
  expect_silent(
    res_cluster <- bHIVE(X = X, 
                         task = "clustering", 
                         affinityFunc = "gaussian", 
                         distFunc = "euclidean", 
                         nAntibodies = 10, 
                         maxIter = 5, 
                         verbose = FALSE)
  )
  expect_type(res_cluster, "list")
  expect_named(res_cluster, c("antibodies", "assignments", "task"))
  expect_equal(length(res_cluster$assignments), nrow(X))
})

test_that("bHIVE handles different initialization methods correctly", {
  data(iris)
  X <- as.matrix(iris[, 1:4])
  
  init_methods <- c("sample", "random", "random_uniform", "kmeans++")
  for (init in init_methods) {
    expect_silent(
      res <- bHIVE(X = X, 
                   task = "clustering", 
                   affinityFunc = "gaussian", 
                   distFunc = "euclidean", 
                   nAntibodies = 10, 
                   maxIter = 5, 
                   initMethod = init, 
                   verbose = FALSE)
    )
    expect_type(res, "list")
    expect_named(res, c("antibodies", "assignments", "task"))
  }
})

test_that("bHIVE returns correct structure and data types", {
  data(iris)
  X <- as.matrix(iris[, 1:4])
  res <- bHIVE(X = X, 
               task = "clustering", 
               nAntibodies = 10, 
               maxIter = 5)
  
  # Check structure
  expect_type(res, "list")
  expect_named(res, c("antibodies", "assignments", "task"))
  
  # Check types of components
  expect_type(res$antibodies, "double")
  expect_type(res$assignments, "numeric")
  expect_equal(res$task, "clustering")
})
