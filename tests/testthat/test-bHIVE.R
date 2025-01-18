# test script for bHIVE.R - testcases are NOT comprehensive!

set.seed(42)

test_that("bHIVE handles clustering correctly", {
  data(iris)
  X <- as.matrix(iris[, 1:4])
  
  # Run bHIVE for clustering
  res <- bHIVE(X = X, 
               task = "clustering", 
               nAntibodies = 20, 
               beta = 5, 
               epsilon = 0.01, 
               maxIter = 10, 
               k = 3,
               verbose = FALSE)
  
  # Check output structure
  expect_named(res, c("antibodies", "assignments", "task"))
  expect_equal(res$task, "clustering")
  expect_true(is.matrix(res$antibodies))
  expect_true(is.integer(res$assignments))
  expect_length(res$assignments, nrow(X))
})

test_that("bHIVE handles classification correctly", {
  data(iris)
  X <- as.matrix(iris[, 1:4])
  y <- iris$Species
  
  # Run bHIVE for classification
  res <- bHIVE(X = X, 
               y = y, 
               task = "classification", 
               nAntibodies = 20, 
               beta = 5, 
               epsilon = 0.01, 
               maxIter = 10, 
               k = 3,
               verbose = FALSE)
  
  # Check output structure
  expect_named(res, c("antibodies", "assignments", "task"))
  expect_equal(res$task, "classification")
  expect_true(is.matrix(res$antibodies))
  expect_length(res$assignments, nrow(X))
})

test_that("bHIVE handles regression correctly", {
  data(iris)
  X <- as.matrix(iris[, 2:4])
  y <- iris$Sepal.Length
  
  # Run bHIVE for regression
  res <- bHIVE(X = X, 
               y = y, 
               task = "regression", 
               nAntibodies = 20, 
               beta = 5, 
               epsilon = 0.01, 
               maxIter = 10, 
               k = 3,
               verbose = FALSE)
  
  # Check output structure
  expect_named(res, c("antibodies", "assignments", "task"))
  expect_equal(res$task, "regression")
  expect_true(is.matrix(res$antibodies))
  expect_true(is.numeric(res$assignments))
  expect_length(res$assignments, nrow(X))
})

test_that("bHIVE throws an error for invalid task", {
  data(iris)
  X <- as.matrix(iris[, 1:4])
  
  # Invalid task
  expect_error(bHIVE(X = X, task = "invalid_task", nAntibodies = 20, beta = 5, epsilon = 0.01))
})
