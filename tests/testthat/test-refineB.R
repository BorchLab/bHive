# test script for refineB.R - testcases are NOT comprehensive!

# Define test cases
test_that("refineB: Input validation for X and A", {
  A <- matrix(runif(10), nrow = 2, ncol = 5)  # 2 antibodies, 5 features each
  X <- matrix(runif(20), nrow = 4, ncol = 5)  # 4 samples, 5 features each
  assignments <- c(1, 2, 1, 2)  # Valid assignments
  
  # Test valid inputs
  expect_silent(refineB(A = A, X = X, assignments = assignments, task = "clustering"))
  
  # Test mismatched dimensions between X and A
  X_invalid <- matrix(runif(24), nrow = 4, ncol = 6)
  expect_error(refineB(A = A, X = X_invalid, assignments = assignments, task = "clustering"),
               "Number of columns in 'X' must match the number of columns in 'A'.")
})

test_that("refineB: Validation for assignments", {
  A <- matrix(runif(10), nrow = 2, ncol = 5)
  X <- matrix(runif(20), nrow = 4, ncol = 5)
  
  # Valid assignments
  assignments <- c(1, 2, 1, 2)
  expect_silent(refineB(A = A, X = X, assignments = assignments, task = "clustering"))
  
  # Assignments out of range
  assignments_out_of_range <- c(1, 3, 1, 2)
  expect_error(refineB(A = A, X = X, assignments = assignments_out_of_range, task = "clustering"),
               "'assignments' values must be between 1..nAb.")
})

test_that("refineB: Classification task validation", {
  A <- matrix(runif(10), nrow = 2, ncol = 5)
  X <- matrix(runif(20), nrow = 4, ncol = 5)
  y <- factor(c("A", "B", "A", "B"))  # Labels for classification
  assignments <- c(1, 2, 1, 2)
  
  # Valid classification task
  expect_silent(refineB(A = A, X = X, y = y, assignments = assignments, task = "classification"))
})

test_that("refineB: Basic clustering functionality", {
  A <- matrix(runif(10), nrow = 2, ncol = 5)  # 2 antibodies, 5 features
  X <- matrix(runif(20), nrow = 4, ncol = 5)  # 4 samples, 5 features
  assignments <- c(1, 2, 1, 2)  # Samples assigned to antibodies
  
  # Check output dimensions
  result <- refineB(A = A, X = X, assignments = assignments, task = "clustering")
  expect_equal(dim(result), dim(A))
})

test_that("refineB: Regression task validation", {
  A <- matrix(runif(10), nrow = 2, ncol = 5)
  X <- matrix(runif(20), nrow = 4, ncol = 5)
  y <- c(1.5, 2.3, 1.8, 2.1)  # Numeric y for regression
  assignments <- c(1, 2, 1, 2)
  
  # Valid regression task
  expect_silent(refineB(A = A, X = X, y = y, assignments = assignments, task = "regression"))
  
  # Invalid y for regression (non-numeric)
  y_invalid <- factor(c("A", "B", "A", "B"))
  expect_error(refineB(A = A, X = X, y = y_invalid, assignments = assignments, task = "regression"),
               "y must be numeric for regression.")
})

test_that("refineB: Edge case - Empty antibody assignment", {
  A <- matrix(runif(10), nrow = 2, ncol = 5)
  X <- matrix(runif(20), nrow = 4, ncol = 5)
  assignments <- c(1, 1, 1, 1)  # All samples assigned to the same antibody
  
  result <- refineB(A = A, X = X, assignments = assignments, task = "clustering")
  expect_equal(dim(result), dim(A))  # Ensure output dimensions are preserved
})