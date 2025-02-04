# test script for refineB.R - testcases are NOT comprehensive!

# Define test cases
test_that("refineB: Input validation for X and A", {
  A <- matrix(runif(10), nrow = 2, ncol = 5)  # 2 antibodies, 5 features each
  X <- matrix(runif(20), nrow = 4, ncol = 5)  # 4 samples, 5 features each
  assignments <- c(1, 2, 1, 2)  # Valid assignments
  
  # Test valid inputs
  expect_silent(refineB(A = A, 
                        X = X, 
                        assignments = assignments, 
                        task = "clustering", 
                        verbose = FALSE))
  
  # Test mismatched dimensions between X and A
  X_invalid <- matrix(runif(24), nrow = 4, ncol = 6)
  expect_error(refineB(A = A, 
                       X = X_invalid, 
                       assignments = assignments, 
                       task = "clustering", 
                       verbose = FALSE),
               "Number of columns in 'X' must match the number of columns in 'A'.")
})

test_that("refineB: Validation for assignments", {
  A <- matrix(runif(10), nrow = 2, ncol = 5)
  X <- matrix(runif(20), nrow = 4, ncol = 5)
  
  # Valid assignments
  assignments <- c(1, 2, 1, 2)
  expect_silent(refineB(A = A, 
                        X = X, 
                        assignments = assignments, 
                        task = "clustering", 
                        verbose = FALSE))
  
  # Assignments out of range
  assignments_out_of_range <- c(1, 3, 1, 2)
  expect_error(refineB(A = A, 
                       X = X, 
                       assignments = assignments_out_of_range, 
                       task = "clustering", 
                       verbose = FALSE),
               "'assignments' contains invalid values: 3. Valid range is 1..2.")
})

test_that("refineB: Classification task validation", {
  A <- matrix(runif(10), nrow = 2, ncol = 5)
  X <- matrix(runif(20), nrow = 4, ncol = 5)
  y <- factor(c("A", "B", "A", "B"))  # Labels for classification
  assignments <- c(1, 2, 1, 2)
  
  # Valid classification task
  expect_silent(refineB(A = A, 
                        X = X, 
                        y = y, 
                        assignments = assignments, 
                        task = "classification", 
                        verbose = FALSE))
})

test_that("refineB: Basic clustering functionality", {
  A <- matrix(runif(10), nrow = 2, ncol = 5)  # 2 antibodies, 5 features
  X <- matrix(runif(20), nrow = 4, ncol = 5)  # 4 samples, 5 features
  assignments <- c(1, 2, 1, 2)  # Samples assigned to antibodies
  
  # Check output dimensions
  result <- refineB(A = A, 
                    X = X, 
                    assignments = assignments, 
                    task = "clustering", 
                    verbose = FALSE)
  expect_equal(dim(result), dim(A))
})

test_that("refineB: Regression task validation", {
  A <- matrix(runif(10), nrow = 2, ncol = 5)
  X <- matrix(runif(20), nrow = 4, ncol = 5)
  y <- c(1.5, 2.3, 1.8, 2.1)  # Numeric y for regression
  assignments <- c(1, 2, 1, 2)
  
  # Valid regression task
  expect_silent(refineB(A = A, 
                        X = X, 
                        y = y, 
                        assignments = assignments, 
                        task = "regression", 
                        verbose = FALSE))
  
  # Invalid y for regression (non-numeric)
  y_invalid <- factor(c("A", "B", "A", "B"))
  expect_error(refineB(A = A, 
                       X = X, 
                       y = y_invalid, 
                       assignments = assignments, 
                       task = "regression",
                       verbose = FALSE),
               "y must be numeric for regression.")
})

test_that("refineB: Edge case - Empty antibody assignment", {
  A <- matrix(runif(10), nrow = 2, ncol = 5)
  X <- matrix(runif(20), nrow = 4, ncol = 5)
  assignments <- c(1, 1, 1, 1)  # All samples assigned to the same antibody
  
  result <- refineB(A = A, 
                    X = X, 
                    assignments = assignments, 
                    task = "clustering", 
                    verbose = FALSE)
  expect_equal(dim(result), dim(A))  # Ensure output dimensions are preserved
})

test_that("refineB: Optimizer variants and hyperparameter customization (classification)", {
  A <- matrix(runif(10), nrow = 2, ncol = 5)
  X <- matrix(runif(20), nrow = 4, ncol = 5)
  assignments <- c(1, 2, 1, 2)
  y <- factor(c("A", "B", "A", "B"))
  
  # Test with basic SGD (default)
  result_sgd <- refineB(A = A, 
                        X = X, 
                        y = y, 
                        assignments = assignments, 
                        task = "classification", 
                        optimizer = "sgd", 
                        verbose = FALSE)
  expect_equal(dim(result_sgd), dim(A))
  
  # Test momentum with custom momentum_coef
  result_momentum <- refineB(A = A, 
                             X = X, 
                             y = y, 
                             assignments = assignments, 
                             task = "classification", 
                             optimizer = "momentum", 
                             momentum_coef = 0.8, 
                             verbose = FALSE)
  expect_equal(dim(result_momentum), dim(A))
  
  # Test adagrad with custom epsilon
  result_adagrad <- refineB(A = A, 
                            X = X, 
                            y = y, 
                            assignments = assignments, 
                            task = "classification", 
                            optimizer = "adagrad", 
                            epsilon = 1e-6, 
                            verbose = FALSE)
  expect_equal(dim(result_adagrad), dim(A))
  
  # Test adam with custom beta1 and beta2
  result_adam <- refineB(A = A, 
                         X = X, 
                         y = y, 
                         assignments = assignments, 
                         task = "classification", 
                         optimizer = "adam", 
                         beta1 = 0.8, 
                         beta2 = 0.95, 
                         verbose = FALSE)
  expect_equal(dim(result_adam), dim(A))
  
  # Test RMSProp with custom rmsprop_decay
  result_rmsprop <- refineB(A = A, 
                            X = X, 
                            y = y, 
                            assignments = assignments, 
                            task = "classification", 
                            optimizer = "rmsprop", 
                            rmsprop_decay = 0.85, 
                            verbose = FALSE)
  expect_equal(dim(result_rmsprop), dim(A))
})

test_that("refineB: Optimizer variants for regression", {
  A <- matrix(runif(10), nrow = 2, ncol = 5)
  X <- matrix(runif(20), nrow = 4, ncol = 5)
  y <- c(1.5, 2.3, 1.8, 2.1)
  assignments <- c(1, 2, 1, 2)
  
  optimizers <- c("sgd", "momentum", "adagrad", "adam", "rmsprop")
  for (opt in optimizers) {
    result <- refineB(A = A, 
                      X = X, 
                      y = y, 
                      assignments = assignments, 
                      task = "regression", 
                      optimizer = opt, 
                      verbose = FALSE)
    expect_equal(dim(result), dim(A))
  }
})
