# test script for refineB.R - testcases are NOT comprehensive!

test_that("refineB Classification: basic usage with categorical_crossentropy, push_away=TRUE", {
  set.seed(1)
  
  # Create a small toy classification dataset
  X <- matrix(c(1,2,3,4,
                2,3,4,5,
                10,11,12,13,
                11,12,13,14), nrow=4, byrow=TRUE)
  y <- factor(c("A","A","B","B"))
  
  # Suppose we have 2 antibodies
  A <- matrix(c(2,3,3,4,
                10,11,12,13), nrow=2, byrow=TRUE)
  
  # Assignments (each row of X => index 1 or 2)
  assignments <- c(1,1,2,2)  # first 2 rows to antibody1, last 2 to antibody2
  
  # refineB
  A_refined <- refineB(
    A = A,
    X = X,
    y = y,
    assignments = assignments,
    loss = "categorical_crossentropy",
    task = "classification",
    steps = 5,
    lr = 0.1,
    push_away = TRUE
  )
  
  # Check dimension is same
  expect_equal(dim(A_refined), dim(A))
  
  # Expect that the refined positions differ from original
  expect_false(all(A_refined == A))
})

test_that("refineB Classification: no push_away, simple mse approach on factor data", {
  set.seed(2)
  
  # Another small dataset
  X <- matrix(runif(20), nrow=5, ncol=4)
  y <- factor(c("X","X","Y","Y","X"))
  
  A <- matrix(runif(8), nrow=2, ncol=4)
  assignments <- c(1,1,2,2,1)  # naive cluster assignment
  
  # Use a classification naive 'mse' approach, push_away=FALSE
  A_refined <- refineB(
    A = A,
    X = X,
    y = y,
    assignments = assignments,
    loss = "mse",
    task = "classification",
    steps = 3,
    lr = 0.05,
    push_away = FALSE
  )
  
  expect_equal(dim(A_refined), c(2,4))
  # We expect some difference but it's possible if random was small
  expect_false(all(A_refined == A))
})

test_that("refineB Regression: MSE refinement changes prototype positions", {
  set.seed(3)
  
  X <- matrix(rnorm(20), nrow=5, ncol=4)
  y <- rnorm(5, mean=10, sd=2)  # numeric
  A <- matrix(rnorm(8), nrow=2, ncol=4)
  
  # Assign points to the 2 antibodies arbitrarily
  assignments <- c(1,1,2,2,1)
  
  A_refined <- refineB(
    A = A,
    X = X,
    y = y,
    assignments = assignments,
    loss = "mse",
    task = "regression",
    steps = 5,
    lr = 0.01
  )
  
  expect_equal(dim(A_refined), c(2,4))
  expect_false(all(A_refined == A))
})

test_that("refineB Regression: MAE with push_away should have no effect (ignored in regression)", {
  set.seed(4)
  
  X <- matrix(rnorm(20), nrow=5)
  y <- rnorm(5, mean=5)
  A <- matrix(rnorm(4), nrow=2)
  assignments <- c(1,1,2,2,1)
  
  A_refined_push <- refineB(
    A = A,
    X = X,
    y = y,
    assignments = assignments,
    loss = "mae",
    task = "regression",
    steps = 3,
    lr = 0.1,
    push_away = TRUE
  )
  A_refined_nopush <- refineB(
    A = A,
    X = X,
    y = y,
    assignments = assignments,
    loss = "mae",
    task = "regression",
    steps = 3,
    lr = 0.1,
    push_away = FALSE
  )
  expect_equal(dim(A_refined_push), c(2,1))
  expect_equal(dim(A_refined_nopush), c(2,1))
  expect_false(all(A_refined_push == A))
  expect_false(all(A_refined_nopush == A))
})

test_that("refineB empty cluster: assigned no points to an antibody doesn't break the function", {
  set.seed(5)
  
  X <- matrix(rnorm(30), nrow=6, ncol=5)
  y <- factor(c("A","A","B","B","B","A"))
  A <- matrix(rnorm(10), nrow=2)
  
  # Suppose antibody 2 gets no points
  assignments <- c(1,1,1,1,1,1)
  
  # classification
  expect_silent({
    A_ref <- refineB(
      A = A,
      X = X,
      y = y,
      assignments = assignments,
      loss = "categorical_crossentropy",
      task = "classification",
      steps = 3,
      lr = 0.1
    )
  })
  # Should not error. Dimension is same
  expect_equal(dim(A_ref), dim(A))
})

test_that("refineB invalid arguments: wrong assignment length or mismatch task/loss", {
  set.seed(6)
  
  X <- matrix(rnorm(20), nrow=5)
  y_fac <- factor(c("A","A","B","B","A"))
  A <- matrix(rnorm(6), nrow=2)
  
  expect_error(
    refineB(
      A = A,
      X = X,
      y = y_fac,
      assignments = c(1,1,2,2), # length mismatch
      loss = "categorical_crossentropy",
      task = "classification"
    ), 
    "must match the number of rows in X"
  )
  
  # Using a regression loss but task=classification
  expect_error(
    refineB(
      A = A,
      X = X,
      y = y_fac,
      assignments = c(1,1,2,2,2),
      loss = "mae",
      task = "classification"
    ),
    "For classification, loss must be one of"
  )
})

test_that("refineB huber loss: modifies prototypes in regression scenario", {
  set.seed(7)
  
  X <- matrix(rnorm(20), nrow=5)
  y <- rnorm(5, mean=10)
  A <- matrix(rnorm(6), nrow=2)
  assignments <- c(1,1,2,2,1)
  
  A_ref <- refineB(
    A = A,
    X = X,
    y = y,
    assignments = assignments,
    loss = "huber",
    task = "regression",
    steps = 5,
    lr = 0.05,
    huber_delta = 1.0
  )
  
  expect_equal(dim(A_ref), c(2,1))
  expect_false(all(A_ref == A))
})

test_that("refineB poisson loss: no errors for regression scenario", {
  set.seed(8)
  
  X <- matrix(rpois(25, lambda=5), nrow=5, ncol=5)  # random counts
  y <- rpois(5, lambda=10)  # regression-like
  A <- matrix(rnorm(10), nrow=2)
  assignments <- c(1,2,1,2,1)
  
  A_ref <- refineB(
    A = A,
    X = X,
    y = y,
    assignments = assignments,
    loss = "poisson",
    task = "regression",
    steps = 3,
    lr = 0.02
  )
  
  expect_equal(dim(A_ref), dim(A))
  expect_false(all(A_ref == A))
})

test_that("refineB cosine classification: dimension check and changes", {
  set.seed(9)
  
  X <- matrix(rnorm(16), nrow=4)
  y <- factor(c("C","C","D","D"))
  A <- matrix(rnorm(8), nrow=2)
  assignments <- c(1,1,2,2)
  
  A_ref <- refineB(
    A = A,
    X = X,
    y = y,
    assignments = assignments,
    loss = "cosine",
    task = "classification",
    steps = 4,
    lr = 0.1,
    push_away = TRUE
  )
  
  expect_equal(dim(A_ref), c(2,2))
  expect_false(all(A_ref == A))
})

