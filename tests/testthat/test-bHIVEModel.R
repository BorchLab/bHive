# test script for bHIVEmodel.R - testcases are NOT comprehensive!

set.seed(42)
library(caret)

test_that("bHIVEmodel integrates with caret for classification", {
  data(iris)
  X <- as.matrix(iris[, 1:4])
  y <- iris$Species
  
  # Train model using caret
  trainControl <- trainControl(method = "cv", number = 3)
  model <- train(
    x = X,
    y = y,
    method = bHIVEmodel,
    trControl = trainControl,
    tuneGrid = expand.grid(
      nAntibodies = c(10, 20),
      beta = c(3, 5),
      epsilon = c(0.01, 0.05)
    ),
    verbose = FALSE
  )
  
  # Check output structure
  expect_s3_class(model, "train")
  expect_true(!is.null(model$results))
})

test_that("bHIVEmodel integrates with caret for regression", {
  data(iris)
  X <- as.matrix(iris[, 2:4])
  y <- iris$Sepal.Length
  
  # Train model using caret
  trainControl <- trainControl(method = "cv", number = 3)
  model <- train(
    x = X,
    y = y,
    method = bHIVEmodel,
    trControl = trainControl,
    tuneGrid = expand.grid(
      nAntibodies = c(10, 20),
      beta = c(3, 5),
      epsilon = c(0.01, 0.05)
    ),
    verbose = FALSE
  )
  
  # Check output structure
  expect_s3_class(model, "train")
  expect_true(!is.null(model$results))
})