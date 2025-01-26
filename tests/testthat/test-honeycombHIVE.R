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

test_that("honeycombHIVE refine=TRUE for classification with cross-entropy", {
  data(iris)
  X <- as.matrix(iris[,1:4])
  y <- iris$Species
  
  # Compare refine=FALSE vs. refine=TRUE
  res_no_ref <- honeycombHIVE(
    X = X,
    y = y,
    task = "classification",
    layers = 2,
    nAntibodies = 10,
    maxIter = 5,
    verbose = FALSE,
    refine = FALSE
  )
  
  res_ref <- honeycombHIVE(
    X = X,
    y = y,
    task = "classification",
    layers = 2,
    nAntibodies = 10,
    maxIter = 5,
    verbose = FALSE,
    refine = TRUE,
    refineLoss = "categorical_crossentropy",
    refineSteps = 3,
    refineLR = 0.05,
    refinePushAway = TRUE
  )
  
  # Check structure
  expect_length(res_ref, 2)
  for (layer_obj in res_ref) {
    expect_named(layer_obj, c("antibodies","assignments","task",
                              "predictions","membership"))
    expect_true(is.data.frame(layer_obj$antibodies))
  }
  
  # Check that the final prototypes differ
  final_no_ref <- res_no_ref[[2]]$antibodies
  final_ref    <- res_ref[[2]]$antibodies
  expect_false(isTRUE(all.equal(final_no_ref, final_ref)))
  
  preds_no_ref <- res_no_ref[[2]]$predictions
  preds_ref    <- res_ref[[2]]$predictions
})

test_that("honeycombHIVE refine=TRUE for regression (MSE)", {
  if (!requireNamespace("MASS", quietly = TRUE)) install.packages("MASS")
  library(MASS)
  data(Boston)
  X <- as.matrix(Boston[, -14])
  y <- Boston$medv
  
  # Without refine
  res_no_ref <- honeycombHIVE(
    X = X,
    y = y,
    task = "regression",
    layers = 2,
    nAntibodies = 10,
    maxIter = 5,
    verbose = FALSE,
    refine = FALSE
  )
  
  # With refine
  res_ref <- honeycombHIVE(
    X = X,
    y = y,
    task = "regression",
    layers = 2,
    nAntibodies = 10,
    maxIter = 5,
    verbose = FALSE,
    refine = TRUE,
    refineLoss = "mse",
    refineSteps = 5,
    refineLR = 0.01
  )
  
  # Check structure
  expect_length(res_ref, 2)
  for (ly in res_ref) {
    expect_named(ly, c("antibodies","assignments","task","predictions","membership"))
    expect_equal(ly$task, "regression")
    expect_true(is.data.frame(ly$antibodies))
  }
  
  # Compare final prototypes
  final_no_ref <- res_no_ref[[2]]$antibodies
  final_ref    <- res_ref[[2]]$antibodies
  expect_false(isTRUE(all.equal(final_no_ref, final_ref)))
})

test_that("honeycombHIVE refineSteps=0 does not change prototypes", {
  data(iris)
  X <- as.matrix(iris[,1:4])
  
  # refine=TRUE but steps=0 => effectively no gradient updates
  res_steps0 <- honeycombHIVE(
    X = X,
    task = "clustering",
    layers = 1,
    nAntibodies = 8,
    maxIter = 3,
    verbose = FALSE,
    refine = TRUE,
    refineLoss = "mae",
    refineSteps = 0,  
    refineLR = 0.1
  )
  
  # normal refineSteps>0
  res_steps5 <- honeycombHIVE(
    X = X,
    task = "clustering",
    layers = 1,
    nAntibodies = 8,
    maxIter = 3,
    verbose = FALSE,
    refine = TRUE,
    refineLoss = "mae",
    refineSteps = 5,
    refineLR = 0.1
  )
  
  # Compare final prototypes
  no_refine <- honeycombHIVE(
    X = X,
    task = "clustering",
    layers = 1,
    nAntibodies = 8,
    maxIter = 3,
    verbose = FALSE,
    refine = FALSE
  )
  
  prot_steps0 <- res_steps0[[1]]$antibodies
  prot_steps5 <- res_steps5[[1]]$antibodies
  prot_no_ref <- no_refine[[1]]$antibodies
  
  diff0 <- sum(abs(as.matrix(prot_steps0) - as.matrix(prot_no_ref)))
  diff5 <- sum(abs(as.matrix(prot_steps5) - as.matrix(prot_no_ref)))
  expect_true(diff0 < diff5,
              info="Refinement steps=0 should yield prototypes closer to no_refine than steps=5.")
})

test_that("honeycombHIVE fails gracefully with invalid refineLoss for given task", {
  data(iris)
  X <- as.matrix(iris[,1:4])
  y <- iris$Species
  
  expect_error({
    honeycombHIVE(
      X = X,
      y = y,
      task = "classification",
      layers = 1,
      refine = TRUE,
      refineLoss = "huber",  
      refineSteps = 3,
      refineLR = 0.01,
      verbose = FALSE
    )
  }, regexp="For classification, loss must be one of|.*huber.*not.*supported")
})

test_that("honeycombHIVE clustering works with refine=TRUE, default refineLoss=mae", {
  data(iris)
  X <- as.matrix(iris[, 1:4])
  
  # refine=FALSE => normal honeycombHIVE
  res_no_ref <- honeycombHIVE(
    X = X,
    task = "clustering",
    layers = 2,
    nAntibodies = 8,
    maxIter = 5,
    refine = FALSE,
    verbose = FALSE
  )
  
  # refine=TRUE => add refine step
  res_ref <- honeycombHIVE(
    X = X,
    task = "clustering",
    layers = 2,
    nAntibodies = 8,
    maxIter = 5,
    refine = TRUE,            
    refineLoss = "mae",       
    refineSteps = 3,
    refineLR = 0.01,
    verbose = FALSE
  )
  
  # Check structure
  expect_length(res_ref, 2)  # 2 layers
  for (layer_obj in res_ref) {
    expect_true("antibodies" %in% names(layer_obj))
    expect_true("assignments" %in% names(layer_obj))
    expect_true("membership"  %in% names(layer_obj))
    expect_equal(layer_obj$task, "clustering")
  }
  
  # Compare final prototypes 
  final_no_ref <- res_no_ref[[2]]$antibodies
  final_ref    <- res_ref[[2]]$antibodies
  # Usually not identical
  expect_false(isTRUE(all.equal(final_no_ref, final_ref)))
})

test_that("honeycombHIVE classification: refine=TRUE with cross_entropy", {
  data(iris)
  X <- as.matrix(iris[, 1:4])
  y <- iris$Species
  
  # Non-refined
  res_no_ref <- honeycombHIVE(
    X = X,
    y = y,
    task = "classification",
    layers = 2,
    nAntibodies = 10,
    maxIter = 5,
    refine = FALSE,
    verbose = FALSE
  )
  
  # Refined with cross-entropy
  res_ref <- honeycombHIVE(
    X = X,
    y = y,
    task = "classification",
    layers = 2,
    nAntibodies = 10,
    maxIter = 5,
    refine = TRUE,
    refineLoss = "categorical_crossentropy",
    refineSteps = 4,
    refineLR = 0.01,
    refinePushAway = TRUE,
    verbose = FALSE
  )
  
  # Structure checks
  expect_length(res_ref, 2)
  for (ly in res_ref) {
    expect_equal(ly$task, "classification")
    expect_true(is.data.frame(ly$antibodies))
    # We also expect "assignments", "predictions", "membership"
    expect_true("assignments" %in% names(ly))
    expect_true("membership"  %in% names(ly))
  }
  
  # Compare final layer prototypes
  final_no_ref <- res_no_ref[[2]]$antibodies
  final_ref    <- res_ref[[2]]$antibodies
  expect_false(isTRUE(all.equal(final_no_ref, final_ref)))
})

test_that("honeycombHIVE regression: refine=TRUE with mse refineLoss", {
  if (!requireNamespace("MASS", quietly=TRUE)) {
    install.packages("MASS")
  }
  library(MASS)
  data(Boston)
  
  X <- as.matrix(Boston[, -14])
  y <- Boston$medv
  
  # Non-refined 
  res_no_ref <- honeycombHIVE(
    X = X,
    y = y,
    task = "regression",
    layers = 2,
    nAntibodies = 10,
    maxIter = 5,
    refine = FALSE,
    verbose = FALSE
  )
  
  # Refined
  res_ref <- honeycombHIVE(
    X = X,
    y = y,
    task = "regression",
    layers = 2,
    nAntibodies = 10,
    maxIter = 5,
    refine = TRUE,
    refineLoss = "mse",
    refineSteps = 3,
    refineLR = 0.005,
    verbose = FALSE
  )
  
  # Check structure
  expect_length(res_ref, 2)
  for (ly in res_ref) {
    expect_true("antibodies"  %in% names(ly))
    expect_true("assignments" %in% names(ly))   
    expect_true("predictions" %in% names(ly))   
    expect_true("membership"  %in% names(ly))
    expect_equal(ly$task, "regression")
  }
  
  # Compare final prototypes
  final_no_ref <- res_no_ref[[2]]$antibodies
  final_ref    <- res_ref[[2]]$antibodies
  expect_false(isTRUE(all.equal(final_no_ref, final_ref)))
  
})

test_that("honeycombHIVE refineSteps=0 yields same result as refine=FALSE", {
  data(iris)
  X <- as.matrix(iris[,1:4])
  
  # refine=FALSE
  res_no_ref <- honeycombHIVE(
    X = X,
    task = "clustering",
    layers = 1,
    nAntibodies = 6,
    maxIter = 3,
    refine = FALSE,
    verbose = FALSE
  )
  
  res_steps0 <- honeycombHIVE(
    X = X,
    task = "clustering",
    layers = 1,
    nAntibodies = 6,
    maxIter = 3,
    refine = TRUE,
    refineLoss = "mse",   
    refineSteps = 0,
    refineLR = 0.01,
    verbose = FALSE
  )
  
  # Compare final prototypes
  final_no_ref <- res_no_ref[[1]]$antibodies
  final_steps0 <- res_steps0[[1]]$antibodies
  
  # If you do random seeds, you can set a seed to ensure identical. If not, do a tolerance check:
  expect_equal(final_no_ref, final_steps0, tolerance=1e-14)
})

test_that("honeycombHIVE refineLoss mismatch: attempts cross_entropy in regression => should error or fallback", {
  if (!requireNamespace("MASS", quietly=TRUE)) {
    install.packages("MASS")
  }
  library(MASS)
  data(Boston)
  X <- as.matrix(Boston[, -14])
  y <- Boston$medv
  
  expect_error(
    honeycombHIVE(
      X = X,
      y = y,
      task = "regression",
      layers = 1,
      nAntibodies = 5,
      maxIter = 2,
      refine = TRUE,
      refineLoss = "categorical_crossentropy",
      refineSteps = 2,
      refineLR = 0.01,
      verbose = FALSE
    ),
    regexp="some.*error.*cross_entropy.*not.*valid.*regression"
  )
})
