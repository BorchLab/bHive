# test script for visualizeHIVE.R - testcases are NOT comprehensive!
library(testthat)
library(ggplot2)

### SETUP FAKE DATA AND RESULT OBJECTS ###
# Create a numeric matrix X with 10 rows and 3 columns
set.seed(123)
X_mat <- matrix(rnorm(30), nrow = 10, ncol = 3)
colnames(X_mat) <- c("Feature1", "Feature2", "Feature3")

fake_result_single <- list(
  antibodies = matrix(rnorm(6), nrow = 2, ncol = 3),
  assignments = factor(sample(1:2, 10, replace = TRUE)),  # used for boxplot/violin/density
  membership = factor(sample(1:3, 10, replace = TRUE))    # used for scatterplots in clustering
)

# Fake result for regression (grouping is forced to "All")
fake_result_regression <- list(
  antibodies = matrix(rnorm(6), nrow = 2, ncol = 3),
  assignments = factor(sample(1:2, 10, replace = TRUE))
)

# Fake multi-layer result: a list of two layers where each element has an antibodies field.
fake_multilayer <- list(
  list(
    antibodies = matrix(rnorm(6), nrow = 2, ncol = 3),
    assignments = factor(sample(1:2, 10, replace = TRUE)),
    membership = factor(sample(1:2, 10, replace = TRUE))
  ),
  list(
    antibodies = matrix(rnorm(6), nrow = 2, ncol = 3),
    assignments = factor(sample(1:2, 10, replace = TRUE)),
    membership = factor(sample(1:2, 10, replace = TRUE))
  )
)


test_that("Scatter plot returns a ggplot object for clustering", {
  p <- visualizeHIVE(result = fake_result_single,
                     X = X_mat,
                     plot_type = "scatter",
                     title = "Test Scatter",
                     layer = 1,
                     task = "clustering",
                     transformation_method = "none",
                     transform = FALSE)
  expect_s3_class(p, "ggplot")
  expect_true(grepl("Test Scatter", p$labels$title))
})

test_that("Scatter plot returns a ggplot object for regression", {
  p <- visualizeHIVE(result = fake_result_regression,
                     X = X_mat,
                     plot_type = "scatter",
                     title = "Regression Scatter",
                     layer = 1,
                     task = "regression",
                     transformation_method = "PCA",
                     transform = TRUE)
  expect_s3_class(p, "ggplot")
  expect_true(grepl("Regression Scatter", p$labels$title))
})

test_that("Boxplot returns a ggplot object for classification", {
  p <- visualizeHIVE(result = fake_result_single,
                     X = as.data.frame(X_mat),
                     plot_type = "boxplot",
                     feature = "Feature2",
                     title = "Boxplot Test",
                     layer = 1,
                     task = "classification")
  expect_s3_class(p, "ggplot")
  expect_true(grepl("Boxplot Test", p$labels$title))
})

test_that("Boxplot returns a ggplot object for regression", {
  p <- visualizeHIVE(result = fake_result_regression,
                     X = as.data.frame(X_mat),
                     plot_type = "boxplot",
                     feature = 2,
                     title = "Boxplot Regression",
                     layer = 1,
                     task = "regression")
  expect_s3_class(p, "ggplot")
})

test_that("Violin plot returns a ggplot object", {
  p <- visualizeHIVE(result = fake_result_single,
                     X = as.data.frame(X_mat),
                     plot_type = "violin",
                     feature = "Feature3",
                     title = "Violin Plot",
                     layer = 1,
                     task = "classification")
  expect_s3_class(p, "ggplot")
})

test_that("Density plot returns a ggplot object", {
  p <- visualizeHIVE(result = fake_result_single,
                     X = as.data.frame(X_mat),
                     plot_type = "density",
                     feature = "Feature1",
                     title = "Density Plot",
                     layer = 1,
                     task = "classification")
  expect_s3_class(p, "ggplot")
})

test_that("Error is thrown if layer is non-numeric", {
  expect_error(visualizeHIVE(result = fake_result_single,
                             X = X_mat,
                             plot_type = "scatter",
                             title = "Layer error",
                             layer = "one",
                             task = "clustering"),
               "'layer' must be numeric")
})


test_that("Multi-layer scatter plot includes facet wrapping", {
  p <- visualizeHIVE(result = fake_multilayer,
                     X = X_mat,
                     plot_type = "scatter",
                     title = "Multi-layer Scatter",
                     layer = c(1, 2),
                     task = "clustering",
                     transformation_method = "none",
                     transform = FALSE)
  expect_s3_class(p, "ggplot")
  expect_true(inherits(p$facet, "FacetWrap"))
})

