# test script for visualizeHIVE.R - testcases are NOT comprehensive!

dummy_class_result <- list(
  list(
    antibodies = matrix(c(5, 3, 1, 
                          6, 4, 2, 
                          7, 5, 3), nrow = 3, byrow = TRUE),
    assignments = c("setosa", "versicolor", "virginica", "setosa", "versicolor"),
    membership = c(1, 2, 3, 1, 2),
    task = "classification",
    predictions = c("setosa", "versicolor", "virginica", "setosa", "versicolor")
  ),
  list(
    antibodies = matrix(c(5.5, 3.2, 1.1, 
                          6.1, 3.6, 2.1, 
                          7.2, 3.8, 2.5), nrow = 3, byrow = TRUE),
    assignments = c("setosa", "versicolor", "virginica", "setosa", "versicolor"),
    membership = c(1, 2, 3, 1, 2),
    task = "classification",
    predictions = c("setosa", "versicolor", "virginica", "setosa", "versicolor")
  )
)

dummy_X <- iris[, 1:4]  # use the iris data for testing

### Dummy Regression Result
set.seed(123)
dummy_reg_result <- list(
  list(
    antibodies = matrix(rnorm(12), nrow = 3, ncol = 4),
    assignments = rep("All", 50),
    membership = sample(1:3, 50, replace = TRUE),
    task = "regression"
  )
)

test_that("visualizeHIVE returns a ggplot for classification scatter plot", {
  p <- visualizeHIVE(
    result = dummy_class_result,
    X = dummy_X,
    plot_type = "scatter",
    layer = 1,
    task = "classification",
    transform = FALSE
  )
  expect_s3_class(p, "ggplot")
})

test_that("visualizeHIVE facets when multiple layers provided", {
  p <- visualizeHIVE(
    result = dummy_class_result,
    X = dummy_X,
    plot_type = "scatter",
    layer = c(1, 2),
    task = "classification",
    transform = FALSE
  )
  # p$facet should be an object of class FacetWrap if facet_wrap was used.
  expect_true("FacetWrap" %in% class(p$facet))
})


test_that("visualizeHIVE errors when X is missing for boxplot/violin/density", {
  expect_error(
    visualizeHIVE(
      result = dummy_class_result,
      X = NULL,
      plot_type = "boxplot",
      feature = "Sepal.Width",
      layer = 1,
      task = "classification"
    ),
    "'data' must be of a vector type, was 'NULL'"
  )
})
