# test script for visualizeHIVE.R - testcases are NOT comprehensive!

# Dummy input feature matrix (6 observations, 5 features)
dummy_X <- matrix(rnorm(30), nrow = 6, ncol = 5)
colnames(dummy_X) <- paste0("V", 1:5)

# Dummy multilayer result for clustering:
dummy_layer1 <- list(
  antibodies = matrix(1:10, nrow = 2, ncol = 5), # 2 prototypes, 5 features
  membership = c(1, 1, 2, 2, 1, 2)
)
dummy_layer2 <- list(
  antibodies = matrix(11:20, nrow = 2, ncol = 5),
  membership = c(2, 2, 1, 1, 2, 1)
)
dummy_result_cluster <- list(dummy_layer1, dummy_layer2)

# Dummy multilayer result for classification:
dummy_layer_class <- list(
  antibodies = matrix(1:10, nrow = 2, ncol = 5),
  membership = c(1, 1, 2, 2, 1, 2),
  predictions = factor(c("A", "A", "B", "B", "A", "B"))
)
dummy_result_class <- list(dummy_layer_class)

# Dummy multilayer result for regression:
dummy_layer_reg <- list(
  antibodies = matrix(1:10, nrow = 2, ncol = 5),
  membership = c(1, 1, 2, 2, 1, 2),
  predictions = c(1.1, 1.2, 2.1, 2.0, 1.3, 2.2)
)
dummy_result_reg <- list(dummy_layer_reg)

# Dummy result with no membership/predictions.
dummy_layer_noGroup <- list(
  antibodies = matrix(1:10, nrow = 2, ncol = 5)
)
dummy_result_noGroup <- list(dummy_layer_noGroup)

test_that("visualizeHIVE returns a ggplot for scatter plots (clustering)", {
  p <- visualizeHIVE(result = dummy_result_cluster,
                     X = dummy_X,
                     plot_type = "scatter",
                     layer = 1,
                     task = "clustering")
  expect_true(inherits(p, "ggplot"))
})

test_that("visualizeHIVE returns a ggplot for scatter plots (classification)", {
  p <- visualizeHIVE(result = dummy_result_class,
                     X = dummy_X,
                     plot_type = "scatter",
                     layer = 1,
                     task = "classification")
  expect_true(inherits(p, "ggplot"))
})

test_that("visualizeHIVE returns a ggplot for scatter plots (regression)", {
  p <- visualizeHIVE(result = dummy_result_reg,
                     X = dummy_X,
                     plot_type = "scatter",
                     layer = 1,
                     task = "regression")
  expect_true(inherits(p, "ggplot"))
})

test_that("visualizeHIVE defaults to layer 1 when layer not specified", {
  p <- visualizeHIVE(result = dummy_result_cluster,
                     X = dummy_X,
                     plot_type = "scatter",
                     task = "clustering")
  expect_true(inherits(p, "ggplot"))
})

test_that("visualizeHIVE facets when multiple layers are specified", {
  p <- visualizeHIVE(result = dummy_result_cluster,
                     X = dummy_X,
                     plot_type = "scatter",
                     layer = c(1,2),
                     task = "clustering")
  # Check that the ggplot object has a facet_wrap (class "FacetWrap")
  expect_true("FacetWrap" %in% class(p$facet))
})

test_that("visualizeHIVE errors for boxplot-type plots when X is missing", {
  expect_error(visualizeHIVE(result = dummy_result_cluster,
                             plot_type = "boxplot",
                             layer = 1,
                             task = "clustering"),
               "X must be provided")
})

test_that("visualizeHIVE returns a ggplot for boxplot (classification)", {
  p <- visualizeHIVE(result = dummy_result_class,
                     X = dummy_X,
                     plot_type = "boxplot",
                     feature = "V3",
                     layer = 1,
                     task = "classification")
  expect_true(inherits(p, "ggplot"))
})

test_that("visualizeHIVE returns a ggplot for violin plot (regression)", {
  p <- visualizeHIVE(result = dummy_result_reg,
                     X = dummy_X,
                     plot_type = "violin",
                     feature = "V3",
                     layer = 1,
                     task = "regression")
  expect_true(inherits(p, "ggplot"))
})

test_that("visualizeHIVE returns a ggplot for density plot (clustering)", {
  p <- visualizeHIVE(result = dummy_result_cluster,
                     X = dummy_X,
                     plot_type = "density",
                     feature = "V4",
                     layer = 1,
                     task = "clustering")
  expect_true(inherits(p, "ggplot"))
})

test_that("visualizeHIVE handles missing membership/predictions gracefully", {
  p <- visualizeHIVE(result = dummy_result_noGroup,
                     X = dummy_X,
                     plot_type = "scatter",
                     layer = 1,
                     task = "clustering")
  expect_true(inherits(p, "ggplot"))
})

test_that("visualizeHIVE errors with invalid plot_type", {
  expect_error(visualizeHIVE(result = dummy_result_cluster,
                             X = dummy_X,
                             plot_type = "invalid",
                             layer = 1,
                             task = "clustering"),
               "must be one of")
})

test_that("visualizeHIVE transformation works with PCA", {
  p <- visualizeHIVE(result = dummy_result_cluster,
                     X = dummy_X,
                     plot_type = "scatter",
                     layer = 1,
                     task = "clustering",
                     transform = TRUE,
                     transformation_method = "PCA")
  expect_true(inherits(p, "ggplot"))
})

test_that("visualizeHIVE transformation with tSNE works for scatter plot", {
  p <- visualizeHIVE(result = dummy_result_cluster,
                     X = dummy_X,
                     plot_type = "scatter",
                     layer = 1,
                     task = "clustering",
                     transform = TRUE,
                     transformation_method = "tSNE")
  expect_true(inherits(p, "ggplot"))
})

test_that("visualizeHIVE Transformation with UMAP works for scatter plot", {
  p <- visualizeHIVE(result = dummy_result_cluster,
                     X = dummy_X,
                     plot_type = "scatter",
                     layer = 1,
                     task = "clustering",
                     transform = TRUE,
                     transformation_method = "UMAP")
  expect_true(inherits(p, "ggplot"))
})

