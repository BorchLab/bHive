# test script for utils.R - testcases are NOT comprehensive!


test_that("X must be a matrix or data frame", {
  # Non-matrix non-data.frame => expect error
  bad_X <- list(a=1, b=2)  # just a list
  expect_error(
    .validate_bHIVE_input(X=bad_X),
    "Input X must be a matrix or data frame."
  )
  
  # Data frame => no error
  df_X <- data.frame(a=1:3, b=2:4)
  expect_silent(
    .validate_bHIVE_input(X=df_X)
  )
  
  # Matrix => no error
  mat_X <- matrix(1:6, nrow=3, ncol=2)
  expect_silent(
    .validate_bHIVE_input(X=mat_X)
  )
})

test_that("y must be factor or numeric if not NULL", {
  df_X <- data.frame(a=1:3, b=2:4)
  
  # y as character => error
  y_char <- c("cat","dog","bird")
  expect_error(
    .validate_bHIVE_input(X=df_X, y=y_char),
    "y must be a factor.*numeric"
  )
  
  # y as factor => OK
  y_factor <- factor(c("A","B","C"))
  expect_silent(
    .validate_bHIVE_input(X=df_X, y=y_factor)
  )
  
  # y as numeric => OK
  y_num <- c(10,20,30)
  expect_silent(
    .validate_bHIVE_input(X=df_X, y=y_num)
  )
})

test_that("X and y must have same number of rows", {
  df_X <- data.frame(a=1:3, b=2:4)
  
  # y with different length => error
  y_mismatch <- c(1,2,3,4)
  expect_error(
    .validate_bHIVE_input(X=df_X, y=y_mismatch),
    "X and y must have the same number of rows."
  )
  
  # matching length => OK
  y_ok <- c(1,2,3)
  expect_silent(
    .validate_bHIVE_input(X=df_X, y=y_ok)
  )
})

test_that("Input X cannot have missing values", {
  df_X <- data.frame(a=c(NA,2,3), b=c(4,5,6))
  
  # X has NA => error
  expect_error(
    .validate_bHIVE_input(X=df_X),
    "Input X contains missing values.*"
  )
  
  # If no missing => OK
  no_na <- data.frame(a=1:3, b=4:6)
  expect_silent(
    .validate_bHIVE_input(X=no_na)
  )
})

test_that("Function returns TRUE invisibly on success", {
  # We can check the return value using expect_invisible
  df_X <- data.frame(x=1:5, y=6:10)
  y_num <- 1:5
  
  expect_invisible(
    .validate_bHIVE_input(X=df_X, y=y_num)
  )
  out <- .validate_bHIVE_input(X=df_X, y=y_num)
  expect_true(isTRUE(out))
})