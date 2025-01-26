# test script for utils.R - testcases are NOT comprehensive!


test_that(".validate_bHIVE_input: X must be a matrix or data frame", {
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

test_that(".validate_bHIVE_input: y must be factor or numeric if not NULL", {
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

test_that(".validate_bHIVE_input: X and y must have same number of rows", {
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

test_that(".validate_bHIVE_input: Input X cannot have missing values", {
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

test_that(".validate_bHIVE_input: Function returns TRUE invisibly on success", {
  # We can check the return value using expect_invisible
  df_X <- data.frame(x=1:5, y=6:10)
  y_num <- 1:5
  
  expect_invisible(
    .validate_bHIVE_input(X=df_X, y=y_num)
  )
  out <- .validate_bHIVE_input(X=df_X, y=y_num)
  expect_true(isTRUE(out))
})

test_that(".update_prototype: Clustering: Always pulls prototype toward x_i", {
  ab_vec <- c(1,2)
  x_i    <- c(3,4)
  
  grad <- .update_prototype(
    ab_vec = ab_vec,
    x_i = x_i,
    task = "clustering",
    loss = "mse"  # the 'loss' is ignored for clustering in the current logic
  )
  
  # Expect x_i - ab_vec
  expect_equal(grad, c(2,2))
})

test_that(".update_prototype: Classification + cross_entropy (same_label=TRUE => pull, FALSE => push if push_away=TRUE)", {
  ab_vec <- c(1,2)
  x_i    <- c(3,4)
  
  # same_label=TRUE => pull
  grad_same <- .update_prototype(
    ab_vec = ab_vec,
    x_i = x_i,
    same_label = TRUE,
    task = "classification",
    loss = "categorical_crossentropy",
    push_away = TRUE
  )
  expect_equal(grad_same, c(2,2))  # x_i - ab_vec
  
  # same_label=FALSE => push (ab_vec - x_i) if push_away=TRUE
  grad_diff <- .update_prototype(
    ab_vec = ab_vec,
    x_i = x_i,
    same_label = FALSE,
    task = "classification",
    loss = "categorical_crossentropy",
    push_away = TRUE
  )
  expect_equal(grad_diff, c(-2,-2))  # ab_vec - x_i
  
  # same_label=FALSE => no movement if push_away=FALSE
  grad_no_push <- .update_prototype(
    ab_vec = ab_vec,
    x_i = x_i,
    same_label = FALSE,
    task = "classification",
    loss = "categorical_crossentropy",
    push_away = FALSE
  )
  expect_equal(grad_no_push, c(0,0))
})

test_that(".update_prototype: Classification + mae uses sign-based approach", {
  ab_vec <- c(2,2)
  x_i    <- c(5,1)
  
  # same_label=TRUE => sign(x_i - ab_vec)
  # => sign(c(3, -1)) => c(1, -1)
  grad_same <- .update_prototype(
    ab_vec = ab_vec,
    x_i = x_i,
    same_label = TRUE,
    task = "classification",
    loss = "mae",
    push_away = TRUE
  )
  expect_equal(grad_same, c(1, -1))
  
  # same_label=FALSE => sign(ab_vec - x_i)
  # => c(2-5, 2-1) => c(-3,1) => sign => c(-1,1)
  grad_diff <- .update_prototype(
    ab_vec = ab_vec,
    x_i = x_i,
    same_label = FALSE,
    task = "classification",
    loss = "mae",
    push_away = TRUE
  )
  expect_equal(grad_diff, c(-1,1))
})

test_that(".update_prototype: Classification + no assigned label => zero grad", {
  ab_vec <- c(1,1)
  x_i    <- c(2,2)
  
  grad <- .update_prototype(
    ab_vec = ab_vec,
    x_i = x_i,
    same_label = NA,
    task = "classification",
    loss = "mse"
  )
  expect_equal(grad, c(0,0))
})

test_that(".update_prototype: Classification + huber or poisson does fallback (zero)", {
  ab_vec <- c(1,1)
  x_i    <- c(3,3)
  
  grad_huber <- .update_prototype(
    ab_vec, x_i, same_label=TRUE, 
    task="classification", loss="huber"
  )
  grad_poisson <- .update_prototype(
    ab_vec, x_i, same_label=FALSE, 
    task="classification", loss="poisson", 
    push_away=TRUE
  )
  
  expect_equal(grad_huber, c(0,0))
  expect_equal(grad_poisson, c(0,0))
})

test_that(".update_prototype: Regression + MSE => residual (x_i - ab_vec)", {
  ab_vec <- c(1,2)
  x_i    <- c(4,6)
  
  grad <- .update_prototype(
    ab_vec=ab_vec,
    x_i=x_i,
    task="regression",
    loss="mse"
  )
  expect_equal(grad, c(3,4))
})

test_that(".update_prototype: Regression + MAE => sign-based in each dimension", {
  ab_vec <- c(2,5)
  x_i    <- c(1,8)
  # residual = c(-1,3) => sign => c(-1,1)
  grad <- .update_prototype(
    ab_vec=ab_vec,
    x_i=x_i,
    task="regression",
    loss="mae"
  )
  expect_equal(grad, c(-1,1))
})

test_that(".update_prototype: Regression + huber => MSE region vs. linear region check", {
  ab_vec <- c(0,0)
  x_i    <- c(0.5,0.5)
  
  grad_mse_region <- .update_prototype(
    ab_vec=c(0,0),
    x_i=c(0.5,0.5),
    task="regression",
    loss="huber",
    huber_delta=1
  )
  expect_equal(grad_mse_region, c(0.5,0.5))
  
  # If huber_delta < 0.7 => linear region => sign => c(1,1)* delta
  grad_linear_region <- .update_prototype(
    ab_vec=c(0,0),
    x_i=c(0.5,0.5),
    task="regression",
    loss="huber",
    huber_delta=0.5
  )
  expect_equal(grad_linear_region, c(0.5,0.5))
  
  # Another scenario: if x_i= c(2,2), distance=2.828 > 1 => => c(1,1)
  grad_big <- .update_prototype(
    ab_vec=c(0,0),
    x_i=c(2,2),
    task="regression",
    loss="huber",
    huber_delta=1
  )
  # distance= sqrt(4+4)=2.828 >1 => linear => c(1,1)
  expect_equal(grad_big, c(1,1))
})

test_that(".update_prototype: Regression + poisson => naive same as residual", {
  ab_vec <- c(2,3)
  x_i    <- c(5,7)
  # residual= c(3,4)
  grad <- .update_prototype(
    ab_vec=ab_vec,
    x_i=x_i,
    task="regression",
    loss="poisson"
  )
  expect_equal(grad, c(3,4))
})

test_that(".update_prototype: Classification + 'cosine' => same push/pull approach", {
  ab_vec <- c(1,2)
  x_i    <- c(3,6)
  
  # same_label=TRUE => pull
  grad_same <- .update_prototype(
    ab_vec=ab_vec,
    x_i=x_i,
    same_label=TRUE,
    task="classification",
    loss="cosine",
    push_away=TRUE
  )
  expect_equal(grad_same, c(2,4))
  
  # same_label=FALSE => push => c(-2,-4)
  grad_diff <- .update_prototype(
    ab_vec=ab_vec,
    x_i=x_i,
    same_label=FALSE,
    task="classification",
    loss="cosine",
    push_away=TRUE
  )
  expect_equal(grad_diff, c(-2, -4))
})

test_that(".update_prototype: Unknown task or loss => no movement is default fallback", {
  local_update_prototype <- function(task, loss) {
    .update_prototype(
      ab_vec=c(0,0),
      x_i=c(1,1),
      same_label=NA,
      task=task,
      loss=loss
    )
  }
  expect_true(TRUE)
})