.validate_bHIVE_input <- function(X, 
                                  y = NULL) {
  if (!is.matrix(X) && !is.data.frame(X)) {
    stop("Input X must be a matrix or data frame.")
  }
  if (!is.null(y)) {
    if (!is.factor(y) & !is.numeric(y)) {
      stop("y must be a factor (classification) or numeric vector (regression).")
    }
    if (nrow(X) != length(y)) {
      stop("X and y must have the same number of rows.")
    }
  }
  if (anyNA(X)) {
    stop("Input X contains missing values. Please handle missing values before running bHIVE.")
  }
  invisible(TRUE)
}

#Determines how to move a single prototype 'ab_vec' in feature space
.update_prototype <- function(ab_vec, x_i, 
                              same_label = NA, 
                              task = c("clustering","classification","regression"),
                              loss = c("categorical_crossentropy", "binary_crossentropy", 
                                       "kullback_leibler", "cosine",
                                       "mse","mae","poisson","huber"),
                              push_away = TRUE,
                              huber_delta = 1.0) {
  
  task <- match.arg(task)
  loss <- match.arg(loss)
  
  # If you want to handle "clustering" directly, do so:
  if (task == "clustering") {
    # Typically always pull the prototype toward x_i
    return(x_i - ab_vec)
  }
  
  # Classification logic
  if (task == "classification") {
    # If no assigned label => no movement
    if (is.na(same_label)) {
      return(rep(0, length(ab_vec)))
    }
    # same_label => TRUE or FALSE
    is_same <- isTRUE(same_label)
    
    # We'll define a generic "pull" = x_i - ab_vec, "push" = ab_vec - x_i
    # Then interpret the loss in a naive push/pull style.
    
    pull_vec <- (x_i - ab_vec)
    push_vec <- (ab_vec - x_i)
    
    if (loss %in% c("categorical_crossentropy","binary_crossentropy","kullback_leibler")) {
      # Typically "cross_entropy" or "kl" => push/pull
      if (is_same) {
        return(pull_vec)  # pull
      } else if (push_away) {
        return(push_vec)  # push
      } else {
        return(rep(0, length(ab_vec)))
      }
    }
    else if (loss == "cosine") {
      if (is_same) {
        return(pull_vec)
      } else if (push_away) {
        return(push_vec)
      } else {
        return(rep(0, length(ab_vec)))
      }
    }
    else if (loss == "mse") {
      if (is_same) {
        return(pull_vec)
      } else if (push_away) {
        return(push_vec)
      } else {
        return(rep(0, length(ab_vec)))
      }
    }
    else if (loss == "mae") {
      # sign-based. If same => sign(pull_vec), else => sign(push_vec) if push_away
      if (is_same) {
        return(sign(pull_vec))
      } else if (push_away) {
        return(sign(push_vec))
      } else {
        return(rep(0, length(ab_vec)))
      }
    }
    return(rep(0, length(ab_vec)))
  }
  
  # Regression logic
  if (task == "regression") {
    residual <- (x_i - ab_vec)
    
    if (loss == "mse") {
      # standard pull
      return(residual)
    }
    else if (loss == "mae") {
      # sign-based pull
      return(sign(residual))
    }
    else if (loss == "poisson") {
      return(residual)
    }
    else if (loss == "huber") {
      # If distance <= huber_delta => MSE region, else => linear region
      dist_val <- sqrt(sum(residual^2))  # L2 norm
      if (dist_val <= huber_delta) {
        # MSE-like
        return(residual)
      } else {
        # linear region
        return(huber_delta * sign(residual))
      }
    }
    return(residual)
  }
  
  # default fallback
  return(rep(0, length(ab_vec)))
}
