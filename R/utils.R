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
.update_prototype <- function(ab_vec, 
                              x_i, 
                              same_label = NA, 
                              task = c("clustering","classification","regression"),
                              loss = c("categorical_crossentropy", "binary_crossentropy", 
                                       "kullback_leibler", "cosine",
                                       "mse","mae","poisson","huber"),
                              push_away = TRUE,
                              huber_delta = 1.0) {
  task <- match.arg(task)
  loss <- match.arg(loss)
  
  # Precompute
  d <- length(ab_vec)
  zero_vec <- numeric(d)          # a zero vector for quick returns
  pull_vec <- x_i - ab_vec        # pulling ab_vec toward x_i
  push_vec <- -pull_vec           # pushing ab_vec away from x_i
  
  # =================
  # 1) CLUSTERING
  # =================
  if (task == "clustering") {
    return(pull_vec)
  }
  
  # =================
  # 2) CLASSIFICATION
  # =================
  if (task == "classification") {
    # If no assigned label => no movement
    if (is.na(same_label)) {
      return(zero_vec)
    }
    is_same <- isTRUE(same_label)
    
    # (A) If 'is_same' => we PULL
    # (B) Else if push_away => we PUSH
    # (C) else => zero
    if (loss %in% c("categorical_crossentropy","binary_crossentropy","kullback_leibler","cosine","mse")) {
      # These classification losses share the same naive push/pull approach:
      if (is_same) {
        return(pull_vec)
      } else if (push_away) {
        return(push_vec)
      } else {
        return(zero_vec)
      }
    }
    
    else if (loss == "mae") {
      # sign-based approach
      if (is_same) {
        return(sign(pull_vec))
      } else if (push_away) {
        return(sign(push_vec))
      } else {
        return(zero_vec)
      }
    }
    
    # If user tries "poisson" or "huber" for classification, we might do a fallback or a no-op:
    return(zero_vec)
  }
  
  # =================
  # 3) REGRESSION
  # =================
  if (task == "regression") {
    # Residual in feature space
    residual <- pull_vec  # i.e., x_i - ab_vec
    
    if (loss == "mse") {
      # standard: partial derivative of 0.5 * ||residual||^2 => residual
      return(residual)
    }
    else if (loss == "mae") {
      # partial derivative => sign of residual
      return(sign(residual))
    }
    else if (loss == "poisson") {
      # In a pure param sense => grad ~ (pred - y)/pred * x_i, but we don't store "pred" or "y" separately.
      return(residual)
    }
    else if (loss == "huber") {
      # Huber: if ||residual|| <= huber_delta => MSE region; else => L1 region
      dist_val <- sqrt(sum(residual^2))
      if (dist_val <= huber_delta) {
        # MSE-like
        return(residual)
      } else {
        # linear region => scale by huber_delta
        # we do sign-based in each dimension, so:
        return(huber_delta * sign(residual))
      }
    }
    # Fallback if user gave a classification-based loss for regression => naive MSE approach
    return(residual)
  }
  
  # If something unexpected
  return(zero_vec)
}
