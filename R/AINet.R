#' @title AINet
#' @description R6 implementation of the Artificial Immune Network algorithm.
#' This is the core bHIVE algorithm using C++ backends for performance-critical
#' operations. Supports composable modules for somatic hypermutation, idiotypic
#' network regulation, germinal center selection, and more.
#'
#' @examples
#' # Clustering with Iris data
#' data(iris)
#' X <- as.matrix(iris[, 1:4])
#' model <- AINet$new(nAntibodies = 15, maxIter = 10, verbose = FALSE)
#' model$fit(X, task = "clustering")
#' table(model$result$assignments)
#'
#' # Classification
#' model2 <- AINet$new(nAntibodies = 20, maxIter = 10, verbose = FALSE)
#' model2$fit(X, iris$Species, task = "classification")
#' mean(model2$result$assignments == as.character(iris$Species))
#'
#' # Predict on new data
#' preds <- model2$predict(X[1:10, ])
#'
#' @importFrom R6 R6Class
#' @importFrom stats rnorm runif sd quantile
#' @export
AINet <- R6::R6Class(
  "AINet",
  inherit = ImmuneAlgorithm,
  public = list(

    #' @description Create a new AINet algorithm instance.
    #' @param nAntibodies Integer. Initial antibody population size.
    #' @param beta Numeric. Clone multiplier.
    #' @param epsilon Numeric. Suppression distance threshold.
    #' @param maxIter Integer. Maximum iterations.
    #' @param k Integer. Top-k antibodies to clone per data point.
    #' @param affinityFunc Character. Affinity function name.
    #' @param distFunc Character. Distance function name.
    #' @param affinityParams List. Parameters for affinity/distance functions.
    #' @param mutationDecay Numeric. Per-iteration mutation rate decay.
    #' @param mutationMin Numeric. Minimum mutation rate.
    #' @param maxClones Numeric. Maximum clones per antibody.
    #' @param stopTolerance Numeric. Early stopping tolerance.
    #' @param noImprovementLimit Integer. Early stopping patience.
    #' @param initMethod Character. Initialization method.
    #' @param shm An SHMEngine instance or NULL for default uniform mutation.
    #' @param init A VDJLibrary instance or NULL for default initialization.
    #' @param activation An ActivationGate instance or NULL.
    #' @param idiotypic An IdiotypicNetwork instance or NULL.
    #' @param germinalCenter A GerminalCenter instance or NULL.
    #' @param microenvironment A Microenvironment instance or NULL.
    #' @param memory A MemoryPool instance or NULL.
    #' @param classSwitcher A ClassSwitcher instance or NULL.
    #' @param consolidate Logical. For clustering, run Lloyd-style
    #'   consolidation (an M-step) after affinity maturation so antibodies are
    #'   pulled onto the data manifold and become true data-space prototypes.
    #'   Has no effect on classification. Default TRUE.
    #' @param consolidationSteps Integer. Maximum consolidation iterations.
    #' @param scale Character. Per-feature input scaling applied at \code{fit()}
    #'   and re-applied to new data at \code{predict()}. One of \code{"none"}
    #'   (default, no transform), \code{"zscore"} (center/SD), \code{"robust"}
    #'   (median / IQR, outlier-tolerant), or \code{"arcsinh"} (inverse
    #'   hyperbolic sine with cofactor \code{scaleCofactor}, the standard
    #'   mass-cytometry transform). Because \code{epsilon}, mutation scale, and
    #'   all distances live in feature units, scaling makes the same defaults
    #'   behave consistently across datasets of different magnitude.
    #' @param scaleCofactor Numeric. Cofactor for \code{scale = "arcsinh"}
    #'   (\code{asinh(x / cofactor)}). Default 5 (CyTOF convention; use ~150 for
    #'   fluorescence flow).
    #' @param targetK Integer or NULL. If set, force the clustering solution to
    #'   exactly \code{targetK} clusters. Affinity maturation still discovers
    #'   where prototypes belong, but the final consolidation seeds a K-means
    #'   (Lloyd) refinement at exactly \code{targetK} centroids: surviving
    #'   antibodies are agglomerated (if more than K) or split with k-means++ (if
    #'   fewer than K) before refinement. This decouples the reported cluster
    #'   count from the emergent suppression dynamics. NULL (default) keeps the
    #'   emergent, self-selected K. Ignored for classification.
    #' @param epsilonQuantile Numeric in (0, 1) or NULL. If set, the suppression
    #'   threshold is recomputed each iteration as this quantile of the pairwise
    #'   distances among the current antibodies, making suppression scale-free
    #'   and adaptive instead of using the fixed \code{epsilon}. NULL (default)
    #'   uses the fixed \code{epsilon}.
    #' @param coverageBoost Logical. Clustering only. After maturation, find
    #'   data points that no surviving antibody covers well (max affinity in the
    #'   bottom \code{coverageQuantile} tail) and seed fresh antibodies there
    #'   with k-means++. Counters the clonal-selection bias toward dense regions,
    #'   which otherwise leaves rare populations unrepresented. Pairs naturally
    #'   with \code{targetK}: the extra seeds give the forced-K refinement
    #'   candidate prototypes for sparse populations. Default FALSE.
    #' @param coverageQuantile Numeric in (0, 1). Affinity-coverage tail that
    #'   defines "poorly covered" points for \code{coverageBoost}. Default 0.05.
    #' @param verbose Logical. Print progress.
    initialize = function(nAntibodies = 20,
                          beta = 5,
                          epsilon = 0.01,
                          maxIter = 50,
                          k = 3,
                          affinityFunc = "gaussian",
                          distFunc = "euclidean",
                          affinityParams = list(alpha = 1, c = 1, p = 2,
                                                Sigma = NULL),
                          mutationDecay = 1.0,
                          mutationMin = 0.01,
                          maxClones = Inf,
                          stopTolerance = 0.0,
                          noImprovementLimit = Inf,
                          initMethod = "sample",
                          consolidate = TRUE,
                          consolidationSteps = 10L,
                          scale = c("none", "zscore", "robust", "arcsinh"),
                          scaleCofactor = 5,
                          targetK = NULL,
                          epsilonQuantile = NULL,
                          coverageBoost = FALSE,
                          coverageQuantile = 0.05,
                          shm = NULL,
                          init = NULL,
                          activation = NULL,
                          idiotypic = NULL,
                          germinalCenter = NULL,
                          microenvironment = NULL,
                          memory = NULL,
                          classSwitcher = NULL,
                          verbose = TRUE) {

      # Validate numeric parameters
      stopifnot(
        "nAntibodies must be a positive integer" = is.numeric(nAntibodies) && nAntibodies >= 1,
        "beta must be positive" = is.numeric(beta) && beta > 0,
        "epsilon must be non-negative" = is.numeric(epsilon) && epsilon >= 0,
        "maxIter must be a positive integer" = is.numeric(maxIter) && maxIter >= 1,
        "k must be a positive integer" = is.numeric(k) && k >= 1,
        "mutationDecay must be in (0, 1]" = is.numeric(mutationDecay) && mutationDecay > 0 && mutationDecay <= 1,
        "mutationMin must be non-negative" = is.numeric(mutationMin) && mutationMin >= 0
      )

      affinityFunc <- match.arg(affinityFunc, c("gaussian", "laplace",
                                                "polynomial", "cosine", "hamming"))
      distFunc     <- match.arg(distFunc, c("euclidean", "manhattan",
                                            "minkowski", "cosine",
                                            "mahalanobis", "hamming"))
      scale        <- match.arg(scale)

      if (!is.null(targetK)) {
        stopifnot("targetK must be a positive integer" =
                    is.numeric(targetK) && length(targetK) == 1L && targetK >= 1)
        targetK <- as.integer(targetK)
      }
      if (!is.null(epsilonQuantile)) {
        stopifnot("epsilonQuantile must be in (0, 1)" =
                    is.numeric(epsilonQuantile) && length(epsilonQuantile) == 1L &&
                    epsilonQuantile > 0 && epsilonQuantile < 1)
      }

      # --- Affinity/distance metric guard ---
      # Clonal selection matures antibodies to maximize `affinityFunc`, but
      # cluster assignment uses `distFunc`. If the two disagree on geometry the
      # model optimizes one space and is read out in another. The worst case is
      # cosine affinity (scores angle only, leaves antibody magnitude a free
      # random walk) paired with a Euclidean-family distance: antibodies drift
      # far off the data manifold and every point collapses onto one of them.
      # Align the distance to the affinity's natural metric and warn rather than
      # silently returning a degenerate clustering.
      natural_dist <- c(gaussian = "euclidean", laplace = "manhattan",
                        polynomial = "euclidean", cosine = "cosine",
                        hamming = "hamming")[[affinityFunc]]
      cosine_mismatch <- (affinityFunc == "cosine") != (distFunc == "cosine")
      if (cosine_mismatch) {
        warning(sprintf(
          "affinityFunc='%s' and distFunc='%s' use inconsistent geometries; %s. Set distFunc='%s' to silence this.",
          affinityFunc, distFunc,
          sprintf("overriding distFunc to '%s'", natural_dist),
          natural_dist), call. = FALSE)
        distFunc <- natural_dist
      }

      config <- list(
        nAntibodies       = as.integer(nAntibodies),
        beta              = beta,
        epsilon           = epsilon,
        maxIter           = maxIter,
        k                 = k,
        affinityFunc      = affinityFunc,
        distFunc          = distFunc,
        affinityParams    = affinityParams,
        mutationDecay     = mutationDecay,
        mutationMin       = mutationMin,
        maxClones         = maxClones,
        stopTolerance     = stopTolerance,
        noImprovementLimit = noImprovementLimit,
        initMethod        = match.arg(initMethod, c("sample", "random",
                                                     "random_uniform", "kmeans++")),
        consolidate       = isTRUE(consolidate),
        consolidationSteps = as.integer(consolidationSteps),
        scale             = scale,
        scaleCofactor     = scaleCofactor,
        targetK           = targetK,
        epsilonQuantile   = epsilonQuantile,
        coverageBoost     = isTRUE(coverageBoost),
        coverageQuantile  = coverageQuantile,
        scaling           = NULL,  # populated at fit() with learned stats
        verbose           = verbose
      )

      modules <- list(
        shm              = shm,
        init             = init,
        activation       = activation,
        idiotypic        = idiotypic,
        germinalCenter   = germinalCenter,
        microenvironment = microenvironment,
        memory           = memory,
        classSwitcher    = classSwitcher
      )
      # Remove NULL modules
      modules <- modules[!vapply(modules, is.null, logical(1))]

      super$initialize(config = config, modules = modules)
    },

    #' @description Fit the AINet algorithm to data.
    #' @param X Numeric matrix or data frame (n x d).
    #' @param y Optional factor target for classification.
    #' @param task Character: "clustering" or "classification".
    #'   Inferred from y if NULL.
    #' @param ... Additional arguments (currently unused).
    #' @return Invisible self, with \code{result} populated.
    fit = function(X, y = NULL, task = NULL, ...) {
      # Validate inputs
      .validate_bHIVE_input(X, y)

      # Infer task
      if (is.null(task)) {
        task <- if (is.null(y)) "clustering" else "classification"
      }
      task <- match.arg(task, c("clustering", "classification"))

      X <- as.matrix(X)

      # ================================
      # 0. Input scaling (C)
      # ================================
      # Learn the transform on this training matrix and stash the stats so
      # predict() applies the identical map to new data. epsilon, mutation
      # scale and every distance are in feature units, so a fixed default
      # only behaves sensibly once features share a comparable scale.
      scaling <- private$.learn_scaling(X, self$config$scale,
                                        self$config$scaleCofactor)
      X <- private$.apply_scaling(X, scaling)
      self$config$scaling <- scaling

      n <- nrow(X)
      d <- ncol(X)
      cfg <- self$config

      # ================================
      # 1. Initialize antibody population
      # ================================
      A <- private$.initialize_antibodies(X, cfg$nAntibodies, cfg$initMethod,
                                          init_lib = self$modules$init)

      # Memory recall (clustering only): merge relevant prior memory cells
      # into the starting repertoire. Disabled for classification because
      # MemoryPool$recall returns cells without their original class labels.
      mem <- self$modules$memory
      if (task == "clustering" && !is.null(mem) && mem$size() > 0L) {
        recalled <- mem$recall(X,
                                affinityFunc = cfg$affinityFunc,
                                affinityParams = list(
                                  alpha = cfg$affinityParams$alpha %||% 1,
                                  c     = cfg$affinityParams$c %||% 1,
                                  p     = cfg$affinityParams$p %||% 2))
        if (nrow(recalled) > 0 && ncol(recalled) == d) {
          A <- rbind(A, recalled)
        }
      }

      self$repertoire <- ImmuneRepertoire$new(A)
      m <- self$repertoire$size()

      # ================================
      # 2. Task-specific setup
      # ================================
      task_int <- switch(task, clustering = 0L, classification = 1L)
      nClasses <- 0L
      classes  <- NULL
      if (task == "classification") {
        classes  <- levels(y)
        nClasses <- length(classes)
        y_num    <- as.numeric(y) - 1  # 0-indexed class
      } else {
        y_num <- rep(0, n)
      }

      # Affinity/distance params (base values; iter_alpha may be modulated
      # by ClassSwitcher below).
      # alpha = "auto" sets the RBF/Laplace bandwidth by the median heuristic:
      # alpha = 1 / median(||x_i - x_j||^2) over a subsample of the (scaled)
      # data, so the kernel resolves structure at the data's own length scale
      # instead of the fixed alpha = 1 that is arbitrary once features are
      # scaled. Only meaningful for distance-based kernels (gaussian/laplace).
      base_alpha <- cfg$affinityParams$alpha %||% 1
      if (is.character(base_alpha) && identical(base_alpha, "auto")) {
        base_alpha <- private$.auto_bandwidth(X, cfg$affinityFunc)
        if (cfg$verbose) {
          cat(sprintf("Auto bandwidth: alpha = %.5g\n", base_alpha))
        }
      }
      c_p   <- cfg$affinityParams$c %||% 1
      p_p   <- cfg$affinityParams$p %||% 2
      iter_alpha <- base_alpha
      Sigma_inv <- if (!is.null(cfg$affinityParams$Sigma)) {
        solve(cfg$affinityParams$Sigma)
      } else {
        matrix(0, 0, 0)
      }

      # SHM dispatch params -- forwarded to clonal_selection_iteration_cpp.
      # When no SHM module is supplied, fall back to "uniform" which
      # reproduces the legacy decay-based mutation behavior.
      shm_engine <- self$modules$shm
      shm_method <- if (is.null(shm_engine)) "uniform" else shm_engine$method
      shm_p <- if (is.null(shm_engine)) {
        list(c_rate = 1, temperature = 0.5, E_0 = 1, base_rate = 0.1,
             beta1 = 0.9, beta2 = 0.999, adam_epsilon = 1e-8)
      } else {
        shm_engine$params
      }

      # Adaptive SHM moment matrices. Only allocated when used.
      use_adaptive <- identical(shm_method, "adaptive")
      m1_state <- if (use_adaptive) matrix(0, m, d) else matrix(0, 0, 0)
      m2_state <- if (use_adaptive) matrix(0, m, d) else matrix(0, 0, 0)
      if (use_adaptive && is.function(shm_engine$init_state)) {
        shm_engine$init_state(m, d)
      }

      # Early stopping state
      noImproveCount <- 0
      prevCount <- m

      # Pre-allocate antibody_classes so gated-aside bookkeeping has a
      # value to index into on iter=1 (the original code computed this
      # only after clonal_selection_iteration_cpp inside the loop).
      antibody_classes <- if (task == "classification") {
        sample(classes, size = m, replace = TRUE)
      } else {
        rep(NA_character_, m)
      }

      # ================================
      # 3. Main iteration loop
      # ================================
      for (iter in seq_len(cfg$maxIter)) {
        A_current <- self$repertoire$as_matrix()
        m <- nrow(A_current)

        # (a0) ActivationGate: gate antibodies in over-dense neighborhoods
        # OUT of this round of clonal selection. They sit aside unchanged
        # while clonal selection runs on the sparse subset, then rejoin
        # the repertoire. Prevents runaway cloning into already-crowded
        # regions of feature space.
        gate <- self$modules$activation
        gated_aside_A      <- NULL
        gated_aside_classes <- NULL
        gated_aside_m1     <- NULL
        gated_aside_m2     <- NULL
        if (!is.null(gate) && m > 4L) {
          Ab_Ab_aff <- compute_affinity_matrix(A_current, A_current,
                                            cfg$affinityFunc,
                                            iter_alpha, c_p, p_p)
          diag(Ab_Ab_aff) <- 0
          density <- rowSums(Ab_Ab_aff)

          # threshold2 in [0, 1] is the density quantile above which an
          # antibody is gated (e.g. 0.75 = top quartile sits out).
          q_cut <- gate$threshold2 %||% 0.75
          q_cut <- max(0.5, min(0.95, q_cut))  # clamp to sensible range
          density_cut <- quantile(density, q_cut, na.rm = TRUE)
          gated_idx   <- which(density > density_cut)

          # Apply Signal 1 (affinity threshold) as an additional gate: an
          # antibody whose max affinity to any data point is below
          # threshold1 is also gated out (it isn't binding anything).
          if (!is.null(gate$threshold1) && gate$threshold1 > 0) {
            Ab_X_aff <- compute_affinity_matrix(X, A_current,
                                             cfg$affinityFunc,
                                             iter_alpha, c_p, p_p)
            max_aff_per_ab <- apply(Ab_X_aff, 2, max)
            low_aff_idx    <- which(max_aff_per_ab < gate$threshold1)
            gated_idx      <- union(gated_idx, low_aff_idx)
          }

          # Ensure at least 2 antibodies still enter clonal selection
          if (length(gated_idx) > 0 && (m - length(gated_idx)) >= 2L) {
            gated_aside_A <- A_current[gated_idx, , drop = FALSE]
            if (task == "classification") {
              gated_aside_classes <- antibody_classes[gated_idx]
            }
            if (use_adaptive) {
              gated_aside_m1 <- m1_state[gated_idx, , drop = FALSE]
              gated_aside_m2 <- m2_state[gated_idx, , drop = FALSE]
              m1_state <- m1_state[-gated_idx, , drop = FALSE]
              m2_state <- m2_state[-gated_idx, , drop = FALSE]
            }
            A_current <- A_current[-gated_idx, , drop = FALSE]
            if (task == "classification") {
              antibody_classes <- antibody_classes[-gated_idx]
            }
          }
        }

        # (a) Clonal selection + SHM-dispatched mutation [C++]
        cs_result <- clonal_selection_iteration_cpp(
          A_current, X, y_num, task_int, cfg$k, cfg$beta,
          cfg$maxClones, cfg$mutationDecay, cfg$mutationMin,
          iter, cfg$affinityFunc, iter_alpha, c_p, p_p, nClasses,
          shm_method,
          shm_p$c_rate, shm_p$temperature, shm_p$E_0, shm_p$base_rate,
          shm_p$beta1, shm_p$beta2, shm_p$adam_epsilon,
          m1_state, m2_state
        )
        if (use_adaptive) {
          m1_state <- cs_result$m1_state
          m2_state <- cs_result$m2_state
        }

        # Update labels (on the selected subset only, then rejoin gated)
        if (task == "classification") {
          new_classes <- apply(cs_result$class_counts, 1, function(row) {
            if (all(row == 0)) classes[sample(nClasses, 1)]
            else classes[which.max(row)]
          })
        }

        # Rejoin gated-aside antibodies to the post-selection repertoire.
        # In classification, refresh their class labels by majority-vote
        # of their nearest data points so stale random labels from the
        # pre-allocation don't poison final predictions.
        if (!is.null(gated_aside_A)) {
          self$repertoire$cells <- rbind(cs_result$A, gated_aside_A)
          if (task == "classification") {
            ga_aff <- compute_affinity_matrix(X, gated_aside_A,
                                                cfg$affinityFunc,
                                                iter_alpha, c_p, p_p)
            refreshed <- vapply(seq_len(nrow(gated_aside_A)), function(j) {
              top <- order(ga_aff[, j], decreasing = TRUE)[
                seq_len(min(cfg$k, nrow(ga_aff)))]
              tab <- table(y[top])
              names(tab)[which.max(tab)]
            }, character(1))
            antibody_classes <- c(new_classes, refreshed)
          }
          if (use_adaptive) {
            m1_state <- rbind(m1_state, gated_aside_m1)
            m2_state <- rbind(m2_state, gated_aside_m2)
          }
        } else {
          self$repertoire$cells <- cs_result$A
          if (task == "classification") {
            antibody_classes <- new_classes
          }
        }

        # (a1) GerminalCenter: Tfh-mediated quality selection. Probabilistic
        # survival weighted by task-aware quality score (clustering: average
        # affinity to assigned points; classification: majority-class purity).
        # Survivor indices are mirrored onto antibody_classes and SHM state.
        gc_mod <- self$modules$germinalCenter
        if (!is.null(gc_mod) && self$repertoire$size() > gc_mod$nTfh) {
          gc_surv <- gc_mod$select(
            self$repertoire, X, y, task,
            affinityFunc   = cfg$affinityFunc,
            affinityParams = list(alpha = iter_alpha, c = c_p, p = p_p))
          if (task == "classification") {
            antibody_classes <- antibody_classes[gc_surv]
          }
          if (use_adaptive) {
            m1_state <- m1_state[gc_surv, , drop = FALSE]
            m2_state <- m2_state[gc_surv, , drop = FALSE]
          }
        }

        # (a2) Microenvironment-aware mutation jitter [optional]
        # Density-dependent perturbation: antibodies in over-dense regions
        # of feature space get small jitter (stabilize / memory-like);
        # antibodies in sparse regions get large jitter (explore / push
        # outward), countering the clonal-selection drift toward the
        # data centroid. Class labels are preserved across the jitter.
        microenv <- self$modules$microenvironment
        env <- NULL
        if (!is.null(microenv) && self$repertoire$size() > 4L) {
          env <- microenv$assess(self$repertoire, X,
                                  affinityFunc = cfg$affinityFunc,
                                  affinityParams = list(alpha = iter_alpha,
                                                          c = c_p,
                                                          p = p_p))
          A_post   <- self$repertoire$as_matrix()
          x_sd     <- apply(X, 2, sd)
          decay    <- cfg$mutationDecay ^ max(iter - 1, 0)
          base_amp <- 0.005 * decay  # ~0.5% of feature SD on iter 1, decaying
          mods     <- env$mutation_modifiers
          d_cols   <- ncol(A_post)
          for (j in seq_len(nrow(A_post))) {
            if (mods[j] <= 0) next
            sigma <- base_amp * mods[j] * x_sd
            A_post[j, ] <- A_post[j, ] + rnorm(d_cols, 0, sigma)
          }
          self$repertoire$cells <- A_post
        }

        # (a3) ClassSwitcher: bind isotype to microenvironment zone and use
        # the population-mean per-isotype alpha for the NEXT iteration's
        # affinity calls. Requires Microenvironment to have run this iter
        # (otherwise we have no zones to switch on). Per-antibody alpha is
        # aggregated to a scalar since the C++ kernels take a scalar alpha.
        cs_mod <- self$modules$classSwitcher
        if (!is.null(cs_mod) && !is.null(env)) {
          alphas <- cs_mod$switch_isotypes(self$repertoire, env$zones)
          iter_alpha <- mean(alphas)
        }

        # (b) Idiotypic regulation [C++, optional]
        # Bell-curve Ab-Ab dynamics cull antibodies in over-crowded niches
        # (over-stimulation -> suppression) and isolated antibodies (under-
        # stimulation -> death), leaving a diversity-preserving repertoire.
        # Runs BEFORE epsilon-ball network suppression so the two operators
        # ablate independently. See IdiotypicNetwork for parameter semantics.
        idi <- self$modules$idiotypic
        if (!is.null(idi)) {
          idi_out <- idiotypic_dynamics_cpp(
            self$repertoire$as_matrix(),
            cfg$affinityFunc, iter_alpha, c_p, p_p,
            idi$theta_low, idi$theta_high,
            idi$source_rate, idi$decay_rate,
            idi$dt, as.integer(idi$timeSteps),
            idi$survival_threshold
          )
          surv_idx <- which(as.logical(idi_out$keep))

          # Safety net: if dynamics would kill every antibody (e.g. ill-tuned
          # thresholds for the current data scale), keep the top-population
          # antibodies so the iteration can continue and a downstream sweep
          # can still penalize this configuration via low Silhouette / kappa.
          if (length(surv_idx) == 0L) {
            pop <- as.numeric(idi_out$population)
            keep_n <- max(1L, floor(0.1 * length(pop)))
            surv_idx <- order(pop, decreasing = TRUE)[seq_len(keep_n)]
          }

          self$repertoire$subset(surv_idx)
          if (task == "classification") {
            antibody_classes <- antibody_classes[surv_idx]
          }
          if (use_adaptive) {
            m1_state <- m1_state[surv_idx, , drop = FALSE]
            m2_state <- m2_state[surv_idx, , drop = FALSE]
          }
        }

        # (c) Network suppression [C++]
        # Removes near-duplicate antibodies within an epsilon-ball in distFunc.
        # With epsilonQuantile set, the threshold tracks the current antibody
        # spread (a low quantile of their pairwise distances) so suppression is
        # scale-free and adapts as the repertoire contracts, rather than using a
        # fixed feature-unit epsilon that means different things per dataset.
        eps_iter <- cfg$epsilon
        if (!is.null(cfg$epsilonQuantile)) {
          eps_iter <- private$.adaptive_epsilon(
            self$repertoire$cells, cfg$distFunc, cfg$epsilonQuantile,
            p_p, Sigma_inv, fallback = cfg$epsilon)
        }
        keep <- network_suppression_cpp(
          self$repertoire$cells, cfg$distFunc, eps_iter,
          p_p, Sigma_inv
        )
        kept_idx <- which(keep)
        self$repertoire$subset(kept_idx)
        m_new <- self$repertoire$size()

        if (task == "classification") {
          antibody_classes <- antibody_classes[kept_idx]
        }
        if (use_adaptive) {
          m1_state <- m1_state[kept_idx, , drop = FALSE]
          m2_state <- m2_state[kept_idx, , drop = FALSE]
        }

        if (m_new == 0) {
          stop("All antibodies were suppressed. Increase nAntibodies or decrease epsilon.")
        }

        # Record iteration history
        self$history[[iter]] <- list(n_antibodies = m_new)

        # Early stopping
        changeCount <- abs(m_new - prevCount)
        if (changeCount <= cfg$stopTolerance) {
          noImproveCount <- noImproveCount + 1
        } else {
          noImproveCount <- 0
        }
        prevCount <- m_new

        if (noImproveCount >= cfg$noImprovementLimit) {
          if (cfg$verbose) {
            cat("Early stopping: no improvement for", noImproveCount, "iterations.\n")
          }
          break
        }

        if (cfg$verbose) {
          cat(sprintf("Iteration %d | #Antibodies: %d | noImproveCount: %d\n",
                      iter, m_new, noImproveCount))
        }
      }

      # ================================
      # 4. Final assignment [C++]
      # ================================
      A_final <- self$repertoire$as_matrix()
      m <- nrow(A_final)

      # ================================
      # 4a. Orphan-antibody pruning
      # ================================
      # Drop antibodies that are not the nearest neighbor to any training
      # point. These are surviving "ghost" antibodies (passed idiotypic
      # and network suppression but bind nothing) that inflate the
      # repertoire without contributing predictions. Pruning them stops
      # them from causing test-time mis-assignments and tightens the
      # cluster-count = effective-antibody-count relationship.
      fa_pre <- final_assignment_cpp(X, A_final, cfg$affinityFunc,
                                       cfg$distFunc,
                                       switch(task, clustering = 0L,
                                                       classification = 1L),
                                       iter_alpha, c_p, p_p, Sigma_inv)
      assigned_to <- if (task == "clustering") {
        as.integer(fa_pre$assignments)
      } else {
        as.integer(fa_pre$best_antibody_idx)
      }
      counts <- tabulate(assigned_to, nbins = m)
      orphan <- counts == 0L
      if (any(orphan) && sum(!orphan) >= 2L) {
        keep_idx <- which(!orphan)
        A_final  <- A_final[keep_idx, , drop = FALSE]
        if (task == "classification") {
          antibody_classes <- antibody_classes[keep_idx]
        }
        self$repertoire$subset(keep_idx)
        m <- nrow(A_final)
      }

      # ================================
      # 4a'. Coverage boost (D)
      # ================================
      # Clonal expansion follows density, so rare populations can end the run
      # with no nearby prototype. Find the worst-covered points (lowest max
      # affinity to any surviving antibody) and drop k-means++ seeds among them
      # so the final assignment / forced-K refinement has candidate prototypes
      # for those regions. Clustering only (new seeds carry no class label).
      if (task == "clustering" && isTRUE(cfg$coverageBoost) && nrow(A_final) >= 1L) {
        cov_aff <- compute_affinity_matrix(X, A_final, cfg$affinityFunc,
                                           iter_alpha, c_p, p_p)
        coverage <- apply(cov_aff, 1, max)
        thr <- stats::quantile(coverage, cfg$coverageQuantile, names = FALSE,
                               na.rm = TRUE)
        under <- which(coverage <= thr)
        if (length(under) >= 2L) {
          n_new <- min(length(under), max(2L, round(0.5 * cfg$nAntibodies)))
          new_seeds <- init_kmeanspp_cpp(X[under, , drop = FALSE], n_new)
          A_final <- rbind(A_final, new_seeds)
          self$repertoire$cells <- A_final
          if (cfg$verbose) {
            cat(sprintf("Coverage boost: +%d antibodies for %d under-covered points\n",
                        n_new, length(under)))
          }
        }
      }

      # ================================
      # 4b. Final assignment [C++]
      # ================================
      if (task == "clustering") {
        if (!is.null(cfg$targetK)) {
          # Target-K mode (A): the matured antibodies supply seeds, but the
          # reported partition is forced to exactly targetK via a seeded Lloyd
          # refinement. Decouples the cluster count from the emergent
          # suppression dynamics, which otherwise under-/over-cluster depending
          # on scale.
          fk <- private$.force_k(X, A_final, cfg$targetK, cfg,
                                 iter_alpha, c_p, p_p, Sigma_inv)
          A_final     <- fk$antibodies
          assignments <- fk$assignments
        } else if (cfg$consolidate && cfg$consolidationSteps > 0L &&
            nrow(A_final) >= 2L) {
          # Consolidation (M-step): pull antibodies onto the data manifold so
          # they are genuine data-space prototypes, not affinity-maximizing
          # directions that may sit far off the data (see metric guard).
          cons <- private$.consolidate_clusters(X, A_final, cfg,
                                                 iter_alpha, c_p, p_p, Sigma_inv)
          A_final     <- cons$antibodies
          assignments <- cons$assignments
          # Note: self$repertoire is left as the matured antibodies (with their
          # isotype/age/lineage metadata and memory state). Consolidation merges
          # clusters, so the centroids have no 1:1 identity with repertoire rows;
          # they are the reported prototypes (result$antibodies), while the
          # repertoire stays the biological population that memory archives.
        } else {
          fa <- final_assignment_cpp(X, A_final, cfg$affinityFunc, cfg$distFunc,
                                     0L, iter_alpha, c_p, p_p, Sigma_inv)
          assignments <- as.numeric(factor(fa$assignments))
        }
        self$result <- list(
          antibodies  = A_final,
          assignments = assignments,
          # Prototypes back-transformed to the original feature space (for
          # interpretation); identical to `antibodies` when scale = "none".
          antibodies_unscaled =
            private$.invert_scaling(A_final, cfg$scaling),
          task        = task
        )
      } else {
        fa <- final_assignment_cpp(X, A_final, cfg$affinityFunc, cfg$distFunc,
                                   1L, iter_alpha, c_p, p_p, Sigma_inv)
        assignments <- antibody_classes[fa$best_antibody_idx]
        self$result <- list(
          antibodies       = A_final,
          assignments      = assignments,
          antibody_classes = antibody_classes,
          antibodies_unscaled =
            private$.invert_scaling(A_final, cfg$scaling),
          task             = task
        )
      }

      # ================================
      # 5. Memory archive (post-training)
      # ================================
      # High-affinity antibodies become long-lived memory cells that
      # persist on the MemoryPool across fit() calls. For classification,
      # carry class labels in repertoire metadata so recall consumers can
      # use them later.
      if (!is.null(mem)) {
        if (task == "classification") {
          self$repertoire$metadata$class_label <- antibody_classes
        }
        mem$archive(self$repertoire, X,
                    affinityFunc   = cfg$affinityFunc,
                    affinityParams = list(alpha = iter_alpha,
                                          c     = c_p,
                                          p     = p_p))
      }

      invisible(self)
    }
  ),

  private = list(

    # --- Input scaling (C) --------------------------------------------------
    # Thin method wrappers around the file-level scaling helpers so the fit()
    # body reads cleanly; predict() (in the base class) calls the helpers
    # directly since it has no access to these private methods.
    .learn_scaling = function(X, method, cofactor) {
      .bhive_learn_scaling(X, method, cofactor)
    },
    .apply_scaling = function(X, scaling) {
      .bhive_apply_scaling(X, scaling)
    },
    .invert_scaling = function(A, scaling) {
      .bhive_invert_scaling(A, scaling)
    },

    # --- Median-heuristic RBF bandwidth (C) ---------------------------------
    # alpha = 1 / median(||x_i - x_j||^2). Distance-kernel only; identity for
    # cosine/polynomial/hamming where a Euclidean length scale is meaningless.
    .auto_bandwidth = function(X, affinityFunc) {
      if (!affinityFunc %in% c("gaussian", "laplace")) return(1)
      n <- nrow(X)
      idx <- if (n > 1000L) sample.int(n, 1000L) else seq_len(n)
      Xs <- X[idx, , drop = FALSE]
      d2 <- as.vector(stats::dist(Xs))^2
      d2 <- d2[d2 > 0]
      if (length(d2) == 0) return(1)
      med <- stats::median(d2)
      if (!is.finite(med) || med <= 0) return(1)
      1 / med
    },

    # --- Adaptive suppression threshold (C) ---------------------------------
    # Return a low quantile of the pairwise antibody distances so the epsilon
    # ball tracks the repertoire's own spread. Falls back to the fixed epsilon
    # when there are too few antibodies to form a distribution.
    .adaptive_epsilon = function(A, distFunc, q, p, Sigma_inv, fallback) {
      m <- nrow(A)
      if (is.null(m) || m < 3L) return(fallback)
      D <- compute_distance_matrix(A, A, distFunc, p, Sigma_inv)
      du <- D[upper.tri(D)]
      du <- du[is.finite(du) & du > 0]
      if (length(du) == 0) return(fallback)
      as.numeric(stats::quantile(du, probs = q, names = FALSE))
    },

    # --- Force exactly K clusters (A) ---------------------------------------
    # Use the matured antibodies as informed seeds, coerce to exactly K
    # centroids, then run Euclidean Lloyd to convergence. If the network kept
    # more than K prototypes, agglomerate the closest pairs (Ward on the
    # antibodies) down to K; if fewer, add k-means++ seeds drawn from the data
    # so under-clustering can be corrected. Euclidean refinement is used because
    # the arithmetic mean is the L2-optimal centroid (consistent, monotone),
    # independent of the training affinity which already drove the search.
    .force_k = function(X, A, K, cfg, alpha, c_p, p_p, Sigma_inv) {
      d <- ncol(X)
      n <- nrow(X)
      K <- min(K, n)  # cannot ask for more clusters than points
      m <- nrow(A)

      if (m > K) {
        # Agglomerate antibodies to K groups, seed = group means.
        hc  <- stats::hclust(stats::dist(A), method = "ward.D2")
        grp <- stats::cutree(hc, k = K)
        cent <- t(vapply(sort(unique(grp)), function(g)
          colMeans(A[grp == g, , drop = FALSE]), numeric(d)))
      } else if (m < K) {
        # Augment with k-means++ seeds from the data to reach K.
        extra <- init_kmeanspp_cpp(X, K)
        cent  <- rbind(A, extra)[seq_len(K), , drop = FALSE]
      } else {
        cent <- A
      }

      # Seed assignment by affinity argmax to the seed prototypes, then Lloyd.
      assign <- as.integer(final_assignment_cpp(
        X, cent, cfg$affinityFunc, cfg$distFunc, 1L,
        alpha, c_p, p_p, Sigma_inv)$best_antibody_idx)
      steps <- max(cfg$consolidationSteps, 10L)
      prev  <- NULL
      for (s in seq_len(steps)) {
        ks <- sort(unique(assign))
        cent <- t(vapply(ks, function(g)
          colMeans(X[assign == g, , drop = FALSE]), numeric(d)))
        new_assign <- as.integer(final_assignment_cpp(
          X, cent, cfg$affinityFunc, "euclidean", 0L,
          alpha, c_p, p_p, Sigma_inv)$assignments)
        if (!is.null(prev) && identical(new_assign, prev)) {
          assign <- new_assign; break
        }
        prev   <- new_assign
        assign <- new_assign
      }
      list(antibodies  = cent,
           assignments = as.numeric(factor(assign)))
    },

    # Lloyd-style consolidation of matured antibodies into data-space
    # prototypes. Affinity maturation (clonal selection + SHM) finds where the
    # prototypes should point, but under magnitude-blind affinities (cosine) it
    # leaves them off the data manifold, and even Euclidean affinities do not
    # guarantee a prototype equals the centroid of the points it wins. This
    # runs a few k-means-style refinement steps to fix that.
    #
    # The seed assignment is by AFFINITY (argmax), not distance: a Euclidean
    # seed over off-manifold antibodies can collapse every point onto a single
    # far-flung antibody, whereas the affinity seed preserves the partition the
    # network actually learned. After the first M-step the prototypes are
    # on-manifold, so subsequent E-steps use the (guard-consistent) distFunc.
    .consolidate_clusters = function(X, A, cfg, alpha, c_p, p_p, Sigma_inv) {
      d <- ncol(X)
      # Seed by AFFINITY (argmax) using the trained affinity: robust to
      # off-manifold antibodies, and it preserves the partition the network
      # learned (an affinity-blind Euclidean seed over far-flung antibodies can
      # collapse every point onto one of them).
      seed <- final_assignment_cpp(X, A, cfg$affinityFunc, cfg$distFunc, 1L,
                                   alpha, c_p, p_p, Sigma_inv)$best_antibody_idx
      assign <- as.integer(seed)
      # The refinement itself is Euclidean Lloyd: the arithmetic mean is the
      # centroid that minimizes within-cluster L2, so reassignment must also be
      # Euclidean for the steps to be consistent and monotone. This is
      # independent of the training affinity/distance -- the affinity drove the
      # search; consolidation commits the prototypes to the data manifold.
      cent <- A
      prev <- NULL
      for (s in seq_len(cfg$consolidationSteps)) {
        ks <- sort(unique(assign))
        # M-step: each prototype is the mean of the points assigned to it.
        cent <- t(vapply(ks, function(g)
          colMeans(X[assign == g, , drop = FALSE]), numeric(d)))
        # E-step: reassign points to the nearest consolidated prototype (L2).
        new_assign <- as.integer(final_assignment_cpp(
          X, cent, cfg$affinityFunc, "euclidean", 0L,
          alpha, c_p, p_p, Sigma_inv)$assignments)
        if (!is.null(prev) && identical(new_assign, prev)) {
          assign <- new_assign
          break
        }
        prev   <- new_assign
        assign <- new_assign
      }
      list(antibodies  = cent,
           assignments = as.numeric(factor(assign)))
    },

    .initialize_antibodies = function(X, nAntibodies, method, init_lib = NULL) {
      # When a VDJLibrary (or any object exposing $generate(n, X)) is supplied
      # via the `init` module, route initialization through V(D)J combinatorial
      # assembly. This produces a structured, diverse starting repertoire that
      # spans the data manifold rather than clumping near the centroid.
      if (!is.null(init_lib) && is.function(init_lib$generate)) {
        return(init_lib$generate(nAntibodies, X))
      }

      n <- nrow(X)
      d <- ncol(X)
      switch(
        method,
        "sample" = X[sample.int(n, size = nAntibodies, replace = TRUE), , drop = FALSE],
        "random" = {
          xMean <- colMeans(X)
          xSd   <- apply(X, 2, sd) + 1e-8
          mat   <- matrix(rnorm(nAntibodies * d), nrow = nAntibodies)
          mat   <- sweep(mat, 2, xSd, `*`)
          sweep(mat, 2, xMean, `+`)
        },
        "random_uniform" = {
          xMin <- apply(X, 2, min)
          xMax <- apply(X, 2, max)
          mat  <- matrix(runif(nAntibodies * d), nrow = nAntibodies)
          sweep(sweep(mat, 2, xMax - xMin, `*`), 2, xMin, `+`)
        },
        "kmeans++" = init_kmeanspp_cpp(X, nAntibodies)
      )
    }
  )
)


# ============================================================================
# Internal scaling helpers (shared by AINet$fit and ImmuneAlgorithm$predict)
# ============================================================================
# Kept as free functions, not R6 methods, so the base-class predict() can apply
# the identical transform to new data without reaching into AINet internals.

#' Learn a per-feature scaling from a training matrix
#' @param X numeric matrix
#' @param method one of "none","zscore","robust","arcsinh"
#' @param cofactor arcsinh cofactor
#' @return list(method, center, scale, cofactor) with stats in original units
#' @keywords internal
#' @noRd
.bhive_learn_scaling <- function(X, method = "none", cofactor = 5) {
  if (is.null(method) || method == "none") {
    return(list(method = "none"))
  }
  if (method == "arcsinh") {
    return(list(method = "arcsinh", cofactor = cofactor))
  }
  if (method == "zscore") {
    ctr <- colMeans(X, na.rm = TRUE)
    scl <- apply(X, 2, stats::sd, na.rm = TRUE)
  } else if (method == "robust") {
    ctr <- apply(X, 2, stats::median, na.rm = TRUE)
    scl <- apply(X, 2, stats::IQR, na.rm = TRUE)
  } else {
    stop("Unknown scaling method: ", method)
  }
  # Guard against zero-variance columns (constant features) -> divide by 1.
  scl[!is.finite(scl) | scl <= 0] <- 1
  list(method = method, center = ctr, scale = scl)
}

#' Apply a learned scaling to a matrix
#' @keywords internal
#' @noRd
.bhive_apply_scaling <- function(X, scaling) {
  if (is.null(scaling) || scaling$method == "none") return(X)
  if (scaling$method == "arcsinh") return(asinh(X / scaling$cofactor))
  sweep(sweep(X, 2, scaling$center, `-`), 2, scaling$scale, `/`)
}

#' Invert a learned scaling (map prototypes back to original units)
#' @keywords internal
#' @noRd
.bhive_invert_scaling <- function(A, scaling) {
  if (is.null(scaling) || scaling$method == "none") return(A)
  if (scaling$method == "arcsinh") return(sinh(A) * scaling$cofactor)
  sweep(sweep(A, 2, scaling$scale, `*`), 2, scaling$center, `+`)
}
