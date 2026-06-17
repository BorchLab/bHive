// [[Rcpp::depends(RcppArmadillo)]]
#include "affinity_distance.h"
#include "shm.h"

// Core clonal selection + mutation loop for one iteration
//
// For each data point:
//   1. Compute affinity to all antibodies (uses pre-computed affinity matrix)
//   2. Identify top-k antibodies
//   3. Clone and mutate top-k via the SHM strategy `shm_method`, keeping
//      improvements. Adaptive method threads m1_state/m2_state in/out.
//   4. Accumulate task-specific statistics
//
// Returns a list with:
//   - A: updated antibody matrix
//   - class_counts: (classification) weighted vote matrix
//   - m1_state, m2_state: updated moment matrices (adaptive SHM)
//
// [[Rcpp::export]]
Rcpp::List clonal_selection_iteration_cpp(
    arma::mat A,                     // passed by value (modified in place)
    const arma::mat& X,
    const arma::vec& y_num,          // 0-indexed class int (classification)
    int task_int,                    // 0=clustering, 1=classification
    int k,
    double beta,
    double maxClones,
    double mutationDecay,
    double mutationMin,
    int iter,
    const std::string& affinity_type,
    double alpha,
    double c_param,
    double p_param,
    int nClasses,
    const std::string& shm_method,   // SHM strategy
    double shm_c_rate,
    double shm_temperature,
    double shm_E_0,
    double shm_base_rate,
    double shm_beta1,
    double shm_beta2,
    double shm_adam_epsilon,
    arma::mat m1_state,              // adaptive: (m x d) first moments
    arma::mat m2_state) {            // adaptive: (m x d) second moments

  const arma::uword n = X.n_rows;
  const arma::uword m = A.n_rows;
  const arma::uword d = X.n_cols;
  const TaskType task = static_cast<TaskType>(task_int);
  const AffinityType aff_type = parse_affinity_type(affinity_type);

  // Adaptive method requires (m x d) state. If caller passed empty
  // matrices (e.g. SHM not adaptive), allocate zero state lazily so the
  // dispatcher does not segfault on empty access.
  const bool is_adaptive = (shm_method == "adaptive");
  if (is_adaptive) {
    if (m1_state.n_rows != m || m1_state.n_cols != d) {
      m1_state = arma::zeros<arma::mat>(m, d);
    }
    if (m2_state.n_rows != m || m2_state.n_cols != d) {
      m2_state = arma::zeros<arma::mat>(m, d);
    }
  }

  // Task-specific accumulators
  arma::mat class_counts;
  if (task == TaskType::CLASSIFICATION) {
    class_counts = arma::zeros<arma::mat>(m, nClasses);
  }

  // Process each data point. We recompute each point's 1 x m affinity row
  // against the current (mutating) antibody set on the fly, instead of
  // precomputing the full n x m matrix once and patching a column on every
  // accepted mutation. The old patch-on-mutation scheme nested an O(n) column
  // recompute inside this per-point loop, so with acceptances scaling in n the
  // routine was O(n^2) in the number of cells -- the dominant cost at
  // CyTOF/scRNA scale. The per-row recompute is the same O(n*m*d) total work
  // and yields identical values (a point always sees every mutation made by
  // earlier points in the pass), but keeps the routine linear in n.
  for (arma::uword i = 0; i < n; ++i) {
    const arma::rowvec aff_values =
        affinity_matrix_cpp(X.rows(i, i), A, aff_type, alpha, c_param, p_param).row(0);

    const double max_aff = aff_values.max();
    if (max_aff <= 0.0 || !std::isfinite(max_aff)) continue;

    // Top-k indices (descending by affinity)
    const arma::uword k2 = std::min(static_cast<arma::uword>(k), m);
    const arma::uvec sorted_idx = arma::sort_index(aff_values, "descend");
    const arma::uvec top_idx = sorted_idx.head(k2);

    // Task-specific accumulation
    if (task == TaskType::CLASSIFICATION) {
      const int class_col = static_cast<int>(y_num(i));
      for (arma::uword ki = 0; ki < k2; ++ki) {
        const arma::uword jj = top_idx(ki);
        class_counts(jj, class_col) += aff_values(jj);
      }
    }

    // Clone and mutate
    const arma::rowvec x_i = X.row(i);
    for (arma::uword ki = 0; ki < k2; ++ki) {
      const arma::uword jj = top_idx(ki);
      const double f_j = aff_values(jj);
      const int nClonesInt = std::min(
        static_cast<int>(maxClones),
        static_cast<int>(std::floor(beta * (f_j / max_aff)))
      );
      if (nClonesInt <= 0) continue;

      for (int clone_id = 0; clone_id < nClonesInt; ++clone_id) {
        // Dispatch mutation by SHM method. mutationDecay and mutationMin
        // are reused for uniform-style decay; other methods use their own
        // parameters via shm_*.
        arma::rowvec mutated = mutate_by_method(
          shm_method, A.row(jj), x_i, f_j, iter,
          mutationDecay, mutationMin,
          shm_c_rate, shm_temperature, shm_E_0, shm_base_rate,
          shm_beta1, shm_beta2, shm_adam_epsilon,
          m1_state, m2_state, jj
        );

        // Use scalar affinity (no matrix allocation in hot path)
        const double f_mutated = affinity_scalar_cpp(
          x_i, mutated, aff_type, alpha, c_param, p_param
        );

        if (std::isfinite(f_mutated) && f_mutated > f_j) {
          // Accept the improved antibody in place. No affinity bookkeeping is
          // needed here: the next point recomputes its own affinity row against
          // this updated A (see the per-point recompute above), so the mutation
          // is seen immediately and exactly, with no O(n) column patch.
          A.row(jj) = mutated;
        }
      }
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("A") = A,
    Rcpp::Named("class_counts") = class_counts,
    Rcpp::Named("m1_state") = m1_state,
    Rcpp::Named("m2_state") = m2_state
  );
}


// Final assignment computation
// For clustering: assign each point to nearest antibody by distance
// For classification: assign each point to antibody with highest affinity
// [[Rcpp::export]]
Rcpp::List final_assignment_cpp(
    const arma::mat& X,
    const arma::mat& A,
    const std::string& affinity_type,
    const std::string& dist_type,
    int task_int,
    double alpha,
    double c_param,
    double p_param,
    const arma::mat& Sigma_inv) {

  const arma::uword n = X.n_rows;
  const TaskType task = static_cast<TaskType>(task_int);

  if (task == TaskType::CLUSTERING) {
    const DistanceType dtype = parse_distance_type(dist_type);
    const arma::mat D = distance_matrix_cpp(X, A, dtype, p_param, Sigma_inv);
    arma::uvec assignments(n);
    for (arma::uword i = 0; i < n; ++i) {
      assignments(i) = D.row(i).index_min();
    }
    const arma::ivec result = arma::conv_to<arma::ivec>::from(assignments) + 1;
    return Rcpp::List::create(Rcpp::Named("assignments") = result);

  } else {
    const AffinityType atype = parse_affinity_type(affinity_type);
    const arma::mat aff = affinity_matrix_cpp(X, A, atype, alpha, c_param, p_param);
    arma::uvec best_idx(n);
    for (arma::uword i = 0; i < n; ++i) {
      best_idx(i) = aff.row(i).index_max();
    }
    const arma::ivec result = arma::conv_to<arma::ivec>::from(best_idx) + 1;
    return Rcpp::List::create(Rcpp::Named("best_antibody_idx") = result);
  }
}


// kmeans++ initialization in C++
// [[Rcpp::export]]
arma::mat init_kmeanspp_cpp(const arma::mat& X, int nCenters) {
  const arma::uword n = X.n_rows;
  const arma::uword d = X.n_cols;
  arma::mat centers(nCenters, d);

  // Choose first center uniformly at random (safe int sampling)
  arma::uword idx = static_cast<arma::uword>(std::floor(R::runif(0.0, static_cast<double>(n))));
  if (idx >= n) idx = n - 1;  // guard against runif returning exactly n
  centers.row(0) = X.row(idx);

  if (nCenters > 1) {
    arma::vec dists(n);

    for (int cId = 1; cId < nCenters; ++cId) {
      // Compute min squared distance to any chosen center
      const arma::mat chosen = centers.rows(0, cId - 1);
      arma::mat D = distance_matrix_cpp(X, chosen, DistanceType::EUCLIDEAN,
                                        2.0, arma::mat());
      D = D % D;  // square distances

      // Min distance to any center for each point
      dists = arma::min(D, 1);

      // Sample proportional to squared distance
      const double total = arma::accu(dists);
      if (total <= 0.0) {
        idx = static_cast<arma::uword>(std::floor(R::runif(0.0, static_cast<double>(n))));
        if (idx >= n) idx = n - 1;
      } else {
        const arma::vec probs = dists / total;
        const double r = R::runif(0.0, 1.0);
        double cumsum = 0.0;
        idx = 0;
        for (arma::uword i = 0; i < n; ++i) {
          cumsum += probs(i);
          if (r <= cumsum) {
            idx = i;
            break;
          }
          idx = i;  // fallback: last element if rounding prevents reaching 1.0
        }
      }
      centers.row(cId) = X.row(idx);
    }
  }

  return centers;
}
