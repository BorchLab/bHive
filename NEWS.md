# bHIVE 0.99.6

## New Behavior
* **Unified engine.** `bHIVE()` is now a thin wrapper over the `AINet` R6
  engine instead of a separate pure-R implementation. `swarmbHIVE()` and
  `honeycombHIVE()` call `bHIVE()`, so they now inherit the C++ clonal-selection
  and suppression backends, Lloyd consolidation, input scaling, target-K, and
  all composable immunology modules. The old module-free R loop is retired. The
  return value gains `antibodies_unscaled` and (from `bHIVE()`) a `model`
  handle; the historical `antibodies` / `assignments` / `task` fields are
  unchanged. `swarmbHIVE()` now forwards every grid column to `bHIVE()`, so any
  `bHIVE`/`AINet` argument can be tuned, not just `nAntibodies`/`beta`/`epsilon`.

* **Input scaling.** `AINet$new()` and `bHIVE()` gain a `scale` argument
  (`"none"`, `"zscore"`, `"robust"`, `"arcsinh"`). The transform is learned at
  `fit()` and re-applied to new data at `predict()`. Because `epsilon`,
  mutation scale, and every distance live in feature units, scaling makes the
  same defaults behave consistently across datasets of different magnitude.
  `arcsinh` (cofactor `scaleCofactor`, default 5) is the standard mass-cytometry
  transform. `affinityParams$alpha = "auto"` sets the RBF/Laplace bandwidth by
  the median heuristic (alpha = 1 / median pairwise squared distance).

* **Target-K clustering.** New `targetK` argument forces the clustering result
  to exactly K clusters. Affinity maturation still discovers where prototypes
  belong, but a seeded K-means refinement coerces the count: surviving
  antibodies are agglomerated (if more than K) or split with k-means++ (if
  fewer). This decouples the reported cluster count from the scale-sensitive
  emergent suppression dynamics.

* **Scale-free suppression.** New `epsilonQuantile` argument recomputes the
  suppression threshold each iteration as a quantile of pairwise antibody
  distances, adapting to the repertoire's own spread instead of a fixed
  feature-unit `epsilon`.

* **Rare-population protection.** New `coverageBoost` (with `coverageQuantile`)
  seeds fresh antibodies among the worst-covered data points after maturation,
  countering the clonal-selection bias toward dense regions that otherwise
  leaves rare populations unrepresented. `honeycombHIVE()` gains
  `smallClusterAction` (`"merge"`, `"keep"`, `"drop"`): small clusters are now
  merged into the nearest prototype (no observations dropped) or kept as rare
  prototypes by default, instead of being silently deleted.

# bHIVE 0.99.5

## Performance
* `clonal_selection_iteration_cpp` no longer recomputes the full n-length
  affinity column on every accepted mutation. That recompute was nested inside
  the per-point loop, so with acceptances scaling in the number of points the
  routine ran in O(n^2) in the cell count, the dominant cost at CyTOF and
  scRNA scale. Each point now computes its own 1 x m affinity row against the
  current antibody set. This is the same O(n*m*d) work per iteration but
  linear in n, and it is numerically identical because a point still sees
  every mutation made by earlier points in the pass. Clustering and
  classification fits drop from roughly 101s to 1s at n=2000, and per-fit cost
  now scales linearly (a 16,000-cell fit runs in about 8s instead of an
  extrapolated 1.8 hours).

## New Behavior
* Affinity/distance metric guard in `AINet$new()`. Clonal selection matures
  antibodies to maximize `affinityFunc`, but cluster assignment uses
  `distFunc`. When the two used different geometries the model optimized one
  space and read out in another. A cosine versus non-cosine mismatch now warns
  and aligns the distance to the affinity's natural metric (cosine to cosine,
  otherwise the Euclidean family). Without it, cosine affinity scores angle
  only and never constrains antibody magnitude, so antibodies drifted far off
  the data manifold and every point collapsed onto a single antibody (one
  cluster, ARI 0).
* Clustering consolidation, controlled by `consolidate` (default `TRUE`) and
  `consolidationSteps` (default 10). After affinity maturation, `AINet$fit()`
  seeds an assignment by affinity, which stays robust to off-manifold
  antibodies, then runs Euclidean Lloyd refinement so prototypes become true
  data-space centroids rather than affinity-maximizing directions that may sit
  far off the data. The repertoire keeps the matured antibodies and their
  metadata, so memory archiving and isotype state are unchanged. Only
  `result$antibodies` reports the consolidated prototypes. Set
  `consolidate = FALSE` to recover the raw affinity-assignment behavior. The
  metric guard and default-on consolidation change clustering output, so
  cached benchmarks should be regenerated.

# bHIVE 0.99.4

## Bug Fixes
* Fix issue with Windows linker failure (Build R-release/devel/oldrel for Windows)
* Fixed R-devel/Linux arm64 test failure (Error in [.default(tb, cl, cl): subscript out of bounds)

# bHIVE 0.99.3

## Module Integration

Modules previously accepted by `AINet$new()` were stored on the algorithm
but never invoked by the fit loop. v0.99.3 wires every module into
`AINet$fit()` so it affects predictions:

* `init` (`VDJLibrary` or any object exposing `$generate(n, X)`) now
  controls initial-repertoire construction. When supplied, it overrides
  `initMethod`. This makes V(D)J combinatorial assembly the actual
  starting state rather than dead metadata.
* `idiotypic` (`IdiotypicNetwork`) is now invoked between clonal selection
  and network suppression. The bell-curve `idiotypic_dynamics_cpp` kernel
  culls antibodies in over-clumped niches (suppression) and isolated
  antibodies (insufficient stimulation). A safety net keeps the top-10%
  population when ill-tuned thresholds would kill every antibody, so the
  iteration can continue and downstream selection can penalize the config.
* `activation` (`ActivationGate`) now pre-screens antibodies in over-dense
  feature-space neighborhoods *out* of each round of clonal selection.
  Gated antibodies sit aside unchanged while selection runs on the sparse
  subset, then rejoin the repertoire. In classification mode their class
  labels are refreshed by majority-vote of their nearest data points so
  pre-allocated random labels do not poison final predictions.
* `microenvironment` (`Microenvironment`) is assessed each iteration after
  clonal selection. Its `mutation_modifiers` scale a per-antibody Gaussian
  jitter applied in feature units (`0.5%` of feature SD on iter 1,
  decaying with `mutationDecay`). Sparse-region antibodies explore
  outward, dense-region antibodies stabilize.
* `germinalCenter` (`GerminalCenter`) now runs Tfh-mediated quality
  selection after clonal-selection rejoin. Survivor indices are composed
  across selection rounds and mirrored onto `antibody_classes` and the
  adaptive-SHM moment matrices so external per-antibody state stays
  aligned. `GerminalCenter$select()` now returns the composed survivor
  vector (relative to the input repertoire) and exposes it as
  `$last_survivors`.
* `classSwitcher` (`ClassSwitcher`) runs after `microenvironment` and
  binds each antibody's isotype (IgM/IgG/IgA) to its zone. The
  population-mean per-isotype alpha is then used as the *next* iteration's
  affinity kernel width, modulating matching breadth based on the current
  density landscape. Requires `microenvironment` to be present.
* `memory` (`MemoryPool`) is invoked at two points. Pre-iteration, in
  clustering mode, relevant memory cells from prior fits are recalled
  (via `recall(X)`) and merged into the starting repertoire so prior
  knowledge seeds the search. Post-iteration, after orphan pruning, the
  surviving high-affinity antibodies are archived back into the pool
  (`archive()`) so they persist across fits. Classification archives also
  attach `class_label` to repertoire metadata so future consumers can read
  them back; recall is suppressed for classification to avoid label loss.
* `shm` (`SHMEngine`) is now dispatched inside `clonal_selection_iteration_cpp`
  rather than ignored. All five strategies — `uniform`, `airs`, `hotspot`,
  `energy`, `adaptive` — produce distinct mutation behavior. The adaptive
  strategy's Adam-style first/second moment matrices are threaded through
  the C++ kernel and re-aligned in R after every subset/reorder operation
  (activation gating, germinal-center selection, idiotypic culling,
  network suppression) so per-antibody moment tracking stays consistent
  across iterations.

## C++ Backend

* `clonal_selection_iteration_cpp` now accepts SHM dispatch parameters
  (`shm_method`, `shm_c_rate`, `shm_temperature`, `shm_E_0`, `shm_base_rate`,
  `shm_beta1`, `shm_beta2`, `shm_adam_epsilon`) plus adaptive moment
  matrices `m1_state`/`m2_state`. The R wrapper provides defaults
  (`shm_method = "uniform"`, empty state matrices) so legacy callers that
  pass the original 15 arguments still get the previous behavior.
* New `src/shm.h` exposes the five SHM strategy helpers and a
  `mutate_by_method` dispatcher shared by `shm_mutate_cpp` and
  `clonal_selection_iteration_cpp`. The helpers are no longer `static`.

## New Behavior

* Orphan-antibody pruning: `AINet$fit()` now drops surviving antibodies
  that are not the nearest neighbor to any training point before computing
  the final assignment. Eliminates ghost antibodies that pass suppression
  but bind nothing, and tightens the relationship between cluster count
  and effective antibody count.
* `antibody_classes` is now pre-allocated before the iteration loop in
  classification mode (uniform random over `levels(y)`), so the new
  pre-selection gate logic has a value to index into on iteration 1
  rather than referencing a not-yet-defined variable.

# bHIVE 0.99.2

## Breaking Changes
* Removed regression task support across the package. The previous regression
  implementation was unreliable and is being redesigned. `bHIVE()`,
  `honeycombHIVE()`, `swarmbHIVE()`, `refineB()`, `AINet`, `GerminalCenter`,
  `bHIVEmodel`, `honeycombHIVEmodel`, and `visualizeHIVE()` no longer accept
  `task = "regression"` or numeric `y`.
* Removed regression-only loss functions (`mse`, `huber`, `poisson`) and
  regression metrics (`rmse`, `mae`, `r2`) from `swarmbHIVE()`.
* Removed `refineHuberDelta` parameter from `honeycombHIVE()` and
  `honeycombHIVEmodel`.
* `final_assignment_cpp()` no longer accepts `antibody_values` or
  `overall_mean` arguments.

# bHIVE 0.99.1

## C++ Backend
* Added RcppArmadillo backend for BLAS-optimized bulk affinity and distance
  matrix computation, replacing per-element R loops
* Scalar affinity function for clone/mutate hot path avoids 1x1 matrix
  allocation overhead
* C++ implementations for clonal selection iteration, network suppression,
  kmeans++ initialization, final assignment, somatic hypermutation (5 methods),
  and idiotypic network dynamics

## R6 Class Architecture
* New `ImmuneRepertoire` class for antibody collections with metadata tracking
  (isotype, state, age, lineage)
* New `ImmuneAlgorithm` abstract base class with fit/predict/summary interface
* New `AINet` class wrapping core algorithm with composable module injection
* Both R6 composition and functional API equally supported

## New Modules
* `SHMEngine` — Five somatic hypermutation strategies: uniform (original
  behavior), airs (affinity-proportional), hotspot (feature-gradient-weighted),
  energy (budget-constrained), and adaptive (per-feature Adam-like moment
  tracking)
* `IdiotypicNetwork` — Antibody-antibody network dynamics with bell-shaped
  activation function replacing epsilon-threshold suppression
* `GerminalCenter` — T-follicular helper mediated selection with task-aware
  quality scoring and resource competition
* `Microenvironment` — Density-dependent zone classification
  (stable/explore/boundary) with chemokine-like gradient computation
* `VDJLibrary` — Combinatorial V(D)J gene library initialization via PCA,
  k-means clustering, or random partition of feature space
* `ActivationGate` — Two-signal activation gate requiring both antigen
  recognition and costimulatory context (density, danger, or entropy)
* `MemoryPool` — Archive high-affinity antibodies as long-lived memory cells
  with threshold-based recall
* `ClassSwitcher` — Isotype class switching (IgM broad, IgG specific, IgA
  boundary) modulating effective kernel width
* `ConvergentSelector` — Cross-repertoire consensus identification of public
  antibodies for ensemble methods

## Documentation
* Complete README rewrite covering functional API, R6 API, module reference
  table, and architecture overview
* New pkgdown website with Bootstrap 5 Flatly theme, organized reference groups,
  and tutorial navigation
* New article: "Composing Immune Modules" — R6 composition patterns for all 9
  modules with worked examples
* New article: "Advanced Tuning & Workflows" — swarmbHIVE grid search,
  honeycombHIVE multilayer refinement, refineB optimizer comparison, caret
  integration, and visualizeHIVE plot types
* New article: "Algorithm & Biological Foundations" — comprehensive mathematical
  reference covering all affinity kernels, distance functions, SHM strategies,
  idiotypic ODE system, germinal center selection, and parameter guidance
* Added roxygen @examples to refineB, bHIVEmodel, and ImmuneAlgorithm (now 90%
  example coverage)
* GitHub Actions workflow for automated pkgdown deployment to borch.dev

## Package Infrastructure
* Created `R/bHIVE-package.R` with roxygen-managed `@useDynLib` and
 `@importFrom Rcpp sourceCpp` directives
* Added `%||%` operator `@name null-coalesce` to avoid illegal characters in Rd
  `\name` field
* Moved tutorial vignettes to `vignettes/articles/` (pkgdown-only, not installed
  with package) to reduce installed size
* Added pkgdown configuration (`_pkgdown.yml`) with 7 reference groups and
  structured article hierarchy

## Testing
* Added comprehensive unit tests for all 12 R6 module classes (ImmuneRepertoire,
  ImmuneAlgorithm, AINet, SHMEngine, IdiotypicNetwork, GerminalCenter,
  Microenvironment, VDJLibrary, ActivationGate, MemoryPool, ClassSwitcher,
  ConvergentSelector) and C++ backend functions
* Test coverage increased from ~26% to ~85% (681 tests, 0 failures)

## BiocCheck Compliance
* Replaced all `sapply()` calls with `vapply()` in bHiVE.R and visualizeHIVE.R
* Replaced all `1:n` patterns with `seq_len()` / `seq_along()`
* Removed `install.packages()` calls from vignettes
* Removed `LazyData: true` from DESCRIPTION
* Updated R dependency to >= 4.5.0
* Updated biocViews to `Software, Clustering, Classification, Regression,
  Network`
* Added class-level `@param` documentation for all R6 initialize() arguments
  across 11 module classes
* Added `@param ... Not used.` to all R6 `print()` methods

## Bug Fixes
* Fixed kmeans++ sampling: corrected cumulative sum fallback index and
  runif-to-integer cast
* Fixed cosine similarity epsilon (1e-12) for numerical stability
* Fixed division-by-zero guard in idiotypic dynamics when theta_low equals
  theta_high
* Fixed VDJLibrary NaN in kmeans for single-allele edge cases
* Fixed VDJLibrary NA subscript when feature dimensions not divisible by 3
* Removed duplicate `Classification` from biocViews, added `Clustering`
* Fixed vignette YAML parsing error from bare `---` horizontal rule

# bHIVE 0.99.0

* Initial submission version with core AIS functionality
* Clonal selection, network suppression, and mutation for clustering,
  classification, and regression
* honeycombHIVE multilayer architecture
* swarmbHIVE hyperparameter tuning via BiocParallel
* caret model integration (bHIVEmodel, honeycombHIVEmodel)
* Visualization utilities via ggplot2
