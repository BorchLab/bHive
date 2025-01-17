# bHIVE

## Artifical Immune Network-based Machine Learning Package

<img align="right" src="https://github.com/ncborcherding/bHive/blob/main/www/bhive_hex.png" width="305" height="352">

### Introduction

**bHIVE** is an R package implementing an Artificial Immune Network (AI-Net) algorithm. 
Inspired by principles of immunology, the bHIVE algorithm evolves a population of 
"antibodies" to model and analyze datasets. It supports three tasks:

- **Clustering**: Groups data points based on similarity.
- **Classification**: Assigns labels to data points based on learned patterns.
- **Regression**: Predicts continuous outcomes based on input features.

The AI-Net algorithm uses clonal selection and mutation, network suppression, 
and affinity/distance metrics to iteratively refine the antibody population.

### Algorithm Basis

The algorithm operates in the following steps:

1. **Affinity Calculation**: For a data point $( x )$ and an antibody $( a )$, the affinity is defined as:
   $$
   \text{Affinity}(x, a) = f(x, a; \text{params})
   $$
   where \( f \) is a similarity kernel such as Gaussian (RBF), Laplace, 
   or cosine similarity.
   
2. **Clonal Expansion and Mutation**: Top-matching antibodies are cloned and mutated:
   $$
   a' = a + \epsilon
   $$
   where $\epsilon$ is sampled from a distribution scaled by affinity and 
   iteration parameters.

3. **Network Suppression**: Antibodies that are too similar (within a threshold $\epsilon$) are suppressed to maintain diversity:
   $$
   \text{Distance}(a_i, a_j) < \epsilon \implies \text{Suppress } a_j
   $$

4. **Assignment**: Data points are assigned to antibodies using affinity (for classification/regression) or distance (for clustering).

## Installation

To install the **bHIVE** from GitHub, use:

```R
# Install bHIVE from GitHub
devtools::install_github("yourusername/bHIVE")
```

## Examples

### Clustering

```
# Load the Iris dataset
data(iris)
X <- as.matrix(iris[, 1:4])  # Use numeric features only

# Run bHIVE for clustering
res <- bHIVE(X = X, 
             task = "clustering", 
             nAntibodies = 30, 
             beta = 5, 
             epsilon = 0.01, 
             maxIter = 20, 
             k = 3)

# View clustering results
table(res$assignments)
```

### Classification

```
# Load the Iris dataset
data(iris)
X <- as.matrix(iris[, 1:4])  # Use numeric features
y <- iris$Species           # Classification labels

# Run bHIVE for classification
res <- bHIVE(X = X, 
             y = y, 
             task = "classification", 
             nAntibodies = 30, beta = 5, 
             epsilon = 0.01, 
             maxIter = 20, 
             k = 3)

# Compare predicted labels with actual labels
table(Predicted = res$assignments, Actual = y)
```

### Regression

```
# Load the Iris dataset
data(iris)
X <- as.matrix(iris[, 2:4])  # Use other features as predictors
y <- iris$Sepal.Length       # Regression target

# Run bHIVE for regression
res <- bHIVE(X = X, 
             y = y, 
             task = "regression",
             nAntibodies = 30, 
             beta = 5, 
             epsilon = 0.01, 
             maxIter = 20, 
             k = 3)

# Compare predicted values with actual values
cor(res$assignments, y)
```

## Bug Reports/New Features

#### If you run into any issues or bugs please submit a [GitHub issue](https://github.com/ncborcherding/bHIVE/issues) with details of the issue.

- If possible please include a [reproducible example](https://reprex.tidyverse.org/). 

##### Any requests for new features or enhancements can also be submitted as [GitHub issues](https://github.com/ncborcherding/bHIVE/issues).

### Contributing

We welcome contributions to the bHIVE project! To contribute:

* Fork the repository.
* Create a feature branch (git checkout -b feature-branch).
* Commit your changes (git commit -m "Add new feature").
* Push to the branch (git push origin feature-branch).
* Open a pull request.