# Bigcluster

This is a slight modification of the `vcovCL` function for cluster-robust standard errors found in the R package [sandwich](https://cran.r-project.org/web/packages/sandwich/index.html) (used version 2.4-0) so that it can be used with `biglm` (and `bigglm`) objects. It is for the admittedly limited use case where data is small enough to fit in memory but too big for models using the built-in `lm` or `glm` procedures with clustered standard errors. So far only the HC0 and HC1 methods are supported (for more information on these different methods, see [MacKinnon and White (1985)](http://qed.econ.queensu.ca/working_papers/papers/qed_wp_537.pdf) and [Zeileis (2004)](https://www.jstatsoft.org/article/view/v011i10)). Only ordinary least-squares (biglm) and Poisson and negative binomial (bigglm) regression are supported for now.

## Getting Started

The function was tested using R version 3.4.2. 

### Prerequisites

The following packages are required for functionality (version numberes in parentheses):
* biglm (0.9-1)
* bigmemory (4.5.19)
* bigtabulate (1.1.5)
* bootSVD (0.5)
* data.table (1.10.4)
* ff (2.2-13)
* magrittr (1.5)
* MASS (7.3-47)
* Matrix (1.2-11)

### Installation

Cloning the repo and installing the above prerequisites into R should be enough to get you up and running.

## Usage

First, load the `BigCluster` function and its helper functions into R via `source`
```
source('bigcluster_sandwich.R')
```

To obtain the cluster-robust standard errors of a `biglm` object

```
bigvcov <- BigCluster(YOUR_BIGLM_OBJECT, YOUR_CLUSTERS, YOUR_DATA, sumMethod = "slow")
```

where `YOUR_BIGLM_OBJECT` is the `biglm` model you have already run/trained; `YOUR_CLUSTERS` have the cluster assignments for all observations; and `YOUR_DATA` contains the response and explanatory variable valus for all observations (in the same order as the `YOUR_CLUSTERS` object). 

The `sumMethod` option can take on the value `"slow"` or `"fast"`. The "slow" method uses less memory but takes much longer, while the "fast" method is much faster but takes up more memory. The empirical estimating functions in `BigCluster` are stored as a sparse matrix. The "fast" method converts the entire sparse matrix into a dense one and applies the `sum` function. The "slow" method does not convert the entire matrix into a dense one. If your data is large enough that you need to use this function and `biglm`, then you are most likely stuck with the "slow" option.

You can also specify a `type` option, which can only take on the values `"HC0"` or `"HC1"`. If the type is not specified, then it will default to HC0-type errors for OLS models and HC1 for everything else. This is the same as the default behavior in `vcovCL` from the sandwich package. 

If you want to conduct hypothesis testing with this cluster-robust covariance matrix, you can use the `coeftest` function from the lmtest package

```
library(lmtest)

coeftest(YOUR_BIGLM_OBJECT, bigvcov)
```

## Tests

Tests comparing the output of `BigCluster` with sandwich's `vcovCL` are included in the test\_bigcluster.R script, which relies on helper functions in prof\_func.R to generate negative-binomial distributed data. The possible options are
* `nRows`: how many observations the simulated data should have
* `nRegressors`: how many non-cluster explanatory variables there should be
* `nClusterVars`: how many cluster variables there should be
* `nClusterCats`: how many possible cluster values each cluster variable can take
* `nSims`: how many instances of simulated data should be generated
* `modelFamily`: whether to run an OLS "ols", Poisson "poisson", or negative binomial "negbin" model on the simulated data

The for loop will generate the data, run `lm` and `biglm` models on the data, and then save the differences between the covariance matrix generated via `vcovCL` and the one generated via `BigCluster`. Values will be printed to demonstrate the differences between each matrix. Values include a random entry in the matrix; the entry where the output differs the most; the mean percentage difference between the two matrices; and the largest magnitude difference between entries. In my own tests, the largest magnitude difference is generally on the order of 10^-5. A plot can be shown that indicates that the highest percentage error occurs in entries that have the lowest magnitudes.

NOTE: Differences are low in magnitude when comparing `vcovCL` and `BigCluster` when the number of cluster categories is less than 10. As of sandwich version 2.4-0, when the number of clusteres is greater than 10, then errors tend to increase. This is due to differences between how `vcovCL` and `BigCluster` identify the intersection of clusters when there are more than one cluster. The `vcovCL` way will conflate some cluster combinations as the same even when they are different.

