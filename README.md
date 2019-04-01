
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mvtraits

The goal of `mvtraits` is to efficiently fit Bayesian (hierarchical)
multivariate models to trait data. A key feature is the ability to
perform multiple imputation of partially missing data.

## Installation

This package is not currently available on CRAN. To install directly
from this repository, do the following:

``` r
devtools::install_github("ashiklom/mvtraits")
```

## Example: The Iris dataset

This is a basic example showing how a multivariate distribution can be
fit to the classic Iris dataset.

``` r
library(mvtraits)
data(iris)
```

### Multivariate fit

The input data must be a matrix, with each trait in its own column.
Below, we drop the `Species` column and convert the result to a matrix.

``` r
iris_mat <- as.matrix(iris[, !grepl("Species", colnames(iris))])
```

Fit a multivariate distribution with `fit_mvnorm`:

``` r
iris_fit <- fit_mvnorm(iris_mat, parallel = FALSE)
#> Running sampler...
#> ===========================================================================
#> ===========================================================================
#> ===========================================================================
#> [1] "All parameters have converged"
#> Calculating correlation matrices...
#> Converting samples to coda mcmc.list object...
#> Preparing summary table...
#> Joining, by = "rowname"
str(iris_fit, max.level = 2)
#> List of 3
#>  $ summary_table:Classes 'tbl_df', 'tbl' and 'data.frame':   20 obs. of  12 variables:
#>   ..$ rowname       : chr [1:20] "mu..Sepal.Length" "mu..Sepal.Width" "mu..Petal.Length" "mu..Petal.Width" ...
#>   ..$ variable      : chr [1:20] "mu" "mu" "mu" "mu" ...
#>   ..$ index         : chr [1:20] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width" ...
#>   ..$ Mean          : num [1:20] 5.783 3.064 3.63 1.146 0.695 ...
#>   ..$ SD            : num [1:20] 0.0678 0.0359 0.1438 0.0624 0.0803 ...
#>   ..$ Naive SE      : num [1:20] 0.000782 0.000415 0.00166 0.000721 0.000928 ...
#>   ..$ Time-series SE: num [1:20] 0.000799 0.000409 0.001645 0.000713 0.000938 ...
#>   ..$ 2.5%          : num [1:20] 5.648 2.993 3.351 1.021 0.553 ...
#>   ..$ 25%           : num [1:20] 5.738 3.04 3.533 1.103 0.638 ...
#>   ..$ 50%           : num [1:20] 5.783 3.064 3.631 1.146 0.689 ...
#>   ..$ 75%           : num [1:20] 5.83 3.088 3.728 1.188 0.745 ...
#>   ..$ 97.5%         : num [1:20] 5.915 3.135 3.906 1.265 0.868 ...
#>  $ stats        :List of 3
#>   ..$ mu   :List of 3
#>   ..$ Sigma:List of 3
#>   ..$ Corr :List of 3
#>  $ samples      :List of 3
#>   ..$ : 'mcmc' num [1:5000, 1:20] 5.8 5.76 5.81 5.73 5.81 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>   .. ..- attr(*, "mcpar")= num [1:3] 1 5000 1
#>   ..$ : 'mcmc' num [1:5000, 1:20] 5.82 5.73 5.78 5.87 5.83 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>   .. ..- attr(*, "mcpar")= num [1:3] 1 5000 1
#>   ..$ : 'mcmc' num [1:5000, 1:20] 5.84 5.8 5.85 5.8 5.89 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>   .. ..- attr(*, "mcpar")= num [1:3] 1 5000 1
#>   ..- attr(*, "class")= chr "mcmc.list"
```

Compare the results to univariate summary statistics.

``` r
iris_fit[[c("stats", "mu", "Mean")]]
#> Sepal.Length  Sepal.Width Petal.Length  Petal.Width 
#>     5.782941     3.064301     3.629673     1.145615
colMeans(iris_mat)
#> Sepal.Length  Sepal.Width Petal.Length  Petal.Width 
#>     5.843333     3.057333     3.758000     1.199333
iris_fit[[c("stats", "Sigma", "Mean")]]
#>              Sepal.Length Sepal.Width Petal.Length Petal.Width
#> Sepal.Length   0.69486901 -0.04274195    1.2799306   0.5188079
#> Sepal.Width   -0.04274195  0.19647519   -0.3298816  -0.1217633
#> Petal.Length   1.27993063 -0.32988162    3.1353441   1.3013037
#> Petal.Width    0.51880794 -0.12176327    1.3013037   0.5903194
cov(iris_mat)
#>              Sepal.Length Sepal.Width Petal.Length Petal.Width
#> Sepal.Length    0.6856935  -0.0424340    1.2743154   0.5162707
#> Sepal.Width    -0.0424340   0.1899794   -0.3296564  -0.1216394
#> Petal.Length    1.2743154  -0.3296564    3.1162779   1.2956094
#> Petal.Width     0.5162707  -0.1216394    1.2956094   0.5810063
```

### Partial missing data

Fitting models also works with partially missing data. Here, we remove
half of the data at random.

``` r
iris_mat2 <- iris_mat
nmat <- length(iris_mat2)
iris_mat2[sample.int(nmat, floor(nmat / 2))] <- NA
head(iris_mat2)
#>      Sepal.Length Sepal.Width Petal.Length Petal.Width
#> [1,]          5.1          NA          1.4          NA
#> [2,]          4.9         3.0           NA          NA
#> [3,]          4.7          NA          1.3         0.2
#> [4,]           NA         3.1          1.5          NA
#> [5,]           NA         3.6          1.4          NA
#> [6,]          5.4         3.9          1.7         0.4
```

Now, re-fit the model.

``` r
iris_fit2 <- fit_mvnorm(iris_mat2, parallel = FALSE)
#> Running sampler...
#> ===========================================================================
#> ===========================================================================
#> ===========================================================================
#> [1] "All parameters have converged"
#> Calculating correlation matrices...
#> Converting samples to coda mcmc.list object...
#> Preparing summary table...
#> Joining, by = "rowname"
iris_fit2[[c("stats", "mu", "Mean")]]
#> Sepal.Length  Sepal.Width Petal.Length  Petal.Width 
#>     5.810113     3.044235     3.676280     1.151628
iris_fit2[[c("stats", "Sigma", "Mean")]]
#>              Sepal.Length Sepal.Width Petal.Length Petal.Width
#> Sepal.Length    0.7070620  -0.1265480    1.3028964   0.5152218
#> Sepal.Width    -0.1265480   0.2334190   -0.3971634  -0.1350961
#> Petal.Length    1.3028964  -0.3971634    3.0563905   1.2090920
#> Petal.Width     0.5152218  -0.1350961    1.2090920   0.5576549
```

### Hierarchical model

`mvtraits` can also fit hierarchical models with the `fit_mvnorm_hier`
function. The arguments are the same as `fit_mvnorm`, but also need a
`group` argument that takes a vector of groups (similar to base R’s
`tapply`).

``` r
iris_fit2h <- fit_mvnorm_hier(iris_mat2, iris[["Species"]], parallel = FALSE)
#> Running sampler...
#> =================================================================================================================================================================================================================================[1] "All parameters have converged"
#> Calculating correlation matrices...
#> Converting samples to coda mcmc.list object...
#> Preparing summary table...
#> Joining, by = "rowname"
str(iris_fit2h, max.level = 2)
#> List of 3
#>  $ summary_table:Classes 'tbl_df', 'tbl' and 'data.frame':   80 obs. of  13 variables:
#>   ..$ rowname       : chr [1:80] "mu..global..Sepal.Length" "mu..global..Sepal.Width" "mu..global..Petal.Length" "mu..global..Petal.Width" ...
#>   ..$ variable      : chr [1:80] "mu" "mu" "mu" "mu" ...
#>   ..$ group         : chr [1:80] "global" "global" "global" "global" ...
#>   ..$ index         : chr [1:80] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width" ...
#>   ..$ Mean          : num [1:80] 1.2991 0.968 -0.0193 -0.0732 21.807 ...
#>   ..$ SD            : num [1:80] 1.078 0.997 0.955 0.868 28.667 ...
#>   ..$ Naive SE      : num [1:80] 0.0124 0.0115 0.011 0.01 0.331 ...
#>   ..$ Time-series SE: num [1:80] 0.0347 0.0356 0.0236 0.0358 0.4216 ...
#>   ..$ 2.5%          : num [1:80] -0.773 -1.037 -1.848 -1.834 2.952 ...
#>   ..$ 25%           : num [1:80] 0.57 0.298 -0.66 -0.626 8.399 ...
#>   ..$ 50%           : num [1:80] 1.3023 0.9872 -0.0296 -0.0476 14.381 ...
#>   ..$ 75%           : num [1:80] 2.028 1.666 0.608 0.487 24.957 ...
#>   ..$ 97.5%         : num [1:80] 3.41 2.84 1.89 1.64 87.15 ...
#>  $ stats        :List of 6
#>   ..$ mu_global   :List of 3
#>   ..$ Sigma_global:List of 3
#>   ..$ Corr_global :List of 3
#>   ..$ mu_group    :List of 3
#>   ..$ Sigma_group :List of 3
#>   ..$ Corr_group  :List of 3
#>  $ samples      :List of 3
#>   ..$ : 'mcmc' num [1:5000, 1:80] 1.94 2.74 2.95 2.68 2.3 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>   .. ..- attr(*, "mcpar")= num [1:3] 1 5000 1
#>   ..$ : 'mcmc' num [1:5000, 1:80] 2.77 2.53 3.31 4.14 4.27 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>   .. ..- attr(*, "mcpar")= num [1:3] 1 5000 1
#>   ..$ : 'mcmc' num [1:5000, 1:80] 1.55 3.22 3.82 4.5 3.69 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>   .. ..- attr(*, "mcpar")= num [1:3] 1 5000 1
#>   ..- attr(*, "class")= chr "mcmc.list"
```
