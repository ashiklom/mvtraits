
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
iris_fit <- fit_mvnorm(iris_mat)
#> Running sampler...
#> All parameters have converged
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
#>   ..$ Mean          : num [1:20] 5.782 3.063 3.63 1.146 0.695 ...
#>   ..$ SD            : num [1:20] 0.0667 0.036 0.1456 0.0634 0.0817 ...
#>   ..$ Naive SE      : num [1:20] 0.000771 0.000415 0.001681 0.000732 0.000943 ...
#>   ..$ Time-series SE: num [1:20] 0.000771 0.000415 0.001747 0.000732 0.000943 ...
#>   ..$ 2.5%          : num [1:20] 5.648 2.994 3.34 1.019 0.554 ...
#>   ..$ 25%           : num [1:20] 5.737 3.038 3.531 1.104 0.637 ...
#>   ..$ 50%           : num [1:20] 5.783 3.063 3.633 1.148 0.688 ...
#>   ..$ 75%           : num [1:20] 5.829 3.087 3.731 1.19 0.746 ...
#>   ..$ 97.5%         : num [1:20] 5.91 3.135 3.908 1.264 0.874 ...
#>  $ stats        :List of 3
#>   ..$ mu   :List of 3
#>   ..$ Sigma:List of 3
#>   ..$ Corr :List of 3
#>  $ samples      :List of 3
#>   ..$ : 'mcmc' num [1:5000, 1:20] 5.82 5.74 5.79 5.82 5.74 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>   .. ..- attr(*, "mcpar")= num [1:3] 1 5000 1
#>   ..$ : 'mcmc' num [1:5000, 1:20] 5.85 5.74 5.78 5.82 5.74 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>   .. ..- attr(*, "mcpar")= num [1:3] 1 5000 1
#>   ..$ : 'mcmc' num [1:5000, 1:20] 5.83 5.75 5.78 5.82 5.73 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>   .. ..- attr(*, "mcpar")= num [1:3] 1 5000 1
#>   ..- attr(*, "class")= chr "mcmc.list"
```

Compare the results to univariate summary statistics.

``` r
iris_fit[[c("stats", "mu", "Mean")]]
#> Sepal.Length  Sepal.Width Petal.Length  Petal.Width 
#>     5.782125     3.062608     3.630155     1.146483
colMeans(iris_mat)
#> Sepal.Length  Sepal.Width Petal.Length  Petal.Width 
#>     5.843333     3.057333     3.758000     1.199333
iris_fit[[c("stats", "Sigma", "Mean")]]
#>              Sepal.Length Sepal.Width Petal.Length Petal.Width
#> Sepal.Length   0.69497210 -0.04291286    1.2810738   0.5189763
#> Sepal.Width   -0.04291286  0.19642316   -0.3308101  -0.1223293
#> Petal.Length   1.28107379 -0.33081015    3.1404326   1.3031064
#> Petal.Width    0.51897633 -0.12232934    1.3031064   0.5909668
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
#> [1,]          5.1         3.5          1.4         0.2
#> [2,]          4.9         3.0           NA         0.2
#> [3,]          4.7         3.2           NA         0.2
#> [4,]          4.6         3.1          1.5         0.2
#> [5,]          5.0          NA          1.4         0.2
#> [6,]           NA          NA          1.7         0.4
```

Now, re-fit the model.

``` r
iris_fit2 <- fit_mvnorm(iris_mat2)
#> Running sampler...
#> All parameters have converged
#> Calculating correlation matrices...
#> Converting samples to coda mcmc.list object...
#> Preparing summary table...
#> Joining, by = "rowname"
iris_fit2[[c("stats", "mu", "Mean")]]
#> Sepal.Length  Sepal.Width Petal.Length  Petal.Width 
#>     5.790697     3.047071     3.632278     1.146587
iris_fit2[[c("stats", "Sigma", "Mean")]]
#>              Sepal.Length Sepal.Width Petal.Length Petal.Width
#> Sepal.Length   0.72814113 -0.07677293    1.3472674   0.5443821
#> Sepal.Width   -0.07677293  0.17147239   -0.2328279  -0.1018946
#> Petal.Length   1.34726742 -0.23282791    3.1963764   1.3270015
#> Petal.Width    0.54438209 -0.10189459    1.3270015   0.6116293
```

### Hierarchical model

`mvtraits` can also fit hierarchical models with the `fit_mvnorm_hier`
function. The arguments are the same as `fit_mvnorm`, but also need a
`group` argument that takes a vector of groups (similar to base R’s
`tapply`).

``` r
iris_fit2h <- fit_mvnorm_hier(iris_mat2, iris[["Species"]])
#> Running sampler...
#> All parameters have converged
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
#>   ..$ Mean          : num [1:80] 1.222 0.944 -0.068 -0.173 22.097 ...
#>   ..$ SD            : num [1:80] 1.029 0.926 0.962 0.813 27.465 ...
#>   ..$ Naive SE      : num [1:80] 0.01188 0.01069 0.0111 0.00939 0.31708 ...
#>   ..$ Time-series SE: num [1:80] 0.0317 0.0347 0.0236 0.0326 0.4221 ...
#>   ..$ 2.5%          : num [1:80] -0.867 -0.957 -1.983 -1.77 3.22 ...
#>   ..$ 25%           : num [1:80] 0.545 0.323 -0.694 -0.714 8.589 ...
#>   ..$ 50%           : num [1:80] 1.234 0.98 -0.069 -0.164 14.442 ...
#>   ..$ 75%           : num [1:80] 1.923 1.594 0.58 0.356 25.147 ...
#>   ..$ 97.5%         : num [1:80] 3.24 2.66 1.82 1.43 88.81 ...
#>  $ stats        :List of 6
#>   ..$ mu_global   :List of 3
#>   ..$ Sigma_global:List of 3
#>   ..$ Corr_global :List of 3
#>   ..$ mu_group    :List of 3
#>   ..$ Sigma_group :List of 3
#>   ..$ Corr_group  :List of 3
#>  $ samples      :List of 3
#>   ..$ : 'mcmc' num [1:5000, 1:80] 2.01 1.67 2.01 3.06 2.51 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>   .. ..- attr(*, "mcpar")= num [1:3] 1 5000 1
#>   ..$ : 'mcmc' num [1:5000, 1:80] 2.93 3.86 3.47 3.89 3.29 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>   .. ..- attr(*, "mcpar")= num [1:3] 1 5000 1
#>   ..$ : 'mcmc' num [1:5000, 1:80] 2.42 3.57 4.16 4.31 4.27 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>   .. ..- attr(*, "mcpar")= num [1:3] 1 5000 1
#>   ..- attr(*, "class")= chr "mcmc.list"
```
