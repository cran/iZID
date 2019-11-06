
# iZID

<!-- badges: start -->
<!-- badges: end -->

iZID computes bootstrapped Monte Carlo estimate of p-value of KS test and likelihood ratio test for zero-inflated count data based on the previous work of Aldirawi et al. (2019). This package also enables user to compute maximum likelihood estimate of data from standard, zero-inflated or hurdle beta binomial, beta negative binomial, negative binomial and Poisson distributions. Besides, user can generate random deviates from the aforementioned distributions.


## Installation

The released version of iZID can be downloaded from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("iZID")
```

## Architecture 

12 functions are exported from this package which can be classified as five classes:

* bb.mle, bnb.mle, nb.mle, poisson.mle: calculate maximum likelihood estimates for general distributions.
* bb.zihmle, bnb.zihmle, nb.zihmle, poisson.zihmle: calculate maximum likelihood estimates for zero-inflated or hurdle distributions.
* dis.kstest: conduct one-sample KS test and output bootstrapped p-value.
* model.lrt: conduct likelihood ratio test to compare two models and output bootstrapped p-value.
* sample.h, sample.zi: simulate random deviates from zero-inflated or hurdle models.

