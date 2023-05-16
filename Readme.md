rocmice
================

# Introduction

`rocmice` is an R package to work with ROC (Receiver Operating Curves)
in multiply imputed data. Using multiple imputation in prediction
modelling can be helpful when missingness patterns don’t follow MCAR and
complete case analysis, single imputation or similar methods could
seriously bias the model. Using pooled ROC and AUC makes evaluation of
such models easier.

# Installation

You can install the package with

``` r
# install.packages("devtools")
devtools::install_github("https://github.com/jonas-schropp/rocmice.git")
```

And then call it like any other package using

``` r
library(rocmice)
```

# Usage

## Use case

To illustrate the use of `rocmice` we use an especially pathological
simulated data set. To learn more how it was created and what makes ROC
analysis on this data set tricky, call

``` r
data(patho)
help("patho")
```

    ## No documentation for 'patho' in specified packages and libraries:
    ## you could try '??patho'

## Plotting the ROC curve

First we pool the ROC for the imputed data sets using `pool_roc_rr`:

``` r
pooled_roc <- pool_roc_rr(
    data = patho$patho_imp,
    score = "score",
    target = "outcome",
    unique_vals = NULL,
    fpr_vals = seq(from=0.001, to=0.999, by=0.001),
    backtransform = TRUE,
    ci.level = 0.95,
    corr = 0.5,
    verbose = TRUE
)
```

    ## 
    ## Calculating TPR for every FPR: 
    ## ================================================================================
    ## Combining and pooling results.

Then we plot the ROC for the complete data set and the complete cases
from the data set with missing values using `plotROC::geom_roc` and add
the pooled ROC created with `pool_roc_rr` using `ggplot::geom_step`:

``` r
ggplot() + 
  geom_roc(
    data = patho$patho_complete, 
    aes(m = score, d = outcome, color = 'complete data')
    ) +
  geom_roc(
    data = patho$patho_amp, 
    aes(m = score, d = outcome, color = 'complete cases')
    ) +
  geom_step(
    data = pooled_roc, 
    aes(fpr, roc, color = 'multiple imputation'), 
    size = 1
    ) +
  labs(x = "FPR", y = "TPR") +
  theme_minimal() +
  scale_color_manual(
    breaks = c(
      'complete data', 
      'complete cases', 
      'multiple imputation'
      ),
    values = c(
      'complete data'='blue', 
      'complete cases'='red', 
      'multiple imputation'='black'
      )
    ) +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom"
    )
```

![](Readme_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

It is easy to see that the ROC created from multiply imputed data
closely tracks the one on a complete data set without missing values,
while the one on complete cases shows an overly optimistic bias.
