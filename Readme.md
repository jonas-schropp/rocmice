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

Notable, it contains a binary variable `outcome` that represents some
kind of status, for example presence of a disease and a `score` that
represents the result of some kind of predictive / diagnostic algorithm
(it was calculated as `invlogit(0.5 * X1 + X2 + 1.5 * X3)`). This score
performs worse when `X4` is 1 then when it is 0.

## Pooling the AUC

We can calculate the AUC for each imputation using `pROC::roc` on each
list element, or, to make it easier, using the provided convenience
function `apply_roc`.

``` r
# So this:
rocs1 <- list()
for (i in 1:length(patho$patho_imp)) {
  rocs1[[i]] <- pROC::roc(
    outcome ~ score, 
    data = patho$patho_imp[[i]],
    direction = ">", levels = c(0, 1)
    )
}

# is pretty much equivalent to this:
rocs2 <- apply_roc(
  milist = patho$patho_imp,
  selection = "all",
  score = "score", 
  target = "outcome"
)
```

And then pool it with `pool_auc_rr`:

``` r
pool_auc_rr(rocs2, ci.level = 0.95, transform = "logit")
```

    ##         auc        ll        ul auc_logit var_total_logit var_between_logit
    ## 1 0.4749906 0.3933517 0.5579887 -0.100121      0.02749194       0.006515622
    ##   var_within_logit
    ## 1       0.01967319

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

    ## Calculation of confidence intervals is experimental. 
    ##           The lower and upper bounds that are currently supplied are simply 
    ##           the lower and upper bounds of the pooled TPR at a given alpha level.

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

![](Readme_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

It is easy to see that the ROC created from multiply imputed data
closely tracks the one on a complete data set without missing values,
while the one on complete cases shows an overly optimistic bias.

You can in theory create confidence intervals for the pooled ROC curve,
but this functionality is experimental and currently simply creates
pooled confidence intervals for each TPR value:

``` r
ggplot(pooled_roc, 
       aes(fpr, roc, ymin = ll_roc, ymax = ul_roc)
       ) + 
  geom_ribbon(alpha = 0.2) +
  geom_step(size = 1) +
  labs(x = "FPR", y = "TPR") +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom"
  ) +
  theme_minimal()
```

![](Readme_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

## Comparing the AUC between two groups

To compare the area under the curve, we rely on a modified version of
Delong’s test that pools the variance of AUC1 - AUC2 within and between
each imputation using Rubin’s Rules.

``` r
mi_roc_test(
  data = patho$patho_imp,
  target = "outcome",
  score = "score",
  group = "X4",
  groups = c(0, 1),
  levels = c(0, 1),
  direction = "<",
  paired = FALSE
) |> 
  round(3)
```

    ##   delta_auc auc1  auc2 t.value p.value var.total var.within var.between  riv
    ## 1     0.438  0.8 0.362   4.663       0     0.094       0.09       0.003 0.04
    ##   lambda   fmi
    ## 1  0.038 0.049

Alternatively, you could rely on `pROC::roc` for the calculation of the
ROC for each group and imputation. This results in more code and
slightly longer computation time, but might be helpful if you need to
access the objects for each imputation created by `pROC` directly.

``` r
g1 <- list()
g2 <- list()
for (i in 1:length(patho$patho_imp)) {
  
  g1[[i]] <- pROC::roc(
    outcome ~ score, 
    data = patho$patho_imp[[i]][patho$patho_imp[[i]]$X4 == 0, ],
    direction = "<", levels = c(0, 1)
    )
  g2[[i]] <- pROC::roc(
    outcome ~ score, 
    data = patho$patho_imp[[i]][patho$patho_imp[[i]]$X4 == 1, ],
    direction = "<", levels = c(0, 1)
    )
}
```

``` r
mi_roc_test(rocs1 = g1, rocs2 = g2, paired = FALSE) |> round(3)
```

    ##   delta_auc auc1  auc2 t.value p.value var.total var.within var.between  riv
    ## 1     0.438  0.8 0.362   4.663       0     0.094       0.09       0.003 0.04
    ##   lambda   fmi
    ## 1  0.038 0.049

## Comparing the AUC for two scores

We could now rightfully decide that our score needs to be updated to
account for the unfair bias against group 1. It’s obvious that the score
performs much better in group 0 than group 1.

``` r
pl1 <- ggplot(patho$patho_imp[[1]], aes(score, outcome, group = factor(X4), color = factor(X4))) + 
  geom_smooth(method = 'loess', formula = 'y ~ x') +
  theme(legend.position = "bottom")

pl2 <- ggplot(patho$patho_imp[[1]], aes(X4, fill = factor(outcome))) + 
  geom_bar() +
  theme(legend.position = "bottom")

ggpubr::ggarrange(pl1, pl2)
```

![](Readme_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

It does appear like group 1 is also several times less likely to suffer
from the outcome condition, so we could incorporate that information
into a new score.

``` r
calc_score2 <- function(d) {
  eta <- 0.5 * d[["X1"]] + d[["X2"]] + 1.5 * d[["X3"]] + 4 * d[["X4"]]
  rocmice:::invlogit(eta)
}

for (i in 1:length(patho$patho_imp)) {
  patho$patho_imp[[i]]$score2 <- calc_score2(patho$patho_imp[[i]])
}
```

If your imputed data sets are in a long format, you can use the same
syntax as above. Alternatively, you can specify two score variables like
below if the data sets are in wide format (like ours here):

``` r
mi_roc_test(
  data = patho$patho_imp,
  target = "outcome",
  score = "score",
  score2 = "score2",
  levels = c(0, 1),
  direction = "<",
  paired = TRUE
) |> 
  round(3)
```

    ##   delta_auc  auc1  auc2     Z p.value var.total var.within var.between   riv
    ## 1      0.12 0.525 0.405 5.517       0         0          0           0 0.138
    ##   lambda   fmi
    ## 1  0.121 0.128
