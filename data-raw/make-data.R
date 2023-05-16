library(MASS)
library(mice)


set.seed(101) 

# SIMULATION -------------------------------------------------------------------

n <- 300
mu <- rep(0, 10)
Sigma <- matrix(0.6, nrow=10, ncol=10)
diag(Sigma) <- 1
df_complete <- data.frame(mvrnorm(n = n, mu = mu, Sigma = Sigma))

df_complete[[4]] <- as.integer(df_complete[[4]] > 0)

# Simulate the outcome variable
eta <- 0.6 * df_complete[[1]] + 
  df_complete[[2]] + 
  1.3 * df_complete[[3]] + 
  -3 * df_complete[[4]] + 
  -2 * df_complete[[3]] * df_complete[[4]] +
  -2 * df_complete[[2]] * df_complete[[4]] +
  rnorm(300, 0, 0.5)

outcome <- rbinom(300, 1, rocmice:::invlogit(eta))


# AMPUTATION -------------------------------------------------------------------

patterns <- matrix(1, 4, 11)
colnames(patterns) <- c(paste0("X", 1:10), "outcome")
patterns[4,1:3] <- 0
diag(patterns[1:3, 1:3]) <- 0
weights <- patterns
weights[,"outcome"] <- 5
weights[4,4] <- -10


df_complete$outcome <- outcome

df_miss <- ampute(
  data = df_complete,
  patterns = patterns,
  prop = 0.8,
  mech = "MAR"
)$amp

# IMPUTATION -------------------------------------------------------------------

imp <- mice(df_miss, 5)
df_imp <- complete(imp, "all")

# ADD THE SCORE ----------------------------------------------------------------

calc_score <- function(x1, x2, x3) {
  eta <- 0.5 * x1 + x2 + 1.5 * x3
  rocmice:::invlogit(eta)
}

df_complete$score <- calc_score(
  df_complete[[1]], df_complete[[2]], 
  df_complete[[3]]
)
df_miss$score <- calc_score(
  df_miss[[1]], df_miss[[2]], 
  df_miss[[3]]
)
for (i in 1:length(df_imp)) {
  df_imp[[i]]$score <- calc_score(
    df_imp[[i]][[1]], df_imp[[i]][[2]], 
    df_imp[[i]][[3]]
  )
}

# PUT TOGETHER AND SAVE --------------------------------------------------------

patho <- list(df_complete, df_miss, df_imp)

names(patho) <- c("patho_complete", "patho_amp", "patho_imp")

usethis::use_data(
  patho,
  overwrite = TRUE
)

