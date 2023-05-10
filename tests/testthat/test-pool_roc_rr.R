r <- c(1, 0, 1, 0, 1, 0, 0, 1)
p <- c(1, 0, 1, 1, 0, 1, 0, 1)

test_that("tpr works", {
  expect_equal(tpr(r, p), 3/4)
})

test_that("fpr works", {
  expect_equal(fpr(r, p), 2/4)
})


test_that("var_prop works", {
  est <- c(0.2, 0.4, 0.6, 0.8)
  n <- 100
  expect_equal(var_prop(est, n), c(0.0016, 0.0024, 0.0024, 0.0016))
})


test_that("cut_off works", {
  score <- c(1.2, 0.5, 0.8, 1.5, 2.1)
  cutoff <- 1
  expect_equal(cut_off(score, cutoff), c(1,0,0,1,1))
})


m <- matrix(
  c(rep(1,8),0,1,rep(0,5)),
  nrow = 5, ncol = 3, byrow=T)

test_that("apply_cut_off works", {
  score <- c(1.2, 0.8, 0.5, 0.3, 0.1)
  unique_vals <- c(0.2, 0.4, 0.6, 0.8)
  expect_equal(
    apply_cut_off(score, unique_vals), 
    m)
})


r <- c(1, 0, 1, 0, 1)
apply_metrics(m, r)
l <- list(
  tpri = c(2/3, 2/3, 1/3),
  fpri = c(1, 1/2, 1/2),
  vartpri = c(0.074, 0.074, 0.074)
)
test_that("apply_metrics works", {
  expect_equal(apply_metrics(m, r), l, tolerance = 3)
})



test_that("rcpp_rowMeans works", {
  expect_equal(rcpp_rowMeans(m), rowMeans(m))
})


test_that("rcpp_var works", {
  expect_equal(rcpp_var(1:10), var(1:10))
})


est = matrix(c(2:5), 2,2)
var = matrix(c(0.2, 0.3, 0.4, 0.2), 2, 2)
test_that("rcpp_var works", {
  pool_rr_m_R <- function (est, var) {
    
    m <- ncol(est)
    mean_est <- rowMeans(est)
    var_w <- rowMeans(var)
    var_b <- apply(est, 1, var) 
    var_b[is.nan(var_b)] <- 0
    var_T <- var_w + (1 + 1 / m) * var_b
    
    return(list(
      roc = mean_est, 
      var_roc = var_T
    ))
  }
  expect_equal(
    pool_rr_m(est, var), pool_rr_m_R(est, var)
  )
})


