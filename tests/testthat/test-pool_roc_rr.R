r <- c(1, 0, 1, 0, 1, 0, 0, 1)
p <- c(1, 0, 1, 1, 0, 1, 0, 1)

test_that("tpr works", {
  expect_equal(tpr(r, p, 0)[1], logittrans(3/4))
})

test_that("fpr works", {
  expect_equal(fpr(r, p, 0)[1], logittrans(2/4))
})


m <- matrix(
  c(rep(1,8),0,1,rep(0,5)),
  nrow = 5, ncol = 3, byrow=T)



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


test_that("filterArray works", {
  x <- matrix(
    data = c(
      0, 1,
      0, 2,
      0, 1,
      1, 2,
      1, 1,
      3, 1,
      4, 0,
      5, 5,
      5, 2
    ),
    nrow=9, ncol=2, byrow = TRUE
  )
  expect_equal(filterArray(x), c(1, 0, 1, 0, 1, 0, 0, 0, 1))
})
