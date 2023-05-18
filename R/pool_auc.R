#' Pool AUC after Multiple Imputation
#' 
#' Pools AUC, it's variance and confidence intervals from a list of
#' roc objects using the methods described by Delong. 
#' This function assumes that we did NOT impute the dependent variable.
#' 
#' @describeIn pool_auc_rr
#' 
#' @param rocs A list of `pROC::roc` objects for each imputed data set
#' @param ci.level CI level for the AUC. Double between 0-1, by default 0.95.
#' @param transform A character string with the transformation if the AUC should 
#' be transformed before pooling or FALSE if the raw AUC should be pooled. 
#' "logit" by default.
#' 
#' @author 
#' Jonas Schropp, code fragments from `pROC::roc.test` and `psfmi::pool_auc`.
#' 
#' @export
#' 
pool_auc_rr <- function(
    rocs, 
    ci.level = 0.95, 
    transform = "logit"
    ) {
  
  # Number of mi data sets
  m <- length(rocs)
  # Variances for each AUC
  S <- double(length = m)
  # Each AUC
  aucs <- double(length = m)
  
  for (i in 1:m) {
    
    nh <- length(rocs[[i]]$controls)
    nc <- length(rocs[[i]]$cases)
    
    V <- delongPlacements(rocs[[i]])
    
    aucs[i] <- V$theta 
    
    SX <- sum((V$X - V$theta) * (V$X - V$theta))/(nc-1)
    SY <- sum((V$Y - V$theta) * (V$Y - V$theta))/(nh-1)
    
    S[i] <- SX/nc + SY/nh   # Variance R identical with var(roc1)
    
  }
  
  SE <- sqrt(S) # Standard Error for each AUC
  
  pool_auc_se(
    AUC = aucs, 
    SE = SE, 
    transform = transform, 
    ci.level = ci.level,
    m = m
    )
  
}



#' Based on `psfmi::pool_auc`, pools auc based on vectors of AUC and SE
#' 
#' @param AUC double vector of AUC values
#' @param SE double vector of SE values
#' @param transform should the AUC be transformed before pooling?
#' @param ci.level double between 0 and 1
#' @param m number of imputations
#' 
#' @keywords Internal
#' @noRd
#' 
pool_auc_se <- function (
    AUC, 
    SE, 
    transform, 
    ci.level,
    m
) {
  
  if (transform == "logit") {
    res <- pool_auc_log(AUC, SE, ci.level, m)
  } else if (isFALSE(transform)) {
    res <- pool_auc_unt(AUC, SE, ci.level, m)
  } else {
    stop("Currently only logit-transformation or pooling of untransformed 
         AUC supported.")
  }
  
  return(res)
}



#' Logit transforms and pools AUC
#' 
#' @param AUC double vector of AUC values
#' @param SE double vector of SE values
#' @param ci.level double between 0 and 1
#' @param m number of imputations
#' 
#' @importFrom stats qlogis
#' 
#' @keywords Internal
#' @noRd
#' 
pool_auc_log <- function(AUC, SE, ci.level, m) {
  
  AUC_log <- qlogis(AUC)
  SE_log <- SE / (AUC * (1 - AUC)) # see miceafter::logit_trans
  
  SET <- rr_se(AUC_log, SE_log, ci.level, m)
  
  AUC_log <- mean(AUC_log)
  
  auc <- invlogit(AUC_log)
  ul <- invlogit(AUC_log + (SET[2] * SET[1]))
  ll <- invlogit(AUC_log - (SET[2] * SET[1]))
  
  res <- data.frame(
    auc = auc, 
    ll = ll, 
    ul = ul, 
    auc_logit = AUC_log, 
    var_logit = sqrt(SET[1])
  )
  
  return(res)
  
}


#' Pools untransformed AUC
#' 
#' @keywords Internal
#' @noRd
#' 
pool_auc_unt <- function(AUC, SE, ci.level, m) {
  
  auc <- mean(AUC)
  SET <- rr_se(AUC, SE, ci.level, m)
  ul <- auc + (SET[2] * SET[1])
  if (ul > 1) ul <- 1
  ll <- auc - (SET[2] * SET[1])
  if (ll < 0) ll <- 0

  res <- data.frame(
    auc = auc, ll = ll, ul = ul, 
    var = sqrt(SET[1])
  )
  
  return(res)
  
}



#' Calculates Total Standard Error and t value based on ci.level
#' 
#' @param est double vector of AUC values
#' @param SE double vector of SE values
#' @param ci.level double between 0 and 1
#' @param m number of imputations
#' 
#' @keywords Internal
#' @noRd
#' 
rr_se <- function(est, SE, ci.level, m) {
  
  VW <- mean(SE^2)
  VB <- var(est)
  VT <- VW + (1 + (1/m)) * VB  # VT = total variance
  SET <- sqrt(VT)  # SE total
  riv <- (1 + 1/m) * (VB / VW)  # riv
  v <- (m - 1) * (1 + (1 / riv)) ^ 2
  
  alpha <- 1 - (1 - ci.level) / 2
  t <- qt(alpha, v)  # t value 
  res <- c(SET, t)
  
  return(res)
  
}