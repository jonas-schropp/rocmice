#' Pool Delong Test after Multiple Imputation
#' 
#' Pools the variance and covariance of two (un-)paired ROC after multiple 
#' imputation to compare their AUC. This function assumes that we did NOT 
#' impute the dependent variable.
#' 
#' @describeIn mi_roc_test
#' 
#' @param rocs1 A list of "roc" objects for each imputed data set
#' @param rocs2 A list of "roc" objects for each imputed data set 
#' to compare to rocs1
#' @param paired Are the roc curves paired?
#' 
#' @author 
#' Jonas Schropp, code fragments from `pROC::roc.test` and `MKmisc::mi.t.test`.
#' 
#' @export
#' 
mi_roc_test <- function(rocs1, rocs2, paired = FALSE) {
  
  if (length(rocs1) != length(rocs2)) {
    stop("rocs1 and rocs2 must be the same length")
  }
  
  if (paired) {
    mi_roc_test.paired(rocs1, rocs2)
  } else if (!paired) {
    mi_roc_test.unpaired(rocs1, rocs2)
  }
  
}



#' @describeIn mi_roc_test
#' 
#' @param rocs1 A list of "roc" objects for each imputed data set
#' @param rocs2 A list of "roc" objects for each imputed data set 
#' to compare to rocs1
#' 
#' @details 
#' Unpaired version
#' 
#' @importFrom stats var
#' @importFrom stats pt
#' @importFrom stats pnorm
#' 
#' @keywords Internal
#' 
mi_roc_test.unpaired <- function(rocs1, rocs2) {
  
  # Number of mi data sets
  m <- length(rocs1)
  # Variances for each AUC
  SR <- double(length = m)
  SS <- double(length = m)
  # Each AUC
  aucsR <- double(length = m)
  aucsS <- double(length = m)
  
  for (i in 1:m) {
    
    nR <- length(rocs1[[i]]$controls)
    mR <- length(rocs1[[i]]$cases)
    
    nS <- length(rocs2[[i]]$controls)
    mS <- length(rocs2[[i]]$cases)
    
    VR <- delongPlacements(rocs1[[i]])
    VS <- delongPlacements(rocs2[[i]])
    
    aucsR[i] <- VR$theta 
    aucsS[i] <- VS$theta
    
    SRX <- sum((VR$X - VR$theta) * (VR$X - VR$theta))/(mR-1)
    SSX <- sum((VS$X - VS$theta) * (VS$X - VS$theta))/(mS-1)
    
    SRY <- sum((VR$Y - VR$theta) * (VR$Y - VR$theta))/(nR-1)
    SSY <- sum((VS$Y - VS$theta) * (VS$Y - VS$theta))/(nS-1)
    
    SR[i] <- SRX/mR + SRY/nR   # Variance R identical with var(roc1)
    SS[i] <- SSX/mS + SSY/nS   # Variance S identical with var(roc2)
    
  }
  
  # Calculate Covariances
  SSR <- sqrt((SR) + (SS))
  
  # Within-imputation variance
  VWR <- mean(SR) # ubar
  VWS <- mean(SS)
  VWSR <- mean(SSR)
  
  # Between-imputation variance
  VBR <- var(aucsR) # b
  VBS <- var(aucsS)
  VBSR <- var(aucsR - aucsS)
  
  # Total Variance
  VTR <- VWR + (1 + 1/m)*VBR # t in mice
  VTS <- VWS + (1 + 1/m)*VBS
  VTSR <- VWSR + (1 + 1/m)*VBSR
  
  # Average AUC
  aucR <- mean(aucsR) # qbar
  aucS <- mean(aucsS)
  thetadiff <- mean(aucsR - aucsS)
  
  # riv: Relative increase in variance due to nonresponse
  rivR <- (1 + 1/m) * VBR/VWR
  rivS <- (1 + 1/m) * VBS/VWS
  rivSR <- (1 + 1/m) * VBSR/VWSR
  
  # lambda: Proportion of total variance due to missingness
  lambdaR <- (1 + 1/m) * VBR/VTR
  lambdaS <- (1 + 1/m) * VBS/VTS
  lambdaSR <- (1 + 1/m) * VBSR/VTSR
  
  # Assumes that we did not impute the dependent variable / reference
  # therefore n are equal in each iteration and we can just use the last one
  ntotR <- nR + mR
  ntotS <- nS + mS
  
  # pooled variance within each imputation
  v0s <- ((SR) + (SS))^2 / 
    (((SR)^2 / (ntotR - 1)) + ((SS)^2 / (ntotS - 1 )))
  
  # inspired by MKmisc::mi.t.test
  vm <- (m-1)*(1 + VWSR / ((1 + 1/m) * VBSR))^2
  v0 <- mean(v0s) 
  dfobs <- (v0 + 1) / (v0 + 3)*v0*(1 - lambdaSR) # barnard.rubin
  
  # Overall df
  df.mod <- 1 / (1/vm + 1/dfobs)
  
  # And t value
  t <- thetadiff / VTSR
  
  # fmi: Fraction of missing information
  fmiSR <- (rivSR + 2 / (df.mod + 3)) / (rivSR + 1)
  
  # p.value
  p <- 2*pt(-abs(t), df = df.mod)#pt(t, df.mod)
  
  
  res <- data.frame(
    delta_auc = thetadiff, 
    t.value = t, 
    #df = df.mod, 
    p.value = p,
    var.total = VTSR, 
    var.within = VWSR, 
    var.between = VBSR, 
    riv = rivSR, 
    lambda = lambdaSR, 
    fmi = fmiSR
    )
  
  return(res)
  
}





#' @describeIn mi_roc_test
#' 
#' @param rocs1 A list of "roc" objects for each imputed data set
#' @param rocs2 A list of "roc" objects for each imputed data set 
#' to compare to rocs1
#' 
#' @details 
#' Paired version
#' 
#' @importFrom stats var
#' @importFrom stats pt
#' 
#' @keywords Internal
#' 
mi_roc_test.paired <- function(rocs1, rocs2) {
  
  # Number of mi data sets
  m <- length(rocs1)
  
  # Common Variance for each imp
  sig <- double(length = m)  # sigma
  SR <- double(length = m)
  SS <- double(length = m)
  
  # Each AUC
  aucsR <- double(length = m)
  aucsS <- double(length = m)
  
  # Extract numbers of cases and controls for both
  # For paired tests, these are identical
  # And since we did not impute the dependent variable, identical for each m
  nn <- length(rocs1[[1]]$controls) 
  mm <- length(rocs1[[1]]$cases)    
  
  for (i in 1:m) {
    
    # C function from pROC, rewrite in C++ for CMD checks???
    VR <- delongPlacements(rocs1[[i]])
    VS <- delongPlacements(rocs2[[i]])
    
    aucsR[i] <- VR$theta
    aucsS[i] <- VS$theta
    
    # Calculate variance and covariance 
    SX <- matrix(NA, ncol = 2, nrow = 2)
    SX[1, 1] <- sum((VR$X - aucsR[i]) * (VR$X - aucsR[i]))/(mm - 1)
    SX[1, 2] <- sum((VR$X - aucsR[i]) * (VS$X - aucsS[i]))/(mm - 1)
    SX[2, 1] <- sum((VS$X - aucsS[i]) * (VR$X - aucsR[i]))/(mm - 1)
    SX[2, 2] <- sum((VS$X - aucsS[i]) * (VS$X - aucsS[i]))/(mm - 1)
    SY <- matrix(NA, ncol = 2, nrow = 2)
    SY[1, 1] <- sum((VR$Y - aucsR[i]) * (VR$Y - aucsR[i]))/(nn - 1)
    SY[1, 2] <- sum((VR$Y - aucsR[i]) * (VS$Y - aucsS[i]))/(nn - 1)
    SY[2, 1] <- sum((VS$Y - aucsS[i]) * (VR$Y - aucsR[i]))/(nn - 1)
    SY[2, 2] <- sum((VS$Y - aucsS[i]) * (VS$Y - aucsS[i]))/(nn - 1)
    
    # Variance - Covariance matrix
    S <- SX/mm + SY/nn
    L <- c(1, -1)
    
    # Variance
    SR[i] <- S[1,1]
    SS[i] <- S[2,2]
    
    # sigma se / sd
    sig[i] <- sqrt(L %*% S %*% L)
    
  }
  
  d <- aucsR - aucsS
  
  # Within-imputation variance
  VWR <- mean(SR) # ubar
  VWS <- mean(SS)
  VWSR <- mean(sig^2)
  
  # Between-imputation variance
  VBR <- var(aucsR) # b
  VBS <- var(aucsS)
  VBSR <- var(d)
  
  # Total Variance
  VTR <- VWR + (1 + 1/m)*VBR # t in mice
  VTS <- VWS + (1 + 1/m)*VBS
  VTSR <- VWSR + (1 + 1/m)*VBSR
  
  # Average AUC
  aucR <- mean(aucsR) # qbar
  aucS <- mean(aucsS)
  thetadiff <- mean(d)
  
  # riv: Relative increase in variance due to nonresponse
  rivR <- (1 + 1/m) * VBR/VWR
  rivS <- (1 + 1/m) * VBS/VWS
  rivSR <- (1 + 1/m) * VBSR/VWSR
  
  # lambda: Proportion of total variance due to missingness
  lambdaR <- (1 + 1/m) * VBR/VTR
  lambdaS <- (1 + 1/m) * VBS/VTS
  lambdaSR <- (1 + 1/m) * VBSR/VTSR
  
  # Assumes that we did not impute the dependent variable / reference
  # therefore n are equal in each iteration and we can just use the last one
  ntot <- nn + mm
  
  # dfs within each imputation
  v0s <- ((SR) + (SS))^2 / ((SR^2 + SS^2)^2 / (ntot - 1))
  
  # inspired by MKmisc::mi.t.test
  vm <- (m-1)*(1 + VWSR / ((1 + 1/m) * VBSR))^2
  v0 <- mean(v0s) 
  dfobs <- (v0 + 1) / (v0 + 3)*v0*(1 - lambdaSR) # barnard.rubin
  
  # Overall df
  df.mod <- 1 / (1/vm + 1/dfobs) # for t
  
  # And Z value
  t <- thetadiff / VTSR # compare to z
  zscore <- thetadiff / sqrt(VTSR)
  if (is.nan(zscore) && thetadiff == 0 && VTSR == 0) zscore <- 0
  
  # fmi: Fraction of missing information
  fmiSR <- (rivSR + 2 / (df.mod + 3)) / (rivSR + 1)
  
  # p.value
  #p <- pt(t, df.mod)
  pval <- 2 * pnorm(-abs(zscore))
  
  #if(conf.level > 0) {
  #  crit_z <- qnorm(1 - ((1 - conf.level)/2))
  #  out <- list()
  #  ul <- with(calcs, thetadiff + crit_z * VTSR)
  #  ll <- with(calcs, thetadiff - crit_z * VTSR)
  #}
  
  res <- data.frame(
    delta_auc = thetadiff, 
    Z = zscore, 
    #df = df.mod, too high? why?
    #p.value.t = p,
    p.value = pval,
    var.total = VTSR, 
    var.within = VWSR, 
    var.between = VBSR, 
    riv = rivSR, 
    lambda = lambdaSR, 
    fmi = fmiSR
  )
  
  return(res)
  
}


#' delongPlacements
#' 
#' Shamelessly stolen from pROC
#' 
#' @param roc pROC::roc object
#' @noRd
#' @keywords Internal
#' 
delongPlacements <- function(roc) {
  
  placements <- delongPlacementsCpp(roc)
  auc <- roc$auc/ifelse(roc$percent, 100, 1)
  if (!isTRUE(all.equal(placements$theta, auc))) {
    sessionInfo <- sessionInfo()
    save(roc, placements, sessionInfo, file = "pROC_bug.RData")
    stop(
      sprintf(
        "pROC: error in calculating DeLong's theta: got %.20f instead of %.20f. Diagnostic data saved in pROC_bug.RData. Please report this bug to <%s>.", 
        placements$theta, auc, utils::packageDescription("pROC")$BugReports))
  }
  
  return(placements)
  
}




