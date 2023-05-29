#' Pool Delong Test after Multiple Imputation
#' 
#' Pools the variance and covariance of two (un-)paired ROC after multiple 
#' imputation to compare their AUC. This function assumes that we did NOT 
#' impute the dependent variable.
#' 
#' @describeIn mi_roc_test
#' 
#' @param data A list of imputed data sets.
#' @param target Character, the name of the outcome variable in data.
#' @param score Character, the name of the predictions variable in data.
#' @param score2 Optional: if you are testing AUC for paired ROC curves and 
#' your data.frames are in wide format, the name of the predictions variable 
#' for your second score in data.
#' @param group Character, the name of the variable that the defines the groups 
#' to compare. Overwritten by `score2` if it is specified.
#' @param groups The levels in `group` to compare. If not specified, the unique 
#' values in `data[[1]]` are used. If `group` has more than two values, or one 
#' of the groups is rare and might not occur in all imputations (which could 
#' bias estimates!), this must be specified. Overwritten by `score2`.
#' @param rocs1 A list of `pROC::roc` objects for each imputed data set. 
#' Overwrites all other arguments. 
#' @param rocs2 A list of `pROC::roc` objects for each imputed data set. 
#' to compare to rocs1. Overwrites all other arguments. 
#' @param paired Logical. Are the ROC curves paired?
#' @param levels Levels of the outcome variable `target`. By default `c(0, 1)`.
#' @param direction In which direction to make the comparison? The default `<` 
#' is different from `pROC` because direction is not automatically checked 
#' (which could lead to weird results after multiple imputation) and means that 
#' predictor values for controls are assumed to be lower than for cases. 
#' 
#' @author 
#' Jonas Schropp, some code from `pROC::roc.test`.
#' 
#' @export
#' 
mi_roc_test <- function(
    data, 
    target, 
    score, 
    score2 = NULL,
    group,
    groups = NULL,
    rocs1 = NULL, 
    rocs2 = NULL, 
    paired = FALSE, 
    levels = c(0, 1),
    direction = '<'
    ) {
  
  if (is.null(direction)) {
    message("Setting direction to '<'.")
    direction <- '<'
  }
  
  if (is.null(levels)) {
    message("Setting levels to c(0, 1).")
    levels <- c(0, 1)
  }
  
  if (!is.null(rocs1)) {
    if (is.null(rocs2)) {
      stop("If rocs1 is supplied, rocs2 must be supplied as well.")
    }
    if (length(rocs1) != length(rocs2)) {
      stop("rocs1 and rocs2 must be the same length")
    }
  }

  
  if (paired) {
    
    if (is.null(rocs1)) {
      
      tmp <- prep_rocs.paired(
        data, target, score, score2, group, groups, levels, direction  
      )
      rocs1 <- tmp[[1]]
      rocs2 <- tmp[[2]]
    }
    
    mi_roc_test.paired(rocs1, rocs2)
    
  } else if (!paired) {
    
    if (is.null(rocs1)) {
      tmp <- prep_rocs.unpaired(data, target, score, group, groups, direction)
      rocs1 <- tmp[[1]]
      rocs2 <- tmp[[2]]
    }
    
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
    auc1 = aucR,
    auc2 = aucS,
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
    
    # C++ function from pROC
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
    auc1 = aucR,
    auc2 = aucS,
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
  
  return(placements)
  
}




#' Extract cases and controls for two groups 
#' 
#' Helper function to prepare data for delongPlacements
#' 
#' @param d a data.frame
#' @param target character, the outcome variable
#' @param score character, the prediction variable
#' @param group character, the group variable
#' @param groups vector with the two group levels to compare
#' @param m imputation number
#' 
#' @keywords Internal
#' 
cases_controls.long <- function(d, target, score, group, groups, m) {
  
  roc1 <- list()
  roc2 <- list()
  
  if (any(groups != sort(unique(d[[group]])))) {
    warning("Different groups in imputation ", m, " detected. 
            Results might be incorrect.")
  }
  
  roc1$cases <- d[d[[group]] == groups[1] & d[[target]] == 1, ][[score]]
  roc1$controls <- d[d[[group]] == groups[1] & d[[target]] == 0, ][[score]]
  roc2$cases <- d[d[[group]] == groups[2] & d[[target]] == 1, ][[score]]
  roc2$controls <- d[d[[group]] == groups[2] & d[[target]] == 0, ][[score]]
  
  return(list(roc1, roc2))
  
}



#' Extract cases and controls for two groups 
#' 
#' Helper function to prepare data for delongPlacements
#' 
#' @param d a data.frame
#' @param target character, the outcome variable
#' @param score character, the prediction variable for group 1
#' @param score2 character, the prediction variable for group 2
#' 
#' @keywords Internal
#' 
cases_controls.wide <- function(d, target, score, score2) {
  
  roc1 <- list()
  roc2 <- list()
  
  roc1$cases <- d[d[[target]] == 1, ][[score]]
  roc1$controls <- d[d[[target]] == 0, ][[score]]
  roc2$cases <- d[d[[target]] == 1, ][[score2]]
  roc2$controls <- d[d[[target]] == 0, ][[score2]]
  
  return(list(roc1, roc2))
  
}




#' Prepares list of ROCs for paired case
#' 
#' @param data see mi_roc_test
#' @param target see mi_roc_test
#' @param score see mi_roc_test
#' @param score2 see mi_roc_test
#' @param group see mi_roc_test
#' @param groups see mi_roc_test
#' @param levels see mi_roc_test
#' @param direction see mi_roc_test
#' 
#' @keywords Internal
prep_rocs.paired <- function(
    data, target, score, score2, group, groups, levels, direction  
  ) {
  
  rocs1 <- list()
  rocs2 <- list()
  
  if (is.null(score2)) {
    if (is.null(groups)) {
      groups <- unique(data[[1]][[group]])
      if (length(groups) != 2) stop("Can only compare exactly 2 groups.")
      message("Setting `groups` to ", groups[1], " and ", groups[2], ".")
    } else {
      if (length(groups) != 2) stop("Can only compare exactly 2 groups.")
    }
    for (m in 1:length(data)){
      
      tmp <- cases_controls.long(data[[m]], target, score, group, groups, m)
      rocs1[[m]] <- tmp[[1]]
      rocs2[[m]] <- tmp[[2]]
      rocs1[[m]]$direction <- direction
      rocs2[[m]]$direction <- direction
      
    }
  } else {
    
    for (m in 1:length(data)){
      
      tmp <- cases_controls.wide(data[[m]], target, score, score2)
      rocs1[[m]] <- tmp[[1]]
      rocs2[[m]] <- tmp[[2]]
      rocs1[[m]]$direction <- direction
      rocs2[[m]]$direction <- direction
      
    }
  }
  
  return(list(rocs1, rocs2))
  
}


#' Prepares list of ROCs for unpaired case
#' 
#' @param data see mi_roc_test
#' @param target see mi_roc_test
#' @param score see mi_roc_test
#' @param group see mi_roc_test
#' @param groups see mi_roc_test
#' @param direction see mi_roc_test
#' 
#' @keywords Internal
prep_rocs.unpaired <- function(data, target, score, group, groups, direction) {
  
  rocs1 <- list()
  rocs2 <- list()
  
  if (is.null(groups)) {
    groups <- unique(data[[1]][[group]])
    if (length(groups) != 2) stop("Can only compare exactly 2 groups.")
    message("Setting `groups` to ", groups[1], " and ", groups[2], ".")
  } else {
    if (length(groups) != 2) stop("Can only compare exactly 2 groups.")
  }
  
  for (m in 1:length(data)){
    
    tmp <- cases_controls.long(data[[m]], target, score, group, groups, m)
    rocs1[[m]] <- tmp[[1]]
    rocs2[[m]] <- tmp[[2]]
    rocs1[[m]]$direction <- direction
    rocs2[[m]]$direction <- direction
    
  }
  
  return(list(rocs1, rocs2))
  
}