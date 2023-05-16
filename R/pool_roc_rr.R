#' Pool ROC curves based on multiple imputation
#' 
#' The function `pool_roc_rr` pools roc curves created from a list of multiply 
#' imputed data sets by applying Rubin's Rules to the logit-transformed 
#' TPR and FPR at each unique value of the supplied score. 
#' The function returns a data.frame of transformed TPR and FPR at each cutoff, 
#' their respective variances and optionally confidence intervals as well as 
#' backtransformed values.
#' 
#' @param data A list of data.frames
#' @param score Character, the name of the score variable that should be present 
#' in each data.frame in data
#' @param target The name of the target variable present in each 
#' data.frame  in data (must be binary).
#' @param unique_vals Optional. Define the unique values that should be used as 
#' cutoffs for the score variable. If not supplied, every unique cutoff is used. 
#' @param fpr_vals The FPR values over which to pool the TPR values. Defaults 
#' to `seq(from=0.001, to=0.999, by=0.001)`
#' @param backtransform Should the results be transformed back?
#' @param ci.level Double between 0-1. NULL if no confidence intervals are 
#' desired (or you are going to bootstrap confidence intervals).
#' @param corr continuity correction. `FALSE` if no continuity correction should 
#' be used instead for 0 values.
#' @param verbose TRUE/FALSE, should information be printed to the screen while 
#' the algorithm is running?
#'  
#'  
#' @return 
#' a data.frame with the elements:
#' \item{tpr_logit}{logit-transformed TPR values}
#' \item{fpr_logit}{logit-transformed FPR values}
#' \item{var_tpr_logit}{description a}
#' \item{tpr}{description a}
#' \item{fpr}{description a}
#' \item{ll_tpr_logit}{description a}
#' \item{ul_tpr_logit}{description a}
#' 
#' @export
#' 
pool_roc_rr <- function(
    data,
    score = "fib4",
    target = "advanced_fibrosis",
    unique_vals = NULL,
    fpr_vals = seq(from=0.001, to=0.999, by=0.001),
    backtransform = FALSE,
    ci.level = NULL,
    corr = 0.1,
    verbose = TRUE
) {
  
  #if (length(target) == 1) {
  #  r <- data[[1]][[target]] 
  #} else if (is.null(target)) {
  #  stop("No target variable supplied.")
  #} 
  
  if (is.null(unique_vals)) unique_vals <- get_unique_vals(data, score)
  fpr_vals <- qlogis(fpr_vals)
  
  
  # Loop over imputations
  cs <- list()
  if(verbose) {
    cat("\nCalculating TPR for every FPR: \n")
    pb <- txtProgressBar(min = 0, max = length(data), initial = 0) 
  }
  
  for (i in 1:length(data)) {
    
    r <- data[[i]][[target]] 
    p <- data[[i]][[score]]
    m <- apply_cut_off(p, unique_vals)
    zero <- qlogis(corr / (sum(r == 1) + 2*corr))
    
    cs[[i]] <- apply_metrics(m, r, corr) |> 
      interpol(fpr_vals, zero)
    
    if (verbose) setTxtProgressBar(pb,i)
  }
  
  if (verbose) cat("\nCombining and pooling results.\n")
  
  # Combine and pool
  l <- combine_metrics(cs) |> pool_metrics()
  
  if (!is.null(ci.level)) l <- add_ci(l, ci.level, target = r)
  
  df <- data.frame(l)
  names(df) <- paste0(names(df), "_logit")
  
  if (backtransform) {
    df <- backtransform_df(df, ci = !is.null(ci.level))
  }
  
  return(df)
  
}



#' Logit Transformation inv
#' @param eta value to transform
#' @importFrom stats make.link
#' @keywords Internal
invlogit <- make.link("logit")$linkinv





#' Interpolate values for tpri 
#' 
#' Unlike the method for meta-analysis by Martínez-Camblor, P. (2016), this 
#' function does not use linear interpolation because we have enough values 
#' for tpr available to simply use the values in the data.
#' 
#' @param l a list with 3 elements: `tpri`, `fpri`, `vartpri` as returned by 
#' `apply_metrics`
#' @param fpr_vals the FPR values to interpolate for 
#' @param zero value to use for 0 cells
#' 
#' @return 
#' a data.frame with 3 elements:
#' \item{tpri}{logit-transformed TPR values}
#' \item{fpri}{logit-transformed FPR values passed on from `fpr_vals`}
#' \item{vartpri}{logit-transformed variance for TPR}
#' 
#' @import dplyr
#' @importFrom tidyr fill
#' @keywords Internal
interpol <- function(l, fpr_vals, zero) {
  
  # For CMD checks, try to remove dependencies
  tpri = NULL; fpri = NULL; vartpri = NULL; varfpri = NULL
  
  res <- data.frame(fpri = fpr_vals) %>%
    full_join(as.data.frame(l), by = "fpri") %>%
    arrange(fpri, tpri) %>%
    fill(tpri, vartpri, varfpri, .direction = "down") %>% 
    mutate(tpri = ifelse(is.na(tpri), zero, tpri)) %>%
    distinct() %>%
    filter(fpri %in% fpr_vals) %>%
    group_by(fpri) %>%
    filter(tpri == max(tpri)) %>%
    ungroup()
    
  res
  
}


#' Adds roc variance to df following the formula in Martinez
#' 
#' varroc = vartpr + varfpr
#' 
#' @param df data.frame with tpr, fpr, vartpr
#' @param r response vector
#' @param p prediction vector
#' @param corr continuity correction
#' 
#' @keywords Internal
add_var_roc <- function(df, r, p, corr) {
  
  neg <- sum(r == 0) + 2*corr
  pos <- sum(r == 1) + 2*corr
  
  # new
  #eneg <- sum(r == 0 & p == 1) + corr
  #epos <- sum(r == 1 & p == 1) + corr
  
  #vartpri <- varlogit2(df$tpri, pos)
  #varfpri <- varlogit2(df$fpri, neg)
  
  vartpri <- var_prop(df$tpri, pos)
  varfpri <- var_prop(df$fpri, neg)
  df$varroci <- vartpri + varfpri
 
  df 
}


#' transform metrics
#' 
#' @param l list containing 4 elements: a vector of `tpr`, a vector of `fpr` and 
#' the respective variances. The length of each vector is the number of 
#' distinct cutoffs (`unique_vals - 1`)
#' 
#' #' @return 
#' a list with 3 elements:
#' \itemize{
#'   \item{tpri}{logit-transformed TPR values}
#'   \item{fpri}{logit-transformed original FPR values}
#'   \item{vartpri}{logit-transformed variance for TPR}
#' }
#' 
#' @keywords Internal
transform_metrics <- function(l) {
  
  # removed new
  l$varroci <- varlogit(l$varroci, l$tpri)
  l$tpri <- qlogis(l$tpri)
  l$fpri <- qlogis(l$fpri)
  
  l
}

#' Transforms a list of lists to a list of three matrices
#' 
#' @param ll A list of lists. For each imputation, one `l` as described in 
#' `transform_metrics`.
#' 
#' #' @return 
#' a list with 3 elements:
#' \itemize{
#'   \item{tprm}{matrix of logit-transformed TPR values, one column for each 
#' imputation}
#'   \item{fprm}{vector of logit-transformed FPR values passed on from `fpr_vals`, 
#' one column for each imputation}
#'   \item{vartprm}{matrix of logit-transformed variance for TPR, one column for 
#' each imputation}
#' }
#' 
#' @keywords Internal
combine_metrics <- function(ll) {
  
  tprm <- matrix(ncol = length(ll), nrow = length(ll[[1]][[1]]))
  varrocm <- matrix(ncol = length(ll), nrow = length(ll[[1]][[1]]))
  
  for (i in 1:length(ll)) {
    tprm[,i] <- ll[[i]]$tpri
    varrocm[,i] <- ll[[i]]$vartpri
    #varfprm[,i] <- ll[[i]]$varfpri
  }
  
  return(list(
    tprm = tprm,
    fprm = ll[[1]]$fpri,
    varrocm = varrocm
  ))
  
}


#' Pool metrics
#' 
#' @param l list containing: (1) tprm, a matrix of tpr by imp, (2) fprm, a 
#' vector of fpr and (3) varrocm, matrix of the variance of tpr
#' 
#' @keywords Internal
#' @noRd
pool_metrics <- function(l) {
  
  tprm <- pool_rr_m(l$tprm, l$varrocm)
  tprm$fpr <- l$fprm
  
  return(tprm)
  
}



#' Gets all unique values of a score in ascending order
#' 
#' @param data list of imputed data sets
#' @param score character, name of score variable in data
#' @param digits FALSE if values should not be rounded, otherwise any number 
#' of digits.
#' 
#' @keywords Internal
get_unique_vals <- function(data, score, digits = 3) {
  
  res <- double()
  
  for (i in 2:length(data)) res <- c(res, data[[i]][[score]])
  
  if (digits) res <- round(res, digits)
  
  sort(unique(res))
  
}




#' Adds backtransformed tpr, fpr and optionally confidence intervals to df
#' 
#' @param df data.frame containing tpr, fpr and their variances
#' @param ci Logical, is there a ci to backtransform
#' 
#' #' @return 
#' a data.frame with the elements:
#' \itemize{
#'   \item{tpr_logit}{logit-transformed TPR values}
#'   \item{fpr_logit}{logit-transformed FPR values}
#'   \item{var_tpr_logit}{logit-transformed variance for TPR}
#'   \item{tpr}{TPR}
#'   \item{fpr}{FPR}
#'   \item{ll_tpr_logit}{lower limit of CI for logit-transformed TPR values}
#'   \item{ul_tpr_logit}{upper limit of CI for logit-transformed TPR values}
#'   \item{ll_tpr}{lower limit of CI for TPR values}
#'   \item{ul_tpr}{upper limit of CI for TPR values}
#' }
#' 
#' @keywords Internal
backtransform_df <- function(df, ci) {
  
  df$roc <- invlogit(df$roc_logit)
  df$fpr <- invlogit(df$fpr_logit)
  
  if (ci) {
    df$ll_roc <- invlogit(df$ll_roc_logit)
    df$ul_roc <- invlogit(df$ul_roc_logit)
  }
  
  df
  
}

#' Adds lower and upper bounds of Confidence Interval to df
#' 
#' @param df data.frame containing:
#' \itemize{
#'   \item{tpr}{logit-transformed TPR values}
#'   \item{fpr}{logit-transformed FPR values}
#'   \item{var_tpr}{logit-transformed variance for TPR}
#' }
#' @param ci.level double between 0-1. 
#' @param target outcome vector
#' 
#' #' @return 
#' a data.frame with the elements:
#' \itemize{
#'   \item{tpr}{logit-transformed TPR values}
#'   \item{fpr}{logit-transformed FPR values}
#'   \item{var_tpr}{logit-transformed variance for TPR}
#'   \item{ll_tpr}{lower limit of CI for logit-transformed TPR values}
#'   \item{ul_tpr}{upper limit of CI for logit-transformed TPR values}
#' }
#' @keywords Internal
add_ci <- function(df, ci.level, target) {
  
  alpha <- 1 - (1 - ci.level) / 2
  
  t <- qt(alpha, length(target) - 1)
  
  df$ll_roc <- df$roc - (t * sqrt(df$var_roc))
  df$ul_roc <- df$roc + (t * sqrt(df$var_roc))
  
  df
  
}








