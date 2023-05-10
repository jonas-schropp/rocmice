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
#' data.frame (must be binary).
#' @param unique_vals Optional. Define the unique values that should be used as 
#' cutoffs for the score variable. If not supplied, every unique cutoff is used. 
#' @param fpr_vals The FPR values over which to pool the TPR values. Defaults 
#' to `seq(from=0.001, to=0.999, by=0.001)`
#' @param backtransform Should the results be transformed back?
#' @param ci.level Double between 0-1. NULL if no confidence intervals are 
#' desired (or you are going to bootstrap confidence intervals).
#' @param tol tolerance to use for 0 and 1 values (due to logit trans). 1e-3 
#' by default. `FALSE` if no tolerance should be used instead for 0 values.
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
    tol = 0.001
) {
  
  #if (length(target) == 1) {
  #  r <- data[[1]][[target]] 
  #} else if (is.null(target)) {
  #  stop("No target variable supplied.")
  #} 
  
  if (is.null(unique_vals)) unique_vals <- get_unique_vals(data, score)
  
  # Loop over imputations
  cs <- list()
  for (i in 1:length(data)) {
    r <- data[[i]][[target]] 
    nit_vals <- data[[i]][[score]]
    m <- apply_cut_off(nit_vals, unique_vals)
    print(paste0("Dataset ", i))
    cs[[i]] <- apply_metrics(m, r) |> 
      interpol(fpr_vals) |>
      add_var_roc(r) |>
      transform_metrics()
  }
  
  # Combine and pool
  l <- combine_metrics(cs) |> pool_metrics()
  
  if (!is.null(ci.level)) l <- add_ci(l, ci.level, target = r)
  
  df <- data.frame(l)
  names(df) <- paste0(names(df), "_logit")
  
  if (backtransform) {
    df <- backtransform_df(df, ci = !is.null(ci.level))
  }
  
  df$var_roc_logit[is.infinite(df$var_roc_logit)] <- tol
  df$roc_logit[is.infinite(df$roc_logit) & df$roc_logit < 0] <- logit(tol)
  df$roc_logit[is.infinite(df$roc_logit) & df$roc_logit > 0] <- logit(1-tol)
  
  return(df)
  
}



#' Logit Transformation
#' @param mu value to transform
#' @importFrom stats make.link
#' 
#' @keywords Internal
logit <- make.link("logit")$linkfun

#' Logit Transformation inv
#' @param mu value to transform
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
interpol <- function(l, fpr_vals) {
  
  # For CMD checks, try to remove dependencies
  tpri = NULL; vartpri = NULL; fpri = NULL
  
  res <- data.frame(fpri = fpr_vals) %>%
    full_join(as.data.frame(l), by = "fpri") %>%
    arrange(fpri, tpri) %>%
    fill(tpri, vartpri, .direction = "down") %>% 
    mutate(tpri = ifelse(is.na(tpri), 0, tpri)) %>%
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
#' 
#' @keywords Internal
add_var_roc <- function(df, r) {
  
  neg <- sum(r == 0)
  
  df$varroci <- df$vartpri + var_prop(df$fpri, neg)
 
  df 
}


#' @param l list containing 4 elements: a vector of `tpr`, a vector of `fpr` and 
#' the respective variances. The length of each vector is the number of 
#' distinct cutoffs (`unique_vals - 1`)
#' 
#' #' @return 
#' a list with 3 elements:
#' \item{tpri}{logit-transformed TPR values}
#' \item{fpri}{logit-transformed original FPR values}
#' \item{vartpri}{logit-transformed variance for TPR}
#' 
#' @keywords Internal
transform_metrics <- function(l) {
  
  l$vartpri <- varlogit(l$vartpri, l$tpri)
  l$varroci <- varlogit(l$varroci, l$tpri)
  l$tpri <- logit(l$tpri)
  l$fpri <- logit(l$fpri)
  
  l
}

#' Transforms a list of lists to a list of three matrices
#' 
#' @param ll A list of lists. For each imputation, one `l` as described in 
#' `transform_metrics`.
#' 
#' #' @return 
#' a list with 3 elements:
#' \item{tprm}{matrix of logit-transformed TPR values, one column for each 
#' imputation}
#' \item{fprm}{vector of logit-transformed FPR values passed on from `fpr_vals`, 
#' one column for each imputation}
#' \item{vartprm}{matrix of logit-transformed variance for TPR, one column for 
#' each imputation}
#' 
#' @keywords Internal
combine_metrics <- function(ll) {
  
  tprm <- matrix(ncol = length(ll), nrow = length(ll[[1]][[1]]))
  varrocm <- matrix(ncol = length(ll), nrow = length(ll[[1]][[1]]))
  
  for (i in 1:length(ll)) {
    tprm[,i] <- ll[[i]]$tpri
    varrocm[,i] <- ll[[i]]$varroci
  }
  
  return(list(
    tprm = tprm,
    fprm = ll[[1]]$fpri,
    varrocm = varrocm
  ))
  
}


#' @param l list containing: (1) tprm, a matrix of tpr by imp, (2) fprm, a 
#' vector of fpr and (3) varrocm, matrix of the variance of tpr
#' 
#' @keywords Internal
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
#' 
#' #' @return 
#' a data.frame with the elements:
#' \item{tpr_logit}{logit-transformed TPR values}
#' \item{fpr_logit}{logit-transformed FPR values}
#' \item{var_tpr_logit}{logit-transformed variance for TPR}
#' \item{tpr}{TPR}
#' \item{fpr}{FPR}
#' \item{ll_tpr_logit}{lower limit of CI for logit-transformed TPR values}
#' \item{ul_tpr_logit}{upper limit of CI for logit-transformed TPR values}
#' \item{ll_tpr}{lower limit of CI for TPR values}
#' \item{ul_tpr}{upper limit of CI for TPR values}
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
#' @param df data.frame containing tpr, fpr and their variances
#' @param ci.level double between 0-1. 
#' @param target outcome vector
#' 
#' #' @return 
#' a data.frame with the elements:
#' \item{tpr}{logit-transformed TPR values}
#' \item{fpr}{logit-transformed FPR values}
#' \item{var_tpr}{logit-transformed variance for TPR}
#' \item{ll_tpr}{lower limit of CI for logit-transformed TPR values}
#' \item{ul_tpr}{upper limit of CI for logit-transformed TPR values}
#' 
#' @keywords Internal
add_ci <- function(df, ci.level, target) {
  
  alpha <- 1 - (1 - ci.level) / 2
  
  t <- qt(alpha, length(target) - 1)
  
  df$ll_roc <- df$roc - (t * (df$var_roc^2))
  df$ul_roc <- df$roc + (t * (df$var_roc^2))
  
  df$ll_roc[is.nan(df$ll_roc)] <- Inf
  df$ul_roc[is.nan(df$ul_roc)] <- -Inf
  
  df
  
}


