#' Calculate ROC objects for a list of data sets
#' 
#' Simple convenience function to apply `pROC::roc` over a list of data sets with 
#' multiple imputations. You could do it using `map` or `apply` or `with`,
#' but this is just a bit easier. It does not allow for smoothed ROC curves.
#' 
#' For more details on the arguments see `pROC::roc`.
#' 
#' @param milist List of data.frames
#' @param selection either "all" if roc object should be returned for every 
#' element in milist, or an integer vector containing the indices of the elements 
#' in milist to use (for example when the first element contains the original 
#' data with missing values)
#' @param score Character, the name of the predictor variable that must be 
#' available in each data.frame in milist.
#' @param target Character, the name of the response variable that must be 
#' available in each data.frame in milist.
#' @param direction either "<" or ">". "auto" is not supported because it might 
#' be different by imputation. By default ">".
#' @param levels the value of the response for controls and cases respectively. 
#' By default, the first two values of levels(as.factor(response)) are taken, 
#' and the remaining levels are ignored. It usually captures two-class factor 
#' data correctly, but will frequently fail for other data types (response 
#' factor with more than 2 levels, or for example if your response is coded 
#' “controls” and “cases”, the levels will be inverted) and must then be 
#' specified here. If your data is coded as 0 and 1 with 0 being the controls, 
#' you can safely omit this argument. With several data frames it is recommended 
#' you set this value explicitly, otherwise it will be determined for each data 
#' set separately.
#' @param algorithm the method used to compute sensitivity and specificity, an 
#' integer of length 1 between 0 and 6. 1: a safe, well-tested, pure-R code that 
#' is efficient when the number of thresholds is low. It goes with O(T*N). 2: 
#' an alternative pure-R algorithm that goes in O(N). Typically faster than 1 
#' when the number of thresholds of the ROC curve is above 1000. Less tested 
#' than 1. 3: a C++ implementation of 1, about 3-5x faster. Typically the 
#' fastest with ROC curves with less than 50-100 thresholds, but has a very 
#' bad worst-case when that number increases. 4 (debug only, slow): runs 
#' algorithms 1 to 3 and makes sure they return the same values. 5: select 2 
#' or 3 based on the number of thresholds. 6 (default): quickly select the 
#' algorithm on the class of the data: 2 for numeric and 3 for ordered. 0: 
#' use microbenchmark to choose between 2 and 3. By default 6.
#' @param quiet set to TRUE to turn off messages when levels is auto-detected.
#' @param auc compute the area under the curve (AUC)? TRUE by default.
#' 
#' @returns A list of `pROC::roc` objects
#' 
#' @importFrom pROC roc_
#' 
#' @export
#' 
apply_roc <- function(
    milist, 
    selection,
    score,
    target,
    direction = ">",
    levels = c(0, 1),
    algorithm = 6,
    quiet = FALSE,
    auc = TRUE
) {
  
  if (typeof(milist) != "list") {
    stop("Milist must be a list of data.frames with multiple imputations.")
  } 
  
  if (selection == "all") {
    ids <- 1:length(milist)
  } else if (!is.numeric(selection)) {
    ids <- as.integer(selection)
    stop("Selection must be either a vector of integers")
  }
  
  if (is.numeric(selection) & isFALSE(all.equal(selection, ids))) {
    stop("Selection appears to contain double-valued ids")
  }
  
  if ((any(!(ids %in% 1:length(milist))))) {
    stop("Some selected indices are not present in milist. Check 'selection'.")
  }
  
  if (direction != "<" & direction != ">") {
    stop("Direction must be either < or >.")
  }
  
  # initialize list
  rocs <- list()
  
  # then simply loop over provided ids
  for (i in ids) {
    
    rocs[[i]] <- pROC::roc_(
      data = milist[[i]],
      response = target,
      predictor = score,
      direction = direction,
      levels = levels,
      algorithm = algorithm,
      quiet = quiet,
      auc = auc
    )
    
  }
  
  return(rocs)
  
  
}