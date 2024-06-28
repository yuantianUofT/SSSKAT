# Updated: 2024-06-28

#' test directly with S for both labeled and unlabeled data
#'
#' @param Y Binary response variable
#' @param X Covariates
#' @param G Genotype data
#' @param S Continuous surrogate
#' @param id.t Row id of labeled data in the dataset.
#' @param wBurden Vector of SNP burden weights.
#' @param wSKAT Vector of SNP SKAT weights.
#' @param wACAT Vector of SNP ACAT weights.
#' @param mac.thresh A threshold of the minor allele count (MAC). The Burden test will be used to aggregate the SNPs with MAC less than this threshold.
#' @param testtype Type of test, default is "all", can also be "SKAT", "Burden", "ACAT".
#' @return Naive pvalues.
#' @export

naive_test <- function(X, G, S, id.t, wBurden = NULL, wSKAT = NULL, wACAT = NULL, mac.thresh = 10, testtype = "all") {
  

  if (! testtype %in% c("all", "SKAT", "Burden", "ACAT")) {
    stop("testtype is not correctly specified")
  }
  
  if (testtype == "all") {
    
    # check weights
    if (is.null(wBurden) | is.null(wSKAT) | is.null(wACAT)) {
      stop("weights are not correctly specified")
    }

    # SKAT and Burden NULL model
    obj.S <- SKAT_Null_Model(S ~ X, out_type="C", Adjustment = F)
    # SKAT results
    SKAT_result <- SKAT(G, obj.S, method="davies", weights=wSKAT)
    # Burden results
    Burden_result <- SKAT(G, obj.S, method="davies", r.corr = 1, weights=wBurden)
    # ACAT NULL model
    G <- Matrix::Matrix(G, sparse = TRUE)
    obj <- ACAT::NULL_Model(S, X)
    ACAT_result <- ACAT::ACAT_V(G, obj, weights = wACAT, mac.thresh = mac.thresh)

    # final pvalues
    pvalue <- c(SKAT_p = SKAT_result$p.value, Burden_p = Burden_result$p.value, ACAT_p = ACAT_result)
    # final weights
    weights <- list(wSKAT = wSKAT, wBurden = wBurden, wACAT = wACAT)

  } else if (testtype == "SKAT") {
    
    # check weights
    if (is.null(wSKAT)) {
      stop("weights are not correctly specified")
    }
    
    # SKAT and Burden NULL model
    obj.S <- SKAT_Null_Model(S ~ X, out_type="C", Adjustment = F)
    # SKAT results
    SKAT_result <- SKAT(G, obj.S, method="davies", weights=wSKAT)
    
    # final pvalues
    pvalue <- c(SKAT_p = SKAT_result$p.value)
    # final weights
    weights <- list(wSKAT = wSKAT)
    
  } else if (testtype == "Burden") {
    
    # check weights
    if (is.null(wBurden)) {
      stop("weights are not correctly specified")
    }
    
    # SKAT and Burden NULL model
    obj.S <- SKAT_Null_Model(S ~ X, out_type="D", Adjustment = F)
    # Burden results
    Burden_result <- SKAT(G, obj.S, method="davies", r.corr = 1, weights=wBurden)
    
    # final pvalues
    pvalue <- c(Burden_p = Burden_result$p.value)
    # final weights
    weights <- list(wBurden = wBurden)
    
  } else if (testtype == "ACAT") {
    
    # check weights
    if (is.null(wACAT)) {
      stop("weights are not correctly specified")
    }
    
    # ACAT NULL model
    G <- Matrix::Matrix(G, sparse = TRUE)
    obj <- ACAT::NULL_Model(S, X)
    ACAT_result <- ACAT::ACAT_V(G, obj, weights = wACAT, mac.thresh = mac.thresh)
    
    # final pvalues
    pvalue <- c(ACAT_p = ACAT_result)
    # final weights
    weights <- list(wACAT = wACAT)
    
  }
  
  
  # final results
  results <- list(nobs = nrow(G), pvalue = pvalue, weights = weights)
  
  return(results)
  
}
