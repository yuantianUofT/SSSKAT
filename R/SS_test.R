# Updated: 2024-06-27

#' Semi-supervised association test function
#'
#' @param Y Binary response variable
#' @param X Covariates
#' @param G Genotype data
#' @param S Continuous surrogate
#' @param id.t Row id of labeled data in the dataset.
#' @param para_results Saved parametric boostrapped individual c value matrix.
#' @param wBurden Vector of SNP burden weights.
#' @param wSKAT Vector of SNP SKAT weights.
#' @param wACAT Vector of SNP ACAT weights.
#' @param weights.beta Weights of the SNPs, default NULL treat the SNPs equally, else take parameters of beta distribution.
#' @param mac.thresh A threshold of the minor allele count (MAC). The Burden test will be used to aggregate the SNPs with MAC less than this threshold.
#' @param full_NR_evaluation Full optimization iteration, default to be TRUE.
#' @param nit Number of iteration for optimization if full_NR_evaluation is FALSE.
#' @param NULL_nlog_like Function to be optimized, predefined, no need to change
#' @param testtype Type of test, default is "all", can also be "SKAT", "Burden", "ACAT".
#' @param boot Whether perform parametric bootstrap, default is TRUE.
#' @param distri Distribution of S|Y, either "normal" or "beta"
#' @param theta Value of estimated parameters, default is NULL.
#' @return Total sample size, parameter estimates, SS scores, SS score variance, SS pvalues, weights used in each test and weights.beta.
#' @export


SS_test <- function(Y, X, G, S, id.t, para_results, 
                    wBurden = NULL, wSKAT = NULL, wACAT = NULL, weights.beta = NULL, mac.thresh = 10,
                    full_NR_evaluation = TRUE, nit = NULL, NULL_nlog_like, 
                    testtype = "all", boot = T, distri, theta = NULL) {
  if (! testtype %in% c("all", "SKAT", "Burden", "ACAT")) {
    stop("testtype is not correctly specified")
  }
  # data preparation
  Z <- data.matrix(cbind(1, X))
  if (testtype %in% c("all", "ACAT")) {
    is.very.rare <- rare_func(G, mac.thresh)
  }
  
  # cvalue estimates
  if (is.null(theta)) {
    theta <- ssl_theta(Y = Y, X = X, S = S, Z = Z, id.t = id.t, distri = distri, weights = NULL, full_eval = TRUE, NULL_nlog_like, nit)$final_est
  }
  cvalue <- c_func(Y, X, S, Z, id.t, theta = theta, distri = distri)
  
  # weights
  if (is.null(wBurden) & is.null(wSKAT) & is.null(wACAT)) {
    wBurden <- SKATBurdenW_func(G, weights.beta)
    wSKAT <- wBurden
    wACAT <- ACATW_func(G, is.very.rare, wBurden, weights.beta)
  } else if ((!is.null(wBurden)) & is.null(wSKAT) & is.null(wACAT)) {
    wSKAT <- wBurden
    wACAT <- ACATW_func(G, is.very.rare, wBurden, weights.beta)
  } else if (is.null(wBurden) & (!is.null(wSKAT)) & is.null(wACAT)) {
    wBurden <- wSKAT
    wACAT <- ACATW_func(G, is.very.rare, wBurden, weights.beta)
  } else if (is.null(wBurden) & is.null(wSKAT) & (!is.null(wACAT))) {
    wBurden <- SKATBurdenW_func(G, weights.beta)
    wSKAT <- wBurden
  } else if (is.null(wBurden) & (!is.null(wSKAT)) & (!is.null(wACAT))) {
    wBurden <- wSKAT
  } else if ((!is.null(wBurden)) & is.null(wSKAT) & (!is.null(wACAT))) {
    wSKAT <- wBurden
  } else if ((!is.null(wBurden)) & (!is.null(wSKAT)) & is.null(wACAT)) {
    wACAT <- ACATW_func(G, is.very.rare, wBurden, weights.beta)
  }
  
  if (testtype == "all") {
    
    # scores
    SKATscore <- SKATQ_fun(G, cvalue, wSKAT)
    Burdenscore <- BurdenQ_fun(G, cvalue, wBurden)
    ACATscore <- ACATsingleQ_fun(G, cvalue)
    scoreQ <- list(SKATscoreQ = SKATscore$Q, BurdenscoreQ = Burdenscore$Q, ACATscoreQs = ACATscore$Q)
    scores <- list(SKATscores = SKATscore$score, Burdenscores = Burdenscore$score, ACATscores = ACATscore$score)
    
    # Scovs
    if (boot == T) {
      Scovs <- Var_boot(G, para_cvalue = para_results, is.very.rare, mac.thresh, wBurden, wSKAT, type=testtype)
      SKAT_Scov <- Scovs$SKATSvar
      Burden_Scov <- Scovs$BurdenSvar
      ACAT_Scov <- Scovs$ACATSvar
    } else {
      
    }
    
    # pvalues
    Liu_result <- Liu(A = diag(ncol(SKAT_Scov)), X_Sigma = SKAT_Scov, X_mu = rep(0, ncol(SKAT_Scov)), X = unlist(SKATscore$score))
    SKAT_p <- unlist(Liu_result$pvalue)
    Burden_p <- Burden(S = Burdenscore$score, sigma = Burden_Scov)
    ACATpvalues <- mapply(Burden, ACATscore$score, ACAT_Scov)
    is.keep <- rep(T, length(ACATpvalues))
    is.keep[which(ACATpvalues == 1)] <- F
    ACAT_p <- ACAT(Pvals=ACATpvalues[is.keep], weights=wACAT[is.keep])
    pvalue <- c(SKAT_p = SKAT_p, Burden_p = Burden_p, ACAT_p = ACAT_p)
    
    # final weights
    weights <- list(wSKAT = wSKAT, wBurden = wBurden, wACAT = wACAT[is.keep])
    
  } else if (testtype == "SKAT") {
    
    # scores
    SKATscore <- SKATQ_fun(G, cvalue, wSKAT)
    scoreQ <- list(SKATscoreQ = SKATscore$Q)
    scores <- list(SKATscores = SKATscore$score)
    
    # Scovs
    if (boot == T) {
      Scovs <- Var_boot(G, para_cvalue = para_results, is.very.rare, mac.thresh, wBurden, wSKAT, type=testtype)
      SKAT_Scov <- Scovs$SKATSvar
    } else {

    }
    
    # pvalues
    Liu_result <- Liu(A = diag(ncol(SKAT_Scov)), X_Sigma = SKAT_Scov, X_mu = rep(0, ncol(SKAT_Scov)), X = unlist(SKATscore$score))
    SKAT_p <- unlist(Liu_result$pvalue)
    pvalue <- c(SKAT_p = SKAT_p)
    
    # final weights
    weights <- list(wSKAT = wSKAT)
    
  } else if (testtype == "Burden") {
    
    # scores
    Burdenscore <- BurdenQ_fun(G, cvalue, wBurden)
    scoreQ <- list(BurdenscoreQ = Burdenscore$Q)
    scores <- list(Burdenscores = Burdenscore$score)
    
    # Scovs
    if (boot == T) {
      Scovs <- Var_boot(G, para_cvalue = para_results, is.very.rare, mac.thresh, wBurden, wSKAT, type=testtype)
      Burden_Scov <- Scovs$BurdenSvar
    } else {

    }
    
    # pvalues
    Burden_p <- Burden(S = Burdenscore$score, sigma = Burden_Scov)
    pvalue <- c(Burden_p = Burden_p)
    
    # final weights
    weights <- list(wBurden = wBurden)
    
  } else if (testtype == "ACAT") {
    
    # scores
    ACATscore <- ACATsingleQ_fun(G, cvalue)
    scoreQ <- list(ACATscoreQs = ACATscore$Q)
    scores <- list(ACATscores = ACATscore$score)
    
    # Scovs
    if (boot == T) {
      Scovs <- Var_boot(G, para_cvalue = para_results, is.very.rare, mac.thresh, wBurden, wSKAT, type=testtype)
      ACAT_Scov <- Scovs$ACATSvar
    } else {

    }
    
    # pvalues
    ACATpvalues <- mapply(Burden, ACATscore$score, ACAT_Scov)
    is.keep <- rep(T, length(ACATpvalues))
    is.keep[which(ACATpvalues == 1)] <- F
    ACAT_p <- ACAT(Pvals=ACATpvalues[is.keep], weights=wACAT[is.keep])
    pvalue <- c(ACAT_p = ACAT_p)
    
    # final weights
    weights <- list(wACAT = wACAT[is.keep])
    
  }
  
  # final results
  results <- list(nobs = nrow(G), theta_est = theta, scoreQ = scoreQ, scores = scores, pvalue = pvalue, 
                  Scovs = Scovs, weights = weights, weights.beta = weights.beta, 
                  ACATpvalues = ACATpvalues, cvalue = cvalue)
  
  return(results)
}
