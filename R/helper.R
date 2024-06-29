# Updated: 2024-06-28

# -----------------------------------------------------------------------------
# Helper Functions
# -----------------------------------------------------------------------------

#' Define logit function
#' @param xx a numeric object
#' @export

logit <- function(xx) {
  log(xx/(1-xx))
}


#' Define expit function
#' @param xx a numeric object
#' @export

g.logit <- function(xx) {
  exp(xx)/(1+exp(xx))
}


#' Expit function under the null model
#' @param a intercept
#' @param b slope
#' @param x covariates
#' @return expit under the null model
#' @export

y.p_NULL <- function(a, b, x) {
  return(1/(1 + exp(-x%*%b-a)))
}


#' Full log-likelihood function under the null
#' @param Y Binary response variable
#' @param X Covariates
#' @param S Continuous surrogate
#' @param Z Covaraites with intercept column design matrix
#' @param id.t Row id of labeled data.
#' @param theta Distribution parameters.
#' @return log-likelihood function under the null
#' @export

NULL_log_like <- function(Y, X, S, Z, id.t, theta) {

  # parameters setup
  alpha <- theta[1:(ncol(X)+1)]
  m1 <- theta[(ncol(X)+2)]
  m0 <- theta[(ncol(X)+3)]
  log.s1 <- theta[(ncol(X)+4)]
  s1 <- exp(log.s1)
  log.s0 <- theta[(ncol(X)+5)]
  s0 <- exp(log.s0)

  # labeled likelihood
  if (length(id.t) == 0) {
    l_ob <- 0
  } else {
    l_ob_y11 <- Y[id.t] * log(g.logit(Z[id.t, ] %*% alpha))
    l_ob_y12 <- Y[id.t] * log(dnorm(S[id.t], mean=m1, sd=s1))
    l_ob_y21 <- (1-Y[id.t]) * log(1-g.logit(Z[id.t,] %*% alpha))
    l_ob_y22 <- (1-Y[id.t]) * log(dnorm(S[id.t], mean=m0, sd=s0))
    l_ob <- sum(l_ob_y11 + l_ob_y12 + l_ob_y21 + l_ob_y22)
  }

  # Unlabeled likelihood
  if (length(id.t) == 0) {
    l_unob_y1 <- (1 - g.logit(Z %*% alpha)) * dnorm(S, mean=m0, sd=s0)
    l_unob_y2 <- g.logit(Z %*% alpha) * dnorm(S, mean=m1, sd=s1)
  } else {
    l_unob_y1 <- (1 - g.logit(Z[-id.t,] %*% alpha)) * dnorm(S[-id.t], mean=m0, sd=s0)
    l_unob_y2 <- g.logit(Z[-id.t,] %*% alpha) * dnorm(S[-id.t], mean=m1, sd=s1)
  }
  l_unob <- sum(log(l_unob_y1 + l_unob_y2))

  # combine likelihoods
  l <- l_ob + l_unob

  # return likelihood
  return(l)
}


#' Negative log-likelihood function under the null
#' @param Y Binary response variable
#' @param X Covariates
#' @param S Continuous surrogate
#' @param Z Covaraites with intercept column design matrix
#' @param id.t Row id of labeled data.
#' @param theta Distribution parameters.
#' @return Negative log-likelihood function under the null
#' @export

NULL_nlog_like <- function(Y, X, S, Z, id.t, theta) {
  return(-NULL_log_like(Y, X, S, Z, id.t, theta))
}


#' Initial parameter estimates under the null
#' @param Y Binary response variable
#' @param S Continuous surrogate
#' @param X Covaraites with intercept column design matrix
#' @param id.t Row id of labeled data.
#' @param weights Weights of the individual in the sample, default treat data equally.
#' @return Initial estimates of the parameters
#' @export

sl_theta = function(Y, S, X, id.t, weights = NULL){

  # prevalence estimate
  if (length(id.t) != 0) {
    prev_est = mean(Y[id.t])
  }
  
  # weights setup
  if (is.null(weights)) {
    weights = rep(1, length(id.t))
  }
  
  # separate labeled_id based on case/control status
  id.1 = intersect(which(Y==1), id.t)
  id.0 = intersect(which(Y==0), id.t)

  # with some labeled data, under the null hypothesis, get supervised estimates
  if (length(id.t) == 0) {
    alpha=lm(S~X)$coef
  } else{
    alpha=glm(Y[id.t] ~ X[id.t, ],family=binomial(link = "logit"), weights = weights)$coef
  }

  # without labeled data, threshold S lower than 20% as labeled controls and greater than 80% as labeled cases
  if (length(id.t) == 0) {
    thre_value1 <- quantile(S, 0.8)
    thre_value2 <- quantile(S, 0.2)
    m1= mean(S[which(S>thre_value1)])
    s1= sd(S[which(S>thre_value1)])
    m0=mean(S[S<thre_value2])
    s0=sd(S[S<thre_value2])
  } else{
    m1= mean(S[id.1])
    s1= sd(S[id.1])
    m0=mean(S[id.0])
    s0=sd(S[id.0])
  }

  # initial estimates of the parameters
  init_est <- c(alpha = alpha, m1 = m1, logs1 = log(s1), m0 = m0, logs0 = log(s0))
  return(init_est)
}


#' Final parameter estimates under the null
#' @param Y Binary response variable
#' @param X Covariates
#' @param S Continuous surrogate
#' @param Z Covaraites with intercept column design matrix
#' @param id.t Row id of labeled data.
#' @param weights Weights of the individual in the sample, default treat data equally.
#' @param full_eval Full optimization iteration, default to be TRUE.
#' @param NULL_nlog_like Function to be optimized, predefined as negative log-likelihood function under the null
#' @param nit Number of iteration for optimization if full_NR_evaluation is FALSE.
#' @return Final estimates of the parameters
#' @export

ssl_theta <- function(Y, X, S, Z, id.t, weights = NULL,
                      full_eval = TRUE, NULL_nlog_like,
                      nit) {


  # initial estimates of the parameters
  init_sl = sl_theta(Y, S, X, id.t, weights)

  # final estimates maximizing likelihood
  if (full_eval){
    optim_temp <- tryCatch(optim(par=init_sl, fn=NULL_nlog_like, Y = Y, X = X, S = S, Z = Z, id.t=id.t,
                                 method="BFGS"), error=function(e) NA)
    if (length(optim_temp) == 1) {
      warning("Switched from BFGS to SANN in optim")
      optim_temp <- tryCatch(optim(par=init_sl, fn=NULL_nlog_like, Y = Y, X = X, S = S, Z = Z, id.t=id.t,
                                   method="SANN"), error=function(e) NA)
    }
    final_est <- optim_temp$par
    converge_steps <- optim_temp$counts[2]
  } else {
    optim_temp <- tryCatch(optim(par=init_sl, fn=NULL_nlog_like, Y = Y, X = X, S = S, Z = Z, id.t=id.t,
                                 method="BFGS", control = list(maxit=nit)), error=function(e) NA)
    final_est <- optim_temp$par
    converge_steps <- optim_temp$counts[2]
  }

  # final parameter estimates
  return(list(final_est = final_est, converge_steps = converge_steps, l_value = optim_temp$value))
}


#' c value calculation to be used in score function
#' @param Y Binary response variable
#' @param X Covariates
#' @param S Continuous surrogate
#' @param Z Covaraites with intercept column design matrix
#' @param id.t Row id of labeled data.
#' @param theta Distribution parameters
#' @return c value to be used in score function
#' @export

c_func <- function(Y, X, S, Z, id.t, theta) {
  # parameters setup
  alpha <- theta[1:(ncol(X)+1)]
  m1 <- theta[(ncol(X)+2)]
  m0 <- theta[(ncol(X)+3)]
  s1 <- exp(theta[(ncol(X)+4)])
  s0 <- exp(theta[(ncol(X)+5)])
  
  # c value in score function
  if (length(id.t) == 0) {
    eta2 = Z %*% alpha
    
    f0 <- dnorm(S, mean=m0, sd=s0)
    f1 <- dnorm(S, mean=m1, sd=s1)
    
    c2 = ((f1-f0)*g.logit(eta2)*(1-g.logit(eta2)))/(g.logit(eta2)*f1+(1-g.logit(eta2))*f0)
    
    c = rep(0, length(S))
    c = c2
  } else {
    eta1 = Z[id.t,] %*% alpha
    eta2 = Z[-id.t,] %*% alpha
    
    f0 <- dnorm(S[-id.t], mean=m0, sd=s0)
    f1 <- dnorm(S[-id.t], mean=m1, sd=s1)
    
    c1 = (Y[id.t] - g.logit(eta1))
    c2 = ((f1-f0)*g.logit(eta2)*(1-g.logit(eta2)))/(g.logit(eta2)*f1+(1-g.logit(eta2))*f0)
    
    c = rep(0, length(S))
    c[id.t] = c1
    c[-id.t] = c2
  }
  
  return(c)
}


#' SS_SKAT score and test statistics function
#' @param G Genotype data
#' @param cvalue Individual c values
#' @param wSKAT Vector of SNP SKAT weights.
#' @return SS_SKAT score and test statistics
#' @export

SKATQ_fun <- function(G, cvalue, wSKAT) {
  
  # SKAT weighted genotype data
  W <- diag(wSKAT)
  GW <- G %*% W

  # Obtain score and test statistics
  score = (t(GW) %*% cvalue)
  Q = sum(score^2)

  return(list(Q = Q, score = score))
}


#' SS_SKAT score function
#' @param GW SKAT weighted genotype data
#' @param cvalue Individual c values
#' @return SS_SKAT score
#' @export

SKATscore_fun <- function(GW, cvalue) {
  
  # Obtain score
  score <- (t(GW) %*% cvalue)
  score <- as.vector(score)
  
  return(score)
}


#' SS_burden score function
#' @param Gw Burden weighted genotype data
#' @param cvalue Individual c values
#' @return SS_burden score
#' @export

Burdenscore_fun <- function(Gw, cvalue) {
  
  # Obtain score
  score <- (t(Gw) %*% cvalue)
  score <- as.vector(score)
  
  return(score)
}


#' SS_burden score and test statistics function
#' @param G Genotype data
#' @param cvalue Individual c values
#' @param wBurden Vector of SNP burden weights.
#' @return SS_burden score and test statistics
#' @export

BurdenQ_fun <- function(G, cvalue, wBurden) {
  
  # Burden weighted genotype data
  Gw <- G %*% wBurden

  # Obtain score and test statistics
  score = (t(Gw) %*% cvalue)
  Q = as.vector(t(Gw) %*% cvalue %*% t(cvalue) %*% Gw)

  return(list(Q = Q, score = score))
}


#' SS_ACAT single SNP score and test statistics function
#' @param G Genotype data
#' @param cvalue Individual c values
#' @return SS_ACAT single SNP scores and test statistics
#' @export

ACATsingleQ_fun <- function(G, cvalue) {

  # Obtain score
  score = (t(G) %*% cvalue)
  Q = score^2

  return(list(Q = as.vector(Q), score = as.vector(score)))
}


#' SS_ACAT single SNP score function
#' @param G Genotype data
#' @param cvalue Individual c values
#' @return SS_ACAT single SNP score
#' @export

ACATsinglescore_fun <- function(G, cvalue) {
  
  score <- (t(G) %*% cvalue)
  score <- as.vector(score)
  
  return(score)
}


#' Rareness filtering function for ACAT
#' @param G Genotype data
#' @param mac.thresh A threshold of the minor allele count (MAC). The Burden test will be used to aggregate the SNPs with MAC less than this threshold.
#' @return If the SNPs are very rare. TRUE for rare.
#' @export

rare_func <- function(G, mac.thresh) {
  mac <- colSums(G)
  is.very.rare <- mac <= mac.thresh
  return(is.very.rare)
}


#' SS_ACAT SNP score and test statistics function
#' @param G Genotype data
#' @param cvalue Individual c values
#' @param mac.thresh A threshold of the minor allele count (MAC). The Burden test will be used to aggregate the SNPs with MAC less than this threshold.
#' @param is.very.rare If the SNPs are very rare. TRUE for rare.
#' @param wBurden Vector of SNP burden weights.
#' @return SS_ACAT scores and test statistics
#' @export

ACATQ_fun <- function(G, cvalue, mac.thresh, is.very.rare, wBurden) {
  
  if (sum(!is.very.rare) == 0) { # All SNPs are very rare, all SNPs are aggregated
    
    # score and test statistics
    BurdenscoreQs <- BurdenQ_fun(G = G, cvalue = cvalue, wBurden = wBurden)
    score <- BurdenscoreQs$score
    Q <- BurdenscoreQs$Q
    
  } else if (sum(is.very.rare) == 0) { # No SNPs are very rare, all SNPs are tested individually
    
    # score and test statistics
    ACATscoreQs <- ACATsingleQ_fun(G, cvalue)
    score <- ACATscoreQs$score
    Q <- ACATscoreQs$Q
    
  } else { # Some SNPs are very rare, the very rare SNPs are aggregated and the others are tested individually
    
    # separate the very rare SNPs and the others
    rareG <- G[, is.very.rare, drop = FALSE]
    denseG <- G[, (!is.very.rare), drop = FALSE]
    
    # score and test statistics
    wrare <- wBurden[is.very.rare]
    BurdenscoreQs <- BurdenQ_fun(rareG, cvalue, wBurden = wrare)
    ACATscoreQs <- ACATsingleQ_fun(denseG, cvalue)
    score <- c(BurdenscoreQs$score, ACATscoreQs$score)
    Q <- c(BurdenscoreQs$Q, ACATscoreQs$Q)
    
  }

  return(list(Q = as.vector(Q), score = as.vector(score)))
}


#' SS_ACAT SNP score function
#' @param Gw Burden weighted genotype data
#' @param G Genotype data
#' @param cvalue Individual c values
#' @param mac.thresh A threshold of the minor allele count (MAC). The Burden test will be used to aggregate the SNPs with MAC less than this threshold.
#' @param is.very.rare If the SNPs are very rare. TRUE for rare.
#' @return SS_ACAT scores
#' @export

ACATscore_fun <- function(Gw, G, cvalue, mac.thresh, is.very.rare) {
  
  if (sum(!is.very.rare) == 0) {
    Burdenscore <- Burdenscore_fun(Gw, cvalue = cvalue)
    score <- Burdenscore
  } else if (sum(is.very.rare) == 0) {
    ACATscore <- ACATsinglescore_fun(G, cvalue)
    score <- ACATscore
  } else {
    rareG <- G[, is.very.rare, drop = FALSE]
    denseG <- G[, (!is.very.rare), drop = FALSE]
    
    Burdenscore <- Burdenscore_fun(Gw[, is.very.rare], cvalue)
    ACATscore <- ACATsinglescore_fun(denseG, cvalue)
    
    score <- c(Burdenscore, ACATscore)
  }
  
  score <- as.vector(score)
  return(score)
}


#' Weights function for ACAT
#' @param G Genotype matrix.
#' @param is.very.rare If the SNPs are very rare. TRUE for rare.
#' @param wBurden Vector of SNP burden weights. When it is NULL, the beta weight with the “weights.beta” parameter is used.
#' @param weights.beta Weights beta parameters.
#' @return ACAT weights
#' @export

ACATW_func <- function(G, is.very.rare, wBurden, weights.beta) {

  mac <- colSums(G)
  n <- nrow(G)
  MAF <- mac/(2 * n)

  if (sum(!is.very.rare) == 0) { # All SNPs are very rare, all SNPs are aggregated
    
    # take mean of the weights for all SNPs
    if (is.null(wBurden)) {
      W <- mean(SKATBurdenW_func(G, weights.beta))
    } else {
      W <- mean(wBurden)
    }
    
  } else if (sum(is.very.rare) == 0) { # No SNPs are very rare, all SNPs are tested individually
    
    # calculate the weights for each SNPs, should be squared of that for Burden
    if (is.null(wBurden)) {
      W <- (dbeta(MAF, weights.beta[1], weights.beta[2])/dbeta(MAF, 0.5, 0.5))^2
    } else {
      W <- wBurden^2
    }
    
  } else { # Some SNPs are very rare, the very rare SNPs are aggregated and the others are tested individually
    
    # separate the very rare SNPs and the others
    # take mean of the weights for very rare SNPs
    # calculate the weights for not very rare SNPs, should be squared of that for Burden
    Grare <- G[, is.very.rare, drop = FALSE]
    if (is.null(wBurden)) {
      mafs <- c(mean(MAF[is.very.rare]), MAF[!is.very.rare])
      W_ACAT <- (dbeta(mafs,weights.beta[1],weights.beta[2])/dbeta(mafs,0.5,0.5))^2
      W <- W_ACAT
    } else {
      weights.sparse <- wBurden[is.very.rare]
      weights.dense <- wBurden[!is.very.rare]
      W <- c(mean(weights.sparse), weights.dense)
    }
  }
  
  return(W)
}


#' Weights function for Burden or SKAT
#' @param G Genotype matrix.
#' @param weights.beta Weights beta parameters.
#' @return Burden or SKAT weights
#' @export

SKATBurdenW_func <- function(G, weights.beta) {
  
  if (is.null(weights.beta)) {
    W <- rep(1, ncol(G))
  } else {
    MAF <- Matrix::colSums(G)/(2 * dim(G)[1])
    W <- dbeta(MAF, weights.beta[1], weights.beta[2])
  }
  return(W)
}


#' ACAT function to combine SNP pvalues
#' @param Pvals Single SNP pvalues to be combined
#' @param weights ACAT weights
#' @param is.check Check pvalues
#' @param is.small.thre Threshold for pvalue to be adjusted
#' @return ACAT combined pvalue
#' @export

ACAT <- function(Pvals, weights = NULL, is.check = T, is.small.thre = 1e-15) {

  is.keep <- rep(T, length(Pvals))
  is.keep[which(Pvals == 1)] <- F

  Pvals <- Pvals[is.keep]
  weights <- weights[is.keep]

  Pvals <- as.matrix(Pvals)
  if (is.check) {
    if (sum(is.na(Pvals)) > 0) {
      stop("Cannot have NAs in the p-values!")
    }
    if ((sum(Pvals < 0) + sum(Pvals > 1)) > 0) {
      stop("P-values must be between 0 and 1!")
    }
    is.zero <- (colSums(Pvals == 0) >= 1)
    is.one <- (colSums(Pvals == 1) >= 1)
    if (sum((is.zero + is.one) == 2) > 0) {
      stop("Cannot have both 0 and 1 p-values in the same column!")
    }
    if (sum(is.zero) > 0) {
      warning("There are p-values that are exactly 0!")
    }
    if (sum(is.one) > 0) {
      warning("There are p-values that are exactly 1!")
    }
  }
  if (is.null(weights)) {
    is.weights.null <- TRUE
  }
  else {
    is.weights.null <- FALSE
    weights <- as.matrix(weights)
    if (sum(dim(weights) != dim(Pvals)) > 0) {
      stop("The dimensions of weights and Pvals must be the same!")
    }
    else if (is.check & (sum(weights < 0) > 0)) {
      stop("All the weights must be nonnegative!")
    }
    else {
      w.sum <- colSums(weights)
      if (sum(w.sum <= 0) > 0) {
        stop("At least one weight should be positive in each column!")
      }
      else {
        for (j in 1:ncol(weights)) {
          weights[, j] <- weights[, j]/w.sum[j]
        }
      }
    }
  }
  is.small <- (Pvals < is.small.thre)
  if (is.weights.null) {
    Pvals[!is.small] <- tan((0.5 - Pvals[!is.small]) * pi)
    Pvals[is.small] <- 1/Pvals[is.small]/pi
    cct.stat <- colMeans(Pvals)
  }
  else {
    Pvals[!is.small] <- weights[!is.small] * tan((0.5 - Pvals[!is.small]) *
                                                   pi)
    Pvals[is.small] <- (weights[is.small]/Pvals[is.small])/pi
    cct.stat <- colSums(Pvals)
  }
  pval <- pcauchy(cct.stat, lower.tail = F)
  return(pval)
}


#' SS_SKAT score variance function
#' @param G Genotype data
#' @param para_cvalue Parametric boostrapped individual c values
#' @param wSKAT Vector of SNP SKAT weights.
#' @return SS_SKAT variance matrix
#' @export

SKATSVar_fun <- function(G, para_cvalue, wSKAT) {

  W <- diag(wSKAT)
  GW <- G %*% W

  # use each column of G times each column of para_cvalue
  para_scores <- as.matrix(t(apply(para_cvalue, 2, function(x) colSums(GW * x))))
  para_scores <- para_scores[complete.cases(para_scores), ]
  
  Scov <- cov(para_scores)
  
  return(Scov)
}


#' SS_Burden score variance function
#' @param G Genotype data
#' @param para_cvalue Parametric boostrapped individual c values
#' @param wBurden Vector of SNP burden weights.
#' @return SS_Burden variance matrix
#' @export

BurdenSVar_fun <- function(G, para_cvalue, wBurden) {
  
  Gw <- t(apply(G, 1, function(x) t(x * wBurden)))
  if (dim(G)[2] == 1) {
    Gw <- t(Gw)
  }
  
  para_scores <- apply(para_cvalue, 2, function(x) sum(Gw * x))
  para_scores <- para_scores[complete.cases(para_scores)]
  
  varb <- var(para_scores)

  return(varb)
}


#' SS_ACAT single SNP score function variances
#' @param G Genotype data
#' @param para_cvalue Parametric boostrapped individual c values
#' @return SS_ACAT single SNP variances
#' @export

ACATSsingleVar_fun <- function(G, para_cvalue) {

  para_scores <- apply(para_cvalue, 2, function(x) sum(G * x))
  para_scores <- para_scores[complete.cases(para_scores)]
  
  var_single <- var(para_scores)

  return(var_single)
}


#' SS_ACAT SNP score variance function
#' @param G Genotype data
#' @param para_cvalue Parametric boostrapped individual c values
#' @param is.very.rare If the SNPs are very rare.
#' @param mac.thresh A threshold of the minor allele count (MAC). The Burden test will be used to aggregate the SNPs with MAC less than this threshold.
#' @param wBurden Vector of SNP burden weights. 
#' @return SS_ACAT variance
#' @export

ACATSVar_fun <- function(G, para_cvalue, is.very.rare, mac.thresh, wBurden) {

  if (sum(!is.very.rare) == 0) {
    Burden_var <- BurdenSVar_fun(G, para_cvalue, wBurden)
    S_vars <- Burden_var
  } else if (sum(is.very.rare) == 0) {
    ACAT_vars <- ACATSsingleVar_fun(G, para_cvalue)
    S_vars <- ACAT_vars
  } else {
    rareG <- G[, is.very.rare, drop = FALSE]
    denseG <- G[, (!is.very.rare), drop = FALSE]
    wrare <- wBurden[is.very.rare]

    Burden_var <- BurdenSVar_fun(rareG, para_cvalue, wrare)
    ACAT_vars <- ACATSsingleVar_fun(denseG, para_cvalue)

    S_vars <- c(Burden_var, ACAT_vars)
  }

  return(S_vars)
}


#' Single run of parametric bootstrap variance function for ACAT, Burden and SKAT
#' @param G Genotype data
#' @param para_cvalue Parametric boostrapped individual c values
#' @param is.very.rare If the SNPs are very rare.
#' @param mac.thresh A threshold of the minor allele count (MAC). The Burden test will be used to aggregate the SNPs with MAC less than this threshold.
#' @param wBurden Vector of SNP burden weights. 
#' @param wSKAT Vector of SNP SKAT weights.
#' @param type Type of SS tests
#' @return Bootstrapped SS score variances
#' @export

Var_boot <- function(G, para_cvalue, is.very.rare, mac.thresh, wBurden, wSKAT, type) {
  
  # prepare dat
  n <- nrow(G)
  W <- diag(wSKAT)
  GW <- G %*% W
  Gw <- G %*% wBurden
  
  if (type == "ACAT") {
    ACATSvar <- ACATSVar_fun(G, para_cvalue, is.very.rare, mac.thresh, wBurden)
    Svar <- list(ACATSvar = ACATSvar)
    
  } else if (type == "SKAT") {
    SKATSvar <- SKATSVar_fun(G, para_cvalue, wSKAT)
    Svar <- list(SKATSvar = SKATSvar)
    
  } else if (type == "Burden") {
    BurdenSvar <- BurdenSVar_fun(G, para_cvalue, wBurden)
    Svar <- list(BurdenSvar = BurdenSvar)
    
  } else {
    ACATSvar <- ACATSVar_fun(G, para_cvalue, is.very.rare, mac.thresh, wBurden)
    SKATSvar <- SKATSVar_fun(G, para_cvalue, wSKAT)
    BurdenSvar <- BurdenSVar_fun(G, para_cvalue, wBurden)
    Svar <- list(ACATSvar = ACATSvar, SKATSvar = SKATSvar, BurdenSvar = BurdenSvar)
  }
  
  return(Svar)
}


#' Liu's non-central chi-square approximation with tail adjust
#' @param A The weight matrix in test statistics Q=S'AS.
#' @param X_Sigma Covariance matrix of the score function S.
#' @param X_mu Mean of the score function.
#' @param X Score vector S.
#' @return pvalue, non-central chi-square df and nc.
#' @export

Liu_adj <- function(A, X_Sigma, X_mu, X) {
  Q_t <- t(X) %*% A %*% X

  c1 <- tr((A %*% X_Sigma)) + t(X_mu) %*% A %*% X_mu
  c2 <- tr((A %*% X_Sigma) %^% 2) + 2 * t(X_mu) %*% A %*% X_Sigma %*% A %*% X_mu
  c3 <- tr((A %*% X_Sigma) %^% 3) + 3 * t(X_mu) %*% ((A %*% X_Sigma) %^% 2) %*% A %*% X_mu
  c4 <- tr((A %*% X_Sigma) %^% 4) + 4 * t(X_mu) %*% ((A %*% X_Sigma) %^% 3) %*% A %*% X_mu
  s1 <- c3/(c2^(3/2))
  s2 <- c4/(c2^2)

  Q_mu <- c1
  Q_sigma <- sqrt(2*c2)
  Q_skew <- sqrt(8)*s1
  Q_kurt <- 12*s2

  if (s1^2 > s2) {
    a = 1/(s1-sqrt(s1^2 - s2))
    delta <- s1*a^3 - a^2
    l <- a^2 - 2*delta
  } else {
    delta <- 0
    l <- 1/(s2)
  }

  Q_t_new <- (Q_t - Q_mu)/Q_sigma
  Chi_mu <- l + delta
  Sigma_chi <- sqrt(2) * sqrt(l + 2*delta)

  pvalue <- pchisq((Q_t_new*Sigma_chi + Chi_mu), df = l, ncp = delta, lower.tail = F)

  result <- list(pvalue = pvalue, df = l, ncp = delta)
  return(result)
}


#' Liu's non-central chi-square approximation for SS_SKAT pvalue calculation
#' @param A The weight matrix in test statistics Q=S'AS.
#' @param X_Sigma Covariance matrix of the score function S.
#' @param X_mu Mean of the score function.
#' @param X Score vector S.
#' @return pvalue, non-central chi-square df and nc.
#' @export

Liu <- function(A, X_Sigma, X_mu, X) {
  Q_t <- t(X) %*% A %*% X

  c1 <- tr((A %*% X_Sigma)) + t(X_mu) %*% A %*% X_mu
  c2 <- tr((A %*% X_Sigma) %^% 2) + 2 * t(X_mu) %*% A %*% X_Sigma %*% A %*% X_mu
  c3 <- tr((A %*% X_Sigma) %^% 3) + 3 * t(X_mu) %*% ((A %*% X_Sigma) %^% 2) %*% A %*% X_mu
  c4 <- tr((A %*% X_Sigma) %^% 4) + 4 * t(X_mu) %*% ((A %*% X_Sigma) %^% 3) %*% A %*% X_mu
  s1 <- c3/(c2^(3/2))
  s2 <- c4/(c2^2)

  Q_mu <- c1
  Q_sigma <- sqrt(2*c2)
  Q_skew <- sqrt(8)*s1
  Q_kurt <- 12*s2

  if (s1^2 > s2) {
    a = 1/(s1-sqrt(s1^2 - s2))
    delta <- s1*a^3 - a^2
    l <- a^2 - 2*delta
  } else {
    delta <- 0
    l <- 1/(s1^2)
  }

  Q_t_new <- (Q_t - Q_mu)/Q_sigma
  Chi_mu <- l + delta
  Sigma_chi <- sqrt(2) * sqrt(l + 2*delta)

  pvalue <- pchisq((Q_t_new*Sigma_chi + Chi_mu), df = l, ncp = delta, lower.tail = F)

  result <- list(pvalue = pvalue, df = l, ncp = delta)
  return(result)
}


#' Chi-square for SS_Burden pvalue calculation
#' @param S S score of burden, Q=S'S.
#' @param sigma Variance of S.
#' @return pvalue.
#' @export

Burden <- function(S, sigma) {
  V <- S/sqrt(sigma)
  Q <- V^2
  pval <- 1 - pchisq(Q, df = 1)
  return(pval)
}


#' MinP calculation
#' @param ps Pvalue to be combined
#' @return pvalue.
#' @export

minP <- function(ps) {
  p1 <- min(ps)
  pvalue <- length(ps) * p1
  return(pvalue)
}


#' SimesP calculation
#' @param ps Pvalue to be combined
#' @return pvalue.
#' @export

SimesP <- function(ps) {
  sortps <- sort(ps)
  newps <- sortps*length(ps)/(1:length(ps))
  pvalue <- min(newps)
  return(pvalue)
}

