# Updated: 2024-06-28

#' Parametric boostrap function for individual c values
#' 
#' @param nn Number of samples
#' @param theta Distribution parameters
#' @param Y Binary response variable
#' @param X Covariates
#' @param S Continuous surrogate
#' @param Z Covaraites with intercept column design matrix
#' @param id.t Row id of labeled data.
#' @param distri Distribution of S|Y, either "normal" or "beta"
#' @return Parametric boostrapped individual c values
#' @export

para_func <- function(nn, theta, Y, X, S, Z, id.t, distri) {
  
  # predict Y
  est_alpha <- theta[-c((length(theta)-3):length(theta))]
  pred_Y <- g.logit(Z %*% est_alpha)
  
  # parametric bootstrap
  if (distri == "normal") {
    m1 <- theta[(length(theta)-3)]
    s1 <- exp(theta[(length(theta)-2)])
    m0 <- theta[(length(theta)-1)]
    s0 <- exp(theta[length(theta)])
  } else {
    m1 <- exp(theta[(length(theta)-3)])
    s1 <- exp(theta[(length(theta)-2)])
    m0 <- exp(theta[(length(theta)-1)])
    s0 <- exp(theta[length(theta)])
  }
  
  para_Y <- rbinom(n = nn, size = 1, prob = pred_Y)
  para_S <- rep(NA, nn)
  if (distri == "normal") {
    para_S[para_Y == 1] <- rnorm(n = sum(para_Y), mean = m1, sd = s1)
    para_S[para_Y == 0] <- rnorm(n = sum(para_Y == 0), mean = m0, sd = s0)
  } else if (distri == "beta") {
    para_S[para_Y == 1] <- rbeta(n = sum(para_Y), shape1 = m1, shape2 = s1)
    para_S[para_Y == 0] <- rbeta(n = sum(para_Y == 0), shape1 = m0, shape2 = s0)
  }
  smalle <- 1e-35
  if (length(which(para_S == 1)) > 0) {
    para_S[which(para_S == 1)] <- max(para_S[para_S != 1])
  } 
  if (length(which(para_S == 0)) > 0) {
    para_S[which(para_S == 0)] <- smalle
  }
  para_S <- unlist(para_S)
  para_id.t <- sample(1:nn, size = length(id.t), replace = FALSE)
  
  # parametric parameter estimation
  para_est <- ssl_theta(Y = para_Y, X = X, S = para_S, Z = Z, 
                        id.t = para_id.t, weights = NULL, full_eval = TRUE, 
                        NULL_nlog_like, nit, distri = distri)
  para_est_theta <- para_est$final_est
  # parametric c-value estimation
  para_cvalue_est <- c_func(Y = para_Y, X = X, S = para_S, Z = Z, id.t = para_id.t, theta=para_est_theta, distri = distri)
  
  outlist <- list(para_est = para_est, para_cvalue_est = para_cvalue_est)
  return(outlist)
}