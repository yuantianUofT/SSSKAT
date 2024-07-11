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
#' @return Parametric boostrapped individual c values
#' @export

para_func <- function(nn, theta, Y, X, S, Z, id.t) {
  
  # predict Y
  est_alpha <- theta[-c((length(theta)-3):length(theta))]
  pred_Y <- g.logit(Z %*% est_alpha)
  
  # parametric bootstrap
  m1 <- theta[(length(theta)-3)]
  sd1 <- exp(theta[(length(theta)-2)])
  m0 <- theta[(length(theta)-1)]
  sd0 <- exp(theta[length(theta)])
  para_Y <- rbinom(n = nn, size = 1, prob = pred_Y)
  para_S <- rep(NA, nn)
  para_S[para_Y == 1] <- rnorm(n = sum(para_Y), mean = m1, sd = sd1)
  para_S[para_Y == 0] <- rnorm(n = sum(para_Y == 0), mean = m0, sd = sd0)
  para_id.t <- sample(1:nn, size = length(id.t), replace = FALSE)
  
  # parametric parameter estimation 
  para_est_theta <- suppressWarnings(ssl_theta(Y = para_Y, X = X, S = para_S, Z = Z, 
                                               id.t = para_id.t, weights = NULL, full_eval = TRUE, 
                                               NULL_nlog_like, nit)$final_est)
  # parametric c-value estimation
  para_cvalue_est <- c_func(Y = para_Y, X = X, S = para_S, Z = Z, id.t = para_id.t, theta=para_est_theta)
  
  return(para_cvalue_est)
}