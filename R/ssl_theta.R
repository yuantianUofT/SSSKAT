#' Final parameter estimates under the null
#' 
#' @param Y Binary response variable
#' @param X Covariates
#' @param S Continuous surrogate
#' @param Z Covaraites with intercept column design matrix
#' @param id.t Row id of labeled data.
#' @param weights Weights of the individual in the sample, default treat data equally.
#' @param full_eval Full optimization iteration, default to be TRUE.
#' @param nit Number of iteration for optimization if full_NR_evaluation is FALSE.
#' @return Final estimates of the parameters
#' @export

ssl_theta <- function(Y, X, S, Z, id.t, weights = NULL,
                      full_eval = TRUE,
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