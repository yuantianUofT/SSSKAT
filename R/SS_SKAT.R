# Updated: 2023-09-19

#' Semi-supervised SKAT function
#'
#' @param dat Basis matrix of data set to be analyzed.
#' @param id.t Row id of labeled data in the dat dataset.
#' @param weights Weights of the data, default treat data equally.
#' @param full_NR_evaluation Full optimization iteration, default to be TRUE.
#' @param nit Number of iteration for optimization if full_NR_evaluation is FALSE.
#' @param NULL_nlog_like Function to be optimized, predefined, no need to change
#' @param prev_est Estimate of prevalence.
#' @param nboot Number of bootstrap of pvalue calculate.
#' @return Parameter estimates, SS SKAT score and SS SKAT pvalue.
#' @examples
#' analysis_data <- example_data$my_dat;
#' labeled_data_id <- example_data$id.t;
#' SS_SKAT_results <- SS_SKAT(dat = analysis_data, id.t = labeled_data_id, weights = NULL, full_NR_evaluation = TRUE, NULL_nlog_like = NULL_nlog_like, prev_est = mean(analysis_data$Y[labeled_data_id]), nboot = 50);
#' @export


SS_SKAT <- function(dat, id.t, weights = NULL,
                    full_NR_evaluation = TRUE, nit = NULL,
                    NULL_nlog_like, prev_est,
                    nboot) {
  # parameter estimates
  theta_est <- ssl_theta(dat=dat, id.t=id.t, weights = weights, full_eval = full_NR_evaluation, nit = nit,
                         NULL_nlog_like = NULL_nlog_like, prev_est = prev_est)$final_est

  # bootstrap
  boot_scores <- sapply(1:nboot, function(x) boot_Q(dat = dat, id.t = id.t, theta = theta_est, x))

  # SS SKAT score
  scoreQ <- score_fun(dat=dat, id.t=id.t, theta=theta_est)$Q

  # pvalue estimate
  pvalue <- length(which(boot_scores > scoreQ))/length(boot_scores)

  results <- list(theta_est = theta_est, scoreQ = scoreQ, pvalue = pvalue)
  return(results)
}
