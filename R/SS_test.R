# Updated: 2023-09-19

#' Semi-supervised association test function
#'
#' @param dat Basis matrix of data set to be analyzed.
#' @param id.t Row id of labeled data in the dat dataset.
#' @param obsweights Weights of the sample, default treat samples equally.
#' @param weights.beta Weights of the SNPs, default NULL treat the SNPs equally, else take parameters of beta distribution.
#' @param full_NR_evaluation Full optimization iteration, default to be TRUE.
#' @param nit Number of iteration for optimization if full_NR_evaluation is FALSE.
#' @param NULL_nlog_like Function to be optimized, predefined, no need to change
#' @param prev_est Estimate of prevalence.
#' @param testtype Type of test, default is "all", can also be "SKAT", "Burden", "ACAT"
#' @param nboot Number of bootstrap of pvalue calculate.
#' @return Parameter estimates, SS SKAT score and SS SKAT pvalue.
#' @examples
#' analysis_data <- example_data$my_dat;
#' labeled_data_id <- example_data$id.t;
#' SS_SKAT_results <- SS_SKAT(dat = analysis_data, id.t = labeled_data_id, weights = NULL, full_NR_evaluation = TRUE, NULL_nlog_like = NULL_nlog_like, prev_est = mean(analysis_data$Y[labeled_data_id]), nboot = 50);
#' @export


SS_test <- function(dat, id.t, obsweights = NULL, weights.beta = NULL,
                    full_NR_evaluation = TRUE, nit = NULL,
                    NULL_nlog_like, prev_est, testtype = "all",
                    nboot) {

  # data preparation
  Y <- dat$Y
  X <- dplyr::select(dat, starts_with('X'))
  G <- dplyr::select(dat, starts_with('G'))
  G <- as.matrix(G)
  S <- dat$S
  Z <- data.matrix(cbind(1, X))
  if (is.null(weights.beta)) {
    ACATweights <- rep(1, ncol(G))
    SKATweights <- rep(1, ncol(G))
    Burdenweights <- rep(1, ncol(G))
  } else {
    ACATweights <- ACATW_func(G, weights.beta)
    SKATweights <- SKATBurdenW_func(G, weights.beta)
    Burdenweights <- SKATBurdenW_func(G, weights.beta)
  }

  # parameter estimates
  theta_est <- ssl_theta(dat=dat, id.t=id.t, weights = obsweights, full_eval = full_NR_evaluation, nit = nit,
                         NULL_nlog_like = NULL_nlog_like, prev_est = prev_est)$final_est

  # bootstrap
  boot_scores <- lapply(1:nboot, function(x) boot_Q(dat = dat, id.t = id.t, theta = theta_est, testtype = testtype,
                                                    ACATweights = ACATweights, SKATweights = SKATweights, Burdenweights = Burdenweights,
                                                    x))

  # SS SKAT score
  if (testtype == "all") {
    SKATscoreQ <- SKATQ_fun(dat=dat, id.t=id.t, theta=theta_est, w = SKATweights)$Q
    BurdenscoreQ <- BurdenQ_fun(dat=dat, id.t=id.t, theta=theta_est, w = Burdenweights)$Q
    ACATscoreQs <- ACATQ_fun(dat=dat, id.t=id.t, theta=theta_est)$Q
    scoreQ <- list(SKATscoreQ = SKATscoreQ, BurdenscoreQ = BurdenscoreQ, ACATscoreQs = ACATscoreQs)
  } else if (testtype == "SKAT") {
    SKATscoreQ <- SKATQ_fun(dat=dat, id.t=id.t, theta=theta_est, w = SKATweights)$Q
    scoreQ <- list(SKATscoreQ = SKATscoreQ)
  } else if (testtype == "Burden") {
    BurdenscoreQ <- BurdenQ_fun(dat=dat, id.t=id.t, theta=theta_est, w = Burdenweights)$Q
    scoreQ <- list(BurdenscoreQ = BurdenscoreQ)
  } else if (testtype == "ACAT") {
    ACATscoreQs <- ACATQ_fun(dat=dat, id.t=id.t, theta=theta_est)$Q
    scoreQ <- list(ACATscoreQs = ACATscoreQs)
  }


  # pvalue estimate
  if (testtype == "all") {
    SKATboot_scores <- unlist(lapply(1:length(boot_scores),function(x) boot_scores[[x]]$newSKATscoreQ))
    Burdenboot_scores <- unlist(lapply(1:length(boot_scores),function(x) boot_scores[[x]]$newBurdenscoreQ))
    ACATboot_scores <- matrix(unlist(lapply(1:length(boot_scores),function(x) boot_scores[[x]]$newACATscoreQs)), nrow = length(boot_scores), byrow = T)

    SKATpvalue <- length(which(SKATboot_scores > scoreQ$SKATscoreQ))/length(SKATboot_scores)
    Burdenpvalue <- length(which(Burdenboot_scores > scoreQ$BurdenscoreQ))/length(Burdenboot_scores)
    ACATpvalues <- 1 - sapply(1:length(scoreQ$ACATscoreQs), function(j) mean(ACATboot_scores[, j] <= scoreQ$ACATscoreQs[j]))
    ACATpvalue <- ACAT(ACATpvalues, ACATweights)

    pvalue <- c(SKATpvalue = SKATpvalue, Burdenpvalue = Burdenpvalue, ACATpvalue = ACATpvalue)
  } else if (testtype == "SKAT") {
    SKATboot_scores <- unlist(lapply(1:length(boot_scores),function(x) boot_scores[[x]]$newSKATscoreQ))
    SKATpvalue <- length(which(SKATboot_scores > scoreQ$SKATscoreQ))/length(SKATboot_scores)
    pvalue <- c(SKATpvalue = SKATpvalue)
  } else if (testtype == "Burden") {
    Burdenboot_scores <- unlist(lapply(1:length(boot_scores),function(x) boot_scores[[x]]$newBurdenscoreQ))
    Burdenpvalue <- length(which(Burdenboot_scores > scoreQ$BurdenscoreQ))/length(Burdenboot_scores)
    pvalue <- c(Burdenpvalue = Burdenpvalue)
  } else if (testtype == "ACAT") {
    ACATboot_scores <- matrix(unlist(lapply(1:length(boot_scores),function(x) boot_scores[[x]]$newACATscoreQs)), nrow = length(boot_scores), byrow = T)
    ACATpvalues <- 1 - sapply(1:length(scoreQ$ACATscoreQs), function(j) mean(ACATboot_scores[, j] <= scoreQ$ACATscoreQs[j]))
    ACATpvalue <- ACAT(ACATpvalues, ACATweights)
    pvalue <- c(ACATpvalue = ACATpvalue)
  }


  results <- list(theta_est = theta_est, scoreQ = scoreQ, pvalue = pvalue)
  return(results)
}
