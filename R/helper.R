# Updated: 2023-09-19

# -----------------------------------------------------------------------------
# Helper Functions
# -----------------------------------------------------------------------------

#' Define logit function
#' @param xx a numeric object
#' @noRd

logit <- function(xx) {
  log(xx/(1-xx))
}

#' Define expit function
#' @param xx a numeric object
#' @noRd

g.logit <- function(xx) {
  exp(xx)/(1+exp(xx))
}

#' Full log-likelihood function
#' @param dat Basis matrix of data set to be analyzed.
#' @param id.t Row id of labeled data in the dat dataset.
#' @param theta Distribution parameters.
#' @param beta Effects of genotype
#' @return Likelihood of labeled and unlabeled parts
#' @noRd

NULL_log_like <- function(dat, id.t, theta, beta=NULL) {
  # data preparation
  Y <- dat$Y
  X <- dplyr::select(dat, starts_with('X'))
  G <- dplyr::select(dat, starts_with('G'))
  S <- dat$S
  Z <- data.matrix(cbind(1, X, G))

  # parameters setup
  alpha <- theta[1:(ncol(X)+1)]
  m1 <- theta[(ncol(X)+2)]
  m0 <- theta[(ncol(X)+3)]
  log.s1 <- theta[(ncol(X)+4)]
  s1 <- exp(log.s1)
  log.s0 <- theta[(ncol(X)+5)]
  s0 <- exp(log.s0)
  beta <- rep(0, ncol(G))

  # labeled likelihood
  if (length(id.t) == 0) {
    l_ob <- 0
  } else {
    l_ob_y11 <- Y[id.t] * log(g.logit(Z[id.t,] %*% c(alpha, beta)))
    l_ob_y12 <- Y[id.t] * log(dnorm(S[id.t], mean=m1, sd=s1))
    l_ob_y21 <- (1-Y[id.t]) * log(1-g.logit(Z[id.t,] %*% c(alpha, beta)))
    l_ob_y22 <- (1-Y[id.t]) * log(dnorm(S[id.t], mean=m0, sd=s0))
    l_ob <- sum(l_ob_y11 + l_ob_y12 + l_ob_y21 + l_ob_y22)
  }

  # Unlabeled likelihood
  if (length(id.t) == 0) {
    l_unob_y1 <- (1 - g.logit(Z %*% c(alpha, beta))) * dnorm(S, mean=m0, sd=s0)
    l_unob_y2 <- g.logit(Z %*% c(alpha, beta)) * dnorm(S, mean=m1, sd=s1)
  } else {
    l_unob_y1 <- (1 - g.logit(Z[-id.t,] %*% c(alpha, beta))) * dnorm(S[-id.t], mean=m0, sd=s0)
    l_unob_y2 <- g.logit(Z[-id.t,] %*% c(alpha, beta)) * dnorm(S[-id.t], mean=m1, sd=s1)
  }
  l_unob <- sum(log(l_unob_y1 + l_unob_y2))

  # combine likelihoods
  l <- l_ob + l_unob

  # return likelihood
  return(l)
}

#' Negative log-likelihood function
#' @param dat Basis matrix of data set to be analyzed.
#' @param id.t Row id of labeled data in the dat dataset.
#' @param theta Distribution parameters.
#' @param beta Effects of genotype
#' @return Negative likelihood of labeled and unlabeled parts
#' @noRd

NULL_nlog_like <- function(dat, id.t, theta, beta = NULL) {
  return(-NULL_log_like(dat, id.t, theta, beta = NULL))
}

#' Initial parameter estimates
#' @param dat Basis matrix of data set to be analyzed.
#' @param id.t Row id of labeled data in the dat dataset.
#' @param weights Weights of the data, default treat data equally.
#' @param prev_est Estimate of prevalence.
#' @return Initial estimates of the parameters
#' @noRd

sl_theta = function(dat, id.t, weights = NULL, prev_est){
  # data preparation
  Y = dat$Y
  S = dat$S
  X = data.matrix((dplyr::select(dat, starts_with('X'))))

  if(is.null(weights)){
    weights = rep(1, length(id.t))
  }

  id.1 = intersect(which(Y==1), id.t)
  id.0 = intersect(which(Y==0), id.t)

  # under the null hypothesis, get supervised estimates
  if (length(id.t) == 0) {
    alpha=lm(S~X)$coef
  } else{
    alpha=glm(Y[id.t]~X[id.t,],family=binomial(link = "logit"), weights = weights)$coef
  }

  if (length(id.t) == 0) {
    thre_value1 <- quantile(S, 1-prev_est)
    thre_value2 <- quantile(S, prev_est)
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

#' Final parameter estimates
#' @param dat Basis matrix of data set to be analyzed.
#' @param id.t Row id of labeled data in the dat dataset.
#' @param weights Weights of the data, default treat data equally.
#' @param full_eval Full optimization iteration, default to be TRUE.
#' @param NULL_nlog_like Function to be optimized, predefined, no need to change
#' @param nit Number of iteration for optimization if full_NR_evaluation is FALSE.
#' @param prev_est Estimate of prevalence.
#' @return Final estimates of the parameters
#' @noRd

ssl_theta <- function(dat, id.t, weights = NULL,
                      full_eval = TRUE, NULL_nlog_like,
                      nit, prev_est) {
  # data preparation
  Y = dat$Y
  X <- dplyr::select(dat, starts_with('X'))

  # initial estimates of the parameters
  init_sl = sl_theta(dat, id.t, weights, prev_est)

  # final estimates maximizing likelihood
  if (full_eval){
    optim_temp <- tryCatch(optim(par=init_sl, fn=NULL_nlog_like, dat=dat, id.t=id.t,
                                 method="BFGS"), error=function(e) NA)
    if (length(optim_temp) == 1) {
      warning("Switched from BFGS to SANN in optim")
      optim_temp <- tryCatch(optim(par=init_sl, fn=NULL_nlog_like, dat=dat, id.t=id.t,
                                   method="SANN"), error=function(e) NA)
    }
    final_est <- optim_temp$par
    converge_steps <- optim_temp$counts[2]
  } else {
    optim_temp <- tryCatch(optim(par=init_sl, fn=NULL_nlog_like, dat=dat, id.t=id.t,
                                 method="BFGS", control = list(maxit=nit)), error=function(e) NA)
    final_est <- optim_temp$par
    converge_steps <- optim_temp$counts[2]
  }

  # final parameter estimates
  return(list(final_est = final_est, converge_steps = converge_steps, l_value = optim_temp$value))
}

#' Expit function under the null model
#' @param a intercept
#' @param b slope
#' @param x covariates
#' @return expit under the null model
#' @noRd

y.p_NULL <- function(a, b, x) {
  return(1/(1 + exp(-x%*%b-a)))
}

#' SS_SKAT score function
#' @param dat Basis matrix of data set to be analyzed.
#' @param id.t Row id of labeled data in the dat dataset.
#' @param theta Distribution parameters.
#' @return SS_SKAT score
#' @noRd

score_fun <- function(dat, id.t, theta) {
  # data preparation
  Y <- dat$Y
  X <- dplyr::select(dat, starts_with('X'))
  G <- dplyr::select(dat, starts_with('G'))
  G <- as.matrix(G)
  S <- dat$S
  Z <- data.matrix(cbind(1, X))

  # parameters setup
  alpha <- theta[1:(ncol(X)+1)]
  m1 <- theta[(ncol(X)+2)]
  m0 <- theta[(ncol(X)+3)]
  s1 <- exp(theta[(ncol(X)+4)])
  s0 <- exp(theta[(ncol(X)+5)])

  # score function
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

  # Obtain score
  score = (t(G) %*% c)
  Q = sum(diag(t(G) %*% c %*% t(c) %*% G))

  return(list(Q = Q, score = score))
}

#' Parametric bootstrap to obtain SS_SKAT score under the null
#' @param dat Basis matrix of data set to be analyzed.
#' @param id.t Row id of labeled data in the dat dataset.
#' @param full_eval Full optimization iteration, default to be TRUE.
#' @param theta Distribution parameters.
#' @return Parametric bootstrap SS_SKAT score
#' @noRd

boot_Q <- function(dat, id.t,
                   full_eval = F, theta = NULL) {
  # data preparation
  X <- dplyr::select(dat, starts_with('X'))
  X <- as.matrix(X)
  G <- dplyr::select(dat, starts_with('G'))

  # parameters setup
  a <- theta[1]
  b <- theta[2:(ncol(X)+1)]
  m1 <- theta[(ncol(X)+2)]
  m0 <- theta[(ncol(X)+3)]
  log.s1 <- theta[(ncol(X)+4)]
  s1 <- exp(log.s1)
  log.s0 <- theta[(ncol(X)+5)]
  s0 <- exp(log.s0)

  # parametric bootstrap
  y.p_value <- sapply(1:nrow(X), function(index) y.p_NULL(a, b, X[index, ]))
  newY <- sapply(1:nrow(X), function(index) rbinom(1, 1, y.p_value[index]))
  newS = rep(NA, nrow(G))
  newS[newY == 1] <- rnorm(n = sum(newY), mean = m1, sd = s1)
  newS[newY == 0] <- rnorm(n = sum(1-newY), mean = m0, sd = s0)

  newdat <- cbind(newY, X, G, newS)
  colnames(newdat)[1] <- "Y"
  colnames(newdat)[ncol(newdat)] <- "S"
  newmy_score <- score_fun(dat=newdat, id.t=id.t, theta=theta)$Q

  # bootstrap SS_SKAT score under the null
  return(newmy_score)
}


