# Updated: 2023-11-19

#' test directly with S for both labeled and unlabeled data
#'
#' @param dat Basis matrix of data set to be analyzed.
#' @param weights.beta Weights of the SNPs, default NULL treat the SNPs equally, else take parameters of beta distribution.
#' @return Naive pvalues.
#' analysis_data <- example_data$my_dat;
#' naive_SKAT_results <- naive_test(dat = analysis_data, weights.beta = NULL);
#' @export

naive_test <- function(dat, weights.beta = NULL) {
  # data preparation
  S <- dat$S
  X <- as.matrix(dat %>% dplyr::select(starts_with('X')))
  G <- as.matrix(dat %>% dplyr::select(starts_with('G')))

  if (is.null(weights.beta)) {
    ACATweights <- rep(1, ncol(G))
    SKATweights <- rep(1, ncol(G))
    Burdenweights <- rep(1, ncol(G))
  } else {
    ACATweights <- ACATW_func(G, weights.beta)
    SKATweights <- SKATBurdenW_func(G, weights.beta)
    Burdenweights <- SKATBurdenW_func(G, weights.beta)
  }

  # SKAT and Burden NULL model
  obj.b <- SKAT_Null_Model(S ~ X, out_type="C", Adjustment = F)
  # SKAT results
  SKAT_result <- SKAT(G, obj.b, method="davies", weights=SKATweights)
  # Burden results
  Burden_result <- SKAT(G, obj.b, method="davies", r.corr = 1, weights=Burdenweights)

  # ACAT NULL model
  G <- Matrix::Matrix(G, sparse = TRUE)
  obj <- ACAT::NULL_Model(S, X)
  ACAT_result <- ACAT::ACAT_V(G, obj)

  result <- list(SKAT_p = SKAT_result$p.value, Burden_p = Burden_result$p.value, ACAT_p = ACAT_result)

  return(result)
}
