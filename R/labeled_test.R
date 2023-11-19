# Updated: 2023-09-19

#' Standard test function using the labeled data only
#'
#' @param dat Basis matrix of data set to be analyzed.
#' @param id.t Row id of labeled data in the dat dataset.
#' @param weights.beta Weights of the SNPs, default NULL treat the SNPs equally, else take parameters of beta distribution.
#' @return Parameter estimates, SS SKAT score and SS SKAT pvalue.
#' @examples
#' analysis_data <- example_data$my_dat;
#' labeled_data_id <- example_data$id.t;
#' labeled_SKAT_result <- labeled_SKAT(dat = analysis_data, id.t = labeled_data_id);
#' @export

labeled_test <- function(dat, id.t, weights.beta = NULL) {
  # prepare data
  y.b <- dat$Y[id.t]
  X <- as.matrix(dat %>% dplyr::select(starts_with('X')))[id.t, ]
  G <- as.matrix(dat %>% dplyr::select(starts_with('G')))[id.t, ]
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
  obj.b <- SKAT_Null_Model(y.b ~ X, out_type="D", Adjustment = F)
  # SKAT results
  SKAT_result <- SKAT(G, obj.b, method="davies", weights=SKATweights)
  # Burden results
  Burden_result <- SKAT(G, obj.b, method="davies", r.corr = 1, weights=Burdenweights)

  # ACAT NULL model
  G <- Matrix::Matrix(G, sparse = TRUE)
  obj <- ACAT::NULL_Model(y.b, X)
  ACAT_result <- ACAT::ACAT_V(G,obj)

  result <- list(SKAT_p = SKAT_result$p.value, Burden_p = Burden_result$p.value, ACAT_p = ACAT_result)
  return(result)
}
