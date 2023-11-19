# Updated: 2023-11-18

#' threshold function using the labeled data and thresholded unlabeled data
#'
#' @param dat Basis matrix of data set to be analyzed.
#' @param id.t Row id of labeled data in the dat dataset.
#' @param thre_value The Threshold value for the surrogate.
#' @param weights.beta Weights of the SNPs, default NULL treat the SNPs equally, else take parameters of beta distribution.
#' @return Thresholded pvalues.
#' @examples
#' analysis_data <- example_data$my_dat;
#' labeled_data_id <- example_data$id.t;
#' Y_threshold <- get_Y_threshold(dat = analysis_data, prev_value = mean(analysis_data$Y[labeled_data_id]));
#' threshold_SKAT_results <- threshold_test(dat = analysis_data, id.t = labeled_data_id, thre_value = Y_threshold, weights.beta = NULL);
#' @export

threshold_test <- function(dat, id.t, thre_value, weights.beta = NULL) {
  # data preparation
  y <- rep(NA, nrow(dat))
  y[id.t] <- dat$Y[id.t]
  S <- dat$S
  X <- as.matrix(dat %>% dplyr::select(starts_with('X')))
  G <- as.matrix(dat %>% dplyr::select(starts_with('G')))

  # threshold y
  y.s <- rep(0, nrow(dat))
  y.s[which(S > thre_value)] <- 1
  y[which(is.na(y))] <- y.s[which(is.na(y))]

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
  obj.b <- SKAT_Null_Model(y ~ X, out_type="D", Adjustment = F)
  # SKAT results
  SKAT_result <- SKAT(G, obj.b, method="davies", weights=SKATweights)
  # Burden results
  Burden_result <- SKAT(G, obj.b, method="davies", r.corr = 1, weights=Burdenweights)

  # ACAT NULL model
  G <- Matrix::Matrix(G, sparse = TRUE)
  obj <- ACAT::NULL_Model(y, X)
  ACAT_result <- ACAT::ACAT_V(G,obj)

  result <- list(SKAT_p = SKAT_result$p.value, Burden_p = Burden_result$p.value, ACAT_p = ACAT_result)
  return(result)

}
