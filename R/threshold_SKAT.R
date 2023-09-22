# Updated: 2023-09-19

#' threshold SKAT function using the labeled data and thresholded unlabeled data
#'
#' @param dat Basis matrix of data set to be analyzed.
#' @param id.t Row id of labeled data in the dat dataset.
#' @param thre_value The Threshold value for the surrogate.
#' @return Parameter estimates, SS SKAT score and SS SKAT pvalue.
#' @examples
#' analysis_data <- example_data$my_dat;
#' labeled_data_id <- example_data$id.t;
#' Y_threshold <- get_Y_threshold(dat = analysis_data, prev_value = mean(analysis_data$Y[labeled_data_id]));
#' threshold_SKAT_results <- threshold_SKAT(dat = analysis_data, id.t = labeled_data_id, thre_value = Y_threshold);
#' @export

threshold_SKAT <- function(dat, id.t, thre_value) {
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

  # SKAT NULL model
  obj.b <- SKAT_Null_Model(y ~ X, out_type="D", Adjustment = F)

  # SKAT results
  result <- SKAT(G, obj.b, method="davies", weights.beta=rep(1, dim(G)[2]))
  return(result)
}
