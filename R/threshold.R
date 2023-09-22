# Updated: 2023-09-19

#' function to obtain threshold
#'
#' @param dat Basis matrix of data set to be analyzed.
#' @param prev_est Estimate of prevalence.
#' @return Parameter estimates, SS SKAT score and SS SKAT pvalue.
#' @examples
#' analysis_data <- example_data$my_dat;
#' labeled_data_id <- example_data$id.t;
#' Y_threshold <- get_Y_threshold(dat = analysis_data, prev_value = mean(analysis_data$Y[labeled_data_id]));
#' threshold_SKAT_results <- threshold_SKAT(dat = analysis_data, id.t = labeled_data_id, thre_value = Y_threshold);
#' @export

get_Y_threshold <- function(dat, prev_value) {
  S <- dat$S
  thre_value <- quantile(S, 1-prev_value)
  return(thre_value)
}
