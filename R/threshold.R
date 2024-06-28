# Updated: 2024-06-28

#' function to obtain threshold
#'
#' @param S Continuous surrogate
#' @param prev_est Estimate of prevalence.
#' @return Parameter estimates, SS SKAT score and SS SKAT pvalue.
#' @export

get_Y_threshold <- function(S, prev_value) {
  thre_value <- quantile(S, 1-prev_value)
  return(thre_value)
}
