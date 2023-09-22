# Updated: 2023-09-19

#' Standard SKAT function using the labeled data only
#'
#' @param dat Basis matrix of data set to be analyzed.
#' @param id.t Row id of labeled data in the dat dataset.
#' @return Parameter estimates, SS SKAT score and SS SKAT pvalue.
#' @examples
#' analysis_data <- example_data$my_dat;
#' labeled_data_id <- example_data$id.t;
#' labeled_SKAT_result <- labeled_SKAT(dat = analysis_data, id.t = labeled_data_id);
#' @export

labeled_SKAT <- function(dat, id.t) {
  # prepare data
  y.b <- dat$Y[id.t]
  X <- as.matrix(dat %>% dplyr::select(starts_with('X')))[id.t, ]
  G <- as.matrix(dat %>% dplyr::select(starts_with('G')))[id.t, ]

  # SKAT NULL model
  obj.b <- SKAT_Null_Model(y.b ~ X, out_type="D", Adjustment = F)

  # SKAT results
  result <- SKAT(G, obj.b, method="davies", weights.beta=c(1, dim(G)[2]))
  return(result)
}
