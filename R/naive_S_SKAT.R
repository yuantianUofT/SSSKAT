# Updated: 2023-09-19

#' SKAT directly with S for both labeled and unlabeled data
#'
#' @param dat Basis matrix of data set to be analyzed.
#' @return Parameter estimates, SS SKAT score and SS SKAT pvalue.
#' analysis_data <- example_data$my_dat;
#' labeled_data_id <- example_data$id.t;
#' naive_SKAT_results <- naive_SKAT(dat = analysis_data);
#' @export

naive_SKAT <- function(dat) {
  # data preparation
  S <- dat$S
  X <- as.matrix(dat %>% dplyr::select(starts_with('X')))
  G <- as.matrix(dat %>% dplyr::select(starts_with('G')))

  # SKAT NULL model
  obj.b <- SKAT_Null_Model(S ~ X, out_type="C", Adjustment = F)

  # SKAT results
  result <- SKAT(G, obj.b, method="davies", weights.beta=rep(1, dim(G)[2]))
  return(result)
}
