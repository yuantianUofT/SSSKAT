# Updated: 2024-06-28

#' STAAR_SS function
#'
#' @param genotype Genotype matrix
#' @param phen Phenotype
#' @param surro Surrogate
#' @param covar Covariate
#' @param id.t Labeled id
#' @param annotation_phred Annotation phred
#' @param rare_maf_cutoff Rare MAF cutoff
#' @param rv_num_cutoff RV number cutoff
#' @param mac.thresh MAC threshold
#' @param boot Parametric Bootstrapping
#' @param para_results Saved parametric Bootstrapping results
#' @param distri Distribution of S|Y, either "normal" or "beta"
#' @return STAAR SS, naive, thesholded, labeled, oracle results
#' @export

STAAR_SS <- function(genotype, phen, surro, covar, id.t, annotation_phred = NULL,
                     rare_maf_cutoff = 0.01, rv_num_cutoff, mac.thresh,
                     boot = T, para_results, distri){

  if(!inherits(genotype, "matrix") && !inherits(genotype, "Matrix")){
    stop("genotype is not a matrix!")
  }

  if(dim(genotype)[2] == 1){
    stop(paste0("Number of rare variant in the set is less than 2!"))
  }

  annotation_phred <- as.data.frame(annotation_phred)
  if(dim(annotation_phred)[1] != 0 & dim(genotype)[2] != dim(annotation_phred)[1]){
    stop(paste0("Dimensions don't match for genotype and annotation!"))
  }

  if(inherits(genotype, "sparseMatrix")){
    genotype <- as.matrix(genotype)
  }
  genotype <- matrix_flip(genotype)
  MAF <- genotype$MAF
  RV_label <- as.vector((MAF<rare_maf_cutoff)&(MAF>0))
  Geno_rare <- genotype$Geno[,RV_label]

  rm(genotype)
  gc()
  annotation_phred <- annotation_phred[RV_label,,drop=FALSE]

  if (sum(RV_label) >= rv_num_cutoff) {
    G <- Geno_rare
    MAF <- MAF[RV_label]
    rm(Geno_rare)
    gc()

    annotation_rank <- 1 - 10^(-annotation_phred/10)

    ## beta(1,25)
    w_1 <- dbeta(MAF,1,25)
    ## beta(1,1)
    w_2 <- dbeta(MAF,1,1)
    if(dim(annotation_phred)[2] == 0){
      ## beta(1,25)
      w_B_1 <- w_S_1 <- w_1
      w_B_2 <- w_S_2 <- w_2
      w_A_1 <- w_1^2/dbeta(MAF,0.5,0.5)^2
      w_A_2 <- w_2^2/dbeta(MAF,0.5,0.5)^2
    }else{
      ## Burden
      w_B_1 <- annotation_rank*w_1
      w_B_2 <- annotation_rank*w_2

      ## SKAT
      w_S_1 <- sqrt(annotation_rank)*w_1
      w_S_2 <- sqrt(annotation_rank)*w_2

      ## ACAT-V
      w_A_1 <- annotation_rank*w_1^2/dbeta(MAF,0.5,0.5)^2
      w_A_2 <- annotation_rank*w_2^2/dbeta(MAF,0.5,0.5)^2
    }
    weights_1 <- list(w_B_1 = as.vector(w_B_1), w_S_1 = as.vector(w_S_1), w_A_1 = as.vector(w_A_1))
    weights_2 <- list(w_B_2 = as.vector(w_B_2), w_S_2 = as.vector(w_S_2), w_A_2 = as.vector(w_A_2))
    
    ## SS
    SS_result1 <- SS_test(Y = phen, para_results = para_results, X = covar, G = G, S = surro, id.t = id.t, wBurden = weights_1$w_B_1, wSKAT = weights_1$w_S_1, wACAT = NULL, 
                          weights.beta = c(1, 25), mac.thresh = 10,
                          full_NR_evaluation = TRUE, nit = NULL, NULL_nlog_like, 
                          testtype = "all", boot = boot, distri = distri) 
    SS_result2 <- SS_test(Y = phen, para_results = para_results, X = covar, G = G, S = surro, id.t = id.t, wBurden = weights_2$w_B_2, wSKAT = weights_2$w_S_2, wACAT = NULL, 
                          weights.beta = c(1, 25), mac.thresh = 10,
                          full_NR_evaluation = TRUE, nit = NULL, NULL_nlog_like, 
                          testtype = "all", boot = boot, distri = distri) 
    ## Native
    naive_result1 <- naive_test(X = covar, G = G, S = surro, id.t = id.t, wBurden = weights_1$w_B_1, wSKAT = weights_1$w_S_1, wACAT = weights_1$w_A_1, 
                                mac.thresh = 10, testtype = "all")
    naive_result2 <- naive_test(X = covar, G = G, S = surro, id.t = id.t, wBurden = weights_2$w_B_2, wSKAT = weights_2$w_S_2, wACAT = weights_2$w_A_2,
                                mac.thresh = 10, testtype = "all")
    # thresholded
    thre_value <- get_Y_threshold(surro, prev_value = mean(phen[id.t]))
    thresholded_result1 <- threshold_test(Y = phen, X = covar, G = G, S = surro, id.t = id.t, thre_value = thre_value, wBurden = weights_1$w_B_1, wSKAT = weights_1$w_S_1, wACAT = weights_1$w_A_1, 
                                            mac.thresh = 10, testtype = "all")
    thresholded_result2 <- threshold_test(Y = phen, X = covar, G = G, S = surro, id.t = id.t, thre_value = thre_value, wBurden = weights_2$w_B_2, wSKAT = weights_2$w_S_2, wACAT = weights_2$w_A_2,
                                            mac.thresh = 10, testtype = "all")
    # oracle
    oracle_result1 <- labeled_test(Y = phen, X = covar, G = G, id.t = 1:length(surro), wBurden = weights_1$w_B_1, wSKAT = weights_1$w_S_1, wACAT = weights_1$w_A_1, 
                                  mac.thresh = 10, testtype = "all")
    oracle_result2 <- labeled_test(Y = phen, X = covar, G = G, id.t = 1:length(surro), wBurden = weights_2$w_B_2, wSKAT = weights_2$w_S_2, wACAT = weights_2$w_A_2,
                                  mac.thresh = 10, testtype = "all")
    # labeled
    labeled_result1 <- suppressWarnings(tryCatch({
       labeled_test(Y = phen, X = covar, G = G, id.t = id.t, wBurden = weights_1$w_B_1, wSKAT = weights_1$w_S_1, wACAT = weights_1$w_A_1, 
                                      mac.thresh = 10, testtype = "all")
      
    }, error = function(e) {
      return(NA)
    }))
    labeled_result2 <- suppressWarnings(tryCatch({
       labeled_test(Y = phen, X = covar, G = G, id.t = id.t, wBurden = weights_2$w_B_2, wSKAT = weights_2$w_S_2, wACAT = weights_2$w_A_2,
                                      mac.thresh = 10, testtype = "all")
      
    }, error = function(e) {
      return(NA)
    }))
    
    
    # final results
    num_variant <- sum(RV_label) #dim(G)[2]
    cMAC <- sum(G)
    num_annotation <- dim(annotation_phred)[2]+1
    
    final_results <- list(num_variant = num_variant,
                          cMAC = cMAC,
                          RV_label = RV_label,
                          num_annotation = num_annotation,
                          results_SS_1_1 = SS_result1,
                          results_SS_1_25 = SS_result2,
                          results_naive_1_1 = naive_result1,
                          results_naive_1_25 = naive_result2,
                          results_thresholded_1_1 = thresholded_result1,
                          results_thresholded_1_25 = thresholded_result2,
                          results_oracle_1_1 = oracle_result1,
                          results_oracle_1_25 = oracle_result2,
                          results_labeled_1_1 = labeled_result1,
                          results_labeled_1_25 = labeled_result2)
    return(final_results)
  }else{
    stop(paste0("Number of rare variant in the set is less than ",rv_num_cutoff,"!"))
  }

}

