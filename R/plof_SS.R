plof_SS <- function(chr, gene_name, genofile, obj_SS, genes,
                    QC_label="annotation/info/QC_label", variant_type=c("SNV","Indel","variant"), geno_missing_imputation=c("mean","minor"),
                    Annotation_dir="annotation/info/FunctionalAnnotation", Annotation_name_catalog,
                    Use_annotation_weights=c(TRUE,FALSE), Annotation_name=NULL,
                    rare_maf_cutoff = 0.01, rv_num_cutoff = 1, mac.thresh = 10,
                    boot = T, para_results){

	## evaluate choices
	variant_type <- match.arg(variant_type)
	geno_missing_imputation <- match.arg(geno_missing_imputation)

	# pheno, covar, surro, theta_est
	phen <- obj_SS$y
	surro <- obj_SS$s
	covar <- obj_SS$covar
	id.t <- obj_SS$id.t
	patient.id <- obj_SS$patient.id
	theta_est <- obj_SS$theta_est

	## get SNV id, position, REF, ALT (whole genome)
	filter <- seqGetData(genofile, QC_label)
	if(variant_type=="variant")
	{
		SNVlist <- filter == "PASS"
	}

	if(variant_type=="SNV")
	{
		SNVlist <- (filter == "PASS") & isSNV(genofile)
	}

	if(variant_type=="Indel")
	{
		SNVlist <- (filter == "PASS") & (!isSNV(genofile))
	}

	position <- as.numeric(seqGetData(genofile, "position"))
	variant.id <- seqGetData(genofile, "variant.id")

	rm(filter)
	gc()

	### Gene
	kk <- which(genes[,1]==gene_name)

	sub_start_loc <- genes[kk,3]
	sub_end_loc <- genes[kk,4]

	is.in <- (SNVlist)&(position>=sub_start_loc)&(position<=sub_end_loc)
	variant.id.gene <- variant.id[is.in]

	seqSetFilter(genofile,variant.id=variant.id.gene,sample.id=patient.id)

	## plof
	## Gencode_Exonic
	GENCODE.EXONIC.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.EXONIC.Category")]))
	## Gencode
	GENCODE.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Category")]))

	variant.id.gene <- seqGetData(genofile, "variant.id")
	lof.in.plof <- (GENCODE.EXONIC.Category=="stopgain")|(GENCODE.EXONIC.Category=="stoploss")|(GENCODE.Category=="splicing")|(GENCODE.Category=="exonic;splicing")|(GENCODE.Category=="ncRNA_splicing")|(GENCODE.Category=="ncRNA_exonic;splicing")
	variant.id.gene <- variant.id.gene[lof.in.plof]

	seqSetFilter(genofile,variant.id=variant.id.gene,sample.id=patient.id)

	## genotype id
	id.genotype <- seqGetData(genofile,"sample.id")
	# id.genotype.match <- rep(0,length(id.genotype))

	id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
	patient.id.merge <- data.frame(patient.id)
	patient.id.merge <- dplyr::left_join(patient.id.merge,id.genotype.merge,by=c("patient.id"="id.genotype"))
	id.genotype.match <- patient.id.merge$index

	## Genotype
	Geno <- seqGetData(genofile, "$dosage")
	Geno <- Geno[id.genotype.match,,drop=FALSE]

	## impute missing
	if(!is.null(dim(Geno)))
	{
		if(dim(Geno)[2]>0)
		{
			if(geno_missing_imputation=="mean")
			{
				Geno <- matrix_flip_mean(Geno)$Geno
			}
			if(geno_missing_imputation=="minor")
			{
				Geno <- matrix_flip_minor(Geno)$Geno
			}
		}
	}

	# ## Annotation
	# Anno.Int.PHRED.sub <- NULL
	# Anno.Int.PHRED.sub.name <- NULL
	# 
	# if(variant_type=="SNV") {
	# 	if(Use_annotation_weights)
	# 	{
	# 		for(k in 1:length(Annotation_name))
	# 		{
	# 			if(Annotation_name[k]%in%Annotation_name_catalog$name)
	# 			{
	# 				Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,Annotation_name[k])
	# 				Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name==Annotation_name[k])]))
	# 
	# 				if(Annotation_name[k]=="CADD")
	# 				{
	# 					Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
	# 				}
	# 
	# 				if(Annotation_name[k]=="aPC.LocalDiversity")
	# 				{
	# 					Annotation.PHRED.2 <- -10*log10(1-10^(-Annotation.PHRED/10))
	# 					Annotation.PHRED <- cbind(Annotation.PHRED,Annotation.PHRED.2)
	# 					Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,paste0(Annotation_name[k],"(-)"))
	# 				}
	# 				Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub,Annotation.PHRED)
	# 			}
	# 		}
	# 
	# 		Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
	# 		colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
	# 	}
	# }

	SS_result <- STAAR_SS(genotype=Geno, para_results = para_results, phen=phen, surro=surro, covar=covar, id.t=id.t, annotation_phred = NULL,
	                      rare_maf_cutoff=rare_maf_cutoff, rv_num_cutoff=rv_num_cutoff, mac.thresh=mac.thresh,
	                      boot=boot)
	SS_p_STAARO <- ACAT(Pvals = c(SS_result$results_SS_1_1$pvalue, SS_result$results_SS_1_25$pvalue))
	naive_p_STAARO <- ACAT(Pvals = c(SS_result$results_naive_1_1$pvalue, SS_result$results_naive_1_25$pvalue))
	thresholded_p_STAARO <- ACAT(Pvals = c(SS_result$results_thresholded_1_1$pvalue, SS_result$results_thresholded_1_25$pvalue))
	oracle_p_STAARO <- ACAT(Pvals = c(SS_result$results_oracle_1_1$pvalue, SS_result$results_oracle_1_25$pvalue))
	if (any(is.na(c(SS_result$results_labeled_1_1, SS_result$results_labeled_1_25)))) {
	  labeled_p_STAARO <- NA
	} else {
	  labeled_p_STAARO <- ACAT(Pvals = c(SS_result$results_labeled_1_1$pvalue, SS_result$results_labeled_1_25$pvalue))
	}
	

	results <- c(NA, 40)
	results[3] <- "plof"
	results[2] <- chr
	results[1] <- as.character(genes[kk,1])
	results[4] <- SS_result$num_variant
	results[5] <- SS_result$cMAC
	results[6:11] <- c(SS_result$results_SS_1_1$pvalue, SS_result$results_SS_1_25$pvalue)
	results[12] <- SS_p_STAARO
	results[13:18] <- c(SS_result$results_naive_1_1$pvalue, SS_result$results_naive_1_25$pvalue)
	results[19] <- naive_p_STAARO
	results[20:25] <- c(SS_result$results_thresholded_1_1$pvalue, SS_result$results_thresholded_1_25$pvalue)
	results[26] <- thresholded_p_STAARO
	results[27:32] <- c(SS_result$results_oracle_1_1$pvalue, SS_result$results_oracle_1_25$pvalue)
	results[33] <- oracle_p_STAARO
	if (any(is.na(c(SS_result$results_labeled_1_1, SS_result$results_labeled_1_25)))) {
	  results[34:39] <- NA
	  results[40] <- NA
	} else {
	  results[34:39] <- c(SS_result$results_labeled_1_1$pvalue, SS_result$results_labeled_1_25$pvalue)
	  results[40] <- labeled_p_STAARO
	}
	
	results <- as.data.frame(t(results))
	
	if(!is.null(results)) {
	  colnames(results)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
	  colnames(results)[6:dim(results)[2]] <- c("SS_SKAT_p_1_1", "SS_Burden_p_1_1", "SS_ACAT_V_p_1_1", 
	                                            "SS_SKAT_p_1_25", "SS_Burden_p_1_25", "SS_ACAT_V_p_1_25", 
	                                            "SS_STAAR-O", 
	                                            "Naive_SKAT_p_1_1", "Naive_Burden_p_1_1", "Naive_ACAT_V_p_1_1",
	                                            "Naive_SKAT_p_1_25", "Naive_Burden_p_1_25", "Naive_ACAT_V_p_1_25",
	                                            "Naive_STAAR-O",
	                                            "Thresholded_SKAT_p_1_1", "Thresholded_Burden_p_1_1", "Thresholded_ACAT_V_p_1_1",
	                                            "Thresholded_SKAT_p_1_25", "Thresholded_Burden_p_1_25", "Thresholded_ACAT_V_p_1_25",
	                                            "Thresholded_STAAR-O",
	                                            "Oracle_SKAT_p_1_1", "Oracle_Burden_p_1_1", "Oracle_ACAT_V_p_1_1",
	                                            "Oracle_SKAT_p_1_25", "Oracle_Burden_p_1_25", "Oracle_ACAT_V_p_1_25",
	                                            "Oracle_STAAR-O",
	                                            "Labeled_SKAT_p_1_1", "Labeled_Burden_p_1_1", "Labeled_ACAT_V_p_1_1",
	                                            "Labeled_SKAT_p_1_25", "Labeled_Burden_p_1_25", "Labeled_ACAT_V_p_1_25",
	                                            "Labeled_STAAR-O")
	}
	results_list <- list(results = results, SS_result = SS_result)

	seqResetFilter(genofile)

	return(results_list)
}

