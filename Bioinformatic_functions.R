###Useful bioinformatic functions


check_sc_enrichment <- function(input_sc_dataset,input_gmts,input_geneset_name){
  seurat_obj <- NULL
  seurat_obj <- try(
    AddModuleScore(
      object   = input_sc_dataset,
      features = list(input_gmts[[input_geneset_name]]),
      name     = "CustomGeneset"
    ),
    silent = TRUE
  )
  if(!is.null(seurat_obj)){
    temp_result <- as.data.frame(seurat_obj@meta.data %>% 
                                   group_by(Cell_type) %>% 
                                   summarize(mean_ModuleScore = mean(CustomGeneset1)))
    temp_min <- min(temp_result$mean_ModuleScore)
    temp_result$norm_Module_score <- temp_result$mean_ModuleScore - temp_min
    temp_max_i <- which.max(temp_result$norm_Module_score)
    temp_max_value <- temp_result$norm_Module_score[temp_max_i]
    temp_max_cell_type <- temp_result$Cell_type[temp_max_i]
    comp_result <- temp_result$norm_Module_score[temp_result$Cell_type != temp_max_cell_type]
    temp_result <- list("max_cell_type" = temp_max_cell_type,
                        "max_module_score" = temp_max_value,
                        "other_module_scores" = comp_result)
    return(temp_result)
  }else{
    temp_result <- list("max_cell_type" = NA,
                        "max_module_score" = 0,
                        "other_module_scores" = c(100,100))
    return(temp_result)
  }
}

find_cellEnrichedPrognostic_scores <- function(input_expression,
                                               input_meta,
                                               input_outcome_time,
                                               input_outcome_result,
                                               input_sc_dataset,
                                               input_gmts,
                                               max_iterations=100,
                                               input_cp_set=c(0.05,0.02),
                                               initial_minsplit_fraction=0.4,
                                               initial_minbucket_fraction=0.4,
                                               minsplit_delta=0.1,
                                               minbucket_delta=0.1,
                                               max_prog_features_per_fit=5,
                                               max_enrichment_iterations=10,
                                               plot_name="initial_survplot"
){
  
  
  if(sum(input_meta[,input_outcome_time] <= 0) > 0){
    print(paste0("Error with survival times of 0 or below. These samples will be removed"))
  }
  temp_meta <- input_meta[input_meta[,input_outcome_time] > 0,]
  
  
  ###Add GSVA scoring here for input expression data rather than meta with GSVA scores added
  cat(paste0("Running GSVA scoring on selected gene sets with included expression data.\n"))
  temp_gsva <- gsva(gsvaParam(input_expression,input_gmts))
  temp_gsva <- t(temp_gsva)[rownames(temp_meta),]
  colnames(temp_gsva) <- paste0("GSVA_",colnames(temp_gsva))
  temp_meta <- cbind(temp_meta,temp_gsva)
  
  ##Create survival formula for assessing all included GSVA scores in relation to overall survival
  temp_formula <- as.formula(paste0("Surv(",input_outcome_time,", ",input_outcome_result,") ~ ",
                                    paste(paste0("GSVA_",names(input_gmts)),collapse = " + ")))
  
  temp_minsplit <- round(nrow(temp_meta)*initial_minsplit_fraction)
  temp_minbucket <- round(nrow(temp_meta)*initial_minbucket_fraction)
  temp_cp_set <- input_cp_set
  
  flag_primary_cutpoint_sc_enrichment <- FALSE
  total_iterations <- 0
  temp_progFeature <- NULL
  temp_progCutoff <- NULL
  filtered_gsva_columns <- colnames(temp_meta)[grepl("GSVA_",colnames(temp_meta))]
  sub_gsva_columns <- filtered_gsva_columns
  cat(paste0("Running survival cutpoint analysis to identify the best prognostic feature and cutpoint with cell type association.\n"))
  while(flag_primary_cutpoint_sc_enrichment == FALSE & total_iterations <= max_iterations){
    if(total_iterations > 0){
      sub_gsva_columns <- filtered_gsva_columns
      temp_minsplit <- round(nrow(temp_meta)*(initial_minsplit_fraction - minsplit_delta))
      temp_minbucket <- round(nrow(temp_meta)*(initial_minsplit_fraction - minbucket_delta))
      temp_formula <- as.formula(paste0("Surv(",input_outcome_time,", ",input_outcome_result,") ~ ",
                                        paste(sub_gsva_columns,collapse = " + ")))
    }
    for(temp_cp_i in 1:length(temp_cp_set)){
      sub_gsva_columns <- filtered_gsva_columns
      temp_formula <- as.formula(paste0("Surv(",input_outcome_time,", ",input_outcome_result,") ~ ",
                                        paste(sub_gsva_columns,collapse = " + ")))
      temp_cp <- temp_cp_set[temp_cp_i]
      temp_fit <- rpart(temp_formula, data = temp_meta, method = "exp", 
                        control = rpart.control(cp=temp_cp,
                                                maxdepth=1,
                                                minsplit = temp_minsplit,
                                                minbucket = temp_minbucket))
      if(length(names(temp_fit$variable.importance)) > 0){
        for(temp_prog_feature_i in 1:min(max_prog_features_per_fit,length(names(temp_fit$variable.importance)))){
          temp_progFeature <- names(temp_fit$variable.importance)[temp_prog_feature_i]
          temp_progCutoff <- temp_fit$splits[temp_progFeature,"index"]
          cat(paste0("Testing cp ",temp_cp," minsplit ",temp_minsplit," and minbucket ",temp_minbucket," with prog feature ",temp_prog_feature_i," : ",temp_progFeature,"\n"))
          temp_meta$test_surv_group <- as.factor(temp_meta[,temp_progFeature] < temp_progCutoff)
          temp_surv_formula <- as.formula(paste0("Surv(",input_outcome_time,", ",input_outcome_result,") ~ test_surv_group"))
          fit_cox <- coxph(temp_surv_formula, data = temp_meta)
          temp_summary <- summary(fit_cox)
          print(temp_summary)
          temp_coxph_pval <- as.numeric(temp_summary$waldtest["pvalue"])
          if(temp_coxph_pval < 0.05){
            input_geneset_name <- str_replace(temp_progFeature,"GSVA_","")
            comp_result <- check_sc_enrichment(input_sc_dataset,input_gmts,input_geneset_name)
            if(sum(comp_result$max_module_score >= comp_result$other_module_scores*2) == length(comp_result$other_module_scores)){
              flag_primary_cutpoint_sc_enrichment <- TRUE
              temp_survFit <- eval(bquote(
                survfit(.(as.formula(
                  paste0("Surv(", input_outcome_time, ", ", input_outcome_result, ") ~ test_surv_group")
                )), data = .(temp_meta))
              ))
              temp_best_surv_plot <- ggsurvplot(temp_survFit, data = temp_meta,
                                                pval = TRUE,            # adds p-value for log-rank test
                                                conf.int = TRUE,        # adds confidence intervals
                                                risk.table = TRUE,      # add risk table at the bottom
                                                xlab = "Time",
                                                ylab = "Survival Probability",
                                                legend.title = paste0(temp_progFeature, " < ",round(temp_progCutoff,3)),
                                                legend.labs = levels(temp_meta$test_surv_group))
              pdf(paste0("E:/Projects/BioApp_1/Example_data/TCGA_PAAD/",plot_name,".pdf"))
              print(temp_best_surv_plot)
              dev.off()
              temp_max_cell_type <- comp_result$max_cell_type
              break
            }else{
              ##Add code here to remove genesets that fail to show cellular enrichment from future iterations
              filtered_gsva_columns <- filtered_gsva_columns[filtered_gsva_columns != temp_progFeature]
            }
          }
          sub_gsva_columns <- sub_gsva_columns[sub_gsva_columns != temp_progFeature]
          total_iterations <- total_iterations + 1
        }
      }
      if(flag_primary_cutpoint_sc_enrichment == TRUE){
        break
      }
      sub_gsva_columns <- sub_gsva_columns[sub_gsva_columns != temp_progFeature]
    }
    sub_gsva_columns <- sub_gsva_columns[sub_gsva_columns != temp_progFeature]
  }
  
  if(flag_primary_cutpoint_sc_enrichment == TRUE){
    
    cat(paste0("Running single cell marker derivation to obtain a marker list between 500 and 1000 genes in size.\n"))
    
    temp_thresholdFC <- 1.5
    temp_scMarker_size <- 0
    temp_scMarkers <- NULL
    enrichment_iterations <- 0
    while((temp_scMarker_size < 500 | temp_scMarker_size > 1000) & enrichment_iterations < max_enrichment_iterations){
      temp_scMarkers <- FindMarkers(
        object        = input_sc_dataset,
        ident.1       = temp_max_cell_type,       # cell type or cluster name
        only.pos      = TRUE,            # only keep genes positively enriched in "T cells"
        logfc.threshold = temp_thresholdFC,          # filter out genes with low logFC
        min.pct       = 0.1              # filter out genes not expressed in at least 10% of T cells
      )
      temp_scMarker_size <- nrow(temp_scMarkers)
      if(temp_scMarker_size > 1000){
        temp_thresholdFC <- temp_thresholdFC * 1.1
      }else if(temp_scMarker_size < 500){
        temp_thresholdFC <- temp_thresholdFC * 0.5
      }
      enrichment_iterations <- enrichment_iterations + 1
    }
    if(enrichment_iterations == max_enrichment_iterations){
      cat(paste0("Failed to identify sufficient numbers of enriched markers for the cell type associating with the identified prognostic feature.\n",
                 "The best obtained marker list will be returned despite being outside the optimal range of approximately 500-1000 genes.\n"))
    }
    final_prog_info <- list(final_meta=temp_meta,
                            final_progFeature=temp_progFeature,
                            final_progCutoff=temp_progCutoff,
                            final_progGroup=as.factor(as.numeric(temp_meta[,temp_progFeature] < temp_progCutoff)),
                            #final_survPlot=temp_best_surv_plot,
                            final_cp=temp_cp,
                            final_minsplit=temp_minsplit,
                            final_minbucket=temp_minbucket,
                            final_enriched_cellType=temp_max_cell_type,
                            final_scMarkers=rownames(temp_scMarkers),
                            final_scMarker_enrichment_threshold=temp_thresholdFC)
    return(final_prog_info)
  }else{
    temp_message <- paste0("Failed to identify meaningful prognostic groups for a score showing cell type enrichment.\n")
    return(temp_message)
  }
}

create_ncvTSP_classifier <- function(input_expression, input_group, input_scMarkers, input_meta=NULL){
  temp_expr <- input_expression[rownames(input_expression) %in% input_scMarkers,]
  if(!is.null(input_meta)){
    temp_expr <- t(t(temp_expr)[rownames(input_meta),])
  }
  
  temp_TSP_result <- SWAP.KTSP.Train(temp_expr,as.factor(input_group),krange = 50,FilterFunc = NULL)
  
  ##Format gene pairs into binary matrix of comparisons for each sample to use in logistic regression
  temp_resultMat <- list()
  for(temp_i in 1:nrow(temp_TSP_result$TSPs)){
    temp_name <- paste0(temp_TSP_result$TSPs[temp_i,1],"_",temp_TSP_result$TSPs[temp_i,2])
    temp_resultMat[[temp_name]] <- as.vector(as.integer(temp_expr[temp_TSP_result$TSPs[temp_i,1],] > temp_expr[temp_TSP_result$TSPs[temp_i,2],]))
  }
  temp_resultMat <- as.data.frame(temp_resultMat,row.names = colnames(temp_expr))
  
  ##Perform logistic regression using top scoring pair result matrix and classification group to derive coefficients for each gene pair
  test_ncv_result <- cv.ncvreg(temp_resultMat,
                               input_group, 
                               alpha=0.5, nfolds = nrow(temp_resultMat))
  
  temp_coefs <- coef(test_ncv_result)
  temp_coefs <- temp_coefs[temp_coefs != 0]
  temp_TSPs <- as.data.frame(list(geneA=unlist(lapply(names(temp_coefs)[-1], function(x) strsplit(x, "_")[[1]][1])),
                                  geneB=unlist(lapply(names(temp_coefs)[-1], function(x) strsplit(x, "_")[[1]][2])),
                                  coef=as.numeric(temp_coefs[-1])))
  
  ##Derive recommended probability cutoff
  all_TSPs <- c(temp_TSPs$geneA,temp_TSPs$geneB)
  sub_expr <- temp_expr[rownames(temp_expr) %in% all_TSPs,]
  temp_TSPs <- temp_TSPs[temp_TSPs$geneA %in% rownames(sub_expr) & temp_TSPs$geneB %in% rownames(sub_expr),]
  
  temp_resultMat <- list()
  for(temp_i in 1:nrow(temp_TSPs)){
    temp_name <- paste0(temp_TSPs$geneA[temp_i],"_",temp_TSPs$geneB[temp_i])
    temp_resultMat[[temp_name]] <- as.vector(as.integer(sub_expr[temp_TSPs$geneA[temp_i],] > sub_expr[temp_TSPs$geneB[temp_i],]))
  }
  temp_resultMat <- as.data.frame(temp_resultMat,row.names = colnames(sub_expr))
  
  ### p = exp(x*coef) / (1 + exp(x*coef))
  temp_resultVector <- rowSums(t(t(temp_resultMat) * temp_TSPs$coef))
  temp_resultProb <- exp(temp_resultVector) / (1 + exp(temp_resultVector))
  all_cutoff_values <- seq(1,99)/100
  min_errors <- length(input_group)
  for(temp_cutoff_i in 1:length(all_cutoff_values)){
    temp_cutoff_value <- all_cutoff_values[temp_cutoff_i]
    temp_table <- table(input_group,temp_resultProb > temp_cutoff_value)
    if(nrow(temp_table) == 2 & ncol(temp_table) == 2){
      temp_errors <- sum(temp_table[1,2] + temp_table[2,1])
      if(temp_errors < min_errors){
        min_errors <- temp_errors
        best_cutoffs <- c(temp_cutoff_value)
      }else if(temp_errors == min_errors){
        best_cutoffs <- c(best_cutoffs,temp_cutoff_value)
      }
    }
  }
  
  ##Return classifier
  temp_TSP_classifier <- list(intercept=as.numeric(temp_coefs[1]),TSPs=temp_TSPs,best_cutoff_probs=best_cutoffs)
  return(temp_TSP_classifier)
}


apply_TSP_classifier <- function(input_expression, input_classifier, use_with_missing=FALSE, best_cutoff_probability=NULL){
  
  if(is.null(best_cutoff_probability)){
    best_cutoff_probability <- input_classifier$best_cutoff_probs[1]
  }
  ##Check for missing classifier genes in input dataset with option to drop those pairs
  all_TSPs <- c(input_classifier$TSPs$geneA,input_classifier$TSPs$geneB)
  if(sum(all_TSPs %in% rownames(input_expression)) < length(all_TSPs)){
    missing_TSPs <- all_TSPs[!(all_TSPs %in% rownames(input_expression))]
    cat(paste0("The following genes were missing from the input dataset:\n",
               paste(missing_TSPs,collapse=", ")))
    cat(paste0("Full input dataset:\n"))
    print(input_classifier$TSPs)
    
    if(use_with_missing==TRUE){
      classifier_rows_to_drop <- c()
      for(temp_i in 1:nrow(input_classifier$TSPs)){
        if(!(input_classifier$TSPs$geneA[temp_i] %in% rownames(input_expression)) | 
           !(input_classifier$TSPs$geneB[temp_i] %in% rownames(input_expression))){
          classifier_rows_to_drop <- c(classifier_rows_to_drop,temp_i)
        }
      }
      input_classifier <- input_classifier$TSPs[-classifier_rows_to_drop,]
    }else{
      return(paste0("Error: Failed due to missing classifier genes in the input data"))
    }
  }
  
  all_TSPs <- c(input_classifier$TSPs$geneA,input_classifier$TSPs$geneB)
  temp_expr <- input_expression[rownames(input_expression) %in% all_TSPs,]
  
  temp_resultMat <- list()
  for(temp_i in 1:nrow(input_classifier$TSPs)){
    temp_name <- paste0(input_classifier$TSPs[temp_i,1],"_",input_classifier$TSPs[temp_i,2])
    temp_resultMat[[temp_name]] <- as.vector(as.integer(temp_expr[input_classifier$TSPs[temp_i,1],] > temp_expr[input_classifier$TSPs[temp_i,2],]))
  }
  temp_resultMat <- as.data.frame(temp_resultMat,row.names = colnames(temp_expr))
  
  ### p = exp(x*coef) / (1 + exp(x*coef))
  temp_resultVector <- rowSums(t(t(temp_resultMat) * input_classifier$TSPs$coef))
  temp_resultProb <- exp(temp_resultVector) / (1 + exp(temp_resultVector))
  temp_resultGroup <- as.factor(as.numeric(temp_resultProb > best_cutoff_probability))
  
  temp_classifier_result <- as.data.frame(list(prob=temp_resultProb,group=temp_resultGroup))
  rownames(temp_classifier_result) <- colnames(input_expression)
  return(temp_classifier_result)
}

count2tpm<- function(rse){
  count_matrix <- rse@assays@data$raw_counts
  gene_length <- rse@rowRanges$bp_length
  reads_per_rpk <- count_matrix/gene_length
  per_mil_scale <- colSums(reads_per_rpk)/1000000
  tpm_matrix <- t(t(reads_per_rpk)/per_mil_scale)
  return(tpm_matrix)
}

getCountMatrix<- function(rse){
  count_matrix <- rse@assays@data$raw_counts
  return(count_matrix)
}

pullDataFromRecount3 <- function(input_project_name, input_project_source){
  
  human_projects <- available_projects()
  if(input_project_name %in% unique(human_projects$project) & input_project_source %in% unique(human_projects$file_source)){
    temp_project_info = subset(
      human_projects,
      project == input_project_name & file_source == input_project_source & project_type == "data_sources"
    )
    temp_info <- map(seq(nrow(temp_project_info)), ~temp_project_info[.x, ])
    ## create the RangedSummarizedExperiment. the create_rse function works on one row a time 
    temp_dataset <- map(temp_info, ~create_rse(.x))
    temp_tpm_data<- map(temp_dataset, count2tpm)
    temp_count_data <- map(temp_dataset, getCountMatrix)
    ## get the metadata column 
    temp_meta_data<- map(temp_dataset, ~.x@colData@listData %>% as.data.frame())
    # bind the data matrix across cancer types together 
    temp_tpms<- purrr::reduce(temp_tpm_data, cbind)
    temp_counts <- purrr::reduce(temp_count_data,cbind)
    ## bind the metadata across cancer types together
    temp_meta<- purrr::reduce(temp_meta_data, rbind)
    rownames(temp_meta) <- temp_meta$external_id
    ##Ensure idential ID naming across datasets and save them
    #rownames(temp_meta) <- temp_meta$tcga.gdc_cases.samples.submitter_id
    #colnames(temp_tpms) <- rownames(temp_meta)
    #colnames(temp_counts) <- rownames(temp_meta)
    temp_result <- list(meta=temp_meta,counts=temp_counts,tpms=temp_tpms)
    return(temp_result)
  }else{
    return(paste0("Error with input"))
  }
}

processRecount3Data_TCGA <- function(temp_result,filter_genes=FALSE){
  #Load expression and metadata files
  input_counts <- test_data$counts
  input_counts <- input_counts[,order(colnames(input_counts))]
  input_tpm <- test_data$tpms
  input_tpm <- input_tpm[,order(colnames(input_tpm))]
  input_meta <- test_data$meta
  input_meta <- input_meta[order(rownames(input_meta)),]
  
  ##Filter non-expressed genes
  if(filter_genes == TRUE){
    keep <- rowSums(input_counts > 0) > 0
    input_counts <- input_counts[keep,]
    input_tpm <- input_tpm[keep,]
  }
  
  ##Load and add survival and sample phenotype data to metadata
  input_surv <- read.table("E:/Projects/Cancer/TCGA_PanCancer_Survival_SupplementalTable_S1_20171025_xena_sp",sep="\t",header=TRUE,row.names=1)
  colnames(input_surv) <- paste0("Outcomes_",colnames(input_surv))
  input_extra <- read.table("E:/Projects/Cancer/TCGA_PanCancer_phenotype_denseDataOnlyDownload.tsv",sep="\t",header=TRUE,row.names=1)
  colnames(input_extra) <- paste0("SamplePhenotype_",colnames(input_extra))
  input_meta$Outcomes_matchingID <- unlist(lapply(input_meta$tcga.tcga_barcode,function(x){
    substr(paste0(str_split(x,"-")[[1]][1:4],collapse = "-"),
           1,nchar(paste0(str_split(x,"-")[[1]][1:4],collapse = "-")) - 1)
  }))
  input_meta <- cbind(input_meta,input_surv[input_meta$Outcomes_matchingID,],input_extra[input_meta$Outcomes_matchingID,])
  ##Reformat and factor primary sample type column used as covariate in modeling
  input_meta$sample_type <- input_meta$SamplePhenotype_sample_type
  input_meta$sample_type[input_meta$sample_type == "Primary Tumor"] <- "Tumor"
  input_meta$sample_type[input_meta$sample_type == "Metastatic"] <- "Met"
  input_meta$sample_type[input_meta$sample_type == "Solid Tissue Normal"] <- "Normal"
  input_meta$sample_type <- factor(input_meta$sample_type,levels=c("Normal","Tumor","Met"))
  ##Get HGNC symbols for remaining genes and replace ENSEMBL rownames
  gene_ids <- rownames(input_counts)
  
  # Remove version suffix if present
  gene_ids_clean <- sub("\\..*", "", gene_ids)
  # Connect to the Ensembl database
  # You can specify an Ensembl release if needed, e.g., useEnsembl(biomart="ensembl", 
  # version=109, dataset="hsapiens_gene_ensembl")
  symbol_map <- mapIds(
    x         = org.Hs.eg.db,
    keys      = gene_ids_clean,
    column    = "SYMBOL",     # We want HGNC gene symbols
    keytype   = "ENSEMBL",    # Our keys are Ensembl gene IDs
    multiVals = "first"       # If multiple symbols map to one ID, take the first
  )
  
  ##Replace ensembl IDs with symbols if not NA and not duplicated
  temp_rownames <- ifelse(
    is.na(symbol_map[gene_ids_clean]),
    gene_ids_clean,
    symbol_map[gene_ids_clean]
  )
  rownames(input_counts) <- ifelse(
    duplicated(temp_rownames),
    gene_ids_clean,
    temp_rownames
  )
  rownames(input_tpm) <- rownames(input_counts)
  
  if(all(colnames(input_counts) %in% rownames(input_meta)) & all(colnames(input_counts) == rownames(input_meta))){
    input_project_dataset <- list(meta=input_meta,log2TPM=log2(input_tpm+1))
    return(input_project_dataset)
  }else{
    return(paste0("Error with metadata and expression column numbers"))
  }
}

count2tpm<- function(rse){
  count_matrix <- rse@assays@data$raw_counts
  gene_length <- rse@rowRanges$bp_length
  reads_per_rpk <- count_matrix/gene_length
  per_mil_scale <- colSums(reads_per_rpk)/1000000
  tpm_matrix <- t(t(reads_per_rpk)/per_mil_scale)
  return(tpm_matrix)
}

getCountMatrix<- function(rse){
  count_matrix <- rse@assays@data$raw_counts
  return(count_matrix)
}

importPDAC_genesets <- function(){
  pdac_genesets <- pdacR::gene_lists
  pdac_genesets$Moffitt.Tumor <- NULL
  pdac_genesets$Puleo.Centroids <- NULL
  pdac_genesets$ICGC.SAM <- NULL
  pdac_genesets$Moffitt.Top5s <- NULL
  pdac_genesets$Moffitt.Top25s <- NULL
  return(pdac_genesets)
}

import_msigdb_genesets <- function(target_collection){
  msig_h <- msigdbdf::msigdbdf(target_species = "HS")
  
  msig_h$gs_collection_name <- paste0(msig_h$gs_collection,"_",msig_h$gs_name)
  
  sub_msig <- msig_h[msig_h$gs_collection == target_collection,]
  msig_list <- sub_msig %>%
    group_by(gs_collection_name) %>%
    summarize(genes = list(unique(db_gene_symbol)))
  names(msig_list$genes) <- msig_list$gs_collection_name
  msig_list <- msig_list$genes
  return(msig_list)
}
