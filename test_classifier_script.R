
##Load necessary packages for classifier testing example
#devtools::install_github('immunogenomics/presto')
library(dplyr)
library(tidyr)
library(TCGAbiolinks)
library(DT)
library(purrr)
library(recount3)
library(DESeq2)
library(edgeR)
library(limma)
library(GSVA)
library(qusage)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(rpart)
library(pdacR)
library(Seurat)
library(switchBox)
library(ncvreg)
library(stringr)
library(survival)
library(survminer)
library(ggplot2)
library(msigdbr)
library(msigdbdf)
#Load the helper file with documented custom functions for classifier workflow
source("E:/Projects/BioApp_1/Bioinformatic_functions.R")

##Pull expression and metadata with outcomes for TCGA PAAD using recount3 and xenabrowser data
test_data <- pullDataFromRecount3("PAAD","tcga")
test_data2 <- processRecount3Data_TCGA(test_data)

##Import relevant genesets from msigdb and pdacR
msigdb_hp_gmt <- import_msigdb_genesets("H")
pdac_gmt <- importPDAC_genesets()
all_gmts <- c(msigdb_hp_gmt,pdac_gmt)

##Load reference PDAC scRNAseq dataset from https://pubmed.ncbi.nlm.nih.gov/35847558/
test_sc_dataset <- readRDS("E:/Projects/Cancer/PDAC Chijimatsu reconstruction/scDatasets/pk_all.rds")

##Run classifier using existing functions with more detailed description in "Bioinformatic_functions.R" file

##Identify prognostic GSVA scores with cell type enrichment by scRNAseq
test_result <- find_cellEnrichedPrognostic_scores(test_data2$log2TPM,
                                                  test_data2$meta,
                                                  "Outcomes_OS.time",
                                                  "Outcomes_OS",
                                                  test_sc_dataset,
                                                  all_gmts)
##Derive an ncvTSP classifier for the prognostic and cell type enriched score result from the "find_cellEnrichedPrognostic_scores" function
test_classifier_result <- create_ncvTSP_classifier(test_data2$log2TPM, 
                                                   test_result$final_progGroup, 
                                                   test_result$final_scMarkers,
                                                   test_result$final_meta)
##Apply the derived classifier result from the "create_ncvTSP_classifier" function to the TCGA PAAD dataset
test_classifier_calls <- apply_TSP_classifier(test_data2$log2TPM,
                                              test_classifier_result)



##Check the classifier result versus the final prognostic group to ensure classifier accuracy is sufficient
table(test_classifier_calls[rownames(test_result$final_meta),"group"],test_result$final_progGroup)

##Run survival analysis on the classifier result to confirm its prognostic value when applied to the dataset
test_result$final_meta$test_surv_group <- test_classifier_calls[rownames(test_result$final_meta),"group"]
surv_object <- Surv(time = test_result$final_meta[,"Outcomes_OS.time"], event = test_result$final_meta[,"Outcomes_OS"])

##Create Kaplan-Meier fit by survival group
fit1 <- survfit(surv_object ~ test_surv_group, data = test_result$final_meta)

##Plot survival curve of binary classifier groups using ggsurvplot
temp_plot1 <- ggsurvplot(fit1, 
                         data = test_result$final_meta,
                         pval = TRUE,            # adds p-value for log-rank test
                         conf.int = TRUE,        # adds confidence intervals
                         risk.table = TRUE,      # add risk table at the bottom
                         xlab = "Time",
                         ylab = "Survival Probability",
                         legend.title = "Group",
                         legend.labs = levels(test_result$final_meta$test_surv_group),
                         title = "Overall survival with scClassifier groups")

print(temp_plot1)