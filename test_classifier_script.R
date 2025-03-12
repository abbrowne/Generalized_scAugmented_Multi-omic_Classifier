
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
source("E:/Projects/BioApp_1/Bioinformatic_functions.R")

##Expression and metadata with outcome
test_data <- pullDataFromRecount3("PAAD","tcga")
test_data2 <- processRecount3Data_TCGA(test_data)

##Genesets
msigdb_hp_gmt <- import_msigdb_genesets("H")
pdac_gmt <- importPDAC_genesets()
all_gmts <- c(msigdb_hp_gmt,pdac_gmt)

##scRNAseq
test_sc_dataset <- readRDS("E:/Projects/Cancer/PDAC Chijimatsu reconstruction/scDatasets/pk_all.rds")

##Streamlined classifier run

test_result <- find_cellEnrichedPrognostic_scores(test_data2$log2TPM,
                                                  test_data2$meta,
                                                  "Outcomes_OS.time",
                                                  "Outcomes_OS",
                                                  test_sc_dataset,
                                                  all_gmts)
test_classifier_result <- create_ncvTSP_classifier(test_data2$log2TPM, 
                                                   test_result$final_progGroup, 
                                                   test_result$final_scMarkers,
                                                   test_result$final_meta)

test_classifier_calls <- apply_TSP_classifier(test_data2$log2TPM,
                                              test_classifier_result)



###Check results
table(test_classifier_calls[rownames(test_result$final_meta),"group"],test_result$final_progGroup)

test_result$final_meta$test_surv_group <- test_classifier_calls[rownames(test_result$final_meta),"group"]

# Create the Surv object
surv_object <- Surv(time = test_result$final_meta[,"Outcomes_OS.time"], event = test_result$final_meta[,"Outcomes_OS"])

# Kaplan-Meier fit by 'group'
fit1 <- survfit(surv_object ~ test_surv_group, data = test_result$final_meta)

# Use ggsurvplot for a survival curve of both classifier groups
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