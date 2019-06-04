

library(flexdashboard)
library(data.table)
library(ComplexHeatmap)
library(tokenizers)
library(dplyr)
library(DT)
library(RSQLite)
library(readxl)

sig_transcripts <- readRDS("data/heatmap_sig_transcripts.rds")
sig <- readRDS("data/heatmap_sig_transcripts_bytime.rds")
LFC_all <- readRDS("data/heatmap_LFC.rds")
heatmap_transcripts_to_genes <- readRDS("data/heatmap_transcripts_to_genes.rds")
overall_LRT_pvalues <- readRDS("data/heatmap_overall_LRT_pvalues.rds")
logtpm <- readRDS("data/heatmap_logtpm.rds")
qval_all <- readRDS("data/qval_all.rds")
## Reorder LFC and qval so that it is in right order
LFC_all <- LFC_all[,c(2,1,3,4)]
qval_all <- qval_all[,c(2,1,3,4)]
colnames(LFC_all) <- strsplit(colnames(LFC_all), split = ".b", fixed= TRUE) %>% unlist()
colnames(qval_all) <- strsplit(colnames(qval_all), split = ".qval", fixed= TRUE) %>% unlist()
## Rerrange LFC_all, qval_all, logtpm, and overall_LRT_pvalues, and heatmap_transcripts_to_genes to have matching order
LFC_all <- LFC_all[match(rownames(logtpm), rownames(LFC_all)),]
qval_all <- qval_all[match(rownames(logtpm), rownames(qval_all)),]
overall_LRT_pvalues <- overall_LRT_pvalues[match(rownames(logtpm), overall_LRT_pvalues$target_id),]
heatmap_transcripts_to_genes <-  
  heatmap_transcripts_to_genes[match(rownames(logtpm), 
                                     heatmap_transcripts_to_genes$target_id),]

## Read in results at a padj = 0.05 threshold (7480 genes)
sig_transcripts_05 <- character()
sig <- vector("list", length = 4)
names(sig) <- c("4dpi-0dpi", "2dpi-0dpi", "7dpi-0dpi",  "12dpi-0dpi")
for(i in c("4dpi-0dpi", "2dpi-0dpi", "7dpi-0dpi",  "12dpi-0dpi")) {
  tmp <- read.csv(paste0("../RNAseq/4_vs0_model/results_", i, ".csv"))
  sig[[i]] <- as.character(tmp[which(tmp$qval < 0.05), "target_id"])
  sig_transcripts_05 <- c(sig_transcripts_05, sig[[i]])
}
sig_transcripts_05 <- unique(sig_transcripts_05)

LFC_05 <- LFC_all[which(rownames(LFC_all) %in% sig_transcripts_05),]
TPM_05 <- logtpm[which(rownames(logtpm) %in% sig_transcripts_05),]
TPM_Z_05 <- data.frame(t(scale(t(TPM_05), center = TRUE, scale = TRUE)), check.names=FALSE)

## Add gene names
LFC_05 <- data.frame(target_id = rownames(LFC_05), LFC_05, check.names=FALSE) %>%
  left_join(., heatmap_transcripts_to_genes, by = "target_id") %>%
  dplyr::select(target_id, ext_gene, everything())
TPM_05 <- data.frame(target_id = rownames(TPM_05), TPM_05, check.names=FALSE) %>%
  left_join(., heatmap_transcripts_to_genes, by = "target_id") %>%
  dplyr::select(target_id, ext_gene, everything())
TPM_Z_05 <- data.frame(target_id = rownames(TPM_Z_05), TPM_Z_05, check.names=FALSE) %>%
  left_join(., heatmap_transcripts_to_genes, by = "target_id") %>%
  dplyr::select(target_id, ext_gene, everything())

writexl::write_xlsx(list(LFC_7480genes_05 = LFC_05, TPM_7480genes_05 = TPM_05, 
                         TPM_Z_7480genes_05 = TPM_Z_05),
                    path = "data_for_7480genes_p05.xlsx")
