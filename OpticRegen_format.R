
library(org.Dr.eg.db)
library(TxDb.Drerio.UCSC.danRer10.refGene)
library(rtracklayer)
library(AnnotationHub)
library(forcats)
library(flexdashboard)
library(data.table)
library(biomaRt)
library(ComplexHeatmap)
library(tokenizers)
library(readxl)
library(writexl)
library(GenomicFeatures)
library(ChIPpeakAnno)
library(topGO)
library(tidyverse)
library(sleuth)
library(readxl)
library(tidyverse)
library(RSQLite)

#------------------------------------------------------------------------------------------------------------------
## Format GO analysis
#------------------------------------------------------------------------------------------------------------------

## Collect gene names from biomaRt
mart <- biomaRt::useMart(biomart = "ensembl", 
                         dataset = "drerio_gene_ensembl")
## Get gene names and GO terms, remove blank entries, convert to list
GTOGO <- biomaRt::getBM(attributes = c("external_gene_name", "go_id"), 
                        mart = mart)
GTOGO <- GTOGO[GTOGO$go_id != '',]
geneID2GO <- by(GTOGO$go_id, GTOGO$external_gene_name, function(x) as.character(x))
saveRDS(geneID2GO, "data/geneID2GO.rds")


ensembl <- useMart(host = 'mar2016.archive.ensembl.org', 
                   biomart='ENSEMBL_MART_ENSEMBL')
ds <- listDatasets(ensembl) %>% arrange(description)
index <- c(14, 27, 37, 56)
ds <- bind_rows(ds[index,], ds[-index,])
saveRDS(ds, "data/ensembl_ds.rds")

#------------------------------------------------------------------------------------------------------------------
## Format significant transcripts for each comparison (heatmap data)
#------------------------------------------------------------------------------------------------------------------

## Read in results
sig_transcripts <- LFC_all <- qval_all <- character()
sig <- LFC <- qval <- vector("list", length = 4)
names(sig) <- names(LFC) <- names(qval) <- c("4dpi-0dpi", "2dpi-0dpi", "7dpi-0dpi",  "12dpi-0dpi")
for(i in c("4dpi-0dpi", "2dpi-0dpi", "7dpi-0dpi",  "12dpi-0dpi")) {
  tmp <- read.csv(paste0("../RNAseq/4_vs0_model/results_", i, ".csv"))
  sig[[i]] <- as.character(tmp[which(tmp$qval < 0.01), "target_id"])
  LFC[[i]] <- tmp[, c("target_id", "b")]
  qval[[i]] <- tmp[, c("target_id", "qval")]
  sig_transcripts <- c(sig_transcripts, sig[[i]])
}
sig_transcripts <- unique(sig_transcripts)
LFC_all <- full_join(LFC[[1]], LFC[[2]], by = "target_id") %>%
  full_join(., LFC[[3]], by = "target_id") %>%
  full_join(., LFC[[4]], by = "target_id")
colnames(LFC_all) <- c("target_id", "4dpi-0dpi.b", "2dpi-0dpi.b", "7dpi-0dpi.b", "12dpi-0dpi.b")
qval_all <- full_join(qval[[1]], qval[[2]], by = "target_id") %>%
  full_join(., qval[[3]], by = "target_id") %>%
  full_join(., qval[[4]], by = "target_id")
colnames(qval_all) <- c("target_id", "4dpi-0dpi.qval", "2dpi-0dpi.qval", "7dpi-0dpi.qval", "12dpi-0dpi.qval")
rownames(LFC_all) <- LFC_all[,1]
rownames(qval_all) <- qval_all[,1]
LFC_all <- LFC_all %>% dplyr::select(-contains("target_id"))
qval_all <- qval_all %>% dplyr::select(-contains("target_id"))

transcripts_to_genes <- read.csv(paste0("../RNAseq/4_vs0_model/results_", i, ".csv")) %>%
   dplyr::select(target_id, ext_gene)

## Load overall p-values from LRT
load("../RNAseq/2_splines_model.RData")
overall_LRT_pvalues <- results_table_factor %>% 
  dplyr::select(target_id, qval, ext_gene)
## Load TPM values
source("../RNAseq/my_sleuth_functions.R")
load("../RNAseq/4_vs0_model.RData")
obj <- so_factor_vs0
tabd_df <- obj$obs_norm
tabd_df <- dplyr::select(tabd_df, target_id, sample, tpm)
tabd_df_wide <- reshape2::dcast(tabd_df, target_id ~ sample, value.var = "tpm")
rownames(tabd_df_wide) <- tabd_df_wide$target_id
tabd_df_wide$target_id <- NULL
time_pt <- strsplit(colnames(tabd_df_wide), split="_", fixed=TRUE) %>%
  lapply(., function(x) x[3]) %>% unlist %>% strsplit(., split="RNA", fixed=TRUE) %>%
  lapply(., function(x) x[1]) %>% unlist %>% as.numeric()
replicate <- strsplit(colnames(tabd_df_wide), split="_", fixed=TRUE) %>%
  lapply(., function(x) x[3]) %>% unlist %>% strsplit(., split="RNA", fixed=TRUE) %>%
  lapply(., function(x) x[2]) %>% unlist %>% as.numeric()
colnames(tabd_df_wide) <- paste0(time_pt, "dpi (rep ", replicate, ")")
o <- order(time_pt)
tabd_df_wide <- log(tabd_df_wide[,o] + 1)
logtpm <- tabd_df_wide

## Save RDS
saveRDS(sig_transcripts, file = "data/heatmap_sig_transcripts.rds")
saveRDS(sig, file = "data/heatmap_sig_transcripts_bytime.rds")
saveRDS(LFC_all, file = "data/heatmap_LFC.rds")
saveRDS(qval_all, file = "data/qval_all.rds")

saveRDS(transcripts_to_genes, file = "data/heatmap_transcripts_to_genes.rds")
saveRDS(overall_LRT_pvalues, file = "data/heatmap_overall_LRT_pvalues.rds")
saveRDS(tabd_df_wide, file = "data/heatmap_logtpm.rds")


## Add human gene symbols here
transcript_to_gene <- readRDS("data/transcripts_to_genes_map.rds")
ensembl <- useMart(host = 'mar2016.archive.ensembl.org', 
                   biomart='ENSEMBL_MART_ENSEMBL')
human <- useDataset(ds[2,1], mart=ensembl)
zebrafish <- useDataset("drerio_gene_ensembl", mart=ensembl)
human_gene_ids <- getLDS(attributes=c("ensembl_gene_id"), 
                         filters="ensembl_gene_id",
                         values=transcript_to_gene$ens_gene,
                         mart=zebrafish, attributesL=c("ensembl_gene_id"), martL=human)
human_gene_ids <- unique(human_gene_ids)
human_symbols <- getBM(filters= "ensembl_gene_id",
      attributes=c("ensembl_gene_id", "hgnc_symbol"),
      values=human_gene_ids[,2], mart=human)
human_symbols <- unique(human_symbols)
colnames(human_gene_ids) <- c("ens_gene", "human_ens_gene")
colnames(human_symbols) <- c("human_ens_gene", "hgnc_symbol")

human_ids <- full_join(human_gene_ids, human_symbols, by = "human_ens_gene")
saveRDS(human_ids, file = "data/human_ids.rds")
# transcript_to_gene <- left_join(transcript_to_gene,  human_ids, by = "ens_gene")
# saveRDS(transcript_to_gene, file = "data/transcripts_to_genes_map_hgnc.rds")

# transcripts_to_genes <- readRDS("data/heatmap_transcripts_to_genes.rds")
# transcripts_to_genes <- left_join(transcripts_to_genes, transcript_to_gene, by = c("target_id", "ext_gene")) 
# saveRDS(transcripts_to_genes, file = "data/heatmap_transcripts_to_genes_hgnc.rds")



#------------------------------------------------------------------------------------------------------------------
  
## Read in peaklet sequences, import the MACS output (skip the transgene), convert chromosome names using IGV alias
alias <- read.table("../ATACseq/danRer10_alias.tab")
colnames(alias) <- c("chr_old", "chr")
sequences <- read.table("../ATACseq/PEAKS_TRANS_PEAKLETS_ALLTIMES/ALLMERGED_ATAC.nodup.unique.macs.peaklets_peaks.pvalsort.narrowPeak_500bp.seq",
                        stringsAsFactors = FALSE)
macs <- data.frame(peak=sequences[seq(from=1, to=nrow(sequences), by=2),], 
                   sequence=sequences[seq(from=2, to=nrow(sequences), by=2),]) %>%
  separate(peak, into=c("chr", "pos"), sep=":") %>%
  separate(chr, into=c("X", "chr_old"), sep=">") %>% dplyr::select(-X) %>%
  separate(pos, into=c("start", "end"), sep="-") %>%
  left_join(., alias, by="chr_old") %>%
  dplyr::select(chr_old, chr, start, end, sequence)
## Save peak sequences in df format
saveRDS(macs, "data/peak_sequences_df.rds")
write.table(macs, "data/peak_sequences_df.txt", row.names = FALSE)

macs_alias <- macs %>%
  dplyr::select(chr, start, end) %>%
  dplyr::mutate(seqnames=chr) %>%
  dplyr::select(seqnames, start, end)
macsOutput <- GRanges(macs_alias)

## Annotation Hub
ah <- AnnotationHub()
ahDb <- query(ah, pattern = c("Danio rerio", "EnsDb", 90))
ahEdb <- ahDb[[1]]
ens90.exons <- exons(ahEdb, columns = listColumns(ahEdb, "gene"))
mcols(ens90.exons) <- mcols(ens90.exons)$gene_id
ens90.exons <- unique(ens90.exons)

## Match up chromosome naming scheme
seqlevelsStyle(ens90.exons) <- "UCSC"
seq_tmp <- seqlevels(macsOutput)
seq_tmp[grep("chrUn", seq_tmp)] <- substr(seq_tmp[grep("chrUn", seq_tmp)], 7, 30)
seq_tmp[grep("v1", seq_tmp)] <- gsub("v", ".", seq_tmp[grep("v1", seq_tmp)])
seqlevels(macsOutput) <- seq_tmp

## Find peaks overlapping exons and TSS
peaks_intersect_exons <- findOverlaps(macsOutput, ens90.exons, 
                                      minoverlap = 50) %>% as.data.frame %>% dplyr::select(queryHits) %>% unique() %>% unlist()
peaks_overlapping_exons <- macsOutput[peaks_intersect_exons]
ens90.refGene <- genes(ahEdb)
ens90.tss <- promoters(ahEdb, upstream = 0, downstream = 0, columns = "gene_id")
seqlevelsStyle(ens90.refGene) <- "UCSC"
seqlevelsStyle(ens90.tss) <- "UCSC"

## Identify proximal and distal peaks
# Proximal peaks = overlapping TSS and/or within 1kb +/- TSS but not overlapping an exon by 50bp
# Distal peaks = in (1kb, 100kb) +/- TSS but not overlapping an exon by 50bp 
#ensembl_ids <- readRDS("data/transcripts_to_genes_map.rds") ## I think this isn't right as it takes ALL transcripts (even not expressed)
transcript_to_gene <- readRDS("data/transcripts_to_genes_map.rds")
ensembl_ids <- data.frame(target_id = rownames(readRDS("data/heatmap_LFC.rds"))) %>%
  left_join(., transcript_to_gene, by = "target_id") %>% na.omit()

ens90.tss_subset <- ens90.tss[which(unlist(ens90.tss$gene_id) %in% ensembl_ids$ens_gene)]
names(ens90.tss_subset) <- ens90.tss_subset$gene_id
ens90.refGene_subset <- ens90.refGene[which(unlist(ens90.refGene$gene_id) %in% ensembl_ids$ens_gene)]
names(ens90.refGene_subset) <- ens90.refGene_subset$gene_id
  
peaks_overlapping_tss_subset <- 
    annotatePeakInBatch(macsOutput, 
                        PeakLocForDistance = "middle",
                        FeatureLocForDistance = "start",
                        AnnotationData=ens90.tss_subset, 
                        output="overlapping", gap=-1L) 
peaks_overlapping_tss_subset <- 
    peaks_overlapping_tss_subset[which(!is.na(peaks_overlapping_tss_subset$feature))]
  
peaks_in_1kb_subset <- 
    annotatePeakInBatch(macsOutput, 
                        PeakLocForDistance = "middle",
                        FeatureLocForDistance = "TSS",
                        AnnotationData=ens90.refGene_subset, 
                        output="overlapping", bindingRegion=c(-1000,1000))
peaks_in_1kb_subset <- 
    peaks_in_1kb_subset[which(!is.na(peaks_in_1kb_subset$feature))]
  
peaks_in_100kb_subset <- 
    annotatePeakInBatch(macsOutput, 
                        PeakLocForDistance = "middle",
                        FeatureLocForDistance = "TSS",
                        AnnotationData=ens90.refGene_subset, 
                        output="overlapping", bindingRegion=c(-100000,100000))
peaks_in_100kb_subset <- 
    peaks_in_100kb_subset[which(!is.na(peaks_in_100kb_subset$feature))]

mcols(peaks_in_1kb_subset) <- DataFrame(ens_gene=mcols(peaks_in_1kb_subset)$feature)
mcols(peaks_overlapping_tss_subset) <- 
  DataFrame(ens_gene=mcols( peaks_overlapping_tss_subset)$feature)
mcols(peaks_in_100kb_subset) <- DataFrame(ens_gene=mcols(peaks_in_100kb_subset)$feature)

## Two options: filter exons or not
for(exonfilter in c(TRUE, FALSE)) {
  if(exonfilter) {
    
    proximal_peaks_subset <- anti_join(data.frame(peaks_in_1kb_subset),
                                       data.frame(peaks_overlapping_exons), 
                                       by = c("seqnames", "start", "end", "width", "strand")) %>%
      bind_rows(., data.frame(peaks_overlapping_tss_subset)) %>%
      dplyr::select(seqnames, start, end, ens_gene) %>%
      unique() %>%
      arrange(seqnames, start, end)
    
    distal_peaks_subset <- anti_join(data.frame(peaks_in_100kb_subset),
                                       data.frame(proximal_peaks_subset, width=500, 
                                                  strand="*"), 
                                       by = c("seqnames", "start", "end", "width", "strand", "ens_gene")) %>%
      anti_join(., data.frame(peaks_overlapping_exons), by = c("seqnames", "start", "end", "width", "strand")) %>%
      dplyr::select(seqnames, start, end, ens_gene) %>%
      unique() %>%
      arrange(seqnames, start, end)
      
  } else {
    proximal_peaks_subset <- 
      bind_rows(data.frame(peaks_in_1kb_subset), data.frame(peaks_overlapping_tss_subset)) %>%
      dplyr::select(seqnames, start, end, ens_gene) %>%
      unique() %>%
      arrange(seqnames, start, end)
    
    distal_peaks_subset <- anti_join(data.frame(peaks_in_100kb_subset),
                                     data.frame(proximal_peaks_subset, width=500, 
                                                strand="*"), 
                                     by = c("seqnames", "start", "end", "width", "strand", "ens_gene")) %>%
      dplyr::select(seqnames, start, end, ens_gene) %>%
      unique() %>%
      arrange(seqnames, start, end)
  }
  
  proximal_peaks_df <- proximal_peaks_subset %>%
    as.data.frame(., row.names=NULL) %>% na.omit() %>% dplyr::select(chr=seqnames, start, end, ens_gene) %>% 
    mutate(start = as.numeric(start), end = as.numeric(end), 
           ens_gene = as.character(ens_gene),
           type = "proximal") %>% arrange(chr, start, end) %>% unique()
  distal_peaks_df <- distal_peaks_subset %>%
    as.data.frame(., row.names=NULL) %>% na.omit() %>% dplyr::select(chr=seqnames, start, end, ens_gene) %>% 
    mutate(start = as.numeric(start), end = as.numeric(end), 
           ens_gene = as.character(ens_gene),
           type = "distal") %>% arrange(chr, start, end) %>% unique()
  
  all_proximal_distal <- bind_rows(proximal_peaks_df, distal_peaks_df) %>%
    arrange(chr, start, end)
  saveRDS(all_proximal_distal, paste0("data/all_proximal_distal", ifelse(exonfilter, "_exonfilter", "_noexonfilter"), ".rds"))
  write.table(all_proximal_distal, paste0("data/all_proximal_distal", ifelse(exonfilter, "_exonfilter", "_noexonfilter"), ".txt"),
              row.names=FALSE)
}

#------------------------------------------------------------------------------------------------------------------
## Make SQLite databases
#------------------------------------------------------------------------------------------------------------------

## 1. peak_sequences_df.sqlite
peak_sequences_df <- fread("data/peak_sequences_df.txt")
mydb <- dbConnect(RSQLite::SQLite(), "data/rosetta.sqlite")
dbWriteTable(mydb, name="peak_sequences_df", value=peak_sequences_df, overwrite=TRUE, append=FALSE)
dbDisconnect(mydb)


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


## Read in ATAC data
all_proximal_distal_exonfilter <- readRDS("data/all_proximal_distal_exonfilter.rds")
all_proximal_distal_noexonfilter <- readRDS("data/all_proximal_distal_noexonfilter.rds")




#------------------------------------------------------------------------------------------------------------------
## Save pre-named cluster gene IDs
#------------------------------------------------------------------------------------------------------------------

prenamed_clusters <- read_excel("data/DEgenes_heatmap_transcript IDs.xls",
                                sheet="heatmap_clusters") %>%
  dplyr::select(ID, cluster_name = `cluster name`)
saveRDS(prenamed_clusters, "data/prenamed_clusters.rds")


#------------------------------------------------------------------------------------------------------------------
## Save additional pre-named clusters
#------------------------------------------------------------------------------------------------------------------

library(biomaRt)
ensembl <- useMart(host = 'mar2016.archive.ensembl.org', 
                   biomart='ENSEMBL_MART_ENSEMBL')
ds <- readRDS("data/ensembl_ds.rds")
human_ids <- readRDS("data/human_ids.rds")

LFC_all <- readRDS("data/heatmap_LFC.rds")
logtpm <- readRDS("data/heatmap_logtpm.rds")
LFC_all <- LFC_all[,c(2,1,3,4)]
colnames(LFC_all) <- strsplit(colnames(LFC_all), split = ".b", fixed= TRUE) %>% unlist()
LFC_all <- LFC_all[match(rownames(logtpm), rownames(LFC_all)),]
available_transcripts <- rownames(LFC_all)
transcript_to_gene <- readRDS("data/transcripts_to_genes_map.rds")

expressed <-  rownames(LFC_all)[which(rowSums(is.na(LFC_all)) < 4)]

final_list <- vector("list", 13)
names(final_list) <- c("Belin", "Chandran", "Crippa_mouse", "Crippa_zebrafish", "Dwaraka", "Herman_brain",
                       #"Herman_brain_TF-encoding", 
                       "Herman_SC", 
                       #"Herman_SC_TF-encoding", 
                       "Li",
                       "MGI_cholesterol_metabolic", "Rabinowitz", "Sauls", "Sifuentes", "Umansky")
for(i in names(final_list)) {
  tmmp <- read_excel("data/rosetta_extra_data.xlsx", sheet = i)
  cat(i, ": ", nrow(tmmp), "inputs \n")
  search <- unlist(unique(tmmp[,1]))
  cat(length(search), "unique input genes \n")

  ## Mouse Ensembl gene IDs
  if(i %in% c("Belin", "Chandran", "Crippa_mouse", "MGI_cholesterol_metabolic", "Sauls", "Umansky")) {
    neworg <- useDataset(ds[which(ds$dataset == "mmusculus_gene_ensembl"),1], mart=ensembl)
    zebrafish <- useDataset("drerio_gene_ensembl", mart=ensembl)
    new_gene_names <- getLDS(attributes=c("ensembl_gene_id"), filters="ensembl_gene_id",
                             values=search,
                             mart=neworg, attributesL=c("ensembl_gene_id"), martL=zebrafish)[,2] 
    tmp <- suppressWarnings(left_join(data.frame(ens_gene = new_gene_names),
                       transcript_to_gene, by = "ens_gene"))
    final <- unique(na.omit(tmp$target_id))
    final_list[[i]] <- data.frame(ID = final[which(final %in% expressed)], stringsAsFactors = FALSE)
    cat(nrow(final_list[[i]]), "unique expressed zebrafish transcripts \n")
    
  ## Zebrafish gene IDs  
  } else if(i %in% c("Crippa_zebrafish", "Rabinowitz", "Sifuentes")) {
    tmp <- suppressWarnings(left_join(data.frame(ens_gene = search),
                                      transcript_to_gene, by = "ens_gene"))
    final <- unique(na.omit(tmp$target_id))
    final_list[[i]] <- data.frame(ID = final[which(final %in% expressed)], stringsAsFactors = FALSE)
    cat(nrow(final_list[[i]]), "unique expressed zebrafish transcripts \n")
    
  ## HGNC symbols  
  } else if(i %in% c("Dwaraka", "Herman_brain", "Herman_brain_TF-encoding", 
                     "Herman_SC", "Herman_SC_TF-encoding")){
    search_tmp <- paste0(unlist(lapply(search, function(x) paste0("^", x, "$", collapse=""))), 
                         collapse = "|") 
    find_ensgene <- human_ids %>%
      dplyr::filter(grepl(search_tmp, hgnc_symbol, ignore.case = TRUE)) %>%
      dplyr::select(ens_gene) %>% unlist
    find_ensdart <- unique(transcript_to_gene[which(transcript_to_gene$ens_gene %in% find_ensgene),1])
    final <- unique(na.omit(find_ensdart))
    final_list[[i]] <- data.frame(ID = final[which(final %in% expressed)], stringsAsFactors = FALSE)
    cat(nrow(final_list[[i]]), "unique expressed zebrafish transcripts \n")  
    
  ## Rat Ensembl gene IDs  
  } else {
    neworg <- useDataset(ds[which(ds$dataset == "rnorvegicus_gene_ensembl"),1], mart=ensembl)
    zebrafish <- useDataset("drerio_gene_ensembl", mart=ensembl)
    new_gene_names <- getLDS(attributes=c("ensembl_gene_id"), filters="ensembl_gene_id",
                             values=search,
                             mart=neworg, attributesL=c("ensembl_gene_id"), martL=zebrafish)[,2] 
    tmp <- suppressWarnings(left_join(data.frame(ens_gene = new_gene_names),
                                      transcript_to_gene, by = "ens_gene"))
    final <- unique(na.omit(tmp$target_id))
    final_list[[i]] <- data.frame(ID = final[which(final %in% expressed)], stringsAsFactors = FALSE)
    cat(nrow(final_list[[i]]), "unique expressed zebrafish transcripts \n")
  }
  cat("*****\n")
}


names(final_list) <- c("(Belin 2015) Mouse ONI Intact v Axotomized RGC",
"(Chandran 2016) Mouse PNS regeneration",
"(Crippa 2016) Mouse Heart regeneration",
"(Crippa 2016) Zebrafish Heart regeneration",
"(Dwaraka 2018) Salamander Limb regeneration",
"(Herman 2018) Lamprey SCI Brain response",
"(Herman 2018) Lamprey SCI SC response",
"(Li 2013) Rat SNI proximal segment response",
"MGI cholesterol metabolic",
"(Rabinowitz 2017) Zebrafish Fin regeneration",
"(Sauls 2018) Mouse embryonic fibroblast cardiac reprogramming",
"(Sifuentes 2016) Zebrafish acute photic lesion Muller glia response",
"(Umansky 2015) Mouse Muscle regeneration")

final_list_df <- bind_rows(final_list, .id = "cluster_name")

final_prenamed_clusters <- bind_rows(prenamed_clusters, final_list_df)
saveRDS(final_prenamed_clusters, "data/prenamed_clusters.rds")

# Belin :  62 inputs 
# 61 unique input genes 
# 178 unique expressed zebrafish transcripts 
# *****
#   Chandran :  62 inputs 
# 60 unique input genes 
# 137 unique expressed zebrafish transcripts 
# *****
#   Crippa_mouse :  6592 inputs 
# 5580 unique input genes 
# 10199 unique expressed zebrafish transcripts 
# *****
#   Crippa_zebrafish :  891 inputs 
# 485 unique input genes 
# 878 unique expressed zebrafish transcripts 
# *****
#   Dwaraka :  405 inputs 
# 405 unique input genes 
# 761 unique expressed zebrafish transcripts 
# *****
#   Herman_brain :  3664 inputs 
# 2325 unique input genes 
# 3712 unique expressed zebrafish transcripts 
# *****
#   Herman_SC :  3998 inputs 
# 2519 unique input genes 
# 4076 unique expressed zebrafish transcripts 
# *****
#   Li :  456 inputs 
# 445 unique input genes 
# 683 unique expressed zebrafish transcripts 
# *****
#   MGI_cholesterol_metabolic :  125 inputs 
# 125 unique input genes 
# 232 unique expressed zebrafish transcripts 
# *****
#   Rabinowitz :  566 inputs 
# 566 unique input genes 
# 927 unique expressed zebrafish transcripts 
# *****
#   Sauls :  3463 inputs 
# 3382 unique input genes 
# 7616 unique expressed zebrafish transcripts 
# *****
#   Sifuentes :  3143 inputs 
# 2685 unique input genes 
# 4303 unique expressed zebrafish transcripts 
# *****
#   Umansky :  636 inputs 
# 623 unique input genes 
# 1055 unique expressed zebrafish transcripts 
# *****
