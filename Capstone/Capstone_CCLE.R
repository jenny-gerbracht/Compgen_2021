library(caret)
library(dplyr)
library(tidyr)
library(tibble)
library(ranger)

memory.limit(9999999999)

setwd("X:/Workshops/Compgen/Capstone")

#####
#gene expression
#####

CCLE_gex <- read.delim(file = "prepared/CCLE/all_intersect_genes/gex.tsv.gz")

#Sample names as rownames
rownames(CCLE_gex) <- CCLE_gex$X
CCLE_gex$X <- NULL

#Take top 1000 most variable predictors
SDs <- apply(CCLE_gex, 2, sd)
topPreds <- order(SDs, decreasing = TRUE)[1:1000]
CCLE_gex_top <- CCLE_gex[,topPreds]

#Center and scale
preproc <- preProcess(CCLE_gex_top, method = c("center", "scale"))
CCLE_gex_top <- predict(preproc, CCLE_gex_top)

#Filter variables that are highly correlated
corrFilt <- preProcess(CCLE_gex_top, method = "corr", cutoff = 0.9)
CCLE_gex_top  <- predict(corrFilt, CCLE_gex_top )

#Check for NA values
anyNA(CCLE_gex_top)

#####
#CNV
#####

CCLE_cnv <- read.delim(file = "prepared/CCLE/all_intersect_genes/cnv.tsv.gz")
rownames(CCLE_cnv) <- CCLE_cnv$X
CCLE_cnv$X <- NULL

descrCor <-  cor(CCLE_cnv)
highlyCorDescr <- findCorrelation(descrCor, cutoff = .9)
CCLE_cnv_cor <- CCLE_cnv[, -highlyCorDescr]

#####
#mutation
#####

CCLE_mut <- read.delim(file = "prepared/CCLE/all_intersect_genes/mut.tsv.gz")

#Count number of non-silent mutations per gene and tumor line

CCLE_mut %>%
  select(Tumor_Sample_Barcode, gene_id) %>% 
  dplyr::count(Tumor_Sample_Barcode, gene_id) %>% 
  pivot_wider(id_cols = Tumor_Sample_Barcode, names_from = gene_id, values_from = n) %>% 
  distinct() -> CCLE_mut_select_counts_wide

CCLE_mut_select_counts_wide[is.na(CCLE_mut_select_counts_wide)] <- 0
CCLE_mut_select_counts_wide <- as.data.frame(CCLE_mut_select_counts_wide)
rownames(CCLE_mut_select_counts_wide) <- CCLE_mut_select_counts_wide$Tumor_Sample_Barcode
CCLE_mut_select_counts_wide$Tumor_Sample_Barcode <- NULL
CCLE_mut_bin <- as.matrix((CCLE_mut_select_counts_wide > 0)+0)
rownames(CCLE_mut_bin) <- rownames(CCLE_mut_select_counts_wide)

#Only keep genes that are mutated in < 10% and > 90% of the samples
mut_filter <- apply(CCLE_mut_bin, 2, mean)
mut_filter_cols <- mut_filter[mut_filter > 0.1 & mut_filter < 0.9]
CCLE_mut_bin_nzv <- CCLE_mut_bin[,names(mut_filter_cols)]

#Rename columns
colnames(CCLE_gex_top) <- paste(colnames(CCLE_gex_top), "gex", sep = "_")
colnames(CCLE_cnv_cor) <- paste(colnames(CCLE_cnv_cor), "cnv", sep = "_")
colnames(CCLE_mut_bin_nzv) <- paste(colnames(CCLE_mut_bin_nzv), "mut", sep = "_")

CCLE_gex_top <- rownames_to_column(CCLE_gex_top, var = "sample")
CCLE_cnv_cor <- rownames_to_column(CCLE_cnv_cor, var = "sample")
CCLE_mut_bin_nzv <- rownames_to_column(as.data.frame(CCLE_mut_bin_nzv), var = "sample")

merged_matrix <- inner_join(CCLE_gex_top, CCLE_cnv_cor, by = "sample")
merged_matrix <- inner_join(merged_matrix, CCLE_mut_bin_nzv, by = "sample")

#rm(list= ls()[!(ls() %in% c("CCLE_gex_top", "CCLE_cnv_cor", "CCLE_mut_bin_nzv", "merged_matrix"))])

###
#drug response
###

filenames <- list.files(path = "prepared/CCLE/drug_response/", pattern = ".tsv.gz", full.names = TRUE)
ldf <- lapply(filenames, read.delim)
combined_dr <- do.call("rbind", ldf)
combined_dr %>%
  select(sample_id, column_name, value) %>% 
  pivot_wider(id_cols = sample_id,
              names_from = column_name,
              values_from = value) -> combined_dr_wide
combined_dr_wide <- dplyr::rename(combined_dr_wide, "sample" = "sample_id")
merged_matrix_dr <- inner_join(combined_dr_wide, merged_matrix, by = "sample")

