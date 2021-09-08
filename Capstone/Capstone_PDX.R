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

PDX_gex <- read.delim(file = "prepared/PDX_all/all_intersect_genes/gex.tsv.gz")

#Sample names as rownames
rownames(PDX_gex) <- PDX_gex$X
PDX_gex$X <- NULL

#Take top 1000 most variable predictors
SDs <- apply(PDX_gex, 2, sd)
topPreds <- order(SDs, decreasing = TRUE)[1:1000]
PDX_gex_top <- PDX_gex[,topPreds]

#Center and scale
preproc <- preProcess(PDX_gex_top, method = c("center", "scale"))
PDX_gex_top <- predict(preproc, PDX_gex_top)

#Filter variables that are highly correlated
corrFilt <- preProcess(PDX_gex_top, method = "corr", cutoff = 0.9)
PDX_gex_top  <- predict(corrFilt, PDX_gex_top )

#Check for NA values
anyNA(PDX_gex_top)

#####
#CNV
#####

PDX_cnv <- read.delim(file = "prepared/PDX_all/all_intersect_genes/cnv.tsv.gz")
rownames(PDX_cnv) <- PDX_cnv$X
PDX_cnv$X <- NULL

descrCor <-  cor(PDX_cnv)
highlyCorDescr <- findCorrelation(descrCor, cutoff = .9)
PDX_cnv_cor <- PDX_cnv[, -highlyCorDescr]

#####
#mutation
#####

PDX_mut <- read.delim(file = "prepared/PDX_all/all_intersect_genes/mut.tsv.gz")

#Count number of non-silent mutations per gene and tumor line

PDX_mut %>%
  select(Sample, Gene) %>% 
  dplyr::count(Sample, Gene) %>% 
  pivot_wider(id_cols = Sample, names_from = Gene, values_from = n) %>% 
  distinct() -> PDX_mut_select_counts_wide

PDX_mut_select_counts_wide[is.na(PDX_mut_select_counts_wide)] <- 0
PDX_mut_select_counts_wide <- as.data.frame(PDX_mut_select_counts_wide)
rownames(PDX_mut_select_counts_wide) <- PDX_mut_select_counts_wide$Sample
PDX_mut_select_counts_wide$Sample <- NULL
PDX_mut_bin <- as.matrix((PDX_mut_select_counts_wide > 0)+0)
rownames(PDX_mut_bin) <- rownames(PDX_mut_select_counts_wide)

#Only keep genes that are mutated in < 10% and > 90% of the samples
mut_filter <- apply(PDX_mut_bin, 2, mean)
mut_filter_cols <- mut_filter[mut_filter > 0.1 & mut_filter < 0.9]
PDX_mut_bin_nzv <- PDX_mut_bin[,names(mut_filter_cols)]

#Rename columns
colnames(PDX_gex_top) <- paste(colnames(PDX_gex_top), "gex", sep = "_")
colnames(PDX_cnv_cor) <- paste(colnames(PDX_cnv_cor), "cnv", sep = "_")
colnames(PDX_mut_bin_nzv) <- paste(colnames(PDX_mut_bin_nzv), "mut", sep = "_")

PDX_gex_top <- rownames_to_column(PDX_gex_top, var = "sample")
PDX_cnv_cor <- rownames_to_column(PDX_cnv_cor, var = "sample")
PDX_mut_bin_nzv <- rownames_to_column(as.data.frame(PDX_mut_bin_nzv), var = "sample")

PDX_merged_matrix <- inner_join(PDX_gex_top, PDX_cnv_cor, by = "sample")
PDX_merged_matrix <- inner_join(PDX_merged_matrix, PDX_mut_bin_nzv, by = "sample")

#rm(list= ls()[!(ls() %in% c("PDX_gex_top", "PDX_cnv_cor", "PDX_mut_bin_nzv", "merged_matrix"))])

###
#drug response
###

filenames <- list.files(path = "prepared/PDX_all/drug_response/", pattern = ".tsv.gz", full.names = TRUE)
ldf <- lapply(filenames, read.delim)
combined_dr <- do.call("rbind", ldf)
combined_dr %>%
  select(sample_id, column_name, value) %>% 
  pivot_wider(id_cols = sample_id,
              names_from = column_name,
              values_from = value) -> combined_dr_wide
combined_dr_wide <- dplyr::rename(combined_dr_wide, "sample" = "sample_id")
PDX_merged_matrix_dr <- inner_join(combined_dr_wide, PDX_merged_matrix, by = "sample")

#data

PDX_data <- list()
drugs <- unique(combined_dr$column_name)

for (i in drugs){
  print(i)
  PDX_merged_matrix_dr %>% 
    drop_na(i) %>% 
    select(i, 8:ncol(PDX_merged_matrix_dr)) -> PDX_data[[i]]
}

