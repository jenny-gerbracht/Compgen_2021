library(caret)
library(dplyr)
library(tidyr)
library(tibble)
library(ranger)

memory.limit(9999999999)

setwd("X:/Workshops/Compgen/Capstone/")

#####
#gene expression
#####

sanger_gex <- read.table(file = "prepared/GDSC_sanger/all_intersect_genes/gex.tsv.gz",
                        sep = ",",
                        header = TRUE)

#Sample names as rownames
rownames(sanger_gex) <- sanger_gex$V1
sanger_gex$V1 <- NULL

#Take top 1000 most variable predictors
SDs <- apply(sanger_gex, 2, sd)
topPreds <- order(SDs, decreasing = TRUE)[1:1000]
sanger_gex_top <- sanger_gex[,topPreds]

#Center and scale
preproc <- preProcess(sanger_gex_top, method = c("center", "scale"))
sanger_gex_top <- predict(preproc, sanger_gex_top)

#Filter variables that are highly correlated
corrFilt <- preProcess(sanger_gex_top, method = "corr", cutoff = 0.9)
sanger_gex_top  <- predict(corrFilt, sanger_gex_top )

#Check for NA values
anyNA(sanger_gex_top)

#####
#CNV
#####

sanger_cnv <- read.delim(file = "prepared/GDSC_sanger/all_intersect_genes/cnv.tsv.gz",
                        sep = ",",
                        header = TRUE)
rownames(sanger_cnv) <- sanger_cnv$V1
sanger_cnv$V1 <- NULL

descrCor <-  cor(sanger_cnv)
highlyCorDescr <- findCorrelation(descrCor, cutoff = .9)
sanger_cnv_cor <- sanger_cnv[, -highlyCorDescr]

#####
#mutation
#####

sanger_mut <- read.table(file = "prepared/GDSC_sanger/all_intersect_genes/mut.tsv.gz",
                        sep = ",",
                        header = TRUE)

#Count number of non-silent mutations per gene and tumor line

sanger_mut %>%
  select(Tumor_Sample_Barcode, gene_id) %>% 
  dplyr::count(Tumor_Sample_Barcode, gene_id) %>% 
  pivot_wider(id_cols = Tumor_Sample_Barcode, names_from = gene_id, values_from = n) %>% 
  distinct() -> sanger_mut_select_counts_wide

sanger_mut_select_counts_wide[is.na(sanger_mut_select_counts_wide)] <- 0
sanger_mut_select_counts_wide <- as.data.frame(sanger_mut_select_counts_wide)
rownames(sanger_mut_select_counts_wide) <- sanger_mut_select_counts_wide$Tumor_Sample_Barcode
sanger_mut_select_counts_wide$Tumor_Sample_Barcode <- NULL
sanger_mut_bin <- as.matrix((sanger_mut_select_counts_wide > 0)+0)
rownames(sanger_mut_bin) <- rownames(sanger_mut_select_counts_wide)

#Only keep genes that are mutated in < 10% and > 90% of the samples
mut_filter <- apply(sanger_mut_bin, 2, mean)
mut_filter_cols <- mut_filter[mut_filter > 0.1 & mut_filter < 0.9]
sanger_mut_bin_nzv <- sanger_mut_bin[,names(mut_filter_cols)]

#Rename columns
colnames(sanger_gex_top) <- paste(colnames(sanger_gex_top), "gex", sep = "_")
colnames(sanger_cnv_cor) <- paste(colnames(sanger_cnv_cor), "cnv", sep = "_")
colnames(sanger_mut_bin_nzv) <- paste(colnames(sanger_mut_bin_nzv), "mut", sep = "_")

sanger_gex_top <- rownames_to_column(sanger_gex_top, var = "sample")
sanger_cnv_cor <- rownames_to_column(sanger_cnv_cor, var = "sample")
sanger_mut_bin_nzv <- rownames_to_column(as.data.frame(sanger_mut_bin_nzv), var = "sample")

sanger_merged_matrix <- inner_join(sanger_gex_top, sanger_cnv_cor, by = "sample")
sanger_merged_matrix <- inner_join(sanger_merged_matrix, sanger_mut_bin_nzv, by = "sample")

#rm(list= ls()[!(ls() %in% c("sanger_gex_top", "sanger_cnv_cor", "sanger_mut_bin_nzv", "sanger_merged_matrix"))])

###
#drug response
###

filenames <- list.files(path = "prepared/GDSC_sanger/drug_response/", pattern = ".tsv.gz", full.names = TRUE)
ldf <- lapply(filenames, read.delim)
combined_dr <- do.call("rbind", ldf)
combined_dr %>%
  select(sample_id, column_name, value) %>% 
  dplyr::group_by(sample_id, column_name) %>% 
  mutate(value_mean = mean(value)) %>% 
  select(!value) %>% 
  distinct() %>% 
  pivot_wider(id_cols = sample_id,
              names_from = column_name,
              values_from = value_mean) -> combined_dr_wide
combined_dr_wide <- dplyr::rename(combined_dr_wide, "sample" = "sample_id")
sanger_merged_matrix_dr <- inner_join(combined_dr_wide, sanger_merged_matrix, by = "sample")

