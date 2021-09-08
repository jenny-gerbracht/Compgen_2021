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

broad_gex <- read.table(file = "prepared/GDSC_broad/all_intersect_genes/gex.tsv.gz",
                        sep = ",",
                        header = TRUE)
head(broad_gex)
#Sample names as rownames
rownames(broad_gex) <- broad_gex$V1
broad_gex$V1 <- NULL

#Take top 1000 most variable predictors
SDs <- apply(broad_gex, 2, sd)
topPreds <- order(SDs, decreasing = TRUE)[1:1000]
broad_gex_top <- broad_gex[,topPreds]

#Center and scale
preproc <- preProcess(broad_gex_top, method = c("center", "scale"))
broad_gex_top <- predict(preproc, broad_gex_top)

#Filter variables that are highly correlated
corrFilt <- preProcess(broad_gex_top, method = "corr", cutoff = 0.9)
broad_gex_top  <- predict(corrFilt, broad_gex_top )

#Check for NA values
anyNA(broad_gex_top)

#####
#CNV
#####

broad_cnv <- read.delim(file = "prepared/GDSC_broad/all_intersect_genes/cnv.tsv.gz",
                        sep = ",",
                        header = TRUE)
rownames(broad_cnv) <- broad_cnv$V1
broad_cnv$V1 <- NULL

descrCor <-  cor(broad_cnv)
highlyCorDescr <- findCorrelation(descrCor, cutoff = .9)
broad_cnv_cor <- broad_cnv[, -highlyCorDescr]

#####
#mutation
#####

broad_mut <- read.table(file = "prepared/GDSC_broad/all_intersect_genes/mut.tsv.gz",
                        sep = ",",
                        header = TRUE)

#Count number of non-silent mutations per gene and tumor line

broad_mut %>%
  select(Tumor_Sample_Barcode, gene_id) %>% 
  dplyr::count(Tumor_Sample_Barcode, gene_id) %>% 
  pivot_wider(id_cols = Tumor_Sample_Barcode, names_from = gene_id, values_from = n) %>% 
  distinct() -> broad_mut_select_counts_wide

broad_mut_select_counts_wide[is.na(broad_mut_select_counts_wide)] <- 0
broad_mut_select_counts_wide <- as.data.frame(broad_mut_select_counts_wide)
rownames(broad_mut_select_counts_wide) <- broad_mut_select_counts_wide$Tumor_Sample_Barcode
broad_mut_select_counts_wide$Tumor_Sample_Barcode <- NULL
broad_mut_bin <- as.matrix((broad_mut_select_counts_wide > 0)+0)
rownames(broad_mut_bin) <- rownames(broad_mut_select_counts_wide)

#Only keep genes that are mutated in < 10% and > 90% of the samples
mut_filter <- apply(broad_mut_bin, 2, mean)
mut_filter_cols <- mut_filter[mut_filter > 0.1 & mut_filter < 0.9]
broad_mut_bin_nzv <- broad_mut_bin[,names(mut_filter_cols)]

#Rename columns
colnames(broad_gex_top) <- paste(colnames(broad_gex_top), "gex", sep = "_")
colnames(broad_cnv_cor) <- paste(colnames(broad_cnv_cor), "cnv", sep = "_")
colnames(broad_mut_bin_nzv) <- paste(colnames(broad_mut_bin_nzv), "mut", sep = "_")

broad_gex_top <- rownames_to_column(broad_gex_top, var = "sample")
broad_cnv_cor <- rownames_to_column(broad_cnv_cor, var = "sample")
broad_mut_bin_nzv <- rownames_to_column(as.data.frame(broad_mut_bin_nzv), var = "sample")

broad_merged_matrix <- inner_join(broad_gex_top, broad_cnv_cor, by = "sample")
broad_merged_matrix <- inner_join(broad_merged_matrix, broad_mut_bin_nzv, by = "sample")

#rm(list= ls()[!(ls() %in% c("broad_gex_top", "broad_cnv_cor", "broad_mut_bin_nzv", "broad_merged_matrix"))])

###
#drug response
###

filenames <- list.files(path = "prepared/GDSC_broad/drug_response/", pattern = ".tsv.gz", full.names = TRUE)
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
broad_merged_matrix_dr <- inner_join(combined_dr_wide, broad_merged_matrix, by = "sample")

