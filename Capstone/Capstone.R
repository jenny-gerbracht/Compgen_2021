library(caretEnsemble)
library(dplyr)

###
#models
###

#add pdx
merged_matrix_con <- (merged_matrix_dr[,colnames(merged_matrix_dr) %in% colnames(broad_merged_matrix_dr) &
                        colnames(merged_matrix_dr) %in% colnames(sanger_merged_matrix_dr) &
                          colnames(merged_matrix_dr) %in% colnames(PDX_merged_matrix_dr)])

broad_merged_matrix_con <- (broad_merged_matrix_dr[,colnames(broad_merged_matrix_dr) %in% colnames(merged_matrix_dr) &
                        colnames(broad_merged_matrix_dr) %in% colnames(sanger_merged_matrix_dr) &
                          colnames(broad_merged_matrix_dr) %in% colnames(PDX_merged_matrix_dr)])

sanger_merged_matrix_con <- (sanger_merged_matrix_dr[,colnames(sanger_merged_matrix_dr) %in% colnames(merged_matrix_dr) &
                              colnames(sanger_merged_matrix_dr) %in% colnames(broad_merged_matrix_dr)&
                                colnames(sanger_merged_matrix_dr) %in% colnames(PDX_merged_matrix_dr)])

PDX_merged_matrix_con <- (PDX_merged_matrix_dr[,colnames(PDX_merged_matrix_dr) %in% colnames(merged_matrix_dr) &
                                                       colnames(PDX_merged_matrix_dr) %in% colnames(broad_merged_matrix_dr)&
                                                       colnames(PDX_merged_matrix_dr) %in% colnames(sanger_merged_matrix_dr)])

###
#models
###
#CCLE
###
#train model
###
data <- list()
drugs <- unique(combined_dr$column_name)
formula <- list()
for (i in c(1:6)){
  formula[[i]] <- formula(paste(drugs[i], "~."))
}

for (i in drugs){
  print(i)
  merged_matrix_con %>% 
    drop_na(i) %>% 
    select(i, 8:ncol(merged_matrix_con)) -> data[[i]]
}
trctrl <- trainControl(method = "cv", 
                       number = 10)
models <- list()
for(i in c(1:6)){
  models[[i]] <- train(formula[[i]], 
                       data = data[[i]], 
                       method = "ranger",
                       trControl = trctrl)
}
###
#broad
###
broad_data <- list()
for (i in drugs){
  print(i)
  broad_merged_matrix_con %>%
    ungroup() %>% 
    drop_na(i) %>% 
    select(i, 8:ncol(broad_merged_matrix_con)) -> broad_data[[i]]
}
broad_models <- list()
for(i in c(1:6)){
  broad_models[[i]] <- train(formula[[i]], 
                       data = broad_data[[i]], 
                       method = "ranger",
                       trControl = trctrl)
}
###
#sanger
###
sanger_data <- list()
for (i in drugs){
  sanger_merged_matrix_con %>%
    ungroup() %>% 
    drop_na(i) %>% 
    select(i, 8:ncol(sanger_merged_matrix_con)) -> sanger_data[[i]]
}
sanger_models <- list()
for(i in c(1:6)){
  sanger_models[[i]] <- train(formula[[i]], 
                             data = sanger_data[[i]], 
                             method = "ranger",
                             trControl = trctrl)
}
###
#PDX
###
PDX_data <- list()
for (i in drugs){
  PDX_merged_matrix_con %>% 
    drop_na(i) %>% 
    select(i, 8:ncol(PDX_merged_matrix_con)) -> PDX_data[[i]]
}

###
#predictions
###

#CCLE data with CCLE models
predictions_CCLE_CCLE <- list()
for (i in c(1:6)){
  predict <- predict(models[[i]], 
                     data[[i]][,-1])
  predictions_CCLE_CCLE[[i]] <- cor(predict, data[[i]][,1])^2
}

#CCLE data with broad models
predictions_CCLE_broad <- list()
for (i in c(1:6)){
  predict <- predict(broad_models[[i]], 
                     data[[i]][,-1])
  predictions_CCLE_broad[[i]] <- cor(predict, data[[i]][,1])^2
}

#CCLE data with sanger models
predictions_CCLE_sanger <- list()
for (i in c(1:6)){
  predict <- predict(sanger_models[[i]], 
                     data[[i]][,-1])
  predictions_CCLE_sanger[[i]] <- cor(predict, data[[i]][,1])^2
}

#broad data with CCLE model
predictions_broad_CCLE <- list()
for (i in c(1:6)){
  predict <- predict(models[[i]], 
                     broad_data[[i]][,-1])
  predictions_broad_CCLE[[i]] <- cor(predict, broad_data[[i]][,1])^2
}

#broad data with broad model
predictions_broad_broad <- list()
for (i in c(1:6)){
  predict <- predict(broad_models[[i]], 
                     broad_data[[i]][,-1])
  predictions_broad_broad[[i]] <- cor(predict, broad_data[[i]][,1])^2
}

#broad data with sanger model
predictions_broad_sanger <- list()
for (i in c(1:6)){
  predict <- predict(sanger_models[[i]], 
                     broad_data[[i]][,-1])
  predictions_broad_sanger[[i]] <- cor(predict, broad_data[[i]][,1])^2
}

#sanger data with CCLE model
predictions_sanger_CCLE <- list()
for (i in c(1:6)){
  predict <- predict(models[[i]], 
                     sanger_data[[i]][,-1])
  predictions_sanger_CCLE[[i]] <- cor(predict, sanger_data[[i]][,1])^2
}

#sanger data with broad model
predictions_sanger_broad <- list()
for (i in c(1:6)){
  predict <- predict(broad_models[[i]], 
                     sanger_data[[i]][,-1])
  predictions_sanger_broad[[i]] <- cor(predict, sanger_data[[i]][,1])^2
}

#sanger data with sanger model
predictions_sanger_sanger <- list()
for (i in c(1:6)){
  predict <- predict(sanger_models[[i]], 
                     sanger_data[[i]][,-1])
  predictions_sanger_sanger[[i]] <- cor(predict, sanger_data[[i]][,1])^2
}

task1 <- data.frame(CCLE = unlist(predictions_CCLE_CCLE),
                    GDSC_Broad = unlist(predictions_broad_broad),
                    GDSC_Sanger = unlist(predictions_sanger_sanger))
rownames(task1) <- drugs

task1_cross <- data.frame(CCLE_CCLE = unlist(predictions_CCLE_CCLE),
                          CCLE_Broad  = unlist(predictions_CCLE_broad),
                          CCLE_Sanger = unlist(predictions_CCLE_sanger),
                          Broad_CCLE = unlist(predictions_broad_CCLE),
                          Broad_Broad = unlist(predictions_broad_broad),
                          Broad_Sanger = unlist(predictions_broad_sanger),
                          Sanger_CCLE = unlist(predictions_sanger_CCLE),
                          Sanger_Broad = unlist(predictions_sanger_broad),
                          Sanger_Sanger = unlist(predictions_sanger_sanger))
rownames(task1_cross) <- drugs

write.table(task1,
            file = "task1.txt",
            quote = FALSE,
            sep = "\t")

write.table(task1_cross,
            file = "task1_cross.txt",
            quote = FALSE,
            sep = "\t")


task1_cross_mean <- apply(task1_cross, 1, mean, na.rm = TRUE)
task1_cross_median <- apply(task1_cross, 1, median, na.rm = TRUE)

task1_cross_summ <- cbind(task1_cross_mean, task1_cross_median)
colnames(task1_cross_summ) <- c("mean", "median")

boxplot(t(task1_cross))

#mean R2 per drug
for (i in c(1:6)){
  print(apply((task1_cross[i, c(4,7)]), 1, mean, na.rm = TRUE))
  print(apply((task1_cross[i, c(2,8)]), 1, mean, na.rm = TRUE))
  print(apply((task1_cross[i, c(3,6)]), 1, mean, na.rm = TRUE))
}

#predict PDX set
predictions_PDX <- list()
for (i in c(1:6)){
  predict <- predict(models[[i]], 
                     PDX_data[[i]][,-1])
  predictions_PDX[[i]] <- cor(predict, PDX_data[[i]][,1])^2
}

task1_PDX <- unlist(predictions_PDX)

write.table(task1_PDX,
            file = "task1_PDX.txt",
            quote = FALSE,
            sep = "\t")

CCLE_gex_test <- CCLE_gex_top
PDX_gex_test <- PDX_gex_top

CCLE_gex_test <- CCLE_gex_test[,colnames(CCLE_gex_test) %in% colnames(PDX_gex_test)]
PDX_gex_test <- PDX_gex_test[,colnames(PDX_gex_test) %in% colnames(CCLE_gex_test)]

CCLE_gex_test$source <- c("CCLE")
PDX_gex_test$source <- c("PDX")

CCLE_PDX_gex <- rbind(CCLE_gex_test, PDX_gex_test)


pr=prcomp(CCLE_PDX_gex[,-c(1,318)])
plot(pr$x[,1],pr$x[,2],col=as.factor(CCLE_PDX_gex$source))
