#                                                                                                                                                   #
#                                                           Imagen: Mental health profiles                                                         #                               
#                                                                                                                                                   #
#####################################################################################################################################################
#written by S.Mareva, Last revised November 2024

# ----------------------------------------------------- Load packages ----------------------------------------------------------------------------####
#manipulation/tables
library(dplyr);library(stargazer); 

#analyses
library(cluster) ; library(psych);  library(MASS); library(smacof);library(uwot); library(lsr); library(caret)
# ----------------------------------------------------- Read data ####
sdqDataT1<-read.csv('imagenSdqT1Age.csv')
sdqDataT2<-read.csv('imagenSdqT2.csv')
# ----------------------------------------------------- Evaluate missingness ####
#check if those missing are different from those that are followed-up
sdqDataT1$followUp<-ifelse(sdqDataT1$PSC2 %in% sdqDataT2$PSC2, 1, 0)

#get a vector with all sdq subscales
sdqScales<-c('semotion','sconduct','shyper','speer','sprosoc')

#compare SDQ
followUpCompTable <-data.frame(); for (i in (c(sdqScales))) {
  followUpCompTable[i,c('N_no_Followup','M_no_Followup','SD_no_Followup')]<-as.data.frame(psych::describe(sdqDataT1[sdqDataT1$followUp == 0,i]))[,c('n','mean','sd')]
  followUpCompTable[i,c('N_followup','M_followup','SD_followup')]<-as.data.frame(psych::describe(sdqDataT1[sdqDataT1$followUp == 1,i]))[,c('n','mean','sd')]
  followUpCompTable[i,c('t','df','p')]<-t.test(sdqDataT1[,i] ~sdqDataT1$followUp, var.equal= FALSE, paired=F)[c('statistic','parameter','p.value')]
  followUpCompTable[i,c('d')]<-cohen.d(sdqDataT1[,i],sdqDataT1$followUp)$cohen.d[2]}
followUpCompTable<-round(followUpCompTable,3) 
# export table
stargazer(followUpCompTable, summary=F, type = 'text', digits =2) 

# ----------------------------------------------------- Participants with two time points ####
#merge
imagenLongi<-merge(sdqDataT1, sdqDataT2, by = 'PSC2',suffixes = c('T1','T2'))
dim(imagenLongi) #1304 final number of participants

#get vectors of SDQ variables
sdqScalesT1<-paste(sdqScales,'T1',sep="")
sdqScalesT2<-paste(sdqScales,'T2',sep="")

#sex
imagenLongiSex<-data.frame(boy= table(imagenLongi$genderT1)[1],  girls= table(imagenLongi$genderT1)[2])
stargazer(imagenLongiSex, summary=F, type = 'text', out='sexBaseline.html', digits =2) 

#age T1
imagenLongiAge<-as.data.frame(describe(as.numeric(imagenLongi$ageFinal)))
stargazer(imagenLongiAge, summary=F, type = 'text', out='ageBasline.html', digits =2) 

#age T2 
imagenLongiAgeT2<-as.data.frame(describe(as.numeric(imagenLongi$ageT2)))
stargazer(imagenLongiAgeT2, summary=F, type = 'text', out='ageFollowup.html', digits =2) 

#sdq descriptives M, SD 
sdqDescrTable<-describe(imagenLongi[,c(sdqScalesT1,sdqScalesT2)])[,c('n','mean','sd')]
stargazer(sdqDescrTable, summary=F, type = 'text', digits =2) 
# ----------------------------------------------------- Subset non-elevated & difficulties groups ####
#screen participants with difficulties
#emotional problems: 0-4 average; >4
#conduct problems: 0-3 average; >3
#hypearcitivity: 0-5 average; >5
#peer: 0-2 average; >2
#prosocial: 7-10 average; <7

#sdq criteria = 1 = selected for further analyes (one of criteria above met)
#sdq criteria = 0 = no difficulites (no criteria above met)
imagenLongi$sdqCriteriaT1<-ifelse((imagenLongi$semotionT1 > 4 | imagenLongi$sconductT1 > 3 | imagenLongi$speerT1 > 2 |
                                    imagenLongi$shyperT1 > 5 | imagenLongi$sprosocT1  < 7),1,0)
imagenLongi$sdqCriteriaT2<-ifelse((imagenLongi$semotionT2 > 4 | imagenLongi$sconductT2 > 3 | imagenLongi$speerT2 > 2 |
                                    imagenLongi$shyperT2 > 5 | imagenLongi$sprosocT2  < 7),1,0)
table(imagenLongi$sdqCriteriaT1) 
table(imagenLongi$sdqCriteriaT2)

#Create a matrix indicating which conditions each participant meets
condition_matrixT1 <- cbind(imagenLongi$semotionT1 > 4,
                            imagenLongi$sconductT1 > 3,
                            imagenLongi$speerT1 > 2,
                            imagenLongi$shyperT1 > 5,
                            imagenLongi$sprosocT1 < 7)

condition_matrixT2 <- cbind(imagenLongi$semotionT2 > 4,
                            imagenLongi$sconductT2 > 3,
                            imagenLongi$speerT2 > 2,
                            imagenLongi$shyperT2 > 5,
                            imagenLongi$sprosocT2 < 7)

# Count the number of conditions each participant meets
imagenLongi$num_conditions_metT1 <- rowSums(condition_matrixT1)
imagenLongi$num_conditions_metT2 <- rowSums(condition_matrixT2)

# Number of participants meeting each number of conditions
table(imagenLongi$num_conditions_metT1)
table(imagenLongi$num_conditions_metT2)

# ----------------------------------------------------- Clustering + cross-validation ####
######################   T1
#T1 umap
set.seed(175); umapSdqT1<-uwot::umap(imagenLongi[imagenLongi$sdqCriteriaT1 ==1 ,sdqScalesT1],metric = 'correlation',
                                         min_dist = 0.001,n_neighbors = 20)
umapSdqT1Df<-data.frame(d1=umapSdqT1[,1], d2 = umapSdqT1[,2])

#T1 kmeans
kmeansT1 <- NbClust::NbClust(umapSdqT1Df, method = "kmeans", index = c("silhouette"))
kmeansT1$All.index
table(kmeansT1$Best.partition)

#T1 cross-validation
set.seed(45); k_values <- 2:15 # Number of clusters to explore (k values from 2 to 15)
n_folds <- 10# Number of folds for cross-validation
t1_clust_data <- umapSdqT1Df[, 1:2] # Data for clustering

# Set up K-Fold Cross-validation using caret package
cv_folds_t1 <- createFolds(1:nrow(t1_clust_data), k = n_folds, list = TRUE, returnTrain = TRUE)
silhouette_metrics_per_k_t1 <- list() # Store clustering metrics for t1

# Cross-validation loop over different values of k for t1
for (k in k_values) {
  fold_metrics_t1 <- data.frame(mean = numeric(n_folds),
                                sd = numeric(n_folds))
  # Data frame to store metrics for each fold
  
  for (i in 1:n_folds) {
    # Split into training and testing sets
    train_index_t1 <- cv_folds_t1[[i]]
    test_index_t1 <- setdiff(1:nrow(t1_clust_data), train_index_t1)
    
    X_train_t1 <- t1_clust_data[train_index_t1, , drop = FALSE]
    X_test_t1 <- t1_clust_data[test_index_t1, , drop = FALSE]
    
    # Perform K-Means clustering using NbClust for the current k
    kmeans_cs_T1 <- NbClust::NbClust(X_train_t1, method = "kmeans", 
                                 index = c("silhouette"), max.nc = k, min.nc = k)
    
    # Extract the "Best Partition" from NbClust results
    best_partitionT1 <-  kmeans_cs_T1$Best.partition
    
    # Evaluate clustering quality on the training set using silhouette score
    dist_train_t1 <- dist(X_train_t1)  # Calculate pairwise distances for silhouette calculation
    silhouette_train_t1 <- silhouette(best_partitionT1, dist_train_t1)  # Compute silhouette score
    
    # Store silhouette scores for the current fold
    fold_metrics_t1$mean[i] <- mean(silhouette_train_t1[, 3])  # Mean silhouette score for the fold
    fold_metrics_t1$sd[i] <- sd(silhouette_train_t1[, 3])      # Standard deviation of silhouette scores for the fold
  }
  
  # Store the metrics for the current k for t1
  silhouette_metrics_per_k_t1[[as.character(k)]] <- fold_metrics_t1}

# Create a summary table of mean and SD of silhouette scores for each k for t1
summary_table_t1 <- data.frame(k = integer(), mean_silhouette = numeric(), sd_silhouette = numeric())

for (k in k_values) {
  # Extract the metrics for the current k for t1
  fold_metrics_t1 <- silhouette_metrics_per_k_t1[[as.character(k)]]
  
  # Calculate mean and SD of silhouette scores across folds for t1
  mean_sil_t1 <- mean(fold_metrics_t1$mean)
  sd_sil_t1 <- mean(fold_metrics_t1$sd)
  
  # Append to the summary table for t1
  summary_table_t1 <- rbind(summary_table_t1, data.frame(k = k, mean_silhouette = mean_sil_t1, sd_silhouette = sd_sil_t1))}

# Print the final summary table for t1
print(summary_table_t1)




######################### T2 
#T2 umap
set.seed(175); umapSdqT2<-uwot::umap(imagenLongi[imagenLongi$sdqCriteriaT2 ==1 ,sdqScalesT2],metric = 'correlation',
                                         min_dist = 0.001,n_neighbors = 20)
umapSdqT2Df<-data.frame(d1=umapSdqT2[,1], d2 = umapSdqT2[,2])


#T2 kmeans 
kmeansT2 <- NbClust::NbClust(umapSdqT2Df, method = "kmeans", index = c("silhouette"))
kmeansT2$All.index
kmeansT2$Best.partition
table(kmeansT2$Best.partition)

#T2 cross validation
set.seed(45)
t2_clust_data <- umapSdqT2Df[, 1:2] # Data for clustering

# Set up K-Fold Cross-validation using caret package
cv_folds_t2 <- createFolds(1:nrow(t2_clust_data), k = n_folds, list = TRUE, returnTrain = TRUE)
silhouette_metrics_per_k_t2 <- list() # Store clustering metrics for t2

# Cross-validation loop over different values of k for t2
for (k in k_values) {
  fold_metrics_t2 <- data.frame(mean = numeric(n_folds),
                                sd = numeric(n_folds))  # Data frame to store metrics for each fold
  
  for (i in 1:n_folds) {
    # Split into training and testing sets
    train_index_t2 <- cv_folds_t2[[i]]
    test_index_t2 <- setdiff(1:nrow(t2_clust_data), train_index_t2)
    
    X_train_t2 <- t2_clust_data[train_index_t2, , drop = FALSE]
    X_test_t2 <- t2_clust_data[test_index_t2, , drop = FALSE]
    
    # Perform K-Means clustering using NbClust for the current k
    kmeans_cs_t2 <- NbClust::NbClust(X_train_t2, method = "kmeans", 
                                 index = c("silhouette"), max.nc = k, min.nc = k)
    
    # Extract the "Best Partition" from NbClust results
    best_partitiont2 <- kmeans_cs_t2$Best.partition
    
    # Evaluate clustering quality on the training set using silhouette score
    dist_train_t2 <- dist(X_train_t2)  # Calculate pairwise distances for silhouette calculation
    silhouette_train_t2 <- silhouette(best_partitiont2, dist_train_t2)  # Compute silhouette score
    
    # Store silhouette scores for the current fold
    fold_metrics_t2$mean[i] <- mean(silhouette_train_t2[, 3])  # Mean silhouette score for the fold
    fold_metrics_t2$sd[i] <- sd(silhouette_train_t2[, 3])      # Standard deviation of silhouette scores for the fold
  }
  
  # Store the metrics for the current k for t2
  silhouette_metrics_per_k_t2[[as.character(k)]] <- fold_metrics_t2
}

# Create a summary table of mean and SD of silhouette scores for each k for t2
summary_table_t2 <- data.frame(k = integer(), mean_silhouette = numeric(), sd_silhouette = numeric())

for (k in k_values) {
  # Extract the metrics for the current k for t2
  fold_metrics_t2 <- silhouette_metrics_per_k_t2[[as.character(k)]]
  
  # Calculate mean and SD of silhouette scores across folds for t2
  mean_sil_t2 <- mean(fold_metrics_t2$mean)
  sd_sil_t2 <- mean(fold_metrics_t2$sd)
  
  # Append to the summary table for t2
  summary_table_t2 <- rbind(summary_table_t2, data.frame(k = k, mean_silhouette = mean_sil_t2, sd_silhouette = sd_sil_t2))
}

# Print the final summary table for t2
print(summary_table_t2)

# ----------------------------------------------------- T1 + T2: clusters description ####
#T1 data with partition
imagenLong$clusterT1<-imagenLong$sdqCriteriaT1
imagenLong[imagenLong$sdqCriteriaT1 ==1 ,'clusterT1']<- kmeansT1$Best.partition
imagenLong$clusterT1Name<-factor(imagenLong$clusterT1,
                                 levels = c(0,1,2,3),
                                 labels = c('ned', 'ext', 'int','social'))


#T2 data with partition
imagenLong$clusterT2<-imagenLong$sdqCriteriaT2
imagenLong[imagenLong$sdqCriteriaT2 ==1 ,'clusterT2']<- kmeansT2$Best.partition
imagenLong$clusterT2Name<-factor(imagenLong$clusterT2,
                                 levels = c(0,1,2,3),
                                 labels = c('ned', 'ext', 'int','social'))


#SDQ-------
#T1 comparison to no difficulties
compareSdqClustersT1<-data.frame()
for (j in 1:3){ for (i in (c(sdqScalesT1))) {
  compareSdqClustersT1[paste(i,j, sep ="_C"),c('N','M','SD')]<-as.data.frame(psych::describe(imagenLong[imagenLong$clusterT1 == j,i]))[,c('n','mean','sd')]
  compareSdqClustersT1[paste(i,j, sep ="_C"),c('N_nd','M_nd','SD_nd')]<-as.data.frame(psych::describe(imagenLong[imagenLong$clusterT1 == 0,i]))[,c('n','mean','sd')]
  compareSdqClustersT1[paste(i,j, sep ="_C"),c('t','df','p')]<-t.test(imagenLong[imagenLong$clusterT1 ==j |imagenLong$clusterT1 == 0,i] ~
                                                                        imagenLong[imagenLong$clusterT1 ==j |imagenLong$clusterT1 == 0,'clusterT1'],
                                                                      data = imagenLong[imagenLong$clusterT1 ==j |imagenLong$clusterT1 == 0,],
                                                                      var.equal= FALSE, paired=F)[c('statistic','parameter','p.value')]
  compareSdqClustersT1[paste(i,j, sep ="_C"),c('low','d','up')]<-cohen.d(imagenLong[imagenLong$clusterT1 ==j |imagenLong$clusterT1 == 0,i],
                                                                         imagenLong[imagenLong$clusterT1 ==j |imagenLong$clusterT1 == 0,'clusterT1'])$cohen.d[,1:3]}}
compareSdqClustersT1<-round(compareSdqClustersT1,3)
compareSdqClustersT1$p.fdr<-p.adjust(compareSdqClustersT1$p, method = 'fdr')

#T2 comparison to no difficulties
compareSdqClustersT2<-data.frame()
for (j in 1:3){ for (i in (c(sdqScalesT2))) {
  compareSdqClustersT2[paste(i,j, sep ="_C"),c('N','M','SD')]<-as.data.frame(psych::describe(imagenLong[imagenLong$clusterT2 == j,i]))[,c('n','mean','sd')]
  compareSdqClustersT2[paste(i,j, sep ="_C"),c('N_nd','M_nd','SD_nd')]<-as.data.frame(psych::describe(imagenLong[imagenLong$clusterT2 == 0,i]))[,c('n','mean','sd')]
  compareSdqClustersT2[paste(i,j, sep ="_C"),c('t','df','p')]<-t.test(imagenLong[imagenLong$clusterT2 ==j |imagenLong$clusterT2 == 0,i] ~
                                                                        imagenLong[imagenLong$clusterT2 ==j |imagenLong$clusterT2 == 0,'clusterT2'],
                                                                      data = imagenLong[imagenLong$clusterT2 ==j |imagenLong$clusterT2 == 0,],
                                                                      var.equal= FALSE, paired=F)[c('statistic','parameter','p.value')]
  compareSdqClustersT2[paste(i,j, sep ="_C"),c('low','d','up')]<-cohen.d(imagenLong[imagenLong$clusterT2 ==j |imagenLong$clusterT2 == 0,i],
                                                                         imagenLong[imagenLong$clusterT2 ==j |imagenLong$clusterT2 == 0,'clusterT2'])$cohen.d[,1:3]}}
compareSdqClustersT2<-round(compareSdqClustersT2,3)
compareSdqClustersT2$p.fdr<-p.adjust(compareSdqClustersT2$p, method = 'fdr')

# ----------------------------------------------------- Age and sex across clusters ####
#T1: sex 
resultsSexT1<-data.frame();  for (i in c('int_t1','ext_t1','social_t1','nd_t1')) {
  resultsSexT1[paste(i),c('boys','girls')]<-chisq.test(x = table(imagenLong$genderT1,imagenLong$clusterT1Factor)[,i], p = prop.table(table(imagenLong$genderT1)))$observed
  resultsSexT1[paste(i),'chi-square']<-chisq.test(x = table(imagenLong$genderT1,imagenLong$clusterT1Factor)[,i], p = prop.table(table(imagenLong$genderT1)))$statistic
  resultsSexT1[paste(i),'p.value']<-chisq.test(x = table(imagenLong$genderT1,imagenLong$clusterT1Factor)[,i], p = prop.table(table(imagenLong$genderT1)))$p.value}
resultsSexT1<-round(resultsSexT1,3)
stargazer(resultsSexT1, summary=F, type = 'text', out='tables/resultsSexT1.html', digits =2) 

clustChi<-chisq.test(x = table(imagenLong$clusterT1Name,imagenLong$clusterT2Name))
clustChi$observed
clustChi$expected

#T1: age
aov(imagenLong$ageFinal ~ imagenLong$clusterT1) %>% summary()

#T2 sex 
resultsSexT2<-data.frame();  for (i in c('int_t2','ext_t2','social_t2','nd_t2')) {
  resultsSexT2[paste(i),c('boys','girls')]<-chisq.test(x = table(imagenLong$genderT2,imagenLong$clusterT2Factor)[,i], p = prop.table(table(imagenLong$genderT2)))$observed
  resultsSexT2[paste(i),'chi-square']<-chisq.test(x = table(imagenLong$genderT2,imagenLong$clusterT2Factor)[,i], p = prop.table(table(imagenLong$genderT2)))$statistic
  resultsSexT2[paste(i),'p.value']<-chisq.test(x = table(imagenLong$genderT2,imagenLong$clusterT2Factor)[,i], p = prop.table(table(imagenLong$genderT2)))$p.value}
round(resultsSexT2,3)
resultsSexT2<-round(resultsSexT2,3)
stargazer(resultsSexT2, summary=F, type = 'text', out='tables/resultsSexT2.html', digits =2) 


#T2 age
aov(imagenLong$ageT2 ~ imagenLong$clusterT2) %>% summary()

# ----------------------------------------------------- Read + merge cognitive data ####
#read data
cogDatT1<-read.csv('cognitive.data.baseline.edit.csv', na.strings = c('NA','N/A','#VALUE!'))
cogDatT2<-read.csv('cognitive.data.T2.csv',na.strings = c('NA','N/A','#VALUE!'))

#vectors with names
cogTasks<-c('SWM_Between_errors','CGT_Delay_aversion',
            'CGT_Risk_adjustment','CGT_Risk_taking')

cogTasksT1<-c('SWM_Between_errorsT1','CGT_Delay_aversionT1',
              'CGT_Risk_adjustmentT1', 'CGT_Risk_takingT1')

cogTasksT2<-c('SWM_Between_errorsT2', 'CGT_Delay_aversionT2', 'CGT_Risk_adjustmentT2','CGT_Risk_takingT2')


#merging cognitive data----
#merge T1 cog
imagenLongCog<-merge(imagenLong, cogDatT1[,-c(2:3)], by = 'PSC2', all.x = T)
imagenLongCogT2<-merge(imagenLong, cogDatT2[,-c(2:3)], by = 'PSC2', all.x = T)
imagenLongCogT2<-merge(imagenLongCog,  imagenLongCogT2[,c('PSC2',cogTasks)],by = 'PSC2',
                       suffixes = c('T1','T2'))
# ----------------------------------------------------- Descriptives Cognitive ####
cogDescrTable<-describe(imagenLongCogT2[,c(cogTasksT1,cogTasksT2)])[,c('n','mean','sd')]
stargazer(cogDescrTable, summary=F, type = 'text', digits =2)
#T1 table comparison to no difficulties----
compareSdqClustersCogT1<-data.frame()
for (j in c('ext','int','social')){ for (i in (c(cogTasksT1))) {
  compareSdqClustersCogT1[paste(i,j, sep =" "),c('N','M','SD')]<-as.data.frame(psych::describe(imagenLongCogT2[imagenLongCogT2$clusterT1Name == j,i]))[,c('n','mean','sd')]
  compareSdqClustersCogT1[paste(i,j, sep =" "),c('N_nd','M_nd','SD_nd')]<-as.data.frame(psych::describe(imagenLongCogT2[imagenLongCogT2$clusterT1Name == 'nd',i]))[,c('n','mean','sd')]
  compareSdqClustersCogT1[paste(i,j, sep =" "),c('t','df','p')]<-t.test(imagenLongCogT2[imagenLongCogT2$clusterT1Name ==j |imagenLongCogT2$clusterT1Name == 'nd',i] ~
                                                                          imagenLongCogT2[imagenLongCogT2$clusterT1Name ==j |imagenLongCogT2$clusterT1Name == 'nd','clusterT1'],
                                                                        data = imagenLongCogT2[imagenLongCogT2$clusterT1Name ==j |imagenLongCogT2$clusterT1Name == 'nd',],
                                                                        var.equal= FALSE, paired=F)[c('statistic','parameter','p.value')]
  compareSdqClustersCogT1[paste(i,j, sep =" "),c('low','d','up')]<-cohen.d(imagenLongCogT2[imagenLongCogT2$clusterT1Name ==j |imagenLongCogT2$clusterT1Name == 'nd',i],
                                                                           as.numeric(as.factor(imagenLongCogT2[imagenLongCogT2$clusterT1Name ==j |
                                                                                                                  imagenLongCogT2$clusterT1Name == 'nd','clusterT1Name'])))$cohen.d[,1:3]}}
compareSdqClustersCogT1<-round(compareSdqClustersCogT1,3)
compareSdqClustersCogT1$p.fdr<-p.adjust(compareSdqClustersCogT1$p, method = 'fdr')
stargazer(compareSdqClustersCogT1, summary=F, type = 'text') 

#T2: table comparison to no difficulties----
compareSdqClustersCogT2<-data.frame()
for (j in c('ext','int','social')){ for (i in (c(cogTasksT2))) {
  compareSdqClustersCogT2[paste(i,j, sep =" "),c('N','M','SD')]<-as.data.frame(psych::describe(imagenLongCogT2[imagenLongCogT2$clusterT2Name == j,i]))[,c('n','mean','sd')]
  compareSdqClustersCogT2[paste(i,j, sep =" "),c('N_nd','M_nd','SD_nd')]<-as.data.frame(psych::describe(imagenLongCogT2[imagenLongCogT2$clusterT2Name == 'nd',i]))[,c('n','mean','sd')]
  compareSdqClustersCogT2[paste(i,j, sep =" "),c('t','df','p')]<-t.test(imagenLongCogT2[imagenLongCogT2$clusterT2Name ==j |imagenLongCogT2$clusterT2Name == 'nd',i] ~
                                                                          imagenLongCogT2[imagenLongCogT2$clusterT2Name ==j |imagenLongCogT2$clusterT2Name == 'nd','clusterT2'],
                                                                        data = imagenLongCogT2[imagenLongCogT2$clusterT2 ==j |imagenLongCogT2$clusterT2Name == 'nd',],
                                                                        var.equal= FALSE, paired=F)[c('statistic','parameter','p.value')]
  compareSdqClustersCogT2[paste(i,j, sep =" "),c('low','d','up')]<-cohen.d(imagenLongCogT2[imagenLongCogT2$clusterT2Name ==j |imagenLongCogT2$clusterT2Name == 'nd',i],
                                                                           imagenLongCogT2[imagenLongCogT2$clusterT2Name ==j |imagenLongCogT2$clusterT2Name == 'nd','clusterT2Name'])$cohen.d[,1:3]}}
compareSdqClustersCogT2<-round(compareSdqClustersCogT2,3)
compareSdqClustersCogT2$p.fdr<-p.adjust(compareSdqClustersCogT2$p, method = 'fdr')
stargazer(compareSdqClustersCogT2, summary=F, type = 'text')

# ----------------------------------------------------- Transitions ####
#get transition variable
imagenLongCogT2$transition<-paste(imagenLongCogT2$clusterT1Name, imagenLongCogT2$clusterT2Name, sep = "->")
transitionNumbers<-as.data.frame(table(imagenLongCogT2$transition))
transitionNumbers$ClusterSize<-c(rep(sum(imagenLong$clusterT1Name == 'ext'),4),
                                 rep(sum(imagenLong$clusterT1Name == 'int'),4),
                                 rep(sum(imagenLong$clusterT1Name == 'ned'),4),
                                 rep(sum(imagenLong$clusterT1Name == 'social'),4))

transitionNumbers$percent<-round((transitionNumbers$Freq*100)/transitionNumbers$ClusterSize,0)

#get p-values
for (i in 1:dim(table4)[1]){ 
  transitionNumbers$p.value[i]<-round(prop.test(x = transitionNumbers$Freq[i],
                                                n = transitionNumbers$ClusterSize[i], 
                                                p = 0.25)$p.value,3)}
transitionNumbers

# ----------------------------------------------------- Transition comparisons ####
meltTransitions<-reshape2::melt(imagenLongCogT2[,c(cogTasksT2,cogTasksT1, 'transition')], id='transition')
meltTransitions$clusterT1Name<-as.factor(meltTransitions$transition)

#all comparisons
transitionsComp<-compare_means(value ~ transition, data = meltTransitions,
                               var.equal = F,group.by = 'variable',
                               method = "t.test") 
