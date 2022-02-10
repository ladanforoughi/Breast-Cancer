if(!require(pacman))install.packages("pacman")
pacman::p_load(
  ggplot2,
  dplyr,
  ggthemes,
  corrplot,
  corrgram,
  caret,
  dslabs,
  matrixStats, # colSdv defined 
  tidyr, # for gather function
  factoextra, # for clusstering plot
  princom, # for PCA
  gam, # for loess model
)
options(digits = 3)

data("brca")
breast.cancer <- data.frame(brca)

################### Dimention and Properties of Data set ############
# Remove "x." from name of columns
colnames(breast.cancer) <- gsub("x.","",colnames(breast.cancer))

# The First Six rows of data set
head(breast.cancer)

# The structure of Data set 
str(breast.cancer)

# Missing Value 
any(is.na(breast.cancer))

# Summary Data set 
summary(breast.cancer)

# The definition of them are below:
# radius. Nucleus radius (mean of distances from center to points on perimeter).
# texture. Nucleus texture (standard deviation of grayscale values).
# perimeter. Nucleus perimeter.
# area. Nucleus area.
# smoothness. Nucleus smoothness (local variation in radius lengths).
# compactness. Nucleus compactness (perimeter^2/area - 1).
# concavity, Nucleus concavity (severity of concave portions of the contour).
# concave_pts. Number of concave portions of the nucleus contour.
# symmetry. Nucleus symmetry.
# fractal_dim. Nucleus fractal dimension ("coastline approximation" -1).

table(breast.cancer$y)

############################### Scaling  ####################
scale_x <- scale(breast.cancer[,1:30])

################################ calculate the Distination ###################
dis.x <- dist(scale_x)

#################### Principle component Analysis (PCA) #####################
pc.breastcancer <- princomp(scale_x, cor = TRUE)

# Information of PC
names(pc.breastcancer)

# Quick summary
summary(pc.breastcancer)

#Eigenvalues and Eigenvectors
Eigenvectors <- pc.breastcancer$loadings
Eigenvalues <- pc.breastcancer$sdev * pc.breastcancer$sdev

round(cor(breast.cancer[1:30],pc.breastcancer$scores),3)
plot(round(cor(breast.cancer[1:30],pc.breastcancer$scores),3))

# scree plot to find out importance of comp
screeplot(pc.breastcancer,type = "l", main = "Screeplot for breast cancer data")
abline(3 , 0, col = "red", lty = 4)

plot(pc.breastcancer$scores[,1:2], xlab = "PC1", ylab ="PC2", type = "n")
points(pc.breastcancer$scores[,1:2], cex = 0.5)
text(pc.breastcancer$scores[,1:2], label = breast.cancer$y, cex = 0.5)


# The Second method
pca <- prcomp(scale_x)
summary(pca)
  # The PC1 and PC2 has the 64% of the variability 

data.frame(pca$x[,1:2], type = breast.cancer$y) %>% 
  ggplot(aes(PC1,PC2,color = type))+
  geom_point(size = 2) + 
  ggtitle("Principal Component Analysis")
  
  # The malignant tends to have a larger value of PC1 compared to benign

data.frame(pca$x[,1:10], type = breast.cancer$y) %>%
  gather(key = "PC", value = "value", -type) %>%
  ggplot(aes(PC, value, fill = type)) +
  geom_boxplot()
  
  # PC1 are significantly different enough by tumor type that there
  # is no overlap in the interquartile ranges (IQRs) for benign and 
  # malignant samples 

###################### Splitting data to train and test ##############
set.seed(123,sample.kind = "Rounding")
test_index <- createDataPartition(breast.cancer$y, times = 1, p = 0.2, list = FALSE)

train_x <- scale_x[-test_index,]
test_x <- scale_x[test_index,]
train_y <- breast.cancer$y[-test_index]
test_y <- breast.cancer$y[test_index]

###################### Logistic regression model #########
train_glm <- train(train_x,train_y, method = "glm")
glm_predict <- predict(train_glm, test_x)
# mean(glm_predict == test_y)
Accuracy_glm <- confusionMatrix(factor(glm_predict), factor(test_y))$overall["Accuracy"]
Sensitivity_glm <- sensitivity(factor(glm_predict), factor(test_y))
Specificity_glm <- specificity(factor(glm_predict), factor(test_y))


Results <- tibble(method = "Logistic Regression Model",
                  Accuracy = Accuracy_glm,
                  Sensitivity = Sensitivity_glm,
                  Specificity = Specificity_glm)

print(Results)
################## Linear Discriminant Analysis (LDA) #######################
train_lda <- train(train_x,train_y, method = "lda")
lda_predict <- predict(train_lda, test_x)
# mean(lda_predict == test_y)
Accuracy_lda <- confusionMatrix(factor(lda_predict), factor(test_y))$overall["Accuracy"]
Sensitivity_lda <- sensitivity(factor(lda_predict), factor(test_y))
Specificity_lda <- specificity(factor(lda_predict), factor(test_y))


Results <- add_row(Results , method = "Linear Discriminant Analysis",
                  Accuracy = Accuracy_lda,
                  Sensitivity = Sensitivity_lda,
                  Specificity = Specificity_lda)

print(Results)

################## Quadratic Discriminant Analysis (QDA) #######################
train_qda <- train(train_x,train_y, method = "qda")
qda_predict <- predict(train_qda, test_x)
# mean(qda_predict == test_y)
Accuracy_qda <- confusionMatrix(factor(qda_predict), factor(test_y))$overall["Accuracy"]
Sensitivity_qda <- sensitivity(factor(qda_predict), factor(test_y))
Specificity_qda <- specificity(factor(qda_predict), factor(test_y))


Results <- add_row(Results , method = "Quadratic Discriminant Analysis",
                   Accuracy = Accuracy_qda,
                   Sensitivity = Sensitivity_qda,
                   Specificity = Specificity_qda)

print(Results)

################################## Loess model #####################
set.seed(5, sample.kind = "Rounding")
train_loess <- train(train_x,train_y, method = "gamLoess")
loess_predict <- predict(train_loess, test_x)
# mean(loess_predict == test_y)
Accuracy_loess <- confusionMatrix(factor(loess_predict), factor(test_y))$overall["Accuracy"]
Sensitivity_loess <- sensitivity(factor(loess_predict), factor(test_y))
Specificity_loess <- specificity(factor(loess_predict), factor(test_y))


Results <- add_row(Results , method = "Loess Model",
                   Accuracy = Accuracy_loess,
                   Sensitivity = Sensitivity_loess,
                   Specificity = Specificity_loess)

print(Results)

############################# K-Nearest Neighbors Model ############
set.seed(7, sample.kind = "Rounding")
# find the best k
train_knn <- train(train_x, train_y ,
                   method = "knn", 
                   tuneGrid = data.frame(k = seq(1,21,2)))

Best.k <- train_knn$bestTune
Best.k

knn_predict <- predict(train_knn, test_x)
# mean(knn_predict == test_y)
Accuracy_knn <- confusionMatrix(factor(knn_predict), factor(test_y))$overall["Accuracy"]
Sensitivity_knn <- sensitivity(factor(knn_predict), factor(test_y))
Specificity_knn <- specificity(factor(knn_predict), factor(test_y))


Results <- add_row(Results , method = "KNN Model",
                   Accuracy = Accuracy_knn,
                   Sensitivity = Sensitivity_knn,
                   Specificity = Specificity_knn)

print(Results)

########################### Random Forest Model ####################
set.seed(9, sample.kind = "Rounding")
train_rf <- train(train_x,train_y, method = "rf",
                  tuneGrid = data.frame(mtry = seq(1,21,2)),
                  importance = TRUE)

train_rf$bestTune

rf_predict <- predict(train_rf, test_x)
# mean(knnrf_predict == test_y)
Accuracy_rf <- confusionMatrix(factor(rf_predict), factor(test_y))$overall["Accuracy"]
Sensitivity_rf <- sensitivity(factor(rf_predict), factor(test_y))
Specificity_rf <- specificity(factor(rf_predict), factor(test_y))


Results <- add_row(Results , method = "Random Forest Model",
                   Accuracy = Accuracy_rf,
                   Sensitivity = Sensitivity_rf,
                   Specificity = Specificity_rf)

print(Results)

#The most important variable in the random forest model
varImp(train_rf)
ggplot(varImp(train_rf))

######################### Ensemble ##############################3
ensemble <- cbind(glm = glm_predict == "B", 
                  lda = lda_predict == "B", 
                  qda = qda_predict == "B", 
                  loess = loess_predict == "B", 
                  rf = rf_predict == "B", 
                  knn = knn_predict == "B")
ensemble_predict <- ifelse(rowMeans(ensemble) > 0.5, "B", "M")

Accuracy_ensemble <- confusionMatrix(factor(ensemble_predict), factor(test_y))$overall["Accuracy"]
Sensitivity_ensemble <- sensitivity(factor(ensemble_predict), factor(test_y))
Specificity_ensemble <- specificity(factor(ensemble_predict), factor(test_y))


Results <- add_row(Results , method = "Ensemble Model",
                   Accuracy = Accuracy_ensemble,
                   Sensitivity = Sensitivity_ensemble,
                   Specificity = Specificity_ensemble)

print(Results)

##################### K-means Clustering ###########################
brca.label <- brca$y
table(brca.label)
brca_new <- brca$x
brca_data <- data.frame(brca_new)

#x_scaled.new <- scale(bcra_new)
brca_newnew.scale <- scale(brca_newnew)

#data.new <- dist(x_scaled.new)
data.newnew <- dist(brca_newnew.scale)

fviz_nbclust(brca_newnew.scale,kmeans ,method = "wss") + labs(subtitle = "Elbow Method")
# fviz_nbclust(x_scaled.new,kmeans ,method = "wss") + labs(subtitle = "Elbow Method")

km_out <- kmeans(brca_newnew.scale, centers = 2,nstart = 100)
print(km_out)

km.cluster <- km_out$cluster
rownames(brca_newnew.scale) <- paste(brca$y,1:dim(brca$x)[1], sep ="_")
fviz_cluster(list(data = brca_newnew.scale,cluster = km.cluster))

table(km.cluster,brca$y)

hc.out <- hclust(data.newnew, method = "complete")
plot(hc.out)

rect.hclust(hc.out,k = 2 , border = 2:3)

brca.cluster <-cutree(hc.out,k=2)
split(names(brca.cluster), brca.cluster)

fviz_cluster(list(data = x_scaled.new,cluster = brca.cluster))
table(brca.cluster,brca$y)


corr_mat <- cor(x_scaled.new[,1:ncol(brca$x)])
corrplot(corr_mat, order = "hclust", tl.cex = 0.5, addrect = 5)

