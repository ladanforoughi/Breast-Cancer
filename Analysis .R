options(digits = 3)

if(!require(pacman))install.packages("pacman")
pacman::p_load(
   matrixStats,
   tidyverse,
   caret,
   dslabs,
   kableExtra,
   GMCM,
   GGally,
   gam
)

data(brca)
str(brca)

# Number of sample in dataset
Number_of_samples <- dim(brca$x)[1]
# Number of predictor in dataset
Number_of_predictors <- dim(brca$x)[2]
# portion of Malignant
portion_of_melignant <- mean(brca$y == "M")
sum(brca$y == "M")
# portion of Benign
Portion_of_benign <- mean(brca$y == "B")
sum(brca$y == "B")


Info_dataset <- cbind(Number_of_samples,
                      Number_of_predictors,
                      Portion_of_benign,
                      portion_of_melignant)

kable(t(Info_dataset),
      "pandoc",
      caption = "Information of Breast cancer dataset",
      align = "c")

# which column of matrix has highest and lowest mean and standard deviation  
which.max(colMeans(brca$x))
which.min(colMeans(brca$x))
which.max(colSds(brca$x))
which.min(colSds(brca$x))

# scaling the matrix
x_centered <- sweep(brca$x,2,colMeans(brca$x))
x_scaled <- sweep(x_centered,2,colSds(brca$x),FUN = '/')

# verified the standardized matrix by column 1 and 30
sd_col1<- sd(x_scaled[,1])
avg_col1 <- mean(x_scaled[,1])
sd_col30<- sd(x_scaled[,30])
avg_col30 <- mean(x_scaled[,30])
verified <- t(cbind(sd_col1,avg_col1,sd_col30,avg_col30))
kable(verified,"pandoc",
      caption = "verifification of Standardized column 1 and 30", 
      align = "c")

# Calculate the distance between all samples in scaled matrix and average of them
d_samples <- dist(x_scaled)
mean(d_samples)

# the average distance between the first sample, 
#which is benign, and other benign samples
dist_BtoB <- as.matrix(d_samples)[1, brca$y == "B"]
mean(dist_BtoB[2:length(dist_BtoB)])

# the average distance between the first sample, 
#which is malignant, and other malignant samples
dist_MtoM <- as.matrix(d_samples)[1, brca$y == "M"]
mean(dist_MtoM[2:length(dist_MtoM)])

# Heatmap of features
d_features <- dist(t(x_scaled))
heatmap(as.matrix(d_features),labRow = NA, labCol = NA)

# Hierarchical clustering
h <- hclust(d_features)
plot(h, cex = 0.75)
groups <- cutree(h, k = 5)
split(names(groups), groups)

# Principal Component Analysis (PCA)
pca <- prcomp(x_scaled)
summary(pca)
# Proportion of variance of the first principal component
# see PC1 Cumulative Proportion
# How many principal components are required to explain at least 90% of the variance?
# # first value of Cumulative Proportion that exceeds 0.9: PC7

#Plot the first two principal components with color representing tumor type (benign/malignant).
data.frame(pca$x[,1:2],type = brca$y) %>% 
   ggplot(aes(PC1,PC2, col = type)) +
   geom_point()
#Make a boxplot of the first 10 PCs grouped by tumor type.
data.frame(type = brca$y, pca$x[,1:10]) %>% 
   gather(key = "PC" , value = "value", -type) %>%
   ggplot(aes(PC,value, fill = type)) +
   geom_boxplot()
# Which PCs are significantly different enough by tumor 
#type that there is no overlap in the interquartile ranges 
#(IQRs) for benign and malignant samples?
#PC1

#Set the seed to 1, then create a data partition splitting 
#brca$y and the scaled version of the brca$x matrix into a 
#20% test set and 80% train using the following code:
set.seed(1, sample.kind = "Rounding")
test_index <- createDataPartition(brca$y, p = 0.2, times = 1, list = FALSE)
train_set_x <- x_scaled[-test_index,]
train_set_y <- brca$y[-test_index]
test_set_x <- x_scaled[test_index,]
test_set_y <- brca$y[test_index]

#What proportion of the training set is benign?
mean(train_set_y == "B")
# What proportion of the test set is benign?
mean(test_set_y == "B")

#K-means Clustering
predict_kmeans <- function(x, k) {
   centers <- k$centers    # extract cluster centers
   # calculate distance to cluster centers
   distances <- sapply(1:nrow(x), function(i){
      apply(centers, 1, function(y) dist(rbind(x[i,], y)))
   })
   max.col(-t(distances))  # select cluster with min distance to center
}

#Set the seed to 3. Perform k-means clustering on the training
# set with 2 centers and assign the output to k. Then use the 
#predict_kmeans() function to make predictions on the test set.
# What is the overall accuracy?
set.seed(3, sample.kind = "Rounding")
k <- kmeans(train_set_x, centers = 2)
k_mean_predic <- ifelse(predict_kmeans(test_set_x,k)==1,"B","M")
mean(k_mean_predic == test_set_y)

#What proportion of benign tumors are correctly identified?
sensitivity(factor(k_mean_predic), test_set_y, positive = "B")
#What proportion of malignant tumors are correctly identified?
sensitivity(factor(k_mean_predic), test_set_y, positive = "M")

# Logistic regression model
#What is the accuracy of the logistic regression model on the test set?
train_glm <- train(train_set_x,train_set_y,method = "glm")
glm_predic <- predict(train_glm, test_set_x)
mean(glm_predic == test_set_y)

# LDA and QDA models
train_lda <- train(train_set_x, train_set_y, method = "lda")
lda_predic <- predict(train_lda,test_set_x)
mean(lda_predic == test_set_y)

train_qda <- train(train_set_x, train_set_y, method = "qda")
qda_predic <- predict(train_qda,test_set_x)
mean(qda_predic == test_set_y)

#Loess model
set.seed(5, sample.kind = "Rounding")
train_loess <- train(train_set_x,train_set_y, method = "gamLoess")
loess_predic <- predict(train_loess, test_set_x)
mean(loess_predic == test_set_y)

# K-nearest neighbors model
set.seed(7, sample.kind = "Rounding")
train_knn <- train(train_set_x, train_set_y, 
                   method = "knn",
                   tuneGrid = data.frame(k= seq(3,21,2)))
train_knn$bestTune
plot(train_knn)

knn_predic <- predict(train_knn, test_set_x)
mean(knn_predic == test_set_y)
varImp(train_knn)
# Random forest model
set.seed(9, sample.kind = "Rounding")
train_rf <- train(train_set_x, train_set_y, 
                  method = "rf",
                  tuneGrid = data.frame(mtry= c(3,5,7,9)))
train_rf$bestTune
plot(train_rf)
rf_predic <- predict(train_rf,test_set_x)
mean(rf_predic == test_set_y)
# The most important variable in the random forest model
varImp(train_rf)
plot(varImp(train_rf))
# Ensemble models 
ensemble_model <- cbind(rf = rf_predic == "B",
                        knn = knn_predic == "B",
                        lda = lda_predic == "B",
                        qda = qda_predic == "B",
                        loess = loess_predic == "B",
                        glm = glm_predic == "B",
                        kmean = k_mean_predic == "B")
ensemble_predic <- ifelse(rowMeans(ensemble_model) > 0.5 , "B", "M")
mean(ensemble_predic == test_set_y)

#Make a table of the accuracies of the 7 models and the accuracy of the ensemble model.
models <- c("glm","lda","qda","k_means","knn","rf","ensemble")
accuracy <- c(mean(glm_predic == test_set_y),
              mean(lda_predic == test_set_y),
              mean(qda_predic == test_set_y),
              mean(k_mean_predic == test_set_y),
              mean(knn_predic == test_set_y),
              mean(rf_predic == test_set_y),
              mean(ensemble_predic == test_set_y))
data.frame(Model = models, Accuracy = accuracy)


