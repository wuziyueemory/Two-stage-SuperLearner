############################################################################
# BOLD
# 1. Stage-1 numerical results
############################################################################
## load packages
library(SuperLearner)
library(e1071)
library(kernlab)
library(xgboost)
library(ROCR)
library(arm)
library(caret)
library(extraTrees)
library(gam)
library(ipred)
library(KernelKnn)
library(class)
library(MASS)
library(rpart)

## Loading BOLD data
train <- get(load("~/Desktop/First paper/HSORM/BOLD/bold.RData"))
rownames(train) <- seq(1,nrow(train))


# All calculations are based on 10-fold cross-validation

############################################################################
# Spine-related total RVU
# total_spine_RVU12mo
############################################################################
# 10-fold cross validation
set.seed(1)
flds <- createFolds(seq(1,nrow(train)), k = 10, list = TRUE, returnTrain = FALSE)

# create empty object to store results
pred_onestage.total <- rep(NA, nrow(train))
cvRisk.total <- NULL
coef.total <- NULL


for (i in 1:10){
  print(i)
  # create train & test set
  flags <- unlist(flds[i])
  train_temp <- train[-flags, ]
  test_temp <- train[flags,]
  
  ## fit one-stage super learner
  ## Outcome is binary (Y>0)
  binary.fit <- SuperLearner(Y=as.numeric(train_temp$total_spine_RVU12mo>0), 
                             X=train_temp[,-c(52:57)], 
                             newX = test_temp[,-c(52:57)], 
                             family=binomial(), 
                             SL.library=c("SL.glm", "SL.glmnet", "SL.ranger", "SL.ksvm", 
                                          "SL.ipredbagg", "SL.xgboost"),
                             verbose=FALSE,
                             method=method.CC_nloglik,
                             obsWeights = NULL,
                             control=list(saveCVFitLibrary=F),
                             cvControl=list(V = 10))
  
  cvRisk.temp <- binary.fit$cvRisk
  cvRisk.total <- cbind(cvRisk.total, cvRisk.temp)
  
  coef.temp <- binary.fit$coef
  coef.total <- cbind(coef.total, coef.temp)
  
  pred_onestage.total[flags] <- binary.fit$SL.predict
  
}

cvRisk.total
apply(cvRisk.total, 1, mean)

coef.total
apply(coef.total, 1, mean)

## Evaluation metrics
# Review AUC - Area Under Curve
y_test.total <- as.numeric(train$total_spine_RVU12mo>0)
pred_rocr.total = ROCR::prediction(pred_onestage.total, y_test.total)
auc.total = ROCR::performance(pred_rocr.total, measure = "auc", x.measure = "cutoff")@y.values[[1]]
auc.total

# Accuracy
acc.total <- performance(pred_rocr.total, "acc")
# max(acc@y.values[[1]])
acc.total@y.values[[1]][max(which(acc.total@x.values[[1]] >= 0.5))]

# Phi/mat
phi.total <- performance(pred_rocr.total, "phi")
# phi.total@y.values[[1]][min(which(phi.total@x.values[[1]] >= 0.5))]
max(phi.total@y.values[[1]], na.rm=T)

# Negative Logloss
nlogloss.total <- -mean(ifelse(y_test.total, log(pred_onestage.total), log(1 - pred_onestage.total)))
nlogloss.total













############################################################################
# Spine-related imaging RVU
# total_image_spine_RVU12mo
############################################################################

# create empty object to store results
pred_onestage.image <- rep(NA, nrow(train))
cvRisk.image <- NULL
coef.image <- NULL


for (i in 1:10){
  print(i)
  # create train & test set
  flags <- unlist(flds[i])
  train_temp <- train[-flags, ]
  test_temp <- train[flags,]
  
  ## fit one-stage super learner
  ## Outcome is binary (Y>0)
  binary.fit <- SuperLearner(Y=as.numeric(train_temp$total_image_spine_RVU12mo>0), 
                             X=train_temp[,-c(52:57)], 
                             newX = test_temp[,-c(52:57)], 
                             family=binomial(), 
                             SL.library=c("SL.glm", "SL.glmnet", "SL.ranger", "SL.ksvm", 
                                          "SL.ipredbagg", "SL.xgboost"),
                             verbose=FALSE,
                             method=method.CC_nloglik,
                             obsWeights = NULL,
                             control=list(saveCVFitLibrary=F),
                             cvControl=list(V = 10))
  
  cvRisk.temp <- binary.fit$cvRisk
  cvRisk.image <- cbind(cvRisk.image, cvRisk.temp)
  
  coef.temp <- binary.fit$coef
  coef.image <- cbind(coef.image, coef.temp)
  
  pred_onestage.image[flags] <- binary.fit$SL.predict
  
}

cvRisk.image
apply(cvRisk.image, 1, mean)
mean(cvRisk.image[4,], na.rm = T)

coef.image
apply(coef.image, 1, mean)

## Evaluation metrics
# Review AUC - Area Under Curve
y_test.image <- as.numeric(train$total_image_spine_RVU12mo>0)
pred_rocr.image = ROCR::prediction(pred_onestage.image, y_test.image)
auc.image = ROCR::performance(pred_rocr.image, measure = "auc", x.measure = "cutoff")@y.values[[1]]
auc.image

# Accuracy
acc.image <- performance(pred_rocr.image, "acc")
# max(acc@y.values[[1]])
acc.image@y.values[[1]][max(which(acc.image@x.values[[1]] >= 0.5))]

# Phi/mat
phi.image <- performance(pred_rocr.image, "phi")
phi.image@y.values[[1]][max(which(phi.image@x.values[[1]] >= 0.5))]

# Negative Logloss
nlogloss.image <- -mean(ifelse(y_test.image, log(pred_onestage.image), log(1 - pred_onestage.image)))
nlogloss.image
















############################################################################
# Spine-related physical therapy RVU
# total_manual_RVU12mo
############################################################################

# create empty object to store results
pred_onestage.manual <- rep(NA, nrow(train))
cvRisk.manual <- NULL
coef.manual <- NULL


for (i in 1:10){
  print(i)
  # create train & test set
  flags <- unlist(flds[i])
  train_temp <- train[-flags, ]
  test_temp <- train[flags,]
  
  ## fit one-stage super learner
  ## Outcome is binary (Y>0)
  binary.fit <- SuperLearner(Y=as.numeric(train_temp$total_manual_RVU12mo>0), 
                             X=train_temp[,-c(52:57)], 
                             newX = test_temp[,-c(52:57)], 
                             family=binomial(), 
                             SL.library=c("SL.glm", "SL.glmnet", "SL.ranger", "SL.ksvm", 
                                          "SL.ipredbagg", "SL.xgboost"),
                             verbose=FALSE,
                             method=method.CC_nloglik,
                             obsWeights = NULL,
                             control=list(saveCVFitLibrary=F),
                             cvControl=list(V = 10))
  
  cvRisk.temp <- binary.fit$cvRisk
  cvRisk.manual <- cbind(cvRisk.manual, cvRisk.temp)
  
  coef.temp <- binary.fit$coef
  coef.manual <- cbind(coef.manual, coef.temp)
  
  pred_onestage.manual[flags] <- binary.fit$SL.predict
  
}

cvRisk.manual
apply(cvRisk.manual, 1, mean)

coef.manual
apply(coef.manual, 1, mean)

## Evaluation metrics
# Review AUC - Area Under Curve
y_test.manual <- as.numeric(train$total_manual_RVU12mo>0)
pred_rocr.manual = ROCR::prediction(pred_onestage.manual, y_test.manual)
auc.manual = ROCR::performance(pred_rocr.manual, measure = "auc", x.measure = "cutoff")@y.values[[1]]
auc.manual

# Accuracy
acc.manual <- performance(pred_rocr.manual, "acc")
# max(acc@y.values[[1]])
acc.manual@y.values[[1]][max(which(acc.manual@x.values[[1]] >= 0.5))]

# Phi/mat
phi.manual <- performance(pred_rocr.manual, "phi")
phi.manual@y.values[[1]][max(which(phi.manual@x.values[[1]] >= 0.5))]

# Negative Logloss
nlogloss.manual <- -mean(ifelse(y_test.manual, log(pred_onestage.manual), log(1 - pred_onestage.manual)))
nlogloss.manual
















############################################################################
# Spine-related injection RVU
# total_percut_RVU12mo
############################################################################

# create empty object to store results
pred_onestage.percut <- rep(NA, nrow(train))
cvRisk.percut <- NULL
coef.percut <- NULL


for (i in 1:10){
  print(i)
  # create train & test set
  flags <- unlist(flds[i])
  train_temp <- train[-flags, ]
  test_temp <- train[flags,]
  
  ## fit one-stage super learner
  ## Outcome is binary (Y>0)
  binary.fit <- SuperLearner(Y=as.numeric(train_temp$total_percut_RVU12mo>0), 
                             X=train_temp[,-c(52:57)], 
                             newX = test_temp[,-c(52:57)], 
                             family=binomial(), 
                             SL.library=c("SL.glm", "SL.glmnet", "SL.ranger", "SL.ksvm", 
                                          "SL.ipredbagg", "SL.xgboost"),
                             verbose=FALSE,
                             method=method.CC_nloglik,
                             obsWeights = NULL,
                             control=list(saveCVFitLibrary=F),
                             cvControl=list(V = 10))
  
  cvRisk.temp <- binary.fit$cvRisk
  cvRisk.percut <- cbind(cvRisk.percut, cvRisk.temp)
  
  coef.temp <- binary.fit$coef
  coef.percut <- cbind(coef.percut, coef.temp)
  
  pred_onestage.percut[flags] <- binary.fit$SL.predict
  
}

cvRisk.percut
apply(cvRisk.percut, 1, mean)

coef.percut
apply(coef.percut, 1, mean)

## Evaluation metrics
# Review AUC - Area Under Curve
y_test.percut <- as.numeric(train$total_percut_RVU12mo>0)
pred_rocr.percut = ROCR::prediction(pred_onestage.percut, y_test.percut)
auc.percut = ROCR::performance(pred_rocr.percut, measure = "auc", x.measure = "cutoff")@y.values[[1]]
auc.percut

# Accuracy
acc.percut <- performance(pred_rocr.percut, "acc")
# max(acc@y.values[[1]])
acc.percut@y.values[[1]][max(which(acc.percut@x.values[[1]] >= 0.5))]

# Phi/mat
phi.percut <- performance(pred_rocr.percut, "phi")
phi.percut@y.values[[1]][max(which(phi.percut@x.values[[1]] >= 0.5))]

# Negative Logloss
nlogloss.percut <- -mean(ifelse(y_test.percut, log(pred_onestage.percut), log(1 - pred_onestage.percut)))
nlogloss.percut
