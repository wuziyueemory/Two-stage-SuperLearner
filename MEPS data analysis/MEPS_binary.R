############################################################################
# MEPS
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

## Loading MEPS data
train <- get(load(paste0("~/Desktop/First paper/HSORM/MEPS/train.RData")))
test <- get(load(paste0("~/Desktop/First paper/HSORM/MEPS/test.RData")))


## fit one-stage super learner
## Outcome is binary (Y>0)

binary.fit <- SuperLearner(Y=as.numeric(train$TOTEXP>0), 
                           X=train[,-c(1,21)], 
                           newX = test[,-c(1,21)], 
                           family=binomial(), 
                           SL.library=c("SL.glm", "SL.glmnet", "SL.ranger", "SL.ksvm", 
                                        "SL.ipredbagg", "SL.xgboost"),
                           verbose=FALSE,
                           method=method.CC_nloglik,
                           obsWeights = NULL,
                           control=list(saveCVFitLibrary=T),
                           cvControl=list(V = 10))


## Weights
binary.fit


# Predict back on the test dataset.
# onlySL is set to TRUE so we don't fit algorithms that had weight = 0, saving computation.
binary.pred = predict(binary.fit, test[,-c(1,21)], onlySL = TRUE)


# Histogram of our predicted values.
library(ggplot2)
qplot(binary.pred$pred[, 1]) + theme_minimal() + xlab("Predicted probability") + ylab("Frequency") + 
  ggtitle("Distribution of predicted probability of Y>0")


# Scatterplot of original values (0, 1) and predicted values.
# Ideally we would use jitter or slight transparency to deal with overlap.
qplot(as.numeric(test$TOTEXP>0), binary.pred$pred[, 1]) + theme_minimal() + xlab("Truth") + ylab("Frequency") + 
  ggtitle("Distribution of true Y>0")




## Evaluation metrics
# Review AUC - Area Under Curve
y_test <- as.numeric(test$TOTEXP>0)
pred_rocr = ROCR::prediction(binary.pred$pred, y_test)
auc = ROCR::performance(pred_rocr, measure = "auc", x.measure = "cutoff")@y.values[[1]]
auc

# Accuracy
acc <- performance(pred_rocr, "acc")
# max(acc@y.values[[1]])
acc@y.values[[1]][max(which(acc@x.values[[1]] >= 0.5))]

# Phi/mat
phi <- performance(pred_rocr, "phi")
phi@y.values[[1]][max(which(phi@x.values[[1]] >= 0.5))]

# Negative Logloss
nlogloss <- -mean(ifelse(y_test, log(binary.pred$pred), log(1-binary.pred$pred)))
nlogloss
