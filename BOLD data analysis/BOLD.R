################################################################################
# BOLD data analysis -- 4 outcomes
################################################################################

# load libraries
library(readxl)
library(mgcv)
library(quadprog)
library(SuperLearner)
library(earth)
library(pscl) 
library(VGAM) 
library(cplm)
library(caret) 
library(randomForest) 
library(e1071) 
library(moments) 
library(flexsurv) 
library(survival) 
library(quantreg) 
library(sandwich)
library(glmnet)
library(foreign)
library(survey)
library(dplyr)
library(nnet)
library(gbm)

# source functions 
source("/home/zwu56/BOLD/twostageSL.R")
source("/home/zwu56/BOLD/Estimators.R")

# load dataset
train <- read.csv("/home/zwu56/BOLD/bold data.csv")


################################################################################################
# Fit two stage Super Learner
################################################################################################


#######################################################################
# Spine-related total RVUs
#######################################################################
# 10-fold cross validation
flds <- createFolds(seq(1,nrow(train)), k = 10, list = TRUE, returnTrain = FALSE)

pred_twostage <- rep(NA,nrow(train))
pred_onestage <- rep(NA,nrow(train))
pred_discrete <- rep(NA,nrow(train))
pred_single <- matrix(NA,nrow=nrow(train),ncol=38)

for (i in 1:10){
  # create train & test set
  flags <- unlist(flds[i])
  train_temp <- train[-flags, ]
  test_temp <- train[flags,]
  # fit two-stage super learner
  twostage.fit <- twostageSL(Y = train_temp$total_spine_RVU12mo, X = train_temp[,c(1:24)], newX = test_temp[,c(1:24)],
                             library.2stage = list(stage1=c("SL.glm","SL.rf.caret1","SL.glmnet"),
                                                   stage2=c("SL.logOLS.smear","SL.manningGLM",
                                                          "SL.gammaIdentityGLM",
                                                          "SL.nnet",
                                                          "SL.gbm.caret1","SL.earth",
                                                          "SL.rpart","SL.ipredbagg",
                                                          "SL.rf.caret1","SL.glmnet")),
                             library.1stage = c("SL.gbm.caret1","SL.lm","SL.rpart","SL.ipredbagg","SL.earth",
                                              "SL.rf.caret1","SL.glmnet","SL.nnet"),
                             twostage = TRUE,
                             family.1 = binomial,
                             family.2 = gaussian,
                             family.single = gaussian,
                             cvControl = list(V = 10))
  
  
  # construct one-stage superlearner
  # extract onestage matrix z1
  z1 <- twostage.fit$Z[,31:38]
  onestagename <- colnames(twostage.fit$library.predict[,31:38])
  # get optimum weights for each algorithm in one-stage
  getCoef <- method.CC_LS.scale()$computeCoef(Z=z1,Y=train_temp$total_spine_RVU12mo,libraryNames=onestagename,
                                              verbose=FALSE)
  coef.onestage <- getCoef$coef
  # Prediction for each algorithm in one-stage superlearner
  predY.onestage = twostage.fit$library.predict[,31:38]
  # Compute onestage superlearner predictions on newX.
  onestage.pred <- twostage.fit$method$computePred(predY = predY.onestage, coef = coef.onestage, 
                                                   control = twostage.fit$control)
  # get discrete two-stage superlearner
  discrete.pred <- twostage.fit$library.predict[,which.min(twostage.fit$cvRisk)]
  
  pred_twostage[flags] <- twostage.fit$SL.predict
  pred_onestage[flags] <- onestage.pred
  pred_discrete[flags] <- discrete.pred
  pred_single[flags,] <- twostage.fit$library.predict
}


## get prediction performance
# MSE
mse_1 <- c(apply(pred_single, 2, function(x) mean((train$total_spine_RVU12mo-x)^2)),
         mean((train$total_spine_RVU12mo - pred_onestage)^2),
         mean((train$total_spine_RVU12mo - pred_twostage)^2),
         mean((train$total_spine_RVU12mo - pred_discrete)^2))

# MAE
mae_1 <- c(apply(pred_single, 2, function(x) mean(abs(train$total_spine_RVU12mo-x))),
         mean(abs(train$total_spine_RVU12mo - pred_onestage)),
         mean(abs(train$total_spine_RVU12mo - pred_twostage)),
         mean(abs(train$total_spine_RVU12mo - pred_discrete)))

# R sqaure
Rsq_1 <-  1 - mse/var(train$total_spine_RVU12mo)




#######################################################################
# Spine-related injection RVUs
#######################################################################

# 10-fold cross validation
flds <- createFolds(seq(1,nrow(train)), k = 10, list = TRUE, returnTrain = FALSE)

pred_twostage <- rep(NA,nrow(train))
pred_onestage <- rep(NA,nrow(train))
pred_discrete <- rep(NA,nrow(train))
pred_single <- matrix(NA,nrow=nrow(train),ncol=38)

for (i in 1:10){
  # create train & test set
  flags <- unlist(flds[i])
  train_temp <- train[-flags, ]
  test_temp <- train[flags,]
  # fit two-stage super learner
  twostage.fit <- twostageSL(Y = train_temp$total_percut_RVU12mo, X = train_temp[,c(1:24)], newX = test_temp[,c(1:24)],
                             library.2stage = list(stage1=c("SL.glm","SL.rf.caret1","SL.glmnet"),
                                                   stage2=c("SL.logOLS.smear","SL.gammaLogGLM",
                                                          "SL.gammaIdentityGLM",
                                                          "SL.nnet",
                                                          "SL.gbm.caret1","SL.earth",
                                                          "SL.rpart","SL.ipredbagg",
                                                          "SL.rf.caret1","SL.glmnet")),
                             library.1stage = c("SL.gbm.caret1","SL.lm","SL.rpart","SL.ipredbagg","SL.earth",
                                              "SL.rf.caret1","SL.glmnet","SL.nnet"),
                             twostage = TRUE,
                             family.1 = binomial,
                             family.2 = gaussian,
                             family.single = gaussian,
                             cvControl = list(V = 10))
  
  
  # construct one-stage superlearner
  # extract onestage matrix z1
  z1 <- twostage.fit$Z[,31:38]
  onestagename <- colnames(twostage.fit$library.predict[,31:38])
  # get optimum weights for each algorithm in one-stage
  getCoef <- method.CC_LS.scale()$computeCoef(Z=z1,Y=train_temp$total_percut_RVU12mo,libraryNames=onestagename,
                                              verbose=FALSE)
  coef.onestage <- getCoef$coef
  # Prediction for each algorithm in one-stage superlearner
  predY.onestage = twostage.fit$library.predict[,31:38]
  # Compute onestage superlearner predictions on newX.
  onestage.pred <- twostage.fit$method$computePred(predY = predY.onestage, coef = coef.onestage, 
                                                   control = twostage.fit$control)
  # get discrete two-stage superlearner
  discrete.pred <- twostage.fit$library.predict[,which.min(twostage.fit$cvRisk)]
  
  pred_twostage[flags] <- twostage.fit$SL.predict
  pred_onestage[flags] <- onestage.pred
  pred_discrete[flags] <- discrete.pred
  pred_single[flags,] <- twostage.fit$library.predict
}


## get prediction performance
# MSE
mse_2 <- c(apply(pred_single, 2, function(x) mean((train$total_percut_RVU12mo-x)^2)),
         mean((train$total_percut_RVU12mo - pred_onestage)^2),
         mean((train$total_percut_RVU12mo - pred_twostage)^2),
         mean((train$total_percut_RVU12mo - pred_discrete)^2))

# MAE
mae_2 <- c(apply(pred_single, 2, function(x) mean(abs(train$total_percut_RVU12mo-x))),
         mean(abs(train$total_percut_RVU12mo - pred_onestage)),
         mean(abs(train$total_percut_RVU12mo - pred_twostage)),
         mean(abs(train$total_percut_RVU12mo - pred_discrete)))

# R sqaure
Rsq_2 <-  1 - mse/var(train$total_percut_RVU12mo)





#######################################################################
# Spine-related imaging RVUs
#######################################################################

# 10-fold cross validation
flds <- createFolds(seq(1,nrow(train)), k = 10, list = TRUE, returnTrain = FALSE)

pred_twostage <- rep(NA,nrow(train))
pred_onestage <- rep(NA,nrow(train))
pred_discrete <- rep(NA,nrow(train))
pred_single <- matrix(NA,nrow=nrow(train),ncol=38)

for (i in 1:10){
  # create train & test set
  flags <- unlist(flds[i])
  train_temp <- train[-flags, ]
  test_temp <- train[flags,]
  # fit two-stage super learner
  twostage.fit <- twostageSL(Y = train_temp$total_image_spine_RVU12mo, X = train_temp[,c(1:24)], newX = test_temp[,c(1:24)],
                             library.2stage = list(stage1=c("SL.glm","SL.rf.caret1","SL.glmnet"),
                                                   stage2=c("SL.logOLS.smear","SL.gammaLogGLM",
                                                          "SL.gammaIdentityGLM",
                                                          "SL.nnet",
                                                          "SL.gbm.caret1","SL.earth",
                                                          "SL.rpart","SL.ipredbagg",
                                                          "SL.rf.caret1","SL.glmnet")),
                             library.1stage = c("SL.gbm.caret1","SL.lm","SL.rpart","SL.ipredbagg","SL.earth",
                                              "SL.rf.caret1","SL.glmnet","SL.nnet"),
                             twostage = TRUE,
                             family.1 = binomial,
                             family.2 = gaussian,
                             family.single = gaussian,
                             cvControl = list(V = 10))
  
  
  # construct one-stage superlearner
  # extract onestage matrix z1
  z1 <- twostage.fit$Z[,31:38]
  onestagename <- colnames(twostage.fit$library.predict[,31:38])
  # get optimum weights for each algorithm in one-stage
  getCoef <- method.CC_LS.scale()$computeCoef(Z=z1,Y=train_temp$total_image_spine_RVU12mo,libraryNames=onestagename,
                                              verbose=FALSE)
  coef.onestage <- getCoef$coef
  # Prediction for each algorithm in one-stage superlearner
  predY.onestage = twostage.fit$library.predict[,31:38]
  # Compute onestage superlearner predictions on newX.
  onestage.pred <- twostage.fit$method$computePred(predY = predY.onestage, coef = coef.onestage, 
                                                   control = twostage.fit$control)
  # get discrete two-stage superlearner
  discrete.pred <- twostage.fit$library.predict[,which.min(twostage.fit$cvRisk)]
  
  pred_twostage[flags] <- twostage.fit$SL.predict
  pred_onestage[flags] <- onestage.pred
  pred_discrete[flags] <- discrete.pred
  pred_single[flags,] <- twostage.fit$library.predict
}


## get prediction performance
# MSE
mse_3 <- c(apply(pred_single, 2, function(x) mean((train$total_image_spine_RVU12mo-x)^2)),
         mean((train$total_image_spine_RVU12mo - pred_onestage)^2),
         mean((train$total_image_spine_RVU12mo - pred_twostage)^2),
         mean((train$total_image_spine_RVU12mo - pred_discrete)^2))

# MAE
mae_3 <- c(apply(pred_single, 2, function(x) mean(abs(train$total_image_spine_RVU12mo-x))),
         mean(abs(train$total_image_spine_RVU12mo - pred_onestage)),
         mean(abs(train$total_image_spine_RVU12mo - pred_twostage)),
         mean(abs(train$total_image_spine_RVU12mo - pred_discrete)))

# R sqaure
Rsq_3 <-  1 - mse/var(train$total_image_spine_RVU12mo)







#######################################################################
# Spine-related physical therapy RVUs
#######################################################################

# 10-fold cross validation
flds <- createFolds(seq(1,nrow(train)), k = 10, list = TRUE, returnTrain = FALSE)

pred_twostage <- rep(NA,nrow(train))
pred_onestage <- rep(NA,nrow(train))
pred_discrete <- rep(NA,nrow(train))
pred_single <- matrix(NA,nrow=nrow(train),ncol=38)

for (i in 1:10){
  # create train & test set
  flags <- unlist(flds[i])
  train_temp <- train[-flags, ]
  test_temp <- train[flags,]
  # fit two-stage super learner
  twostage.fit <- twostageSL(Y = train_temp$total_manual_RVU12mo, X = train_temp[,c(1:24)], newX = test_temp[,c(1:24)],
                             library.2stage = list(stage1=c("SL.glm","SL.rf.caret1","SL.glmnet"),
                                                   stage2=c("SL.logOLS.smear","SL.gammaLogGLM",
                                                          "SL.gammaIdentityGLM",
                                                          "SL.nnet",
                                                          "SL.gbm.caret1","SL.earth",
                                                          "SL.rpart","SL.ipredbagg",
                                                          "SL.rf.caret1","SL.glmnet")),
                             library.1stage = c("SL.gbm.caret1","SL.lm","SL.rpart","SL.ipredbagg","SL.earth",
                                              "SL.rf.caret1","SL.glmnet","SL.nnet"),
                             twostage = TRUE,
                             family.1 = binomial,
                             family.2 = gaussian,
                             family.single = gaussian,
                             cvControl = list(V = 10))
  
  
  # construct one-stage superlearner
  # extract onestage matrix z1
  z1 <- twostage.fit$Z[,31:38]
  onestagename <- colnames(twostage.fit$library.predict[,31:38])
  # get optimum weights for each algorithm in one-stage
  getCoef <- method.CC_LS.scale()$computeCoef(Z=z1,Y=train_temp$total_manual_RVU12mo,libraryNames=onestagename,
                                              verbose=FALSE)
  coef.onestage <- getCoef$coef
  # Prediction for each algorithm in one-stage superlearner
  predY.onestage = twostage.fit$library.predict[,31:38]
  # Compute onestage superlearner predictions on newX.
  onestage.pred <- twostage.fit$method$computePred(predY = predY.onestage, coef = coef.onestage, 
                                                   control = twostage.fit$control)
  # get discrete two-stage superlearner
  discrete.pred <- twostage.fit$library.predict[,which.min(twostage.fit$cvRisk)]
  
  pred_twostage[flags] <- twostage.fit$SL.predict
  pred_onestage[flags] <- onestage.pred
  pred_discrete[flags] <- discrete.pred
  pred_single[flags,] <- twostage.fit$library.predict
}


## get prediction performance
# MSE
mse_4 <- c(apply(pred_single, 2, function(x) mean((train$total_manual_RVU12mo-x)^2)),
         mean((train$total_manual_RVU12mo - pred_onestage)^2),
         mean((train$total_manual_RVU12mo - pred_twostage)^2),
         mean((train$total_manual_RVU12mo - pred_discrete)^2))

# MAE
mae_4 <- c(apply(pred_single, 2, function(x) mean(abs(train$total_manual_RVU12mo-x))),
         mean(abs(train$total_manual_RVU12mo - pred_onestage)),
         mean(abs(train$total_manual_RVU12mo - pred_twostage)),
         mean(abs(train$total_manual_RVU12mo - pred_discrete)))

# R sqaure
Rsq_4 <-  1 - mse/var(train$total_manual_RVU12mo)





########################### Combine algorithm names with evaluation metrics ###########################

# algorithm name
algo.name <- c( "S1: SL.glm + S2: SL.logOLS.smear",           
                "S1: SL.glm + S2: SL.gammaLogGLM",            
                "S1: SL.glm + S2: SL.gammaIdentityGLM",       
                "S1: SL.glm + S2: SL.nnet",             
                "S1: SL.glm + S2: SL.gbm",               
                "S1: SL.glm + S2: SL.earth",                  
                "S1: SL.glm + S2: SL.rpart",               
                "S1: SL.glm + S2: SL.ipredbagg",              
                "S1: SL.glm + S2: SL.rf.caret1",             
                "S1: SL.glm + S2: SL.glmnet",             
                "S1: SL.rf.caret1 + S2: SL.logOLS.smear",     
                "S1: SL.rf.caret1 + S2: SL.gammaLogGLM",     
                "S1: SL.rf.caret1 + S2: SL.gammaIdentityGLM", 
                "S1: SL.rf.caret1 + S2: SL.nnet",       
                "S1: SL.rf.caret1 + S2: SL.gbm",         
                "S1: SL.rf.caret1 + S2: SL.earth",            
                "S1: SL.rf.caret1 + S2: SL.rpart",         
                "S1: SL.rf.caret1 + S2: SL.ipredbagg",        
                "S1: SL.rf.caret1 + S2: SL.rf.caret1",        
                "S1: SL.rf.caret1 + S2: SL.glmnet",       
                "S1: SL.glmnet + S2: SL.logOLS.smear",    
                "S1: SL.glmnet + S2: SL.gammaLogGLM",     
                "S1: SL.glmnet + S2: SL.gammaIdentityGLM",
                "S1: SL.glmnet + S2: SL.nnet",      
                "S1: SL.glmnet + S2: SL.gbm",        
                "S1: SL.glmnet + S2: SL.earth",           
                "S1: SL.glmnet + S2: SL.rpart",        
                "S1: SL.glmnet + S2: SL.ipredbagg",       
                "S1: SL.glmnet + S2: SL.rf.caret1",       
                "S1: SL.glmnet + S2: SL.glmnet",      
                "Single: SL.gbm",                                
                "Single: SL.lm",                                  
                "Single: SL.rpart",                                 
                "Single: SL.ipredbagg", 
                "Single: SL.earth",                
                "Single: SL.rf.caret1",                           
                "Single: SL.glmnet",
                "Single: SL.nnet",
                "One-stage SuperLearner",
                "Two-stage SuperLearner",
                "Discrete SuperLearner")
algo.num = length(algo.name)

# MSE
mse_spine_total <- data.frame("Algorithm"=algo.name,"MSE"=mse_1)
mse_spine_injection <- data.frame("Algorithm"=algo.name,"MSE"=mse_2)
mse_spine_imaging <- data.frame("Algorithm"=algo.name,"MSE"=mse_3)
mse_spine_manual <- data.frame("Algorithm"=algo.name,"MSE"=mse_4)


# MAE
mae_spine_total <- data.frame("Algorithm"=algo.name,"MAE"=mae_1)
mae_spine_injection <- data.frame("Algorithm"=algo.name,"MAE"=mae_2)
mae_spine_imaging <- data.frame("Algorithm"=algo.name,"MAE"=mae_3)
mae_spine_manual <- data.frame("Algorithm"=algo.name,"MAE"=mae_4)

# R square
Rsq_spine_total <- data.frame("Algorithm"=algo.name,"R_square"=Rsq_1)
Rsq_spine_injection <- data.frame("Algorithm"=algo.name,"R_square"=Rsq_2)
Rsq_spine_imaging <- data.frame("Algorithm"=algo.name,"R_square"=Rsq_3)
Rsq_spine_manual <- data.frame("Algorithm"=algo.name,"R_square"=Rsq_4)



# save output
save(mse_spine_total, file=paste0("/home/zwu56/BOLD/MSE_spine_total",".RData"))
save(mse_spine_injection, file=paste0("/home/zwu56/BOLD/MSE_spine_injection",".RData"))
save(mse_spine_imaging, file=paste0("/home/zwu56/BOLD/MSE_spine_imaging",".RData"))
save(mse_spine_manual, file=paste0("/home/zwu56/BOLD/MSE_spine_manual",".RData"))

save(mae_spine_total, file=paste0("/home/zwu56/BOLD/MAE_spine_total",".RData"))
save(mae_spine_injection, file=paste0("/home/zwu56/BOLD/MAE_spine_injection",".RData"))
save(mae_spine_imaging, file=paste0("/home/zwu56/BOLD/MAE_spine_imaging",".RData"))
save(mae_spine_manual, file=paste0("/home/zwu56/BOLD/MAE_spine_manual",".RData"))

save(Rsq_spine_total, file=paste0("/home/zwu56/BOLD/Rsq_spine_total",".RData"))
save(Rsq_spine_injection, file=paste0("/home/zwu56/BOLD/Rsq_spine_injection",".RData"))
save(Rsq_spine_imaging, file=paste0("/home/zwu56/BOLD/Rsq_spine_imaging",".RData"))
save(Rsq_spine_manual, file=paste0("/home/zwu56/BOLD/Rsq_spine_manual",".RData"))
