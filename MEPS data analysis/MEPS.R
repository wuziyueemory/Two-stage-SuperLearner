#==================================================================================#
# MEPS data analysis
#==================================================================================#

# load libraries
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
library(gbm) 
library(moments) 
library(flexsurv) 
library(survival) 
library(quantreg) 
library(sandwich)
library(glmnet)
library(foreign)
library(survey)
library(dplyr)
library(ggplot2)

# source functions 
source("/home/zwu56/MEPS/twostageSL.R")
source("/home/zwu56/MEPS/Estimators.R")

# Loading MEPS train_data
# train
train <- read.csv("/home/zwu56/MEPS/train.csv")
# test
test <- read.csv("/home/zwu56/MEPS/test.csv")

#==================================================================================#
# Fit two stage Super Learner
#==================================================================================#

# fit two-stage superlearner 
twostage.fit <- twostageSL(Y = train$TOTEXP, X = train[,-c(1,21)], newX = test[,-c(1,21)],
                           library.2stage = list(stage1=c("SL.glm","SL.rf.caret1","SL.glmnet"),
                                                 stage2=c("SL.logOLS.smear","SL.gammaLogGLM",
                                                          "SL.gammaIdentityGLM",
                                                          "SL.manningGLM",
                                                          "SL.gengamma","SL.coxph",
                                                          "SL.wangZhou","SL.gilleskie",
                                                          "SL.rf.caret1","SL.glmnet")),
                           library.1stage = c("SL.mean","SL.lm","SL.zip","SL.zinb","SL.tobit",
                                              "SL.tweedie","SL.rf.caret1","SL.glmnet"),
                           twostage = TRUE,
                           family.1 = binomial,
                           family.2 = gaussian,
                           family.single = gaussian,
                           cvControl = list(V = 2))

# construct one-stage superlearner
# extract onestage matrix z1
z1 <- twostage.fit$Z[,31:38]
onestagename <- colnames(twostage.fit$library.predict[,31:38])
# get optimum weights for each algorithm in one-stage
getCoef <- method.CC_LS.scale()$computeCoef(Z=z1,Y=train$TOTEXP,libraryNames=onestagename,
                                            verbose=FALSE)
coef.onestage <- getCoef$coef
# Prediction for each algorithm in one-stage superlearner
predY.onestage = twostage.fit$library.predict[,31:38]
# Compute onestage superlearner predictions on newX.
onestage.pred <- twostage.fit$method$computePred(predY = predY.onestage, coef = coef.onestage, 
                                                 control = twostage.fit$control)
# get discrete two-stage superlearner
discrete.pred <- twostage.fit$library.predict[,which.min(twostage.fit$cvRisk)]

## get prediction performance
# MSE
mse <- c(apply(twostage.fit$library.predict, 2, function(x) mean((test$TOTEXP-x)^2)),
         mean((test$TOTEXP - onestage.pred)^2),
         mean((test$TOTEXP - twostage.fit$SL.predict)^2),
         mean((test$TOTEXP - discrete.pred)^2))
# MAE
mae <- c(apply(twostage.fit$library.predict, 2, function(x) mean(abs(test$TOTEXP-x))),
         mean(abs(test$TOTEXP - onestage.pred)),
         mean(abs(test$TOTEXP - twostage.fit$SL.predict)),
         mean(abs(test$TOTEXP - discrete.pred)))
# R^2
Rsq <-  1 - mse/var(test$TOTEXP)

# algorithm name
algo.name <- c( "S1: SL.glm + S2: SL.logOLS.smear",           
                "S1: SL.glm + S2: SL.gammaLogGLM",            
                "S1: SL.glm + S2: SL.gammaIdentityGLM",       
                "S1: SL.glm + S2: SL.manningGLM",             
                "S1: SL.glm + S2: SL.gengamma",               
                "S1: SL.glm + S2: SL.coxph",                  
                "S1: SL.glm + S2: SL.wangZhou",               
                "S1: SL.glm + S2: SL.gilleskie",              
                "S1: SL.glm + S2: SL.rf.caret1",             
                "S1: SL.glm + S2: SL.glmnet",             
                "S1: SL.rf.caret1 + S2: SL.logOLS.smear",     
                "S1: SL.rf.caret1 + S2: SL.gammaLogGLM",     
                "S1: SL.rf.caret1 + S2: SL.gammaIdentityGLM", 
                "S1: SL.rf.caret1 + S2: SL.manningGLM",       
                "S1: SL.rf.caret1 + S2: SL.gengamma",         
                "S1: SL.rf.caret1 + S2: SL.coxph",            
                "S1: SL.rf.caret1 + S2: SL.wangZhou",         
                "S1: SL.rf.caret1 + S2: SL.gilleskie",        
                "S1: SL.rf.caret1 + S2: SL.rf.caret1",        
                "S1: SL.rf.caret1 + S2: SL.glmnet",       
                "S1: SL.glmnet + S2: SL.logOLS.smear",    
                "S1: SL.glmnet + S2: SL.gammaLogGLM",     
                "S1: SL.glmnet + S2: SL.gammaIdentityGLM",
                "S1: SL.glmnet + S2: SL.manningGLM",      
                "S1: SL.glmnet + S2: SL.gengamma",        
                "S1: SL.glmnet + S2: SL.coxph",           
                "S1: SL.glmnet + S2: SL.wangZhou",        
                "S1: SL.glmnet + S2: SL.gilleskie",       
                "S1: SL.glmnet + S2: SL.rf.caret1",       
                "S1: SL.glmnet + S2: SL.glmnet",      
                "Single: SL.mean",                                
                "Single: SL.lm",                                  
                "Single: SL.zip",                                 
                "Single: SL.zinb",                                
                "Single: SL.tobit",                               
                "Single: SL.tweedie",                             
                "Single: SL.rf.caret1",                           
                "Single: SL.glmnet",
                "One-stage SuperLearner",
                "Two-stage SuperLearner",
                "Discrete SuperLearner")
algo.num = length(algo.name)

# Combine algorithm name with evaluation metrics
# MSE
mse_result <- data.frame("Algorithm"=algo.name,"MSE"=mse)
# MAE
mae_result <- data.frame("Algorithm"=algo.name,"MSE"=mae)
# R square
Rsq_result <- data.frame("Algorithm"=algo.name,"MSE"=Rsq)

# Save output
save(mse_result, file = paste0("/home.zwu56/MEPS/MSE.RDATA"))
         
save(mae_result, file = paste0("/home.zwu56/MEPS/MAE.RDATA"))
         
save(Rsq_result, file = paste0("/home.zwu56/MEPS/Rsq.RDATA"))

   
