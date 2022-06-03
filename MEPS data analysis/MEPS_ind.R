# get environment variables
# if you want to see commands in output file
options(echo=TRUE) 
# import args
args=(commandArgs(TRUE))
iter <- eval(parse(text = args[1]))
print(iter)


parameter_grid <- expand.grid(
seed = 1:1
)


############################################################################
# MEPS
# 2. Final prediction results
############################################################################

# load libraries
library(mgcv)
library(quadprog)
library(SuperLearner)
library(nloptr)
library(pscl) 
library(VGAM) 
library(cplm)
library(caret) 
library(randomForest)
library(ranger)
library(kernlab)
library(xgboost)
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

## source functions 
source("/projects/dbenkes/ziyue/sim_1/MEPS/1/new_twosl.R")


## Loading MEPS data
train <- get(load(paste0("/projects/dbenkes/ziyue/sim_1/MEPS/1/train.RData")))
test <- get(load(paste0("/projects/dbenkes/ziyue/sim_1/MEPS/1/test.RData")))

rownames(train) <- seq(1,nrow(train))
rownames(test) <- seq(1,nrow(test))



## SL library
SL.lib.stage1 <- c("SL.glm", "SL.glmnet", "SL.ranger", "SL.ksvm", "SL.ipredbagg", "SL.xgboost")

SL.lib.stage2 <- c("SL.logOLS.smear", "SL.gammaLogGLM", "SL.gammaIdentityGLM", "SL.manningGLM",
                   "SL.gengamma", "SL.coxph", "SL.wangZhou","SL.gilleskie", "SL.rf.caret1",
                   "SL.glmnet", "SL.ksvm", "SL.ipredbagg", "SL.xgboost")

SL.lib.stage.single <- c("SL.mean", "SL.lm", "SL.zip", "SL.zinb", "SL.tobit", "SL.tweedie", "SL.rf.caret1",
                         "SL.glmnet", "SL.ksvm", "SL.ipredbagg", "SL.xgboost")


###################################################################################
# Fit two-stage Super Learner with indicators
###################################################################################

twostage.fit <- twostageSL.indicator(Y = train$TOTEXP, 
                           X = train[,-c(1,21)], 
                           newX = test[,-c(1,21)],
                           library.2stage = list(stage1=SL.lib.stage1,
                                                 stage2=SL.lib.stage2),
                           library.1stage = SL.lib.stage.single,
                           twostage = TRUE,
                           family.1 = binomial,
                           family.2 = gaussian,
                           family.single = gaussian,
                           cvControl = list(V = 10),
                           # threshold value
                           c = c(0.25, 0.5))

# Save coefficients
coef.twostage.indicator <- data.frame("Algorithm" = twostage.fit$libraryNames,
                                      "CV.Risk" = twostage.fit$cvRisk,
                                      "Coef" = twostage.fit$coef)


###################################################################################
# Fit original two-stage Super Learner (no indicator)
###################################################################################
# extract twostage matrix z.twostage
z.twostage <- twostage.fit$Z[,1:89]
twostage.name <- colnames(twostage.fit$library.predict[,1:89])

# get optimum weights for each algorithm in two-stage
getCoef.twostage <- method.CC_LS.scale()$computeCoef(Z = z.twostage, 
                                                     Y = train$TOTEXP,
                                                     libraryNames = twostage.name,
                                                     verbose = FALSE)

# Save coefficients
coef.twostage <- data.frame("Algorithm" = twostage.name,
                            "CV.Risk" = getCoef.twostage$cvRisk,
                            "Coef" = getCoef.twostage$coef)

# Prediction for each algorithm in two-stage superlearner
predY.twostage = twostage.fit$library.predict[,1:89]

# Compute twostage superlearner predictions on newX
twostage.pred <- twostage.fit$method$computePred(predY = predY.twostage, coef = coef.twostage$Coef, 
                                                 control = twostage.fit$control)


###################################################################################
# Fit one-stage Super Learner
###################################################################################
# extract onestage matrix z1
z.onestage <- twostage.fit$Z[,79:89]
onestage.name <- colnames(twostage.fit$library.predict[,79:89])

# get optimum weights for each algorithm in one-stage
getCoef.onestage <- method.CC_LS.scale()$computeCoef(Z = z.onestage,
                                                     Y = train$TOTEXP,
                                                     libraryNames = onestage.name,
                                                     verbose=FALSE)
# Save coefficients
coef.onestage <- data.frame("Algorithm" = onestage.name,
                            "CV.Risk" = getCoef.onestage$cvRisk,
                            "Coef" = getCoef.onestage$coef)

# Prediction for each algorithm in one-stage superlearner
predY.onestage = twostage.fit$library.predict[,79:89]

# Compute onestage superlearner predictions on newX.
onestage.pred <- twostage.fit$method$computePred(predY = predY.onestage, coef = coef.onestage$Coef, 
                                                 control = twostage.fit$control)


####################################################################################
# get discrete two-stage superlearner with indicators
####################################################################################
discrete.pred.twostage.indicator <- twostage.fit$library.predict[,which.min(twostage.fit$cvRisk)]


###################################################################################
# get discrete two-stage superlearner without indicators
###################################################################################
discrete.pred.twostage <- twostage.fit$library.predict[,1:89][,which.min(twostage.fit$cvRisk[1:89])]



###################################################################################
# Store results
###################################################################################
## get prediction performance
# MSE
mse.value <- c(# Individual algorithm
               apply(twostage.fit$library.predict, 2, function(x) mean((test$TOTEXP-x)^2)),
               # One-stage SL
               mean((test$TOTEXP - onestage.pred)^2),
               # Two-stage SL
               mean((test$TOTEXP - twostage.pred)^2),
               # Two-stage SL with indicator
               mean((test$TOTEXP - twostage.fit$SL.predict)^2),
               # Discrete two-stage SL
               mean((test$TOTEXP - discrete.pred.twostage)^2),
               # Discfege two-stage SL with indicators
               mean((test$TOTEXP - discrete.pred.twostage.indicator)^2))

mse <- data.frame("Algo" = c(twostage.fit$libraryNames, "One-stage SL", "Two-stage SL", 
                             "Two-stage SL + Ind", "Discrete two-stage", "Discrete two-stage + Ind"),
                  "MSE" = mse.value)


# MAE
mae.value <- c(# Individual algorithm
               apply(twostage.fit$library.predict, 2, function(x) mean(abs(test$TOTEXP-x))),
               # One-stage SL
               mean(abs(test$TOTEXP - onestage.pred)),
               # Two-stage SL
               mean(abs(test$TOTEXP - twostage.pred)),
               # Two-stage SL with indicator
               mean(abs(test$TOTEXP - twostage.fit$SL.predict)),
               # Discrete two-stage SL
               mean(abs(test$TOTEXP - discrete.pred.twostage)),
               # Discfege two-stage SL with indicators
               mean(abs(test$TOTEXP - discrete.pred.twostage.indicator)))

mae <- data.frame("Algo" = c(twostage.fit$libraryNames, "One-stage SL", "Two-stage SL", 
                             "Two-stage SL + Ind", "Discrete two-stage", "Discrete two-stage + Ind"),
                  "MAE" = mae.value)

# R sqaure
Rsq.value <-  1 - mse.value/var(test$TOTEXP)

Rsq <- data.frame("Algo" = c(twostage.fit$libraryNames, "One-stage SL", "Two-stage SL", 
                             "Two-stage SL + Ind", "Discrete two-stage", "Discrete two-stage + Ind"),
                  "Rsq" = Rsq.value)



###################################################################################
# Save results
###################################################################################
# save output

save(coef.twostage.indicator, file=paste0("/projects/dbenkes/ziyue/sim_1/MEPS/1/results/coef_two_stage_ind_seed=", parameter_grid$seed[iter],".RData"))

save(coef.twostage, file=paste0("/projects/dbenkes/ziyue/sim_1/MEPS/1/results/coef_two_stage_seed=", parameter_grid$seed[iter],".RData"))

save(coef.onestage, file=paste0("/projects/dbenkes/ziyue/sim_1/MEPS/1/results/coef_one_stage_seed=", parameter_grid$seed[iter],".RData"))

save(mse, file=paste0("/projects/dbenkes/ziyue/sim_1/MEPS/1/results/MSE_seed=", parameter_grid$seed[iter],".RData"))

save(mae, file=paste0("/projects/dbenkes/ziyue/sim_1/MEPS/1/results/MAE_seed=", parameter_grid$seed[iter],".RData"))

save(Rsq, file=paste0("/projects/dbenkes/ziyue/sim_1/MEPS/1/results/Rsq_seed=", parameter_grid$seed[iter],".RData"))
