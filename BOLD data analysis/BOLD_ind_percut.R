# get environment variables
# if you want to see commands in output file
options(echo=TRUE) 
# import args
args=(commandArgs(TRUE))
iter <- eval(parse(text = args[1]))
print(iter)



parameter_grid <- expand.grid(
  seed = 1:10
)


############################################################################
# BOLD
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
library(matrixcalc)
library(Matrix)
library(arm)

## source functions 
source("/projects/dbenkes/ziyue/sim_1/BOLD/new_twosl.R")


## Loading MEPS data
# train <- get(load("/projects/dbenkes/ziyue/sim_1/BOLD/bold.RData"))
train <- get(load("~/Desktop/First paper/HSORM/BOLD/bold.RData"))
rownames(train) <- seq(1,nrow(train))


## SL library
SL.lib.stage1 <- c("SL.glm", "SL.glmnet", "SL.ranger", "SL.ksvm", "SL.ipredbagg", "SL.xgboost")

SL.lib.stage2 <- c("SL.logOLS.smear", "SL.gammaLogGLM", "SL.gammaIdentityGLM", "SL.bayesglm",
                         "SL.rpart", "SL.coxph", "SL.glmnet", "SL.rf.caret1",
                         "SL.ksvm", "SL.ipredbagg", "SL.xgboost")

SL.lib.stage.single.percut <- c("SL.mean", "SL.lm", "SL.glmnet", "SL.zip", "SL.zinb", "SL.tweedie", 
                         "SL.ksvm","SL.rf.caret1", "SL.ipredbagg", "SL.xgboost", "SL.nnet")


###############################################################################################################
# Fit two-stage Super Learner with indicators
###############################################################################################################
# 10-fold cross validation

# flds <- readRDS("/projects/dbenkes/ziyue/sim_1/BOLD/folds.RData")

flds <- readRDS("~/Desktop/First paper/HSORM/BOLD/folds.RData")

############################################################################
# Spine-related injection RVU
# total_percut_RVU12mo
############################################################################

# create empty object to store results
pred_twostage_ind <- rep(NA, nrow(train))
pred_twostage <- rep(NA, nrow(train))
pred_onestage <- rep(NA, nrow(train))
pred_discrete_twostage <- rep(NA, nrow(train))
pred_discrete_twostage_ind <- rep(NA, nrow(train))
pred_single <- matrix(NA,nrow=nrow(train),ncol = 209)


  # Starting fitting
  i <- parameter_grid$seed[iter]
  print(i)
  # create train & test set
  flags <- unlist(flds[i])
  length(flags)
  train_temp <- train[-flags, ]
  test_temp <- train[flags,]
  dim(test_temp)

  rownames(train_temp) <- seq(1,nrow(train_temp))
  rownames(test_temp) <- seq(1,nrow(test_temp))
  
  # Fit two-stage SL with indicator
  twostage.fit <- twostageSL.indicator(Y = train_temp$total_percut_RVU12mo, 
                                       X = train_temp[,-c(52:57)], 
                                       newX = test_temp[,-c(52:57)],
                                       library.2stage = list(stage1=SL.lib.stage1,
                                                             stage2=SL.lib.stage2),
                                       library.1stage = SL.lib.stage.single.percut,
                                       twostage = TRUE,
                                       family.1 = binomial,
                                       family.2 = gaussian,
                                       family.single = gaussian,
                                       cvControl = list(V = 5),
                                       # threshold value
                                       c = c(0.01, 0.05))
  
  # Save cvRisk
  cvRisk_twostage_ind <- twostage.fit$cvRisk
  
  # Save coefficient
  coef_twostage_ind <- twostage.fit$coef
  
  # Compute twostage superlearner predictions on newX
  pred_twostage_ind[flags] <- twostage.fit$SL.predict

  # Compute individual algorithm predictions on newX
  colnames(pred_single) <- twostage.fit$libraryNames
  pred_single[flags,] <- twostage.fit$library.predict
  
  
##################### Fit original two-stage Super Learner (no indicator) #######################

  # extract twostage matrix z.twostage
  z.twostage <- twostage.fit$Z[,1:77]
  twostage.name <- colnames(twostage.fit$library.predict[,1:77])

  # get optimum weights for each algorithm in two-stage
  getCoef.twostage <- method.CC_LS.scale()$computeCoef(Z = z.twostage, 
                                                       Y = train_temp$total_percut_RVU12mo,
                                                       libraryNames = twostage.name,
                                                       verbose = FALSE)

  # Save cvRisk
  cvRisk_twostage <- getCoef.twostage$cvRisk

  # Save coefficients
  coef_twostage <- getCoef.twostage$coef

  # Prediction for each algorithm in two-stage superlearner
  predY.twostage = twostage.fit$library.predict[,1:77]

  # Compute twostage superlearner predictions on newX
  pred_twostage[flags] <- twostage.fit$method$computePred(predY = predY.twostage, coef = coef_twostage, 
                                                   control = twostage.fit$control)


################################ Fit one-stage Super Learner #################################

  # extract onestage matrix z1
  z.onestage <- twostage.fit$Z[,67:77]
  onestage.name <- colnames(twostage.fit$library.predict[,67:77])

  # get optimum weights for each algorithm in one-stage
  getCoef.onestage <- method.CC_LS.scale()$computeCoef(Z = z.onestage,
                                                       Y = train_temp$total_percut_RVU12mo,
                                                       libraryNames = onestage.name,
                                                       verbose=FALSE)

  # Save cvRisk
  cvRisk_onestage <- getCoef.onestage$cvRisk
  
  # Save coefficients
  coef_onestage <- getCoef.onestage$coef

  # Prediction for each algorithm in one-stage superlearner
  predY.onestage = twostage.fit$library.predict[,67:77]

  # Compute onestage superlearner predictions on newX.
  pred_onestage[flags] <- twostage.fit$method$computePred(predY = predY.onestage, coef = coef_onestage, 
                                                                control = twostage.fit$control)



#################### get discrete two-stage superlearner with indicators ######################
  pred_discrete_twostage_ind[flags] <- twostage.fit$library.predict[,which.min(twostage.fit$cvRisk)]


################## get discrete two-stage superlearner without indicators #####################
  pred_discrete_twostage[flags] <- twostage.fit$library.predict[,1:77][,which.min(twostage.fit$cvRisk[1:77])]





###################################################################################
# Save results
###################################################################################
# save output
# cvRisk
save(cvRisk_twostage_ind, file=paste0("/projects/dbenkes/ziyue/sim_1/BOLD/percut/results/percut_cvRisk_two_stage_ind_seed=", parameter_grid$seed[iter],".RData"))

save(cvRisk_twostage, file=paste0("/projects/dbenkes/ziyue/sim_1/BOLD/percut/results/percut_cvRisk_two_stage_seed=", parameter_grid$seed[iter],".RData"))

save(cvRisk_onestage, file=paste0("/projects/dbenkes/ziyue/sim_1/BOLD/percut/results/percut_cvRisk_one_stage_seed=", parameter_grid$seed[iter],".RData"))


# Coef
save(coef_twostage_ind, file=paste0("/projects/dbenkes/ziyue/sim_1/BOLD/percut/results/percut_coef_two_stage_ind_seed=", parameter_grid$seed[iter],".RData"))

save(coef_twostage, file=paste0("/projects/dbenkes/ziyue/sim_1/BOLD/percut/results/percut_coef_two_stage_seed=", parameter_grid$seed[iter],".RData"))

save(coef_onestage, file=paste0("/projects/dbenkes/ziyue/sim_1/BOLD/percut/results/percut_coef_one_stage_seed=", parameter_grid$seed[iter],".RData"))


# Data
save(pred_twostage_ind, file=paste0("/projects/dbenkes/ziyue/sim_1/BOLD/percut/results/percut_pred_twostage_ind_seed=", parameter_grid$seed[iter],".RData"))

save(pred_twostage, file=paste0("/projects/dbenkes/ziyue/sim_1/BOLD/percut/results/percut_pred_twostage_seed=", parameter_grid$seed[iter],".RData"))

save(pred_onestage, file=paste0("/projects/dbenkes/ziyue/sim_1/BOLD/percut/results/percut_pred_onestage_seed=", parameter_grid$seed[iter],".RData"))

save(pred_discrete_twostage, file=paste0("/projects/dbenkes/ziyue/sim_1/BOLD/percut/results/percut_pred_discrete_twostage_seed=", parameter_grid$seed[iter],".RData"))

save(pred_discrete_twostage_ind, file=paste0("/projects/dbenkes/ziyue/sim_1/BOLD/percut/results/percut_pred_discrete_twostage_ind_seed=", parameter_grid$seed[iter],".RData"))

save(pred_single, file=paste0("/projects/dbenkes/ziyue/sim_1/BOLD/percut/results/percut_pred_single_seed=", parameter_grid$seed[iter],".RData"))
