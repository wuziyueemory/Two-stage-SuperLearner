#! /usr/bin/env Rscript


# This file was used to submit the simulation files to a cluster system.
# Using the submit_jobs.sh shell script one can submit each simulation 
# in sequence. First, the created simulated data files are analyzed in 
# the 'run' execution. Then the results are collated in the 'merge' 
# execution

## get environment variables
myscratch <- Sys.getenv('myscratch')
resultdir <- Sys.getenv('resultdir')
job_control <- Sys.getenv("job_control")
iter <- as.numeric(Sys.getenv("iter"))

## simulation parameters
parameter_grid <- expand.grid(
  seed = 1:1000,
  zero_inflate = c(0.05, 0.7),
  dist = c("lognormal", "gamma", "tweedie", "mixture"),
  interaction = c(TRUE, FALSE),
  linear = c(TRUE) ,
  sample_size = c(500, 2000)
)

# source in functions 
source("/home/zwu56/sim/twostageSL.R")
source("/home/zwu56/sim/createData.R")
source("/home/zwu56/sim/Estimators.R")

#################### execute job #######################
if(job_control == "run"){

  # do your simulation for row iter of parameter_grid
  # set seed according to parameter_grid$seeds[iter]
    set.seed(parameter_grid$seed[iter])
    # make training & testing data based on parameter_grid[iter, ]
    train <- createData(n = parameter_grid$sample_size[iter],
                        zero = parameter_grid$zero_inflate[iter],
                        nonzero = parameter_grid$dist[iter],
                        interaction = parameter_grid$interaction[iter],
                        linear = parameter_grid$linear[iter]
                        )
    
    test <- createData(n = parameter_grid$sample_size[iter],
                        zero = parameter_grid$zero_inflate[iter],
                        nonzero = parameter_grid$dist[iter],
                        interaction = parameter_grid$interaction[iter],
                        linear = parameter_grid$linear[iter]
    )

    ## load libraries
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
    
    ## fit two-stage superlearner    
    twostage.fit <- twostageSL(Y = train$y, X = train[,-c(1,12)], newX = test[,-c(1,12)],
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
                               cvControl = list(V = 5))
    
    ###################### construct one-stage superlearner #######################
    # extract onestage matrix z1
    z1 <- twostage.fit$Z[,31:38]
    onestagename <- colnames(twostage.fit$library.predict[,31:38])
    # get optimum weights for each algorithm in one-stage
    getCoef <- method.CC_LS.scale()$computeCoef(Z=z1,Y=train$y,libraryNames=onestagename,
                                                verbose=FALSE)
    coef.onestage <- getCoef$coef
    # Prediction for each algorithm in one-stage superlearner
    predY.onestage = twostage.fit$library.predict[,31:38]
    # Compute onestage superlearner predictions on newX.
    onestage.pred <- twostage.fit$method$computePred(predY = predY.onestage, coef = coef.onestage, 
                                                     control = twostage.fit$control)
    
    # get discrete two-stage superlearner
    discrete.pred <- twostage.fit$library.predict[,which.min(twostage.fit$cvRisk)]
    
    ########################## get prediction metrics ########################
    # MSE
    mse <- c(apply(twostage.fit$library.predict, 2, function(x) mean((test$y-x)^2)),
             mean((test$y - onestage.pred)^2),
             mean((test$y - twostage.fit$SL.predict)^2),
             mean((test$y - discrete.pred)^2))
    
    # MAE
    mae <- c(apply(twostage.fit$library.predict, 2, function(x) mean(abs(test$y-x))),
             mean(abs(test$y - onestage.pred)),
             mean(abs(test$y - twostage.fit$SL.predict)),
             mean(abs(test$y - discrete.pred)))
    
    # R sqaure
    Rsq <-  1 - mse/var(test$y)
    
    ######################### save output ###########################
    save(mse, file=paste0("/home/zwu56/sim/scratch/MSE_n=",parameter_grid$sample_size[iter],
                          "_zero=",parameter_grid$zero_inflate[iter],
                          "_nonzero=",parameter_grid$dist[iter],
                          "_interaction=",parameter_grid$interaction[iter],
                          "_linear=",parameter_grid$linear[iter],
                          "_seed=", parameter_grid$seed[iter],".RData"))
    
    save(mae, file=paste0("/home/zwu56/sim/scratch/MAE_n=",parameter_grid$sample_size[iter],
                          "_zero=",parameter_grid$zero_inflate[iter],
                          "_nonzero=",parameter_grid$dist[iter],
                          "_interaction=",parameter_grid$interaction[iter],
                          "_linear=",parameter_grid$linear[iter],
                          "_seed=", parameter_grid$seed[iter],".RData"))
    
    save(Rsq, file=paste0("/home/zwu56/sim/scratch/Rsq_n=",parameter_grid$sample_size[iter],
                          "_zero=",parameter_grid$zero_inflate[iter],
                          "_nonzero=",parameter_grid$dist[iter],
                          "_interaction=",parameter_grid$interaction[iter],
                          "_linear=",parameter_grid$linear[iter],
                          "_seed=", parameter_grid$seed[iter],".RData"))
    }


############################## merge job ################################### 
if(job_control == "merge"){
  
  # algorithm name
  algo.name <- c( "S1: SL.glm_All + S2: SL.logOLS.smear_All",           
                  "S1: SL.glm_All + S2: SL.gammaLogGLM_All",            
                  "S1: SL.glm_All + S2: SL.gammaIdentityGLM_All",       
                  "S1: SL.glm_All + S2: SL.manningGLM_All",             
                  "S1: SL.glm_All + S2: SL.gengamma_All",               
                  "S1: SL.glm_All + S2: SL.coxph_All",                  
                  "S1: SL.glm_All + S2: SL.wangZhou_All",               
                  "S1: SL.glm_All + S2: SL.gilleskie_All",              
                  "S1: SL.glm_All + S2: SL.rf.caret1_All",             
                  "S1: SL.glm_All + S2: SL.glmnet_All",             
                  "S1: SL.rf.caret1_All + S2: SL.logOLS.smear_All",     
                  "S1: SL.rf.caret1_All + S2: SL.gammaLogGLM_All",     
                  "S1: SL.rf.caret1_All + S2: SL.gammaIdentityGLM_All", 
                  "S1: SL.rf.caret1_All + S2: SL.manningGLM_All",       
                  "S1: SL.rf.caret1_All + S2: SL.gengamma_All",         
                  "S1: SL.rf.caret1_All + S2: SL.coxph_All",            
                  "S1: SL.rf.caret1_All + S2: SL.wangZhou_All",         
                  "S1: SL.rf.caret1_All + S2: SL.gilleskie_All",        
                  "S1: SL.rf.caret1_All + S2: SL.rf.caret1_All",        
                  "S1: SL.rf.caret1_All + S2: SL.glmnet_All",       
                  "S1: SL.glmnet_All + S2: SL.logOLS.smear_All",    
                  "S1: SL.glmnet_All + S2: SL.gammaLogGLM_All",     
                  "S1: SL.glmnet_All + S2: SL.gammaIdentityGLM_All",
                  "S1: SL.glmnet_All + S2: SL.manningGLM_All",      
                  "S1: SL.glmnet_All + S2: SL.gengamma_All",        
                  "S1: SL.glmnet_All + S2: SL.coxph_All",           
                  "S1: SL.glmnet_All + S2: SL.wangZhou_All",        
                  "S1: SL.glmnet_All + S2: SL.gilleskie_All",       
                  "S1: SL.glmnet_All + S2: SL.rf.caret1_All",       
                  "S1: SL.glmnet_All + S2: SL.glmnet_All",      
                  "Single: SL.mean_All",                                
                  "Single: SL.lm_All",                                  
                  "Single: SL.zip_All",                                 
                  "Single: SL.zinb_All",                                
                  "Single: SL.tobit_All",                               
                  "Single: SL.tweedie_All",                             
                  "Single: SL.rf.caret1_All",                           
                  "Single: SL.glmnet_All",
                  "One-stage SuperLearner",
                  "Two-stage SuperLearner",
                  "Discrete SuperLearner")
  algo.num <- length(algo.name)
  
  ## MSE
	mse_result <- matrix(NA, nrow = nrow(parameter_grid), ncol = algo.num)
	for(i in c(1:nrow(parameter_grid))){
		tmp <- get(load(paste0("/home/zwu56/sim/scratch/MSE_n=",parameter_grid$sample_size[i],
		                       "_zero=",parameter_grid$zero_inflate[i],
		                       "_nonzero=",parameter_grid$dist[i],
		                       "_interaction=",parameter_grid$interaction[i],
		                       "_linear=",parameter_grid$linear[i],
		                       "_seed=", parameter_grid$seed[i],".RData")))
		mse_result[i,] <- tmp
	}
	mse_result <- cbind(parameter_grid,mse_result)
	colnames(mse_result)[7:47] <- algo.name
	
  ## MAE
	mae_result <- matrix(NA, nrow = nrow(parameter_grid), ncol = algo.num)
	for(i in c(1:nrow(parameter_grid))){
	  tmp <- get(load(paste0("/home/zwu56/sim/scratch/MAE_n=",parameter_grid$sample_size[i],
	                         "_zero=",parameter_grid$zero_inflate[i],
	                         "_nonzero=",parameter_grid$dist[i],
	                         "_interaction=",parameter_grid$interaction[i],
	                         "_linear=",parameter_grid$linear[i],
	                         "_seed=", parameter_grid$seed[i],".RData")))
	  mae_result[i,] <- tmp
	}
	mae_result <- cbind(parameter_grid,mae_result)
	colnames(mae_result)[7:47] <- algo.name
	
  ## R square
	Rsq_result <- matrix(NA, nrow = nrow(parameter_grid), ncol = algo.num)
	for(i in c(1:nrow(parameter_grid))){
	  tmp <- get(load(paste0("/home/zwu56/sim/scratch/Rsq_n=",parameter_grid$sample_size[i],
	                         "_zero=",parameter_grid$zero_inflate[i],
	                         "_nonzero=",parameter_grid$dist[i],
	                         "_interaction=",parameter_grid$interaction[i],
	                         "_linear=",parameter_grid$linear[i],
	                         "_seed=", parameter_grid$seed[i],".RData")))
	  Rsq_result[i,] <- tmp
	}
	Rsq_result <- cbind(parameter_grid,Rsq_result)
	colnames(Rsq_result)[7:47] <- algo.name
	
  ## save result object
	save(mse_result, 
	     file=paste0("/home/zwu56/sim/MSE_result",".RData"))
	save(mae_result, 
	     file=paste0("/home/zwu56/sim/MAE_result",".RData"))
	save(Rsq_result, 
	     file=paste0("/home/zwu56/sim/Rsq_result",".RData"))
}
