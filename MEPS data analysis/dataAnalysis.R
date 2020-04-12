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
source(/home/zwu56/sim/twostageSL.R)
source(/home/zwu56/sim/Estimators.R)

# Loading MEPS train_data
h202 = read.xport("/home/zwu56/read data analysis/h202.ssp")

# Restrict to only adults (age > 18 years)
data <- h202[h202$AGE1X >= 18,]

# Prepare for training dataset
# Include a small part of variables
# Age, Gender, Race, Education, Income as percent poverty, Health insurance variables, Region of the country,
# Diabetes, Hypertension, Cancer, Heart disease, Total expenditures in 2016
train_dat = data %>%
  select(DUPERSID, AGE1X, SEX, RACETHX, EDUCYR, POVLEVY1, INSCOVY1, MCREVY1, MCDEVY1, OPAEVY1, OPBEVY1, REGIONY1,
         DIABDXY1,  HIBPDXY1,  CANCERY1, MIDXY1, OHRTDXY1, CHDDXY1, TOTEXPY1)

# Continuous variables: Age, education, Income as percent poverty, Total expenditures in 2016
# Discard records with continuous variables <0
train_dat <- train_dat[train_dat$EDUCYR>=0,]
train_dat <- train_dat[train_dat$POVLEVY1>=0,]
train_dat <- train_dat[train_dat$REGIONY1>=0,]
train_dat <- train_dat[train_dat$HIBPDXY1>=0,]
train_dat <- train_dat[train_dat$CANCERY1>=0,]

# create dummy  variables for categorical variables
# SEX (2 levels)
train_dat$SEX_dum = train_dat$SEX - 1
# RACE (5 levels): using WHITE as the reference
train_dat$RACETHX_HISPANIC = as.numeric(train_dat$RACETHX == 1)
train_dat$RACETHX_BLACK = as.numeric(train_dat$RACETHX == 3)
train_dat$RACETHX_ASIAN = as.numeric(train_dat$RACETHX == 4)
train_dat$RACETHX_OTHER = as.numeric(train_dat$RACETHX == 5)
# REGION (4 levels): using SOUTH as the reference
train_dat$REGION_NE = as.numeric(train_dat$REGIONY1 == 1)
train_dat$REGION_MW = as.numeric(train_dat$REGIONY1 == 2)
train_dat$REGION_W = as.numeric(train_dat$REGIONY1 == 4)
# DIABETES (2 levels)
train_dat$DIABDX_dum <- 2-train_dat$DIABDXY1
# HYPERTENSION (2 levels)
train_dat$HIBPDX_dum <- 2-train_dat$HIBPDXY1
# CANCER (2 levels)
train_dat$CANCER_dum <- 2-train_dat$CANCERY1

# create new variables for train_dataset
# health insurance variables
# private
train_dat$insure1 = 0
train_dat$insure1[train_dat$INSCOVY1==1] = 1
# medicare
train_dat$insure2 = 0
train_dat$insure2[train_dat$MCREVY1==1] = 1
# public
train_dat$insure3 = 0
train_dat$insure3[train_dat$MCDEVY1==1 | train_dat$OPAEVY1==1 | train_dat$OPBEVY1==1] = 1
# uninsured
train_dat$uninsured = 0
train_dat$uninsured[train_dat$INSCOVY1==3]=1
# heart disease variable
train_dat$heart_disease = 0
train_dat$heart_disease[train_dat$MIDXY1==1 | train_dat$OHRTDXY1==1 | train_dat$CHDDXY1==1] = 1

# create train data subset to include target predictors 
train <- train_dat[,c(1,2,5,6,20:35,19)]
colnames(train)[4] <- 'POVLEV'
colnames(train)[21] <- 'TOTEXP'


# Create test data in a similar way
test_dat = data %>%
  select(DUPERSID, AGE1X, SEX, RACETHX, EDUCYR, POVLEVY2, INSCOVY2, MCREVY2, MCDEVY2, OPAEVY2, OPBEVY2, REGIONY2,
         DIABDXY2,  HIBPDXY2,  CANCERY2, MIDXY2, OHRTDXY2, CHDDXY2, TOTEXPY2)

# Discard records with continuous variables <0
test_dat <- test_dat[test_dat$EDUCYR>=0,]
test_dat <- test_dat[test_dat$POVLEVY2>=0,]
test_dat <- test_dat[test_dat$REGIONY2>=0,]
test_dat <- test_dat[test_dat$HIBPDXY2>=0,]
test_dat <- test_dat[test_dat$CANCERY2>=0,]
test_dat <- test_dat[test_dat$TOTEXPY2>=0,]

# create dummy  variables for categorical variables: SEX, RACETHX
# SEX (2 levels)
test_dat$SEX_dum = test_dat$SEX - 1
# RACE (5 levels): using WHITE as the reference
test_dat$RACETHX_HISPANIC = as.numeric(test_dat$RACETHX == 1)
test_dat$RACETHX_BLACK = as.numeric(test_dat$RACETHX == 3)
test_dat$RACETHX_ASIAN = as.numeric(test_dat$RACETHX == 4)
test_dat$RACETHX_OTHER = as.numeric(test_dat$RACETHX == 5)
# REGION (4 levels): using SOUTH as the reference
test_dat$REGION_NE = as.numeric(test_dat$REGIONY2 == 1)
test_dat$REGION_MW = as.numeric(test_dat$REGIONY2 == 2)
test_dat$REGION_W = as.numeric(test_dat$REGIONY2 == 4)
# DIABETES (2 levels)
test_dat$DIABDX_dum <- 2-test_dat$DIABDXY2
# HYPERTENSION (2 levels)
test_dat$HIBPDX_dum <- 2-test_dat$HIBPDXY2
# CANCER (2 levels)
test_dat$CANCER_dum <- 2-test_dat$CANCERY2

# create new variables for test_dataset
# health insurance variables
# private
test_dat$insure1 = 0
test_dat$insure1[test_dat$INSCOVY2==1] = 1
# medicare
test_dat$insure2 = 0
test_dat$insure2[test_dat$MCREVY2==1] = 1
# public
test_dat$insure3 = 0
test_dat$insure3[test_dat$MCDEVY2==1 | test_dat$OPAEVY2==1 | test_dat$OPBEVY2==1] = 1
# uninsured
test_dat$uninsured = 0
test_dat$uninsured[test_dat$INSCOVY2==3]=1
# heart disease variable
test_dat$heart_disease = 0
test_dat$heart_disease[test_dat$MIDXY2==1 | test_dat$OHRTDXY2==1 | test_dat$CHDDXY2==1] = 1

# create train data subset to include target predictors 
test <- test_dat[,c(1,2,5,6,20:35,19)]
colnames(test)[4] <- 'POVLEV'
colnames(test)[21] <- 'TOTEXP'

#==================================================================================#
# Exploratory data analysis
#==================================================================================#
# Outcome distribution
# train
quantile(train$TOTEXP, c(.05, .25, .50, .75, .95))
mean(train$TOTEXP)
sqrt(var(train$TOTEXP))
skewness((train$TOTEXP))
# test
quantile(test$TOTEXP, c(.05, .25, .50, .75, .95))
mean(test$TOTEXP)
sqrt(var(test$TOTEXP))
skewness((test$TOTEXP))

# descriptive statistics of covariates
# train
mean(train$AGE1X)
sqrt(var(train$AGE1X))
table(train_dat$SEX)/nrow(train_dat)*100
table(train_dat$RACETHX)/nrow(train_dat)*100
table(train_dat$REGIONY1)/nrow(train_dat)*100
mean(train$EDUCYR)
mean(train$POVLEV)
mean(train$insure1)*100
mean(train$insure2)*100
mean(train$insure3)*100
mean(train$uninsured)*100
mean(train$DIABDX_dum)*100
mean(train$HIBPDX_dum)*100
mean(train$CANCER_dum)*100
mean(train$heart_disease)*100
# test
mean(test$AGE1X)
sqrt(var(test$AGE1X))
table(test_dat$SEX)/nrow(test_dat)*100
table(test_dat$RACETHX)/nrow(test_dat)*100
table(test_dat$REGIONY2)/nrow(test_dat)*100
mean(test$EDUCYR)
mean(test$POVLEV)
mean(test$insure1)*100
mean(test$insure2)*100
mean(test$insure3)*100
mean(test$uninsured)*100
mean(test$DIABDX_dum)*100
mean(test$HIBPDX_dum)*100
mean(test$CANCER_dum)*100
mean(test$heart_disease)*100
# zero percentage
sum(train$TOTEXP==0)/nrow(train)*100
sum(test$TOTEXP==0)/nrow(test)*100
# extreme percentage
sum(train$TOTEXP>50000)/nrow(train)*100
sum(test$TOTEXP>50000)/nrow(test)*100
max(train$TOTEXP)
max(test$TOTEXP)

# histogram
par(mfrow=c(1,2))
hist(train$TOTEXP,breaks = 200,lty="blank",main = "Total healthcare expenditures - 2016",
     xlab = "Total healthcare expenditures (dollars)",ylab="Percentage",yaxt='n',
     col=rgb(0,0,1,1/2),ylim = c(0,7000))
axis(side=2, at=seq(0,7000,1000),labels=seq(0,0.7,0.1))

hist(test$TOTEXP,breaks = 200,lty="blank",main = "Total healthcare expenditures - 2017",
     xlab = "Total healthcare expenditures (dollars)",ylab="Percentage",yaxt='n',
     col=rgb(0,0,1,1/2),ylim = c(0,7000))
axis(side=2, at=seq(0,7000,1000),labels=seq(0,0.7,0.1))

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

# MSE
mse_result <- data.frame("Algorithm"=algo.name,"MSE"=mse)
# MAE
mae_result <- data.frame("Algorithm"=algo.name,"MSE"=mae)
# R square
Rsq_result <- data.frame("Algorithm"=algo.name,"MSE"=Rsq)

# order MSE ascendingly
mse_reorder <- mse_result[order(mse_result$MSE),]
rownames(mse_reorder) <- seq(1,41,1)
mse_result$relativeMSE <- mse_result$MSE/mse_result$MSE[2]
# order MAE ascendingly
mae_reorder <- mae_result[order(mae_result$MSE),]
# order Rsq ascendingly
Rsq_result$RE <- Rsq_result$R_square/Rsq_result$R_square[40]
Rsq_reorder <- Rsq_result[order(Rsq_result$R_square),]
rownames(Rsq_reorder) <- seq(1,41,1)

# plot For MSE
mse_reorder$Algorithm <- factor(mse_reorder$Algorithm, levels = mse_reorder$Algorithm[order(mse_reorder$MSE)])
mse_plot <- ggplot(mse_reorder, aes(x=reorder(Algorithm, -MSE), y=MSE))+
  geom_point()+
  ggtitle('MSE - 2017 MEPS data')+
  labs(x="Algorithm",y="MSE")+
  scale_y_continuous(breaks=seq(200000000,260000000,5000000))+
  theme(legend.position = 'none',
        # Remove panel border
        panel.border = element_blank(),  
        # Remove panel grid lines
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Remove panel background
        panel.background = element_blank(),
        # Add axis line
        axis.line = element_line(colour = "grey"))
mse_plot + coord_flip() 

# plot For relative MSE
mse_reorder$relative_MSE <- mse_reorder$MSE / mse_reorder$MSE[15]
mse_relative_plot <- ggplot(mse_reorder, aes(x=reorder(Algorithm, -relative_MSE), y=relative_MSE))+
  geom_point()+
  ggtitle('Relative MSE - 2017 MEPS data')+
  labs(x="Algorithm",y="Relative MSE")+
  scale_y_continuous(breaks=seq(0.94,1.14,0.02))+
  theme(legend.position = 'none',
        # Remove panel border
        panel.border = element_blank(),  
        # Remove panel grid lines
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Remove panel background
        panel.background = element_blank(),
        # Add axis line
        axis.line = element_line(colour = "grey"))
mse_relative_plot + coord_flip()

# plot For Rsq & RE
Rsq_reorder$Algorithm <- factor(Rsq_reorder$Algorithm, levels = Rsq_reorder$Algorithm[order(Rsq_reorder$R_square)])
Rsq_1 <- Rsq_reorder[,c(1,3)]
colnames(Rsq_1) <- c("Algorithm","Value")
Rsq_2 <- Rsq_reorder[,c(1,4)]
colnames(Rsq_2) <- c("Algorithm","Value")
Rsq <- rbind(Rsq_1,Rsq_2)
rownames(Rsq) <- seq(1,82,1)
Rsq$metric <- c(rep("R^2",41),rep("RE",41))
mycolors <- c("R^2"="blue", "RE"="red")
tit <- expression(paste("R"^"2"," & RE - 2017 MEPS data"))
Rsq_plot <- ggplot(Rsq, aes(x=Algorithm, y=Value, group=metric,color=metric))+
  geom_point()+
  ggtitle(tit)+
  scale_y_continuous(name=bquote(R^2), sec.axis = sec_axis(~ 1*., name="RE"))+
  scale_color_manual(name="Metric", values = mycolors) +
  theme(legend.position = c(0.95, 0.2),
        legend.justification = c("right", "top"),
        # Remove panel border
        panel.border = element_blank(),  
        # Remove panel grid lines
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Remove panel background
        panel.background = element_blank(),
        # Add axis line
        axis.line = element_line(colour = "grey"),
        axis.title.y = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.title.y.right = element_text(color = "black"),
        axis.text.y.right = element_text(color = "black"))
Rsq_plot + coord_flip()
