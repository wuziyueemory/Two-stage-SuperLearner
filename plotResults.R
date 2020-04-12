#==============================================================#
# Code used to plot simulation results
#==============================================================#

# load libraries
library(gplots)
library(ggplot2)

# parameter grid
parameter_grid <- expand.grid(
  seed = 1:1000,
  zero_inflate = c(0.05, 0.7),
  dist = c("lognormal", "gamma", "tweedie", "mixture"),
  interaction = c(TRUE, FALSE),
  linear = c(TRUE),
  sample_size = c(500, 2000)
)

# load whole dataset
mse <- get(load(paste0("D:/Emory University/David Benkeser/collaboration/simulation/results_1000/MSE",".RData")))
mae <- get(load(paste0("D:/Emory University/David Benkeser/collaboration/simulation/results_1000/MAE",".RData")))
Rsq <- get(load(paste0("D:/Emory University/David Benkeser/collaboration/simulation/results_1000/Rsq",".RData")))

rownames(mse) <- seq(1,32000,1)
rownames(mae) <- seq(1,32000,1)
rownames(Rsq) <- seq(1,32000,1)

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

colnames(mse)[7:47] <- algo.name
colnames(mae)[7:47] <- algo.name
colnames(Rsq)[7:47] <- algo.name

#======================================================================#
# Summary table of MSE, MAE, R^2: mean & SD
#======================================================================#
mean_grid <- expand.grid(
  zero_inflate = c(0.05, 0.7),
  dist = c("lognormal", "gamma", "tweedie", "mixture"),
  interaction = c(TRUE, FALSE),
  linear = c(TRUE),
  sample_size = c(500, 2000)
)

mse.mean <- NULL
mae.mean <- NULL
Rsq.mean <- NULL

for (i in 1:32){
  mse.tmp = mse[seq(1000*(i-1)+1,1000*i),7:47]
  mae.tmp = mae[seq(1000*(i-1)+1,1000*i),7:47]
  Rsq.tmp = Rsq[seq(1000*(i-1)+1,1000*i),7:47]  
  mse.mean.tmp <- apply(mse.tmp,2,mean)
  mse.mean <- rbind(mse.mean,mse.mean.tmp)
  mae.mean.tmp <- apply(mae.tmp,2,mean)
  mae.mean <- rbind(mae.mean,mae.mean.tmp)
  Rsq.mean.tmp <- apply(Rsq.tmp,2,mean)
  Rsq.mean <- rbind(Rsq.mean,Rsq.mean.tmp)
}

mse_mean <- apply(mse.mean,2,mean)
mae_mean <- apply(mae.mean,2,mean)
Rsq_mean <- apply(Rsq.mean,2,mean)
mse_se <- apply(mse.mean,2,sd)
mae_se <- apply(mae.mean,2,sd)
Rsq_se <- apply(Rsq.mean,2,sd)

# Summary
summary <- data.frame("Algorithm"=algo.name,"MSE_mean"=mse_mean,"MSE_se"=mse_se,
                      "MAE_mean"=mae_mean,"MAE_se"=mae_se,
                      "Rsq_mean"=Rsq_mean,"Rsq_se"=Rsq_se)
summary$relativeMSE <- summary$MSE_mean/summary$MSE_mean[2]

#==================================================================================#
# Plot results - Aggregate 32 settings into one boxplot
#==================================================================================#
# MSE
mse.rank <- NULL
for (i in 1:32){
  mse.tmp = mse[seq(1000*(i-1)+1,1000*i),7:47]
  mse.mean.tmp <- apply(mse.tmp,2,mean)
  mse.rank.tmp <- rank(mse.mean.tmp)
  mse.rank <- rbind(mse.rank,mse.rank.tmp)
}
mse.rank <- cbind(mse.mean[,-c(47:49)],mse.rank)
# Extract MSE & Rank from dataset
mse.rank.mse <- mse.rank[,6:46]
mse.rank.rank <- mse.rank[,47:87]
# Add average MSE / Rank row for each dataset
mse.rank.mse <- rbind(mse.rank.mse,apply(mse.rank[,6:46],2,mean))
mse.rank.rank <- rbind(mse.rank.rank,apply(mse.rank[,47:87],2,mean))
# Reorder according to average MSE
mse.rank.mse.reorder <- mse.rank.mse[,order(mse.rank.mse[33,])]
mse.rank.rank.reorder <- mse.rank.rank[,order(mse.rank.rank[33,])]
# combine together
mse.rank.reorder <- cbind(mean_grid,mse.rank.mse.reorder[-33,],mse.rank.rank.reorder[-33,])

# Box plot - MSE rank
m.0 = as.matrix(mse.rank.reorder[,c(47:87)])
q.0 = c()
for (i in seq(1:ncol(m.0))){
  q.0 = c(q.0, m.0[,i])
}

mserank <- data.frame(
  Algorithm = factor(rep(1:41, each = 32), labels = colnames(mse.rank.reorder)[47:87]),
  MSE_rank = q.0
  )

bp.0 <- ggplot(mserank, aes(x=Algorithm, y=MSE_rank))+
  geom_boxplot()+
  ggtitle("MSE rank for 32 data settings")+
  labs(x="Algorithm",y="MSE rank")+
  scale_y_continuous(breaks=seq(0,40,5))+
  stat_summary(fun.y=mean, geom="point", shape=23, size=2)+
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
# Horizontal box plot
bp.0 + coord_flip()

# Boxplot - Raw MSE
m.1 = as.matrix(mse.rank.reorder[,c(6:46)])
q.1 = c()
for (i in seq(1:ncol(m.1))){
  q.1 = c(q.1, m.1[,i])
}

mseraw <- data.frame(
  Algorithm = factor(rep(1:41, each = 32), labels = colnames(mse.rank.reorder)[6:46]),
  MSE_raw = q.1
)

bp.1 <- ggplot(mseraw, aes(x=Algorithm, y=MSE_raw))+
  geom_boxplot()+
  ggtitle("Raw MSE for 32 data settings")+
  labs(x="Algorithm",y="MSE")+
  stat_summary(fun.y=mean, geom="point", shape=23, size=2)+
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
# Horizontal box plot
bp.1 + coord_flip()

# Relative MSE - log scale (S1: glm + S2: gammaLogGLM as reference)
for (i in 1:41){
  mse.rank.reorder[,87+i] = mse.rank.reorder[,5+i]/mse.rank.reorder[,13]
}
colnames(mse.rank.reorder)[88:128] <- colnames(mse.rank.reorder)[6:46]
# order algorithm by their geometric mean
sort(apply(mse.rank.reorder[,c(88:128)],2,function(x) exp(mean(log(x)))))
# Extract relative MSE from dataset
mse.rank.mserelative <- mse.rank.reorder[,88:128]
# Add geometric mean MSE relative row for each dataset
mse.rank.mserelative <- rbind(mse.rank.mserelative,
                              apply(mse.rank.mserelative,2,function(x) exp(mean(log(x)))))
# Reorder according to average MSE
mse.rank.mserelative.reorder <- mse.rank.mserelative[,order(mse.rank.mserelative[33,])]
mse.rank.mserelative.reorder <- mse.rank.mserelative.reorder[-33,]

# Box plot - Relative MSE (log scale)
m.2 = as.matrix(mse.rank.mserelative.reorder)
q.2 = c()
for (i in seq(1:ncol(m.2))){
  q.2 = c(q.2, log(m.2[,i]))
}

mserelative <- data.frame(
  Algorithm = factor(rep(1:41, each = 32), labels = colnames(mse.rank.mserelative.reorder)),
  Log_MSE_relative = q.2
)
                                  
bp.2 <- ggplot(mserelative, aes(x=Algorithm, y=Log_MSE_relative))+
  geom_boxplot()+
  ggtitle("Log relative MSE for 32 data settings")+
  labs(x="Algorithm",y="Log relative MSE ")+
  scale_y_continuous(breaks=seq(-1,4,1))+
  stat_summary(fun.y=mean, geom="point", shape=23, size=2)+
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
# Horizontal box plot
bp.2 + coord_flip()

# R^2
Rsq.mean <- NULL
for (i in 1:32){
  Rsq.tmp = Rsq[seq(1000*(i-1)+1,1000*i),7:47]
  Rsq.mean.tmp <- apply(Rsq.tmp,2,mean)
  Rsq.mean <- rbind(Rsq.mean,Rsq.mean.tmp)
}
# Add average Rsq row for each dataset
Rsq.mean <- rbind(Rsq.mean,apply(Rsq.mean,2,mean))
# Reorder according to average MSE
Rsq.mean.reorder <- Rsq.mean[,order(Rsq.mean[33,],decreasing = T)]
Rsq.mean.reorder <- Rsq.mean.reorder[-33,]

# Box plot - R^2
m.3 = as.matrix(Rsq.mean.reorder)
q.3 = c()
for (i in seq(1:ncol(m.3))){
  q.3 = c(q.3, m.3[,i])
}

Rsqmean <- data.frame(
  Algorithm = factor(rep(1:41, each = 32), labels = colnames(Rsq.mean.reorder)),
  Rsq = q.3
)

bp.3 <- ggplot(Rsqmean, aes(x=Algorithm, y=Rsq))+
  geom_boxplot()+
  ggtitle("Boxplot of MSE rank for 32 data settings")+
  scale_y_continuous(breaks=seq(-1.5,1,0.1))+
  stat_summary(fun.y=mean, geom="point", shape=23, size=2)+
  theme(legend.position = 'none',
        # Remove panel border
        panel.border = element_blank(),  
        # Remove panel grid lines
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Remove panel background
        panel.background = element_blank(),
        # Add axis line
        axis.line = element_line(colour = "grey"))+
  ylim(-6,2)
# Horizontal box plot
bp.3 + coord_flip()

#=======================================================================================#
# Plot results - Sub-Boxplots grouped by 3 different data generating parameters
#=======================================================================================#

# Boxplot - small sample size (n=500)
# Raw MSE
mse.small <- mse.rank.reorder[mse.rank.reorder$sample_size==500,c(6:46)]
# Add average MSE row 
mse.small <- rbind(mse.small,apply(mse.small,2,mean))
# Reorder according to average MSE
mse.small <- mse.small[,order(mse.small[17,])]
mse.small <- mse.small[-17,]

m.500 = as.matrix(mse.small)
q.500 = c()
for (i in seq(1:ncol(m.500))){
  q.500 = c(q.500, m.500[,i])
}

mse_small <- data.frame(
  Algorithm = factor(rep(1:41, each = 16), labels = colnames(mse.small)),
  MSE = q.500
)

bp.500 <- ggplot(mse_small, aes(x=Algorithm, y=MSE))+
  geom_boxplot()+
  ggtitle("(a) MSE - n=500")+
  stat_summary(fun.y=mean, geom="point", shape=23, size=2)+
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
# Horizontal box plot
bp.500 <- bp.500 + coord_flip()

# Boxplot - large sample size (n=2000)
# Raw MSE
mse.large <- mse.rank.reorder[mse.rank.reorder$sample_size==2000,c(6:46)]
# Add average MSE row 
mse.large <- rbind(mse.large,apply(mse.large,2,mean))
# Reorder according to average MSE
mse.large <- mse.large[,order(mse.large[17,])]
mse.large <- mse.large[-17,]

m.2000 = as.matrix(mse.large)
q.2000 = c()
for (i in seq(1:ncol(m.2000))){
  q.2000 = c(q.2000, m.2000[,i])
}

mse_large <- data.frame(
  Algorithm = factor(rep(1:41, each = 16), labels = colnames(mse.large)),
  MSE = q.2000
)

bp.2000 <- ggplot(mse_large, aes(x=Algorithm, y=MSE))+
  geom_boxplot()+
  ggtitle("(b) MSE - n=2000")+
  stat_summary(fun.y=mean, geom="point", shape=23, size=2)+
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
# Horizontal box plot
bp.2000 <- bp.2000 + coord_flip()

library(gridExtra)
grid.arrange(bp.500, bp.2000, ncol=2)

# Boxplot - low zero percent (5%)
# Raw MSE
mse.low <- mse.rank.reorder[mse.rank.reorder$zero_inflate==0.05,c(6:46)]
# Add average MSE row 
mse.low <- rbind(mse.low,apply(mse.low,2,mean))
# Reorder according to average MSE
mse.low <- mse.low[,order(mse.low[17,])]
mse.low <- mse.low[-17,]

m.low = as.matrix(mse.low)
q.low = c()
for (i in seq(1:ncol(m.low))){
  q.low = c(q.low, m.low[,i])
}

mse_low <- data.frame(
  Algorithm = factor(rep(1:41, each = 16), labels = colnames(mse.low)),
  MSE = q.low
)

bp.low <- ggplot(mse_low, aes(x=Algorithm, y=MSE))+
  geom_boxplot()+
  ggtitle("(c) MSE - 5% zero")+
  stat_summary(fun.y=mean, geom="point", shape=23, size=2)+
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
# Horizontal box plot
bp.low <- bp.low + coord_flip()

# Boxplot - high zero percentage (n=70%)
# Raw MSE
mse.high <- mse.rank.reorder[mse.rank.reorder==0.7,c(6:46)]
# Add average MSE row 
mse.high <- rbind(mse.high,apply(mse.high,2,mean))
# Reorder according to average MSE
mse.high <- mse.high[,order(mse.high[17,])]
mse.high <- mse.high[-17,]

m.high = as.matrix(mse.high)
q.high = c()
for (i in seq(1:ncol(m.high))){
  q.high = c(q.high, m.high[,i])
}

mse_high <- data.frame(
  Algorithm = factor(rep(1:41, each = 16), labels = colnames(mse.high)),
  MSE = q.high
)

# Box plot
bp.high <- ggplot(mse_high, aes(x=Algorithm, y=MSE))+
  geom_boxplot()+
  ggtitle("(d) MSE - 70% zero")+
  stat_summary(fun.y=mean, geom="point", shape=23, size=2)+
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
# Horizontal box plot
bp.high <- bp.high + coord_flip()

library(gridExtra)
grid.arrange(bp.low, bp.high, ncol=2)

# Boxplot - lognormal
# Raw MSE
mse.lognormal <- mse.rank.reorder[mse.rank.reorder$dist=="lognormal",c(6:46)]
# Add average MSE row 
mse.lognormal <- rbind(mse.lognormal,apply(mse.lognormal,2,mean))
# Reorder according to average MSE
mse.lognormal <- mse.lognormal[,order(mse.lognormal[9,])]
mse.lognormal <- mse.lognormal[-9,]

m.lognormal = as.matrix(mse.lognormal)
q.lognormal = c()
for (i in seq(1:ncol(m.lognormal))){
  q.lognormal = c(q.lognormal, m.lognormal[,i])
}

mse_lognormal <- data.frame(
  Algorithm = factor(rep(1:41, each = 8), labels = colnames(mse.lognormal)),
  MSE = q.lognormal
)
# restrict to top 10 algorithm
mse_lognormal.subset <- mse_lognormal[c(1:80),]

bp.lognormal <- ggplot(mse_lognormal.subset, aes(x=Algorithm, y=MSE))+
  geom_boxplot()+
  ggtitle("(a) MSE - Lognormal")+
  stat_summary(fun.y=mean, geom="point", shape=23, size=2)+
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
# Horizontal box plot
bp.lognormal <- bp.lognormal + coord_flip()

# Boxplot - gamma
# Raw MSE
mse.gamma <- mse.rank.reorder[mse.rank.reorder$dist=="gamma",c(6:46)]
# Add average MSE row 
mse.gamma <- rbind(mse.gamma,apply(mse.gamma,2,mean))
# Reorder according to average MSE
mse.gamma <- mse.gamma[,order(mse.gamma[9,])]
mse.gamma <- mse.gamma[-9,]

m.gamma = as.matrix(mse.gamma)
q.gamma = c()
for (i in seq(1:ncol(m.gamma))){
  q.gamma = c(q.gamma, m.gamma[,i])
}

mse_gamma <- data.frame(
  Algorithm = factor(rep(1:41, each = 8), labels = colnames(mse.gamma)),
  MSE = q.gamma
)
# restrict to top 10 algorithm
mse_gamma.subset <- mse_gamma[c(1:80),]

bp.gamma <- ggplot(mse_gamma.subset, aes(x=Algorithm, y=MSE))+
  geom_boxplot()+
  ggtitle("(b) MSE - Gamma")+
  stat_summary(fun.y=mean, geom="point", shape=23, size=2)+
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
# Horizontal box plot
bp.gamma <- bp.gamma + coord_flip()

# Boxplot - Tweedie
# Raw MSE
mse.tweedie <- mse.rank.reorder[mse.rank.reorder$dist=="tweedie",c(6:46)]
# Add average MSE row 
mse.tweedie <- rbind(mse.tweedie,apply(mse.tweedie,2,mean))
# Reorder according to average MSE
mse.tweedie <- mse.tweedie[,order(mse.tweedie[9,])]
mse.tweedie <- mse.tweedie[-9,]

m.tweedie = as.matrix(mse.tweedie)
q.tweedie = c()
for (i in seq(1:ncol(m.tweedie))){
  q.tweedie = c(q.tweedie, m.tweedie[,i])
}

mse_tweedie <- data.frame(
  Algorithm = factor(rep(1:41, each = 8), labels = colnames(mse.tweedie)),
  MSE = q.tweedie
)
# restrict to top 10 algorithm + tweedie
mse_tweedie.subset <- mse_tweedie[c(1:80,137:144),]

# Box plot
bp.tweedie <- ggplot(mse_tweedie.subset, aes(x=Algorithm, y=MSE))+
  geom_boxplot()+
  ggtitle("(c) MSE - Tweedie")+
  stat_summary(fun.y=mean, geom="point", shape=23, size=2)+
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
# Horizontal box plot
bp.tweedie <- bp.tweedie + coord_flip()

# Boxplot - mixture
# Raw MSE
mse.mixture <- mse.rank.reorder[mse.rank.reorder$dist=="mixture",c(6:46)]
# Add average MSE row 
mse.mixture <- rbind(mse.mixture,apply(mse.mixture,2,mean))
# Reorder according to average MSE
mse.mixture <- mse.mixture[,order(mse.mixture[9,])]
mse.mixture <- mse.mixture[-9,]

m.mixture = as.matrix(mse.mixture)
q.mixture = c()
for (i in seq(1:ncol(m.mixture))){
  q.mixture = c(q.mixture, m.mixture[,i])
}

mse_mixture <- data.frame(
  Algorithm = factor(rep(1:41, each = 8), labels = colnames(mse.mixture)),
  MSE = q.mixture
)
# restrict to top 10 algorithm
mse_mixture.subset <- mse_mixture[c(1:80),]

bp.mixture <- ggplot(mse_mixture.subset, aes(x=Algorithm, y=MSE))+
  geom_boxplot()+
  ggtitle("(d) MSE - Mixture")+
  stat_summary(fun.y=mean, geom="point", shape=23, size=2)+
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
# Horizontal box plot
bp.mixture <- bp.mixture + coord_flip()

library(gridExtra)
grid.arrange(bp.lognormal, bp.gamma, bp.tweedie, bp.mixture,nrow=2, ncol=2)

                                   
