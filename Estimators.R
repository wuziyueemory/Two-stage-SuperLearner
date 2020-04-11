#===================================================================================#
# Estimators used in the simulation studies and MEPS data analysis
#===================================================================================#

# library relevant packages
library(SuperLearner)
library(earth)
library(pscl) # zero-inflated-poisson
library(VGAM) # tobit regression
library(cplm) # tweedie regression
library(caret) # random forest 
library(randomForest) # random forest
library(e1071) # random forest 
library(moments) # Adative GLM
library(flexsurv) # Accelerated Failure Time Models (AFT)
library(survival) # Accelerated Failure Time Models (AFT) 
library(quantreg) # Quantile regression
library(sandwich) # Adaptive Hazard method (Gilleskie & Mroz)


##################################### One-part model #######################################################

#==================================================#
#	Zero-Inflated Poisson Model
#==================================================#

SL.zip <- function(Y, X, newX, family, obsWeights, ...){
  if(family$family=="gaussian"){
    # round outcome Y to be interger
    Y.int <- round(Y)
    suppressWarnings(
    fit.zip <- zeroinfl(Y.int ~ . | ., data=X,weights = obsWeights)
    )
    pred <- predict(fit.zip, newdata=newX, type="response")
    fit <- list(object = fit.zip)
    class(fit) <- "SL.glm" # can use predict.SL.glm
    out <- list(pred=pred, fit=fit)
    return(out)
  }else{
    stop("SL.zip not written for binomial family")
  }
}

#==================================================#
# Zero-Inflated Negative Binomial Model
#==================================================#

SL.zinb <- function(Y, X, newX, family, obsWeights, ...){
  if(family$family=="gaussian"){
    # round outcome Y to be interger
    Y.int <- round(Y)
    suppressWarnings(
    fit.zinb <- zeroinfl(Y.int ~ . | ., data=X,weights = obsWeights,dist = "negbin")
    )
    pred <- predict(fit.zinb, newdata=newX, type="response")
    fit <- list(object = fit.zinb)
    class(fit) <- "SL.glm" # can use predict.SL.glm
    out <- list(pred=pred, fit=fit)
    return(out)
  }else{
    stop("SL.zinb not written for binomial family")
  }
}

#==================================================#
#	Tobit Model
#==================================================#

SL.tobit <- function(Y, X, newX, family, obsWeights, ...){
  if(family$family=="gaussian"){
    suppressWarnings(
    fit.tobit <- vglm(Y ~., tobit(Lower = 0,type.fitted = "censored"),data=X,maxit=100)
    )
    pred <- predict(fit.tobit, newdata=newX, type="response")
    # in case generate negative prediction
    pred[pred<0]=0
    fit <- list(object = fit.tobit)
    class(fit) <- "SL.tobit"
    out <- list(pred=pred, fit=fit)
    return(out)
  }else{
    stop("SL.tobit not written for binomial family")
  }
}

predict.SL.tobit <- function(object, newdata, ...) {
  # newdata must be a dataframe, not a matrix.
  if (is.matrix(newdata)) {
    newdata = as.data.frame(newdata)
  }
  pred <- predict(object = object$object, newdata = newdata, type = "response")
  # in case generate negative prediction
  pred[pred<0]=0
  pred
}

#==================================================#
#	Tweedie Model
#==================================================#

SL.tweedie <- function(Y, X, newX, family, obsWeights, ...){
  if(family$family=="gaussian"){
    # using optimizer bobyqa
    suppressWarnings(
    fit.tweedie <-  cpglm(Y~.,data=X,optimizer = "bobyqa")
    )
    pred <- predict(fit.tweedie, newdata=newX, type="response")
    fit <- list(object = fit.tweedie)
    class(fit) <- "SL.tweedie" 
    out <- list(pred=pred, fit=fit)
    return(out)
  }else{
    stop("SL.tweedie not written for binomial family")
  }
}

predict.SL.tweedie <- function(object, newdata, ...) {
  # newdata must be a dataframe, not a matrix.
  if (is.matrix(newdata)) {
    newdata = as.data.frame(newdata)
  }
  pred <- predict(object = object$object, newdata = newdata, type = "response")
  pred
}

#======================================================================#
# Modified version of SL.caret that prints less annoying GBM output
#======================================================================#
# cv.number = 10

SL.caret1 <- function(Y, X, newX, family, obsWeights, method = "rf", tuneLength = 3, 
                      trControl=trainControl(method = "cv", number = 10, verboseIter = FALSE),
                      metric,...) 
{
  if (length(unique(Y))>2){
    if(is.matrix(Y)) Y <- as.numeric(Y)
    metric <- "RMSE"
    if(method=="gbm"){
      suppressWarnings(
        fit.train <- caret::train(x = X, y = Y, weights = obsWeights, 
                                  metric = metric, method = method, 
                                  tuneLength = tuneLength, 
                                  trControl = trControl,verbose=FALSE)
      )
    }else{
      suppressWarnings(
        fit.train <- caret::train(x = X, y = Y, weights = obsWeights, 
                                  metric = metric, method = method, 
                                  tuneLength = tuneLength, 
                                  trControl = trControl)
      )
    }
    pred <- predict(fit.train, newdata = newX, type = "raw")
  }
  if (length(unique(Y))<=2) {
    metric <- "Accuracy"
    Y.f <- as.factor(Y)
    levels(Y.f) <- c("A0", "A1")
    fit.train <- caret::train(x = X, y = Y.f, weights = obsWeights,
                              metric = metric, method = method, 
                              tuneLength = tuneLength, 
                              trControl = trControl)
    pred <- predict(fit.train, newdata = newX, type = "prob")[,2]
  }
  fit <- list(object = fit.train)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.caret")
  return(out)
}

#==================================================#
# CV-Random Forest
#==================================================#
# tuneLength=3

SL.rf.caret1 <- function(...,method="rf",tuneLength=3){
  SL.caret1(...,method=method,tuneLength=tuneLength)
}


######################################### Two-part model #######################################################

# Part 2:

#=========================================================#
# log-OLS: OLS on ln(y) + smear retransformation
# GLM with Gaussian family and id link on log(Y) + Duan (1983) correction
#=========================================================#

SL.logOLS.smear <- function(Y, X, newX, family, obsWeights, ...){
  if(family$family=="gaussian"){
    logY <- log(Y)
    fit.logGLM <- glm(logY ~ ., data=X, family=family, weights=obsWeights)
    mu <- predict(fit.logGLM, type="response", newdata=X)
    resid <- logY - mu
    pred <- exp(predict(fit.logGLM, type="response",newdata=newX))*mean(exp(resid))
    fit <- list(object=fit.logGLM, mean(exp(resid)))
    class(fit) <- "SL.logOLS.smear"
  }else{
    stop("SL.logGLM.smear not written for binomial family")
  }
  out <- list(fit=fit, pred=pred)
  return(out)
}

# predict function for SL.logOLS.smear
predict.SL.logOLS.smear <- function(object, newdata, ...){
  mu <- predict(object$object, newdata=newdata, type="response")
  correction <- object[[2]]
  return(exp(mu)*correction) 
}

#=========================================================#
# GLM (Gamma distribution + log-link)
#=========================================================#

SL.gammaLogGLM <- function(Y, X, newX, family, obsWeights, ...){
  if(family$family=="gaussian"){
    fit.glm <- glm(Y ~ ., data=X, family=Gamma(link='log'), weights=obsWeights,
                   control=list(maxit=100))
    pred <- predict(fit.glm, newdata=newX, type="response")
    fit <- list(object = fit.glm)
    class(fit) <- "SL.glm" # can use predict.SL.glm
    out <- list(pred=pred, fit=fit)
    return(out)
  }else{
    stop("SL.logGLM not written for binomial family")
  }
}

#=========================================================#
# GLM (Gamma distribution + Identity-link)
#=========================================================#

SL.gammaIdentityGLM <- function(Y, X, newX, family, obsWeights,...){
  if(family$family=="gaussian"){
    fit.glm <- glm(Y ~ ., data=X, family=Gamma(link='identity'), 
                   weights=obsWeights,
                   control=list(maxit=100), start=c(mean(Y),rep(0,ncol(X))))
    pred <- predict(fit.glm, newdata=newX, type="response")
    fit <- list(object = fit.glm)
    class(fit) <- "SL.glm"
    out <- list(pred=pred, fit=fit)
    return(out)
  }else{
    stop("SL.gammaIdentityGLM not written for binomial family")
  }
}

#=========================================================#
# Adaptive GLM algorithm of Manning (2001)
#=========================================================#

SL.manningGLM <- function(Y, X, newX, family, obsWeights, 
                          kCut = 3, # kurtosis cutpoint
                          lambdaCut = c(0.5,1.5,2.5), # skew cutpoint
                          startNLS=0, # starting values for NLS?
                          ...){
  if(family$family=="gaussian"){
    require(moments)
    # first do ols on log scale
    logY <- log(Y)
    fit.logGLM <- glm(logY ~ ., data=X, family=family, weights=obsWeights)
    mu <- predict(fit.logGLM, type="response", newdata=X)
    resid <- logY - mu
    # check kurtosis of residuals
    k <- kurtosis(resid)
    # by default use these methods
    # some of the other GLMs are unstable and if they fail, this 
    # algorithm returns log OLS + smearning estimate
    pred <- exp(predict(fit.logGLM, type="response",newdata=newX))*mean(exp(resid))
    fit <- list(object=fit.logGLM, mean(exp(resid)))
    class(fit) <- "SL.logOLS.smear"
    try({
      if(k < kCut){
        # park test
        fit.initGLM <- glm(Y ~ ., data=X, weights=obsWeights, family="gaussian")
        muPark <- predict(fit.initGLM, type="response", newdata=X)
        resid2Park <- (Y - muPark)^2
        fit.parkGLM <- glm(resid2Park ~ muPark, family="gaussian")
        lambda1 <- fit.parkGLM$coef[2]
        # use nls
        if(lambda1 < lambdaCut[1]){
          xNames <- colnames(X)
          d <- length(xNames)
          bNames <- paste0("b",1:d)
          form <- apply(matrix(1:d), 1, function(i){
            paste(c(bNames[i],xNames[i]),collapse="*")
          })
          formula <- paste(form,collapse=" + ")
          try({
            fit.nls <- nls(as.formula(paste0("Y ~ exp(b0 +",formula,")")), data= data.frame(Y, X),
                           start=eval(parse(text=paste0(
                             "list(b0=0.5,",paste(paste0(bNames, "=", startNLS),collapse=","),")"))))
          })
          pred <- predict(fit.nls, newdata=newX)
          fit <- list(object=fit.nls)
          class(fit) <- "SL.manningGLM"
        }else if(lambda1 < lambdaCut[2] & lambda1 >= lambdaCut[1]){
          # use poisson glm
          fit.poisGLM <- suppressWarnings(
            glm(Y ~ ., data=X, weights=obsWeights, family="poisson",control=list(maxit=100))
          )
          pred <- predict(fit.poisGLM, newdata=newX, type="response")
          fit <- list(object=fit.poisGLM)
          class(fit) <- "SL.manningGLM"
        }else if(lambda1 < lambdaCut[3] & lambda1 >= lambdaCut[2]){
          # use gamma glm
          fit.gammaGLM <- glm(Y ~ ., data=X, weights=obsWeights, family=Gamma(link='log'),control=list(maxit=100))
          pred <- predict(fit.gammaGLM, newdata=newX, type="response")
          fit <- list(object=fit.gammaGLM)
          class(fit) <- "SL.manningGLM"
        }else if(lambda1 > lambdaCut[3]){
          # use inverse gaussian glm -- not very stable
          fit.invGaussianGLM <- glm(Y ~ ., data=X,weights=obsWeights, family=inverse.gaussian(link="log"),control=list(maxit=100))
          pred <- predict(fit.invGaussianGLM, newdata=newX, type="response")
          fit <- list(object=fit.invGaussianGLM)
          class(fit) <- "SL.manningGLM"
        }
      }
    }, silent=TRUE)
  }else{
    stop("SL.manningGLM doesn't work with binomial family.")
  }
  out <- list(pred = pred, fit=fit)
  return(out)
}

# predict function
predict.SL.manningGLM <- function(object, newdata,...){
  if(!is.list(object$object)){
    pred <- predict(object=object$object, newdata=newdata, type="response")
  }else{
    pred <- predict(object=object$object, newdata=newdata, type="response")
  }
  pred
}

#=========================================================#
# Accelerated Failure Time Models (AFT)
#=========================================================#

SL.flexsurvreg <- function(Y, X, newX, family, obsWeights,
                           dist="gengamma",...){
  require(flexsurv)
  if(family$family=="gaussian"){
    fit.flexSurv <- flexsurvreg(
      as.formula(paste0("Surv(Y, rep(1, length(Y))) ~", paste(colnames(X),collapse="+"))) ,
      data=X, dist=dist
    )
    pred <- predict.SL.flexsurvreg(object=list(object=fit.flexSurv), newdata=newX, type="mean")
    fit <- list(object=fit.flexSurv)
    class(fit) <- "SL.flexsurvreg"
    out <- list(fit=fit, pred=pred)
  }else{
    stop("SL.genGamma not implemented for binominal family")
  }
  out
}

#prediction
predict.SL.flexsurvreg <- function(object, newdata, type="mean", ...){
  # function to return survival probability based on flexsurv object
  .getSurv <- function(x, fit, thisnewdata){
    summary(fit, t=x, B=0, newdata=thisnewdata)[[1]][,2]
  }
  pred <- as.numeric(apply(matrix(1:nrow(newdata)), 1, function(i){
    upper <- Inf
    out <- NA; class(out) <- "try-error"
    # integrate can be finnicky, so for stability, we first try to integrate with
    # upper limit = Inf, but if that fails move to 1e8, which sometimes is able to 
    # provide a sane answer when upper limit=Inf fails. Keep trying smaller and smaller
    # values, but don't go smaller than 1e6. If you try, then it just returns a random 
    # number between 0 and 1e6, which prevents Super Learner from crashing. 
    while(class(out)=="try-error" & upper > 1e6){
      out <- try(integrate(.getSurv, fit=object$object, thisnewdata=newdata[i,],lower=0,upper=upper)$value, silent=TRUE)
      if(upper==Inf){
        upper <- 1e8
      }else{
        upper <- upper/2
      }
    }
    if(class(out)=="try-error"){
      warning("Unable to integrate survival function. Returning random number between 0 and 100k")
      out <- runif(1,0,100000)
    }
    out
  }))
  pred
}

#=================================================================#
# Accelerated Failure Time Models (AFT): Genrealized Gamma
#=================================================================#

SL.gengamma <- function(..., dist="gengamma"){
  SL.flexsurvreg(...,dist=dist)
}

#=================================================================#
# Cox Proportional Hazard 
#=================================================================#

SL.coxph  <- function(Y, X, newX, family, obsWeights,
                      dist="gengamma",...){
  if(family$family=="gaussian"){
    library(survival)
    fit.coxph <- coxph(Surv(Y,rep(1,length(Y)))~., data=X)
    fit <- list(object=fit.coxph)
    class(fit) <- "SL.coxph"
    pred <- predict.SL.coxph(object=list(object=fit.coxph), newdata=newX)
  }else{
    stop("SL.coxph not implemented for binominal family")
  }
  return(list(fit=fit,pred=pred))
}

# prediction
predict.SL.coxph <- function(object,newdata,type="mean",...){
  # use surv.fit to get survival estimate and because by default it uses
  # nelson-aalen hazard, easy to convert back to an estimate of the mean
  surv.fit <- survfit(object$object, newdata=newdata)
  pred <- colSums(
    diff(c(0,surv.fit$time))*rbind(
      rep(1,dim(surv.fit$surv)[2]),
      surv.fit$surv[1:(dim(surv.fit$surv)[1]-1),]
    )
  )
  pred
}

#=================================================================#
# Quantile regression method of Wang and Zhou (2009)
#=================================================================#

SL.wangZhou <- function(Y, X, newX, family, obsWeights, 
                        g="log", # transformation of Y
                        m=length(Y), # number of quantiles
                        c=0.2, # for calculating truncated mean
                        b=0.05,# for calculating truncated mean
                        ...){
  require(quantreg)
  if(family$family=="gaussian"){
    n <- length(Y)
    # calculate alpha_n for calculating truncated mean
    alpha <- c*n^(-1/(1+4*b))
    tau <- seq(alpha, 1-alpha, length=m)
    # transform Y
    if(g=="log"){
      thisY <- log(Y)
      ginv <- function(x){ exp(x) }
    }else{
      stop("SL.wangZhou only implemented for log transforms")
    }
    # get quantile regressions
    suppressWarnings(
      fm <- rq(formula=as.formula("thisY~."), tau=tau, weights=obsWeights,
               data=data.frame(thisY,X))
    )
    QhatList <- predict(fm, newdata=newX, stepfun=TRUE, type="Qhat")
    QhatRearrange <- lapply(QhatList, rearrange)
    # transform to means
    pred <- unlist(lapply(QhatRearrange, FUN=function(Q){
      Qw <- ginv(environment(Q)$y[-which(duplicated(environment(Q)$x))])
      1/(1-2*alpha) * sum(Qw * diff(c(0,tau)))
    }))
  }else{
    stop("SL.wangZhou not written for binomial family")
  }
  fit <- list(object=fm, alpha=alpha, ginv=ginv)
  class(fit) <- "SL.wangZhou"
  out <- list(pred=pred, fit=fit)
  return(out)
}

# predict function for SL.wangZhou
predict.SL.wangZhou <- function(object, newdata, ...){
  require(quantreg)
  QhatList <- predict(object$object, newdata=newdata, stepfun=TRUE, type="Qhat")
  QhatRearrange <- lapply(QhatList, rearrange)
  pred <- mapply(Q=QhatRearrange, dt=diff(c(0,object$object$tau)), function(Q,dt){
    Qw <- do.call(object$ginv,args=list(x=environment(Q)$y[-which(duplicated(environment(Q)$x))]))
    1/(1-2*object$alpha) * sum(Qw * dt)
  })    
  pred
}

#==========================================================================================#
# Discrete conditional density estimator /  Adaptive Hazard method (Gilleskie & Mroz)
#==========================================================================================#

SL.gilleskie <- function(Y, X, newX, family, obsWeights,
                         kValues=c(5,15,25), # number of intervals
                         yBoundaries, # boundaries on the y-variable
                         maxPoly=2, # maximum polynomial in hazard regressions
                         ...){
  # need the sandwich package for covariate selection algorithm
  library(sandwich)
  # maxPoly describes the polynomial used for the partition variable
  # in the hazard regression. Choosing the number of partitions to be
  # less than this value leads to rank definiciency in glm()
  if(any(kValues < maxPoly)){ 
    warning("kValue specified that is less than maxPoly. These kValues will be ignored")
    kValues <- kValues[kValues>maxPoly]
  }
  
  #====================================================
  # get hazard fit over different partitions of data
  #====================================================
  
  outList <- lapply(split(kValues, 1:length(kValues)), FUN=function(K){
    # break up Y into K+1 partitions
    Ytilde <- cut(Y, breaks=quantile(Y, p=seq(0,1,length=K+1)), labels=FALSE,
                  include.lowest=TRUE)
    # make a long versions data set
    longDat <- data.frame(Ytilde,X,id=1:length(Y))[rep(1:length(Y),Ytilde),]
    # assign parition number variable
    row.names(longDat)[row.names(longDat) %in% paste(row.names(data.frame(Ytilde,X)))] <- paste(row.names(data.frame(Ytilde,X)),".0",sep="")  
    longDat$k <- as.numeric(paste(unlist(strsplit(row.names(longDat),".",fixed=T))[seq(2,nrow(longDat)*2,2)]))+1
    # indicator of falling in a particular partition
    longDat$indY <- as.numeric(longDat$k==longDat$Ytilde)
    # loop to do covariate selection
    pVal <- Inf
    d <- maxPoly
    while(pVal > 0.05 & d>=1){
      # generate the regression equation
      rhs <- NULL
      for(i in 1:(ncol(X)-1)){
        rhs <- c(rhs, paste0("poly(",colnames(X)[i],",",ifelse(length(unique(X[,i]))>d, d, length(unique(X[,i]))-1),")*poly(k,",d,")*",colnames(X)[(i+1):(ncol(X))],collapse="+"))
      }
      rhs <- c(rhs, paste0("poly(",colnames(X)[ncol(X)],",",ifelse(length(unique(X[,i]))>d, d, length(unique(X[,i]))-1),")*poly(k,",d,")"))
      # fit the hazard regression
      suppressWarnings(
        fm <- glm(as.formula(paste0("indY ~ ",paste0(rhs,collapse="+"))),
                  data=longDat, family="binomial")
      )
      # get coefficients of degree d
      dropNum <- NULL
      for(cn in colnames(X)){
        dropNum <- c(dropNum, grep(paste0(cn,", ",d,")",d), names(fm$coef[!is.na(fm$coef)])))
      }
      dropCoef <- fm$coef[!is.na(fm$coef)][dropNum]
      # get covariance matrix for all thos ecoefficients
      fullCov <- vcovHC(fm,type="HC0")
      dropCov <- fullCov[dropNum,dropNum]
      # test significance of those coefficients
      chi2Stat <- tryCatch({
        t(dropCoef)%*%solve(dropCov)%*%dropCoef
      },error=function(e){ return(0) }
      )
      pVal <- pchisq(chi2Stat, lower.tail=FALSE, df=length(dropCoef))
      d <- d-1
    }
    # after done dropping polynomial terms, get hazard predictions
    suppressWarnings(
      longDat$haz <- predict(fm, newdata=longDat,type="response")
    )
    # calculate likelihood
    tmp <- by(longDat, factor(longDat$id), FUN=function(x){
      prod((c(1,cumprod(1-x$haz[x$k < x$Ytilde])) * x$haz)^x$indY)
    })
    LKR <- sum(log(as.numeric(tmp))) + length(Y)*log(K)
    # return the likelihood ratio and estimate hazard regression
    return(list(LKR=LKR, fm=fm))
  })
  # figure out which one had highest likelihood ratio
  LKRs <- unlist(lapply(outList, function(x){x[[1]]}))
  maxLKR <- which(LKRs==max(LKRs))
  maxK <- kValues[maxLKR]
  thisOut <- outList[[maxLKR]]
  # get mean in each partition for transforming back to mean-scale
  Ytilde <- cut(Y, breaks=quantile(Y, p=seq(0,1,length=maxK+1)), labels=FALSE,
                include.lowest=TRUE)
  hK <- apply(matrix(1:maxK), 1, function(x){
    mean(Y[Ytilde==x])
  })
  # calculate mean by calculating density of each partition 
  pred <- apply(matrix(1:length(newX[,1])), 1, FUN=function(x){
    suppressWarnings(
      haz <- predict(thisOut$fm, 
                     newdata=data.frame(newX[x,],k=1:maxK),
                     type="response")
    )
    dens <- c(1,cumprod(1-haz[1:(maxK-1)])) * haz
    sum(hK*dens)
  })
  fit <- list(object=thisOut$fm, maxK=maxK, hK=hK)
  class(fit) <- "SL.gilleskie"
  out <- list(pred=pred, fit=fit)
  out
}

# predict function for SL.gilleskie
predict.SL.gilleskie <- function(object, newdata,...){
  pred <- apply(matrix(1:length(newdata[,1])), 1, FUN=function(x){
    suppressWarnings(
      haz <- predict(object$object,newdata=data.frame(newdata[rep(x,object$maxK),],k=1:object$maxK),
                     type="response")
    )
    dens <- c(1,cumprod(1-haz[1:(object$maxK-1)])) * haz
    sum(object$hK*dens)
  })
  pred
}
