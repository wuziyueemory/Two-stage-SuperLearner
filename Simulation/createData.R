# Functions for data generating process
# Options: 
# 1) sample size
# 2) zero-inflation percentage 
# 3) non-zero distribution 
# 4) model complexity - with/without interaction 

# Predictors: 10 variables (5 used, 5 non-used)

library(mgcv) # for generating tweedie distribution

# define functions for creating mixture distribution of log-normal & gamma
rlnormgammamix <- function(n,prob,shape,mu) {
  ifelse(runif(n)<prob,
         rgamma(n,shape = shape, scale = 1.5),
         exp(mu + rnorm(n,0,0.2)))
}

# Functions for generating simulation data
createData <- function(n,zero,nonzero,linear,interaction){
  # generate covariates (10)
  # w1 & w6 have bernoulli distribution
  # w2 & w7 have Uniform distribution
  # w3 & w8 have a normal distribution
  # w4 & w9 have a gamma distribution
  # w5 & w10 have a Poisson distribution
  w1 <- rbinom(n,1,0.5)
  w2 <- runif(n,-2,2)
  w3 <- rnorm(n,0,1)
  w4 <- rgamma(n,1,0.5)
  w5 <- rpois(n,1)
  w6 <- rbinom(n,1,0.2)
  w7 <- runif(n,0,1)
  w8 <- rnorm(n,0,3)
  w9 <- rgamma(n,0.5,1)
  w10 <- rpois(n,2)
  y <- rep(NA,n)
  
  # whether interaction
  if (interaction){
    # with interaction
    # whether linear
    if (linear){
      # linear 
      # linear regression for positive (non-zero) part
      main <- 1 + 0.3*w1 - w2 + 0.2*w3 - 0.15*w4 + 0.1*w5 - 0.01*w1*w2 + 
        0.03*w1*w3 - 0.1*w2*w3 + 0.02*w3*w4 - 0.04*w4*w5
      
      # zero proportion
      if (zero==0.05) {
        # probability of y=0
        prob <- plogis(4.25 + w1 + w2 + w3 + w4 + w5 + w1*w2 + w1*w3 + w1*w4 + w1*w5 + 
                         w2*w3 + w2*w4 + w2*w5 + w3*w4 + w3*w5 + w4*w5)
        g <- rbinom(n,1,prob)
        # assign g=0 costs
        ind <- g==0
        y[ind] <- 0
        
        # nonzero dist.
        # assign g=1 costs
        ind <- g==1
        if (nonzero=="lognormal") { #log-normal
          y[ind] <- exp(9 + main[ind] + rnorm(sum(ind),0,0.3))
        } else if (nonzero=="gamma") { #gamma
          y[ind] <- rgamma(sum(ind),shape = exp(9 + main[ind]), scale = 1.3)
        } else if (nonzero=="tweedie") { #Tweedie
          y[ind] <- rTweedie(exp(9 + main[ind]),p = 1.5,phi = 1.8)
        } else if (nonzero=="mixture") { #mixture
          y[ind] <- rlnormgammamix(sum(ind),0.5,shape=exp(9 + main[ind]),mu=9+main[ind])
        }
      } else if (zero==0.7){
        # probability of y=0
        prob <- plogis(-1.65 + 0.1*(w1 + w2 + w3 + w4 + w5 + w1*w2 + w1*w3 + w1*w4 + w1*w5 + 
                                      w2*w3 + w2*w4 + w2*w5 + w3*w4 + w3*w5 + w4*w5))
        g <- rbinom(n,1,prob)
        # assign g=0 costs
        ind <- g==0
        y[ind] <- 0
        
        # nonzero dist.
        # assign g=1 costs
        ind <- g==1
        if (nonzero=="lognormal") { #log-normal
          y[ind] <- exp(9 + main[ind] + rnorm(sum(ind),0,0.3))
        } else if (nonzero=="gamma") { #gamma
          y[ind] <- rgamma(sum(ind),shape = exp(9 + main[ind]), scale = 1.3)
        } else if (nonzero=="tweedie") { #Tweedie
          y[ind] <- rTweedie(exp(9 + main[ind]),p = 1.5,phi = 1.8)
        } else if (nonzero=="mixture") { #mixture
          y[ind] <- rlnormgammamix(sum(ind),0.5,shape=exp(9 + main[ind]),mu=9+main[ind])
        }
      }
    } else {
      # nonlinear
      # non-linear regression for positive (non-zero) part
      main <- 1 - 0.5*w1 - 1.2*w2 - 0.1*w3 + 0.2*w4 - 0.3*w5 - 0.15*w1^2 + 0.2*w2^2 - 
        0.05*w3^2 - 0.03*w4^2 - 0.01*w5^2 - 0.04*w1*w3 + 0.02*w1*w5 - 0.06*w2*w3 + 
        0.07*w2*w4 - 0.05*w3*w4
      
      # zero proportion
      if (zero==0.05) {
        # probability of y=0
        prob <- plogis(27.432+ w1 + 2*w2 + 3*w3 + 0.5*w4 + 1.5*w5 + w1^2 - w2^2 + w3^2 - 
                         w4^2 + w5^2 + w1*w2 - w1*w3 + w2*w4 - w3*w5 + w4*w5)
        g <- rbinom(n,1,prob)
        # assign g=0 costs
        ind <- g==0
        y[ind] <- 0
        
        # nonzero dist.
        # assign g=1 costs
        ind <- g==1
        if (nonzero=="lognormal") { #log-normal
          y[ind] <- exp(9 + main[ind] + rnorm(sum(ind),0,0.3))
        } else if (nonzero=="gamma") { #gamma
          y[ind] <- rgamma(sum(ind),shape = exp(9 + main[ind]), scale = 1.3)
        } else if (nonzero=="tweedie") { #Tweedie
          y[ind] <- rTweedie(exp(9 + main[ind]),p = 1.5,phi = 1.8)
        } else if (nonzero=="mixture") { #mixture
          y[ind] <- rlnormgammamix(sum(ind),0.5,shape=exp(9 + main[ind]),mu=9+main[ind])
        }
      } else if (zero==0.7){
        # probability of y=0
        prob <- plogis(-6.3 + w1 + 2*w2 + 3*w3 + 0.5*w4 + 1.5*w5 + w1^2 - w2^2 + w3^2 - 
                         w4^2 + w5^2 + w1*w2 - w1*w3 + w2*w4 - w3*w5 + w4*w5)
        g <- rbinom(n,1,prob)
        # assign g=0 costs
        ind <- g==0
        y[ind] <- 0
        
        # nonzero dist.
        # assign g=1 costs
        ind <- g==1
        if (nonzero=="lognormal") { #log-normal
          y[ind] <- exp(9 + main[ind] + rnorm(sum(ind),0,0.3))
        } else if (nonzero=="gamma") { #gamma
          y[ind] <- rgamma(sum(ind),shape = exp(9 + main[ind]), scale = 1.3)
        } else if (nonzero=="tweedie") { #Tweedie
          y[ind] <- rTweedie(exp(9 + main[ind]),p = 1.5,phi = 1.8)
        } else if (nonzero=="mixture") { #mixture
          y[ind] <- rlnormgammamix(sum(ind),0.5,shape=exp(9 + main[ind]),mu=9+main[ind])
        }
      }
    }
  } else {
    # no interaction
    # whether linear
    if (linear){
      # linear regression for positive (non-zero) part
      main <- 0.1 + 0.5*w1 - 1.5*w2 + 0.4*w3 - 0.3*w4 - 0.2*w5
      
      # zero proportion
      if (zero==0.05) {
        # probability of y=0
        prob <- plogis(1.459 + w1 + w2 + w3 + w4 + w5)
        g <- rbinom(n,1,prob)
        # assign g=0 costs
        ind <- g==0
        y[ind] <- 0
        
        # nonzero dist.
        # assign g=1 costs
        ind <- g==1
        if (nonzero=="lognormal") { #log-normal
          y[ind] <- exp(9 + main[ind] + rnorm(sum(ind),0,0.3))
        } else if (nonzero=="gamma") { #gamma
          y[ind] <- rgamma(sum(ind),shape = exp(9 + main[ind]), scale = 1.3)
        } else if (nonzero=="tweedie") { #Tweedie
          y[ind] <- rTweedie(exp(9 + main[ind]),p = 1.5,phi = 1.8)
        } else if (nonzero=="mixture") { #mixture
          y[ind] <- rlnormgammamix(sum(ind),0.5,shape=exp(9 + main[ind]),mu=9+main[ind])
        }
      } else if (zero==0.7){
        # probability of y=0
        prob <- plogis(-1.21 + 0.1*(w1 + w2 + w3 + w4 + w5))
        g <- rbinom(n,1,prob)
        # assign g=0 costs
        ind <- g==0
        y[ind] <- 0
        
        # nonzero dist.
        # assign g=1 costs
        ind <- g==1
        if (nonzero=="lognormal") { #log-normal
          y[ind] <- exp(9 + main[ind] + rnorm(sum(ind),0,0.3))
        } else if (nonzero=="gamma") { #gamma
          y[ind] <- rgamma(sum(ind),shape = exp(9 + main[ind]), scale = 1.3)
        } else if (nonzero=="tweedie") { #Tweedie
          y[ind] <- rTweedie(exp(9 + main[ind]),p = 1.5,phi = 1.8)
        } else if (nonzero=="mixture") { #mixture
          y[ind] <- rlnormgammamix(sum(ind),0.5,shape=exp(9 + main[ind]),mu=9+main[ind])
        }
      }
    } else {
      # non-linear regression for positive (non-zero) part
      main <- 0.1 - 0.5*w1 + 1.3*w2 - 0.2*w3 + 0.15*w4 - 0.3*w5 + 0.3*w1^2 - 0.4*w2^2 - 
        0.2*w3^2 + 0.07*w4^2 - 0.05*w5^2
      
      # zero proportion
      if (zero==0.05) {
        # probability of y=0
        prob <- plogis(42.88 + w1 + 2*w2 + 0.8*w3 + 1.2*w4 + 0.5*w5 + 0.5*w1^2 - 
                         0.5*w2^2 + 1.5*w3^2 - 1.5*w4^2 + w5^2)
        g <- rbinom(n,1,prob)
        # assign g=0 costs
        ind <- g==0
        y[ind] <- 0
        
        # nonzero dist.
        # assign g=1 costs
        ind <- g==1
        if (nonzero=="lognormal") { #log-normal
          y[ind] <- exp(9 + main[ind] + rnorm(sum(ind),0,0.3))
        } else if (nonzero=="gamma") { #gamma
          y[ind] <- rgamma(sum(ind),shape = exp(9 + main[ind]), scale = 1.3)
        } else if (nonzero=="tweedie") { #Tweedie
          y[ind] <- rTweedie(exp(9 + main[ind]),p = 1.5,phi = 1.8)
        } else if (nonzero=="mixture") { #mixture
          y[ind] <- rlnormgammamix(sum(ind),0.5,shape=exp(9 + main[ind]),mu=9+main[ind])
        }
      } else if (zero==0.7){
        # probability of y=0
        prob <- plogis(-3.75 + w1 + 2*w2 + 0.8*w3 + 1.2*w4 + 0.5*w5 + 0.5*w1^2 - 
                         0.5*w2^2 + 1.5*w3^2 - 1.5*w4^2 + w5^2)
        g <- rbinom(n,1,prob)
        # assign g=0 costs
        ind <- g==0
        y[ind] <- 0
        
        # nonzero dist.
        # assign g=1 costs
        ind <- g==1
        if (nonzero=="lognormal") { #log-normal
          y[ind] <- exp(9 + main[ind] + rnorm(sum(ind),0,0.3))
        } else if (nonzero=="gamma") { #gamma
          y[ind] <- rgamma(sum(ind),shape = exp(9 + main[ind]), scale = 1.3)
        } else if (nonzero=="tweedie") { #Tweedie
          y[ind] <- rTweedie(exp(9 + main[ind]),p = 1.5,phi = 1.8)
        } else if (nonzero=="mixture") { #mixture
          y[ind] <- rlnormgammamix(sum(ind),0.5,shape=exp(9 + main[ind]),mu=9+main[ind])
        }
      }
    }
  }
  # generate a dataframe
  return(data.frame(pid=1:n,w1=w1,w2=w2,w3=w3,w4=w4,w5=w5,w6=w6,
                    w7=w7,w8=w8,w9=w9,w10=w10,y=y))
}
