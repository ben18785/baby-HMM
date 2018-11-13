# initial stuff ---------------
rm(list=ls())
setwd("~/Desktop/baby-HMM/Analysis")
library(rstan)
library(ggplot2)
library(reshape2)
options(mc.cores=parallel::detectCores())

# get data for participant 3 for block 1 ---------------
aParticipantNumber <- 5
lID <- unique(aDF$ID)
aDF <- read.csv('babydata.csv')
aDF <- aDF[order(aDF$ID,aDF$BLOCK,aDF$TRIAL,aDF$BIN),]
aDF1 <- aDF[aDF$ID==lID[aParticipantNumber],]
aDF1 <- aDF1[aDF1$BLOCK==1,]

# estimate model -------------
## compile model 
bModel <- stan_model('hmm.stan')
## function to generate initial parameter values (not always needed but here for some reason
## it is!)
initFn1 <- function(){
  list(mu1 = c(40,60,110),
       sigma1 = rnorm(K,40,0.1),
       mu2 = c(40,60,110),
       sigma2 = rnorm(K,40,0.1))}
## number of states
K <- 3
## estimate model using Stan's HMC
aFit <- sampling(bModel, data=list(N=nrow(aDF1),simple=aDF1$simple,complex=aDF1$complex,K=3),
                 iter=200,chains=4,init=initFn1)
## check MCMC diagnostics
print(aFit)

# graph the results -----------------------
## extrate mean state
lState <- apply(extract(aFit,'state')[[1]], 2, median)
aGraph_DF <- data.frame(time=aDF1$BIN, simple=aDF1$simple/120,
                        complex=aDF1$complex/120,neither=aDF1$neither/120,
                        state=lState)
aGraph_DF <- melt(aGraph_DF,id.vars='time')
ggplot(aGraph_DF, aes(x=time, y=value, colour=as.factor(variable))) +
  geom_path()


# fit all of the participants' data using the same model --------------
## estimate model for all participants
lFits <- vector(length = length(lID),mode = 'list')
for(i in seq_along(lID)){
  print(i)
  aDF1 <- aDF[aDF$ID==lID[i],]
  aDF1 <- aDF1[aDF1$BLOCK==1,]
  lFits[[i]] <- sampling(bModel,data=list(N=nrow(aDF1),simple=aDF1$simple,complex=aDF1$complex,K=3),iter=200,chains=4,init=initFn1)
}

## extract states for all the individuals and combine in a data frame that can be used to graph
lStates <- vector(mode='list',length = length(lID))
for(i in seq_along(lID)){
  print(i)
  aDF1 <- aDF[aDF$ID==lID[i],]
  aDF1 <- aDF1[aDF1$BLOCK==1,]
  aFit <- lFits[[i]]
  aState <- apply(extract(aFit,'state')[[1]], 2, median)  
  lStates[[i]] <- data.frame(state=aState,simple=aDF1$simple/120,complex=aDF1$complex/120,neither=aDF1$neither/120,time=aDF1$BIN)
}

# graph a specific participant's fit
aNumber <- 1
aDF_short <- lStates[[aNumber]]
aDF_short <- melt(aDF_short,id.vars = 'time')
ggplot(aDF_short,aes(time,value,colour=as.factor(variable))) + geom_path()
print(lFits[[aNumber]])
