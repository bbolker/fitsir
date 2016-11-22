library(fitsir)
library(bbmle) ## needed for coef() ...
library(splines) ## for ns()

##' Stochastic simulation (observation/uncorrelated error only)
##' 
##' @param pars parameters (on "original" constrained scale)
##' @param tmax max times
##' @param dt time step
##' @param rfun function for generating random (observation) noise
##' @param rmean name of mean for \code{rfun}
##' @param rpars additional arguments for \code{rfun}
simfun <- function(pars=c(beta=0.2,gamma=0.1,N=1000,i0=0.01),
                   tmax=100,dt=1,
                   rfun=rnorm,
                   rmean="mean",
                   rpars=list(sd = 3)
                   ) {
    tvec <- seq(0,tmax,by=dt)
    ss <- SIR.detsim(tvec,pars)
    noiseArgs <- c(setNames(list(length(ss),ss),c("n",rmean)),
                   rpars)
    count <- do.call(rfun,noiseArgs)
    return(data.frame(tvec,count))
}

set.seed(101)
s0 <- simfun()

## get starting values and trajectory based on them
ss0 <- startfun(s0)
ss2 <- SIR.detsim(s0$tvec,trans.pars(ss0))

## fit and corresponding trajectory
t1 <- system.time(f1 <- fitsir(s0,start=ss0))
ss3 <- SIR.detsim(s0$tvec,trans.pars(coef(f1)))

## GAM fit: match number of degrees of freedom
## glm() lists df as df+2 
##          == length(coef())+1
##   counting for intercept and residual var(?)
m1 <- glm(count~ns(tvec,df=3),
          family=gaussian(link="log"),data=s0)

plot(count~tvec,data=s0)
lines(s0$tvec,ss2)  ## incidence/prevalence mismatch?
lines(s0$tvec,ss3,col=2)  ## decent fit anyway
lines(s0$tvec,predict(m1,type="response"),col=4)
-logLik(m1)-(-logLik(f1))

## takes about 1 minute per sim ...
fitfun <- function(data) {
    t1 <- system.time(f1 <- fitsir(data,start=startfun(auto=TRUE,data=data)))
    m1 <- glm(count~ns(tvec,df=3),family=gaussian(link="log"),data=data)
    res <- c(t=unname(t1["elapsed"]),coef(f1),
      nll.SIR=c(-logLik(f1)),nll.gam=c(-logLik(m1)))
    return(res)
}

load("stochsim.rda")
fitfun(s0)

library(plyr)
set.seed(101)
res1 <- raply(20,fitfun(simfun(rpars=list(size=10))),
              .progress="text")
res1 <- as.data.frame(res1)
save("res1",file="stochsim.rda")
with(res1,hist(nll.SIR-nll.gam,col="gray"))
