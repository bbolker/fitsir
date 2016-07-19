##' Stochastic simulation (observation/uncorrelated error only)
##' 
##' should this be called SIR.stochsim instead?
##' allow facility for (discrete and/or continuous time) process-error
##'   sims?
##'
##' @param pars parameters (on "original" constrained scale)
##' @param tmax max times
##' @param dt time step
##' @param rfun function for generating random (observation) noise
##' @param rmean name of mean for \code{rfun}
##' @param rpars additional arguments for \code{rfun}
simfun <- function(pars=c(beta=0.2,gamma=0.1,N=1000,i0=0.01),
                   tmax=100,dt=1,
                   rfun=rnbinom,
                   rmean="mu",
                   rpars=list(size=1),
                   seed=NULL,
                   drop.zeros=TRUE
                   ) {
    if (!is.null(seed)) set.seed(seed)
    tvec <- seq(0,tmax,by=dt)
    ss <- SIR.detsim(tvec,pars)
    noiseArgs <- c(setNames(list(length(ss),ss),c("n",rmean)),
                   rpars)
    count <- do.call(rfun,noiseArgs)
    ##
    lastpos <- tail(which(count>0),1)
    return(data.frame(tvec,count)[1:lastpos,])
}

fitfun <- function(data) {
    t1 <- system.time(f1 <- fitsir(data))
                                   ## start=startfun(auto=TRUE,data=data)))
    m1 <- glm(count~ns(tvec,df=3),family=gaussian(link="log"),data=data)
    res <- c(t=unname(t1["elapsed"]),coef(f1),coef(m1),
      nll.SIR=c(-logLik(f1)),nll.gam=c(-logLik(m1)))
    return(res)
}

fitfun2 <- function(data,
                    start_method=c("auto","lhs","true","single"),
                    truepars=NULL,
                    plot.it=FALSE,...) {
    require(bbmle)  ## for coef, logLik ... should patch up ...
    require(splines) ## ns()
    start_method <- match.arg(start_method)
    ss0 <- switch(start_method,
                  auto=startfun(data=data,auto=TRUE),
                  true=truepars,
                  single=startfun(), ## default starting params
                  lhs=NA)
    if (start_method!="lhs") {
        t1 <- system.time(f1 <- fitsir(data,start=ss0))
    } else {
        stop("lhs not implemented yet")
    }
    fcoef <- bbmle::coef(f1)
    ## m1 <- glm(count~ns(tvec,df=3),family=gaussian(link="log"),data=data)
    ## m1 <- lm(log(count+1)~ns(tvec,df=3),data=data)
    nzdat <- subset(data,count>0)
    m1 <- smooth.spline(nzdat$tvec,log(nzdat$count),nknots=4)
    m2 <- lm(log(count)~ns(tvec,df=6),data=nzdat)
    spred1 <- exp(fitted(m1))
    spred2 <- exp(fitted(m2))
    fpred <- SIR.detsim(nzdat$tvec,trans.pars(fcoef))
    mse <- c(mse.fitsir=mean((1-fpred/nzdat$count)^2),
             mse.spline1=mean((1-spred1/nzdat$count)^2),
             mse.spline2=mean((1-spred2/nzdat$count)^2))
    if (plot.it) {
        plot(data$tvec,data$count,xlab="time",ylab="count",type="l",...)
        matpoints(nzdat$tvec,cbind(fpred,spred1),col=c(2,4),pch=1:2)
        legend("topright",
               c("data","fitsir","spline"),
               col=c(1,2,4),lty=1)
    }
    res <- c(time=unname(t1["elapsed"]),fcoef,
             nll.SIR=c(-logLik(f1)),
             nll.spline2=c(-logLik(m2)), ## FIXME: adjust for transformation
             ## (+ sum(1/y) ?)
             mse)
    return(res)
}

nmvec2 <- c("time", "log.beta", "log.gamma", "log.N", "logit.i", "nll.SIR", 
            "nll.spline2", "mse.fitsir", "mse.spline1", "mse.spline2")
