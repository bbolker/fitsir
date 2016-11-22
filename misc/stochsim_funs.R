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
                   rfun=rnorm,
                   rmean="mean",
                   rpars=list(sd=10),
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
    lastpos <- head(which(count <= 0),1)
    return(data.frame(tvec,count)[1:(lastpos-1),])
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
                    start_method=c("auto","true"),
                    ## other possibilities: "single", "lhs"
                    truepars=NULL,
                    plot.it=FALSE,
                    spline.method="ss",
                    spline.var="log",
                    spline.df=6,
                    ...) {
    require(bbmle)  ## for coef, logLik ... should patch up ...
    require(splines) ## ns()
    fitres <- list()
    nzdat <- subset(data,count>0)
    for (s in start_method) {
        ss0 <- switch(s,
                      auto=startfun(data=data,auto=TRUE),
                      true=truepars,
                      single=startfun(), ## default starting params
                      lhs=NA)
        if (s!="lhs") {
            t1 <- system.time(f1 <- fitsir(data,start=ss0))
        } else {
            stop("lhs not implemented yet")
        }
        fcoef <- coef(f1)
        fpred <- SIR.detsim(nzdat$tvec,trans.pars(fcoef))
        mse <- mean((1-fpred/nzdat$count)^2)
        fitres[[s]] <- list(time=unname(t1["elapsed"]),fit=f1,coef=fcoef,
                            pred=fpred,mse=mse,nll=c(-logLik(f1)))
    }
    nspline <- max(length(spline.method),length(spline.df),
                   length(spline.var))
    spline.method <- rep(spline.method,length.out=nspline)
    spline.var <- rep(spline.var,length.out=nspline)
    spline.df <- rep(spline.df,length.out=nspline)
    splineres <- list()
    for (i in 1:nspline) {
        meth <- spline.method[i]
        cm1 <- setNames(rep(NA,spline.df[i]),
                        c("intercept",paste0("beta",1:(spline.df[i]-1))))
        if (meth=="ss") {
            t1 <- system.time(m1 <- smooth.spline(nzdat$tvec,log(nzdat$count),
                                nknots=spline.df[i]-2))
            spred <- exp(fitted(m1))
            nll <- NA

        } else {  ## method == "ns"
            if (spline.var[i]=="log") {
                t1 <- system.time(m1 <- lm(log(count)~ns(tvec,df=spline.df[i]-1),data=nzdat))
                nll <- c(-logLik(m1))+sum(1/nzdat$count)
            } else {
                t1 <- system.time(m1 <- glm(count~ns(tvec,df=spline.df[i]-1),
                                            data=nzdat,
                                            family=gaussian(link="log")))
                nll <- c(-logLik(m1))
            }
            cm1[] <- coef(m1) ## leave names alone, reset values
        }
        spred <- exp(fitted(m1))
        splineres[[i]] <- list(fit=m1,
                               ## hack coeff names to be nicer ...
                               coef=cm1,
                               pred=spred,
                               mse=mean((1-spred/nzdat$count)^2))
    }
    ## if (plot.it) {
    ##     plot(data$tvec,data$count,xlab="time",ylab="count",type="l",...)
    ##     matpoints(nzdat$tvec,cbind(fpred,spred1),col=c(2,4),pch=1:2)
    ##     legend("topright",
    ##            c("data","fitsir","spline"),
    ##            col=c(1,2,4),lty=1)
    ## }
    ## collect
    res <- numeric(0)
    sumfun <- function(x,nm) {
        r <- unlist(x[c("time","coef","mse","nll")])
        names(r) <- paste(nm,names(r),sep="_")
        return(r)
    }
    for (i in seq_along(fitres))
        res <- c(res,sumfun(fitres[[i]],paste("fitsir",i,sep=".")))
    for (i in seq_along(splineres))
        res <- c(res,sumfun(splineres[[i]],paste("spline",i,sep=".")))
    return(res)
}

