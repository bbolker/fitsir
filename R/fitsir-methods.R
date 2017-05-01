##' @import methods
NULL

##' @aliases plot,fitsir-class
##' @describeIn fitsir plot deterministic trajectory
##' @importFrom bbmle plot
##' @examples
##' harbin2 <- setNames(harbin, c("times", "count"))
##' ff <- fitsir(harbin2, type="death", method="BFGS")
##' plot(ff)
##' 
##' ff2 <- fitsir(harbin2, type="death", dist="nbinom", method="BFGS")
##' ff3 <- fitsir(harbin2, type="death", dist="quasipoisson", method="BFGS")
##' plot(ff2, level=0.95, col.traj="red", main="Negative binomial error vs. Quasipoisson error CIs")
##' plot(ff3, add=TRUE, level=0.95, col.traj="blue", col.conf="blue")
##' legend(2, 270, legend = c("NB2", "Quasipoisson"), col=c("red", "blue"), lty=1)
##' 
setMethod("plot", signature(x="fitsir", y="missing"),
    function(x, level,
             method=c("delta", "mvrnorm", "wmvrnorm", "sobol"),
             main, xlim, ylim, xlab, ylab, add=FALSE,
             col.traj="black",lty.traj=1,
             col.conf="red",lty.conf=4,
             ...){
        method <- match.arg(method)
        count <- x@data$count
        pred <- predict(x,level,method=method)
        times <- pred[["times"]]
        i.hat <- pred[["mean"]]
        
        if (missing(main)) main <- "fitsir result"
        if (missing(xlab)) xlab <- "time"
        if (missing(ylab)) ylab <- "count"
        if (missing(ylim)) {
            ymin <- min(i.hat, count)
            ymax <- 1.1 * max(i.hat, count)
            ylim <- c(ymin, ymax)
        }
        if (missing(xlim)) xlim <- c(min(times), max(times))
        
        if (!add) plot(times, count, xlim=xlim, ylim=ylim,
                       xlab=xlab, ylab=ylab, main=main, ...)
        
        
        lines(times, i.hat, col=col.traj, lty=lty.traj)
        
        if (!missing(level)) {
            matlines(times, pred[,3:4], col=col.conf, lty=lty.conf)
        }
        
        invisible()
    }
)

##' @aliases coef,fitsir-class
##' @describeIn fitsir extract coefficients
##' @importFrom bbmle coef
##' @examples 
##' coef(ff)
##' coef(ff,"trans")
##' coef(ff,"summary")
##' 
setMethod("coef", "fitsir", 
    function(object,type=c("raw","trans","summary")){
        type <- match.arg(type)
        cc <- object@coef
        switch(type,
            raw=cc,
            trans=trans.pars(cc),
            summary=summarize.pars(cc)
        )
    }
)

##' @aliases predict,fitsir-class
##' @describeIn fitsir predict deterministic trajectory
##' @importFrom bbmle predict
##' @importFrom MASS mvrnorm
##' @importFrom mvtnorm dmvnorm
##' @examples
##' predict(ff, level=0.95)
##'
##' 
setMethod("predict", "fitsir",
    function(object,
             level,times,
             method=c("delta", "mvrnorm", "wmvrnorm", "sobol"),
             debug=FALSE){
        if(missing(times)) times <- object@data$times
        method <- match.arg(method)
        type <- object@data$type
        dist <- object@data$dist
        pars <- coef(object, "trans")
        i.hat <- SIR.detsim(times, pars, type)
        df <- data.frame(times=times,mean=i.hat)
        
        ## wquant from King et al.
        wquant <- function (x, weights, probs = c(0.025, 0.975)) {
            idx <- order(x)
            x <- x[idx]
            weights <- weights[idx]
            w <- cumsum(weights)/sum(weights)
            rval <- approx(w,x,probs,rule=1)
            rval$y
        }
        
        if (!missing(level)) {
            ll <- (1-level)/2
            
            nsim <- 2000
            
            if (method != "delta") {
                simtraj <- matrix(NA,nrow=length(times),ncol=nsim)
                if (method == "sobol") {
                    ranges <- confint(object)
                    simpars <- sobolDesign(lower=ranges[,1],
                                           upper=ranges[,2],
                                           nseq=nsim)
                } else {
                    simpars <- mvrnorm(nsim,mu=coef(object),
                                       Sigma=vcov(object))
                }
                
                for (i in 1:nsim) {
                    simtraj[,i] <- SIR.detsim(times, 
                                              trans.pars(simpars[i,]),
                                              type)
                }
            }
            
            cmat <- switch(method,
                delta={
                    nmat <- as.matrix(SIR.detsim(times, pars, type, grad=TRUE)[,-1])
                    xvcov <- object@vcov[1:4,1:4]
                    if(any(diag(xvcov < 0)))
                        warning("At least one entries in diag(vcov) is negative. Confidence interval may not be accurate.")
                    
                    ivcov <- nmat %*% xvcov %*% t(nmat)
                    ierr <- sqrt(diag(ivcov))
                    z <- -qnorm(ll)
                    cmat <- data.frame(i.hat - z * ierr, i.hat + z * ierr)
                    cmat
                },
                mvrnorm={
                    cmat <- t(apply(simtraj,1,quantile,c(ll,1-ll)))
                    cmat
                },
                wmvrnorm={
                    traj.logLik <- -apply(simpars, 1, SIR.logLik, count=object@data$count, times, dist=dist, type=type)
                    ##FIXME: vcov not symmetric for low tolerance?
                    i <- 10
                    while(!isSymmetric(round(vcov(object), i))) i <- i - 1
                    sample.logLik <- dmvnorm(simpars, coef(object), round(vcov(object), i), log=TRUE)
                    ww <- exp(traj.logLik-sample.logLik)
                    cmat <- t(apply(simtraj, 1, wquant, weights=ww, probs=c(ll, 1-ll)))
                    cmat
                },
                sobol={
                    traj.logLik <- -apply(simpars, 1, SIR.logLik, count=object@data$count, times, dist=dist, type=type)
                    cc <- which(traj.logLik > max(traj.logLik) - 100)
                    traj.logLik <- traj.logLik[cc]
                    ww <- exp(traj.logLik-mean(traj.logLik))
                    cmat <- t(apply(simtraj[,cc], 1, wquant, weights=ww, probs=c(ll, 1-ll)))
                    cmat
                })
            
            cmat <- setNames(cmat, c(paste(100*ll, "%"), paste(100*(1-ll), "%")))
            
            if (debug & method != "delta") {
                matplot(times, simtraj, type="l",col=adjustcolor("black", alpha=0.1), lty=1)
                matlines(times, cmat, col=2, lty=1, lwd=2)
            }
            
            df <- cbind(df, cmat)
        }
        df
    }
)


##' @aliases residuals,fitsir-class
##' @describeIn fitsir calculate residuals between data and deterministic trajectory
##' @importFrom bbmle residuals
##' @examples
##' residuals(ff)
##' 
setMethod("residuals", "fitsir",
    function(object,type=c("pearson", "raw")){
        type <- match.arg(type)
        pred <- predict(object)
        mean <- pred$mean
        count <- object@data$count
        dd <- mean-count
        var <- dd^2
        
        switch(type,
            pearson=var/mean,
            raw=dd
        )
    }
)

##' @exportMethod sigma
setGeneric("sigma", function(object, ...) standardGeneric("sigma"))
setMethod("sigma", "fitsir",
    function(object,dist=c("quasipoisson", "nbinom", "nbinom1")){
        dist <- match.arg(dist)
        pred <- predict(object)
        mean <- pred$mean
        count <- object@data$count
        n <- length(count)
        var <- (mean - count)^2
        if(dist != "quasipoisson")
            log.dsp <- unname(mledsp(count,mean,dist))
          
        switch(dist,
            quasipoisson=sum(var/mean)/(n-1),
            nbinom=exp(log.dsp),
            nbinom1=exp(log.dsp)
        )
    }
)

##' @aliases summary,fitsir-class
##' @describeIn fitsir summarize fit using meaningful parameters
##' @importFrom bbmle summary
setMethod("summary", "fitsir",
    function(object,...){
        cc <- object@coef
        ss <- callNextMethod()
        m <- summarize.pars.jacobian(cc)
        cc.vcov <- m %*% object@vcov[1:4,1:4] %*% t(m)
        smat <- rbind(
            Estimate=summarize.pars(cc),
            `Std. Error` = sqrt(diag(cc.vcov))
        )
        m2logL <- 2*object@min
        new("summary.fitsir", call=object@call.orig, coef=ss@coef, summary=smat, m2logL=m2logL)
    }
)

##' @aliases show,summary.fitsir-class
##' @describeIn summary.fitsir pretty-prints \code{object}
setMethod("show", "summary.fitsir",
    function(object){
        cat("Maximum likelihood estimation\n\nCall:\n")
        print(object@call)
        cat("\nCoefficients:\n")
        printCoefmat(object@summary)
        cat("\n-2 log L:", object@m2logL, "\n")
    }
)
