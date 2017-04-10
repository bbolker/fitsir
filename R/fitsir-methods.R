##' @import methods
NULL

##' S4 plot method for fitsir objects
##' 
##' @name plot.fitsir
##' @rdname plot.fitsir
##' @param x a fitsir object
##' @param level the confidence level required
##' @param main main title
##' @param xlab x label
##' @param ylab y label
##' @param add (logical) add trajectories to an existing figure
##' @param col.traj color for main trajectory
##' @param lty.traj line type for main trajectory
##' @param col.conf color for trajectories based on confidence intervals
##' @param lty.conf line type for trajectories based on confidence intervals
##' @param ... additional arguments
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
setMethod("plot", signature(x="fitsir", y="missing"),
    function(x, level,
             method=c("delta", "sample"),
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

##' S4 coef method for fitsir objects
##' @name coef.fitsir
##' @rdname coef.fitsir
##' @param object a fitsir object
##' @param type types of returned parameters
##' @importFrom bbmle coef
##' @examples 
##' ff <- fitsir(harbin, type="death", method="BFGS", tcol="week", icol="Deaths")
##' coef(ff)
##' coef(ff,"trans")
##' coef(ff,"summary")
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

##' S4 predict method for fitsir objects
##' @name predict.fitsir
##' @rdname predict.fitsir
##' @param object a fitsir object
##' @param level the confidence level required
##' @param times new time vector
##' @importFrom bbmle predict
##' @importFrom MASS mvrnorm
##' @examples
##' ff <- fitsir(harbin, type="death", method="BFGS", tcol="week", icol="Deaths")
##' predict(ff, level=0.95)
setMethod("predict", "fitsir",
    function(object,
             level,times,
             method=c("delta", "sample"),
             debug=FALSE){
        if(missing(times)) times <- object@data$times
        method <- match.arg(method)
        type <- object@data$type
        pars <- coef(object, "trans")
        i.hat <- SIR.detsim(times, pars, type)
        df <- data.frame(times=times,mean=i.hat)
        
        if (!missing(level)) {
            ll <- (1-level)/2
            
            if (method=="delta") {
                nmat <- as.matrix(SIR.detsim(times, pars, type, grad=TRUE)[,-1])
                xvcov <- object@vcov[1:4,1:4]
                if(any(diag(xvcov < 0)))
                    warning("At least one entries in diag(vcov) is negative. Confidence interval may not be accurate.")
                
                ivcov <- nmat %*% xvcov %*% t(nmat)
                ierr <- sqrt(diag(ivcov))
                z <- -qnorm(ll)
                cmat <- data.frame(i.hat - z * ierr, i.hat + z * ierr)
                cmat <- setNames(cmat, c(paste(100*ll, "%"), paste(100*(1-ll), "%")))
            } else {
                nsim <- 1000
                simtraj <- matrix(NA,nrow=length(times),ncol=nsim)
                simpars <- mvrnorm(nsim,mu=coef(object),
                                   Sigma=vcov(object))
                for (i in 1:nsim) {
                    simtraj[,i] <- SIR.detsim(times, 
                                              trans.pars(simpars[i,]),
                                              type)
                }
                cmat <- t(apply(simtraj,1,quantile,
                                c(ll,1-ll)))
                if (debug) {
                    matplot(times, simtraj, type="l",col=adjustcolor("black", alpha=0.1), lty=1)
                    matlines(times, cmat, col=2, lty=1, lwd=2)
                }
            }
            df <- cbind(df, cmat)
        }
        df
    }
)

##' S4 residuals method for fitsir objects
##' @name residuals.fitsir
##' @rdname residuals.fitsir
##' @param object a fitsir object
##' @param type type of residuals
##' @importFrom bbmle residuals
##' @examples
##' ff <- fitsir(harbin, type="death", method="BFGS", tcol="week", icol="Deaths")
##' residuals(ff)
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
            nbinom=1/exp(log.dsp),
            nbinom1=exp(log.dsp)
        )
    }
)

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

setMethod("show", "summary.fitsir",
    function(object){
        cat("Maximum likelihood estimation\n\nCall:\n")
        print(object@call)
        cat("\nCoefficients:\n")
        printCoefmat(object@summary)
        cat("\n-2 log L:", object@m2logL, "\n")
    }
)
