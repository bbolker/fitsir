##' S4 plot method for fitsir objects
##' 
##' @name plot.fitsir
##' @rdname plot.fitsir
##' @param x An object of class \code{fitsir}
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
##' @examples
##' harbin2 <- setNames(harbin, c("times", "count"))
##' ff <- fitsir(harbin2, type="death")
##' plot(ff)
##' 
##' ff2 <- fitsir(harbin2, type="death", dist="nbinom")
##' ff3 <- fitsir(harbin2, type="death", dist="quasipoisson")
##' plot(ff2, level=0.95, col.traj="red", main="Negative binomial error vs. Quasipoisson error CIs")
##' plot(ff3, add=TRUE, level=0.95, col.traj="blue", col.conf="blue")
##' legend(2, 270, legend = c("NB2", "Quasipoisson"), col=c("red", "blue"), lty=1)
setMethod("plot", signature(x="fitsir", y="missing"),
    function(x, level,
             main, xlim, ylim, xlab, ylab, add=FALSE,
             col.traj="black",lty.traj=1,
             col.conf="red",lty.conf=4,
             ...){
        count <- x@data$count
        pred <- predict(x,level)
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
        
        if (!add) {
            plot(times, i.hat, type="l",
                 col=col.traj, lty=lty.traj, xlim=xlim, ylim=ylim,
                 xlab=xlab, ylab=ylab, main=main, ...)
            points(times, count)
        } else {
            lines(times, i.hat, col=col.traj, lty=lty.traj, ...)
        }
        
        if (!missing(level)) {
            matlines(times, pred[,3:4], col=col.conf, lty=lty.conf)
        }
        
        invisible()
    }
)

##' Coef method for fitsir objects
##' @name coef.fitsir
##' @rdname coef.fitsir
##' @param object An object of class \code{fitsir}
##' @param type types of returned parameters
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

setMethod("predict", "fitsir",
    function(object,level,times){
        if(missing(times)) times <- object@data$times
        type <- object@data$type
        pars <- coef(object, "trans")
        i.hat <- SIR.detsim(times, pars, type)
        df <- data.frame(times=times,mean=i.hat)
        
        if(!missing(level)) {
            nmat <- as.matrix(SIR.detsim(times, pars, type, grad=TRUE)[,-1])
            xvcov <- object@vcov
            if(any(diag(xvcov < 0)))
                warning("At least one entries in diag(vcov) is negative. Confidence interval may not be accurate.")
            
            ivcov <- nmat %*% xvcov %*% t(nmat)
            ierr <- sqrt(diag(ivcov))
            ll <- (1-level)/2
            z <- -qnorm(ll)
            cmat <- data.frame(i.hat + z * ierr, i.hat - z * ierr)
            cmat <- setNames(cmat, c(paste(100*ll, "%"), paste(100*(1-ll), "%")))
            df <- cbind(df, cmat)
        }
        df
    }
)

setMethod("summary", "fitsir",
    function(object,...){
        cc <- object@coef
        ss <- callNextMethod()
        m <- summarize.pars.jacobian(cc)
        cc.vcov <- m %*% object@vcov %*% t(m)
        smat <- rbind(
            Estimate=summarize.pars(cc),
            `Std. Error` = sqrt(diag(cc.vcov))
        )
        m2logL <- 2*object@min
        new("summary-fitsir", call=object@call.orig, coef=ss@coef, summary=smat, m2logL=m2logL)
    }
)

setMethod("show", "summary-fitsir",
    function(object){
        cat("Maximum likelihood estimation\n\nCall:\n")
        print(object@call)
        cat("\nCoefficients:\n")
        printCoefmat(object@summary)
        cat("\n-2 log L:", object@m2logL, "\n")
    }
)
