##' Plotting method for fitsir objects
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
##' plot(ff, level=0.95, log="y")
##' ff2 <- fitsir(harbin2, type="death", method="BFGS")
##' plot(ff2, add=TRUE, col.traj="blue")
##' 
##' plot(ff2, level=0.8, col.traj="red", log="y")
##' ff3 <- fitsir(harbin2, type="death", dist="pois", method="BFGS")
##' plot(ff3, add=TRUE, level=0.8, col.traj="blue", col.conf="blue")
##' 
##' ff4 <- fitsir(harbin2, type="death", dist="nbinom")
##' plot(ff4, level=0.8, log="y")
setMethod("plot", signature(x="fitsir", y="missing"),
    function(x, level,
             main, xlim, ylim, xlab, ylab, add=FALSE,
             col.traj="black",lty.traj=1,
             col.conf="red",lty.conf=4,
             ...){
        times <- x@data$times
        count <- x@data$count
        pars <- coef(x)
        type <- x@data$type
        i.hat <- SIR.detsim(times, trans.pars(pars), type)
        
        if (missing(main)) main <- paste("fitsir result:", x@data$dist)
        if (missing(xlab)) xlab <- "time"
        if (missing(ylab)) ylab <- paste(type, "count")
        if (missing(ylim)) {
            ymin <- min(i.hat, count)
            ymax <- max(i.hat, count)
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
            cc <- confint(x, level=level)
            if (class(cc) == "mle2") {
                message("Plotting the better fit instead\n")
                bpars <- coef(cc)
                ## FIXME: I think this is an error with bbmle
                names(bpars) <- c("log.beta", "log.gamma", "log.N", "logit.i")
                bfit <- SIR.detsim(times, trans.pars(bpars), type)
                lines(times, bfit, col=col.conf, lty=lty.conf)
            } else {
                conf.traj <- apply(cc, 2, function(x) SIR.detsim(times, trans.pars(x), type))
                matlines(times, conf.traj, col=col.conf, lty=lty.conf)
            }
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

setMethod("summary", "fitsir",
    function(object,...){
        cc <- object@coef
        m <- unname(do.call(rbind, summarize.pars.deriv(cc)))
        cc.vcov <- m %*% object@vcov %*% t(m)
        smat <- rbind(
            Estimate=summarize.pars(cc),
            `Std. Error` = sqrt(diag(cc.vcov))
        )
        m2logL <- 2*object@min
        new("summary-fitsir", call=object@call.orig, coef=cc, summary=smat, m2logL=m2logL)
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
