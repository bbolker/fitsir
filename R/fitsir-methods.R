##' @import methods
##' @import stats
##' @import graphics
NULL

##' Plot a fitsir object
##' @aliases plot,fitsir-method
##' @param x fitsir object
##' @param level the confidence level required
##' @param method confidence interval method
##' @param main an overall title for the plot
##' @param xlim the x limit of the plot
##' @param ylim the y limit of the plot
##' @param xlab a label for the x axis
##' @param ylab a label for the y axis
##' @param add (logical) add to an existing plot?
##' @param col.traj colour of the estimated trajectory
##' @param lty.traj line type of the estimated trajectory
##' @param col.conf colour of the confidence intervals
##' @param lty.conf line type of the confidence intervals
##' @importFrom bbmle plot
##' @examples
##' harbin2 <- setNames(harbin, c("times", "count"))
##' ss <- startfun(harbin2, type="death")
##' ff <- fitsir(harbin2, start=ss, type="death", method="BFGS")
##' plot(ff)
##' 
##' ff2 <- fitsir(harbin2, start=c(ss, ll.k=5), type="death", dist="nbinom", method="BFGS")
##' ff3 <- fitsir(harbin2, start=ss, type="death", dist="quasipoisson", method="BFGS")
##' plot(ff2, level=0.95, col.traj="red", main="Negative binomial error vs. Quasipoisson error CIs")
##' plot(ff3, add=TRUE, level=0.95, col.traj="blue", col.conf="blue")
##' legend(2, 270, legend = c("NB2", "Quasipoisson"), col=c("red", "blue"), lty=1)
##' @docType methods
##' @exportMethod plot
setMethod("plot", signature(x="fitsir", y="missing"),
    function(x, level,
             method=c("delta", "mvrnorm", "wmvrnorm"),
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

##' Extract parameter of a fit
##' @param object fitsir object
##' @param type type of parameter to be returned
##' @importFrom bbmle coef
##' @importFrom bbmle vcov
##' @details \code{raw} returns unconstrained parameters; \code{trans} returns constrained parameters; and
##' \code{summary} returns summarized parameters.
##' @examples 
##' coef(ff)
##' coef(ff,"trans")
##' coef(ff,"summary")
##' @docType methods
##' @exportMethod coef
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

##' Forecast from an SIR fit and find confidence interval
##' @param object fitsir object
##' @param level the confidence level required
##' @param times time vector to predict over. Default is set to the time frame of the data.
##' @param method confidence interval method. Default is set to Delta method.
##' @param debug print debugging output?
##' @details
##' See vignette for different methods: \code{vignette("details", package="fitsir")}
##' @importFrom bbmle predict
##' @importFrom bbmle confint
##' @importFrom MASS mvrnorm
##' @importFrom grDevices adjustcolor
##' @examples
##' predict(ff, level=0.95)
##' @docType methods
##' @exportMethod predict
setMethod("predict", "fitsir",
    function(object,
             level,times,
             method=c("delta", "mvrnorm", "wmvrnorm"),
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
                simpars <- mvrnorm(nsim,mu=coef(object),
                                       Sigma=vcov(object))
                
                for (i in 1:nsim) {
                    simtraj[,i] <- SIR.detsim(times, 
                                              trans.pars(simpars[i,]),
                                              type)
                }
            }
            
            cmat <- switch(method,
                delta={
                    nmat <- as.matrix(SIR.detsim(times, pars, type, grad=TRUE)[,-1])
                    nmat <- with(as.list(pars), {
                        logSens <- c(beta, gamma, N, i0*(1-i0))
                        sweep(nmat, 2, logSens, "*")
                    })
                    
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
                    traj.logLik <- -apply(simpars, 1, SIR.logLik, count=object@data$count, times, model=select_model(dist), type=type)
                    ##FIXME: vcov not symmetric for low tolerance?
                    i <- 10
                    while(!isSymmetric(round(vcov(object), i))) i <- i - 1
                    sample.logLik <- mvtnorm::dmvnorm(simpars, coef(object), round(vcov(object), i), log=TRUE)
                    ww <- exp(traj.logLik-sample.logLik)
                    cmat <- t(apply(simtraj, 1, wquant, weights=ww, probs=c(ll, 1-ll)))
                    cmat
                })
            
            cmat <- setNames(cmat, c(paste(100*ll, "%"), paste(100*(1-ll), "%")))
            
            if (debug & method != "delta") {
                matplot(times, simtraj, type="l",col=adjustcolor("black", alpha.f=0.1), lty=1)
                matlines(times, cmat, col=2, lty=1, lwd=2)
            }
            
            df <- cbind(df, cmat)
        }
        df
    }
)


##' Find residuals between the fit and the data
##' @param object fitsir object
##' @param type type of residuals. Default is set to pearson residuals.
##' @details 
##' \code{pearson} returns \eqn{(X_i - \mu_i)^2/\mu_i} and \code{raw} retrns \eqn{X_i -\mu_i}, 
##' where \eqn{X_i} is the observed counts and \eqn{\mu_i} is the expected coutns.
##' @importFrom bbmle residuals
##' @examples
##' residuals(ff)
##' @docType methods
##' @exportMethod residuals
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

setGeneric("dispersion", function(object, ...) standardGeneric("dispersion"))

##' Find dispersion parameter
##' @param object fitsir object
##' @param dist distribution
##' @details
##' \code{quasipoisson} returns the sum of pearson residuals divided by the degrees of freedom.
##' \code{nbinom} assumes quadratic mean-variance relation (var=mu+mu^2/k) and estimates k based on maximum likelihood.
##' \code{nbinom1} assumes linear mean-variance relation (var=(1+phi)mu) and estimates phi based on maximum likelihood.
##' @exportMethod dispersion
##' @docType methods
##' @exportMethod dispersion
setMethod("dispersion", "fitsir",
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

##' Summarize the fit
##' @param object fitsir object
##' @importFrom bbmle summary
##' @docType methods
##' @exportMethod summary
setMethod("summary", "fitsir",
    function(object){
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

##' Show summary of a fit
##' @param summary.fitsir object
##' @docType methods
##' @exportMethod show
setMethod("show", "summary.fitsir",
    function(object){
        cat("Maximum likelihood estimation\n\nCall:\n")
        print(object@call)
        cat("\nCoefficients:\n")
        printCoefmat(object@summary)
        cat("\n-2 log L:", object@m2logL, "\n")
    }
)
